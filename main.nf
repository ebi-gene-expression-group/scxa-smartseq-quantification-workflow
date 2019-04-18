#!/usr/bin/env nextflow

sdrfFile = params.sdrf
resultsRoot = params.resultsRoot
referenceFasta = params.referenceFasta
contaminationIndex = params.contaminationIndex

// Read ENA_RUN column from an SDRF

Channel
    .fromPath(sdrfFile, checkIfExists: true)
    .splitCsv(header:true, sep:"\t")
    .filter{ row -> (! row.containsKey(params.fields.quality)) || ( row["${params.fields.quality}"].toLowerCase() != 'not ok') }
    .into {
        SDRF_FOR_FASTQS
        SDRF_FOR_STRAND
        SDRF_FOR_TECHREP
        SDRF_FOR_COUNT
        SDRF_FOR_TECHREP_COUNT
    }

SDRF_FOR_FASTQS
    .map{ row-> tuple(row["${params.fields.run}"], row["${params.fields.fastq}"]) }
    .set { FASTQ_RUNS }

REFERENCE_FASTA = Channel.fromPath( referenceFasta, checkIfExists: true )

SDRF_FOR_COUNT
    .map{ row-> tuple(row["${params.fields.run}"]) }
    .unique()
    .count()
    .set { RUN_COUNT }

// Get the file names from the URLs

process get_download_filename {
    
    executor 'local'

    input:
        set runId, runURI from FASTQ_RUNS
    
    output:
        set val(runId), val(runURI), stdout into FASTQ_RUNS_FILES

    """
        basename $runURI | tr -d \'\\n\'
    """
}

// Call the download script to retrieve run fastqs

process download_fastqs {
    
    conda "${baseDir}/envs/atlas-fastq-provider.yml"
    
    maxForks params.maxConcurrentDownloads
    time { 1.hour * task.attempt }

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }   
    maxRetries 3
 
    input:
        set runId, runURI, runFastq from FASTQ_RUNS_FILES

    output:
        set val(runId), file("${runFastq}") into DOWNLOADED_FASTQS

    """
        confPart=''
        if [ -e "$NXF_TEMP/atlas-fastq-provider/download_config.sh" ]; then
            confPart=" -c $NXF_TEMP/atlas-fastq-provider/download_config.sh"
        fi 
        fetchFastq.sh -f ${runURI} -t ${runFastq} -m ${params.downloadMethod} \$confPart
    """
}

// Copy channels to allow different operations

DOWNLOADED_FASTQS.into {
    DOWNLOADED_FASTQS_COUNTS
    DOWNLOADED_FASTQS_FILTERING
    DOWNLOADED_FASTQS_FASTQC
}

// Run fastQC on raw reads

process raw_fastqc {
   
    conda "${baseDir}/envs/fastqc.yml"
   
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }   
    maxRetries 3
    
    memory { 2.GB * task.attempt }
 
    publishDir "$resultsRoot/qc/fastqc/raw", mode: 'copy', overwrite: true
    
    input:
        set val(runId), file(runFastq) from DOWNLOADED_FASTQS_FASTQC
    output:
        set val(runId), file("${runFastq.simpleName}_fastqc.zip") into RAW_FASTQC

    """
        fastqc -t 8 --noextract ${runFastq}
    """
}

// Filter based on qualities

process quality_filter {
    
    conda "${baseDir}/envs/fastx_toolkit.yml"

    errorStrategy {  task.attempt < 3 ? 'retry' : 'ignore' }   
    maxRetries 3
    
    input:
        set val(runId), file(runFastq) from DOWNLOADED_FASTQS_FILTERING
    output:
        set val(runId), file("qfilt/${runFastq}") into QFILT_FASTQS

    beforeScript 'mkdir -p qfilt'

    """
        zcat ${runFastq} | fastq_quality_filter -o qfilt/${runFastq} -v \
            -Q ${params.fastq_quality_filter.Q} -p ${params.fastq_quality_filter.p} \
            -q ${params.fastq_quality_filter.q}
   """
}

// Trim based on qualities

process quality_trim {
    
    conda "${baseDir}/envs/fastx_toolkit.yml"
    
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }   
    maxRetries 3

    input:
        set val(runId), file(runFastq) from QFILT_FASTQS
    output:
        set val(runId), file("qtrim/${runFastq}") into QTRIM_FASTQS

    beforeScript 'mkdir -p qtrim'
    
    """
        cat ${runFastq} | fastq_quality_trimmer -v -Q ${params.fastq_quality_trimmer.Q} \
            -t ${params.fastq_quality_trimmer.t}  -l ${params.fastq_quality_trimmer.l} \
            -o qtrim/${runFastq}
    """
}

// Trim poly-ats (Nuno's fastq_utils) 

process quality_polya {

    conda "${baseDir}/envs/fastq_utils.yml"
    
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }   
    maxRetries 3
    
    input:
        set val(runId), file(runFastq) from QTRIM_FASTQS
    output:
        set val(runId), file("polyatrim/${runFastq}") into POLYATRIM_FASTQS

    beforeScript 'mkdir -p polyatrim'
    
    """
        cat ${runFastq} | fastq_trim_poly_at --min_len ${params.fastq_trim_poly_at.min_len} \
            --min_poly_at_len ${params.fastq_trim_poly_at.min_poly_at_len} \
            --file - --outfile polyatrim/${runFastq}
    """
}

// Run the fastx_artifacts filter

process quality_artifacts {
    
    conda "${baseDir}/envs/fastx_toolkit.yml"
    
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }   
    maxRetries 3
    
    input:
        set val(runId), file(runFastq) from POLYATRIM_FASTQS
    output:
        set val(runId), file("qualart/${runFastq}") into ARTEFACTS_FASTQS

    beforeScript 'mkdir -p qualart'
    
    """
        zcat ${runFastq} | fastx_artifacts_filter -v -Q 33 | \
            gzip -c - > "qualart/${runFastq}"
    """
}

ARTEFACTS_FASTQS.into {
    ARTEFACTS_FASTQS_CONT
    ARTEFACTS_FASTQS_COUNTS
}

// Do alignments agains a contamination index and filter out matches

process quality_contamination {
    
    conda 'bowtie2 samtools r-data.table'
    
    memory { 5.GB * task.attempt }
    time { 3.hour * task.attempt }
    cpus 8

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'ignore' }
    maxRetries 3

    input:
        set val(runId), file(runFastq) from ARTEFACTS_FASTQS_CONT
    output:
        set val(runId), file("contfilt/${runFastq}") into CONT_FASTQS
        set val("${runFastq.simpleName}"), file("contfilt/${runId}.cont.bam") into CONT_BAMS
    
    beforeScript 'mkdir -p contfilt'

    """
        bowtie2 -p 8 --very-fast --un contfilt/${runFastq} --fast-local --phred33 \
            -x $contaminationIndex -U ${runFastq} -S /dev/stdout | \
            samtools view -S -b -F 4 - > contfilt/${runId}.cont.bam
    """
}

CONT_FASTQS.into {
    CONT_FASTQS_UNCALLED
    CONT_FASTQS_COUNTS
    CONT_FASTQS_CHAR
}

// Filter uncalled bases

process quality_uncalled {
    
    conda "${baseDir}/envs/fastq_utils.yml"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }   
    maxRetries 3

    input:
        set val(runId), file(runFastq) from CONT_FASTQS_UNCALLED
    output:
        set val(runId), file("uncalled/${runFastq}") into FILTERED_FASTQS

    beforeScript 'mkdir -p uncalled'
    
    """
        fastq_filter_n -n ${params.fastq_filter_n.n} ${runFastq} | gzip -c - > uncalled/${runFastq}
    """
}

// Filter out near-0 sizes files. Zipped empty files seem to have a size of 20.

FILTERED_FASTQS
    .into{
        FILTERED_FASTQS_FOR_DOWNSTREAM
        FILTERED_FASTQS_FOR_COUNT
    }

FILTERED_FASTQS_FOR_COUNT
    .map{ tuple(it[0]) }
    .unique()
    .count()
    .set { FILTERED_FASTQ_RUN_COUNT }

FILTERED_FASTQS_FOR_DOWNSTREAM.filter{ it.get(1).size()>40 }.into {
    FILTERED_NONEMPTY_FASTQS_COUNTS
    FILTERED_NONEMPTY_FASTQS_QUANT
    FILTERED_NONEMPTY_FASTQS_FASTQC
    FILTERED_NONEMPTY_FASTQS_FOR_COUNT
}

FILTERED_NONEMPTY_FASTQS_FOR_COUNT
    .map{ tuple(it[0]) }
    .unique()
    .count()
    .set { FILTERED_NONEMPTY_FASTQ_RUN_COUNT }

// Run fastQC on filtered reads

process filtered_fastqc {
   
    conda "${baseDir}/envs/fastqc.yml"
 
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }   
    maxRetries 3
    memory { 2.GB * task.attempt }

    publishDir "$resultsRoot/qc/fastqc/filtered", mode: 'copy', overwrite: true
    
    input:
        set val(runId), file(runFastq) from FILTERED_NONEMPTY_FASTQS_FASTQC
    output:
        set val(runId), file("${runFastq.simpleName}_fastqc.zip") into FILTERED_FASTQC

    """
        fastqc -t 8 --noextract ${runFastq}
    """
}

// Extract TSVs from the filtered fastq

process filtered_fastqc_tsv {
   
    conda 'irap-components'

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }   
    maxRetries 3
    
    publishDir "$resultsRoot/qc/fastqc/filtered", mode: 'copy', overwrite: true
    
    input:
        set val(runId), file(runFilteredFastqc) from FILTERED_FASTQC
    output:
        set val(runId), file("${runId}.filtered_fastqc.tsv") into FILTERED_FASTQC_TSV

    """
        irap_fastqc2tsv ${runFilteredFastqc} > ${runId}.filtered_fastqc.tsv
    """
}

// Generate a summary of counts at different filtering steps
// First collect the related outputs re-keyed by file name 


DOWNLOADED_FASTQS_COUNTS
    .map{ tuple(it[1].simpleName, it[1]) }
    .join( 
        ARTEFACTS_FASTQS_COUNTS
        .map{ tuple(it[1].simpleName, it[1]) }
    )
    .join( 
        CONT_FASTQS_COUNTS
        .map{ tuple(it[1].simpleName, it[1]) }
    )
    .join( FILTERED_NONEMPTY_FASTQS_COUNTS.map{ tuple(it[1].simpleName, it[1]) }, remainder: true)
    .set {
        FASTQS_FOR_COUNTING_BY_FILENAME
    }


// Summarise counts at each filtering step


process count_reads {
    
    conda 'irap-components'

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }   
    maxRetries 3
    
    input:
        set val(fileName), file(runFastq), file('art.fastq.gz'), file('cont.fastq.gz'), file('filt.fastq.gz') from FASTQS_FOR_COUNTING_BY_FILENAME    
    output:
        stdout  FASTQ_COUNTS
        
    """
        echo ${runFastq.simpleName},`num_reads.sh ${runFastq}`,`num_reads.sh art.fastq.gz`,`num_reads.sh cont.fastq.gz`,`num_reads.sh filt.fastq.gz`
    """
}

// Collect the count lines and add a header

FASTQ_COUNTS
    .collectFile(name: 'fastq_counts_nohead.csv', sort: true)
    .set{
        MERGED_FASTQ_COUNTS
    }

process head_counts {

    publishDir "$resultsRoot/qc/counts", mode: 'copy', overwrite: true

    input:
       file(mergedCounts) from MERGED_FASTQ_COUNTS

    output:
       file('fastq_counts.csv') into FINAL_COUNTS

    """
        echo "name,raw,artefacts_removed,contamination_filtered,uncalled_filtered" > fastq_counts.csv
        cat $mergedCounts >> fastq_counts.csv
    """ 
}


// Group read files by run name with strandedness

SDRF_FOR_STRAND
    .map{ row-> tuple(row["${params.fields.run}"], params.fields.containsKey('strand') && row.containsKey(params.fields.strand) ? row["${params.fields.strand}"] : 'not applicable', row["${params.fields.layout}"]) }
    .set {
        RUN_META
    }

FILTERED_NONEMPTY_FASTQS_QUANT
    .groupTuple(sort: true)
    .set{ GROUPED_FILTERED_FASTQS }

// Note: following the below these tuples will now be: <run id> <strandedness> <layout> <fastq files>

RUN_META.join( GROUPED_FILTERED_FASTQS ).set { GROUPED_FILTERED_FASTQS_WITH_META }

// Group read files by run name, or by technical replicate group if specified

if ( params.fields.containsKey('techrep')){

    // If technical replicates are present, create a channel containing that info 

    SDRF_FOR_TECHREP
        .map{ row-> tuple(row["${params.fields.run}"], row["${params.fields.techrep}"]) }
        .groupTuple()
        .map{ row-> tuple( row[0], row[1][0]) }
        .set{ TECHREPS }

    // The target set of results will now be the technical replicate group number

    SDRF_FOR_COUNT
        .map{ row-> tuple(row["${params.fields.techrep}"]) }
        .unique()
        .count()
        .into {
            TARGET_COUNT 
        }
    
    // Now add the tech rep group to the run info, group by it, and create a
    // tuple of files keyed by techrep group

    TECHREPS.join( GROUPED_FILTERED_FASTQS_WITH_META )
        .groupTuple(by: 1)
        .map{ row-> tuple( row[1], row[2][0], row[3][0], row[4].flatten()) }
        .set{
            TECHREP_GROUPED_FILTERED_FASTQS_WITH_META
        }

    // Now collapse all of the FASTQs prior to quantification

    process merge_techrep_fastqs {

        input:
            set val(groupId), val(strand), val(layout), file('*') from TECHREP_GROUPED_FILTERED_FASTQS_WITH_META
        
        output:
            set val(groupId), val(strand), val(layout), file("${groupId}*.fastq.gz") into FINAL_GROUPED_FASTQS

        """
            cat *_1.fastq.gz 2>/dev/null > ${groupId}_1.fastq.gz || :
            cat *_2.fastq.gz 2>/dev/null > ${groupId}_2.fastq.gz || :
            find . -name '*.fastq.gz' ! -name '*_1.fastq.gz' ! -name '*_2.fastq.gz' -exec cat {} \\; > ${groupId}.fastq.gz
            find . -empty -type f -delete 
        """
    } 
}else{
    GROUPED_FILTERED_FASTQS_WITH_META.set{ FINAL_GROUPED_FASTQS }

    RUN_COUNT
        .set{
            TARGET_COUNT
        }
}

// Validate the number of read files by layout

process validate_layout {

    errorStrategy 'ignore'

    input:
       set val(runId), val(strand), val(layout), file('*') from FINAL_GROUPED_FASTQS 

    output:
       set val(runId), val(strand), val(layout), file('validated/*.fastq.gz') into FINAL_VALIDATED_GROUPED_FASTQS 

    """
        if [ "$layout" == 'PAIRED' ]; then
            if [ ! -e "${runId}_1.fastq.gz" ] || [ ! -e "${runId}_2.fastq.gz" ]; then
                echo "One or more paired end read files not found for ${runId}" 1>&2
                exit 1
            fi
        elif [ ! -e "${runId}.fastq.gz" ]; then
            echo "Single-end read file not found for ${runId}"        
            exit 1
        fi  
        
        mkdir -p validated
        cp -P *.fastq.gz validated
    """
}

// Separate paired and unpaired runs

PAIRED = Channel.create()
UNPAIRED = Channel.create()

FINAL_VALIDATED_GROUPED_FASTQS.choice( UNPAIRED, PAIRED ) {a -> 
    a[2] == 'PAIRED' ? 1 : 0
}

// Synchronise paired-end read files 

process synchronise_pairs {
  
    conda "${baseDir}/envs/fastq_pair.yml"
  
    memory { 5.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 3
    
    input:
        set val(runId), val(strand), val(layout), file('*') from PAIRED

    output:
        set val(runId), val(strand), file( "matched/${runId}_1.fastq.gz" ), file("matched/${runId}_2.fastq.gz") into MATCHED_PAIRED_FASTQS
        set val(runId), val(strand), file( "unmatched/${runId}_1.fastq.gz" ), file( "unmatched/${runId}_2.fastq.gz" ) into UNMATCHED_PAIRED_FASTQS

    beforeScript 'mkdir -p matched && mkdir -p unmatched'

    """
        zcat ${runId}_1.fastq.gz > ${runId}_1.fastq
        zcat ${runId}_2.fastq.gz > ${runId}_2.fastq

        fastq_pair ${runId}_1.fastq ${runId}_2.fastq

        gzip ${runId}_1.fastq.single.fq && mv ${runId}_1.fastq.single.fq.gz unmatched/${runId}_1.fastq.gz
        gzip ${runId}_2.fastq.single.fq && mv ${runId}_2.fastq.single.fq.gz unmatched/${runId}_2.fastq.gz
        gzip ${runId}_1.fastq.paired.fq && mv ${runId}_1.fastq.paired.fq.gz matched/${runId}_1.fastq.gz
        gzip ${runId}_2.fastq.paired.fq && mv ${runId}_2.fastq.paired.fq.gz matched/${runId}_2.fastq.gz

        rm -f ${runId}_1.fastq ${runId}_2.fastq
    """          
}

// Generate a Kallisto index 

process kallisto_index {

    conda "${baseDir}/envs/kallisto.yml"

    // Cache based on path and file size
    cache 'deep'
    
    memory { 5.GB * task.attempt }
    time { 3.hour * task.attempt }
    cpus 8

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 5 ? 'retry' : 'ignore' }
    maxRetries 5
    
    input:
         file(referenceFastq) from REFERENCE_FASTA

    output:
        file ("${referenceFastq}.index") into KALLISTO_INDEX
        
    """
        kallisto index --kmer-size ${params.kallisto.index.kmer_size} \
            -i ${referenceFastq}.index ${referenceFastq}
    """
}

// Re-use index for both single- and paired-end reads

KALLISTO_INDEX.into{
    KALLISTO_INDEX_SINGLE
    KALLISTO_INDEX_PAIRED
}

// Run Kallisto for each single- read file
process kallisto_single {

    conda "${baseDir}/envs/kallisto.yml"
    
    publishDir "$resultsRoot/kallisto", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    time { 3.hour * task.attempt }
    cpus 8
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 5 ? 'retry' : 'ignore' }
    maxRetries 5

    input:
        file(kallistoIndex) from KALLISTO_INDEX_SINGLE.first()
        set val(runId), val(strand), val(layout), file(runFastq) from UNPAIRED

    output:
        file "${runId}" into KALLISTO_SINGLE

    script:
        
        def strandedness = ''

        if ( strand == 'first strand' ){
            strandedness = '--fr-stranded'
        }else if ( strand == 'second strand' ){
            strandedness = '--rf-stranded'
        }

        """
            kallisto quant $strandedness -i ${kallistoIndex} --single \
                -l ${params.kallisto.quant.se.l} -s ${params.kallisto.quant.se.s} \
                -t ${task.cpus} -o ${runId} ${runFastq}          
        """
}

// Run Kallisto for each synchronised paired-end read set

process kallisto_paired {

    conda "${baseDir}/envs/kallisto.yml"
    
    publishDir "$resultsRoot/kallisto", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    time { 3.hour * task.attempt }
    cpus 8
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 5? 'retry' : 'ignore' }
    maxRetries 5

    input:
        file(kallistoIndex) from KALLISTO_INDEX_PAIRED.first()
        set val(runId), val(strand), file(read1), file(read2) from MATCHED_PAIRED_FASTQS

    output:
        file "${runId}" into KALLISTO_PAIRED

    script:
        
        def strandedness = ''

        if ( strand == 'first strand' ){
            strandedness = '--fr-stranded'
        }else if ( strand == 'second strand' ){
            strandedness = '--rf-stranded'
        }

        """
            kallisto quant ${strandedness} -i ${kallistoIndex} -t ${task.cpus} -o ${runId} ${read1} ${read2}          
        """
}

// Check the total number of runs we have 

KALLISTO_SINGLE
    .concat(KALLISTO_PAIRED)
    .count()
    .set{ KALLISTO_RESULTS_COUNT } 

process validate_results {
    
    executor 'local'
    
    errorStrategy 'finish'  
    
    input:
        val(kallistoResultCount) from KALLISTO_RESULTS_COUNT 
        val(runCount) from RUN_COUNT
        val(targetCount) from TARGET_RESULT_COUNT
        val(filteredFastqRunCount) from FILTERED_FASTQ_RUN_COUNT
        val(filteredNonemptyFastqRunCount) from FILTERED_NONEMPTY_FASTQ_RUN_COUNT
        file(finalCounts) from FINAL_COUNTS

    output:
        stdout DONE

    """
    if [ "$kallistoResultCount" -ne "$targetCount" ]; then
        
        echo "Kallisto results count of $kallistoResultCount does not match expected results number ($targetCount)" 1>&2

        if [ "$filteredFastqRunCount" -ne "$runCount" ]; then
            echo "... filtering failed" 1>&2
            exit 1
        elif [ "$filteredNonemptyFastqRunCount" -ne "$filteredFastqRunCount" ]; then
            echo "... this is just because filtering removed some runs" 1>&2
            exit 0
        else
            echo "... unrecoverable errors may have occurred at quantification" 1>&2
            exit 1 
        fi
    else
        echo "Kallisto results count of $kallistoResultCount matches expected results number ($targetCount)"
    fi
    """
}   


