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
    }

SDRF_FOR_FASTQS
    .map{ row-> tuple(row["${params.fields.run}"], row["${params.fields.fastq}"]) }
    .set { FASTQ_RUNS1 }

FASTQ_RUNS1
    .into { 
        PRINT_FASTQ_RUNS 
        FASTQ_RUNS
    }

REFERENCE_FASTA = Channel.fromPath( referenceFasta, checkIfExists: true )

// Call the download script to retrieve run fastqs

process download_fastqs {
    
    maxForks 12

    errorStrategy { task.attempt<=10 ? 'retry' : 'finish' } 
    
    input:
        set runId, runFastq from FASTQ_RUNS
    
    output:
        set val(runId), file("${runId}_1.fastq.gz") optional true into DOWNLOADED_FASTQS_R1
        set val(runId), file("${runId}_2.fastq.gz") optional true into DOWNLOADED_FASTQS_R2
        set val(runId), file("${runId}.fastq.gz") optional true into DOWNLOADED_FASTQS_UNPAIRED

    """
        fetchEnaFastqFtpUriViaSsh.sh ${runFastq}
    """
}

// Merge paired and unpaired FASTQS for processing

DOWNLOADED_FASTQS_UNPAIRED
    .concat(DOWNLOADED_FASTQS_R1)
    .concat(DOWNLOADED_FASTQS_R2)
    .set {
        DOWNLOADED_FASTQS
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
   
    errorStrategy { task.attempt<=10 ? 'retry' : 'finish' } 
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

    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' } 
    
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
    
    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' } 

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
    
    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' } 
    
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
    
    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' } 
    
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

    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 20

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

    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' } 
    
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

FILTERED_FASTQS.filter{ it.get(1).size()>40 }.into {
    FILTERED_FASTQS_COUNTS
    FILTERED_FASTQS_QUANT
    FILTERED_FASTQS_FASTQC
}

// Run fastQC on filtered reads

process filtered_fastqc {
   
    conda "${baseDir}/envs/fastqc.yml"
 
    errorStrategy { task.attempt<=10 ? 'retry' : 'finish' } 
    memory { 2.GB * task.attempt }

    publishDir "$resultsRoot/qc/fastqc/filtered", mode: 'copy', overwrite: true
    
    input:
        set val(runId), file(runFastq) from FILTERED_FASTQS_FASTQC
    output:
        set val(runId), file("${runFastq.simpleName}_fastqc.zip") into FILTERED_FASTQC

    """
        fastqc -t 8 --noextract ${runFastq}
    """
}

// Extract TSVs from the filtered fastq

process filtered_fastqc_tsv {
   
    conda 'irap-components'

    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' } 
    
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
    .join( 
        FILTERED_FASTQS_COUNTS
        .map{ tuple(it[1].simpleName, it[1]) }
    )
    .set {
        FASTQS_FOR_COUNTING_BY_FILENAME
    }


// Summarise counts at each filtering step

process count_reads {
    
    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' } 
    
    input:
        set val(fileName), file(runFastq), file(artFastq), file(contFastq), file(filtFastq) from FASTQS_FOR_COUNTING_BY_FILENAME    
    output:
        set val("${fileName}"), file("${fileName}.counts.tsv") into FASTQ_COUNTS
        
    """
        echo ${runFastq.simpleName},${runFastq.countFastq()},${artFastq.countFastq()},${contFastq.countFastq()},${filtFastq.countFastq()} > ${fileName}.counts.tsv
    """
}

// Collect the count lines and add a header

countsHeader = file('counts_header.txt')
countsHeader.text = "name,raw,artefacts_removed,contamination_filtered,uncalled_filtered"

Channel.value( countsHeader )
    .concat(FASTQ_COUNTS)
    .collectFile(name: 'fastq_counts.csv', storeDir: "$resultsRoot/qc/counts")
    .set{
        COLLECTED_FASTQ_COUNTS
    }

// Group read files by run name with strandedness

SDRF_FOR_STRAND
    .map{ row-> tuple(row["${params.fields.run}"], params.fields.containsKey('strand') && row.containsKey(params.fields.strand) ? row["${params.fields.strand}"] : 'not applicable') }
    .set {
        RUN_STRANDEDNESS
    }

FILTERED_FASTQS_QUANT
    .groupTuple(sort: true)
    .set{ GROUPED_FILTERED_FASTQS }

// Note: following the below these tuples will now be: <run id> <strandedness> <fastq files>

RUN_STRANDEDNESS.join( GROUPED_FILTERED_FASTQS ).set { GROUPED_FILTERED_FASTQS_WITH_STRAND }

// Group read files by run name, or by technical replicate group if specified

if ( params.fields.containsKey('techrep')){

    // If technical replicates are present, create a channel containing that info 

    SDRF_FOR_TECHREP
        .map{ row-> tuple(row["${params.fields.run}"], row["${params.fields.techrep}"]) }
        .groupTuple()
        .map{ row-> tuple( row[0], row[1][0]) }
        .set{ TECHREPS }

    // Now add the tech rep group to the run info, group by it, and create a
    // tuple of files keyed by techrep group

    TECHREPS.join( GROUPED_FILTERED_FASTQS_WITH_STRAND )
        .groupTuple(by: 1)
        .map{ row-> tuple( row[1], row[2][0], row[3].flatten()) }
        .set{
            TECHREP_GROUPED_FILTERED_FASTQS_WITH_STRAND
        }

    // Now collapse all of the FASTQs prior to quantification

    process merge_techrep_fastqs {

        input:
            set val(groupId), val(strand), file('*') from TECHREP_GROUPED_FILTERED_FASTQS_WITH_STRAND
        
        output:
            set val(groupId), val(strand), file("${groupId}*.fastq.gz") into FINAL_GROUPED_FASTQS

        """
            cat *_1.fastq.gz 2>/dev/null > ${groupId}_1.fastq.gz || :
            cat *_2.fastq.gz 2>/dev/null > ${groupId}_2.fastq.gz || :
            find . -name '*.fastq.gz' ! -name '*_1.fastq.gz' ! -name '*_2.fastq.gz' -exec cat {} \\; > ${groupId}.fastq.gz
            find . -empty -type f -delete 
        """
    } 
}else{
    GROUPED_FILTERED_FASTQS_WITH_STRAND.set{ FINAL_GROUPED_FASTQS }
}

// Separate paired and unpaired runs

PAIRED = Channel.create()
UNPAIRED = Channel.create()

FINAL_GROUPED_FASTQS.choice( UNPAIRED, PAIRED ) {a -> 
    a[2].size() == 2 ? 1 : 0
}

// Synchronise paired-end read files 

process synchronise_pairs {
  
    conda "${baseDir}/envs/fastq_utils.yml"
  
    memory { 5.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 3
    
    input:
        set val(runId), val(strand), file('*') from PAIRED

    output:
        set val(runId), val(strand), file( "matched/${runId}_1.fastq.gz" ), file("matched/${runId}_2.fastq.gz") into MATCHED_PAIRED_FASTQS
        set val(runId), val(strand), file( "unmatched/${runId}.fastq.gz" ) into UNMATCHED_PAIRED_FASTQS

    beforeScript 'mkdir -p matched && mkdir -p unmatched'

    """
        fastq_filterpair ${runId}_1.fastq.gz ${runId}_2.fastq.gz matched/${runId}_1.fastq.gz matched/${runId}_2.fastq.gz unmatched/${runId}.fastq.gz sorted 
    """          
}

// Generate a Kallisto index 

process kallisto_index {

    conda "${baseDir}/envs/kallisto.yml"
    
    memory { 5.GB * task.attempt }
    time { 3.hour * task.attempt }
    cpus 8

    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 20
    
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
    
    publishDir "$resultsRoot/kallisto", mode: 'move', overwrite: true
    
    memory { 4.GB * task.attempt }
    time { 3.hour * task.attempt }
    cpus 8
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 20

    input:
        file(kallistoIndex) from KALLISTO_INDEX_SINGLE.first()
        set val(runId), val(strand), file(runFastq) from UNPAIRED

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
    
    publishDir "$resultsRoot/kallisto", mode: 'move', overwrite: true
    
    memory { 4.GB * task.attempt }
    time { 3.hour * task.attempt }
    cpus 8
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 20

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

