import java.util.regex.Pattern

class WorkflowParamValidator {
    private static final Pattern TOKEN = Pattern.compile(/[A-Za-z0-9][A-Za-z0-9._+-]*/)
    private static final Pattern INTEGER = Pattern.compile(/[0-9]+/)
    private static final Pattern DECIMAL = Pattern.compile(/[0-9]+(\.[0-9]+)?/)
    private static final Pattern PATH_VALUE = Pattern.compile(/[A-Za-z0-9._+,:=@%\/-]+/)
    private static final Pattern HEADER_VALUE = Pattern.compile(/[^\p{Cntrl}$`"';|&<>]+/)

    static void validate(def params) {
        requirePath(params, 'sdrf')
        requirePath(params, 'resultsRoot')
        requirePath(params, 'transcriptomeIndex')
        requireToken(params, 'downloadMethod')
        requireInteger(params, 'maxConcurrentDownloads')

        optionalPath(params, 'manualDownloadFolder')
        optionalPath(params, 'fastqProviderConfig')
        optionalPath(params, 'contaminationIndex')

        requireFields(params.fields, ['run', 'fastq', 'layout'])
        optionalFields(params.fields, ['quality', 'controlled_access', 'strand', 'techrep'])

        requireNestedIntegerMap(params.fastq_quality_filter, 'params.fastq_quality_filter', ['Q', 'p', 'q'])
        requireNestedIntegerMap(params.fastq_quality_trimmer, 'params.fastq_quality_trimmer', ['Q', 't', 'l'])
        requireNestedIntegerMap(params.fastq_trim_poly_at, 'params.fastq_trim_poly_at', ['min_len', 'min_poly_at_len'])
        requireNestedIntegerMap(params.fastq_filter_n, 'params.fastq_filter_n', ['n'])
        requireNestedDecimal(params.kallisto?.quant?.se, 'params.kallisto.quant.se', 'l')
        requireNestedDecimal(params.kallisto?.quant?.se, 'params.kallisto.quant.se', 's')
    }

    static String safeToken(value, String fieldName) {
        def text = value == null ? '' : value.toString()
        if (!(TOKEN.matcher(text).matches())) {
            throw new IllegalArgumentException("Unsafe SDRF value for ${fieldName}: '${text}'")
        }
        text
    }

    static String safeUri(value, String fieldName) {
        def text = value == null ? '' : value.toString()
        if (!(text ==~ /[^\p{Cntrl}\s]+/)) {
            throw new IllegalArgumentException("Unsafe SDRF URI value for ${fieldName}: '${text}'")
        }
        text
    }

    static String safeControlledAccess(value) {
        def text = value == null ? 'no' : value.toString().toLowerCase()
        if (!(text in ['yes', 'no'])) {
            throw new IllegalArgumentException("Unsafe SDRF controlled access value: '${value}'")
        }
        text
    }

    static String safeLayout(value, String fieldName) {
        def text = value == null ? '' : value.toString()
        if (!(text in ['SINGLE', 'PAIRED'])) {
            throw new IllegalArgumentException("Unsafe SDRF layout value for ${fieldName}: '${text}'")
        }
        text
    }

    static String shellQuote(value) {
        "'" + value.toString().replace("'", "'\"'\"'") + "'"
    }

    private static void requireFields(def fields, List required) {
        if (fields == null) {
            throw new IllegalArgumentException('Missing workflow params.fields settings')
        }
        required.each { requireNestedHeader(fields, 'params.fields', it) }
    }

    private static void optionalFields(def fields, List optional) {
        optional.each { optionalNestedHeader(fields, 'params.fields', it) }
    }

    private static void requireNestedIntegerMap(def params, String scope, List keys) {
        if (params == null) {
            throw new IllegalArgumentException("Missing workflow ${scope} settings")
        }
        keys.each { requireNestedInteger(params, scope, it) }
    }

    private static void requirePath(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", PATH_VALUE)
    }

    private static void optionalPath(def params, String name) {
        if (has(params, name) && params.get(name) != null && params.get(name).toString() != '') {
            assertPattern(params.get(name), "params.${name}", PATH_VALUE)
        }
    }

    private static void requireToken(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", TOKEN)
    }

    private static void requireInteger(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", INTEGER)
    }

    private static void requireNestedInteger(def params, String scope, String name) {
        requireValue(params, "${scope}.${name}", name)
        assertPattern(params.get(name), "${scope}.${name}", INTEGER)
    }

    private static void requireNestedDecimal(def params, String scope, String name) {
        requireValue(params, "${scope}.${name}", name)
        assertPattern(params.get(name), "${scope}.${name}", DECIMAL)
    }

    private static void requireNestedHeader(def params, String scope, String name) {
        requireValue(params, "${scope}.${name}", name)
        assertPattern(params.get(name), "${scope}.${name}", HEADER_VALUE)
    }

    private static void optionalNestedHeader(def params, String scope, String name) {
        if (has(params, name) && params.get(name) != null && params.get(name).toString() != '') {
            assertPattern(params.get(name), "${scope}.${name}", HEADER_VALUE)
        }
    }

    private static void requireValue(def params, String label, String key) {
        if (params == null || !has(params, key) || params.get(key) == null || params.get(key).toString() == '') {
            throw new IllegalArgumentException("Missing required workflow parameter ${label}")
        }
    }

    private static void assertPattern(value, String label, Pattern pattern) {
        def text = value == null ? '' : value.toString()
        if (!pattern.matcher(text).matches()) {
            throw new IllegalArgumentException("Invalid workflow parameter ${label}: '${text}'")
        }
    }

    private static boolean has(def params, String key) {
        params != null && params.containsKey(key)
    }
}
