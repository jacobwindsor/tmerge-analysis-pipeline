nextflow.preview.dsl=2

tmerge2_path = "/users/rg/jwindsor/tmerge/2.0"
sirvs_path = "/users/rg/jwindsor/tests/tmerge/runner/data/SIRVome.gff"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"

sensitivy_threshold = 96.4
precision_threshold = 39.4

include gffToBED from './utils'

process runTmerge2 {
    input:
    path input

    output:
        path 'output.gff'

    shell:
    '''
    module load Python
    source !{venv}

    python !{tmerge2_path}/tmerge.py -i !{input} -o output.gff -t 2 -f 4
    '''
}

process runGFFCompare {
    input:
    path input

    output:
    stdout()

    shell:
    '''
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -o !{input.baseName}.gffcompare !{input} -r !{sirvs_path} -T
    cat output.gffcompare
    '''
}

process checkSensitivityPrecision{
    input:
    val input

    executor "local"

    exec:
    (transcript_level_row, precision, sensitivity) = (input =~ /.*Transcript level:\s*(\d*\.\d*)\s*\|\s*(\d*.\d*)/)[0]

    if ((Float.parseFloat(precision) < precision_threshold) || (Float.parseFloat(sensitivity) < sensitivy_threshold)) {
        throw new Exception("Sensitivity and/or precision are not above the required threshold. Sensitivity: $sensitivity, Precision: $precision")
    }
}

workflow.onComplete {
    if (workflow.success) {
        println "Pipeline completed in: $workflow.duration"
    }
    else {
        println "Oops!"
    }
}

workflow {
    Channel.fromPath("/users/rg/jwindsor/tests/tmerge/runner/data/sirv_reads.gff") | runTmerge2 | runGFFCompare | checkSensitivityPrecision
}