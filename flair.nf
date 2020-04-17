nextflow.preview.dsl=2

/*
* Configuration options. 
*/
// Options that are only used in this context
tmerge2_path = "/users/rg/jwindsor/tmerge/2.0"
flair_path = "/users/rg/jwindsor/flair"
read_mapping_path = "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/readMapping/"
sirvs_path = "/users/rg/jlagarde/genomes/lexogen_SIRVs/SIRV_Set1_Lot00141_Sequences_181206a/SIRVome_isoforms_Lot00141_C_181206a.gtf.unix.corrected.gene_types.gtf"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"

// Options used by included files and this context
params.output_dir = "/users/rg/jwindsor/tests/tmerge/results/flair"
params.julien_utils_path = "/users/rg/jlagarde/julien_utils/"

include processForSIRVs as processTmerge2SIRVs from './utils'
include processForSIRVs as processFLAIRSIRVs from './utils'
include gffToBED from './utils'

inputFiles = Channel.fromList([
    [
        "standard",
        "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/pacBio:Cshl:Smarter:Corr0_HpreCap_0+_Brain01Rep1.strandedHCGMs.gff.gz"
    ]
])

no_redundant_choices = ["none", "longest", "best_only"]
support_choices = [0, 1, 2, 5, 8, 10]


process copyInputGFF {
    input:
    val x

    publishDir "$params.output_dir"

    output:
    path '*.input.gff'

    shell:
    '''
    zcat "!{x[1]}" > !{x[0]}.input.gff
    '''    
}

process bedToGFF {
    input:
    path inputBED

    publishDir "$params.output_dir"
    errorStrategy "ignore"

    output:
    path '*.gff'

    shell:
    '''
    PATH="$PATH:!{params.julien_utils_path}"

    cat !{inputBED} | sortgff | awk -f !{params.julien_utils_path}bed12fields2gff.awk > !{inputBED.baseName}.gff
    '''
}


process runTmerge2 {
    input:
    path input

    memory '30 GB'
    publishDir "$params.output_dir"

    output:
        path '*.output.tmerge2.gff'

    shell:
    '''
    module load Python
    source !{venv}

    python !{tmerge2_path}/tmerge.py -i !{input} -o !{input.simpleName}.output.tmerge2.gff
    '''
}

process runFLAIR {
    input:
    val input

    memory '30GB'
    publishDir "$params.output_dir"

    output:
        path '*.output.flair.*.bed'

    shell:
        '''
        python !{flair_path}/bin/collapse_isoforms_precise.py -q !{input[0]} -s !{input[2]} -n !{input[1]} -o !{input[0].simpleName}.output.flair.support_!{input[2]}.tolerance_!{input[1]}.bed
        '''
}


process runGFFCompare {
    input:
    path input

    publishDir "$params.output_dir"

    output:
    path '*.gffcompare*'

    shell:
    '''
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -o !{input.baseName}.gffcompare !{input} -r !{sirvs_path} -T
    '''
}

workflow runTools {
    take: inputs
    main:
        // TMERGE 1 and 2
        tmerge2 = runTmerge2(inputs)
        tmerge2_sirvs = tmerge2 | processTmerge2SIRVs

        // FLAIR
        flair = gffToBED(inputs) | combine(no_redundant_choices) | combine(support_choices) | runFLAIR
        flairGFF = bedToGFF(flair) | processFLAIRSIRVs | mix(tmerge2_sirvs) | runGFFCompare

    emit:
        tmerge2.mix(flair)
}

workflow convertToBed {
    take: inputs
    main:
        gffToBED(inputs)
}

workflow {
    outputs = copyInputGFF(inputFiles) | runTools
    outputs.filter( ~/.*\.gff$/ ) | convertToBed
}