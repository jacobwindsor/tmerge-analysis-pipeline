nextflow.preview.dsl=2


include processForSIRVs from './utils'
include gffToBED from './utils'

tmerge1_path = "/users/rg/jlagarde/julien_utils"
tmerge2_path = "/users/rg/jwindsor/tmerge/2.0"
sirvs_path = "/users/rg/jlagarde/genomes/lexogen_SIRVs/SIRV_Set1_Lot00141_Sequences_181206a/SIRVome_isoforms_Lot00141_C_181206a.gtf.unix.corrected.gene_types.gtf"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"

params.julien_utils_path = "/users/rg/jlagarde/julien_utils/"
params.output_dir = "/users/rg/jwindsor/tests/tmerge/results/tolerances"


inputFiles = Channel.from("/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/pacBio:Cshl:Smarter:Corr0_HpreCap_0+_Brain01Rep1.strandedHCGMs.gff.gz")
tolerances = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
process copyInputFiles {
    input:
    val x

    publishDir "$params.output_dir"

    output:
    path 'input.gff'

    shell:
    '''
    zcat "!{x}" > input.gff
    '''    
}

process runTmerge2 {
    input:
    val input_tolerance_pair

    publishDir "$params.output_dir"

    output:
        path 'output.tmerge2.tolerance_*.gff'

    shell:
    '''
    module load Python
    source !{venv}

    python !{tmerge2_path}/tmerge.py -i !{input_tolerance_pair[0]} -o output.tmerge2.tolerance_!{input_tolerance_pair[1]}.gff -t !{input_tolerance_pair[1]}
    '''
}

process runTmerge1 {
    input:
    val input_tolerance_pair

    publishDir "$params.output_dir"

    output:
        path 'output.tmerge1.tolerance_*.gff'

    shell:
    '''
    PATH="$PATH:!{params.julien_utils_path}"

    perl !{tmerge1_path}/tmerge --exonOverhangTolerance=!{input_tolerance_pair[1]} !{input_tolerance_pair[0]} | sortgff > output.tmerge1.tolerance_!{input_tolerance_pair[1]}.gff
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
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -o !{input.baseName}.gffcompare !{input} -r !{sirvs_path}
    '''
}

workflow {
    sirv = copyInputFiles(inputFiles) | processForSIRVs
    output = sirv | combine(tolerances) | (runTmerge2 & runTmerge1) | mix
    sirv | mix(output) | gffToBED
    output | runGFFCompare
}