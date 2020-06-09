nextflow.preview.dsl=2

tmerge1_path = "/users/rg/jlagarde/julien_utils"
sirvome_path = "/users/rg/jwindsor/genomes/SIRVome.gff"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"

params.julien_utils_path = "/users/rg/jlagarde/julien_utils/"
params.output_dir = "/users/rg/jwindsor/tests/tmerge/results/tolerances"
params.results_path = params.output_dir + "/results.csv"

// Remove time stats first
new File(params.results_path).delete()  

include processForSIRVs from './utils'
include gffToBED from './utils'

// Just run with PacBio high confidence reads
inputFiles = Channel.fromPath(
    // "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/pacBio:Cshl:Smarter:Corr0_HpreCap_0+_Brain01Rep1.strandedHCGMs.gff.gz"
    "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/pacBio*.gff.gz"
).map { file -> tuple(file.simpleName, file) }

tolerances = [0,1,2,3,5]
min_reads = [0,2,3,4,5]
end_fuzz = [0,1,2,4,6]

process getGFF {
    input:
    tuple ID, val(fileName)
    
    output:
    path '*.input.gff'

    shell:
    '''
    zcat !{fileName} > !{ID}.input.gff
    '''    
}

process runTmerge2 {
    input:
    tuple path(gff), tolerance, support, fuzz

    output:
        tuple path("*.output.gff"), tolerance, support, fuzz

    shell:
    '''
    module load Python
    source !{venv}

    tmerge -i !{gff} -o !{gff.simpleName}.output.gff --tolerance=!{tolerance} --min_read_support=!{support} --end_fuzz=!{fuzz}
    '''
}

process runTmerge1 {
    input:
    tuple path(gff), tolerance, support, fuzz

    memory "10G"

    output:
        tuple path("*.output.gff"), tolerance, support, fuzz

    shell:
    '''
    PATH="$PATH:!{params.julien_utils_path}"

    perl !{tmerge1_path}/tmerge --exonOverhangTolerance=!{tolerance} --minReadSupport=!{support} --endFuzz=!{fuzz} !{gff} | sortgff > !{gff.simpleName}.output.gff
    '''
}

process addGFFCompare {
    input:
    tuple type, path(gff), tolerance, support, fuzz
    
    output:
        tuple type, path(gff), tolerance, support, fuzz, path("sensitivity.txt"), path("precision.txt")

    shell:
    '''
    # Retrieve only SIRVs
    cat !{gff} | awk -F"\t" '$1=="SIRVome_isoforms"' > sirvs.gff

    # Run GFF compare
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -o tmp.gffcompare sirvs.gff  -r !{sirvome_path}

    cat tmp.gffcompare

    # Sensitivity
    cat tmp.gffcompare | grep "Transcript level" | grep -Eo '[0-9]*\\.[0-9]*' | awk '{i++}i==1' > sensitivity.txt
    # Precision
    cat tmp.gffcompare | grep "Transcript level" | grep -Eo '[0-9]*\\.[0-9]*' | awk '{i++}i==2' > precision.txt
    '''
}

process recordResults {
    // Append to CSV
    // Columns:
    // Type, dataset name, tolerance, support, fuzz, sensitivity (%), precision (%)
    input:
    tuple type, path(gff), tolerance, support, fuzz, path(sensitivity), path(precision)
    
    shell:
    '''
    touch !{params.results_path}

    PREC="$(cat !{precision} | tr -d '\n')"
    SENS="$(cat !{sensitivity} | tr -d '\n')"

    echo "!{type},!{gff.simpleName},!{tolerance},!{support},!{fuzz},$SENS,$PREC" >> !{params.results_path}
    '''
}

workflow {
    sirv = getGFF(inputFiles) | processForSIRVs
    combined = sirv | combine(tolerances) | combine(min_reads) | combine(end_fuzz)
    runTmerge1(combined).map { x -> tuple("tmerge1", x[0], x[1], x[2], x[3]) }
        .mix(runTmerge2(combined).map {x -> tuple("tmerge2",  x[0], x[1], x[2], x[3]) }) | addGFFCompare | recordResults
}