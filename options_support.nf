nextflow.preview.dsl=2

tmerge1_path = "/users/rg/jlagarde/julien_utils"
sirvome_path = "/users/rg/jwindsor/genomes/SIRVome.gff"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"

params.julien_utils_path = "/users/rg/jlagarde/julien_utils/"
params.output_dir = "/users/rg/jwindsor/tests/tmerge/results/options_support"

// Just run with PacBio high confidence reads
inputFiles = Channel.fromPath([
    "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/*.gff.gz",
    "/users/rg/jlagarde/projects/encode/scaling/whole_genome/202006_lrgasp/mappings/highConfidenceReads/*.gff.gz"
]).map { file -> tuple(file.simpleName, file) }

process getGFF {
    input:
    tuple ID, val(fileName)
    
    output:
    path '*.input.gff'

    shell:
    '''
    # Retrieve only SIRVs
    zcat !{fileName} | awk -F"\t" '$1=="SIRVome_isoforms"' > !{ID}.input.gff
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

process addGFFCompare {
    input:
    tuple path(gff), tolerance, support, fuzz
    
    output:
        tuple path(gff), tolerance, support, fuzz, path("sensitivity.txt"), path("precision.txt")

    shell:
    '''
    # Run GFF compare
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -o tmp.gffcompare !{gff}  -r !{sirvome_path}

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
    tuple path(gff), tolerance, support, fuzz, path(sensitivity), path(precision)

    output:
        path("tmp.csv")
    
    shell:
    '''
    touch !{params.results_path}

    PREC="$(cat !{precision} | tr -d '\n')"
    SENS="$(cat !{sensitivity} | tr -d '\n')"

    echo "!{gff.simpleName},!{tolerance},!{support},!{fuzz},$SENS,$PREC" > tmp.csv
    '''
}

workflow {
    tolerances = [0, 2, 4, 6, 8, 10, 12, 14, 16, 20]
    support = [1,2,3,4,5,6,8,10,12,14,16,18,20]
    fuzzs = [0, 1, 2, 4, 6, 8, 10]

    combined = getGFF(inputFiles)  | combine(tolerances) | combine(support) | combine(fuzzs)
    runTmerge2(combined) | addGFFCompare | recordResults | collectFile(name: "results.csv", newLine: true, storeDir: params.output_dir)
}