nextflow.preview.dsl=2

sirvome_path = "/users/rg/jwindsor/annotations/gencode.v34.annotation+SIRVs.exons.gtf"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"

params.julien_utils_path = "/users/rg/jlagarde/julien_utils/"
params.output_dir = "/users/rg/jwindsor/tests/tmerge/results/options"


// Just run with high confidence reads
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
    # If the file has SIRVs, use that for comparison for speed
    # Otherwise, compare to GENCODE v34
    if grep -q "SIRVome_isoforms" !{fileName}; then
        zcat !{fileName} | awk -F"\t" '$1=="SIRVome_isoforms"' > !{ID}_SIRVs.input.gff
    else
        zcat !{fileName} > !{ID}.input.gff
    fi
    '''    
}

process runTmerge2_isoform_fraction {
    input:
    tuple path(gff), min_length, tolerance, fuzz, isoform_fraction

    output:
        tuple path("*.output.gff"), val("isoform_fraction"), min_length, tolerance, fuzz, isoform_fraction

    shell:
        '''
        module load Python
        source !{venv}

        tmerge -i !{gff} -o !{gff.simpleName}.output.gff --tolerance=!{tolerance} --min_isoform_fraction=!{isoform_fraction} --end_fuzz=!{fuzz} --min_length=!{min_length}
        '''
}

process runTmerge2_support {
    input:
    tuple path(gff), min_length, tolerance, fuzz, support

    output:
        tuple path("*.output.gff"), val("support"), min_length, tolerance, fuzz, support

    shell:
        '''
        module load Python
        source !{venv}

        tmerge -i !{gff} -o !{gff.simpleName}.output.gff --tolerance=!{tolerance} --min_read_support=!{support} --end_fuzz=!{fuzz} --min_length=!{min_length}
        '''
}

process addGFFCompare {
    input:
    tuple path(gff), type, min_length, tolerance, fuzz, isoform_fraction_or_support
    
    output:
        tuple path(gff), type, min_length, tolerance, fuzz, isoform_fraction_or_support, path("sensitivity.txt"), path("precision.txt")

    shell:
    '''
    # Run GFF compare
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -o tmp.gffcompare !{gff}  -r !{sirvome_path} -R

    cat tmp.gffcompare

    # Sensitivity
    cat tmp.gffcompare | grep "Transcript level" | grep -Eo '[0-9]*\\.[0-9]*' | awk '{i++}i==1' > sensitivity.txt
    # Precision
    cat tmp.gffcompare | grep "Transcript level" | grep -Eo '[0-9]*\\.[0-9]*' | awk '{i++}i==2' > precision.txt
    '''
}

process recordResults {
    // Append to CSV
    input:
    tuple path(gff), type, min_length, tolerance, fuzz, isoform_fraction_or_support, path(sensitivity), path(precision)

    output:
        path("tmp.csv")
    
    shell:
    '''
    PREC="$(cat !{precision} | tr -d '\n')"
    SENS="$(cat !{sensitivity} | tr -d '\n')"

    echo "!{gff.simpleName},!{type},!{min_length},!{tolerance},!{fuzz},!{isoform_fraction_or_support},$SENS,$PREC" > tmp.csv
    '''
}

workflow {
    min_lengths = [0]
    tolerances = [0, 2, 4, 6, 8, 10, 12, 16, 20]
    isoform_fractions = [0, 0.01, 0.02, 0.04, 0.06, 0.09, 0.1, 0.2]
    fuzzs = [0, 1, 2, 4, 6, 8, 10]
    supports = [1, 2, 4, 6, 8, 10, 14, 20]

    // min_lengths = [0]
    // tolerances = [0]
    // isoform_fractions = [0]
    // fuzzs = [0]
    // supports = [1]

    combined = getGFF(inputFiles) | combine(min_lengths) | combine(tolerances) | combine(fuzzs)

    // Isoform fraction
    isoform_fraction = combined | combine(isoform_fractions) | runTmerge2_isoform_fraction
    // support
    support = combined | combine(supports) | runTmerge2_support

    isoform_fraction.mix(support) | addGFFCompare | recordResults \
        | collectFile(name: "results.csv", newLine: false, storeDir: params.output_dir)
}