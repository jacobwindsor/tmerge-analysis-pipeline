nextflow.preview.dsl=2

/*
* tmerge 2.0 testing pipeline
*
* ALL FILE NAMES MUST FOLLOW THIS CONVENTION:
* [nickname].[output|input|time].[tool].[sirv?].[ext]
*/

/*
* Configuration options. 
*/
// Options that are only used in this context
sirvome_path = "/users/rg/jlagarde/genomes/lexogen_SIRVs/SIRV_Set1_Lot00141_Sequences_181206a/SIRVome_isoforms_Lot00141_C_181206a.gtf.unix.corrected.gene_types.gtf"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"

// Options used by included files and this context
params.output_dir = "/users/rg/jwindsor/tests/tmerge/results/splice_scoring"
params.julien_utils_path = "/users/rg/jlagarde/julien_utils/"
params.results_path = params.output_dir + "/results.csv"

// Remove time stats first
new File(params.results_path).delete()  

// Run on all high confidence
inputFiles = Channel.fromPath([
    "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/*.gff.gz",
    "/users/rg/jlagarde/projects/encode/scaling/whole_genome/202006_lrgasp/mappings/highConfidenceReads/*.gff.gz"
]
).map { file -> tuple(file.simpleName, file) }

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
    path inputGFF

    memory '30 GB'

    output:
        path('*.output.gff')

    shell:
    '''
    module load Python
    source !{venv}

    tmerge -i=!{inputGFF} -o=!{inputGFF.simpleName}.output.gff --splice_scoring --acceptor_path=/users/rg/jlagarde/projects/splice_sites/yeo_burge/Hsap.acceptor.mecoef --donor_path=/users/rg/jlagarde/projects/splice_sites/yeo_burge/Hsap.donor.mecoef --fasta_path=/users/rg/jwindsor/genomes/hg38.fa
    '''
}


process recordResults {
    // Append to CSV
    // Columns:
    input:
    tuple  path(output), path(sensitivity), path(precision)
    
    shell:
    '''
    touch !{params.results_path}

    PREC="$(cat !{precision} | tr -d '\n')"
    SENS="$(cat !{sensitivity} | tr -d '\n')"

    echo "!{output.simpleName},$SENS,$PREC" >> !{params.results_path}
    '''
}

process addGFFCompare {
    input:
    path output
    
    output:
    tuple path(output), path("sensitivity.txt"), path("precision.txt")

    shell:
    '''
    # Run GFF compare
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -o tmp.gffcompare !{output}  -r /users/rg/jwindsor/annotations/gencode/gencode.v34.annotation.exons.gtf

    cat tmp.gffcompare

    # Sensitivity
    cat tmp.gffcompare | grep "Transcript level" | grep -Eo '[0-9]*\\.[0-9]*' | awk '{i++}i==1' > sensitivity.txt
    # Precision
    cat tmp.gffcompare | grep "Transcript level" | grep -Eo '[0-9]*\\.[0-9]*' | awk '{i++}i==2' > precision.txt
    '''
}


workflow {
    getGFF(inputFiles) | runTmerge2 | addGFFCompare | recordResults
}