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
tmerge1_path = "/users/rg/jlagarde/julien_utils"
flair_path = "/users/rg/jwindsor/flair"
stringtie2_path = "/users/rg/jwindsor/stringtie2"
read_mapping_path = "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/readMapping/"
sirvome_path = "/users/rg/jlagarde/genomes/lexogen_SIRVs/SIRV_Set1_Lot00141_Sequences_181206a/SIRVome_isoforms_Lot00141_C_181206a.gtf.unix.corrected.gene_types.gtf"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"

// Options used by included files and this context
params.output_dir = "/users/rg/jwindsor/tests/tmerge/results/compare"
params.julien_utils_path = "/users/rg/jlagarde/julien_utils/"
params.results_path = params.output_dir + "/results.csv"

// Remove time stats first
new File(params.results_path).delete()  

inputFiles = Channel.fromPath([
    // Run on all high confidence
    "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/*.gff.gz",
    // Include the two raw pacbio sets with SIRVs
    // "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/strandGffs/pacBio:Cshl:Smarter*.gff.gz"
    // "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/pacBio:Cshl:Smarter:Corr0_HpreCap_0+_Brain01Rep1.strandedHCGMs.gff.gz"
]).map { file -> tuple(file.simpleName, file) }

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

/* 
STEP 1a:
 */
process getBED {
    input:
    tuple ID, val(fileName)

    memory "30G"

    output:
    path '*.input.bed'

    shell:
    '''
    PATH="$PATH:!{params.julien_utils_path}"

    zcat !{fileName} | gff2bed_full.pl -| sortbed > !{ID}.input.bed
    '''
}

process getBAM {
    input:
    tuple ID, val(fileName)

    memory "30G"

    output:
    path '*.input.bam'
    
    shell:
        '''
        PATH="$PATH:!{params.julien_utils_path}"
        module load SAMtools/1.5-foss-2016b;

        # select all read IDs (transcript_id's) from GTF:
        zcat !{fileName} | extractGffAttributeValue.pl transcript_id | sort|uniq > tmp.ids
        
        # extract SAM header from the BAM file (if you don't, the subsequent SAM->BAM conversion will not work):
        samtools view -H !{read_mapping_path}!{ID}.bam > tmp.sam

        # extract SAM alignment records and filter them based on your list of read IDs.
        samtools view "!{read_mapping_path}!{ID}.bam" | fgrep -w -f tmp.ids >> tmp.sam

        # convert SAM to BAM:
        samtools view -b tmp.sam > !{ID}.input.bam
        '''
}

process runTmerge1 {
    input:
    path inputGFF

    memory '30 GB'
    time '10h'
    errorStrategy 'ignore'

    output:
        tuple path('*.output.tmerge1.gff'), path('*.time.tmerge1.txt') // time file outputed as time(s),peak_mem(kb)
    
    shell:
    '''
    PATH="$PATH:!{params.julien_utils_path}"
    /usr/bin/time -f "%e,%M" -o !{inputGFF.simpleName}.time.tmerge1.txt perl !{tmerge1_path}/tmerge --exonOverhangTolerance=4 --minReadSupport=5 --endFuzz=1 !{inputGFF} | sortgff > !{inputGFF.simpleName}.output.tmerge1.gff
    '''
}

process runTmerge2 {
    input:
    path inputGFF

    memory '30 GB'
    errorStrategy 'ignore'

    output:
        tuple path('*.output.tmerge2.gff'), path('*.time.tmerge2.txt') // time file outputed as time(s),peak_mem(kb)

    shell:
    '''
    module load Python
    source !{venv}

    /usr/bin/time -f "%e,%M" -o !{inputGFF.simpleName}.time.tmerge2.txt tmerge -i !{inputGFF} -o !{inputGFF.simpleName}.output.tmerge2.gff --tolerance=7 --end_fuzz=2 --min_read_support=1 --min_abundance=0.086
    '''
}

process runTmerge2_splice_scoring {
    input:
    path inputGFF

    memory '30 GB'
    errorStrategy 'ignore'

    output:
        tuple path('*.output.tmerge2.gff'), path('*.time.tmerge2.txt') // time file outputed as time(s),peak_mem(kb)

    shell:
    '''
    module load Python
    source !{venv}
    
    /usr/bin/time -f "%e,%M" -o !{inputGFF.simpleName}.time.tmerge2.txt tmerge -i !{inputGFF} -o !{inputGFF.simpleName}.output.tmerge2.gff --tolerance=7 --end_fuzz=0 --min_read_support=1 --min_abundance=0 --splice_scoring --acceptor_path=/users/rg/jlagarde/projects/splice_sites/yeo_burge/Hsap.acceptor.mecoef --donor_path=/users/rg/jlagarde/projects/splice_sites/yeo_burge/Hsap.donor.mecoef --fasta_path=/users/rg/jwindsor/genomes/hg38.fa
    '''
}

process runFLAIR {
    // default params (Tested on pacBio:Cshl:Smarter:Corr0_HpreCap_0+_Brain01Rep1.input.gff (standard) to get best results)
    input:
    path inputBED

    memory '30GB'
    errorStrategy 'ignore'

    output:
        tuple path('*.output.flair.gff'), path('*.time.flair.txt') // time file outputed as time(s),peak_mem(kb)

    shell:
        '''
        /usr/bin/time -f "%e,%M" -o !{inputBED.simpleName}.time.flair.txt time python !{flair_path}/bin/collapse_isoforms_precise.py -s 0.8 -q !{inputBED} -o tmp.bed
        
        # Convert to GFF
        PATH="$PATH:!{params.julien_utils_path}"

        cat tmp.bed | sortgff | awk -f !{params.julien_utils_path}bed12fields2gff.awk > !{inputBED.simpleName}.output.flair.gff
        '''
}

process runStringTie2 {
    // Options used. Tested on pacBio:Cshl:Smarter:Corr0_HpreCap_0+_Brain01Rep1.input.bam (standard) to get best results
    // -L (long read mode)
    // -t (disable trimming)
    input:
    path inputBAM

    errorStrategy 'ignore'

    output:
        tuple path('*.output.stringtie2.gff'), path('*.time.stringtie2.txt') // time file outputed as time(s),peak_mem(kb)

    shell:
    '''
    PATH="$PATH:!{params.julien_utils_path}"
    /usr/bin/time -f "%e,%M" -o !{inputBAM.simpleName}.time.stringtie2.txt time !{stringtie2_path}/stringtie -L !{inputBAM} -o !{inputBAM.simpleName}.output.stringtie2.tmp.gff
    cat !{inputBAM.simpleName}.output.stringtie2.tmp.gff | sortgff > !{inputBAM.simpleName}.output.stringtie2.gff
    '''
}

process addGFFCompare {
    input:
    tuple type, path(output), path(time)
    
    output:
    tuple type, path(output), path(time), path("sensitivity.txt"), path("precision.txt")

    shell:
    '''
    # Retrieve only SIRVs
    cat !{output} | awk -F"\t" '$1=="SIRVome_isoforms"' > sirvs.gff

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
    // Type, dataset name, time (s), peak mem (kb), sensitivity (%), precision (%)
    input:
    tuple type, path(output), path(time), path(sensitivity), path(precision)
    
    shell:
    '''
    touch !{params.results_path}

    TIME="$(cat !{time} | tr -d '\n')"
    PREC="$(cat !{precision} | tr -d '\n')"
    SENS="$(cat !{sensitivity} | tr -d '\n')"

    echo "!{type},!{output.simpleName},$TIME,$SENS,$PREC\n" >> !{params.results_path}
    '''
}

workflow runTools {
    take: inputFiles
    main:
        inputGFFs = inputFiles.filter(~/.*\.gff$/ )
        inputBAMs = inputFiles.filter( ~/.*\.bam$/ )
        inputBEDs = inputFiles.filter( ~/.*\.bed$/ )
        
        tmerge1 = runTmerge1(inputGFFs)
        tmerge2 = runTmerge2(inputGFFs)
        tmerge2_splice_scoring = runTmerge2_splice_scoring(inputGFFs)

        flair = runFLAIR(inputBEDs)

        stringtie2 = runStringTie2(inputBAMs)
    emit:
        tmerge1.map { out -> tuple("tmerge1", out[0], out[1]) }
            .mix(tmerge2.map { out -> tuple("tmerge2", out[0], out[1]) })
            .mix(tmerge2_splice_scoring.map { out -> tuple("tmerge2_splice_scoring", out[0], out[1]) })
            .mix(flair.map { out -> tuple("FLAIR", out[0], out[1]) })
            .mix(stringtie2.map{ out -> tuple("StringTie2", out[0], out[1]) })
}

workflow {
    inputs = getGFF(inputFiles)
        .mix(getBAM(inputFiles))
        .mix(getBED(inputFiles))
    
    runTools(inputs) | addGFFCompare | recordResults
}