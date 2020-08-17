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
sirvome_path = "/users/rg/jwindsor/annotations/SIRVome_isoforms_Lot00141_C_181206a.gtf.unix.corrected.gtf"
gencode_path = "/users/rg/jwindsor/annotations/gencode/gencode.v34.annotation.gtf"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"

// Options used by included files and this context
params.output_dir = "/users/rg/jwindsor/tests/tmerge/results/compare"
params.julien_utils_path = "/users/rg/jlagarde/julien_utils/"

inputFiles = Channel.fromPath([
    // Run on all high confidence
    "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/*.gff.gz",
    // "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/pacBio:Cshl:Smarter:Corr0_HpreCap_0+_Heart01Rep1.strandedHCGMs.gff.gz"
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
    
    if [[ !{inputGFF.simpleName} == *"pacBio"* ]]; then
        /usr/bin/time -f "%e,%M" -o !{inputGFF.simpleName}.time.tmerge2.txt tmerge -i !{inputGFF} -o !{inputGFF.simpleName}.output.tmerge2.gff --pacbio
    else
        /usr/bin/time -f "%e,%M" -o !{inputGFF.simpleName}.time.tmerge2.txt tmerge -i !{inputGFF} -o !{inputGFF.simpleName}.output.tmerge2.gff --ont
    fi
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
    
    if [[ !{inputGFF.simpleName} == *"pacBio"* ]]; then
        /usr/bin/time -f "%e,%M" -o !{inputGFF.simpleName}.time.tmerge2.txt tmerge -i !{inputGFF} -o !{inputGFF.simpleName}.output.tmerge2.gff --pacbio --splice_scoring --acceptor_path=/users/rg/jlagarde/projects/splice_sites/yeo_burge/Hsap.acceptor.mecoef --donor_path=/users/rg/jlagarde/projects/splice_sites/yeo_burge/Hsap.donor.mecoef --fasta_path=/users/rg/jwindsor/genomes/hg38.fa
    else
        /usr/bin/time -f "%e,%M" -o !{inputGFF.simpleName}.time.tmerge2.txt tmerge -i !{inputGFF} -o !{inputGFF.simpleName}.output.tmerge2.gff --ont --splice_scoring --acceptor_path=/users/rg/jlagarde/projects/splice_sites/yeo_burge/Hsap.acceptor.mecoef --donor_path=/users/rg/jlagarde/projects/splice_sites/yeo_burge/Hsap.donor.mecoef --fasta_path=/users/rg/jwindsor/genomes/hg38.fa
    fi
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
        /usr/bin/time -f "%e,%M" -o !{inputBED.simpleName}.time.flair.txt time python !{flair_path}/bin/collapse_isoforms_precise.py -q !{inputBED} -o tmp.bed
        
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

process addGFFCompare_sirvs {
    input:
    tuple type, path(output), path(time)
    
    output:
    tuple type, path(output), path(time), path("sensitivity_sirvs.txt"), path("precision_sirvs.txt")

    shell:
    '''
    # extract only SIRVs
    cat !{output} | awk -F"\t" '$1=="SIRVome_isoforms"' > sirvs.gff 

    # Run GFF compare
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -R -o tmp.gffcompare sirvs.gff -r !{sirvome_path}

    # Sensitivity
    cat tmp.gffcompare | grep "Transcript level" | grep -Eo '[0-9]*\\.[0-9]*' | awk '{i++}i==1' > sensitivity_sirvs.txt
    # Precision
    cat tmp.gffcompare | grep "Transcript level" | grep -Eo '[0-9]*\\.[0-9]*' | awk '{i++}i==2' > precision_sirvs.txt
    '''
}

process addGFFCompare_gencode {
    input:
    tuple type, path(output), path(time), path(sensitivity_sirvs), path(precision_sirvs)
    
    output:
    tuple type, path(output), path(time), path(sensitivity_sirvs), path(precision_sirvs), path("sensitivity_gencode.txt"), path("precision_gencode.txt")

    shell:
    '''
    # Run GFF compare
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -R -o tmp.gffcompare !{output} -r !{gencode_path}

    # Sensitivity
    cat tmp.gffcompare | grep "Transcript level" | grep -Eo '[0-9]*\\.[0-9]*' | awk '{i++}i==1' > sensitivity_gencode.txt
    # Precision
    cat tmp.gffcompare | grep "Transcript level" | grep -Eo '[0-9]*\\.[0-9]*' | awk '{i++}i==2' > precision_gencode.txt
    '''
}



process recordResults {
    // Append to CSV
    input:
    tuple type, path(output), path(time), path(sensitivity_sirvs), path(precision_sirvs), path(sensitivity_gencode), path(precision_gencode)
    
    output:
    path 'tmp.csv'

    shell:
    '''
    TIME="$(cat !{time} | tr -d '\n')"
    PREC_SIRV="$(cat !{precision_sirvs} | tr -d '\n')"
    SENS_SIRV="$(cat !{sensitivity_sirvs} | tr -d '\n')"
    PREC_GENCODE="$(cat !{precision_gencode} | tr -d '\n')"
    SENS_GENCODE="$(cat !{sensitivity_gencode} | tr -d '\n')"

    # If the file has SIRVs, flag it
    SIRVS="FALSE"
    if grep -q "SIRVome_isoforms" !{output} ; then
        SIRVS="TRUE"
    fi

    echo "!{type},!{output.simpleName},$SIRVS,$TIME,$SENS_SIRV,$PREC_SIRV,$SENS_GENCODE,$PREC_GENCODE" > tmp.csv
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
    
    runTools(inputs) | addGFFCompare_sirvs | addGFFCompare_gencode | recordResults \
        | collectFile(name: "results.csv", newLine: false, storeDir: params.output_dir, sort: false)
}