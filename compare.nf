nextflow.preview.dsl=2


include processForSIRVs as processTmerge2SIRVs from './utils'
include processForSIRVs as processFLAIRSIRVs from './utils'
include processForSIRVs as processStringtie2SIRVs from './utils'
include runGFFCompare as oneVsTwo from './utils'
include runGFFCompare as FLAIRCompare from './utils'
include runGFFCompare as StringTie2Compare from './utils'
include runGFFCompare as tmerge2Compare from './utils'




/*
* tmerge 2.0 testing pipeline
*
* ALL FILE NAMES MUST FOLLOW THIS CONVENTION:
* [nickname].[output|input|time].[tool].[sirv?].[ext]
*/

/*
* Configuration options
*/
tmerge2_path = "/users/rg/jwindsor/tmerge/2.0"
tmerge1_path = "/users/rg/jlagarde/julien_utils"
flair_path = "/users/rg/jwindsor/flair"
stringtie2_path = "/users/rg/jwindsor/stringtie2"
julien_utils_path = "/users/rg/jlagarde/julien_utils/"
read_mapping_path = "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/readMapping/"
sirvs_path = "/users/rg/jlagarde/genomes/lexogen_SIRVs/SIRV_Set1_Lot00141_Sequences_181206a/SIRVome_isoforms_Lot00141_C_181206a.gtf.unix.corrected.gene_types.gtf"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"
output_dir = "/users/rg/jwindsor/tests/tmerge/results"
cache_dir = "/users/rg/jwindsor/tests/tmerge/cache"


inputFiles = Channel.fromList([
    [
        "nickname": "standard",
        "path": "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/pacBio:Cshl:Smarter:Corr0_HpreCap_0+_Brain01Rep1.strandedHCGMs.gff.gz"
    ],
    [
        "nickname": "tricky",
        "path": "/users/project/gencode_006070_no_backup/jlagarde/lncRNACapture_phase3/mappings/strandGffs/ont:Cshl:CapTrap:Corr0_HpreCap_0+_Heart01Rep1.stranded.gff.gz"
    ]
])

process copyInputFiles {
    input:
    val x

    publishDir "$output_dir"

    output:
    path '*.input.gff'

    shell:
    '''
    echo !{x.nickname}
    zcat "!{x.path}" > !{x.nickname}.input.gff
    '''    
}

process getBAM {
    input:
    val x

    publishDir "$output_dir"
    memory "30G"

    output:
    path '*.input.bam'
    
    shell:
    if (x.nickname == "standard")
        '''
        PATH="$PATH:!{julien_utils_path}"
        module load SAMtools/1.5-foss-2016b;

        # select all read IDs (transcript_id's) from GTF:
        zcat !{x.path} | extractGffAttributeValue.pl transcript_id | sort|uniq > tmp.ids
        
        # extract SAM header from the BAM file (if you don't, the subsequent SAM->BAM conversion will not work):
        samtools view -H !{read_mapping_path}$(basename !{x.path} .strandedHCGMs.gff.gz).bam > tmp.sam

        # extract SAM alignment records and filter them based on your list of read IDs.
        samtools view "!{read_mapping_path}$(basename !{x.path} .strandedHCGMs.gff.gz).bam" | fgrep -w -f tmp.ids >> tmp.sam

        # convert SAM to BAM:
        samtools view -b tmp.sam > !{x.nickname}.input.bam
        '''
}

process runTmerge1 {
    input:
    path inputGFF

    memory '30 GB'
    time '25h'
    publishDir "$output_dir"

    output:
        path '*.output.tmerge1.gff', emit: output
        path '*.time.tmerge1.txt', emit: time
    
    shell:
    '''
    PATH="$PATH:!{julien_utils_path}"
    /usr/bin/time -f "%e" -o !{inputGFF.simpleName}.time.tmerge1.txt perl !{tmerge1_path}/tmerge !{inputGFF} | sortgff > !{inputGFF.simpleName}.output.tmerge1.gff
    '''
}

process runTmerge2 {
    input:
    path inputGFF

    memory '30 GB'
    publishDir "$output_dir"

    output:
        path '*.output.tmerge2.gff', emit: output
        path '*.time.tmerge2.txt', emit: time

    shell:
    '''
    module load Python
    source !{venv}

    /usr/bin/time -f "%e" -o !{inputGFF.simpleName}.time.tmerge2.txt python !{tmerge2_path}/tmerge.py -i !{inputGFF} -o !{inputGFF.simpleName}.output.tmerge2.gff
    '''
}

process runFLAIR {
    input:
    path inputBED

    memory '30GB'
    publishDir "$output_dir"

    output:
        path '*.output.flair.bed', emit: output
        path '*.time.flair.txt', emit: time

    shell:
        '''
        /usr/bin/time -f "%e" -o !{inputBED.simpleName}.time.flair.txt time python !{flair_path}/bin/collapse_isoforms_precise.py -q !{inputBED} -o !{inputBED.simpleName}.output.flair.bed
        '''
}

process runStringTie2 {
    input:
    path inputBAM

    publishDir "$output_dir"

    output:
        path '*.output.stringtie2.gff', emit: output
        path '*.time.stringtie2.txt', emit: time

    shell:
    '''
    PATH="$PATH:!{julien_utils_path}"
    /usr/bin/time -f "%e" -o !{inputBAM.simpleName}.time.stringtie2.txt time !{stringtie2_path}/stringtie -L -f 0 -a 1 !{inputBAM} -o !{inputBAM.simpleName}.output.stringtie2.tmp.gff
    cat !{inputBAM.simpleName}.output.stringtie2.tmp.gff | sortgff > !{inputBAM.simpleName}.output.stringtie2.gff
    '''
}

process gffToBED {
    input:
    path inputGFF

    publishDir "$output_dir"

    output:
    path '*.bed'

    shell:
    '''
    PATH="$PATH:!{julien_utils_path}"

    cat !{inputGFF} | gff2bed_full.pl -| sortbed > !{inputGFF.baseName}.bed
    '''
}

process bedToGFF {
    input:
    path inputBED

    publishDir "$output_dir"

    output:
    path '*.gff'

    shell:
    '''
    PATH="$PATH:!{julien_utils_path}"

    cat !{inputBED} | sortgff | awk -f !{julien_utils_path}bed12fields2gff.awk > !{inputBED.baseName}.gff
    '''
}

process recordTime {
    input:
    path tmerge2_time

    shell:
    '''
    LAST_COMMIT_ID="$(cd !{tmerge2_path} && git log --format="%H" -n 1)"

    TMERGE2_TIME="$(cat !{tmerge2_time} | tr -d '\n')"

    echo "$LAST_COMMIT_ID,!{tmerge2_time.simpleName},$TMERGE2_TIME" >> !{output_dir}/tmerge2_time_stats.csv
    '''
}

workflow runTools {
    take: inputFiles
    main:
        inputGFFs = inputFiles.filter(~/.*\.gff$/ )
        inputBAMs = inputFiles.filter( ~/.*\.bam$/ )
        
        tmerge1 = runTmerge1(inputGFFs)
        tmerge2 = runTmerge2(inputGFFs)
        tmerge2_sirvs = tmerge2.output | processTmerge2SIRVs
        recordTime(tmerge2.time)
        // NOTE: All comparisons are done against SIRVs except tmerge1 vs tmerge2
        oneVsTwo(tmerge2.output, tmerge1.output)
        tmerge2Compare(tmerge2_sirvs, sirvs_path)

        flair = gffToBED(inputGFFs) | runFLAIR
        flairGFF = bedToGFF(flair.output) | processFLAIRSIRVs
        FLAIRCompare(flairGFF, sirvs_path)

        stringtie2 = runStringTie2(inputBAMs)
        stringtie2GFF = stringtie2.output | processStringtie2SIRVs
        StringTie2Compare(stringtie2GFF, sirvs_path)
    emit:
        tmerge1.output.mix(tmerge2.output).mix(flair.output).mix(stringtie2.output)
}

workflow convertToBed {
    take: inputs
    main:
        gffToBED(inputs)
}

workflow {
    inputs = copyInputFiles(inputFiles)
        .mix(getBAM(inputFiles.filter { x -> x.nickname == "standard" })) // Only run for standard. Tricky doesnt work for getting BAM
    outputs = runTools(inputs)
    outputs.filter( ~/.*\.gff$/ ) | convertToBed
}