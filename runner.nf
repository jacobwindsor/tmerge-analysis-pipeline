nextflow.preview.dsl=2

include runGFFCompare as oneVsTwo from './gffcompare'
include runGFFCompare as twoVsFLAIR from './gffcompare'
include runGFFCompare as twoVsStringTie2 from './gffcompare'



/*
* tmerge 2.0 testing pipeline
*
* ALL FILE NAMES MUST FOLLOW THIS CONVENTION:
* [nickname].[output|input|time].[tool].[ext]
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
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"
output_dir = "/users/rg/jwindsor/tests/tmerge/results_mock"
cache_dir = "/users/rg/jwindsor/tests/tmerge/cache"


inputFiles = Channel.fromList([
    [
        "nickname": "standard",
        "path": "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/pacBio:Cshl:Smarter:Corr0_HpreCap_0+_Brain01Rep1.strandedHCGMs.gff.gz"
    ],
    // [
    //     "nickname": "tricky",
    //     "path": "/users/project/gencode_006070_no_backup/jlagarde/lncRNACapture_phase3/mappings/strandGffs/ont:Cshl:CapTrap:Corr0_HpreCap_0+_Heart01Rep1.stranded.gff.gz"
    // ]
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

    output:
    path '*.input.bam'
    
    shell:
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
    /usr/bin/time -f "%e" -o !{inputGFF.simpleName}.time.tmerge1.txt perl !{tmerge1_path}/tmerge !{inputGFF} > !{inputGFF.simpleName}.output.tmerge1.gff
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
    /usr/bin/time -f "%e" -o !{inputBAM.simpleName}.time.stringtie2.txt time !{stringtie2_path}/stringtie -L -f 0 !{inputBAM} -o !{inputBAM.simpleName}.output.stringtie2.gff
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
    cat !{inputBED} | awk -f !{julien_utils_path}bed12fields2gff.awk > !{inputBED.baseName}.gff
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
        recordTime(tmerge2.time)
        oneVsTwo(tmerge2.output, tmerge1.output)

        flair = gffToBED(inputGFFs) | runFLAIR
        flairGFF = bedToGFF(flair.output)
        twoVsFLAIR(tmerge2.output, flairGFF)

        stringtie2 = runStringTie2(inputBAMs)
        twoVsStringTie2(tmerge2.output, stringtie2.output)
    emit:
        tmerge1.output.mix(tmerge2.output).mix(flair.output)
}

workflow tmerge1 {
    take: inputGFFs
    main:
        tmerge1 = runTmerge1(inputGFFs)
        tmerge2 = runTmerge2(inputGFFs)    
        recordTime(tmerge2.time)
        oneVsTwo(tmerge2.output, tmerge1.output)
}

workflow convertToBed {
    take: inputs
    main:
        gffToBED(inputs)
}

workflow {
    inputs = copyInputFiles(inputFiles).mix(getBAM(inputFiles))
    outputs = runTools(inputs)
    outputs.filter( ~/.*\.gff$/ ) | convertToBed
}