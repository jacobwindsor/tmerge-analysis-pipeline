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
tmerge2_path = "/users/rg/jwindsor/tmerge/2.0"
tmerge1_path = "/users/rg/jlagarde/julien_utils"
flair_path = "/users/rg/jwindsor/flair"
stringtie2_path = "/users/rg/jwindsor/stringtie2"
read_mapping_path = "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/readMapping/"
sirvs_path = "/users/rg/jlagarde/genomes/lexogen_SIRVs/SIRV_Set1_Lot00141_Sequences_181206a/SIRVome_isoforms_Lot00141_C_181206a.gtf.unix.corrected.gene_types.gtf"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"
cupcake_venv = "/users/rg/jwindsor/venvs/cupcake/bin/activate"

// Options used by included files and this context
params.output_dir = "/users/rg/jwindsor/tests/tmerge/results/compare"
params.julien_utils_path = "/users/rg/jlagarde/julien_utils/"

include processForSIRVs as processTmerge2SIRVs from './utils'
include processForSIRVs as processFLAIRSIRVs from './utils'
include processForSIRVs as processStringtie2SIRVs from './utils'
include processForSIRVs as processCupcakeSIRVs from './utils'
include runGFFCompare as oneVsTwo from './utils'
include runGFFCompare as FLAIRCompare from './utils'
include runGFFCompare as StringTie2Compare from './utils'
include runGFFCompare as tmerge2Compare from './utils'
include runGFFCompare as CupcakeCompare from './utils'
include gffToBED from './utils'

inputFiles2 = Channel.fromList([
    // [
    //     "tricky",
    //     "/users/project/gencode_006070_no_backup/jlagarde/lncRNACapture_phase3/mappings/strandGffs/ont:Cshl:CapTrap:Corr0_HpreCap_0+_Heart01Rep1.stranded.gff.gz"
    // ],
    [
        "standard",
        "/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/pacBio:Cshl:Smarter:Corr0_HpreCap_0+_Brain01Rep1.strandedHCGMs.gff.gz"
    ]
])

// inputFiles2 = Channel
//     .fromPath("/users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/mappings/highConfidenceReads/*01Rep1.strandedHCGMs.gff.gz")
//     .map { x -> [ x.simpleName, x ] }
//     .mix(inputFiles)


process copyInputGFF {
    input:
    val x

    publishDir "$params.output_dir"

    output:
    path '*.input.gff'

    shell:
    '''
    echo !{x[0]}
    zcat "!{x[1]}" > !{x[0]}.input.gff
    '''    
}

process getSAM {
    input:
    val x

    publishDir "$params.output_dir"
    memory "30G"
    errorStrategy "ignore"

    output:
    path '*.input.sam'
    
    shell:
    if (x[0] == "tricky")
        '''
        PATH="$PATH:!{params.julien_utils_path}"
        module load SAMtools/1.5-foss-2016b;

        # select all read IDs (transcript_id's) from GTF:
        zcat !{x[1]} | extractGffAttributeValue.pl transcript_id | sort|uniq > tmp.ids
        
        # extract SAM header from the BAM file (if you don't, the subsequent SAM->BAM conversion will not work):
        samtools view -H !{read_mapping_path}$(basename !{x[1]} .stranded.gff.gz).bam > !{x[0]}.input.sam

        # extract SAM alignment records and filter them based on your list of read IDs.
        samtools view "!{read_mapping_path}$(basename !{x[1]} .stranded.gff.gz).bam" | fgrep -w -f tmp.ids >> !{x[0]}.input.sam
        '''
    else
        '''
        PATH="$PATH:!{params.julien_utils_path}"
        module load SAMtools/1.5-foss-2016b;

        # select all read IDs (transcript_id's) from GTF:
        zcat !{x[1]} | extractGffAttributeValue.pl transcript_id | sort|uniq > tmp.ids
        
        # extract SAM header from the BAM file (if you don't, the subsequent SAM->BAM conversion will not work):
        samtools view -H !{read_mapping_path}$(basename !{x[1]} .strandedHCGMs.gff.gz).bam > !{x[0]}.input.sam

        # extract SAM alignment records and filter them based on your list of read IDs.
        samtools view "!{read_mapping_path}$(basename !{x[1]} .strandedHCGMs.gff.gz).bam" | fgrep -w -f tmp.ids >> !{x[0]}.input.sam
        '''
}

process getFQ {
    input:
    val x

    output:
    path '*.input.fq'

    shell:
    if (x[0] == "tricky")
        '''
        FQNAME=$((basename !{x[1]} .stranded.gff.gz) | sed "s/Corr0//" | sed "s/_0+_/_0+./").fastq.gz

        zcat /users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/fastqs/$FQNAME > !{x[0]}.input.fq
        '''
    else
        '''
        FQNAME=$((basename !{x[1]} .strandedHCGMs.gff.gz) | sed "s/Corr0//" | sed "s/_0+_/_0+./").fastq.gz

        zcat /users/rg/jlagarde/projects/encode/scaling/whole_genome/lncRNACapture_phase3/fastqs/$FQNAME > !{x[0]}.input.fq
        '''
}

process SAMToBAM {
    input:
    val samFile

    output:
    path '*.input.bam'

    shell:
    '''
    samtools view -Sb !{samFile} > !{samFile.simpleName}.input.bam
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

process runTmerge1 {
    input:
    val input

    memory '30 GB'
    time '25h'
    publishDir "$params.output_dir"

    output:
        path '*.output.tmerge1.gff', emit: output
        path '*.time.tmerge1.txt', emit: time
    
    shell:
    '''
    PATH="$PATH:!{params.julien_utils_path}"
    /usr/bin/time -f "%e" -o !{input[0]}.time.tmerge1.txt perl !{tmerge1_path}/tmerge !{input[1]} | sortgff > !{input[0]}.output.tmerge1.gff
    '''
}

process runTmerge2 {
    input:
    val input

    memory '30 GB'
    publishDir "$params.output_dir"

    output:
        path '*.output.tmerge2.gff', emit: output
        path '*.time.tmerge2.txt', emit: time

    shell:
    '''
    module load Python
    source !{venv}

    /usr/bin/time -f "%e" -o !{input[0]}.time.tmerge2.txt python !{tmerge2_path}/tmerge.py -i !{input[1]} -o !{input[0]}.output.tmerge2.gff
    '''
}

process runFLAIR {
    input:
    val input

    memory '30GB'
    publishDir "$params.output_dir"
    errorStrategy "ignore"

    output:
        path '*.output.flair.bed', emit: output
        path '*.time.flair.txt', emit: time

    shell:
        '''
        /usr/bin/time -f "%e" -o !{input[0]}.time.flair.txt time python !{flair_path}/bin/collapse_isoforms_precise.py -q !{input[5]} -o !{input[0]}.output.flair.bed
        '''
}

process runStringTie2 {
    input:
    val input

    publishDir "$params.output_dir"
    errorStrategy "ignore"

    output:
        path '*.output.stringtie2.gff', emit: output
        path '*.time.stringtie2.txt', emit: time

    shell:
    '''
    PATH="$PATH:!{params.julien_utils_path}"
    /usr/bin/time -f "%e" -o !{input[0]}.time.stringtie2.txt time !{stringtie2_path}/stringtie -L -f 0 -a 1 !{input[4]} -o !{input[0]}.output.stringtie2.tmp.gff
    cat !{input[0]}.output.stringtie2.tmp.gff | sortgff > !{input[0]}.output.stringtie2.gff
    '''
}

process runCupcake {
    input:
    val input

    publishDir "$params.output_dir"

    output:
        path '*.output.cupcake.collapsed.gff', emit: output
        path '*.time.cupcake.txt', emit: time

    shell:
    '''
    module load Python
    source !{cupcake_venv}
    
    /usr/bin/time -f "%e" -o !{input[0]}.time.cupcake.txt collapse_isoforms_by_sam.py --input !{input[3]} --fq -s !{input[2]} -o !{input[0]}.output.cupcake
    '''
}

process recordTime {
    input:
    path tmerge2_time

    shell:
    '''
    LAST_COMMIT_ID="$(cd !{tmerge2_path} && git log --format="%H" -n 1)"

    TMERGE2_TIME="$(cat !{tmerge2_time} | tr -d '\n')"

    echo "$LAST_COMMIT_ID,!{tmerge2_time.simpleName},$TMERGE2_TIME" >> !{params.output_dir}/tmerge2_time_stats.csv
    '''
}

workflow runTools {
    take: inputs // [[ nickname, gff, sam, fq, bam, bed ]]
    main:
        // TMERGE 1 and 2
        tmerge1 = runTmerge1(inputs)
        tmerge2 = runTmerge2(inputs)
        tmerge2_sirvs = tmerge2.output | processTmerge2SIRVs
        recordTime(tmerge2.time)
        // NOTE: All comparisons are done against SIRVs except tmerge1 vs tmerge2
        oneVsTwo(tmerge2.output, tmerge1.output)
        tmerge2Compare(tmerge2_sirvs, sirvs_path)

        // FLAIR
        flair = runFLAIR(inputs)
        flairGFF = bedToGFF(flair.output) | processFLAIRSIRVs
        FLAIRCompare(flairGFF, sirvs_path)

        // StringTie2
        stringtie2 = runStringTie2(inputs)
        stringtie2GFF = stringtie2.output | processStringtie2SIRVs
        StringTie2Compare(stringtie2GFF, sirvs_path)

        // Cupcake
        cupcake = runCupcake(inputs)
        cupcakeGFF = cupcake.output | processCupcakeSIRVs
        CupcakeCompare(cupcakeGFF, sirvs_path)
    emit:
        tmerge1.output.mix(tmerge2.output).mix(flair.output).mix(stringtie2.output).mix(cupcake.output)
}

workflow convertToBed {
    take: inputs
    main:
        gffToBED(inputs)
}

workflow {
    gff = copyInputGFF(inputFiles2)| map { x -> [ x.simpleName, x ] }
    sam = getSAM(inputFiles2)| map { x -> [ x.simpleName, x ] }
    fq = getFQ(inputFiles2)| map { x -> [ x.simpleName, x ] }
    bam = sam | map { x -> x[1] } | SAMToBAM | map { x -> [ x.simpleName, x ] }
    bed = gff | map { x -> x[1] } | gffToBED | map { x -> [ x.simpleName, x ] }

    // [[ nickname, gff, sam, fq, bam, bed ]]
    outputs = gff.join(sam).join(fq).join(bam).join(bed) | runTools

    outputs.filter( ~/.*\.gff$/ ) | convertToBed
}