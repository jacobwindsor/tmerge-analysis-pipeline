nextflow.preview.dsl=2

/*
* tmerge 2.0 testing pipeline
*/

/*
* Configuration options
*/
tmerge2_path = "/users/rg/jwindsor/tmerge/2.0"
tmerge1_path = "/users/rg/jlagarde/julien_utils"
julien_utils_path = "/users/rg/jlagarde/julien_utils/"
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
    zcat "!{x.path}" !{x.nickname}.input.gff
    '''    
}

process runTmerge1 {
    input:
    path inputGFF

    memory '30 GB'
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

process outputToBED {
    input:
    path mergedOutput

    publishDir "$output_dir"

    output:
    path '*.bed'

    shell:
    '''
    PATH="$PATH:!{julien_utils_path}"

    cat !{mergedOutput} | gff2bed_full.pl -| sortbed > !{mergedOutput.baseName}.bed
    '''
}

process oneVsTwo {
    input:
    path tmerge1Output
    path tmerge2Output

    publishDir "$output_dir"

    output:
    path '*.gffcompare*'


    shell:
    '''
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -o !{tmerge1Output.simpleName}.tmerge2_vs_tmerge1.gffcompare !{tmerge2Output} -r !{tmerge1Output}
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

workflow {
    inputs = copyInputFiles(inputFiles)
    
    tmerge1 = runTmerge1(inputs)
    tmerge2 = runTmerge2(inputs)
    recordTime(tmerge2.time)

    oneVsTwo(tmerge1.output, tmerge2.output)
    
    tmerge1.output.mix(tmerge2.output) | outputToBED
}