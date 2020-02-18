/*
* tmerge 2.0 testing pipeline
*/

/*
* Configuration options
*/
tmerge2_path = "/users/rg/jwindsor/tmerge/2.0"
tmerge1_path = "/users/rg/jwindsor/tmerge"
julien_utils_path = "/users/rg/jlagarde/julien_utils/"
venv = "/users/rg/jwindsor/venvs/tmerge2/bin/activate"
output_dir = "/users/rg/jwindsor/tests/tmerge/results"
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
    val x from inputFiles

    storeDir "$cache_dir"

    output:
    path '*.input.gff' into inputGFF

    shell:
    '''
    zcat "!{x.path}" > !{x.nickname}.input.gff
    '''
}

inputGFF.into { tmerge1Input; tmerge2Input }

process runTmerge2 {
    input:
    path tmerge2Input

    output:
    stdout tmerge2Speed
    path '*.tmerge2.output.gff' into tmerge2Output

    shell:
    '''
    module load Python
    source !{venv}

    time python !{tmerge2_path}/tmerge.py -i !{tmerge2Input} -o !{tmerge2Input.simpleName}.tmerge2.output.gff
    '''
}

process runTmerge1 {
    input:
    path tmerge1Input

    output:
    stdout tmerge1Speed
    path '*.tmerge1.output.gff' into tmerge1Output

    shell:
    '''
    time perl !{tmerge1_path}/tmerge !{tmerge1Input} > !{tmerge1Input.simpleName}.tmerge1.output.gff
    '''
}

tmerge1Output.into { tmerge1ToBEDInput; tmerge1OneVsTwoInput }
tmerge2Output.into { tmerge2ToBEDInput; tmerge2OneVsTwoInput }
mergedOutput = tmerge2ToBEDInput.mix(tmerge1ToBEDInput)

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

process OneVsTwo {
    input:
    path tmerge1OneVsTwoInput
    path tmerge2OneVsTwoInput

    publishDir "$output_dir"

    output:
    path '*.gffcompare*' into oneVsTwoOutput


    shell:
    '''
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -o !{tmerge1OneVsTwoInput.simpleName}.tmerge2_vs_tmerge1.gffcompare !{tmerge2OneVsTwoInput} -r !{tmerge1OneVsTwoInput}
    '''
}


tmerge2Speed.subscribe { x -> print x }
tmerge1Speed.subscribe { x -> print x }