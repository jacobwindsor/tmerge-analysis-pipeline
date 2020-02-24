output_dir = "/users/rg/jwindsor/tests/tmerge/results"

process runGFFCompare {
    input:
    path input
    path reference

    publishDir "$output_dir"

    output:
    path '*.gffcompare*'

    shell:
    '''
    gffcompare --strict-match --no-merge -e 0 -d 0 --debug -o !{input.simpleName}.!{input.baseName.split("\\\\.")[2]}_vs_!{reference.baseName.split("\\\\.")[2]}.gffcompare !{input} -r !{reference}
    '''
}