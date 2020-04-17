params.output_dir = "/users/rg/jwindsor/tests/tmerge/results"
params.julien_utils_path = "/users/rg/jlagarde/julien_utils/"


process runGFFCompare {
    input:
    path input
    path reference

    publishDir "$params.output_dir"

    output:
    path '*.gffcompare*'

    shell:
    '''
    gffcompare --strict-match --no-merge -T -e 0 -d 0 --debug -o !{input.simpleName}.!{input.baseName.split("\\\\.")[2]}_vs_!{reference.baseName.split("\\\\.")[2]}.gffcompare !{input} -r !{reference}
    '''
}

process processForSIRVs {
    input:
    path inputGFF

    publishDir "$params.output_dir"

    output:
    path '*.gff'

    shell:
    '''
    cat !{inputGFF} | awk -F"\t" '$1=="SIRVome_isoforms"' > !{inputGFF.baseName}.sirv.gff
    '''
}

process gffToBED {
    input:
    path inputGFF

    publishDir "$params.output_dir"

    output:
    path '*.bed'

    shell:
    '''
    PATH="$PATH:!{params.julien_utils_path}"

    cat !{inputGFF} | gff2bed_full.pl -| sortbed > !{inputGFF.baseName}.bed
    '''
}