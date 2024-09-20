#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
translate all the words in the first field of input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the second column

discard untranslated rows
*/
process translate1col {

    input:
        path 'input/file.tsv'
        path 'input/dict.tsv'

    output:
        path 'translated_file.tsv'

    script:
    """
    translator.py \
        input/dict.tsv \
        input/file.tsv \
        1 \
        2 \
        1 \
        0 \
        > translated_file.tsv
    """

}

/*
translate all the words in the 3rd field of input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the second column

do not discard untranslated rows
*/
process translate3col_keep_untranslated {

    memory '64G'
    time '24h'

    input:
        path 'input/file.tsv'
        path 'input/dict.tsv'

    output:
        path 'translated_file.tsv'

    script:
    """
    translator.py \
        input/dict.tsv \
        input/file.tsv \
        1 \
        2 \
        3 \
        1 \
        > translated_file.tsv
    """

}

/*
translate all the words in the first field of input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the second column

don't discard untranslated rows
*/
process translate1col_keep_untranslated {

    input:
        path 'input/file.tsv'
        path 'input/dict.tsv'

    output:
        path 'translated_file.tsv'

    script:
    """
    translator.py \
        input/dict.tsv \
        input/file.tsv \
        1 \
        2 \
        1 \
        1 \
        > translated_file.tsv
    """

}

/*
publish a file giving it an arbitrary path
*/
process phony {

    publishDir "${out_dir}",
                pattern: "${outputFile}",
                mode: 'copy'

    input:
        path 'input/file'
        val outputFile

    output:
        path "${outputFile}"

    script:
    """
    mkdir -p \$(dirname ${outputFile})
    
    cp input/file ${outputFile}
    """

}

/*
split a single input file into chunks of n rows,
each named with a unique identifier

typically, each chunk in the output channel is mapped to its id with
.flatten().map{ file -> tuple( file.baseName, file ) }
*/
process split {

    input:
        path 'input/file'
        val n

    output:
        path "output/*"

    script:
    """
    mkdir -p output

    split -l ${n} -a 16 -x input/file output/
    """

}

/*
concatenate input files removing empty lines
*/
process concatenate {

    input:
        path 'input/*'
    
    output:
        path 'output/file'
    
    script:
    """
    mkdir -p output

    cat input/* | awk 'NF' > output/file
    """

}

/*
sort and uniq a file, except the first row
*/
process sort_uniq_w_header {

    input:
        path 'input/file'
    
    output:
        path 'output/file'
    
    script:
    """
    mkdir -p output

    cat input/file | sed -n '1p' > header

    cat input/file | sed '1d' | sort | uniq | awk 'NF' > body

    cat header body > output/file
    """

}