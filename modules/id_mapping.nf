#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
.META:
1. UniProtKB-AC 
2. Gene Name
*/
process dl_human_ref_proteome_ids {

    publishDir "${out_dir}",
                pattern: 'uniprot/human_ref_proteome_ids.tsv',
                mode: 'copy'

    output:
        path 'uniprot/human_ref_proteome_ids.tsv'

    shell:
    """
    mkdir -p uniprot

    wget \
        -O uniprot/human_ref_proteome_ids.tsv \
        "${params.url_uniprot_reference_proteome_tsv}"
    """
}

/*
.META:
1. Gene Name
2. UniProtKB-AC
*/
process gene_name2uniprot_ac_dict {

    publishDir "${out_dir}",
                pattern: 'uniprot/gene_name2uniprot_ac_dict.tsv',
                mode: 'copy'

    input:
        path 'input/human_ref_proteome_ids.tsv'

    output:
        path 'uniprot/gene_name2uniprot_ac_dict.tsv'

    shell:
    """
    mkdir -p uniprot

    cat input/human_ref_proteome_ids.tsv \
        | sed '1d' \
        | awk '{print \$2"\\t"\$1}' \
        | sort \
        | uniq \
        >  uniprot/gene_name2uniprot_ac_dict.tsv
    """
}

/*
.META:
1. UniProtKB-AC
2. Gene Name
*/
process uniprot_ac2gene_name_dict {

    publishDir "${out_dir}",
                pattern: 'uniprot/uniprot_ac2gene_name_dict.tsv',
                mode: 'copy'

    input:
        path 'input/HUMAN_9606_idmapping.dat.gz'

    output:
        path 'uniprot/uniprot_ac2gene_name_dict.tsv'

    shell:
    """
    mkdir -p uniprot

    zcat input/HUMAN_9606_idmapping.dat.gz \
        | awk '\$2=="Gene_Name"{print \$1"\\t"\$3}' \
        > uniprot/uniprot_ac2gene_name_dict.tsv
    """
}

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
.META:
1. UniProtKB-AC 
2. ID_type 
3. ID
*/
process dl_human {

    publishDir "${out_dir}",
                pattern: 'uniprot/HUMAN_9606_idmapping.dat.gz',
                mode: 'copy'

    output:
        path 'uniprot/HUMAN_9606_idmapping.dat.gz'

    shell:
    """
    mkdir -p uniprot

    wget -P uniprot ${params.human_url}
    """

}

/*
3-fields tab-delimited file

provides mapping in this order:

IDa --->  UniProt AC ---> IDb

all ids of nomenclature IDa are mapped to the corresponding UniProt ACs,
and the UniProt ACs are then mapped to the corresponding ids of nomenclature
IDb. For example, all ENSP ids are mapped to their UniProt AC. Then, each of
the UniProt AC identifiers is mapped to HGNC ids.

.META:
1. IDa
2. UniProt AC referred to IDa
3. IDb referred to UniProt AC
*/
process IDa2uniprot2IDb {

    publishDir "${out_dir}",
                pattern: "uniprot/${IDa}2uniprot2${IDb}.tsv",
                mode: 'copy'

    input:
        path 'input/mapping.tsv.gz'
        val IDa
        val IDb

    output:
        path "uniprot/${IDa}2uniprot2${IDb}.tsv"

    script:
    """
    mkdir -p uniprot

    zcat input/mapping.tsv.gz \
        | grep -w -f <(echo -e "${IDa}\n${IDb}") \
        | gzip > mapping.tsv.gz

    IDa2uniprot2IDb.py \
        mapping.tsv.gz \
        ${IDa} \
        ${IDb} \
        > uniprot/${IDa}2uniprot2${IDb}.tsv
    """

}