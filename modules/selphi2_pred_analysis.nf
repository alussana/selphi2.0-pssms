#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
[...]
*/
process kinase_dendrogram_all {

    publishDir "${out_dir}", 
               pattern: "selphi2_kinase_dendrogram/*tsv",
               mode: 'copy'
    publishDir "${out_dir}", 
               pattern: "selphi2_kinase_dendrogram/*pdf",
               mode: 'copy'

    output:
        path 'selphi2_kinase_dendrogram/*.tsv'
        path 'selphi2_kinase_dendrogram/*.pdf'

    script:
    """
    mkdir -p selphi2_kinase_dendrogram

    ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_high_conf     high_conf

    cat high_conf | awk -F"," '{print \$2"\\t"\$4"\\t"\$5"\\t"\$7"_"\$11}' > k_p_associations.tsv

    #selphi2_kinase_dendrogram.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/all.pdf \
        selphi2_kinase_dendrogram/all.tsv \
        all

    selphi2_kinase_radial_dendrogram.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/radial_all.pdf \
        selphi2_kinase_dendrogram/radial_all.tsv \
        all
    """

}


/*
[...]
*/
process kinase_dendrogram_tyr {

    publishDir "${out_dir}", 
               pattern: "selphi2_kinase_dendrogram/*tsv",
               mode: 'copy'
    publishDir "${out_dir}", 
               pattern: "selphi2_kinase_dendrogram/*pdf",
               mode: 'copy'

    output:
        path 'selphi2_kinase_dendrogram/*.tsv'
        path 'selphi2_kinase_dendrogram/*.pdf'

    script:
    """
    mkdir -p selphi2_kinase_dendrogram

    ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_high_conf     high_conf

    cat high_conf | awk -F"," '{print \$2"\\t"\$4"\\t"\$5"\\t"\$7"_"\$11}' > k_p_associations.tsv

    #selphi2_kinase_dendrogram.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/tyr.pdf \
        selphi2_kinase_dendrogram/tyr.tsv \
        tyr

    selphi2_kinase_radial_dendrogram.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/radial_tyr.pdf \
        selphi2_kinase_dendrogram/radial_tyr.tsv \
        tyr
    """

}


/*
[...]
*/
process kinase_dendrogram_ser_thr {

    publishDir "${out_dir}", 
               pattern: "selphi2_kinase_dendrogram/*tsv",
               mode: 'copy'
    publishDir "${out_dir}", 
               pattern: "selphi2_kinase_dendrogram/*pdf",
               mode: 'copy'

    output:
        path 'selphi2_kinase_dendrogram/*.tsv'
        path 'selphi2_kinase_dendrogram/*.pdf'

    script:
    """
    mkdir -p selphi2_kinase_dendrogram

    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_AGC.csv       AGC.csv
    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_Atypical.csv  Atypical.csv
    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_CAMK.csv      CAMK.csv
    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_CK1.csv       CK1.csv
    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_CMGC.csv      CMGC.csv
    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix.csv           prediction_matrix.csv
    ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_high_conf     high_conf
    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_NA.csv        NA.csv
    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_Other.csv     Other.csv
    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_overlap       overlap
    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_STE.csv       STE.csv
    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_TK.csv        TK.csv
    #ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_TKL.csv       TKL.csv

    cat high_conf | awk -F"," '{print \$2"\\t"\$4"\\t"\$5"\\t"\$7"_"\$11}' > k_p_associations.tsv

    #selphi2_kinase_dendrogram.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/ser_thr.pdf \
        selphi2_kinase_dendrogram/ser_thr.tsv \
        ser/thr

    selphi2_kinase_radial_dendrogram.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/radial_ser_thr.pdf \
        selphi2_kinase_dendrogram/radial_ser_thr.tsv \
        ser/thr
    """

}


/*
[...]
*/
process go_terms_dendrogram_all {

    publishDir "${out_dir}", 
               pattern: "selphi2_kinase_dendrogram/*pdf",
               mode: 'copy'

    output:
        path 'selphi2_kinase_dendrogram/*.pdf'

    script:
    """
    mkdir -p selphi2_kinase_dendrogram

    ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_high_conf     high_conf

    cat high_conf | awk -F"," '{print \$2"\\t"\$4"\\t"\$5"\\t"\$7"_"\$11}' > k_p_associations.tsv

    selphi2_go_terms_radial_dendrogram.py \
        ${kinase_substrates_go_leaves_enrich} \
        selphi2_kinase_dendrogram/radial_terms_all.pdf \
        k_p_associations.tsv \
        all
    """

}


/*
[...]
*/
process go_terms_dendrogram_ser_thr {

    publishDir "${out_dir}", 
               pattern: "selphi2_kinase_dendrogram/*pdf",
               mode: 'copy'

    output:
        path 'selphi2_kinase_dendrogram/*.pdf'

    script:
    """
    mkdir -p selphi2_kinase_dendrogram

    ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_high_conf     high_conf

    cat high_conf | awk -F"," '{print \$2"\\t"\$4"\\t"\$5"\\t"\$7"_"\$11}' > k_p_associations.tsv

    selphi2_go_terms_radial_dendrogram.py \
        ${kinase_substrates_go_leaves_enrich} \
        selphi2_kinase_dendrogram/radial_terms_ser_thr.pdf \
        k_p_associations.tsv \
        ser/thr
    """

}


/*
[...]
*/
process go_terms_dendrogram_tyr {

    publishDir "${out_dir}", 
               pattern: "selphi2_kinase_dendrogram/*pdf",
               mode: 'copy'

    output:
        path 'selphi2_kinase_dendrogram/*.pdf'

    script:
    """
    mkdir -p selphi2_kinase_dendrogram

    ln -s ${selphi2_prediction_matrix_dir}/prediction_matrix_high_conf     high_conf

    cat high_conf | awk -F"," '{print \$2"\\t"\$4"\\t"\$5"\\t"\$7"_"\$11}' > k_p_associations.tsv

    selphi2_go_terms_radial_dendrogram.py \
        ${kinase_substrates_go_leaves_enrich} \
        selphi2_kinase_dendrogram/radial_terms_tyr.pdf \
        k_p_associations.tsv \
        tyr
    """

}
