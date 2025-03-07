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

    #selphi2_kinase_radial_dendrogram.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/radial_all.pdf \
        selphi2_kinase_dendrogram/radial_all.tsv \
        all

    selphi2_kinase_radial_dendrogram_minphos.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/radial_all.pdf \
        selphi2_kinase_dendrogram/radial_all.tsv \
        all \
        ${params.minphos} \
        selphi2_kinase_dendrogram/distance_matrix_all.tsv \
        ${kinase_labels_csv}
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

    #selphi2_kinase_radial_dendrogram.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/radial_tyr.pdf \
        selphi2_kinase_dendrogram/radial_tyr.tsv \
        tyr

    selphi2_kinase_radial_dendrogram_minphos.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/radial_tyr.pdf \
        selphi2_kinase_dendrogram/radial_tyr.tsv \
        tyr \
        ${params.minphos} \
        selphi2_kinase_dendrogram/distance_matrix_tyr.tsv \
        ${kinase_labels_csv}
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

    #selphi2_kinase_radial_dendrogram.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/radial_ser_thr.pdf \
        selphi2_kinase_dendrogram/radial_ser_thr.tsv \
        ser/thr

    selphi2_kinase_radial_dendrogram_minphos.py \
        k_p_associations.tsv \
        selphi2_kinase_dendrogram/radial_ser_thr.pdf \
        selphi2_kinase_dendrogram/radial_ser_thr.tsv \
        ser/thr \
        ${params.minphos} \
        selphi2_kinase_dendrogram/distance_matrix_ser_thr.tsv \
        ${kinase_labels_csv}
    """

}


/*
[...]
*/
process go_terms_dendrogram_all {

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

    cat high_conf | awk -F"," '{print \$2"\\t"\$4"\\t"\$5}' > k_p_associations.tsv

    #selphi2_go_terms_radial_dendrogram.py \
        ${kinase_substrates_go_leaves_enrich} \
        selphi2_kinase_dendrogram/radial_terms_all.pdf \
        k_p_associations.tsv \
        all

    selphi2_go_terms_radial_dendrogram_minphos.py \
        ${kinase_substrates_go_leaves_enrich} \
        selphi2_kinase_dendrogram/radial_terms_all.pdf \
        k_p_associations.tsv \
        all \
        ${params.minphos} \
        selphi2_kinase_dendrogram/distance_matrix_terms_all.tsv \
        ${kinase_labels_csv}
    """

}


/*
[...]
*/
process go_terms_dendrogram_ser_thr {

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

    cat high_conf | awk -F"," '{print \$2"\\t"\$4"\\t"\$5}' > k_p_associations.tsv

    #selphi2_go_terms_radial_dendrogram.py \
        ${kinase_substrates_go_leaves_enrich} \
        selphi2_kinase_dendrogram/radial_terms_ser_thr.pdf \
        k_p_associations.tsv \
        ser/thr

    selphi2_go_terms_radial_dendrogram_minphos.py \
        ${kinase_substrates_go_leaves_enrich} \
        selphi2_kinase_dendrogram/radial_terms_ser_thr.pdf \
        k_p_associations.tsv \
        ser/thr \
        ${params.minphos} \
        selphi2_kinase_dendrogram/distance_matrix_terms_ser_thr.tsv \
        ${kinase_labels_csv}
    """

}


/*
[...]
*/
process go_terms_dendrogram_tyr {

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

    cat high_conf | awk -F"," '{print \$2"\\t"\$4"\\t"\$5}' > k_p_associations.tsv

    #selphi2_go_terms_radial_dendrogram.py \
        ${kinase_substrates_go_leaves_enrich} \
        selphi2_kinase_dendrogram/radial_terms_tyr.pdf \
        k_p_associations.tsv \
        tyr

    selphi2_go_terms_radial_dendrogram_minphos.py \
        ${kinase_substrates_go_leaves_enrich} \
        selphi2_kinase_dendrogram/radial_terms_tyr.pdf \
        k_p_associations.tsv \
        tyr \
        ${params.minphos} \
        selphi2_kinase_dendrogram/distance_matrix_terms_tyr.tsv \
        ${kinase_labels_csv}
    """

}


/*
[...]
*/
process selphi_eval_classifier_w_random_neg_set {

    memory '16G'

    input:
        tuple val(kin_fam),
              val(id),
              file('input/pos_set.tsv')

    output:
        tuple val(kin_fam),
              file("selphi2_100_rand_neg_sets/roc_points/${id}_roc_points.tsv"), emit: roc_points
        tuple val(kin_fam),
              file("selphi2_100_rand_neg_sets/pr_points/${id}_pr_points.tsv"), emit: pr_points
        tuple val(kin_fam),
              file("selphi2_100_rand_neg_sets/roc_points/${id}_roc_auc.txt"), emit: roc_auc
        tuple val(kin_fam),
              file("selphi2_100_rand_neg_sets/pr_points/${id}_pr_auc.txt"), emit: pr_auc

    script:
    """
    mkdir -p selphi2_100_rand_neg_sets/roc_points/
    mkdir -p selphi2_100_rand_neg_sets/pr_points/

    cat ${selphi2_prediction_matrix_dir}/prediction_matrix.csv  \
        | awk -F"," '\$5=="${kin_fam}"' \
        > prediction_matrix.csv

    cat input/pos_set.tsv \
        | awk '{print \$1"_"\$2}' \
        | sed '1d' \
        > pos_set_filter.txt

    cat prediction_matrix.csv | grep -f pos_set_filter.txt \
        | awk -F"," '{print \$2"_"\$6"_"\$9"\\t"\$NF}' \
        > pos_set.tsv

    n_pos=\$(cat pos_set.tsv | wc -l)
    n_neg=\$(echo \$((\${n_pos} * 1)))

    awk \
        -v n=\${n_neg} \
        -v seed=\$RANDOM \
        'BEGIN {srand(seed)} \
        {a[NR]=\$0} END {for (i=1; i<=n; i++) print a[int(rand()*NR)+1]}' \
        prediction_matrix.csv \
        | awk -F"," '{print \$2"_"\$6"_"\$9"\\t"\$NF}' \
        > neg_set.tsv

    cat pos_set.tsv neg_set.tsv \
        > data.tsv

    compute_roc_curve_points.py \
        ${id} \
        data.tsv \
        pos_set_filter.txt \
        selphi2_100_rand_neg_sets/roc_points/

    compute_pr_curve_points.py \
        ${id} \
        data.tsv \
        pos_set_filter.txt \
        selphi2_100_rand_neg_sets/pr_points/
    """

}


/*
given the points for multiple ROC curves, plot their mean, min, and max at each point
*/
process draw_roc_curves_per_kin_fam {

    publishDir "${out_dir}", pattern: "selphi2_100_rand_neg_sets/*.pdf", mode: 'copy'

    input:
        tuple val(kin_fam),
              file('input/*.tsv')
    
    output:
        path "selphi2_100_rand_neg_sets/${kin_fam}_roc_curves.pdf"
    
    script:
    """
    mkdir -p selphi2_100_rand_neg_sets

    ls input/ > curves.txt
    sed -i 's/^/input\\//' curves.txt

    draw_average_curve_from_points.py \
        FPR \
        TPR \
        curves.txt \
        selphi2_100_rand_neg_sets/${kin_fam}_roc_curves.pdf \
        ${kin_fam}
    """

}


/*
[...]
*/
process average_auroc_kin_fam {

    publishDir "${out_dir}", pattern: "selphi2_100_rand_neg_sets/*.txt", mode: 'copy'

    input:
        tuple val(kin_fam),
              file('input/*.txt')
    
    output:
        path "selphi2_100_rand_neg_sets/${kin_fam}_avg_auroc.txt"
    
    script:
    """
    mkdir -p selphi2_100_rand_neg_sets

    cat input/*txt > curves.txt

    average.py \
        curves.txt \
        selphi2_100_rand_neg_sets/${kin_fam}_avg_auroc.txt
    """

}


/*
[...]
*/
process average_aupr_kin_fam {

    publishDir "${out_dir}", pattern: "selphi2_100_rand_neg_sets/*.txt", mode: 'copy'

    input:
        tuple val(kin_fam),
              file('input/*.txt')
    
    output:
        path "selphi2_100_rand_neg_sets/${kin_fam}_avg_aupr.txt"
    
    script:
    """
    mkdir -p selphi2_100_rand_neg_sets

    cat input/*txt > curves.txt

    average.py \
        curves.txt \
        selphi2_100_rand_neg_sets/${kin_fam}_avg_aupr.txt
    """

}


/*
given the points for multiple PR curves, plot their mean, min, and max at each point
*/
process draw_pr_curves_per_kin_fam {

    publishDir "${out_dir}", pattern: "selphi2_100_rand_neg_sets/*.pdf", mode: 'copy'

    input:
        tuple val(kin_fam),
              file('input/*.tsv')

    output:
        path "selphi2_100_rand_neg_sets/${kin_fam}_pr_curves.pdf"
    
    script:
    """
    mkdir -p selphi2_100_rand_neg_sets

    ls input/ > curves.txt
    sed -i 's/^/input\\//' curves.txt

    draw_average_curve_from_points.py \
        Recall \
        Precision \
        curves.txt \
        selphi2_100_rand_neg_sets/${kin_fam}_pr_curves.pdf \
        ${kin_fam}
    """

}


/*
[...]
*/
process selphi_eval_classifier_w_sugiyama_random_neg_set {

    memory '16G'

    input:
        tuple val(kin_fam),
              val(id)

    output:
        tuple val(kin_fam),
              file("selphi2_sugiyama_100_rand_neg_sets/roc_points/${id}_roc_points.tsv"), emit: roc_points
        tuple val(kin_fam),
              file("selphi2_sugiyama_100_rand_neg_sets/pr_points/${id}_pr_points.tsv"), emit: pr_points
        tuple val(kin_fam),
              file("selphi2_sugiyama_100_rand_neg_sets/roc_points/${id}_roc_auc.txt"), emit: roc_auc
        tuple val(kin_fam),
              file("selphi2_sugiyama_100_rand_neg_sets/pr_points/${id}_pr_auc.txt"), emit: pr_auc

    script:
    """
    mkdir -p selphi2_sugiyama_100_rand_neg_sets/roc_points/
    mkdir -p selphi2_sugiyama_100_rand_neg_sets/pr_points/

    cat ${selphi2_prediction_matrix_dir}/prediction_matrix.csv  \
        | awk -F"," '\$5=="${kin_fam}"' \
        > prediction_matrix.csv


    cat prediction_matrix.csv \
        | awk -F"," '\$16=="TRUE"{print \$2"_"\$6"_"\$9"\\t"\$NF}' \
        > pos_set.tsv

    cat pos_set.tsv | cut -f1 > pos_set_filter.txt

    n_pos=\$(cat pos_set.tsv | wc -l)
    n_neg=\$(echo \$((\${n_pos} * 1)))

    awk \
        -v n=\${n_neg} \
        -v seed=\$RANDOM \
        'BEGIN {srand(seed)} \
        {a[NR]=\$0} END {for (i=1; i<=n; i++) print a[int(rand()*NR)+1]}' \
        prediction_matrix.csv \
        | awk -F"," '{print \$2"_"\$6"_"\$9"\\t"\$NF}' \
        > neg_set.tsv

    cat pos_set.tsv neg_set.tsv \
        > data.tsv

    compute_roc_curve_points.py \
        ${id} \
        data.tsv \
        pos_set_filter.txt \
        selphi2_sugiyama_100_rand_neg_sets/roc_points/

    compute_pr_curve_points.py \
        ${id} \
        data.tsv \
        pos_set_filter.txt \
        selphi2_sugiyama_100_rand_neg_sets/pr_points/
    """

}


/*
[...]
*/
process average_auroc_kin_fam_sugiyama {

    publishDir "${out_dir}", pattern: "selphi2_sugiyama_100_rand_neg_sets/*.txt", mode: 'copy'

    input:
        tuple val(kin_fam),
              file('input/*.txt')
    
    output:
        tuple val(kin_fam),
              file("selphi2_sugiyama_100_rand_neg_sets/${kin_fam}_avg_auroc.txt")
    
    script:
    """
    mkdir -p selphi2_sugiyama_100_rand_neg_sets

    cat input/*txt > curves.txt

    average.py \
        curves.txt \
        selphi2_sugiyama_100_rand_neg_sets/${kin_fam}_avg_auroc.txt
    """

}


/*
[...]
*/
process average_aupr_kin_fam_sugiyama {

    publishDir "${out_dir}", pattern: "selphi2_sugiyama_100_rand_neg_sets/*.txt", mode: 'copy'

    input:
        tuple val(kin_fam),
              file('input/*.txt')
    
    output:
        tuple val(kin_fam),
              file("selphi2_sugiyama_100_rand_neg_sets/${kin_fam}_avg_aupr.txt")
    
    script:
    """
    mkdir -p selphi2_sugiyama_100_rand_neg_sets

    cat input/*txt > curves.txt

    average.py \
        curves.txt \
        selphi2_sugiyama_100_rand_neg_sets/${kin_fam}_avg_aupr.txt
    """

}


/*
given the points for multiple ROC curves, plot their mean, min, and max at each point
*/
process draw_roc_curves_per_kin_fam_sugiyama {

    publishDir "${out_dir}", pattern: "selphi2_sugiyama_100_rand_neg_sets/*.pdf", mode: 'copy'

    input:
        tuple val(kin_fam),
              file('input/*.tsv')
    
    output:
        path "selphi2_sugiyama_100_rand_neg_sets/${kin_fam}_roc_curves.pdf"
    
    script:
    """
    mkdir -p selphi2_sugiyama_100_rand_neg_sets

    ls input/ > curves.txt
    sed -i 's/^/input\\//' curves.txt

    draw_average_curve_from_points.py \
        FPR \
        TPR \
        curves.txt \
        selphi2_sugiyama_100_rand_neg_sets/${kin_fam}_roc_curves.pdf \
        ${kin_fam}
    """

}


/*
given the points for multiple PR curves, plot their mean, min, and max at each point
*/
process draw_pr_curves_per_kin_fam_sugiyama {

    publishDir "${out_dir}", pattern: "selphi2_sugiyama_100_rand_neg_sets/*.pdf", mode: 'copy'

    input:
        tuple val(kin_fam),
              file('input/*.tsv')

    output:
        path "selphi2_sugiyama_100_rand_neg_sets/${kin_fam}_pr_curves.pdf"
    
    script:
    """
    mkdir -p selphi2_sugiyama_100_rand_neg_sets

    ls input/ > curves.txt
    sed -i 's/^/input\\//' curves.txt

    draw_average_curve_from_points.py \
        Recall \
        Precision \
        curves.txt \
        selphi2_sugiyama_100_rand_neg_sets/${kin_fam}_pr_curves.pdf \
        ${kin_fam}
    """

}


/*
[...]
*/
process selphi_eval_classifier_w_hijazi_random_neg_set {

    memory '16G'

    input:
        tuple val(kin_fam),
              val(id)

    output:
        tuple val(kin_fam),
              file("selphi2_hijazi_100_rand_neg_sets/roc_points/${id}_roc_points.tsv"), emit: roc_points
        tuple val(kin_fam),
              file("selphi2_hijazi_100_rand_neg_sets/pr_points/${id}_pr_points.tsv"), emit: pr_points
        tuple val(kin_fam),
              file("selphi2_hijazi_100_rand_neg_sets/roc_points/${id}_roc_auc.txt"), emit: roc_auc
        tuple val(kin_fam),
              file("selphi2_hijazi_100_rand_neg_sets/pr_points/${id}_pr_auc.txt"), emit: pr_auc

    script:
    """
    mkdir -p selphi2_hijazi_100_rand_neg_sets/roc_points/
    mkdir -p selphi2_hijazi_100_rand_neg_sets/pr_points/

    cat ${selphi2_prediction_matrix_dir}/prediction_matrix.csv  \
        | awk -F"," '\$5=="${kin_fam}"' \
        > prediction_matrix.csv


    cat prediction_matrix.csv \
        | awk -F"," '\$15=="TRUE"{print \$2"_"\$6"_"\$9"\\t"\$NF}' \
        > pos_set.tsv

    cat pos_set.tsv | cut -f1 > pos_set_filter.txt

    n_pos=\$(cat pos_set.tsv | wc -l)
    n_neg=\$(echo \$((\${n_pos} * 1)))

    awk \
        -v n=\${n_neg} \
        -v seed=\$RANDOM \
        'BEGIN {srand(seed)} \
        {a[NR]=\$0} END {for (i=1; i<=n; i++) print a[int(rand()*NR)+1]}' \
        prediction_matrix.csv \
        | awk -F"," '{print \$2"_"\$6"_"\$9"\\t"\$NF}' \
        > neg_set.tsv

    cat pos_set.tsv neg_set.tsv \
        > data.tsv

    compute_roc_curve_points.py \
        ${id} \
        data.tsv \
        pos_set_filter.txt \
        selphi2_hijazi_100_rand_neg_sets/roc_points/

    compute_pr_curve_points.py \
        ${id} \
        data.tsv \
        pos_set_filter.txt \
        selphi2_hijazi_100_rand_neg_sets/pr_points/
    """

}


/*
given the points for multiple ROC curves, plot their mean, min, and max at each point
*/
process draw_roc_curves_per_kin_fam_hijazi {

    publishDir "${out_dir}", pattern: "selphi2_hijazi_100_rand_neg_sets/*.pdf", mode: 'copy'

    input:
        tuple val(kin_fam),
              file('input/*.tsv')
    
    output:
        path "selphi2_hijazi_100_rand_neg_sets/${kin_fam}_roc_curves.pdf"
    
    script:
    """
    mkdir -p selphi2_hijazi_100_rand_neg_sets

    ls input/ > curves.txt
    sed -i 's/^/input\\//' curves.txt

    draw_average_curve_from_points.py \
        FPR \
        TPR \
        curves.txt \
        selphi2_hijazi_100_rand_neg_sets/${kin_fam}_roc_curves.pdf \
        ${kin_fam}
    """

}


/*
given the points for multiple PR curves, plot their mean, min, and max at each point
*/
process draw_pr_curves_per_kin_fam_hijazi {

    publishDir "${out_dir}", pattern: "selphi2_hijazi_100_rand_neg_sets/*.pdf", mode: 'copy'

    input:
        tuple val(kin_fam),
              file('input/*.tsv')

    output:
        path "selphi2_hijazi_100_rand_neg_sets/${kin_fam}_pr_curves.pdf"
    
    script:
    """
    mkdir -p selphi2_hijazi_100_rand_neg_sets

    ls input/ > curves.txt
    sed -i 's/^/input\\//' curves.txt

    draw_average_curve_from_points.py \
        Recall \
        Precision \
        curves.txt \
        selphi2_hijazi_100_rand_neg_sets/${kin_fam}_pr_curves.pdf \
        ${kin_fam}
    """

}


/*
[...]
*/
process selphi_psp_stats {

    publishDir "${out_dir}", pattern: "selphi2_100_rand_neg_sets/*.txt", mode: 'copy'
    publishDir "${out_dir}", pattern: "selphi2_100_rand_neg_sets/*.tsv", mode: 'copy'

    input:
        tuple val(kin_fam),
              file("input/pos_set.tsv")

    output:
        path "selphi2_100_rand_neg_sets/*.txt"
        path "selphi2_100_rand_neg_sets/*.tsv"

    script:
    """
    mkdir -p selphi2_100_rand_neg_sets/roc_points/
    mkdir -p selphi2_100_rand_neg_sets/pr_points/

    cat ${selphi2_prediction_matrix_dir}/prediction_matrix.csv  \
        | awk -F"," '\$5=="${kin_fam}"' \
        > prediction_matrix.csv

    cat input/pos_set.tsv \
        | awk '{print \$1"_"\$2}' \
        | sed '1d' \
        > pos_set_filter.txt

    cat prediction_matrix.csv | grep -f pos_set_filter.txt \
        | awk -F"," '{print \$2"_"\$6"_"\$9"\\t"\$NF}' \
        > selphi2_100_rand_neg_sets/pos_set_${kin_fam}.tsv

    n_pos=\$(cat selphi2_100_rand_neg_sets/pos_set_${kin_fam}.tsv | wc -l)

    echo "\${n_pos}" > selphi2_100_rand_neg_sets/n_${kin_fam}.txt
    """

}


/*
[...]
*/
process selphi_sugiyama_stats {

    publishDir "${out_dir}", pattern: "selphi2_sugiyama_100_rand_neg_sets/*.txt", mode: 'copy'
    publishDir "${out_dir}", pattern: "selphi2_sugiyama_100_rand_neg_sets/*.tsv", mode: 'copy'

    input:
        val kin_fam

    output:
        path "selphi2_sugiyama_100_rand_neg_sets/*.txt"
        path "selphi2_sugiyama_100_rand_neg_sets/*.tsv"

    script:
    """
    mkdir -p selphi2_sugiyama_100_rand_neg_sets/roc_points/
    mkdir -p selphi2_sugiyama_100_rand_neg_sets/pr_points/

    cat ${selphi2_prediction_matrix_dir}/prediction_matrix.csv  \
        | awk -F"," '\$5=="${kin_fam}"' \
        > prediction_matrix.csv

    cat prediction_matrix.csv \
        | awk -F"," '\$16=="TRUE"{print \$2"_"\$6"_"\$9"\\t"\$NF}' \
        > selphi2_sugiyama_100_rand_neg_sets/pos_set_${kin_fam}.tsv

    n_pos=\$(cat selphi2_sugiyama_100_rand_neg_sets/pos_set_${kin_fam}.tsv | wc -l)

    echo "\${n_pos}" > selphi2_sugiyama_100_rand_neg_sets/n_${kin_fam}.txt
    """

}


/*
[...]
*/
process selphi_hijazi_stats {

    publishDir "${out_dir}", pattern: "selphi2_hijazi_100_rand_neg_sets/*.txt", mode: 'copy'
    publishDir "${out_dir}", pattern: "selphi2_hijazi_100_rand_neg_sets/*.tsv", mode: 'copy'

    input:
        val kin_fam

    output:
        path "selphi2_hijazi_100_rand_neg_sets/*.txt"
        path "selphi2_hijazi_100_rand_neg_sets/*.tsv"

    script:
    """
    mkdir -p selphi2_hijazi_100_rand_neg_sets/roc_points/
    mkdir -p selphi2_hijazi_100_rand_neg_sets/pr_points/

    cat ${selphi2_prediction_matrix_dir}/prediction_matrix.csv  \
        | awk -F"," '\$5=="${kin_fam}"' \
        > prediction_matrix.csv

    cat prediction_matrix.csv \
        | awk -F"," '\$15=="TRUE"{print \$2"_"\$6"_"\$9"\\t"\$NF}' \
        > selphi2_hijazi_100_rand_neg_sets/pos_set_${kin_fam}.tsv

    n_pos=\$(cat selphi2_hijazi_100_rand_neg_sets/pos_set_${kin_fam}.tsv | wc -l)

    echo "\${n_pos}" > selphi2_hijazi_100_rand_neg_sets/n_${kin_fam}.txt
    """

}