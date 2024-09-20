#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Get kinases having a PSSM
*/
process get_kinases_having_pssm {

    input:
        path 'input/pssm_dict.h5'

    output:
        path 'pssm_model_100_rand_neg_sets/kinases_having_pssm.txt'


    script:
    """
    mkdir -p pssm_model_100_rand_neg_sets/

    print_kinases_having_pssm.py \
        input/pssm_dict.h5 \
        > pssm_model_100_rand_neg_sets/kinases_having_pssm.txt
    """

}


/*
Filter positive set to only include kinases that have a PSSM
*/
process filter_set_for_kinases_having_pssm {

    input:
        path 'input/kinases_having_pssm.txt'
        path 'input/pos_set.tsv'

    output:
        path 'pssm_model_100_rand_neg_sets/pos_set.tsv'


    script:
    """
    mkdir -p pssm_model_100_rand_neg_sets/

    cat input/pos_set.tsv \
        | filter1col.py input/kinases_having_pssm.txt 1 \
        > pssm_model_100_rand_neg_sets/pos_set.tsv
    """

}

/*
Filter dataset to only include kinases that have a PSSM
*/
process filter_scores_for_kinases_having_pssm {

    input:
        path 'input/kinases_having_pssm.txt'
        path 'input/k_p_pssm_scores.tsv'

    output:
        path 'pssm_model_100_rand_neg_sets/k_p_pssm_scores.tsv'


    script:
    """
    mkdir -p pssm_model_100_rand_neg_sets/

    cat input/k_p_pssm_scores.tsv \
        | filter1col.py input/kinases_having_pssm.txt 1 \
        > pssm_model_100_rand_neg_sets/k_p_pssm_scores.tsv
    """

}

/*
Sample random kinase-phosphosite pairs to build a negative set of examples that
is equally sized as the positive set
Use such dataset to compute ROC curve points and PR curve points based on the
PSSM score predictor
Curves are interpolated using np.interp to report 200 representative points
*/
process eval_classifier_w_random_neg_set {

    memory '16G'

    input:
        tuple val(id),
              file('input/pos_set.tsv'),
              file('input/k_p_pssm_scores.tsv')

    output:
        path "pssm_model_100_rand_neg_sets/roc_points/${id}_roc_points.tsv", emit: roc_points
        path "pssm_model_100_rand_neg_sets/roc_points/${id}_roc_auc.txt"
        path "pssm_model_100_rand_neg_sets/pr_points/${id}_pr_points.tsv", emit: pr_points
        path "pssm_model_100_rand_neg_sets/pr_points/${id}_pr_auc.txt"


    script:
    """
    mkdir -p pssm_model_100_rand_neg_sets/roc_points/
    mkdir -p pssm_model_100_rand_neg_sets/pr_points/

    cat input/pos_set.tsv \
        | awk 'NF' \
        | awk '{print \$1"_"\$2}' \
        > pos_set.txt

    cat input/k_p_pssm_scores.tsv \
        | awk '{print \$1"_"\$2"_"\$4"\\t"\$7}' \
        > k_p_pssm_scores.tsv

    cat k_p_pssm_scores.tsv \
        | grep -f pos_set.txt \
        > k_p_pssm_scores_pos_set.tsv

    n_pos=\$(cat k_p_pssm_scores_pos_set.tsv | wc -l)
    n_neg=\$(echo \$((\${n_pos} * 1)))
    #n_neg=\${n_pos}

    awk \
        -v n=\${n_neg} \
        -v seed=\$RANDOM \
        'BEGIN {srand(seed)} \
        {a[NR]=\$0} END {for (i=1; i<=n; i++) print a[int(rand()*NR)+1]}' \
        input/k_p_pssm_scores.tsv \
        | awk '{print \$1"_"\$2"_"\$4"\\t"\$7}' \
        > k_p_pssm_scores_neg_set.tsv

    cat k_p_pssm_scores_pos_set.tsv k_p_pssm_scores_neg_set.tsv \
        > k_p_pssm_scores.tsv

    compute_roc_curve_points.py \
        ${id} \
        k_p_pssm_scores.tsv \
        pos_set.txt \
        pssm_model_100_rand_neg_sets/roc_points/

    compute_pr_curve_points.py \
        ${id} \
        k_p_pssm_scores.tsv \
        pos_set.txt \
        pssm_model_100_rand_neg_sets/pr_points/
    """

}

/*
given the points for multiple ROC curves, plot their mean, min, and max at each point
*/
process draw_roc_curves {

    publishDir "${out_dir}", pattern: "pssm_model_100_rand_neg_sets/*.pdf", mode: 'copy'

    input:
        path 'input/*.tsv'
        val id

    output:
        path "pssm_model_100_rand_neg_sets/${id}_roc_curves.pdf"
    
    script:
    """
    mkdir -p pssm_model_100_rand_neg_sets

    ls input/ > curves.txt
    sed -i 's/^/input\\//' curves.txt

    draw_average_curve_from_points.py \
        FPR \
        TPR \
        curves.txt \
        pssm_model_100_rand_neg_sets/${id}_roc_curves.pdf
    """

}

/*
given the points for multiple PR curves, plot their mean, min, and max at each point
*/
process draw_pr_curves {

    publishDir "${out_dir}", pattern: "pssm_model_100_rand_neg_sets/*.pdf", mode: 'copy'

    input:
        path 'input/*.tsv'
        val id

    output:
        path "pssm_model_100_rand_neg_sets/${id}_pr_curves.pdf"
    
    script:
    """
    mkdir -p pssm_model_100_rand_neg_sets

    ls input/ > curves.txt
    sed -i 's/^/input\\//' curves.txt

    draw_average_curve_from_points.py \
        Recall \
        Precision \
        curves.txt \
        pssm_model_100_rand_neg_sets/${id}_pr_curves.pdf
    """

}