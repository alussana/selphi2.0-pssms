#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
[...]
*/
process get_p_sequences {

    label 'phosformer'

    publishDir "${out_dir}", pattern: "feature_table/k_p_seq.tsv", mode: 'copy'

    input:
        path 'input/kinases_table.tsv'
        path 'input/k_p_table.tsv'

    output:
        path 'feature_table/k_p_seq.tsv'


    script:
    """
    mkdir -p feature_table

    cat input/k_p_table.tsv \
        | awk '{print \$1"_"\$2"_"\$4"\t"\$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5}' \
        > k_p_table.tsv

    eval_phosformer_prepare_dataset.py \
        input/kinases_table.tsv \
        k_p_table.tsv \
        > feature_table/k_p_seq.tsv
    """

}


/*
[...]
*/
process get_s_t_phosformer_kinases {

    label 'phosformer'

    publishDir "${out_dir}", pattern: "phosformer_model_100_rand_neg_sets/*tsv", mode: 'copy'

    output:
        path 'phosformer_model_100_rand_neg_sets/s_t_phosformer_kinases.tsv'


    script:
    """
    mkdir -p phosformer_model_100_rand_neg_sets/

    get_phosformer_kinases.py \
        ${projectDir}/src/phosformer \
        S_T \
        > phosformer_model_100_rand_neg_sets/s_t_phosformer_kinases.tsv
    """

}


/*
[...]
*/
process get_y_phosformer_kinases {

    label 'phosformer'

    publishDir "${out_dir}", pattern: "phosformer_model_100_rand_neg_sets/*tsv", mode: 'copy'

    output:
        path 'phosformer_model_100_rand_neg_sets/y_phosformer_kinases.tsv'


    script:
    """
    mkdir -p phosformer_model_100_rand_neg_sets/

    get_phosformer_kinases.py \
        ${projectDir}/src/phosformer \
        Y \
        > phosformer_model_100_rand_neg_sets/y_phosformer_kinases.tsv
    """

}


/*
[...]
*/
process get_1col {

    input:
        path 'input/table.tsv'

    output:
        path 'phosformer_model_100_rand_neg_sets/phosformer_kinases.txt'

    script:
    """
    mkdir -p phosformer_model_100_rand_neg_sets/

    cat input/table.tsv \
        | cut -f1 | sort | uniq | awk 'NF' \
        > phosformer_model_100_rand_neg_sets/phosformer_kinases.txt
    """

}


/*
[...]
*/
process translate_k_p_pos_set {

    input:
        path 'input/file.tsv'
        path 'input/dict.tsv'

    output:
        path 'translated_file.tsv'

    script:
    """
    cat input/file.tsv | tr '_' '\\t' > file.tsv

    translator.py \
        input/dict.tsv \
        file.tsv \
        1 \
        3 \
        1,2 \
        1 \
    | awk '{print \$1"\\t"\$2"\\t"\$3}' \
        > translated_file.tsv
    """

}


/*
Filter positive set to only include given kinases
*/
process filter_set_for_phosformer_kinases {

    label 'phosformer'

    input:
        path 'input/kinases.txt'
        path 'input/pos_set.tsv'

    output:
        path 'phosformer_model_100_rand_neg_sets/pos_set.tsv'


    script:
    """
    mkdir -p phosformer_model_100_rand_neg_sets/

    cat input/pos_set.tsv \
        | filter1col.py input/kinases.txt 1 \
        > phosformer_model_100_rand_neg_sets/pos_set.tsv
    """

}


/*
Filter dataset to only include kinases that have a phosformer
*/
process filter_k_p_for_phosformer_kinases {

    input:
        path 'input/phosformer_kinases.txt'
        path 'input/k_p_phosformer_scores.tsv'

    output:
        path 'phosformer_model_100_rand_neg_sets/k_p_phosformer_kinases.tsv'


    script:
    """
    mkdir -p phosformer_model_100_rand_neg_sets/

    cat input/k_p_phosformer_scores.tsv \
        | filter1col.py input/phosformer_kinases.txt 1 \
        > phosformer_model_100_rand_neg_sets/k_p_phosformer_kinases.tsv
    """

}


/*
Sample random kinase-phosphosite pairs to build a negative set of examples that
is equally sized as the positive set
Use such dataset to compute ROC curve points and PR curve points based on the
phosformer score predictor
Curves are interpolated using np.interp to report 200 representative points
*/
process eval_phosformer_w_random_neg_set {

    label 'phosformer'

    memory '16G'

    input:
        tuple val(id),
              file('input/pos_set.tsv'),
              file('input/k_p.tsv'),
              file('input/kinases_table.tsv')

    output:
        path "phosformer_model_100_rand_neg_sets/roc_points/${id}_roc_points.tsv", emit: roc_points
        path "phosformer_model_100_rand_neg_sets/roc_points/${id}_roc_auc.txt"
        path "phosformer_model_100_rand_neg_sets/pr_points/${id}_pr_points.tsv", emit: pr_points
        path "phosformer_model_100_rand_neg_sets/pr_points/${id}_pr_auc.txt"
        path "phosformer_model_100_rand_neg_sets/k_p_eval_score.tsv"


    script:
    """
    mkdir -p phosformer_model_100_rand_neg_sets/roc_points/
    mkdir -p phosformer_model_100_rand_neg_sets/pr_points/

    cat input/pos_set.tsv \
        | awk 'NF' \
        | awk '{print \$1"_"\$2"_"\$3}' \
        > pos_set.txt

    cat input/k_p.tsv \
        | awk -F "," '{print \$1"_"\$2"_"\$4"\\t"\$1"\\t"\$2"\\t"\$5"\\t"\$4"\\t"\$3}' \
        > k_p.tsv

    cat k_p.tsv \
        | grep -f pos_set.txt \
        > k_p_pos_set.tsv

    eval_phosformer_make_dataset_and_predict.py \
        input/kinases_table.tsv \
        k_p_pos_set.tsv \
        k_p.tsv \
        1 \
        ${projectDir}/src/phosformer \
        > phosformer_model_100_rand_neg_sets/k_p_eval_score.tsv

    phosformer_compute_roc_curve_points.py \
        ${id} \
        phosformer_model_100_rand_neg_sets/k_p_eval_score.tsv \
        phosformer_model_100_rand_neg_sets/roc_points/

    phosformer_compute_pr_curve_points.py \
        ${id} \
        phosformer_model_100_rand_neg_sets/k_p_eval_score.tsv \
        phosformer_model_100_rand_neg_sets/pr_points/
    """

}

/*
given the points for multiple ROC curves, plot their mean, min, and max at each point
*/
process draw_roc_curves_phosformer {

    publishDir "${out_dir}", pattern: "phosformer_model_100_rand_neg_sets/*.pdf", mode: 'copy'

    input:
        path 'input/*.tsv'
        val id

    output:
        path "phosformer_model_100_rand_neg_sets/${id}_roc_curves.pdf"
    
    script:
    """
    mkdir -p phosformer_model_100_rand_neg_sets

    ls input/ > curves.txt
    sed -i 's/^/input\\//' curves.txt

    draw_average_curve_from_points.py \
        FPR \
        TPR \
        curves.txt \
        phosformer_model_100_rand_neg_sets/${id}_roc_curves.pdf
    """

}

/*
given the points for multiple PR curves, plot their mean, min, and max at each point
*/
process draw_pr_curves_phosformer {

    publishDir "${out_dir}", pattern: "phosformer_model_100_rand_neg_sets/*.pdf", mode: 'copy'

    input:
        path 'input/*.tsv'
        val id

    output:
        path "phosformer_model_100_rand_neg_sets/${id}_pr_curves.pdf"
    
    script:
    """
    mkdir -p phosformer_model_100_rand_neg_sets

    ls input/ > curves.txt
    sed -i 's/^/input\\//' curves.txt

    draw_average_curve_from_points.py \
        Recall \
        Precision \
        curves.txt \
        phosformer_model_100_rand_neg_sets/${id}_pr_curves.pdf
    """

}