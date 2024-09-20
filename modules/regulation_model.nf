#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
extract phosphosite only-specific features, make unique phosphosite entries
*/
process make_all_phosphosites_features {

    input:
        tuple val(id),
              file("input/${id}.tsv")

    output:
        path "output/${id}.tsv"

    script:
    """
    mkdir -p output

    cat input/${id}.tsv \
        | sed 's/ /\\t/g' \
        | cut -f2,3,6,15,16,17,18,19,20,21,22,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42,43 \
        | uniq \
        > output/${id}.tsv
    """

}

/*
extract phosphosite only-specific features and labels, make unique phosphosite entries

label has column name "status"; 0 is [...], 1 is [...]
*/
process make_regulation_train_dataset {

    input:
        path "input/regulation_features_table.tsv"

    output:
        path "phosphosite_regulation_sign/regulation_train_dataset.tsv"

    script:
    """
    mkdir -p phosphosite_regulation_sign

    cat input/regulation_features_table.tsv \
        | sed 's/ /_/g' \
        | cut -f2,3,4,7,16,17,18,19,20,21,22,23,26,27,28,29,30,32,33,34,35,36,37,38,39,40,41,42,43,44 \
        > phosphosite_regulation_sign/regulation_train_dataset.tsv
    """

}

/*
train a random forest classifier of phosphosite regulation sign
perform a grid search and run a 10-fold cross-validation
*/
process train_regulation_model {

    publishDir "${out_dir}", pattern: "phosphosite_regulation_sign/*", mode: 'copy'

    input:
        path 'input/regulation_train_dataset.tsv'

    output:
        path 'phosphosite_regulation_sign/model.pkl', emit: model
        path 'phosphosite_regulation_sign/*pdf'

    script:
    """
    mkdir -p phosphosite_regulation_sign

    #train_regulation_model.py \
    #    input/regulation_train_dataset.tsv \
    #    ${params.n_jobs}
    #    phosphosite_regulation_sign

    touch phosphosite_regulation_sign/model.pkl
    touch phosphosite_regulation_sign/tmp.pdf
    """

}

/*
predict regulation sign with trained model for all available phosphosites
*/
process predict_regulation_sign {

    publishDir "${out_dir}", pattern: "phosphosite_regulation_sign/*", mode: 'copy'

    input:
        path 'input/phosphosites_features.tsv'

    output:
        path 'phosphosite_regulation_sign/all_phosphosites_regulation_sign.tsv', emit: tsv
        path 'phosphosite_regulation_sign/*pdf'

    script:
    """
    mkdir -p phosphosite_regulation_sign

    #predict_regulation_sign.py \
    #    input/phosphosites_features.tsv \
    #    phosphosite_regulation_sign \
    #    > phosphosite_regulation_sign/all_phosphosites_regulation_sign.tsv
    
    touch phosphosite_regulation_sign/all_phosphosites_regulation_sign.tsv
    touch phosphosite_regulation_sign/tmp.pdf
    
    """

}