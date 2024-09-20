#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
.META: ${selphi_2_features_table}
     1  prot1
     2  psites
     3  residue
     4  score
     5  coreg293
     6  functional_score
     7  gtex
     8  rna.tiss
     9  rna.cell
    10  sel.1
    11  sel.2
    12  coreg.ntera
    13  coreg.hl
    14  coreg.mcf7
    15  is_disopred
    16  disopred_score
    17  exp3d_ala_ddG_effect
    18  exp3d_acid_ddG_effect
    19  log10_hotspot_pval_min
    20  isHotspot
    21  isInterface
    22  adj_ptms_w21
    23  netpho_max_all
    24  netpho_max_KIN
    25  paxdb_abundance_log10
    26  w0_mya
    27  w3_mya
    28  quant_top1
    29  quant_top5
    30  PWM_max_mss
    31  ACCpro
    32  SSpro
    33  SSpro8
    34  sift_min_score
    35  sift_mean_score
    36  sift_ala_score
    37  sift_acid_score
    38  isProteinDomain
    39  isProteinKinaseDomain
    40  isUniprotRegion
    41  isCytoplasmic
    42  isELMLinearMotif
    43  isEV_ala_prediction_epistatic5
    44  isKinaseCoreg

.META: pssm_score/k_p.tsv
     1  Kinase (gene name)
     2  Phosphorylated target (gene name)
     3  Phosphorylated target (gene name)
     4  Position (int)
     5  Residue (amino acid symbol)
*/
process get_k_p {

    publishDir "${out_dir}", pattern: "pssm_score/*tsv", mode: 'copy'

    output:
        path 'pssm_score/*tsv'

    script:
    """
    mkdir -p pssm_score

    ln -s ${selphi_2_features_table} selphi_2_features_table.tsv

    cat selphi_2_features_table.tsv \
        | sed '1d' | cut -f1-3 | sed 's/ /\\t/g' \
        > k_p.tsv

    cat k_p.tsv \
        | awk '{print \$1"\\t"\$2"\\t"\$2"\\t"\$3"\\t"\$4}' \
        > pssm_score/k_p.tsv
    """

}


/*
.META: pssm_score/k_p_combinations.csv
    1   Kinase (gene name)
    2   target (gene name)
    3   Residue (amino acid symbol)target (UniProt AC)
    4   Position (int)
    5   target (UniProt AC)
    6   Sequence
*/
process get_k_p_combinations {

    publishDir "${out_dir}", pattern: "pssm_score/k_p_combinations.csv", mode: 'copy'

    output:
        path 'pssm_score/k_p_combinations.csv'

    script:
    """
    mkdir -p pssm_score

    cat ${selphi_2_k_p_info} \
        | sed '1d' \
        | awk -F"," '{ \
            \$6 = toupper(\$6); \
            \$6 = substr(\$6, 3, length(\$6) - 5); \
            print \$1","\$2","\$3","\$4","\$5","\$6 \
        }' > pssm_score/k_p_combinations.csv
    """

}


/*
*/
process get_k_p_functional {

    publishDir "${out_dir}", pattern: "pssm_score/*tsv", mode: 'copy'

    output:
        path 'pssm_score/*tsv'

    script:
    """
    mkdir -p pssm_score

    ln -s ${selphi_2_functional_features_table} selphi_2_functional_features_table.tsv

    cat selphi_2_functional_features_table.tsv \
        | sed '1d' | cut -f1-3 | sed 's/ /\\t/g' \
        > k_p.tsv

    cat k_p.tsv \
        | awk '{print \$1"\\t"\$2"\\t"\$2"\\t"\$3"\\t"\$4}' \
        > pssm_score/k_p_functional.tsv
    """

}

/*
*/
process get_k_p_regulation {

    publishDir "${out_dir}", pattern: "pssm_score/*tsv", mode: 'copy'

    output:
        path 'pssm_score/*tsv'

    script:
    """
    mkdir -p pssm_score

    ln -s ${selphi_2_regulation_features_table} selphi_2_regulation_features_table.tsv

    cat selphi_2_regulation_features_table.tsv \
        | sed '1d' | cut -f1,2,4 | sed 's/ /\\t/g' \
        > k_p.tsv

    cat k_p.tsv \
        | awk '{print \$1"\\t"\$2"\\t"\$2"\\t"\$3"\\t"\$4}' \
        > pssm_score/k_p_regulation.tsv
    """

}

/*
.META: ${selphi_2_features_table}
     1  prot1
     2  psites
     3  residue
     4  score
     5  coreg293
     6  functional_score
     7  gtex
     8  rna.tiss
     9  rna.cell
    10  sel.1
    11  sel.2
    12  coreg.ntera
    13  coreg.hl
    14  coreg.mcf7
    15  is_disopred
    16  disopred_score
    17  exp3d_ala_ddG_effect
    18  exp3d_acid_ddG_effect
    19  log10_hotspot_pval_min
    20  isHotspot
    21  isInterface
    22  adj_ptms_w21
    23  netpho_max_all
    24  netpho_max_KIN
    25  paxdb_abundance_log10
    26  w0_mya
    27  w3_mya
    28  quant_top1
    29  quant_top5
    30  PWM_max_mss
    31  ACCpro
    32  SSpro
    33  SSpro8
    34  sift_min_score
    35  sift_mean_score
    36  sift_ala_score
    37  sift_acid_score
    38  isProteinDomain
    39  isProteinKinaseDomain
    40  isUniprotRegion
    41  isCytoplasmic
    42  isELMLinearMotif
    43  isEV_ala_prediction_epistatic5
    44  isKinaseCoreg
    45  SerThrPSSMscore

.META: pssm_score/k_p.tsv
     1  Kinase (gene name)
     2  Phosphorylated target (gene name)
     3  Phosphorylated target (gene name)
     4  Position (int)
     5  Residue (amino acid symbol)
*/
process get_k_p_feature_table {

    publishDir "${out_dir}", pattern: "pssm_score/*tsv", mode: 'copy'

    input:
        path 'input/feature_table.tsv'

    output:
        path 'pssm_score/*tsv'

    script:
    """
    mkdir -p pssm_score

    ln -s ${selphi_2_features_table} selphi_2_features_table.tsv

    cat selphi_2_features_table.tsv \
        | sed '1d' | cut -f1-3 | sed 's/ /\\t/g' \
        > k_p.tsv

    cat k_p.tsv \
        | awk '{print \$1"\\t"\$2"\\t"\$2"\\t"\$3"\\t"\$4}' \
        > pssm_score/k_p_features_table.tsv
    """

}

/*
*/
process get_k_p_translated_feature_table {

    publishDir "${out_dir}", pattern: "pssm_score/*tsv", mode: 'copy'

    input:
        path 'input/feature_table.tsv'

    output:
        path 'pssm_score/*tsv'

    script:
    """
    mkdir -p pssm_score

    cat input/feature_table.tsv \
        | sed '1d' | cut -f1-3 | sed 's/ /\\t/g' \
        > k_p.tsv

    cat k_p.tsv \
        | awk '{print \$1"\\t"\$2"\\t"\$2"\\t"\$3"\\t"\$4}' \
        > pssm_score/k_p_features_table.tsv
    """

}

/*
*/
process compute_ser_thr_pssm_score {

    memory '64G'

    input:
        path 'input/k_p.tsv'
        path 'input/pssm_dict.h5'

    output:
        path 'pssm_score/*tsv'

    script:
    """
    mkdir -p pssm_score

    ln -s ${ser_thr_pssm_bg_scores} input/pssm_background_scores.tsv.gz

    uuid=\$(basename \$(readlink -f input/k_p.tsv))
    #uuid=\$(cat /proc/sys/kernel/random/uuid)

    compute_pssm_scores.py \
        input/k_p.tsv \
        input/pssm_dict.h5 \
        input/pssm_background_scores.tsv.gz \
        > pssm_score/k_p_pssm_scores_\${uuid}.tsv
    """

}


/*
*/
process compute_ser_thr_pssm_score_from_seq {

    memory '64G'

    input:
        tuple file('input/k_p.tsv'),
              file('input/pssm_dict.h5')

    output:
        path 'pssm_score/*tsv'

    script:
    """
    mkdir -p pssm_score

    ln -s ${ser_thr_pssm_bg_scores} input/pssm_background_scores.tsv.gz

    uuid=\$(basename \$(readlink -f input/k_p.tsv))
    #uuid=\$(cat /proc/sys/kernel/random/uuid)

    compute_pssm_scores_from_seq.py \
        input/k_p.tsv \
        input/pssm_dict.h5 \
        input/pssm_background_scores.tsv.gz \
        > pssm_score/k_p_pssm_scores_\${uuid}.tsv
    """

}


/*
*/
process compute_tyr_pssm_score {

    memory '64G'

    input:
        path 'input/k_p.tsv'
        path 'input/pssm_dict.h5'

    output:
        path 'pssm_score/*tsv'

    script:
    """
    mkdir -p pssm_score

    ln -s ${tyr_pssm_bg_scores} input/pssm_background_scores.tsv.gz

    uuid=\$(basename \$(readlink -f input/k_p.tsv))
    #uuid=\$(cat /proc/sys/kernel/random/uuid)

    compute_pssm_scores.py \
        input/k_p.tsv \
        input/pssm_dict.h5 \
        input/pssm_background_scores.tsv.gz \
        > pssm_score/k_p_pssm_scores_\${uuid}.tsv
    """

}


/*
*/
process compute_tyr_pssm_score_from_seq {

    memory '64G'

    input:
        tuple file('input/k_p.tsv'),
              file('input/pssm_dict.h5')

    output:
        path 'pssm_score/*tsv'

    script:
    """
    mkdir -p pssm_score

    ln -s ${tyr_pssm_bg_scores} input/pssm_background_scores.tsv.gz

    uuid=\$(basename \$(readlink -f input/k_p.tsv))
    #uuid=\$(cat /proc/sys/kernel/random/uuid)

    compute_pssm_scores_from_seq.py \
        input/k_p.tsv \
        input/pssm_dict.h5 \
        input/pssm_background_scores.tsv.gz \
        > pssm_score/k_p_pssm_scores_\${uuid}.tsv
    """

}


/*
*/
process paste_features_tables {

    memory '32G'

    input:
        path 'input/k_p_pssm_scores.tsv'

    output:
        path 'features_table/k_p_features.tsv'

    script:
    """
    mkdir -p features_table

    old_header=\$(cat ${selphi_2_features_table} | sed -n '1p')

    header=\$(echo -e "\${old_header}\\tSerThrPSSMscore")

    cat \
        <(echo \${header}) \
        <(paste \
            <(cat ${selphi_2_features_table} \
                | sed '1d' \
            ) \
            <(cat input/k_p_pssm_scores.tsv \
                | cut -f7 \
            ) \
        ) \
    > features_table/k_p_features.tsv

    """

}

/*
*/
process paste_complete_features_tables {

    publishDir "${out_dir}", pattern: "features_table/*tsv", mode: 'copy'

    memory '32G'

    input:
        path 'input/features_table.tsv'
        path 'input/k_p_pssm_scores.tsv'

    output:
        path 'features_table/k_p_features.tsv'

    script:
    """
    mkdir -p features_table

    old_header=\$(cat input/features_table.tsv | sed -n '1p')

    header=\$(echo -e "\${old_header}\\tSerThrPSSMscore\\tTyrPSSMscore")

    cat \
        <(echo \${header}) \
        <(paste \
            <(cat input/features_table.tsv \
                | sed '1d' \
            ) \
            <(cat input/k_p_pssm_scores.tsv \
                | cut -f7 \
            ) \
        ) \
    > features_table/k_p_features.tsv

    """

}

/*
*/
process paste_functional_features_tables {

    publishDir "${out_dir}", pattern: "features_table/*tsv", mode: 'copy'

    memory '32G'

    input:
        path 'input/k_p_pssm_scores.tsv'

    output:
        path 'features_table/k_p_functional_features.tsv'

    script:
    """
    mkdir -p features_table

    old_header=\$(cat ${selphi_2_functional_features_table} | sed -n '1p')

    header=\$(echo -e "\${old_header}\\tSerThrPSSMscore")

    cat \
        <(echo \${header} | sed 's/ /\\t/g') \
        <(paste \
            <(cat ${selphi_2_functional_features_table} \
                | sed '1d' \
            ) \
            <(cat input/k_p_pssm_scores.tsv \
                | cut -f7 \
            ) \
        ) \
    > features_table/k_p_functional_features.tsv

    """

}

/*
*/
process paste_regulation_features_tables {

    publishDir "${out_dir}", pattern: "features_table/*tsv", mode: 'copy'

    memory '32G'

    input:
        path 'input/k_p_pssm_scores.tsv'

    output:
        path 'features_table/k_p_regulation_features.tsv'

    script:
    """
    mkdir -p features_table

    old_header=\$(cat ${selphi_2_regulation_features_table} | sed -n '1p')

    header=\$(echo -e "\${old_header}\\tSerThrPSSMscore")

    cat \
        <(echo \${header} | sed 's/ /\\t/g') \
        <(paste \
            <(cat ${selphi_2_regulation_features_table} \
                | sed '1d' \
            ) \
            <(cat input/k_p_pssm_scores.tsv \
                | cut -f7 \
            ) \
        ) \
    > features_table/k_p_regulation_features.tsv

    """

}