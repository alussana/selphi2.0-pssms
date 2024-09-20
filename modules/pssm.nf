#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Download Supplementary Table 2 from

Johnson, J.L., Yaron, T.M., Huntsman, E.M. et al.
An atlas of substrate specificities for the human serine/threonine kinome.
Nature 613, 759–766 (2023).

<https://doi.org/10.1038/s41586-022-05575-3>
*/
process get_ser_thr_kinome_2023_suppl_table_2 {

    publishDir "${out_dir}", pattern: "datasets/41586_2022_5575_MOESM4_ESM.xlsx", mode: 'copy'

    output:
        path 'datasets/41586_2022_5575_MOESM4_ESM.xlsx'

    script:
    """
    mkdir -p datasets
    
    wget -P datasets/ "${params.url_ser_thr_kinome_2023_suppl_table_2}"
    """

}

/*
Parse output from get_ser_thr_kinome_2023_suppl_table_2()

Save HDF5 datasets with Ser/Thr kinases PSSMs based on the 
"ser_thr_all_norm_scaled_matrice" sheet of 
"input/kinome_2023_suppl_table_2.xlsx" spreadsheet from [1]

Also save each PSSM in text format (pssm/*.tsv) and plot the
sequence logos (pssm/logos/*.svg)

[1]
Johnson, J.L., Yaron, T.M., Huntsman, E.M. et al.
An atlas of substrate specificities for the human serine/threonine kinome.
Nature 613, 759–766 (2023).
<https://doi.org/10.1038/s41586-022-05575-3>
*/
process get_ser_thr_kinase_pssm {

    publishDir "${out_dir}", pattern: "pssm/*tsv", mode: 'copy'
    publishDir "${out_dir}", pattern: "pssm/pssm_dict.h5", mode: 'copy'

    input:
        path 'input/kinome_2023_suppl_table_2.xlsx'

    output:
        path 'pssm/pssm_dict.h5', emit: pssm_dict_h5
        path 'pssm/*.tsv'

    script:
    """
    mkdir -p pssm/logos

    get_pssm_ser_thr_kinome_2023_suppl_table_2.py
    """

}

/*
Download Supplementary Table 2 from

Johnson, J.L., Yaron, T.M., Huntsman, E.M. et al.
An atlas of substrate specificities for the human serine/threonine kinome.
Nature 613, 759–766 (2023).

<https://doi.org/10.1038/s41586-022-05575-3>
*/
process get_tyr_kinome_2024_suppl_table_2 {

    publishDir "${out_dir}", pattern: "datasets/suppl_table_2.xlsx", mode: 'copy'

    output:
        path 'datasets/suppl_table_2.xlsx'

    script:
    """
    mkdir -p datasets
    
    wget -O datasets/suppl_table_2.xlsx "${params.url_tyr_kinome_2024_suppl_table_2}"
    """

}

/*
Parse output from get_tyr_kinome_2024_suppl_table_2()

Translate kinase name with gene synonym --> gene name id mapping

NOTE: some kinases, that are commonly described as Ser/Thr kinases but 
      show activity on Tyr residues too, are here reported as
      <KINASE_NAME>_TYR. These names are kept as they are

NOTE: PAK1 is only a synonym, not a gene name, according to the dictionary
      PKN1 is a gene name according to the dictionary
      But in the PSSM source files, both PAK1 and PKN1 are found and they
      correspond to two distinct PSSMs/kinases. Therefore the translation
      of PAK1 into PKN1 is not performed, otherwise two different PSSM would
      result having the same name, i.e. PKN1
      Same for:
      PDHK1 --> PDK1 translation
      ETK --> EPHA3 translation

Save HDF5 datasets with Tyr kinases PSSMs based on the 
"tyr_all_norm_scaled_matrice" sheet of 
"input/kinome_2024_suppl_table_2.xlsx" spreadsheet from [1]

Also save each PSSM in text format (pssm/*.tsv) and plot the
sequence logos (pssm/logos/*.svg)

[1]
Yaron-Barir, T.M., Joughin, B.A., Huntsman, E.M. et al. 
The intrinsic substrate specificity of the human tyrosine kinome. 
Nature (2024). 
<https://doi.org/10.1038/s41586-024-07407-y>
*/
process get_tyr_kinase_pssm {

    publishDir "${out_dir}", pattern: "pssm/*tsv", mode: 'copy'
    publishDir "${out_dir}", pattern: "pssm/tyr_pssm_dict.h5", mode: 'copy'

    input:
        path 'input/kinome_2024_suppl_table_2.xlsx'

    output:
        path 'pssm/tyr_pssm_dict.h5', emit: pssm_dict_h5
        path 'pssm/*.tsv'

    script:
    """
    mkdir -p pssm/logos

    cp \$(readlink input/kinome_2024_suppl_table_2.xlsx) kinome_2024_suppl_table_2.xlsx

    get_pssm_tyr_kinome_2024_suppl_table_2.py
    """

}