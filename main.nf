#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { translate3col_keep_untranslated } from './modules/utils'
include { translate1col_keep_untranslated } from './modules/utils'
include { phony } from './modules/utils'
include { phony as phony2 } from './modules/utils'
include { split } from './modules/utils'
include { concatenate } from './modules/utils'
include { sort_uniq_w_header } from './modules/utils'
include { sort_uniq_w_header as sort_uniq_w_header2 } from './modules/utils'

include { get_ser_thr_kinome_2023_suppl_table_2 } from './modules/pssm'
include { get_ser_thr_kinase_pssm } from './modules/pssm'
include { get_tyr_kinome_2024_suppl_table_2 } from './modules/pssm'
include { get_tyr_kinase_pssm } from './modules/pssm'

include { dl_human_ref_proteome_ids } from './modules/id_mapping'
include { gene_name2uniprot_ac_dict } from './modules/id_mapping'
include { dl_human } from './modules/id_mapping'
include { IDa2uniprot2IDb } from './modules/id_mapping'
include { uniprot_ac2gene_name_dict } from './modules/id_mapping'

include { get_k_p } from './modules/pssm_score'
include { get_k_p_combinations } from './modules/pssm_score'
include { get_k_p_functional } from './modules/pssm_score'
include { get_k_p_regulation } from './modules/pssm_score'
include { get_k_p_feature_table } from './modules/pssm_score'
include { get_k_p_translated_feature_table } from './modules/pssm_score'
include { compute_ser_thr_pssm_score } from './modules/pssm_score'
include { compute_ser_thr_pssm_score_from_seq } from './modules/pssm_score'
include { compute_tyr_pssm_score } from './modules/pssm_score'
include { compute_tyr_pssm_score_from_seq } from './modules/pssm_score'
include { paste_features_tables } from './modules/pssm_score'
include { paste_complete_features_tables } from './modules/pssm_score'
include { paste_functional_features_tables } from './modules/pssm_score'
include { paste_regulation_features_tables } from './modules/pssm_score'

include { get_kinases_having_pssm } from './modules/pssm_model'
include { filter_set_for_kinases_having_pssm } from './modules/pssm_model'
include { filter_scores_for_kinases_having_pssm } from './modules/pssm_model'
include { eval_classifier_w_random_neg_set } from './modules/pssm_model'
include { draw_roc_curves } from './modules/pssm_model'
include { draw_pr_curves } from './modules/pssm_model'

include { get_p_sequences } from './modules/phosformer_model'
include { get_s_t_phosformer_kinases } from './modules/phosformer_model'
include { get_y_phosformer_kinases } from './modules/phosformer_model'
include { filter_set_for_phosformer_kinases } from './modules/phosformer_model'
include { filter_k_p_for_phosformer_kinases } from './modules/phosformer_model'
include { eval_phosformer_w_random_neg_set } from './modules/phosformer_model'
include { draw_roc_curves_phosformer } from './modules/phosformer_model'
include { draw_pr_curves_phosformer } from './modules/phosformer_model'
include { get_1col } from './modules/phosformer_model'
include { translate_k_p_pos_set } from './modules/phosformer_model'

include { make_all_phosphosites_features } from './modules/regulation_model'
include { make_regulation_train_dataset } from './modules/regulation_model'
include { train_regulation_model } from './modules/regulation_model'
include { predict_regulation_sign } from './modules/regulation_model'

// ===== //

workflow SER_THR_PSSMS {

    // NOT USED

    main:
        kinome_2023_suppl_table_2 = get_ser_thr_kinome_2023_suppl_table_2()
        pssm_out = get_ser_thr_kinase_pssm( kinome_2023_suppl_table_2 )
        pssm = pssm_out.pssm_dict_h5

    emit:
        pssm

}

workflow TYR_PSSMS {

    // NOT USED

    main:
        kinome_2024_suppl_table_2 = get_tyr_kinome_2024_suppl_table_2()
        pssm_out = get_tyr_kinase_pssm( kinome_2024_suppl_table_2 )
        pssm = pssm_out.pssm_dict_h5

    emit:
        pssm

}

workflow GENE_2_AC_ID_DICT {

    main:
        human_ref_ids = dl_human_ref_proteome_ids()

        gene_to_ac_dict = gene_name2uniprot_ac_dict( human_ref_ids )

        /*
        uniprot_id_dict = dl_human()

        synonym_to_name_dict = IDa2uniprot2IDb(
            uniprot_id_dict,
            'Gene_Synonym', 
            'Gene_Name' )
        */        

        dict = gene_to_ac_dict

    emit:
        dict
}

workflow AC_2_REF_GENENAME_DICT {

    main:
        id_mapping_table = dl_human()
        dict = uniprot_ac2gene_name_dict( id_mapping_table )

    emit:
        dict

}

workflow GENESYNONYM_2_GENENAME {

    main:
        uniprot_id_dict = dl_human()
        dict = IDa2uniprot2IDb( uniprot_id_dict,
                                'Gene_Synonym',
                                'Gene_Name' )

    emit:
        dict
}

workflow SER_THR_PSSM_SCORES {

    take:
        ser_thr_pssm_dict_h5
        gene_to_ac_dict

    main:
        //phosphosites = get_k_p()
        phosphosites = get_k_p_combinations()

        translated_phosphosites = translate3col_keep_untranslated(
             phosphosites,
             gene_to_ac_dict
        )

        phony(
            translated_phosphosites,
            'pssm_score/phosphosites_ac_pos_res.tsv'
        )

        translated_phosphosites_chunks = split(
            translated_phosphosites,
            50000
        ).flatten().map{ file -> file }

        phosphosites_pssm_scores_chunks = compute_ser_thr_pssm_score(
            translated_phosphosites_chunks,
            ser_thr_pssm_dict_h5
        ).collect()

        phosphosites_pssm_scores = concatenate(
            phosphosites_pssm_scores_chunks
        )

        phony2(
            phosphosites_pssm_scores,
            'pssm_score/k_p_ser_thr_pssm_scores.tsv'
        )

        features_table = paste_features_tables(
            phosphosites_pssm_scores
        )

    emit:
        translated_phosphosites
        phosphosites_pssm_scores
        features_table

}

workflow SER_THR_PSSM_SCORES_FROM_SEQ {

    take:
        ser_thr_pssm_dict_h5

    main:
        phosphosites = get_k_p_combinations()

        phosphosites_chunks = split(
            phosphosites,
            50000
        ).flatten().map{ file -> file }

        pssm_score_input_ch = phosphosites_chunks.combine(ser_thr_pssm_dict_h5)

        phosphosites_pssm_scores_chunks = compute_ser_thr_pssm_score_from_seq(
            pssm_score_input_ch
        ).collect()

        phosphosites_pssm_scores = concatenate(
            phosphosites_pssm_scores_chunks
        )

        phony2(
            phosphosites_pssm_scores,
            'pssm_score/k_p_ser_thr_pssm_scores.tsv'
        )

        features_table = paste_features_tables(
            phosphosites_pssm_scores
        )

    emit:
        phosphosites
        phosphosites_pssm_scores
        features_table

}

workflow TYR_PSSM_SCORES {

    take:
        ser_thr_features_table
        tyr_pssm_dict_h5
        gene_to_ac_dict

    main:
        phosphosites = get_k_p_translated_feature_table( ser_thr_features_table )

        translated_phosphosites = translate3col_keep_untranslated(
             phosphosites,
             gene_to_ac_dict
        )

        translated_phosphosites_chunks = split(
            translated_phosphosites,
            50000
        ).flatten().map{ file -> file }

        phosphosites_pssm_scores_chunks = compute_tyr_pssm_score(
            translated_phosphosites_chunks,
            tyr_pssm_dict_h5
        ).collect()

        phosphosites_pssm_scores = concatenate(
            phosphosites_pssm_scores_chunks
        )

        phony(
            phosphosites_pssm_scores,
            'pssm_score/k_p_tyr_pssm_scores.tsv'
        )

        features_table = paste_complete_features_tables(
            ser_thr_features_table,
            phosphosites_pssm_scores
        )

    emit:
        phosphosites_pssm_scores
        features_table

}

workflow TYR_PSSM_SCORES_FROM_SEQ {

    take:
        ser_thr_features_table
        phosphosites
        tyr_pssm_dict_h5

    main:
        phosphosites_chunks = split(
            phosphosites,
            50000
        ).flatten().map{ file -> file }

        pssm_score_input_ch = phosphosites_chunks.combine(tyr_pssm_dict_h5)

        phosphosites_pssm_scores_chunks = compute_tyr_pssm_score_from_seq(
            pssm_score_input_ch
        ).collect()

        phosphosites_pssm_scores = concatenate(
            phosphosites_pssm_scores_chunks
        )

        phony(
            phosphosites_pssm_scores,
            'pssm_score/k_p_tyr_pssm_scores.tsv'
        )

        features_table = paste_complete_features_tables(
            ser_thr_features_table,
            phosphosites_pssm_scores
        )

    emit:
        phosphosites_pssm_scores
        features_table

}

workflow FUNCTIONAL_PSSM_SCORES {

    take:
        ser_thr_pssm_dict_h5
        gene_to_ac_dict

    main:
        phosphosites = get_k_p_functional()

        translated_phosphosites = translate3col_keep_untranslated(
             phosphosites,
             gene_to_ac_dict
        )

        phony(
            translated_phosphosites,
            'pssm_score/phosphosites_functional_ac_pos_res.tsv'
        )

        translated_phosphosites_chunks = split(
            translated_phosphosites,
            50000
        ).flatten().map{ file -> file }

        phosphosites_pssm_scores_chunks = compute_pssm_score(
            translated_phosphosites_chunks,
            ser_thr_pssm_dict_h5
        ).collect()

        phosphosites_pssm_scores = concatenate(
            phosphosites_pssm_scores_chunks
        )

        phony2(
            phosphosites_pssm_scores,
            'pssm_score/k_p_functional_pssm_scores.tsv'
        )

        features_table = paste_functional_features_tables(
            phosphosites_pssm_scores
        )

    emit:
        phosphosites_pssm_scores
        features_table

}

workflow REGULATION_PSSM_SCORES {
    
    take:
        ser_thr_pssm_dict_h5
        gene_to_ac_dict

    main:
        phosphosites = get_k_p_regulation()

        translated_phosphosites = translate3col_keep_untranslated(
             phosphosites,
             gene_to_ac_dict
        )

        phony(
            translated_phosphosites,
            'pssm_score/phosphosites_regulation_ac_pos_res.tsv'
        )

        translated_phosphosites_chunks = split(
            translated_phosphosites,
            10000
        ).flatten().map{ file -> file }

        phosphosites_pssm_scores_chunks = compute_pssm_score(
            translated_phosphosites_chunks,
            ser_thr_pssm_dict_h5
        ).collect()

        phosphosites_pssm_scores = concatenate(
            phosphosites_pssm_scores_chunks
        )

        phony2(
            phosphosites_pssm_scores,
            'pssm_score/k_p_regulation_pssm_scores.tsv'
        )

        features_table = paste_regulation_features_tables(
            phosphosites_pssm_scores
        )

    emit:
        phosphosites_pssm_scores
        features_table

}

workflow SER_THR_PSSM_MODEL_ROC_PR_100_RAND_NEG_SETS {

    take:
        k_p_pssm_scores
        ser_thr_pssm_dict_h5

    main:
        kinases_having_pssm = get_kinases_having_pssm( ser_thr_pssm_dict_h5 )

        k_p_pos_set = Channel.fromPath( selphi_2_k_p_pos_set )

        k_p_pos_set_filtered = filter_set_for_kinases_having_pssm( kinases_having_pssm, k_p_pos_set )

        /*k_p_pssm_scores_chunks = split( k_p_pssm_scores, 500000 )
        k_p_pssm_scores_chunks_filtered = filter_scores_for_kinases_having_pssm( kinases_having_pssm, k_p_pssm_scores_chunks )
        k_p_pssm_scores_filtered = concatenate( k_p_pssm_scores_chunks_filtered.collect() )*/
        k_p_pssm_scores_filtered = filter_scores_for_kinases_having_pssm( kinases_having_pssm, k_p_pssm_scores )

        id = Channel.of( (1..100).toList() ).flatten()
        
        combined_ch = id.combine( k_p_pos_set_filtered ).combine( k_p_pssm_scores_filtered )

        eval_results = eval_classifier_w_random_neg_set( combined_ch )
        
        draw_roc_curves( eval_results.roc_points.collect(), 'S_T_PSSM' )

        draw_pr_curves( eval_results.pr_points.collect(), 'S_T_PSSM' )

}

workflow TYR_PSSM_MODEL_ROC_PR_100_RAND_NEG_SETS {

    take:
        k_p_pssm_scores
        tyr_pssm_dict_h5

    main:
        kinases_having_pssm = get_kinases_having_pssm( tyr_pssm_dict_h5 )

        k_p_pos_set = Channel.fromPath( selphi_2_k_p_pos_set )

        k_p_pos_set_filtered = filter_set_for_kinases_having_pssm( kinases_having_pssm, k_p_pos_set )

        /*k_p_pssm_scores_chunks = split( k_p_pssm_scores, 500000 )
        k_p_pssm_scores_chunks_filtered = filter_scores_for_kinases_having_pssm( kinases_having_pssm, k_p_pssm_scores_chunks )
        k_p_pssm_scores_filtered = concatenate( k_p_pssm_scores_chunks_filtered.collect() )*/
        k_p_pssm_scores_filtered = filter_scores_for_kinases_having_pssm( kinases_having_pssm, k_p_pssm_scores )

        id = Channel.of( (1..100).toList() ).flatten()
        
        combined_ch = id.combine( k_p_pos_set_filtered ).combine( k_p_pssm_scores_filtered )

        eval_results = eval_classifier_w_random_neg_set( combined_ch )
        
        draw_roc_curves( eval_results.roc_points.collect(), 'Y_PSSM' )

        draw_pr_curves( eval_results.pr_points.collect(), 'Y_PSSM' )

}

workflow SER_THR_PHOSFORMER_MODEL_ROC_PR_100_RAND_NEG_SETS {

    take:
        k_p
        ac_2_gene_dict
        genesynonym_2_genename_dict

    main:
        s_t_kinases_table = get_s_t_phosformer_kinases()

        s_t_kinases_table_translated = translate1col_keep_untranslated( s_t_kinases_table, ac_2_gene_dict )

        phony( s_t_kinases_table_translated,
               'phosformer_model_100_rand_neg_sets/s_t_phosformer_kinases_translated.tsv' )

        k_p_pos_set = Channel.fromPath( selphi_2_k_p_pos_set )

        k_p_pos_set_translated = translate_k_p_pos_set( k_p_pos_set, genesynonym_2_genename_dict)

        s_t_kinases = get_1col( s_t_kinases_table_translated )

        k_p_pos_set_filtered = filter_set_for_phosformer_kinases( s_t_kinases, k_p_pos_set_translated )

        k_p_filtered = filter_k_p_for_phosformer_kinases( s_t_kinases, k_p )

        id = Channel.of( (1..100).toList() ).flatten()

        combined_ch = id.combine( k_p_pos_set_filtered )
                        .combine( k_p_filtered )
                        .combine( s_t_kinases_table_translated )

        eval_results = eval_phosformer_w_random_neg_set( combined_ch )

        draw_roc_curves_phosformer( eval_results.roc_points.collect(), 'S_T_Phosformer' )

        draw_pr_curves_phosformer( eval_results.pr_points.collect(), 'S_T_Phosformer' )
}

workflow TYR_PHOSFORMER_MODEL_ROC_PR_100_RAND_NEG_SETS {

    take:
        k_p
        ac_2_gene_dict
        genesynonym_2_genename_dict

    main:
        y_kinases_table = get_y_phosformer_kinases()

        y_kinases_table_translated = translate1col_keep_untranslated( y_kinases_table, ac_2_gene_dict )

        phony( y_kinases_table_translated,
               'phosformer_model_100_rand_neg_sets/y_phosformer_kinases_translated.tsv' )

        k_p_pos_set = Channel.fromPath( selphi_2_k_p_pos_set )

        k_p_pos_set_translated = translate_k_p_pos_set( k_p_pos_set, genesynonym_2_genename_dict)

        y_kinases = get_1col( y_kinases_table_translated )

        k_p_pos_set_filtered = filter_set_for_phosformer_kinases( y_kinases, k_p_pos_set_translated )

        k_p_filtered = filter_k_p_for_phosformer_kinases( y_kinases, k_p )

        id = Channel.of( (1..100).toList() ).flatten()

        combined_ch = id.combine( k_p_pos_set_filtered )
                        .combine( k_p_filtered )
                        .combine( y_kinases_table_translated )

        eval_results = eval_phosformer_w_random_neg_set( combined_ch )

        draw_roc_curves_phosformer( eval_results.roc_points.collect(), 'Y_Phosformer' )

        draw_pr_curves_phosformer( eval_results.pr_points.collect(), 'Y_Phosformer' )

}

workflow REGULATION_MODEL {

    take:
        regulation_features_table
        k_p_features_table

    main:
        k_p_features_table_chunks = split( k_p_features_table, 131072 )
                                        .flatten().map{ file -> tuple( file.baseName, file ) }
        all_phosphosites_tsv_chunks = make_all_phosphosites_features( k_p_features_table_chunks )
        all_phosphosites_w_duplicates_tsv = concatenate( all_phosphosites_tsv_chunks.collect() )
        all_phosphosites_tsv = sort_uniq_w_header( all_phosphosites_w_duplicates_tsv )
        phony( all_phosphosites_tsv, "phosphosite_regulation_sign/phosphosites_features.tsv")

        train_dataset_w_duplicates_tsv = make_regulation_train_dataset( regulation_features_table )
        train_dataset_tsv = sort_uniq_w_header2( train_dataset_w_duplicates_tsv )
        phony2( train_dataset_tsv, 'phosphosite_regulation_sign/regulation_train_dataset.tsv' )
      
        model = train_regulation_model( train_dataset_tsv ).model
        
        //phosphosites_signed = predict_regulation_sign( model, all_phosphosites_tsv )

    emit:
        //model
        //phosphosites_signed
        all_phosphosites_tsv

}

workflow {

    // get pssms
    //ser_thr_pssm_dict_h5 = SER_THR_PSSMS().pssm
    ser_thr_pssm_dict_h5 = Channel.fromPath("${projectDir}/data/S_T_PSSMs.h5")
    //tyr_pssm_dict_h5 = TYR_PSSMS().pssm
    tyr_pssm_dict_h5 = Channel.fromPath("${projectDir}/data/Y_PSSMs.h5")

    // generate dictionary to map Gene Name to UniProt AC
    //id_dict = GENE_2_AC_ID_DICT()

    /*ac_2_gene_dict = AC_2_REF_GENENAME_DICT()

    genesynonym_2_genename_dict = GENESYNONYM_2_GENENAME()*/
    
    // compute pssm scores on ${selphi_2_features_table} phosphosites
    k_p_ser_thr_pssm_scores = SER_THR_PSSM_SCORES_FROM_SEQ( ser_thr_pssm_dict_h5 )
    k_p_pssm_scores = TYR_PSSM_SCORES_FROM_SEQ( k_p_ser_thr_pssm_scores.features_table,
                                                k_p_ser_thr_pssm_scores.phosphosites,
                                                tyr_pssm_dict_h5 )

    /*
    // compute pssm scores on ${selphi_2_regulation_features_table} phosphosites
    regulation_features_table = REGULATION_PSSM_SCORES( ser_thr_pssm_dict_h5,
                                                        id_dict ).features_table

    // compute pssm scores on ${selphi_2_functional_features_table} phosphosites
    functional_features_table = FUNCTIONAL_PSSM_SCORES( ser_thr_pssm_dict_h5,
                                                        id_dict ).features_table
    */

    // run 100x 10-fold cross validations of a linear classifier
    /*SER_THR_PSSM_MODEL_ROC_PR_100_RAND_NEG_SETS( k_p_ser_thr_pssm_scores.phosphosites_pssm_scores,
                                                 ser_thr_pssm_dict_h5 )

    TYR_PSSM_MODEL_ROC_PR_100_RAND_NEG_SETS( k_p_pssm_scores.phosphosites_pssm_scores,
                                             tyr_pssm_dict_h5 )

   
    // run 100x validation run with Phosformer using random negative sets
    SER_THR_PHOSFORMER_MODEL_ROC_PR_100_RAND_NEG_SETS( k_p_ser_thr_pssm_scores.phosphosites,
                                                       ac_2_gene_dict,
                                                       genesynonym_2_genename_dict )

    TYR_PHOSFORMER_MODEL_ROC_PR_100_RAND_NEG_SETS( k_p_ser_thr_pssm_scores.phosphosites,
                                                   ac_2_gene_dict,
                                                   genesynonym_2_genename_dict )


    // train classifier of protein activity sign regulated by a phosphosite and make phosphoproteome-wide predictions
    /*regulation_model = REGULATION_MODEL( regulation_features_table,
                                         k_p_ser_thr_pssm_scores.features_table )
    */

}