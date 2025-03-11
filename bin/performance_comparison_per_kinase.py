#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
from itertools import product


def quantile_normalize(df: pd.DataFrame):
    """
    Perform quantile normalization on a pandas DataFrame.

    Args:
    df (pandas.DataFrame): DataFrame to be quantile normalized.

    Returns:
    pandas.DataFrame: Quantile normalized DataFrame.
    """
    # Rank the values in each column
    ranked = df.stack().groupby(df.rank(method="first").stack().astype(int)).mean()
    # Sort the ranks and use them as indices to get the quantiles
    quantiles = df.rank(method="min").stack().astype(int).map(ranked).unstack()
    # Sort the indices of the original DataFrame
    quantiles = quantiles.sort_index()
    return quantiles


def scale_01(df):
    # Subtract the minimum in each column
    scaled_df = df.apply(lambda x: x - x.min())
    # Divide by the maximum in each column
    scaled_df = scaled_df.apply(lambda x: x / x.max())
    return scaled_df


def top5percentKinases(
    kinase_activity_df, regulation_df, metadata_experiments_set
):
    df = pd.DataFrame()
    for experiment in metadata_experiments_set:
        upreg_kinases = list(
            regulation_df["Kinase"].loc[regulation_df["Experiment"] == experiment,]
        )
        if len(upreg_kinases) == 0:
            next
        else:
            activities = kinase_activity_df.loc[
                kinase_activity_df["Experiment"] == experiment
            ]
            # uniq_vals = activities['Kinase activity change'].sort_values(ascending=False).unique()
            uniq_vals = (
                activities["Kinase activity change"]
                .sort_values(ascending=True)
                .unique()
            )
            ranks = np.array([i for i in range(len(uniq_vals))])
            vals2ranks = dict(zip(uniq_vals, ranks))
            activities["Rank"] = [
                vals2ranks[i] for i in activities["Kinase activity change"]
            ]
            n_ranks = activities["Rank"].max()
            for kinase in upreg_kinases:
                rank_k = activities["Rank"].loc[activities["Kinase"] == kinase]
                if len(rank_k) == 0:
                    next
                else:
                    if rank_k.values[0] / n_ranks >= 0.95:
                        rank_k.iloc[0] = 1
                    else:
                        rank_k.iloc[0] = 0
                    df = pd.concat([df, rank_k])
    return df


def bottom5percentKinases(
    kinase_activity_df, regulation_df, metadata_experiments_set
):
    df = pd.DataFrame()
    for experiment in metadata_experiments_set:
        downreg_kinases = list(
            regulation_df["Kinase"].loc[regulation_df["Experiment"] == experiment,]
        )
        if len(downreg_kinases) == 0:
            next
        else:
            activities = kinase_activity_df.loc[
                kinase_activity_df["Experiment"] == experiment
            ]
            # uniq_vals = activities['Kinase activity change'].sort_values(ascending=False).unique()
            uniq_vals = (
                activities["Kinase activity change"]
                .sort_values(ascending=True)
                .unique()
            )
            ranks = np.array([i for i in range(len(uniq_vals))])
            vals2ranks = dict(zip(uniq_vals, ranks))
            activities["Rank"] = [
                vals2ranks[i] for i in activities["Kinase activity change"]
            ]
            n_ranks = activities["Rank"].max()
            for kinase in downreg_kinases:
                rank_k = activities["Rank"].loc[activities["Kinase"] == kinase]
                if len(rank_k) == 0:
                    next
                else:
                    if rank_k.values[0] / n_ranks <= 0.05:
                        rank_k.iloc[0] = 1
                    else:
                        rank_k.iloc[0] = 0
                    df = pd.concat([df, rank_k])
    return df


def make_kinase_activity_df(
    method_name:str,
    input_list_txt:str,
    kinase_activity_metric:str,
    out_prefix:str,
):
    # build dictionary (experiment id --> kinase activity table path)
    data_path_dict = {}
    with open(input_list_txt, "r") as input_list_fh:
        for line in input_list_fh:
            data_path_str = line.strip()
            data_id = os.path.basename(data_path_str)[:-4]
            data_path_dict[data_id] = data_path_str
    # build dataframe of kinase activities
    kinase_activity_df = pd.DataFrame()
    for data_id in data_path_dict.keys():
        df = pd.read_csv(data_path_dict[data_id], sep="\t", index_col=0)
        if len(df.columns) == 0:
            df = pd.read_csv(data_path_dict[data_id], sep=",", index_col=0)
        series = df[kinase_activity_metric]
        series.name = data_id
        kinase_activity_df = kinase_activity_df.join(series, how="outer")
    # normalise and scale
    kinase_activity_df = quantile_normalize(kinase_activity_df)
    kinase_activity_df = scale_01(kinase_activity_df)
    # melt
    kinase_activity_df = pd.melt(
        kinase_activity_df.assign(index=kinase_activity_df.index),
        id_vars=["index"],
    )
    kinase_activity_df.columns = [
        "Kinase",
        "Experiment",
        "Kinase activity change",
    ]
    # remove NAs due to missing phosphosite sequence
    kinase_activity_df = kinase_activity_df.dropna()
    # make index
    kinase_activity_df.index = kinase_activity_df.apply(
        lambda x: f'{x["Experiment"]}__{x["Kinase"]}', axis=1
    )
    kinase_activity_df.to_csv(
        f"{out_prefix}hernandez2017_kinase_activity_{method_name}.tsv",
        header=True,
        index=True,
        sep="\t",
    )
    return kinase_activity_df


def pairwise_comparison(
    method_1_name:str,
    kinase_activity_method_1_df:pd.DataFrame,
    method_2_name:str,
    kinase_activity_method_2_df:pd.DataFrame,
    metadata_df:pd.DataFrame,
    out_prefix:str,
):
    # Take only instances for which a kinase activity could be computed by both methods
    intersection_index = list(
        set(kinase_activity_method_1_df.index)
        .intersection(set(kinase_activity_method_2_df.index))
    )
    kinase_activity_method_1_df = kinase_activity_method_1_df.loc[intersection_index,]
    kinase_activity_method_2_df = kinase_activity_method_2_df.loc[intersection_index,]
    # make sets of indexes for positive and negative examples in upregulation and downregulation, separately
    upregulation_df = metadata_df.loc[metadata_df["Regulation"] == 1]
    downregulation_df = metadata_df.loc[metadata_df["Regulation"] == -1]
    upregulation_df.to_csv(
        f"{out_prefix}hernandez2017_{method_1_name}_{method_2_name}_upregulated.tsv",
        header=True,
        index=True,
        sep="\t"
    )
    downregulation_df.to_csv(
        f"{out_prefix}hernandez2017_{method_1_name}_{method_2_name}_downregulated.tsv",
        header=True,
        index=True,
        sep="\t",
    )
    metadata_kinases_set = set(metadata_df["Kinase"].unique())
    metadata_experiments_set = set(metadata_df["Experiment"].unique())
    possible_indexes_list = list(
        product(metadata_kinases_set, metadata_experiments_set)
    )
    possible_indexes_set = set([f"{i[1]}__{i[0]}" for i in possible_indexes_list])
    upregulation_negative_indexes_list = list(
        possible_indexes_set.difference(set(upregulation_df.index))
    )
    downregulation_negative_indexes_list = list(
        possible_indexes_set.difference(set(downregulation_df.index))
    )
    # top 5% kinases true positives performance per-kinase count
    method_1_top5percentKinases_df = (
        top5percentKinases(
            kinase_activity_method_1_df, upregulation_df, metadata_experiments_set
        )
    )
    method_2_top5percentKinases_df = (
        top5percentKinases(
            kinase_activity_method_2_df, upregulation_df, metadata_experiments_set
        )
    )
    method_1_bottom5percentKinases_df = (
        bottom5percentKinases(
            kinase_activity_method_1_df, downregulation_df, metadata_experiments_set
        )
    )
    method_2_bottom5percentKinases_df = (
        bottom5percentKinases(
            kinase_activity_method_2_df, downregulation_df, metadata_experiments_set
        )
    )
    method_1_joined5percentKinases_df = pd.concat(
        [method_1_top5percentKinases_df, method_1_bottom5percentKinases_df]
    )
    method_2_joined5percentKinases_df = pd.concat(
        [method_2_top5percentKinases_df, method_2_bottom5percentKinases_df]
    )
    method_1_joined5percentKinases_df["Kinase"] = method_1_joined5percentKinases_df.index
    method_1_joined5percentKinases_df["Kinase"] = method_1_joined5percentKinases_df.apply(
        lambda x: x["Kinase"].split("__")[1], axis=1
    )
    method_2_joined5percentKinases_df["Kinase"] = method_2_joined5percentKinases_df.index
    method_2_joined5percentKinases_df["Kinase"] = method_2_joined5percentKinases_df.apply(
        lambda x: x["Kinase"].split("__")[1], axis=1
    )
    method_1_hit_count_df = method_1_joined5percentKinases_df.groupby(by="Kinase").agg('sum')
    method_2_hit_count_df = method_2_joined5percentKinases_df.groupby(by="Kinase").agg('sum')
    method_1_hit_count_df.columns = [method_1_name]
    method_2_hit_count_df.columns = [method_2_name]
    hit_count_df = pd.concat([method_1_hit_count_df, method_2_hit_count_df], axis=1)
    hit_count_df = hit_count_df.replace(np.nan, 0)
    # radial plots
    plt.close()
    labels = list(hit_count_df.index)
    values_1 = list(hit_count_df[method_1_name])
    values_2 = list(hit_count_df[method_2_name])
    num_vars = len(values_1)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    values_1 += values_1[:1]
    values_2 += values_2[:1]
    angles += angles[:1]
    fig, ax = plt.subplots(figsize=(4, 4), subplot_kw={'projection': 'polar'})
    ax.plot(angles, values_1, marker='o', label=method_1_name)
    ax.plot(angles, values_2, marker='x', linestyle='--', label=method_2_name)
    ax.set_xticklabels([])
    for angle, label in zip(angles, labels):
        ax.text(angle, max(max(values_1), max(values_2)) + 2, label, ha='center', va='center')
    fig.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)
    ax.legend(loc='upper right', frameon=False, bbox_to_anchor=(1.3, 1.3))
    plt.savefig(f"{out_prefix}hernandez2017_{method_1_name}_{method_2_name}_topNpercentKinases_radial.pdf")


def main():
    input_list_phosx_txt = sys.argv[1]
    input_list_gsea_txt = sys.argv[2]
    input_list_kinex_txt = sys.argv[3]
    input_list_kstar_txt = sys.argv[4]
    input_list_ptmsea_txt = sys.argv[5]
    input_list_zscore_txt = sys.argv[6]
    metadata_tsv = sys.argv[7]
    kinase_activity_metric_str = sys.argv[8]
    out_prefix = sys.argv[9]
    """
    input_list_phosx_txt = 'input_files_phosx.txt'
    input_list_gsea_txt = 'input_files_gsea.txt'
    input_list_kinex_txt = 'input_files_kinex.txt'
    input_list_kstar_txt = 'input_files_kstar.txt'
    input_list_ptmsea_txt = 'input_files_ptmsea.txt'
    input_list_zscore_txt = 'input_files_zscore.txt'
    metadata_tsv = 'input/metadata.tsv'
    kinase_activity_metric_str = 'Activity Score'
    out_prefix = 'kinase_activity_benchmark/hernandez2017/pairwise/'
    """
    
    
    # metadata - ground truth kinase regulation
    metadata_df = pd.read_csv(metadata_tsv, sep="\t", index_col=None, header=None)
    metadata_df.columns = ["Experiment", "Kinase", "Regulation"]
    metadata_df.index = metadata_df.apply(
        lambda x: f'{x["Experiment"]}__{x["Kinase"]}', axis=1
    )
    

    gsea_kinase_activity_metric_str = "NES"


    # PhosX kinase activity
    kinase_activity_phosx_df = make_kinase_activity_df(
        "PhosX",
        input_list_phosx_txt,
        kinase_activity_metric_str,
        out_prefix,
    )


    # GSEApy kinase activity
    kinase_activity_gsea_df = make_kinase_activity_df(
        "GSEA",  
        input_list_gsea_txt,
        gsea_kinase_activity_metric_str,
        out_prefix,
    )
   

    # Kinex kinase activity
    kinase_activity_kinex_df = make_kinase_activity_df(
        "Kinex",
        input_list_kinex_txt,
        kinase_activity_metric_str,
        out_prefix,
    )

    # KSTAR kinase activity
    kinase_activity_kstar_df = make_kinase_activity_df(
        "KSTAR",
        input_list_kstar_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    

    # PTM-SEA kinase activity
    kinase_activity_ptmsea_df = make_kinase_activity_df(
        "PTM-SEA",
        input_list_ptmsea_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    

    # Z-score kinase activity
    kinase_activity_zscore_df = make_kinase_activity_df(
        "Z-score",
        input_list_zscore_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    
    
    # Pairwise comparisons: PhosX vs <method>
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "Kinex",
        kinase_activity_kinex_df,
        metadata_df,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "GSEApy",
        kinase_activity_gsea_df,
        metadata_df,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "KSTAR",
        kinase_activity_kstar_df,
        metadata_df,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "PTM-SEA",
        kinase_activity_ptmsea_df,
        metadata_df,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "Z-score",
        kinase_activity_zscore_df,
        metadata_df,
        out_prefix,
    )
    

if __name__ == "__main__":
    main()
