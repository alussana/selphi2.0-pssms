#!/usr/bin/env python3

import sys
import pandas as pd
import scipy.cluster.hierarchy as sch
import numpy as np
import radialtree as rt
import matplotlib.pyplot as plt


kinase_family_colors = {
    "Other": "#F8766D",
    "TK": "#D39200",
    "TKL": "#93AA00",
    "AGC": "#00BA38",
    "CMGC": "#00C19F",
    "Atypical": "#00B9E3",
    "CAMK": "#619CFF",
    "CK1": "#DB72FB",
    "STE": "#FF61C3",
    "NA": "#BEBEBE",
}

kinase_specificity_colors = {
    "SerThr": "#ADD8E6",
    "Tyr": "#E6DAA6",
    "Dual": "#000000",
}


def jaccard_similarity(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0

def jaccard_distance(set1, set2):
    return 1 - jaccard_similarity(set1, set2)


def main():
    selphi_k_p_tsv = sys.argv[1]
    out_pdf = sys.argv[2]
    out_tsv = sys.argv[3]
    kinase_family = sys.argv[4]
    """
    selphi_k_p_tsv = "k_p_associations.tsv"
    out_pdf = "circular_all.pdf"
    kinase_family = "all"
    """
   
    selphi_k_p = pd.read_csv(selphi_k_p_tsv, sep="\t", index_col=None)
    selphi_k_p.columns = ["Kinase","Specificity","Family","Phosphosite"]


    if kinase_family == "tyr":
        selphi_k_p = selphi_k_p.loc[selphi_k_p["Family"]=="TK"]
    elif kinase_family == "ser/thr":
        selphi_k_p = selphi_k_p.loc[selphi_k_p["Family"]!="TK"]


    # lookup table of kinase family memberships
    kinase_family_memberships = selphi_k_p[["Kinase", "Family"]].drop_duplicates().reset_index().drop("index", axis=1)
    kinase_specificity_class = selphi_k_p[["Kinase", "Specificity"]].drop_duplicates().reset_index().drop("index", axis=1)
    selphi_k_p.drop(["Specificity","Family"], axis=1, inplace=True)


    # make dictionary of phosphosite sets for each kinase
    kinase_phospho_sets = selphi_k_p.groupby('Kinase')['Phosphosite'].apply(set).to_dict()


    # compute distance matrix
    kinases = list(kinase_phospho_sets.keys())
    n = len(kinases)
    distance_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):  # only upper triangle needed since it's symmetric
            dist = jaccard_distance(
                kinase_phospho_sets[kinases[i]],
                kinase_phospho_sets[kinases[j]]
            )
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist  # mirror to lower triangle
    distance_df = pd.DataFrame(distance_matrix, index=kinases, columns=kinases)


    # create color list for leaves in the order of distance_df index
    kinase_family_memberships['Color code'] = kinase_family_memberships['Family'].map(kinase_family_colors)
    kinase_family_memberships["Family"].replace(np.nan, "NA", inplace=True)
    kinase_family_memberships["Color code"].replace(np.nan, "#BEBEBE", inplace=True)
    kinase_specificity_class['Color code'] = kinase_specificity_class['Specificity'].map(kinase_specificity_colors)
    kinase_specificity_class["Specificity"].replace(np.nan, "NA", inplace=True)
    kinase_specificity_class["Color code"].replace(np.nan, "#BEBEBE", inplace=True)
    kinases = distance_df.index.tolist()
    kinase_colors_family = kinase_family_memberships[["Kinase", "Color code"]]
    kinase_colors_specificity = kinase_specificity_class[["Kinase", "Color code"]]
    kinase_colors_family = dict(zip(
        kinase_family_memberships['Kinase'],
        kinase_family_memberships['Color code']
    ))
    kinase_colors_specificity = dict(zip(
        kinase_specificity_class['Kinase'],
        kinase_specificity_class['Color code']
    ))
    colors_dict = {
        "Kinase Family":[kinase_colors_family.get(kinase, '#BEBEBE') for kinase in kinases],
        "Kinase Specificity":[kinase_colors_specificity.get(kinase, '#BEBEBE') for kinase in kinases]
    }
    # '#BEBEBE' (gray) is a default color for any kinases not found in kinase_family_memberships
    colors_legends={
        "Kinase Family": {
            "colors": [v for k,v in kinase_family_colors.items()],
            "labels": [k for k,v in kinase_family_colors.items()]
        },
        "Kinase Specificity": {
            "colors": [v for k,v in kinase_specificity_colors.items()],
            "labels": [k for k,v in kinase_specificity_colors.items()]
        }
    }


    # convert the distance matrix to a condensed form
    condensed_dist = sch.distance.squareform(distance_df.values)


    # perform hierarchical clustering
    linkage_matrix = sch.linkage(
        condensed_dist,
        method='average',  # Options: 'single', 'complete', 'average', 'ward', etc.
        optimal_ordering=False  # Attempts to optimize the leaf ordering
    )    


    # compute the dendrogram
    dendrogram = sch.dendrogram(
        linkage_matrix,
        labels=distance_df.index,
        no_plot=True,
        #color_threshold=0
    )

    # plot circular dendrogram
    if kinase_family == "all":
        figsize = (14,6)
        fontsize = 3
    elif kinase_family == "tyr":
        figsize = (12,6)
        fontsize = 10
    elif kinase_family == "ser/thr":
        figsize = (14,6)
        fontsize = 4
    rt.plot(
        dendrogram,
        colorlabels=colors_dict,
        colorlabels_legend=colors_legends,
        figsize=figsize,
        fontsize=fontsize,
    )
    plt.tight_layout()
    plt.savefig(out_pdf)


    # get cluster assignments at a specific threshold
    threshold = 0.7
    clusters = sch.fcluster(linkage_matrix, threshold, criterion='distance')
    cluster_df = pd.DataFrame({
        'Kinase': distance_df.index,
        'Cluster': clusters
    })
    cluster_df.to_csv(out_tsv, sep="\t", index=False)
    

if __name__ == "__main__":
    main()
