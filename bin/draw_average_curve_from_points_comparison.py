#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


plt.rcParams['axes.titlesize'] = 8        
plt.rcParams['axes.labelsize'] = 6        
plt.rcParams['xtick.labelsize'] = 6       
plt.rcParams['ytick.labelsize'] = 6       
plt.rcParams['legend.fontsize'] = 6       
plt.rcParams['figure.titlesize'] = 8   
plt.rcParams['font.size'] = 6 


def main():
    x_str = sys.argv[1]
    y_str = sys.argv[2]
    pssm_file_list_txt = sys.argv[3]
    phosformer_file_list_txt = sys.argv[4]
    out_pdf = sys.argv[5]
    if len(sys.argv) > 6:
        title = sys.argv[6]

    
    pssm_x_df = pd.DataFrame()
    pssm_y_df = pd.DataFrame()
    phosformer_x_df = pd.DataFrame()
    phosformer_y_df = pd.DataFrame()
    

    # read curves
    with open(pssm_file_list_txt) as file_list_fh:
        i = 0
        for line in file_list_fh:
            file_str = line.strip()
            df = pd.read_csv(file_str, sep='\t')
            pssm_x_df[i] = df[x_str]
            pssm_y_df[i] = df[y_str]
            i = i + 1
    with open(phosformer_file_list_txt) as file_list_fh:
        i = 0
        for line in file_list_fh:
            file_str = line.strip()
            df = pd.read_csv(file_str, sep='\t')
            phosformer_x_df[i] = df[x_str]
            phosformer_y_df[i] = df[y_str]
            i = i + 1
   

    # mean
    pssm_mean_x_series = pssm_x_df.apply(
        lambda x: x.mean(),
        axis=1
    )
    pssm_mean_y_series = pssm_y_df.apply(
        lambda x: x.mean(),
        axis=1
    )
    phosformer_mean_x_series = phosformer_x_df.apply(
        lambda x: x.mean(),
        axis=1
    )
    phosformer_mean_y_series = phosformer_y_df.apply(
        lambda x: x.mean(),
        axis=1
    )


    # AUC
    pssm_auc = round(np.trapezoid(pssm_mean_y_series, pssm_mean_x_series), 3)
    phosformer_auc = round(np.trapezoid(phosformer_mean_y_series, phosformer_mean_x_series), 3)
    
    
    # min and max
    pssm_max_y_series = pssm_y_df.apply(
        lambda x: x.max(),
        axis=1
    )
    
    pssm_min_y_series = pssm_y_df.apply(
        lambda x: x.min(),
        axis=1
    )
    phosformer_max_y_series = phosformer_y_df.apply(
        lambda x: x.max(),
        axis=1
    )
    phosformer_min_y_series = phosformer_y_df.apply(
        lambda x: x.min(),
        axis=1
    )


    # plot
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    plt.plot(phosformer_mean_x_series, phosformer_mean_y_series, 'k-', label=f'Phosformer (AUC = {phosformer_auc})')
    plt.fill_between(x=phosformer_mean_x_series, y1=phosformer_min_y_series, y2=phosformer_max_y_series, alpha=0.5, color='none', facecolor='black')
    plt.plot(pssm_mean_x_series, pssm_mean_y_series, 'r-', label=f'PSSM (AUC = {pssm_auc})')
    plt.fill_between(x=pssm_mean_x_series, y1=pssm_min_y_series, y2=pssm_max_y_series, alpha=0.5, color='none', facecolor='red')
    ax.set(
        xlabel=f'{x_str}',
        ylabel=f'{y_str}'
    )
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_aspect('equal')
    plt.legend(loc='lower right', frameon=False)
    if 'title' in locals():
        plt.title(title)
    sns.despine()
    plt.tight_layout()
    plt.savefig(out_pdf)
    
   
if __name__ == '__main__':
    main()