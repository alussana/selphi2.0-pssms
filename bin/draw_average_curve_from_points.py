#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt
import seaborn as sns

"""
x_str = 'FPR'
y_str = 'TPR'
file_list_txt = 'curves.txt'
out_pdf = 'pssm_model_100_rand_neg_sets/roc_curves.pdf'
"""

def main():
    x_str = sys.argv[1]
    y_str = sys.argv[2]
    file_list_txt = sys.argv[3]
    out_pdf = sys.argv[4]
    
    x_df = pd.DataFrame()
    y_df = pd.DataFrame()
    
    with open(file_list_txt) as file_list_fh:
        i = 0
        for line in file_list_fh:
            file_str = line.strip()
            df = pd.read_csv(file_str, sep='\t')
            x_df[i] = df[x_str]
            y_df[i] = df[y_str]
            i = i + 1
   
    # mean
    mean_x_series = x_df.apply(
        lambda x: x.mean(),
        axis=1
    )
    
    mean_y_series = y_df.apply(
        lambda x: x.mean(),
        axis=1
    )
    
    # confidence intervals
    #lower_ci_x_series = x_df.apply(
    #    lambda x: st.t.interval(confidence=0.95, df=len(x)-1, loc=np.mean(x), scale=st.sem(x))[0],
    #    axis=1
    #)
    #upper_ci_x_series = x_df.apply(
    #    lambda x: st.t.interval(confidence=0.95, df=len(x)-1, loc=np.mean(x), scale=st.sem(x))[1],
    #    axis=1
    #)
    #lower_ci_y_series = y_df.apply(
    #    lambda x: st.t.interval(confidence=0.95, df=len(x)-1, loc=np.mean(x), scale=st.sem(x))[0],
    #    axis=1
    #)
    #upper_ci_y_series = y_df.apply(
    #    lambda x: st.t.interval(confidence=0.95, df=len(x)-1, loc=np.mean(x), scale=st.sem(x))[1],
    #    axis=1
    #)
    
    # min and max
    max_y_series = y_df.apply(
        lambda x: x.max(),
        axis=1
    )
    
    min_y_series = y_df.apply(
        lambda x: x.min(),
        axis=1
    )
    
    
    fig, ax = plt.subplots(figsize=(4, 4))
    plt.plot(mean_x_series, mean_y_series, 'k-')
    plt.fill_between(x=mean_x_series, y1=min_y_series, y2=max_y_series, alpha=0.5, color='none', facecolor='black')
    ax.set(
        xlabel=f'{x_str}',
        ylabel=f'{y_str}'
    )
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    sns.despine()
    plt.tight_layout()
    plt.savefig(out_pdf)
    
   
if __name__ == '__main__':
    main()