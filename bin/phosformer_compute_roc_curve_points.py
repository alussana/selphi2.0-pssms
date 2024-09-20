#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
from sklearn.metrics import auc


def compute_roc(labels, scores, n_points):
    fpr, tpr, thresholds = roc_curve(labels, scores)
    
    interp_fpr = np.linspace(0, 1, n_points)
    interp_tpr = np.interp(x=interp_fpr, xp=fpr, fp=tpr)
    
    roc_auc = auc(fpr, tpr)
    
    return(interp_fpr, interp_tpr, roc_auc)
    

def main():
    id_int = int(sys.argv[1])
    scores_tsv = sys.argv[2]
    out_prefix = sys.argv[3]
    
    scores_df = pd.read_csv(scores_tsv, sep='\t', index_col=0)
    scores_df.index.name = 'example'
    
    scores_df['score'].loc[scores_df['score'].isna()] = 0
    
    labels = scores_df['true_interaction'] 
    scores = scores_df['score']

    n_points = 200
    
    fpr, tpr, roc_auc = compute_roc(labels, scores, n_points)
    
    roc_points_df = pd.DataFrame({'FPR': fpr, 'TPR': tpr})
    
    roc_points_df.to_csv(f'{out_prefix}{id_int}_roc_points.tsv', sep='\t', index=False)
    
    with open(f'{out_prefix}{id_int}_roc_auc.txt', 'w') as roc_auc_fh:
        print(roc_auc, file=roc_auc_fh)
    
   
if __name__ == '__main__':
    main()