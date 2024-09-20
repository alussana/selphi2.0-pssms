#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve 
from sklearn.metrics import auc


def compute_pr(labels, scores, n_points):
    precision, recall, thresholds = precision_recall_curve(labels, scores)
    
    precision = np.flip(precision)
    recall = np.flip(recall)
    
    interp_recall = np.linspace(0, 1, n_points)
    interp_precision = np.interp(x=interp_recall, xp=recall, fp=precision)
    
    pr_auc = auc(recall, precision)
    
    return(interp_recall, interp_precision, pr_auc)


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
    
    recall, precision, pr_auc = compute_pr(labels, scores, n_points)
    
    pr_points_df = pd.DataFrame({'Recall': recall, 'Precision': precision})
    
    pr_points_df.to_csv(f'{out_prefix}{id_int}_pr_points.tsv', sep='\t', index=False)
        
    with open(f'{out_prefix}{id_int}_pr_auc.txt', 'w') as pr_auc_fh:
        print(pr_auc, file=pr_auc_fh)

   
if __name__ == '__main__':
    main()