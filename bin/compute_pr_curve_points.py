#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve 
from sklearn.metrics import auc

"""
id_int = 0
pssm_scores_tsv = 'k_p_pssm_scores.tsv'
pos_examples_txt = 'pos_set.txt'
out_prefix = 'pssm_model_100_10fold_cv/pr_points/'
"""

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
    pssm_scores_tsv = sys.argv[2]
    pos_examples_txt = sys.argv[3]
    out_prefix = sys.argv[4]
    
    pssm_scores_df = pd.read_csv(pssm_scores_tsv, sep='\t', header=None, index_col=0)
    pssm_scores_df.columns = ['score']
    pssm_scores_df.index.name = 'example'
    
    pos_examples = []
    with open(pos_examples_txt) as pos_examples_fh:
        for line in pos_examples_fh:
            pos_examples.append(line.strip())
    
    pssm_scores_df['true_interaction'] = pssm_scores_df.apply(
        lambda x: x.name in pos_examples,
        axis=1
    )
    
    pssm_scores_df['score'].loc[pssm_scores_df['score'].isna()] = 0
    
    labels = pssm_scores_df['true_interaction'] 
    scores = pssm_scores_df['score']

    n_points = 200
    
    recall, precision, pr_auc = compute_pr(labels, scores, n_points)
    
    pr_points_df = pd.DataFrame({'Recall': recall, 'Precision': precision})
    
    pr_points_df.to_csv(f'{out_prefix}{id_int}_pr_points.tsv', sep='\t', index=False)
        
    with open(f'{out_prefix}{id_int}_pr_auc.txt', 'w') as pr_auc_fh:
        print(pr_auc, file=pr_auc_fh)

   
if __name__ == '__main__':
    main()