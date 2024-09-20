#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import sklearn
import pickle as pkl
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV

def main():
    dataset_tsv = sys.argv[1]
    n_jobs = sys.argv[2]
    out_dir = sys.argv[3]
    
    X = pd.read_csv(dataset_tsv, sep='\t', index_col=0)
    y = X['status']
    X = X.drop('status', axis=1)
    X['is_Y'] = 0
    X['is_S'] = 0
    X['is_T'] = 0
    X.loc[X['residue']=='Y', 'is_Y'] = 1
    X.loc[X['residue']=='S', 'is_S'] = 1
    X.loc[X['residue']=='T', 'is_T'] = 1
    X = X.drop('residue', axis=1)
    
    param_grid = {
        'bootstrap': [True],
        'max_depth': [8, 16, 32],
        'min_samples_split': [4, 8, 16],
        'n_estimators': [128, 256, 512, 1024]
    }
    model = RandomForestClassifier()
    grid_search = GridSearchCV(estimator = model,
                               param_grid = param_grid,
                               cv = 3,
                               verbose = 2,
                               scoring = "roc_auc",
                               n_jobs=n_jobs)
    grid_search.fit(X, y)
    model = grid_search.best_estimator_
    
    skf = StratifiedKFold(n_splits=10)
    scores = cross_val_score(model, np.array(X), np.array(y), scoring = "roc_auc", cv = skf)
    scores = cross_val_score(model, X, y, scoring = "roc_auc", cv = skf)
    cross_probs = cross_val_predict(model, np.array(X), np.array(y), cv = skf, method="predict_proba")
    
   
if __name__ == '__main__':
    main()