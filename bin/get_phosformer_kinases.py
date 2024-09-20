#!/usr/bin/env python

import pandas as pd
import sys
    

def get_kinase_sequences(phosformer_dir: str, kinase_type: str):
    
    # Load the included csv file containing kinase domain sequences
    
    kinase_csv = pd.read_csv(f'{phosformer_dir}/data/reference_human_kinases.csv')
    
    # Select only Ser/Thr or Tyr kinases
    
    if kinase_type == 'S_T':
        kinase_df = kinase_csv.loc[kinase_csv['group'] != 'TK', ]
    elif kinase_type == 'Y':
        kinase_df = kinase_csv.loc[kinase_csv['group'] == 'TK', ]
        
    # Duplicate UniProt AC column to be the index
    
    kinase_df.index = kinase_df['uniprot']
        
    return kinase_df
    

def main():
    phosformer_dir = sys.argv[1]
    kinase_type = sys.argv[2]
    
    sys.path.insert(0, phosformer_dir)
    import Phosformer
    
    kinase_df = get_kinase_sequences(phosformer_dir, kinase_type)
    
    print(kinase_df.to_csv(sep='\t', header=True, index=True))


if __name__ == '__main__':
    main()