#!/usr/bin/env python3

import sys
from functools import lru_cache
import warnings
import requests
import pandas as pd
from sklearn.metrics import roc_curve
from sklearn.metrics import auc

    
@lru_cache
def get_fasta(uniprot_ac:str):
    url_prefix = 'https://rest.uniprot.org/uniprotkb/'
    r = requests.get(f'{url_prefix}{uniprot_ac}.fasta')
    fasta = r.content.decode()
    lines = fasta.split('\n')
    lines.pop(0)
    sequence = ''.join(lines)
    return sequence
    
    
def get_phosphosite(uniprot_ac:str, residue:str, pos:int, before:int, after:int):
    if before < 0 or after < 0:
        raise Exception(
            """
            'before' and 'after' have to be positive integers.
            """
        )
    sequence = get_fasta(uniprot_ac)
    seq_length = len(sequence)
    before_padding = '-' * before
    after_padding = '-' * after
    sequence = before_padding + sequence + after_padding
    try:
        phosphosite = sequence[pos-1:pos+before+after]
        if phosphosite[before] != residue:
            warnings.warn(
                f"""
                Phosphosite residue does not match the UniProt sequence position.
                uniprot_ac: {uniprot_ac}
                pos: {pos}
                expected residue: {residue}
                matched residue: {phosphosite[before]}
                """
            )
        else:
            return phosphosite
    except:
        warnings.warn(
            f"""
            Phosphosite position exceedes sequence length.
            pos: {pos}
            sequence length: {seq_length}
            """
        )
        
        
def make_k_p_entry(series:pd.Series) -> pd.Series:
    uniprot_ac = series['uniprot']
    pos = series['position']
    residue = series['residue']
    sequence = get_phosphosite(
        uniprot_ac=uniprot_ac,
        residue=residue,
        pos=pos,
        before=5,
        after=5
    )
    return sequence
        

def main():
    
    # import kinase sequences
    
    kinases_tsv = sys.argv[1]
    kinases_df = pd.read_csv(kinases_tsv, sep='\t', index_col=0) 
    kinases_df.index.name = 'gene_name'
    kinases_df = kinases_df.rename(columns={'uniprot.1':'uniprot'})

    
    # get entries for positive k-p examples
    
    k_p_pos_set_tsv = sys.argv[2]
    k_p_pos_df = pd.read_csv(k_p_pos_set_tsv, sep='\t', index_col=0, header=None)
    k_p_pos_df.columns = ['k_gene', 'p_gene', 'uniprot', 'position', 'residue']
    k_p_pos_df.index.name = 'index'
    
            
    # get all kinase-phosphosite pairs
    
    k_p_tsv = sys.argv[3]
    k_p_df = pd.read_csv(k_p_tsv, sep='\t', index_col=0, header=None)
    k_p_df.columns = ['k_gene', 'p_gene', 'uniprot', 'position', 'residue']
    k_p_df.index.name = 'index'
    
    
    # get kinase sequences for positive examples
    
    k_p_pos_df['k_sequence'] = k_p_pos_df.apply(lambda x: kinases_df.loc[x['k_gene'], 'sequence'] if x['k_gene'] in list(kinases_df.index) else None, axis=1)

    
    # get phosphosite sequences for positive examples

    k_p_pos_df['p_sequence'] = k_p_pos_df.apply(make_k_p_entry, axis=1)
    
    
    # remove entries that do not have sequences
    
    k_p_pos_df = k_p_pos_df.dropna()
    

    # generate entries for negative examples
    
    neg_set_size_factor = float(sys.argv[4])
    n_neg = len(k_p_pos_df) * neg_set_size_factor
    
    k_p_df = k_p_df.drop(k_p_pos_df.index)
    
    i = 0
    k_p_neg_df = pd.DataFrame(columns=k_p_pos_df.columns)
    while i < n_neg:
        entry = k_p_df.sample(n=1)
        entry['k_sequence'] = entry.apply(lambda x: kinases_df.loc[x['k_gene'], 'sequence'] if x['k_gene'] in list(kinases_df.index) else None, axis=1)
        entry['p_sequence'] = entry.apply(make_k_p_entry, axis=1)
        entry = entry.dropna()
        if len(entry) == 0:
            next
        elif str(entry.index[0]) in list(k_p_neg_df.index):
            next
        else:
            k_p_neg_df = pd.concat([k_p_neg_df, entry])
            i = i + 1


    # concatenate positive and negative examples
    
    k_p_pos_df['true_interaction'] = True
    k_p_neg_df['true_interaction'] = False
    k_p_eval_df = pd.concat([k_p_pos_df, k_p_neg_df], axis=0)
    
    # delete big objects from memory
    
    del k_p_df
    
    
    # import Phosphormer library
    
    phosformer_dir = sys.argv[5]
    
    sys.path.insert(0, phosformer_dir)
    import Phosformer
    
    
    # load Phosformer model and tokenizer
    
    model = Phosformer.RobertaForSequenceClassification.from_pretrained('waylandy/phosformer')
    tokenizer = Phosformer.RobertaTokenizer.from_pretrained('waylandy/phosformer')
    
    
    # disables dropout for deterministic results
    
    model.eval()
    
    
    # make predictions
    
    def make_prediction(x):
        kinase_sequence = x['k_sequence']
        peptide_sequence = x['p_sequence']
        prediction = Phosformer.predict_one(kinase_sequence, peptide_sequence, model=model, tokenizer=tokenizer)
        return prediction    
    
    k_p_eval_df['score'] = k_p_eval_df.apply(make_prediction, axis=1)
    
    
    # export dataset with predictions
    
    print(k_p_eval_df.to_csv(sep='\t', header=True, index=True))
        
    

if __name__ == '__main__':
    main()