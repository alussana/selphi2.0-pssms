#!/usr/bin/env python3

import sys
from functools import lru_cache
import warnings
import requests
import pandas as pd
    

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
    before_padding = '_' * before
    after_padding = '_' * after
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
    series['sequence'] = sequence
    return series
        

def main():
    
    # import kinase sequences
    
    kinases_tsv = sys.argv[1]
    kinases_df = pd.read_csv(kinases_tsv, sep='\t', index_col=0) 
    kinases_df.index.name = 'gene_name'
    kinases_df = kinases_df.rename(columns={'uniprot.1':'uniprot'})
    
    
    # get all kinase-phosphosite pairs
    
    k_p_tsv = sys.argv[2]
    k_p_df = pd.read_csv(k_p_tsv, sep='\t', index_col=0, header=None)
    k_p_df.columns = ['k_gene', 'p_gene', 'uniprot', 'position', 'residue']
    k_p_df.index.name = 'index'
    
    
    # fetch phosphosite sequences, filter out those with no sequence
    
    k_p_df = k_p_df.apply(make_k_p_entry, axis=1)
    k_p_df = k_p_df.dropna()
    
    
    # export dataset
    
    print(k_p_df.to_csv(sep='\t', header=True, index=True))    


if __name__ == '__main__':
    main()