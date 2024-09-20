#!/usr/bin/env python3

import sys
import h5py
import pandas as pd
import requests
import warnings
from functools import lru_cache

AA_LIST = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T','s','t','y']
POSITIONS_LIST = list(range(-5,5))

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

def score_sequence(seq_str:str, pssm_df:pd.DataFrame):
    if len(seq_str) != len(pssm_df):
        raise Exception('Sequence length cannot be different from pssm length')
    """
    # transform pssm into probabilities, score sequence by multiplying values and scale by prob of random peptide
    pssm_df = pssm_df.apply(lambda x: x / x.sum(), axis=1)
    n_aa = len(AA_LIST)
    n_pos = len(seq_str)
    n_res = 0
    p = 1
    for i in range(n_pos):
        if seq_str[i] != '_':
            n_res = n_res + 1
            pos = list(pssm_df.index)[i]
            p = p * pssm_df.loc[pos, seq_str[i]] * 20
    """
    # score sequence by multiplying values and scale by prob of random peptide as described in Supplementary Note 2 of https://doi.org/10.1038/s41586-022-05575-3
    n_pos = len(seq_str)
    p = 1
    for i in range(n_pos):
        if seq_str[i] != '_':
            pos = list(pssm_df.index)[i]
            p = p * pssm_df.loc[pos, seq_str[i]]
    """
    # add scores instead of multiplying
    p = 0
    if pssm_df.loc[0, seq_str[5]] != 0:
        for i in range(len(seq_str)):
            if seq_str[i] != '_':
                pos = list(pssm_df.index)[i]
                p = p + pssm_df.loc[pos, seq_str[i]]
    """
    return p

def pssm_scoring(record: pd.Series, pssm_df_dict: dict, pssm_bg_scores: pd.DataFrame):
    try:
        kinase = record['Kinase']
        seq = record['Sequence'] 
        p = score_sequence(seq, pssm_df_dict[kinase])
        bg = pssm_bg_scores[kinase]
        den = len(bg)
        num = (bg < p).sum()
        record['PSSM Score'] = num / den
    except:
        record['PSSM Score'] = 0
    return record

def make_seqrnk_entry(series:pd.Series) -> pd.Series:
    uniprot_ac = series['UniProt AC']
    pos = series['Position']
    residue = series['Residue']
    sequence = get_phosphosite(
        uniprot_ac=uniprot_ac,
        residue=residue,
        pos=pos,
        before=5,
        after=4
    )
    score = series[-1]
    seqrnk_entry_series = pd.Series(
        {'Sequence': sequence, 'Score':score},
        name=series.name
    )
    return seqrnk_entry_series

def main():
    k_p_tsv = sys.argv[1]
    pssms_h5_file = sys.argv[2]
    pssm_bg_scores_tsv = sys.argv[3]
    """
    k_p_tsv = 'input/k_p.tsv'
    pssms_h5_file = 'input/pssm_dict.h5'
    pssm_bg_scores_tsv = 'input/pssm_background_scores.tsv.gz'
    """
    
    # read input data
    pssms_h5 = h5py.File(pssms_h5_file, 'r')
    pssm_df_dict = {}
    for kinase in pssms_h5.keys():
        pssm_df_dict[kinase] = pd.DataFrame(pssms_h5[kinase])
        pssm_df_dict[kinase].columns = AA_LIST
        pssm_df_dict[kinase].index = POSITIONS_LIST
    
    k_p = pd.read_csv(k_p_tsv, sep='\t', header=None, index_col=None)
    k_p.columns = ['Kinase', 'Target', 'UniProt AC', 'Position', 'Residue']
    #k_p = k_p.dropna()
    k_p['Position'] = k_p['Position'].astype(int)
    
    pssm_bg_scores = pd.read_csv(pssm_bg_scores_tsv, sep='\t')
    
    # get phosphosite sequences
    targets = k_p[['UniProt AC','Position','Residue']].drop_duplicates()
    seqrnk = targets.apply(make_seqrnk_entry, axis=1)
    targets['Sequence'] = seqrnk['Sequence']
    
    # add sequences to original dataframe
    k_p = pd.merge(k_p, targets, on=['UniProt AC','Position','Residue'], how='left')
    
    # compute PSSM scores
    k_p_pssm_scores = k_p.apply(pssm_scoring, args=[pssm_df_dict, pssm_bg_scores], axis=1)

    print(k_p_pssm_scores.to_csv(sep='\t', index=False, header=False))
    
if __name__ == '__main__':
    main()