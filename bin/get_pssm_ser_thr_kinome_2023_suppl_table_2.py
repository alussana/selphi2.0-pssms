#!/usr/bin/env python3

import numpy as np
import pandas as pd

def get_kinase_pssm(kinase_str:str, matrix_df:pd.DataFrame, densitometry_df:pd.DataFrame):
    aa_list = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T','s','t','y']
    pos_list = [-5, -4, -3, -2, -1, 1, 2, 3, 4]
    kinase_series = matrix_df.loc[kinase_str,]
    densitometry_series = densitometry_df.loc[kinase_str,]
    #field_list = list(kinase_series.index)
    #pssm_values_list = np.array([[0 for aa in aa_list] for pos in pos_list], dtype=float)
    pssm_values_list = np.empty([len(pos_list), len(aa_list)], dtype = float)
    for i in range(len(pos_list)):
        pos = pos_list[i]
        for j in range(len(aa_list)):
            aa = aa_list[j]
            pssm_values_list[i][j] = kinase_series[f'{pos}{aa}']
    pssm = pd.DataFrame(pssm_values_list, index=pos_list, columns=aa_list)
    S_0 = np.empty(len(pos_list), dtype = float)
    T_0 = np.empty(len(pos_list), dtype = float)
    for i in range(len(pos_list)):
        pos = pos_list[i]
        S_0[i] = densitometry_series[f'{pos}S']
        T_0[i] = densitometry_series[f'{pos}T']
    S_S = S_0.sum()
    S_T = T_0.sum()
    S_ctrl = 0.75 * S_S - 0.25 * S_T
    T_ctrl = 0.75 * S_T - 0.25 * S_S
    S_0 = S_ctrl / max(S_ctrl, T_ctrl)
    T_0 = T_ctrl / max(S_ctrl, T_ctrl)
    pos_0_df = pd.DataFrame(index=[0], columns=aa_list)
    pos_0_df['S'] = S_0
    pos_0_df['T'] = T_0
    pos_0_df = pos_0_df.fillna(0)
    pssm = pd.concat([pssm, pos_0_df])
    pssm = pssm.sort_index()    
    #pssm = pssm.apply(lambda x: x / x.sum(), axis=1)
    return pssm

def save_dict_to_hdf5(dictionary:dict, output_file:str):
    """
    Save a Python dictionary to an HDF5 file.

    Parameters:
        pssm_dict (dict): The Python dictionary to be saved.
        output_file (str): The name of the HDF5 file to be created.

    Returns:
        None
    """
    import h5py
    try:
        with h5py.File(output_file, 'w') as h5file:
            # Traverse the dictionary and save each key-value pair as a dataset in the HDF5 file
            for key, value in dictionary.items():
                h5file.create_dataset(str(key), data=value)
        print("HDF5 file successfully created and data saved.")
    except Exception as e:
        print(f"An error occurred while saving the data: {e}")

def main():
    norm_scaled_pssm_df = pd.read_excel(
        'input/kinome_2023_suppl_table_2.xlsx',
        sheet_name='ser_thr_all_norm_scaled_matrice',
        index_col=0
    )
    kinase_list = list(norm_scaled_pssm_df.index)
    densitometry_df = pd.read_excel(
        'input/kinome_2023_suppl_table_2.xlsx',
        sheet_name='ser_thr_all_raw_matrices',
        index_col=0
    )
    """
    norm_pssm_df = pd.read_excel(
        '41586_2022_5575_MOESM4_ESM.xlsx',
        sheet_name='ser_thr_all_norm_matrices',
        index_col=0
    )
    """
    
    pssm_dict = {}
    for kinase_str in kinase_list:
        pssm_dict[kinase_str] = get_kinase_pssm(kinase_str, norm_scaled_pssm_df, densitometry_df)
        pssm_dict[kinase_str].to_csv(
            path_or_buf=f'pssm/{kinase_str}.tsv',
            sep='\t',
            index=True,
            header=True
        )
    save_dict_to_hdf5(pssm_dict, 'pssm/pssm_dict.h5')

if __name__ == '__main__':
    main()