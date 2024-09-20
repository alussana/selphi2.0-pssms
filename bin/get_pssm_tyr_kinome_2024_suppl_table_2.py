#!/usr/bin/env python3

import numpy as np
import pandas as pd

def get_kinase_pssm(kinase_str:str, matrix_df:pd.DataFrame):
    aa_list = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T','s','t','y']
    pos_list = [-5, -4, -3, -2, -1, 1, 2, 3, 4]
    kinase_series = matrix_df.loc[kinase_str,]
    pssm_values_list = np.empty([len(pos_list), len(aa_list)], dtype = float)
    for i in range(len(pos_list)):
        pos = pos_list[i]
        for j in range(len(aa_list)):
            aa = aa_list[j]
            pssm_values_list[i][j] = kinase_series[f'{pos}{aa}']
    pssm = pd.DataFrame(pssm_values_list, index=pos_list, columns=aa_list)
    pos_0_df = pd.DataFrame(index=[0], columns=aa_list)
    pos_0_df['Y'] = 1
    pos_0_df = pos_0_df.fillna(0)
    pssm = pd.concat([pssm, pos_0_df])
    pssm = pssm.sort_index()    
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
        'input/kinome_2024_suppl_table_2.xlsx',
        sheet_name='tyrosine_all_norm_scaled_matric',
        index_col=0
    )
    kinase_list = list(norm_scaled_pssm_df.index)

    pssm_dict = {}
    for kinase_str in kinase_list:
        pssm_dict[kinase_str] = get_kinase_pssm(kinase_str, norm_scaled_pssm_df)
        pssm_dict[kinase_str].to_csv(
            path_or_buf=f'pssm/{kinase_str}.tsv',
            sep='\t',
            index=True,
            header=True
        )
    save_dict_to_hdf5(pssm_dict, 'pssm/tyr_pssm_dict.h5')

if __name__ == '__main__':
    main()