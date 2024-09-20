#!/usr/bin/env python3

import sys
import h5py

def main():
    pssm_dict_h5 = sys.argv[1]
    
    # get kinases for which a PSSM exists
    pssms_h5 = h5py.File(pssm_dict_h5, 'r')
    for kinase in pssms_h5.keys():
        print(kinase)
    
   
if __name__ == '__main__':
    main()