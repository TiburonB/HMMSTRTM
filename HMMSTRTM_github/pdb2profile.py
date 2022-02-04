# pdb2profile.py
# given a pdb file, create a msa profile that can be read by HMMSTR.

import argparse
import os
from os import listdir
from os import isdir, isfile

def parse_args():
    parser = argparse.ArgumentParser(description = 'pdb2profile.py')
    parser.add_argument('-pf', dest = 'pdb_file', default = "./HMMS/HMMSTR.hmm")
    parser.add_argument('-nr', dest = 'nr_dataset_path', default = "../../Warehouse/Base/Externals/NR/")
    parser.add_argument('-rsync', dest = 'rsync_path', default = "./rsyncpdb.bin")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    pdb_file = args.pdb_file 
    nr_dataset_path = args.nr_dataset_path
    rsync_path = args.rsync_path

    if i

        

