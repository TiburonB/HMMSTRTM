# max_gamma.py
# TLB 1/19/21
# This program adds the gamma values across all entries in the db and s

import argparse
import os
from Atom import *
from os.path import exists
from os.path import isdir
def get_state_motif_9mer(dat, out_dir):
    if isdir(out_dir) == False:
        os.system('mkdir ' + out_dir)
        print('made directory ' + out_dir)
    for datum in dat:
        try:
            n = int(datum[0])
            codechain = datum[1]
            ind = int(datum[2])
        except:
            continue
        pdb_file = '../../Base/Warehouse/Protein/'+codechain[:4].upper()+'/PDB/'+codechain[:4].lower()+'.pdb'
        if exists(pdb_file)== False:
            continue
        AS = AtomSet.from_file(pdb_file)
        select_atoms = []
        for atom in AS.atoms:
            if atom.chain_id == codechain[4] and atom.res_num > ind -5 and atom.res_num < ind + 5:
                select_atoms.append(atom)
        new_AS = AtomSet(select_atoms)
        for i in range(0,len(new_AS.atoms)):
            if new_AS.atoms[i].res_num == ind:
                new_AS.atoms[i].temp_factor = 50.0
            else:
                new_AS.atoms[i].temp_factor = 20.0
        new_pdb = out_dir + 'HMMSTR_'+str(n)+'_'+codechain[:4].lower()+'_'+str(ind)+'.pdb'
        f = open(new_pdb, 'w')
        in_f = open(pdb_file, 'r')
        for l in in_f:
            if l.startswith('ATOM') == False and \
                    l.startswith('ANISOU') == False and \
                l.startswith('HETATM') == False:
                f.write(l)
        f.write(str(new_AS))
        print("wrote to " + str(new_pdb))
    return

def parse_args():
    parser = argparse.ArgumentParser(description = 'max_gamma.py')
    parser.add_argument('-mf', dest = 'model_file', default = "./HMMSTR/HMMSTR.hmm")
    parser.add_argument('-df', dest = 'gamma_drct', default = "./HMMSTR/HMMSTR_  1GAMMA.drct")
    parser.add_argument('-of', dest = 'out_file', default = "./max_gam.txt")
    parser.add_argument('-od', dest = 'out_dir', default = "./state_motifs/")
    args = parser.parse_args()
    return args

# 1/19 CLI
if __name__ == '__main__':
    args = parse_args()
    model_file = args.model_file
    gamma_drct = args.gamma_drct
    out_file = args.out_file
    out_dir = args.out_dir
    os.system("gfortran -c drct2profile.f95 gamma_io.f95 hmm_io.f95 max_gamma.f95")
    os.system("gfortran drct2profile.o gamma_io.o hmm_io.o max_gamma.o")
    s = "./a.out '" + model_file + "' '" + gamma_drct + "' '" + out_file + "'"
    print(s)
    os.system(s)
    data = [[l.strip() for l in k.split()] for k in open(out_file, 'r') ]
    get_state_motif_9mer(data, out_dir)
