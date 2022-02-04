# gamma_9mer.py

in_f = "m_gam.txt"
dat = [[l.strip() for l in k.split()] for k in open(in_f,'r')]
from Atom import *
out_dir = './state_motifs/'
for datum in dat:
    n = int(datum[0])
    codechain = datum[1]
    ind = int(datum[2])

    pdb_file = '../../Base/Warehouse/Protein/'+codechain[:4].upper()+'/PDB/'+codechain[:4].lower()+'.pdb'
    AS = AtomSet.from_file(pdb_file)
    select_atoms = []
    for atom in AS.atoms:
        if atom.chain_id == codechain[4] and atom.res_num > ind -5 and atom.res_num < ind + 5:
            select_atoms.append(atom)
    new_AS = AtomSet(select_atoms)
    new_pdb = out_dir + 'HMMSTR codechain[:4].lower()+'.pdb'
