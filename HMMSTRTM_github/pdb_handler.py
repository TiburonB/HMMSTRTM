# pdb_handler.py

AALIST = ["ALA", "ARG", "PHE", "PRO", "TYR", "TRP", "GLY", "LEU", "ILE", "GLU", "GLN", "ASP", "ASN", "SER", "CYS", "MET", "THR", "LYS", "VAL", "HIS", "MSE"]
DSSP_AA1LIST = ["A", "R", "F", "P", "Y", "W", "G", "L", "I", "E", "Q", "D", "N", "S", "C", "M", "T", "K", "V", "H", "d", "e", "f", "c", "b", "a"]
import sys
from os import listdir
from Bio.Blast import NCBIWWW, NCBIXML  
from os.path import isfile, join, exists
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
import subprocess
import shlex
import os
import re
import numpy as np
from blaster import *
from os.path import exists
from Atom import *
# used in get_fasta:
AA3to1={
   'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
   'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
   'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
   'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
# NCAA's   
   'MSE':'M', 'MHS':'H', 'AGM':'R', 'GL3':'G', 'OCS':'C', 'SMC':'C', 'MLY':'K', 'DMG':'G', 'SEP':'S', 'PCA':'Q',
   '4NT':'S', 'HYP':'P', 'NH2':'_', 'TRO':'W', 'KCX':'K', 'CME':'C', 'ALY':'K', 'CSD':'C', 'TYS':'Y', 'MLZ':'K',
   'M3L':'K', 'ACE':'_', 'STA':'_', 'MSO':'M', 'SME':'M', 'UNK':'_', 'CSS':'C', 'SNN':'N', 'ACY':'_', 'ORP':'_',
   '5JP':'S', 'OMT':'M', '4BF':'F', 'GLZ':'_', 'TPO':'T', 'PTR':'Y', 'NLE':'L', '2ML':'L', 'DPN':'F', 'SEE':'S',
   '2TL':'T', 'O12':'_', 'LMT':'_', 'CSO':'C', '1MA':'A', 'OMU':'_', 'OMG':'_', 'UR3':'_', 'PSU':'_', 'ACA':'_', 
   'BFD':'D', 'MEQ':'Q', 'LLP':'K', 'SAR':'_', 'MVA':'V', 'DSN':'S', 'N2C':'C', 'NCY':'C', 'SAH':'C', 'PGE':'_',
   'QUI':'_', 'ACT':'_', '15P':'_', 'C37':'_', '5CM':'_', 'CAF':'C', '8AN':'_', 'CRO':'G', '2RA':'A', 'FGL':'_',
   'DVA':'V', 'R2T':'Q', '2TL':'T', '0QZ':'S', 'MB8':'_', 'KFP':'K', 'LP6':'K', 'AYE':'_', '5IU':'_', 'CY3':'C',
   'OXX':'D', '2MR':'R', 'GYS':'_', 'LYR':'K', '73C':'_', '48V':'_', 'PO4':'_', 'NEP':'H', 'CGU':'E', '3DR':'_',
   'HOH':'_', 'TYN':'Y', 'KYQ':'K', 'TPQ':'Y', 'KCR':'K', 'BE2':'_', 'CU1':'_', 'CAS':'C', '9QV':'_', 'AIB':'L',
   'PHD':'D', '8OG':'_', 'YCM':'C', 'DAL':'A', 'ACB':'D', 'FGA':'E', 'DAM':'A', '1ZN':'_', 'SEB':'S', 'DSE':'S',
   'M12':'_', '5PG':'_', 'DLY':'K', 'MED':'M', 'DGN':'Q', 'DAR':'R', 'DPR':'P', 'DPN':'F', 'DCY':'C', 'DGL':'E',
   'DIL':'I', 'DTY':'Y', 'DSG':'N', 'DTH':'T', 'DAS':'D', 'DHI':'H', 'NAG':'_', 'SO4':'_', 'CIT':'_', 'DLE':'L',
   'DTR':'W', 'PHI':'F', 'CSX':'C', 'CRQ':'_', 'FC0':'_', 'OAR':'R', 'BIF':'F', 'TH5':'T', 'MA7':'_', '2YR':'_',
   'HPD':'_', 'AHB':'N', 'ALS':'A', '6V1':'_', 'IML':'I', '6VO':'_', 'C12':'G', 'MGN':'Q', 'DYA':'D', 'CSU':'C',
   'P5P':'A', 'RGL':'R', 'GTP':'_', '3PA':'_', 'GGL':'E', '4BA':'_', 'MLE':'L', 'BMT':'T', 'ABA':'_', 'DOC':'_',
   'BRU':'_', 'FME':'M', 'HIC':'H', 'PHQ':'_', 'CF0':'_', 'ASJ':'_', 'PYL':'K', 'CR2':'_', 'MHO':'M', 'TA4':'_', 
   '9AT':'T', 'FHU':'_', '5CS':'C', 'AYA':'A', 'CXM':'M', 'NRQ':'_', '01B':'_', 'ASA':'D', 'PYR':'_', 'A23':'_', 
   'MYR':'_', 'QBT':'_', '2L6':'_', 'AME':'M', 'SNC':'C', 'EDA':'_', 'PPU':'_', 'HFA':'_', 'BTN':'_', 'ALO':'T',
   '2X0':'_', '4L8':'_', 'KYN':'W', '2DT':'_', 'NMM':'_', 'B74':'_', 'APD':'F', 'FTR':'W', '0A1':'Y', '3WX':'P',
   'GOL':'_', 'LYK':'K', '0QE':'_', 'OAS':'S', 'SCY':'C', 'SCS':'C', 'MDO':'_', 'SCH':'C', 'UF0':'_', 'MPT':'_',
   'U2X':'Y', 'HSO':'H', 'MYK':'K', 'CKC':'_', 'IAS':'D', 'JJJ':'C', '5AA':'_', 'DKA':'_', 'AR7':'_', 'MEA':'F',
   'ALN':'A', 'CIR':'R', '4LT':'_', 'ATP':'_', 'TED':'_', '3X9':'A', '9R1':'_', 'CYG':'C', 'DHL':'_', 'BHD':'D', 
   'QCS':'C', 'NIY':'Y', 'GMA':'_', '1AP':'_', 'PXU':'P', 'OMC':'_', '2CO':'C', 'HCS':'_', 'ALC':'_', 'SAC':'S',
   'JJK':'C', 'OCY':'C', 'DBB':'_', 'DHA':'A', 'TEE':'_', 'MLL':'L', '3CN':'_', 'TRQ':'W', 'APK':'K', 'BOC':'_', 
   'DDZ':'C', '2ZC':'S', '5R0':'_', 'ELY':'K', '5R5':'_', 'ORN':'_', 'ZAL':'A', 'CR8':'_', 'TGP':'_', 'HS8':'H',
   'DAH':'F', 'TY2':'Y', 'DCZ':'_', 'PO2':'_', 'PED':'_', 'PG1':'S', 'C38':'_', 'KPI':'_', '6KM':'C', 'QPA':'C',
   'MEN':'N', 'UD5':'_', 'SEC':'S', 'SVY':'S', 'CS4':'C', 'CCC':'_', '5CY':'_', 'TS6':'_', 'TTD':'_', 'CYW':'C',
   'CSR':'C', '1CC':'_', 'HRG':'R', '2OP':'_', 'MSU':'_', 'ALV':'A', 'SIC':'_', 'FBP':'_', 'PHL':'_', 'OSE':'S',
   'BAL':'A', '4WQ':'_', 'DA2':'R', 'AAR':'R', 'SE7':'_', 'ZPO':'_', 'XZA':'_', '2PR':'_', '64T':'_', '5PY':'_',
   'UMA':'A', '5BU':'_', 'TY5':'Y', 'ABN':'_', 'RDG':'_', '2L5':'F', 'R0E':'_', 'MOH':'_', 'MUM':'_', '1WA':'_',
   'IGU':'_', 'JSP':'_', '1W5':'_', 'HZP':'P', 'CAR':'_', 'MK8':'L', 'CCS':'C', 'OHI':'A', 'SNM':'S', 'HSK':'H',
   'LED':'L', '6FL':'L', 'NRI':'_', 'KWS':'_', 'MP8':'P', 'FP9':'P', '9KK':'L', 'FOX':'_', '6CW':'W', 'IYR':'Y',
   'XYG':'_', 'TYC':'Y', 'NAL':'_', 'TYI':'Y', 'OIC':'_', '200':'A', 'FGP':'_', 'MIR':'S', 'HSE':'S', '2MG':'_',
   'H2U':'_', 'M2G':'_', '6IA':'_', '5MC':'_', '5MU':'_', '1HB':'_', 'ASL':'_', '9FZ':'_', '9G2':'_', 'B3K':'K',
   'B3D':'D', 'B3E':'E', 'B3Y':'Y', 'B3A':'A', 'B3S':'S', 'B3X':'N', 'BEZ':'_', 'M9P':'_', '0AF':'W', 'CMT':'C',
   '6MA':'_', 'TRW':'W', 'MAA':'A', 'CMH':'C', 'C34':'_', 'PFX':'_', 'GLX':'G', 'CSP':'C', 'ONL':'L', 'CBR':'_', 
   'PSW':'A', 'LYZ':'K', 'TSE':'_', '0QL':'A', 'MLU':'L', 'OMZ':'Y', 'GHP':'_', 'OMY':'Y', '3FG':'_', 'TYD':'_',
   'BGC':'_', '6OG':'_', '3DA':'_', 'FAK':'K', 'BTK':'K', 'DAB':'_', 'SLZ':'K', '6MI':'_', 'IIL':'I', '5BV':'_',
   '8AZ':'_', '54C':'_', 'WLU':'L', 'WPA':'F', 'WVL':'_', 'ASX':'X', 'M2L':'_', 'APC':'_', 'OTZ':'_', 'BP4':'_', 
   'PEA':'_', 'N7P':'P', 'BH2':'N', 'GDP':'_', 'LIG':'_', '6MZ':'_', 'YCP':'_', 'TLX':'_', '4SU':'_', 'QUO':'_',
   'G7M':'_', 'HOX':'F', 'IVA':'_', 'DDE':'_', 'NEH':'_', '7MG':'_', 'PHA':'_', 'GLK':'_', 'E1X':'_', 'TTQ':'W',
   'GME':'E', 'MHW':'_', 'MHV':'_', '004':'_', 'GVL':'S', 'DPP':'_', 'VAD':'_', 'DVA':'V', 'FMU':'_', 'CR7':'_',
   'C49':'_', 'VKJ':'_', 'LPD':'_', 'XPC':'_', 'GM8':'_', 'LCK':'K', 'BF9':'_', 'PIV':'_', 'LPL':'L', 'YRR':'_',
   '1QQ':'_', 'R4K':'_', 'AKK':'_', '143':'C', 'TCP':'_', '6FK':'_', 'PYO':'_', 'CGN':'_', 'HL2':'_', 'OTH':'T',
   'THC':'_', 'HF2':'_', '5CT':'K', 'DBZ':'A', 'UBD':'_', 'NFA':'F', 'RC7':'_', 'CZ2':'C', 'OPR':'R', 'ECC':'_',          
   '4L0':'P', '56A':'H', 'XCP':'_', '2KT':'_', 'DBU':'_', 'NA8':'A', 'TH6':'T', 'PRS':'P', 'CPI':'_', '5OM':'_',
   '2LF':'_', '1OL':'_', 'E1H':'_', 'CWR':'_', 'ARO':'R', 'CSZ':'C', 'LYX':'K', 'RSQ':'_', '23G':'_', 'HPH':'_', 
   'DYJ':'P', 'CYD':'_', 'HSL':'_', 'N7X':'_', 'T39':'_', 'C5L':'_', 'A2M':'_', 'PST':'_', '6OO':'_', 'RFJ':'_', 
   '6NW':'_', 'CSA':'C', '0CS':'A', '02J':'_', 'PJE':'_', '010':'_', 'SUI':'_', 'YOF':'Y', 'CSK':'C', '00C':'C',
   '1MG':'_', 'CSJ':'C', 'MME':'M', '6L3':'_', '2LT':'Y', 'KOR':'C', 'PSA':'_', 'GTA':'_', 'UPS':'_', 'GVE':'_', 
   'GCP':'_', 'KI2':'_', 'GPL':'K', '5HC':'_', 'MX5':'_', 'PR4':'P', 'CRK':'_', 'VLM':'_', 'MG1':'_', 'DIV':'V',
   'LEN':'_', '0YG':'G', 'END':'',  'DYG':'_', 'PVO':'_', 'R1A':'_', '4HH':'S', 'PFF':'F', '3MY':'Y', 'PPN':'F',
   'XCY':'_', 'PAE':'_', 'FOE':'_', 'HQA':'A', 'WRP':'W', 'NVA':'_', 'P1L':'C', 'AC5':'_', '23F':'_', '4M9':'_', 
   'LEF':'L', '4H0':'_', 'CRF':'_', 'CYQ':'C', 'PN2':'_',     
 }


class pdb_handler:
    def __init__(self, pdbfile, chain = None):
        self.f = pdbfile
        self.name = pdbfile[pdbfile.rfind('/')+1:pdbfile.rfind('/')+5].lower()
        self.chain = pdbfile[pdbfile.rfind('/')+6]
        self.chain = chain # if chain is None then load all chains
        # unspecified to begin with, can be loaded in at will:
        self.structure = None # get_bpstruct
        self.Atoms = None # get_atoms
        self.Atom_Set = None
        self.recurse = False

    def get_bpstruct(self):
        try:
            #print(self.f)
            parser = PDBParser(QUIET = True)
            self.structure = parser.get_structure(self.name, self.f)
            return self.structure
        except:
            print("ERROR: failed to read file.", file = sys.stderr)
            return None

    # ippds = inter proton-proton (pairwise) distances
    # ( used to gather corpus 6/19 ippds range = 5 Angstroms )
    def get_ippds(self):
        atoms = AtomSet.from_file(self.f)
        self.iads = atoms.get_iads_1d(chain = self.chain)
        return self.iads

    def get_seq(self, chain = None):
        return self.get_sequence(chain = chain).strip()

    def get_sequence(self, chain = None):
        if chain is not None:
            try:
                from Bio.PDB.DSSP import dssp_dict_from_pdb_file
                dssp_tuple = dssp_dict_from_pdb_file(self.f)
                dssp_dict = dssp_tuple[0]
                dict_keys = dssp_tuple[1]
            except:
                print("Unexpected error:", sys.exc_info()[0])
                print(self.f)
                if self.recurse:
                    return None
                self.rsyncpdb(rerun = True)
                self.recurse = True
                return self.get_sequence(chain = chain)
            seq = ""
            for key in dict_keys:
                if key[0] == chain:
                    #print(key)
                    dssp_obj = dssp_dict[key]
                    seq += dssp_obj[0]
            self.seqence = seq
            return seq


    # get a dictionary for sequence
    #   { CHAIN : "SEQUENCE" }
    def get_fasta_rec(self):
        # requires fasta

        if self.fasta_files is None: # if we have no fasta files, then get them
            self.get_fasta()
        try:
            from Bio.PDB.DSSP import dssp_dict_from_pdb_file
            dssp_tuple = dssp_dict_from_pdb_file(self.f)
            dssp_dict = dssp_tuple[0]
            dict_keys = dssp_tuple[1]
        except:
            print("Unexpected error:", sys.exc_info()[0])
            print(self.f)
            if self.recurse:
                return None
            self.rsyncpdb(rerun = True)
            self.recurse = True
            return self.get_dssp(chain = chain)
        seq = {}
        for key in dict_keys:
            if key[0] not in list(seq):
                seq[key[0]] = ""
            dssp_obj = dssp_dict[key]
            seq[key[0]] += dssp_obj[0]
        self.fasta_rec = seq
        return self.fasta_rec
    
    # check if we've got the pdb file, otherwise generate it
    # check if we've already generated the fasta files, otherwise generate them
    # set fast
    def get_fasta(self, rerun = False, this_chain = None):
       
        try:
            from Bio.PDB.DSSP import dssp_dict_from_pdb_file
            dssp_tuple = dssp_dict_from_pdb_file(self.f)
            dssp_dict = dssp_tuple[0]
            dict_keys = dssp_tuple[1]
        except:
            print("Unexpected error:", sys.exc_info()[0])
            print(self.f)
            return 
        self.seq_1 = {}
        for key in dict_keys:
            if key[0] not in list(self.seq_1):
                self.seq_1[key[0]] = ""
            dssp_obj = dssp_dict[key]
            self.seq_1[key[0]] += dssp_obj[0]
        s = ""
        for k in list(self.seq_1):
            if this_chain is None or (this_chain is not None and k == this_chain):
                s += "> " + self.name + ":" + k + "\n" + self.seq_1[k] + "\n"
            #print("> " + self.name + ":" + k + "\n" + self.seq_1[k])
        if this_chain is None:
            fasta_file = fdir + self.name + ".fasta"
        else:
            try:
                self.sequence = self.seq_1[this_chain]
            except: # try again
                self.rsyncpdb(rerun = True)
                return self.get_fasta(this_chain = this_chain, rerun = True)
            #fasta_file = fdir + self.name + this_chain +".fasta"
        #self.fasta_files.append(fasta_file)
        #f = open(fasta_file, "w")
        #f.write(s)
        #f.close()
        #this_chain = 0
        #self.entry.log("Wrote FASTA file.") 
        if this_chain is not None:
            self.sequence = self.seq_1[this_chain]
            return self.seq_1[this_chain]
        else:
            if len(list(self.seq_1)) == 1:
                self.sequence = self.seq_1[list(self.seq_1)[0]]
                return self.seq_1[list(self.seq_1)[0]]
            else:
                print(self.seq_1)
                print("SPECIFY CHAIN TO RETURN FROM get_fasta() in pdb_handler.")


    def get_dssp(self, chain = None): # if chain is specified then look in dssp_dir


        try:
            s = self.structure
        except:
            self.get_bpstruct()
        if self.structure is None:
            self.get_bpstruct()
            if self.structure is None:
                print("Problem generating structure.")
                return
        dssp = ""
        if chain is not None:
            try:
                from Bio.PDB.DSSP import dssp_dict_from_pdb_file
                dssp_tuple = dssp_dict_from_pdb_file(self.f)
                dssp_dict = dssp_tuple[0]
                dict_keys = dssp_tuple[1]
            except:
                print("Unexpected error:", sys.exc_info()[0])
                print(self.f)
                if self.recurse:
                    return None
                self.rsyncpdb(rerun = True)
                self.recurse = True
                return self.get_dssp(chain = chain)
            seq = ""
            for key in dict_keys:
                if key[0] == chain:
                    #print(key)
                    dssp_obj = dssp_dict[key]
        #            print(dssp_obj)
                    #if dssp_obj[0] not in DSSP_AA1LIST: # Found either a weird gap, or a NCAA.
                    #    seq += dssp_obj[0]
                    #    #print(dssp_obj[0])
                    #    continue
                    dssp_char = dssp_obj[1]
                    if dssp == ' ' :
                        dssp += '-'
                    seq += dssp_obj[0]
                    dssp += dssp_char
            #print("seq  = " + seq)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("dssp reported sequence  = " + seq)
        print("dssp SS sequence        = " + dssp)
        self.sequence = self.get_seq(chain = chain).strip()
        print("fasta reported sequence = " + self.sequence)
        from Bio import pairwise2
        from Bio.pairwise2 import format_alignment
        align = pairwise2.align.globalxx(self.sequence, seq)
        gap_seq = align[0].seqB
        if len(gap_seq) > len(self.sequence):
            gap_seq = gap_seq[:len(self.sequence)]
        print("aligned AA seq          = " + gap_seq)
        aligned_dssp = ''
        ind = 0
        for i in range(0, len(gap_seq)):
            if ind >= len(dssp) or gap_seq[i] == '-':
                aligned_dssp += '-'
            else:
                aligned_dssp += dssp[ind]
                ind += 1
        dssp = aligned_dssp
        print("aligned dssp seq        = " + dssp)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        while len(dssp) < len(self.sequence):
            dssp += '-'
        #raise
        self.dssp = dssp
        self.DSSP_AAseq = seq
        #f = dssp_dir + self.name + chain
        #f_out = open(f, "w")
        #f_out.write(self.dssp)
        #f_out.close()
        return self.dssp

    def get_umask(self, chain = None): # get mask of residues with disordered structure
        self.get_pdb_files()
        f = self.base_pdb_file
        seq = self.get_fasta(this_chain = chain).strip()
        ATOM_SET = AtomSet.from_file(f)
        indeces = []
        for atom in ATOM_SET.atoms:
            #print(atom.chain_id)
            if atom.chain_id == chain and atom.res_num not in indeces:
                indeces.append(atom.res_num)
        UMASK = ""
        for k in range(0, len(seq)):
            if k+1 not in indeces:
                UMASK += "1"
            else:
                UMASK += "0"
        return UMASK

