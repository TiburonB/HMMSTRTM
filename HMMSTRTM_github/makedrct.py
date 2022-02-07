# makedrct.py

from config import *
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import argparse
from pdb_handler import *

def parse_args():
    parser = argparse.ArgumentParser(description = 'makedrct.py')
    parser.add_argument('-q', dest = 'quiet', default = "True")
    parser.add_argument('-p', dest = 'pdbfile', default = "")
    parser.add_argument('-f', dest = 'profilefile', default = "")
    parser.add_argument('-c', dest = 'chain', default = "")
    args = parser.parse_args()
    return args

def str2bool(st):
    if st in ["True", "T", "t"]:
        return True
    return False


if __name__ == '__main__':
    dropped = []
    args = parse_args()
    QUIET = str2bool(args.quiet)
    pdbfile = args.pdbfile
    profilefile = args.profilefile
    chain = args.chain
    name = pdbfile[pdbfile.rfind('/')+1:pdbfile.rfind('/')+5] 
    drct_dir = pdbfile[:pdbfile.rfind('/')+1] 
    name_chain = name + chain
    ##############################################
    if not(QUIET):
        print("getting PH object.")
    PH = None 
    PH = pdb_handler(pdbfile,chain=chain)
    ##############################################
    if not(QUIET):
        print("getting sequence.")
    AAseq = None
    try:
        PH.get_fasta(this_chain = chain, rerun = True)
        AAseq = PH.sequence
    except:
        dropped.append(name)
        raise
    ##############################################
    if not(QUIET):
        print("getting SS sequence.")
    DSSP_SSseq = None
    DSSP_AAseq = None
    try:
        PH.get_dssp(chain = chain)
        DSSP_SSseq = PH.dssp
        DSSP_AAseq = PH.DSSP_AAseq 
    except:
        raise
        dropped.append(name)
    ##############################################
    if not(QUIET):
        print("getting PDBTM sequence, defaulting to globular.")
    PDBTM_TMseq = "1"*len(AAseq)
    PDBTM_AAseq = AAseq
    ##############################################
    if not(QUIET):
        print("getting PROFILE sequence")
    if not(QUIET):
        ##############################################
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("FASTA AA SEQUENCE = " + AAseq)
        print("DSSP AA SEQUENCE  = " + DSSP_AAseq)
        print("DSSP SS SEQUENCE  = " + DSSP_SSseq)
        print("PDBTM AA SEQUENCE = " + PDBTM_AAseq)
        print("PDBTM TM SEQUENCE = " + PDBTM_TMseq)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        ##############################################

    
    # aligns compseq to master seq and applies gap seq to retseq, returns retseq
    def align(masterseq, compseq, retseq, gap_char = '-'):
        aligned = pairwise2.align.globalxx(masterseq, compseq)
        gap_seq = aligned[0].seqB
        other_seq = aligned[0].seqA
        #print(gap_seq)
        #print(other_seq)
        if len(gap_seq) > len(masterseq):
            gap_seq = gap_seq[:len(masterseq)]
        aligned_ret = ""
        ind = 0
        for i in range(0, len(gap_seq)):
            #if gap_seq[i] == "-" and other_seq[i+1] == "-" and gap_seq[i+1] != "-":
            #    continue
            if other_seq[i] == "-" and i == 0:
                ind += 1
            elif other_seq[i] == "-":# and gap_seq[i-1] != "-":
                aligned_ret = aligned_ret[:len(aligned_ret)-1]
                if ind < len(retseq):
                    aligned_ret += retseq[ind]
                    ind += 1
            if ind >= len(retseq) or gap_seq[i] == "-":# and i != len(gap_seq) -1 and other_seq[i+1] != "-"):
                aligned_ret += gap_char
            else:
                aligned_ret += retseq[ind]
                ind += 1
        while len(aligned_ret) < len(masterseq):
            aligned_ret += gap_char
#            print(aligned_ret)
#            print(masterseq)
#            raise
        return aligned_ret
    try:
        DSSP_SSseq = align(AAseq, DSSP_AAseq, DSSP_SSseq, gap_char = '-')
        PDBTM_TMseq = align(AAseq, PDBTM_AAseq, PDBTM_TMseq, gap_char = 'U')
    except:
        dropped.append(name)
        raise
    seq = AAseq
    SSseq = DSSP_SSseq
    TMseq = PDBTM_TMseq
    ##################################
    print("PRINTING SEQUENCES")
    print("AA = " + str(seq))
    print("SS = " + str(SSseq))
    print("TM = " + str(TMseq))
    ##################################
    if seq is None or SSseq is None or TMseq is None:
        print("One of our sequences is None.")
        dropped.append(name)
        raise
    if len(seq) != len(SSseq) or len(SSseq) != len(TMseq):
        print("Sequences are of nonequal length.")
        raise 
    #################################
    if not(QUIET):
        print("Make drct.")
    drct_file = drct_dir + name_chain + ".drct"
    profile_dir = drct_dir
    s = "./xmaster2drct " + seq + " " + SSseq + " " + TMseq + " " + profile_dir + name_chain + '.profile' + " " + name_chain + " " + drct_file
    print(s)
    os.system(s)
    #################################
    if not(QUIET):
        print("Add ca coordinates to the drct.")
    s = "./drct3d_ca.pl " + drct_file + " " + pdbfile + " " + chain + " " + name
    os.system(s)
    #################################
    if not(QUIET):
        print("Add ramachandran emissions to the drct.")
    temp_file =  drct_file[:len(drct_file)-5] + "_1.drct"
    s = "./xrama2drct " + drct_file + " " + temp_file
    print(s)
    os.system(s)
    s = "mv " + temp_file + " " + drct_file
    os.system(s)
