# blaster.py
# TLB 2/27/2021
# blaster.py help me run and parse command-line blast stuff
import subprocess
import shlex
import os
from Bio.Blast import NCBIWWW, NCBIXML 
from tqdm import tqdm
from os.path import isfile, join, exists
from config import * 
import argparse
def ff2name(f): # "../../fastas/1ubq_A.fasta" --> "1ubq_A"
    li = -1
    while True:
        try:
            li = f.index('/')
        except:
            break
        f = f[li+1:]
    # f == '1ubq_A.fasta'
    f = f[:len(f)-6]
    return f

program = "blastp"
seq_files = None 
out_dir = None 
outfmt = "4"
path = '~/Warehouse/Base/Externals/NR/nr'

def parse(blastfile, msa_dir):
    
    ind = blastfile.rfind('/')+1
    name = blastfile[ind:ind+4]
#               print("parsing " + f)
#               print("name = " + name)
    try:
        f_in = open(blastfile)
    except:
        print("failed to open blast file in parser.")
        return []
    first = True
    seqs = {}
    next_block = True
    block_size = 0
    block = 0
    q_len = 0
    m_len = 0
    chain = None
    msa_file = msa_dir + name + '.fa'
    f = open(msa_file, 'w')
    f_in = open(blastfile, 'r')
    for line in f_in:
        args = line.split()
        if len(args) > 0 and args[0] == "Query=":
            if len(list(seqs)) != 0 and chain is not None and next_block:
                out_f = msa_dir + name + chain +'.fa'
                f =  open(out_f, "w")
                for l in list(seqs):
#                    print(l)
#                    print(seqs[l])
                    f.write(">" + l + "\n")
                    f.write(seqs[l] + "\n")
#                           print("wrote to " + out_f )
            chain = line.strip()[len(line.strip()) -1 ]
            queries = []
            seqs = {}
        if first or next_block:
            if line[:len("Query_")] == "Query_":
    #            print(line)
                q_len = len(args[2])
                query = args[0]
                if query not in list(seqs):
                    first = True
                if first:
                    block_size = q_len
                    seqs[query] = args[2]
                else:
                    seqs[query] += args[2]
                first = False
                next_block = False
        else:
            if len(args) != 4:
                next_block = True
                block += 1
            else:
                if args[0] in list(seqs): # already added something to this entry
                    seqs[args[0]] += args[2]
                else:
                    m_len = len(args[2])
                    s = "-" * (block * block_size + (q_len) - (m_len))
                    s += args[2]
                    seqs[args[0]] = s
    if len(list(seqs)) != 0:
        out_f = msa_dir + name + chain + '.fa'
        f =  open(out_f, "w")

        for l in list(seqs):
#            print(l)
#            print(seqs[l])
            f.write(">" + l + "\n")
            f.write(seqs[l] + "\n")
        print("wrote to " + out_f)
        #os.system( "cat " + out_f)
    else:
        print("no queries found for " + name)
        return []

# 4 = flat query-anchored, no identities,        ---> Custom Parser
# 5 = XML Blast output,                          ---> Use NCBIXML
def run(seq_file, outfmt, program, out_dir):
    ff = seq_file
    this_name = ff2name(ff)
    save_file = out_dir + this_name + "." + program
    if isfile(save_file):
        print("file already exists, rewriting.")
    print("running " + program + " on fasta file: " + ff)
    print("wait 1~3 minutes")
    if outfmt is None:                            # BLASTP_DATABASE_DIRECTORY
        s = program + " -query " + ff + " -db " + path + " -out " + save_file + " -evalue 0.001"
    else:                                        # BLASTP_DATABASE_DIRECTORY
        s = program + " -query " + ff + " -db " + path + " -out " + save_file + " -outfmt " + str(outfmt) + " -evalue 0.001"
    os.system(s)
    print("saved " + self.program + " file: " + str(save_file))

def parse_args():
    parser = argparse.ArgumentParser(description = 'diagnostics.py')
    parser.add_argument('-if', dest = 'in_file', default = "None")
    parser.add_argument('-f', dest = 'format', default = "4")
    parser.add_argument('-p', dest = 'program', default = "blastp")
    parser.add_argument('-od', dest = 'outdir', default = "./tmp/")
    parser.add_argument('-e', dest = 'exe', default = "run")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    in_file = args.in_file # if exe = run then this is .fasta, 
    format_ = args.format           # if exe = parse then this is .blastp
    program = args.program
    outdir = args.outdir
    exe = args.exe
    if exe == 'run':
        run(in_file, format_, program, outdir) # read .fasta, write .blastp
    else:
        parse(in_file, outdir) # read .blastp, write .fa

    
