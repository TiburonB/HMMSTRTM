# gamma_heat.py
# TLB 7/15/21
# This program adds the gamma values across all entries in the db and s

from gamma import *
import argparse
from os import listdir
from tqdm import tqdm
from hmm_viz import * 
from hmm import * 

def parse_args():
    parser = argparse.ArgumentParser(description = 'gamma_heat.py')
    parser.add_argument('-gf', dest = 'gamma_file', default = "")
    parser.add_argument('-mf', dest = 'model_file', default = "")
    parser.add_argument('-of', dest = 'out_file', default = "")
    parser.add_argument('-s', dest = 'server', default = "False")
    args = parser.parse_args()
    return args

def str2bool(s):
    if s in ["True", 'T', "t", 'true']:
        return True
    return False

if __name__ == '__main__':
    args = parse_args()
    ftitr = args.gamma_file
    modelfile = args.model_file
    of = args.out_file
    server = str2bool(args.server)
    f = open(ftitr)
    gamma = []

    for l in f:
        gamma = [ float(x) for x in l.split() ]
    sum_gamma = sum(gamma)
    for i in range(0,len(gamma)):
        gamma[i] = gamma[i] / sum_gamma * 100
    
    HMM = hmm(modelfile)
    HMM.read()
    if server:
        HMM_GAMMA_IMAGE(HMM, gamma, './static/' + of)
    else:
        HMM_GAMMA_IMAGE(HMM, gamma, of)



