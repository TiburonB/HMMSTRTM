# gamma_heat.py
# TLB 7/15/21
# This program adds the gamma values across all entries in the db and s

from gamma import *
import argparse
from os import listdir
from tqdm import tqdm

def gini(dat):
    G = 0
    for i in range(0, dat.shape[0]):
        for j in range(0, dat.shape[0]):
            G += abs(dat[i] - dat[j])
    G = G / ( 2 * dat.shape[0] * sum(dat))
    return G

def get_gamma_heat(gamma_dir):
    gamma_files = [ f for f in listdir(gamma_dir) ]
    heat = np.zeros(283)
    for f in tqdm(gamma_files):
        gamma = read(f, gamma_dir)
        for i in range(0, gamma.shape[0]):
            for j in range(0, gamma.shape[1]):
                heat[j] = heat[j] + gamma[i][j]
    return heat


def parse_args():
    parser = argparse.ArgumentParser(description = 'confusion.py')
   # parser.add_argument('-m', dest = 'model_file', default = "./HMMS/HMMSTR_ 60.hmm")
    parser.add_argument('-d', dest = 'gamma_dir', default = "./HMMS/HMMSTR_GAMMA/")
    args = parser.parse_args()
    return args

# 7/15 CLI
if __name__ == '__main__':
    args = parse_args()
   # model_file = args.model_file
    gamma_dir = args.gamma_dir
    heat = get_gamma_heat(gamma_dir) 
    print(heat)
    print("GINI COEF = " + str(gini(heat)))
