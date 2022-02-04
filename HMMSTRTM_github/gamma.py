# gamma.py
import os
import sys
from os.path import isfile
import numpy as np

# GAMMA I/O

def read(filename, directory):
    n = -1
    t = -1
    gamma = None
    print("READING GAMMAS FROM " + directory + filename)
    if isfile(directory + filename) == False:
        print("COULDN'T LOCATE GAMMA FILE", file = sys.stderr)
    else:
        for line in open(directory+filename):
            if n != -1 and t != -1 and gamma is None: # initialize gamma array
                gamma = np.zeros((t,n))
            comps = line.split()
            if comps[0] == "N":
                n = int(comps[1]) 
            elif comps[0] == "T":
                t = int(comps[1]) 
            else:
                titr = int(comps[0])-1
                nitr = int(comps[1])-1
                gamma[titr][nitr] = float(comps[2])
    print("GAMMA MATRIX (nxt) = (" + str(n) + "x" + str(t) + ").")
    print(gamma)
    return gamma

def write_gamma_profile(filename, directory, gamma_heat, model_file, drct_file):
    write_gamma_heat(filename, directory, gamma_heat, model_file, drct_file)
def write_gamma_heat(filename, directory, gamma_heat, model_file, drct_file):
    print("writing gammas to " + directory + filename)
    f = open(directory + filename, "w")
    f.write(model_file)
    f.write(drct_file)
    for l in gamma_heat:
        f.write(l)
    f.close()
    return

def read_gamma_profile(filename, directory):
    print('reading gamma profile from ' + directory + filename)
    f = open(directory + filename, 'r')
    lines = [ l for l in f ]
    gamma_prof = []
    for litr in range(0, len(lines)):
        if litr == 0:
            print('model_file = ' + lines[litr])
        elif litr == 1:
            print('drct_file = ' + lines[litr])
        else:
            gamma_prof.append(float(lines[litr].strip()))
    return gamma_prof

if __name__ == '__main__':
    read('7odc', './HMMSTR_50_GAMMA/')

