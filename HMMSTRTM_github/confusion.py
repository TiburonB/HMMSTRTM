from os import listdir
from tqdm import tqdm
import os
from os.path import isfile, isdir
import math
import argparse
SS_bins = [ "S/T", "G/H", "E" ]
RAMA_bins = [ "G/H", "L/l", "B/b", "E/e/d", "c", "x" ]



def write_to_file(filename, var):
    f = open(filename, 'w')
    s = "       "
    SS = False
    RAMA = False
    if len( var ) == 11:
        RAMA = True
    elif len(var) == 3:
        SS = True
    for i in range(0, len(var)):
        if SS:
            s += "{v:>10}".format(v=SS_bins[i])
        elif RAMA:
            s += "{v:>10}".format(v=RAMA_bins[i])
    s += "\n"
    for i in range(0, len(var)):
        if SS:
            s += "{v:>7}".format(v=SS_bins[i])
        elif RAMA:
            s += "{v:>7}".format(v=RAMA_bins[i])
        for j in range(0, len(var[i])):
            s += "{v:>10n}".format(v =var[i][j])
        s += "\n"
    f.write(s)
    #print(s)


def write_one(name, var):
    print(name)
    print("X = PREDICTION")
    print("Y = GROUNDTRUTH")
    s = "       "
    SS = False
    RAMA = False
    if name[:2] == "SS":
        SS = True
    elif name[:4] == "RAMA":
        RAMA = True
    for i in range(0, len(var)):
        if SS:
            s += "{v:>10}".format(v=SS_bins[i]) 
        elif RAMA:
            s += "{v:>10}".format(v=RAMA_bins[i]) 
    s += "\n"
    for i in range(0, len(var)):
        if SS:
            s += "{v:>7}".format(v=SS_bins[i])
        elif RAMA:
            s += "{v:>7}".format(v=RAMA_bins[i])
        for j in range(0, len(var[i])):
            s += "{v:>10n}".format(v =var[i][j])  
        s += "\n"
    print(s)


def writeout(AA_TEST, AA_TRAINING, SS_TEST, SS_TRAINING, RAMA_TEST, RAMA_TRAINING) : 
#    write_one("AA_TEST", AA_TEST)
#    write_one("AA_TRAINING", AA_TRAINING)
    write_one("SS_TEST", SS_TEST)
    write_one("SS_TRAINING", SS_TRAINING)
    write_one("RAMA_TEST", RAMA_TEST)
    write_one("RAMA_TRAINING", RAMA_TRAINING)


# 7/6/21
# Create Confusion matricies on training/test sets 
def confusion(model_file, gamma_dir, testing = False):
    #gamma_dir = model_file[:len(model_file)-4] + "_GAMMA/"
    if not isdir(gamma_dir):
        os.system("mkdir " + gamma_dir)
    #model_file = './HMMS/HMMSTR_ 30.hmm'
    #gamma_dir = './HMMS/HMMSTR_30_GAMMA/' # GAMMA DIR SHOULD CONTAIN ALL GAMMAS FROM final.drct (both training + test)
    gamma_files = [ f for f in listdir(gamma_dir) if isfile(gamma_dir + f) ]
    AA_TRAINING = (0,1)
    AA_TEST = (0,1)
    SS_TRAINING = [[0,0,0],[0,0,0],[0,0,0]]
    SS_TEST = [[0,0,0],[0,0,0],[0,0,0]]
    RAMA_TRAINING = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
    RAMA_TEST = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
    # FIRST... make sure the GAMMAS WERE PROPERLY CREATED...
    #os.system("gfortran -c hmm_io.f95 drct2profile.f95 hmmstrtm.f95 hmm_gamma.f95 hmmstr_gamma.f95")
    #os.system("gfortran hmm_io.o drct2profile.o  hmmstrtm.o hmm_gamma.o hmmstr_gamma.o")
    #os.system("./a.out \'" + model_file + "\' ./testing.drct \'" + gamma_dir + "\'")
    #os.system("./a.out \'" + model_file + "\' ./training.drct \'" + gamma_dir + "\'")
    #raise 
    gamma_files = [ f for f in listdir(gamma_dir) if isfile(gamma_dir + f) ]
    if len(gamma_files) == 0:
        print("NO GAMMA FILES FOUND, PLEASE GENERATE USING diagnostics.py")
        raise

    # LOAD THE CORRECT PROGRAM:
    os.system("gfortran -c hmm_io.f95 drct2profile.f95 hmmstrtm.f95 hmm_gamma.f95 confusion.f95")
    os.system("gfortran hmm_io.o drct2profile.o hmmstrtm.o hmm_gamma.o confusion.o")

    # GET LISTS OF TRAINING/TEST FILES:
    training_chains = [line[:4] for line in open('./training_chains.txt', "r")]
    test_chains = [line[:4] for line in open('./test_chains.txt', "r")]
    i = 0
    for f in tqdm(gamma_files):
        if i % 100 == 0 :
            writeout(AA_TEST, AA_TRAINING, SS_TEST, SS_TRAINING, RAMA_TEST, RAMA_TRAINING)

        if testing:
            s = "./a.out \'" + gamma_dir + '\' \'' + model_file + '\' testing.drct \'' + f + '\''
        else:
            s = "./a.out \'" + gamma_dir + '\' \'' + model_file + '\' final.drct \'' + f + '\''
        #print(s)
        data = os.popen(s).read().split()
        print(data)
        #print(len(data))
        #raise
        if len(data) != 45:
            continue
        if f in training_chains:
            for j in range(0, len(data)):
                if j < 9: # SS CONFUSION
                  i_pos = math.floor(j / len(SS_bins))
                  j_pos = j % len(SS_bins) 
                  SS_TRAINING[i_pos][j_pos] += int(data[j])
                else: # RAMA CONFUSION
                  i_pos = math.floor(( j - 9 ) / len(RAMA_bins))
                  j_pos = (j-9) % len(RAMA_bins)
                  RAMA_TRAINING[i_pos][j_pos] += int(data[j])
        elif f in test_chains:
            for j in range(0, len(data)):
                if j < 9: # SS CONFUSION
                  i_pos = math.floor(j / len(SS_bins))
                  j_pos = j % len(SS_bins) 
                  SS_TEST[i_pos][j_pos] += int(data[j])
                else: # RAMA CONFUSION
                  i_pos = math.floor(( j - 9 ) / len(RAMA_bins))
                  j_pos = (j-9) % len(RAMA_bins)
                  RAMA_TEST[i_pos][j_pos] += int(data[j])
        i += 1

    return writeout(AA_TEST, AA_TRAINING, SS_TEST, SS_TRAINING, RAMA_TEST, RAMA_TRAINING)
        

def parse_args():
    parser = argparse.ArgumentParser(description = 'confusion.py')
    parser.add_argument('-m', dest = 'model_file', default = "./HMMS/HMMSTR_  2.hmm")
    parser.add_argument('-d', dest = 'gamma_dir', default = "./HMMS/HMMSTR_2_GAMMA/")
    args = parser.parse_args()
    return args

# 7/14 CLI
if __name__ == '__main__':
    args = parse_args()
    model_file = args.model_file
    gamma_dir = args.gamma_dir
    confusion(model_file, gamma_dir)
