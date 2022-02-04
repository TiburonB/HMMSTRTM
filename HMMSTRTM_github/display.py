# display_test.py
# TLB 1/9/21 Testing new graphviz display function

from hmm import *
from hmm_viz import *
import argparse
def parse_args():
    parser = argparse.ArgumentParser(description = 'diagnostics.py')
    parser.add_argument('-mf', dest = 'model_file', default = "./HMMS/HMMSTR.hmm")
    parser.add_argument('-df', dest = 'drct_file', default = "None")
    parser.add_argument('-c', dest = 'cutoff', default = "0.1")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    model_file = args.model_file
    cutoff = float(args.cutoff)
    model_name = model_file[model_file.rfind('/')+1:]
    model_dir = model_file[:model_file.rfind('/')+1]
    HMMSTR = hmm(model_name, directory = model_dir)
    HMMSTR.read()
    print('displaying with cutoff = ' + str(cutoff))
    display(HMMSTR, cutoff = cutoff)



