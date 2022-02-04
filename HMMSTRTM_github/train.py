# train.py
# TLB 1/9/21 easy hmm-trainer.

from hmm import *
from hmm_tools import train
import argparse
def parse_args():
    parser = argparse.ArgumentParser(description = 'train.py')
    parser.add_argument('-mf', dest = 'model_file', default = "./HMMS/HMMSTR.hmm")
    parser.add_argument('-df', dest = 'drct_file', default = 'training.drct')
    parser.add_argument('-e', dest = 'epochs', default = 99)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    model_file = args.model_file
    drct_file = args.drct_file
    epochs = int(args.epochs)
    model_name = model_file[model_file.rfind('/')+1:]
    model_dir = model_file[:model_file.rfind('/')+1]
    HMMSTR = hmm(model_name, directory = model_dir)
    HMMSTR.read()
    train(HMMSTR, drct_file, epochs)
    print("completed training " + str(epochs) + " epochs of " + model_name + ".")






