# visualizer.py
# TLB 9/7/2021
# Create a visualization GUI type thing for HMMSTR using Graphviz

import argparse
from hmm import *
from hmm_viz import * 

def parse_args():
    parser = argparse.ArgumentParser(description = 'diagnostics.py')
    parser.add_argument('-mf', dest = 'model_file', default = "./HMMS/HMMSTR.hmm")
    parser.add_argument('-q', dest = 'quiet', default = "True")
    parser.add_argument('-s', dest = 'start', default = "-1")
    parser.add_argument('-e', dest = 'end', default = "-1")
    parser.add_argument('-c', dest = 'cutoff', default = "0.1")
    parser.add_argument('-g', dest = 'greed', default = "False")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    model_file = args.model_file
    greed = args.greed
    start = int(args.start)
    end = int(args.end)
    cutoff = float(args.cutoff)
    if greed == "True":
        greed = True
    else:
        greed = False
    model_dir  = model_file[:model_file.rfind('/')+1]
    model_name =  model_file[model_file.rfind('/')+1:model_file.rfind('.')].replace(" ", "")
    file_name =  model_file[model_file.rfind('/')+1:]
    try:
        model_itr = int(model_file[model_file.rfind('_')+1:model_file.rfind('.')])
    except:
        model_itr = 1
    diagnostics_dir = model_dir + model_name + "_DIAGNOSTICS/"
    gap_dir = model_dir + model_name + "_GAPS/"
    gamma_dir = model_dir + model_name + "_GAMMA/"

    HMMSTR = hmm(file_name, model_dir)
    HMMSTR.read()

    if start == -1 or end == -1:
        try:
            # Create a colored plot by argmax SS/RAMA/TM
            display_by_RAMA_AA(HMMSTR, cutoff = cutoff)
            display_by_TM(HMMSTR, cutoff = cutoff)
        except:
            display_states(HMMSTR, cutoff = cutoff)
    else:
        display_states_by_RAMA_AA(HMMSTR, [ l for l in range(start, end)], cutoff=cutoff)
        display_states_by_TM(HMMSTR, [ l for l in range(start, end)], cutoff= cutoff)
