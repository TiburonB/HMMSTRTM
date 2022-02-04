# path_view.py
# TLB 9/7/2021

import numpy as np
from graphviz import Digraph
import pydot
import os
import random
from tqdm import tqdm
import matplotlib
import matplotlib.pyplot as plt
from gamma import read as read_gamma
import math
from hmm import *
from os import listdir


def read(file_name):
    i_state = []
    j_state = []
    trans = []
    for l in open(file_name):
        args = l.split(' ')
        args = [ l.strip() for l in args if l != '']
        try:
            i_state.append(int(args[0]))
            j_state.append(int(args[1]))
            trans.append(float(args[2]))
        except:
            continue
    n = max([max(i_state), max(j_state)])
    return i_state, j_state, n, trans

def view(file_name):
    i_,j_,n,a = read(file_name)
    QUIET = True
    dot = Digraph(comment="NEW_TOP")

    for i in range(0, n):
        dot.node(str(i), str(i), shape = 'circle', color='red', style='filled')
    for i in range(0,len(i_)):
        dot.edge(str(i_[i]),str(j_[i]), label=str(a[i]))

    dot.render('path.gv', view=not(QUIET), format = 'png')



view('./1cisA00001.dag.path')

