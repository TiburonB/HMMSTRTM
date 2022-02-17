# hmm_viz.py
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

ss_2_color = { "H": 'blue', "E": 'cyan', "G": 'purple', "S": 'green', "T": 'red', "_": 'yellow'}
label2_color = { 'in': 'blue', 'ohelix': 'green', 'ihelix' : 'purple', 'outglob': 'pink', 'out': 'orange'}
    # Glycine = Green
    # Proline = Pink
    # Cysteine = Yellow
    # Hydrophobic = Grey
    # Aromatic = Purple
    # Polar (uncharged) = Teal
    # Charged (+) = Red
    # Charge (-) = Blue
aa_2_color = {'G': 'green',
        'P': 'pink',
        'C': 'yellow',
        'A': 'grey', 'V': 'grey', 'I': 'grey', 'L' : 'grey', 'M': 'grey',
        'F': 'purple', 'W': 'purple', 'Y': 'purple',
        'S': 'cyan', 'T': 'cyan', 'N': 'cyan', 'Q': 'cyan',
        'R': 'red', 'H' : 'red', 'K': 'red',
        'D': 'blue', 'E': 'blue' }
rama_2_shape = {"H": 'triangle', "G": 'triangle',
                "E": 'square', "e": 'square', "d": 'square',
                "L": 'star', "l": 'star',
                "x": 'polygon',
                "b": 'circle', "B": 'circle',
                "c": 'trapezium'}
TM_2_color = { "1": 'blue', "2": 'blue', "B": 'cyan', "H" : 'purple', "C": 'green', "I": 'red', "L": 'yellow', "F": 'orange', "U": 'pink' } 

def display_top_changes(HMM, new_nodes, new_trans, f_name = None):
    QUIET = True
    dot = Digraph(comment="NEW_TOP")

    for i in range(0, HMM.n):
        if i < HMM.u:
            dot.node(str(i), "tie node " + str(i), shape = 'doublecircle',color='blue')
        elif i < HMM.k + HMM.u:
            dot.node(str(i), "naught node " + str(i), shape = 'doublecircle')
        else:
            TM_type = HMM.emissions["TM"][HMM.node_data[i]["TM"].index(max(HMM.node_data[i]["TM"]))].strip()
            #print(rama_type) #HGBEdbeLlxc
            if i in new_nodes:
                dot.node(str(i), str(i), shape = 'circle', color='red', style='filled')
            else:
                dot.node(str(i), str(i), shape = 'circle', color='blue', style='filled')
    for i in range(0, HMM.n):
        for j in range(0, HMM.n):
            if HMM.a[i][j] > 0.0000:
                if (i,j) in new_trans:
                    dot.edge(str(i),str(j), color = 'red')
                else:
                    dot.edge(str(i),str(j))

    if f_name is not None:
        dot.render(f_name + '.gv', view=not(QUIET), format = 'png')
    else:
        dot.render('HMMSTR_NEW_TOP.gv', view=not(QUIET), format = 'png')

def highlight_states(HMM, states):
    dot = Digraph(comment = "highlighting")
    for i in range(0, HMM.n):
        if i == 0:
            dot.node(str(i), "naught node", shape = 'doublecircle')
        else:
            if i in states:
                dot.node(str(i), str(i), shape = 'circle', color='red', style='filled')
            else:
                dot.node(str(i), str(i), shape = 'circle', color='blue', style='filled')
    for i in range(0, HMM.n):
        for j in range(0, HMM.n):
            if HMM.a[i][j] > 0.0000:
                if (i,j) in new_trans:
                    dot.edge(str(i),str(j), color = 'red')
                else:
                    dot.edge(str(i),str(j))
    dot.render('HMMSTR_NEW_TOP.gv', view=not(True), format = 'png')
   
# visualize paradigm seq through looping over time interval in gamma...
def HMM_GAMMA_IMAGE(HMM, gamma, out_f):
    dot = Digraph(comment=HMM.name+'gamma')
    # find the max number of times a node is used, the max heat on any node
    max_num_times = 0
    max_heat = 0
    heats = []

    percentiles = [ float(i+1/10) for i in range(0, 9) ]
    for i in range(0, HMM.n):
        if False:#i == 0 :
            dot.node(str(i), "naught node", shape = 'doublecircle')
        else:
            
            color = 1
            for k,p in enumerate(percentiles):
                if gamma[i] > p:
                    color = k
            #color_map = {0: "gray0", 1: "blue", 2: "dodgerblue", 3: "aqua", 4: "cyan", 5: "darkolivegreen1", 6: "chartreuse", 7:"gold", 8:"chocolate1", 9:"firebrick1", 10:"darkred"}

            dot.node(str(i),  str(i),
                    shape = 'circle',
                    colorscheme = "purd9",
                    color = str(color+1), style = 'filled',
                    )

    for i in range(0, HMM.n):
        for j in range(0, HMM.n):
            if HMM.a[i][j] > 0.1000:
                dot.edge(str(i),str(j))
    
    print("RENDERING GAMMA GRAPH to " + out_f+'.jpg')
    dot.render(out_f , format = 'jpg')

def show(HMM):
    display_by_TM(HMM, QUIET = False)

def HMM_VIZ(HMM): # get graphviz drawings of current HMM and save to directory.
    display_by_SS(HMM)
    display_by_RAMA(HMM)
    display_by_TM(HMM)
    return

def display_by_SS(HMM, QUIET = True):

    dot = Digraph(comment="HMMSTR") #graph_type="digraph")

    for i in range(0, HMM.n):
        if i < HMM.u:
            dot.node(str(i), "tie node " + str(i), shape = 'doublecircle',color='blue')
        elif i < HMM.k + HMM.u:
            dot.node(str(i), "naught node " + str(i), shape = 'doublecircle')
        else:
            ss_type = HMM.emissions["SS"][HMM.node_data[i]["SS"].index(max(HMM.node_data[i]["SS"]))].strip()
            #print(ss_type)
            try:
                dot.node(str(i), str(i) + " " + ss_type, shape = 'circle', color=ss_2_color[ss_type], style='filled')
            except:
                dot.node(str(i), str(i) + " " + ss_type, shape = 'circle', color='black', style='filled')
    for i in range(0, HMM.n):
        for j in range(0, HMM.n):
            if HMM.a[i][j] > 0.0000:
                dot.edge(str(i),str(j))
    if HMM.directory is not None:
        dot.render(HMM.directory + HMM.name + '_SS.gv', view=not(QUIET), format = 'png')
    else:
        dot.render('HMMSTR_SS.gv', view=not(QUIET), format = 'png')

def display_states(HMM, cutoff = 0.1, QUIET = True):

    dot = Digraph(comment="HMMSTR") #graph_type="digraph")

    for i in range(0, HMM.n):
        tlabel = HMM.node_data[i]["label"][0]
        short_label = ""
        if tlabel.startswith('in'):
            short_label = 'in'
        elif tlabel.startswith('outglob') and tlabel != 'outglobSHORT':
            short_label = 'outglob'
        elif tlabel.startswith('out'):
            short_label = 'out'
        elif tlabel.startswith('ohelix'):
            short_label = 'ohelix'
        elif tlabel.startswith('ihelix'):
            short_label = 'ihelix'
        if i < HMM.u:
            dot.node(str(i), "tie node " + str(i), shape = 'doublecircle',color='blue')
        elif i < HMM.k + HMM.u:
            dot.node(str(i), "naught node " + str(i), shape = 'doublecircle')
        else:
            dot.node(str(HMM.node_data[i]["label"][0]), shape = 'circle', color = label2_color[short_label], style='filled')

    for i in range(0, HMM.n):
        for j in range(0, HMM.n):
            if HMM.a[i][j] > 0.0000:
                label_i = HMM.node_data[i]["label"][0]
                label_j = HMM.node_data[j]["label"][0]
                dot.edge(label_i,label_j)
    if HMM.directory is not None:
        dot.render(HMM.directory + HMM.name + '_RAMA.gv', view=not(QUIET), format = 'png')
    else:
        dot.render('HMMSTR_RAMA.gv', view=not(QUIET), format = 'png')


def display_by_RAMA(HMM, QUIET = True):

    dot = Digraph(comment="HMMSTR_RAMA") #graph_type="digraph")

    rama_2_color = {"H": 'blue', "G": 'blue', "E": 'cyan', "e": 'cyan', "d": 'cyan', "L": 'purple', "l": 'purple', "x": 'green', "b": 'red', "B": 'red', "c": 'yellow'}
    aa_2_shape = {"G": 'diamond', 'P': 'square'}

    for i in range(0, HMM.n):
        if i < HMM.u:
            dot.node(str(i), "tie node " + str(i), shape = 'doublecircle',color='blue')
        elif i < HMM.k + HMM.u:
            dot.node(str(i), "naught node " + str(i), shape = 'doublecircle')
        else:
            rama_type = HMM.emissions["RAMA"][HMM.node_data[i]["RAMA"].index(max(HMM.node_data[i]["RAMA"]))].strip()
            aa_type = HMM.emissions["AA"][HMM.node_data[i]["AA"].index(max(HMM.node_data[i]["AA"]))].strip()
            #print(rama_type) #HGBEdbeLlxc
            try:
                    dot.node(str(i), rama_type, shape = aa_2_shape[aa_type], color=rama_2_color[rama_type], style='filled')
            except:
                try:
                    dot.node(str(i), rama_type, shape = 'circle', color=rama_2_color[rama_type], style='filled')
                except:
                    dot.node(str(i), rama_type, shape = 'circle', color='black', style='filled')

    for i in range(0, HMM.n):
        for j in range(0, HMM.n):
            if HMM.a[i][j] > 0.0000:
                dot.edge(str(i),str(j))
    if HMM.directory is not None:
        dot.render(HMM.directory + HMM.name + '_RAMA.gv', view=not(QUIET), format = 'png')
    else:
        dot.render('HMMSTR_RAMA.gv', view=not(QUIET), format = 'png')
def display_by_RAMA_AA(HMM, cutoff = 0.1, QUIET = True):
    display_states_by_RAMA_AA(HMM, [ l for l in range(0, HMM.n)], cutoff = cutoff, QUIET = QUIET)
def display_states_by_RAMA_AA(HMM, states, cutoff = 0.1, QUIET = True):
    display_states(HMM, states, cutoff = cutoff, QUIET = QUIET, RAMA = True, AA = True)

def display_2_by_RAMA(HMM, HMM2, QUIET = True):
    dot = Digraph(comment="HMMSTR_RAMA") #graph_type="digraph")

    rama_2_color = {"H": 'blue', "G": 'blue', "E": 'cyan', "e": 'cyan', "d": 'cyan', "L": 'purple', "l": 'purple', "x": 'green', "b": 'red', "B": 'red', "c": 'yellow'}
    aa_2_shape = {"G": 'diamond', 'P': 'square'}
    subgraphs = {}
    for i in range(0, HMM.n):
        if i < HMM.u:
            dot.node(str(i), "tie node " + str(i), shape = 'doublecircle',color='blue')
        elif i < HMM.k + HMM.u:
            dot.node(str(i), "naught node " + str(i), shape = 'doublecircle')
        else:
            rama_type = HMM.emissions["RAMA"][HMM.node_data[i]["RAMA"].index(max(HMM.node_data[i]["RAMA"]))].strip()
            aa_type = HMM.emissions["AA"][HMM.node_data[i]["AA"].index(max(HMM.node_data[i]["AA"]))].strip()
            rama_type_2 = HMM2.emissions["RAMA"][HMM2.node_data[i]["RAMA"].index(max(HMM2.node_data[i]["RAMA"]))].strip()
            aa_type_2 = HMM2.emissions["AA"][HMM2.node_data[i]["AA"].index(max(HMM2.node_data[i]["AA"]))].strip()
            #print(rama_type) #HGBEdbeLlxc
            if rama_2_shape[rama_type] == rama_2_shape[rama_type_2] and aa_2_color[aa_type] == aa_2_color[aa_type_2]:
                try:
                        #dot.node(str(i), str(i), shape = aa_2_shape[aa_type], color=rama_2_color[rama_type], style='filled')
                        dot.node(str(i), str(i), shape = rama_2_shape[rama_type], color=aa_2_color[aa_type], style='filled')
                except:
                    try:
                        #dot.node(str(i), str(i), shape = 'circle', color=rama_2_color[rama_type], style='filled')
                        dot.node(str(i), str(i), shape = 'circle', color=aa_2_color[aa_type], style='filled')
                    except:
                        dot.node(str(i), str(i), shape = 'circle', color='black', style='filled')
            else:
                this_subgraph_name = 'cluster'+str(len(list(subgraphs)))
                subgraphs[str(i)] = this_subgraph_name
                with dot.subgraph(name = this_subgraph_name) as c:
                    c.attr(color='black')
                    try:
                        #c.node(str(i)+"_1", str(i)+"_1", shape = aa_2_shape[aa_type], color=rama_2_color[rama_type], style='filled')
                        c.node(str(i)+"_1", str(i)+"_1", shape = rama_2_shape[rama_type], color=aa_2_color[aa_type], style='filled')
                    except:
                        try:
                            #c.node(str(i)+"_1", str(i)+"_1", shape = 'circle', color=rama_2_color[rama_type], style='filled')
                            c.node(str(i)+"_1", str(i)+"_1", shape = 'circle', color=aa_2_color[aa_type], style='filled')
                        except:
                            c.node(str(i)+"_1", str(i)+"_1", shape = 'circle', color='black', style='filled')
                    try:
                        #c.node(str(i)+"_2", str(i)+"_2", shape = aa_2_shape[aa_type_2], color=rama_2_color[rama_type_2], style='filled')
                        c.node(str(i)+"_2", str(i)+"_2", shape = rama_2_shape[rama_type_2], color=aa_2_color[aa_type_2], style='filled')
                    except:
                        try:
                            #c.node(str(i)+"_2", str(i)+"_2", shape = 'circle', color=rama_2_color[rama_type_2], style='filled')
                            c.node(str(i)+"_2", str(i)+"_2", shape = 'circle', color=aa_2_color[aa_type_2], style='filled')
                        except:
                            c.node(str(i)+"_2", str(i)+"_2", shape = 'circle', color='black', style='filled')
                    c.edges([(str(i)+"_1", str(i)+"_2")])
                    c.attr(label='')
                    c.node_attr.update(style = 'filled')
    for i in range(0, HMM.n):
        for j in range(0, HMM.n):
            if HMM.a[i][j] > 0.1000 and HMM2.a[i][j] > 0.1000:
                n_i = str(i)
                if n_i in list(subgraphs):
                    n_i = n_i+"_1"#subgraphs[n_i]
                n_j= str(j)
                if n_j in list(subgraphs):
                    n_j = n_j+"_1"#subgraphs[n_j]
                dot.edge(str(n_i),str(n_j))
    if HMM.directory is not None:
        dot.render(HMM.directory + HMM.name + '_RAMA_COMP.gv', view=not(QUIET), format = 'png')
    else:
        dot.render('HMMSTR_RAMA.gv', view=not(QUIET), format = 'png')



def display_by_TM(HMM, cutoff = 0.1, QUIET = True, plot_p = True):
    display_states_by_TM(HMM, [l for l in range(0, HMM.n)], cutoff = cutoff, QUIET = QUIET, plot_p = plot_p)

def display_states_by_TM(HMM, states, cutoff = 0.1, QUIET = True, plot_p = True):
    display_states(HMM, states, cutoff = cutoff, QUIET = QUIET, TM = True)

# TLB 1/9/21
# Create a single display function that captures all others. Fuck This.


def display(HMM, QUIET = True, TM = True, AA = True, RAMA = True, cutoff = 0.1):
    display_states(HMM, [l for l in range(0, HMM.n)], QUIET = QUIET, TM = TM, AA = True, RAMA = True, cutoff = 0.1)

def display_states(HMM, states, QUIET = True, TM = False, AA = False, RAMA = False, cutoff = 0.1, 
        state_colors = None,  # dictionary of { state_number :  
        trans_features = None ):
    dot = Digraph(comment="HMSTR_RAMA")
    color_freqs = {}
    for A in list(aa_2_color):
        color_freqs[A] = 0
    shape_freqs = {}
    for A in list(rama_2_shape):
        shape_freqs[A] = 0
    for i in range(0, HMM.n):
        if i in states:
            if i < HMM.u:
                dot.node(str(i), "tie node " + str(i), shape = 'doublecircle',color='blue')
            elif i < HMM.k + HMM.u:
                dot.node(str(i), "naught node " + str(i), shape = 'doublecircle')
            else:
                TM_type = HMM.emissions["TM"][HMM.node_data[i]["TM"].index(max(HMM.node_data[i]["TM"]))].strip()
                rama_type = HMM.emissions["RAMA"][HMM.node_data[i]["RAMA"].index(max(HMM.node_data[i]["RAMA"]))].strip()
                aa_type = HMM.emissions["AA"][HMM.node_data[i]["AA"].index(max(HMM.node_data[i]["AA"]))].strip()
                if TM and not(RAMA) and not(AA):
                    label = HMM.node_data[i]["label"][0]
                else:
                    label = str(i)
                #try:
                #    label = HMM.node_data[i]["label"][0]
                #except:
                #    label = str(i)

                try:
                        #dot.node(str(i), str(i), shape = aa_2_shape[aa_type], color=rama_2_color[rama_type], style='filled')
                        if TM and AA and RAMA:
                            dot.node(str(i), label, shape = rama_2_shape[rama_type], color = TM_2_color[TM_type], fillcolor=aa_2_color[aa_type], style='filled')
                        elif AA and RAMA and not(TM):
                            dot.node(str(i), label, shape = rama_2_shape[rama_type], fillcolor=aa_2_color[aa_type], style='filled')
                        elif TM and not(RAMA) and not(AA):
                            dot.node(str(i), label, fillcolor=TM_2_color[TM_type], style='filled')
                        color_freqs[aa_type] +=1
                        shape_freqs[rama_type] +=1
                except:
                    raise
                    print("HERE")
                    try:
                        #dot.node(str(i), str(i), shape = 'circle', color=rama_2_color[rama_type], style='filled')
                        dot.node(str(i), label, shape = 'circle', color=aa_2_color[aa_type], style='filled')
                        shape_freqs["B"] +=1
                        color_freqs[aa_type] +=1
                    except:
                        dot.node(str(i), label, shape = 'circle', color='black', style='filled')
                        shape_freqs["B"] +=1
    for i in range(0, HMM.n):
        if i in states:
            for j in range(0, HMM.n):
                if j in states:
                    if HMM.a[i][j] > cutoff:
                        dot.edge(str(i),str(j))
    if HMM.directory is not None:
        dot.render(HMM.directory + HMM.name + '.gv', view=not(QUIET), format = 'png')
    else:
        dot.render('HMMSTR.gv', view=not(QUIET), format = 'png')

# gif display test
def display_state_motif_gif(HMM):
    cutoff = 0.1
    dot = Digraph(comment = 'html-gif-test')

    gifdir = './state_motifs/'
    gifs = [ gifdir + f for f in listdir(gifdir) if f.endswith('.gif')]
    color_freqs = {}
    for A in list(aa_2_color):
        color_freqs[A] = 0
    shape_freqs = {}
    for A in list(rama_2_shape):
        shape_freqs[A] = 0
    for i in range(0, HMM.n):
        if i < HMM.u:
            dot.node(str(i), "tie node " + str(i), shape = 'doublecircle',color='blue')
        elif i < HMM.k + HMM.u:
            dot.node(str(i), "naught node " + str(i), shape = 'doublecircle')
        else:
            TM_type = HMM.emissions["TM"][HMM.node_data[i]["TM"].index(max(HMM.node_data[i]["TM"]))].strip()
            rama_type = HMM.emissions["RAMA"][HMM.node_data[i]["RAMA"].index(max(HMM.node_data[i]["RAMA"]))].strip()
            aa_type = HMM.emissions["AA"][HMM.node_data[i]["AA"].index(max(HMM.node_data[i]["AA"]))].strip()
            label = str(i)
            #label = HMM.node_data[i]["label"][0]
            try:
                for f in gifs:
                    if f.split('_')[2] == str(i):
                        label = '<<TABLE> <TR><TD><IMG SRC= "'+f+'"/></TD></TR></TABLE>>'
            except:
                x = 0 #continue, label is the same as before

            try:
                    #dot.node(str(i), str(i), shape = aa_2_shape[aa_type], color=rama_2_color[rama_type], style='filled')
                    dot.node(str(i), label, shape = rama_2_shape[rama_type], color = TM_2_color[TM_type], fillcolor=aa_2_color[aa_type], style='filled')
            except:
                raise
                print("HERE")
                try:
                    #dot.node(str(i), str(i), shape = 'circle', color=rama_2_color[rama_type], style='filled')
                    dot.node(str(i), label, shape = 'circle', color=aa_2_color[aa_type], style='filled')
                    shape_freqs["B"] +=1
                    color_freqs[aa_type] +=1
                except:
                    dot.node(str(i), label, shape = 'circle', color='black', style='filled')
                    shape_freqs["B"] +=1
    for i in range(0, HMM.n):
        for j in range(0, HMM.n):
            if HMM.a[i][j] > cutoff:
                dot.edge(str(i),str(j))

    dot.render( 'test.gv', view=True, format = 'jpg')
    raise


