# hmm_tools.py
# build a set of tools that can be used to craft HMM topology algorithmically

from hmm import * 
from datamine import * 
import os
from tplot import *
############# HMMSTR INTERACTIONS ################

def train(CANON, database, epochs):
    # TRAIN A VERSION OF HMMSTRTM:
    os.system("gfortran -c -fdefault-double-8 hmm_io.f95 drct2profile.f95 gamma_io.f95 hmm.f95 hmm_train.f95 -fcheck=all")
    os.system("gfortran hmm_io.o drct2profile.o gamma_io.o hmm.o hmm_train.o")
    os.system("./a.out '" + CANON.file + "' " + database + " " +  str(epochs))
    return True 

def scratch(CANON, database, code):
    # SCRATCH(RUN) A VERSION OF HMMSTRTM:
    os.system("gfortran -c -fdefault-double-8 hmm_io.f95 drct2profile.f95 gamma_io.f95 hmm.f95 hmm_run.f95 -fcheck=all")
    os.system("gfortran hmm_io.o drct2profile.o gamma_io.o hmm.o hmm_run.o")
    os.system("./a.out '" +  CANON.file + "' " + database + " " + code)
    return True

def examine(CANON, database, code = None):
    if code is not None: 
        os.system("python3 comprehensive.py -mf '" + CANON.file + "' -df " + database + " -code " + code )
    else:
        os.system("python3 comprehensive.py -mf '" + CANON.file + "' -df " + database )

# evaluate the model, gamma-weighted but run f/b at all start/end time pairs.
def exhaust(CANON, database, code = None):
    os.system("gfortran -c -fdefault-double-8 hmm_io.f95 drct2profile.f95 gamma_io.f95 hmm.f95 hmm_exhaust.f95")
    os.system("gfortran hmm_io.o drct2profile.o gamma_io.o hmm.o hmm_exhaust.o")
    os.system("./a.out '" +  CANON.file + "' " + database + " " + code)
    return True

###################################################

# TOPOLOGY BUILDING ALGORITHMS ###########################
# All topologies should have an additional datastructure which describes topology in terms of:
#    ... A dictionary   { Parent : [(Child, outtrans)] }
#    ... A table   [ 1, 2, 3, 4, 5]
#      + trans     .a[1][2] = x
#    pyramid architecture is represented by a 2d table where 
#       row x column y is the label of the state in level x diagonal y in the pyramid
#   This pyramid has width = 3, states are represented by '#', no info exists in 'X'.
#  [[ # # # ]
#   [ X # # ]
#   [ X X # ]]

def get_pyramid(width, _name_, sink):
    # given diagonal and level, return label
    def DL2l(D,L):
        return "L_d"+_name_+str(D)+"_"+str(L)

    # pyramid_state_table
    pyramid_states = [["0" for d in range(0,width)] for l in range(0,width)]
    for i,x in enumerate(pyramid_states):
        for j,y in enumerate(x):
            if i - j <= 0:
                pyramid_states[i][j] = DL2l(j,i)
    #print(pyramid_states)
    
    # pyramid_state_trans_dict
    pyramid_trans = {}
    for x in pyramid_states:
        for y in x:
            if y != "0":
                pyramid_trans[y] = []
    
    # given the label of the state in pyramid, 
    #   return the Level of that state
    def l2L(l):
        s = l.find(_name_)+len(_name_)
        s = l.index("_",s)+1 # skip first field (Diagonal)
        try:
            ret = int(l[s:l.index("_",s)])
        except:
            ret = int(l[s:len(l)])
        return ret
    # given the label of the state in pyramid, 
    #   return the Diagonal of that state
    def l2D(l):
        s = l.find(_name_)+len(_name_)
        return int(l[s:l.index("_",s)])

    for y in list(pyramid_trans):
        L = l2L(y)
        D = l2D(y)
        # add transition L, D --> L + 1, D + 1 
        j = DL2l(D+1, L+1)
        if j in list(pyramid_trans):
            pyramid_trans[y].append((j, 1))
        # if L is even, add transition L, D -->  L, D + 1
        if L % 2 == 0:
            j = DL2l(D+1, L)
            if j in list(pyramid_trans):
                pyramid_trans[y].append((j, 1))
        # if we haven't added a transition yet, then add transition L, D --> Sink
        if len(pyramid_trans[y]) == 0:
            pyramid_trans[y].append((sink,1))
    #print(pyramid_trans)
    return pyramid_trans


def add_pyramid_cascade(CANON, node_data, width, _name_, E_name):
    pyramid_trans = get_pyramid(width,_name_,E_name)
    pyramid_states = [ l for l in list(pyramid_trans) ]
    for state in pyramid_states:
        node = {}
        for l in list(node_data):
            node[l] = [ k for k in node_data[l]]
        node['label'] = [state]
        CANON.add_node(node, prior = 0)
    for state in list(pyramid_trans):
        i = CANON.get_node_by_label(state)
        for t in pyramid_trans[state]:
            j = CANON.get_node_by_label(t[0])
            CANON.a[i][j] = t[1]


# make l states 
# l states are linearly connected states
# l states = linear connected states, 1 goes to 2 goes to 3, before states, after states
def make_l_states(CANON, l, node_data, before_states, after_states, _name_, label_offset = 0, bstate_trans = -1 ):
    for i in range(0, l):
        node = { }
        for k in list(node_data):
            node[k] = [ j for j in node_data[k]]
        node["label"] = ["l"+_name_+str(i+1+label_offset)]
        print(node["label"])
        CANON.add_node(node, prior = 0)
        if i == 0 :
            if bstate_trans == -1:
                for t in before_states:
                    CANON.a[t][CANON.n-1] = 1 / len(before_states)
            else:
                for t in before_states:
                    print(CANON.node_data[t]['label'][0] + " --> " + CANON.node_data[CANON.n-1]['label'][0] + " = " + str(bstate_trans) + ".")
                    CANON.a[t][CANON.n-1] = bstate_trans
        if i == l - 1:
            for j in after_states:
                CANON.a[CANON.n-1][j] = 1 / len(after_states)
        if i > 0:
            CANON.a[CANON.n-2][CANON.n-1] = 1
    return CANON


# ARCHIVED

# make l states 
# l states are linearly connected states
# l states = linear connected states, 1 goes to 2 goes to 3, before states, after states
#def make_l_states(CANON, l, node_data, before_states, after_states, _name_, label_offset = 0, bstate_trans = -1 ):
#    for i in range(0, l):
#        node = { }
#        for k in list(node_data):
#            node[k] = [ j for j in node_data[k]]
#        node["label"] = ["l"+_name_+str(i+1+label_offset)]
#        print(node["label"])
#        CANON.add_node(node, prior = 0)
#        if i == 0 :
#            if bstate_trans == -1:
#                for t in before_states:
#                    CANON.a[t][CANON.n-1] = 1 / len(before_states)
#            else:
#                for t in before_states:
#                    print(CANON.node_data[t]['label'][0] + " --> " + CANON.node_data[CANON.n-1]['label'][0] + " = " + str(bstate_trans) + ".")
#                    CANON.a[t][CANON.n-1] = bstate_trans
#        if i == l - 1:
#            for j in after_states:
#                CANON.a[CANON.n-1][j] = 1 / len(after_states)
#        if i > 0:
#            #CANON.a[CANON.n-2][CANON.n-1] = 1
#    return CANON


# m states, match a histogram length distribution by transitioning to a discrete number of l states with set probability. 
def make_m_state(CANON, m_ind, before_states, after_states, distribution, _name_, node_data):
    # distribution is a list of length probabiliites index i is the percent chance of m transitioning to l state i 
    rnge = len(distribution)
    make_l_states(CANON, rnge, node_data, [], after_states, "_L"+_name_+"1_")
    for i in range(0, rnge-1):
        l_state = CANON.get_node_by_label('l_L'+_name_+"1_"+str(i))
        CANON.a[m_ind][l_state] = distribution[i]
    CANON.a[before_states[0]][m_ind] = 1
    CANON.a[m_ind][after_states[0]] = 0.01

#def pyramid_cascade(node_data, layer_1_states, _name_, E_name):
#    polar_number = 0
#    for layer in range(0, len(layer_1_states)):
#        if layer % 2 == 0: # disconnected, polar layer
#            polar_number += 1
#            # get ppd:
#            ppd = get_tm_polar_res(data_dict, 30, polar_number) # 10 ... 29 are the important ones
#            #print("POLAR PROBABILITY DISTRIBUTION = " + str(ppd))
#            node = {}
#            for k in list(node_data):
#                node[k] = [ l for l in node_data[k]]
#            for i,a in enumerate(node["AA"]): # create a polar-only distribution
#                node["AA"][i] = 0
#                if alpha[i] in ["K", "R", "H", "Q", "E", "D"]:
#                    node["AA"][i] = 1
#            for Litr in range(layer+1,len(layer_1_states)):
#                before_state = [CANON.get_node_by_label("l_L"+_name_+str(layer+1)+"_"+str(Litr))]
#                after_state = []
#                if Litr == len(layer_1_states)-1:
#                    after_state = [CANON.get_node_by_label("l"+E_name+"1")]
#                # get state transition probability
#                ppd_index = Litr + 9
#                bstate_trans = ppd[ppd_index]
#                #print("BSTATE_TRANS = " + str(bstate_trans))
#                make_l_states(1,node, before_state,after_state,"_L"+_name_+str(layer+2)+"_",label_offset = Litr,bstate_trans = bstate_trans)
#        else:
#            after_state = [CANON.get_node_by_label("l"+E_name+"1")]
#            make_l_states(len(layer_1_states)-1-layer,node_data,[],after_state,"_L"+_name_+str(layer+2)+"_",label_offset = layer + 1)
#            for Litr in range(layer+1, len(layer_1_states)):
#                bstate_name = "l_L"+_name_+str(layer+1)+"_"+str(Litr)
#                astate_name = "l_L"+_name_+str(layer+2)+"_"+str(Litr+1)
#                before_state = CANON.get_node_by_label(bstate_name)
#                current_state = CANON.get_node_by_label(astate_name)
#                print(bstate_name + "    " + astate_name)
#                CANON.a[before_state][current_state] = 1


def make_TM_region(CANON, TM_region_type, node_data, conseq_length_distribution ):
    GLOB_1_entering_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[CANON.u][i] > 0]
    GLOB_1_leaving_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[i][CANON.u] > 0]
    adjusted_distribution = [ l - 15 for l in conseq_length_distribution ]
    length_probabilities = [0 for i in range(0, 20)]
    for l in adjusted_distribution:
        if l > 0 and l < 20:
            length_probabilities[l] += 1
    length_probabilities.reverse()
    s = sum(length_probabilities)
    for i in range(0, len(length_probabilities)):
        length_probabilities[i] /= s

    alpha = CANON.emissions["AA"]

    core_chain_length = 5
    pyra_width = 20
    scap_pyra_width = 4
    ecap_pyra_width = 4

    def DL2l(D,L,name):
        return "L_d"+name+str(D)+"_"+str(L)
    # FIRST CAP STATE
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['S_TM'+TM_region_type+"_0"]
    CANON.add_node(node,prior=0)
    SCAP_node = CANON.get_node_by_label('S_TM'+TM_region_type+"_0")
    # add 5 helix core linear
    make_l_states(CANON,core_chain_length,node_data,[],[],'_TM'+TM_region_type+"CORE_")
       # add SCAP_0 incoming trans
    for node in GLOB_1_leaving_states:
        CANON.a[node][SCAP_node] = 1
       # create SCAP pyramid
    scap_pyra_name = '_TM'+TM_region_type+"SCAP_"
    add_pyramid_cascade(CANON,node_data,scap_pyra_width,scap_pyra_name,"l_TM"+TM_region_type+"CORE_1")
       # add SCAP_0 outgoing trans to SCAP pyramid
    for d in range(0,scap_pyra_width):
        s = DL2l(d,0,scap_pyra_name)
        CANON.a[SCAP_node][CANON.get_node_by_label(s)] = 1
    
    # Create Helix Core Cascade
    #add_pyramid_cascade(CANON,node_data,7,'_TM'+TM_region_type+"ECAP_","l_TM"+TM_region_type+"_ECAP_1")
    before_m = CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_2")
    m = CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_3")
    after_m = CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_4")
    pyra_name = '_TM'+TM_region_type+"CORE_"
    add_pyramid_cascade(CANON,node_data,pyra_width,pyra_name,"l_TM"+TM_region_type+"CORE_4")
    for i in range(0,pyra_width):
        L1state = CANON.get_node_by_label(DL2l(i,0,'_TM'+TM_region_type+'CORE_'))
        CANON.a[m][L1state] = length_probabilities[i]

    # LAST CAP STATE + ECAP pyra
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['E_TM'+TM_region_type+"_0"]
    CANON.add_node(node, prior = 0)
    ECAP_node = CANON.get_node_by_label('E_TM'+TM_region_type+"_0")
    # create ECAP pyramid
    ecap_pyra_name = '_TM'+TM_region_type+"ECAP_"
    add_pyramid_cascade(CANON,node_data,ecap_pyra_width,ecap_pyra_name,"E_TM"+TM_region_type+"_0")
       # add CORE_5 outgoing trans to SCAP pyramid
    for d in range(0,ecap_pyra_width):
        s = DL2l(d,0,ecap_pyra_name)
        CANON.a[CANON.get_node_by_label('l_TM'+TM_region_type+"CORE_5")][CANON.get_node_by_label(s)] = 1
       # add ECAP_0 outgoing trans
    for node in GLOB_1_entering_states:
        CANON.a[ECAP_node][node] = 1

    
    # CORE HELIX B-TIES ###################################
    # make a B-tie of all odd layer states in the pyramid
    # make a B-tie of all even layer states in the pyramid
    core_odd_layers = []
    core_even_layers = []
    for l in range(0,pyra_width):
        if l % 2 == 1:
            for d in range(l,pyra_width):
                core_odd_layers.append(DL2l(d,l,pyra_name))
        else:
            for d in range(l, pyra_width):
                core_even_layers.append(DL2l(d,l,pyra_name))
    core_odd_layer_nodes = [ CANON.get_node_by_label(l) for l in core_odd_layers]
    core_even_layer_nodes = [ CANON.get_node_by_label(l) for l in core_even_layers]

    # add core chain to B-tie with odd layers
    for k in range(0, core_chain_length):
        core_odd_layers.append(CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_"+str(k+1)))
    
    CANON.make_tie_state(core_even_layer_nodes, "B")
    CANON.make_tie_state(core_odd_layer_nodes, "B")
    ######################################################
    # SCAP HELIX B-TIES ###########################
    ######################################################
    scap_odd_layers = []
    scap_even_layers = []
    for l in range(0, scap_pyra_width):
        if l % 2 == 1:
            for d in range(l, scap_pyra_width):
                scap_odd_layers.append(DL2l(d,l,scap_pyra_name))
        else:
            for d in range(l, scap_pyra_width):
                scap_even_layers.append(DL2l(d,l,scap_pyra_name))
    scap_odd_layer_nodes = [ CANON.get_node_by_label(l) for l in scap_odd_layers]
    scap_even_layer_nodes = [ CANON.get_node_by_label(l) for l in scap_even_layers]
    CANON.make_tie_state(scap_even_layer_nodes, "B")
    CANON.make_tie_state(scap_odd_layer_nodes, "B")
    ######################################################
    # ECAP HELIX B-TIES ###########################
    ecap_odd_layers = []
    ecap_even_layers = []
    for l in range(0, ecap_pyra_width):
        if l % 2 == 1:
            for d in range(l, ecap_pyra_width):
                ecap_odd_layers.append(DL2l(d,l,ecap_pyra_name))
        else:
            for d in range(l, scap_pyra_width):
                ecap_even_layers.append(DL2l(d,l,ecap_pyra_name))
    ecap_odd_layer_nodes = [ CANON.get_node_by_label(l) for l in ecap_odd_layers]
    ecap_even_layer_nodes = [ CANON.get_node_by_label(l) for l in ecap_even_layers]
    CANON.make_tie_state(ecap_even_layer_nodes, "B")
    CANON.make_tie_state(ecap_odd_layer_nodes, "B")
    ######################################################
    # make every odd layer in core pyramid polar
    for l in range(0, 20):
        if l % 2 == 1:
            for d in range(l, 20):
                node = CANON.get_node_by_label(DL2l(d,l,pyra_name))
                node_aa = [dat for dat in CANON.node_data[node]["AA"]]
                for i,a in enumerate(node_aa):
                    if alpha[i] in ["K", "R", "H", "Q", "E", "D", 'S', 'T']:
                        node_aa[i] = 1
                    else: 
                        node_aa[i] = 0.1
                CANON.node_data[node]["AA"] = [ dat for dat in node_aa ] 
    
    # make every even layer in SCAP and ECAP pyramid amphipathic
    for l in range(0,4):
        if l % 2 == 0:
            for d in range(l,4):
                scap_node = CANON.get_node_by_label(DL2l(d,l,scap_pyra_name))
                ecap_node = CANON.get_node_by_label(DL2l(d,l,ecap_pyra_name))
                node_aa = [ dat for dat in CANON.node_data[scap_node]["AA"]]
                for i, a in enumerate(node_aa):
                    if alpha[i] in ['W', 'Y', 'K', 'R', 'D', 'S', 'T']:
                        node_aa[i] = 1
                    else: 
                        node_aa[i] = 0.1
                CANON.node_data[scap_node]["AA"] = [ dat for dat in node_aa]
                CANON.node_data[ecap_node]["AA"] = [ dat for dat in node_aa]


    # add transitions from each even numbered core pyramid layer(d,l) to each(d+1,l+1) equal to
    #   probability of position (end-7-20+d+1) being the l+1th polar residue in a TMH. ( use datamine )
    for l in range(0,pyra_width-1):
        if l % 2 == 0:
            # get distribution of l+1th polar residue (counted from end TMH res)
            distribution = polar_res(l+1, 27, tm_type = TM_region_type)
            np_distro = []
            for i,j in enumerate(distribution):
                for k in range(0, j):
                    np_distro.append(i)
            k = int(l/2)+1
            sk = ""
            if k == 1:
                sk = "1st"
            if k == 2:
                sk = "2nd"
            if k == 3:
                sk == "3rd"
            if sk == "":
                sk = str(k) + "th"
            view_np(np_distro, 'TMH ' + sk + " POLAR RESIDUE POSITION", './', 
                    nf = max(np_distro)-min(np_distro), QUIET = False, x = "Position from end of TMH.")
            #view_prob_line_plot(distribution, 'TMH_POLAR_RES_' + str(l/2) + "_DISTRIBUTION",
            #        [ str(i) for i in range(0, len(distribution))], "TMH_POLAR_RES_"+str(l/2), './')  
            l1rescount = sum(distribution) # the total number of TMH's with l+1 polar residues (denominator) 
            for d in range(l,pyra_width):
                try:
                    p = distribution[26-d] / l1rescount
                except:
                    p = 0.001
                even_state_s = DL2l(d,l,pyra_name) 
                even_state = CANON.get_node_by_label(even_state_s)
                adjacent_even_state = CANON.get_node_by_label(DL2l(d+1, l, pyra_name))
                odd_state_s = DL2l(d+1,l+1, pyra_name)
                odd_state = CANON.get_node_by_label(odd_state_s)
                CANON.a[even_state][odd_state] = p
                print(even_state_s + " --> " + odd_state_s + ", " + str(p))
                CANON.a[even_state][adjacent_even_state] = 1-p
    raise
                
# TM_B has 1 start cap state like TMH region, 
# start cap state transitions into either a polar or nonpolar state based on % TMB residue 2 to be polar
# a TMB core helix stems from the polar and nonpolar route with alternating 
def make_TM_B_region(CANON, node_data, distribution): 
    #print(distribution)
    #view_prob_line_plot(distribution, "TM_B_periodicitiy", [str(i) for i in range(0,len(distribution))], path='./', prediction_name='TM_B_periodicity')

    # make source state
    GLOB_1_entering_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[CANON.u][i] > 0]
    GLOB_1_leaving_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[i][CANON.u] > 0]
    adjusted_distribution = [ l - 3 for l in distribution ]
    length_probabilities = [0 for i in range(0, 15)]
    for l in adjusted_distribution:
        if l >  0 and l < 15:
            length_probabilities[l] += 1
    s = sum(length_probabilities)
    for i in range(0, len(length_probabilities)):
        length_probabilities[i] /= s
    length_probabilities.reverse()
    
    pyra_width = 15
    alpha = CANON.emissions["AA"]
    
    # get probability of second residue in a TMB region being polar
    data_dict = datamine('./tm_training.drct')
    p_np_res_2 = [0,0]
    for datum in data_dict:
        if datum["TMSEQ"].find("B") != -1:
            for i in range(0, len(datum["TMSEQ"])):
                if datum["TMSEQ"][i] == "B" and datum["TMSEQ"][i-2] != "B":
                    if datum["AASEQ"][i] in ["R", "K", "Q", "E", "D", 'S', 'T']:
                        p_np_res_2[0] += 1
                    else:
                        p_np_res_2[1] += 1

    print("Probability polar begin TMB = " + str(p_np_res_2[0] / sum(p_np_res_2)))
    print("Probability nonpolar begin TMB = " + str(p_np_res_2[1] / sum(p_np_res_2)))
    def DL2l(D,L,name):
        return "L_d"+name+str(D)+"_"+str(L)
    # FIRST CAP STATE
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['S_TMB_0']
    CANON.add_node(node,prior=0)
    SCAP_node = CANON.get_node_by_label('S_TMB_0')
    # Polar M 
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['M_TMB_P']
    CANON.add_node(node,prior=0)
    PM_node = CANON.get_node_by_label('M_TMB_P')
    # Non-Polar M 
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['M_TMB_NP']
    CANON.add_node(node,prior=0)
    NPM_node = CANON.get_node_by_label('M_TMB_NP')
    # add polar or non-polar transitions
    CANON.a[SCAP_node][PM_node] = p_np_res_2[0] / sum(p_np_res_2)
    CANON.a[SCAP_node][NPM_node] = p_np_res_2[1] / sum(p_np_res_2)
    # ADD SCAP incoming trans
    for node in GLOB_1_leaving_states:
        CANON.a[node][SCAP_node] = 1
    # create ECAP state
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['E_TMB_0']
    CANON.add_node(node,prior=0)
    ECAP_node = CANON.get_node_by_label('E_TMB_0')
    # ADD ECAP outgoing trans
    for node in GLOB_1_entering_states:
        CANON.a[ECAP_node][node] = 1

    # make the two pyramid cascades 
    polar_pyra_name = "_TMBPOLARCORE_"
    add_pyramid_cascade(CANON,node_data,pyra_width,polar_pyra_name,"E_TMB_0")
    for i in range(0,pyra_width):
        L1state = CANON.get_node_by_label(DL2l(i,0,polar_pyra_name))
        CANON.a[PM_node][L1state] = length_probabilities[i]
    # nonpolar
    nonpolar_pyra_name = '_TMBNONPOLARCORE_'
    add_pyramid_cascade(CANON,node_data,pyra_width,nonpolar_pyra_name,"E_TMB_0")
    for i in range(0,pyra_width):
        L1state = CANON.get_node_by_label(DL2l(i,0,nonpolar_pyra_name))
        CANON.a[NPM_node][L1state] = length_probabilities[i]

    # create polar residues within pyramid
    for l in range(0, pyra_width):
        if l % 2 == 1:  # make each odd layer in polar pyramid polar
            for d in range(l, pyra_width):
                s = DL2l(d, l, polar_pyra_name)
                node_data_aa = [ i for i in CANON.node_data[CANON.get_node_by_label(s)]["AA"]]
                for i in range(0, len(node_data_aa)):
                    if alpha[i] in ["R", "K", "Q", "E", "D", 'S', 'T']:
                        node_data_aa[i] = 1
                    else: 
                        node_data_aa[i] = 0.1
                CANON.node_data[CANON.get_node_by_label(s)]["AA"] = [ i for i in node_data_aa]
        else:   # make each even layer in nonpolar pyramid polar
            for d in range(l, pyra_width):
                s = DL2l(d, l, nonpolar_pyra_name)
                node_data_aa = [ i for i in CANON.node_data[CANON.get_node_by_label(s)]["AA"]]
                for i in range(0, len(node_data_aa)):
                    if alpha[i] in ["R", "K", "Q", "E", "D", 'S', 'T']:
                        node_data_aa[i] = 1
                    else: 
                        node_data_aa[i] = 0.1
                CANON.node_data[CANON.get_node_by_label(s)]["AA"] = [ i for i in node_data_aa]
         
    # remove all transitions d,l --> d+1,l
    for l in range(0, pyra_width):
        for d in range(l, pyra_width):
            s = CANON.get_node_by_label(DL2l(d, l, polar_pyra_name))
            e = CANON.get_node_by_label(DL2l(d+1, l, polar_pyra_name))
            CANON.a[s][e] = 0
            s = CANON.get_node_by_label(DL2l(d, l, nonpolar_pyra_name))
            e = CANON.get_node_by_label(DL2l(d+1, l, nonpolar_pyra_name))
            CANON.a[s][e] = 0

            
    # B-tie all odd layer polarcore + all even layer nonpolar core
    B_polar = [PM_node]
    for l in range(0, pyra_width):
        for d in range(l, pyra_width):
            if l % 2 == 0:
                B_polar.append(CANON.get_node_by_label(DL2l(d,l,nonpolar_pyra_name)))
            else:
                B_polar.append(CANON.get_node_by_label(DL2l(d,l,polar_pyra_name)))
    CANON.make_tie_state(B_polar, "B")
    # B-tie all even layer polarcore + all odd layer nonpolar core
    B_nonpolar = [NPM_node]
    for l in range(0, pyra_width):
        for d in range(l, pyra_width):
            if l % 2 == 1:
                B_polar.append(CANON.get_node_by_label(DL2l(d,l,nonpolar_pyra_name)))
            else:
                B_polar.append(CANON.get_node_by_label(DL2l(d,l,polar_pyra_name)))
    CANON.make_tie_state(B_nonpolar, "B")


def make_TM_C_region(CANON, node_data, conseq_length_distribution ):
    TM_region_type = 'C'
    GLOB_1_entering_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[CANON.u][i] > 0]
    GLOB_1_leaving_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[i][CANON.u] > 0]
    adjusted_distribution = [ l - 15 for l in conseq_length_distribution ]
    length_probabilities = [0 for i in range(0, 20)]
    for l in adjusted_distribution:
        if l > 0 and l < 20:
            length_probabilities[l] += 1
    length_probabilities.reverse()
    s = sum(length_probabilities)
    for i in range(0, len(length_probabilities)):
        length_probabilities[i] /= s

    alpha = CANON.emissions["AA"]

    core_chain_length = 5
    pyra_width = 20
    scap_pyra_width = 4
    ecap_pyra_width = 4

    def DL2l(D,L,name):
        return "L_d"+name+str(D)+"_"+str(L)
    # FIRST CAP STATE
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['S_TM'+TM_region_type+"_0"]
    CANON.add_node(node,prior=0)
    SCAP_node = CANON.get_node_by_label('S_TM'+TM_region_type+"_0")
    # add 5 helix core linear
    make_l_states(CANON,core_chain_length,node_data,[],[],'_TM'+TM_region_type+"CORE_")
       # add SCAP_0 incoming trans
    for node in GLOB_1_leaving_states:
        CANON.a[node][SCAP_node] = 1
       # create SCAP pyramid
    scap_pyra_name = '_TM'+TM_region_type+"SCAP_"
    add_pyramid_cascade(CANON,node_data,scap_pyra_width,scap_pyra_name,"l_TM"+TM_region_type+"CORE_1")
       # add SCAP_0 outgoing trans to SCAP pyramid
    for d in range(0,scap_pyra_width):
        s = DL2l(d,0,scap_pyra_name)
        CANON.a[SCAP_node][CANON.get_node_by_label(s)] = 1
    
    # Create Helix Core Chain
    before_m = [CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_2")]
    m = CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_3")
    after_m = [CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_4")]
    pyra_name = '_TM'+TM_region_type+"CORE_"
    make_m_state(CANON, m, before_m, after_m, length_probabilities, pyra_name, node_data)

    # LAST CAP STATE + ECAP pyra
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['E_TM'+TM_region_type+"_0"]
    CANON.add_node(node, prior = 0)
    ECAP_node = CANON.get_node_by_label('E_TM'+TM_region_type+"_0")
    # create ECAP pyramid
    ecap_pyra_name = '_TM'+TM_region_type+"ECAP_"
    add_pyramid_cascade(CANON,node_data,ecap_pyra_width,ecap_pyra_name,"E_TM"+TM_region_type+"_0")
       # add CORE_5 outgoing trans to SCAP pyramid
    for d in range(0,ecap_pyra_width):
        s = DL2l(d,0,ecap_pyra_name)
        CANON.a[CANON.get_node_by_label('l_TM'+TM_region_type+"CORE_5")][CANON.get_node_by_label(s)] = 1
       # add ECAP_0 outgoing trans
    for node in GLOB_1_entering_states:
        CANON.a[ECAP_node][node] = 1

    
    # CORE HELIX B-TIES ###################################
    # make a B-tie of all odd layer states in the pyramid
    # make a B-tie of all even layer states in the pyramid
    core_nodes = []
    for d in range(0, len(length_probabilities)):
        core_nodes.append(CANON.get_node_by_label('l_L_TM'+TM_region_type+"CORE_1_"+str(d))) 
    # add core chain to B-tie with odd layers
    for k in range(0, core_chain_length):
        core_nodes.append(CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_"+str(k+1)))
    CANON.make_tie_state(core_nodes, "B")
    ######################################################
    # SCAP HELIX B-TIES ###########################
    ######################################################
    scap_odd_layers = []
    scap_even_layers = []
    for l in range(0, scap_pyra_width):
        if l % 2 == 1:
            for d in range(l, scap_pyra_width):
                scap_odd_layers.append(DL2l(d,l,scap_pyra_name))
        else:
            for d in range(l, scap_pyra_width):
                scap_even_layers.append(DL2l(d,l,scap_pyra_name))
    scap_odd_layer_nodes = [ CANON.get_node_by_label(l) for l in scap_odd_layers]
    scap_even_layer_nodes = [ CANON.get_node_by_label(l) for l in scap_even_layers]
    CANON.make_tie_state(scap_even_layer_nodes, "B")
    CANON.make_tie_state(scap_odd_layer_nodes, "B")
    ######################################################
    # ECAP HELIX B-TIES ###########################
    ecap_odd_layers = []
    ecap_even_layers = []
    for l in range(0, ecap_pyra_width):
        if l % 2 == 1:
            for d in range(l, ecap_pyra_width):
                ecap_odd_layers.append(DL2l(d,l,ecap_pyra_name))
        else:
            for d in range(l, scap_pyra_width):
                ecap_even_layers.append(DL2l(d,l,ecap_pyra_name))
    ecap_odd_layer_nodes = [ CANON.get_node_by_label(l) for l in ecap_odd_layers]
    ecap_even_layer_nodes = [ CANON.get_node_by_label(l) for l in ecap_even_layers]
    CANON.make_tie_state(ecap_even_layer_nodes, "B")
    CANON.make_tie_state(ecap_odd_layer_nodes, "B")
    ######################################################
    
    # make every even layer in SCAP and ECAP pyramid amphipathic
    for l in range(0,4):
        if l % 2 == 0:
            for d in range(l,4):
                scap_node = CANON.get_node_by_label(DL2l(d,l,scap_pyra_name))
                ecap_node = CANON.get_node_by_label(DL2l(d,l,ecap_pyra_name))
                node_aa = [ dat for dat in CANON.node_data[scap_node]["AA"]]
                for i, a in enumerate(node_aa):
                    if alpha[i] in ['W', 'Y', 'K', 'R', 'D', 'S', 'T']:
                        node_aa[i] = 1
                    else: 
                        node_aa[i] = 0.1
                CANON.node_data[scap_node]["AA"] = [ dat for dat in node_aa]
                CANON.node_data[ecap_node]["AA"] = [ dat for dat in node_aa]

def make_TM_I_region(CANON, node_data, conseq_length_distribution ):
    TM_region_type = 'I'
    GLOB_1_entering_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[CANON.u][i] > 0]
    GLOB_1_leaving_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[i][CANON.u] > 0]
    adjusted_distribution = [ l - 15 for l in conseq_length_distribution ]
    length_probabilities = [0 for i in range(0, 25)]
    for l in adjusted_distribution:
        if l > 0 and l < 25:
            length_probabilities[l] += 1
    length_probabilities.reverse()
    s = sum(length_probabilities)
    for i in range(0, len(length_probabilities)):
        length_probabilities[i] /= s

    alpha = CANON.emissions["AA"]

    core_chain_length = 5
    pyra_width = 20
    scap_pyra_width = 4
    ecap_pyra_width = 4

    def DL2l(D,L,name):
        return "L_d"+name+str(D)+"_"+str(L)
    # FIRST CAP STATE
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['S_TM'+TM_region_type+"_0"]
    CANON.add_node(node,prior=0)
    SCAP_node = CANON.get_node_by_label('S_TM'+TM_region_type+"_0")
    # add 5 helix core linear
    make_l_states(CANON,core_chain_length,node_data,[],[],'_TM'+TM_region_type+"CORE_")
       # add SCAP_0 incoming trans
    for node in GLOB_1_leaving_states:
        CANON.a[node][SCAP_node] = 1
       # create SCAP pyramid
    scap_pyra_name = '_TM'+TM_region_type+"SCAP_"
    add_pyramid_cascade(CANON,node_data,scap_pyra_width,scap_pyra_name,"l_TM"+TM_region_type+"CORE_1")
       # add SCAP_0 outgoing trans to SCAP pyramid
    for d in range(0,scap_pyra_width):
        s = DL2l(d,0,scap_pyra_name)
        CANON.a[SCAP_node][CANON.get_node_by_label(s)] = 1
    
    # Create Helix Core Chain
    before_m = [CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_2")]
    m = CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_3")
    after_m = [CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_4")]
    pyra_name = '_TM'+TM_region_type+"CORE_"
    make_m_state(CANON, m, before_m, after_m, length_probabilities, pyra_name, node_data)

    # LAST CAP STATE + ECAP pyra
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['E_TM'+TM_region_type+"_0"]
    CANON.add_node(node, prior = 0)
    ECAP_node = CANON.get_node_by_label('E_TM'+TM_region_type+"_0")
    # create ECAP pyramid
    ecap_pyra_name = '_TM'+TM_region_type+"ECAP_"
    add_pyramid_cascade(CANON,node_data,ecap_pyra_width,ecap_pyra_name,"E_TM"+TM_region_type+"_0")
       # add CORE_5 outgoing trans to SCAP pyramid
    for d in range(0,ecap_pyra_width):
        s = DL2l(d,0,ecap_pyra_name)
        CANON.a[CANON.get_node_by_label('l_TM'+TM_region_type+"CORE_5")][CANON.get_node_by_label(s)] = 1
       # add ECAP_0 outgoing trans
    for node in GLOB_1_entering_states:
        CANON.a[ECAP_node][node] = 1

    
    # CORE HELIX B-TIES ###################################
    # make a B-tie of all odd layer states in the pyramid
    # make a B-tie of all even layer states in the pyramid
    core_nodes = []
    for d in range(0, len(length_probabilities)):
        core_nodes.append(CANON.get_node_by_label('l_L_TM'+TM_region_type+"CORE_1_"+str(d))) 
    # add core chain to B-tie with odd layers
    for k in range(0, core_chain_length):
        core_nodes.append(CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_"+str(k+1)))
    CANON.make_tie_state(core_nodes, "B")
    ######################################################
    # SCAP HELIX B-TIES ###########################
    ######################################################
    scap_odd_layers = []
    scap_even_layers = []
    for l in range(0, scap_pyra_width):
        if l % 2 == 1:
            for d in range(l, scap_pyra_width):
                scap_odd_layers.append(DL2l(d,l,scap_pyra_name))
        else:
            for d in range(l, scap_pyra_width):
                scap_even_layers.append(DL2l(d,l,scap_pyra_name))
    scap_odd_layer_nodes = [ CANON.get_node_by_label(l) for l in scap_odd_layers]
    scap_even_layer_nodes = [ CANON.get_node_by_label(l) for l in scap_even_layers]
    CANON.make_tie_state(scap_even_layer_nodes, "B")
    CANON.make_tie_state(scap_odd_layer_nodes, "B")
    ######################################################
    # ECAP HELIX B-TIES ###########################
    ecap_odd_layers = []
    ecap_even_layers = []
    for l in range(0, ecap_pyra_width):
        if l % 2 == 1:
            for d in range(l, ecap_pyra_width):
                ecap_odd_layers.append(DL2l(d,l,ecap_pyra_name))
        else:
            for d in range(l, scap_pyra_width):
                ecap_even_layers.append(DL2l(d,l,ecap_pyra_name))
    ecap_odd_layer_nodes = [ CANON.get_node_by_label(l) for l in ecap_odd_layers]
    ecap_even_layer_nodes = [ CANON.get_node_by_label(l) for l in ecap_even_layers]
    CANON.make_tie_state(ecap_even_layer_nodes, "B")
    CANON.make_tie_state(ecap_odd_layer_nodes, "B")
    ######################################################
    
    # make every even layer in SCAP and ECAP pyramid amphipathic
    for l in range(0,4):
        if l % 2 == 0:
            for d in range(l,4):
                scap_node = CANON.get_node_by_label(DL2l(d,l,scap_pyra_name))
                ecap_node = CANON.get_node_by_label(DL2l(d,l,ecap_pyra_name))
                node_aa = [ dat for dat in CANON.node_data[scap_node]["AA"]]
                for i, a in enumerate(node_aa):
                    if alpha[i] in ['W', 'Y', 'K', 'R', 'D', 'S', 'T']:
                        node_aa[i] = 1
                    else: 
                        node_aa[i] = 0.1
                CANON.node_data[scap_node]["AA"] = [ dat for dat in node_aa]
                CANON.node_data[ecap_node]["AA"] = [ dat for dat in node_aa]


def make_TM_L_region(CANON, node_data, conseq_length_distribution ):
    TM_region_type = 'L'
    GLOB_1_entering_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[CANON.u][i] > 0]
    GLOB_1_leaving_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[i][CANON.u] > 0]
    adjusted_distribution = [ l - 15 for l in conseq_length_distribution ]
    length_probabilities = [0 for i in range(0, 20)]
    for l in adjusted_distribution:
        if l > 0 and l < 20:
            length_probabilities[l] += 1
    length_probabilities.reverse()
    s = sum(length_probabilities)
    for i in range(0, len(length_probabilities)):
        length_probabilities[i] /= s

    alpha = CANON.emissions["AA"]

    core_chain_length = 5
    pyra_width = 20
    scap_pyra_width = 4
    ecap_pyra_width = 4

    def DL2l(D,L,name):
        return "L_d"+name+str(D)+"_"+str(L)
    # FIRST CAP STATE
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['S_TM'+TM_region_type+"_0"]
    CANON.add_node(node,prior=0)
    SCAP_node = CANON.get_node_by_label('S_TM'+TM_region_type+"_0")
    # add 5 helix core linear
    make_l_states(CANON,core_chain_length,node_data,[],[],'_TM'+TM_region_type+"CORE_")
       # add SCAP_0 incoming trans
    for node in GLOB_1_leaving_states:
        CANON.a[node][SCAP_node] = 1
       # create SCAP pyramid
    scap_pyra_name = '_TM'+TM_region_type+"SCAP_"
    add_pyramid_cascade(CANON,node_data,scap_pyra_width,scap_pyra_name,"l_TM"+TM_region_type+"CORE_1")
       # add SCAP_0 outgoing trans to SCAP pyramid
    for d in range(0,scap_pyra_width):
        s = DL2l(d,0,scap_pyra_name)
        CANON.a[SCAP_node][CANON.get_node_by_label(s)] = 1
    
    # Create Helix Core Chain
    before_m = [CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_2")]
    m = CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_3")
    after_m = [CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_4")]
    pyra_name = '_TM'+TM_region_type+"CORE_"
    make_m_state(CANON, m, before_m, after_m, length_probabilities, pyra_name, node_data)

    # LAST CAP STATE + ECAP pyra
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['E_TM'+TM_region_type+"_0"]
    CANON.add_node(node, prior = 0)
    ECAP_node = CANON.get_node_by_label('E_TM'+TM_region_type+"_0")
    # create ECAP pyramid
    ecap_pyra_name = '_TM'+TM_region_type+"ECAP_"
    add_pyramid_cascade(CANON,node_data,ecap_pyra_width,ecap_pyra_name,"E_TM"+TM_region_type+"_0")
       # add CORE_5 outgoing trans to SCAP pyramid
    for d in range(0,ecap_pyra_width):
        s = DL2l(d,0,ecap_pyra_name)
        CANON.a[CANON.get_node_by_label('l_TM'+TM_region_type+"CORE_5")][CANON.get_node_by_label(s)] = 1
       # add ECAP_0 outgoing trans
    for node in GLOB_1_entering_states:
        CANON.a[ECAP_node][node] = 1

    
    # CORE HELIX B-TIES ###################################
    # make a B-tie of all odd layer states in the pyramid
    # make a B-tie of all even layer states in the pyramid
    core_nodes = []
    for d in range(0, len(length_probabilities)):
        core_nodes.append(CANON.get_node_by_label('l_L_TM'+TM_region_type+"CORE_1_"+str(d))) 
    # add core chain to B-tie with odd layers
    for k in range(0, core_chain_length):
        core_nodes.append(CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_"+str(k+1)))
    CANON.make_tie_state(core_nodes, "B")
    ######################################################
    # SCAP HELIX B-TIES ###########################
    ######################################################
    scap_odd_layers = []
    scap_even_layers = []
    for l in range(0, scap_pyra_width):
        if l % 2 == 1:
            for d in range(l, scap_pyra_width):
                scap_odd_layers.append(DL2l(d,l,scap_pyra_name))
        else:
            for d in range(l, scap_pyra_width):
                scap_even_layers.append(DL2l(d,l,scap_pyra_name))
    scap_odd_layer_nodes = [ CANON.get_node_by_label(l) for l in scap_odd_layers]
    scap_even_layer_nodes = [ CANON.get_node_by_label(l) for l in scap_even_layers]
    CANON.make_tie_state(scap_even_layer_nodes, "B")
    CANON.make_tie_state(scap_odd_layer_nodes, "B")
    ######################################################
    # ECAP HELIX B-TIES ###########################
    ecap_odd_layers = []
    ecap_even_layers = []
    for l in range(0, ecap_pyra_width):
        if l % 2 == 1:
            for d in range(l, ecap_pyra_width):
                ecap_odd_layers.append(DL2l(d,l,ecap_pyra_name))
        else:
            for d in range(l, scap_pyra_width):
                ecap_even_layers.append(DL2l(d,l,ecap_pyra_name))
    ecap_odd_layer_nodes = [ CANON.get_node_by_label(l) for l in ecap_odd_layers]
    ecap_even_layer_nodes = [ CANON.get_node_by_label(l) for l in ecap_even_layers]
    CANON.make_tie_state(ecap_even_layer_nodes, "B")
    CANON.make_tie_state(ecap_odd_layer_nodes, "B")
    ######################################################
    
    # make every even layer in SCAP and ECAP pyramid amphipathic
    for l in range(0,4):
        if l % 2 == 0:
            for d in range(l,4):
                scap_node = CANON.get_node_by_label(DL2l(d,l,scap_pyra_name))
                ecap_node = CANON.get_node_by_label(DL2l(d,l,ecap_pyra_name))
                node_aa = [ dat for dat in CANON.node_data[scap_node]["AA"]]
                for i, a in enumerate(node_aa):
                    if alpha[i] in ['W', 'Y', 'K', 'R', 'D', 'S', 'T']:
                        node_aa[i] = 1
                    else: 
                        node_aa[i] = 0.1
                CANON.node_data[scap_node]["AA"] = [ dat for dat in node_aa]
                CANON.node_data[ecap_node]["AA"] = [ dat for dat in node_aa]

def make_TM_F_region(CANON, node_data, conseq_length_distribution ):
    TM_region_type = 'F'
    GLOB_1_entering_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[CANON.u][i] > 0]
    GLOB_1_leaving_states = [i for i in range(0,CANON.n) if CANON.node_data[i]["label"][0].startswith("GLOBULAR_SIDE_1") and CANON.a[i][CANON.u] > 0]
    adjusted_distribution = [ l - 15 for l in conseq_length_distribution ]
    length_probabilities = [0 for i in range(0, 20)]
    for l in adjusted_distribution:
        if l > 0 and l < 20:
            length_probabilities[l] += 1
    length_probabilities.reverse()
    s = sum(length_probabilities)
    for i in range(0, len(length_probabilities)):
        length_probabilities[i] /= s

    alpha = CANON.emissions["AA"]

    core_chain_length = 5
    pyra_width = 20
    scap_pyra_width = 4
    ecap_pyra_width = 4

    def DL2l(D,L,name):
        return "L_d"+name+str(D)+"_"+str(L)
    # FIRST CAP STATE
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['S_TM'+TM_region_type+"_0"]
    CANON.add_node(node,prior=0)
    SCAP_node = CANON.get_node_by_label('S_TM'+TM_region_type+"_0")
    # add 5 helix core linear
    make_l_states(CANON,core_chain_length,node_data,[],[],'_TM'+TM_region_type+"CORE_")
       # add SCAP_0 incoming trans
    for node in GLOB_1_leaving_states:
        CANON.a[node][SCAP_node] = 1
       # create SCAP pyramid
    
    # Create Helix Core Chain
    before_m = [CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_2")]
    m = CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_3")
    after_m = [CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_4")]
    pyra_name = '_TM'+TM_region_type+"CORE_"
    make_m_state(CANON, m, before_m, after_m, length_probabilities, pyra_name, node_data)

    CANON.a[SCAP_node, CANON.get_node_by_label('l_TM'+TM_region_type+"CORE_1")] = 1
    # LAST CAP STATE + ECAP pyra
    node = {}
    for k in list(node_data):
        node[k] = [ l for l in node_data[k]]
    node['label'] = ['E_TM'+TM_region_type+"_0"]
    CANON.add_node(node, prior = 0)
    ECAP_node = CANON.get_node_by_label('E_TM'+TM_region_type+"_0")

    CANON.a[CANON.get_node_by_label('l_TM'+TM_region_type+"CORE_5")][ECAP_node] = 1
    for node in GLOB_1_entering_states:
        CANON.a[ECAP_node][node] = 1
    
    # CORE HELIX B-TIES ###################################
    # make a B-tie of all odd layer states in the pyramid
    # make a B-tie of all even layer states in the pyramid
    core_nodes = []
    for d in range(0, len(length_probabilities)):
        core_nodes.append(CANON.get_node_by_label('l_L_TM'+TM_region_type+"CORE_1_"+str(d))) 
    # add core chain to B-tie with odd layers
    for k in range(0, core_chain_length):
        core_nodes.append(CANON.get_node_by_label("l_TM"+TM_region_type+"CORE_"+str(k+1)))
    CANON.make_tie_state(core_nodes, "B")
    ######################################################
    
def copy_TM_region(CANON, TM_region_type, model_file):
    model_name = model_file[model_file.rfind('/')+1:]
    model_dir = model_file[:model_file.rfind('/')+1]
    COPY_CANON = hmm(model_name, directory = model_dir)
    COPY_CANON.read()
    old_As = [] # old node indeces\ PARALLEL LISTS
    new_As = [] # new node indeces/
    for i in range(0, COPY_CANON.n):
        try:
            if COPY_CANON.node_data[i]["label"][0].find("TM"+TM_region_type) != -1:
                CANON.add_node(COPY_CANON.node_data[i], prior = 0)
                old_As.append(i)
                new_As.append(CANON.n-1)
        except:
            print(model_file)
            print(COPY_CANON.node_data[i])
            raise
    for i in range(COPY_CANON.u, COPY_CANON.n):
        for j in range(COPY_CANON.u, COPY_CANON.n):
            ti = i
            tj = j
            if i in old_As:
                ti = new_As[old_As.index(i)]
            if j in old_As:
                tj = new_As[old_As.index(j)]
            if ti != i or tj != j: # at least one of i or j is a new node
                CANON.a[ti][tj] = COPY_CANON.a[i][j]

    
    for i in range(0, COPY_CANON.u):
        # get all B-ties within this u state
        B_states = []
        for j in range(0, len(COPY_CANON.Bties)):
            print(COPY_CANON.Bties)
            if COPY_CANON.Bties[j][1] == i:
                B_states.append(new_As[old_As.index(COPY_CANON.Bties[j][0])]+i)
        print(B_states)
        CANON.make_tie_state(B_states, "B")


def create_TM_regions(CANON, TM_region_type, model_file, glob_n):
    print("CREATING TM" + TM_region_type + " REGION.")
    COPY_CANON = hmm(model_file)
    COPY_CANON.read()
    old_As_s1 = [] # old node indeces\ PARALLEL LISTS
    new_As_s1 = [] # new node indeces/
    
    print("creating s_1 nodes")
    for i in range(0, COPY_CANON.n): # append new nodes
        try:
            if COPY_CANON.node_data[i]["label"][0].find("TM"+TM_region_type) != -1:
                # add new node, change name to s_+name+_1
                node = {}
                for j in list(COPY_CANON.node_data[i]):
                    node[j] = [l for l in COPY_CANON.node_data[i][j] ]
                node["label"] = ["s_" + node["label"][0] + "_1"]
                CANON.add_node(node, prior = 0)
                old_As_s1.append(i)
                new_As_s1.append(CANON.n-1)
        except:
            print(model_file)
            print(COPY_CANON.node_data[i])
            raise
    for i in range(COPY_CANON.u, COPY_CANON.n):
        for j in range(COPY_CANON.u, COPY_CANON.n):
            ti = i
            tj = j
            mark = False
            if i in old_As_s1:
                ti = new_As_s1[old_As_s1.index(i)]
                mark = True
            else:
                ti = i - COPY_CANON.u + 1 + CANON.u
            if j in old_As_s1:
                tj = new_As_s1[old_As_s1.index(j)]
                mark = True
            else:
                tj = j - COPY_CANON.u + 1 + CANON.u
            if mark: # at least one of i or j is a new node
                CANON.a[ti][tj] = COPY_CANON.a[i][j]

    print("copying s_1 b_ties")
    for i in range(0, COPY_CANON.u):
        # get all B-ties within this u state
        B_states = []
        for j in range(0, len(COPY_CANON.Bties)):
        #    print(COPY_CANON.Bties)
            if COPY_CANON.Bties[j][1] == i:
                try:
                    B_states.append(new_As_s1[old_As_s1.index(COPY_CANON.Bties[j][0])]+i)
                except:
                    continue
        #print(B_states)
        CANON.make_tie_state(B_states, "B")
    print("COPIED " + str(COPY_CANON.u) + " tie states.")
    print("Creating s_2 nodes:")
    old_As_s2 = []
    new_As_s2 = []
    for i in range(0, COPY_CANON.n): # append new nodes
        try:
            if COPY_CANON.node_data[i]["label"][0].find("TM"+TM_region_type) != -1:
                # add new node change name to s_+name+_2
                node = {}
                for j in list(COPY_CANON.node_data[i]):
                    node[j] = [l for l in COPY_CANON.node_data[i][j] ]
                node["label"] = ["s_" + node["label"][0] + "_2"]
                CANON.add_node(node, prior=0)
                old_As_s2.append(i)        
                new_As_s2.append(CANON.n-1)
        except:
            print(model_file)
            print(COPY_CANON.node_data[i])
            raise
    for i in range(COPY_CANON.u, COPY_CANON.n):
        for j in range(COPY_CANON.u, COPY_CANON.n):
            ti = i
            tj = j
            mark = False
            if i in old_As_s2:
                ti = new_As_s2[old_As_s2.index(i)]
                mark = True
            else:
                ti = i - COPY_CANON.u + 1 + CANON.u
            if j in old_As_s2:
                tj = new_As_s2[old_As_s2.index(j)]
                mark = True
            else:
                tj = j - COPY_CANON.u + 1 + CANON.u
            if mark: # at least one of i or j is new node
                CANON.a[ti][tj] = COPY_CANON.a[i][j]
    print("Copying s_2 B_ties")
    for i in range(0, COPY_CANON.u):
        # get all B-ties within this u states
        B_states = []
        for j in range(0, len(COPY_CANON.Bties)):
            if COPY_CANON.Bties[j][1] == i:
                try:
                    B_states.append(new_As_s2[old_As_s2.index(COPY_CANON.Bties[j][0])]+1)
                except:
                    continue
        #print(B_states)
        CANON.make_tie_state(B_states, "B")
    print("COPIED " + str(COPY_CANON.u) + " TIE STATES.")

    if TM_region_type not in ["L", "F"]:

        print("make transitions s_E_TMX_0_1 --> glob_2")
        se01 = CANON.get_node_by_label("s_E_TM"+TM_region_type+"_0_1")
        for j in range(0, CANON.n):
            if CANON.node_data[j]['label'][0] == 's_E_TM'+TM_region_type+'_0_1':
                print(j)
        print("SE01 = " + str(se01) + " label = " + str(CANON.node_data[se01]['label'][0]))
        for j in range(0, CANON.n):
            if CANON.a[se01][j] > 0 and CANON.node_data[j]["label"][0] == "GLOBULAR_SIDE_1":
                t = CANON.a[se01][j]
                CANON.a[se01][j+glob_n] = t
                print(str(CANON.node_data[se01]['label'][0]) + " --> " + str(CANON.node_data[j+glob_n]['label'][0]) + " _" + str(j+glob_n))
                CANON.a[se01][j] = 0
            elif CANON.a[se01][j] > 0:
                print(str(CANON.node_data[se01]['label'][0]) + " --> " + str(CANON.node_data[j]['label'][0]) + " " +str(j))


        print("make transitions glob_2 --> s_S_TMX_0_2")
        ss02 = CANON.get_node_by_label('s_S_TM'+TM_region_type+"_0_2")
        print("SS02 = " + str(ss02) + " label = " + str(CANON.node_data[ss02]['label'][0]))
        for i in range(0, CANON.n):
            if CANON.a[i][ss02] > 0 and CANON.node_data[i]["label"][0] == "GLOBULAR_SIDE_1":
                t = CANON.a[i][ss02]
                CANON.a[i+glob_n][ss02] = t
                print(str(CANON.node_data[i+glob_n]['label'][0]) + " _"+str(i+glob_n)+ " --> " + str(CANON.node_data[ss02]['label'][0]))
                CANON.a[i][ss02] = 0 
            elif CANON.a[se01][j] > 0:
                print(str(CANON.node_data[i]['label'][0]) + " "+str(i)+ " --> " + str(CANON.node_data[ss02]['label'][0]))
    else:
        print("make side 2 tranisitons from side 2 to side 2")
        ss02 = CANON.get_node_by_label('s_S_TM'+TM_region_type+"_0_2")
        se02 = CANON.get_node_by_label('s_E_TM'+TM_region_type+"_0_2")
        for j in range(0, CANON.n):
            if CANON.a[se02][j] > 0 and CANON.node_data[j]['label'][0] == "GLOBULAR_SIDE_1":
                t = CANON.a[se02][j]
                CANON.a[se02][j+glob_n] = t
                CANON.a[se02][j] = 0
            if CANON.a[j][ss02] > 0 and CANON.node_data[j]['label'][0] == 'GLOBULAR_SIDE_1':
                t = CANON.a[j][ss02]
                CANON.a[j+glob_n][ss02] = t
                CANON.a[j][ss02] = 0

def copy_labels(mf1, mf2):
    # copy node_data[label] from mf1 to mf2 and write mf2
    mn1 = mf1[mf1.rfind('/')+1:]
    md1 = mf1[:mf1.rfind('/')+1]
    HMM1 = hmm(mn1, directory = md1)
    HMM1.read()
    mn2 = mf2[mf2.rfind('/')+1:]
    md2 = mf2[:mf2.rfind('/')+1]
    HMM2 = hmm(mn2, directory = md2)
    HMM2.read()

    for i in range(0, HMM1.n):
        label = HMM1.node_data[i]["label"]
        HMM2.node_data[i]["label"] = label
    HMM2.write(mf2)

def copy_ties(mf1, mf2):
    HMM1 = hmm(mf1)
    HMM1.read()
    HMM2 = hmm(mf2)
    HMM2.read()

    # Check HMM2's ties
    HMM2.nta = len(HMM2.Aties)
    HMM2.ntb = len(HMM2.Bties)
    HMM2.ntab = len(HMM2.ABties)
    # Check HMM2's tie states
    if HMM2.u > 0:
        found_u_states = [] # look in Bties for u state indeces
        for i in range(0, HMM2.ntb):
            if HMM2.Bties[i][1] not in found_u_states:
                found_u_states.append(HMM2.Bties[i][1])
        remove_u_states = []
        for i in range(0, HMM2.u):
            if i not in found_u_states:
                remove_u_states.append(i)
        remove_u_states.sort()
        remove_u_states.reverse()
        for i in remove_u_states:
            HMM2.remove_node(i)
            print("removed node " + str(i))
    
    # copy over new ties
    tie_groups = []
    for i in range(0, HMM1.u):
        tie_groups.append([])
    for tie in HMM1.Bties:
        tie_groups[tie[1]].append(tie[0]-HMM1.u)
    for i,group in enumerate(tie_groups):
        group = [ g + i for g in group ]
        HMM2.make_tie_state(group, "B")

    HMM2.write(HMM2.file)




def clear_ties(mf):
    HMM = hmm(mf)
    HMM.read()
    for i in range(HMM.u,0,-1):
        HMM.remove_node(i)
    HMM.u = 0
    HMM.ntab = 0
    HMM.nta = 0
    HMM.ntb = 0
    HMM.Bties = []
    HMM.Aties = []
    HMM.ABties = []
    HMM.write(HMM.file)
            
