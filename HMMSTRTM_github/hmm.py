# hmm.py
# read in a hmm-object. do stuff with it.

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

def rgb(minimum, maximum, value):
    minimum, maximum = float(minimum), float(maximum)
    ratio = 2 * (value-minimum) / (maximum - minimum)
    b = int(max(0, 255*(1 - ratio)))
    r = int(max(0, 255*(ratio - 1)))
    g = 255 - b - r
    return r, g, b

def to_percentile(heats, step = 10):
    heats.sort()
    while heats[0] == 0.00:
        heats = heats[1:]
    percentile = [heats[0]]
    step_increment = 100 / step
    last_step = 0
    for i in range(0, len(heats)):
        if i / len(heats) > (last_step+1.0) / step_increment:
            percentile.append(heats[i])
            last_step += 1
    return percentile
def normalize(dat):
    s = sum(dat)
    new_dat = [f/s for f in dat]
    return new_dat

class hmm:
    def __init__(self, f, directory = None):
        self.directory = directory
        if f.rfind('/') != -1:
            self.name = f[f.rfind('/')+1:len(f)-4]
            if self.directory is None:
                self.directory = f[:f.rfind('/')+1]
        else:
            self.name = f[:f.index(".")]
        if self.directory is None:
            self.directory = ""
        self.file = self.directory + self.name + '.hmm'
        self.gamma = None
        self.rama_reweight = [ 1.0 for i in range(1, 7) ]
        self.ss_reweight = [ 1.0 for i in range(1, 4)]
        self.tm_reweight = [ 1.0 for i in range(1, 9)]
        self.Aties = []
        self.ABties = []
        self.Bties = []
        self.nta = 0
        self.ntab = 0
        self.ntb = 0
        self.AA_EPS = 0
        self.SS_EPS = 0
        self.TM_EPS = 0
        self.RAMA_EPS = 0
        

    #def file_reader(self):
    #    for line in open(self.f, "r"):
    #        yield line
########################################################################################################################################
####################### HMMSTR ... HMM I/O #############################################################################################
########################################################################################################################################

    def read_file(self):
        return self.read()
    def read(self):
        try:
            return self.read_file_format_2()
        except:
            return self.read_file_format_0() # try old format reader

    
    def read_file_format_0(self): # The Bystroff-lab grown HMM format.
        self.bg_freqs = []
        self.n = -1
        self.k = -1
        i = 0 
        f = open(self.file, "r")
        self.sumb = 0
        self.sumd = 0
        self.comments = None
        self.u = 0
        self.node_data = []
        lines = (line for line in open(self.file, "r"))
        while True:
            comps = next(lines).split()
            if comps[0] == "default_bg_freqs":
                self.AA_ground = [ float(x) for x in comps[1:] ]
            elif comps[0] == "num_nodes":
                self.n = int(comps[1]) + 1 
                self.priors = np.empty(self.n)   # (n)
                self.dih = np.empty((self.n,3))  # (n,3)
       #        self.ss =  np.empty(self.n)      # (n)
                self.a = np.empty((self.n,self.n))    # transitions (n,n)
                self.b = np.empty((self.n,20))   # (n,m)
                self.r = np.empty((self.n,11))   # (n,mrama)
                self.d = np.empty((self.n,6))    # (n, mdssp)
       #        self.logb = np.empty((self.n,20))# (n,m)
            elif comps[0] == "k_nodes":
                self.k = int(comps[1])
            if self.n == -1:
                continue
            elif comps[0] == "node" or comps[0] == "unk_node":
                uk = comps[0] == "unk_node"
                if uk:
                    self.priors[0] = float(comps[2])
                    self.dih[0] = [-75, -15, 180]
                else:
                    self.priors[i] = float(comps[4])
                    self.dih[i] = [float(s) for s in comps[5:]]
                # AA distribution
                comps = next(lines).split()
                self.b[i][:10] = [float(s) for s in comps]
                comps = next(lines).split()
                self.b[i][10:] = [float(s) for s in comps]
                self.sumb += sum(self.b[i][:])
                # SS
                comps = next(lines).split()
                if uk:
                    comps = comps[2:]
                    self.d[i] = [float(s) for s in comps ]
                else:
                    comps = comps[4:]
                    self.d[i] = [float(s) for s in comps ]
                # rama
                comps = next(lines).split()
                if uk:
                    comps = comps[2:]
                    self.r[i] = [ float(s) for s in comps ]
                else:
                    comps = comps[4:]
                    self.r[i] = [float(s) for s in comps]
                self.node_data.append({'AA' : [ f for f in self.b[i]] , 'SS' : [f for f in self.d[i]], 'TM': [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                    'RAMA': [f for f in self.r[i]], 'DIH': [ f for f in self.dih[i]]})
                #print(self.b[i])
                #print(self.r[i])
                #print(self.d[i])
                i += 1 # node counter
            else:
                if comps[0] == "transit_freqs":
                    for x in range(0, self.n):
                        comps = next(lines).split()
                        self.a[x] = [float(s) for s in comps]
                    # EXPERIMENTAL, TRANSPOSING THE MATRIX
                    #new_a = []
                    #for j in range(0, self.n):
                    #    new_a.append([0 for f in self.a[j]])
                    #    for k in range(0, self.n):
                    #        new_a[j][k] = self.a[k][j]
                    #self.a = new_a
                    ######################################
                    break
            if i > self.n:
                break
        if self.k == -1:
            print("Number of 0-nodes unspecified.")
        if self.n == -1:
            print("Number of nodes unspecified, terminating.", file = sys.stderr)
            raise Exception
        print(len(self.a))
        print(len(self.a[0]))
    
    def read_file_format_2(self): # 2nd iteration

        f = open(self.file, "r")
        self.bg_freqs = []
        self.n = 0
        self.k = 0
        self.u = 0
        i = 0
        self.emissions = {}
        self.node_data = []
        self.priors = []
        self.a = None # use numpy to unsparsify this.
        lines = (line for line in f)
        self.comments = ""
        while True:
            try:
                line = next(lines)
            except StopIteration:
                break
            comps = line.split()
            if len(comps) == 0:
                continue
            if comps[0] == "COMMENT":
                self.comments += line
            elif comps[0] == "AA_EPS":
                self.AA_EPS = float(comps[1])
            elif comps[0] == "SS_EPS":
                self.SS_EPS = float(comps[1])
            elif comps[0] == "TM_EPS":
                self.TM_EPS = float(comps[1])
            elif comps[0] == "RAMA_EPS":
                self.RAMA_EPS = float(comps[1])
            elif comps[0] == "n_a_ties":
                self.nta = int(comps[1])
            elif comps[0] == "n_b_ties":
                self.ntb = int(comps[1])
            elif comps[0] == "n_ab_ties" or comps[0] == 'n_ties':
                self.ntab = int(comps[1])
            elif comps[0] == "n_nodes":
                self.n = int(comps[1])
            elif comps[0] == "k_nodes":
                self.k = int(comps[1])
            elif comps[0] == "u_nodes":
                self.u = int(comps[1])
            elif comps[0] == "emissions":
                for v in comps[1:]:
                    self.emissions[v] = []
                for emission in list(self.emissions):
                    comps = next(lines).split()
                    self.emissions[emission] = comps[0]
            elif comps[0] == "AA_ground" or comps[0] == 'default_bg_freqs':
                self.AA_ground = [float(v) for v in comps[1:]]
            elif comps[0] == "SS_ground":
                self.SS_ground = [float(v) for v in comps[1:]]
            elif comps[0] == "RAMA_ground":
                self.RAMA_ground = [float(v) for v in comps[1:]]
            elif comps[0] == "TM_ground":
                self.TM_ground = [float(v) for v in comps[1:]]
            elif comps[0] == "RAMA_reweight":
                self.rama_reweight = [float(v) for v in comps[1:]]
            elif comps[0] == "SS_reweight":
                self.ss_reweight = [float(v) for v in comps[1:]]
            elif comps[0] == "TM_reweight":
                self.tm_reweight = [float(v) for v in comps[1:]]
            elif comps[0] == "node":
                this_node_data = {}
                this_node_number = int(comps[1])
                self.priors.append(float(comps[2]))
                comps = next(lines).split()
                while len(comps) > 1 and int(comps[1]) == this_node_number:
                    if comps[2] in list(self.emissions):
                        this_node_data[comps[2]] = [float(f) for f in comps[4:]]
                    else:
                        this_node_data[comps[2]] = [f for f in comps[3:]]
                    comps = next(lines).split()
                self.node_data.append(this_node_data)
            elif comps[0] == "START_OUTTRANS":
                self.a = np.zeros((self.n, self.n))
                comps = next(lines).split()
                while comps[0] != "END_OUTTRANS":
                    self.a[int(comps[0])][int(comps[1])] = float(comps[2])
                    try:
                        comps = next(lines).split()
                    except StopIteration:
                        break
            elif comps[0] == "START_TIES":
                self.Aties = []
                self.Bties = []
                self.ABties = []
                comps = next(lines).split()
                tie_type = [ 0,0,0]
                while comps[0] != "END_TIES":
                    if sum(tie_type) == 1:
                        try:
                            if tie_type[0] == 1:
                                self.Aties.append((int(comps[0]),int(comps[1])))
                            elif tie_type[1] == 1:
                                self.Bties.append((int(comps[0]),int(comps[1])))
                            elif tie_type[2] == 1:
                                self.ABties.append((int(comps[0]),int(comps[1])))
                        except:
                            if comps[0] == "A_TIES":
                                tie_type = [1,0,0]
                            elif comps[0] == "B_TIES":
                                tie_type = [0,1,0]
                            elif comps[0] == "AB_TIES":
                                tie_type = [ 0,0,1]
                    if comps[0] == "A_TIES":
                        tie_type = [1,0,0]
                    elif comps[0] == "B_TIES":
                        tie_type = [0,1,0]
                    elif comps[0] == "AB_TIES":
                        tie_type = [ 0,0,1]
                    try:
                        comps = next(lines).split()
                    except StopIteration:
                        break
            elif comps[0] == "EOF":
                print("Successfully Parsed " + self.file +" in hmm.py .")
                
                # print all nodes w/ no intrans or outrans
                x =[]
                for i in range(0, self.n):
                    if sum(self.a[:][i]) == 0 or sum(self.a[i][:]) == 0:
                        x.append(i) #print(i)
                print("nodes with no intrans or outrans:")
                print(x)
                x =[]

                #print(self.node_data[0])
                #print(self.priors)
                break
            else:
                print(comps)
                print("unhandled keyword experienced.")
                raise
    def write(self, filename):
        self.write_file_format_2(filename)

    def write_file_format_2(self, filename): # write itr2
        # TODO: Update transition matrix into sparse representation
        #       Switch to `delimiting keywords` "node 1" rather than "node 0"
        #       keyword template at the top as a commented section. Comment can be a keyword.
        
        out = open(filename, "w")
        
        if self.comments is None: 
            for x in range(0, 10):
                out.write("COMMENT\n")
        else:
            out.write(self.comments)

        out.write("AA_EPS ")
        out.write("{:.15f}".format(self.AA_EPS) + "\n")
        out.write("SS_EPS ")
        out.write("{:.15f}".format(self.SS_EPS) + "\n")
        out.write("TM_EPS ")
        out.write("{:.15f}".format(self.TM_EPS) + "\n")
        out.write("RAMA_EPS ")
        out.write("{:.15f}".format(self.RAMA_EPS) + "\n")
        out.write("n_a_ties " + str(self.nta) + "\n")
        out.write("n_b_ties " + str(self.ntb) + "\n")
        out.write("n_ab_ties " + str(self.ntab) + "\n")
        out.write("n_nodes " + str(self.n) + "\n")
        out.write("k_nodes " + str(self.k) + "\n")
        out.write("u_nodes " + str(self.u) + "\n")
        out.write("emissions ")
        for f in list(self.emissions):
            out.write(f + " ")
        out.write("\n")
        for f in list(self.emissions):
            out.write(self.emissions[f] + "\n")
        out.write("AA_ground ")
        for f in self.AA_ground:
            out.write("{:.5f}".format(f) + " ")
        out.write("\n")
        out.write("SS_ground ")
        for f in self.SS_ground:
            out.write("{:.5f}".format(f) + " ")
        out.write("\n")
        out.write("TM_ground ")
        for f in self.TM_ground:
            out.write("{:.5f}".format(f) + " ")
        out.write("\n")
        out.write("RAMA_ground ")
        for f in self.RAMA_ground:
            out.write("{:.5f}".format(f) + " ")
        out.write("\n")
        out.write("RAMA_reweight ")
        for f in self.rama_reweight:
            out.write("{:.5f}".format(f) + " ")
        out.write("\n")
        out.write("SS_reweight ")
        for f in self.ss_reweight:
            out.write("{:.5f}".format(f) + " ")
        out.write("\n")
        out.write("TM_reweight ")
        for f in self.tm_reweight:
            out.write("{:.5f}".format(f) + " ")
        out.write("\n")
        for i in range(0, self.n):
            out.write("node " + str(i) + " " + "{:.5f}".format(self.priors[i]) + "\n")
            #for emission in list(self.emissions):
            for datum in list(self.node_data[i]):# [emission]:
                out.write("node " + str(i) + " " + datum + " " )
                if datum in self.emissions:
                    this_max_emitter = self.emissions[datum][self.node_data[i][datum].index(max(self.node_data[i][datum]))]
                    out.write(str(this_max_emitter) + " ")
        #        print(self.node_data[i][datum])
        #        print(datum)
                for f in self.node_data[i][datum]:
                    if isinstance(f, float):
                        out.write("{:.5f}".format(f) + " ")
                    else:
                        out.write(str(f) + " ")
                out.write("\n")
            out.write("endnode\n")
        out.write("START_OUTTRANS\n")
        self.normalize_a()
        for i in range(0, self.n):
            for j in range(0, self.n):
                if self.a[i][j] > 0.000:
                    out.write(str(i) + " " + str(j) + " " + "{:.5f}".format(self.a[i][j]) + "\n")
        out.write("END_OUTTRANS\n")
        out.write("START_TIES\n")
        out.write("A_TIES\n")
        for i in self.Aties:
            out.write(str(i[0]) + " " + str(i[1]) + "\n")
        out.write("B_TIES\n")
        for i in self.Bties:
            out.write(str(i[0]) + " " + str(i[1]) + "\n")
        out.write("AB_TIES\n")
        for i in self.ABties:
            out.write(str(i[0]) + " " + str(i[1]) + "\n")
        out.write("END_TIES\n")
        out.write("EOF\n")
        out.close()


        def normalize_a(self):
            for i in range(0, self.n):
                s = 0.0
                for j in range(0, self.n):
                    s += self.a[i][j]
                for j in range(0, self.n):
                    self.a[i][j] /= s
            return

########################################################################################################################################
####################### HMMSTR ... Rabiner Magic #######################################################################################
########################################################################################################################################


    # use A to get random state sequence of defined length
    def stochastic_search(self, length):
        import random
        state_sequence = []
        last_state = 0
        for i in range(0, length):
            rand_trans = random.random()
            this_trans = 0
            for nitr in range(0, self.n):
                this_trans += self.a[last_state][nitr]
                if this_trans > rand_trans:
                    state_sequence.append(nitr)
                    last_state = nitr
                    break
        print(state_sequence)
        self.gamma = np.zeros((self.n, length))
        for i in range(0, len(state_sequence)):
            self.gamma[state_sequence[i]][i] = 1
    
    def train(self, drct_file, epochs): # OKAY SO TURNS OUT WE CAN'T TRAIN B/C READING FORTRAN .drct IN PYTHON IS HARD.
        x = None
        return x
    
    def update_HMMSTR(self): # end of training epoch. Make updates to A,B
        x = None
        return x

    def get_alpha(self):
        x = None
        return x

    def get_beta(self):
        x = None
        return x

    def get_trans(self):
        x = None
        return x

    def get_gamma(self):
        x = None
        return x

    def get_A(self):
        x = None
        return x
    
    def get_B(self):
        x = None
        return x

    def add_A_ties(self, states):
        for state_pair in states:
            self.set_A_tie(state_pair[0], state_pair[1])
    
    def add_AB_ties(self, states):
        for state_pair in states:
            self.set_AB_tie(state_pair[0], state_pair[1])

    def add_B_ties(self, states):
        self.u += 1
        node_data = { "AA": self.AA_ground, "RAMA": self.RAMA_ground, 
                      "TM": self.TM_ground, "SS": self.SS_ground, 'label': ["TIE_STATE"+str(self.u)]}
        self.insert_node(node_data, 0)
        for i,tie in enumerate(self.Bties):
            new_tie = (tie[0], tie[1]+1)
            self.Bties[i] = new_tie
        for state in states:
            self.set_B_tie(state+1, 0)

    def make_tie_state(self, states, tie_type):
        if tie_type.lower() == "a":
            self.add_A_ties(states)
        elif tie_type.lower() == "b":
            self.add_B_ties(states)
        elif tie_type.lower() == "ab":
            self.add_AB_ties(states)
        else:
            print("BAD tie type, accepted tie types are 'A', 'B', or 'AB'.")

    def set_B_tie(self, i, j):
        self.Bties.append((i,j))
        self.ntb += 1
        return
    def set_A_tie(self, i, j):
        self.Aties.append((i,j))
        self.nta += 1
        return
    def set_AB_tie(self, i, j):
        self.ABties.append((i,j))
        #print(str(i) + " " + str(j))
        self.ntab += 1
        return

########################################################################################################################################
####################### HMMSTR ... HMM_GAMMA_FUCNTIONS #################################################################################
########################################################################################################################################

    def merge_state(self, state_i, state_j, weight_i, weight_j):

        profs = [ "AA", "SS", "RAMA", "TM" ]
        new_profs = {} # get new, weighted profile for single state
        for p in profs:
            i_prof = self.node_data[state_i][p]
            j_prof = self.node_data[state_j][p]
            new_prof = [ 0 for l in i_prof]
            for i in range(0, len(new_prof)):
                new_prof[i] = i_prof[i]*weight_i + j_prof[i]*weight_j
            new_profs[p] = [ l for l in new_prof]

        # get new (weighted) incoming and outoging transitions
        new_in_trans = [0 for l in range(0, self.n)]
        new_out_trans = [0 for l in range(0, self.n)]
        for i in range(0, self.n):
            if self.a[i][state_i] > 0 and self.a[i][state_j] > 0:
                new_in_trans[i] = self.a[i][state_i] * weight_i + self.a[i][state_j] * weight_j
            if self.a[state_i][i] > 0 and self.a[state_j][i] > 0:
                new_out_trans[i] = self.a[state_i][i] * weight_i + self.a[state_j][i] * weight_j

        # get new prior.
        new_prior = self.priors[state_i] * weight_i + self.priors[state_j] * weight_j

        if state_i < state : # remove old nodes
            self.remove(state_j) # remove state_j first, 
            self.remove(state_i)
            # then add the new state at state_i's index
            self.insert_node(new_profs, state_i)
        else:
            self.remove(state_i) # else remove state_i first.
            self.remove(state_j)
            self.insert_node(new_profs, state_j)

        # fix transitions ( there will be a -1 offset after i or j is > min(state_i,state_j))
        for i in range(0, self.n):
            nitr = i
            if i > min(state_i, state_j):
                nitr = i - 1
            if new_in_trans[nitr] > 0:
                self.a[nitr][min(state_i, state_j)] = new_in_trans[nitr]
            if new_out_trans[nitr] > 0:
                self.a[min(state_i, state_j)][nitr] = new_in_trans[nitr]
        

    def get_gamma(self, gamma_dir, gamma_file):
        self.gamma = read_gamma(gamma_file, gamma_dir)

    # 6/22/21
    def GAMMA_PREDICT(self, gamma):
        if self.gamma is None:
            self.get_gamma('./HMMSTR_50_GAMMA/', gamma)
        gamma = self.gamma

        # get the paradigm, gamma-weighted profiles for AA, Q3, RAMA_BIN
        t = len(gamma[0])
        AA_paradigm = np.zeros((t, 20))
        SS_paradigm = np.zeros((t, 6))
        TM_paradigm = np.zeros((t, 10))
        RAMA_paradigm = np.zeros((t, 11))
        for titr in range(0, t):
            this_AA_prof   = [0 for t in range(0,20)]
            this_SS_prof   = [0 for t in range(0,6)]
            this_TM_prof   = [0 for t in range(0,10)]
            this_RAMA_prof = [0 for t in range(0,11)]
            for nitr in range(0, self.n):
                for b in range(0,20):
                    this_AA_prof[b] += gamma[titr][nitr] * self.node_data[nitr]["AA"][b]
                    if b < 11:
                        this_RAMA_prof[b] += gamma[titr][nitr] * self.node_data[nitr]["RAMA"][b]
                        if b < 10:
                            this_TM_prof[b] += gamma[titr][nitr] * self.node_data[nitr]["TM"][b]
                            if b < 6:
                                this_SS_prof[b] += gamma[titr][nitr] * self.node_data[nitr]["SS"][b]
            for b in range(0,20):
                AA_paradigm[titr][b] = this_AA_prof[b]
                if b < 11:
                    RAMA_paradigm[titr][b] = this_RAMA_prof[b]
                    if b < 10:
                        TM_paradigm[titr][b] = this_TM_prof[b]
                        if b < 6:
                            SS_paradigm[titr][b] = this_SS_prof[b]
        # print argmax entry in the seq at each time point...
        AA_seq = ""
        SS_seq = ""
        TM_seq = ""
        RAMA_seq = ""
        AA_argmax = np.argmax(AA_paradigm, axis=1)
        SS_argmax = np.argmax(SS_paradigm, axis=1)
        TM_argmax = np.argmax(TM_paradigm, axis=1)
        RAMA_argmax = np.argmax(RAMA_paradigm, axis=1)
        for titr in range(0,t):
            AA_seq += self.emissions["AA"][AA_argmax[titr]]
            SS_seq += self.emissions["SS"][SS_argmax[titr]]
            TM_seq += self.emissions["TM"][TM_argmax[titr]]
            RAMA_seq += self.emissions["RAMA"][RAMA_argmax[titr]]
        
        import logomaker
        import pandas as pd
        if False: # Log transform.. (broken, I think)
            for i in range(0,len(AA_paradigm)):
                for j in range(0,len(AA_paradigm[i])):
                    try:
                        AA_paradigm[i][j] = math.log(AA_paradigm[i][j])
                    except:
                        continue
        if False: # argmax..
            for i in range(0, len(AA_paradigm)):
                this_paradigm = np.zeros(20)
                this_paradigm[AA_argmax[i]] = AA_paradigm[i][AA_argmax[i]]
                AA_paradigm[i] = this_paradigm
        df = pd.DataFrame(AA_paradigm, columns = [f for f in self.emissions["AA"]])
        print("AA PARADIGM MATRIX : ")
        print(df)

        fig,ax = plt.subplots(2,1,figsize=(len(df),3))
        for i in [0,1]:
            if i == 0:
                this_data = df[0:int(len(df)/2)][:]
            else:
                this_data = df[int(len(df)/2):len(df)][:]
            ww_logo = logomaker.Logo(this_data, width = 0.8, ax=ax[i], vpad =0.05, color_scheme = 'dmslogo_charge', font_name = 'Luxi Mono',stack_order='big_on_top'
                    # ,vsep=0.05
                    )
            ww_logo.style_spines(visible=True)
            ww_logo.style_spines(spines=['left','right'],visible=False)
            ww_logo.style_xticks(anchor = 0, spacing = 25, rotation = 45)
            #ww_logo.ax.set_xlim([-1,len(df)*2])
        ww_logo.fig.savefig('./test.png')
        ww_logo.fig.show()
        
        print("AA : " + AA_seq)
        print("SS : " + SS_seq)
        print("TM : " + TM_seq)
        print("RAMA : " + RAMA_seq)
                            

    def gather_accuracy(self, af):
        #print(self.file)
        self.accuracy_file = af
        f = open(self.accuracy_file, "r")
        self.accuracies = [ float(line.split()[0]) for line in f]

    def visualize_accuracies(self):
        plt.style.use('seaborn-white')
        bins = 100 #CONSTRAINT_STYLES[name]['bins']
        rnge = (min(self.accuracies),max(self.accuracies))# min(self.accuracies), max(self.accuracies)) #CONSTRAINT_STYLES[name]['range']
        density = True # CONSTRAINT_STYLES[name]['density']

        plt.hist(self.accuracies, bins, density = density, alpha = 0.75, range = rnge)
        plt.xlabel("MODEL ACCURACIES (P(O| model))")
        if density:
            plt.ylabel('% occurence')
        else:
            plt.ylabel('# of times occured')
        plt.title('Histogram of model accuracies')
        #print(len(data))
        plt.grid(True)
        plt.savefig(self.accuracy_file[:len(self.accuracy_file)-4] + '.png')
        plt.show()

########################################################################################################################################
####################### HMMSTR ... HMM_VISUALIZER ######################################################################################
########################################################################################################################################

    def HMM_CONSISTENCE(self, A_min = 0.7):
        x = 0
    
    def display_top_changes(self, new_nodes, new_trans, f_name = None):
        from hmm_viz import display_top_changes
        display_top_changes(self, new_nodes, new_trans, f_name = f_name)

    # visualize paradigm seq through looping over time interval in gamma...
    def HMM_GAMMA_IMAGE(self, gamma_dir, gamma_file):
        from hmm_viz import HMM_GAMMA_IMAGE
        HMM_GAMMA_IMAGE(self, gamma_dir, gamma_file)
        
    def HMM_VIZ(self): # get graphviz drawings of current HMM and save to directory. 
        from hmm_viz import HMM_VIZ
        HMM_VIZ(self)

    def show(self):
        from hmm_viz import show
        show(self)

########################################################################################################################################
####################### HMMSTR ... ADD/REMOVE NODES ####################################################################################
########################################################################################################################################


    def init_HMMSTR(self): 
        self.bg_freqs = [0.08279, 0.01937, 0.05855, 0.05992, 0.04014, 0.08089, 0.02275, 0.05552, 0.05959, 0.0802, 0.0207, 0.04729, 0.04599, 0.03728, 0.0464, 0.06246, 0.05888, 0.06866, 0.01507, 0.03756]
        self.n = 0
        self.k = 0
        self.u = 0
        i = 0
        self.emissions = {'AA': 'ACDEFGHIKLMNPQRSTVWY', 'SS': 'HEGST_', 'TM': '12BHCILFU_', 'RAMA': 'HGBEdbeLlxc'}
        self.node_data = []
        #[{'AA': [0.08279, 0.01937, 0.05855, 0.05992, 0.04014, 0.08089, 0.02275, 0.05552, 0.05959, 0.0802, 0.0207, 0.04729, 0.04599, 0.03728, 0.0464, 0.06246, 0.05888, 0.06866, 0.01507, 0.03756], 'SS': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'TM': [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'RAMA': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'DIH': ['-75.00000', '-15.00000', '180.00000']}] # THIS IS THE NAUGHT NODE INFORMATION
        self.priors = []
        self.a = None
        self.comments = ""
        self.AA_ground = [0.08279, 0.01937, 0.05855, 0.05992, 0.04014, 0.08089, 0.02275, 0.05552, 0.05959, 0.08020, 0.02070, 0.04729, 0.04599, 0.03728, 0.04640, 0.06246, 0.05888, 0.06866, 0.01507, 0.03756] 
        self.SS_ground = [0.40672, 0.28659, 0.05248, 0.10787, 0.14634, 0.00000] 
        self.TM_ground = [0.99382, 0.00386, 0.00014, 0.00107, 0.00006, 0.00014, 0.00003, 0.00011, 0.00077, 0.00000] 
        self.RAMA_ground = [0.41779, 0.08362, 0.17662, 0.18220, 0.01939, 0.02236, 0.03537, 0.02220, 0.02446, 0.01223, 0.00378]  

    def add_node_(self, prior = 0, 
            AA_ground = [0.08279, 0.01937, 0.05855, 0.05992, 0.04014, 0.08089, 
                0.02275, 0.05552, 0.05959, 0.08020, 0.02070, 0.04729, 0.04599,
                0.03728, 0.04640, 0.06246, 0.05888, 0.06866, 0.01507, 0.03756],
            SS_ground = [0.40672, 0.28659, 0.05248, 0.10787, 0.14634, 0.00000],
            TM_ground = [0.99382, 0.00386, 0.00014, 0.00107, 0.00006, 0.00014, 0.00003, 0.00011, 0.00077, 0.00000],
            RAMA_ground = [0.41779, 0.08362, 0.17662, 0.18220, 0.01939, 0.02236, 0.03537, 0.02220, 0.02446, 0.01223, 0.00378],
            DIH = [-75.00000, -15.00000, 180.00000]):
        dat = {"AA": AA_ground, "SS": SS_ground, "TM":TM_ground, "RAMA":RAMA_ground,"DIH":DIH}
        self.add_node(dat, prior = prior)
    def add_node(self, node_data, prior = 0):
        self.n += 1
        self.priors.append(prior)
       # data = { "AA" : AA_ground, "SS": SS_ground, "RAMA": RAMA_ground, "TM": TM_ground, "DIH": DIH } 
        self.node_data.append(node_data)
        new_A = np.zeros((self.n, self.n))
        for i in range(0, self.n-1):
            for j in range(0, self.n-1):
                if self.a[i][j] > 0 :
                    new_A[i][j] = self.a[i][j]
        self.a = new_A
    def get_node_by_label(self, label):
        for i in range(0, self.n):
            try:
                if self.node_data[i]["label"][0] == label:
                    return i
            except:
                print(self.node_data)
                #print("ERROR, bad label. label = " + self.node_data[i]["label"])
        return -1

    # inserts new naught node after the last naught node. 
    def add_naught_node(self, prior = 0, 
            AA_ground = [0.08279, 0.01937, 0.05855, 0.05992, 0.04014, 0.08089, 0.02275, 0.05552, 0.05959, 0.08020, 0.02070, 0.04729, 0.04599, 0.03728, 0.04640, 0.06246, 0.05888, 0.06866, 0.01507, 0.03756],
            SS_ground = [0.40672, 0.28659, 0.05248, 0.10787, 0.14634, 0.00000],
            TM_ground = [0.99382, 0.00386, 0.00014, 0.00107, 0.00006, 0.00014, 0.00003, 0.00011, 0.00077, 0.00000],
            RAMA_ground = [0.41779, 0.08362, 0.17662, 0.18220, 0.01939, 0.02236, 0.03537, 0.02220, 0.02446, 0.01223, 0.00378],
            DIH = [-75.00000, -15.00000, 180.00000]):
        self.n += 1
        self.k += 1
        self.priors.insert(self.k-1, prior)
        data = { "AA": AA_ground, "SS": SS_ground, "RAMA": RAMA_ground, "TM": TM_ground, "DIH": DIH }
        self.node_data.insert(self.k-1, data)
        new_A = np.zeros((self.n, self.n))
        for i in range(0, self.n-1):
            for j in range(0, self.n-1):
                if self.a[i][j] > 0:
                    ti = i
                    if i >= self.k-1:
                        ti = i + 1
                    tj = j 
                    if j >= self.k-1:
                        tj = j + 1
                    new_A[ti][tj] = self.a[i][j]
        self.a = new_A
        print("NEW N = " + str(self.n))


    def insert_node(self, node_data, nitr): # nitr is the new index of the node being added
        self.n += 1
        self.priors.insert(nitr, 0)
        self.node_data.insert(nitr, node_data)
        new_A = np.zeros((self.n, self.n))
        for i in range(0, self.n-1):
            for j in range(0, self.n-1):
                if self.a[i][j] > 0:
                    ti = i
                    if i >= nitr:
                        ti = i + 1
                    tj = j
                    if j >= nitr:
                        tj = j + 1
                    new_A[ti][tj] = self.a[i][j]
        self.a = new_A

    def add_transition(self, i, j, outtrans_prob) : 
        self.a[i][j] = outtrans_prob

    def init_transition_matrix(self):
        self.a = np.zeros((self.n, self.n))

    def normalize_a(self):
        for i in range(0, self.n):
            s = 0
            for j in range(0, self.n):
                s += self.a[i][j]
            for j in range(0, self.n):
                self.a[i][j] = self.a[i][j] / s
        
    def remove_node(self, nitr):
        new_A = np.zeros((self.n-1, self.n-1))
        for i in range(0, self.n):
            ti = i 
            if i > nitr:
                ti = i - 1
            for j in range(0, self.n):
                tj = j 
                if j > nitr: 
                    tj = j - 1
                if i == nitr or j == nitr: 
                    continue
                new_A[ti][tj] = self.a[i][j]
        self.a = new_A
        self.node_data = [ dat for i,dat in enumerate(self.node_data) if i != nitr]
        self.priors = [ dat for i,dat in enumerate(self.priors) if i != nitr]
        self.n = self.n -1

    def merge(self, HMM2, naught_prob):
        old_n = self.n
        for i in range(1, len(HMM2.node_data)):
            state = HMM2.node_data[i]
            if state["TM"][0] > 0.005:
                continue
            self.add_node(0, state["AA"], state["SS"], state["RAMA"], state["TM"], state["DIH"])
        for i in range(1, HMM2.n): # don't add transitions to naught state
            for j in range(1, HMM2.n):
                if HMM2.a[i][j] > 0 :
                    self.add_transition(i-1 + old_n, j-1 + old_n, HMM2.a[i][j])
        # find states which are transitioned to from the naught node:
        for i in range(1, HMM2.n):
            if HMM2.a[0][i] > 0 :
                self.a[0][i+old_n-1] = naught_prob
                num_normalize = 0
                for j in range(0, self.n):
                    if self.a[0][j] > 0:
                        num_normalize += 1
                for j in range(0, self.n):
                    if self.a[0][j] > 0:
                        self.a[0][j] = self.a[0][j] - naught_prob / num_normalize
        # find states which transition back to naught node:
        for i in range(1, HMM2.n):
            if HMM2.a[i][0] > 0:
                self.a[i+old_n-1][0] = HMM2.a[i][0]

########################################################################################################################################
####################### TM TESTING #####################################################################################################
########################################################################################################################################



    def concurrency_histogram(self, name = "TEST"): # run 10000 out of naught.. create histogram of length before returning to naught.
        
        lengths = []
        for itr in tqdm(range(0,10000)):
            i = 0
            node = self.next_node()
            while node != 0:
                i += 1
                #print("NODE " + str(node))
                node = self.next_node(node)
            lengths.append(i)
        
        plt.style.use('seaborn-white')
        bins = 30 #CONSTRAINT_STYLES[name]['bins']
        rnge = (0,30) #CONSTRAINT_STYLES[name]['range']
        density = True # CONSTRAINT_STYLES[name]['density']
        
        plt.hist(lengths, bins, density = density, alpha = 0.75, range = rnge)
        plt.xlabel(str(name))
        if density:
            plt.ylabel('% occurence')
        else:
            plt.ylabel('# of times occured')
        plt.title('Histogram of ' + str(name))
        #print(len(data))
        plt.grid(True)
        plt.savefig(self.directory + name + '.png')
        plt.show()

    def next_node(self, current_node = 0): # start @ naught by default
        # get all possible transitions:
        transitions = self.a[current_node]
        # get random number [0,1]
        i = random.random()
        l = 0
        #print("i = " + str(i))
        for j,x in enumerate(transitions):
            l += x
            #print(l)
            if l > i:
                return j
        # if for some reason we didn't find it, return most likely next node.
        return np.argmax(transitions) 

    def clone_and_perturb(self, clone_index, new_trans):
        # update new_trans:
        new_new_trans = []
        for trans in new_trans:
            i = trans[0]
            j = trans[1]
            if i > clone_index:
                i += 1
            if j > clone_index:
                j += 1
            new_new_trans.append((i,j))

        new_a = np.zeros((self.n+1, self.n+1))
        # split all transitions to and from the cloned states:
        for itr in range(0, self.n):
            if itr > clone_index:
                i = itr + 1
            else:
                i = itr
            if i == clone_index:
                for jitr in range(0, self.n):
                    if jitr > clone_index:
                        j = jitr + 1
                    else:
                        j = jitr
                    if self.a[clone_index][jitr] > 0:
                        new_a[clone_index][j] = self.a[clone_index][jitr] * 0.5
                        new_a[clone_index+1][j] = self.a[clone_index][jitr] * 0.5
                        new_new_trans.append((i+1,j))
            else:
                for jitr in range(0, self.n):
                    if jitr > clone_index:
                        j = jitr + 1
                    else:
                        j = jitr
                    if self.a[itr][clone_index] > 0 and clone_index == j:
                        new_a[i][j] = self.a[itr][clone_index]*0.5
                        new_a[i][clone_index+1] = self.a[itr][clone_index]*0.5
                        new_new_trans.append((i,clone_index+1))
                    else:
                        new_a[i][j] = self.a[itr][jitr]

        new_node_data = []
        new_priors = []
        # add the clone to the node data list
        for i in range(0, self.n):
            if i == clone_index: # if it's the clone, add another copy
                new_priors.append(self.priors[i] * 0.5 )
                l = self.priors[i]
                new_priors.append(l * 0.5)
                new_node_data.append(self.node_data[i])
                copy_node_data = {}
                for f in list(self.node_data[i]):
                    copy_node_data[f] = [l for l in self.node_data[i][f]]
                new_node_data.append(copy_node_data)
            else:
                new_priors.append(self.priors[i])
                new_node_data.append(self.node_data[i])

        # update class variables
        self.a = new_a
        self.node_data = new_node_data
        self.priors = new_priors
        self.n = self.n + 1

        # greedily perturb the clone and it's copy in opposite directions (RAMA emission)
        original_emission = self.node_data[clone_index]["RAMA"]
        clone_emission = self.node_data[clone_index+1]["RAMA"]
        original_emission[0] = 0
        oe = original_emission
        # increase original's highest rama emisison:
        oe_RAMA_bins = [ oe[0] + oe[1], oe[2] + oe[5], oe[3] + oe[4] + + oe[6], oe[7] + oe[8], oe[9], oe[10]]
        sorted_oe_RAMA_bins = [ f for f in oe_RAMA_bins ]
        sorted_oe_RAMA_bins.sort()
        sorted_oe_RAMA_bins.reverse()
        print("ORIGINAL EMISSION \n" + str(oe))
        highest_bin = oe_RAMA_bins.index(max(oe_RAMA_bins))
        second_highest_bin = oe_RAMA_bins.index(sorted_oe_RAMA_bins[1])
        # PERTURB ( GENTLER )
        if False:
            if highest_bin == 0:
                original_emission[0] += 1
                original_emission[1] += 1
            elif highest_bin == 1:
                original_emission[2] += 1
                original_emission[5] += 1
            elif highest_bin == 2:
                original_emission[3] += 1
                original_emission[4] += 1
                original_emission[6] += 1
            elif highest_bin == 3:
                original_emission[7] += 1
                original_emission[8] += 1
            elif highest_bin == 4:
                original_emission[9] += 1
            elif highest_bin == 5:
                original_emission[10] += 1
            # increase clone's second highest rama emssion: 
            if second_highest_bin == 0:
                clone_emission[0] += 1
                clone_emission[1] += 1
            elif second_highest_bin == 1:
                clone_emission[2] += 1
                clone_emission[5] += 1
            elif second_highest_bin == 2:
                clone_emission[3] += 1
                clone_emission[4] += 1
                clone_emission[6] += 1
            elif second_highest_bin == 3:
                clone_emission[7] += 1
                clone_emission[8] += 1
            elif second_highest_bin == 4:
                clone_emission[9] += 1
            elif second_highest_bin == 5:
                clone_emission[10] += 1
        # PROD ( STRICTER ) 
        if True:
            if highest_bin == 0:
                clone_emission[0] = 0
                clone_emission[1] = 0
            elif highest_bin == 1:
                clone_emission[2] = 0
                clone_emission[5] = 0
            elif highest_bin == 2:
                clone_emission[3] = 0
                clone_emission[4] = 0
                clone_emission[6] = 0
            elif highest_bin == 3:
                clone_emission[7] = 0
                clone_emission[8] = 0
            elif highest_bin == 4:
                clone_emission[9] = 0
            elif highest_bin == 5:
                clone_emission[10] = 0
            if second_highest_bin == 0:
                original_emission[0] = 0
                original_emission[1] = 0
            elif second_highest_bin == 1:
                original_emission[2] = 0
                original_emission[5] = 0
            elif second_highest_bin == 2:
                original_emission[3] = 0
                original_emission[4] = 0
                original_emission[6] = 0
            elif second_highest_bin == 3:
                original_emission[7] = 0
                original_emission[8] = 0
            elif second_highest_bin == 4:
                original_emission[9] = 0
            elif second_highest_bin == 5:
                original_emission[10] = 0
        clone_emission = normalize(clone_emission)
        original_emission = normalize(original_emission)
        print("ORIGINAL EMISSION BIN = " + str(highest_bin) + " " + str(original_emission))
        print("CLONE EMISSION BIN = " + str(second_highest_bin) + " " + str(clone_emission))
        self.node_data[clone_index]["RAMA"] = original_emission
        self.node_data[clone_index+1]["RAMA"] = clone_emission
        return new_new_trans

