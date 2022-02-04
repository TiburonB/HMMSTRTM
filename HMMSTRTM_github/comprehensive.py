# comphrehensive.py
# TLB 10/12/2021
# Run all important diagnostics on a version of HMMSTR

import argparse
import pickle
from tplot import *
from hmm import *
from hmm_viz import *
import pickle
from tqdm import tqdm
from os import listdir
from os.path import isfile, isdir
import os
import numpy as np
from gamma import read as gamma_read
from gamma import write_gamma_heat
from confusion import write_to_file 
import math
SS_bins = [  "G/H","S/T/_", "E" ]
RAMA_bins = [ "G/H", "B/b", "E/e/d", "L/l", "x", "c" ]
TM_bins = [ "1/2", "B", "H", "C", "I", "L", "F", "U" ]

def is_TM(model_file):
    try:
        model_itr = int(model_file[model_file.rfind('_')+1:model_file.rfind('.')])
    except:
        model_itr = 1
    TM = False
    if model_itr == 1:
        ind = model_file.rfind('.')
        TM =  model_file[ind-2:ind] == "TM"
    else:
        ind = model_file.rfind('_')
        TM =  model_file[ind-2:ind] == "TM"
    return TM

def compile_get_call_str(model_file, code, gamma_run, test = True, drct = None):
    model_dir  = model_file[:model_file.rfind('/')+1]
    model_name =  model_file[model_file.rfind('/')+1:model_file.rfind('.')].replace(" ", "")
    try:
        model_itr = int(model_file[model_file.rfind('_')+1:model_file.rfind('.')])
    except:
        model_itr = 1
    TM = is_TM(model_file)
    # make a call to comprehensive.f95
    os.system("gfortran -c hmm_io.f95 drct2profile.f95 gamma_io.f95 hmm.f95 comprehensive.f95 -fcheck=all")
    os.system("gfortran hmm_io.o drct2profile.o gamma_io.o hmm.o comprehensive.o")
    if drct is not None:
        x = 0
    else:
        drct = ""
        if TM:
            if test:
                drct = "tm_testing.drct"
            else:
                drct = "tm_training.drct"
        else:
            if test:
                drct = "testing.drct"
            else:
                drct = "training.drct"
    if code is None:
        s = './a.out \'' + model_file + '\' ' + drct
    else:
        s = './a.out \'' + model_file + '\' ' + drct + ' ' + code
    if gamma_run:
        s += " T"
    else:
        s += " F" # false flag to not printout gamma matrix to terminal.
    print(s)
    return s

def parse(data):
    # retrieve data parsed as
    # START < code > 
    # < t > 
    # < mda accuracy > 
    # SS_CONFUSION ( 3 x 3 ) 
    # RAMA_CONFUSION ( 6 x 6 ) 
    # TM_CONFUSION   ( 8 x 8 )
    # RAMA_PROFILE   ( t x 11 ) 
    # TM_PROFILE   ( t x 10 ) 
    # RAMA_SEQ ( t )
    # TM_SEQ ( t )
    gen_data = ( k for k in data )
    codes = []
    ts = []
    mdas = []
    mdavs = []
    SS_confusion = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    RAMA_confusion = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
    TM_confusion =  [[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    RAMA_Profile = []
    TM_Profile = []
    AA_seq = []
    RAMA_seq = []
    RAMA_seq_pred = []
    TM_seq = []
    TM_seq_pred = []
    while True:
        try:
            k = next(gen_data)
        except:
            break
        if k == "START":
            code = next(gen_data)
            codes.append(code)
            t = int(next(gen_data))
            ts.append(t)
            mda = float(next(gen_data))
            mdavs.append(mda)
            mda = next(gen_data)
            mdas.append(mda)
            for i in range(0, 4): # get dssp (Q3) confusion
                for j in range(0,4):
                    SS_confusion[i][j] += int(next(gen_data))
            for i in range(0, 6): # get RAMA confusion
                for j in range(0,6):
                    RAMA_confusion[i][j] += int(next(gen_data))
            for i in range(0, 8): # get TM confusion
                for j in range(0,8):
                    TM_confusion[i][j] += int(next(gen_data))
            this_RAMA_Profile = []
            for titr in range(0, t):
                this_profile = []
                for j in range(0, 11):
                    this_profile.append(float(next(gen_data)))
                this_RAMA_Profile.append(this_profile)
            RAMA_Profile.append(this_RAMA_Profile)
            this_TM_Profile = []
            for titr in range(0, t):
                this_profile = []
                for j in range(0, 10):
                    this_profile.append(float(next(gen_data)))
                this_TM_Profile.append(this_profile)
            TM_Profile.append(this_TM_Profile)
            AA_seq.append(next(gen_data))
            RAMA_seq.append(next(gen_data))
            RAMA_seq_pred.append(next(gen_data))
            TM_seq.append(next(gen_data))
            TM_seq_pred.append(next(gen_data))

    return codes, ts, mdavs,mdas, SS_confusion, RAMA_confusion, TM_confusion, RAMA_Profile, TM_Profile, RAMA_seq, RAMA_seq_pred, TM_seq, TM_seq_pred


from Bio.PDB import *
from os.path import exists
from Atom import *
def structure_viz(model_file, drct_file, code, chain):

    s = compile_get_call_str(model_file, code, False, test = False, drct = drct_file) 
    data = os.popen(s).read().split()
    model_dir  = model_file[:model_file.rfind('/')+1]
    model_name =  model_file[model_file.rfind('/')+1:model_file.rfind('.')].replace(" ", "")
    TM = is_TM(model_file)
    codes, ts, mdavs, mdas, SS_confusion, RAMA_confusion, TM_confusion, RAMA_Profile, TM_Profile, RAMA_seq, RAMA_seq_pred, TM_seq, TM_seq_pred = parse(data)

    try:
        if ( exists('./' + code.lower() + ".pdb") == False):
            path = '../../Base/Warehouse/Protein/'+code.upper()+'/PDB/'+code.lower()+'.pdb'
            s = 'cp ' + path +  ' ./'
            print(s)
            os.system(s)
    except:
        print("COULDN't FIND PDB FILE.")
    print(RAMA_seq[0])
    print(RAMA_seq_pred[0])
    print(mdas[0])
    AS = AtomSet.from_file('./'+code.lower()+'.pdb')
    AS.restrict_to_atoms(["CA"])
    ind = 0
    #print(AS)
    PTM = get_TMPRED_temp(TM_seq, TM_Profile, model_dir, model_name, code)
    #print(PTM)
    for i,A in enumerate(AS.atoms):
        if A.chain_id == "A":
            try:
                if mdas[0][ind] == "1": #RAMA_seq[0][ind] == RAMA_seq_pred[0][ind]:
                    AS.atoms[i].temp_factor = 50.00
                else:
                    AS.atoms[i].temp_factor = 20.00
            except:
                AS.atoms[i].temp_factor = 20.00
            AS.atoms[i].temp_factor = PTM[ind] * 30 + 20.00
            ind += 1
        else:
            AS.atoms[i].temp_factor = 20.0
    if chain == '':
        AS.write('./out.pdb', 'A')
    else:
        AS.write('./out.pdb', chain)
    raise

def get_comprehensive(model_file, code, testing = False, drct_file = None):
    if drct_file is not None and drct_file == "None":
        drct_file = None
    s = compile_get_call_str(model_file, code, False, test = testing, drct = drct_file) 
    data = os.popen(s).read().split()
    model_dir  = model_file[:model_file.rfind('/')+1]
    model_name =  model_file[model_file.rfind('/')+1:model_file.rfind('.')].replace(" ", "")
    TM = is_TM(model_file)
    codes, ts, mdas, mdaseq, SS_confusion, RAMA_confusion, TM_confusion, RAMA_Profile, TM_Profile, RAMA_seq, RAMA_seq_pred, TM_seq, TM_seq_pred = parse(data)
    # Get TM F+-/T+- Matrix
    small_TM_confusion = TM_small_confusion(TM_confusion)
    # PRINT CONFUSION MATRICES
    printout(TM_confusion, TM_bins, TM_bins)
    printout(small_TM_confusion, ["+", "-"], ["T", "F"])
    printout(RAMA_confusion, RAMA_bins, RAMA_bins)
    printout(SS_confusion, SS_bins, SS_bins)
    # GET MDA
    MDA = get_MDA(ts, mdas)
    # GET Q3 / RAMA / TM % acc
    Q3 = get_MAT_acc(SS_confusion)
    print("Q3 = {:>12.5f}".format(Q3))
    RAMA_acc = get_MAT_acc(RAMA_confusion)
    print("RAMAacc = {:>7.5f}".format(RAMA_acc))
    TM_acc = get_MAT_acc(TM_confusion)
    print("TMacc = {:>9.5f}".format(TM_acc))
    if code is not None: 
        # Get RAMAPRED
        #get_RAMAPRED(RAMA_seq, RAMA_Profile, model_dir, model_name, code)
        if TM:# Get TMPRED
            print(TM_seq)
            get_TMPRED(TM_seq, TM_Profile, model_dir, model_name, code)
    # Get TMROC
    if TM:
        get_TM_ROC_x(TM_seq, TM_Profile, model_name, model_dir, 25)

def get_MDA(ts, mdas):
    t_sum = sum(ts)
    weighted_mdas = []
    for i in range(0, len(ts)):
        weighted_mdas.append(ts[i] * mdas[i])
    mda = sum(weighted_mdas) / t_sum
    print("database mda = " + str(mda))


def get_MAT_acc(MAT):
    acc = 0
    s = 0
    for i in range(0, len(MAT)):
        s += sum(MAT[i])
        acc += MAT[i][i]
    acc = acc / s
    return acc

def TM_small_confusion(TM_confusion):
    # TM_bins = [ "1/2", "B", "H", "C", "I", "L", "F", "U" ]
    Tpos = 0
    Tneg = 0
    Fpos = 0
    Fneg = 0
    # True + = (GT, PRED) = ({ B, H, C, I, L, F },  { B, H, C, I, L, F })
    for i in range(1, 7):
        for j in range(1, 7):
            Tpos += TM_confusion[i][j]
    # True - = (GT, PRED) = ({ 1, 2 },  { 1, 2 })
    Tneg = TM_confusion[0][0]
    # False + = (GT, PRED) = ({ 1, 2 },  { B, H, C, I, L, F })
    for i in range(1, 7):
        Fpos += TM_confusion[i][0]
    # False - = (GT, PRED) = ({ B, H, C, I, L },  { 1, 2 })
    for j in range(1, 7):
        Fneg += TM_confusion[0][j]
    small_TM_confusion = [[Tpos, Tneg],[Fpos,Fneg]]
    return small_TM_confusion

def get_RAMAPRED(RAMASEQ, RAMAPROF, model_dir, model_name, code):
    RAMASEQ = RAMASEQ[0]
    RAMAPROF = RAMAPROF[0]
    PHG = []
    PBb = []
    PEed = []
    Pc = []
    Px = []
    PlL = []
    for k in RAMAPROF:
        if sum(k[:]) == 0:
            PHG.append(0)
            PBb.append(0)
            PEed.append(0)
            Pc.append(0)
            Px.append(0)
            PlL.append(0)
        else:
            PHG.append(min(1,max((k[0]+k[1]) / sum(k[:]),0)))
            PBb.append(min(1,max((k[2] + k[5])/sum(k[:]),0)))
            PEed.append(min(1,max((k[3] + k[4]+k[6])/sum(k[:]),0)))
            Pc.append(min(1,max(k[10]/sum(k[:]),0)))
            Px.append(min(1,max(k[9]/sum(k[:]),0)))
            PlL.append(min(1,max((k[7]+k[8])/sum(k[:]),0)))
    # ADVANCED PLOT T 10/11/21
    RAMASEQ = [ f for f in RAMASEQ ]
    view_RAMAPRED_plot(PHG,PBb,PEed,Pc,Px,PlL, model_dir, model_name +'RAMAPRED: '+code, RAMASEQ)

def get_TMPRED_temp(TMSEQ, TMPROF, model_dir, model_name, code):
    TMSEQ = TMSEQ[0]
    TMPROF = TMPROF[0]
    PTM_B = []
    PTM_H = []
    PTM_C = []
    PTM_I = []
    PTM_L = []
    PTM_F = []
    P_U = []
    PTM = []
    PINSIDE = [] # to be filled later
    POUTSIDE = [] # to be filled later
    INSIDE_LAST = True

    
    for i,f in enumerate(TMPROF):
        PTM_B.append(f[2])
        PTM_H.append(f[3])
        PTM_C.append(f[4])
        PTM_I.append(f[5])
        PTM_L.append(f[6])
        PTM_F.append(f[7])
        P_U.append(f[8])
        PTM.append(sum(f[1:]))
        if i == 0:
            PINSIDE.append(1)
            POUTSIDE.append(0)
        else:
            PTM_dif = PTM[i] - PTM[i-1] # if > 0, PTM increased
            if PTM_dif > 0:
                if PINSIDE[i-1] > POUTSIDE[i-1]:
                    PINSIDE.append(PINSIDE[i-1] - abs(PTM[i] - PTM[i-1]))
                    POUTSIDE.append(POUTSIDE[i-1])
                    if PTM[i-1] != max(PTM[i-1], PINSIDE[i-1], POUTSIDE[i-1]):
                        INSIDE_LAST = True
                else:
                    POUTSIDE.append(POUTSIDE[i-1] - abs(PTM[i] - PTM[i-1]))
                    PINSIDE.append(PINSIDE[i-1])
                    if PTM[i-1] != max(PTM[i-1], PINSIDE[i-1], POUTSIDE[i-1]):
                        INSIDE_LAST = False
            else: # PTM decreasing
                if INSIDE_LAST:
                    if PINSIDE[i-1] == max(PTM[i-1], PINSIDE[i-1], POUTSIDE[i-1]):
                        PINSIDE.append(PINSIDE[i-1] + abs(PTM[i] - PTM[i-1]))
                        POUTSIDE.append(POUTSIDE[i-1])
                    else:
                        POUTSIDE.append(POUTSIDE[i-1] + abs(PTM[i] - PTM[i-1]))
                        PINSIDE.append(PINSIDE[i-1])
                else:
                    if POUTSIDE[i-1] == max(PTM[i-1], PINSIDE[i-1], POUTSIDE[i-1]):
                        POUTSIDE.append(POUTSIDE[i-1] + abs(PTM[i] - PTM[i-1]))
                        PINSIDE.append(PINSIDE[i-1])
                    else:
                        PINSIDE.append(PINSIDE[i-1] + abs(PTM[i] - PTM[i-1]))
                        POUTSIDE.append(POUTSIDE[i-1])
    # ADVANCED PLOT T 10/11/21
    TM_c = [ f for f in TMSEQ ]
    return PTM


def get_TMPRED(TMSEQ, TMPROF, model_dir, model_name, code):
    TMSEQ = TMSEQ[0]
    TMPROF = TMPROF[0]
    PTM_B = []
    PTM_H = []
    PTM_C = []
    PTM_I = []
    PTM_L = []
    PTM_F = []
    P_U = []
    PTM = []
    PINSIDE = [] # to be filled later
    POUTSIDE = [] # to be filled later
    INSIDE_LAST = True

    
    for i,f in enumerate(TMPROF):
        PTM_B.append(f[2])
        PTM_H.append(f[3])
        PTM_C.append(f[4])
        PTM_I.append(f[5])
        PTM_L.append(f[6])
        PTM_F.append(f[7])
        P_U.append(f[8])
        PTM.append(sum(f[1:]))
        if i == 0:
            PINSIDE.append(1)
            POUTSIDE.append(0)
        else:
            PTM_dif = PTM[i] - PTM[i-1] # if > 0, PTM increased
            if PTM_dif > 0:
                if PINSIDE[i-1] > POUTSIDE[i-1]:
                    PINSIDE.append(PINSIDE[i-1] - abs(PTM[i] - PTM[i-1]))
                    POUTSIDE.append(POUTSIDE[i-1])
                    if PTM[i-1] != max(PTM[i-1], PINSIDE[i-1], POUTSIDE[i-1]):
                        INSIDE_LAST = True
                else:
                    POUTSIDE.append(POUTSIDE[i-1] - abs(PTM[i] - PTM[i-1]))
                    PINSIDE.append(PINSIDE[i-1])
                    if PTM[i-1] != max(PTM[i-1], PINSIDE[i-1], POUTSIDE[i-1]):
                        INSIDE_LAST = False
            else: # PTM decreasing
                if INSIDE_LAST:
                    if PINSIDE[i-1] == max(PTM[i-1], PINSIDE[i-1], POUTSIDE[i-1]):
                        PINSIDE.append(PINSIDE[i-1] + abs(PTM[i] - PTM[i-1]))
                        POUTSIDE.append(POUTSIDE[i-1])
                    else:
                        POUTSIDE.append(POUTSIDE[i-1] + abs(PTM[i] - PTM[i-1]))
                        PINSIDE.append(PINSIDE[i-1])
                else:
                    if POUTSIDE[i-1] == max(PTM[i-1], PINSIDE[i-1], POUTSIDE[i-1]):
                        POUTSIDE.append(POUTSIDE[i-1] + abs(PTM[i] - PTM[i-1]))
                        PINSIDE.append(PINSIDE[i-1])
                    else:
                        PINSIDE.append(PINSIDE[i-1] + abs(PTM[i] - PTM[i-1]))
                        POUTSIDE.append(POUTSIDE[i-1])
    # ADVANCED PLOT T 10/11/21
    TM_c = [ f for f in TMSEQ ]
    #view_TMPRED_plot(PTM, PINSIDE, POUTSIDE, model_dir, model_name + 'TMPRED: '+code, TM_c)
    view_TMPRED_plot_more(PTM_B, PTM_H, PTM_C, PTM_I, PTM_L, PTM_F, P_U, PINSIDE, POUTSIDE, model_dir, model_name + 'TMPRED: '+code, TM_c)

import math
def get_TM_ROC_x(TMSEQ, TMPROF, model_name, model_dir, x):
    ROCS = [ [] for i in range(0, x) ]
    for i in range(0, len(TMSEQ)):
        for j in range(0, len(TMSEQ[i])):
            TMness = sum(TMPROF[i][j][1:])
            confidence = (max(TMness, 1 - TMness) - 0.5) * 2.0
            TM = (TMSEQ[i][j] in ["H", "B", "I", "L", "C", "F", "U" ])
            correct = (( TMness < 0.5 and not(TM)) or ( TMness >= 0.5 and TM))
            ind = math.floor(x - confidence * x)
            ROCS[ind].append((confidence, correct, TM))
    nf = 0
    nt = 0
    roc = 0
    roc_dat = []
    bin_splits = []
    for i,ret_ROC in enumerate(ROCS):
        for dat in ret_ROC:
            #print(dat)
            if dat[1]:
                nt = nt + 1
            else:
                roc = roc + nt
                nf = nf + 1
        bin_splits.append(int(len(ret_ROC)/50))
#        print(sum(bin_splits))
        for j in range(0, sum(bin_splits)):
            roc_dat.append(i+1)
    roc = roc / ( nt * nf )
    print("TM ROC = " + str(roc))
    bin_splits = [ i for i in range(1,11)]
    try:
        view_TM_ROC(roc_dat, model_name+"TM_ROC",model_dir, x, bin_splits)
    except:
        print("Failed to visualized TMROC.")

# TM ROC unsplit
def get_TM_ROC(TMSEQ, TMPROF, model_name, model_dir):
    n = 0
    for SEQ in TMSEQ:
        for c in SEQ:
            n += 1
    print(n)
    if n < 1000:
        get_TM_ROC_x(TMSEQ, TMPROF, model_name, model_dir, 100)
    elif n < 100000:
        get_TM_ROC_x(TMSEQ, TMPROF, model_name, model_dir, 50)
    elif n < 10000000:
        get_TM_ROC_x(TMSEQ, TMPROF, model_name, model_dir, 25)

# printout a matrix given the data and axis titles
def printout( MAT, xaxis, yaxis ):
    s = "     "
    for f in xaxis:
        s += '{:>9s}'.format(f)
    s += "\n"
    for i in range(0, len(yaxis)):
        s += '{:>9s}'.format(yaxis[i])
        for j in range(0, len(xaxis)):
            s += '{:>9}'.format(MAT[i][j])
        s += "\n"
    print(s)

def gamma_heat(model_file, code, testing = False):
    if ( code is not None ) :  # if code is specitified
        if ( isfile('./gamma/'+code) ): # if gamma file already generated...             
            gamma = gamma_read(code, './gamma/')
            model_dir  = model_file[:model_file.rfind('/')+1]
            model_name =  model_file[model_file.rfind('/')+1:model_file.rfind('.')].replace(" ", "")
            HMM = hmm(model_file[model_file.rfind('/')+1:], model_dir)
            HMM.read()
            naught_t_step(HMM, gamma)
        else:
            get_comprehensive(model_file, code) # if gamma_file not found, generate & recurse
            gamma_heat(model_file, code)
            return
    else: # code isn't specified... we want gamma across the database.
        s = compile_get_call_str(model_file,code, True, test = testing)
        data = os.popen(s).read().split()

        def parse(data):
            # retrieve data parsed as
            # START < code > 
            # < t > 
            # < n >
            # < GAMMA MATRIX ( t x n ) 
            gen_data = ( k for k in data )
            t = -1
            n = -1
            gamma = []
            while True:
                try:
                    k = next(gen_data)
                except:
                    break
                if k == "START":
                    code = next(gen_data)
                    t = int(next(gen_data))
                    weight = float(next(gen_data))
                    n = int(next(gen_data))
                    # initialize gamma array
                    if len(gamma) == 0:
                        for nitr in range(0, n):
                            gamma.append(0)
                    for titr in range(0, t):
                        for nitr in range(0, n):
                            gamma[nitr] += float(next(gen_data)) * weight
            return gamma
        s = s.split()
        if code is None:
            drct_file = s[len(s)-2]
        else: 
            drct_file = s[len(s) - 3]
        print('check drct file = ' + drct_file)
        gamma = parse(data)
        model_name =  model_file[model_file.rfind('/')+1:model_file.rfind('.')].replace(" ", "")
        write_gamma_heat(model_name+ '.gamma', './gamma/', gamma, model_file, drct_file)

def confusion_reweight(model_file, testing = False, drct_file = None): 
    # iteratively find the reweighting values for RAMA and SS confusion matrices. 
    model_dir  = model_file[:model_file.rfind('/')+1]
    model_name =  model_file[model_file.rfind('/')+1:model_file.rfind('.')].replace(" ", "")
    mf = model_file[model_file.rfind('/')+1:]
    HMMSTR = hmm(mf, model_dir)
    HMMSTR.read()
    HMMSTR.write(model_file)
    print("TM REWEIGHT " + str(HMMSTR.tm_reweight))
    print("RAMA REWEIGHT " + str(HMMSTR.rama_reweight))
    print("SS REWEIGHT " + str(HMMSTR.ss_reweight))
    s = compile_get_call_str(model_file, None, False, test = testing, drct = drct_file)
    data = os.popen(s).read().split()

    def parse(data):
        # retrieve data parsed as
        # START < code > 
        # < t > 
        # < mda accuracy > 
        # SS_CONFUSION ( 3 x 3 ) 
        # RAMA_CONFUSION ( 6 x 6 ) 
        # TM_CONFUSION   ( 8 x 8 )
        # RAMA_PROFILE   ( t x 11 ) 
        # TM_PROFILE   ( t x 10 ) 
        # RAMA_SEQ ( t )
        # TM_SEQ ( t )
        gen_data = ( k for k in data )
        codes = []
        ts = []
        mdas = []
        mdavs = []
        SS_confusion = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
        RAMA_confusion = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
        TM_confusion =  [[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
        RAMA_Profile = []
        TM_Profile = []
        AA_seq = []
        RAMA_seq = []
        RAMA_seq_pred = []
        TM_seq = []
        TM_seq_pred = []
        while True:
            try:
                k = next(gen_data)
            except:
                break
            if k == "START":
                code = next(gen_data)
                codes.append(code)
                t = int(next(gen_data))
                ts.append(t)
                mda = float(next(gen_data))
                mdavs.append(mda)
                mda = next(gen_data)
                mdas.append(mda)
                for i in range(0, 4): # get dssp (Q3) confusion
                    for j in range(0,4):
                        SS_confusion[i][j] += int(next(gen_data))
                for i in range(0, 6): # get RAMA confusion
                    for j in range(0,6):
                        RAMA_confusion[i][j] += int(next(gen_data))
                for i in range(0, 8): # get TM confusion
                    for j in range(0,8):
                        TM_confusion[i][j] += int(next(gen_data))
                this_RAMA_Profile = []
                for titr in range(0, t):
                    this_profile = []
                    for j in range(0, 11):
                        this_profile.append(float(next(gen_data)))
                    this_RAMA_Profile.append(this_profile)
                RAMA_Profile.append(this_RAMA_Profile)
                this_TM_Profile = []
                for titr in range(0, t):
                    this_profile = []
                    for j in range(0, 10):
                        this_profile.append(float(next(gen_data)))
                    this_TM_Profile.append(this_profile)
                TM_Profile.append(this_TM_Profile)
                AA_seq.append(next(gen_data))
                RAMA_seq.append(next(gen_data))
                RAMA_seq_pred.append(next(gen_data))
                TM_seq.append(next(gen_data))
                TM_seq_pred.append(next(gen_data))
        return SS_confusion, RAMA_confusion, TM_confusion, AA_seq, TM_seq, TM_seq_pred, codes

    SS_confusion, RAMA_confusion, TM_confusion, AA_seq, TM_seq, TM_seq_pred, codes = parse(data)
    printout(TM_confusion, TM_bins, TM_bins)
    printout(RAMA_confusion, RAMA_bins, RAMA_bins)
    printout(SS_confusion, SS_bins, SS_bins)


    if False:
        TM_H_segment_ranges = [] 
        # 3d tuple. 1st index is chain #, 2nd and 3rd are start and end indeces within the chain.
        for i,TM in enumerate(TM_seq):
            end = -1
            start = -1
            for j,c in enumerate(TM):
                if c == "H" and TM[j-1] != "H":
                    start = j
                elif c == "H" and TM[j+1] != "H":
                    end = j
                if end != -1 and start != -1:
                    TM_H_segment_ranges.append((i,start,end))
                    start = -1
                    end = -1
        print(TM_H_segment_ranges)
        import pickle
        #pickle.dump(TM_H_segment_ranges, open('./TMH_ranges.pickle', 'wb'))

        missed_TM_H_segments = []
        missed_codes = []
        missed = 0
        hit = 0
        for i,r in enumerate(TM_H_segment_ranges):
            prediction = TM_seq_pred[r[0]][r[1]:r[2]]
            num_H = prediction.count('H')
            pmatch = num_H / len(prediction)
            if pmatch < .1:
                missed += 1
                missed_codes.append(codes[r[0]])
                AA_ = AA_seq[r[0]][r[1]:r[2]]
                missed_TM_H_segments.append(AA_)
            else:
                hit += 1
        print("NUMBER TMH hits = " + str(hit))
        print("NUMBER TMH misses = " + str(missed))
        #pickle.dump(missed_TM_H_segments, open('./missed_TMAAs.pickle', 'wb'))
        #pickle.dump(missed_codes, open('./missed_codes.pickle', 'wb'))
    
        max_length =  max([len(x) for x in missed_TM_H_segments])
        avg_prof = [ [0 for i in range(0,20)] for j in range(0, max_length) ] 
        alpha = HMMSTR.emissions["AA"]
        for AA_ in missed_TM_H_segments:
            for i,c in enumerate(AA_):
                j = alpha.index(c)
                avg_prof[i][j] += 1
        print(avg_prof)
        raise
                
    
    if True:
        write_to_file(model_dir + model_name + "SS_CONFUSION.txt", SS_confusion)
        write_to_file(model_dir + model_name + "RAMA_CONFUSION.txt", RAMA_confusion)
        # print SS reweight, RAMA reweight 1-d matrices
        ss_converge = True
        rama_converge = True
        tm_converge = True
        ss_reweight = []
        for i in range(0 , len(SS_confusion)-1):
            num = 0.0
            denom = 0.0 
            for j in range(0, len(SS_confusion[i])-1):
                num += SS_confusion[i][j]
                denom += SS_confusion[j][i]
            if denom != 0:
                reweight = 1+math.log10(num/denom)
                print(SS_bins[i] + " : " + str(num/denom) + " ... " + str(reweight)) 
                ss_reweight.append(reweight)
                if reweight > 1.05 or reweight < 0.95:
                    ss_converge = False
            else:
                ss_reweight.append(4)
                ss_converge = False
                print(SS_bins[i] + " : INF")
        # print SS reweight, RAMA reweight 1-d matrices
        rama_reweight = []
        for i in range(0 , len(RAMA_confusion)):
            num = 0.0
            denom = 0.0
            for j in range(0, len(RAMA_confusion[i])):
                num += RAMA_confusion[i][j]
                denom += RAMA_confusion[j][i]
            if denom != 0:
                reweight = 1+math.log10(num/denom)
                print(RAMA_bins[i] + " : " + str(num/denom) + " ... " + str(reweight))
                rama_reweight.append(reweight)
                if reweight > 1.05 or reweight < 0.95:
                    rama_converge = False
            else:
                rama_converge = False
                rama_reweight.append(4)
                print(RAMA_bins[i] + " : INF")
        tm_reweight = []
        for i in range(0 , len(TM_confusion)):
            num = 0.0
            denom = 0.0
            for j in range(0, len(TM_confusion[i])):
                num += TM_confusion[i][j]
                denom += TM_confusion[j][i]
            if denom != 0:
                reweight = 1+math.log10(num/denom)
                print(TM_bins[i] + " : " + str(num/denom) + " ... " + str(reweight))
                tm_reweight.append(reweight)
                if reweight > 1.05 or reweight < 0.95:
                    rama_converge = False
            else:
                tm_converge = False
                tm_reweight.append(4)

        print("SS REWEIGHT =  " + str(ss_reweight))
        print("RAMA REWEIGHT = " + str(rama_reweight))
        print("TM REWEIGHT = " + str(tm_reweight))
        printout(RAMA_confusion, RAMA_bins, RAMA_bins)
        printout(SS_confusion, SS_bins, SS_bins)
        printout(TM_confusion, TM_bins, TM_bins)

    if rama_converge == False or ss_converge == False or tm_converge == False:
        for j in range(0, len(RAMA_confusion)):
            HMMSTR.rama_reweight[j] = HMMSTR.rama_reweight[j] *  rama_reweight[j]
        for j in range(0, len(SS_confusion)-1):
            HMMSTR.ss_reweight[j] = HMMSTR.ss_reweight[j] * ss_reweight[j]
        for j in range(0, len(TM_confusion)-1):
            HMMSTR.tm_reweight[j] = HMMSTR.tm_reweight[j] * tm_reweight[j]
        HMMSTR.write(model_file)
        confusion_reweight(model_file, testing = testing)
    else:
        print("CONVERGED.")
    

def parse_args():
    parser = argparse.ArgumentParser(description = 'diagnostics.py')
    parser.add_argument('-mf', dest = 'model_file', default = "./HMMS/HMMSTR.hmm")
    parser.add_argument('-df', dest = 'drct_file', default = "None")
    parser.add_argument('-code', dest = 'code', default = "None")
    parser.add_argument('-GH', dest = 'GH', default = "False")
    parser.add_argument('-T', dest = 'TEST', default = "False")
    parser.add_argument('-c', dest = 'CONFUSION', default = "False")
    parser.add_argument('-ch', dest = 'chain', default = "")
    parser.add_argument('-s', dest = 'STRUCTURE', default = "False")
    args = parser.parse_args()
    return args

def str_b(v):
    if v in [ "True", "T", "true", "t", "Y", "y" ]:
        return True
    else:
        return False

if __name__ == '__main__':
    args = parse_args()
    model_file = args.model_file
    GH = str_b(args.GH)
    code = args.code
    TEST = str_b(args.TEST)
    c = str_b(args.CONFUSION)
    chain = args.chain
    s = str_b(args.STRUCTURE)
    drct_file = args.drct_file
    if drct_file == "None":
        drct_file = None
    if s:
        structure_viz(model_file, drct_file, code, chain)
        raise
    if code == "None":
        code = None
    if c:
        confusion_reweight(model_file, testing = TEST, drct_file = drct_file)
    if GH:
        gamma_heat(model_file, code, testing = TEST)
    else:
        get_comprehensive(model_file, code, testing = TEST, drct_file = drct_file)
