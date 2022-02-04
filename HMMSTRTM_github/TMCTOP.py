# TMCTOP.py
#  1/7/2022. Culminative.
from hmm import *
import os
import pickle
from datamine import * 
glob = hmm('HMMSTR_150.hmm', directory = './HMMSTR/') 
glob.read()
# view the globular model.
from hmm_viz import display_by_RAMA_AA, display_by_TM, display
    #display_by_RAMA_AA(glob,QUIET=False)
from hmm_tools import * 

# Retrieve a globular model ( HMMSTR2K20 )
GLOB_1 = hmm('HMMSTR_150.hmm', directory = './HMMSTR/')
GLOB_1.read()

CANON = hmm('HMMSTRTM.hmm', directory = './TMHTOP/')
CANON.init_HMMSTR()
old_canon = hmm('HMMSTRTM.hmm', directory = './')
old_canon.read()
# copy over important hyperparameters. 
CANON.AA_ground = old_canon.AA_ground
CANON.SS_ground = old_canon.SS_ground
CANON.TM_ground = old_canon.TM_ground
CANON.RAMA_ground = old_canon.RAMA_ground
CANON.AA_EPS = old_canon.AA_EPS
CANON.SS_EPS = old_canon.SS_EPS
CANON.TM_EPS = old_canon.TM_EPS
CANON.RAMA_EPS = old_canon.RAMA_EPS

alpha = CANON.emissions["AA"]
tmalpha = CANON.emissions["TM"]

for i in range(0, GLOB_1.n):
    GLOB_1.node_data[i]["label"] = [ "GLOBULAR_SIDE_1" ]
    for j in range(0, len(GLOB_1.node_data[i]["TM"])):
        if j == 0: # replace TM emission 1 w/ 100% prob
            GLOB_1.node_data[i]["TM"][j] = 1.0
        else:
            GLOB_1.node_data[i]["TM"][j] = 0.1 # not exactly 100%

for i,node in enumerate(GLOB_1.node_data):
    if i == 0: 
        CANON.add_naught_node(GLOB_1.priors[i], node["AA"], node["SS"], node["TM"], node["RAMA"])
        CANON.node_data[0]["label"] = ["GLOBULAR_SIDE_1_NAUGHT"]
    else:
        CANON.add_node(node, prior = GLOB_1.priors[i] )
        CANON.node_data[CANON.n-1]["label"] = ["GLOBULAR_SIDE_1"]

for i in range(0, GLOB_1.n):
    for j in range(0, GLOB_1.n):
        t = GLOB_1.a[i][j]
        if t > 0:
            CANON.a[i][j] = t
Helix_emission =  [ 0.1, 0.01, 0.01, .01, 0.83, 0.01, 0.01, 0.01, 0.01, 0 ]
node_data = { "AA": CANON.AA_ground, "RAMA": CANON.RAMA_ground, "SS": CANON.SS_ground, "TM": Helix_emission }
data_dict = datamine('./tm_training.drct', code = 'null') # grab data_dict of the database. 

#distribution = pickle.load(open('./TMH_periods.pickle', 'rb'))
distribution = get_period(data_dict, 'TMSEQ', 'C')
#view_np(distribution,'TMH_periodicity', './', QUIET = False)
make_TM_C_region(CANON, node_data, distribution)

CANON.AA_EPS = 0.00001
CANON.SS_EPS = 0.00001
CANON.TM_EPS = 0.00000000001
CANON.RAMA_EPS = 0.00000001
CANON.write('./TMCTOP/HMMSTRTM.hmm')
CANON = hmm('HMMSTRTM.hmm', directory = './TMCTOP/')
CANON.read()
display_by_TM(CANON)

raise
