# HMM2TOP.py
#  1/17/2022. Culminative.
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
GLOB_2 = hmm('HMMSTR_150.hmm', directory = './HMMSTR/')
GLOB_2.read()

CANON = hmm('HMMSTRTM.hmm', directory = './HMMSTRTMTOP/')
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

# glob_1 node_data
for i in range(0, GLOB_1.n):
    GLOB_1.node_data[i]["label"] = [ "GLOBULAR_SIDE_1" ]
    for j in range(0, len(GLOB_1.node_data[i]["TM"])):
        if j == 0: # replace TM emission 1 w/ 100% prob
            GLOB_1.node_data[i]["TM"][j] = 1.0
        else:
            GLOB_1.node_data[i]["TM"][j] = 0.0
for i in range(0, GLOB_2.n):
    GLOB_1.node_data[i]["label"] = [ "GLOBULAR_SIDE_2" ]
    for j in range(0, len(GLOB_2.node_data[i]["TM"])):
        if j == 0: # replace TM emission 2 w/ 100% prob
            GLOB_2.node_data[i]["TM"][j] = 1.0
        else:
            GLOB_2.node_data[i]["TM"][j] = 0.0

# glob_2 node data
for i,node in enumerate(GLOB_1.node_data):
    if i == 0:
        CANON.add_naught_node(GLOB_1.priors[i], node["AA"], node["SS"], node["TM"], node["RAMA"])
        CANON.node_data[0]["label"] = ["GLOBULAR_SIDE_1_NAUGHT"]
    else:
        CANON.add_node(node, prior = GLOB_1.priors[i] )
        CANON.node_data[CANON.n-1]["label"] = ["GLOBULAR_SIDE_1"]
for i,node in enumerate(GLOB_2.node_data):
    if i == 0:
        CANON.add_naught_node(0, node["AA"], node["SS"], node["TM"], node["RAMA"])
        CANON.node_data[1]["label"] = ["GLOBULAR_SIDE_2_NAUGHT"]
    else:
        CANON.add_node(node, prior = 0)
        CANON.node_data[CANON.n-1]["label"] = ["GLOBULAR_SIDE_2"]

# add glob_1 transitions
for i in range(0, GLOB_1.n):
    ti = i
    if i > 0:
        ti = i + 1
    for j in range(0, GLOB_1.n):
        tj = j
        if j > 0:
            tj = j + 1
        t = GLOB_1.a[i][j]
        if t > 0:
            CANON.a[ti][tj] = t

# add glbo_2 node_data
for i in range(0, GLOB_2.n):
    if i > 0: 
        ti = i + GLOB_1.n
    else:
        ti = 1
    for j in range(0, GLOB_1.n):
        if j > 0:
            tj = j + GLOB_1.n
        else:
            tj = 1
        t = GLOB_1.a[i][j]
        if t > 0:
            CANON.a[ti][tj] = t

#Helix_emission =  [ 0.1, 0.01, 0.01, .83, 0.01, 0.01, 0.01, 0.01, 0.01, 0 ]
#node_data = { "AA": CANON.AA_ground, "RAMA": CANON.RAMA_ground, "SS": CANON.SS_ground, "TM": Helix_emission }
#data_dict = datamine('./tm_training.drct', code = 'null') # grab data_dict of the database. 

# rerun the other TMTOP python scripts
import os
os.system('python3 TMHTOP.py')
os.system('python3 TMITOP.py')
os.system('python3 TMLTOP.py')
os.system('python3 TMFTOP.py')
os.system('python3 TMCTOP.py')
os.system('python3 TMBTOP.py')

#copy_labels('./TMHTOP/HMMSTRTM.hmm', './TMHTOP/HMMSTRTM_ 56.hmm')
create_TM_regions(CANON, "H", './TMHTOP/HMMSTRTM.hmm', GLOB_2.n-1)
#copy_labels('./TMITOP/HMMSTRTM.hmm', './TMITOP/HMMSTRTM_ 92.hmm')
if True:
    create_TM_regions(CANON, "I", './TMITOP/HMMSTRTM.hmm', GLOB_2.n-1)
    #copy_labels('./TMLTOP/HMMSTRTM.hmm', './TMLTOP/HMMSTRTM_ 93.hmm')
    create_TM_regions(CANON, "L", './TMLTOP/HMMSTRTM.hmm', GLOB_2.n-1)
    #copy_labels('./TMFTOP/HMMSTRTM.hmm', './TMFTOP/HMMSTRTM_ 89.hmm')
    create_TM_regions(CANON, "F", './TMFTOP/HMMSTRTM.hmm', GLOB_2.n-1)
    #copy_labels('./TMCTOP/HMMSTRTM.hmm', './TMCTOP/HMMSTRTM_ 91.hmm')
    create_TM_regions(CANON, "C", './TMCTOP/HMMSTRTM.hmm', GLOB_2.n-1)
    #copy_labels('./TMBTOP/HMMSTRTM.hmm', './TMBTOP/HMMSTRTM_ 99.hmm')
    create_TM_regions(CANON, "B", './TMBTOP/HMMSTRTM.hmm', GLOB_2.n-1)

    # add the following transitions:
    # 1. fcoresink -> hcoresource
    CANON.a[CANON.get_node_by_label('s_l_TMFCORE_4_1')][CANON.get_node_by_label('s_l_TMHCORE_3_1')] = 0.5
    CANON.a[CANON.get_node_by_label('s_l_TMFCORE_4_2')][CANON.get_node_by_label('s_l_TMHCORE_3_2')] = 0.5
    # 2. bcoresink -> icoresource
    CANON.a[CANON.get_node_by_label('s_E_TMB_0_1')][CANON.get_node_by_label('s_l_TMICORE_3_1')] = 0.5
    CANON.a[CANON.get_node_by_label('s_E_TMB_0_2')][CANON.get_node_by_label('s_l_TMICORE_3_2')] = 0.5
    # 3. icoresink -> bcoresource2
    CANON.a[CANON.get_node_by_label('s_l_TMICORE_4_1')][CANON.get_node_by_label('s_S_TMB_0_2')] = 0.5
    CANON.a[CANON.get_node_by_label('s_l_TMICORE_4_2')][CANON.get_node_by_label('s_S_TMB_0_1')] = 0.5
    # 4. hcoresink -> fcoresource
    CANON.a[CANON.get_node_by_label('s_l_TMHCORE_4_1')][CANON.get_node_by_label('s_l_TMFCORE_3_2')] = 0.5
    CANON.a[CANON.get_node_by_label('s_l_TMHCORE_4_2')][CANON.get_node_by_label('s_l_TMFCORE_3_1')] = 0.5

print(CANON.Bties)

CANON.AA_EPS = 0.00001
CANON.SS_EPS = 0.00001
CANON.TM_EPS = 0.00000000001
CANON.RAMA_EPS = 0.00000001
CANON.write('./HMMSTRTM2/HMMSTRTM.hmm')
CANON = hmm('HMMSTRTM.hmm', directory = './HMMSTRTM2/')
CANON.read()
display(CANON)

raise


