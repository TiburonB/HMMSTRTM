# TMHTOP.py
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

Helix_emission =  [ 0.1, 0.01, 0.01, .01, 0.01, 0.01, 0.83, 0.01, 0.01, 0 ]
node_data = { "AA": CANON.AA_ground, "RAMA": CANON.RAMA_ground, "SS": CANON.SS_ground, "TM": Helix_emission }
data_dict = datamine('./tm_training.drct', code = 'null') # grab data_dict of the database. 

distribution = pickle.load(open('./TMH_periods.pickle', 'rb'))
distribution = get_period(data_dict, 'TMSEQ', 'L')
#view_np(distribution,'TMH_periodicity', './', QUIET = False)
make_TM_L_region(CANON, node_data, distribution)

CANON.AA_EPS = 0.00001
CANON.SS_EPS = 0.00001
CANON.TM_EPS = 0.00000000001
CANON.RAMA_EPS = 0.00000001
CANON.write('./TMLTOP/HMMSTRTM.hmm')
CANON = hmm('HMMSTRTM.hmm', directory = './TMLTOP/')
CANON.read()
display(CANON)

raise
add_pyramid_cascade(CANON, node_data, 5, '_TMH_SCAP_',)
make_l_states(CANON,7,node_data,)
pyramid_trans = get_pyramid(3,"_test_", "SINK")
pyramid_states = [ l for l in list(pyramid_trans) ] 
for state in pyramid_states:
    node = {}
    for l in list(node_data):
        node[l] = [ k for k in node_data[l]]
    node['label'] = [state]
    CANON.add_node(node, prior = 0)
sink_node = {}
for l in list(node_data):
    sink_node[l] = [ k for k in node_data[l]]
sink_node['label'] = ['SINK']
CANON.add_node(sink_node, prior = 0)
for state in list(pyramid_trans):
    i = CANON.get_node_by_label(state)
    for t in pyramid_trans[state]:
        j = CANON.get_node_by_label(t[0])
        CANON.a[i][j] = t[1]
display_by_TM(CANON, QUIET = False)  

raise


Globular_emission =  [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
Helix_emission =  [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ]
Coil_emission =  [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ]
Inside_emission =  [ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 ]
Beta_emission =  [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ]
Loop_emission =  [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ]
IF_emission =  [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ]

# make TM-L
L_cap_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Loop_emission, "RAMA": CANON.RAMA_ground, "label": ["L_CAP_STATE_1"]}
L_m_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Loop_emission, "RAMA": CANON.RAMA_ground, "label": ["L_M_1_STATE"]}
make_lml_states(GLOB_1_leaving_states, 7, 7, 15, L_cap_dat, L_m_dat, GLOB_2_entering_states)
L_cap_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Loop_emission, "RAMA": CANON.RAMA_ground, "label": ["L_CAP_STATE_2"]}
L_m_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Loop_emission, "RAMA": CANON.RAMA_ground, "label": ["L_M_2_STATE"]}
make_lml_states(GLOB_2_leaving_states, 7, 7, 15, L_cap_dat, L_m_dat, GLOB_1_entering_states)
# make TM-B
B_cap_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Beta_emission, "RAMA": CANON.RAMA_ground, "label": ["B_CAP_STATE_1"]}
B_m_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Beta_emission, "RAMA": CANON.RAMA_ground, "label": ["B_M_1_STATE"]}
make_lml_states(GLOB_1_leaving_states, 7, 7, 15, B_cap_dat, B_m_dat, GLOB_2_entering_states)
B_cap_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Beta_emission, "RAMA": CANON.RAMA_ground, "label": ["B_CAP_STATE_2"]}
B_m_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Beta_emission, "RAMA": CANON.RAMA_ground, "label": ["B_M_2_STATE"]}
make_lml_states(GLOB_2_leaving_states, 7, 7, 15, B_cap_dat, B_m_dat, GLOB_1_entering_states)
# make TM-C
C_cap_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Coil_emission, "RAMA": CANON.RAMA_ground, "label": ["C_CAP_STATE_1"]}
C_m_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Coil_emission, "RAMA": CANON.RAMA_ground, "label": ["C_M_1_STATE"]}
make_lml_states(GLOB_1_leaving_states, 7, 7, 15, C_cap_dat, C_m_dat, GLOB_2_entering_states)
C_cap_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Coil_emission, "RAMA": CANON.RAMA_ground, "label": ["C_CAP_STATE_2"]}
C_m_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Coil_emission, "RAMA": CANON.RAMA_ground, "label": ["C_M_2_STATE"]}
make_lml_states(GLOB_2_leaving_states, 7, 7, 15, C_cap_dat, C_m_dat, GLOB_1_entering_states)
# make TM-I
I_cap_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Inside_emission, "RAMA": CANON.RAMA_ground, "label": ["I_CAP_STATE_1"]}
I_m_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Inside_emission, "RAMA": CANON.RAMA_ground, "label": ["I_M_1_STATE"]}
make_lml_states(GLOB_1_leaving_states, 7, 7, 15, I_cap_dat, I_m_dat, GLOB_2_entering_states)
I_cap_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Inside_emission, "RAMA": CANON.RAMA_ground, "label": ["I_CAP_STATE_2"]}
I_m_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": Inside_emission, "RAMA": CANON.RAMA_ground, "label": ["I_M_2_STATE"]}
make_lml_states(GLOB_2_leaving_states, 7, 7, 15, I_cap_dat, I_m_dat, GLOB_1_entering_states)
# make IF
IF_cap_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": IF_emission, "RAMA": CANON.RAMA_ground, "label": ["IF_CAP_STATE_1"]}
IF_m_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": IF_emission, "RAMA": CANON.RAMA_ground, "label": ["IF_M_1_STATE"]}
make_lml_states(GLOB_1_leaving_states, 4, 4, 20, IF_cap_dat, IF_m_dat, GLOB_1_entering_states)
IF_cap_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": IF_emission, "RAMA": CANON.RAMA_ground, "label": ["IF_CAP_STATE_2"]}
IF_m_dat = { "AA": CANON.AA_ground, "SS": CANON.SS_ground, "TM": IF_emission, "RAMA": CANON.RAMA_ground, "label": ["IF_M_2_STATE"]}
make_lml_states(GLOB_2_leaving_states, 4, 4, 20, IF_cap_dat, IF_m_dat, GLOB_2_entering_states)

CANON.write('./JANCANON/JANCANONTM.hmm')
raise


