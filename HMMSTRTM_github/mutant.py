# mutant.py 1/7/2022 warm-up script.
# reset globular TM emissions to equal prob of all characters, retrain.

from hmm import * 
from hmm_tools import * 
drct = './tm_training.drct'

if False:
    HMMSTR = hmm('HMMSTR_150.hmm', directory = './HMMSTR/')
    HMMSTR.read()
    for i in range(0, HMMSTR.n):
        temp_tm = [ p for p in HMMSTR.node_data[i]["TM"]]
        for j in range(0, len(temp_tm)):
            temp_tm[j] = 1 # set all TM-emissions equal.
        HMMSTR.node_data[i]["TM"] = [ p for p in temp_tm ]

    HMMSTR.write('./MUTTS/mutant.hmm')

    print("MUHUHUAUAHAUHAHUA")
    mutt = hmm('mutant.hmm' , directory = './MUTTS/')
    mutt.read()

    train(mutt, drct, 10) # approximately 90 minute training time.
if True:
    HMMSTR_mutt_10 = hmm('mutant_ 10.hmm', directory = './MUTTS/')
    HMMSTR_mutt_10.read()
    examine(HMMSTR_mutt_10, drct)
    scratch(HMMSTR_mutt_10, drct, '6xyt')
    scratch(HMMSTR_mutt_10, drct, '7cal')
