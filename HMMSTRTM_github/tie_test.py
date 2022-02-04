# tie_test.py TLB 1/11/22

from hmm import *
hmm_file = './JAN3/HMMSTRTM.hmm'
hmm_file_2 = './JAN4/HMMSTRTM.hmm'
hmm_file_e1 = './JAN3/HMMSTRTM_  2.hmm'
out_file = './tie_test.hmm'

# read an hmm with ties and try to print out the ds, 
#   check that ds is properly printed on file write.
HMM = hmm(hmm_file)
HMM.read()
print(HMM.Bties)
if False:
    HMM.write(out_file)
    new_HMM = hmm(out_file)
    new_HMM.read()
    print(new_HMM.Bties)

# Check Fortran i/o is working properly.
if False:
    from hmm_tools import train
    train(HMM, './just_tm.drct', 1)
    trained_HMM = hmm(hmm_file_e1)
    trained_HMM.read()
    print(trained_HMM.Bties)


# Check that we can copy ties from one file to another
if True:
    from hmm_tools import clear_ties, copy_ties
    clear_ties(hmm_file_2) # reset all ties
    copy_ties(hmm_file, hmm_file_2)
    HMM2 = hmm(hmm_file_2)
    HMM2.read()
    print(HMM2.Bties)

