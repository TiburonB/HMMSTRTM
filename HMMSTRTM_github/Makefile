#| ======================== HMMSTRTM ====================================
#| HMMSTRTM is a model used to predict facets of protein local structure,
#| including amino acid profiles, secondary structure, ramachandran 
#| (dihedral angle) identifiers, and membrane localization.
#| This model is built from:
#| Bystroff, Christopher, Vesteinn Thorsson, and David Baker. "HMMSTR: 
#|     a hidden Markov model for local sequence-structure correlations 
#|     in proteins." Journal of molecular biology 301.1 (2000): 173-190.
#|
#|-------------------------------------------------------------------
#| To install HMMSTR programs, type:
#|
#| make clean
#| make all
#|
#| To test, use 
#| make test 
#|
#| To run:
#|  ./RUNHMMSTRTM.csh PDBFILE chaincode
#|
#| FORTRAN compiler:
#|  set FF to the fortran compiler on your system.
#|  For best results use the gfortran compiler 
#|
#| T.Benavides  20 JAN 2022
#|=======================================================================

FF = gfortran -O2
TESTPDB = 7bpo
TESTCHAIN = A

help :
	more INSTALL

bin32 :
	make all

all : xtrain xscratch xgetchain xinject xrenumber_one x3to1 xmsa2profile xmaster2drct xrama2drct clean  

# TRAIN A VERSION OF HMMSTRTM:
xtrain : hmm_io.f95 profile_io.f95 gamma_io.f95 hmm.f95 hmm_train.f95
	$(FF) -c -fdefault-double-8 hmm_io.f95 profile_io.f95 gamma_io.f95 hmm.f95 hmm_train.f95 -fcheck=all
	$(FF) -o xtrain hmm_io.o profile_io.o gamma_io.o hmm.o hmm_train.o
	# "xtrain HMM.file database.file epochs

# RUN a single pdb through HMMSTRTM (SCRATCH)
xscratch : hmm_io.f95 profile_io.f95 gamma_io.f95 hmm.f95 hmm_run.f95
	$(FF) -c -fdefault-double-8 hmm_io.f95 profile_io.f95 gamma_io.f95 hmm.f95 hmm_run.f95 -fcheck=all
	$(FF) -o xscratch hmm_io.o profile_io.o gamma_io.o hmm.o hmm_run.o
	# "xscratch  HMM.file database.file code

# extract a single chain from a pdb file
xgetchain : getchain.f95
	$(FF) -o xgetchain getchain.f95

xrenumber_one : renumber_one.f95
	$(FF) -o xrenumber_one renumber_one.f95

xmaster2drct : master2drct.f95
	$(FF) -o xmaster2drct master2drct.f95

xrama2drct : rama2drct.f95
	$(FF) -o xrama2drct rama2drct.f95

x3to1 : 3to1.f95
	$(FF) -o x3to1 3to1.f95

# RUN a single pdb through HMMSTRTM (inject)
xinject : profile_io.f95 hmm_inject.f95 
	$(FF) -c -fdefault-double-8 profile_io.f95 hmm_inject.f95
	$(FF) -o xinject profile_io.o hmm_inject.o
	# "xinject sequence"

xmakeprofile: profile_io.f95 make_profile.f95
	$(FF) -c -fdefault-double-8 profile_io.f95 make_profile.f95 -fcheck=all
	$(FF) -o xmakeprofile profile_io.o make_profile.o
	# "xmakeprofile database code"

xmsa2profile: msa2profile.f95
	$(FF) -o xmsa2profile msa2profile.f95

clean:
	rm *.mod
	rm *.o
