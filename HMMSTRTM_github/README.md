To configure HMMSTR:
make all

To run a single pdb that exists in a direct access binary file (\*.drct):
xscratch HMMfile drctfile code

To run an amino acid sequence from the command line:
xinject HMMfile SEQUENCE <code> <chain>

If NCBI blast is installed, you can run a pdb with this command:
./RUNPDB.csh HMMfile code chain


To make Geeky Graphic:
python3 seq2jpg.py -sf './HMMSTR/HMMSTR_100<code>.seq' -od './HMMSTR/HMMSTR_100GAMMA/<code>dir/'
(manually run output scripts on terminal)
python3 gamma_gif.py -mf './HMMSTR/HMMSTR_100.hmm' -gf './HMMSTR/HMMSTR_100GAMMA/<code>'

# some dependencies for python files
pip3 install Bio
pip3 install scipy
pip3 install matplotlib
pip3 install numpy
pip3 install tqdm
pip3 install rmsd

