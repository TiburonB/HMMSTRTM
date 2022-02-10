#!/bin/csh
# RUNPDB
# This script takes as an argument a pdb code and chain
#  and tries to retrieve the pdb from rcsb, seperate the desired chain
#  and run blast, parse to profile and finally run HMMSTR on the profile file.
if ($#argv < 2) then
   echo "RUNPDB.csh < HMM file > < pdb code > < chain > < run blast > "
   exit 1
endif

setenv MODELFILE $1
setenv PDBCODE $2
#echo PDBCODE
setenv PDBCHAIN $3
setenv RUNBLAST $4

#echo PDBCHAIN
setenv THISDIR `pwd`
setenv GDIR $THISDIR
setenv DBDIR '~/Warehouse/Base/Externals/NR/nr'
setenv TEMPDIR "$GDIR/tmp/"
setenv PDBFILE "$TEMPDIR$PDBCODE.pdb"
setenv CIFFILE "$TEMPDIR$PDBCODE.cif"

# collect pdb file if not found
if !(-e $PDBFILE) then
   echo "$PDBFILE not found, collecting using rsync.bin"
   ./rsyncpdb.bin $PDBCODE
   cp "$PDBCODE.pdb" $PDBFILE
   rm "$PDBCODE.cif" 
endif

if !(-e $PDBFILE) then
   echo "rsyncpdb.bin failed, manually download pdb file into tmp directory."
   exit 1
endif

# get desired chain
if !(-e "$TEMPDIR$PDBCODE$PDBCHAIN.pdb") then
	( $GDIR/xgetchain $PDBCHAIN < $PDBFILE > "$TEMPDIR$PDBCODE$PDBCHAIN.pdb" \
		 ||  ( echo "Error in GETCHAIN" ; exit ) ) 
endif
# renumber chain
( $GDIR/xrenumber_one "${TEMPDIR}$PDBCODE$PDBCHAIN.pdb" "${TEMPDIR}$PDBCODE$PDBCHAIN.tmp" \
         ||  ( echo "Error in renmber" ; exit ))
mv "${TEMPDIR}$PDBCODE$PDBCHAIN.tmp" "${TEMPDIR}$PDBCODE$PDBCHAIN.pdb"

# get .seq file
if !(-e "${TEMPDIR}$PDBCODE$PDBCHAIN.fasta") then 
	( $GDIR/x3to1 "+" < "${TEMPDIR}$PDBCODE$PDBCHAIN.pdb" > "${TEMPDIR}$PDBCODE$PDBCHAIN.fasta" \
		 ||  ( echo "Error in 3to1" ; exit ))
endif

if ( $RUNBLAST ) then
	# run blast
	#python3 blaster.py -if "${TEMPDIR}$PDBCODE$PDBCHAIN.fasta" 
	if !(-e "${TEMPDIR}$PDBCODE$PDBCHAIN.blastp") then
		echo "RUNNING BLAST .... " 
		blastp -query "${TEMPDIR}$PDBCODE$PDBCHAIN.fasta" -db  $DBDIR -out  "${TEMPDIR}$PDBCODE$PDBCHAIN.blastp" -outfmt 4 -evalue 0.001
	endif
	# parse blast
#	if !(-e '${TEMPDIR}$PDBCODE$PDBCHAIN.fa') then
		echo "PARSING BLAST .... "
		echo "python3 blaster.py -if '${TEMPDIR}$PDBCODE$PDBCHAIN.blastp' -e parse"
		python3 blaster.py -if "${TEMPDIR}$PDBCODE$PDBCHAIN.blastp" -e parse
#	endif
else
     # copy .fasta to .fa
     echo "copying .fasta to .fa"
     cp $TEMPDIR$PDBCODE$PDBCHAIN.fasta $TEMPDIR$PDBCODE$PDBCHAIN.fa
endif

# make profile 
#if !(-e "${TEMPDIR}$PDBCODE$PDBCHAIN.profile") then
	echo "Making profile file... "
	echo "$GDIR/xmsa2profile '${TEMPDIR}$PDBCODE$PDBCHAIN.fa' '${TEMPDIR}$PDBCODE$PDBCHAIN.profile'"
	$GDIR/xmsa2profile "${TEMPDIR}$PDBCODE$PDBCHAIN.fa" "${TEMPDIR}$PDBCODE$PDBCHAIN.profile"
#endif
#if !(-e ${TEMPDIR}$PDBCODE$PDBCHAIN.drct) then
	echo "Making drct file ...."
	echo "python3 makedrct.py -p '$TEMPDIR$PDBCODE.pdb' -f '${TEMPDIR}$PDBCODE$PDBCHAIN.profile' -c $PDBCHAIN"
	python3 makedrct.py -p "$TEMPDIR$PDBCODE.pdb" -f "${TEMPDIR}$PDBCODE$PDBCHAIN.profile" -c $PDBCHAIN
#endif
# run profile on HMMSTR
echo "Scratching HMMSTR ..."
echo "$GDIR/xscratch $MODELFILE ${TEMPDIR}$PDBCODE$PDBCHAIN.drct" $PDBCODE
$GDIR/xscratch $MODELFILE "${TEMPDIR}$PDBCODE$PDBCHAIN.drct" $PDBCODE





