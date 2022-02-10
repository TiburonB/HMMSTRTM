#!/bin/csh
if ($#argv < 2) then
   echo "RUNPDB.csh < HMM file > < pdb code > < chain > "
   exit 1
endif

setenv MODELFILE $1
setenv PDBCODE $2
#echo PDBCODE
setenv PDBCHAIN $3
#echo PDBCHAIN
setenv THISDIR `pwd`
setenv GDIR $THISDIR
setenv DBDIR '~/Warehouse/Base/Externals/NR/nr'
setenv TEMPDIR "$GDIR/tmp/"

# run blast
#python3 blaster.py -if "${TEMPDIR}$PDBCODE$PDBCHAIN.fasta" 
if !(-e "${TEMPDIR}$PDBCODE$PDBCHAIN.blastp") then
        echo "RUNNING BLAST .... "
        blastp -query "${TEMPDIR}$PDBCODE$PDBCHAIN.fasta" -db  $DBDIR -out  "${TEMPDIR}$PDBCODE$PDBCHAIN.blastp" -outfmt 4 -evalue 0.001
endif
# parse blast
#if !(-e '${TEMPDIR}$PDBCODE$PDBCHAIN.fa') then
        echo "PARSING BLAST .... "
        echo "python3 blaster.py -if '${TEMPDIR}$PDBCODE$PDBCHAIN.blastp' -e parse"
        python3 blaster.py -if "${TEMPDIR}$PDBCODE$PDBCHAIN.blastp" -e parse
#endif
# make profile 
#if !(-e "${TEMPDIR}$PDBCODE$PDBCHAIN.profile") then
        echo "Making profile file... "
        echo "$GDIR/xmsa2profile '${TEMPDIR}$PDBCODE$PDBCHAIN.fa' '${TEMPDIR}$PDBCODE$PDBCHAIN.profile'"
        $GDIR/xmsa2profile "${TEMPDIR}$PDBCODE$PDBCHAIN.fa" "${TEMPDIR}$PDBCODE$PDBCHAIN.profile"
#endif

#if !(-e ${TEMPDIR}$PDBCODE$PDBCHAIN.drct) then
        echo "Making drct file ...."
        echo "$GDIR/xprofile2drct ${TEMPDIR}$PDBCODE$PDBCHAIN.profile $PDBCODE$PDBCHAIN ${TEMPDIR}$PDBCODE$PDBCHAIN.drct"
        $GDIR/xprofile2drct "${TEMPDIR}$PDBCODE$PDBCHAIN.profile $PDBCODE$PDBCHAIN ${TEMPDIR}$PDBCODE$PDBCHAIN.drct"       
#endif
# run profile on HMMSTR
echo "Scratching HMMSTR ..."
echo "$GDIR/xscratch $MODELFILE ${TEMPDIR}$PDBCODE$PDBCHAIN.drct" $PDBCODE
$GDIR/xscratch $MODELFILE "${TEMPDIR}$PDBCODE$PDBCHAIN.drct" $PDBCODE

