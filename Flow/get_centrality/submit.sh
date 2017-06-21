#!/bin/bash

FILE_LIST=$1
DCA_FILE=$2
MOD=1

cd /lustre/nyx/hades/user/parfenov/get-multiplicity/

i=0
while read INFILE; do
  #i=$((i+1))
  #MOD=$((i%100))
  #if [ $MOD == 0 ]
  #then
  #  sleep 10m
  #fi
  OUTFILE=/lustre/nyx/hades/user/parfenov/get-multiplicity/TMP/${INFILE##*/}_multiplicity.root
  sbatch run_mult.sh $INFILE $OUTFILE $DCA_FILE
  #sleep 1
done <${FILE_LIST}
