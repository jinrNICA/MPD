#!/bin/bash

FILE_LIST=$1
MOD=1

cd /lustre/nyx/hades/user/parfenov/get-dca/

i=0
while read INFILE; do
  #i=$((i+1))
  #MOD=$((i%100))
  #if [ $MOD == 0 ]
  #then
  #  sleep 10m
  #fi
  OUTFILE=/lustre/nyx/hades/user/parfenov/get-dca/TMP/${INFILE##*/}_dca.root
  sbatch run_dca.sh $INFILE $OUTFILE
  #sleep 1
done <${FILE_LIST}
