#!/bin/bash

FILE_LIST=$1
MOD=1

cd /lustre/nyx/hades/user/parfenov/restore-dca/

source /lustre/nyx/hades/user/parfenov/Soft/MPDRoot/build/config.sh

i=0
while read INFILE; do
  #~ i=$((i+1))
  #~ MOD=$((i%100))
  #~ if [ $MOD == 0 ]
  #~ then
    #~ sleep 30m
  #~ fi
  OUTFILE=/lustre/nyx/hades/user/parfenov/mpd_data/11gev/restored-dca-dst/${INFILE##*/}
  sbatch run.sh $INFILE $OUTFILE 
  #sleep 90
done <${FILE_LIST}
