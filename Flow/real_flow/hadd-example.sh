#!/bin/bash

#SBATCH --time=0:15:00
#SBATCH -D /tmp

FILE_LIST=$1

. /lustre/nyx/hades/user/svintsov/MPDRoot/build/config.sh

OUTDIR=`dirname "$FILE_LIST"`

hadd ${FILE_LIST}.root `cat "$FILE_LIST"`
