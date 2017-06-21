#!/bin/bash

INFILE=$1
OUTFILE=$2
RES_FIT_FILE=$3
DCAFILE=$4

#SBATCH --time=0:15:00
#SBATCH -D /tmp

PROJECT_DIR=/lustre/nyx/hades/user/parfenov/real-flow/

source /cvmfs/hades.gsi.de/install/5.34.34/hydra2-4.9n/defall.sh
#. /lustre/nyx/hades/user/parfenov/Soft/MPDRoot/build/config.sh

cd $PROJECT_DIR

OUTDIR=${OUTFILE%/*}
AFTERSLASHNAME=${INFILE##*/}
BASENAME=${AFTERSLASHNAME%.*}

root -l -b -q "main_flow.C(\"${INFILE}\",\"${OUTFILE}\",\"${RES_FIT_FILE}\",\"${DCAFILE}\")" 1>> ${OUTDIR}/${BASENAME}_flow_calc.OUT 2>> ${OUTDIR}/${BASENAME}_flow_calc.OUT
