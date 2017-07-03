#!/bin/bash

INFILE=$1
OUTFILE=$2
DCAFILE=$3

#SBATCH --time=0:15:00
#SBATCH -D /tmp

PROJECT_DIR=$PWD #/lustre/nyx/hades/user/parfenov/real-flow/

#. /lustre/nyx/hades/user/parfenov/Soft/MPDRoot/build/config.sh
source /cvmfs/hades.gsi.de/install/5.34.34/hydra2-4.9n/defall.sh

cd $PROJECT_DIR

OUTDIR=${OUTFILE%/*}
AFTERSLASHNAME=${INFILE##*/}
BASENAME=${AFTERSLASHNAME%.*}

root -l -b -q "main_resolutions.C(\"${INFILE}\",\"${OUTFILE}\",\"${DCAFILE}\")" 1>> ${OUTDIR}/${BASENAME}_rec_calc.OUT 2>> ${OUTDIR}/${BASENAME}_rec_calc.ERR


mv BASENAME_res_out.root OUTDIR/.
