#!/bin/bash

#SBATCH -D /tmp
#SBATCH --time=0:25:00

INFILE=$1
OUTFILE=$2

LOG=${OUTFILE%.*}.log

cd /lustre/nyx/hades/user/svintsov/mpdroot_batch/restore-dca/

source /lustre/nyx/hades/user/svintsov/MPDRoot/build/config.sh

root -l rootlogon.C 1>>$LOG 2>>$LOG <<EOF
gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")
mpdloadlibs(kTRUE,kTRUE)
.L restore_dca_c.so
restore_dca("$INFILE","$OUTFILE")
EOF
