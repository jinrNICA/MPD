#!/bin/bash

#SBATCH -D /tmp
#SBATCH --time=0:15:00

INFILE=$1
OUTFILE=$2

LOG=${OUTFILE%.*}.log

cd /lustre/nyx/hades/user/parfenov/get-dca/

. /lustre/nyx/hades/user/parfenov/Soft/MPDRoot/build/config.sh

root -l rootlogon.C 1>>$LOG 2>>$LOG <<EOF
gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")
mpdloadlibs(kTRUE,kTRUE)
.L get_dca_cxx.so
get_dca("$INFILE","$OUTFILE")
EOF
