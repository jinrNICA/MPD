#!/bin/bash

#SBATCH -D /tmp
#SBATCH --time=0:15:00

INFILEHIST=$1
INFILE=$2
OUTFILE=$3
DCAFILE=$4

LOG=${OUTFILE%.*}.log

cd /lustre/nyx/hades/user/parfenov/check-dca/

. /lustre/nyx/hades/user/parfenov/Soft/MPDRoot/build/config.sh

root -l rootlogon.C 1>>$LOG 2>>$LOG <<EOF
gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")
mpdloadlibs(kTRUE,kTRUE)
.L reducedTreeCreator_C.so
reducedTreeCreator rtc = reducedTreeCreator("$INFILEHIST", "$INFILE", "$OUTFILE", "$DCAFILE")
rtc.CreateReducedTree()
EOF
