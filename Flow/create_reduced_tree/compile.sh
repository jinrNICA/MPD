#!/bin/bash

source /lustre/nyx/hades/user/parfenov/Soft/MPDRoot/build/config.sh

root -l rootlogon.C <<EOF
gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")
mpdloadlibs(kTRUE,kTRUE)
.L reducedTreeCreator.C+
EOF
