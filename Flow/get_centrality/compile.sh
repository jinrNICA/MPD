#!/bin/bash

source /lustre/nyx/hades/user/parfenov/Soft/MPDRoot/build/config.sh

root.exe -l rootlogon.C <<EOF
gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")
mpdloadlibs(kTRUE,kTRUE)
.L get_multiplicity.cxx+
EOF
