#!/bin/bash

root -l rootlogon.C <<EOF
gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")
mpdloadlibs(kTRUE,kTRUE)
.L restore_dca.c+
EOF
