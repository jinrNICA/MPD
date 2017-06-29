#!/bin/bash

source /cvmfs/hades.gsi.de/install/5.34.34/hydra2-4.9n/defall.sh
#. /lustre/nyx/hades/user/parfenov/Soft/MPDRoot/build/config.sh

root -l <<EOF
gSystem->Load("libMathMore")
.L ../Utilities/utility.C+
.L MpdCalculator.C+
EOF
