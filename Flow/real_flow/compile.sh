#!/bin/bash

source /cvmfs/hades.gsi.de/install/5.34.34/hydra2-4.9n/defall.sh
. /lustre/nyx/hades/user/parfenov/Soft/MPDRoot/build/config.sh

root -l <<EOF
gSystem->Load("libMathMore")
.L ../Utilities/BinningData.cxx+
.L ../Utilities/utility.cxx+
.L MpdCalculator.cxx+
EOF
