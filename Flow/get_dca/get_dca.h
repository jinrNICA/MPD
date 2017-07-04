#ifndef GET_DCA_H

#define GET_DCA_H

//Including necessary classes & constants
#include "../Utilities/utility.h"

#include <iostream>

//Main function
//
//USAGE:
//root -l rootlogon.C
//  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")
//  mpdloadlibs(kTRUE,kTRUE)
//  .L get_dca.cxx+
//  get_dca("inFileName","outFileName")
//
//ARGUMENTS:
// inFileName - input cbmsim tree with corrected dca values (see restore_dca.h)
// outFileName - output: standard root file with TH1* histograms needed for the further DCA cuts
//
void get_dca(TString inFileName , TString outFileName);

#endif