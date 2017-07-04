#ifndef GET_MULTIPLICITY_H

#define GET_MULTIPLICITY_H

#include <iostream>

//Including necessary classes & constants
#include "../Utilities/utility.h"


//Main function
//
//USAGE:
//root -l rootlogon.C
//  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")
//  mpdloadlibs(kTRUE,kTRUE)
//  .L get_multiplicity.cxx+
//  get_multiplicity("inFileName","outFileName","dcaFileName")
//
//ARGUMENTS:
// inFileName - input cbmsim tree with corrected dca values (see restore_dca.h)
// outFileName - output: standard root file with TH1* histograms needed for the further DCA cuts
// dcaFileName - Second iteration of the dca fitting containing sigma_pt_fit TF1* functions (output MakeFitDCA(...))
//
void get_multiplicity(TString inFileName , TString outFileName , TString dcaFileName);

#endif