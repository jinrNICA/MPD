#ifndef GET_FIT_H

#define GET_FIT_H

#include <iostream>

//Including necessary classes & constants
#include "../Utilities/utility.h"


//Main function - can be executed via standard root
//
//First iteration of the dca fit procedure
//
//USAGE:
//root -l 
//  .L get_fit.cxx+
//  get_fit("inFileName","outFileName")
//
//ARGUMENTS:
// inFileName - input: standard root file with TH1* histograms needed for the further DCA cuts (output file from get_dca(...))
// outFileName - output: standard root file with TH1* histograms and TF1* fitted function needed for the further DCA cuts
//
void get_fit(TString inFileName , TString outFileName);

#endif