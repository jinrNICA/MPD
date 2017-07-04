#ifndef MAKE_FIT_DCA_H

#define MAKE_FIT_DCA_H

//Including necessary classes & constants
#include "../Utilities/utility.h"


//Main function - can be executed via standard root
//
//Second (and final) iteration of the dca fit procedure
//
//USAGE:
//root -l 
//  .L MakeFitDCA.cxx+
//  MakeFitDCA("inFileName","outFileName")
//
//ARGUMENTS:
// inFileName - input: standard root file with TH1* histograms and TF1* fitted function needed for the further DCA cuts (output from get_fit(...))
// outFileName - output: standard root file with TH1* histograms and TF1* fitted function needed for the further DCA cuts (pt_sigma_fit) - for improving pt efficiency
//
void MakeFitDCA(TString inFileName, TString outFileName);

#endif