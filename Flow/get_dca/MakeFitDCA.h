#ifndef MAKE_FIT_DCA_H

#define MAKE_FIT_DCA_H

#include <TMath.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "TROOT.h"
#include <TH1.h>

#include <iostream>

#include "../Utilities/utility.h"

void MakeFitDCA(TString inFileName, TString outFileName);

#endif