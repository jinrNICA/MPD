#ifndef GET_FIT_H

#define GET_FIT_H

#include <TMath.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "TROOT.h"
#include <TH1.h>

#include <iostream>

#include "../Utilities/utility.h"

void get_fit(TString inFileName , TString outFileName);

#endif