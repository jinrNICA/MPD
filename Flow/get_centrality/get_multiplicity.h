#ifndef GET_MULTIPLICITY_H

#define GET_MULTIPLICITY_H

#include <TMath.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "TROOT.h"
#include <FairMCEventHeader.h>
#include <MpdEvent.h>
#include <TClonesArray.h>
#include <MpdTrack.h>
#include <FairMCTrack.h>
#include <TH1.h>
#include <TF1.h>

#include <iostream>

#include "../Utilities/utility.h"

void get_multiplicity(TString inFileName , TString outFileName , TString dcaFileName);

#endif