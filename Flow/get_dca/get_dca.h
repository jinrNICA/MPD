#ifndef GET_DCA_H

#define GET_DCA_H

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

#include "../Utilities/utility.h"

#include <iostream>

void get_dca(TString inFileName , TString outFileName);

#endif