#ifndef GET_RES_H

#define GET_RES_H

#include <cstdio>
//Including necessary classes & constants
#include "../Utilities/utility.h"
#include "./MpdCalculator.cxx"

//Main function - can be executed via standard root
//
//Calculating adn fitting resolution correction factor
//
//USAGE:
//root -l 
//  gSystem->Load("libMathMore")
//  .L get_res.cxx+
//  get_res("inFileName","outFileName")
//
//ARGUMENTS:
// inFileName - input: standard root file with TProfiles and histograms containing data for resolution correction factor (output from MpdCalculator::CalculateResolutions(0))
// outFileName - output: standard root file with TH1* histograms and TF1* fitted function needed for the azimuthal flow calculation (MpdCalculator::CalculateFlow(...))
//
void get_res(TString inFileName , TString outFileName);

#endif