#ifndef GET_RES_H

#define GET_RES_H

#include <cstdio>

#include <TProfile.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TF1.h>

const int _N_HARM = 2;
const int _N_METHOD = 2;

const float centralityBinsFlow[] = {10.,20.,40.,50};
const int NcentralityBinsFlow = 3;

const float centralityBinsRes[] = {0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100};
const int NcentralityBinsRes = 20;

const float ptBins[] = {0.,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.};
const int NptBins = 12;

const float etaBins[] = {-1.5,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.5};
const int NetaBins = 14;

const float rapidityBins[] = {-2., -1.8, -1.6, -1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.};
const int NrapidityBins = 20;

const TString methods_names[_N_METHOD] = {TString("TPC") , TString("FHCal")};

using TMath::Sqrt;

void get_res(TString inFileName , TString outFileName);

#endif