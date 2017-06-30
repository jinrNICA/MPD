#ifndef UTILITY_H

#define UTILITY_H

#define _MAX_TRACKS 200000
#define _N_ARM 2 // 3 IS FULL DETECTOR
#define _N_HARM 2
#define _N_SORTS 4
#define _N_MODULES_TOTAL 90
#define _N_METHOD 2 // 0 - TPC , 1 - ZDC
#define _N_QCOMP 2

#include <TMath.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TChain.h>
#include <TFile.h>
#include <TStyle.h>
#include <TApplication.h>
#include "SpecFuncMathMore.h"
#include "BinningData.cxx"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>

using TMath::Abs;
using TMath::Cos;
using TMath::Sin;
using TMath::ATan2;
using TMath::Sqrt;
using TMath::LocMin;
using TMath::Pi;

using std::cout;
using std::endl;

const Float_t Cut_Pt_Min = 0.;
const Float_t Cut_Eta_Min = 0.7;
const Float_t Cut_Eta_Max = 1.5;
const Int_t Cut_No_Of_hits_min = 32;

const int Ndim = 3;

const float centralityBinsFlow[] = {10.,20.,40.,50};
const int NcentralityBinsFlow = 3;

const float centralityBinsRes[] = {0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100};
//const float centralityBinsRes[] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,46.,47.,48.,49.,50.,51.,52.,53.,54.,55.,56.,57.,58.,59.,60.,61.,62.,63.,64.,65.,66.,67.,68.,69.,70.,71.,72.,73.,74.,75.,76.,77.,78.,79.,80.,81.,82.,83.,84.,85.,86.,87.,88.,89.,90.,91.,92.,93.,94.,95.,96.,97.,98.,99.,100.};
const int NcentralityBinsRes = 20;

const float ptBins[] = {0.,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.};
const int NptBins = 12;

const float etaBins[] = {-1.5,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.5};
const int NetaBins = 14;

const float rapidityBins[] = {-2., -1.8, -1.6, -1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.};
const int NrapidityBins = 20;

void FormKinematicBins(BinningData* bins)
{
    bins->SetPtBins(ptBins,NptBins);
    bins->SetEtaBins(etaBins,NetaBins);
    bins->SetRapidityBins(rapidityBins,NrapidityBins);
}

const TString arm_names[_N_ARM] = {TString("R"),TString("L")};
const TString sorts_of_particles[4] = {TString("all sorts") , TString("pions (211)") , TString("protons (2212)") , TString("kaons (321)")};
const TString methods_names[_N_METHOD] = {TString("TPC") , TString("FHCal")};

class FlowParticle
{
public:

    double Eta;
    double Pt;
    double Phi;
    double Rapidity;

    FlowParticle();

    FlowParticle(double Eta, double Pt, double Phi, double Rapidity);
};

class EPParticle
{
public:

    double Eta;
    double Pt;
    double Phi;

    EPParticle();

    EPParticle(double Eta, double Pt, double Phi);
};

Double_t ResEventPlane(Double_t chi, Int_t harm); //harm = 1 or 2 for our case

Double_t Chi(Double_t res, Int_t harm); //harm = 1 or 2 for our case

Double_t* GetAngles();

Float_t Unfold(Float_t phiEP_mc, Float_t psi_N_FULL, Int_t harm);

#endif
