#include <TMath.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "TROOT.h"
#include <TH1.h>

#include <iostream>

int GetPtBin(Float_t pt);
int GetEtaBin(Float_t eta);

const float ptBins[] = {0.,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.};
const int NptBins = 12;

const float etaBins[] = {-1.5,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.5};
const int NetaBins = 14;

const int Ndim = 3;

void get_fit(TString inFileName , TString outFileName)
{
	TH1F* h_dca[Ndim][NptBins][NetaBins];
	TF1 *dca_fit[Ndim][NptBins][NetaBins];
	
	char name[200];
	char title[200];
	
	TFile *inFile = new TFile(inFileName.Data(),"READ");
	TFile  *outFile = new TFile(outFileName.Data(),"RECREATE");
	
	for (int dim = 0; dim < Ndim; ++dim)
	{
		for (int ptbin = 0; ptbin < NptBins; ++ptbin)
		{
			for (int etabin = 0; etabin < NetaBins; ++etabin)
			{
				sprintf(name,"h_dca[%i][%i][%i]",dim,ptbin,etabin);
				h_dca[dim][ptbin][etabin] = (TH1F*)inFile->Get(name);
				
				sprintf(name,"dca_fit[%i][%i][%i]",dim,ptbin,etabin);
				dca_fit[dim][ptbin][etabin] = new TF1(name,"gaus",-0.2,0.2);
				h_dca[dim][ptbin][etabin]->Fit(name,"R");
				h_dca[dim][ptbin][etabin]->Write();
				dca_fit[dim][ptbin][etabin]->Write();
			}
		}
	}
}
