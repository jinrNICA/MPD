#include "get_fit.h"

void get_fit(TString inFileName , TString outFileName)
{
	BinningData* bins = new BinningData;
	FormKinematicBins(bins);

	const int NptBins = bins->GetPtBinSize();
	const int NetaBins = bins->GetEtaBinSize();

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
