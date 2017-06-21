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

#include <iostream>

int GetPtBin(Float_t pt);
int GetEtaBin(Float_t eta);

const float ptBins[] = {0.,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.};
const int NptBins = 12;

const float etaBins[] = {-1.5,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.5};
const int NetaBins = 14;

const int Ndim = 3;

void get_dca(TString inFileName , TString outFileName)
{	
	TH1F* h_dca[Ndim][NptBins][NetaBins];
	
	char name[200];
	char title[200];
	
	for (int dim = 0; dim < Ndim; ++dim)
	{
		for (int ptbin = 0; ptbin < NptBins; ++ptbin)
		{
			for (int etabin = 0; etabin < NetaBins; ++etabin)
			{
				sprintf(name,"h_dca[%i][%i][%i]",dim,ptbin,etabin);
				sprintf(title,"DCA distribution for %.2f < p_{t} < %.2f and %.2f < #eta < %.2f", ptBins[ptbin], ptBins[ptbin+1],etaBins[etabin],etaBins[etabin+1]);
				h_dca[dim][ptbin][etabin] = new TH1F(name,title,4000,-50.,50);
			}
		}
	}
	
	TFile *inFile = new TFile(inFileName.Data(),"READ");
	TTree *inTree = (TTree*) inFile->Get("cbmsim");
	
	TFile  *outFile = new TFile(outFileName.Data(),"RECREATE");
	
	MpdEvent *MPDEvent=0;
	inTree->SetBranchAddress("MPDEvent.", &MPDEvent);
    TClonesArray *MpdGlobalTracks=0;
    
    Int_t n_events = inTree->GetEntries();
    //Int_t n_events = 50;
    for (int event = 0; event < n_events; ++event)
    {
		cout << "EVENT N "<< event <<endl;
		inTree->GetEntry(event);
	    MpdGlobalTracks = MPDEvent->GetGlobalTracks();
	    for (int track = 0; track < MpdGlobalTracks->GetEntriesFast(); ++track)
	    {
			MpdTrack *mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(track); 
			Int_t pt_bin = GetPtBin(TMath::Abs(mpdtrack->GetPt()));
			Int_t eta_bin = GetEtaBin(mpdtrack->GetEta());
			if ( (eta_bin == -1) || (pt_bin == -1 )) continue;
			h_dca[0][pt_bin][eta_bin]->Fill(mpdtrack->GetDCAX());
			h_dca[1][pt_bin][eta_bin]->Fill(mpdtrack->GetDCAY());
			h_dca[2][pt_bin][eta_bin]->Fill(mpdtrack->GetDCAZ());
		} 
	}
	
	outFile->cd();
	for (int dim = 0; dim < Ndim; ++dim)
	{
		for (int ptbin = 0; ptbin < NptBins; ++ptbin)
		{
			for (int etabin = 0; etabin < NetaBins; ++etabin)
			{
				h_dca[dim][ptbin][etabin]->Write();
			}
		}
	}
}

int GetPtBin(Float_t pt)
{
	int pt_bin = -1;
	for (int i = 0; i < NptBins; ++i)
	if ((pt > ptBins[i]) && (pt <= ptBins[i + 1])) 
		pt_bin = i;
	return pt_bin;
}

int GetEtaBin(Float_t eta)
{
	int eta_bin = -1;
	for (int i = 0; i < NetaBins; ++i)
	if ((eta > etaBins[i]) && (eta <= etaBins[i + 1])) 
		eta_bin = i;
	return eta_bin;
}
