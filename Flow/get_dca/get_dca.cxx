#include "get_dca.h"

void get_dca(TString inFileName , TString outFileName)
{	

	BinningData* bins = new BinningData;
	FormKinematicBins(bins);

	const int NptBins = bins->GetPtBinSize();
	const int NetaBins = bins->GetEtaBinSize();

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
				sprintf(title,"DCA distribution for %.2f < p_{t} < %.2f and %.2f < #eta < %.2f", bins->GetPtBinContent(ptbin), bins->GetPtBinContent(ptbin+1),bins->GetEtaBinContent(etabin),bins->GetEtaBinContent(etabin+1));
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
			Int_t pt_bin = bins->GetPtBin(TMath::Abs(mpdtrack->GetPt()));
			Int_t eta_bin = bins->GetEtaBin(mpdtrack->GetEta());
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
