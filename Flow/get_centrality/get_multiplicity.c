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

void get_multiplicity(TString inFileName , TString outFileName , TString dcaFileName)
{

	const float pt_bins[] = {0.,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.};
	const int n_pt_bin = 12;

	const float eta_bins[] = {-1.5,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.5};
	const int n_eta_bin = 14;

	const Int_t   n_proj    = 3;
	const Int_t   n_hit_cut = 32;

	TH1F* h_multiplicity = new TH1F("h_multiplicity","multiplicity in TPC after cuts;N_{tracks} in TPC;N_{events};"           ,1600,0.,1600.);
	TH1F* h_mult_compare = new TH1F("h_mult_compare","Normalized multiplicity in TPC after cuts;N_{tracks} in TPC;N_{events};",1600,0.,1600.);
	TH1F* h_mult_old     = new TH1F("h_mult_old"    ,"multiplicity in TPC before cuts;N_{tracks} in TPC;N_{events};"          ,1600,0.,1600.);

	
	TFile *inFile = new TFile(inFileName.Data(),"READ");
	TTree *inTree = (TTree*) inFile->Get("cbmsim");
	
	TFile  *outFile = new TFile(outFileName.Data(),"RECREATE");

	//TFile  *dcaFile = new TFile("/lustre/nyx/hades/user/parfenov/dca_out_file_TDR.root","READ");
	TFile  *dcaFile = new TFile(dcaFileName.Data(),"READ");
	
	MpdEvent *MPDEvent=0;
	inTree->SetBranchAddress("MPDEvent.", &MPDEvent);
    TClonesArray *MpdGlobalTracks=0;
    TClonesArray *MCTracks=0;
	inTree->SetBranchAddress("MCTrack", &MCTracks);

	//TF1* f_dca[n_proj][n_pt_bin][n_eta_bin];
	TF1* f_pt_fit[n_proj][n_eta_bin];

	for (Int_t i_proj=0;i_proj<n_proj;i_proj++){
	    //for (Int_t i_pt=0;i_pt<n_pt_bin;i_pt++){
		for (Int_t i_eta=0;i_eta<n_eta_bin;i_eta++){
		    //f_dca[i_proj][i_pt][i_eta] = (TF1*) dcaFile->Get(Form("dca_fit[%i][%i][%i]",i_proj,i_pt,i_eta));
		    f_pt_fit[i_proj][i_eta] = (TF1*) dcaFile->Get(Form("sigma_pt_fit%i%i",i_proj,i_eta));
		}
	    //}
	}
  
    Int_t n_events = inTree->GetEntries();
    //Int_t n_events = 50;
    Int_t pt_bin;
    Int_t eta_bin;
    for (int event = 0; event < n_events; ++event)
    {
		cout << "EVENT N "<< event <<endl;
		inTree->GetEntry(event);
	    MpdGlobalTracks = MPDEvent->GetGlobalTracks();
	   
	    h_mult_old->Fill(MpdGlobalTracks->GetEntriesFast());

	    Long_t n_tracks_mpd=0;
	    for (Int_t track = 0;track<MpdGlobalTracks->GetEntriesFast();track++){
			MpdTrack* mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(track);
			
			
			pt_bin=-1;
			eta_bin=-1;
			if (mpdtrack->GetNofHits()<n_hit_cut) continue;

			for (Int_t i_pt=0;i_pt<n_pt_bin;i_pt++)
				if (TMath::Abs(mpdtrack->GetPt())>=pt_bins[i_pt] && TMath::Abs(mpdtrack->GetPt())<=pt_bins[i_pt+1]) pt_bin=i_pt;
			for (Int_t i_eta=0;i_eta<n_eta_bin;i_eta++)
				if (mpdtrack->GetEta()>eta_bins[i_eta] && mpdtrack->GetEta()<=eta_bins[i_eta+1]) eta_bin=i_eta;

			if (pt_bin==-1) continue;
			if (eta_bin==-1) continue;
			Double_t Pt = TMath::Abs(mpdtrack->GetPt());
			TF1 sigma_fit_X = *f_pt_fit[0][eta_bin];
			TF1 sigma_fit_Y = *f_pt_fit[1][eta_bin];
			TF1 sigma_fit_Z = *f_pt_fit[2][eta_bin];

			if (TMath::Abs(mpdtrack->GetDCAX()) >= sigma_fit_X(Pt)*2) continue;
			if (TMath::Abs(mpdtrack->GetDCAY()) >= sigma_fit_Y(Pt)*2) continue;
			if (TMath::Abs(mpdtrack->GetDCAZ()) >= sigma_fit_Z(Pt)*2) continue;
					

			if (TMath::Abs(mpdtrack->GetEta())>1.5) continue;
			if (TMath::Abs(mpdtrack->GetPt())<0 && TMath::Abs(mpdtrack->GetPt())>3) continue;
			
			n_tracks_mpd++;
	    }
		
	    if (n_tracks_mpd>0){ 
		h_multiplicity->Fill(n_tracks_mpd);
		h_mult_compare->Fill((Float_t) (n_tracks_mpd*1.5));

	    }
	}
	
	outFile->cd();
	h_multiplicity->Write();
	h_mult_compare->Write();
	h_mult_old    ->Write();
}
