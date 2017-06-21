#define _MAX_TRACKS 5000
#define _N_ARM 2 // 3 IS FULL DETECTOR
#define _N_HARM 2
#define _N_SORTS 4
#define _N_MODULES_TOTAL 90
#define _N_METHOD 2 // 0 - TPC , 1 - ZDC
#define _N_QCOMP 2

#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TChain.h>
#include <TFile.h>

void fill_multiplicity(TString inFileName , TString outFileName)
{
	TH1F* h_multiplicity_before = new TH1F("h_multiplicity_before","multiplicity in TPCbefore cuts;N_{tracks} in TPC;N_{events};",1600,0.,1600.);
	TProfile* p_b_vs_multiplicity = new TProfile("p_b_vs_multiplicity","Impact parameter vs TPC multiplicity",800,0,1600);
	TH1F* centrality = new TH1F("h_centrality","h_centrality",100,0,100);
	
	inChain = new TChain("cbmsim_reduced");
	inChain->Add(inFileName.Data());

	outFile = new TFile(outFileName.Data(), "RECREATE");
	
	Float_t b_mc;
	Long64_t n_tracks_mc;
	Long64_t n_tracks_mpd;
    Float_t signed_pt_mpd[_MAX_TRACKS];
    Float_t pt_mc[_MAX_TRACKS];
	Float_t eta_mpd[_MAX_TRACKS];
	Int_t n_hits_mpd[_MAX_TRACKS];
	Float_t phi_mpd[_MAX_TRACKS];
	Float_t theta_mpd[_MAX_TRACKS];
	Float_t tof_beta_mpd[_MAX_TRACKS];
    Long64_t id_from_mc_mpd[_MAX_TRACKS];
    Int_t mother_ID_mc[_MAX_TRACKS];
    Int_t PDG_code_mc[_MAX_TRACKS];
    Float_t phiEP_mc;
    Float_t ZDC_energy_mpd[_N_MODULES_TOTAL];
    Float_t chi2_mpd[_MAX_TRACKS];
    Float_t eta_mc[_MAX_TRACKS];
    Float_t px_mc[_MAX_TRACKS];
    Float_t py_mc[_MAX_TRACKS];
    Float_t pz_mc[_MAX_TRACKS];
    Float_t energy_mc[_MAX_TRACKS];
    Float_t q_vectors_mpd[_N_QCOMP][_N_HARM][_N_METHOD][_N_ARM];
    Float_t q_norm_mpd[_N_HARM][_N_METHOD][_N_ARM];
    Float_t eta_match_mc_mpd[_MAX_TRACKS];
    Float_t phi_match_mc_mpd[_MAX_TRACKS];
    Float_t signed_pt_match_mc_mpd[_MAX_TRACKS];
    Int_t centrality_tpc_mpd;	

	inChain->SetBranchAddress("b_mc",&b_mc);
	inChain->SetBranchAddress("n_tracks_mc",&n_tracks_mc);
	inChain->SetBranchAddress("n_tracks_mpd",&n_tracks_mpd);
	inChain->SetBranchAddress("signed_pt_mpd",signed_pt_mpd);
	inChain->SetBranchAddress("eta_mpd",eta_mpd);
	inChain->SetBranchAddress("n_hits_mpd",n_hits_mpd);
	inChain->SetBranchAddress("phi_mpd",phi_mpd);
	inChain->SetBranchAddress("theta_mpd",theta_mpd);
	inChain->SetBranchAddress("tof_beta_mpd",tof_beta_mpd);
	inChain->SetBranchAddress("id_from_mc_mpd",id_from_mc_mpd);
	inChain->SetBranchAddress("mother_ID_mc",mother_ID_mc);
	inChain->SetBranchAddress("PDG_code_mc",PDG_code_mc);
	inChain->SetBranchAddress("phiEP_mc",&phiEP_mc);
	inChain->SetBranchAddress("ZDC_energy_mpd",ZDC_energy_mpd);
	inChain->SetBranchAddress("centrality_tpc_mpd",&centrality_tpc_mpd);

    inChain->SetBranchAddress("pt_mc",pt_mc);
    inChain->SetBranchAddress("chi2_mpd",chi2_mpd);
    inChain->SetBranchAddress("eta_mc",eta_mc);
    inChain->SetBranchAddress("px_mc",px_mc);
    inChain->SetBranchAddress("py_mc",py_mc);
    inChain->SetBranchAddress("pz_mc",pz_mc);
    inChain->SetBranchAddress("energy_mc",energy_mc);
    inChain->SetBranchAddress("q_vectors_mpd",q_vectors_mpd);
    inChain->SetBranchAddress("q_norm_mpd",q_norm_mpd);
    inChain->SetBranchAddress("eta_match_mc_mpd",eta_match_mc_mpd);
    inChain->SetBranchAddress("phi_match_mc_mpd",phi_match_mc_mpd);
    inChain->SetBranchAddress("signed_pt_match_mc_mpd",signed_pt_match_mc_mpd);
    
    Int_t nevents = inChain->GetEntries();
	for (Int_t event = 0; event < nevents; ++event)
	{
		inChain->GetEntry(event); //reading all the branches for current event
		cout << "EVENT # "<<event << endl;
		centrality->Fill(centrality_tpc_mpd); 
		Int_t multiplicity = GetMultiplicityTPC(id_from_mc_mpd, n_tracks_mpd);
		if (multiplicity == 0) continue;
		p_b_vs_multiplicity->Fill(multiplicity,b_mc);
		h_multiplicity_before->Fill(multiplicity);
	}
	
	outFile->cd();
	h_multiplicity_before->Write();
	p_b_vs_multiplicity->Write();
	centrality->Write();
	outFile->Close();
}

Int_t GetMultiplicityTPC(Long64_t* id_from_mc_mpd, Long64_t n_tracks_mpd) //should be called in loop over events
{
	Int_t multiplicity = 0;
	for (Long64_t track = 0; track < n_tracks_mpd; ++track)//loop over mpdtracks, cut on mother ID will be everywhere
	{
		if (id_from_mc_mpd[track] == -1) continue; //equivalent to mother id cut	
		multiplicity++; //multiplicity in TPC before eta, nhits and pt cuts, but after mother ID cut
	}
	return multiplicity;
}
