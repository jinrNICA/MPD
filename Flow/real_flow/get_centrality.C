#define _MAX_TRACKS 5000
#define _N_ARM 2 // 3 IS FULL DETECTOR
#define _N_HARM 2
#define _N_SORTS 4
#define _N_MODULES_TOTAL 90
#define _N_METHOD 2 // 0 - TPC , 1 - ZDC
#define _N_QCOMP 2


#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom1.h>

const Int_t NmultiplicityBins = 100;

void get_centrality(TString inFileHistName, TString inFileTreeName, TString outFileName)
{
	CentralityCalc cc = CentralityCalc(inFileHistName,inFileTreeName,outFileName);
	return;
}

class CentralityCalc
{
	public:
		CentralityCalc(TString inFileHistName, TString inFileTreeName, TString outFileName);
		Int_t GetCentrality(Int_t multiplicity);
		Int_t GetMultiplicityTPC();
	private:
		Float_t b_mc;
		Float_t x_vertex_mc;
		Float_t y_vertex_mc;
		Float_t z_vertex_mc;
		Long64_t n_tracks_mpd;
		Long64_t n_tracks_mc;
		Float_t eta_mpd[_MAX_TRACKS];
		Float_t eta_mc[_MAX_TRACKS];
		Int_t n_hits_mpd[_MAX_TRACKS];
		Int_t n_hits_poss_mpd[_MAX_TRACKS];
		Float_t pid_tpc_prob_electron_mpd[_MAX_TRACKS];
		Float_t pid_tpc_prob_pion_mpd[_MAX_TRACKS];
		Float_t pid_tpc_prob_kaon_mpd[_MAX_TRACKS];
		Float_t pid_tpc_prob_proton_mpd[_MAX_TRACKS];
		Float_t pid_tof_prob_electron_mpd[_MAX_TRACKS];
		Float_t pid_tof_prob_pion_mpd[_MAX_TRACKS];
		Float_t pid_tof_prob_kaon_mpd[_MAX_TRACKS];
		Float_t pid_tof_prob_proton_mpd[_MAX_TRACKS];
		Float_t tof_beta_mpd[_MAX_TRACKS];
		Float_t tof_mass2_mpd[_MAX_TRACKS];
		Float_t dEdx_tpc_mpd[_MAX_TRACKS];
		Float_t chi2_mpd[_MAX_TRACKS];
		Float_t pt_error_mpd[_MAX_TRACKS];
		Float_t theta_error_mpd[_MAX_TRACKS];
		Float_t phi_error_mpd[_MAX_TRACKS];
		Float_t DCA_x_mpd[_MAX_TRACKS];
		Float_t DCA_y_mpd[_MAX_TRACKS];
		Float_t DCA_z_mpd[_MAX_TRACKS];
		Float_t DCA_global_x_mpd[_MAX_TRACKS];
		Float_t DCA_global_y_mpd[_MAX_TRACKS];
		Float_t DCA_global_z_mpd[_MAX_TRACKS];
		Float_t signed_pt_mpd[_MAX_TRACKS];
		Float_t p_mpd[_MAX_TRACKS];
		Float_t pt_mc[_MAX_TRACKS];
		Int_t mother_ID_mc[_MAX_TRACKS];
		Int_t PDG_code_mc[_MAX_TRACKS];
		Float_t px_mc[_MAX_TRACKS];
		Float_t py_mc[_MAX_TRACKS];
		Float_t pz_mc[_MAX_TRACKS];
		Float_t start_x_mc[_MAX_TRACKS];
		Float_t start_y_mc[_MAX_TRACKS];
		Float_t start_z_mc[_MAX_TRACKS];
		Float_t mass_mc[_MAX_TRACKS];
		Float_t energy_mc[_MAX_TRACKS];
		Float_t phi_mpd[_MAX_TRACKS];
		Float_t theta_mpd[_MAX_TRACKS];
		Int_t TOF_flag_mpd[_MAX_TRACKS];
		Float_t ZDC_energy_mpd[_N_MODULES_TOTAL];
		Long64_t mpd_side[_MAX_TRACKS];
		Float_t phiEP_mc;
		Int_t centrality_tpc_mpd;
		
		TRandom *RNG;
		
		Float_t multiplicity_bins[NmultiplicityBins+1];
};

CentralityCalc::CentralityCalc(TString inFileHistName, TString inFileTreeName, TString outFileName)
{	
	RNG = new TRandom();

	TFile *inFileHist = new TFile(inFileHistName.Data());
	TFile *inFileTree = new TFile(inFileTreeName.Data());
	TFile *outFile = new TFile(outFileName.Data(),"RECREATE");
	
	TTree *inTree = (TTree*)inFileTree->Get("cbmsim_reduced");
	TTree *outTree = new TTree("cbmsim_reduced","cbmsim_reduced");
	
	for (int i = 0; i <= NmultiplicityBins; ++i) multiplicity_bins[i] = 1;
	
	inTree->SetBranchAddress("b_mc",&b_mc);
	inTree->SetBranchAddress("x_vertex_mc",&x_vertex_mc);
	inTree->SetBranchAddress("y_vertex_mc",&y_vertex_mc);
	inTree->SetBranchAddress("z_vertex_mc",&z_vertex_mc);
	inTree->SetBranchAddress("n_tracks_mpd",&n_tracks_mpd);
	inTree->SetBranchAddress("n_tracks_mc",&n_tracks_mc);
	inTree->SetBranchAddress("eta_mpd",eta_mpd);
	inTree->SetBranchAddress("eta_mc",eta_mc);
	inTree->SetBranchAddress("n_hits_mpd",n_hits_mpd);
	inTree->SetBranchAddress("n_hits_poss_mpd",n_hits_poss_mpd);
	inTree->SetBranchAddress("pid_tpc_prob_electron_mpd",pid_tpc_prob_electron_mpd);
	inTree->SetBranchAddress("pid_tpc_prob_pion_mpd",pid_tpc_prob_pion_mpd);
	inTree->SetBranchAddress("pid_tpc_prob_kaon_mpd",pid_tpc_prob_kaon_mpd);
	inTree->SetBranchAddress("pid_tpc_prob_proton_mpd",pid_tpc_prob_proton_mpd);
	inTree->SetBranchAddress("pid_tof_prob_electron_mpd",pid_tof_prob_electron_mpd);
	inTree->SetBranchAddress("pid_tof_prob_pion_mpd",pid_tof_prob_pion_mpd);
	inTree->SetBranchAddress("pid_tof_prob_kaon_mpd",pid_tof_prob_kaon_mpd);
	inTree->SetBranchAddress("pid_tof_prob_proton_mpd",pid_tof_prob_proton_mpd);
	inTree->SetBranchAddress("tof_beta_mpd",tof_beta_mpd);
	inTree->SetBranchAddress("tof_mass2_mpd",tof_mass2_mpd);
	inTree->SetBranchAddress("dEdx_tpc_mpd",dEdx_tpc_mpd);
	inTree->SetBranchAddress("chi2_mpd",chi2_mpd);
	inTree->SetBranchAddress("pt_error_mpd",pt_error_mpd);
	inTree->SetBranchAddress("theta_error_mpd",theta_error_mpd);
	inTree->SetBranchAddress("phi_error_mpd",phi_error_mpd);
	inTree->SetBranchAddress("DCA_x_mpd",DCA_x_mpd);
	inTree->SetBranchAddress("DCA_y_mpd",DCA_y_mpd);
	inTree->SetBranchAddress("DCA_z_mpd",DCA_z_mpd);
	inTree->SetBranchAddress("DCA_global_x_mpd",DCA_global_x_mpd);
	inTree->SetBranchAddress("DCA_global_y_mpd",DCA_global_y_mpd);
	inTree->SetBranchAddress("DCA_global_z_mpd",DCA_global_z_mpd);
	inTree->SetBranchAddress("signed_pt_mpd",signed_pt_mpd);
	inTree->SetBranchAddress("pt_mc",pt_mc);
	inTree->SetBranchAddress("mother_ID_mc",mother_ID_mc);
	inTree->SetBranchAddress("PDG_code_mc",PDG_code_mc);
	inTree->SetBranchAddress("px_mc",px_mc);
	inTree->SetBranchAddress("py_mc",py_mc);
	inTree->SetBranchAddress("pz_mc",pz_mc);
	inTree->SetBranchAddress("start_x_mc",start_x_mc);
	inTree->SetBranchAddress("start_y_mc",start_y_mc);
	inTree->SetBranchAddress("start_z_mc",start_z_mc);
	inTree->SetBranchAddress("mass_mc",mass_mc);
	inTree->SetBranchAddress("energy_mc",energy_mc);
	inTree->SetBranchAddress("phi_mpd",phi_mpd);
	inTree->SetBranchAddress("theta_mpd",theta_mpd);
	inTree->SetBranchAddress("TOF_flag_mpd",TOF_flag_mpd);
	inTree->SetBranchAddress("ZDC_energy_mpd",ZDC_energy_mpd);
	inTree->SetBranchAddress("id_from_mc_mpd",mpd_side);
	inTree->SetBranchAddress("phiEP_mc",&phiEP_mc);
	
	outTree->Branch("b_mc",&b_mc,"b_mc/F");
	outTree->Branch("x_vertex_mc",&x_vertex_mc,"x_vertex_mc/F");
	outTree->Branch("y_vertex_mc",&y_vertex_mc,"y_vertex_mc/F");
	outTree->Branch("z_vertex_mc",&z_vertex_mc,"z_vertex_mc/F");
	outTree->Branch("n_tracks_mpd",&n_tracks_mpd,"n_tracks_mpd/L");
	outTree->Branch("n_tracks_mc",&n_tracks_mc,"n_tracks_mc/L");
	outTree->Branch("eta_mpd",eta_mpd,"eta_mpd[n_tracks_mpd]/F");
	outTree->Branch("eta_mc",eta_mc,"eta_mc[n_tracks_mc]/F");
	outTree->Branch("n_hits_mpd",n_hits_mpd,"n_hits_mpd[n_tracks_mpd]/I");
	outTree->Branch("n_hits_poss_mpd",n_hits_poss_mpd,"n_hits_poss_mpd[n_tracks_mpd]/I");
	outTree->Branch("pid_tpc_prob_electron_mpd",pid_tpc_prob_electron_mpd,"pid_tpc_prob_electron_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tpc_prob_pion_mpd",pid_tpc_prob_pion_mpd,"pid_tpc_prob_pion_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tpc_prob_kaon_mpd",pid_tpc_prob_kaon_mpd,"pid_tpc_prob_kaon_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tpc_prob_proton_mpd",pid_tpc_prob_proton_mpd,"pid_tpc_prob_proton_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tof_prob_electron_mpd",pid_tof_prob_electron_mpd,"pid_tof_prob_electron_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tof_prob_pion_mpd",pid_tof_prob_pion_mpd,"pid_tof_prob_pion_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tof_prob_kaon_mpd",pid_tof_prob_kaon_mpd,"pid_tof_prob_kaon_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tof_prob_proton_mpd",pid_tof_prob_proton_mpd,"pid_tof_prob_proton_mpd[n_tracks_mpd]/F");
	outTree->Branch("tof_beta_mpd",tof_beta_mpd,"tof_beta_mpd[n_tracks_mpd]/F");
	outTree->Branch("tof_mass2_mpd",tof_mass2_mpd,"tof_mass2_mpd[n_tracks_mpd]/F");
	outTree->Branch("dEdx_tpc_mpd",dEdx_tpc_mpd,"dEdx_tpc_mpd[n_tracks_mpd]/F");
	outTree->Branch("chi2_mpd",chi2_mpd,"chi2_mpd[n_tracks_mpd]/F");
	outTree->Branch("pt_error_mpd",pt_error_mpd,"pt_error_mpd[n_tracks_mpd]/F");
	outTree->Branch("theta_error_mpd",theta_error_mpd,"theta_error_mpd[n_tracks_mpd]/F");
	outTree->Branch("phi_error_mpd",phi_error_mpd,"phi_error_mpd[n_tracks_mpd]/F");
	outTree->Branch("DCA_x_mpd",DCA_x_mpd,"DCA_x_mpd[n_tracks_mpd]/F");
	outTree->Branch("DCA_y_mpd",DCA_y_mpd,"DCA_y_mpd[n_tracks_mpd]/F");
	outTree->Branch("DCA_z_mpd",DCA_z_mpd,"DCA_z_mpd[n_tracks_mpd]/F");
	outTree->Branch("DCA_global_x_mpd",DCA_global_x_mpd,"DCA_global_x_mpd[n_tracks_mpd]/F");
	outTree->Branch("DCA_global_y_mpd",DCA_global_y_mpd,"DCA_global_y_mpd[n_tracks_mpd]/F");
	outTree->Branch("DCA_global_z_mpd",DCA_global_z_mpd,"DCA_global_z_mpd[n_tracks_mpd]/F");
	outTree->Branch("signed_pt_mpd",signed_pt_mpd,"signed_pt_mpd[n_tracks_mpd]/F");
    outTree->Branch("p_mpd",p_mpd,"p_mpd[n_tracks_mpd]/F");
	outTree->Branch("pt_mc",pt_mc,"pt_mc[n_tracks_mc]/F");
	outTree->Branch("mother_ID_mc",mother_ID_mc,"mother_ID_mc[n_tracks_mc]/I");
	outTree->Branch("PDG_code_mc",PDG_code_mc,"PDG_code_mc[n_tracks_mc]/I");
	outTree->Branch("px_mc",px_mc,"px_mc[n_tracks_mc]");
	outTree->Branch("py_mc",py_mc,"py_mc[n_tracks_mc]");
	outTree->Branch("pz_mc",pz_mc,"pz_mc[n_tracks_mc]");
	outTree->Branch("start_x_mc",start_x_mc,"start_x_mc[n_tracks_mc]");
	outTree->Branch("start_y_mc",start_y_mc,"start_y_mc[n_tracks_mc]");
	outTree->Branch("start_z_mc",start_z_mc,"start_z_mc[n_tracks_mc]");
	outTree->Branch("mass_mc",mass_mc,"mass_mc[n_tracks_mc]");
	outTree->Branch("energy_mc",energy_mc,"energy_mc[n_tracks_mc]");
	outTree->Branch("phi_mpd",phi_mpd,"phi_mpd[n_tracks_mpd]/F");
	outTree->Branch("theta_mpd",theta_mpd,"theta_mpd[n_tracks_mpd]/F");
	outTree->Branch("TOF_flag_mpd",TOF_flag_mpd,"TOF_flag_mpd[n_tracks_mpd]/I");
	outTree->Branch("ZDC_energy_mpd",ZDC_energy_mpd,"ZDC_energy_mpd[90]/F"); //////////////////////////////////SHIIIIIIIIIIIIIEEEEEEEEET
	outTree->Branch("id_from_mc_mpd",mpd_side,"id_from_mc_mpd[n_tracks_mpd]/L");
	outTree->Branch("phiEP_mc",&phiEP_mc,"phiEP_mc/F");
	outTree->Branch("centrality_tpc_mpd",&centrality_tpc_mpd,"centrality_tpc_mpd/I");
	
	TH1F *h_multiplicity_before = (TH1F*)inFileHist->Get("h_multiplicity_before");
	Float_t one_tenth = h_multiplicity_before->Integral("WIDTH") / NmultiplicityBins; //the fraction of events equal to on centrality bin
	
	Int_t n_mult_bins = h_multiplicity_before->GetNbinsX(); //total number of bins in multiplicity histogram
	multiplicity_bins[0] = 0.; // the first bin is always zero
	for (int i = 1; i <= NmultiplicityBins; ++i)
	{
		Float_t sum = one_tenth*i;
		multiplicity_bins[i] = integrate(h_multiplicity_before, n_mult_bins, sum);
	}
	//the last bin is always the last
	multiplicity_bins[NmultiplicityBins] = h_multiplicity_before->GetBinLowEdge(n_mult_bins) + h_multiplicity_before->GetBinWidth(1);
	
	for (int i = 0; i <= NmultiplicityBins; ++i) cout << "multiplicity bin = " << multiplicity_bins[i] << endl;
	
	Int_t nevents = inTree->GetEntries();
	//nevents =
	for (Int_t event = 0; event < nevents; ++event) //loop over events
	{
		cout << "EVENT # "<<event << endl;
		inTree->GetEntry(event); //reading all branches for a given event
		
		Int_t multiplicity = GetMultiplicityTPC();
		if (multiplicity == 0) continue;
		
		centrality_tpc_mpd = GetCentrality(multiplicity);
		
		for (long int j = 0; j < n_tracks_mpd; ++j)
	    {
			p_mpd[j] = TMath::Abs(signed_pt_mpd[j])/TMath::Sin(theta_mpd[j]);
		} 
		
		outTree->Fill(); 
	}
	
	outFile->cd();
	outFile->Write();
	outFile->Close();
}

Int_t integrate(TH1F *h, Int_t max_bin, Float_t sum)
{
	for (Int_t mult_bin = 1; mult_bin <= max_bin; ++mult_bin)
	{
		if (h->Integral(1,mult_bin,"WIDTH") >= sum)
		{
			//return floating-point multiplicity value
			return h->GetBinLowEdge(mult_bin) + 1 - ( (h->Integral(1,mult_bin,"WIDTH") - sum) / h->GetBinContent(mult_bin));
		}
	}
	return max_bin;
}

Int_t CentralityCalc::GetMultiplicityTPC() //should be called in loop over events
{
	Int_t multiplicity = 0;
	for (Long64_t track = 0; track < n_tracks_mpd; ++track)//loop over mpdtracks, cut on mother ID will be everywhere
	{
		if (mpd_side[track] == -1) continue; //equivalent to mother id cut	
		multiplicity++; //multiplicity in TPC before eta, nhits and pt cuts, but after mother ID cut
	}
	return multiplicity;
}

Int_t CentralityCalc::GetCentrality(Int_t multiplicity)
{
	int centrality_bin = -1;
	for (int multiplicityBin = 0; multiplicityBin < NmultiplicityBins; ++multiplicityBin)
	{
		if ((multiplicity > (Int_t)multiplicity_bins[multiplicityBin]) && (multiplicity <= (Int_t)multiplicity_bins[multiplicityBin+1])) //inclusive borders
		{
			if ((multiplicity == (Int_t)multiplicity_bins[multiplicityBin+1]))
			{
				if ((multiplicity == (Int_t)multiplicity_bins[multiplicityBin+2]))
				{
					
					Double_t random_number = RNG->Uniform();
					if (((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) > multiplicity_bins[multiplicityBin+1]) &&((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) < multiplicity_bins[multiplicityBin+2]))
						return multiplicityBin+1;
					else if (((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) > multiplicity_bins[multiplicityBin+2]) /*&&((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) < multiplicity_bins[multiplicityBin+3])*/)
						return multiplicityBin+2;
					else 
						return multiplicityBin;
				}
				if ((RNG->Uniform() +(Int_t)multiplicity_bins[multiplicityBin+1]) > multiplicity_bins[multiplicityBin+1]) 
					return multiplicityBin+1; 
				else 
					return multiplicityBin;
				//~ Int_t n = 0;
				//~ while (multiplicity == (Int_t)multiplicity_bins[multiplicityBin + 1 + n]) ++n;
				//~ cout << "multiplicity = "<< multiplicity <<"; n = " << n << endl;
				//~ Float_t *array = new Float_t[n+2];
				//~ array[0] = (Int_t)multiplicity_bins[multiplicityBin + 1];
				//~ for (int k = 1; k < n+2; ++k) array[k] = multiplicity_bins[k];
				//~ Double_t random_number = RNG->Uniform();
				//~ for (int j = 0; j < n + 1; ++j) if ( (((Int_t)multiplicity_bins[multiplicityBin+1] + random_number) > array[j]) &&
					//~ (((Int_t)multiplicity_bins[multiplicityBin+1] + random_number) < array[j+1]) ) return  multiplicityBin + j;
				//~ delete[] array;
			}
			return multiplicityBin;
		}
	}
	return centrality_bin;
}
//ASS
