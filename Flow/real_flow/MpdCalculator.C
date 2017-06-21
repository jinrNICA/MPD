//!!!!!!!!!!!!!!!!!!!!!!
#include "MpdCalculator.h"

MpdCalculator::MpdCalculator(TString inFileName , TString outFileName, TString dcaFileName)
{
	
	this->inFileName = inFileName;
	this->outFileName = outFileName;

	inChain = new TChain("cbmsim_reduced");
	inChain->Add(inFileName.Data());

	outFile = new TFile(outFileName.Data(), "RECREATE");
	
	dcaFile = new TFile(dcaFileName.Data(),"READ");

	inChain->SetBranchAddress("b_mc",&b_mc);
	inChain->SetBranchAddress("centrality_tpc_mpd",&centrality_tpc_mpd);
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
	inChain->SetBranchAddress("ZDC_energy_mpd",ZDC_energy_mpd);;
    inChain->SetBranchAddress("pt_mc",pt_mc);
    inChain->SetBranchAddress("chi2_mpd",chi2_mpd);
    inChain->SetBranchAddress("eta_mc",eta_mc);
    inChain->SetBranchAddress("px_mc",px_mc);
    inChain->SetBranchAddress("py_mc",py_mc);
    inChain->SetBranchAddress("pz_mc",pz_mc);
    inChain->SetBranchAddress("energy_mc",energy_mc);
    inChain->SetBranchAddress("DCA_x_mpd",DCA_x_mpd);
    inChain->SetBranchAddress("DCA_y_mpd",DCA_y_mpd);
    inChain->SetBranchAddress("DCA_z_mpd",DCA_z_mpd);

	h_nhits_TPC = new TH1F("h_nhits_TPC","number of hits in TPC distribution;N_{hits};N_{tracks};",60,0.,60.);
	h_nhits_TPC_after = new TH1F("h_nhits_TPC_after","number of hits in TPC distribution after cuts;N_{hits};N_{tracks};",60,0.,60.);
	h_multiplicity = new TH1F("h_multiplicity","multiplicity in TPC;N_{tracks} in TPC;N_{events};",1600,0.,1600.);
	h_multiplicity_before = new TH1F("h_multiplicity_before","multiplicity in TPCbefore cuts;N_{tracks} in TPC;N_{events};",1600,0.,1600.);
	h_impact_parameter = new TH1F("h_impact_parameter","impact parameter distribution;b,fm;N_{events};",100,0.,20.);
	h2_energyZDC_L_vs_energy_ZDC_R = new TH2F("h2_energyZDC_L_vs_energy_ZDC_R","FHCal energy left vs right;E_{left};E_{right};",100,0.,22.,100,0.,22.);
	p_b_vs_multiplicity = new TProfile("p_b_vs_multiplicity","Impact parameter vs TPC multiplicity",800,0,1600);
	p_b_vs_multiplicity->Sumw2();
	p_b_vs_energy = new TProfile("p_b_vs_energy","Impact parameter vs FHCal energy",200,0,70);
	p_b_vs_energy->Sumw2();
	h2_b_vs_centrality = new TH2F("h2_b_vs_centrality","h2_b_vs_centrality",100,0,100,200,0.,20.);
	
	char name[200];
	char title[200];
	
	for (Int_t dim = 0; dim < Ndim; dim++)
	{	
		for (Int_t eta_bin=0; eta_bin<NetaBins; eta_bin++){
			f_pt_fit[dim][eta_bin] = (TF1*) dcaFile->Get(Form("sigma_pt_fit%i%i",dim,eta_bin));
		}
		for (Int_t pt_bin = 0; pt_bin < NptBins; pt_bin++)
		{
			for (Int_t eta_bin = 0; eta_bin < NetaBins ; eta_bin++ )
			{
				//f_dca[dim][pt_bin][eta_bin] = (TF1*) dcaFile->Get(Form("dca_fit[%i][%i][%i]",dim,pt_bin,eta_bin));				

				sprintf(name,"h_DCA_primary[%i][%i][%i]",dim,pt_bin,eta_bin);
				sprintf(title,"DCA_{%i} primary for %.2f < p_{t} < %.2f and %.2f < #eta < %.2f ", dim,ptBins[pt_bin], ptBins[pt_bin+1],etaBins[eta_bin],etaBins[eta_bin+1]);
				h_DCA_primary[dim][pt_bin][eta_bin] = new TH1F(name,title,200,-50.,50.);
				
				sprintf(name,"h_DCA_secondary[%i][%i][%i]",dim,pt_bin,eta_bin);
				sprintf(title,"DCA_{%i} secondary for %.2f < p_{t} < %.2f and %.2f < #eta < %.2f", dim,ptBins[pt_bin], ptBins[pt_bin+1],etaBins[eta_bin],etaBins[eta_bin+1]);
				h_DCA_secondary[dim][pt_bin][eta_bin] = new TH1F(name,title,200,-50.,50.);
				
				sprintf(name,"h_DCA_all[%i][%i][%i]",dim,pt_bin,eta_bin);
				sprintf(title,"DCA_{%i} all for %.2f < p_{t} < %.2f and %.2f < #eta < %.2f", dim,ptBins[pt_bin], ptBins[pt_bin+1],etaBins[eta_bin],etaBins[eta_bin+1]);
				h_DCA_all[dim][pt_bin][eta_bin] = new TH1F(name,title,200,-50.,50.);
			}
		}
	}
	
	for (int sort = 0; sort < _N_SORTS; ++sort)
	{
		sprintf(name,"p_momenta_resolution[%i]",sort);
		sprintf(title,"#Delta p_{t} / p_{t} for %s", sorts_of_particles[sort].Data());
		p_momenta_resolution[sort] = new TProfile(name,title,NptBins,ptBins);
	}

	for (int harm = 0; harm < _N_HARM; ++harm)
	{
		sprintf(name,"p_flow_wrt_RP_vs_centrality[%i]",harm);
		sprintf(title,"v_{%i} wrt RP for ;centrailty, %;",harm+1);
		p_flow_wrt_RP_vs_centrality[harm] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
		
		for (int cenralityBin = 0; cenralityBin < NcentralityBinsFlow; ++cenralityBin)
		{	
			sprintf(name,"p_flow_wrt_RP_vs_pt[%i][%i]",cenralityBin,harm);
			sprintf(title,"v_{%i} wrt RP for %.2f < b < %.2f;p_{t};",harm+1,centralityBinsFlow[cenralityBin],centralityBinsFlow[cenralityBin+1]);
			p_flow_wrt_RP_vs_pt[cenralityBin][harm] = new TProfile(name,title,NptBins,ptBins);

			sprintf(name,"p_flow_wrt_RP_vs_eta[%i][%i]",cenralityBin,harm);
			sprintf(title,"v_{%i} wrt RP for %.2f < b < %.2f;#eta;",harm+1,centralityBinsFlow[cenralityBin],centralityBinsFlow[cenralityBin+1]);
			p_flow_wrt_RP_vs_eta[cenralityBin][harm] = new TProfile(name,title,NetaBins,etaBins);

			sprintf(name,"p_flow_wrt_RP_vs_rapidity[%i][%i]",cenralityBin,harm);
			sprintf(title,"v_{%i} wrt RP for %.2f < b < %.2f;#y;",harm+1,centralityBinsFlow[cenralityBin],centralityBinsFlow[cenralityBin+1]);
			p_flow_wrt_RP_vs_rapidity[cenralityBin][harm] = new TProfile(name, title, NrapidityBins, rapidityBins);
		}
		
		for (int _harm = 0; _harm < _N_HARM; ++ _harm)
		{
			for (int method = 0; method < _N_METHOD; ++method)
			{
				sprintf(name,"p_true_Res_vs_b[%i][%i][%i]",harm,_harm,method);
				sprintf(title,"%s%i%s%i%s%s%s","<cos",harm+1,"(#Psi_{",_harm + 1,"}^{FULL} - #Psi_{RP})> using ",
						methods_names[method].Data(),";b,fm;");
				p_true_Res_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);

				sprintf(name,"%s%i%s%i%s%i%s","p_Res2Psi_vs_b[",harm,"][",_harm,"][",method,"]");
				sprintf(title,"%s%i%s%i%s%s%s%i%s%s%s","cos(",harm + 1,"(#Psi_{",_harm + 1,",",methods_names[method].Data(),"}^{R} - #Psi_{",
						_harm + 1,",",methods_names[method].Data(),"}^{L});b,fm;");
				p_Res2Psi_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
				
				sprintf(name,"p_flow_wrt_full_vs_centrality[%i][%i][%i]",harm,_harm,method);
				sprintf(title,"v_{%i} wrt #Psi_{%i,%s}^{FULL} ;centrality, %;", harm+1,_harm+1,methods_names[method].Data());
				p_flow_wrt_full_vs_centrality[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);

				sprintf(name,"p_flow_wrt_full_vs_centrality_divided[%i][%i][%i]",harm,_harm,method);
				sprintf(title,"v_{%i} wrt #Psi_{%i,%s}^{FULL} divided ;centrality, %;",harm+1,_harm+1,methods_names[method].Data());
				p_flow_wrt_full_vs_centrality_divided[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
				
				for (int cenralityBin = 0; cenralityBin < NcentralityBinsFlow; ++cenralityBin)
				{	
					sprintf(name,"p_flow_wrt_full_vs_pt[%i][%i][%i][%i]",cenralityBin,harm,_harm,method);
					sprintf(title,"v_{%i} wrt #Psi_{%i,%s}^{FULL} for %.2f < b < %.2f;p_{t};", harm+1,_harm+1,methods_names[method].Data(),
							centralityBinsFlow[cenralityBin], centralityBinsFlow[cenralityBin+1]);
					p_flow_wrt_full_vs_pt[cenralityBin][harm][_harm][method] = new TProfile(name,title,NptBins,ptBins);

					sprintf(name,"p_flow_wrt_full_vs_pt_divided[%i][%i][%i][%i]",cenralityBin,harm,_harm,method);
					sprintf(title,"v_{%i} wrt #Psi_{%i,%s}^{FULL} divided for %.2f < b < %.2f;p_{t};"
							,harm+1,_harm+1,methods_names[method].Data(),centralityBinsFlow[cenralityBin], centralityBinsFlow[cenralityBin+1]);
					p_flow_wrt_full_vs_pt_divided[cenralityBin][harm][_harm][method] = new TProfile(name,title,NptBins,ptBins);

					sprintf(name,"p_flow_wrt_full_vs_eta[%i][%i][%i][%i]",cenralityBin,harm,_harm,method);
					sprintf(title,"v_{%i} wrt #Psi_{%i,%s}^{FULL} for %.2f < b < %.2f;#eta;",harm+1,_harm+1,methods_names[method].Data(),
							centralityBinsFlow[cenralityBin], centralityBinsFlow[cenralityBin+1]);
					p_flow_wrt_full_vs_eta[cenralityBin][harm][_harm][method] = new TProfile(name,title,NetaBins,etaBins);

					sprintf(name,"p_flow_wrt_full_vs_eta_divided[%i][%i][%i][%i]",cenralityBin,harm,_harm,method);
					sprintf(title,"v_{%i} wrt #Psi_{%i,%s}^{FULL} divided for %.2f < b < %.2f;#eta;",harm+1,_harm+1,methods_names[method].Data(),
							centralityBinsFlow[cenralityBin], centralityBinsFlow[cenralityBin+1]);
					p_flow_wrt_full_vs_eta_divided[cenralityBin][harm][_harm][method] = new TProfile(name,title,NetaBins,etaBins);

					sprintf(name,"p_flow_wrt_full_vs_rapidity[%i][%i][%i][%i]",cenralityBin,harm,_harm,method);
					sprintf(title,"v_{%i} wrt #Psi_{%i,%s}^{FULL} for %.2f < b < %.2f;#y;",harm+1,_harm+1,methods_names[method].Data(),
							centralityBinsFlow[cenralityBin], centralityBinsFlow[cenralityBin+1]);
					p_flow_wrt_full_vs_rapidity[cenralityBin][harm][_harm][method] = new TProfile(name,title,NrapidityBins,rapidityBins);

					sprintf(name,"p_flow_wrt_full_vs_rapidity_divided[%i][%i][%i][%i]",cenralityBin,harm,_harm,method);
					sprintf(title,"v_{%i} wrt #Psi_{%i,%s}^{FULL} divided for %.2f < b < %.2f;y;",harm+1,_harm+1,methods_names[method].Data(),
							centralityBinsFlow[cenralityBin], centralityBinsFlow[cenralityBin+1]);
					p_flow_wrt_full_vs_rapidity_divided[cenralityBin][harm][_harm][method] = new TProfile(name,title,NrapidityBins,rapidityBins);
				}
			}		
		}
	}

	for (int arm = 0; arm < _N_ARM; ++arm)
	{

		sprintf(name,"%s%i%s","h_energy_ZDC_total[",arm,"]");
		sprintf(title,"%s%s%s","FHCal energy distribution for ",arm_names[arm].Data()," arm;E_{total};N_{events};");
		h_energy_ZDC_total[arm] = new TH1F(name,title,100,0.,22.);

		sprintf(name,"%s%i%s","h2_energy_ZDC_vs_multiplicity[",arm,"]");
		sprintf(title,"%s%s%s","Total FHCal energy in ",arm_names[arm].Data()," arm vs TPC multiplicity;E_{total};multiplicity;");
		h2_energy_ZDC_vs_multiplicity[arm] = new TH2F(name,title,100,0.,22.,100,0.,600.);

		for (int harm = 0; harm < _N_HARM; ++harm)
		{
			for (int method = 0; method < _N_METHOD; ++method)
			{
				sprintf(name,"%s%i%s%i%s%i%s","p_qx_vs_b[",arm,"][",harm,"][",method,"]");
				sprintf(title,"%s%i%s%s%s%s","Q_{x}^{",harm+1,",",arm_names[arm].Data(),"} using ",methods_names[method].Data());
				p_qx_vs_b[arm][harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);

				sprintf(name,"%s%i%s%i%s%i%s","p_qy_vs_b[",arm,"][",harm,"][",method,"]");
				sprintf(title,"%s%i%s%s%s%s","Q_{y}^{",harm+1,",",arm_names[arm].Data(),"} using ",methods_names[method].Data());
				p_qy_vs_b[arm][harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);

				for(int _harm = 0; _harm < _N_HARM; ++ _harm)
				{
					sprintf(name,"p_true_Res_half_vs_b[%i][%i][%i][%i]",arm,harm,_harm,method);
					sprintf(title,"<cos(%i(#Psi_{%i,%s}^{%s} - #Psi_{RP} ))> ;b,fm;",harm+1,_harm+1,methods_names[method].Data(),arm_names[arm].Data());
					p_true_Res_half_vs_b[arm][harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
				}
			}
		}
	}

	for (int centralityBin = 0; centralityBin < NcentralityBinsRes; ++centralityBin)
	{
		for (int sort = 0; sort < _N_SORTS; ++sort)
		{
			sprintf(name,"%s%i%s%i%s","h_pt[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","P_{t} distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";p_{t};N_{tracks};");
			h_pt[centralityBin][sort] = new TH1F(name,title,100,0.,3.5);

			sprintf(name,"%s%i%s%i%s","h_eta[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","#eta distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";#eta;N_{tracks};");
			h_eta[centralityBin][sort] = new TH1F(name,title,100,-2.,2.);

			sprintf(name,"%s%i%s%i%s","h_phi[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","#phi distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";#phi;N_{tracks};");
			h_phi[centralityBin][sort] = new TH1F(name,title,100,-4.,4.);

			sprintf(name,"%s%i%s%i%s","h2_pt_vs_eta[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","p_{t} vs #eta distribution for ",sorts_of_particles[sort].Data()," for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";p_{t};#eta;");
			h2_pt_vs_eta[centralityBin][sort] = new TH2F(name,title,100,0.,3.5,100.,-1.5,1.5);

			sprintf(name,"%s%i%s%i%s","h2_pt_vs_phi[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","p_{t} vs #phi distribution for ",sorts_of_particles[sort].Data()," for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";p_{t};#phi;");
			h2_pt_vs_phi[centralityBin][sort] = new TH2F(name,title,100,0.,3.5,100.,-4.,4.);

			sprintf(name,"%s%i%s%i%s","h2_phi_vs_eta[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","#phi vs #eta distribution for ",sorts_of_particles[sort].Data()," for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";#phi;#eta;");
			h2_phi_vs_eta[centralityBin][sort] = new TH2F(name,title,100,-4.,4.,100.,-1.5,1.5);//
			
			sprintf(name,"%s%i%s%i%s","h_pt_after[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","P_{t} distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";p_{t};N_{tracks};");
			h_pt_after[centralityBin][sort] = new TH1F(name,title,100,0.,3.5);

			sprintf(name,"%s%i%s%i%s","h_eta_after[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","#eta distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";#eta;N_{tracks};");
			h_eta_after[centralityBin][sort] = new TH1F(name,title,100,-2.,2.);

			sprintf(name,"%s%i%s%i%s","h_phi_after[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","#phi distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";#phi;N_{tracks};");
			h_phi_after[centralityBin][sort] = new TH1F(name,title,100,-4.,4.);

			sprintf(name,"%s%i%s%i%s","h2_pt_vs_eta_after[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","p_{t} vs #eta distribution for ",sorts_of_particles[sort].Data()," for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";p_{t};#eta;");
			h2_pt_vs_eta_after[centralityBin][sort] = new TH2F(name,title,100,0.,3.5,100.,-1.5,1.5);

			sprintf(name,"%s%i%s%i%s","h2_pt_vs_phi_after[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","p_{t} vs #phi distribution for ",sorts_of_particles[sort].Data()," for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";p_{t};#phi;");
			h2_pt_vs_phi_after[centralityBin][sort] = new TH2F(name,title,100,0.,3.5,100.,-4.,4.);

			sprintf(name,"%s%i%s%i%s","h2_phi_vs_eta_after[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","#phi vs #eta distribution for ",sorts_of_particles[sort].Data()," for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";#phi;#eta;");
			h2_phi_vs_eta_after[centralityBin][sort] = new TH2F(name,title,100,-4.,4.,100.,-1.5,1.5);//
			
			sprintf(name,"%s%i%s%i%s","h_pt_mc[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","P_{t} distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";p_{t};N_{tracks};");
			h_pt_mc[centralityBin][sort] = new TH1F(name,title,100,0.,3.5);

			sprintf(name,"%s%i%s%i%s","h_eta_mc[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","#eta distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";#eta;N_{tracks};");
			h_eta_mc[centralityBin][sort] = new TH1F(name,title,100,-2.,2.);

			sprintf(name,"%s%i%s%i%s","h_phi_mc[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","#phi distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";#phi;N_{tracks};");
			h_phi_mc[centralityBin][sort] = new TH1F(name,title,100,-4.,4.);
			
			sprintf(name,"%s%i%s%i%s","h_pt_mc_after[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","P_{t} distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";p_{t};N_{tracks};");
			h_pt_mc_after[centralityBin][sort] = new TH1F(name,title,100,0.,3.5);

			sprintf(name,"%s%i%s%i%s","h_eta_mc_after[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","#eta distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";#eta;N_{tracks};");
			h_eta_mc_after[centralityBin][sort] = new TH1F(name,title,100,-2.,2.);

			sprintf(name,"%s%i%s%i%s","h_phi_mc_after[",centralityBin,"][",sort,"]");
			sprintf(title,"%s%s%s%.2f%s%.2f%s","#phi distribution for ",sorts_of_particles[sort].Data(),
					" for ",centralityBinsRes[centralityBin]," < b < ",centralityBinsRes[centralityBin+1],";#phi;N_{tracks};");
			h_phi_mc_after[centralityBin][sort] = new TH1F(name,title,100,-4.,4.);
		}
	}
}

void MpdCalculator::FillTPC(Int_t cenrality_bin, EPParticle *event_plane_buffer, Int_t ep_particle_count)
{

	Double_t qx, qy;
	Double_t psi_N_HALF;
	phiEP_mc = ATan2(Sin(phiEP_mc), Cos(phiEP_mc));

	for (Int_t harm = 0; harm < _N_HARM; ++harm)
	{
		Double_t psi_N_R = GetPsiHalfTpc(event_plane_buffer,ep_particle_count, 1, harm+1, qx, qy);
		Double_t psi_N_L = GetPsiHalfTpc(event_plane_buffer,ep_particle_count, -1, harm+1, qx, qy);

		Double_t psi_N_FULL = GetPsiFullTpc(event_plane_buffer,ep_particle_count, harm+1);
		for (Int_t _harm = 0; _harm < _N_HARM; ++_harm)
		{
			Double_t psi_N_R = GetPsiHalfTpc(event_plane_buffer,ep_particle_count, 1, _harm+1, qx, qy);
			Double_t psi_N_L = GetPsiHalfTpc(event_plane_buffer,ep_particle_count, -1, _harm+1, qx, qy);
			// << "Filling with centrality = " << centralityBinsRes[cenrality_bin]+0.1 << "Res2 = " << Cos((harm+1)*(psi_N_R - psi_N_L)) << endl;
			p_Res2Psi_vs_b[harm][_harm][0]->Fill(centralityBinsRes[cenrality_bin]+0.1, Cos((harm+1)*(psi_N_R - psi_N_L)));

			Double_t psi_N_FULL = GetPsiFullTpc(event_plane_buffer,ep_particle_count, harm+1);
			p_true_Res_vs_b[harm][_harm][0]->Fill(centralityBinsRes[cenrality_bin]+0.1,Cos((harm+1)*(psi_N_FULL - phiEP_mc)));

			for (Int_t arm = 0; arm < _N_ARM; ++arm)
			{
				Double_t psi_N_HALF = GetPsiHalfTpc(event_plane_buffer,ep_particle_count, -2*arm + 1, _harm+1, qx, qy);
				p_true_Res_half_vs_b[arm][harm][_harm][0]->Fill(centralityBinsRes[cenrality_bin]+0.1,Cos((harm+1)*(psi_N_HALF - phiEP_mc)));
			}
		}

		for (Int_t arm = 0; arm < _N_ARM; ++arm)
		{
			psi_N_HALF = GetPsiHalfTpc(event_plane_buffer,ep_particle_count, -2*arm + 1, harm+1, qx, qy);
			p_qx_vs_b[arm][harm][0]->Fill(centralityBinsRes[cenrality_bin] + 0.1, qx);
			p_qy_vs_b[arm][harm][0]->Fill(centralityBinsRes[cenrality_bin] + 0.1, qy);
		}
	}
}

void MpdCalculator::FillZDC(Int_t cenrality_bin, Float_t *ZDC_energy_mpd)
{
	Float_t values[10] , absvalues[10];
	Double_t qx,qy;
	Double_t psi_N_HALF;
	phiEP_mc = ATan2(Sin(phiEP_mc), Cos(phiEP_mc));

	for (Int_t harm = 0; harm < _N_HARM; ++harm)
	{
		for (Int_t _harm = 0; _harm < _N_HARM; ++_harm)
		{
			Double_t psi_N_R = GetPsiHalfZdc(ZDC_energy_mpd, 0, _harm + 1, qx, qy);
			Double_t psi_N_L = GetPsiHalfZdc(ZDC_energy_mpd, 1, _harm + 1, qx, qy);
			//cout << "Filling with centrality = " << centralityBinsRes[cenrality_bin]+0.1 << "Res2 = " << Cos((harm+1)*(psi_N_R - psi_N_L)) << endl;
			p_Res2Psi_vs_b[harm][_harm][1]->Fill(centralityBinsRes[cenrality_bin]+0.1, Cos((harm+1)*(psi_N_R - psi_N_L)));

			Double_t psi_N_FULL = GetPsiFullZdc(ZDC_energy_mpd, _harm+1);
			p_true_Res_vs_b[harm][_harm][1]->Fill(centralityBinsRes[cenrality_bin]+0.1,Cos((harm+1)*(psi_N_FULL - phiEP_mc)));

			for (Int_t arm = 0; arm < _N_ARM; ++arm)
			{
				Double_t psi_N_HALF = GetPsiHalfZdc(ZDC_energy_mpd, arm, _harm + 1, qx, qy);
				p_true_Res_half_vs_b[arm][harm][_harm][1]->Fill(centralityBinsRes[cenrality_bin]+0.1,Cos((harm+1)*(psi_N_HALF - phiEP_mc)));
			}
		}

		for (Int_t arm = 0; arm < _N_ARM; ++arm)
		{
			psi_N_HALF = GetPsiHalfZdc(ZDC_energy_mpd, arm, harm + 1, qx, qy);
			p_qx_vs_b[arm][harm][1]->Fill(centralityBinsRes[cenrality_bin] + 0.1, qx);
			p_qy_vs_b[arm][harm][1]->Fill(centralityBinsRes[cenrality_bin], qy);
		}
	}
}

void MpdCalculator::CalculateResolutions(Int_t nevents)
{

	EPParticle *event_plane_buffer = new EPParticle[10000];
	if (nevents == 0) nevents = inChain->GetEntries();
	for (Int_t event = 0; event < nevents; ++event)
	{
		inChain->GetEntry(event); //reading all the branches for current event
		cout << "EVENT # "<<event << endl;

		centrality_tpc_mpd = 100 - centrality_tpc_mpd; //FIXED CENTRALITY!!!!!!!!!!!!!!!!!!!!!!
	
		h_impact_parameter->Fill(b_mc); //filling impact parameter hitsogram with no cuts on events
		//Int_t multiplicity = GetMultiplicityTPC();
		if (n_tracks_mpd == 0) continue;//skipping empty events
		h_multiplicity_before->Fill(n_tracks_mpd); //filling multiplicity TPC before cuts hitsogram
		p_b_vs_multiplicity->Fill(n_tracks_mpd,b_mc);
		h2_b_vs_centrality->Fill(centrality_tpc_mpd,b_mc);
		
		int cenrality_bin_res = GetCentralityBinRes(centrality_tpc_mpd);//getting the resolution bin in whic the event is 
		if (cenrality_bin_res == -1) continue; // just in case...
		
		for (Long64_t track = 0; track < n_tracks_mc; ++track)
		{
				 //~ _      _  _                                              _             __                                 _        
	 //~ | |    (_)| |                                            | |           / _|                               | |       
	 //~ | |__   _ | |_  ___    __ _  _ __  __ _  _ __ ___   ___  | |__    ___ | |_  ___   _ __  ___    ___  _   _ | |_  ___ 
	 //~ | '_ \ | || __|/ _ \  / _` || '__|/ _` || '_ ` _ \ / __| | '_ \  / _ \|  _|/ _ \ | '__|/ _ \  / __|| | | || __|/ __|
	 //~ | | | || || |_| (_) || (_| || |  | (_| || | | | | |\__ \ | |_) ||  __/| | | (_) || |  |  __/ | (__ | |_| || |_ \__ \
	 //~ |_| |_||_| \__|\___/  \__, ||_|   \__,_||_| |_| |_||___/ |_.__/  \___||_|  \___/ |_|   \___|  \___| \__,_| \__||___/
							//~ __/ |                                                                                        
						   //~ |___/    
			int sort = 0;
			h_phi_mc[cenrality_bin_res][0]->Fill(ATan2(px_mc[track],py_mc[track]));
			h_eta_mc[cenrality_bin_res][0]->Fill(eta_mc[track]);
			h_pt_mc[cenrality_bin_res][0]->Fill(pt_mc[track]);
			
			switch (PDG_code_mc[id_from_mc_mpd[track]]) {
			case 211:
				sort = 1;
				break;
			case 2212:
				sort = 2;
				break;
			case 321:
				sort = 3;
				break;
			default:
				break;
			}
			
			if ((sort == 1) || (sort == 2) || (sort == 3))
			{
				h_phi_mc[cenrality_bin_res][sort]->Fill(ATan2(px_mc[track],py_mc[track]));
				h_eta_mc[cenrality_bin_res][sort]->Fill(eta_mc[track]);
				h_pt_mc[cenrality_bin_res][sort]->Fill(pt_mc[track]);
			}
			
								  //~ | |       
		   //~ ___  _   _ | |_  ___ 
		  //~ / __|| | | || __|/ __|
		 //~ | (__ | |_| || |_ \__ \
		  //~ \___| \__,_| \__||___/
			
			Int_t pt_bin = GetPtBin(pt_mc[track]);
			Int_t eta_bin = GetEtaBin(eta_mc[track]);
			if ( (eta_bin == -1) || (pt_bin == -1 )) continue;
			if (mother_ID_mc[track] > -1) continue;
			
						//~ | |    (_)     | |                                                    / _|| |                           | |       
 //~ | |__   _  ___ | |_  ___    __ _  _ __  __ _  _ __ ___   ___    __ _ | |_ | |_  ___  _ __    ___  _   _ | |_  ___ 
 //~ | '_ \ | |/ __|| __|/ _ \  / _` || '__|/ _` || '_ ` _ \ / __|  / _` ||  _|| __|/ _ \| '__|  / __|| | | || __|/ __|
 //~ | | | || |\__ \| |_| (_) || (_| || |  | (_| || | | | | |\__ \ | (_| || |  | |_|  __/| |    | (__ | |_| || |_ \__ \
 //~ |_| |_||_||___/ \__|\___/  \__, ||_|   \__,_||_| |_| |_||___/  \__,_||_|   \__|\___||_|     \___| \__,_| \__||___/
                             //~ __/ |                                                                                 
                            //~ |___/  
            sort = 0;
			h_phi_mc_after[cenrality_bin_res][0]->Fill(ATan2(px_mc[track],py_mc[track]));
			h_eta_mc_after[cenrality_bin_res][0]->Fill(eta_mc[track]);
			h_pt_mc_after[cenrality_bin_res][0]->Fill(pt_mc[track]);
			
			switch (PDG_code_mc[track]) {
			case 211:
				sort = 1;
				break;
			case 2212:
				sort = 2;
				break;
			case 321:
				sort = 3;
				break;
			default:
				break;
			}
			
			if ((sort == 1) || (sort == 2) || (sort == 3))
			{
				h_phi_mc_after[cenrality_bin_res][sort]->Fill(ATan2(px_mc[track],py_mc[track]));
				h_eta_mc_after[cenrality_bin_res][sort]->Fill(eta_mc[track]);
				h_pt_mc_after[cenrality_bin_res][sort]->Fill(pt_mc[track]);
			}
		}
		
		int ep_particle_count = 0, multiplicity_tpc_after = 0 ;
		for (Long64_t track = 0; track < n_tracks_mpd; ++track)//loop over mpdtracks
		{
	 //~ _      _  _                                              _             __                                 _        
	 //~ | |    (_)| |                                            | |           / _|                               | |       
	 //~ | |__   _ | |_  ___    __ _  _ __  __ _  _ __ ___   ___  | |__    ___ | |_  ___   _ __  ___    ___  _   _ | |_  ___ 
	 //~ | '_ \ | || __|/ _ \  / _` || '__|/ _` || '_ ` _ \ / __| | '_ \  / _ \|  _|/ _ \ | '__|/ _ \  / __|| | | || __|/ __|
	 //~ | | | || || |_| (_) || (_| || |  | (_| || | | | | |\__ \ | |_) ||  __/| | | (_) || |  |  __/ | (__ | |_| || |_ \__ \
	 //~ |_| |_||_| \__|\___/  \__, ||_|   \__,_||_| |_| |_||___/ |_.__/  \___||_|  \___/ |_|   \___|  \___| \__,_| \__||___/
							//~ __/ |                                                                                        
						   //~ |___/  
			Int_t pt_bin = GetPtBin(Abs(signed_pt_mpd[track]));
			Int_t eta_bin = GetEtaBin(eta_mpd[track]); 
			
			h_DCA_all[0][pt_bin][eta_bin]->Fill(DCA_x_mpd[track]);
			h_DCA_all[1][pt_bin][eta_bin]->Fill(DCA_y_mpd[track]);
			h_DCA_all[2][pt_bin][eta_bin]->Fill(DCA_z_mpd[track]);
			if (mother_ID_mc[id_from_mc_mpd[track]] == -1)
			{    
				h_DCA_primary[0][pt_bin][eta_bin]->Fill(DCA_x_mpd[track]);
				h_DCA_primary[1][pt_bin][eta_bin]->Fill(DCA_y_mpd[track]);
				h_DCA_primary[2][pt_bin][eta_bin]->Fill(DCA_z_mpd[track]);
			}
			else if (mother_ID_mc[id_from_mc_mpd[track]] > -1)
			{
				h_DCA_secondary[0][pt_bin][eta_bin]->Fill(DCA_x_mpd[track]);
				h_DCA_secondary[1][pt_bin][eta_bin]->Fill(DCA_y_mpd[track]);
				h_DCA_secondary[2][pt_bin][eta_bin]->Fill(DCA_z_mpd[track]);
			}
			                                                                                 
			h_nhits_TPC->Fill(n_hits_mpd[track]); //number of hits in TPC per track, before eta, nhits and pt cuts, but after mother ID cut
			int sort = 0;
			h_phi[cenrality_bin_res][0]->Fill(ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])));
			h_eta[cenrality_bin_res][0]->Fill(eta_mpd[track]);
			h_pt[cenrality_bin_res][0]->Fill(Abs(signed_pt_mpd[track]));
			h2_pt_vs_eta[cenrality_bin_res][0]->Fill(Abs(signed_pt_mpd[track]),eta_mpd[track]);
			h2_pt_vs_phi[cenrality_bin_res][0]->Fill(Abs(signed_pt_mpd[track]),ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])));
			h2_phi_vs_eta[cenrality_bin_res][0]->Fill(ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])),eta_mpd[track]);
			p_momenta_resolution[0]->Fill(Abs(signed_pt_mpd[track]),Abs(Abs(signed_pt_mpd[track]) - pt_mc[id_from_mc_mpd[track]])
			/Abs(signed_pt_mpd[track]));
				
			switch (PDG_code_mc[id_from_mc_mpd[track]]) {
			case 211:
				sort = 1;
				break;
			case 2212:
				sort = 2;
				break;
			case 321:
				sort = 3;
				break;
			default:
				break;
			}
			
			if ((sort == 1) || (sort == 2) || (sort == 3))
			{
				h_phi[cenrality_bin_res][sort]->Fill(ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])));
				h_eta[cenrality_bin_res][sort]->Fill(eta_mpd[track]);
				h_pt[cenrality_bin_res][sort]->Fill(Abs(signed_pt_mpd[track]));
				h2_pt_vs_eta[cenrality_bin_res][sort]->Fill(Abs(signed_pt_mpd[track]),eta_mpd[track]);
				h2_pt_vs_phi[cenrality_bin_res][sort]->Fill(Abs(signed_pt_mpd[track]),ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])));
				h2_phi_vs_eta[cenrality_bin_res][sort]->Fill(ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])),eta_mpd[track]);
				p_momenta_resolution[sort]->Fill(Abs(signed_pt_mpd[track]),Abs(Abs(signed_pt_mpd[track]) - pt_mc[id_from_mc_mpd[track]])
				/Abs(signed_pt_mpd[track]));
			}
			     
					  //~ | |       
		   //~ ___  _   _ | |_  ___ 
		  //~ / __|| | | || __|/ __|
		 //~ | (__ | |_| || |_ \__ \
		  //~ \___| \__,_| \__||___/
			
			//if (id_from_mc_mpd[track] == -1) continue; //equivalent to mother id cut
			pt_bin = GetPtBin(Abs(signed_pt_mpd[track]));
			eta_bin = GetEtaBin(eta_mpd[track]);
			if ( (eta_bin == -1) || (pt_bin == -1 )) continue;
			if (n_hits_mpd[track] < Cut_No_Of_hits_min) continue; //n hits in TPC per track cut
			
			//if (TMath::Abs(DCA_x_mpd[track]) >= f_dca[0][pt_bin][eta_bin]->GetParameter(2)*2) continue;
			//if (TMath::Abs(DCA_y_mpd[track]) >= f_dca[1][pt_bin][eta_bin]->GetParameter(2)*2) continue;
			//if (TMath::Abs(DCA_z_mpd[track]) >= f_dca[2][pt_bin][eta_bin]->GetParameter(2)*2) continue;
			TF1 sigma_fit_res_X = *f_pt_fit[0][eta_bin];
			TF1 sigma_fit_res_Y = *f_pt_fit[1][eta_bin];
			TF1 sigma_fit_res_Z = *f_pt_fit[2][eta_bin];
		

			if (TMath::Abs(DCA_x_mpd[track]) >= sigma_fit_res_X(Abs(signed_pt_mpd[track]))*2) continue;
			if (TMath::Abs(DCA_y_mpd[track]) >= sigma_fit_res_Y(Abs(signed_pt_mpd[track]))*2) continue;
			if (TMath::Abs(DCA_z_mpd[track]) >= sigma_fit_res_Z(Abs(signed_pt_mpd[track]))*2) continue;
			//if (mother_ID_mc[id_from_mc_mpd[track]] > -1) continue;
			
			//~ | |    (_)     | |                                                    / _|| |                           | |       
 //~ | |__   _  ___ | |_  ___    __ _  _ __  __ _  _ __ ___   ___    __ _ | |_ | |_  ___  _ __    ___  _   _ | |_  ___ 
 //~ | '_ \ | |/ __|| __|/ _ \  / _` || '__|/ _` || '_ ` _ \ / __|  / _` ||  _|| __|/ _ \| '__|  / __|| | | || __|/ __|
 //~ | | | || |\__ \| |_| (_) || (_| || |  | (_| || | | | | |\__ \ | (_| || |  | |_|  __/| |    | (__ | |_| || |_ \__ \
 //~ |_| |_||_||___/ \__|\___/  \__, ||_|   \__,_||_| |_| |_||___/  \__,_||_|   \__|\___||_|     \___| \__,_| \__||___/
                             //~ __/ |                                                                                 
                            //~ |___/                                                                                  
			
			h_nhits_TPC_after->Fill(n_hits_mpd[track]); //n hits in TPC after eta, nhits, pt and mother ID cuts
			multiplicity_tpc_after++; //multiplicity in TPC after eta, nhits, pt and mother ID cuts
			sort = 0;
			//filling some distributions after track cuts, depending on sort of the particle from MC data
			h_phi_after[cenrality_bin_res][0]->Fill(ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])));
			h_eta_after[cenrality_bin_res][0]->Fill(eta_mpd[track]);
			h_pt_after[cenrality_bin_res][0]->Fill(Abs(signed_pt_mpd[track]));
			h2_pt_vs_eta_after[cenrality_bin_res][0]->Fill(Abs(signed_pt_mpd[track]),eta_mpd[track]);
			h2_pt_vs_phi_after[cenrality_bin_res][0]->Fill(Abs(signed_pt_mpd[track]),ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])));
			h2_phi_vs_eta_after[cenrality_bin_res][0]->Fill(ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])),eta_mpd[track]);
			//~ p_momenta_resolution[0]->Fill(Abs(signed_pt_mpd[track]),Abs(Abs(signed_pt_mpd[track]) - pt_mc[id_from_mc_mpd[track]])
			//~ /Abs(signed_pt_mpd[track]));
				
			switch (PDG_code_mc[id_from_mc_mpd[track]]) {
			case 211:
				sort = 1;
				break;
			case 2212:
				sort = 2;
				break;
			case 321:
				sort = 3;
				break;
			default:
				break;
			}
			
			if ((sort == 1) || (sort == 2) || (sort == 3))
			{
				h_phi_after[cenrality_bin_res][sort]->Fill(ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])));
				h_eta_after[cenrality_bin_res][sort]->Fill(eta_mpd[track]);
				h_pt_after[cenrality_bin_res][sort]->Fill(Abs(signed_pt_mpd[track]));
				h2_pt_vs_eta_after[cenrality_bin_res][sort]->Fill(Abs(signed_pt_mpd[track]),eta_mpd[track]);
				h2_pt_vs_phi_after[cenrality_bin_res][sort]->Fill(Abs(signed_pt_mpd[track]),ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])));
				h2_phi_vs_eta_after[cenrality_bin_res][sort]->Fill(ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track])),eta_mpd[track]);
				//~ p_momenta_resolution[sort]->Fill(Abs(signed_pt_mpd[track]),Abs(Abs(signed_pt_mpd[track]) - pt_mc[id_from_mc_mpd[track]])
				//~ /Abs(signed_pt_mpd[track]));
			}
			
			//~ if (Abs(eta_mpd[track]) > Cut_Eta_Min) //and finally, the resolution gap eta cut, then filling the buffer of EP particles:
			//~ {
				event_plane_buffer[ep_particle_count].Pt = Abs(signed_pt_mpd[track]);
				event_plane_buffer[ep_particle_count].Eta = eta_mpd[track];
				event_plane_buffer[ep_particle_count].Phi = ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track]));
				ep_particle_count++;
			//~ }
		}
		h_multiplicity->Fill(multiplicity_tpc_after); //filling multiplicity TPC after cuts hitsogram, now that the loop over mpd tracks is over
		
		for (int arm = 0; arm < _N_ARM; ++arm)
		{
			h_energy_ZDC_total[arm]->Fill(GetTotalEnergy(ZDC_energy_mpd,arm)); //filling the total ZDC energy for given event
			h2_energy_ZDC_vs_multiplicity[arm]->Fill(GetTotalEnergy(ZDC_energy_mpd,arm),multiplicity_tpc_after);//filling total ZDC energy vs TPC 
			//multiplicity after cuts for a given event
		}
		p_b_vs_energy->Fill(GetTotalEnergy(ZDC_energy_mpd,0)+GetTotalEnergy(ZDC_energy_mpd,1),b_mc);
		h2_energyZDC_L_vs_energy_ZDC_R->Fill(GetTotalEnergy(ZDC_energy_mpd,0) , GetTotalEnergy(ZDC_energy_mpd,1)); //filling ZDC energy L vs R for a given event

		int multiplicity_R = 0 , multiplicity_L = 0;
		for (int ep_particle = 0; ep_particle < ep_particle_count;++ep_particle) //TPC subevent multiplicity after cuts
		{
			if (event_plane_buffer[ep_particle].Eta > 0) multiplicity_R++;
			else if (event_plane_buffer[ep_particle].Eta < 0) multiplicity_L++;
		}
		
		//calculate the EP resolution using TPC if subevents have enough multipllicity
		if ((multiplicity_L >= 4) && (multiplicity_R >= 4)) FillTPC(cenrality_bin_res,event_plane_buffer,ep_particle_count);
		
		//if both ZDCs have some signal then calculate EP resolutions using ZDC
		if ((GetTotalEnergy(ZDC_energy_mpd, 0) != 0)&&(GetTotalEnergy(ZDC_energy_mpd, 1 != 0))) FillZDC(cenrality_bin_res,ZDC_energy_mpd);
	}

	return;
}

   //~ __  _                 
  //~ / _|| |                
 //~ | |_ | |  ___ __      __
 //~ |  _|| | / _ \\ \ /\ / /
 //~ | |  | || (_) |\ V  V / 
 //~ |_|  |_| \___/  \_/\_/  
                         

void MpdCalculator::CalculateFlow(Int_t nevents, TString fitFile)
{
	TFile *f = new TFile(fitFile.Data()); //opening file with fit data
	TF1 *resolution_fit[_N_HARM][_N_HARM][_N_METHOD]; //declaring the fit-function array
	for (int harm = 0; harm < _N_HARM; ++harm)
	{
		for (int _harm = 0; _harm < _N_HARM; ++_harm)
		{
			for (int method = 0; method < _N_METHOD; ++method)
			{
				char name[200];
				sprintf(name,"resolution_fit[%i][%i][%i]",harm,_harm,method);
				resolution_fit[harm][_harm][method] = (TF1*)f->Get(name); //reading fit-functions from file
			}
		}
	}

	if (nevents == 0) nevents = inChain->GetEntries();
	for (Int_t event = 0; event < nevents; ++event) //loop over events
	{
		inChain->GetEntry(event); //reading all branches for a given event
		cout << "FLOW EVENT # "<<event << endl;
		
		//Int_t multiplicity = GetMultiplicityTPC();
		centrality_tpc_mpd = 100 - centrality_tpc_mpd;
		
		int cenrality_bin_flow = GetCentralityBinFlow(centrality_tpc_mpd); //getting b bin of a given event for flow measurements
		//if (cenrality_bin_flow == -1) continue;
		int cenrality_bin_res = GetCentralityBinRes(centrality_tpc_mpd);//getting the resolution bin in whic the event is 
		if (cenrality_bin_res == -1) continue; // just in case...
		//skip the event if the ZDC signal is zero
		if ((GetTotalEnergy(ZDC_energy_mpd, 0) == 0) || (GetTotalEnergy(ZDC_energy_mpd, 1) == 0)) continue;
		
		Double_t phi_EP = ATan2(Sin(phiEP_mc),Cos(phiEP_mc)); //unfold the generated event plane
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//RECONSTRUCTED
		const Float_t p_mass2[]={0.88035,0.2437,0.01948};
		const Int_t pdg_codes[]={2212,321,211}; //0 - PROTON, 1 - KAON, 2 -PION
		Int_t sort = 0;
		for (Long64_t track = 0; track < n_tracks_mpd; ++track)//loop over reconstructed tracks to fill flow buffer
		{	
			//if (id_from_mc_mpd[track] == -1) continue; //equivalent to mother ID cut
			if (n_hits_mpd[track] < Cut_No_Of_hits_min) continue; //n hits in TPC cut
			if (PDG_code_mc[id_from_mc_mpd[track]] != pdg_codes[sort]) continue;
			
			Float_t Pt = Abs(signed_pt_mpd[track]);
			Float_t Eta = eta_mpd[track];
			Float_t Phi = ATan2(Sin(phi_mpd[track]),Cos(phi_mpd[track]));
			//~ Float_t Rapidity = 0.5*TMath::Log((energy_mc[id_from_mc_mpd[track]] + pz_mc[id_from_mc_mpd[track]])
				//~ /(energy_mc[id_from_mc_mpd[track]] - pz_mc[id_from_mc_mpd[track]]));
			Float_t Rapidity = TMath::Log( (TMath::Sqrt(p_mass2[sort]+Pt*Pt*TMath::CosH(Eta)*TMath::CosH(Eta))+Pt*TMath::SinH(Eta))/( TMath::Sqrt(p_mass2[sort]+Pt*Pt)) );
			
			Int_t pt_bin = GetPtBin(Pt);
			Int_t eta_bin = GetEtaBin(Eta);
			if ( (eta_bin == -1) || (pt_bin == -1 )) continue;
			
			//if (TMath::Abs(DCA_x_mpd[track]) >= f_dca[0][pt_bin][eta_bin]->GetParameter(2)*2) continue;
			//if (TMath::Abs(DCA_y_mpd[track]) >= f_dca[1][pt_bin][eta_bin]->GetParameter(2)*2) continue;
			//if (TMath::Abs(DCA_z_mpd[track]) >= f_dca[2][pt_bin][eta_bin]->GetParameter(2)*2) continue;

			TF1 sigma_fit_X = *f_pt_fit[0][eta_bin];
			TF1 sigma_fit_Y = *f_pt_fit[1][eta_bin];
			TF1 sigma_fit_Z = *f_pt_fit[2][eta_bin];
		

			if (TMath::Abs(DCA_x_mpd[track]) >= sigma_fit_X(Pt)*2) continue;
			if (TMath::Abs(DCA_y_mpd[track]) >= sigma_fit_Y(Pt)*2) continue;
			if (TMath::Abs(DCA_z_mpd[track]) >= sigma_fit_Z(Pt)*2) continue;

			
			Double_t Psi_1_FULL_ZDC = GetPsiFullZdc(ZDC_energy_mpd,1); //getting event plane using ZDC for 1st harmonic
		
			int sign;
			if (Eta < 0) sign = -1;
			else sign = 1;
			TF1 res11 = *resolution_fit[0][0][1];
			Float_t mid_bin_cent = (centralityBinsRes[cenrality_bin_res]+centralityBinsRes[cenrality_bin_res+1]) / 2;
			if (res11(mid_bin_cent) != 0)
			{	
				if (cenrality_bin_flow != -1)
				{
					if ((Abs(Eta) > 0.2) && (Abs(Eta) < Cut_Eta_Max)){
						p_flow_wrt_full_vs_pt_divided[cenrality_bin_flow][0][0][1]->Fill(Pt,sign*Cos(Psi_1_FULL_ZDC - Phi)/res11(mid_bin_cent));
					}
					//if (Pt > Cut_Pt_Min)
					p_flow_wrt_full_vs_eta_divided[cenrality_bin_flow][0][0][1]->Fill(Eta,Cos(Psi_1_FULL_ZDC - Phi)/res11(mid_bin_cent));
					//if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min))
					p_flow_wrt_full_vs_rapidity_divided[cenrality_bin_flow][0][0][1]->Fill(Rapidity,Cos(Psi_1_FULL_ZDC - Phi)/res11(mid_bin_cent));
				}
				//cout << "filling with centralityBinsRes[cenrality_bin_res] = " << centralityBinsRes[cenrality_bin_res] << " Cos(Psi_1_FULL_ZDC - Phi)/res11(mid_bin_cent) = " << Cos(Psi_1_FULL_ZDC - Phi)/res11(mid_bin_cent) << endl;
				//if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min))
				p_flow_wrt_full_vs_centrality_divided[0][0][1]->Fill(centralityBinsRes[cenrality_bin_res] + 0.1,sign*Cos(Psi_1_FULL_ZDC - Phi)/res11(mid_bin_cent));
			}
			TF1 res21 = *resolution_fit[1][0][1];
			if (res21(mid_bin_cent) != 0)
			{
				if (cenrality_bin_flow != -1)
				{
					//if (Abs(Eta) < Cut_Eta_Max)
					p_flow_wrt_full_vs_pt_divided[cenrality_bin_flow][1][0][1]->Fill(Pt,Cos(2.*(Psi_1_FULL_ZDC - Phi))/res21(mid_bin_cent));
					//if (Pt > Cut_Pt_Min)	
					p_flow_wrt_full_vs_eta_divided[cenrality_bin_flow][1][0][1]->Fill(Eta,Cos(2.*(Psi_1_FULL_ZDC - Phi))/res21(mid_bin_cent));
					//if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min))
					p_flow_wrt_full_vs_rapidity_divided[cenrality_bin_flow][1][0][1]->Fill(Rapidity,Cos(2.*(Psi_1_FULL_ZDC - Phi))/res21(mid_bin_cent));
				}
				//if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min))
				p_flow_wrt_full_vs_centrality_divided[1][0][1]->Fill(centralityBinsRes[cenrality_bin_res] + 0.1,Cos(2.*(Psi_1_FULL_ZDC - Phi))/res21(mid_bin_cent));
			}
		
		}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//GENERATED
		for (Long64_t track = 0; track < n_tracks_mc; ++track)
		{
			if (PDG_code_mc[track] != pdg_codes[sort]) continue;
			
			Float_t Pt = Abs(pt_mc[track]);
			Float_t Eta = eta_mc[track];
			Float_t Phi = ATan2(py_mc[track],px_mc[track]);
			Float_t Rapidity = .5*TMath::Log((energy_mc[track] + pz_mc[track])/(energy_mc[track] - pz_mc[track]));
			
			Int_t pt_bin = GetPtBin(Pt);
			Int_t eta_bin = GetEtaBin(Eta);
			if ( (eta_bin == -1) || (pt_bin == -1 )) continue;
			if (mother_ID_mc[track] > -1 ) continue;
			
			int sign;
			if (Eta < 0) sign = -1;
			else sign = 1;
			
			//if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min)) 
			//{
				p_flow_wrt_RP_vs_centrality[0]->Fill(centralityBinsRes[cenrality_bin_res] + 0.1, sign*Cos(phi_EP - Phi));
				//cout << "Filling with centralityBinsRes[cenrality_bin_res] + 0.1 = " << centralityBinsRes[cenrality_bin_res] + 0.1 << "Cos(phi_EP - Phi) = " << Cos(phi_EP - Phi) << endl;
				p_flow_wrt_RP_vs_centrality[1]->Fill(centralityBinsRes[cenrality_bin_res] + 0.1, Cos(2.*(phi_EP - Phi)));
			//}
			
			if (cenrality_bin_flow != -1)
			{
				if ((Abs(Eta) > 0.2) && (Abs(Eta) < Cut_Eta_Max))
				p_flow_wrt_RP_vs_pt[cenrality_bin_flow][0]->Fill(Pt , sign*Cos(phi_EP - Phi));
				//if (Abs(Eta) < Cut_Eta_Max)	
				p_flow_wrt_RP_vs_pt[cenrality_bin_flow][1]->Fill(Pt , Cos(2.*(phi_EP-Phi)));
				
				//if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min))
				p_flow_wrt_RP_vs_rapidity[cenrality_bin_flow][0]->Fill(Rapidity , Cos(phi_EP	- Phi));
				//if ((Abs(Eta) < Cut_Eta_Max) && (Pt > Cut_Pt_Min))	
				p_flow_wrt_RP_vs_rapidity[cenrality_bin_flow][1]->Fill(Rapidity , Cos(2.*(phi_EP - Phi)));
				//if (Pt > Cut_Pt_Min)	
				p_flow_wrt_RP_vs_eta[cenrality_bin_flow][0]->Fill(Eta , Cos(phi_EP - Phi));
				//if (Pt > Cut_Pt_Min)	
				p_flow_wrt_RP_vs_eta[cenrality_bin_flow][1]->Fill(Eta , Cos(2.*(phi_EP - Phi)));
			}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}
}

void MpdCalculator::GetQsTpc(EPParticle* event_plane_buffer, Int_t buffer_size, Int_t sign, Int_t harm, Double_t &Qx, Double_t &Qy)
{
	Double_t Qcos=0., Qsin=0. , total_momenta = 0;

	total_momenta = GetTotalMomenta(event_plane_buffer,buffer_size, sign);

	Int_t w_sign = 0;
	if (harm == 2) w_sign = 1;
	else w_sign = sign;

	for (int ep_particle = 0; ep_particle < buffer_size; ++ep_particle)
	{
		if (!(event_plane_buffer[ep_particle].Eta*sign > 0)) continue;
		Qcos += w_sign*event_plane_buffer[ep_particle].Pt / total_momenta * Cos(harm*event_plane_buffer[ep_particle].Phi);
		Qsin += w_sign*event_plane_buffer[ep_particle].Pt / total_momenta * Sin(harm*event_plane_buffer[ep_particle].Phi);
	}

	Qx = Qcos;
	Qy = Qsin;
}

void MpdCalculator::GetQsZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t harm, Double_t &Qx, Double_t &Qy)
{
	Double_t *phi_angle_of_modules = GetAngles();
	Double_t Qcos = 0 , Qsin = 0;
	Double_t w_sign = 0;
	Double_t total_energy = GetTotalEnergy(zdc_energy,zdc_ID);

	if (harm == 2) w_sign = 1;
	else if (harm == 1)
	{
		if (zdc_ID == 0) w_sign = 1;
		else if (zdc_ID == 1) w_sign = -1;
	}


	for (int module = _N_MODULES_TOTAL/2 * zdc_ID; module < _N_MODULES_TOTAL/2 * (zdc_ID + 1); ++module)
	{
		if ((module==22) || (module==67)) continue;
//		if ((i ==15)||(i ==21)||(i==23)||(i ==29)||(i==60)||
//				(i==66)||(i==68)||(i==74)) continue;
		Qcos += w_sign*zdc_energy[module] / total_energy * Cos(harm*phi_angle_of_modules[module]);
        Qsin += w_sign*zdc_energy[module] / total_energy * Sin(harm*phi_angle_of_modules[module]);
	}

	Qx = Qcos; Qy = Qsin;
}

Double_t MpdCalculator::GetPsiHalfTpc(EPParticle* event_plane_buffer, Int_t buffer_size, Int_t sign, Int_t harm, Double_t &Qx, Double_t &Qy) //if sign == 1 then its positive pseudorapidity and all weights are positive,
{
	GetQsTpc(event_plane_buffer,buffer_size,sign,harm,Qx,Qy);
	Double_t PsiEP = (1/(Double_t)harm) * ATan2(Qx,Qy);
	return PsiEP;
}

Double_t MpdCalculator::GetPsiHalfZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t n, Double_t &qx, Double_t &qy)
{
	Double_t Qcos=0., Qsin=0.;
	GetQsZdc(zdc_energy,zdc_ID,n,Qcos,Qsin);

	qx = Qcos;
	qy = Qsin;
	Double_t PsiEP = (1/(Double_t)n) * ATan2(Qsin,Qcos);
	return PsiEP;
}

Double_t MpdCalculator::GetPsiFullZdc(Float_t* zdc_energy, Int_t n)
{
	Double_t QcosR=0., QsinR=0.;
	GetQsZdc(zdc_energy,0,n,QcosR, QsinR);

	Double_t QcosL=0., QsinL=0.;
	GetQsZdc(zdc_energy,1,n,QcosL,QsinL);

	Double_t psiEP = ATan2(QsinR + QsinL,QcosR + QcosL)/(Double_t) n; // (-pi/n,pi/n]
    return psiEP;
}

Double_t MpdCalculator::GetTotalMomenta(EPParticle* event_plane_buffer, Int_t buffer_size, Int_t sign)
{
	Double_t total_momenta = 0.;
	for (int ep_particle = 0; ep_particle < buffer_size; ++ep_particle)
	{
		if (!(event_plane_buffer[ep_particle].Eta*sign > 0)) continue;
		total_momenta += event_plane_buffer[ep_particle].Pt;
	}
	return total_momenta;
}

Double_t MpdCalculator::GetTotalEnergy(Float_t* zdc_energy, Int_t zdc_ID)
{
	Double_t total_energy = 0.;
	for (int i = _N_MODULES_TOTAL/2 * zdc_ID; i < _N_MODULES_TOTAL/2 * (zdc_ID + 1); ++i)
	{
		if ((i==22) || (i==67)) continue;
//		if ((i ==15)||(i ==21)||(i==23)||(i ==29)||(i==60)||
//				(i==66)||(i==68)||(i==74)) continue;
		total_energy += zdc_energy[i];
	}
	return total_energy;
}

Double_t MpdCalculator::GetPsiFullTpc(EPParticle* event_plane_buffer, Int_t buffer_size, Int_t harm)
{
	Double_t QcosR=0., QsinR=0.;
	GetQsTpc(event_plane_buffer,buffer_size,1,harm,QcosR,QsinR);

	Double_t QcosL=0., QsinL=0.;
	GetQsTpc(event_plane_buffer,buffer_size,-1,harm,QcosL,QsinL);

	Double_t psiEP = ATan2(QsinR + QsinL,QcosR + QcosL)/(Double_t) harm; // (-pi/n,pi/n]
	return psiEP;
}

                  //~ _  _        
                 //~ (_)| |       
 //~ __      __ _ __  _ | |_  ___ 
 //~ \ \ /\ / /| '__|| || __|/ _ \
  //~ \ V  V / | |   | || |_|  __/
   //~ \_/\_/  |_|   |_| \__|\___|

void MpdCalculator::Write()
{
	outFile->cd();
	
	h_nhits_TPC->Write();
	h_nhits_TPC_after->Write();
	h_impact_parameter->Write();
	h_multiplicity->Write();
	h_multiplicity_before->Write();
	h2_energyZDC_L_vs_energy_ZDC_R->Write();
	p_b_vs_multiplicity->Write();
	p_b_vs_energy->Write();
	h2_b_vs_centrality->Write();
	
	for (Int_t dim = 0; dim < Ndim; dim++)
	{	
		for (Int_t pt_bin = 0; pt_bin < NptBins; pt_bin++)
		{
			for (Int_t eta_bin = 0; eta_bin < NetaBins ; eta_bin++ )
			{
				h_DCA_all[dim][pt_bin][eta_bin]->Write();
				h_DCA_primary[dim][pt_bin][eta_bin]->Write();
				h_DCA_secondary[dim][pt_bin][eta_bin]->Write();
			}
		}
	}
	
	for (int sort = 0; sort < _N_SORTS; ++sort)
	{
		p_momenta_resolution[sort]->Write();
	}

	for (int harm = 0; harm < _N_HARM; ++harm)
	{
		p_flow_wrt_RP_vs_centrality[harm]->Write();
		for (int centralityBin = 0; centralityBin < NcentralityBinsFlow; ++centralityBin)
		{
			p_flow_wrt_RP_vs_pt[centralityBin][harm]->Write();
			p_flow_wrt_RP_vs_eta[centralityBin][harm]->Write();
			p_flow_wrt_RP_vs_rapidity[centralityBin][harm]->Write();
		}

		for (int method = 0; method < _N_METHOD; ++method)
		{
			for (int _harm = 0; _harm < _N_HARM; ++_harm)
			{
				p_flow_wrt_full_vs_centrality[harm][_harm][method]->Write();
				p_flow_wrt_full_vs_centrality_divided[harm][_harm][method]->Write();
				p_Res2Psi_vs_b[harm][_harm][method]->Write();
				p_true_Res_vs_b[harm][_harm][method]->Write();
				for (int centralityBin = 0; centralityBin < NcentralityBinsFlow; ++centralityBin)
				{
					p_flow_wrt_full_vs_pt[centralityBin][harm][_harm][method]->Write();
					p_flow_wrt_full_vs_eta[centralityBin][harm][_harm][method]->Write();
					p_flow_wrt_full_vs_rapidity[centralityBin][harm][_harm][method]->Write();
					p_flow_wrt_full_vs_pt_divided[centralityBin][harm][_harm][method]->Write();
					p_flow_wrt_full_vs_eta_divided[centralityBin][harm][_harm][method]->Write();
					p_flow_wrt_full_vs_rapidity_divided[centralityBin][harm][_harm][method]->Write();
				}
			}
		}
	}

	for (int arm = 0; arm < _N_ARM; ++arm)
	{
		h_energy_ZDC_total[arm]->Write();
		h2_energy_ZDC_vs_multiplicity[arm]->Write();

		for (int harm = 0; harm < _N_HARM; ++harm)
		{
			for (int method = 0; method < _N_METHOD; ++method)
			{
				p_qx_vs_b[arm][harm][method]->Write();
				p_qy_vs_b[arm][harm][method]->Write();
				for (int _harm = 0; _harm < _N_HARM; ++_harm)
				{
					p_true_Res_half_vs_b[arm][harm][_harm][method]->Write();
				}
			}
		}
	}

	for (int centralityBin = 0; centralityBin < NcentralityBinsRes; ++centralityBin)
	{
		for (int sort = 0; sort < _N_SORTS; ++sort)
		{
			h_phi[centralityBin][sort]->Write();
			h_pt[centralityBin][sort]->Write();
			h_eta[centralityBin][sort]->Write();
			h2_pt_vs_eta[centralityBin][sort]->Write();
			h2_pt_vs_phi[centralityBin][sort]->Write();
			h2_phi_vs_eta[centralityBin][sort]->Write();
			h_phi_after[centralityBin][sort]->Write();
			h_pt_after[centralityBin][sort]->Write();
			h_eta_after[centralityBin][sort]->Write();
			h2_pt_vs_eta_after[centralityBin][sort]->Write();
			h2_pt_vs_phi_after[centralityBin][sort]->Write();
			h2_phi_vs_eta_after[centralityBin][sort]->Write();
			
			h_phi_mc[centralityBin][sort]->Write();
			h_pt_mc[centralityBin][sort]->Write();
			h_eta_mc[centralityBin][sort]->Write();
			h_phi_mc_after[centralityBin][sort]->Write();
			h_pt_mc_after[centralityBin][sort]->Write();
			h_eta_mc_after[centralityBin][sort]->Write();
		}
	}
	outFile->Close();
}

Int_t MpdCalculator::GetCentralityBinRes(Int_t centrality)
{
	int centrality_bin = -1;
	for (int c_bin = 0; c_bin < NcentralityBinsRes; ++c_bin){
	if ((centrality > centralityBinsRes[c_bin]) && (centrality <= centralityBinsRes[c_bin+1])) 
		centrality_bin = c_bin;
	}
	return centrality_bin;
}

Int_t MpdCalculator::GetCentralityBinFlow(Int_t centrality)
{
	int centrality_bin = -1;
	for (int c_bin = 0; c_bin < NcentralityBinsFlow; ++c_bin)
	if ((centrality > centralityBinsFlow[c_bin]) && (centrality <= centralityBinsFlow[c_bin+1])) 
		centrality_bin = c_bin;
	if (centrality_bin == 1) centrality_bin = -1;
	else if (centrality_bin == 2) centrality_bin = 1; ////////////////////PODGON!!!!
	return centrality_bin;
}

Int_t MpdCalculator::GetMultiplicityTPC() //should be called in loop over events
{
	Int_t multiplicity = 0;
	for (Long64_t track = 0; track < n_tracks_mpd; ++track)//loop over mpdtracks, cut on mother ID will be everywhere
	{
		if (id_from_mc_mpd[track] == -1) continue; //equivalent to mother id cut	
		multiplicity++; //multiplicity in TPC before eta, nhits and pt cuts, but after mother ID cut
	}
	return multiplicity;
}

int MpdCalculator::GetPtBin(Float_t pt)
{
	int pt_bin = -1;
	for (int i = 0; i < NptBins; ++i)
	if ((pt > ptBins[i]) && (pt <= ptBins[i + 1])) 
		pt_bin = i;
	return pt_bin;
}

int MpdCalculator::GetEtaBin(Float_t eta)
{
	int eta_bin = -1;
	for (int i = 0; i < NetaBins; ++i)
	if ((eta > etaBins[i]) && (eta <= etaBins[i + 1])) 
		eta_bin = i;
	return eta_bin;
}
