//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Class MpdCalculator contains the required methods for the flow analysis. There're two main procedures:   			            //
// 	1. Calculate resolution factor correction (Requires to run get_res.cxx afterwards!!!)                      		            //
//  2. Calculate azimuthal flow v1 and v2 using Event Plane method.                                                             //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//USAGE:    can be executed via standard root																		            //
//root -l               																							            //
//  gSystem->Load("libMathMore")                                                                                                //
//  .L ../Utilities/BinningData.cxx+                    															            //
//  .L ../Utilities/utility.cxx+																					            //
//  .L MpdCalculator.cxx+   																						            //
///////CALCULATE RESOLUTION FACTOR CORRECTION/////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                              //
//  MpdCalculator mpd = MpdCalculator(inFileName,outFileName,dcaFileName)                                         	            //
//	mpd.CalculateResolutions(0) 																					            //
//  mpd.Write()                                                                                                                 //
//  !!!!!!!!!!!!RUN get_res.cxx!!!!!!!!!!!                                                                                      //
//ARGUMENTS:																										            //
// inFileName - input standard root file with TTree cbmsim_reduced for the further analysis (output from reducedTreeCreator)    //
// outFileName - output standard root file with TProfiles and histograms containing data for resolution correction factor       //
// dcaFileName - Second iteration of the dca fitting containing sigma_pt_fit TF1* functions (output MakeFitDCA(...))            //
//                                                                                                                              //
///////CALCULATE AZIMUTHAL FLOW///////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                              //
//  MpdCalculator mpd = MpdCalculator(inFileName,outFileName,dcaFileName)                                         	            //
//	mpd.CalculateFlow(0, resFitFile.Data())			    															            //
//ARGUMENTS:																										            //
// inFileName - input standard root file with TTree cbmsim_reduced for the further analysis (output from reducedTreeCreator)    //
// outFileName - output standard root file with TProfiles and histograms containing data with azimuthal flow (v1 and v2)        //
// dcaFileName - Second iteration of the dca fitting containing sigma_pt_fit TF1* functions (output MakeFitDCA(...))            //
// resFitFile - input standard root file with TH1* histograms and TF1* fitted function of resolution (output from get_res(...)) //
//                                                                                                                              //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MPD_CALCULATOR_H

#define MPD_CALCULATOR_H

//Including necessary classes & constants
#include "../Utilities/utility.cxx"

class MpdCalculator
{
private:

    //Used varables
	TString inFileName;
	TString outFileName;

	TChain *inChain;
	TFile *outFile;
	
	TFile  *dcaFile;
    //Histograms
    TH1F	 *h_pt[NcentralityBinsRes][_N_SORTS], *h_eta[NcentralityBinsRes][_N_SORTS], *h_phi[NcentralityBinsRes][_N_SORTS];//
    TH1F	 *h_pt_mc[NcentralityBinsRes][_N_SORTS], *h_eta_mc[NcentralityBinsRes][_N_SORTS], *h_phi_mc[NcentralityBinsRes][_N_SORTS];//
    TH1F	 *h_pt_after[NcentralityBinsRes][_N_SORTS], *h_eta_after[NcentralityBinsRes][_N_SORTS], *h_phi_after[NcentralityBinsRes][_N_SORTS];//
    TH1F	 *h_pt_mc_after[NcentralityBinsRes][_N_SORTS], *h_eta_mc_after[NcentralityBinsRes][_N_SORTS], *h_phi_mc_after[NcentralityBinsRes][_N_SORTS];//
    TH2F	 *h2_pt_vs_eta[NcentralityBinsRes][_N_SORTS], *h2_pt_vs_phi[NcentralityBinsRes][_N_SORTS], *h2_phi_vs_eta[NcentralityBinsRes][_N_SORTS];//
    TH2F	 *h2_pt_vs_eta_after[NcentralityBinsRes][_N_SORTS], *h2_pt_vs_phi_after[NcentralityBinsRes][_N_SORTS], *h2_phi_vs_eta_after[NcentralityBinsRes][_N_SORTS];//
    TH1F 	 *h_nhits_TPC, *h_nhits_TPC_after, *h_multiplicity , *h_multiplicity_before, *h_impact_parameter;//
    TH2F	 *h2_energy_ZDC_vs_multiplicity[_N_ARM];//
    TH2F	 *h2_energyZDC_L_vs_energy_ZDC_R;//
    TH2F     *h2_b_vs_centrality;
    TH1F	 *h_energy_ZDC_total[_N_ARM];//
    TH1F     *h_DCA_primary[Ndim][NptBins][NetaBins], *h_DCA_secondary[Ndim][NptBins][NetaBins], *h_DCA_all[Ndim][NptBins][NetaBins];
    
    //TF1 	 *f_dca[Ndim][NptBins][NetaBins];
    TF1         *f_pt_fit[Ndim][NetaBins];
    //Profiles
    TProfile *p_Res2Psi_vs_b[_N_HARM][_N_HARM][_N_METHOD];//
    TProfile *p_true_Res_vs_b[_N_HARM][_N_HARM][_N_METHOD];//
    TProfile *p_true_Res_half_vs_b[_N_ARM][_N_HARM][_N_HARM][_N_METHOD];//
    TProfile *p_qx_vs_b[_N_ARM][_N_HARM][_N_METHOD];//
    TProfile *p_qy_vs_b[_N_ARM][_N_HARM][_N_METHOD];//
    TProfile *p_flow_wrt_full_vs_centrality[_N_HARM][_N_HARM][_N_METHOD], *p_flow_wrt_full_vs_pt[NcentralityBinsFlow][_N_HARM][_N_HARM][_N_METHOD], *p_flow_wrt_full_vs_eta[NcentralityBinsFlow][_N_HARM][_N_HARM][_N_METHOD], *p_flow_wrt_full_vs_rapidity[NcentralityBinsFlow][_N_HARM][_N_HARM][_N_METHOD];
    //////////         RESULTED AZIMUTHAL FLOW            //////////
    //Reconstructed (realistic values) 
    TProfile *p_flow_wrt_full_vs_centrality_divided[_N_HARM][_N_HARM][_N_METHOD], *p_flow_wrt_full_vs_pt_divided[NcentralityBinsFlow][_N_HARM][_N_HARM][_N_METHOD], *p_flow_wrt_full_vs_eta_divided[NcentralityBinsFlow][_N_HARM][_N_HARM][_N_METHOD], *p_flow_wrt_full_vs_rapidity_divided[NcentralityBinsFlow][_N_HARM][_N_HARM][_N_METHOD];
    //Generated (idealisic values)
    TProfile *p_flow_wrt_RP_vs_centrality[_N_HARM], *p_flow_wrt_RP_vs_pt[NcentralityBinsFlow][_N_HARM], *p_flow_wrt_RP_vs_eta[NcentralityBinsFlow][_N_HARM], *p_flow_wrt_RP_vs_rapidity[NcentralityBinsFlow][_N_HARM];
    ////////////////////////////////////////////////////////////////
    TProfile *p_momenta_resolution[_N_SORTS];
    TProfile *p_b_vs_multiplicity, *p_b_vs_energy;

	Float_t b_mc;
	Int_t centrality_tpc_mpd;
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
	Float_t DCA_x_mpd[_MAX_TRACKS];
	Float_t DCA_y_mpd[_MAX_TRACKS];
	Float_t DCA_z_mpd[_MAX_TRACKS];
    BinningData* bins;

    //Calculate Event Plane  in TPC subevent (half TPC) harm = 1,2
	Double_t GetPsiHalfTpc(EPParticle* event_plane_buffer, Int_t buffer_size, Int_t sign, Int_t harm, Double_t &Qx, Double_t &Qy); //if sign == 1 then its positive pseudorapidity and all weights are positive
	//Calculate Event Plane  in FHCla subevent (half FHCal) harm = 1,2
    Double_t GetPsiHalfZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t n, Double_t &qx, Double_t &qy);
    //Calculate Event plane in all TPC harm = 1,2
	Double_t GetPsiFullTpc(EPParticle* event_plane_buffer, Int_t buffer_size, Int_t harm); //harm = 
    //Calculate Event plane in all FHCal n = 1,2
	Double_t GetPsiFullZdc(Float_t* zdc_energy, Int_t n);
    //Normalize factor for Q vector calculated in TPC
	Double_t GetTotalMomenta(EPParticle* event_plane_buffer, Int_t buffer_size, Int_t sign);
    //Calculate Q vector in TPC harm = 1,2 
	void GetQsTpc(EPParticle* event_plane_buffer, Int_t buffer_size, Int_t sign, Int_t harm, Double_t &Qx, Double_t &Qy);
    //Calculate Q vector in FHCal harm = 1,2
	void GetQsZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t harm, Double_t &Qx, Double_t &Qy);
    //Normalize factor for Q vector calculated in FHCal
	Double_t GetTotalEnergy(Float_t* zdc_energy, Int_t zdc_ID);
    //Cuts for flow analysis
	bool FlowCut(Float_t eta, Float_t pt, Int_t pdg, Int_t goal);

public:

    //Constructor of the MpdCalculator class
	MpdCalculator(TString inFileName, TString outFileName, TString dcaFileName);
    //Fill and calculate Event plane related values from TPC detector
	void FillTPC(Int_t centrality_bin, EPParticle *event_plane_buffer, Int_t ep_particle_count);
    //Fill and calculate Event plane related values from FHCal detector
	void FillZDC(Int_t centrality_bin, Float_t *ZDC_energy_mpd);
    ////////////////////////////////EXTERNAL USE/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Calculate values for get_res(...) for resolution correction factor
	void CalculateResolutions(Int_t nevents = 0);
    //Calculate azimuthal flow with resolution correction factor for reconstructed (close to experimental) and generated (idealisitc values from the model - no detector effects)
	void CalculateFlow(Int_t nevents, TString fitFile);
    ///////////////////////////////EXTERNAL USE//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Write result in output standard root file outFileName
	void Write();
    //Centrality binning for resolution correction factor
	Int_t GetCentralityBinRes(Int_t multiplicity);
    //Centrality binning for azimuthal flow (v1 and v2)
	Int_t GetCentralityBinFlow(Int_t multiplicity);
    //Additional method for work with multiplicity distribution
	Int_t GetMultiplicityTPC();
};

#endif
