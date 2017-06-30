#include "get_res.h"

void get_res(TString inFileName , TString outFileName)
{
	
	gSystem->Load("./utility_C.so");
	gSystem->Load("./MpdCalculator_C.so");
		
	TFile *inFile = new TFile(inFileName.Data());
	TFile *outFile = new TFile(outFileName.Data(),"RECREATE");
	
	TProfile *p_Res2Psi_vs_b[_N_HARM][_N_HARM][_N_METHOD];
	TProfile *p_ResPsi_vs_b[_N_HARM][_N_HARM][_N_METHOD];
	TF1 *resolution_fit[_N_HARM][_N_HARM][_N_METHOD];
	
	for (int harm = 0; harm < _N_HARM; ++harm)
	{
		for (int _harm = 0; _harm < _N_HARM; ++_harm)
		{
			for (int method = 0; method < _N_METHOD; ++method)
			{
				
				char name[200];
				char title[200];
				
				sprintf(name,"p_Res2Psi_vs_b[%i][%i][%i]",harm,_harm,method);
				p_Res2Psi_vs_b[harm][_harm][method] = (TProfile*)inFile->Get(name);
				
				sprintf(name,"p_ResPsi_vs_b[%i][%i][%i]",harm,_harm,method);
				sprintf(title,"Res(#Psi_{%i,%s}) for v_{%i};b,fm;",_harm+1,methods_names[method].Data(),harm+1);
				p_ResPsi_vs_b[harm][_harm][method] = new TProfile(name,title,NcentralityBinsRes,centralityBinsRes);
			}
		}
	}
	
	outFile->cd();
	for (int harm = 0; harm < _N_HARM; ++harm)
	{
		for (int _harm = 0; _harm < _N_HARM; ++_harm)
		{
			for (int method = 0; method < _N_METHOD; ++method)
			{
				for (int centralitybin = 0; centralitybin < NcentralityBinsRes; ++centralitybin)
				{
					Double_t res2 = p_Res2Psi_vs_b[harm][_harm][method]->GetBinContent(centralitybin + 1);
					if (res2 < 0)res2 = 0;
					Double_t chi_half = Chi(TMath::Sqrt(res2),harm+1);
					p_ResPsi_vs_b[harm][_harm][method]->Fill(centralityBinsRes[centralitybin] + 0.1,ResEventPlane(TMath::Sqrt(2.) * chi_half,harm+1));
				}
				
				char name[200];
				sprintf(name,"resolution_fit[%i][%i][%i]",harm,_harm,method);
				resolution_fit[harm][_harm][method] = new TF1(name,"pol8",0.01,100);
				p_ResPsi_vs_b[harm][_harm][method]->Fit(name,"W");
				p_ResPsi_vs_b[harm][_harm][method]->Write();
				resolution_fit[harm][_harm][method]->Write();
			}
		}
	}
	
	outFile->Close();
	inFile->Close();
}
