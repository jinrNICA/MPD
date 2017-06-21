
#include <cstdio>

#include <TProfile.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TF1.h>

const int _N_HARM = 2;
const int _N_METHOD = 2;

const float centralityBinsFlow[] = {10.,20.,40.,50};
const int NcentralityBinsFlow = 3;

const float centralityBinsRes[] = {0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100};
const int NcentralityBinsRes = 20;

const float ptBins[] = {0.,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.};
const int NptBins = 12;

const float etaBins[] = {-1.5,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.5};
const int NetaBins = 14;

const float rapidityBins[] = {-2., -1.8, -1.6, -1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.};
const int NrapidityBins = 20;

const TString methods_names[_N_METHOD] = {TString("TPC") , TString("FHCal")};

using TMath::Sqrt;

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
