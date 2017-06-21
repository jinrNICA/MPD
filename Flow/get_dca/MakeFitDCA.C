void MakeFitDCA(TString inFileName, TString outFileName)
{
	
	const float ptBins[] = {0.,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.};
	const int NptBins = 12;

	const float etaBins[] = {-1.5,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.5};
	const int NetaBins = 14;

	const int Ndim = 3;

	const float par[] = {0.,0.,0.,0.,0.};

	TF1*   func[Ndim][NptBins][NetaBins];
	TFile* inFile = new TFile(inFileName.Data(),"read");

	TH2F* h_integral[Ndim];
	TH2F* h_mean[Ndim];
	TH2F* h_sigma[Ndim]; 
	TF1*  sigma_pt_fit[Ndim][NetaBins];
	TH1F* h_sigma_X[Ndim][NetaBins];

	for (Int_t i_dim=0;i_dim<Ndim;i_dim++){
		h_integral[i_dim] = new TH2F(Form("h_integral%i",i_dim),Form("Integral of the fit function %i;p_{T};#eta;",i_dim),NptBins,ptBins,NetaBins,etaBins);
		h_mean[i_dim] 	 = new TH2F(Form("h_mean%i",i_dim)   ,Form("Mean of the fit function %i;p_{T};#eta;",i_dim)   ,NptBins,ptBins,NetaBins,etaBins);
		h_sigma[i_dim] 	 = new TH2F(Form("h_sigma%i",i_dim)   ,Form("Sigma of the fit function %i;p_{T};#eta;",i_dim)   ,NptBins,ptBins,NetaBins,etaBins);
	}

	for (Int_t i_dim=0;i_dim<Ndim;i_dim++){
		for (Int_t i_eta=0;i_eta<NetaBins;i_eta++){
			sigma_pt_fit[i_dim][i_eta] = new TF1(Form("sigma_pt_fit%i%i",i_dim,i_eta),"pol 7",0.,3.);
			//sigma_pt_fit[i_dim][i_eta]->SetParameters(par);
		}
	}

	for (Int_t i_dim=0;i_dim<Ndim;i_dim++){
		for (Int_t i_pt=0; i_pt<NptBins;i_pt++){
			for (Int_t i_eta=0;i_eta<NetaBins;i_eta++){
				cout << Form("dca_fit[%i][%i][%i]",i_dim,i_pt,i_eta) << endl;
				
				func[i_dim][i_pt][i_eta] = (TF1*) inFile->Get(Form("dca_fit[%i][%i][%i]",i_dim,i_pt,i_eta));

				h_integral[i_dim]->SetBinContent(i_pt+1,i_eta+1,func[i_dim][i_pt][i_eta]->GetParameter(0));
				h_integral[i_dim]->SetBinError(  i_pt+1,i_eta+1,func[i_dim][i_pt][i_eta]->GetParError(0));
				h_mean[i_dim]   ->SetBinContent(i_pt+1,i_eta+1,func[i_dim][i_pt][i_eta]->GetParameter(1));
				h_mean[i_dim]   ->SetBinError(  i_pt+1,i_eta+1,func[i_dim][i_pt][i_eta]->GetParError(1));
				h_sigma[i_dim]   ->SetBinContent(i_pt+1,i_eta+1,func[i_dim][i_pt][i_eta]->GetParameter(2));
				h_sigma[i_dim]   ->SetBinError(  i_pt+1,i_eta+1,func[i_dim][i_pt][i_eta]->GetParError(2));
			}
		}
	}

	for (Int_t i_dim=0;i_dim<Ndim;i_dim++){
		for (Int_t i_eta=0;i_eta<NetaBins;i_eta++){
			h_sigma_X[i_dim][i_eta] = (TH1F*) h_sigma[i_dim]->ProjectionX(Form("h_sigma_X%i%i",i_dim,i_eta),i_eta+1,i_eta+2,"e");
		}
	}

	for (Int_t i_dim=0;i_dim<Ndim;i_dim++){
		for (Int_t i_eta=0;i_eta<NetaBins;i_eta++){
			h_sigma_X[i_dim][i_eta]->Fit(sigma_pt_fit[i_dim][i_eta],"R");
		}
	}


	TFile* outFile = new TFile(outFileName.Data(),"recreate");
	outFile->cd();
	for (Int_t i_dim=0;i_dim<Ndim;i_dim++){
		h_integral[i_dim]->Write();
		h_mean[i_dim]   ->Write();
		h_sigma[i_dim]   ->Write();
		for (Int_t i_eta=0; i_eta<NetaBins; i_eta++){
			h_sigma_X[i_dim][i_eta]->Write();
			sigma_pt_fit[i_dim][i_eta]->Write();
		}
	}
}
