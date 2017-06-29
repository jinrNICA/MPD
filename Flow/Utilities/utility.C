#include "utility.h"

FlowParticle::FlowParticle()
{
	Eta = 0. , Pt = 0., Phi = 0., Rapidity = 0.;
};

FlowParticle::FlowParticle(double Eta, double Pt, double Phi, double Rapidity)
{
	this->Eta = Eta, this->Pt = Pt, this->Phi = Phi, this->Rapidity = Rapidity;
};

EPParticle::EPParticle()
{
	Eta = 0., Pt = 0., Phi = 0.;
}

EPParticle::EPParticle(double Eta, double Pt, double Phi)
{
	this->Eta = Eta, this->Pt = Pt, this->Phi = Phi;
}

Double_t ResEventPlane(Double_t chi, Int_t harm) //harm = 1 or 2 for our case
{
  Double_t con = TMath::Sqrt(TMath::Pi()/2)/2 ;
  Double_t arg = chi * chi / 4.;
  if (harm == 1) { Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg)); return res; }
  else if (harm == 2) { Double_t res = con * chi * exp(-arg) * (ROOT::Math::cyl_bessel_i(0.5,arg) + ROOT::Math::cyl_bessel_i(1.5,arg)); return res; }
  else cout << "Wrong harmonic in resolution extrapolation! " <<endl;
}

Double_t Chi(Double_t res, Int_t harm) //harm = 1 or 2 for our case
{

  Double_t chi   = 2.0;
  Double_t delta = 1.0;
  for(int i = 0; i < 50; i++)
  {
   if(ResEventPlane(chi, harm) < res) { chi = chi + delta ; }
   else                         { chi = chi - delta ; }
   delta = delta / 2.;
  }

  return chi ;
}

Double_t* GetAngles()
{
	Double_t *phi_angle_of_module = new Double_t[_N_MODULES_TOTAL];

	for (int i = 0; i < _N_ARM; ++i)
	{
		Int_t x_axis_switch;
		if (i == 0) x_axis_switch = 1;
		else if (i == 1) x_axis_switch = -1;

		for (Int_t j = 0; j < _N_MODULES_TOTAL/2; ++j)
		{
			Double_t y = 0, x = 0;

			if ((j>=0) && (j<=4))
			{
				y = 45., x = (j-2)*15.;
				phi_angle_of_module[j + i*_N_MODULES_TOTAL/2] = ATan2(y,x_axis_switch*x);
			}
			else if ((j>=5) && (j<=39))
			{
				y = (3-(j+2)/7)*15, x = (3-(j+2)%7)*15;
				phi_angle_of_module[j + i*_N_MODULES_TOTAL/2] = ATan2(y,x_axis_switch*x);
			}
			else if ((j>=40) && (j<=44))
			{
				y = -45. , x = (j-42)*15.;
				phi_angle_of_module[j + i*_N_MODULES_TOTAL/2] = ATan2(y,x_axis_switch*x);
			}
		}
	}

	return phi_angle_of_module;
}

Float_t Unfold(Float_t phiEP_mc, Float_t psi_N_FULL, Int_t harm)
{
	Float_t values[10] , absvalues[10];

	values[0] = phiEP_mc - psi_N_FULL;
	absvalues[0] = Abs(values[0]);
	values[1] = phiEP_mc + 2. * Pi() / harm - psi_N_FULL;
	absvalues[1] = Abs(values[1]);
	values[2] = phiEP_mc - 2. * Pi() / harm - psi_N_FULL;
	absvalues[2] = Abs(values[2]);
    return values[LocMin(3, absvalues)];
}
