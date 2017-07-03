int main_resolutions(TString inFileName , TString outFileName, TString dcaFileName)
{
	gSystem->Load("../Utilities/utility_cxx.so");
	cout << "loading libraries" << endl;
	//gSystem->Load("/lustre/nyx/hades/user/parfenov/real-flow/utility_C.so");
	gSystem->Load("./MpdCalculator_cxx.so");
	
	MpdCalculator mpd = MpdCalculator(inFileName,outFileName,dcaFileName);
	mpd.CalculateResolutions(0);
	mpd.Write();
}
