int main_flow(TString inFileName , TString outFileName, TString resFitFile, TString dcaFile)
{
	gSystem->Load("../Utilities/utility_cxx.so");
	gSystem->Load("./MpdCalculator_cxx.so");
	
	MpdCalculator mpd = MpdCalculator(inFileName,outFileName,dcaFile);
	mpd.CalculateFlow(0, resFitFile.Data());
	mpd.Write();
}
