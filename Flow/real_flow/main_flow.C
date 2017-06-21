int main_flow(TString inFileName , TString outFileName, TString resFitFile, TString dcaFile)
{
	gSystem->Load("./utility_C.so");
	gSystem->Load("./MpdCalculator_C.so");
	
	MpdCalculator mpd = MpdCalculator(inFileName,outFileName,dcaFile);
	mpd.CalculateFlow(0, resFitFile.Data());
	mpd.Write();
}
