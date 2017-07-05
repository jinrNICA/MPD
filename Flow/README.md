#Quick information with step-by-step instruction.

1. Initial simulated files via MPDRoot and choosen generator (UrQMD, LAQGSM etc.) have information about distance of closest approach (DCA) according to z-dependent calcuations on the 2D transverse plain. There are more precise method based on 3D helisity fitting implemented in MPDRoot. To rewrite correct DCA values into cbmsim tree, use restore_dca:  
    root -l rootlogon.C    
    -------------------in root session---------------------  
    gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")  
    mpdloadlibs(kTRUE,kTRUE)  
    .L restore_dca.c+  
    restore_dca("inFileName","outFileName")  
    -------------------root session end--------------------  
Where arguments:  
    inFileName - input cbmsim tree with uncorrect dca values.  
    outFileName - output cbmsim tree with correct dca values. Everything else is the same as it was in inFileName.  
Macro restore_data.c copies cbmsim files and replaces incorrected dca values in the cbmsim trees.  
2. For the further analysis one does need the information about dca distributions to distinguish primary particles from secondary ones. While the former carry the needed information about initial geometry of the colliding system, the latter deteriorate this signal. Such process includes 3 main steps:  
    --Get dca distributions and store them into calibration file;  
    --Fit dca distributions via gaus function to make primary particles selection in the terms of n-sigma;  
    --Fit pt dependence of the dca distributions via polinomial function to reduce pt efficiency loss due to the dca distributions are splitted into discrete pt bins.  
To get light calibration file containing only histograms with dca distributions get_dca.cxx is used:  
    root -l rootlogon.C  
    -------------------in root session---------------------  
    gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")  
    mpdloadlibs(kTRUE,kTRUE)  
    .L get_dca.cxx+  
    get_dca("inFileName","outFileName")  
    -------------------root session end--------------------  
Where arguments:  
    inFileName - input cbmsim tree with corrected dca values (see restore_dca.h).  
    outFileName - output: standard root file with TH1* histograms needed for the further DCA cuts.  
Resulting file contains TH1* histograms of the dca distributions.  
!!! If get_dca() was used for several files, use hadd to merge outFileName files into 1 merged file. !!!  
Next, get_fit.cxx is used for 1-st iteration fitting procedure. It fits dca distributions with gaus functions:  
    root -l  
    -------------------in root session---------------------  
    .L get_fit.cxx+  
    get_fit("inFileName","outFileName")  
    -------------------root session end--------------------  
Where arguments:  
    inFileName - input: standard root file with TH1* histograms needed for the further DCA cuts (output file from get_dca(...)).  
    outFileName - output: standard root file with TH1* histograms and TF1* fitted function needed for the further DCA cuts.  
Resulting file contains TH1* histograms from the get_dca.cxx output file and their TF1* fitted functions.  
Finaly, to be able to distinguish primary particles without pt efficiency loss due to pt-dependence of the dca distributions, 2-nd iteration of the fitting procedure is used:  
    root -l  
    -------------------in root session---------------------  
    .L MakeFitDCA.cxx+  
    MakeFitDCA("inFileName","outFileName")  
    -------------------root session end--------------------  
Where arguments:  
    inFileName - input: standard root file with TH1* histograms and TF1* fitted function needed for the further DCA cuts (output from get_fit(...)).  
    outFileName - output: standard root file with TH1* histograms and TF1* fitted function needed for the further DCA cuts (pt_sigma_fit) - for improving pt efficiency.  
Resulting file contains sigma_pt_fit - TF1* functions that will be using furhter.  
3. To get centrality values from multiplicity in TPC one should get the inforamtion about multiplicity in all used statistics. Thus, one should use calibration file with multiplicity along with dca fit files in the previous step. To get multiplicity calibration file, use get_multiplicity.cxx in get_centrality directory:  
    root -l rootlogon.C  
    -------------------in root session---------------------  
    gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")  
    mpdloadlibs(kTRUE,kTRUE)  
    .L get_multiplicity.cxx+  
    get_multiplicity("inFileName","outFileName","dcaFileName")  
    -------------------root session end--------------------  
Where arguments:  
    inFileName - input cbmsim tree with corrected dca values (see restore_dca.h).  
    outFileName - output: standard root file with TH1* histograms needed for the further DCA cuts.  
    dcaFileName - Second iteration of the dca fitting containing sigma_pt_fit TF1* functions (output MakeFitDCA(...)).  
Resulting file contains TH1* histogram of multiplicity n TPC.  
!!! If get_multiplicity() was used for several files, use hadd to merge outFileName files into 1 merged file. !!!  
4. For the further analysis faster data construction is needed since flow measurements requires large statistics (100k~1M events). To convert cbmsim into cbmsim_reduced files (which are standard TTrees) use reducedTreeCreator.C in create_reduced_tree directory:  
    root -l rootlogon.C  
    -------------------in root session---------------------  
    gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")  
    mpdloadlibs(kTRUE,kTRUE)  
    reducedTreeCreator rtc = reducedTreeCreator("inFileHistName", "inFileTreeName", "outFileName", "dcaFileName")  
    rtc.CreateReducedTree()  
    -------------------root session end--------------------  
Where arguments:  
    inFileHistName - input standard root file with TH1* histograms of multiplicity (output get_multiplicity(...)).  
    inFileTreeName - input cbmsim tree with corrected dca values (output restore_dca(...)).  
    outFileName    - output standard root file with TTree cbmsim_reduced for the further analysis.  
    dcaFileName - Second iteration of the dca fitting containing sigma_pt_fit TF1* functions (output MakeFitDCA(...)).  
Resulting file contains standard TTree cbmsim_reduced with all needed information. This TTrees works ~1000 times faster than cbmsim.  
5. To exclude effects of the detector, one should calculate resolution correction factor before measuring the azimuthal flow. MpdCalculator.cxx in real-flow directory allows to do that:  
    root -l  
    -------------------in root session---------------------  
    gSystem->Load("libMathMore")  
    .L ../Utilities/BinningData.cxx+  
    .L ../Utilities/utility.cxx+  
    .L MpdCalculator.cxx+  
    MpdCalculator mpd = MpdCalculator(inFileName,outFileName,dcaFileName)  
    mpd.CalculateResolutions(0)  
    mpd.Write()  
    -------------------root session end--------------------  
Where arguments:  
    inFileName - input standard root file with TTree cbmsim_reduced for the further analysis (output from reducedTreeCreator).  
    outFileName - output standard root file with TProfiles and histograms containing data for resolution correction factor.  
    dcaFileName - Second iteration of the dca fitting containing sigma_pt_fit TF1* functions (output MakeFitDCA(...)).  
Resulting file contains TProfiles of cosine functions of the subevents angles.  
!!! If MpdCalculator::CalculateResolutions(0) was used for several files, use hadd to merge outFileName files into 1 merged file. !!!  
To calculate resolution correction factor from cosine of the difference between 2 subevent angles, averaged over events < cos(Psi_A-Psi_B) >, get_res.cxx is used:  
    root -l  
    -------------------in root session---------------------  
    gSystem->Load("libMathMore")  
    .L get_res.cxx+  
    get_res("inFileName","outFileName")  
    -------------------root session end--------------------  
Where arguments:  
    inFileName - input: standard root file with TProfiles and histograms containing data for resolution correction factor (output from MpdCalculator::CalculateResolutions(0)).  
    outFileName - output: standard root file with TH1* histograms and TF1* fitted function needed for the azimuthal flow calculation (MpdCalculator::CalculateFlow(...)).  
Resulting file contains TH1* histograms and TF1* functions of the resolution correction factor.  
6. To calculate azimuthal flow use MpdCalculator again:  
    root -l  
    -------------------in root session---------------------  
    gSystem->Load("libMathMore")  
    .L ../Utilities/BinningData.cxx+  
    .L ../Utilities/utility.cxx+  
    .L MpdCalculator.cxx+  
    MpdCalculator mpd = MpdCalculator(inFileName,outFileName,dcaFileName)  
    mpd.CalculateFlow(0, resFitFile.Data())  
    mpd.Write()  
    -------------------root session end--------------------  
Where arguments:  
    inFileName - input standard root file with TTree cbmsim_reduced for the further analysis (output from reducedTreeCreator).  
    outFileName - output standard root file with TProfiles and histograms containing data for resolution correction factor.  
    dcaFileName - Second iteration of the dca fitting containing sigma_pt_fit TF1* functions (output MakeFitDCA(...)).  
    resFitFile - input standard root file with TH1* histograms and TF1* fitted function of resolution (output from get_res(...)).  
OUTPUT:  
Resulting file contains TProfiles of v1 and v2 as a function of:  
    centrality (reconstructed: p_flow_wrt_full_vs_centrality_divided, generated: p_flow_wrt_RP_vs_centrality),  
    transverse momentum (reconstructed: p_flow_wrt_full_vs_pt_divided, generated: p_flow_wrt_RP_vs_pt),  
    pseudorapidity (reconstructed: p_flow_wrt_full_vs_eta_divided, generated: p_flow_wrt_RP_vs_eta),  
    rapidity (reconstructed: p_flow_wrt_full_vs_rapidity_divided, generated: p_flow_wrt_RP_vs_rapidity).  
!!! By the default, azimuthal flow is calculated for protons. !!!  
To change it to pions or kaons one should change the line 700 in MpdCalculator.cxx:  
from  
    Int_t sort = 0; - proton  
to  
    Int_t sort = 1; - kaon  
or  
    Int_t sort = 2; - pion  
and recompile the MpdCalculator.cxx.  
