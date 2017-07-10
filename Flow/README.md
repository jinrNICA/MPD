# Quick information with step-by-step instruction.
***
## Simulation of the data

To simulate data one needs to install correct FHCal geometry, install generator (UrQMD) and run MPDRoot simulation code.

### 1. Installation of the correct FHCal geometry in the MPDRoot framework

<u><big>1.1</big></u> Use branch marina_070417 in GIT repository:
<center><a href="https://git.jinr.ru/nica/mpdroot" target="_blank">https://git.jinr.ru/nica/mpdroot</a></center>  
One can download files from GIT using:

        git clone git@git.jinr.ru:nica/mpdroot.git
        git checkout marina_070417

<u><big>1.2</big></u> To create corresponding root file with correct FHCal geometry for reconstruction, run:

        root -l mpdroot/macro/mpd/geometry/create_rootgeom_zdc_oldnames_7sect_noSlots_v1.C

Such root file will be stored in the directory `mpdroot/geometry/zdc_oldnames_7sect_noSlots_v1.root`. 

Then include root file name into mpd geometry file `mpdroot/macro/mpd/geometry_stage1.C` (line 51):

        Zdc->SetGeometryFileName("zdc_oldnames_7sect_noSlots_v1.root");

For the only PSD performance we keep information about hits in the way:
- FHCal number
- module number
- slice number (1 - 42)
- deposited energy in slice (in the event)

The number of hits in each event is not more than 2x45x42 which allows to
have rather small output root files.
In this case information about MC tracks is lost.

<u><big>1.3</big></u> For the further analysis, information  about MC tracks is needed. To have such information, one should do:

        cd mpdroot/zdc
        cp MpdZdc.cxx.ProcHitsTracks_newgeomZDC  MpdZdc.cxx

and recompile MPDRoot.

<u><big>1.4</big></u> One needs to change `mpdroot/macro/mpd/runMC.C`:
line 49: 
        
        #define URQMD
Optionally, change `auau.09gev.mbias.98k.ftn14` to `test.f14` everywhere in the code since UrQMD will be used as a particle generator.

<u><big>1.5</big></u> There are 2 tracking algorithm that can be used:
- Idealistic, called hit producer (by default);
- Realistic, called cluster finder.

To change from hit producer to cluster finder, find in `mpdroot/macro/mpd/reco.C` lines:

        MpdTpcHitProducer* hitPr = new MpdTpcHitProducer();
        hitPr->SetModular(0);
        fRun->AddTask(hitPr);
and change them into:

        MpdTpcDigitizerAZ* tpcDigitizer = new MpdTpcDigitizerAZ();
        tpcDigitizer->SetPersistence(kFALSE);
        fRun->AddTask(tpcDigitizer);

Also, emc tracking was not used, so lines

        FairTask *emcHP = new MpdEmcHitProducer();
        fRun->AddTask(emcHP);
can be deleted.


### 2. Installation of the UrQMD generator

<u><big>2.1</big></u> Download urqmd.tar.gz:  
    <center><a href="http://urqmd.org/download/src/urqmd-3.4.tar.gz" target="_blank">http://urqmd.org/download/src/urqmd-3.4.tar.gz</a></center>  
and place it in the directory.

<u><big>2.2</big></u> Untar file  

        tar xzvvf urqmd-3.4.tar.gz  

<u><big>2.3</big></u> Compile UrQMD

        cd urqmd-3.4/  
        make  

<u><big>2.4</big></u> Set UrQMD to generate 250 events of Au-Au collision at 11 GeV energy for impact parameter range (0,20) fm. Change urqmd-3.4/inputfile into  

        pro 197 79  
        tar 197 79  

        nev 250  
        imp -20  

        ecm 11  
        tim 200 200  
        rsd randomnumber  

        f13  
        #f14  
        f15  
        f16  
        f19  
        f20  

        xxx  

In order to set unique set of random number sequence, change `randomnumber` value every time the UrQMD generator is used. For example, one can do in bash:  

        sed -e "s|randomnumber|$RANDOM|" -i inputfile

### 3. Simulation of the data

Simulation is contain 3 main parts:  
-- Particle generation (UrQMD);  
-- Monte Carlo simulation (GEANT4);  
-- Reconstruction procedure (Tracking, etc.).

<u><big>3.1</big></u> First of all, one should set FairSoft and MPDRoot variables:  

        source MPDRoot/build/config.sh  

or source it to the place where `config.sh` is stored.  

        export SIMPATH=<path to FairSoft>  
        export ROOTSYS=$SIMPATH  
        export PATH=$SIMPATH/bin:$PATH  
        export LD_LIBRARY_PATH=$SIMPATH/lib:$SIMPATH/lib/root:$LD_LIBRARY_PATH  
        source geant4.sh  
        platform=(root-config --arch)  

<u><big>3.2</big></u> Simulate particles using UrQMD generator. Run `runqmd.bash` in urqmd-3.4 directory:

        . runqmd.bash

Resulting file will be `test.f14`.  

<u><big>3.3</big></u> Run Monte Carlo simulation to calculate particle emission inside the MPD detector. In the MPDRoot/macro/mpd/ directory run `runMC.C`:  

        root -l runMC.C(inFile, outFile, nStartEvent, nEvents, flag_store_FairRadLenPoint, FieldSwitcher)

Where the arguments are:  
`inFile` - input `test.f14` file from generator (`"./test.f14"` by the default).  
`outFile` - output file from GEANT4 (`"./evetest.root"` by the default).  
`nStartEvent` - starting GEANT4 from this event (`0` by default).  
`nEvents` - number of all events that was in `test.f14` (`250` by the default).  
`flag_store_FairRadLenPoint` should be set by the default (`kFALSE`).  
`FieldSwitcher` should be set by the default (`0`).  

If `test.f14` contains 250 events, the `runMC.C` macro can have only 1 argument:  

        root -l runMC.C("<path to test.f14>")

<u><big>3.4</big></u> Run reconstruction procedure to get final output with simulated data. In the MPDRoot/macro/mpd/ directory run `reco.C`:

        root -l reco.C(inFile, outFile, nStartEvent, nEvents, run_type)

Where the arguments are:  
`inFile` - input `evetest.root` file from GEANT4 simulation (`"./evetest.root"` by the default).  
`outFile` - output file from reconstruction procedure (`"./mpddst.root"` by the default).  
`nStartEvent` - starting reconstruction from this event (`0` by default).  
`nEvents` - number of all events that was in `evetest.root` (`250` by the default).  
`run_type` should be set by default (`"local"`).

If `evetest.root` contains 250 events and was stored in the same direction, the `reco.C` macro can be used as:  

        root -l reco.C


<u><big>3.5</big></u> The final step is to correct dca values.
This procedure is realized in the Flow analysis framework:
    <center><a href="https://github.com/jinrNICA/MPD" target="_blank">https://github.com/jinrNICA/MPD</a></center>  
    
It can be installed by

        git clone git@github.com:jinrNICA/MPD.git
Initial simulated files via MPDRoot and choosen generator (UrQMD) have information about distance of closest approach (DCA) according to z-dependent calcuations on the 2D transverse plain. There are more precise method based on 3D helisity fitting implemented in MPDRoot. To rewrite correct DCA values into cbmsim tree, use MPD/Flow/restore_dca:  

        root -l rootlogon.C  
        -------------------in root session---------------------  
        gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")  
        mpdloadlibs(kTRUE,kTRUE)  
        .L restore_dca.c+  
        restore_dca("inFileName","outFileName")  
        -------------------root session end--------------------  
Where arguments:  
    `inFileName` - input cbmsim tree with uncorrect dca values.  
    `outFileName` - output cbmsim tree with correct dca values. Everything else is the same as it was in inFileName.  
Macro `restore_data.c` copies cbmsim files and replaces incorrected dca values in the cbmsim trees.  

Resulting `mpddst.root` files contain results from both simulation via GEANT4 (generated data) and reconstruction procedure (reconstructed data).
***
## Creating light data for the analysis

Since azimuthal anisotropic flow analysis requires large statistics, simulated data should be converted into light and fast accessible data structure that allows to reduce CPU time. Convertion procedure requires 2 calibration files: dca and multiplicity. This procedure is realized in the Flow analysis framework (see <u><big>3.5</big></u>).

### 4. DCA calibration file

For the further analysis one does need the information about dca distributions to distinguish primary particles from secondary ones. While the former carry the needed information about initial geometry of the colliding system, the latter deteriorate this signal. Such process includes 3 main steps:  
    --Get dca distributions and store them into calibration file;  
    --Fit dca distributions via gaus function to make primary particles selection in the terms of n-sigma;  
    --Fit pt dependence of the dca distributions via polinomial function to reduce pt efficiency loss due to the dca distributions are splitted into discrete pt bins.  
    
<u><big>4.1</big></u>  To get light calibration file containing only histograms with dca distributions `get_dca.cxx` is used in the MPD/Flow/get_dca/ directory:  
    
        root -l rootlogon.C  
        -------------------in root session---------------------  
        gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")  
        mpdloadlibs(kTRUE,kTRUE)  
        .L get_dca.cxx+  
        get_dca("inFileName","outFileName")  
        -------------------root session end--------------------  
Where the arguments are:  
    `inFileName` - input cbmsim tree with corrected dca values (see `restore_dca.h`).  
    `outFileName` - output: standard root file with `TH1*` histograms needed for the further DCA cuts.  
If compilation failed due to `FairMCEventHeader.h` not found, delete this line.
Resulting file contains `TH1*` histograms of the dca distributions.  

<center>!!! If `get_dca()` was used for several files, use `hadd` to merge `outFileName` files into 1 merged file. !!!</center>  

<u><big>4.2</big></u> Next, `get_fit.cxx` is used for 1-st iteration fitting procedure. It fits dca distributions with gaus functions:  
        
        root -l  
        -------------------in root session---------------------  
        .L get_fit.cxx+  
        get_fit("inFileName","outFileName")  
        -------------------root session end--------------------  
Where the arguments are:  
    `inFileName` - input: standard root file with `TH1*` histograms needed for the further DCA cuts (output file from `get_dca(...)`).  
    `outFileName` - output: standard root file with `TH1*` histograms and `TF1*` fitted function needed for the further DCA cuts.  
    
Resulting file contains `TH1*` histograms from the `get_dca.cxx` output file and their `TF1*` fitted functions.  

<u><big>4.3</big></u> Finaly, to be able to distinguish primary particles without pt efficiency loss due to pt-dependence of the dca distributions, 2-nd iteration of the fitting procedure is used:  

        root -l  
        -------------------in root session---------------------  
        .L MakeFitDCA.cxx+  
        MakeFitDCA("inFileName","outFileName")  
        -------------------root session end--------------------  
Where the arguments are:  
    `inFileName` - input: standard root file with `TH1*` histograms and `TF1*` fitted function needed for the further DCA cuts (output from `get_fit(...)`).  
    `outFileName` - output: standard root file with `TH1*` histograms and `TF1*` fitted function needed for the further DCA cuts (`pt_sigma_fit`) - for improving pt efficiency.  
Resulting file contains `sigma_pt_fit` - `TF1*` functions that will be using furhter.  

### 5. Multiplicity calibration file

<u><big>5.1</big></u> To get centrality values from multiplicity in TPC one should get the inforamtion about multiplicity in all used statistics. Thus, one should use calibration file with multiplicity along with dca fit files in the previous step. To get multiplicity calibration file, use `get_multiplicity.cxx` in MPD/Flow/get_centrality directory:  

        root -l rootlogon.C  
        -------------------in root session---------------------  
        gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")  
        mpdloadlibs(kTRUE,kTRUE)  
        .L get_multiplicity.cxx+  
        get_multiplicity("inFileName","outFileName","dcaFileName")  
        -------------------root session end--------------------  
Where the arguments are:  
    `inFileName` - input cbmsim tree with corrected dca values (see `restore_dca.h`).  
    `outFileName` - output: standard root file with `TH1*` histograms needed for the further DCA cuts.  
    `dcaFileName` - Second iteration of the dca fitting containing `sigma_pt_fit TF1*` functions (output `MakeFitDCA(...)`).  
    
Resulting file contains `TH1*` histogram of multiplicity in TPC.  

<center>!!! If `get_multiplicity()` was used for several files, use hadd to merge `outFileName` files into 1 merged file. !!!</center>  

### 6. Conversion into light data files

<u><big>6.1</big></u> To convert cbmsim into cbmsim_reduced files (which are standard TTrees) use `reducedTreeCreator.C` in MPD/Flow/create_reduced_tree directory:  
    
        root -l rootlogon.C  
        -------------------in root session---------------------  
        gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")  
        mpdloadlibs(kTRUE,kTRUE)  
        reducedTreeCreator.C+
        reducedTreeCreator rtc = reducedTreeCreator("inFileHistName", "inFileTreeName", "outFileName", "dcaFileName")  
        rtc.CreateReducedTree()  
        -------------------root session end--------------------  
Where the arguments are:  
    `inFileHistName` - input standard root file with `TH1*` histograms of multiplicity (output from `get_multiplicity(...)`).  
    `inFileTreeName` - input cbmsim tree with corrected dca values (output from `restore_dca(...)`).  
    `outFileName` - output standard root file with TTree cbmsim_reduced for the further analysis.  
    `dcaFileName` - Second iteration of the dca fitting containing `sigma_pt_fit TF1*` functions (output from `MakeFitDCA(...)`).  
    
Resulting file contains standard TTree cbmsim_reduced with all needed information. This TTrees works ~1000 times faster than cbmsim.  
***
## Azimuthal anisotropic flow measurement

Azimuthal anisotropic flow measurements is implemented in Flow analysis framework (see <u><big>3.5</big></u>). It works with converted light data files (see <u><big>4-6</big></u>).

Before measuring direct and elliptic flow (v1 and v2) one should calculate the resolution correction factor.

### 7. Resolution correction factor

<u><big>7.1</big></u> To exclude effects of the detector, one should calculate resolution correction factor before measuring the azimuthal flow. `MpdCalculator.cxx` in MPD/Flow/real-flow directory allows to do that:  
    
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
Where the arguments are:  
    `inFileName` - input standard root file with TTree cbmsim_reduced for the further analysis (output from `reducedTreeCreator` class).  
    `outFileName` - output standard root file with TProfiles and histograms containing data for resolution correction factor.  
    `dcaFileName` - Second iteration of the dca fitting containing `sigma_pt_fit TF1*` functions (output from `MakeFitDCA(...)`).  
Resulting file contains TProfiles of cosine functions of the subevents angles.  

<center>!!! If `MpdCalculator::CalculateResolutions(0)` was used for several files, use hadd to merge `outFileName` files into 1 merged file. !!!</center>  

To calculate resolution correction factor from cosine of the difference between 2 subevent angles, averaged over events < cos(Psi_A-Psi_B) >, `get_res.cxx` is used:  
    
        root -l  
        -------------------in root session---------------------  
        gSystem->Load("libMathMore")  
        .L get_res.cxx+  
        get_res("inFileName","outFileName")  
        -------------------root session end--------------------  
Where the arguments are:  
    `inFileName` - input: standard root file with TProfiles and histograms containing data for resolution correction factor (output from `MpdCalculator::CalculateResolutions(0)`).  
    `outFileName` - output: standard root file with `TH1*` histograms and `TF1*` fitted function needed for the azimuthal flow calculation (`MpdCalculator::CalculateFlow(...)`).  
Resulting file contains `TH1*` histograms and `TF1*` functions of the resolution correction factor.  

### 8. Azimuthal anisotropic flow calculation

<u><big>8.1</big></u> To calculate azimuthal flow use `MpdCalculator.cxx` again:  
    
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
Where the arguments are:  
    `inFileName` - input standard root file with TTree cbmsim_reduced for the further analysis (output from `reducedTreeCreator`).  
    `outFileName` - output standard root file with TProfiles and histograms containing data for resolution correction factor.  
    `dcaFileName` - Second iteration of the dca fitting containing `sigma_pt_fit TF1*` functions (output `MakeFitDCA(...)`).  
    `resFitFile` - input standard root file with `TH1*` histograms and `TF1*` fitted function of resolution (output from `get_res(...)`).  

### OUTPUT:  
Resulting file contains TProfiles of v1 and v2 as a function of:  
    - centrality (reconstructed: p_flow_wrt_full_vs_centrality_divided, generated: p_flow_wrt_RP_vs_centrality),  
    - transverse momentum (reconstructed: p_flow_wrt_full_vs_pt_divided, generated: p_flow_wrt_RP_vs_pt),  
    - pseudorapidity (reconstructed: p_flow_wrt_full_vs_eta_divided, generated: p_flow_wrt_RP_vs_eta),  
    - rapidity (reconstructed: p_flow_wrt_full_vs_rapidity_divided, generated: p_flow_wrt_RP_vs_rapidity).  
!!! By the default, azimuthal flow is calculated for protons. !!!  
To change it to pions or kaons one should change the line 700 in `MpdCalculator.cxx`:  
from  
    `Int_t sort = 0;` - proton  
to  
    `Int_t sort = 1;` - kaon  
or  
    `Int_t sort = 2;` - pion  
and recompile the `MpdCalculator.cxx`.  
***