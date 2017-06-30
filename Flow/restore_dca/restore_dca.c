#include "MpdFillDstTask.h"
#include "MpdKalmanTrack.h"
#include "MpdEctKalmanTrack.h"
#include "MpdTofMatching.h"
#include "MpdEtofMatching.h"
#include "MpdVertex.h"
#include "MpdHelix.h"
#include "MpdParticleIdentification.h"

#include "FairMCTrack.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include <TMath.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>

#include <iostream>

#include "restore_dca.h"

#define UNDEFINED_DCA -9999

const Double_t F_CUR0 = 0.3 * 0.01 * 5 / 10; // 5kG

MpdHelix MakeHelix(const MpdKalmanTrack *tr) 
{
    Double_t r = tr->GetPosNew();
    Double_t phi = tr->GetParam(0) / r;
    Double_t x = r * TMath::Cos(phi);
    Double_t y = r * TMath::Sin(phi);
    Double_t dip = tr->GetParam(3);
    Double_t cur = F_CUR0 * TMath::Abs (tr->GetParam(4));
    TVector3 o(x, y, tr->GetParam(1));
    Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
    MpdHelix helix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
    return helix;
}

void fill_zeroes(TClonesArray *MpdGlobalTracks)
{
    Int_t nGlobalTracks = MpdGlobalTracks->GetEntriesFast();
    for (Long_t track = 0; track < nGlobalTracks; ++track)
    {
        MpdTrack *mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(track); 
        mpdtrack->SetDCAX(UNDEFINED_DCA);
        mpdtrack->SetDCAY(UNDEFINED_DCA);
        mpdtrack->SetDCAZ(UNDEFINED_DCA);
    }
}

void restore_dca(TTree *inTree, TTree *outTree,
                 TClonesArray *mpdKalmanTracks, TClonesArray *MCTracks,
                 MpdEvent *MPDEvent, TClonesArray *vertexes)
{
    if (!inTree)
    {
        cout << "restore_dca(TTree *inTree, TTree *outTree, TClonesArray *mpdKalmanTracks, TClonesArray *MCTracks, MpdEvent *MPDEvent, TClonesArray *vertexes) Error - inTree is NULL" << endl;
        return;
    }
    if (!outTree)
    {
        cout << "restore_dca(TTree *inTree, TTree *outTree, TClonesArray *mpdKalmanTracks, TClonesArray *MCTracks, MpdEvent *MPDEvent, TClonesArray *vertexes) Error - outTree is NULL" << endl;
        return;
    }
    if (!mpdKalmanTracks)
    {
        cout << "restore_dca(TTree *inTree, TTree *outTree, TClonesArray *mpdKalmanTracks, TClonesArray *MCTracks, MpdEvent *MPDEvent, TClonesArray *vertexes) Error - mpdKalmanTracks is NULL" << endl;
        return;
    }
    if (!MCTracks)
    {
        cout << "restore_dca(TTree *inTree, TTree *outTree, TClonesArray *mpdKalmanTracks, TClonesArray *MCTracks, MpdEvent *MPDEvent, TClonesArray *vertexes) Error - MCTracks is NULL" << endl;
        return;
    }
    if (!MPDEvent)
    {
        cout << "restore_dca(TTree *inTree, TTree *outTree, TClonesArray *mpdKalmanTracks, TClonesArray *MCTracks, MpdEvent *MPDEvent, TClonesArray *vertexes) Error - MPDEvent is NULL" << endl;
        return;
    }
    if (!vertexes)
    {
        cout << "restore_dca(TTree *inTree, TTree *outTree, TClonesArray *mpdKalmanTracks, TClonesArray *MCTracks, MpdEvent *MPDEvent, TClonesArray *vertexes) Error - vertexes is NULL" << endl;
        return;
    }
    
    TClonesArray *MpdGlobalTracks=0;
    Int_t n_entries = inTree->GetEntries();

    for (int event = 0; event < n_entries; ++event)
    {
		cout << "EVENT N "<< event <<endl;
		inTree->GetEntry(event);
		
        Int_t nKalmanTracks = mpdKalmanTracks->GetEntriesFast();
        MpdGlobalTracks = MPDEvent->GetGlobalTracks();
        Int_t nGlobalTracks = MpdGlobalTracks->GetEntriesFast();
        MpdVertex *vertex = (MpdVertex*) vertexes->First();
		TVector3 primaryVertex;
        vertex->Position(primaryVertex);
        
        if (nGlobalTracks != nKalmanTracks)
        {
            cout << "no bijection..." << endl;
            fill_zeroes(MpdGlobalTracks);
            continue;
        }
        
        for (Long_t track = 0; track < nKalmanTracks; ++track)
		{
            MpdKalmanTrack *kalmanTrack = (MpdKalmanTrack*) mpdKalmanTracks->UncheckedAt(track);
			MpdTrack *mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(track); 
            if ((Int_t) mpdtrack->GetID() != (Int_t) kalmanTrack->GetTrackID())
            {
                cout << "tracks' ids are messed up..." << endl;
                fill_zeroes(MpdGlobalTracks);
                break;
            }
            
            MpdHelix helix = MakeHelix(kalmanTrack);
			Double_t pathLength = helix.pathLength(primaryVertex);
			TVector3 pca;
			pca = helix.at(pathLength);
			pca -= primaryVertex;
			mpdtrack->SetDCAX(pca.X());
			mpdtrack->SetDCAY(pca.Y());
			mpdtrack->SetDCAZ(pca.Z());
        }
        outTree->Fill();
	}
    outTree->AutoSave();
}

void restore_dca(const TString &inFileName, const TString &inFileTreeName, const TString &outFileName, const TString &mpdKalmanTracksBranchName,
                 const TString &MCTracksBranchName, const TString &MPDEventName, const TString &vertexesBranchName)
{
    TFile *inFile = new TFile(inFileName.Data(),"READ");
	TTree *inTree = (TTree*) inFile->Get(inFileTreeName);
	
    TFile *outFile = new TFile(outFileName.Data(),"RECREATE");
    TTree *outTree = inTree->CloneTree(0);
    
    TClonesArray *mpdKalmanTracks = (TClonesArray*) inFile->FindObjectAny(mpdKalmanTracksBranchName);
	inTree->SetBranchAddress(mpdKalmanTracksBranchName, &mpdKalmanTracks);
	TClonesArray *MCTracks=0;
	inTree->SetBranchAddress(MCTracksBranchName, &MCTracks);
	MpdEvent *MPDEvent=0;
	inTree->SetBranchAddress(MPDEventName, &MPDEvent);
    TClonesArray *vertexes = (TClonesArray*) inFile->FindObjectAny(vertexesBranchName);
	inTree->SetBranchAddress(vertexesBranchName, &vertexes);
    
    restore_dca(inTree, outTree, mpdKalmanTracks, MCTracks, MPDEvent, vertexes);
    
    inFile->Close();
    outFile->Close();
	return;
}

void restore_dca(const TString &inFileName, const TString &outFileName)
{
    restore_dca(inFileName, DCA_INPUT_TREE_NAME, outFileName, MPD_KALMAN_TRACKS_BRANCH_NAME, MC_TRACK_BRANCH_NAME, MPD_EVENT_TRACKS_BRANCH_NAME, VERTEXES_BRANCH_NAME);
}