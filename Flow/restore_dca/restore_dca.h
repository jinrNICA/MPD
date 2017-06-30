#ifndef RESTORE_DCA_H
#define RESTORE_DCA_H

#define UNDEFINED_DCA -9999

#include <TString.h>

const TString DCA_INPUT_TREE_NAME = "cbmsim";
const TString MPD_KALMAN_TRACKS_BRANCH_NAME = "TpcKalmanTrack";
const TString MC_TRACK_BRANCH_NAME = "MCTrack";
const TString MPD_EVENT_TRACKS_BRANCH_NAME = "MPDEvent.";
const TString VERTEXES_BRANCH_NAME = "Vertex";

class MpdGlobalTracks;
class MpdHelix;
class MpdKalmanTrack;
class TTree;
class TClonesArray;

MpdHelix MakeHelix(const MpdKalmanTrack *tr);
void fill_zeroes(TClonesArray *MpdGlobalTracks);

void restore_dca(const TString &inFileName, const TString &outFileName);

void restore_dca(TTree *inTree, TTree *outTree,
                 TClonesArray *mpdKalmanTracks, TClonesArray *MCTracks,
                 MpdEvent *MPDEvent, TClonesArray *vertexes);
                 
void restore_dca(const TString &inFileName, const TString &inFileTreeName, const TString &outFileName, const TString &mpdKalmanTracksBranchName,
                 const TString &MCTracksBranchName, const TString &MPDEventName, const TString &vertexesBranchName);

#endif // RESTORE_DCA_H