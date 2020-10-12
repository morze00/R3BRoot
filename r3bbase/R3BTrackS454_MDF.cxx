/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// ------------------------------------------------------------
// -----                  R3BTrackS454_MDF                -----
// -----          Created April 13th 2016 by M.Heil       -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BCalifaMappedData.h"
#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
#include "R3BLosMappedData.h"

#include "R3BBeamMonitorMappedData.h"

#include "R3BTrackS454_MDF.h"

#include "R3BSci8CalData.h"
#include "R3BSci8MappedData.h"

#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdMappedData.h"

#include "R3BRoluCalData.h"
#include "R3BRoluMappedData.h"

#include "R3BPaddleCalData.h"

#include "R3BPspxCalData.h"
#include "R3BPspxMappedData.h"

#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"

#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberHitData.h"
#include "R3BBunchedFiberMappedData.h"

#include "R3BMCTrack.h"
#include "R3BTrack.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TCutG.h"
#include "R3BMDFWrapper.h"

#include "TClonesArray.h"
#include "TFile.h"
#include "TMath.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

R3BTrackS454_MDF::R3BTrackS454_MDF()
    : R3BTrackS454_MDF("Track", 1)
{
}

R3BTrackS454_MDF::R3BTrackS454_MDF(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fCuts(0)
    , fGraphCuts(0)
    , fGhost(0)
    , fPairs(0)
    , fB(-1672)
    , fSimu(0)
      , fNEvents(0)
    , fTrackItems(new TClonesArray("R3BTrack"))
      , fNofTrackItems()
{
}

R3BTrackS454_MDF::~R3BTrackS454_MDF()
{
    delete fTrackItems;
}

InitStatus R3BTrackS454_MDF::Init()
{
    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(INFO) << "R3BTrackS454_MDF::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");
    FairRunOnline* run = FairRunOnline::Instance();

    // Get objects for detectors on all levels
    fMCTrack = (TClonesArray*)mgr->GetObject("MCTrack");
    if (fMCTrack)
        mgr->Register("MCTrack", "Monte Carlo Tracks", fMCTrack, kTRUE);

    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    printf("Have %d fiber detectors.\n", NOF_FIB_DET);
    fMappedItems.push_back((TClonesArray*)mgr->GetObject(Form("%sMappedData", fDetectorNames[0])));
    if (NULL == fMappedItems.at(0))
    {
        printf("Could not find mapped data for '%s'.\n", fDetectorNames[0]);
    }
    fCalItems.push_back((TClonesArray*)mgr->GetObject(Form("%sCrystalCalData", fDetectorNames[0])));
    if (NULL == fCalItems.at(0))
    {
        printf("Could not find Cal data for '%s'.\n", fDetectorNames[0]);
    }
    fHitItems.push_back((TClonesArray*)mgr->GetObject(Form("%sHit", fDetectorNames[0])));
    if (NULL == fHitItems.at(0))
    {
        printf("Could not find hit data for '%s'.\n", fDetectorNames[0]);
    }
    for (int det = 1; det < DET_MAX; det++)
    {
        fMappedItems.push_back((TClonesArray*)mgr->GetObject(Form("%sMapped", fDetectorNames[det])));
        if (NULL == fMappedItems.at(det))
        {
            printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
        }
        if (det == 9)
            maxevent = mgr->CheckMaxEventNo();
        fCalItems.push_back((TClonesArray*)mgr->GetObject(Form("%sCal", fDetectorNames[det])));
        if (NULL == fCalItems.at(det))
        {
            printf("Could not find Cal data for '%s'.\n", fDetectorNames[det]);
        }
        fHitItems.push_back((TClonesArray*)mgr->GetObject(Form("%sHit", fDetectorNames[det])));
        if (NULL == fHitItems.at(det))
        {
            printf("Could not find hit data for '%s'.\n", fDetectorNames[det]);
        }
    }

    mgr->Register("Track", "Land", fTrackItems, kTRUE);

    //------------------------------------------------------------------------
    // graphical cuts
    //------------------------------------------------------------------------
    if (fGraphCuts)
    {
        cut_Fi13vsTofd = NULL;
        cut_Fi10vsTofd = NULL;

        TFile* f = TFile::Open("myCuts.root", "read");
        cut_Fi13vsTofd = dynamic_cast<TCutG*>(f->Get("Fi13vsTofd"));
        cut_Fi10vsTofd = dynamic_cast<TCutG*>(f->Get("Fi10vsTofd"));
    }
    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    //-----------------------------------------------------------------------
    // BeamMonitor

    // get the theoretical calib factors for SEETRAM
    Double_t fexp = float(fsens_SEE + 9);
    Double_t fpow = float(pow(10., fexp));
    calib_SEE = 135641.7786 * fpow;
    LOG(DEBUG) << fsens_SEE << ", " << fexp << ", " << fpow << ", " << calib_SEE << endl;

    // MDF tracker.
    if (tracker)
    {
        //Reconstruction of the target X coordinate
        //MDF_X0 = new R3BMDFWrapper();
        //MDF_X0->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/X0/MDF_X0_20200622.txt");
        //MDF_X0->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/X0/PCA_X0_20200622.txt");

        //New 
        MDF_X0 = new R3BMDFWrapper();
        MDF_X0->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/X0/MDF_X0_20200824.txt");
        MDF_X0->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/X0/PCA_X0_20200824.txt");


        //Reconstruction of the target TX0 angle using posiion X0 at the target
        MDF_TX0_targ = new R3BMDFWrapper();
        MDF_TX0_targ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/TX0/MDF_TX0_20200622.txt");
        MDF_TX0_targ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/TX0/PCA_TX0_20200622.txt");

        //Reconstruction of the target TX0 angle using positions in the fiber3 
        MDF_TX0_f3 = new R3BMDFWrapper();
        MDF_TX0_f3->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/TX0/MDF_TX0_20200824.txt");
        MDF_TX0_f3->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/TX0/PCA_TX0_20200824.txt");


        MDF_PoQ = new R3BMDFWrapper();
        //Reconstruction of the PoQ (produciton runs)

        //Reconstruction of the PoQ (production runa, I=1672A)
        //if(fB == -1672)
        //{
        MDF_PoQ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/MDF_PoQ_1672A_20200622.txt");
        MDF_PoQ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/PCA_PoQ_1672A_20200622.txt");
        //}

        //Reconstruction of the PoQ (run 395, I=1102A)
        //if(fB == -1102)
        //{
        //    MDF_PoQ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/MDF_PoQ_1102A_20200713.txt");
        //    MDF_PoQ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/PCA_PoQ_1102A_20200713.txt");
        //}

        ////Reconstruction of the PoQ (run 391, I=1292A)
        //if(fB == -1292)
        //{
        //    MDF_PoQ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/MDF_PoQ_1292A_20200713.txt");
        //    MDF_PoQ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/PCA_PoQ_1292A_20200713.txt");
        //}

        ////Reconstruction of the PoQ (run 389, I=1482A)
        //if(fB == -1482)
        //{
        //    MDF_PoQ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/MDF_PoQ_1482A_20200713.txt");
        //    MDF_PoQ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/PCA_PoQ_1482A_20200713.txt");
        //}
        //
        ////Reconstruction of the PoQ (run 399, I=1862A)
        //if(fB == -1862)
        //{
        //    MDF_PoQ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/MDF_PoQ_1862A_20200803.txt");
        //    MDF_PoQ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/PCA_PoQ_1862A_20200803.txt");
        //}

        ////Reconstruction of the PoQ (run 401, I=2052)
        //if(fB == -2052)
        //{
        //    MDF_PoQ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/MDF_PoQ_2052A_20200803.txt");
        //    MDF_PoQ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/PCA_PoQ_2052A_20200803.txt");
        //}
        //
        ////Reconstruction of the PoQ (run 405, I=2242)
        //if(fB == -2242)
        //{
        //    MDF_PoQ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/MDF_PoQ_2242A_20200803.txt");
        //    MDF_PoQ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/PCA_PoQ_2242A_20200803.txt");
        //}


    }

    tree_out.SetName("tree_out");

    tree_out.Branch("N_glob_tracks_left",     &N_glob_tracks_left,  "N_glob_tracks_left/i");
    tree_out.Branch("N_glob_tracks_right",     &N_glob_tracks_right,  "N_glob_tracks_right/i");

    //============= Left arm detectors ================

    tree_out.Branch("f10_X",     f10_X,  "f10_X[N_glob_tracks_left]/F");
    tree_out.Branch("f10_Y",     f10_Y,  "f10_Y[N_glob_tracks_left]/F");
    tree_out.Branch("f10_Z",     f10_Z,  "f10_Z[N_glob_tracks_left]/F");
    tree_out.Branch("f10_Q",     f10_Q,  "f10_Q[N_glob_tracks_left]/F");
    tree_out.Branch("f10_T",     f10_T,  "f10_T[N_glob_tracks_left]/F");

    tree_out.Branch("f12_X",     f12_X,  "f12_X[N_glob_tracks_left]/F");
    tree_out.Branch("f12_Y",     f12_Y,  "f12_Y[N_glob_tracks_left]/F");
    tree_out.Branch("f12_Z",     f12_Z,  "f12_Z[N_glob_tracks_left]/F");
    tree_out.Branch("f12_Q",     f12_Q,  "f12_Q[N_glob_tracks_left]/F");
    tree_out.Branch("f12_T",     f12_T,  "f12_T[N_glob_tracks_left]/F");

    tree_out.Branch("f3b_X",     f3b_X,  "f3b_X[N_glob_tracks_left]/F");
    tree_out.Branch("f3b_Y",     f3b_Y,  "f3b_Y[N_glob_tracks_left]/F");
    tree_out.Branch("f3b_Z",     f3b_Z,  "f3b_Z[N_glob_tracks_left]/F");
    tree_out.Branch("f3b_Q",     f3b_Q,  "f3b_Q[N_glob_tracks_left]/F");
    tree_out.Branch("f3b_T",     f3b_T,  "f3b_T[N_glob_tracks_left]/F");
    
    tree_out.Branch("tofd_left_X",     tofd_left_X,  "tofd_left_X[N_glob_tracks_left]/F");
    tree_out.Branch("tofd_left_Y",     tofd_left_Y,  "tofd_left_Y[N_glob_tracks_left]/F");
    tree_out.Branch("tofd_left_Z",     tofd_left_Z,  "tofd_left_Z[N_glob_tracks_left]/F");
    tree_out.Branch("tofd_left_Q",     tofd_left_Q,  "tofd_left_Q[N_glob_tracks_left]/F");
    tree_out.Branch("tofd_left_T",     tofd_left_T,  "tofd_left_T[N_glob_tracks_left]/F");

    //================ Right arm detectors ==============

    tree_out.Branch("f11_X",     f11_X,  "f11_X[N_glob_tracks_right]/F");
    tree_out.Branch("f11_Y",     f11_Y,  "f11_Y[N_glob_tracks_right]/F");
    tree_out.Branch("f11_Z",     f11_Z,  "f11_Z[N_glob_tracks_right]/F");
    tree_out.Branch("f11_Q",     f11_Q,  "f11_Q[N_glob_tracks_right]/F");
    tree_out.Branch("f11_T",     f11_T,  "f11_T[N_glob_tracks_right]/F");

    tree_out.Branch("f13_X",     f13_X,  "f13_X[N_glob_tracks_right]/F");
    tree_out.Branch("f13_Y",     f13_Y,  "f13_Y[N_glob_tracks_right]/F");
    tree_out.Branch("f13_Z",     f13_Z,  "f13_Z[N_glob_tracks_right]/F");
    tree_out.Branch("f13_Q",     f13_Q,  "f13_Q[N_glob_tracks_right]/F");
    tree_out.Branch("f13_T",     f13_T,  "f13_T[N_glob_tracks_right]/F");

    tree_out.Branch("f3a_X",     f3a_X,  "f3a_X[N_glob_tracks_right]/F");
    tree_out.Branch("f3a_Y",     f3a_Y,  "f3a_Y[N_glob_tracks_right]/F");
    tree_out.Branch("f3a_Z",     f3a_Z,  "f3a_Z[N_glob_tracks_right]/F");
    tree_out.Branch("f3a_Q",     f3a_Q,  "f3a_Q[N_glob_tracks_right]/F");
    tree_out.Branch("f3a_T",     f3a_T,  "f3a_T[N_glob_tracks_right]/F");
    
    tree_out.Branch("tofd_right_X",     tofd_right_X,  "tofd_right_X[N_glob_tracks_right]/F");
    tree_out.Branch("tofd_right_Y",     tofd_right_Y,  "tofd_right_Y[N_glob_tracks_right]/F");
    tree_out.Branch("tofd_right_Z",     tofd_right_Z,  "tofd_right_Z[N_glob_tracks_right]/F");
    tree_out.Branch("tofd_right_Q",     tofd_right_Q,  "tofd_right_Q[N_glob_tracks_right]/F");
    tree_out.Branch("tofd_right_T",     tofd_right_T,  "tofd_right_T[N_glob_tracks_right]/F");

    //======== Tracking variables left arm ============= 

    tree_out.Branch("TX0_mdf_left",       TX0_mdf_left,      "TX0_mdf_left[N_glob_tracks_left]/F");
    tree_out.Branch("TX0_f3_mdf_left",    TX0_f3_mdf_left,   "TX0_f3_mdf_left[N_glob_tracks_left]/F");
    tree_out.Branch("PoQ_mdf_left",       PoQ_mdf_left,      "PoQ_mdf_left[N_glob_tracks_left]/F");
    tree_out.Branch("X0_residual_left",   X0_residual_left,  "X0_residual_left[N_glob_tracks_left]/F");
    tree_out.Branch("TX0_residual_left",  TX0_residual_left, "TX0_residual_left[N_glob_tracks_left]/F");

    tree_out.Branch("TX0_residual_f3_left",   TX0_residual_f3_left,  "TX0_residual_f3_left[N_glob_tracks_left]/F");
    tree_out.Branch("TX0_residual_mdf_left",   TX0_residual_mdf_left,  "TX0_residual_mdf_left[N_glob_tracks_left]/F");
    tree_out.Branch("X0_proj_by_f3_left",     X0_proj_by_f3_left,    "X0_proj_by_f3_left[N_glob_tracks_left]/F");
    tree_out.Branch("TX0_proj_by_f3_left",    TX0_proj_by_f3_left,   "TX0_proj_by_f3_left[N_glob_tracks_left]/F");
    tree_out.Branch("Xoffset_f3b",        Xoffset_f3b,       "Xoffset_f3b[N_glob_tracks_left]/F");
    
    //======== Tracking variables right arm ============= 

    tree_out.Branch("TX0_mdf_right",       TX0_mdf_right,      "TX0_mdf_right[N_glob_tracks_right]/F");
    tree_out.Branch("TX0_f3_mdf_right",    TX0_f3_mdf_right,   "TX0_f3_mdf_right[N_glob_tracks_right]/F");
    tree_out.Branch("PoQ_mdf_right",       PoQ_mdf_right,      "PoQ_mdf_right[N_glob_tracks_right]/F");
    tree_out.Branch("X0_residual_right",   X0_residual_right,  "X0_residual_right[N_glob_tracks_right]/F");
    tree_out.Branch("TX0_residual_right",  TX0_residual_right, "TX0_residual_right[N_glob_tracks_right]/F");
    tree_out.Branch("TX0_residual_f3_right",   TX0_residual_f3_right,  "TX0_residual_f3_right[N_glob_tracks_right]/F");
    tree_out.Branch("TX0_residual_mdf_right",   TX0_residual_mdf_right,  "TX0_residual_mdf_right[N_glob_tracks_right]/F");
    tree_out.Branch("X0_proj_by_f3_right",     X0_proj_by_f3_right,    "X0_proj_by_f3_right[N_glob_tracks_right]/F");
    tree_out.Branch("TX0_proj_by_f3_right",    TX0_proj_by_f3_right,   "TX0_proj_by_f3_right[N_glob_tracks_right]/F");
    tree_out.Branch("Xoffset_f3a",        Xoffset_f3a,       "Xoffset_f3a[N_glob_tracks_right]/F");

    //======== Tracking variables right arm ============= 
    
    tree_out.Branch("Nhits_f3a",   &Nhits_f3a,  "Nhits_f3a/i");
    tree_out.Branch("Nhits_f3b",   &Nhits_f3b,  "Nhits_f3b/i");
    tree_out.Branch("Nhits_f10",   &Nhits_f10,  "Nhits_f10/i");
    tree_out.Branch("Nhits_f11",   &Nhits_f11,  "Nhits_f11/i");
    tree_out.Branch("Nhits_f12",   &Nhits_f12,  "Nhits_f12/i");
    tree_out.Branch("Nhits_f13",   &Nhits_f13,  "Nhits_f13/i");
    tree_out.Branch("Nhits_tofd",   &Nhits_tofd,  "Nhits_tofd/i");

    return kSUCCESS;
}

void R3BTrackS454_MDF::Exec(Option_t* option)
{
    if (fNEvents / 10000. == (int)fNEvents / 10000)
        std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
            << " %) "
            << " Tofd: " << counterTofd << " tracked: " << counter2 << " chix: " << counter3
            << " chiy: " << counter4 << std::flush;
    fNEvents += 1;

    //cout << "\n\nNew event ******************************" << endl;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    //===========================================================================

    if (header)
    {
        time = header->GetTimeStamp();

        if (time_start == 0 && time > 0)
        {
            time_start = time;
            fNEvents_start = fNEvents;
        }

        if (header->GetTrigger() == 12)
        {
            // spill start in nsec
            // cout << "spill start" << endl;
            num_spills++;
        }
        if (header->GetTrigger() == 13)
        {
            // spill end  in nsec
            // cout << "spill stop" << endl;
        }

        //   check for requested trigger (Todo: should be done globablly / somewhere else)
        if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        {
            counterWrongTrigger++;
            return;
        }

        // fTpat = 1-16; fTpat_bit = 0-15
        Int_t fTpat_bit = fTpat - 1;
        Int_t itpat;
        Int_t tpatvalue;
        if (fTpat_bit >= 0)
        {
            itpat = header->GetTpat();
            tpatvalue = (itpat && (1 << fTpat_bit)) >> fTpat_bit;
            if (tpatvalue == 0)
            {
                counterWrongTpat++;
                return;
            }
        }
    }//end if(header)

    //===========================================================================

    if (fMappedItems.at(DET_BMON))
    {
        unsigned long IC;
        unsigned long SEETRAM_raw;
        Double_t SEETRAM;
        unsigned long TOFDOR;

        auto detBmon = fMappedItems.at(DET_BMON);
        Int_t nHitsbm = detBmon->GetEntriesFast();
        // cout<<"Bmon hits: "<<nHitsbm<<endl;

        for (Int_t ihit = 0; ihit < nHitsbm; ihit++)
        {
            R3BBeamMonitorMappedData* hit = (R3BBeamMonitorMappedData*)detBmon->At(ihit);
            if (!hit)
                continue;

            IC = hit->GetIC(); // negative values if offset not high enough
            counts_IC += (double)IC;

            SEETRAM_raw = hit->GetSEETRAM();           // raw counts
            SEETRAM = (double)SEETRAM_raw * calib_SEE; // calibrated SEETRAM counts
            counts_SEE += SEETRAM;

            TOFDOR = hit->GetTOFDOR(); // only positive values possible
            counts_TofD += TOFDOR;

            if (fNEvents == fNEvents_start)
            {
                see_start = SEETRAM;
                ic_start = IC;
                tofdor_start = TOFDOR;
            }
        }
    }//end if (fMappedItems.at(DET_BMON))

    //===========================================================================

    Bool_t RoluCut = false;
    if (fMappedItems.at(DET_ROLU))
    {
        // rolu
        auto detRolu = fMappedItems.at(DET_ROLU);
        Int_t nHitsRolu = detRolu->GetEntriesFast();
        // cout<<"ROLU hits: "<<nHitsRolu<<endl;

        for (Int_t ihit = 0; ihit < nHitsRolu; ihit++)
        {
            R3BRoluMappedData* hitRolu = (R3BRoluMappedData*)detRolu->At(ihit);
            if (!hitRolu)
                continue;

            // channel numbers are stored 1-based (1..n)
            Int_t iDet = hitRolu->GetDetector(); // 1..
            Int_t iCha = hitRolu->GetChannel();  // 1..
            RoluCut = true;
        }
    }

    if (RoluCut)
    {
        counterRolu++;
        return;
    }

    //===========================================================================

    Bool_t CalifaHit = false;
    if (fMappedItems.at(DET_CALIFA))
    {
        // CALIFA
        auto detCalifa = fMappedItems.at(DET_CALIFA);
        Int_t nHitsCalifa = detCalifa->GetEntriesFast();
        // cout<<"Califa hits: "<<nHitsCalifa<<endl;

        for (Int_t ihit = 0; ihit < nHitsCalifa; ihit++)
        {
            R3BCalifaMappedData* hitCalifa = (R3BCalifaMappedData*)detCalifa->At(ihit);
            if (!hitCalifa)
                continue;

            Int_t Crystal = hitCalifa->GetCrystalId();
            Int_t Energy = hitCalifa->GetEnergy();
            // cout << "Califa: " << Crystal << " Energy: " << Energy << endl;
            if (Energy > 0)
            {
                CalifaHit = true;
            }
        }
    }
    if (CalifaHit)
    {
        counterCalifa++;
        //return;
    }

    //variables to work with
    //===========================================================================

    Int_t max = 10000;
    Int_t detector[max];
    Double_t xdet[max];
    Double_t ydet[max];
    Double_t zdet[max];
    Double_t tdet[max];
    Double_t qdet[max];

    countdet = 0;

    Double_t track[12]; // x,y,z, px, py, pz
    Double_t chi[6];    // x,y,z, px, py, pz

    Int_t n_det = 10;

    Double_t x[n_det];
    Double_t y[n_det];
    Double_t z[n_det];
    Double_t q[n_det];
    Double_t t[n_det];
    Double_t x1[n_det];
    Double_t y1[n_det];
    Double_t z1[n_det];
    Double_t q1[n_det];
    Double_t t1[n_det];
    Double_t x2[n_det];
    Double_t y2[n_det];
    Double_t z2[n_det];
    Double_t q2[n_det];
    Double_t t2[n_det];
    Double_t xMax[n_det];
    Double_t yMax[n_det];
    Double_t zMax[n_det];
    Double_t qMax[n_det];
    Double_t tMax[n_det];

    Double_t xTrack[n_det];
    Double_t yTrack[n_det];
    Double_t zTrack[n_det];
    Double_t qTrack[n_det];

    Bool_t pat1[2 * n_det];
    Bool_t pat2[2 * n_det];

    Int_t id, id1, id2;

    Int_t det = 0;
    Int_t det1 = 0;
    Int_t det2 = 0;

    Double_t track1_best[6]; // x,y,z, px, py, pz
    Double_t track2_best[6]; // x,y,z, px, py, pz
    Double_t chi_best[6];    // chi2, chi2_red

    Double_t chi2;

    for (int i = 0; i < n_det; i++)
    {
        x[i] = 0.;
        y[i] = 0.;
        z[i] = 0.;
        q[i] = 0.;
        t[i] = -1000.;

        x1[i] = 0.;
        y1[i] = 0.;
        z1[i] = 0.;
        q1[i] = 0.;
        t1[i] = -1000.;

        x2[i] = 0.;
        y2[i] = 0.;
        z2[i] = 0.;
        q2[i] = 0.;
        t2[i] = -1000.;

        xMax[i] = -1000.;
        yMax[i] = -1000.;
        zMax[i] = -1000.;
        qMax[i] = -1000.;
        tMax[i] = -1000.;
    }

    // is also number of ifibcount
    Int_t fi3a = 0;
    Int_t fi3b = 1;
    Int_t fi10 = 2;
    Int_t fi11 = 3;
    Int_t fi12 = 4;
    Int_t fi13 = 5;
    Int_t tofd1r = 6;
    Int_t tofd1l = 7;
    Int_t tofd2r = 8;
    Int_t tofd2l = 9;
    Int_t ghost = 10;

    Double_t tof = 0.;
    Bool_t single = false;
    Double_t tStart = -1e9;
    Bool_t first = true;
    Bool_t alpha = false;
    Bool_t carbon = false;

    for (Int_t i = 0; i < 10; i++)
    {
        tPrev[i] = -1000.;
        detPrev[i] = -1;
    }


    auto detTofd = fHitItems.at(DET_TOFD);
    Int_t nHits = detTofd->GetEntriesFast();
    // cout << "ToFD hits: " << nHits << endl;

    if (nHits > 0)
    {
        counterTofd++;
    }

    if (nHits > 100)
        return;

    Int_t multTofd = 0;

    //===========================================================================
    // loop over ToFD
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        R3BTofdHitData* hitTofd = (R3BTofdHitData*)detTofd->At(ihit);

        if (IS_NAN(hitTofd->GetTime()))
            continue;

        Double_t ttt = hitTofd->GetTime();
        //if (fCuts && (ttt < -100. || ttt > 100.) && !fSimu) // trigger window -1500, 1500
        if (fCuts && (ttt < -20. || ttt > 20.) && !fSimu) // trigger window -1500, 1500
            continue;

        Double_t qqq = hitTofd->GetEloss();
        Double_t xxx = hitTofd->GetX();
        Double_t yyy = hitTofd->GetY();
        Double_t y_corr = 0.;

        // first looking for the right charge
        if (fB == -1102 && !fSimu)
        {
            //            if (qqq < 10. || qqq > 14.)
            //                continue;
            if (xxx >= 46. && xxx < 47.)
                y_corr = 0.08;
            if (xxx >= 47. && xxx < 48.)
                y_corr = 0.04;
            if (xxx >= 49. && xxx < 50.)
                y_corr = 0.05;
        }
        if (fB == -1292 && !fSimu)
        {
            //            if (qqq < 8.5 || qqq > 10.)
            //                continue;
            if (xxx > 30. && xxx < 31.)
                y_corr = 0.03;
            if (xxx >= 32. && xxx < 33.)
                y_corr = -0.16;
            if (xxx >= 33. && xxx < 34.)
                y_corr = 0.0;
            if (xxx >= 35. && xxx < 36.)
                y_corr = 0.12;
        }
        if (fB == -1482 && !fSimu)
        {
            if (qqq < 1. || qqq > 9.)
                continue;
            if (xxx >= 15. && xxx < 16.)
                y_corr = -0.01;
            if (xxx >= 16. && xxx < 17.)
                y_corr = 0.00;
            if (xxx >= 18. && xxx < 19.)
                y_corr = -0.04;
        }
        if (fB == -1862 && !fSimu)
        {
            //            if (qqq < 6.5 || qqq > 8.)
            //                continue;
            if (xxx >= -17. && xxx < -16.)
                y_corr = 0.12;
            if (xxx >= -16. && xxx < -15.)
                y_corr = 0.08;
            if (xxx >= -15. && xxx < -14.)
                y_corr = 0.0;
        }
        if (fB == -2052 && !fSimu)
        {
            //            if (qqq < 6. || qqq > 7.)
            //                continue;
            if (xxx >= -34. && xxx < -33.)
                y_corr = 0.07;
            if (xxx >= -33. && xxx < -32.)
                y_corr = 0.07;
            if (xxx >= -32. && xxx < -31.)
                y_corr = -0.06;
            if (xxx >= -31. && xxx < -30.)
                y_corr = -0.05;
        }
        if (fB == -2242 && !fSimu)
        {
            //            if (qqq < 5.8 || qqq > 6.6)
            //                continue;
            if (xxx >= -49. && xxx < -48.)
                y_corr = 0.05;
            if (xxx >= -48. && xxx < -47.)
                y_corr = -0.035;
            if (xxx >= -47. && xxx < -46.)
                y_corr = 0.04;
        }

        if (fB == -1672)
        {

            //if (!fPairs && (qqq < 7.5 || qqq > 8.5))
            //    continue;

            //if (fPairs && !(qqq > 1.5 && qqq < 2.5) && !(qqq > 5.5 && qqq < 6.5))
            //    continue;

            y_corr = 0.0;
        }

        //======= Determine which TOFD is hit ========
        id2 = hitTofd->GetDetId();
        if (hitTofd->GetX() / 100. <= 0.)
        {
            // tof rechts
            if (id2 == 1)
            {
                det2 = tofd1r;
            }
            else if (id2 == 2)
            {
                det2 = tofd2r;
            }
        }
        else
        {
            // tof links
            if (id2 == 1)
            {
                det2 = tofd1l;
            }
            else if (id2 == 2)
            {
                det2 = tofd2l;
            }
        }

        x2[det2] = hitTofd->GetX() / 100.;
        y2[det2] = hitTofd->GetY() / 100. + y_corr;
        z2[det2] = 0.;
        q2[det2] = hitTofd->GetEloss();
        t2[det2] = hitTofd->GetTime();

        // Achtung, ändern
        //if (!fPairs)
        //    q2[det2] = 8.;

        // register hits for tracker as long a time is in the coincidence window
        if ((fabs(t2[det2] - t1[det1]) < 2.) || first)
        {
            // register point for tracker
            detector[countdet] = det2;
            xdet[countdet] = x2[det2];
            ydet[countdet] = y2[det2];
            zdet[countdet] = z2[det2];
            //qdet[countdet] = (int)(q2[det2] + 0.5);
            qdet[countdet] = q2[det2];
            tdet[countdet] = t2[det2];

            countdet++;

            //if(det2==7)
            tStart = t2[det2];

            if (fabs(qdet[countdet] - 2.) < 0.5)
                alpha = true;
            if (fabs(qdet[countdet] - 6.) < 0.5)
                carbon = true;

            single = true;
            first = false;



            det1 = det2;
            x1[det1] = x2[det2];
            y1[det1] = y2[det2];
            z1[det1] = 0.;
            q1[det1] = q2[det2];
            t1[det1] = t2[det2];
            id1 = id2;

            // since we had a coincidence, continue with next event, if not last event.
            if (ihit < nHits - 1)
                continue;
        }

        if (!single)
            continue;

        if (fPairs && !(alpha && carbon))
            continue;

        alpha = false;
        carbon = false;
        single = false;
        first = true;

        if (ihit < nHits - 1)
            ihit--;

        det1 = det2;
        x1[det1] = x2[det2];
        y1[det1] = y2[det2];
        z1[det1] = 0.;
        q1[det1] = q2[det2];
        t1[det1] = t2[det2];
        id1 = id2;

        multTofd++;

        counterTofdMulti++;

    }//=========================== end loop on TofD hits ========================

    if(multTofd==0) return;

    // cut in ToT for Fibers
    Bool_t maxWerte = false;

    Double_t cutQ = 0.;

    if (!fPairs || fB != -1672)
        cutQ = 4.;


    //===========================================================================
    // loop over fiber 13
    auto detHit13 = fHitItems.at(DET_FI13);
    Int_t nHits13 = detHit13->GetEntriesFast();
    LOG(DEBUG) << "Fi13 hits: " << nHits13 << endl;

    Int_t mult13 = 0;
    for (Int_t ihit13 = 0; ihit13 < nHits13; ihit13++)
    {
        det = fi13;
        R3BBunchedFiberHitData* hit13 = (R3BBunchedFiberHitData*)detHit13->At(ihit13);
        x1[det] = hit13->GetX() / 100.;
        y1[det] = hit13->GetY() / 100.;
        z1[det] = 0.;
        q1[det] = hit13->GetEloss();
        t1[det] = hit13->GetTime();
        tof = tStart - t1[det];

        LOG(DEBUG2) << "Fi13 bc: " << ihit13 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
            << " t1: " << t1[det]
            << " tof " << tof << endl;

        //-------------------- Set Cuts -----------------------

        if (fCuts && x1[det] * 100. < -24.4)
            continue;

        if (fCuts && (tof < -20 || tof > 20))
            continue;

        if (fCuts && y1[det] < -100.)
            continue;

        //if (fCuts && !fPairs && q1[det] < cutQ)
        //    continue;
        //if (fGraphCuts && !cut_Fi13vsTofd->IsInside(x1[tofd1r] * 100., x1[det] * 100.))
        //    continue;

        mult13++;

        if (q1[det] > qMax[det])
        {
            qMax[det] = q1[det];
            xMax[det] = x1[det];
            yMax[det] = y1[det];
            zMax[det] = z1[det];
            //tMax[det] = t1[det];
            tMax[det] = tof;
        }

        if (!maxWerte)
        {
            detector[countdet] = det;
            xdet[countdet] = x1[det];// + gRandom->Uniform(-0.00025, 0.00025);
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = q1[det];
            //tdet[countdet] = t1[det];
            tdet[countdet] = tof;
            countdet++;
        }
    }

    if (maxWerte && mult13>0)
    {
        detector[countdet] = det;
        xdet[countdet] = xMax[det];// + gRandom->Uniform(-0.00025, 0.00025);
        ydet[countdet] = yMax[det];
        zdet[countdet] = zMax[det];
        qdet[countdet] = qMax[det];
        tdet[countdet] = tMax[det];
        countdet++;
    }


    //===========================================================================
    // loop over fiber 11
    auto detHit11 = fHitItems.at(DET_FI11);
    Int_t nHits11 = detHit11->GetEntriesFast();
    LOG(DEBUG) << "Fi11 hits: " << nHits11 << endl;
    Int_t mult11 = 0;
    for (Int_t ihit11 = 0; ihit11 < nHits11; ihit11++)
    {
        det = fi11;
        R3BBunchedFiberHitData* hit11 = (R3BBunchedFiberHitData*)detHit11->At(ihit11);
        x1[det] = hit11->GetX() / 100.;
        y1[det] = hit11->GetY() / 100.;
        z1[det] = 0.;
        q1[det] = hit11->GetEloss();
        t1[det] = hit11->GetTime();
        tof = tStart - t1[det];

        LOG(DEBUG2) << "Fi11 bc: " << ihit11 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
            << " t1: " << t1[det] << endl;

        //-------------------- Set Cuts -----------------------

        if (fCuts && x1[det] * 100. < -24.4)
            continue;

        if (fCuts && (tof < -20 || tof > 20))
            continue;

        if (fCuts && y1[det] < -100.)
            continue;

        //if (fCuts && !fPairs && q1[det] < cutQ)
        //    continue;
        //if (fGraphCuts && !cut_Fi13vsTofd->IsInside(x1[tofd1r] * 100., x1[det] * 100.))
        //    continue;

        mult11++;

        if (q1[det] > qMax[det])
        {
            qMax[det] = q1[det];
            xMax[det] = x1[det];
            yMax[det] = y1[det];
            zMax[det] = z1[det];
            //tMax[det] = t1[det];
            tMax[det] = tof;
        }

        LOG(DEBUG2) << "Fi11: " << ihit11 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
            << " t1: " << t1[det] << endl;

        if (!maxWerte)
        {
            detector[countdet] = det;
            xdet[countdet] = x1[det];// + gRandom->Uniform(-0.00025, 0.00025);
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = q1[det];
            //tdet[countdet] = t1[det];
            tdet[countdet] = tof;
            countdet++;
        }
    }
    if (maxWerte && mult11>0)
    {
        detector[countdet] = det;
        xdet[countdet] = xMax[det];// + gRandom->Uniform(-0.00025, 0.00025);
        ydet[countdet] = yMax[det];
        zdet[countdet] = zMax[det];
        qdet[countdet] = qMax[det];
        tdet[countdet] = tMax[det];
        countdet++;
    }


    //===========================================================================
    // loop over fiber 10
    auto detHit10 = fHitItems.at(DET_FI10);
    Int_t nHits10 = detHit10->GetEntriesFast();
    LOG(DEBUG) << "Fi10 hits: " << nHits10 << endl;
    Int_t mult10 = 0;
    for (Int_t ihit10 = 0; ihit10 < nHits10; ihit10++)
    {
        det = fi10;
        R3BBunchedFiberHitData* hit10 = (R3BBunchedFiberHitData*)detHit10->At(ihit10);
        x1[det] = hit10->GetX() / 100.;
        y1[det] = hit10->GetY() / 100.;
        z1[det] = 0.;
        q1[det] = hit10->GetEloss();
        t1[det] = hit10->GetTime();
        tof = tStart - t1[det];

        LOG(DEBUG2) << "Fi10: " << ihit10 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
            << " t1: " << t1[det] << endl;

        //-------------------- Set Cuts -----------------------

        if (fCuts && x1[det] * 100. < -24.4)
            continue;

        if (fCuts && (tof < -20 || tof > 20))
            continue;

        if (fCuts && y1[det] < -100.)
            continue;

        //if (fCuts && !fPairs && q1[det] < cutQ)
        //    continue;
        //if (fGraphCuts && !cut_Fi13vsTofd->IsInside(x1[tofd1r] * 100., x1[det] * 100.))
        //    continue;

        mult10++;

        if (q1[det] > qMax[det])
        {
            qMax[det] = q1[det];
            xMax[det] = x1[det];
            yMax[det] = y1[det];
            zMax[det] = z1[det];
            //tMax[det] = t1[det];
            tMax[det] = tof;
        }

        if (!maxWerte)
        {
            detector[countdet] = det;
            xdet[countdet] = x1[det];//+ gRandom->Uniform(-0.0005, 0.0005);
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = q1[det];
            //tdet[countdet] = t1[det];
            tdet[countdet] = tof;
            countdet++;
        }
    }
    if (maxWerte && mult10>0)
    {
        detector[countdet] = det;
        xdet[countdet] = xMax[det];//+ gRandom->Uniform(-0.0005, 0.0005);
        ydet[countdet] = yMax[det];
        zdet[countdet] = zMax[det];
        qdet[countdet] = qMax[det];
        tdet[countdet] = tMax[det];
        countdet++;
    }


    //===========================================================================
    // loop over fiber 12
    auto detHit12 = fHitItems.at(DET_FI12);
    Int_t nHits12 = detHit12->GetEntriesFast();
    LOG(DEBUG) << "Fi12 hits: " << nHits12 << endl;
    Int_t mult12 = 0;
    for (Int_t ihit12 = 0; ihit12 < nHits12; ihit12++)
    {
        det = fi12;
        R3BBunchedFiberHitData* hit12 = (R3BBunchedFiberHitData*)detHit12->At(ihit12);
        x1[det] = hit12->GetX() / 100.;
        y1[det] = hit12->GetY() / 100.;
        z1[det] = 0.;
        q1[det] = hit12->GetEloss();
        t1[det] = hit12->GetTime();
        tof = tStart - t1[det];

        LOG(DEBUG2) << "Fi12: " << ihit12 << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
            << " t1: " << t1[det] << endl;


        //-------------------- Set Cuts -----------------------

        if (fCuts && x1[det] * 100. < -24.4)
            continue;

        if (fCuts && (tof < -20 || tof > 20))
            continue;

        if (fCuts && y1[det] < -100.)
            continue;

        //if (fCuts && !fPairs && q1[det] < cutQ)
        //    continue;
        //if (fGraphCuts && !cut_Fi13vsTofd->IsInside(x1[tofd1r] * 100., x1[det] * 100.))
        //    continue;

        mult12++;

        if (q1[det] > qMax[det])
        {
            qMax[det] = q1[det];
            xMax[det] = x1[det];
            yMax[det] = y1[det];
            zMax[det] = z1[det];
            //tMax[det] = t1[det];
            tMax[det] = tof;
        }

        if (!maxWerte)
        {
            detector[countdet] = det;
            xdet[countdet] = x1[det];// + gRandom->Uniform(-0.0005, 0.0005);
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = q1[det];
            //tdet[countdet] = t1[det];
            tdet[countdet] = tof;
            countdet++;
        }
    }
    if (maxWerte && mult12>0)
    {
        detector[countdet] = det;
        xdet[countdet] = xMax[det];//+ gRandom->Uniform(-0.0005, 0.0005);
        ydet[countdet] = yMax[det];
        zdet[countdet] = zMax[det];
        qdet[countdet] = qMax[det];
        tdet[countdet] = tMax[det];
        countdet++;
    }


    //===========================================================================
    // loop over fiber 3a
    auto detHit3a = fHitItems.at(DET_FI3A);
    Int_t nHits3a = detHit3a->GetEntriesFast();
    LOG(DEBUG) << "Fi3a hits: " << nHits3a << endl;
    Int_t mult3a = 0;
    for (Int_t ihit3a = 0; ihit3a < nHits3a; ihit3a++)
    {
        det = fi3a;
        R3BBunchedFiberHitData* hit3a = (R3BBunchedFiberHitData*)detHit3a->At(ihit3a);
        x1[det] = hit3a->GetX() / 100.;
        y1[det] = hit3a->GetY() / 100.;
        z1[det] = 0.;
        q1[det] = hit3a->GetEloss();
        t1[det] = hit3a->GetTime();

        //if (q1[det] > 5.)
        //    q1[det] = 8.;
        //else if (q1[det] > 2.)
        //    q1[det] = 6.;
        //else if (q1[det] > 0.)
        //    q1[det] = 2.;

        tof = tStart - t1[det];

        // Cuts on Fi3a
        //            if (fCuts && y1[det] * 100. > 50.)
        //                continue;
        //            if (fCuts && y1[det] *100. < -50.)
        //                continue;
        //if (fCuts && !fPairs && q1[det] < cutQ)
        //    continue;
        //if (fCuts && (tof < -40 || tof > 20))
        if (fCuts && (tof < -20 || tof > 20))
            continue;

        mult3a++;
        //if (mult3a > 10)
        //    continue;

        if (q1[det] > qMax[det])
        {
            qMax[det] = q1[det];
            xMax[det] = x1[det];
            yMax[det] = y1[det];
            zMax[det] = z1[det];
            tMax[det] = t1[det];
        }

        LOG(DEBUG2) << "Fi3a " << ihit3a << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
            << " t1: " << tof << endl;

        if (!maxWerte)
        {
            detector[countdet] = det;
            xdet[countdet] = x1[det];
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = q1[det];
            //tdet[countdet] = t1[det];
            tdet[countdet] = tof;
            countdet++;
        }
    }
    if (maxWerte && mult3a>0)
    {
        detector[countdet] = det;
        xdet[countdet] = xMax[det];
        ydet[countdet] = yMax[det];
        zdet[countdet] = zMax[det];
        qdet[countdet] = qMax[det];
        tdet[countdet] = tMax[det];
        countdet++;
    }


    //===========================================================================
    // loop over fiber 3b
    auto detHit3b = fHitItems.at(DET_FI3B);
    Int_t nHits3b = detHit3b->GetEntriesFast();
    LOG(DEBUG) << "Fi3b hits: " << nHits3b << endl;
    Int_t mult3b = 0;

    for (Int_t ihit3b = 0; ihit3b < nHits3b; ihit3b++)
    {
        det = fi3b;
        R3BBunchedFiberHitData* hit3b = (R3BBunchedFiberHitData*)detHit3b->At(ihit3b);
        x1[det] = hit3b->GetX() / 100.;
        y1[det] = hit3b->GetY() / 100.;
        z1[det] = 0.;
        q1[det] = hit3b->GetEloss();
        t1[det] = hit3b->GetTime();

        //if (q1[det] > 5.)
        //    q1[det] = 8.;
        //else if (q1[det] > 2.)
        //    q1[det] = 6.;
        //else if (q1[det] > 0.)
        //    q1[det] = 2.;

        tof = tStart - t1[det];

        // Cuts on Fi3b
        //            if (fCuts && y1[det] * 100. > 50.)
        //                continue;
        //            if (fCuts && y1[det] * 100. < -50.)
        //                continue;
        //if (fCuts && !fPairs && q1[det] < cutQ)
        //    continue;
        //if (fCuts && (tof < -40 || tof > 20))
        //    continue;
        if (fCuts && (tof < -20 || tof > 20))
            //if (fCuts && (tof < 0 || tof > 7))
            continue;



        mult3b++;
        //if (mult3b > 10)
        //    continue;

        if (q1[det] > qMax[det])
        {
            qMax[det] = q1[det];

            xMax[det] = x1[det];
            yMax[det] = y1[det];
            zMax[det] = z1[det];
            tMax[det] = t1[det];
        }

        LOG(DEBUG2) << "Fi3b " << ihit3b << " x1: " << x1[det] << " y1: " << y1[det] << " q1: " << q1[det]
            << " t1: " << tof << endl;

        if (!maxWerte)
        {
            detector[countdet] = det;
            xdet[countdet] = x1[det];
            ydet[countdet] = y1[det];
            zdet[countdet] = z1[det];
            qdet[countdet] = q1[det];
            //tdet[countdet] = t1[det];
            tdet[countdet] = tof;
            countdet++;
        }
    }
    if (maxWerte && mult3b>0)
    {
        detector[countdet] = det;
        xdet[countdet] = xMax[det];
        ydet[countdet] = yMax[det];
        zdet[countdet] = zMax[det];
        qdet[countdet] = qMax[det];
        tdet[countdet] = tMax[det];
        countdet++;
    }

    //========== TRACKER ==============================
    Current =0;
    N_glob_tracks_left = 0;
    N_glob_tracks_right = 0;
    bool is_used;

    if (tracker && 
             mult10>0 && mult12>0 && multTofd>0 && mult3b>0
              //(mult11>0 && mult13>0 && multTofd>0 && mult3a>0) )
       )
    {
        counter2++;

        Nhits_f3a  = mult3a; 
        Nhits_f3b  = mult3b;
        Nhits_f10  = mult10;
        Nhits_f11  = mult11;
        Nhits_f12  = mult12;
        Nhits_f13  = mult13;
        Nhits_tofd = multTofd;

        Make_Lab_Coordinates(countdet, detector, qdet, tdet, xdet, ydet, zdet);

        cout << "\n============= Tracking Left Arm ============== ";
        //============== Tracking left arm ===================
        PoQ_data.clear();
        PoQ_data_filtered.clear();
        Track_MDF(1);//1 - left arm (F10+F12+TOFD+F3B)
        for(auto & _poq : PoQ_data) //Filtering collected data
        {
            if( !(_poq.is_left_arm) ) continue;
            is_used = false;
            for(auto & _poq_filter : PoQ_data_filtered) //check if already used
            {
                if( Used_Track_Left(_poq, _poq_filter) )
                {   is_used = true;  break;  }
            }
            if(is_used) continue;
            //If good unused track, add it to the filtered list
            if( Good_Track_Left(_poq) ){  PoQ_data_filtered.push_back(_poq); }
            else continue;

            //Loop remaining data to see if there is anything beter
            for(auto & _poq_check : PoQ_data)
            {
                if( !(_poq_check.is_left_arm) || Same_Track_Left(_poq, _poq_check)
                        || !Good_Track_Left(_poq_check) ) continue;

                if( Used_Track_Left(_poq, _poq_check) && 
                        sqrt(pow(_poq_check.tx0_res, 2) + pow(_poq_check.tx0_f3_res,2) ) <
                        sqrt(pow(_poq.tx0_res, 2) + pow(_poq.tx0_f3_res,2)) )
                {
                    cout << "\nReplacing by a better track...";
                    PoQ_data_filtered.pop_back();
                    PoQ_data_filtered.push_back(_poq_check);
                }
            }
        }//end for loop filtering left arm
            
        cout << "\n Number of al potential tracks = " << PoQ_data.size() ;
        cout << "\n NUmber tracks after filtering = " << PoQ_data_filtered.size();


        //============== Collecting tracked data from the left arm===================
        for(auto & _poq_data : PoQ_data_filtered)
        {
            if(!_poq_data.is_left_arm) continue;
            f3b_X[N_glob_tracks_left] = _poq_data.f3b_hit.X;
            f3b_Y[N_glob_tracks_left] = _poq_data.f3b_hit.Y;
            f3b_Z[N_glob_tracks_left] = _poq_data.f3b_hit.Z;
            f3b_Q[N_glob_tracks_left] = _poq_data.f3b_hit.Q;
            f3b_T[N_glob_tracks_left] = _poq_data.f3b_hit.T;

            f10_X[N_glob_tracks_left] = _poq_data.f10_hit.X;
            f10_Y[N_glob_tracks_left] = _poq_data.f10_hit.Y;
            f10_Z[N_glob_tracks_left] = _poq_data.f10_hit.Z;
            f10_Q[N_glob_tracks_left] = _poq_data.f10_hit.Q;
            f10_T[N_glob_tracks_left] = _poq_data.f10_hit.T;

            f12_X[N_glob_tracks_left] = _poq_data.f12_hit.X;
            f12_Y[N_glob_tracks_left] = _poq_data.f12_hit.Y;
            f12_Z[N_glob_tracks_left] = _poq_data.f12_hit.Z;
            f12_Q[N_glob_tracks_left] = _poq_data.f12_hit.Q;
            f12_T[N_glob_tracks_left] = _poq_data.f12_hit.T;

            tofd_left_X[N_glob_tracks_left] = _poq_data.dtof_hit.X;
            tofd_left_Y[N_glob_tracks_left] = _poq_data.dtof_hit.Y;
            tofd_left_Z[N_glob_tracks_left] = _poq_data.dtof_hit.Z;
            tofd_left_Q[N_glob_tracks_left] = _poq_data.dtof_hit.Q;
            tofd_left_T[N_glob_tracks_left] = _poq_data.dtof_hit.T;

            TX0_mdf_left[N_glob_tracks_left]    = _poq_data.value_tx0;
            TX0_f3_mdf_left[N_glob_tracks_left] = _poq_data.value_tx0_f3;
            PoQ_mdf_left[N_glob_tracks_left]    = _poq_data.value_poq;

            X0_residual_left[N_glob_tracks_left]     = _poq_data.x0_res;
            TX0_residual_left[N_glob_tracks_left]    = _poq_data.tx0_res;
            TX0_residual_f3_left[N_glob_tracks_left] = _poq_data.tx0_f3_res;
            TX0_residual_mdf_left[N_glob_tracks_left] = _poq_data.tx0_mdf_res;
            X0_proj_by_f3_left[N_glob_tracks_left]  = _poq_data.x0_proj_by_f3b;
            TX0_proj_by_f3_left[N_glob_tracks_left] = _poq_data.tx0_f3;
            Xoffset_f3b[N_glob_tracks_left] = 
                _poq_data.f3b_hit.Z*(TX0_mdf_left[N_glob_tracks_left]- TX0_proj_by_f3_left[N_glob_tracks_left]) ;

            N_glob_tracks_left++;

        }//end filling tree variables from filtered data of left arm

        ////============== Tracking right arm ===================
        //PoQ_data.clear();
        //PoQ_data_filtered.clear();
        //Track_MDF(2);//2 = right arm (F11+F13+TOFD+F3a)
        //for(auto & _poq : PoQ_data) //Filtering collected data
        //{
        //    if( !(_poq.is_right_arm) ) continue;
        //    is_used = false;
        //    for(auto & _poq_filter : PoQ_data_filtered) //check if already used
        //    {
        //        if( Used_Track_Right(_poq, _poq_filter) )
        //        {   is_used = true;  break;  }
        //    }
        //    if(is_used) continue;
        //    //If good unused track, add it to the filtered list
        //    if( Good_Track_Right(_poq) ){  PoQ_data_filtered.push_back(_poq); }
        //    else continue;

        //    //Loop remaining data to see if there is anything beter
        //    for(auto & _poq_check : PoQ_data)
        //    {
        //        if( !(_poq_check.is_right_arm) || Same_Track_Right(_poq, _poq_check)
        //                || !Good_Track_Right(_poq_check) ) continue;

        //        if( Used_Track_Right(_poq, _poq_check) && 
        //                sqrt(pow(_poq_check.tx0_res, 2) + pow(_poq_check.tx0_f3_res,2) ) <
        //                sqrt(pow(_poq.tx0_res, 2) + pow(_poq.tx0_f3_res,2)) )
        //        {
        //            PoQ_data_filtered.pop_back();
        //            PoQ_data_filtered.push_back(_poq_check);
        //        }
        //    }
        //}//end for loop filtering right arm

        ////============== Collecting tracked data from the right arm===================
        //for(auto & _poq_data : PoQ_data_filtered)
        //{
        //    if(!_poq_data.is_right_arm) continue;
        //    f3a_X[N_glob_tracks_right] = _poq_data.f3a_hit.X;
        //    f3a_Y[N_glob_tracks_right] = _poq_data.f3a_hit.Y;
        //    f3a_Z[N_glob_tracks_right] = _poq_data.f3a_hit.Z;
        //    f3a_Q[N_glob_tracks_right] = _poq_data.f3a_hit.Q;
        //    f3a_T[N_glob_tracks_right] = _poq_data.f3a_hit.T;

        //    f11_X[N_glob_tracks_right] = _poq_data.f11_hit.X;
        //    f11_Y[N_glob_tracks_right] = _poq_data.f11_hit.Y;
        //    f11_Z[N_glob_tracks_right] = _poq_data.f11_hit.Z;
        //    f11_Q[N_glob_tracks_right] = _poq_data.f11_hit.Q;
        //    f11_T[N_glob_tracks_right] = _poq_data.f11_hit.T;

        //    f13_X[N_glob_tracks_right] = _poq_data.f13_hit.X;
        //    f13_Y[N_glob_tracks_right] = _poq_data.f13_hit.Y;
        //    f13_Z[N_glob_tracks_right] = _poq_data.f13_hit.Z;
        //    f13_Q[N_glob_tracks_right] = _poq_data.f13_hit.Q;
        //    f13_T[N_glob_tracks_right] = _poq_data.f13_hit.T;

        //    tofd_right_X[N_glob_tracks_right] = _poq_data.dtof_hit.X;
        //    tofd_right_Y[N_glob_tracks_right] = _poq_data.dtof_hit.Y;
        //    tofd_right_Z[N_glob_tracks_right] = _poq_data.dtof_hit.Z;
        //    tofd_right_Q[N_glob_tracks_right] = _poq_data.dtof_hit.Q;
        //    tofd_right_T[N_glob_tracks_right] = _poq_data.dtof_hit.T;

        //    TX0_mdf_right[N_glob_tracks_right]    = _poq_data.value_tx0;
        //    TX0_f3_mdf_right[N_glob_tracks_right] = _poq_data.value_tx0_f3;
        //    PoQ_mdf_right[N_glob_tracks_right]    = _poq_data.value_poq;

        //    X0_residual_right[N_glob_tracks_right]     = _poq_data.x0_res;

        //    TX0_residual_right[N_glob_tracks_right]    = _poq_data.tx0_res;
        //    TX0_residual_f3_right[N_glob_tracks_right] = _poq_data.tx0_f3_res;
        //    TX0_residual_mdf_right[N_glob_tracks_right] = _poq_data.tx0_mdf_res;

        //    X0_proj_by_f3_right[N_glob_tracks_right]  = _poq_data.x0_proj_by_f3a;
        //    TX0_proj_by_f3_right[N_glob_tracks_right] = _poq_data.tx0_f3;
        //    Xoffset_f3a[N_glob_tracks_right] =
        //        _poq_data.f3a_hit.Z*(TX0_mdf_right[N_glob_tracks_right]- TX0_proj_by_f3_right[N_glob_tracks_right]) ;

        //    N_glob_tracks_right++;
        //}//end filling tree variables from filtered data of right arm

        Current = fB;
        tree_out.Fill();
    }//end tracker
    return;
}

void R3BTrackS454_MDF::Output1(Double_t track[12], Double_t chi[6])
{
}

void R3BTrackS454_MDF::Output2(Double_t track_parameter[12], Double_t chi_single_parameter[6])
{
}


void R3BTrackS454_MDF::FinishEvent()
{


    fTrackItems->Clear();
    fNofTrackItems = 0;

    for (Int_t det = 0; det < DET_MAX; det++)
    {
        if (fMappedItems.at(det))
        {
            fMappedItems.at(det)->Clear();
        }
        if (fCalItems.at(det))
        {
            fCalItems.at(det)->Clear();
        }
        if (fHitItems.at(det))
        {
            fHitItems.at(det)->Clear();
        }
    }




}

void R3BTrackS454_MDF::FinishTask()
{

    cout << "Statistics:" << endl;
    cout << "Events: " << fNEvents << endl;
    cout << "Wrong Trigger: " << counterWrongTrigger << endl;
    cout << "Wrong Tpat: " << counterWrongTpat << endl;
    cout << "ROLU veto: " << counterRolu << endl;
    cout << "Califa veto: " << counterCalifa << endl;
    cout << "TofD: " << counterTofd << endl;
    cout << "TofD multi: " << counterTofdMulti << endl;
    cout << "Tracker: " << counterTracker << endl;

    tree_out.Write();
}

void R3BTrackS454_MDF::Detector_Hit::Set_XYZQ(Double_t _x, Double_t _y, Double_t _z, Double_t _q, Double_t _t, Int_t _det)
{
    X =_x; Y=_y; Z=_z; Q=_q; T=_t; Detector = _det;
    return;
}

void R3BTrackS454_MDF::Detector_Hit::Set_XYZ_det(Double_t _x, Double_t _y, Double_t _z)
{
    Xdet =_x; Ydet=_y; Zdet=_z;
    return;
}

bool R3BTrackS454_MDF::Same_Track_Right(MDF_Data_PoQ p1, MDF_Data_PoQ p2)
{
    bool decision = false;
    if(p1.f11_hit.X      == p2.f11_hit.X 
            &&  p1.f13_hit.X  == p2.f13_hit.X
            &&  p1.dtof_hit.X == p2.dtof_hit.X
            &&  p1.f3a_hit.X  == p2.f3a_hit.X)
        decision = true;
    return decision;
}

bool R3BTrackS454_MDF::Used_Track_Right(MDF_Data_PoQ p1, MDF_Data_PoQ p2)
{
    bool decision = false;
    if(p1.f11_hit.X      == p2.f11_hit.X 
            ||  p1.f13_hit.X  == p2.f13_hit.X
            ||  p1.dtof_hit.X == p2.dtof_hit.X
            ||  p1.f3a_hit.X  == p2.f3a_hit.X)
        decision = true;
    return decision;
}

bool R3BTrackS454_MDF::Good_Track_Right(MDF_Data_PoQ p1)
{
    bool decision = false;
    if( 
            fabs(p1.tx0_mdf_res)   < Res_Limit
            //&& fabs(p1.tx0_res)    < Res_Limit
            //&& fabs(p1.tx0_f3_res) < Res_Limit
            && p1.value_tx0 < 0
            && p1.value_tx0_f3 < 0
            && p1.tx0_f3 < 0   
            && p1.x0_proj_by_f3a > -0.4 && p1.x0_proj_by_f3a <0.4
            && p1.dtof_hit.Q < 8.5
      )
        decision = true;
    return decision;
}

bool R3BTrackS454_MDF::Same_Track_Left(MDF_Data_PoQ p1, MDF_Data_PoQ p2)
{
    bool decision = false;
    if(p1.f10_hit.X      == p2.f10_hit.X 
            &&  p1.f12_hit.X  == p2.f12_hit.X
            &&  p1.dtof_hit.X == p2.dtof_hit.X
            &&  p1.f3b_hit.X  == p2.f3b_hit.X)
        decision = true;
    return decision;
}

bool R3BTrackS454_MDF::Used_Track_Left(MDF_Data_PoQ p1, MDF_Data_PoQ p2)
{
    bool decision = false;
    if(p1.f10_hit.X      == p2.f10_hit.X 
            ||  p1.f12_hit.X  == p2.f12_hit.X
            ||  p1.dtof_hit.X == p2.dtof_hit.X
            ||  p1.f3b_hit.X  == p2.f3b_hit.X)
        decision = true;
    return decision;
}


bool R3BTrackS454_MDF::Good_Track_Left(MDF_Data_PoQ p1)
{
    bool decision = false;
    if( 
            fabs(p1.tx0_mdf_res)   < Res_Limit
            && fabs(p1.tx0_res)    < Res_Limit
            && fabs(p1.tx0_f3_res) < Res_Limit
            && p1.value_tx0 > 0
            && p1.value_tx0_f3 > 0
            && p1.tx0_f3 > 0
            && p1.x0_proj_by_f3b > -0.4 && p1.x0_proj_by_f3b <0.4
            && p1.dtof_hit.Q < 8.5
      )
        decision = true;
    return decision;
}

void R3BTrackS454_MDF::Make_Lab_Coordinates(Int_t detcount, Int_t* det, Double_t* qd, Double_t* td, Double_t * xd, Double_t * yd, Double_t * zd)
{
    //Make sure that hit containers are empty
    f3a_hits.clear();
    f3b_hits.clear();
    f10_hits.clear();
    f11_hits.clear();
    f12_hits.clear();
    f13_hits.clear();
    dtof_hits.clear();

    //Temporary containers
    Detector_Hit det_hit;

    //At first convert each detector hit to the lab coordinates (in cm as r3broot)
    for(int i=0; i<detcount; i++)
    {
        if(det[i] == 0) //=========== fi3a
        {
            det_hit.Set_XYZQ(
                    xd[i]*(-100.) - 0.25 -  f3ab_halfwidth + 0.3724 + F1113_X0par + F3b_dX,// 5 mm slit, flipped
                    yd[i]*100.,
                    zd[i]+80.5,
                    qd[i], td[i], det[i]);
            f3a_hits.push_back(det_hit);
        }
        else if(det[i] == 1) //========== fi3b
        {
            det_hit.Set_XYZQ(
                    xd[i]*100. + 0.25 + f3ab_halfwidth + 0.3724 + F1012_X0par + F3b_dX, //5mm slit
                    yd[i]*100.,
                    zd[i]+80.5,
                    qd[i], td[i], det[i]);
            det_hit.Set_XYZ_det(xd[i]*100., yd[i]*100., zd[i]*100);//saving original det coordiantes
            f3b_hits.push_back(det_hit);
        }
        else if(det[i] == 2) //============ fi10
        {
            det_hit.Set_XYZQ(
                    xd[i]*100.*TMath::Cos(Angle + F10_TX1par) - 119.3 + fib_halfwidth * TMath::Cos(Angle + F10_TX1par) + F1012_Xpar,
                    yd[i]*100.,
                    0. + xd[i]*100. * TMath::Sin(Angle + F10_TX1par) + 660.2 + fib_halfwidth * TMath::Sin(Angle+F10_TX1par) + F1012_Zpar ,
                    qd[i], td[i], det[i]);
            det_hit.Set_XYZ_det(xd[i]*100., yd[i]*100., zd[i]*100);//saving original det coordiantes
            f10_hits.push_back(det_hit);
        } 
        else if(det[i] == 3) //============= fi11
        {
            det_hit.Set_XYZQ(
                    xd[i]*(-100.)*TMath::Cos(Angle + F11_TX1par) - 103.3 - fib_halfwidth * TMath::Cos(Angle + F11_TX1par) + F1113_Xpar,
                    yd[i]*100.,
                    xd[i]*(-100.) * TMath::Sin(Angle + F11_TX1par) + 595.6 - fib_halfwidth * TMath::Sin(Angle + F11_TX1par) + F1113_Zpar,
                    qd[i], td[i], det[i]);
            det_hit.Set_XYZ_det(xd[i]*100., yd[i]*100., zd[i]*100);//saving original det coordiantes
            f11_hits.push_back(det_hit);
        }
        else if(det[i] == 4) //================ fi12 -->Similar to Fib10
        {
            det_hit.Set_XYZQ(
                    xd[i]*100.*TMath::Cos(Angle + F12_TX1par) - 92.4 + fib_halfwidth * TMath::Cos(Angle + F12_TX1par) + F1012_Xpar,
                    yd[i]*100.,
                    0. + xd[i] * 100. * TMath::Sin(Angle + F12_TX1par) + 575.5 + fib_halfwidth * TMath::Sin(Angle + F12_TX1par) + F1012_Zpar ,
                    qd[i], td[i], det[i]);
            det_hit.Set_XYZ_det(xd[i]*100., yd[i]*100., zd[i]*100);//saving original det coordiantes
            f12_hits.push_back(det_hit);
        }
        else if(det[i] == 5) //=============== fi13 --> Similar to f11
        {
            det_hit.Set_XYZQ(
                    xd[i]*(-100.)*TMath::Cos(Angle + F13_TX1par) - 131.9 - fib_halfwidth * TMath::Cos(Angle+F13_TX1par) + F1113_Xpar,
                    yd[i]*100.,
                    xd[i]*(-100.) * TMath::Sin(Angle + F13_TX1par) + 682.6 - fib_halfwidth * TMath::Sin(Angle+F13_TX1par) + F1113_Zpar ,
                    qd[i], td[i], det[i]);
            det_hit.Set_XYZ_det(xd[i]*100., yd[i]*100., zd[i]*100);//saving original det coordiantes
            f13_hits.push_back(det_hit);
        }

        else if(det[i] == 6 || det[i] ==7) //dTOF 6 and 7
        {
            det_hit.Set_XYZQ(
                    xd[i]*100.*TMath::Cos(Angle) - 136.05,
                    yd[i]*100.,
                    0. + xd[i]*100. * TMath::Sin(Angle) + 703.95,
                    qd[i], td[i], det[i]);
            det_hit.Set_XYZ_det(xd[i]*100., yd[i]*100., zd[i]*100);//saving original det coordiantes
            dtof_hits.push_back(det_hit);
        }

        else if(det[i] == 8 || det[i] == 9) //dTOF 8, 9
        {
            det_hit.Set_XYZQ(
                    xd[i]*100.*TMath::Cos(Angle) - 137.55,
                    yd[i]*100.,
                    0. + xd[i]*100. * TMath::Sin(Angle) + 708.65,
                    qd[i], td[i], det[i]);
            det_hit.Set_XYZ_det(xd[i]*100., yd[i]*100., zd[i]*100);//saving original det coordiantes
            dtof_hits.push_back(det_hit);
        }
        else continue;
    }//end of transformation to the lab system

    return;
}

void R3BTrackS454_MDF::Track_MDF(int Tracking_Arm)
{
    //MDF data vectors
    PoQ_data.clear();
    MDF_Data_PoQ PoQ;
    double Xproj=0;

    //============== Tracking in the left arm ================== 
    if(Tracking_Arm==1 && f10_hits.size()>0 && f12_hits.size()>0 && dtof_hits.size()>0 && f3b_hits.size()>0)
    {
        for(auto & _hit_f12 : f12_hits)
        {
            for(auto & _hit_f10 : f10_hits)
            {
                for(auto & _hit_tofd : dtof_hits)
                {
                    //Project onto dtof from F10 (last fiber)
                    Xproj = _hit_f10.X + (_hit_tofd.Z - _hit_f10.Z) * (_hit_f10.X -_hit_f12.X)/(_hit_f10.Z - _hit_f12.Z);
                    if(_hit_tofd.Detector==6 || _hit_tofd.Detector==8) continue;
                    if(_hit_tofd.Q>10) continue;
                    if((Xproj -  _hit_tofd.X)>2. || (Xproj -  _hit_tofd.X)<(-2.5)) continue;

                    //=============== Next look for coincidences in F3b ==============
                    for(auto & _hit_f3b : f3b_hits)
                    {
                        //========== Reset PoQ object and collect all data ==========
                        PoQ.is_right_arm    = false;
                        PoQ.is_left_arm    = false;

                        //============ Saving data from each detector into the PoQ object
                        PoQ.f10_hit.Set_XYZQ( _hit_f10.X, _hit_f10.Y, _hit_f10.Z,
                                _hit_f10.Q, _hit_f10.T, _hit_f10.Detector);

                        PoQ.f12_hit.Set_XYZQ( _hit_f12.X, _hit_f12.Y, _hit_f12.Z,
                                _hit_f12.Q, _hit_f12.T, _hit_f12.Detector);

                        PoQ.dtof_hit.Set_XYZQ( _hit_tofd.X, _hit_tofd.Y, _hit_tofd.Z,
                                _hit_tofd.Q, _hit_tofd.T, _hit_tofd.Detector);

                        PoQ.f3b_hit.Set_XYZQ( _hit_f3b.X, _hit_f3b.Y, _hit_f3b.Z,
                                _hit_f3b.Q, _hit_f3b.T, _hit_f3b.Detector);

                        //===== Separately saving detector coordinates of the hits =========
                        PoQ.f10_hit.Set_XYZ_det(_hit_f10.Xdet,  _hit_f10.Ydet,  _hit_f10.Zdet );
                        PoQ.f12_hit.Set_XYZ_det(_hit_f12.Xdet,  _hit_f12.Ydet,  _hit_f12.Zdet );
                        PoQ.f3b_hit.Set_XYZ_det(_hit_f3b.Xdet,  _hit_f3b.Ydet,  _hit_f3b.Zdet );

                        //==== Tracking TX0 angle from the target ========
                        PoQ.edata_tx0[0] = F1012_X0par;//X0 nominal
                        PoQ.edata_tx0[1] = 0;//Z0
                        PoQ.edata_tx0[2] = _hit_f12.X;//X1
                        PoQ.edata_tx0[3] = _hit_f12.Z;//Z1
                        PoQ.edata_tx0[4] = (_hit_f10.X -_hit_f12.X)/(_hit_f10.Z - _hit_f12.Z);//TX1
                        MDF_TX0_targ->X2P(PoQ.edata_tx0, PoQ.pdata_tx0);
                        PoQ.value_tx0 = MDF_TX0_targ->MDF(PoQ.pdata_tx0);

                        //==== Tracking TX0 angle from f3b ========
                        PoQ.edata_tx0_f3[0] = _hit_f3b.X;//X from f3b
                        PoQ.edata_tx0_f3[1] = _hit_f3b.Z;//Z from f3b
                        PoQ.edata_tx0_f3[2] = _hit_f12.X;//X1
                        PoQ.edata_tx0_f3[3] = _hit_f12.Z;//Z1
                        PoQ.edata_tx0_f3[4] = (_hit_f10.X -_hit_f12.X)/(_hit_f10.Z - _hit_f12.Z);//TX1
                        MDF_TX0_f3->X2P(PoQ.edata_tx0_f3, PoQ.pdata_tx0_f3);
                        PoQ.value_tx0_f3 = MDF_TX0_f3->MDF(PoQ.pdata_tx0_f3);

                        //==== Tracking PoQ from the target ========
                        PoQ.edata_poq[0] = F1012_X0par;//X0 
                        PoQ.edata_poq[1] = 0.;//Z0
                        PoQ.edata_poq[2] = _hit_f12.X;
                        PoQ.edata_poq[3] = _hit_f12.Z;
                        PoQ.edata_poq[4] = (_hit_f10.X -_hit_f12.X)/(_hit_f10.Z - _hit_f12.Z);
                        MDF_PoQ->X2P(PoQ.edata_poq, PoQ.pdata_poq);
                        PoQ.value_poq = MDF_PoQ->MDF(PoQ.pdata_poq) * 1672./1672.;

                        PoQ.tx0_f3  =(PoQ.f3b_hit.X -  F1012_X0par) / PoQ.f3b_hit.Z;
                        PoQ.tx0_res = PoQ.value_tx0 - PoQ.tx0_f3;
                        PoQ.tx0_f3_res = PoQ.value_tx0_f3 - PoQ.tx0_f3;
                        PoQ.tx0_mdf_res = PoQ.value_tx0_f3 - PoQ.value_tx0;
                        PoQ.x0_proj_by_f3b =  PoQ.f3b_hit.X - (PoQ.value_tx0_f3 * PoQ.f3b_hit.Z);
                        PoQ.x0_res =  PoQ.x0_proj_by_f3b - F1012_X0par;

                        PoQ.is_left_arm = true;

                        if(fabs(PoQ.tx0_mdf_res)<0.1 && PoQ.tx0_f3 < 0.1 && PoQ.tx0_f3>0)
                        {
                            PoQ_data.push_back(PoQ);
                        }
                    }//end f3b
                    break; //using only one coincident hit with TOFD
                }//end tofd
            }//end f10
        }//end f12
    }//endif left arm

    //============== Second Tracking in the right arm ================== 
    if(Tracking_Arm==2 && f11_hits.size()>0 && f13_hits.size()>0 && dtof_hits.size()>0 && f3a_hits.size()>0)
    {
        for(auto & _hit_f11 : f11_hits)
        {
            for(auto & _hit_f13 : f13_hits)
            {
                for(auto & _hit_tofd : dtof_hits)
                {
                    //Project onto dtof from F13 (last fiber)
                    if(_hit_tofd.Detector==7 || _hit_tofd.Detector==9) continue;
                    if(_hit_tofd.Q>10) continue;
                    Xproj = _hit_f13.X + (_hit_tofd.Z - _hit_f13.Z) * (_hit_f13.X -_hit_f11.X)/(_hit_f13.Z - _hit_f11.Z);
                    if((Xproj -  _hit_tofd.X)>2. || (Xproj -  _hit_tofd.X)<(-2.5)) continue;

                    //=============== Next look for coincidences in F3a ==============
                    for(auto & _hit_f3a : f3a_hits)
                    {
                        //========== Reset PoQ object and collect all data ==========
                        PoQ.is_right_arm  = false;
                        PoQ.is_left_arm  = false;

                        //============ Saving data from each detector into the PoQ object
                        PoQ.f11_hit.Set_XYZQ( _hit_f11.X, _hit_f11.Y, _hit_f11.Z,
                                _hit_f11.Q, _hit_f11.T, _hit_f11.Detector);

                        PoQ.f13_hit.Set_XYZQ( _hit_f13.X, _hit_f13.Y, _hit_f13.Z,
                                _hit_f13.Q, _hit_f13.T, _hit_f13.Detector);

                        PoQ.dtof_hit.Set_XYZQ( _hit_tofd.X, _hit_tofd.Y, _hit_tofd.Z,
                                _hit_tofd.Q, _hit_tofd.T, _hit_tofd.Detector);

                        PoQ.f3a_hit.Set_XYZQ( _hit_f3a.X, _hit_f3a.Y, _hit_f3a.Z,
                                _hit_f3a.Q, _hit_f3a.T, _hit_f3a.Detector);

                        //===== Separately saving detector coordinates of the hits =========
                        PoQ.f11_hit.Set_XYZ_det(_hit_f11.Xdet,  _hit_f11.Ydet,  _hit_f11.Zdet );
                        PoQ.f13_hit.Set_XYZ_det(_hit_f13.Xdet,  _hit_f13.Ydet,  _hit_f13.Zdet );
                        PoQ.f3a_hit.Set_XYZ_det(_hit_f3a.Xdet,  _hit_f3a.Ydet,  _hit_f3a.Zdet );

                        //==== Tracking TX0 angle from the target ========
                        PoQ.edata_tx0[0] = F1113_X0par;//X0 nominal
                        PoQ.edata_tx0[1] = 0;//Z0
                        PoQ.edata_tx0[2] = _hit_f11.X;//X1
                        PoQ.edata_tx0[3] = _hit_f11.Z;//Z1
                        PoQ.edata_tx0[4] = (_hit_f13.X -_hit_f11.X)/(_hit_f13.Z - _hit_f11.Z);//TX1
                        MDF_TX0_targ->X2P(PoQ.edata_tx0, PoQ.pdata_tx0);
                        PoQ.value_tx0 = MDF_TX0_targ->MDF(PoQ.pdata_tx0);

                        //==== Tracking TX0 angle from f3b ========
                        PoQ.edata_tx0_f3[0] = _hit_f3a.X;//X from f3b
                        PoQ.edata_tx0_f3[1] = _hit_f3a.Z;//Z from f3b
                        PoQ.edata_tx0_f3[2] = _hit_f11.X;//X1
                        PoQ.edata_tx0_f3[3] = _hit_f11.Z;//Z1
                        PoQ.edata_tx0_f3[4] = (_hit_f13.X -_hit_f11.X)/(_hit_f13.Z - _hit_f11.Z);//TX1
                        MDF_TX0_f3->X2P(PoQ.edata_tx0_f3, PoQ.pdata_tx0_f3);
                        PoQ.value_tx0_f3 = MDF_TX0_f3->MDF(PoQ.pdata_tx0_f3);

                        //==== Tracking PoQ from the target ========
                        PoQ.edata_poq[0] = F1113_X0par;//X0 
                        PoQ.edata_poq[1] = 0.;//Z0
                        PoQ.edata_poq[2] = _hit_f11.X;
                        PoQ.edata_poq[3] = _hit_f11.Z;
                        PoQ.edata_poq[4] = (_hit_f13.X -_hit_f11.X)/(_hit_f13.Z - _hit_f11.Z);
                        MDF_PoQ->X2P(PoQ.edata_poq, PoQ.pdata_poq);
                        PoQ.value_poq = MDF_PoQ->MDF(PoQ.pdata_poq) * 1672./1672.;

                        PoQ.tx0_f3  =(PoQ.f3a_hit.X -  F1113_X0par) / PoQ.f3a_hit.Z;
                        PoQ.tx0_res = PoQ.value_tx0 - PoQ.tx0_f3;
                        PoQ.tx0_f3_res = PoQ.value_tx0_f3 - PoQ.tx0_f3;
                        PoQ.tx0_mdf_res = PoQ.value_tx0_f3 - PoQ.value_tx0;
                        PoQ.x0_proj_by_f3a =  PoQ.f3a_hit.X - (PoQ.value_tx0_f3 * PoQ.f3a_hit.Z);
                        PoQ.x0_res =  PoQ.x0_proj_by_f3a - F1113_X0par;

                        PoQ.is_right_arm = true;

                        if(fabs(PoQ.tx0_mdf_res)<0.1 && PoQ.tx0_f3>(-0.1) && PoQ.tx0_f3<0)
                        {
                            PoQ_data.push_back(PoQ);
                        }
                    }//end f3a
                    break; //using only one coincident hit with TOFD
                }//end tofd
            }//end f13
        }//end f11
    }//endif right arm
    return;
}

ClassImp(R3BTrackS454_MDF)
