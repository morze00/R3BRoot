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
        MDF_X0 = new R3BMDFWrapper();
        MDF_X0->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/X0/MDF_X0_20200622.txt");
        MDF_X0->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/X0/PCA_X0_20200622.txt");

        //Reconstruction of the target X coordinate
        MDF_TX0 = new R3BMDFWrapper();
        MDF_TX0->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/TX0/MDF_TX0_20200622.txt");
        MDF_TX0->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/TX0/PCA_TX0_20200622.txt");

        MDF_PoQ = new R3BMDFWrapper();
        //Reconstruction of the PoQ (produciton runs)
        
        //Reconstruction of the PoQ (production runa, I=1672A)
        if(fB == -1672)
        {
            MDF_PoQ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/MDF_PoQ_1672A_20200622.txt");
            MDF_PoQ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/PCA_PoQ_1672A_20200622.txt");
        }

        //Reconstruction of the PoQ (run 395, I=1102A)
        if(fB == -1102)
        {
            MDF_PoQ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/MDF_PoQ_1102A_20200713.txt");
            MDF_PoQ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/PCA_PoQ_1102A_20200713.txt");
        }

        //Reconstruction of the PoQ (run 391, I=1292A)
        if(fB == -1292)
        {
            MDF_PoQ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/MDF_PoQ_1292A_20200713.txt");
            MDF_PoQ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/PCA_PoQ_1292A_20200713.txt");
        }

        //Reconstruction of the PoQ (run 389, I=1482A)
        if(fB == -1482)
        {
            MDF_PoQ->InitMDF("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/MDF_PoQ_1482A_20200713.txt");
            MDF_PoQ->InitPCA("/Users/vpanin/r3broot/my_codes/s454_ana/MDF_params/PoQ/PCA_PoQ_1482A_20200713.txt");
        }

    }

    h_X0_mdf = new TH1F("h_X0_mdf","h_X0_mdf",1000,-20,20);
    h_TX0_mdf = new TH1F("h_TX0_mdf","h_TX0_mdf",4000,-0.01,0.01);
    h_PoQ_mdf = new TH1F("h_PoQ_mdf","h_PoQ_mdf",10000,2,2.5);

    h_lab_xz = new TH2F("h_lab_xz","h_lab_xz",10000,0,800,10000,-500,500);
    h_glob_track_mul = new TH1F("h_glob_track_mul","h_glob_track_mul",100,0,100);
    h_glob_track_mul_corr = new TH2F("h_glob_track_mul_corr","h_glob_track_mul_corr",100,0,100,100,0,100);

    h_fib10_tof = new TH1F("h_fib10_tof","h_fib10_tof",1000,-20,20);
    h_fib12_tof = new TH1F("h_fib12_tof","h_fib12_tof",1000,-20,20);
    h_fib12_10_tof = new TH1F("h_fib12_10_tof","h_fib12_10_tof",1000,-40,40);

    h_fib10_q = new TH1F("h_fib10_q","h_fib10_q",1000,0,20);
    h_fib12_q = new TH1F("h_fib12_q","h_fib12_q",1000,0,20);

    h_fib10_x = new TH1F("h_fib10_x","h_fib10_x",2000,-120,0);
    h_fib12_x = new TH1F("h_fib12_x","h_fib12_x",2000,-120,0);

    h_tofd7_q = new TH1F("h_tofd7_q","h_tofd7_q",1000,0,20);
    h_tofd9_q = new TH1F("h_tofd9_q","h_tofd9_q",1000,0,20);

    is_init_out =false;
    tree_out.SetName("tree_out");
    tree_out.Branch("f10_X",     &f10_X,  "f10_X/D");
    tree_out.Branch("f12_X",     &f12_X,  "f12_X/D");
    tree_out.Branch("Current",     &Current,  "Current/I");

    return kSUCCESS;
}

void R3BTrackS454_MDF::Exec(Option_t* option)
{
    if(!is_init_out)
    {
        //tree_out = new TTree("tree_out","tree_out"); 
        is_init_out = true;
    }

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
        if (fCuts && (ttt < -100. || ttt > 100.) && !fSimu) // trigger window -1500, 1500
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

            if (!fPairs && (qqq < 7.5 || qqq > 8.5))
                continue;

            if (fPairs && !(qqq > 1.5 && qqq < 2.5) && !(qqq > 5.5 && qqq < 6.5))
                continue;

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

            if(det2==7) tStart = t2[det2];

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
        y1[det] = 0.; // hit3a->GetY() / 100.;
        z1[det] = 0.;
        q1[det] = hit3a->GetEloss();
        t1[det] = hit3a->GetTime();

        if (q1[det] > 5.)
            q1[det] = 8.;
        else if (q1[det] > 2.)
            q1[det] = 6.;
        else if (q1[det] > 0.)
            q1[det] = 2.;

        tof = tStart - t1[det];

        // Cuts on Fi3a
        //            if (fCuts && y1[det] * 100. > 50.)
        //                continue;
        //            if (fCuts && y1[det] *100. < -50.)
        //                continue;
        //if (fCuts && !fPairs && q1[det] < cutQ)
        //    continue;
        //if (fCuts && (tof < -40 || tof > 20))
        //    continue;

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
            tdet[countdet] = t1[det];
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
        y1[det] = 0.; // hit3b->GetY() / 100.;
        z1[det] = 0.;
        q1[det] = hit3b->GetEloss();
        t1[det] = hit3b->GetTime();

        if (q1[det] > 5.)
            q1[det] = 8.;
        else if (q1[det] > 2.)
            q1[det] = 6.;
        else if (q1[det] > 0.)
            q1[det] = 2.;

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
            tdet[countdet] = t1[det];
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

    //===============================================================
    // Tracker
    chi2 = 1.E100;
    counter2++;

    //Dummy values fo for the output tree
    f10_X = -9999.;
    f12_X = -9999.;
    Current =0;

    if (tracker)
    {
        Double_t test[9];
        if(mult10>0 //&& mult10<20
                && mult12>0// && mult12<20
                //&& mult11 ==0 && mult13==0
                //&& mult3a ==0 && mult3b==0
          )
        {
            TrackMDF(countdet, detector, qdet, tdet, xdet, ydet, zdet);
            if(PoQ_data.size()==1)
            {
                for(auto & _poq_data : PoQ_data)
                {
                    h_PoQ_mdf->Fill( _poq_data.value);
                    h_lab_xz->Fill(_poq_data.f10_hit.Z, _poq_data.f10_hit.X);
                    h_lab_xz->Fill(_poq_data.f12_hit.Z, _poq_data.f12_hit.X);

                    h_fib10_tof->Fill(_poq_data.f10_hit.T);
                    h_fib12_tof->Fill(_poq_data.f12_hit.T);
                    h_fib12_10_tof->Fill(_poq_data.f12_hit.T - _poq_data.f10_hit.T);

                    h_fib10_q->Fill(_poq_data.f10_hit.Q);
                    h_fib12_q->Fill(_poq_data.f12_hit.Q);
                    
                    h_fib10_x->Fill(_poq_data.f10_hit.X);
                    h_fib12_x->Fill(_poq_data.f12_hit.X);

                    //Fill output tree for alignment
                    f10_X = _poq_data.f10_hit.Xdet;
                    f12_X = _poq_data.f12_hit.Xdet;
                    Current = fB;
                    tree_out.Fill();

                }
                for(auto & _x0_data : X0_data)
                    h_X0_mdf->Fill( _x0_data.value);

                for(auto & _tx0_data : TX0_data)
                    h_TX0_mdf->Fill( _tx0_data.value);
            }
        }

        if(f10_hits.size()>0 
                && f12_hits.size()>0 
                && dtof_hits.size()>0
                && f3a_hits.size()==0
                && f3b_hits.size()==0
          )
        {

            h_glob_track_mul->Fill(X0_data.size());
            h_glob_track_mul_corr->Fill(X0_data.size(), (f12_hits.size() * f10_hits.size()));
        }



        Bool_t det_coord = true;

        Bool_t st = false;

        //cout << "\n\nNew event ******************************" << endl;
        //cout << "\ncountdet = " << countdet;
        //cout << "\n\nNhits in fi3a " << f3a_hits.size();
        //cout << "\nNhits in fi3b " << f3b_hits.size();
        //cout << "\nNhits in fi10 " << f10_hits.size();
        //cout << "\nNhits in fi11 " << f11_hits.size();
        //cout << "\nNhits in fi12 " << f12_hits.size();
        //cout << "\nNhits in fi13 " << f13_hits.size();
        //cout << "\nNhits in dTof " << dtof_hits.size();
        //    
        //cout << "\n\ndTof hits:";
        //for(auto & _hit : dtof_hits)
        //{
        //    cout << "\ndTof.X = " << _hit.X <<  "\tdTof.Y = " << _hit.Y
        //        << "\tdTof.Z = " << _hit.Z << "\tdTof.Q = " << _hit.Q 
        //        << "\tdTof.Det = " << _hit.Detector; 
        //}

        //cout << "\n\nf10 hits:";
        //for(auto & _hit : f10_hits)
        //{
        //    cout << "\nf10.X = " << _hit.X <<  "\tf10.Y = " << _hit.Y
        //        << "\tf10.Z = " << _hit.Z << "\tf10.Q = " << _hit.Q 
        //        << "\tf10.Det = " << _hit.Detector; 
        //}

        //cout << "\n\nf12 hits:";
        //for(auto & _hit : f12_hits)
        //{
        //    cout << "\nf12.X = " << _hit.X <<  "\tf12.Y = " << _hit.Y
        //        << "\tf12.Z = " << _hit.Z << "\tf12.Q = " << _hit.Q 
        //        << "\tf12.Det = " << _hit.Detector; 
        //}


        chi2 = chi[0] + chi[1];

        if (chi[0] < 1.e10)
            counter3++;
        if (chi[1] < 1.e10)
            counter4++;
        if (chi[0] < 1.e10 && chi[1] < 1.e10)
        {
            // fill histograms
            Output2(track, chi);
        }

        if (chi[0] < 1.e10 && chi[1] < 1.e10)
        {
            counterTracker++;
            LOG(DEBUG) << "track1: " << track[0] << "  " << track[1] << "  " << track[2] << endl;
            LOG(DEBUG) << "track1: " << track[3] << "  " << track[4] << "  " << track[5] << endl;
            LOG(DEBUG) << "chi: " << chi[0] << "  " << chi[1] << "  " << chi[2] << "  " << chi[3] << "  "
                << chi[4] << "  " << chi[5] << endl;

            // we have a hit
            for (Int_t i = 0; i < ndet; i++)
            {
                xTrack[i] = -1000.;
                yTrack[i] = -1000.;
                zTrack[i] = -1000.;
                qTrack[i] = -1000.;
            }
            Int_t charge = 0;
            LOG(DEBUG2) << "# of points back" << countdet << endl;
            for (Int_t i = 0; i < countdet; i++)
            {

                LOG(DEBUG2) << "back #" << i << " Det: " << detector[i] << " x: " << xdet[i]
                    << " y: " << ydet[i] << " q: " << qdet[i] << endl;
                xTrack[detector[i]] = xdet[i];
                yTrack[detector[i]] = ydet[i];
                zTrack[detector[i]] = zdet[i];
                qTrack[detector[i]] = qdet[i];
                if (qdet[i] > charge)
                    charge = qdet[i];
            }

            // store hits in track level
            new ((*fTrackItems)[fNofTrackItems++]) R3BTrack(
                    track[0], track[1], track[2], track[3], track[4], track[5], charge, 2, chi[0], chi[1], 0);
        }
    }

    for (Int_t i = 0; i < 12; i++)
    {
        track[i] = 0.;
    }
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

    h_X0_mdf->Write();
    h_TX0_mdf->Write();
    h_PoQ_mdf->Write();
    h_lab_xz->Write();
    h_glob_track_mul->Write();
    h_glob_track_mul_corr->Write();

    h_fib10_tof->Write();
    h_fib12_tof->Write();
    h_fib12_10_tof->Write();
    h_fib12_q->Write();
    h_fib10_q->Write();

    h_fib12_x->Write();
    h_fib10_x->Write();

    h_tofd7_q->Write();
    h_tofd9_q->Write();

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

void R3BTrackS454_MDF::TrackMDF(Int_t detcount, Int_t* det, Double_t* qd, Double_t* td, Double_t * xd, Double_t * yd, Double_t * zd)
{
    //Make sure that hit containers are empty
    f3a_hits.clear();
    f3b_hits.clear();
    f10_hits.clear();
    f11_hits.clear();
    f12_hits.clear();
    f13_hits.clear();
    dtof_hits.clear();

    f10_hits_det.clear();
    f12_hits_det.clear();


    //MDF data vectors
    X0_data.clear();
    TX0_data.clear();
    PoQ_data.clear();

    //Temporary containers
    Detector_Hit det_hit;
    MDF_Data_X0 X0;
    MDF_Data_TX0 TX0;
    MDF_Data_PoQ PoQ;

    double Angle = 17.7 * TMath::Pi()/180.;//Central turning angle from Daniel
    double Z0 = 277.7; //cm from target middle to the central turning point in GLAD
    double f3ab_halfwidth =  0.021*256;//210 um pitch
    double fib_halfwidth =  512.*0.05;//500um pitch


    bool good_tofd7=false;
    bool good_tofd9=false;

    //At first convert each detector hit to the lab coordinates (in cm as r3broot)
    for(int i=0; i<detcount; i++)
    {
        if(det[i] == 0) //fi3a
        {
            det_hit.Set_XYZQ(
                    xd[i]*(-100.) - 0.25 -  f3ab_halfwidth,// 5 mm slit, flipped
                    yd[i],
                    zd[i]+80.5,
                    qd[i], td[i], det[i]);

            f3a_hits.push_back(det_hit);
        }

        else if(det[i] == 1) //fi3b
        {
            det_hit.Set_XYZQ(
                    xd[i]*100. + 0.25 + f3ab_halfwidth, //5mm slit
                    yd[i],
                    zd[i]+80.5,
                    qd[i], td[i], det[i]);

            f3b_hits.push_back(det_hit);
        }

        else if(det[i] == 2) //fi10
        {
            det_hit.Set_XYZQ(
                    xd[i]*100.*TMath::Cos(Angle) - 119.3 + fib_halfwidth * TMath::Cos(Angle),
                    yd[i],
                    0. + xd[i]*100. * TMath::Sin(Angle) + 660.2 + fib_halfwidth * TMath::Sin(Angle) ,
                    qd[i], td[i], det[i]);
                
            det_hit.Set_XYZ_det(xd[i]*100., yd[i]*100., zd[i]*100);//saving original det coordiantes

            f10_hits.push_back(det_hit);
            
        } 

        else if(det[i] == 3) //fi11
        {
            det_hit.Set_XYZQ(xd[i]*100,
                    yd[i],
                    zd[i],
                    qd[i], td[i], det[i]);
            f11_hits.push_back(det_hit);
        }

        else if(det[i] == 4) //fi12 -->Similar to Fib10
        {
            det_hit.Set_XYZQ(
                    xd[i]*100.*TMath::Cos(Angle) - 92.4 + fib_halfwidth * TMath::Cos(Angle),
                    yd[i],
                    0. + xd[i] * 100. * TMath::Sin(Angle) + 575.5 + fib_halfwidth * TMath::Sin(Angle) ,
                    qd[i], td[i], det[i]);

            det_hit.Set_XYZ_det(xd[i]*100., yd[i]*100., zd[i]*100);//saving original det coordiantes

            f12_hits.push_back(det_hit);
        }

        else if(det[i] == 5) //fi13
        {
            det_hit.Set_XYZQ(
                    xd[i]*100.,
                    yd[i],
                    zd[i],
                    qd[i], td[i], det[i]);
            f13_hits.push_back(det_hit);
        }

        else if(det[i] == 6 || det[i] ==7 || det[i] ==8 || det[i]==9) //dTOF
        {
            det_hit.Set_XYZQ(
                    xd[i]*100.,
                    yd[i]*100.,
                    zd[i],
                    qd[i], td[i], det[i]);
            dtof_hits.push_back(det_hit);

                   }

        else continue;
    }//end of transformation to the lab system

    //Now actual tracking
    if(f12_hits.size()==0 || f10_hits.size()==0 
            //|| f12_hits.size()>4 || f10_hits.size()>4
            || dtof_hits.size()==0
            || f3a_hits.size()!=0 || f3b_hits.size()!=0) return;


    bool is_used_combo;

    //cout << "\n\nNew event ******************************" << endl;
    for(auto & _hit_f12 : f12_hits)
    {

        if(_hit_f12.Q<6) continue;

        for(auto & _hit_f10 : f10_hits)
        {

            if(_hit_f10.Q<6) continue;

            //if((_hit_f12.T-_hit_f10.T)<(-1.2) ||
            //       (_hit_f12.T-_hit_f10.T)>0.4 ) continue;


            //==== Tracking X0 coordiante (target X position) ========
            X0.edata[0] = 0.; // TX0 value for unreacted beam
            X0.edata[1] = _hit_f12.X;
            X0.edata[2] = _hit_f12.Z;
            X0.edata[3] = (_hit_f10.X -_hit_f12.X)/(_hit_f10.Z - _hit_f12.Z);
            MDF_X0->X2P(X0.edata, X0.pdata);
            X0.value = MDF_X0->MDF(X0.pdata);

            //==== Tracking TX0 angle from the target ========
            TX0.edata[0] = 0;//X0
            TX0.edata[1] = 0.;//Z0
            TX0.edata[2] = _hit_f12.X;
            TX0.edata[3] = _hit_f12.Z;
            TX0.edata[4] = (_hit_f10.X -_hit_f12.X)/(_hit_f10.Z - _hit_f12.Z);
            MDF_TX0->X2P(TX0.edata, TX0.pdata);
            TX0.value = MDF_TX0->MDF(TX0.pdata);

            //==== Repeat X0 after tracking ========
            //X0.edata[0] = TX0.value; // TX0 value for unreacted beam
            //MDF_X0->X2P(X0.edata, X0.pdata);
            //X0.value = MDF_X0->MDF(X0.pdata);

            //==== Repeat TX0 angle from the target ========
            //TX0.edata[0] = X0.value;//X0
            //MDF_TX0->X2P(TX0.edata, TX0.pdata);
            //TX0.value = MDF_TX0->MDF(TX0.pdata);

            //// //==== Repeat X0 after tracking ========
            //X0.edata[0] = TX0.value; // TX0 value for unreacted beam
            //MDF_X0->X2P(X0.edata, X0.pdata);
            //X0.value = MDF_X0->MDF(X0.pdata);

            //==== Tracking PoQ from the target ========
            PoQ.edata[0] = X0.value;//X0
            PoQ.edata[1] = 0.;//Z0
            PoQ.edata[2] = _hit_f12.X;
            PoQ.edata[3] = _hit_f12.Z;
            PoQ.edata[4] = (_hit_f10.X -_hit_f12.X)/(_hit_f10.Z - _hit_f12.Z);
            MDF_PoQ->X2P(PoQ.edata, PoQ.pdata);
            //PoQ.value = MDF_PoQ->MDF(PoQ.pdata) * 1482./ 1672.;
            PoQ.value = MDF_PoQ->MDF(PoQ.pdata);

            if(     
                    //run 395
                    //X0.value < -1.5 || 
                    //X0.value > 0.2 || 
                    //TX0.value < -0.005 || 
                    //TX0.value > 0.

                    //run 391
                    //X0.value < -1. || 
                    //X0.value > 1. || 
                    //TX0.value < -0.003 || 
                    //TX0.value > 0.003

                    //run 389
                    X0.value < 0. || 
                    X0.value > 1.3 || 
                    TX0.value < 0 || 
                    TX0.value > 0.005 
              ) continue; 

            //cout << "\n\nF12 hit X = " <<  _hit_f12.X;
            //cout << "\nF10 hit X = " << _hit_f10.X;

            PoQ.f10_hit.Set_XYZQ(
                    _hit_f10.X,
                    _hit_f10.Y,
                    _hit_f10.Z,
                    _hit_f10.Q,
                    _hit_f10.T,
                    _hit_f10.Detector);

            PoQ.f12_hit.Set_XYZQ(
                    _hit_f12.X,
                    _hit_f12.Y,
                    _hit_f12.Z,
                    _hit_f12.Q,
                    _hit_f12.T,
                    _hit_f12.Detector);

            //save also det coordinates
            PoQ.f10_hit.Set_XYZ_det(
                    _hit_f10.Xdet,
                    _hit_f10.Ydet,
                    _hit_f10.Zdet
                    );

          //save also det coordinates
            PoQ.f12_hit.Set_XYZ_det(
                    _hit_f12.Xdet,
                    _hit_f12.Ydet,
                    _hit_f12.Zdet
                    );


            for(auto & _hit_tofd : dtof_hits)
            {
                if(_hit_tofd.Detector==7
                        //&& _hit_tofd.Q>7
                  )
                {
                    h_tofd7_q->Fill(_hit_tofd.Q);
                    good_tofd7=true;
                }

                if(_hit_tofd.Detector==9 
                        //&& _hit_tofd.Q>7
                  )
                {
                    h_tofd9_q->Fill(_hit_tofd.Q);
                    good_tofd9=true;
                }

            }

            if(good_tofd7 && good_tofd9)
            {
                TX0_data.push_back(TX0);
                X0_data.push_back(X0);
                PoQ_data.push_back(PoQ);
            }
            //cout << "\n\tReconstructed X0 = " << X0.value;
        }
    }
    return;
}

ClassImp(R3BTrackS454_MDF)
