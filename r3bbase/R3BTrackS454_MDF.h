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
// -----            Created 13-04-2016 by M.Heil          -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#ifndef R3BTRACKS454_MDF
#define R3BTRACKS454_MDF
#define N_PLANE_MAX_TOFD 4
#define N_PADDLE_MAX_TOFD 50
#define N_PADDLE_MAX_PTOF 100
#define N_PSPX_S454 4

#include "FairTask.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include "TCutG.h"

#include "TClonesArray.h"
#include "TMath.h"
#include <cstdlib>

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;
class R3BMDFWrapper;

/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BTrackS454_MDF : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BTrackS454_MDF();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BTrackS454_MDF(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTrackS454_MDF();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    virtual void Output1(Double_t tracker[6], Double_t chi2[2]);
    virtual void Output2(Double_t tracker[6], Double_t chi2[2]);
    
    //Main tracking function
    void TrackMDF(Int_t detcount, Int_t* det, Double_t* qd, Double_t* td, Double_t * xd, Double_t * yd, Double_t * zd);

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

    /**
     * Methods for setting reset and readout times for Bmon
     */
    inline void SetBmon(Int_t sens_SEE, Int_t sens_IC)
    {
        fsens_SEE = sens_SEE;
        fsens_IC = sens_IC;
    }

    /**
     * Methods for setting cuts
     */
    inline void SetCuts(Bool_t cuts)
    {
        fCuts = cuts;
    }

    inline void SetGraphicalCuts(Bool_t graphCuts)
    {
        fGraphCuts = graphCuts;
    }

    inline void SetGhost(Bool_t ghost)
    {
        fGhost = ghost;
    }

    inline void SetPairs(Bool_t p)
    {
        fPairs = p;
    }

    inline void SetBfield(Int_t B)
    {
        fB = B;
    }
    inline void SetSimu(Int_t simu)
    {
        fSimu = simu;
    }

  private:
    std::vector<TClonesArray*> fMappedItems;
    std::vector<TClonesArray*> fCalItems;
    std::vector<TClonesArray*> fHitItems;
    TClonesArray* fMCTrack;
    TClonesArray* fTrackItems;
    Int_t fNofTrackItems;


    enum DetectorInstances
    {
        DET_CALIFA,
        DET_BMON,
        DET_ROLU,
        DET_FI_FIRST,
        DET_FI3A = DET_FI_FIRST,
        DET_FI3B,
        DET_FI10,
        DET_FI11,
        DET_FI12,
        DET_FI13,
        DET_FI_LAST = DET_FI13,
        DET_TOFD,
        DET_MAX
    };

#define NOF_FIB_DET (DET_FI_LAST - DET_FI_FIRST + 1)

    const char* fDetectorNames[DET_MAX + 1] = { "Califa", "BeamMonitor", "Rolu", "Fi3a", "Fi3b", "Fi10",
                                                "Fi11",   "Fi12",        "Fi13", "Tofd", NULL };

    // If FiberI is present or not:
    Int_t ifibdet;
    // Number of fibers per detector
    Double_t n_fiber[NOF_FIB_DET] = { 512., 512., 2048., 2048., 2048., 2048. };

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat;
	Bool_t fCuts;
	Bool_t fGhost;
	Bool_t fPairs;
	Bool_t fGraphCuts;
	Bool_t fSimu;
	Int_t fB;
	Bool_t tracker = true;

	TCutG *cut_Fi10vsTofd;
	TCutG *cut_Fi13vsTofd;
	
    unsigned long long time_start = 0, time = 0;
    unsigned long ic_start = 0, see_start = 0, tofdor_start = 0;
    unsigned long fNEvents = 0, fNEvents_start = 0; /**< Event counter. */

    Int_t maxevent;

    Int_t fsens_SEE, fsens_IC; // SEETRAM and IC sensitivity, between -4 and -10
    Double_t calib_SEE = 1.;   // SEETRAM calibration factor
    Double_t see_offset = 7.1; // SEETRAM offset in kHz
    Double_t counts_SEE = 0;
    Double_t counts_IC = 0;
    Double_t counts_TofD = 0;

	Double_t XHes, YHes, ZHes, XCs, YCs, ZCs, THes, TCs;
	Double_t pHexs, pHeys, pHezs, pCxs, pCys, pCzs, pHes, pCs;
	
	Double_t amu = 931.49410242;
	Double_t pHex, pHey, pHez, pCx, pCy, pCz;
	Double_t Pxf, Pyf, Pzf, Xf, Yf, Zf, Pf_tot;
//	Double_t mHe = 4.00260325413*amu;
//	Double_t mC = 12. * amu;
	Double_t mHe = 3727.409;
	Double_t mC = 11174.950;
	
	Int_t Q = 0;
	Double_t tPrev[10];
	Int_t detPrev[10];
	
	Int_t counter1 = 0;
	Int_t counter2 = 0;
	Int_t counter3 = 0;
	Int_t counter4 = 0;
	Int_t counterTofd = 0;
	Int_t counterTofdMulti = 0;
	Int_t counterCalifa = 0;
	Int_t counterWrongTpat = 0;
	Int_t counterWrongTrigger = 0;
	Int_t counterRolu = 0;
	Int_t counterTracker = 0;
	Int_t countdet;

    Bool_t writeFile = false;

    UInt_t num_spills = 0;

	Int_t ndet = 10;


    //Define what is needed for MDF tracking

    TH1F * h_X0_mdf;
    TH1F * h_TX0_mdf;
    TH1F * h_PoQ_mdf;
    TH1F * h_glob_track_mul;
    TH2F * h_glob_track_mul_corr;
    TH2F * h_lab_xz;

    TH1F * h_fib10_tof;
    TH1F * h_fib12_tof;
    TH1F * h_fib12_10_tof;
    
    TH1F * h_fib10_q;
    TH1F * h_fib12_q;

    TH1F * h_fib10_x;
    TH1F * h_fib12_x;

    TH1F * h_tofd9_q;
    TH1F * h_tofd7_q;

    R3BMDFWrapper * MDF_X0;
    R3BMDFWrapper * MDF_TX0;
    R3BMDFWrapper * MDF_PoQ;

  public:

    struct Detector_Hit
    {
        Double_t X;
        Double_t Xdet;//internal X coordiante
        Double_t Y;
        Double_t Ydet;//internal Y coordiante
        Double_t Z;
        Double_t Zdet;//internal Z coordinate
        Double_t Q;
        Double_t T;
        Int_t Detector;
        void Set_XYZQ(Double_t _x, Double_t _y, Double_t _z, Double_t _q, Double_t _t, Int_t _det);
        void Set_XYZ_det(Double_t _x, Double_t _y, Double_t _z);//saving internal det coordinates (for alignemnt)
    };


    std::vector<Detector_Hit> f3a_hits;
    std::vector<Detector_Hit> f3b_hits;
    std::vector<Detector_Hit> f10_hits;
    std::vector<Detector_Hit> f11_hits;
    std::vector<Detector_Hit> f12_hits;
    std::vector<Detector_Hit> f13_hits;
    std::vector<Detector_Hit> dtof_hits;

    //Internal det coordinate storage
    std::vector<Detector_Hit> f10_hits_det;
    std::vector<Detector_Hit> f12_hits_det;

    //MDF data containers
    struct MDF_Data_X0 //TX0, X1, Z1, TX1
    {
        Double_t edata[4];
        Double_t pdata[4];
        Double_t value;
    };

    struct MDF_Data_TX0 //X0, Z0, X1, Z1, TX1
    {
        Double_t edata[5];
        Double_t pdata[5];
        Double_t value;
    };

    struct MDF_Data_PoQ //X0, Z0, X1, Z1, TX1
    {
        Double_t edata[5];
        Double_t pdata[5];
        Double_t value;
        Detector_Hit f10_hit;
        Detector_Hit f12_hit;
    };


    std::vector<MDF_Data_X0> X0_data;
    std::vector<MDF_Data_TX0> TX0_data;
    std::vector<MDF_Data_PoQ> PoQ_data;

  public:
    //saving track data for tracker alignment
    TTree tree_out;
    Double_t f10_X;
    Double_t f12_X;
    Int_t Current;
    bool is_init_out;

    

  public:
    ClassDef(R3BTrackS454_MDF, 1)
};

#endif
