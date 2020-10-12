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

        //iConverting data to lab coordinates
        void Make_Lab_Coordinates(Int_t detcount, Int_t* det, Double_t* qd, Double_t* td, Double_t * xd, Double_t * yd, Double_t * zd);
        void Track_MDF(int Tracking_Arm);

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


        R3BMDFWrapper * MDF_X0;
        R3BMDFWrapper * MDF_TX0_targ;
        R3BMDFWrapper * MDF_TX0_f3;
        R3BMDFWrapper * MDF_PoQ;

    public:

        struct Detector_Hit
        {
            Double_t X;
            Double_t Y;
            Double_t Z;
            Double_t Q;
            Double_t T;
            Int_t Detector;
            Double_t Xdet;//internal X coordiante
            Double_t Ydet;//internal Y coordiante
            Double_t Zdet;//internal Z coordinate

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

        struct MDF_Data_PoQ //X0, Z0, X1, Z1, TX1
        {
            Double_t edata_poq[5];//X0, Z0, X1, Z1, TX1
            Double_t pdata_poq[5];

            Double_t edata_tx0[5];//X0, Z0, X1, Z1, TX1
            Double_t pdata_tx0[5];
        
            Double_t edata_tx0_f3[5];//X0, Z0, X1, Z1, TX1
            Double_t pdata_tx0_f3[5];

            //MDF values
            Double_t value_poq;
            Double_t value_tx0;
            Double_t value_tx0_f3;

            Double_t tx0_f3;
            Double_t tx0_mdf_res;
            Double_t tx0_f3_res;
            Double_t tx0_res;

            Double_t x0_res;
            Double_t x0_proj_by_f3b;
            Double_t x0_proj_by_f3a;

            Detector_Hit f3b_hit;
            Detector_Hit f3a_hit;
            Detector_Hit f10_hit;
            Detector_Hit f12_hit;
            Detector_Hit f11_hit;
            Detector_Hit f13_hit;
            Detector_Hit dtof_hit;

            bool is_right_arm;
            bool is_left_arm;
        };

        std::vector<MDF_Data_PoQ> PoQ_data;
        std::vector<MDF_Data_PoQ> PoQ_data_filtered;

        bool Same_Track_Right(MDF_Data_PoQ p1, MDF_Data_PoQ p2);
        bool Used_Track_Right(MDF_Data_PoQ p1, MDF_Data_PoQ p2);
        bool Good_Track_Right(MDF_Data_PoQ p1);

        bool Same_Track_Left(MDF_Data_PoQ p1, MDF_Data_PoQ p2);
        bool Used_Track_Left(MDF_Data_PoQ p1, MDF_Data_PoQ p2);
        bool Good_Track_Left(MDF_Data_PoQ p1);



    public:
        //saving track data for tracker alignment
        TTree tree_out;

        UInt_t N_glob_tracks_left;
        UInt_t N_glob_tracks_right;

        Float_t f10_X[400];
        Float_t f10_Y[400];
        Float_t f10_Z[400];
        Float_t f10_Q[400];
        Float_t f10_T[400];
    
        Float_t f11_X[400];
        Float_t f11_Y[400];
        Float_t f11_Z[400];
        Float_t f11_Q[400];
        Float_t f11_T[400];

        Float_t f12_X[400];
        Float_t f12_Y[400];
        Float_t f12_Z[400];
        Float_t f12_Q[400];
        Float_t f12_T[400];

        Float_t f13_X[400];
        Float_t f13_Y[400];
        Float_t f13_Z[400];
        Float_t f13_Q[400];
        Float_t f13_T[400];

        Float_t f3b_X[400];
        Float_t f3b_Y[400];
        Float_t f3b_Z[400];
        Float_t f3b_Q[400];
        Float_t f3b_T[400];

        Float_t f3a_X[400];
        Float_t f3a_Y[400];
        Float_t f3a_Z[400];
        Float_t f3a_Q[400];
        Float_t f3a_T[400];

        Float_t tofd_left_X[400];
        Float_t tofd_left_Y[400];
        Float_t tofd_left_Z[400];
        Float_t tofd_left_Q[400];
        Float_t tofd_left_T[400];
        
        Float_t tofd_right_X[400];
        Float_t tofd_right_Y[400];
        Float_t tofd_right_Z[400];
        Float_t tofd_right_Q[400];
        Float_t tofd_right_T[400];

        Float_t TX0_mdf_left[400];
        Float_t TX0_f3_mdf_left[400];
        Float_t PoQ_mdf_left[400];
        Float_t X0_residual_left[400];
        Float_t TX0_residual_left[400];
        Float_t TX0_residual_f3_left[400];
        Float_t TX0_residual_mdf_left[400];
        Float_t X0_proj_by_f3_left[400];
        Float_t TX0_proj_by_f3_left[400];
        Float_t Xoffset_f3b[400];
 
        Float_t TX0_mdf_right[400];
        Float_t TX0_f3_mdf_right[400];
        Float_t PoQ_mdf_right[400];
        Float_t X0_residual_right[400];
        Float_t TX0_residual_right[400];
        Float_t TX0_residual_mdf_right[400];
        Float_t TX0_residual_f3_right[400];
        Float_t X0_proj_by_f3_right[400];
        Float_t TX0_proj_by_f3_right[400];
        Float_t Xoffset_f3a[400];

        UInt_t Nhits_f3a;
        UInt_t Nhits_f3b;
        UInt_t Nhits_f10;
        UInt_t Nhits_f11;
        UInt_t Nhits_f12;
        UInt_t Nhits_f13;
        UInt_t Nhits_tofd;

        Int_t Current;

        //Alignment constants
        const double Angle = 17.7 * TMath::Pi()/180.;//Central turning angle from Daniel
        const double Z0 = 277.7; //cm from target middle to the central turning point in GLAD
        const double f3ab_halfwidth =  0.021*256;//210 um pitch
        const double fib_halfwidth =  512.*0.05;//500um pitch

        //Alignment offsets for F10,12
        const double F10_TX1par  =  0.04862;
        const double F12_TX1par  = -0.01913;
        const double F1012_Xpar  = -0.5151;
        const double F1012_Zpar  =  0.8256;
        const double F1012_X0par =  0.2449;

        //Alignment offsets for F11,13
        const double F11_TX1par  = -0.02797;
        const double F13_TX1par  = -0.04109;
        const double F1113_Xpar  = -0.4749;
        const double F1113_Zpar  =  1.321;
        const double F1113_X0par =  0.23;

        //Alignment offset of f3b
        const double F3b_dX  =  -0.253;

        const double Res_Limit = 0.003;

    public:
        ClassDef(R3BTrackS454_MDF, 1)
};

#endif
