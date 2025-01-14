/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2025 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

/********************************************************
 *
 * Structure for ext_data_fetch_event() filling.
 *
 * Do not edit - automatically generated.
 */

#ifndef __GUARD_H101_FIBFOURTY_EXT_H101_FIBFOURTY_H__
#define __GUARD_H101_FIBFOURTY_EXT_H101_FIBFOURTY_H__

#ifndef __CINT__
# include <stdint.h>
#else
/* For CINT (old version trouble with stdint.h): */
# ifndef uint32_t
typedef unsigned int uint32_t;
typedef          int  int32_t;
# endif
#endif
#ifndef EXT_STRUCT_CTRL
# define EXT_STRUCT_CTRL(x)
#endif

/********************************************************
 *
 * Plain structure (layout as ntuple/root file):
 */

typedef struct EXT_STR_h101_FIBFOURTY_t
{
  /* RAW */
  uint32_t FIBFOURTY_TBLCM /* [1,1032] */;
  uint32_t FIBFOURTY_TBLCMI[1032 EXT_STRUCT_CTRL(FIBFOURTY_TBLCM)] /* [1,1032] */;
  uint32_t FIBFOURTY_TBLCME[1032 EXT_STRUCT_CTRL(FIBFOURTY_TBLCM)] /* [1,66048] */;
  uint32_t FIBFOURTY_TBLC /* [0,66048] */;
  uint32_t FIBFOURTY_TBLCv[66048 EXT_STRUCT_CTRL(FIBFOURTY_TBLC)] /* [0,65535] */;
  uint32_t FIBFOURTY_TBLFM /* [1,1032] */;
  uint32_t FIBFOURTY_TBLFMI[1032 EXT_STRUCT_CTRL(FIBFOURTY_TBLFM)] /* [1,1032] */;
  uint32_t FIBFOURTY_TBLFME[1032 EXT_STRUCT_CTRL(FIBFOURTY_TBLFM)] /* [1,66048] */;
  uint32_t FIBFOURTY_TBLF /* [0,66048] */;
  uint32_t FIBFOURTY_TBLFv[66048 EXT_STRUCT_CTRL(FIBFOURTY_TBLF)] /* [0,65535] */;
  uint32_t FIBFOURTY_TBTCM /* [1,1032] */;
  uint32_t FIBFOURTY_TBTCMI[1032 EXT_STRUCT_CTRL(FIBFOURTY_TBTCM)] /* [1,1032] */;
  uint32_t FIBFOURTY_TBTCME[1032 EXT_STRUCT_CTRL(FIBFOURTY_TBTCM)] /* [1,66048] */;
  uint32_t FIBFOURTY_TBTC /* [0,66048] */;
  uint32_t FIBFOURTY_TBTCv[66048 EXT_STRUCT_CTRL(FIBFOURTY_TBTC)] /* [0,65535] */;
  uint32_t FIBFOURTY_TBTFM /* [1,1032] */;
  uint32_t FIBFOURTY_TBTFMI[1032 EXT_STRUCT_CTRL(FIBFOURTY_TBTFM)] /* [1,1032] */;
  uint32_t FIBFOURTY_TBTFME[1032 EXT_STRUCT_CTRL(FIBFOURTY_TBTFM)] /* [1,66048] */;
  uint32_t FIBFOURTY_TBTF /* [0,66048] */;
  uint32_t FIBFOURTY_TBTFv[66048 EXT_STRUCT_CTRL(FIBFOURTY_TBTF)] /* [0,65535] */;
  uint32_t FIBFOURTY_TTLCM /* [1,1032] */;
  uint32_t FIBFOURTY_TTLCMI[1032 EXT_STRUCT_CTRL(FIBFOURTY_TTLCM)] /* [1,1032] */;
  uint32_t FIBFOURTY_TTLCME[1032 EXT_STRUCT_CTRL(FIBFOURTY_TTLCM)] /* [1,66048] */;
  uint32_t FIBFOURTY_TTLC /* [0,66048] */;
  uint32_t FIBFOURTY_TTLCv[66048 EXT_STRUCT_CTRL(FIBFOURTY_TTLC)] /* [0,65535] */;
  uint32_t FIBFOURTY_TTLFM /* [1,1032] */;
  uint32_t FIBFOURTY_TTLFMI[1032 EXT_STRUCT_CTRL(FIBFOURTY_TTLFM)] /* [1,1032] */;
  uint32_t FIBFOURTY_TTLFME[1032 EXT_STRUCT_CTRL(FIBFOURTY_TTLFM)] /* [1,66048] */;
  uint32_t FIBFOURTY_TTLF /* [0,66048] */;
  uint32_t FIBFOURTY_TTLFv[66048 EXT_STRUCT_CTRL(FIBFOURTY_TTLF)] /* [0,65535] */;
  uint32_t FIBFOURTY_TTTCM /* [1,1032] */;
  uint32_t FIBFOURTY_TTTCMI[1032 EXT_STRUCT_CTRL(FIBFOURTY_TTTCM)] /* [1,1032] */;
  uint32_t FIBFOURTY_TTTCME[1032 EXT_STRUCT_CTRL(FIBFOURTY_TTTCM)] /* [1,66048] */;
  uint32_t FIBFOURTY_TTTC /* [0,66048] */;
  uint32_t FIBFOURTY_TTTCv[66048 EXT_STRUCT_CTRL(FIBFOURTY_TTTC)] /* [0,65535] */;
  uint32_t FIBFOURTY_TTTFM /* [1,1032] */;
  uint32_t FIBFOURTY_TTTFMI[1032 EXT_STRUCT_CTRL(FIBFOURTY_TTTFM)] /* [1,1032] */;
  uint32_t FIBFOURTY_TTTFME[1032 EXT_STRUCT_CTRL(FIBFOURTY_TTTFM)] /* [1,66048] */;
  uint32_t FIBFOURTY_TTTF /* [0,66048] */;
  uint32_t FIBFOURTY_TTTFv[66048 EXT_STRUCT_CTRL(FIBFOURTY_TTTF)] /* [0,65535] */;
  uint32_t FIBFOURTY_TRIGLCM /* [1,16] */;
  uint32_t FIBFOURTY_TRIGLCMI[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGLCM)] /* [1,16] */;
  uint32_t FIBFOURTY_TRIGLCME[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGLCM)] /* [1,16] */;
  uint32_t FIBFOURTY_TRIGLC /* [0,16] */;
  uint32_t FIBFOURTY_TRIGLCv[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGLC)] /* [0,65535] */;
  uint32_t FIBFOURTY_TRIGLFM /* [1,16] */;
  uint32_t FIBFOURTY_TRIGLFMI[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGLFM)] /* [1,16] */;
  uint32_t FIBFOURTY_TRIGLFME[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGLFM)] /* [1,16] */;
  uint32_t FIBFOURTY_TRIGLF /* [0,16] */;
  uint32_t FIBFOURTY_TRIGLFv[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGLF)] /* [0,65535] */;
  uint32_t FIBFOURTY_TRIGTCM /* [1,16] */;
  uint32_t FIBFOURTY_TRIGTCMI[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGTCM)] /* [1,16] */;
  uint32_t FIBFOURTY_TRIGTCME[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGTCM)] /* [1,16] */;
  uint32_t FIBFOURTY_TRIGTC /* [0,16] */;
  uint32_t FIBFOURTY_TRIGTCv[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGTC)] /* [0,65535] */;
  uint32_t FIBFOURTY_TRIGTFM /* [1,16] */;
  uint32_t FIBFOURTY_TRIGTFMI[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGTFM)] /* [1,16] */;
  uint32_t FIBFOURTY_TRIGTFME[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGTFM)] /* [1,16] */;
  uint32_t FIBFOURTY_TRIGTF /* [0,16] */;
  uint32_t FIBFOURTY_TRIGTFv[16 EXT_STRUCT_CTRL(FIBFOURTY_TRIGTF)] /* [0,65535] */;

} EXT_STR_h101_FIBFOURTY;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBFOURTY_onion_t
{
  /* RAW */
  uint32_t FIBFOURTY_TBLCM;
  uint32_t FIBFOURTY_TBLCMI[1032 /* FIBFOURTY_TBLCM */];
  uint32_t FIBFOURTY_TBLCME[1032 /* FIBFOURTY_TBLCM */];
  uint32_t FIBFOURTY_TBLC;
  uint32_t FIBFOURTY_TBLCv[66048 /* FIBFOURTY_TBLC */];
  uint32_t FIBFOURTY_TBLFM;
  uint32_t FIBFOURTY_TBLFMI[1032 /* FIBFOURTY_TBLFM */];
  uint32_t FIBFOURTY_TBLFME[1032 /* FIBFOURTY_TBLFM */];
  uint32_t FIBFOURTY_TBLF;
  uint32_t FIBFOURTY_TBLFv[66048 /* FIBFOURTY_TBLF */];
  uint32_t FIBFOURTY_TBTCM;
  uint32_t FIBFOURTY_TBTCMI[1032 /* FIBFOURTY_TBTCM */];
  uint32_t FIBFOURTY_TBTCME[1032 /* FIBFOURTY_TBTCM */];
  uint32_t FIBFOURTY_TBTC;
  uint32_t FIBFOURTY_TBTCv[66048 /* FIBFOURTY_TBTC */];
  uint32_t FIBFOURTY_TBTFM;
  uint32_t FIBFOURTY_TBTFMI[1032 /* FIBFOURTY_TBTFM */];
  uint32_t FIBFOURTY_TBTFME[1032 /* FIBFOURTY_TBTFM */];
  uint32_t FIBFOURTY_TBTF;
  uint32_t FIBFOURTY_TBTFv[66048 /* FIBFOURTY_TBTF */];
  uint32_t FIBFOURTY_TTLCM;
  uint32_t FIBFOURTY_TTLCMI[1032 /* FIBFOURTY_TTLCM */];
  uint32_t FIBFOURTY_TTLCME[1032 /* FIBFOURTY_TTLCM */];
  uint32_t FIBFOURTY_TTLC;
  uint32_t FIBFOURTY_TTLCv[66048 /* FIBFOURTY_TTLC */];
  uint32_t FIBFOURTY_TTLFM;
  uint32_t FIBFOURTY_TTLFMI[1032 /* FIBFOURTY_TTLFM */];
  uint32_t FIBFOURTY_TTLFME[1032 /* FIBFOURTY_TTLFM */];
  uint32_t FIBFOURTY_TTLF;
  uint32_t FIBFOURTY_TTLFv[66048 /* FIBFOURTY_TTLF */];
  uint32_t FIBFOURTY_TTTCM;
  uint32_t FIBFOURTY_TTTCMI[1032 /* FIBFOURTY_TTTCM */];
  uint32_t FIBFOURTY_TTTCME[1032 /* FIBFOURTY_TTTCM */];
  uint32_t FIBFOURTY_TTTC;
  uint32_t FIBFOURTY_TTTCv[66048 /* FIBFOURTY_TTTC */];
  uint32_t FIBFOURTY_TTTFM;
  uint32_t FIBFOURTY_TTTFMI[1032 /* FIBFOURTY_TTTFM */];
  uint32_t FIBFOURTY_TTTFME[1032 /* FIBFOURTY_TTTFM */];
  uint32_t FIBFOURTY_TTTF;
  uint32_t FIBFOURTY_TTTFv[66048 /* FIBFOURTY_TTTF */];
  uint32_t FIBFOURTY_TRIGLCM;
  uint32_t FIBFOURTY_TRIGLCMI[16 /* FIBFOURTY_TRIGLCM */];
  uint32_t FIBFOURTY_TRIGLCME[16 /* FIBFOURTY_TRIGLCM */];
  uint32_t FIBFOURTY_TRIGLC;
  uint32_t FIBFOURTY_TRIGLCv[16 /* FIBFOURTY_TRIGLC */];
  uint32_t FIBFOURTY_TRIGLFM;
  uint32_t FIBFOURTY_TRIGLFMI[16 /* FIBFOURTY_TRIGLFM */];
  uint32_t FIBFOURTY_TRIGLFME[16 /* FIBFOURTY_TRIGLFM */];
  uint32_t FIBFOURTY_TRIGLF;
  uint32_t FIBFOURTY_TRIGLFv[16 /* FIBFOURTY_TRIGLF */];
  uint32_t FIBFOURTY_TRIGTCM;
  uint32_t FIBFOURTY_TRIGTCMI[16 /* FIBFOURTY_TRIGTCM */];
  uint32_t FIBFOURTY_TRIGTCME[16 /* FIBFOURTY_TRIGTCM */];
  uint32_t FIBFOURTY_TRIGTC;
  uint32_t FIBFOURTY_TRIGTCv[16 /* FIBFOURTY_TRIGTC */];
  uint32_t FIBFOURTY_TRIGTFM;
  uint32_t FIBFOURTY_TRIGTFMI[16 /* FIBFOURTY_TRIGTFM */];
  uint32_t FIBFOURTY_TRIGTFME[16 /* FIBFOURTY_TRIGTFM */];
  uint32_t FIBFOURTY_TRIGTF;
  uint32_t FIBFOURTY_TRIGTFv[16 /* FIBFOURTY_TRIGTF */];

} EXT_STR_h101_FIBFOURTY_onion;

/*******************************************************/

#define EXT_STR_h101_FIBFOURTY_ITEMS_INFO(ok,si,offset,struct_t,printerr) do { \
  ok = 1; \
  /* RAW */ \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBLCM,                 UINT32,\
                    "FIBFOURTY_TBLCM",1032,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBLCMI,                UINT32,\
                    "FIBFOURTY_TBLCMI",                "FIBFOURTY_TBLCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBLCME,                UINT32,\
                    "FIBFOURTY_TBLCME",                "FIBFOURTY_TBLCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBLC,                  UINT32,\
                    "FIBFOURTY_TBLC",66048,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBLCv,                 UINT32,\
                    "FIBFOURTY_TBLCv",                 "FIBFOURTY_TBLC",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBLFM,                 UINT32,\
                    "FIBFOURTY_TBLFM",1032,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBLFMI,                UINT32,\
                    "FIBFOURTY_TBLFMI",                "FIBFOURTY_TBLFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBLFME,                UINT32,\
                    "FIBFOURTY_TBLFME",                "FIBFOURTY_TBLFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBLF,                  UINT32,\
                    "FIBFOURTY_TBLF",66048,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBLFv,                 UINT32,\
                    "FIBFOURTY_TBLFv",                 "FIBFOURTY_TBLF",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBTCM,                 UINT32,\
                    "FIBFOURTY_TBTCM",1032,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBTCMI,                UINT32,\
                    "FIBFOURTY_TBTCMI",                "FIBFOURTY_TBTCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBTCME,                UINT32,\
                    "FIBFOURTY_TBTCME",                "FIBFOURTY_TBTCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBTC,                  UINT32,\
                    "FIBFOURTY_TBTC",66048,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBTCv,                 UINT32,\
                    "FIBFOURTY_TBTCv",                 "FIBFOURTY_TBTC",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBTFM,                 UINT32,\
                    "FIBFOURTY_TBTFM",1032,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBTFMI,                UINT32,\
                    "FIBFOURTY_TBTFMI",                "FIBFOURTY_TBTFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBTFME,                UINT32,\
                    "FIBFOURTY_TBTFME",                "FIBFOURTY_TBTFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBTF,                  UINT32,\
                    "FIBFOURTY_TBTF",66048,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TBTFv,                 UINT32,\
                    "FIBFOURTY_TBTFv",                 "FIBFOURTY_TBTF",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTLCM,                 UINT32,\
                    "FIBFOURTY_TTLCM",1032,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTLCMI,                UINT32,\
                    "FIBFOURTY_TTLCMI",                "FIBFOURTY_TTLCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTLCME,                UINT32,\
                    "FIBFOURTY_TTLCME",                "FIBFOURTY_TTLCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTLC,                  UINT32,\
                    "FIBFOURTY_TTLC",66048,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTLCv,                 UINT32,\
                    "FIBFOURTY_TTLCv",                 "FIBFOURTY_TTLC",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTLFM,                 UINT32,\
                    "FIBFOURTY_TTLFM",1032,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTLFMI,                UINT32,\
                    "FIBFOURTY_TTLFMI",                "FIBFOURTY_TTLFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTLFME,                UINT32,\
                    "FIBFOURTY_TTLFME",                "FIBFOURTY_TTLFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTLF,                  UINT32,\
                    "FIBFOURTY_TTLF",66048,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTLFv,                 UINT32,\
                    "FIBFOURTY_TTLFv",                 "FIBFOURTY_TTLF",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTTCM,                 UINT32,\
                    "FIBFOURTY_TTTCM",1032,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTTCMI,                UINT32,\
                    "FIBFOURTY_TTTCMI",                "FIBFOURTY_TTTCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTTCME,                UINT32,\
                    "FIBFOURTY_TTTCME",                "FIBFOURTY_TTTCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTTC,                  UINT32,\
                    "FIBFOURTY_TTTC",66048,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTTCv,                 UINT32,\
                    "FIBFOURTY_TTTCv",                 "FIBFOURTY_TTTC",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTTFM,                 UINT32,\
                    "FIBFOURTY_TTTFM",1032,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTTFMI,                UINT32,\
                    "FIBFOURTY_TTTFMI",                "FIBFOURTY_TTTFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTTFME,                UINT32,\
                    "FIBFOURTY_TTTFME",                "FIBFOURTY_TTTFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTTF,                  UINT32,\
                    "FIBFOURTY_TTTF",66048,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TTTFv,                 UINT32,\
                    "FIBFOURTY_TTTFv",                 "FIBFOURTY_TTTF",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGLCM,               UINT32,\
                    "FIBFOURTY_TRIGLCM",16,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGLCMI,              UINT32,\
                    "FIBFOURTY_TRIGLCMI",              "FIBFOURTY_TRIGLCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGLCME,              UINT32,\
                    "FIBFOURTY_TRIGLCME",              "FIBFOURTY_TRIGLCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGLC,                UINT32,\
                    "FIBFOURTY_TRIGLC",16,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGLCv,               UINT32,\
                    "FIBFOURTY_TRIGLCv",               "FIBFOURTY_TRIGLC",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGLFM,               UINT32,\
                    "FIBFOURTY_TRIGLFM",16,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGLFMI,              UINT32,\
                    "FIBFOURTY_TRIGLFMI",              "FIBFOURTY_TRIGLFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGLFME,              UINT32,\
                    "FIBFOURTY_TRIGLFME",              "FIBFOURTY_TRIGLFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGLF,                UINT32,\
                    "FIBFOURTY_TRIGLF",16,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGLFv,               UINT32,\
                    "FIBFOURTY_TRIGLFv",               "FIBFOURTY_TRIGLF",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGTCM,               UINT32,\
                    "FIBFOURTY_TRIGTCM",16,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGTCMI,              UINT32,\
                    "FIBFOURTY_TRIGTCMI",              "FIBFOURTY_TRIGTCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGTCME,              UINT32,\
                    "FIBFOURTY_TRIGTCME",              "FIBFOURTY_TRIGTCM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGTC,                UINT32,\
                    "FIBFOURTY_TRIGTC",16,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGTCv,               UINT32,\
                    "FIBFOURTY_TRIGTCv",               "FIBFOURTY_TRIGTC",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGTFM,               UINT32,\
                    "FIBFOURTY_TRIGTFM",16,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGTFMI,              UINT32,\
                    "FIBFOURTY_TRIGTFMI",              "FIBFOURTY_TRIGTFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGTFME,              UINT32,\
                    "FIBFOURTY_TRIGTFME",              "FIBFOURTY_TRIGTFM",0/*flags*/); \
  EXT_STR_ITEM_INFO2_LIM(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGTF,                UINT32,\
                    "FIBFOURTY_TRIGTF",16,0/*flags*/); \
  EXT_STR_ITEM_INFO2_ZZP(ok,si,offset,struct_t,printerr,\
                     FIBFOURTY_TRIGTFv,               UINT32,\
                    "FIBFOURTY_TRIGTFv",               "FIBFOURTY_TRIGTF",0/*flags*/); \
  \
} while (0);
#endif/*__GUARD_H101_FIBFOURTY_EXT_H101_FIBFOURTY_H__*/

/*******************************************************/
