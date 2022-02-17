/* ===========================================================================

         FILENAME:     swp_h_ycinstr.h

     PURPOSE:      A yield curve can be thought of as composed
                       of sets of instruments of various types
                       (swap      , cash      , futs);
                                         here they are.

   ===========================================================================
 */

#ifndef SWP_H_YCINSTR_H
#define SWP_H_YCINSTR_H

#include "swp_h_swapdp.h"
#include "swp_h_zc_types.h"

typedef struct {
  BasisCode basis_code;
  Date start;
  Date end;
  Date settle;
} FutDP;
typedef FutDP FraDP;

typedef struct {
  BasisCode basis_code;
  Date start;
  Date end;
} CashDP;
typedef CashDP ToyDP;

typedef union {
  FutDP futdates;
  SwapDP swapdates;
  CashDP cashdates;
  FraDP fradates;
} DateParam;

typedef enum YCInstrType_ {
  SWAPINSTR,
  CASHINSTR,
  FUT1MINSTR,
  FUT3MINSTR,
  LASTYCINSTRTYPE
} YCInstrType;

typedef struct {
  YCInstrType type; /* one of Future      , Swap      , Cash */
  double rate;
  DateParam dp;
} YCInstr;

Err string_to_YCInstr(String iname, YCInstr *instr, double *rate,
                      SrtCcyParam *ccy_param, Date spot);
Err YC_get_Instr(YC_Obj *yc, YCInstr *instr, int i);
Err YC_set_Instr(YC_Obj *yc, YCInstr *instr, int i);
int YC_index_Instr(YC_Obj *yc, YCInstr *in1);
int YCInstr_comp(YCInstr *in1, YCInstr *in2);

#endif
