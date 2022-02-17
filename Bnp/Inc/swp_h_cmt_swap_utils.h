#ifndef SWP_H_CMT_SWAP_UTILS_H
#define SWP_H_CMT_SWAP_UTILS_H

/* ====
Please remember that  ,
- in utdtsched.h:

typedef struct{
        Date *date;
        int len;
        SwapDateDir dir;
        DateListType type;      	BROKEN or NOTBROKEN
        Date prev;				date before 2nd date  , if
BROKEN }DateList  ,SrtDateList;

- in srt_h_zc_utils :

typedef struct{
                double 		*d;
                int 		len;
                }Dlist;

Dlist time_list(DateList  , Date);
Dlist cvg_list(DateList  , BasisCode);

- in srt_h_gen_constant :

typedef enum StructType {BOND_OPTION  ,SWAPTION  ,CAPFLOOR  ,
                        SWAP  ,BOND  ,
                        LASTSTRUCTTYPE}StructType;
   ==== */

#define DEFAULT_GEN_NOTIONAL 1.0e8

typedef enum CMRateType { SWAP_RATE, TREASURY_RATE } CMRateType;

typedef struct gen_swap_leg {
  LegType leg_type;
  SrtReceiverType rec_pay;
  SwapDP sdp;
  DateList dl;
  Dlist time;
  Dlist temp_fwd;
  Dlist fwd;
  Dlist vol;
  Dlist cpn; /* == convexity adjusted forward or fixed or fra...*/
  double spread;
  Dlist cvg;
  double initial_exchange;
  double final_exchange;
  Dlist payment;
  Dlist df;
  double leg_value;
  int leg_length;
  long spot_date;
  long today;
} GenSwapLeg;

typedef struct gen_full_swap {
  SwapType swap_type;
  GenSwapLeg leg[2];
  long today;
  long spotdate;
  long enddate;
  double notional;
  double swap_value;
} GenFullSwap;

Dlist const_list(DateList dl, double constant);

Dlist fra_list(DateList dl, SrtBasisCode basiscode, SrtCurvePtr yc_crv);

Dlist fwd_treas_list(DateList dl, SrtCurvePtr cmt_crv);

Dlist fwd_swap_list(DateList dl, SrtCurvePtr cmt_crv);

Dlist swap_vol_list(DateList dl, SrtCurvePtr cmt_crv);

Dlist treas_vol_list(Dlist swap_vol, Dlist swap_rate, Dlist treas_rate,
                     SrtCurvePtr cmt_crv);

Dlist conv_fwd_list(long today, DateList dl, Dlist fwd_rate, Dlist vol,
                    SrtCurvePtr cmt_crv, CMRateType rate_type);

#endif
