/* =============================================================================

   FILENAME:    swp_h_swap_simple.h

   PURPOSE:     Provide a few function for simple swap calculations:
                                - PV
                                - Fwd Rate
                                - IRR
                                - Duration
                                - Convexity
   ============================================================================= */
#ifndef SWP_H_SWAP_SIMPLE_H
#define SWP_H_SWAP_SIMPLE_H

typedef struct simple_swap_obj
{
    SwapDP   sdp;
    DateList dl;
    Dlist    t;
    Dlist    cpn;
} SimpleSwap;

/* Swaps calculation functions based on swapdp */

Err zcswap(
    SwapDP*    sdp,
    double     strike,
    double     ini,
    double     fin,
    SwapOutput m,
    String     ycname,
    Date       value_date,
    double*    answer);

Err zcsens(String tenor, SrtCurvePtr crv, double* answer);

Err bond_modified_duration(SwapDP* sdp, double coupon, double pv, double* duration);

Err bond_modified_convexity(SwapDP* sdp, double coupon, double pv, double* convexity);

SimpleSwap make_SimpleSwap(
    SwapDP* sdp, double strike, double ini, double fin, Date today, StructType t);

int free_inSimpleSwap(SimpleSwap* s);

double value_SimpleSwap(SimpleSwap* s, SrtCurvePtr crv);

double frate_SimpleSwap(SimpleSwap* s, SrtCurvePtr crv);

Err srt_f_SimpleSwap(
    long    start,
    long    end_nfp,
    String  compStr,
    String  basisStr,
    double  coupon,
    double  initial,
    double  final,
    String  ycname,
    String  strMessage,
    long    value_date,
    double* answer);

#endif
