/**
        File_Name:	swp_h_swap_compute.h
        Author:		G. Amblard

**/
#ifndef SWP_H_SWAP_COMPUTE_H
#define SWP_H_SWAP_COMPUTE_H

typedef enum SwapOutput {
  COMPUTE_PV,
  COMPUTE_FWD_PV,
  COMPUTE_FWD_RATE,
  COMPUTE_LEVEL,
  COMPUTE_MARGIN,
  COMPUTE_IRR,
  COMPUTE_FWD_IRR,
  COMPUTE_DURATION,
  COMPUTE_CONVEXITY,
  COMPUTE_THIRDMOMENT,
  COMPUTE_MODIFIED_DURATION,
  COMPUTE_MODIFIED_CONVEXITY,
  COMPUTE_MODIFIED_MATCHING_RATIO,
  COMPUTE_SPREAD,
  COMPUTE_FWD_SPREAD
} SwapOutput;

enum swaption_types { COMPUTE_STANDARD_BOND_OPTION, COMPUTE_STANDARD_SWAPTION };

enum capfloor_types { COMPUTE_STANDARD_CAPFLOOR, COMPUTE_QUANTO_CAPFLOOR };

#define GREEKS(val)                                                            \
  ((val) == PREMIUM || (val) == FWD_DELTA_S || (val) == SPOT_DELTA_S ||        \
   (val) == FWD_DELTA_K || (val) == SPOT_DELTA_K || (val) == FWD_GAMMA ||      \
   (val) == SPOT_GAMMA || (val) == THETA || (val) == THETA_1D ||               \
   (val) == VEGA)

double get_greek(struct greek_struct greeks, SrtGreekType greek);

/* ------------------- Function to interpret the swap info required
 * --------------------- */
Err swp_f_interp_swap_message(String messageStr, SwapOutput *output);

Err SWAP_compute(SrtCurvePtr crv, Message mess, Arg_Obj *arg);

Err swap_unwind_compute(SrtCurvePtr crv, Arg_Obj *arg, SRT_Boolean b, int eod);

double black_scholes_tmp(double fwd_price, double strike, double vol,
                         double mat, double disc, SrtReceiverType rec_or_pay,
                         struct greek_struct *greeks);

#endif
