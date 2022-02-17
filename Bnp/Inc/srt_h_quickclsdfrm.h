/* ===================================================================================
   FILENAME:      srt_h_quickclsdfrm.h

   PURPOSE:       Compute quick prices for saptions  , cap floors for a generic
   Cheyette Beta model the formula is a complete approximation and offers
   absolutely no accuracy. It is just used in the calibration with FIXED POINT
   as a way of getting a Beta sensitivity...
   ===================================================================================
 */

#ifndef SRT_H_QUICKCLSDFRM_H
#define SRT_H_QUICKCLSDFRM_H

double srt_f_quick_beta_bond_option(int n, double bond_strike, TermStruct *ts,
                                    double fixing_time, double *coupon,
                                    double *pay_time, double *df,
                                    SrtReceiverType rec_pay,
                                    SrtMdlType mdl_type, SrtMdlDim eModelDim,
                                    String szYieldCurveName);

/* ----------------------------------------------------------------------------
 */

double srt_f_quick_beta_capfloor(
    TermStruct *ts, double *fixing_time,
    double *period_time,         /* period_time[0] is the strike payment time */
    double *df, double *payment, /* This is cvg * ( cash_fwd + spread ) */
    int num_caplets, SrtReceiverType rec_pay, SrtMdlType mdl_type,
    SrtMdlDim eModelDim, String yc_name);

/* ------------------------------------------------------------------------------
 */

#endif
