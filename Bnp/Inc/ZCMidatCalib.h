
#ifndef __ZCMIDAT_CALIB_H
#define __ZCMIDAT_CALIB_H

#include "AmortMidatCalib.h"
#include "CPDCalib.h"

/*	Calibrate lgm on ZC swaptions : main function */
Err ZCMidat_cpd_calib_diagonal(
    char* yc_name,               /*	Name of the yield curve */
    char* vol_curve_name,        /*	Name of the market vol curve */
    char* default_ref_rate_name, /*	Name of the reference rate */
    Err (*get_cash_vol)(         /*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    double vol_shift,
    int    shift_type,    /*	0:	Additive
                                          1:	Multiplicative */
                          /*	If ex_date is NULL,
                          exercise dates will be generated 2bd before start */
    int num_ex_dates,     /*	Exercise dates,
                                                          all supposed to be on or after today */
    long* ex_date,        /*	Supposed to be sorted
                                                          NULL = 2bd before each coupon */
    long    start_date,   /*	Start date of the amortised swap */
    long    end_date,     /*	End date for the amortised swap */
    double* long_strike,  /*	Diagonal swaption strikes
                                                                  NULL = ATM */
    double* short_strike, /*	Short swaption strikes
                                                                  NULL = ATM */
    int strike_type,      /*	0: ATM
                                                  1: CASH
                                                  2: SWAP
                                                  3: STD */
    double* diag_prices,  /* Diagonal Prices */
    double  max_std_long,
    double  max_std_short,
    char*   ref_rate_name, /*	Name of the reference rate */
    char*   swaption_freq, /*	Frequency, basis and ref. rate of underlying swaptions */
    char*   swaption_basis,
    int     fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                to diagonal */
    int one_f_equi,     /*	1F equivalent flag:
                                                if set to 1, then 2F lambda will calibrate
                                                to the cap priced within calibrated 1F
                                                with the given lambda */
    int skip_last,      /*	If 1, the last option is disregarded
                                                and the forward volatility is flat from option
                                                n-1 */
    double* lambda,     /*	Lambda: may be changed in the process */
    int     one2F,      /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double   alpha,
    double   gamma,
    double   rho,
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data); /*	NULL = don't save calibration instrument data */

#endif