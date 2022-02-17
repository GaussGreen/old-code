#ifndef __CTS_QTO_PROD_STRUCT_H
#define __CTS_QTO_PROD_STRUCT_H

#include "FundLegProdStruct.h"
#include "LGMQuantoUnd.h"
#include "RangeAccrualProdStruct.h"

#define CTS_QTO_NCPN 512
#define CTS_QTO_MAX_NUM_OBS 512

/*	------------------------------------------------------ */
/*	Structures and functions for callable Time Swap Quanto */
/*	------------------------------------------------------ */

/*	Structures and functions for the calls */

/*	Call */
typedef struct {
  int pay_rec; /*	0: rec pd upon exercise      , 1: pay pd upon exercise
                */
  /*	Specs */
  long ex_date;     /*	Exercise date */
  double ex_time;   /*	Exercise time */
  int fund_idx;     /*	Index of the first funding coupon to be called */
  int num_fund_cpn; /*	Number of funding coupons called (including redemption)
                     */
  int ra_idx;       /*	Index of the first ra coupon to be called */
  int num_ra_cpn;   /*	Number of cf coupons called (excluding redemption) */
  /*	Fee upon exercise */
  long set_date;   /*	Fee payment date */
  double set_time; /*	Fee payment time */
  double fee;      /*	Amount of the fee */
} ctsqto_call, *CTSQTO_CALL;

/*	Callable Time Swap Quanto */
typedef struct {
  RangeAccrualStruct *ra_leg;
  funding_leg *fund_leg;
  int num_calls;
  ctsqto_call *call;
} ctsqto_str, *CTSQTO_STR;

Err ctsqto_fill_calls(
    /*	Exercises before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag,           /*	0: I      , 1: E */
    int ncall, int pay_rec, /*	0: rec pd      , 1: pay pd */
    long *ex_date, long *set_date, double *fee, CTSQTO_STR ctsqto);

/*	Check dates consistency */
Err ctsqto_check_calls(CTSQTO_STR ctsqto);

/*	Free */
Err ctsqto_free_calls(CTSQTO_STR ctsqto);

Err ctsqto_free_all_struct(CTSQTO_STR ctsqto, LGMQTO_UND und, int call_feat,
                           LGMQTO_ADI_ARG adi_arg);

Err ctsqto_calc_mdl_iv_fwd(CTSQTO_STR ctsqto, LGMQTO_UND und,
                           LGMQTO_ADI_ARG adi_arg,

                           int num_hermite,
                           /*	Result */
                           double *premium);

Err ctsqto_launch_adi(CTSQTO_STR cstqto, LGMQTO_ADI_ARG adi_arg,
                      /*	Result */
                      double *prem);

Err ctsqto_fill_check_all_struct(
    /*	Today's date */
    long today,
    /*	The underlying */
    int use_calib, /*	0: use fx3dund      , 1: calibrate */
    /*		if calib */
    char *dom_yc,         /*	yc */
    char *dom_vc,         /*	vc (only if calib) */
    char *dom_ref,        /*	ref rate (only if calib) */
    char *dom_swap_freq,  /*	swap freq (only if calib) */
    char *dom_swap_basis, /*	swap basis (only if calib) */
    double dom_lambda,    /*	lambda if unique */

    char *for_yc,         /*	yc */
    char *for_vc,         /*	vc */
    char *for_ref,        /*	ref rate */
    char *for_swap_freq,  /*	swap freq */
    char *for_swap_basis, /*	swap basis */
    double for_lambda,    /*	lambda if unique */

    int forcalib, /*	0 : RA Und      , 1 : Diag */

    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*		if no calilb */
    char *fx3dund,
    /*	The structure */

    /*		funding */
    double fund_not, int fund_ncpn, long *fund_fix, long *fund_start,
    long *fund_pay, char **fund_basis, double *fund_spr, double *fund_mrg,
    /*		ra */
    double ra_not, int ra_cpn_type, int n_fxvol_dates, long *fxvol_dates,
    double *fxvol, double *qtocorrel,

    int typeVol,

    int n_periods, long *dates, /* n_periods + 1 */
    double *cpns, char *basis, int *ra_nfixings, long **ra_fixingdates,
    double **ra_fixings, /*	Past coupon fixing if relevant */

    // RA floating coupons
    int ra_float_refrate_is_dom_for, char *ra_float_refrate,
    long *ra_float_fixl, double *ra_float_past_fixings,
    double *ra_float_gearings,

    double *upper_barr, double *lower_barr, char *recpay, char *ra_refrate,
    double c_spread, int obs_freq, double rho_df,

    // Params for floating coupon
    double correl_start, double correl_end, int float_adj_strike,

    /*		calls */
    int ncall, int pay_rec, /*	0: rec pd      , 1: pay pd */
    long *ex_date, long *set_date, double *fee,
    /*	Numerical params */
    int req_stp, int req_stpx,
    /*	Calib params */
    int dom_force_atm, /*	force atm calib domestic und */
    int for_force_atm, /*	force atm calib foreign und */
    double max_std_long, double max_std_short,
    int fix_lambda, /*	0: calib lambda to cap      , 1: fix lambda calib
                                                to diagonal */
    int one_f_equi, /*	1F equivalent flag:
                                                if set to 1      , then 2F
                   lambda will calibrate to the cap priced within calibrated
                   1F with the given lambda */
    int skip_last,  /*	If 1      , the last option is disregarded
                                                and the forward volatility is
                   flat from option  n-1 */
    long *fx_mkt_vol_date, double *fx_mkt_vol, int num_fx_mkt_vol,

    double *corr_times, double *correl_dom_for, double *correl_dom_fx,
    double *correl_for_fx, int corr_n_times,
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I      , 1: E */
    int eod_ex_flag,  /*	0: I      , 1: E */

    //------Spot Lag-------------
    int fund_spot_lag, int ra_spot_lag,

    /*	Results */
    CTSQTO_STR ctsqto, LGMQTO_UND und,
    int *call_feat, /*	0: No callable feature to be valued
                        1: Callable feature to be valued through adi */
    LGMQTO_ADI_ARG adi_arg);

#endif