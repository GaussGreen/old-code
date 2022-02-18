
#ifndef __CCF_CALLER_H
#define __CCF_CALLER_H

#include "CCFProdStruct.h"

/*	Caller for callable yield curve options */
/*	------------------------------- */

Err ccf_caller(
    /*	Today's date */
    long today,
    /*	The underlying */
    int use_calib, /*	0: use lgm2fund, 1: calibrate */
    /*		if calib */
    char*   yc,         /*	yc */
    char*   vc,         /*	vc */
    char*   ref,        /*	ref rate (only if calib) */
    char*   swap_freq,  /*	swap freq (only if calib) */
    char*   swap_basis, /*	swap basis (only if calib) */
    int     lam_ts,     /*	0: use unique lambda, 1: use ts */
    double  lambda,     /*	lambda if unique */
    int     tsnlam,     /*	number of lambdas if ts */
    double* tslamtime,  /*	lambda times i.e. (date - today) / 365 if ts */
    double* tslam,      /*	corresponding lambdas if ts */
    double  alpha,      /*	alpha */
    double  gamma,      /*	gamma */
    double  rho,        /*	rho */
    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    /*		if no calilb */
    char* lgm2dund,
    /*	The structure */
    long start_date, /*	Date at which initial notional exchange occurs */
    /*		funding */
    int     fund_ccy,    /*	0: domestic, 1: other */
    double  fund_not,    /*	If different from domestic or foreign (fund_ccy = 1) */
    char*   fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 1) */
    double  fx_fund_dom, /*	If different from domestic or foreign (fund_ccy = 1) 2 bd fwd */
    long    fx_fund_dom_spot_date,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    double* fund_fix_cpn, /*	Past coupon fixing if relevant,
                                          includes spr, but not mrg, cvg and notional */
    /*		cf */
    double  cf_not,
    int     cf_ncpn,
    long*   cf_fix,
    long*   cf_start,
    long*   cf_pay,
    char**  cf_basis,
    char**  cf_cms_tenor1,
    char*   cf_cms_freq1,
    char*   cf_cms_basis1,
    double* cf_cms_spread1,
    char**  cf_cms_tenor2,
    char*   cf_cms_freq2,
    char*   cf_cms_basis2,
    double* cf_cms_spread2,
    long    spread_vol_n,
    double* spread_vol_time,
    double* spread_vol_floor,
    double* spread_vol_cap,

    double* spread_slbeta1, /* Shifted log beta on the CMS1 */
    double* spread_slbeta2, /* Shifted log beta on the CMS2 */

    double* cf_alpha,
    double* cf_beta,
    double* cf_gamma,
    int*    cf_floored,
    double* cf_floor,
    int*    cf_capped,
    double* cf_cap,
    double* cf_fix_cpn, /*	Past coupon fixing if relevant */

    int      cf_nopt,      /* Number of spread options */
    double** cf_notopt,    /* Notional of the spread options */
    double** cf_strikeopt, /* spread option strikes */
    int**    cf_typeopt,   /* spread option type 0 call 1 put */

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,
    /*	Numerical params */
    int  req_stp,
    int  req_stpx,
    long req_paths,
    /*	Calib params */
    int    force_atm,
    double max_std_long,
    double max_std_short,
    int    fix_lambda,         /*	0: calib lambda to cap, 1: fix lambda calib
                                                               to diagonal */
    int    cal_vol_shift_type, /*	vol shift type for volatility */
    double cal_vol_shift,      /*	vol shift */
    double cal_lambda_shift,   /*	shift on lambda after calibration */
    int    one_f_equi,         /*	1F equivalent flag:
                                                               if set to 1, then 2F lambda will calibrate
                                                               to the cap priced within calibrated 1F
                                                               with the given lambda */
    int skip_last, /*	If 1, the last option is disregarded and the forward volatility is flat from
                      option n-1 */
    double long_prec,  /*	Precision on primary instruments */
    double short_prec, /*	Precision on secondary instruments */
    double min_fact,   /*	Maximum down jump on variance */
    double max_fact,   /*	Maximum up jump on variance */
    int    use_jumps,  /*	Allow vol term structure to jump */
    int    proba_weight,

    int calc_fwdiv, /*	output the model and market fwdiv  */
    int adjust_fee, /*	adjust the fee to price correctly the fwd iv */

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_pay_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */
    /*	CMS can be valued with smile */
    int cms_adj,                /*	1: adjust for CMS effect
                                                0: don't */
    int cms_for_iv_only,        /*	1: adjust for CMS effect only IV
                                                        0: adjust both IV and call */
    int use_cms_smile,          /*	1: value CMS with ATM
                                                        0: value CMS with smile */
    int cms_vol_adj,            /*	1: adjust for CMS vol effect
                                                        0: don't */
    double           cms_beta1, /*	How adjustment varies with rates */
    double           cms_beta2,
    int              num_strikes_in_vol, /*	Array of strikes in vol matrix */
    double*          strikes_in_vol,
    SrtDiffusionType vol_type, /*	Type of vol in matrix, SRT_NORMAL or SRT_LOGNORMAL */
    int              cash_vol, /*	1: matrix is a cash vol
                                               0: matrix is a swap vol */
    int cvg_sv,                /*	1: coverging model on spread vol/correlation
                                               0: sliding model on spread vol/correlation */
    int is_corr,               /*	1: spread vols are correlations
                                               0: spread vols are vols */

    int use_SL,   /*  1: use Shifted Log
                              0: don't use				*/
    int calib_SL, /*  1: Calibrate the Shifted Log models on CMS1 and on CMS2 *
                                  0: use the Shifted Log beta given */

    double NbStdcalib_SL, /*  Nb of std for the calibration of to the skew */

    int calib_correl_SL, /*	1: Calibrate the correlation between the two SL to get the same ATM
                                                         normal spread vol
                                                 0: use the normal spread correl for the sl correl
                          */

    int use_cfoptions, /*  1: Use the spread options and don't take into account
                                                       the floor and cap in the cf coupons
                                               0: take into account the floor and cap in the cf
                          coupons */

    /*	Exercised flag */
    int    exercised,   /*	Flag */
    long   ex_date_ex,  /*	Date when exercised */
    long   ex_date_set, /*	Corresponding settlement date */
    double ex_fee,      /*	Corresponding fee */
    /*	Results */
    double* fund_val,  /*	Value of the funding leg */
    double* cf_val,    /*	Value of the Power Dual leg */
    double* call_val,  /*	Value of the callable feature */
    int     export_ts, /*	1: Export TS, 0: don't */
    CCF_UND und_exp);

Err ccf_calc_mdl_iv_fwd(
    CCF_STR     ccf,
    CCF_UND     und,
    CCF_ADI_ARG adi_arg,

    int num_hermite,
    /*	Result */
    double* premium);

Err ccf_caller_LambdaShift(
    /*	Today's date */
    long today,
    /*	The underlying */
    int use_calib, /*	0: use lgm2fund, 1: calibrate */
    /*		if calib */
    char*   yc,         /*	yc */
    char*   vc,         /*	vc */
    char*   ref,        /*	ref rate (only if calib) */
    char*   swap_freq,  /*	swap freq (only if calib) */
    char*   swap_basis, /*	swap basis (only if calib) */
    int     lam_ts,     /*	0: use unique lambda, 1: use ts */
    double  lambda,     /*	lambda if unique */
    int     tsnlam,     /*	number of lambdas if ts */
    double* tslamtime,  /*	lambda times i.e. (date - today) / 365 if ts */
    double* tslam,      /*	corresponding lambdas if ts */
    double  alpha,      /*	alpha */
    double  gamma,      /*	gamma */
    double  rho,        /*	rho */
    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    /*		if no calilb */
    char* lgm2dund,
    /*	The structure */
    long start_date, /*	Date at which initial notional exchange occurs */
    /*		funding */
    int     fund_ccy,    /*	0: domestic, 1: other */
    double  fund_not,    /*	If different from domestic or foreign (fund_ccy = 1) */
    char*   fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 1) */
    double  fx_fund_dom, /*	If different from domestic or foreign (fund_ccy = 1) 2 bd fwd */
    long    fx_fund_dom_spot_date,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    double* fund_fix_cpn, /*	Past coupon fixing if relevant,
                                          includes spr, but not mrg, cvg and notional */
    /*		cf */
    double  cf_not,
    int     cf_ncpn,
    long*   cf_fix,
    long*   cf_start,
    long*   cf_pay,
    char**  cf_basis,
    char**  cf_cms_tenor1,
    char*   cf_cms_freq1,
    char*   cf_cms_basis1,
    double* cf_cms_spread1,
    char**  cf_cms_tenor2,
    char*   cf_cms_freq2,
    char*   cf_cms_basis2,
    double* cf_cms_spread2,
    long    spread_vol_n,
    double* spread_vol_time,
    double* spread_vol_floor,
    double* spread_vol_cap,

    double* spread_slbeta1, /* Shifted log beta on the CMS1 */
    double* spread_slbeta2, /* Shifted log beta on the CMS2 */

    double* cf_alpha,
    double* cf_beta,
    double* cf_gamma,
    int*    cf_floored,
    double* cf_floor,
    int*    cf_capped,
    double* cf_cap,
    double* cf_fix_cpn, /*	Past coupon fixing if relevant */

    int      cf_nopt,      /* Number of spread options */
    double** cf_notopt,    /* Notional of the spread options */
    double** cf_strikeopt, /* spread option strikes */
    int**    cf_typeopt,   /* spread option type 0 call 1 put */

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,
    /*	Numerical params */
    int  req_stp,
    int  req_stpx,
    long req_paths,
    /*	Calib params */
    int    force_atm,
    double max_std_long,
    double max_std_short,
    int    fix_lambda,         /*	0: calib lambda to cap, 1: fix lambda calib
                                                               to diagonal */
    int    cal_vol_shift_type, /*	vol shift type for volatility */
    double cal_vol_shift,      /*	vol shift */
    double cal_lambda_shift,   /*	shift on lambda after calibration */
    int    one_f_equi,         /*	1F equivalent flag:
                                                               if set to 1, then 2F lambda will calibrate
                                                               to the cap priced within calibrated 1F
                                                               with the given lambda */
    int skip_last, /*	If 1, the last option is disregarded and the forward volatility is flat from
                      option n-1 */
    double long_prec,  /*	Precision on primary instruments */
    double short_prec, /*	Precision on secondary instruments */
    double min_fact,   /*	Maximum down jump on variance */
    double max_fact,   /*	Maximum up jump on variance */
    int    use_jumps,  /*	Allow vol term structure to jump */
    int    proba_weight,

    int calc_fwdiv, /*	output the model and market fwdiv  */
    int adjust_fee, /*	adjust the fee to price correctly the fwd iv */

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_pay_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */
    /*	CMS can be valued with smile */
    int cms_adj,                /*	1: adjust for CMS effect
                                                0: don't */
    int cms_for_iv_only,        /*	1: adjust for CMS effect only IV
                                                        0: adjust both IV and call */
    int use_cms_smile,          /*	1: value CMS with ATM
                                                        0: value CMS with smile */
    int cms_vol_adj,            /*	1: adjust for CMS vol effect
                                                        0: don't */
    double           cms_beta1, /*	How adjustment varies with rates */
    double           cms_beta2,
    int              num_strikes_in_vol, /*	Array of strikes in vol matrix */
    double*          strikes_in_vol,
    SrtDiffusionType vol_type, /*	Type of vol in matrix, SRT_NORMAL or SRT_LOGNORMAL */
    int              cash_vol, /*	1: matrix is a cash vol
                                               0: matrix is a swap vol */
    int cvg_sv,                /*	1: coverging model on spread vol/correlation
                                               0: sliding model on spread vol/correlation */
    int is_corr,               /*	1: spread vols are correlations
                                               0: spread vols are vols */

    int use_SL,          /*  1: use Shifted Log
                                     0: don't use				*/
    int calib_SL,        /*  1: Calibrate the Shifted Log models on CMS1 and on CMS2 *
                                         0: use the Shifted Log beta given */
    int calib_correl_SL, /*	1: Calibrate the correlation between the two SL to get the same ATM
                                                         normal spread vol
                                                 0: use the normal spread correl for the sl correl
                          */

    int use_cfoptions, /*  1: Use the spread options and don't take into account
                                                       the floor and cap in the cf coupons
                                               0: take into account the floor and cap in the cf
                          coupons */

    /*	Exercised flag */
    int    exercised,   /*	Flag */
    long   ex_date_ex,  /*	Date when exercised */
    long   ex_date_set, /*	Corresponding settlement date */
    double ex_fee,      /*	Corresponding fee */
    double dLambdaShift,
    /*	Results */
    double* fund_val,  /*	Value of the funding leg */
    double* cf_val,    /*	Value of the Power Dual leg */
    double* call_val,  /*	Value of the callable feature */
    int     export_ts, /*	1: Export TS, 0: don't */
    CCF_UND und_exp);

#endif