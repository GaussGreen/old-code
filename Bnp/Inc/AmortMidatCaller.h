
#ifndef __AMORT_MIDAT_CALLER_H
#define __AMORT_MIDAT_CALLER_H

#include "amortmidatprodstruct.h"

/*	------------------------------- */
/*	---------Amortized Midat------- */
/*	------------------------------- */

Err am_caller(
    /*	Today's date */
    long today,

    /*	The underlying */
    int use_calib, /*	0: use lgm2fund  , 1: calibrate */

    /*		if calib */
    char *yc,         /*	yc */
    char *vc,         /*	vc */
    char *ref,        /*	ref rate (only if calib) */
    char *swap_freq,  /*	swap freq (only if calib) */
    char *swap_basis, /*	swap basis (only if calib) */
    double lambda,    /*	lambda if unique */
    double alpha,     /*	alpha */
    double gamma,     /*	gamma */
    double rho,       /*	rho */

    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),

    /*		if no calilb */
    char *lgm2dund,

    /*	The structure */
    long start_date,  /*	Date at which initial notional exchange occurs */
    long theoEndDate, /* End Date of the structure */

    /*		funding */
    char *fund_ref, int fund_ccy, /*	0: domestic  , 1: other */
    double
        *fund_not, /*	If different from domestic or foreign (fund_ccy = 1) */
    char *
        fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 1) */
    double fx_fund_dom, /*	If different from domestic or foreign (fund_ccy
                           = 1) 2 bd fwd */
    long fx_fund_dom_spot_date, int fund_ncpn, long *fund_fix, long *fund_start,
    long *fund_end, long *fund_pay, char **fund_basis, double *fund_spr,
    double *fund_mrg,
    double *
        fund_fix_cpn, /*	Past coupon fixing if relevant  ,
                              includes spr  , but not mrg  , cvg and notional */
    /*		fix */
    double *fix_not, int fix_ncpn,
    //			long		*fix_fix  ,
    long *fix_start, long *fix_end, long *fix_pay, char **fix_basis,
    double *fix_rate, double *fix_fee, /*	Exercise Fee */
    //			double		*fix_fix_cpn  ,			/*	Past coupon fixing if relevant
    //*/

    /*		calls */
    int ncall, int pay_rec, /*	0: rec pd  , 1: pay pd */
    long *ex_date, long *set_date, double *fee,

    /*	Numerical params */
    int req_stp, int req_stpx,

    /*	Calib params */
    double mintime, double mininterval, int notperiod, int one2F, int use_jump,
    double max_var_jump, int strike_type, int european_model,

    /*	Function to get correlation stored from a vol cube */
    Err (*get_correl)(char *correl_cube_name, double start_date,
                      double end_date, double strike, double *vol),
    char *CorrelName,
    /*
                            Err (*GetVolForBadr)( Date  , Date  , double  ,
       SRT_Boolean  , double *)  , char *cVolType  ,
    */
    double max_std_short, int fix_lambda, /*	0: calib lambda to cap  , 1: fix
                                             lambda calib to diagonal */
    int one_f_equi,                       /*	1F equivalent flag:
                                                                          if set to 1  , then 2F
                                             lambda will calibrate                       to the cap priced within calibrated
                                             1F                       with the given lambda */
    int skip_last, /*	If 1  , the last option is disregarded
                                                   and the forward volatility is
                      flat from option n-1 */

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I  , 1: E */
    int eod_pay_flag, /*	0: I  , 1: E */
    int eod_ex_flag,  /*	0: I  , 1: E */

    /*	Exercised flag */
    int exercised,    /*	Flag */
    long ex_date_ex,  /*	Date when exercised */
    long ex_date_set, /*	Corresponding settlement date */
    double ex_fee,    /*	Corresponding fee */

    /*	Results */
    double *fund_val, /*	Value of the funding leg */
    double *fix_val,  /*	Value of the Power Dual leg */
    double *call_val, /*	Value of the callable feature */
    int export_ts,    /*	1: Export TS  , 0: don't */
    AM_UND und_exp);

Err am_calc_mdl_iv_fwd(AM_STR am, AM_UND und, AM_ADI_ARG adi_arg,

                       int num_hermite,
                       /*	Result */
                       double *premium);

Err am_caller_ts(
    /*	Today's date */
    long today,

    /*	The underlying */
    int use_calib, /*	0: use lgm2fund  , 1: calibrate */

    /*		if calib */
    char *yc,              /*	yc */
    char *vc,              /*	vc */
    char *ref,             /*	ref rate (only if calib) */
    char *swap_freq,       /*	swap freq (only if calib) */
    char *swap_basis,      /*	swap basis (only if calib) */
    double *pdLambdaValue, /*	lambda if unique */
    double *pdLambTime, int nLambdaSize, double alpha, /*	alpha */
    double gamma,                                      /*	gamma */
    double rho,                                        /*	rho */

    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),

    /*		if no calilb */
    char *lgm2dund,

    /*	The structure */
    long start_date,  /*	Date at which initial notional exchange occurs */
    long theoEndDate, /* End Date of the structure */

    /*		funding */
    char *fund_ref, int fund_ccy, /*	0: domestic  , 1: other */
    double
        *fund_not, /*	If different from domestic or foreign (fund_ccy = 1) */
    char *
        fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 1) */
    double fx_fund_dom, /*	If different from domestic or foreign (fund_ccy
                           = 1) 2 bd fwd */
    long fx_fund_dom_spot_date, int fund_ncpn, long *fund_fix, long *fund_start,
    long *fund_end, long *fund_pay, char **fund_basis, double *fund_spr,
    double *fund_mrg,
    double *
        fund_fix_cpn, /*	Past coupon fixing if relevant  ,
                              includes spr  , but not mrg  , cvg and notional */
    /*		fix */
    double *fix_not, int fix_ncpn,
    //			long		*fix_fix  ,
    long *fix_start, long *fix_end, long *fix_pay, char **fix_basis,
    double *fix_rate, double *fix_fee, /*	Exercise Fee */
    //			double		*fix_fix_cpn  ,			/*	Past coupon fixing if relevant
    //*/

    /*		calls */
    int ncall, int pay_rec, /*	0: rec pd  , 1: pay pd */
    long *ex_date, long *set_date, double *fee,

    /*	Numerical params */
    int req_stp, int req_stpx,

    /*	Calib params */
    double mintime, double mininterval, int notperiod, int one2F, int use_jump,
    double max_var_jump, int strike_type, int european_model,

    /*	Function to get correlation stored from a vol cube */
    Err (*get_correl)(char *correl_cube_name, double start_date,
                      double end_date, double strike, double *vol),
    char *CorrelName,
    /*
                            Err (*GetVolForBadr)( Date  , Date  , double  ,
       SRT_Boolean  , double *)  , char *cVolType  ,
    */
    double max_std_short, int fix_lambda, /*	0: calib lambda to cap  , 1: fix
                                             lambda calib to diagonal */
    int one_f_equi,                       /*	1F equivalent flag:
                                                                          if set to 1  , then 2F
                                             lambda will calibrate                       to the cap priced within calibrated
                                             1F                       with the given lambda */
    int skip_last, /*	If 1  , the last option is disregarded
                                                   and the forward volatility is
                      flat from option n-1 */

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I  , 1: E */
    int eod_pay_flag, /*	0: I  , 1: E */
    int eod_ex_flag,  /*	0: I  , 1: E */

    /*	Exercised flag */
    int exercised,    /*	Flag */
    long ex_date_ex,  /*	Date when exercised */
    long ex_date_set, /*	Corresponding settlement date */
    double ex_fee,    /*	Corresponding fee */

    /*	Results */
    double *fund_val, /*	Value of the funding leg */
    double *fix_val,  /*	Value of the Power Dual leg */
    double *call_val, /*	Value of the callable feature */
    int export_ts,    /*	1: Export TS  , 0: don't */
    AM_UND und_exp);

#endif