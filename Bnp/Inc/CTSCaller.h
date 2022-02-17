
#ifndef __CTS_CALLER_H
#define __CTS_CALLER_H

/* Set the default values for all the CTS autocal parameters */
Err cts_set_default_params(
    int *accrue_on_barrier, int *value_zero, int *use_cmsopt,
    double *correl_start, double *correl_end, int *float_adjust_type,
    int *iv_call_sell, double *iv_call_spread, int *call_call_sell,
    double *call_call_spread, double *numer_call_spread,

    int *iv_trim_type, int *iv_max_fix, double *iv_min_fix_time,
    int *call_trim_type, int *call_max_fix, double *call_min_fix_time,

    int *iNbX, double *iNbSigmaXLeft, double *iNbSigmaXRight,
    double *dIntegParam, int *iIntegMethod, double *dVolLimit, int *iCalibLGM,
    double *dMinStd, double *dMaxStd, double *numer_tstar, double *precision,
    int *nbIterMax, int *keep_first,

    int *pde_or_mc, int *nstpt, int *nstpx, int *nstpvol, int *nstpphi,
    double *mc_mintime, long *npaths, double *integ_mintime,

    int *nb_factor, double *lgm_alpha, double *lgm_gamma, double *lgm_rho,
    int *calib_strategy, int *fix_lambda, int *fix_smile,
    int *smile_calib_months, LGMSV_CalibParams *lgmsv_calib_params,
    int *force_atm, int *strike_flag_long, int *strike_flag_short,
    double *max_std_long, double *max_std_short, double *vol_shift_long,
    DIAGCALIB_VOLTYPE *vol_type_long, DIAGCALIB_SHIFTTYPE *vol_shift_type_long,
    double *vol_shift_short, DIAGCALIB_VOLTYPE *vol_type_short,
    DIAGCALIB_SHIFTTYPE *vol_shift_type_short, double *lambda_shift,
    double *min_time, int *skip_last, double *min_fact, double *max_fact,
    int *use_jumps, int ref_months, char *short_tenor,

    int *calc_fwdiv, int *adjust_fee,

    int *do_one_time, int *one_time_index, int *compute_reserve,
    int *reserve_method, int *lgm_reserve, int *lgm_nstpt, int *lgm_nstpx,
    int *midat_reserve, double *one_time_vega, double *lambda_reserve,
    int *recalib_european, int *recalc_one_factor, int *euro_nb_iter,

    double *tstar,

    int *compatibility_flag);

/*	Caller for callable time swaps */
/*	------------------------------- */

Err cts_caller(
    /*	Market */
    long today, char *yc, /*	yc */
    char *vc,             /*	vc */
    char *ref,            /*	ref rate (only if calib) */
    char *swap_freq,      /*	swap freq (only if calib) */
    char *swap_basis,     /*	swap basis (only if calib) */
    char *lib_freq,       /*	libor freq (only if calib) */
    char *lib_basis,      /*	libor basis (only if calib) */
    Err (*get_cash_vol)(  /*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	The underlying */
    double tstar, int use_calib, /*	0: use lgmsvund  , 1: calibrate */
    /*		if calib */
    int nb_factor, double lambda, /*	LGM lambda */
    double lgm_alpha, double lgm_gamma, double lgm_rho, int nsmilepar,
    double *smilepartime, double *alphaeps, /*	alpha */
    double *rhoeps,                         /*	rho */
    double *ldaeps,                         /*	ldaeps */
    double *rho2eps,                        /*	rho2eps */
    /*	End of calib params */
    char *lgmsvund,
    /*	The structure */
    long start_date, /*	Date at which initial notional exchange occurs */
    /*		funding */
    int fund_ccy, /*	0: domestic  , 1: other */
    double *fund_not,
    char *
        fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 1) */
    double fx_fund_dom, /*	If different from domestic or foreign (fund_ccy
                           = 1) 2 bd fwd */
    long fx_fund_dom_spot_date, int fund_ncpn, long *fund_fix, long *fund_start,
    long *fund_end, long *fund_pay, char **fund_basis, double *fund_spr,
    double *fund_mrg,
    double *fund_fix_cpn, /*	Past coupon fixing if relevant  ,
                                                  includes spr  , but not mrg  ,
                             cvg and notional */
    double *cpn_not,
    /*			number of coupons */
    int ncpn,
    /*			coupon description */
    /*				cpn */
    int *cpn_type, /*	0:	fixed
                                   1:	libor fixed at start
                                   2:	libor fixed at end
                                   3:	CIF coupon */
    long *cpn_start_date, long *cpn_end_date, long *cpn_pay_date,
    double *cpn_coupon, char **cpn_basis,
    /*				pay libor */
    long *pay_fix_date, /*	fixing date of the Libor */
    int pay_months, char *pay_freq, char *pay_basis, double *pay_fwd_spread,
    double *pay_gearing,
    /*				fix libor */
    /*			general Libor properties */
    int fix_lag_bd,
    /*				fixings */
    int *nfix, double **weights, long **fix_dates, char **fix_tenor,
    char *fix_freq, char *fix_basis, double **fix_fwd_spreads,
    /*				profiles */
    int tstype, /*	0: generic  , 1: range */
    /*				profile 0 */
    double *alpha, double *beta, int *nstr, double **str, double **nbopt,
    /*				profile 1 */
    int iv_buy_sell,   /*	1: BNPP buys  , -1: BNPP sells */
    int call_buy_sell, /*	1: BNPP buys  , -1: BNPP sells */
    double iv_call_spread, double call_call_spread, double numer_call_spread,
    int accrue_on_barrier, /*	0: on barrier no accrue  , 1: accrue on barrier
                            */
    int value_zero,        /*	0: Do not value low strike options  , 1: do */
    double *lb,            /*	Lower bounds of ranges */
    double *ub,            /*	Upper bounds of ranges */
    double *payoff,        /*	Payoff in ranges */

    /*				fixings	*/
    double **fixed_ref, /*	Past coupon reference fixing if relevant */
    double *fixed_pay,  /*	Past coupon payment libor fixing if relevant */

    /*	Calls */
    int ncall, int pay_rec, /*	1: rec pd  , -1: pay pd */
    long *ex_date, long *set_date, double *fee,

    /*		Trimming */
    int iv_trim_type, /*	0: no trim
                                      1: x fixings max
                                      2: x time min between two fixings */
    int iv_max_fix, double iv_min_fix_time,
    int call_trim_type, /*	0: no trim
                                1: x fixings max
                                2: x time min between two fixings */
    int call_max_fix, double call_min_fix_time,
    /*	Extra model parameters*/
    int use_cmsopt, /*	Use CmsOption to value fix option  , use BS on CmsRate
                       otherwise */
    double correl_start, /*	Correl between libor fixing and libor paid at
                            start */
    double
        correl_end, /*	Correl between libor fixing and libor paid at start */
    int float_adjust_type, /*	type of adjustment for the floating coupon  ,
                                                   0: ATM vol  , 1: Strike Vol
                            */
    /*	Numerical params */
    /*		CF */
    int iNbX, double iNbSigmaXGridLeft, double iNbSigmaXGridRight,
    double dIntegParam, int iIntegMethod, double dVolLimit, int iCalibLGM,
    double dMinStd, double dMaxStd, double numer_tstar,
    /*		ADI */
    int pde_or_mc, int req_stp, int req_stppsi, int req_stpx, int req_stpz,
    double req_mintime, long req_paths, double integ_mintime,

    /*		Calib */
    char *cal_tenor, char *cal_ref, char *cal_freq, char *cal_basis,
    int force_atm, /*	force atm calib */
    double max_std_long, double max_std_short, double vol_shift_long,
    DIAGCALIB_VOLTYPE vol_type_long, DIAGCALIB_SHIFTTYPE vol_shift_type_long,
    double vol_shift_short, DIAGCALIB_VOLTYPE vol_type_short,
    DIAGCALIB_SHIFTTYPE vol_shift_type_short, double lambda_shift,
    int calib_strategy, /*	0: autocal  , 1: swaptions / cap  , 2: cap /
                           swaptions */
    int fix_lambda,     /*	0: calib lambda to cap  , 1: fix lambda calib
                                        to diagonal */
    char *short_tenor, char *short_refrate, char *short_freq, char *short_basis,
    int fix_smile,          /*	0: calib smile parameters to market smile */
    int smile_calib_months, /* 0: co-terminal swaption  , otherwise underlyings
                               with required nb months */
    LGMSV_CalibParams *lgmsv_calib_params, double min_time,
    int skip_last,   /*	If 1  , the last option is disregarded and the forward
                        volatility is flat from option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps,   /*	1: we allow jumps on vol  , 0: we don't */
    double prec, int maxiter, int keep_first,
    /*	Strike choice */
    int long_strike_flag,  /*	0: ATM
                                                   1: Coupon
                                                   2: Eq (PV/Lvl) */
    int short_strike_flag, /*	0: ATM  ,
                                                   1: implied digital caplet
                              strike 2: same number of std */
    /*		IV calculation */
    int calc_fwd_iv, int adj_fee,

    /*	Flag for extra calculation / adjustment of one-time callable */
    int do_one_time,       /*	1: calc the one time */
    int one_time_index,    /*	0: choose automatically the index  , >0: index
                              provided by user */
    int compute_reserve,   /*	1: adjust the price from the reserve */
    int reserve_method,    /*	0: old method  , 1: new method */
    int lgm_reserve,       /*	1: calculate the reserve in the one factor */
    int lgm_nstpt,         /*	number of time steps for LGM midat */
    int lgm_nstpx,         /*	number of space steps for LGM midat */
    int midat_reserve,     /*	1: reserve calculated using the
                              midat-replication */
    double one_time_vega,  /*	Vega to be applied one the one time callable if
                              required */
    double lambda_reserve, /*	Percent of the switch to be reserved */
    int recalib_european,  /*	If set to 1  , then european are recalibrated
                              when changing lambda */
    int recalc_one_factor, /*	Recalculate MCEB in one factor for reserve */
    int euro_nb_iter,      /*	Number of Levenberg iteration for calibration of
                              european */

    /*	vol Matrix parameters */
    int num_strikes_in_vol, /*	Array of strikes in vol matrix */
    double *strikes_in_vol,
    SrtDiffusionType
        vol_type, /*	Type of vol in matrix  , SRT_NORMAL or SRT_LOGNORMAL */
    int cash_vol, /*	1: matrix is a cash vol
                                  0: matrix is a swap vol */
    /*	Flags */
    /*		EOD */
    int eod_fix_flag, /*	0: I  , 1: E */
    int eod_pay_flag, /*	0: I  , 1: E */
    int eod_ex_flag,  /*	0: I  , 1: E */
    /*	Exercised flag */
    int exercised,    /*	Flag */
    long ex_date_ex,  /*	Date when exercised */
    long ex_date_set, /*	Corresponding settlement date */
    double ex_fee,    /*	Corresponding fee */

    /* compatibility flag for old version of Callable TS */
    int compatibility_flag,

    /*	Results */
    double *fund_val,        /*	Value of the funding leg */
    double *cpn_val,         /*	Value of the Power Dual leg */
    double *call_val,        /*	Value of the multi-callable feature */
    double *onetime_val,     /*	Value of the one-time callable feature */
    double *onetime_reserve, /*	One Time Callable Reserve */
    double *switch_reserve,  /*	Switch reserve */
    /*	Feedback */
    int export_und, CTS_UND und_exp, int save_inst_data,
    cpd_calib_inst_data *inst_data, int save_fwdiv, CTS_IV fwd_iv_info,
    int save_extra_infos, CTS_EXTRA_INFOS extra_infos);

#endif