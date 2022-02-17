
#ifndef __CPD_CALLER_H
#define __CPD_CALLER_H

/*
Convert the funding leg into a domestic one
                -	spreads and margins
                -	past fixings
                -	final notional exchange
                -	initial notional exchange
*/
Err convert_funding_to_domestic(
    /*	Inputs */
    long today,                 /*	Today */
    long not_ex_date,           /*	Date at which the
                                                                    initial notional
                                                                    exchange takes
                               place           (or has taken place) */
    int eod_fix_flag,           /*	0: I      , 1: E */
    int eod_pay_flag,           /*	0: I      , 1: E */
    double fx_fund_dom,         /*	Fx fund/dom      , 2bd fwd */
    long fx_fund_dom_spot_date, /*	Spot date for Fx */
    double dom_not,             /*	Domestic notional */
    char *dom_yc,               /*	Domestic discount curve */
    int fund_ncpn,              /*	Number of coupon */
    long *fund_fix,             /*	Fixing dates */
    long *fund_start,           /*	Start dates */
    long *fund_pay,             /*	Pay dates */
    char **fund_basis,          /*	Basis */
    /*	The following are modified */
    char *fund_yc,        /*	Funding discount curve      ,
                                              changed to domestic */
    double *fund_not,     /*	Funding notional
                                                      in funding ccy      ,
                                              converted to domestic */
    double *fund_spr,     /*	Spread in funding ccy      ,
                                              put to 0 */
    double *fund_mrg,     /*	Margin in funding ccy      ,
                                              converted to margin over
                                              cash libor in domestic currency */
    double *fund_fix_cpn, /*	Fixing: contains spread
                                                      but not margin      ,
                                              converted to equivalent
                                              domestic cash-flows */
    /*	The following are returned */
    long *fund_start_date, /*	Start date of the funding */
    double *eq_final_ex,   /*	Domestic cash-flow equivalent
                                                               to final
                          exchange   (to be delivered   at funding start date) */
    double
        *eq_init_ex); /*	Domestic cash-flow equivalent
                                                              to initial
                         exchange (to be delivered at initial exchange date) */

Err transform_interp_coupon(int npts, double *pts, double *cpn, int linxtr_l,
                            int linxtr_r, double *cst, double *wspot,
                            int *first_pt, int *last_pt, double *weights);

void transform_oldspec_into_newspec(double alpha, double beta, int floored,
                                    double floor, int capped, double cap,
                                    int *nstrikes, double *wcst, double *wspot,
                                    double *strikes, double *weights);

/*	Caller for callable power duals */
/*	------------------------------- */

Err cpd_caller(
    /*	Today's date */
    long today,
    /*	The underlying */
    int use_calib, /*	0: use fx3dund      , 1: calibrate */
    /*		if calib */
    double fx_spot,                   /*	2bd fwd */
    long fx_spot_date, int dom_calib, /*	Calibrate domestic underlying */
    char *dom_und,        /*	If no      , domestic underlying to be used */
    char *dom_yc,         /*	Domestic yc */
    char *dom_vc,         /*	Domestic vc (only if calib) */
    char *dom_ref,        /*	Domestic ref rate (only if calib) */
    char *dom_swap_freq,  /*	Domestic swap freq (only if calib) */
    char *dom_swap_basis, /*	Domestic swap basis (only if calib) */
    double dom_lam,       /*	Domestic lambda */
    int for_calib,        /*	Same for foreign */
    char *for_und, char *for_yc, char *for_vc, char *for_ref,
    char *for_swap_freq, char *for_swap_basis, double for_lam,

    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps,   /*	Allow vol term structure to jump */

    double *corr_times, double *correl_dom_for, /*	Correlations */
    double *correl_dom_fx, double *correl_for_fx, long corr_n_times,
    CPDBETADLMPARAMS cpd_dlm_params,
    Err (*get_ir_cash_vol)(/*	Function to get IR cash vol from the markets */
                           char *vol_curve_name, double start_date,
                           double end_date, double cash_strike, int zero,
                           char *ref_rate_name, double *vol, double *power),
    /*	Fx vol from the market */
    long *fx_mkt_vol_date, double *fx_mkt_vol, int num_fx_mkt_vol,
    /*	Fx SABR parameters from the market */
    double *fx_mkt_smile_alpha, double *fx_mkt_smile_beta,
    double *fx_mkt_smile_rho, double *fx_mkt_smile_pi,
    int use_sabr,      /*	0: no smile adj-t      , 1: only for underlying      ,
                      2: fees      adj-t for the call/KO */
    int use_3F_interp, /*	0: interpolates the ATM market vols linearly ,
                      1: uses the 3F for the interpolation */
    int smile_spec_type, //	0: lognormal vol + SABR params      , 1:
                         //sigma-beta
                         //+ SABR params      , 2: BMM (not yet supported)
    /*		if no calilb */
    char *fx3dund,
    /*	The structure */
    long start_date, /*	Date at which initial notional exchange occurs */
    /*		funding */
    double fund_not, /*	Notional */
    int fund_ccy,    /*	0: domestic      , 1: foreign 2: other */
    char *
        fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 2) */
    double fx_fund_dom, /*	If different from domestic or foreign (fund_ccy
                       = 2) 2 bd fwd */
    long fx_fund_dom_spot_date, int fund_ncpn, /*	Number of coupons */
    long *fund_fix,                            /*	Fixing dates */
    long *fund_start,                          /*	Start dates */
    long *fund_pay,                            /*	Pay dates */
    char **fund_basis,                         /*	Basis */
    double *fund_spr,                          /*	Forward spreads */
    double *fund_mrg,                          /*	Margins */
    double *fund_fix_cpn, /*	Past coupon fixing if relevant      ,
                                                      includes spr      , but
                         not mrg      , cvg and notional */
    /*		pd */
    double pd_not,   /*	Notional */
    int pd_ncpn,     /*	Number of coupons */
    long *pd_fix,    /*	Fx fixing dates */
    long *pd_start,  /*	Start dates */
    long *pd_pay,    /*	Pay dates */
    char **pd_basis, /*	Basis */
    double
        *pd_alpha, /*	Coupon = alpha + beta * fx [capped      , floored] */
    double *pd_beta, int *pd_floored, double *pd_floor, int *pd_capped,
    double *pd_cap,
    /*		pd interp coupon specification */
    int *pd_nfxpts,         /*	Number of coupon interpolation points */
    double **pd_fxpts,      /*	Coupon interpolation points */
    double **pd_cpn_at_pts, /*	Coupon values at interpolation points
                           (interpolation is linear) */
    int *pd_lin_xtrpl_l,    /*	0 = flat extrapolation to the left      , 1 =
                           linear    w first 2 pts slope */
    int *pd_lin_xtrpl_r,    /*	0 = flat extrapolation to the right      , 1 =
                           linear w last 2 pts slope */
    double *pd_fix_fx,      /*	Past Fx fixing if relevant */

    /*		pd not refund */
    long *pd_not_ref_fix,    //	fx fixing dates FOR EACH CALL DATE + in the end
                             // if relevant
    double pd_not_ref_alpha, /*	Final notional on PD leg */
    double pd_not_ref_beta, int pd_not_ref_floored, double pd_not_ref_floor,
    int pd_not_ref_capped, double pd_not_ref_cap,
    /*		pd not interp specification */
    int *pd_not_ref_nfxpts, /*	Number of notional interpolation points FOR EACH
                           CALL DATE + in the end */
    double **pd_not_ref_fxpts, double **pd_not_ref_cpn_at_pts,
    int *pd_not_ref_lin_xtrpl_l, /*	0 = flat extrapolation to the left , 1
                                = linear w first 2 pts slope */
    int *pd_not_ref_lin_xtrpl_r, /*	0 = flat extrapolation to the right , 1
                                = linear w last 2 pts slope */
    double *pd_not_ref_fix_fx,   //	fx fixings FOR EACH CALL DATE + in the
                                 // end if relevant
    /*		calls */
    int *call_type,    /*	0: call      , 1: KO */
    int ncall,         /*	Number of calls */
    int pay_rec,       /*	0: rec pd      , 1: pay pd */
    long *ex_date,     /*	Call dates */
    long *set_date,    /*	Settlement dates */
    double *barrier,   /*	in case of a pure KO or a Callable KO */
    int *bar_type,     /*	0: up and in      , 1: down and in */
    double *fees,      /*  fees if deal is called in domestic currency */
    int TARN_Do,       /*	Is it a Powerdual TARN */
    double TARN_Floor, /*  Floor of the TARN level */
    /*	Numerical params */
    long req_stp,      /*	Number of time steps in the tree */
    long req_pth,      /*	Number of paths in the MC */
    double bar_smooth, /*	Smoothing factor for barriers */
    int do_pecs,       /*	Do PECS in the MC */
    int forcetree,     /*	If equal to 1 then the valuation is done in a tree
                        */
    int do_optim,      /*	If equal to 1 then the call are replaced by optimal KO
                        */
    int force_optim, /*	If equal to 1 then all call will be replaced by optimal
                    KO	*/
    int fx_bound,    /*	If equal to 1 then optimisation on the Fx      , on the
                    IV    otherwise	*/
    int use_bound,   /*	If equal to 1 then prices the call as UO on the Fx using
                    a provided boundary	*/
    int do_infos,    /*	infos on callable right */
    /*	EOD Fixing Flag */
    int eod_fix_flag, /*	0: I      , 1: E */
    /*	EOD Payment Flag */
    int eod_pay_flag, /*	0: I      , 1: E */
    /*	EOD Exercise Flag */
    int eod_ex_flag, /*	0: I      , 1: E */
    /*	Vega */
    double dom_vol_shift, double for_vol_shift, double fx_vol_shift,
    /*	Exercised flag */
    int exercised,    /*	Flag */
    long ex_date_ex,  /*	Date when exercised */
    long ex_date_set, /*	Corresponding settlement date */
    /*	Prune calls */
    int prune_calls,    /*	Flag: 0 no prune      , 1 do prune */
    int no_prune_years, /*	Number of years from the next call to the first
                       prune */
    int prune_factor,   /*	Ex: 2 -> one call out of two */

    /* Disable call dates */
    int disable_calls_method, /* Name of the disable call method used : 0 None
                                 , 1 All in range 2 list of call ( Parameter1
                             and Parameter 2 not used */
    int first_call_index,     /* Parameter 1 of the disable method */
    int last_call_index,      /* Parameter 2 of the disable method */

    int iNbIndex_List, /* Number of erased call in the list */
    int *Index_List,   /* Array of index of calls to be erased */

    int erasing_call_done, // always to 0      , switch to 1 after haveng
                           // disabled calls
    int *erased_call_list, // list of 0 & 1 used only if erasing_call_done = 1
                           // and if use_GMA = TRUE

    /* Alpha Beta Smile Model */
    int use_smile, /* Flag: 0 no use alphabeta model      , 1 use of the
                  alphabeta model */
    double alpha,  /* Alpha      , only used if use_smile =1 */
    double beta,   /* Beta      ,  only used if use_smile =1 */
    /*For the smile impact*/
    int use_GMA,
    int use_3f_optim_barrier, // Optimize the boundary with the 3F and then
                              // price the smile adjustment
    double SSLstd, int SSLniter, int CopNsimul, int CopDoPecs,
    int CummulNpoints, int CummulNstd, double CummulPrecision, int CummulLinear,
    int PayoffFunction, double fwdSmileVisuNstd, int smileOtc, int smileFee,
    int fwdVolMethod, int smileModel, double BMpi, int otc, int FundingSpeedUp,
    /* Do not calculate all OTC*/
    int nStart, int oneOutOfN,
    /*Fast MC*/
    int FMC_do, double FMC_precision, int FMC_min_paths,
    /*	Results */
    double *fund_val,         /*	Value of the funding leg */
    double *pd_val,           /*	Value of the Power Dual leg */
    double *call_val,         /*	Value of the callable feature */
    double *call_stdev,       /*	Standard deviation of the call if applicable */
    double *smile_adjustment, /*	GMA smile adjustment */
    double ***optim_bar,      /*	Contains the value of the optimal KO for
                             corresponding calls */
    double **GMA_Results,     /*	Contains the OTC/KO info */
    int export_ts,            /*	1: Export TS      , 0: don't */
    CPD_UND und_exp);         /*	TS to be exported */

#endif