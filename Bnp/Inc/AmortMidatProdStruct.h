
#ifndef __AMORTMIDAT_PROD_STRUCT_H
#define __AMORTMIDAT_PROD_STRUCT_H

#include "cpdcalib.h"
#include "utdates.h"

#define AM_NCPN 512

/*	Structures and functions for callable cms floaters */
/*	-------------------------------------------------- */

/*	Structures and functions for the funding leg */

/*	Funding cpn */
typedef struct
{
    long   start_date; /*	Coupon start date */
    double start_time; /*	Coupon start time */
    long   end_date;   /*	Coupon end date */
    double end_time;   /*	Coupon end time */
    long   pay_date;   /*	Coupon pay date */
    double pay_time;   /*	Coupon pay time */
    double cvg;        /*	Cvg */
    double cpn;        /*	Notional * (fwd_spr + margin) * cvg */
} am_fund_cpn, *AM_FUND_CPN;

/*	Funding leg */
typedef struct
{
    double* notional; /*	Notional */
    double* margin;   /*	Margins */
    double* spread;   /*	Spreads */
    int     num_cpn;  /*	Number of coupons */
    /*	0..num_cpn-1 */
    am_fund_cpn* cpn;
} am_fund_leg, *AM_FUND_LEG;

char* am_fill_fund_leg(
    /*	Coupons that started before today are disregarded */
    long today,
    /*	EOD Flag */
    int         eod_flag, /*	0: I, 1: E */
    double*     fund_not,
    int         fund_ncpn,
    long*       fund_fix,
    long*       fund_start,
    long*       fund_end,
    long*       fund_pay,
    char**      fund_basis,
    double*     fund_spr,
    double*     fund_mrg,
    AM_FUND_LEG fund_leg);

/*	Check dates consistency */
char* am_check_fund_leg(AM_FUND_LEG fund_leg);

/*	Free */
char* am_free_fund_leg(AM_FUND_LEG fund_leg);

/*	Structures and functions for the fixed leg */

/*	CMS Float cpn */
typedef struct
{
    long   start_date; /*	Coupon start date */
    double start_time; /*	Coupon start time */
    long   end_date;   /*	Coupon end date */
    double end_time;   /*	Coupon end time */
    long   pay_date;   /*	Coupon pay date */
    double pay_time;   /*	Coupon pay time */
    double cpn;        /*	Fixed Rate * coverage * notional */
    double cvg;        /*	Fixed Rate * coverage * notional */

} am_fix_cpn, *AM_FIX_CPN;

/*	Fixed leg */
typedef struct
{
    double* notional;
    double* rate;
    double* fee;
    int     num_cpn; /*	Number of coupons */
    /*	0..num_cpn-1 */
    am_fix_cpn* cpn;

} am_fix_leg, *AM_FIX_LEG;

char* am_fill_fix_leg(
    /*	Coupons that fixed before today are disregarded */
    long today,
    /*	EOD Flag */
    int     eod_flag, /*	0: I, 1: E */
    double* fix_not,
    int     fix_ncpn,
    // long		*fix_fix,
    long*      fix_start,
    long*      fix_end,
    long*      fix_pay,
    char**     fix_basis,
    double*    fix_rate,
    double*    fix_fee,
    AM_FIX_LEG fix_leg);

/*	Check dates consistency */
char* am_check_fix_leg(AM_FIX_LEG fix_leg);

/*	Free */
char* am_free_fix_leg(AM_FIX_LEG fix_leg);

/*	Structures and functions for the calls */

/*	Call */
typedef struct
{
    int pay_rec; /*	0: rec fix upon exercise, 1: pay fix upon exercise */
    /*	Specs */
    long   ex_date;      /*	Exercise date */
    double ex_time;      /*	Exercise time */
    int    fund_idx;     /*	Index of the first funding coupon to be called */
    int    num_fund_cpn; /*	Number of funding coupons called (including redemption) */
    int    fix_idx;      /*	Index of the first fix coupon to be called */
    int    num_fix_cpn;  /*	Number of fix coupons called (excluding redemption) */
    /*	Fee upon exercise */
    long   set_date; /*	Fee payment date */
    double set_time; /*	Fee payment time */
    double fee;      /*	Amount of the fee */
} am_call, *AM_CALL;

/*	Amortized Midat */
typedef struct
{
    am_fix_leg*  fix_leg;
    am_fund_leg* fund_leg;
    int          num_calls;
    long         theoEndDate;
    am_call*     call;

} am_str, *AM_STR;

char* am_fill_calls(
    /*	Exercises before today are disregarded */
    long today,
    /*	EOD Flag */
    int     eod_flag, /*	0: I, 1: E */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,
    AM_STR  am);

/*	Check dates consistency */
char* am_check_calls(AM_STR am);

/*	Free */
char* am_free_calls(AM_STR am);

/*	Structures and functions for the underlying and its term structures */

/*	Underlying term structs */
typedef struct
{
    char    name[256];
    long    today;
    char    yc[256];
    char    vc[256];
    char    ref[256]; /*	Reference rate used for getting the vol */
    char    swap_freq[256];
    char    swap_basis[256];
    int     cvg_sv;
    int     sigma_n;
    double* sigma_date;
    double* sigma_time;
    double* sigma; /*	Sigma1 */
    double  lambda;
    double  alpha; /*	Sigma2 = alpha * Sigma1 */
    double  gamma; /*	Lambda2 = Lambda1 + gamma */
    double  rho;

    /*	Calibration instrument data */
    int                 has_inst_data;
    cpd_calib_inst_data inst_data;

    /*	Value of the fwd amortized swap */
    int     num_prices;
    double* exercise_dates;
    double* mkt_prices;
    double* mdl_prices;

    /* added lambda term structure by Albert Wang 12/15/03
    ideally single lambda variable above should be eliminated.
     will come back later .. */
    long    lambda_n;
    double* pdlambda_date;
    double* pdlambda_time;
    double* pdlambda;

} am_und, *AM_UND;

/*	Fill underlying structure from a predefined underlying */
char* am_fill_und(
    char* lgm2dund, char* vc, char* ref, char* swap_freq, char* swap_basis, AM_UND und);

/*	Fill underlying structure from calibration instruments */
char* am_calib_und(
    long today,
    /*	EOD Flag */
    int    eod_flag,   /*	0: I, 1: E */
    char*  yc,         /*	yc */
    char*  vc,         /*	vc (only if calib) */
    char*  ref,        /*	ref rate (only if calib) */
    char*  swap_freq,  /*	swap freq (only if calib) */
    char*  swap_basis, /*	swap basis (only if calib) */
    double lambda,     /*	lambda if unique */
    double alpha,      /*	alpha */
    double gamma,      /*	gamma */
    double rho,        /*	rho */
    /*	Calib params */
    double mintime,
    double mininterval,

    int    notperiod,
    double max_std_short,
    int    one2F,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                       to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */
    double max_var_jump,

    int strike_type,
    int european_model,

    char* (*get_correl)(
        char* correl_cube_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,

    char* (*GetVolForBadr)(Date, Date, double, SRT_Boolean, double*),
    char* cVolType,

    /*	End of calib params */
    AM_STR am,            /*	structure */
    char* (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                          char*   vol_curve_name,
                          double  start_date,
                          double  end_date,
                          double  cash_strike,
                          int     zero,
                          char*   ref_rate_name,
                          double* vol,
                          double* power),
    AM_UND und);

void am_copy_und(AM_UND src, AM_UND dest);

char* am_free_und(AM_UND und);

/*	Constants used for reconstruction and evaluation  */
/*	------------------------------------------------- */

#define MAXDF 1000
typedef struct
{
    /*	Discount Factors */
    /*	Reconstruction is exp ( - alpha - beta * r1 - gamma * r2 ) */

    /*	All df */
    int    num_df;          /*	Number of df required */
    double df_mat[MAXDF];   /*	Residual maturity of df */
    double df_alpha[MAXDF]; /*	For reconstruction */
    double df_beta[MAXDF];
    double df_gamma[MAXDF];
    double df[MAXDF]; /*	Space to store the df */

    /*	Funding */
    int do_fund;
    int num_fund_cpn;
    int start_idx;
    int fund_idx[AM_NCPN];

    /*	Fix discount */
    int do_fix_disc;
    int num_fix_cpn;
    int fix_disc_idx[AM_NCPN];

    /*	Fix coupons */
    int do_fix_fwd;
    int fix_idx[AM_NCPN];

    /*	Fee */
    int fee_idx;

} am_eval_const, *AM_EVAL_CONST;

/*	Fill structure */
char* am_fill_eval_const(
    AM_UND und,
    AM_STR am,
    /*	Index of the current call */
    int           call_idx,
    AM_EVAL_CONST eval_const);

char* am_fill_eval_const_ts(
    int     nLambda,
    double* pdLambdaValue,
    double* pdLambdaTime,
    AM_UND  und,
    AM_STR  am,
    /*	Index of the current call */
    int           call_idx,
    AM_EVAL_CONST eval_const);

/*	Arguments to all payoff evaluation functions */
/*	-------------------------------------------- */

typedef struct
{
    /*	Underlying */
    am_und* und;

    /*	Structure */
    am_str* am;

    /*	Constants used for reconstruction and evaluation */
    am_eval_const eval_const;

    /*	Index of the current call */
    int call_idx;

} am_pay_arg, *AM_PAY_ARG;

/*	Arguments to the adi function */
/*	----------------------------- */

typedef struct
{
    int     nstp;
    double* time;
    double* date;
    int     nstpx;
    double* sig_time;
    double* sig1;
    int     nb_sig;
    double  lam;
    double  alpha;
    double  gamma;
    double  rho;
    void**  void_prm;
    int*    is_event;
    double* ifr;
    char    yc[256];
} am_adi_arg, *AM_ADI_ARG;

char* am_fill_adi_arg(
    AM_UND und,
    AM_STR am,
    char* (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                          char*   vol_curve_name,
                          double  start_date,
                          double  end_date,
                          double  cash_strike,
                          int     zero,
                          char*   ref_rate_name,
                          double* vol,
                          double* power),
    /*	Required number of steps */
    int        req_stp,
    int        req_stpx,
    AM_ADI_ARG am_arg);

char* am_fill_adi_arg_ts(

    int     nLambda,
    double* pdLambdaValue,
    double* pdLambdaTime,

    AM_UND und,
    AM_STR am,
    char* (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                          char*   vol_curve_name,
                          double  start_date,
                          double  end_date,
                          double  cash_strike,
                          int     zero,
                          char*   ref_rate_name,
                          double* vol,
                          double* power),
    /*	Required number of steps */
    int        req_stp,
    int        req_stpx,
    AM_ADI_ARG am_arg);

char* am_free_adi_arg(AM_ADI_ARG am_arg);

/*	Main function to be called in order to fill and check all structures */
/*	==================================================================== */

/*	Fill and check all the relevant structures */
char* am_fill_check_all_struct(
    /*	Today's date */
    long today,
    long theoEndDate,

    /*	The underlying */
    int use_calib, /*	0: use lgm2dund, 1: calibrate */

    /*		if calib */
    char*  yc,         /*	yc */
    char*  vc,         /*	vc (only if calib) */
    char*  ref,        /*	ref rate (only if calib) */
    char*  swap_freq,  /*	swap freq (only if calib) */
    char*  swap_basis, /*	swap basis (only if calib) */
    double lambda,     /*	lambda if unique */
    double alpha,      /*	alpha */
    double gamma,      /*	gamma */
    double rho,        /*	rho */
    /*	End of calib params */

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */

    /*		funding */
    char*   fund_ref,
    double* fund_not,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,

    /*		cf */
    double* fix_not,
    int     fix_ncpn,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee,

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    char* (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                          char*   vol_curve_name,
                          double  start_date,
                          double  end_date,
                          double  cash_strike,
                          int     zero,
                          char*   ref_rate_name,
                          double* vol,
                          double* power),

    double mintime,
    double mininterval,

    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    char* (*get_correl)(
        char* vol_curve_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,

    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                               to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Results */
    AM_STR am,
    AM_UND und,
    int*   call_feat, /*	0: No callable feature to be valued
                              1: Callable feature to be valued through adi */
    AM_ADI_ARG adi_arg);

/*	Fill and check all the relevant structures */
char* am_fill_check_all_struct_ts(
    /*	Today's date */
    long today,
    long theoEndDate,

    /*	The underlying */
    int use_calib, /*	0: use lgm2dund, 1: calibrate */

    /*		if calib */
    char* yc,         /*	yc */
    char* vc,         /*	vc (only if calib) */
    char* ref,        /*	ref rate (only if calib) */
    char* swap_freq,  /*	swap freq (only if calib) */
    char* swap_basis, /*	swap basis (only if calib) */

    int     nlambda,
    double* pdlambda_time,
    double* pdlambda,

    double alpha, /*	alpha */
    double gamma, /*	gamma */
    double rho,   /*	rho */
    /*	End of calib params */

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */

    /*		funding */
    char*   fund_ref,
    double* fund_not,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,

    /*		cf */
    double* fix_not,
    int     fix_ncpn,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee,

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    char* (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                          char*   vol_curve_name,
                          double  start_date,
                          double  end_date,
                          double  cash_strike,
                          int     zero,
                          char*   ref_rate_name,
                          double* vol,
                          double* power),

    double mintime,
    double mininterval,

    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    char* (*get_correl)(
        char* vol_curve_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,

    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                               to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Results */
    AM_STR am,
    AM_UND und,
    int*   call_feat, /*	0: No callable feature to be valued
                              1: Callable feature to be valued through adi */
    AM_ADI_ARG adi_arg);

char* am_free_all_struct(AM_STR am, AM_UND und, int call_feat, AM_ADI_ARG adi_arg);

/*	Payoff function for adi (callable) */
/*	-----------------------------------	*/

char* am_payoff_4_lgm2f_adi(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    void*   yc,
    double* lam,
    double* ts_time,
    int     nb_ts,
    double  gamma,
    double  rho,
    double  phi1,
    double  phi2,
    double  phi12,
    /* Nodes data */
    int      l1,
    int      u1,
    int      l2,
    int      u2,
    double*  r1,
    double** r2,
    int      nprod,
    /* Vector of results to be updated */
    double*** prod_val);

/*	Main pricing function for adi */

/*	Launch the adi */
char* am_launch_adi(
    AM_STR     am,
    AM_UND     und,
    AM_ADI_ARG adi_arg,
    /*	Result */
    double* prem);

char* am_launch_adi_ts(
    double*    pdLambdaValue,
    double*    pdLambdaTime,
    int        nLambdaSize,
    AM_STR     am,
    AM_UND     und,
    AM_ADI_ARG adi_arg,
    /*	Result */
    double* prem);

#endif
