
#ifndef __CPD_PROD_STRUCT_H
#define __CPD_PROD_STRUCT_H

#include "Fx3FBetaDLMCalculations.h"
#include "Fx3FBetaDLMCalibration.h"
#include "Fx3FBetaDLMUtil.h"

/*	Structure for BetaDLM model */
typedef struct
{
    int use_beta_dlm;
    int use_beta_dlm_for_iv;
    int use_tree_iv_for_call;

    /* Model informations */
    double tstar;
    double B0;
    double C0;
    double alpha;
    double lambda;

    /* Smile informations */
    int    calib_smile;
    double calib_smile_std_up;
    double calib_smile_std_down;
    double calib_smile_maturity;
    double calib_smile_shift;
    long   smile_settlmt_date;
    double smile_strikes[2];
    double smile_bs_vols[2];

    /* Numerical parameters */
    FxBetaDLM_OptNumerParams*  NumParams;
    FxBetaDLM_GRFNNumerParams* GrfnParams;

} cpd_beta_dlm_params, *CPDBETADLMPARAMS;

void cpd_init_beta_dlm_params(cpd_beta_dlm_params* cpd_dlm_params);

Err cpd_free_beta_dlm_params(cpd_beta_dlm_params* cpd_dlm_params);

/*	Structures and functions for callable power duals */
/*	------------------------------------------------- */

/*	Structures and functions for the funding leg */

/*	Funding cpn */
typedef struct
{
    long   start_date; /*	Coupon start date */
    double start_time; /*	Coupon start time */
    long   pay_date;   /*	Coupon pay date */
    double pay_time;   /*	Coupon pay time */
    double cpn;        /*	notional * (fwd_spr + margin) * cvg */
} pd_fund_cpn, *PD_FUND_CPN;

/*	Funding leg */
typedef struct
{
    int    dom_for;  /*	0: funding domestic, 1: funding foreign */
    double notional; /*	Notional */
    int    num_cpn;  /*	Number of coupons, not including notional refund */
    /*	0..num_cpn-1 */
    pd_fund_cpn* cpn;
} pd_fund_leg, *PD_FUND_LEG;

Err cpd_fill_fund_leg(
    /*	Coupons that started before today are disregarded */
    long today,
    /*	EOD Flag */
    int         eod_flag, /*	0: I, 1: E */
    double      fund_not,
    int         fund_ccy, /*	0: domestic, 1: foreign */
    int         fund_ncpn,
    long*       fund_fix,
    long*       fund_start,
    long*       fund_pay,
    char**      fund_basis,
    double*     fund_spr,
    double*     fund_mrg,
    PD_FUND_LEG fund_leg);

/*	Check dates consistency */
Err cpd_check_fund_leg(PD_FUND_LEG fund_leg);

/*	Free */
Err cpd_free_fund_leg(PD_FUND_LEG fund_leg);

/*	Structures and functions for the exotic leg */

/*	Power Dual cpn */
typedef struct
{
    long   start_date;  /*	Coupon start date */
    double start_time;  /*	Coupon start time */
    long   pay_date;    /*	Coupon pay date */
    double pay_time;    /*	Coupon pay time */
    long   fx_fix_date; /*	Fx fixing dates for coupon */
    double fx_fix_time; /*	Fx fixing times for coupon */
    long   fx_val_date; /*	Fx value dates for coupon */
    double fx_val_time; /*	Fx value times for coupon */
    /*	Coupon is of the form alpha + beta * Fx [floored, capped] */
    /*	Notional and coverage are already included in alpha, beta, floor, cap */
    double alpha;
    double beta;
    int    floored;
    double floor;
    int    capped;
    double cap;

    // Interp coupon specification: coupon is in the form C + W0*Fx + sum[Wi * max(Fx-Ki,0)]
    // Notional and coverage are already included in C and Wi
    int     use_opt_str;
    int     nstrikes;
    double  wcst;
    double  wspot;
    double* strikes;
    double* weights;

} pd_exo_cpn, *PD_EXO_CPN;

/*	Power Dual leg */
typedef struct
{
    /*	Always paid in domestic */
    int num_cpn; /*	Number of coupons, not including notional refund */
    /*	0..num_cpn-1 */
    pd_exo_cpn* cpn;
    /*	Notional refund */
    pd_exo_cpn not_ref; /*	Notional refund amount */
} pd_exo_leg, *PD_EXO_LEG;

Err cpd_fill_exo_leg(
    /*	Coupons that fixed before today are disregarded */
    long today,
    /*	EOD Flag */
    int     eod_flag, /*	0: I, 1: E */
    double  pd_not,
    int     pd_ncpn,
    long*   pd_fix,
    long*   pd_start,
    long*   pd_pay,
    char**  pd_basis,
    double* pd_alpha,
    double* pd_beta,
    int*    pd_floored,
    double* pd_floor,
    int*    pd_capped,
    double* pd_cap,

    //		pd interp coupon:
    int      use_cpn_opt_str,
    int*     pd_num_strikes,
    double*  pd_wcst,
    double*  pd_wspot,
    double** pd_strikes,
    double** pd_weights,

    /*		pd not refund */
    long   pd_not_ref_fix,
    double pd_not_ref_alpha,
    double pd_not_ref_beta,
    int    pd_not_ref_floored,
    double pd_not_ref_floor,
    int    pd_not_ref_capped,
    double pd_not_ref_cap,

    //		pd interp notional:
    int     use_not_opt_str,
    int     pd_not_num_strikes,
    double  pd_not_wcst,
    double  pd_not_wspot,
    double* pd_not_strikes,
    double* pd_not_weights,

    PD_EXO_LEG exo_leg);

/*	Check dates consistency */
Err cpd_check_exo_leg(PD_EXO_LEG exo_leg);

/*	Free */
Err cpd_free_exo_leg(PD_EXO_LEG exo_leg);

/*	Structures and functions for the calls */

/*	Call */
typedef struct
{
    int    call_type; /*	0: Call, 1: KO */
    double barrier;   /*	KO only */
    double smooth;    /*	Smoothing, barrier only */
    int    bar_type;  /*	0: up and in, 1: down and in */
    double fee;       /*  fee if deal is called in domestic currency */
    double extra_fee; /*	fee adjustment if use_sabr = 2 */
    double orig_fee;  /*	original fee saved here if use_sabr = 2 */
    int    pay_rec;   /*	0: rec pd upon exercise, 1: pay pd upon exercise */
    /*	Specs */
    long   ex_date;      /*	Exercise date */
    long   ex_date2bd;   //	Exercise date + 2bd
    double ex_time;      /*	Exercise time */
    double ex_time2bd;   //	Exercise time + 2bd
    int    fund_idx;     /*	Index of the first funding coupon to be called */
    int    num_fund_cpn; /*	Number of funding coupons called (including redemption) */
    int    pd_idx;       /*	Index of the first pd coupon to be called */
    int    num_pd_cpn;   /*	Number of pd coupons called (excluding redemption) */
    /*	Notional exchange upon exercise */
    long   set_date;     /*	Notional exchange date */
    double set_time;     /*	Notional exchange time */
    double fund_not_amt; /*	Amount of the funding leg notional refund */
    double pd_not_amt;   /*	Amount of the pd leg notional refund */

    int fx_bound; /*	do we optimise on the Fx or on the IV */

    // Exotic early redemption specification: redemption is in the form C + W0*Fx + sum[Wi *
    // max(Fx-Ki,0)] Notional is already included in C and Wi
    int     use_opt_str;
    int     nstrikes;
    double  wcst;
    double  wspot;
    double* strikes;
    double* weights;

    long   fx_fix_date; /*	Fx fixing date in case of exotic redemption */
    double fx_fix_time; /*	Fx fixing time in case of exotic redemption */
    long   fx_val_date; /*	Fx value date in case of exotic redemption */
    double fx_val_time; /*	Fx value time in case of exotic redemption */

    // Target note params
    int    TARN_Do;
    double orig_barrier;
} pd_call, *PD_CALL;

/*	Callable Power Dual */
typedef struct
{
    pd_exo_leg*  pd_leg;
    pd_fund_leg* fund_leg;
    int          type; /*	-1: Call KO, 0: Call only, 1: KO only */
    int          num_calls;
    pd_call*     call;
} cpd_str, *CPD_STR;

Err cpd_fill_calls(
    /*	Exercises before today are disregarded */
    long today,
    /*	EOD Flag */
    int   eod_flag, /*	0: I, 1: E */
    int   ncall,
    int*  type,    /*	0: call, 1: KO */
    int   pay_rec, /*	0: rec pd, 1: pay pd */
    long* ex_date,
    long* set_date,
    long* pd_not_ref_fix,  // Fx fixing dates in case of exotic redemption at call date
    /*	Notionals on both legs in their own currencies,
                    to be refunded upon call
                    in order to determine the Bond Strike */
    double* fund_not,
    double* pd_not,
    //		pd interp notional:
    int      use_not_opt_str,
    int*     pd_not_num_strikes,
    double*  pd_not_wcst,
    double*  pd_not_wspot,
    double** pd_not_strikes,
    double** pd_not_weights,

    double* barrier,  /*	KO only */
    int*    bar_type, /*	0: up and in, 1: down and in */
    double* fees,     /*  fees if deal is called in domestic currency */
    int     TARN_Do,  /*  1: Target note powerdual */
    double  smooth,   /*	Smoothing factor */
    int     fx_bound,
    int     use_bound,
    int     force_optim,
    CPD_STR cpd);

/*	Check dates consistency */
Err cpd_check_calls(CPD_STR cpd);

/*	Free */
Err cpd_free_calls(CPD_STR cpd);

/*	Structures and functions for the underlying and its term structures */

/*	Underlying term structs */
typedef struct
{
    char    name[256];
    long    today;
    double  spot_fx;
    char    dom_yc[256];
    char    for_yc[256];
    long    sigma_n_rates;
    double* sigma_date_rates;
    double* sigma_time_rates;
    double* sigma_dom;
    double  lda_dom;
    double* sigma_for;
    double  lda_for;
    long    sigma_n_fx;
    double* sigma_date_fx;
    double* sigma_time_fx;
    double* sigma_fx;
    double* corr_times;
    double* correl_dom_for;
    double* correl_dom_fx;
    double* correl_for_fx;
    long    corr_n_times;

    /* DLM model if needed */
    int                        use_beta_dlm;
    FxBetaDLM_model*           model;
    FxBetaDLM_Hermite*         hermite;
    FxBetaDLM_OptNumerParams*  num_params;
    FxBetaDLM_GRFNNumerParams* grfn_params;

    /* Extra Fees if required */
    int     nb_fees;
    double* fees_dates;
    double* fees;

} cpd_und, *CPD_UND;

/*	Fill underlying structure from a predefined underlying */
Err cpd_fill_und(
    char*            fx3dund,
    CPD_UND          und,
    CPDBETADLMPARAMS cpd_dlm_params,
    double           dom_vol_shift,
    double           for_vol_shift,
    double           fx_vol_shift);

/*	Fill underlying structure from calibration instruments */
Err cpd_calib_und(
    long today,
    /*	EOD Flag */
    int    eod_flag, /*	0: I, 1: E */
    double fx_spot,
    long   fx_spot_date,
    int    dom_calib,      /*	Calibrate domestic underlying */
    char*  dom_und,        /*	If no, domestic underlying to be used */
    char*  dom_yc,         /*	Domestic yc */
    char*  dom_vc,         /*	Domestic vc (only if calib) */
    char*  dom_ref,        /*	Domestic ref rate (only if calib) */
    char*  dom_swap_freq,  /*	Domestic swap freq (only if calib) */
    char*  dom_swap_basis, /*	Domestic swap basis (only if calib) */
    double dom_lam,        /*	Domestic lambda */
    int    for_calib,      /*	Same for foreign */
    char*  for_und,
    char*  for_yc,
    char*  for_vc,
    char*  for_ref,
    char*  for_swap_freq,
    char*  for_swap_basis,
    double for_lam,

    double min_fact,  /*	Maximum down jump on variance */
    double max_fact,  /*	Maximum up jump on variance */
    int    use_jumps, /*	Allow vol term structure to jump */

    double*          corr_times,
    double*          correl_dom_for, /*	Correlations */
    double*          correl_dom_fx,
    double*          correl_for_fx,
    long             corr_n_times,
    CPDBETADLMPARAMS cpd_dlm_params,
    CPD_STR          cpd,  /*	Structure */
    Err (*get_ir_cash_vol)(/*	Function to get IR cash vol from the markets */
                           char*   vol_curve_name,
                           double  start_date,
                           double  end_date,
                           double  cash_strike,
                           int     zero,
                           char*   ref_rate_name,
                           double* vol,
                           double* power),
    /*	Fx vol from the market */
    long*   fx_mkt_vol_date,
    double* fx_mkt_vol,
    int     num_fx_mkt_vol,
    CPD_UND und,
    double  dom_vol_shift,
    double  for_vol_shift,
    double  fx_vol_shift);

void cpd_copy_und(CPD_UND src, CPD_UND dest);

Err cpd_AlphaFudgeUpdateUnd(CPD_UND und, int UpOrDown); /* Up = 0 and Down = 1 */

Err cpd_free_und(CPD_UND und);

/*	Constants used for reconstruction and evaluation  */
/*	------------------------------------------------- */

#define MAX_NCPN 512
typedef struct
{
    int call_idx;

    /*	Discount Factor reconstruction information */

    /*	Spot fx */
    long   fx_spot_date;
    double fx_spot_time;
    double spot_fx_dff_log_ratio;
    double spot_fx_dom_gam;
    double spot_fx_cvx;
    double spot_fx_for_gam;

    /*	Funding */
    int    do_fund;
    int    num_fund_cpn;
    int    fund_var;
    double start_dff;
    double start_gam;
    double start_gam2;
    double fund_dff[MAX_NCPN];
    double fund_gam[MAX_NCPN];
    double fund_gam2[MAX_NCPN];

    double fund_I[MAX_NCPN];
    double fund_J[MAX_NCPN];
    double fund_K[MAX_NCPN];

    double fund_DI[MAX_NCPN];
    double fund_DJ[MAX_NCPN];
    double fund_DK[MAX_NCPN];

    /*	PD discount */
    int    do_pd_disc;
    int    num_pd_cpn;
    double pd_disc_dff[MAX_NCPN];
    double pd_disc_gam[MAX_NCPN];
    double pd_disc_gam2[MAX_NCPN];

    double pd_disc_I[MAX_NCPN];
    double pd_disc_J[MAX_NCPN];
    double pd_disc_K[MAX_NCPN];

    double pd_disc_DI[MAX_NCPN];
    double pd_disc_DJ[MAX_NCPN];
    double pd_disc_DK[MAX_NCPN];

    /*	PD coupons */
    int    do_pd_fwd;
    double pd_fwd_dff_ratio[MAX_NCPN];
    double pd_fwd_dom_gam[MAX_NCPN];
    double pd_fwd_cvx[MAX_NCPN];
    double pd_fwd_for_gam[MAX_NCPN];
    double pd_std[MAX_NCPN];
    double pd_half_std[MAX_NCPN];
    double pd_abs_beta[MAX_NCPN];
    int    pd_floor_type[MAX_NCPN];
    double pd_floor_str[MAX_NCPN];
    int    pd_cap_type[MAX_NCPN];
    double pd_cap_str[MAX_NCPN];

    int index_fwd[MAX_NCPN];
    int index_cap[MAX_NCPN];

    double pd_log_fwd_I[MAX_NCPN];
    double pd_log_fwd_J[MAX_NCPN];
    double pd_log_fwd_K[MAX_NCPN];

    double pd_log_fwd_DI[MAX_NCPN];
    double pd_log_fwd_DJ[MAX_NCPN];
    double pd_log_fwd_DK[MAX_NCPN];

    double pd_fwd_I[MAX_NCPN];
    double pd_fwd_J[MAX_NCPN];
    double pd_fwd_K[MAX_NCPN];

    double pd_fwd_DI[MAX_NCPN];
    double pd_fwd_DJ[MAX_NCPN];
    double pd_fwd_DK[MAX_NCPN];

    /*	Beta DLM extra informations */
    int                    nb_pd_inst;
    FxBetaDLM_FxOptInst*   pd_inst;
    FxBetaDLM_InstPrecalc* pd_beta_precalc;

    /*	DLM X reconstruction */
    double pd_X_const;
    double pd_X_dom;
    double pd_X_for;

    /*	DLM Spot reconstruction */
    double pd_S_const;
    double pd_S_lin;
    double pd_S_quad;
    double pd_S_dom;
    double pd_S_for;

    /*	Initial notional exchange */
    double in_not_fund_dff;
    double in_not_fund_gam;
    double in_not_fund_gam2;
    double in_not_pd_dff;
    double in_not_pd_gam;
    double in_not_pd_gam2;
    double fee;
    double extra_fee; /*	fee adjustment if use_sabr = 2 */
    double orig_fee;  /*	original fee saved here if use_sabr = 2 */

    //	Initial notional exchange in case of option string specification:
    double in_not_pd_fwd_dff_ratio;
    double in_not_pd_fwd_dom_gam;
    double in_not_pd_fwd_cvx;
    double in_not_pd_fwd_for_gam;
    double in_not_pd_std;
    double in_not_pd_half_std;

    /*	Next call settlement date notional exchange */
    double next_fund_dff;
    double next_fund_gam;
    double next_fund_gam2;
    double next_pd_dff;
    double next_pd_gam;
    double next_pd_gam2;

    /*	Next call start date notional exchange */
    double next_start_fund_dff;
    double next_start_fund_gam;
    double next_start_fund_gam2;

    /*	Final notional exchange */
    /*	Disc */
    double fin_not_disc_dff;
    double fin_not_disc_gam;
    double fin_not_disc_gam2;
    /*	Cpn */
    double fin_not_fwd_dff_ratio;
    double fin_not_fwd_dom_gam;
    double fin_not_fwd_cvx;
    double fin_not_fwd_for_gam;
    double fin_not_std;
    double fin_not_half_std;
    double fin_not_abs_beta;
    int    fin_not_floor_type;
    double fin_not_floor_str;
    int    fin_not_cap_type;
    double fin_not_cap_str;

    int fin_not_index_fwd;
    int fin_not_index_cap;

    /*	Beta DLM extra informations */
    FxBetaDLM_FxOptInst   fin_not_inst;
    FxBetaDLM_InstPrecalc fin_not_beta_precalc;

    /* Saved infos for Beta DLM */
    double beta_dom;
    double beta_for;
    double phi_dom;
    double phi_for;
    double varFFX;

} cpd_eval_const, *CPD_EVAL_CONST;

/*	Fill structure */
Err cpd_fill_eval_const(
    CPD_UND und,
    CPD_STR cpd,
    /*	Index of the current call/KO */
    int            call_idx,
    CPD_EVAL_CONST eval_const);

Err cpd_fill_eval_const_dlm(
    CPD_UND und,
    CPD_STR cpd,
    /*	Index of the current call/KO */
    int            call_idx,
    int            is_mc,
    CPD_EVAL_CONST eval_const);

void cpd_free_eval_const(CPD_UND und, CPD_STR cpd, CPD_EVAL_CONST eval_const);

/*	Arguments to all payoff evaluation functions */
/*	-------------------------------------------- */

typedef struct
{
    /*	Underlying */
    cpd_und* und;

    /*	Structure */
    cpd_str* cpd;

    /*	Constants used for reconstruction and evaluation */
    cpd_eval_const eval_const;

    /*	Index of the current call/KO */
    int call_idx;

} cpd_pay_arg, *CPD_PAY_ARG;

/*	Arguments to the tree function */
/*	------------------------------ */

typedef struct
{
    long    nstp;
    double* time;
    double* date;
    int*    vol_change;
    double* sig_dom;
    double* sig_for;
    double* sig_fx;
    double* dom_ifr;
    double* dom_fwd;
    double* dom_var;
    double* for_ifr;
    double* for_fwd;
    double* for_var;
    double* fx_fwd;
    double* fx_var;
    void**  void_prm;
    int*    is_event;
    double* bar_lvl;
    int*    bar_cl;
    int*    is_bar;
    double  dom_lam;
    double  for_lam;
    double* corr_dom_for;
    double* corr_dom_fx;
    double* corr_for_fx;
    double  spot_fx;
    char    dom_yc[256];
    char    for_yc[256];

    /* For Beta DLM */
    double* dr_const_dom;
    double* dr_coef_dom;
    double* dr_const_for;
    double* dr_coef1_for;
    double* dr_coef2_for;
    double* dr_coef3_for;
    double* dr_const_fx;
    double* dr_coef1_fx;
    double* dr_coef2_fx;
    double* dr_coef3_fx;

    double* glob_corr_dom_for;
    double* glob_corr_dom_fx;
    double* glob_corr_for_fx;

} cpd_tree_arg, *CPD_TREE_ARG;

Err cpd_fill_tree_arg(
    CPD_UND und,
    CPD_STR cpd,
    /*	Required number of time steps */
    long         req_stp,
    CPD_TREE_ARG tree_arg);

Err cpd_free_tree_arg(CPD_TREE_ARG tree_arg);

/*	Arguments to the MC function */
/*	---------------------------- */

typedef struct
{
    long    npaths;
    int     num_col;
    double* time;
    double* date;
    long    nb_dates;
    double* dom_ifr;
    double* dom_fwd;
    double* dom_std;
    double* dom_phi;
    double* dom_beta;
    double* dom_bond_pay;
    double* dom_beta_pay;
    double* for_ifr;
    double* for_fwd;
    double* for_std;
    double* for_phi;
    double* for_beta;
    double* fx_fwd;
    double* fx_std;
    double* dom_for_cov;
    double* dom_fx_cov;
    double* for_fx_cov;
    void**  void_prm;
    double  dom_lam;
    double  for_lam;
    double* corr_times;
    double* corr_dom_for;
    double* corr_dom_fx;
    double* corr_for_fx;
    long    corr_n_times;
    double  spot_fx;
    char    dom_yc[256];
    char    for_yc[256];
    int     do_pecs;
    double  smooth;

    /* For Beta DLM */
    double* for_fwd_const;
    double* for_fwd_lin;
    double* ffx_var;

} cpd_mc_arg, *CPD_MC_ARG;

Err cpd_fill_mc_arg(
    CPD_UND und,
    CPD_STR cpd,
    /*	Required number of paths */
    long req_pth,
    /*	Do PECS */
    int        do_pecs,
    CPD_MC_ARG mc_arg);

Err cpd_free_mc_arg(CPD_MC_ARG mc_arg);

/*	Main function to be called in order to fill and check all structures */
/*	==================================================================== */

/*	Fill and check all the relevant structures */
Err cpd_fill_check_all_struct(
    /*	Today's date */
    long today,
    /*	The underlying */
    int use_calib, /*	0: use fx3dund, 1: calibrate */
    /*		if calib */
    double fx_spot,
    long   fx_spot_date,
    int    dom_calib,      /*	Calibrate domestic underlying */
    char*  dom_und,        /*	If no, domestic underlying to be used */
    char*  dom_yc,         /*	Domestic yc */
    char*  dom_vc,         /*	Domestic vc (only if calib) */
    char*  dom_ref,        /*	Domestic ref rate (only if calib) */
    char*  dom_swap_freq,  /*	Domestic swap freq (only if calib) */
    char*  dom_swap_basis, /*	Domestic swap basis (only if calib) */
    double dom_lam,        /*	Domestic lambda */
    int    for_calib,      /*	Same for foreign */
    char*  for_und,
    char*  for_yc,
    char*  for_vc,
    char*  for_ref,
    char*  for_swap_freq,
    char*  for_swap_basis,
    double for_lam,

    double min_fact,  /*	Maximum down jump on variance */
    double max_fact,  /*	Maximum up jump on variance */
    int    use_jumps, /*	Allow vol term structure to jump */

    double*          corr_times,
    double*          correl_dom_for, /*	Correlations */
    double*          correl_dom_fx,
    double*          correl_for_fx,
    long             corr_n_times,
    CPDBETADLMPARAMS cpd_dlm_params,
    Err (*get_ir_cash_vol)(/*	Function to get IR cash vol from the markets */
                           char*   vol_curve_name,
                           double  start_date,
                           double  end_date,
                           double  cash_strike,
                           int     zero,
                           char*   ref_rate_name,
                           double* vol,
                           double* power),
    /*	Fx vol from the market */
    long*   fx_mkt_vol_date,
    double* fx_mkt_vol,
    int     num_fx_mkt_vol,
    /*		if no calilb */
    char* fx3dund,
    /*	The structure */
    /*		funding */
    double  fund_not,
    int     fund_ccy, /*	0: domestic, 1: foreign */
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    /*		pd */
    double  pd_not,
    int     pd_ncpn,
    long*   pd_fix,
    long*   pd_start,
    long*   pd_pay,
    char**  pd_basis,
    double* pd_alpha,
    double* pd_beta,
    int*    pd_floored,
    double* pd_floor,
    int*    pd_capped,
    double* pd_cap,

    //		pd interp coupon:
    int      use_cpn_opt_str,
    int*     pd_num_strikes,
    double*  pd_wcst,
    double*  pd_wspot,
    double** pd_strikes,
    double** pd_weights,

    /*		pd not refund */
    long*  pd_not_ref_fix,
    double pd_not_ref_alpha,
    double pd_not_ref_beta,
    int    pd_not_ref_floored,
    double pd_not_ref_floor,
    int    pd_not_ref_capped,
    double pd_not_ref_cap,

    //		pd interp notional:
    int      use_not_opt_str,
    int*     pd_not_num_strikes,
    double*  pd_not_wcst,
    double*  pd_not_wspot,
    double** pd_not_strikes,
    double** pd_not_weights,

    /*		calls */
    int*    call_type, /*	0: call, 1: KO */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* barrier,  /*	KO only */
    int*    bar_type, /*	0: up and in, 1: down and in */
    double* fees,     /*  fees if deal is called in domestic currency */
    int     TARN_Do,

    /*	Numerical params */
    long   req_stp,
    long   req_pth,
    int    do_pecs,
    int    forcetree,
    int    do_optim,    /*	If equal to 1 then the call are replaced by optimal KO	*/
    int    force_optim, /*	If equal to 1 then all call will be replaced by optimal KO	*/
    int    fx_bound,    /*	If equal to 1 then optimisation on the Fx, on the IV otherwise	*/
    int    use_bound,
    double smooth,
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */
    /*	Results */
    CPD_STR cpd,
    CPD_UND und,
    int*    call_feat, /*	0: No callable feature to be valued
                                       1: Callable feature to be valued through tree
                                       2: Callable feature to be valued through MC	*/

    CPD_TREE_ARG tree_arg,
    CPD_MC_ARG   mc_arg,
    double       dom_vol_shift,
    double       for_vol_shift,
    double       fx_vol_shift,
    int          skip_fill); /* to skip tree_arg and mc_arg fill */

/*	Free all structures */
Err cpd_free_all_struct(CPD_STR cpd, CPD_UND und, CPD_TREE_ARG tree_arg, CPD_MC_ARG mc_arg);

/*	Payoff function for tree (callable) */
/*	---------------------------------	*/

Err cpd_payoff_4_3dfx_tree(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double spot_fx,
    void*  dom_yc,
    double dom_lam,
    double dom_phi,
    void*  for_yc,
    double for_lam,
    double for_phi,
    /* Nodes data */
    long n1,
    long n2,
    long n3,
    /* i: d1, j: d2, k: d3, l = {0: xDom, 1: xFor, 2: log (Fx/Fx0)} */
    double**** sv,
    /* Vector of results to be updated */
    long       nprod,
    double**** prod_val,
    /* Barrier details */
    int    is_bar,
    int    bar_k,
    int**  bar_idx,
    int    bar_col,
    double bar_lvl);

/*	Main pricing function for tree */

/*	Launch the tree */
Err cpd_launch_tree(
    CPD_STR      cpd,
    CPD_UND      und,
    CPD_TREE_ARG tree_arg,
    /*	Result */
    double* prem,
    double* IV);

/*	Launch the tree ko */
Err cpd_launch_tree_ko(
    CPD_STR      cpd,
    CPD_UND      und,
    CPD_TREE_ARG tree_arg,
    /*	Result */
    double* prem);

/*	Payoff function for mc (KO) */
/*	---------------------------	*/

Err cpd_fut_intr_val(
    /* State variables */
    double Xdomfor[], /*	[0]: Xdom, [1]: Xfor */
    double Zfx,
    double spot_fx,
    /* Structure */
    CPD_STR cpd,
    /*	Call */
    int    fund_idx,
    int    num_fund_cpn,
    double fund_not_amt,
    int    pd_idx,
    int    num_pd_cpn,
    double pd_not_amt,
    /* Precalculated constants */
    CPD_EVAL_CONST eval_const,
    /* Results */
    double* fund_leg,
    double* pd_leg);

Err cpd_current_pd_cpn_val(
    /* State variables */
    double Xdomfor[], /*	[0]: Xdom, [1]: Xfor */
    double Zfx,
    double spot_fx,
    /* Structure */
    CPD_STR cpd,
    /*	Call */
    int    pd_idx,
    int    num_pd_cpn,
    double pd_not_amt,
    /* Precalculated constants */
    CPD_EVAL_CONST eval_const,
    /* Results */
    double* pd_leg);

void cpd_init_4_3dfx_mc(void);

Err cpd_payoff_4_3dfx_mc_dlm(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    double Xdom,
    double Yfor,
    double Zfx,
    /* Results */
    int     num_col,
    double* res,
    int*    stop_path);

// Target Note
Err cpd_payoff_4_3dfx_tarn(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double spot_fx,
    void*  dom_yc,
    double dom_lam,
    double dom_phi,
    void*  for_yc,
    double for_lam,
    double for_phi,
    double Xdom,
    double Yfor,
    double Zfx,
    /* Results */
    int     num_col,
    double* res,
    int*    stop_path);

/* MC but on total swap ! */
Err cpd_payoff_total_4_3dfx_mc_dlm(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    double Xdom,
    double Yfor,
    double Zfx,
    /* Results */
    int     num_col,
    double* res,
    int*    stop_path);

Err cpd_payoff_4_3dfx_mc(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double spot_fx,
    void*  dom_yc,
    double dom_lam,
    double dom_phi,
    void*  for_yc,
    double for_lam,
    double for_phi,
    double Xdom,
    double Yfor,
    double Zfx,
    /* Results */
    int     num_col,
    double* res,
    int*    stop_path);

Err cpd_payoff_4_3dfx_optiFx_mc(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double spot_fx,
    void*  dom_yc,
    double dom_lam,
    double dom_phi,
    void*  for_yc,
    double for_lam,
    double for_phi,
    double Xdom,
    double Yfor,
    double Zfx,
    /* Results */
    int     num_col,
    double* res,
    int*    stop_path);

/*	Main pricing function for mc */

/*	Launch the mc */
Err cpd_launch_mc(
    CPD_STR    cpd,
    CPD_UND    und,
    CPD_MC_ARG mc_arg,
    /*	Result */
    double* prem,
    double* std);

Err cpd_launch_opti_mc(
    CPD_STR    cpd,
    CPD_UND    und,
    CPD_MC_ARG mc_arg,
    /*	Result */
    double*   prem,
    double*   std,
    int       do_infos,
    double*** optim_bar);

/*	Functions used to deal with the Fx smile using a GRFN tableau generation */

/*	A reduced version of cpd_fill_check_all_struct that does only the cpd (product) structure */
Err cpd_fill_check_prod_struct(
    /*	Today's date */
    long today,
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */
    /*	The structure */
    /*		funding */
    double  fund_not,
    int     fund_ccy, /*	0: domestic, 1: foreign */
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    /*		pd */
    double  pd_not,
    int     pd_ncpn,
    long*   pd_fix,
    long*   pd_start,
    long*   pd_pay,
    char**  pd_basis,
    double* pd_alpha,
    double* pd_beta,
    int*    pd_floored,
    double* pd_floor,
    int*    pd_capped,
    double* pd_cap,
    /*		pd not refund */
    long   pd_not_ref_fix,
    double pd_not_ref_alpha,
    double pd_not_ref_beta,
    int    pd_not_ref_floored,
    double pd_not_ref_floor,
    int    pd_not_ref_capped,
    double pd_not_ref_cap,
    /*		calls */
    int*    call_type, /*	0: Call, 1: KO */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* barrier,  /*	KO only */
    int*    bar_type, /*	0: up and in, 1: down and in */
    double* fees,     /*  fees if deal is called in domestic currency */
    double  smooth,
    /*	Result */
    CPD_STR cpd);

/*	Free product structure */
Err cpd_free_prod_struct(CPD_STR cpd);

/*	Function that produces a GRFN tableau for a CPD deal */

#define MAX_NEVT 256
#define MAX_COL 25
#define MAX_STR_LEN 256
#define TOTAL_SIZE 1638400
#define NO_FLOOR -99999999999
#define NO_CAP 99999999999

/*	Function that produces a GRFN tableau for a CPD deal */
Err cpd_make_grfn_tableau(
    /*	The structure */
    CPD_STR cpd,
    /*	Extra cash-flows coming from past coupons and initial exchange */
    int     fund_ncf,
    long*   fund_cfdtes,
    double* fund_cf,
    int     pd_ncf,
    long*   pd_cfdtes,
    double* pd_cf,
    /*	The underlyings */
    long  today,
    char* dom_name,
    char* for_name,
    char* fx_name,
    /*	The tableau and auxiliaries
            Must be allocated with maximum values
                    and initialised with zeros */
    int*     num_evt_dtes,
    int*     num_cols,
    int*     aux_rows,
    int*     aux_cols,
    long*    evt_dtes,
    char***  evt,
    int**    mask,
    double** aux);

/*	Function that produces a GRFN tableau for a CPD deal */
Err cpd_make_grfn_tableau_credit(
    /*	The structure */
    CPD_STR cpd,
    /*	Extra cash-flows coming from past coupons and initial exchange */
    int     fund_ncf,
    long*   fund_cfdtes,
    double* fund_cf,
    int     pd_ncf,
    long*   pd_cfdtes,
    double* pd_cf,
    /*	The underlyings */
    long  today,
    char* dom_name,
    char* for_name,
    char* fx_name,

    /*	Credit inputs */
    char* risky_yc,
    long  credit_freq,
    /*	The tableau and auxiliaries
            Must be allocated with maximum values
                    and initialised with zeros */
    int*     num_evt_dtes,
    int*     num_cols,
    int*     aux_rows,
    int*     aux_cols,
    long*    evt_dtes,
    char***  evt,
    int**    mask,
    double** aux);

/*	Structures and functions for the SMILEd underlying and its term structures */

/*	Underlying term structs */
typedef struct
{
    cpd_und und;

    double alpha;
    double beta;
    double fx2bdfwd;

} cpd_smile_und, *CPD_SMILE_UND;

/*	Fill underlying structure from a predefined underlying */
Err cpd_fill_smile_und(
    char*         fx3dund,
    double        alpha,
    double        beta,
    CPD_SMILE_UND und,
    double        dom_vol_shift,
    double        for_vol_shift,
    double        fx_vol_shift);

/*	Fill underlying structure from calibration instruments */
Err cpd_calib_smile_und(
    long today,
    /*	EOD Flag */
    int     eod_flag, /*	0: I, 1: E */
    double  fx_spot,
    long    fx_spot_date,
    int     dom_calib,      /*	Calibrate domestic underlying */
    char*   dom_und,        /*	If no, domestic underlying to be used */
    char*   dom_yc,         /*	Domestic yc */
    char*   dom_vc,         /*	Domestic vc (only if calib) */
    char*   dom_ref,        /*	Domestic ref rate (only if calib) */
    char*   dom_swap_freq,  /*	Domestic swap freq (only if calib) */
    char*   dom_swap_basis, /*	Domestic swap basis (only if calib) */
    double  dom_lam,        /*	Domestic lambda */
    int     for_calib,      /*	Same for foreign */
    char*   for_und,
    char*   for_yc,
    char*   for_vc,
    char*   for_ref,
    char*   for_swap_freq,
    char*   for_swap_basis,
    double  for_lam,
    double* corr_times,
    double* correl_dom_for, /*	Correlations */
    double* correl_dom_fx,
    double* correl_for_fx,
    long    corr_n_times,
    double  alpha,
    double  beta,
    CPD_STR cpd,           /*	Structure */
    Err (*get_ir_cash_vol)(/*	Function to get IR cash vol from the markets */
                           char*   vol_curve_name,
                           double  start_date,
                           double  end_date,
                           double  cash_strike,
                           int     zero,
                           char*   ref_rate_name,
                           double* vol,
                           double* power),
    /*	Fx vol from the market */
    long*         fx_mkt_vol_date,
    double*       fx_mkt_vol,
    int           num_fx_mkt_vol,
    CPD_SMILE_UND und,
    double        dom_vol_shift,
    double        for_vol_shift,
    double        fx_vol_shift);

void cpd_copy_smile_und(CPD_SMILE_UND src, CPD_SMILE_UND dest);

Err cpd_free_smile_und(CPD_SMILE_UND und);

/*	Launch the smile tree */
Err cpd_launch_smile_tree(
    CPD_STR       cpd,
    CPD_SMILE_UND und,
    /*	Number of steps */
    int num_stp,
    /*	GRFN tableau structures */
    int*     num_evt_dtes,
    int*     num_cols,
    int*     aux_rows,
    int*     aux_cols,
    long*    evt_dtes,
    char***  evt,
    int**    mask,
    double** aux,
    /*	Results */
    double* fund_leg,
    double* pd_leg,
    double* call);

Err call_GRFN_alphabetaModel(
    /* Input */
    CPD_UND cpd_und, /* The CPD underlying */
    CPD_STR cpd_str,

    int   use_calib, /*	0: use fx3dund, 1: calibrate */
    char* fx_name,
    int   dom_calib, /*	Calibrate domestic underlying */
    char* dom_name,  /*	If no, domestic underlying to be used */
    char* dom_yc,    /*	Domestic yc */
    int   for_calib, /*	Same for foreign */
    char* for_name,
    char* for_yc,

    /*	The structure */
    long start_date, /*	Date at which initial notional exchange occurs */

    /*		funding */
    double  fund_not,    /*	Notional */
    int     for_fund,    /*	1: if funding has been transform into domestic */
    double  eq_final_ex, /*	equivalent final exchange use only if for_fund = 1 */
    double  eq_init_ex,  /*	equivalent init exchange use only if for_fund = 1 */
    long    fund_start_date,
    int     fund_ncpn,    /*	Number of coupons */
    long*   fund_fix,     /*	Fixing dates */
    long*   fund_start,   /*	Start dates */
    long*   fund_pay,     /*	Pay dates */
    char**  fund_basis,   /*	Basis */
    double* fund_mrg,     /*	Margins */
    double* fund_fix_cpn, /*	Past coupon fixing if relevant,
                                                                  includes spr, but not mrg, cvg and
                             notional */
    /*		pd */
    double  pd_not,   /*	Notional */
    int     pd_ncpn,  /*	Number of coupons */
    long*   pd_fix,   /*	Fx fixing dates */
    long*   pd_start, /*	Start dates */
    long*   pd_pay,   /*	Pay dates */
    char**  pd_basis, /*	Basis */
    double* pd_alpha, /*	Coupon = alpha + beta * fx [capped, floored] */
    double* pd_beta,
    int*    pd_floored,
    double* pd_floor,
    int*    pd_capped,
    double* pd_cap,
    double* pd_fix_fx, /*	Past Fx fixing if relevant */

    /* Call */
    int     ncall,   /*	Number of calls */
    long*   ex_date, /*	Call dates */
    double* barrier, /*	in case of a pure KO or a Callable KO */

    /*	EOD Fixing Flag */
    int eod_fix_flag, /*	0: I, 1: E */
    /*	EOD Payment Flag */
    int eod_pay_flag, /*	0: I, 1: E */

    /* Model Parameters */
    double alpha,
    double beta,

    /* Numerical parameter */
    int nb_step, /* Number of step for the tree */

    /* Fx vol from the market */
    long*   fx_mkt_vol_date,
    double* fx_mkt_vol,

    int    num_fx_mkt_vol,
    double fx_spot, /* Spot Fx */

    /* OutPut */
    double* fund_val, /*	Value of the funding leg */
    double* pd_val,   /*	Value of the Power Dual leg */
    double* call_val, /*	Value of the callable feature */
    double* call_stdev /*	Standard deviation of the call if applicable */);

Err cpd_payoff_4_3dfx_BetaDLM_tree(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double spot_fx,
    void*  dom_yc,
    void*  for_yc,
    /* Nodes data */
    long n1,
    long n2,
    long n3,
    /* i: d1, j: d2, k: d3, l = {0: xDom, 1: xFor, 2: log (Fx/Fx0)} */
    double**** sv,
    /* Vector of results to be updated */
    long       nprod,
    double**** prod_val,
    /* Barrier details */
    int    is_bar,
    int    bar_k,
    int**  bar_idx,
    int    bar_col,
    double bar_lvl);

#endif