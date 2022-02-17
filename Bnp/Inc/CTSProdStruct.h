
#ifndef __CTS_PROD_STRUCT_H
#define __CTS_PROD_STRUCT_H

#define CTS_MINBARRIER -0.01
#define CTS_MAXBARRIER 0.50
#define CTS_NA -9999.99
#define CTS_EPS 1.0e-16
#define CTS_IS_NA(X) (((X)-CTS_NA) < CTS_EPS)
#define MAXDF 1000
#define MAXTS 7500
#define CTS_MAX_STR 10

#include "DiagCalibDLM.h"
#include "LGMSVCalibApprox.h"
#include "MCEBOptimisation.h"

/*	Market structure */

typedef struct {
  long today;

  char yc[256];
  char vc[256];
  char ref[256]; /*	Reference rate used for getting the vol */
  char swap_freq[256];
  char swap_basis[256];
  SrtCompounding swap_srt_freq;
  SrtBasisCode swap_srt_basis;
  char lib_freq[256];
  char lib_basis[256];
  SrtCompounding lib_srt_freq;
  SrtBasisCode lib_srt_basis;

  GET_CASH_VOL_FUNC get_cash_vol;
  GET_CASH_VOL_CONVERT_FUNC get_cash_vol_convert;
  SrtDiffusionType vol_type;
  int cash_vol;

  /*	For DRS adjustment      , strikes in the market volatility matrix for
   * replication */
  int num_strikes_in_vol;
  double *strikes_in_vol;

} cts_mkt, *CTS_MKT, lgmsv_mkt, *LGMSV_MKT;

/*	Init */
Err cts_init_mkt(long today, char *yc, char *vc, char *ref, char *swap_freq,
                 char *swap_basis, char *lib_freq, char *lib_basis,
                 GET_CASH_VOL_FUNC get_cash_vol,
                 GET_CASH_VOL_CONVERT_FUNC get_cash_vol_convert,
                 SrtDiffusionType vol_type, int cash_vol,
                 int num_strikes_in_vol, double *strikes_in_vol, CTS_MKT mkt);

/*	Copy */
Err cts_copy_mkt(CTS_MKT src, CTS_MKT dest);

/*	Free */
void cts_free_mkt(CTS_MKT mkt);

/*	Structures and functions for callable cts floaters */
/*	-------------------------------------------------- */

/*	Structures and functions for the funding leg */

/*	Funding cpn */
typedef struct {
  double not ;                /*	Notional */
  long start_date;            /*	Coupon start date */
  double start_time;          /*	Coupon start time */
  long end_date;              /*	Coupon end date */
  double end_time;            /*	Coupon end time */
  long pay_date;              /*	Coupon pay date */
  double pay_time;            /*	Coupon pay time */
  SrtBasisCode basis;         /*	coupon basis */
  double cvg;                 /*	Cvg */
  double cpn;                 /*	Notional * (fwd_spr + margin) * cvg */
  double cpn_plus_ex;         /*  Coupon including notional repay/increase */
  double cpn_plus_ex_partial; /*  Coupon including notional repay/increase for
                                 partial */

  /*	PV today */
  double mkt_val;
} cts_fund_cpn, *CTS_FUND_CPN;

/*	Funding leg */
typedef struct {
  int num_cpn; /*	Number of coupons */
  /*	0..num_cpn-1 */
  cts_fund_cpn *cpn;
} cts_fund_leg, *CTS_FUND_LEG;

Err cts_fill_fund_leg(
    CTS_MKT mkt,
    /*	Coupons that started before today are disregarded */
    /*	EOD Flag */
    int eod_flag, /*	0: I      , 1: E */
    int fund_ncpn, double *fund_not, long *fund_fix, long *fund_start,
    long *fund_end, long *fund_pay, char **fund_basis, double *fund_spr,
    double *fund_mrg, CTS_FUND_LEG fund_leg);

/*	Check dates consistency */
Err cts_check_fund_leg(CTS_FUND_LEG fund_leg);

/*	Free */
void cts_free_fund_leg(CTS_FUND_LEG fund_leg);

/*	Structures and functions for the exotic leg */

/*	Fixing */
typedef struct {
  /*	Fixing date */
  long ref_fix_date;
  double ref_fix_time;
  long ref_theo_end_date;

  /*	Reference CMS */
  CalibCpnScheduleDLM schedule;

  double ref_sum_df;
  double ref_level;
  double ref_swp_cash;
  double ref_fwd_spread;

  /*	Params for CMS adjustment in market value */
  double ref_cpnd;
  double ref_conv;

} cts_fix, *CTS_FIX;

/*	Fixing Initialisation */
Err cts_init_fixing(long today, char *yc, int ref_fix_lag, long ref_start_date,
                    long ref_theo_end_date, SrtCompounding srt_freq,
                    SrtBasisCode srt_basis, double ref_fwd_spread,
                    CTS_FIX fixing);

/*	Check */
Err cts_check_fixing(CTS_FIX fixing);

/*	Coupon */
typedef struct {
  /*	Libor paid if relevant */
  long pay_fix_date;
  double pay_fix_time;
  long pay_start_date;
  double pay_start_time;
  long pay_theo_end_date;
  long pay_end_date;
  double pay_end_time;
  double pay_cvg;
  double pay_cpnd;
  double pay_conv;
  double pay_fwd_spread;
  double pay_gearing;

  /*	Coupon payment */
  long cpn_start_date;
  double cpn_start_time;
  long cpn_end_date;
  double cpn_end_time;
  long cpn_pay_date;
  double cpn_pay_time;
  double cpn_coupon;
  double cpn_cvg;
  double cpn_not;
  SrtBasisCode cpn_basis;

  /*	Option structure */
  /*	alpha + beta * fwd + string of options on forward */
  double alpha;
  double beta;
  int nstr;
  double *str;
  double *nbopt;

  /*	Original fixings */
  int nfix;
  double *weights;
  cts_fix *fix;

  /*	Used fixings */
  int used_nfix;
  int *used_fixidx;
  double *used_fixweights;

  int type;               /*	0:	fixed
                                          1:	libor fixed at start
                                          2:	libor fixed at end */
  double swaption_strike; /*	Implied swaption strike or CTS_NA */
  double caplet_strike;   /*	Implied caplet strike or CTS_NA */

  /*	Market value */
  double mkt_val;
  double mkt_fixed_part;
  double mkt_float_part;

  double mdl_val;

} cts_cpn, *CTS_CPN;

/* Initialisation of the payoff function */

/* General case */
Err cts_init_cpn_payoff(double alpha, double beta, int nstr, double *str,
                        double *nbopt, CTS_CPN cpn);

/* For a range */
Err cts_init_cpn_payoff_for_range(
    int buy_sell,   /*	1: BNPP buys      , -1: BNPP sells */
    int value_zero, /*	0: Do not value low strike options      , 1: do */
    double min_barrier, double max_barrier,
    double lb,     /*	Lower bound of range */
    double ub,     /*	Upper bound of range */
    double payoff, /*	Payoff in range */
    double call_spread, CTS_CPN cpn);

/* payoff = (gearing * PayLibor + margin) floored and capped */
Err cts_init_cpn_payoff_for_cif(int value_zero, double min_barrier,
                                double max_barrier, double gearing,
                                double margin, double floor, double cap,
                                CTS_CPN cpn);

/*	Adjust call spread values */
Err cts_adjust_call_spread(
    int buy_sell,   /*	1: BNPP buys      , -1: BNPP sells */
    int value_zero, /*	0: Do not value low strike options      , 1: do */
    double min_barrier, double max_barrier,
    double lb,     /*	Lower bound of range */
    double ub,     /*	Upper bound of range */
    double payoff, /*	Payoff in range */
    CTS_CPN cpn, CTS_MKT mkt, double call_spread, double numer_spread,
    double precision, int maxiter);

/*	Init */
Err cts_init_cpn(
    /*	Market */
    CTS_MKT mkt,
    /*	Coupon details */
    long cpn_start_date, long cpn_end_date, long cpn_pay_date,
    double cpn_coupon, char *cpn_basis, double cpn_not,
    /*	General libor fixing and basis properties */
    int fix_lag_bd,
    /*	Coupon type */
    int cpn_type, /*	0:	fixed
                              1:	libor fixed at start
                              2:	libor fixed at end */
    /*	Paid libor details */
    long pay_fix_date, /*	fixing date of the Libor */
    int pay_months,    /*	Length of the paid libor in months */
    char *pay_freq, char *pay_basis, double pay_fwd_spread, double pay_margin,
    /*	Fixing details */
    int nfix, double *weights,
    /*		ref libor */
    long *ref_fix_dates, /*	Fixing dates */
    char *fix_tenor,     /*	Fixing tenors */
    char *fix_freq, char *fix_basis,
    double *ref_fwd_spreads, /*	Reference libor spreads corresponding to fixing
                            dates */
    /*	Profile details */
    int tstype, /*	0: generic      , 1: range */
    /*		Profile 0 */
    double alpha, double beta, int nstr, double *str, double *nbopt,
    /*		Profile 1 */
    int buy_sell,   /*	1: BNPP buys      , -1: BNPP sells */
    int value_zero, /*	0: Do not value low strike options      , 1: do */
    double lb,      /*	Lower bound of range */
    double ub,      /*	Upper bound of range */
    double payoff,  /*	Payoff in range */
    double call_spread, double numer_spread,
    /*		Trimming */
    int trim_type, /*	0: no trim
                               1: x fixings max
                               2: x time min between two fixings */
    int max_fix, double min_fix_time, int calc_mkt_iv,
    int use_cmsopt, /*	Use CmsOption to value fix option      , use BS on
                   CmsRate otherwise */
    /*		Correl */
    double correl_start, /*	Correl between libor fixing and libor paid at
                        start */
    double
        correl_end, /*	Correl between libor fixing and libor paid at start */
    int float_adjust_type, /*	type of adjustment for the floating coupon ,
                               0: ATM vol      , 1: Strike Vol */
    CTS_CPN cpn);

/*	Trim fixings */
Err cts_trim_fixings(
    CTS_CPN cpn, int trim_type, /*	0: no trim
                                            1: x fixings max
                                            2: x time min between two fixings */
    int max_fix, double min_fix_time);

Err cts_calc_fixed_cpn_mkt_value(
    /*	Market */
    CTS_MKT mkt,
    /*	Coupon details */
    long cpn_start_date, long cpn_end_date, long cpn_pay_date,
    double cpn_coupon, char *cpn_basis, double cpn_not,
    /*	General libor fixing and basis properties */
    int fix_lag_bd,
    /*	Coupon type */
    int cpn_type, /*	0:	fixed
                              1:	libor fixed at start
                              2:	libor fixed at end
                              3:	midat */
    /*	Paid libor details */
    long pay_fix_date, /*	fixing date of the Libor */
    int pay_months,    /*	Length of the paid libor in months */
    char *pay_freq, char *pay_basis, double pay_fwd_spread, double pay_gearing,
    double fixed_pay,
    /*	Fixing details */
    int nfix, double *weights,
    /*		ref libor */
    double *fixed_ref, /*	Fixing dates */
    /*	Profile details */
    int tstype, /*	0: generic      , 1: range */
    /*		Profile 0 */
    double alpha, double beta, int nstr, double *str, double *nbopt,
    /*		Profile 1 */
    double lb,                                  /*	Lower bound of range */
    double ub,                                  /*	Upper bound of range */
    double payoff,                              /*	Payoff in range */
    double accrue_on_barrier, int eod_fix_flag, /*	0: I      , 1: E */
    int eod_pay_flag,                           /*	0: I      , 1: E */
    int eod_ex_flag,                            /*	0: I      , 1: E */
    double *mkt_value);

/*	Calc market values */
Err cts_calc_cpn_mkt_value(
    CTS_CPN cpn, CTS_MKT mkt,
    int use_cmsopt,      /*	Use CmsOption to value fix option      , use BS on
                        CmsRate      otherwise */
    double correl_start, /*	Correl between libor fixing and libor paid at
                        start */
    double
        correl_end, /*	Correl between libor fixing and libor paid at start */
    int float_adjust_type); /*	type of adjustment for the floating coupon ,
                                            0: ATM vol      , 1: Strike Vol */

/*	Check */
Err cts_check_coupon(CTS_CPN cpn);

/*	Free */
void cts_free_coupon(CTS_CPN cpn);

/*	Exotic leg */
typedef struct {
  int num_cpn; /*	Number of coupons */
  /*	0..num_cpn-1 */
  cts_cpn *cpn;

} cts_exo_leg, *CTS_EXO_LEG;

/*	Init */
Err cts_fill_cpn_leg(
    CTS_MKT mkt,
    /*	Coupons that fixed before today are disregarded */
    /*	EOD Flag */
    int eod_flag, /*	0: I      , 1: E */
    /*	General Libor properties */
    int fix_lag_bd,
    /*	Number of coupons */
    int ncpn,
    /*	Coupon description */
    /*		cpn */
    int *cpn_type, /*	0:	fixed
                               1:	libor fixed at start
                               2:	libor fixed at end */
    long *cpn_start_date, long *cpn_end_date, long *cpn_pay_date,
    double *cpn_coupon, char **cpn_basis, double *cpn_not,
    /*		pay libor */
    long *pay_fix_date, /*	fixing date of the Libor */
    int pay_months, char *pay_freq, char *pay_basis, double *pay_fwd_spread,
    double *pay_margin,
    /*		fix libor */
    /*	fixings */
    int *nfix, double **weights, long **ref_fix_dates,
    char **fix_tenor, /*	Fixing tenors */
    char *fix_freq, char *fix_basis, double **ref_fwd_spreads,
    /*	Profiles */
    int tstype, /*	0: generic      , 1: range */
    /*	profile 0 */
    double *alpha, double *beta, int *nstr, double **str, double **nbopt,
    /*	profile 1 */
    int buy_sell,   /*	1: BNPP buys      , -1: BNPP sells */
    int value_zero, /*	0: Do not value low strike options      , 1: do */
    double *lb,     /*	Lower bounds of ranges */
    double *ub,     /*	Upper bounds of ranges */
    double *payoff, /*	Payoff in ranges */
    double call_spread, double numer_call_spread,
    /*		Trimming */
    int trim_type, /*	0: no trim
                               1: x fixings max
                               2: x time min between two fixings */
    int max_fix, double min_fix_time,
    /*	Extra model parameters*/
    int calc_mkt_iv, int use_cmsopt, /*	Use CmsOption to value fix option      ,
                                use BS on CmsRate otherwise */
    double correl_start, /*	Correl between libor fixing and libor paid at
                        start */
    double
        correl_end, /*	Correl between libor fixing and libor paid at start */
    int float_adjust_type, /*	type of adjustment for the floating coupon , 0:
                              ATM vol      , 1: Strike Vol
                        */
    CTS_EXO_LEG exo_leg);

/*	Check */
Err cts_check_exo_leg(CTS_EXO_LEG exo_leg);

/*	Free */
void cts_free_exo_leg(CTS_EXO_LEG exo_leg);

/*	Structures and functions for the calls */

typedef struct {
  int pay_rec; /*	0: rec exo leg upon exercise      , 1: pay exo leg upon
                  exercise */
  /*	Specs */
  int index;       /*	Index of the call */
  long ex_date;    /*	Exercise date */
  double ex_time;  /*	Exercise time */
  int exo_idx;     /*	Index of the first exotic coupon to be called */
  int num_exo_cpn; /*	Number of exotic coupons called */
  int num_partial_exo_cpn;
  int fund_idx;     /*	Index of the first funding coupon to be called */
  int num_fund_cpn; /*	Number of funding coupons called (including redemption)
                     */
  int num_partial_fund_cpn;

  /*	Fee upon exercise */
  long set_date;    /*	Fee payment date */
  double set_time;  /*	Fee payment time */
  double fee;       /*	Amount of the fee */
  double extra_fee; /*	Extra fee to adjust forward IV */
  double total_fee; /*	= fee + extra_fee */

  double implied_strike;
  int is_midat;

} cts_call, *CTS_CALL;

/*	CTS */
typedef struct {
  cts_exo_leg *exo_leg;
  cts_fund_leg *fund_leg;
  int num_calls;
  cts_call *call;
} cts, *CTS;

Err cts_fill_calls(
    /*	Exercises before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag,           /*	0: I      , 1: E */
    int ncall, int pay_rec, /*	1: rec pd      , -1: pay pd */
    long *ex_date, long *set_date, double *fee, CTS cts);

/* Adjustment for fixing in floating case */
Err cts_adjust_floating_fixing(CTS cts, int pde_or_mc);

/*	Check dates consistency */
Err cts_check_calls(CTS cts);

/*	Free */
void cts_free_calls(CTS cts);

/* Find equivalent cash strikes */
Err cts_calc_equivalent_strikes(CTS cts_iv, CTS cts_call, CTS_MKT mkt);

/*	Structures and functions for the underlying and its term structures */

/*	Underlying term structs */
typedef struct {
  char name[256];
  CTS_MKT mkt;
  LGMSV_model model;

} cts_und, *CTS_UND, lgmsv_und, *LGMSV_UND;

/*	Fill underlying structure from a predefined underlying */
Err cts_fill_und(CTS_MKT mkt, char *lgmsvund, double tstar, CTS_UND und);

Err cts_implied_alpha_approx(double maturity,
                             double param[], /* first param = Alpha^2      ,
                                                second param = 2.0 * LamEps */
                             double *price, double *gradient, int nb_param);

Err cts_get_mkt_sabr_beta_param(CTS_MKT mkt, long exe_date, long start_date,
                                long theo_end_date, double fixed_beta,
                                int use_levenberg, double *alpha, double *rho,
                                double *fitting_error);

/*	Calibration of alpha      , rho and lameps to the market smile */
Err cts_calib_const_alpha_and_rho(
    CTS_MKT mkt, CTS cts,
    int calib_months, /* 0: co-terminal swaption      , otherwise underlyings
                     with required nb months */
    double fixed_beta, double min_time, double *alpha, double *rho,
    double *lameps, int save_inst_data, cpd_calib_inst_data *inst_data);

/*	Fill underlying structure from calibration instruments */
Err cts_calib_und(
    CTS_MKT mkt,
    /*	EOD Flag */
    int eod_flag,  /*	0: I      , 1: E */
    double lambda, /*	lambda */
    int nb_factor, double lgm_alpha, double lgm_gamma, double lgm_rho,
    int nsmilepar, double *smilepartime, double *alpha, /*	alpha */
    double *rho,                                        /*	rho */
    double *ldaeps,                                     /*	ldaeps */
    double *rho2,                                       /*	rho2eps */
    double tstar,
    /*	Numerical CF params */
    LGMSV_NUMERPARAMS NumerParams,
    /*	Calib params */
    char *cal_tenor, char *cal_ref, char *cal_freq, char *cal_basis,
    int force_atm, /*	force atm calib */
    double max_std_long, double max_std_short, double vol_shift_long,
    DIAGCALIB_VOLTYPE vol_type_long, DIAGCALIB_SHIFTTYPE vol_shift_type_long,
    double vol_shift_short, DIAGCALIB_VOLTYPE vol_type_short,
    DIAGCALIB_SHIFTTYPE vol_shift_type_short, double lambda_shift,

    int calib_strategy, /*	-1: autocal      , 0: swaptions / cap      , 1:
                       cap / swaptions */
    int fix_lambda,     /*	0: calib lambda to cap      , 1: fix lambda calib to
                       diagonal */
    char *short_tenor, char *short_refrate, char *short_freq, char *short_basis,
    int fix_smile,          /*	0: calib smile parameters to market smile */
    int smile_calib_months, /* 0: co-terminal swaption      , otherwise
                           underlyings with required nb months */
    LGMSV_CalibParams *calib_params, double min_time,
    int skip_last,   /*	If 1      , the last option is disregarded and the
                    forward   volatility is flat from option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps,   /*	1: we allow jumps on vol      , 0: we don't */
    double prec, int maxiter, int keep_first,
    /*	Strike choice */
    int long_strike_flag,  /*	0: ATM
                                               1: Coupon
                                               2: Eq (PV/Lvl) */
    int short_strike_flag, /*	0: ATM      ,
                                               1: implied digital caplet
                          strike 2: same number of std */
    /*	End of calib params */
    CTS cts, /*	structure */
    CTS_UND und, int save_inst_data, cpd_calib_inst_data *inst_data);

Err cts_copy_und(CTS_UND src, CTS_UND dest);

void cts_free_und(CTS_UND und);

/*	IV info for feedback */
typedef struct {
  int ncall;
  long *call_date;
  double *market_fwdiv;
  double *model_fwdiv;
  double *extra_fee;

  int has_one_time;
  double *one_time_call;
  double *one_time_strike;
  double *one_time_norm_vol;
  double *one_time_log_vol;
  double *one_time_coupon;
  double *one_time_notional;

} cts_iv, *CTS_IV;

/*	Initialise the IV infos */
void cts_init_cts_iv(CTS_IV iv_infos);

/* Free the IV infos */
void cts_free_cts_iv(CTS_IV iv_infos);

/*	IV info for feedback */
typedef struct {
  int nrow;
  int ncol;
  double **extra_infos;
  double call1;
  double call2;

  double algo_fwd_iv1;
  double algo_fwd_iv2;

} cts_extra_infos, *CTS_EXTRA_INFOS;

/* Free the Extra infos */
void cts_free_extra_infos(CTS_EXTRA_INFOS extra_infos);

/*	Calculate forward and market IVs */
Err cts_calc_model_fwd_iv(CTS_UND und, CTS cts, int calc_fwd_iv, int adj_fee,
                          /*	Numerical CF params */
                          LGMSV_NUMERPARAMS NumerParams);

/*	Adjust model to match market IVs */
Err cts_adjust_model_fwd_iv(CTS_UND und, CTS market_cts, CTS model_cts,
                            int for_fund, int calc_fwd_iv, int adj_fee,
                            int pde_or_mc,
                            /*	Feedback */
                            int save_fwdiv, CTS_IV fwd_iv_info);

/*	Constants used for reconstruction and evaluation  */
/*	------------------------------------------------- */

typedef struct {
  /*	DF (Tpay) / DF (Tstar) */
  /*	Reconstruction is exp ( - alpha - beta * f(t      ,T*) - gamma * Psit )
   */
  double tpay_tstar_alpha;
  double tpay_tstar_beta;
  double tpay_tstar_gamma;

  /* 2 Factor extra */
  double tpay_tstar_beta2;
  double tpay_tstar_gamma2;
  double tpay_tstar_gamma12;

  /*	DF (Ts) / DF (Te) */
  /*	Reconstruction is exp ( - alpha - beta * f(t      ,T*) - gamma * Psit )
   */
  double ts_te_alpha;
  double ts_te_beta;
  double ts_te_gamma;

  /* 2 Factor extra */
  double ts_te_beta2;
  double ts_te_gamma2;
  double ts_te_gamma12;

  /*	Special flag for initialisation */
  int do_init_col2;

  /*	Special option precalculations */
  /*	Alternative notation */
  /*	a + b * X + string of options on X */
  /*	X = 1 + cvg * cash_libor */
  double a;
  double b;
  double s[CTS_MAX_STR];
  double n[CTS_MAX_STR];

} cts_eval_const_fix, *CTS_EVAL_CONST_FIX;

typedef struct {
  /*	DF (Tpay) / DF (Tstar) */
  /*	Reconstruction is exp ( - alpha - beta * f(t      ,T*) - gamma * Psit )
   */
  double tpay_tstar_alpha;
  double tpay_tstar_beta;
  double tpay_tstar_gamma;

  /*	DF (Ti) / DF (Tstar) */
  /*	Reconstruction is exp ( - alpha - beta * f(t      ,T*) - gamma * Psit )
   */
  double ti_tstar_alpha[MAX_CPN];
  double ti_tstar_beta[MAX_CPN];
  double ti_tstar_gamma[MAX_CPN];

  /* 2 Factor extra */
  double tpay_tstar_beta2;
  double tpay_tstar_gamma2;
  double tpay_tstar_gamma12;

  double ti_tstar_beta2[MAX_CPN];
  double ti_tstar_gamma2[MAX_CPN];
  double ti_tstar_gamma12[MAX_CPN];

  /*	Special flag for initialisation */
  int do_init_col2;

  /*	Special option precalculations */
  /*	Alternative notation */
  /*	a + b * X + string of options on X */
  /*	X = 1 + cvg * cash_libor */
  double a;
  double b;
  double s[CTS_MAX_STR];
  double n[CTS_MAX_STR];

} cts_eval_const_fix_cms, *CTS_EVAL_CONST_FIX_CMS;

typedef struct {
  /*	DF (Tpay) / DF (Tstar) */
  /*	Reconstruction is exp ( - alpha - beta * f(t      ,T*) - gamma * Psit )
   */
  double tpay_tstar_alpha;
  double tpay_tstar_beta;
  double tpay_tstar_gamma;

  /*	DF (Ts) / DF (Te) */
  /*	Reconstruction is exp ( - alpha - beta * f(t      ,T*) - gamma * Psit )
   */
  double ts_te_alpha;
  double ts_te_beta;
  double ts_te_gamma;

  /* 2 Factor extra */
  double tpay_tstar_beta2;
  double tpay_tstar_gamma2;
  double tpay_tstar_gamma12;

  double ts_te_beta2;
  double ts_te_gamma2;
  double ts_te_gamma12;

  /*	Reconstruction of the geared Libor */
  double shift;  /* gearing * (-1/cvg + spread) + margin */
  double multip; /* gearing / cvg */

} cts_eval_const_cpn, *CTS_EVAL_CONST_CPN;

typedef struct {
  /*	DF (Tfund) / DF (Tstar) */
  /*	Reconstruction is exp ( - alpha - beta * f(t      ,T*) - gamma * Psit )
   */
  double tpay_tstar_alpha[MAXDF];
  double tpay_tstar_beta[MAXDF];
  double tpay_tstar_gamma[MAXDF];

  /*	DF (Tset) / DF (Tstar) */
  /*	Reconstruction is exp ( - alpha - beta * f(t      ,T*) - gamma * Psit )
   */
  double tset_tstar_alpha;
  double tset_tstar_beta;
  double tset_tstar_gamma;

  /*	DF (Tpayfixed) / DF (Tstar) */
  /*	Reconstruction is exp ( - alpha - beta * f(t      ,T*) - gamma * Psit )
   */
  double tpay_tstar_alpha2[MAXDF];
  double tpay_tstar_beta2[MAXDF];
  double tpay_tstar_gamma2[MAXDF];

  /* 2 Factor extra */
  double tpay_tstar_beta_2[MAXDF];
  double tpay_tstar_gamma_2[MAXDF];
  double tpay_tstar_gamma_12[MAXDF];

  double tset_tstar_beta2;
  double tset_tstar_gamma2;
  double tset_tstar_gamma12;

  double tpay_tstar_beta2_2[MAXDF];
  double tpay_tstar_gamma2_2[MAXDF];
  double tpay_tstar_gamma2_12[MAXDF];

  /* flag for 1-time callable valuation */
  int eval_one_time;

} cts_eval_const_call, *CTS_EVAL_CONST_CALL;

/*	Fill fixing structure */
Err cts_fill_eval_const_fix(CTS_UND und, CTS cts,
                            /*	Index of the current cpn */
                            int cpn_idx,
                            /*	Index of the current fixing */
                            int fix_idx, CTS_EVAL_CONST_FIX eval_const_fix);

/*	Fill fixing structure in the case of a CMS */
Err cts_fill_eval_const_fix_cms(CTS_UND und, CTS cts,
                                /*	Index of the current cpn */
                                int cpn_idx,
                                /*	Index of the current fixing */
                                int fix_idx,
                                CTS_EVAL_CONST_FIX_CMS eval_const_fix_cms);

/*	Fill coupon structure */
Err cts_fill_eval_const_cpn(CTS_UND und, CTS cts,
                            /*	Index of the current cpn */
                            int cpn_idx, CTS_EVAL_CONST_CPN eval_const_cpn);

/*	Fill call structure */
Err cts_fill_eval_const_call(CTS_UND und, CTS cts,
                             /*	Index of the current call */
                             int call_idx, CTS_EVAL_CONST_CALL eval_const_call);

/*	Arguments to all payoff evaluation functions */
/*	-------------------------------------------- */

typedef struct {
  /*	Underlying */
  CTS_UND und;

  /*	Structure */
  CTS cts;

  /*	Payoff type */
  int is_fixing;
  int is_fixing_cms;
  int is_cpn;
  int is_call;

  /*	Index of the current fixing      , coupon and call */
  CTS_FIX fix;
  CTS_CPN cpn;
  CTS_CALL call;

  /*	Eval consts */
  cts_eval_const_fix eval_const_fix;
  cts_eval_const_fix_cms eval_const_fix_cms;
  cts_eval_const_cpn eval_const_cpn;
  cts_eval_const_call eval_const_call;

  /* Extra path dependent infos */
  double *dPathInfos;

} cts_pay_arg, *CTS_PAY_ARG;

/*	Arguments to the numerical functions */
/*	------------------------------------ */

typedef struct {
  /* Algo Choice */
  int pde_or_mc;

  int nstp;
  double *time;
  double *date;
  int initnstp;
  int nstppsi;
  int nstpx;
  int nstpz;
  double integ_mintime;
  int has_one_time;
  void **void_prm;
  int *is_event;
  double *ifr;
  LGMSVPARAM params;

  /* For MC */
  int nb_event;
  double mc_mintime;
  long nsteps;
  long npaths;

  MCEBPARAMS mcebparams;
  int *optimise;

  /* For 1F MC */
  double *dSigma;
  double *dAlpha;
  double *dLambdaEps;
  double *dLvlEps;
  double *dRho;

  /* For 2F MC */
  double *dRho2;
  double *dLGMAlpha;
  double *dLGMRho;

  /* For Floating MC */
  double *dPathInfos;

} cts_adi_arg, *CTS_ADI_ARG;

Err cts_find_one_time_index(CTS_UND und, CTS cts, int *one_time_index);

Err cts_fill_algo_arg(
    CTS_UND und, CTS cts, int pde_or_mc, /* 0: PDE      , 1: MC */
    /*	Required number of steps for PDE */
    int req_stp, int req_stppsi, int req_stpx, int req_stpz,
    /*	Required number of steps for MC */
    double req_mintime, long req_paths,
    /*	Extra numerical parameter */
    double integ_mintime,
    /*	Flag for extra calculation / adjustment of one-time callable */
    int do_one_time,    /*	1: calc the one time */
    int one_time_index, /*	0: choose automatically the index      , >0:
                       index provided by user */
    /*	T-star */
    CTS_ADI_ARG adi_arg);

/* Free */
void cts_free_adi_arg(CTS_ADI_ARG adi_arg);

/*	Find Numer TStar */
void cts_get_numer_tstar(long today, CTS cts, double *numer_tstar);

/*	Main function to be called in order to fill and check all structures */
/*	==================================================================== */

/*	Fill and check all the relevant structures */
Err cts_fill_check_all_struct(
    /*	Market */
    CTS_MKT mkt,
    /*	The underlying */
    double tstar, int use_calib, /*	0: use lgmsvund      , 1: calibrate */
    /*		if calib */
    double lambda, /*	LGM lambda */
    int nb_factor, double lgm_alpha, double lgm_gamma, double lgm_rho,
    int nsmilepar, double *smilepartime, double *alphaeps, /*	alpha */
    double *rhoeps,                                        /*	rho */
    double *ldaeps,                                        /*	ldaeps */
    double *rho2eps,
    /*	End of calib params */
    char *lgmsvund,
    /*	The structure */
    /*		funding */
    int fund_ncpn, double *fund_not, long *fund_fix, long *fund_start,
    long *fund_end, long *fund_pay, char **fund_basis, double *fund_spr,
    double *fund_mrg,
    /*		ts */
    /*			general Libor properties */
    int fix_lag_bd,
    /*			number of coupons */
    int ncpn,
    /*			coupon description */
    /*				cpn */
    int *cpn_type, /*	0:	fixed
                               1:	libor fixed at start
                               2:	libor fixed at end */
    long *cpn_start_date, long *cpn_end_date, long *cpn_pay_date,
    double *cpn_coupon, char **cpn_basis, double *cpn_not,
    /*				pay libor */
    long *pay_fix_date, /*	fixing date of the Libor */
    int pay_months, char *pay_freq, char *pay_basis, double *pay_fwd_spread,
    double *pay_margin,
    /*				fix libor */
    /*				fixings */
    int *nfix, double **weights, long **fix_dates,
    char **fix_tenor, /*	Fixing tenors */
    char *fix_freq, char *fix_basis, double **fix_fwd_spreads,
    /*				profiles */
    int tstype, /*	0: generic      , 1: range */
    /*				profile 0 */
    double *alpha, double *beta, int *nstr, double **str, double **nbopt,
    /*				profile 1 */
    int buy_sell,   /*	1: BNPP buys      , -1: BNPP sells */
    int value_zero, /*	0: Do not value low strike options      , 1: do */
    double *lb,     /*	Lower bounds of ranges */
    double *ub,     /*	Upper bounds of ranges */
    double *payoff, /*	Payoff in ranges */
    double call_spread, double numer_spread,
    /*		Trimming */
    int trim_type, /*	0: no trim
                               1: x fixings max
                               2: x time min between two fixings */
    int max_fix, double min_fix_time,
    /*	Extra model parameters*/
    int use_cmsopt,      /*	Use CmsOption to value fix option      , use BS on
                        CmsRate      otherwise */
    double correl_start, /*	Correl between libor fixing and libor paid at
                        start */
    double
        correl_end, /*	Correl between libor fixing and libor paid at start */
    int float_adjust_type, /*	type of adjustment for the floating coupon , 0:
                              ATM vol      , 1: Strike Vol
                        */
    /*	Calls */
    int ncall, int pay_rec, /*	1: rec pd      , -1: pay pd */
    long *ex_date, long *set_date, double *fee, int adj_fee,

    /*	Flag for extra calculation / adjustment of one-time callable */
    int do_one_time,    /*	1: calc the one time */
    int one_time_index, /*	0: choose automatically the index      , >0:
                       index provided by user */
    int just_recalib,   /*	Just recalibrate the underlying and skip all the
                       other parts */

    /*	Numerical params */
    /*		CF */
    LGMSV_NUMERPARAMS NumerParams,
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
    int calib_strategy, /*	-1: autocal      , 0: swaptions / cap      , 1:
                       cap / swaptions */
    int fix_lambda,     /*	0: calib lambda to cap      , 1: fix lambda calib
                                    to diagonal */
    char *short_tenor, char *short_refrate, char *short_freq, char *short_basis,
    int fix_smile,          /*	0: calib smile parameters to market smile */
    int smile_calib_months, /* 0: co-terminal swaption      , otherwise
                           underlyings with required nb months */
    LGMSV_CalibParams *calib_params, double min_time,
    int skip_last,   /*	If 1      , the last option is disregarded and the
                    forward   volatility is flat from option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps,   /*	1: we allow jumps on vol      , 0: we don't */
    double numer_tstar, double prec, int maxiter, int keep_first,
    /*	Strike choice */
    int long_strike_flag,  /*	0: ATM
                                               1: Coupon
                                               2: Eq (PV/Lvl) */
    int short_strike_flag, /*	0: ATM      ,
                                               1: implied digital caplet
                          strike 2: same number of std */
    CTS cts_iv,            /*	Needed to compute implied strike */
    /*	Flags */
    /*		IV calculation */
    int calc_mkt_iv, int calc_fwd_iv,
    /*		EOD */
    int eod_fix_flag, /*	0: I      , 1: E */
    int eod_ex_flag,  /*	0: I      , 1: E */
    /*	Results */
    CTS cts, CTS_UND und, int *call_feat, /*	0: No callable feature to be
                                 valued 1: Callable feature to be
                                 valued through adi */
    CTS_ADI_ARG adi_arg,
    /*	Feedback */
    int save_inst_data, cpd_calib_inst_data *inst_data, int save_fwdiv,
    CTS_IV fwd_iv_info);

/*	Free all structures */
void cts_free_all_struct(CTS cts, CTS_UND und,
                         int call_feat, /*	0: No callable feature to be
                                       valued 1: Callable feature to be
                                       valued through adi */
                         CTS_ADI_ARG adi_arg);

/*	Payoff function for adi (callable) */
/*	-----------------------------------	*/

Err cts_payoff_4_lgmsv_adi(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    void *yc, double lam, double tstar, double alpha, double rho, double sigma,
    /* Nodes data */
    int lpsi, int upsi, int lx, int ux, int lz, int uz, double *psi, double *x,
    double *z, int *nprod,
    /* Vector of results to be updated */
    double ****prod_val);

/*	Payoff function for MC (callable with optimise KO) */
/*	-------------------------------------------------- */

Err cts_payoff_4_lgmsv_mc(
    /* Event */
    long path_index, double evt_date, double evt_time, void *func_parm,
    double ft, double psi, double v, int nprod,
    /* Vector of results to be updated */
    double *prod_val, int *stop_path);

Err cts_payoff_4_lgmsv2F_mc(
    /* Event */
    long path_index, double evt_date, double evt_time, void *func_parm,
    double ft1, double ft2, double psi1, double psi2, double psi12, double v,
    int nprod,
    /* Vector of results to be updated */
    double *prod_val, int *stop_path);

/*	Main pricing function for adi */

/*	Launch the adi */
Err cts_launch_algo(CTS cts, CTS_UND und, CTS_ADI_ARG adi_arg,

                    /*	Result */
                    double *multi_pv, double *onetime_pv, double *fwd_iv);

#endif