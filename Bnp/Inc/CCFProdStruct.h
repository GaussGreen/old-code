
#ifndef __CCF_PROD_STRUCT_H
#define __CCF_PROD_STRUCT_H

#define CCF_NCPN 512
#define CCF_NOPT 10

/*	Structures and functions for callable cms floaters */
/*	-------------------------------------------------- */

/*	Structures and functions for the funding leg */

/*	Funding cpn */
typedef struct {
  long start_date;   /*	Coupon start date */
  double start_time; /*	Coupon start time */
  long pay_date;     /*	Coupon pay date */
  double pay_time;   /*	Coupon pay time */
  double cvg;        /*	Cvg */
  double cpn;        /*	Notional * (fwd_spr + margin) * cvg */
} cf_fund_cpn, *CF_FUND_CPN;

/*	Funding leg */
typedef struct {
  double notional; /*	Notional */
  int num_cpn;     /*	Number of coupons */
  /*	0..num_cpn-1 */
  cf_fund_cpn *cpn;
} cf_fund_leg, *CF_FUND_LEG;

Err ccf_fill_fund_leg(
    /*	Coupons that started before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I      , 1: E */
    double fund_not, int fund_ncpn, long *fund_fix, long *fund_start,
    long *fund_pay, char **fund_basis, double *fund_spr, double *fund_mrg,
    CF_FUND_LEG fund_leg);

/*	Check dates consistency */
Err ccf_check_fund_leg(CF_FUND_LEG fund_leg);

/*	Free */
Err ccf_free_fund_leg(CF_FUND_LEG fund_leg);

/*	Structures and functions for the exotic leg */

/*	CMS description */
typedef struct {
  /*	Underlying swap */
  long start_date;
  double start_time;
  long end_date;
  double end_time;
  /*	Coupons */
  int num_cpn;
  long cpn_pay_date[CCF_NCPN];
  double cpn_pay_time[CCF_NCPN];
  double cvg[CCF_NCPN];

  /*	For CMS calculation	*/
  int cpnd;
  double nfp;
  double delay;
  double fix_time;

  /*	Forward (cash) */
  double fwd;
  double fwd_cms;

  /*	Spread from a cash swap */
  double swap_spread;

  /*	Spread from libor */
  double lib_spread;

  /*	Normal variance */
  double atmvar;   /*	ATM: for CMS adjsutment */
  double floorvar; /*	Floor: for floor */
  double capvar;   /*	Cap: for cap */

  double optvar[CCF_NOPT]; /* opt normal */

  /* Shifted Log parameter */
  double slbeta;
  double slshift;
  double slvol;

  /*	CMS Adjustment */
  int adj_type; /*	0: no adjustment      , 1: atm normal adjustment      ,
                   2: smile adjustment */
  /*	If smile is used      , CMS = F+lambda*F^(2*beta)*mat
          beta is input      , lambda = (CMS-F)/(F^(2*beta)*mat)

          if shifted log is used      , if smile is used      , ignores the
     input beta , CMS = F + lambda*(F+slshift)^2*mat
          (calibrated to hit today's adjustment) */
  double lambda;
  double beta;

} cf_cms_desc, *CF_CMS_DESC;

/*	CMS Float cpn */
typedef struct {
  long start_date;     /*	Coupon start date */
  double start_time;   /*	Coupon start time */
  long pay_date;       /*	Coupon pay date */
  double pay_time;     /*	Coupon pay time */
  long cms_fix_date;   /*	CMS fixing dates for coupon */
  double cms_fix_time; /*	CMS fixing times for coupon */

  double cvg; /*	Coupon coverage * notional */

  int ncms; /*	0 for type 0      , 1 for types 1 and 2      , 2 for type 3 */
  cf_cms_desc cms[2];

  /*	Coupon type
  0:	Midat			(only gamma is used)
  1:	CIF				(only alpha and cms1 are used)
  2:	CMS				(only alpha and cms1 are used)
  3:	CMS Spread */
  int type;

  /*	Coupon is of the form alpha * CMS1 + beta * CMS2 + gamma [floored      ,
   * capped] */
  double alphabeta[2]; /*	0: alpha      , 1: beta */
  double gamma;
  int floored;
  double floor;
  int capped;
  double cap;

  int use_SL; /* 1: use shifted log model to value the options
                         0: don't
               */

  int use_cfoptions;          /* 1: use the string of options  0: dont */
  int nopt;                   /* Number of spread options */
  double notopt[CCF_NOPT];    /* Notional of the spread options */
  double strikeopt[CCF_NOPT]; /* spread option strikes */
  int typeopt[CCF_NOPT];      /* spread option type 0 call 1 put */

} cf_exo_cpn, *CF_EXO_CPN;

/*	CMS Float leg */
typedef struct {
  double notional;

  /*	Leg type
  0:	Midat			(all coupons are Midat)
  1:	CIF				(all coupons are Midat and at least 1 is
  CIF)
  2:	CMS				(all coupons are Midat or CIF and at
  least 1 is CMS) 3:	CMS Spread		(at least 1 coupon is CMS
  spread) */
  int type;

  int num_cpn; /*	Number of coupons */
  /*	0..num_cpn-1 */
  cf_exo_cpn *cpn;

} cf_exo_leg, *CF_EXO_LEG;

/*	ATTENTION: cms lambdas are not filled and must be handled separetely */
Err ccf_fill_exo_leg(
    /*	Coupons that fixed before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I      , 1: E */
    double cf_not, int cf_ncpn, long *cf_fix, long *cf_start, long *cf_pay,
    char **cf_basis, char **cf_cms_tenor1, char *cf_cms_freq1,
    char *cf_cms_basis1, double *cf_cms_spread1, char **cf_cms_tenor2,
    char *cf_cms_freq2, char *cf_cms_basis2, double *cf_cms_spread2,
    double *cf_alpha, double *cf_beta, double *cf_gamma, int *cf_floored,
    double *cf_floor, int *cf_capped, double *cf_cap,

    long spread_vol_n, double *spread_vol_time,
    double *spread_slbeta1, /* Shifted log beta on the CMS1 */
    double *spread_slbeta2, /* Shifted log beta on the CMS2 */

    int cf_nopt,           /* Number of spread options */
    double **cf_notopt,    /* Notional of the spread options */
    double **cf_strikeopt, /* spread option strikes */
    int **cf_typeopt,      /* spread option type 0 call 1 put */

    int use_SL,   /*  1: use Shifted Log
                  0: don't use				*/
    int calib_SL, /*  1: Calibrate the Shifted Log models on CMS1 and on CMS2 *
                              0: use the Shifted Log beta given */
    double NbStdcalib_SL, /*  Nb of std for the calibration of to the skew */
    int calib_correl_SL,  /*	1: Calibrate the correlation between the two SL
                         to get the same ATM  normal spread vol  0: use the
                         normal spread correl for the sl correl */
    int use_cfoptions, /*  1: Use the spread options and don't take into account
                                           the floor and cap in the cf
                      coupons
                                   0: take into account the floor and cap in
                      the cf coupons */

    int cms_adj, int cms_smile, int cms_vol_adj, /*	1: adjust for CMS vol
                                        effect 0: don't */
    double cms_beta1, double cms_beta2,
    int num_strikes_in_vol, /*	Array of strikes in vol matrix */
    double *strikes_in_vol,
    SrtDiffusionType vol_type, /*	Type of vol in matrix      , SRT_NORMAL
                                  or SRT_LOGNORMAL */
    int cash_vol,              /*	1: matrix is a cash vol
     /*	Needed to calculate ATM (converging) normal volatility associated to
         CMS coupons */
    char *yc,                  /*	yc */
    char *vc,                  /*	vc */
    char *ref,                 /*	ref rate */
    char *swap_freq,           /*	swap freq */
    char *swap_basis,          /*	swap basis */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    CF_EXO_LEG exo_leg);

/*	Check dates consistency */
Err ccf_check_exo_leg(CF_EXO_LEG exo_leg);

/*	Free */
Err ccf_free_exo_leg(CF_EXO_LEG exo_leg);

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
  int cf_idx;       /*	Index of the first cf coupon to be called */
  int num_cf_cpn;   /*	Number of cf coupons called (excluding redemption) */
  /*	Fee upon exercise */
  long set_date;   /*	Fee payment date */
  double set_time; /*	Fee payment time */
  double fee;      /*	Amount of the fee */
} cf_call, *CF_CALL;

/*	Callable CMS Floater */
typedef struct {
  cf_exo_leg *cf_leg;
  cf_fund_leg *fund_leg;
  int num_calls;
  cf_call *call;
} ccf_str, *CCF_STR;

Err ccf_fill_calls(
    /*	Exercises before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag,           /*	0: I      , 1: E */
    int ncall, int pay_rec, /*	0: rec pd      , 1: pay pd */
    long *ex_date, long *set_date, double *fee, CCF_STR ccf);

/*	Check dates consistency */
Err ccf_check_calls(CCF_STR ccf);

/*	Free */
Err ccf_free_calls(CCF_STR ccf);

/*	Structures and functions for the underlying and its term structures */

/*	Underlying term structs */
typedef struct {
  char name[256];
  long today;
  char yc[256];
  char vc[256];
  char ref[256]; /*	Reference rate used for getting the vol */
  char swap_freq[256];
  char swap_basis[256];
  long spread_vol_n;
  double *spread_vol_time;
  double *spread_vol_floor;
  double *spread_vol_cap;
  int cvg_sv;
  int is_corr;
  long sigma_n;
  double *sigma_date;
  double *sigma_time;
  double *sigma; /*	Sigma1 */
  long lambda_n;
  double *lambda_date;
  double *lambda_time;
  double *lambda;
  double alpha; /*	Sigma2 = alpha * Sigma1 */
  double gamma; /*	Lambda2 = Lambda1 + gamma */
  double rho;

  /*	Calibration instrument data */
  int has_inst_data;
  cpd_calib_inst_data inst_data;

  /*	Value of the fwd receiver CCF */
  int has_fwd_iv;
  int nb_fwdiv;
  double *exercise_date;
  double *market_fwdiv;
  double *model_fwdiv;
  double *extra_fees;

} ccf_und, lgm2f_und, *CCF_UND, *LGM2F_UND;

/*	Fill underlying structure from a predefined underlying */
Err ccf_fill_und(char *lgm2dund, char *vc, char *ref, char *swap_freq,
                 char *swap_basis, long spread_vol_n, double *spread_vol_time,
                 double *spread_vol_floor, double *spread_vol_cap, int cvg_sv,
                 int is_corr, CCF_UND und);

/*	Fill underlying structure from calibration instruments */
Err ccf_calib_und(
    long today,
    /*	EOD Flag */
    int eod_flag,      /*	0: I      , 1: E */
    char *yc,          /*	yc */
    char *vc,          /*	vc (only if calib) */
    char *ref,         /*	ref rate (only if calib) */
    char *swap_freq,   /*	swap freq (only if calib) */
    char *swap_basis,  /*	swap basis (only if calib) */
    int lam_ts,        /*	0: use unique lambda      , 1: use ts */
    double lambda,     /*	lambda if unique */
    int tsnlam,        /*	number of lambdas if ts */
    double *tslamtime, /*	lambda times i.e. (date - today) / 365 if ts */
    double *tslam,     /*	corresponding lambdas if ts */
    double alpha,      /*	alpha */
    double gamma,      /*	gamma */
    double rho,        /*	rho */
    /*	Calib params */
    int force_atm, /*	force atm calib */
    double max_std_long, double max_std_short,
    int fix_lambda,         /*	0: calib lambda to cap      , 1: fix lambda calib
                                                        to diagonal */
    int cal_vol_shift_type, /*	vol shift type for volatility */
    double cal_vol_shift,   /*	vol shift */
    double cal_lambda_shift, /*	shift on lambda after calibration */
    int one_f_equi,          /*	1F equivalent flag:
                                                         if set to 1      , then 2F
                            lambda will calibrate          to the cap priced within
                            calibrated          1F          with the given lambda */
    int skip_last,     /*	If 1      , the last option is disregarded and the
                      forward     volatility is flat from option n-1 */
    double long_prec,  /*	Precision on primary instruments */
    double short_prec, /*	Precision on secondary instruments */
    double min_fact,   /*	Maximum down jump on variance */
    double max_fact,   /*	Maximum up jump on variance */
    int use_jumps,     /*	Allow vol term structure to jump */
    int proba_weight, int use_exe_bound, double *proba, double *exe_bound,
    /*	End of calib params */
    CCF_STR ccf,        /*	structure */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    CCF_UND und, long spread_vol_n, double *spread_vol_time,
    double *spread_vol_floor, double *spread_vol_cap, int cvg_sv, int is_corr);

void ccf_copy_und(CCF_UND src, CCF_UND dest);

Err ccf_free_und(CCF_UND und);

/*	Constants used for reconstruction and evaluation  */
/*	------------------------------------------------- */

#define MAXDF 1000
typedef struct {
  /*	Discount Factors */
  /*	Reconstruction is exp ( - alpha - beta * r1 - gamma * r2 ) */

  /*	All df */
  int num_df;             /*	Number of df required */
  double df_mat[MAXDF];   /*	Residual maturity of df */
  double df_alpha[MAXDF]; /*	For reconstruction */
  double df_beta[MAXDF];
  double df_gamma[MAXDF];
  double df[MAXDF]; /*	Space to store the df */

  /*	Funding */
  int do_fund;
  int num_fund_cpn;
  int start_idx;
  int fund_idx[CCF_NCPN];

  /*	CF discount */
  int do_cf_disc;
  int num_cf_cpn;
  int cf_disc_idx[CCF_NCPN];

  /*	CF coupons */
  int do_cf_fwd;
  int cf_fwd_start_idx[2][CCF_NCPN];
  int cf_fwd_idx[2][CCF_NCPN][CCF_NCPN];

  /*	Fee */
  int fee_idx;

  /*	Volatility */

  /*	CMS effect */
  double cf_cms_vol[2][CCF_NCPN];

  /*	Spread options */
  int cf_spr_floor_type[CCF_NCPN];
  double cf_spr_floor_str[CCF_NCPN];
  double cf_spr_floor_std[CCF_NCPN];

  int cf_spr_cap_type[CCF_NCPN];
  double cf_spr_cap_str[CCF_NCPN];
  double cf_spr_cap_std[CCF_NCPN];

  double cf_spr_opt_std[CCF_NCPN]
                       [CCF_NOPT]; /* Normal variance of the spread option */
  double cf_spr_slcorr[CCF_NCPN];  /* shifted log correlation */

} ccf_eval_const, *CCF_EVAL_CONST;

/*	Fill structure */
Err ccf_fill_eval_const(CCF_UND und, CCF_STR ccf,
                        /*	Index of the current call */
                        int call_idx, CCF_EVAL_CONST eval_const);

/*	Arguments to all payoff evaluation functions */
/*	-------------------------------------------- */

typedef struct {
  /*	Underlying */
  ccf_und *und;

  /*	Structure */
  ccf_str *ccf;

  /*	Constants used for reconstruction and evaluation */
  ccf_eval_const eval_const;

  /*	Index of the current call */
  int call_idx;

} ccf_pay_arg, *CCF_PAY_ARG;

/*	Arguments to the adi function */
/*	----------------------------- */

typedef struct {
  int nstp;
  double *time;
  double *date;
  int nstpx;
  double *sig_time;
  double *sig1;
  int nb_sig;
  double *lam_time;
  double *lam;
  int nb_lam;
  double alpha;
  double gamma;
  double rho;
  void **void_prm;
  int *is_event;
  double *ifr;
  char yc[256];
} ccf_adi_arg, *CCF_ADI_ARG;

Err ccf_fill_adi_arg(
    CCF_UND und, CCF_STR ccf,
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Required number of steps */
    int req_stp, int req_stpx, CCF_ADI_ARG adi_arg);

Err ccf_free_adi_arg(CCF_ADI_ARG adi_arg);

typedef struct {
  long npaths;
  double *time;
  double *date;
  int nb_dates;
  int jumping_num;
  long pay_date;
  double pay_time;
  double *dom_fwd1;
  double *dom_fwd2;
  double *dom_exp1;
  double *dom_exp2;
  double *dom_phi1;
  double *dom_phi2;
  double *dom_phi12;
  double *dom_gam1_fwd;
  double *dom_gam2_fwd;
  double *dom_bond_pay;
  double *dom_gam1_pay;
  double *dom_gam2_pay;
  double ***covar;

  void **void_prm;

} ccf_mc_arg, *CCF_MC_ARG;

Err ccf_fill_mc_arg(
    CCF_UND und, CCF_STR ccf,
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Required number of steps */
    long req_paths, int jumping_num, CCF_MC_ARG mc_arg);

Err ccf_free_mc_arg(CCF_MC_ARG mc_arg);

/*	Main function to be called in order to fill and check all structures */
/*	==================================================================== */

/*	Fill and check all the relevant structures */
Err ccf_fill_check_all_struct(
    /*	Today's date */
    long today,
    /*	The underlying */
    int use_calib, /*	0: use fx3dund      , 1: calibrate */
    /*		if calib */
    char *yc,          /*	yc */
    char *vc,          /*	vc (only if calib) */
    char *ref,         /*	ref rate (only if calib) */
    char *swap_freq,   /*	swap freq (only if calib) */
    char *swap_basis,  /*	swap basis (only if calib) */
    int lam_ts,        /*	0: use unique lambda      , 1: use ts */
    double lambda,     /*	lambda if unique */
    int tsnlam,        /*	number of lambdas if ts */
    double *tslamtime, /*	lambda times i.e. (date - today) / 365 if ts */
    double *tslam,     /*	corresponding lambdas if ts */
    double alpha,      /*	alpha */
    double gamma,      /*	gamma */
    double rho,        /*	rho */
    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*		if no calilb */
    char *lgm2dund,
    /*	The structure */
    /*		funding */
    double fund_not, int fund_ncpn, long *fund_fix, long *fund_start,
    long *fund_pay, char **fund_basis, double *fund_spr, double *fund_mrg,
    /*		cf */
    double cf_not, int cf_ncpn, long *cf_fix, long *cf_start, long *cf_pay,
    char **cf_basis, char **cf_cms_tenor1, char *cf_cms_freq1,
    char *cf_cms_basis1, double *cf_cms_spread1, char **cf_cms_tenor2,
    char *cf_cms_freq2, char *cf_cms_basis2, double *cf_cms_spread2,
    long spread_vol_n, double *spread_vol_time, double *spread_vol_floor,
    double *spread_vol_cap,

    double *spread_slbeta1, /* Shifted log beta on the CMS1 */
    double *spread_slbeta2, /* Shifted log beta on the CMS2 */

    int cvg_sv, int is_corr, double *cf_alpha, double *cf_beta,
    double *cf_gamma, int *cf_floored, double *cf_floor, int *cf_capped,
    double *cf_cap,

    int cf_nopt,           /* Number of spread options */
    double **cf_notopt,    /* Notional of the spread options */
    double **cf_strikeopt, /* spread option strikes */
    int **cf_typeopt,      /* spread option type 0 call 1 put */

    int use_SL,   /*  1: use Shifted Log
                  0: don't use				*/
    int calib_SL, /*  1: Calibrate the Shifted Log models on CMS1 and on CMS2 *
                              0: use the Shifted Log beta given */
    double NbStdcalib_SL, /*  Nb of std for the calibration of to the skew */
    int calib_correl_SL,  /*	1: Calibrate the correlation between the two SL
                         to get the same ATM  normal spread vol  0: use the
                         normal spread correl for the sl correl */
    int use_cfoptions, /*  1: Use the spread options and don't take into account
                                           the floor and cap in the cf
                      coupons
                                   0: take into account the floor and cap in
                      the cf coupons */

    int cms_adj, int cms_smile, int cms_vol_adj, /*	1: adjust for CMS vol
                                        effect 0: don't */
    double cms_beta1, double cms_beta2,
    int num_strikes_in_vol, /*	Array of strikes in vol matrix */
    double *strikes_in_vol,
    SrtDiffusionType vol_type, /*	Type of vol in matrix      , SRT_NORMAL
                                  or SRT_LOGNORMAL */
    int cash_vol,              /*	1: matrix is a cash vol
         /*		calls */
    int ncall, int pay_rec,    /*	0: rec pd      , 1: pay pd */
    long *ex_date, long *set_date, double *fee,
    /*	Numerical params */
    int req_stp, int req_stpx, long req_paths,
    /*	Calib params */
    int force_atm, double max_std_long, double max_std_short,
    int fix_lambda,         /*	0: calib lambda to cap      , 1: fix lambda calib
                                                                to diagonal */
    int cal_vol_shift_type, /*	vol shift type for volatility */
    double cal_vol_shift,   /*	vol shift */
    double cal_lambda_shift, /*	shift on the calibrated lambda */
    int one_f_equi,          /*	1F equivalent flag:
                                                         if set to 1      , then 2F
                            lambda will calibrate          to the cap priced within
                            calibrated          1F          with the given lambda */
    int skip_last,     /*	If 1      , the last option is disregarded and the
                      forward     volatility is flat from option n-1 */
    double long_prec,  /*	Precision on primary instruments */
    double short_prec, /*	Precision on secondary instruments */
    double min_fact,   /*	Maximum down jump on variance */
    double max_fact,   /*	Maximum up jump on variance */
    int use_jumps,     /*	Allow vol term structure to jump */
    int proba_weight,
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I      , 1: E */
    int eod_ex_flag,  /*	0: I      , 1: E */
    /*	Results */
    CCF_STR ccf, CCF_UND und,
    int *call_feat, /*	0: No callable feature to be valued
                                1: Callable feature to be valued through adi
                 */
    CCF_ADI_ARG adi_arg);

/*	Free all structures */
Err ccf_free_all_struct(CCF_STR ccf, CCF_UND und, int call_feat,
                        CCF_ADI_ARG adi_arg);

/*	Payoff function for adi (callable) */
/*	-----------------------------------	*/

Err ccf_payoff_4_3dfx_adi(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    void *yc, double *lam, double *ts_time, int nb_ts, double gamma, double rho,
    double phi1, double phi2, double phi12,
    /* Nodes data */
    int l1, int u1, int l2, int u2, double *r1, double **r2, int nprod,
    /* Vector of results to be updated */
    double ***prod_val);

/*	Payoff function for adi (callable) with shifted log or string of options
 */
/*	-------------------------------------------------------------------------*/
Err ccf_payoff_4_3dfx_adi_sl(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    void *yc, double *lam, double *ts_time, int nb_ts, double gamma, double rho,
    double phi1, double phi2, double phi12,
    /* Nodes data */
    int l1, int u1, int l2, int u2, double *r1, double **r2, int nprod,
    /* Vector of results to be updated */
    double ***prod_val);

/*	Main pricing function for adi */

/*	Launch the adi */
Err ccf_launch_adi(CCF_STR ccf, CCF_UND und, CCF_ADI_ARG adi_arg,
                   /*	Result */
                   double *prem);

/*	Launch the MC */
Err ccf_launch_mc(CCF_STR ccf, CCF_UND und, CCF_MC_ARG mc_arg,
                  /*	Result */
                  double *prem, double *probas, double *exe_bound);

#endif