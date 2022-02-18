
#ifndef Fx3FCalib_h
#define Fx3FCalib_h

#include "uterror.h"
#include "srt_h_und_struct.h"

/* Integral of exp(-x*(T-u)) between s and t */
double Phi_Func(double x, double T, double s, double t);

/* Integral of (1.0 - exp(-x*(T-u))) / x between s and t */
double Etha_Func(double x, double T, double s, double t);

/* Integral of (1.0 - exp(-x*(T-u))) / x * (1.0 - exp(-y*(T-u))) / x between s and t */
double Psi_Func(double x, double y, double T, double s, double t);

/* Integral of (1.0 - exp(-x*(T-u))) / x * exp(-y*(T-u)) between s and t */
double Gamma_Func(double x, double y, double T, double s, double t);

/* Integral of (1.0 - exp(-x*(Tx-u))) / x * (1.0 - exp(-y*(Ty-u))) / y between s and t */
double Psi2_Func(double x, double y, double Tx, double Ty, double s, double t);

/* Integral of exp(-x(T-u)) * (1.0 - exp(-y*(T-u))) / y between s and t */
double Zeta_Func(double x, double y, double T, double s, double t);

/* Integral of (1.0 - exp(-x*(Tx-u))) / x * exp(-y*(Ty-u)) between s and t */
double Gamma2_Func(double x, double y, double Tx, double Ty, double s, double t);

/* Defines the long term Fx volatility using the 3F model
                                5 strategies
                                1.)	Calibrate all
                                2.)	Calibrate all, optimise fvol stationarity in slr
                                3.) Calibrate partial, extrapolate with slr
                                4.) Calibrate partial, best fit rest with slr
                                5.) Optimise fvol stationarity under bid/offer constraint */

Err Fx3DDefImpVol(/*	Total number of options */
                  int num_opt,
                  /*	Number of options to calibrate to */
                  int num_calib,
                  /*	Implied volatility */
                  double* opt_mat,
                  double* ivol,
                  double* bid,
                  double* offer,
                  /*	Model parameters and bounds */
                  double* fvol,
                  double* sig_dom,
                  double  sig_dom_min,
                  double  sig_dom_max,
                  double* sig_for,
                  double  sig_for_min,
                  double  sig_for_max,
                  double* lam_dom,
                  double  lam_dom_min,
                  double  lam_dom_max,
                  double* lam_for,
                  double  lam_for_min,
                  double  lam_for_max,
                  double* rho_dom_for,
                  double  rho_dom_for_min,
                  double  rho_dom_for_max,
                  double* rho_dom_fx,
                  double  rho_dom_fx_min,
                  double  rho_dom_fx_max,
                  double* rho_for_fx,
                  double  rho_for_fx_min,
                  double  rho_for_fx_max,
                  int     ext,
                  double  exta,
                  double  extb,
                  /*	Number of iterations */
                  int num_iter,
                  /*	What to optimise on */
                  int opt_sigma,
                  int opt_lambda,
                  int opt_correl,
                  int force_eq_sig,
                  int force_eq_lam,
                  /*	Method 1 to 5 */
                  int method);

Err Fx3DDefImpVol2(/*	Total number of options */
                   int num_opt,

                   /*	Number of options to calibrate to */
                   int num_calib,

                   /*	Implied volatility */
                   double* opt_exe,
                   double* opt_mat,
                   double* ivol,

                   /*	Market parameters */
                   long  today,
                   long  spot_date_dom,
                   long  spot_date_for,
                   long  maturity_date,
                   char* yc_dom,
                   char* vol_dom,
                   char* ref_rate_dom,
                   char* basis_dom,
                   char* cpd_dom,
                   char* ccy_dom,
                   char* yc_for,
                   char* vol_for,
                   char* ref_rate_for,
                   char* basis_for,
                   char* cpd_for,
                   char* ccy_for,

                   /* Model parameters */
                   double lam_dom,
                   double lam_for,
                   double rho_dom_for,
                   double rho_dom_fx,
                   double rho_for_fx,

                   /* Extrapolation parameters */
                   double vol_lim,
                   double vol_lam,

                   /* Calibration parameters */
                   char*  call_freq,
                   double dom_vol_shift,
                   double for_vol_shift,

                   /* Vol Function */
                   Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                                       char*   vol_curve_name,
                                       double  start_date,
                                       double  end_date,
                                       double  cash_strike,
                                       int     zero,
                                       char*   ref_rate_name,
                                       double* vol,
                                       double* power),

                   /* Answer */
                   double* vol);

Err Fx3DDefImpVol2_corr(/*	Total number of options */
                        int num_opt,

                        /*	Number of options to calibrate to */
                        int num_calib,

                        /*	Implied volatility */
                        double* opt_exe,
                        double* opt_mat,
                        double* ivol,

                        /*	Market parameters */
                        long  today,
                        long  spot_date_dom,
                        long  spot_date_for,
                        long  maturity_date,
                        char* yc_dom,
                        char* vol_dom,
                        char* ref_rate_dom,
                        char* basis_dom,
                        char* cpd_dom,
                        char* ccy_dom,
                        char* yc_for,
                        char* vol_for,
                        char* ref_rate_for,
                        char* basis_for,
                        char* cpd_for,
                        char* ccy_for,

                        /* Model parameters */
                        double  lam_dom,
                        double  lam_for,
                        double* rho_times,
                        double* rho_dom_for_ts,
                        double* rho_dom_fx_ts,
                        double* rho_for_fx_ts,
                        long    nb_rho,

                        /* Extrapolation parameters */
                        double vol_lim,
                        double vol_lam,

                        /* Calibration parameters */
                        char*  call_freq,
                        double dom_vol_shift,
                        double for_vol_shift,

                        /* Vol Function */
                        Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                                            char*   vol_curve_name,
                                            double  start_date,
                                            double  end_date,
                                            double  cash_strike,
                                            int     zero,
                                            char*   ref_rate_name,
                                            double* vol,
                                            double* power),

                        /* Answer */
                        double* vol);

/* Fx3D functions that can be directly used from the structure */

/*	Check that there is no term structure of tau and get lambda */
Err get_unique_lambda(double* tau, int ntau, double* lda);
/*	Merge rates term structures */
Err merge_rates_ts(
    double*  sig_dom_mat,
    double*  sig_dom,
    int      sig_n_dom,
    double*  sig_for_mat,
    double*  sig_for,
    int      sig_n_for,
    double** sig_merge_mat,
    double** sig_merge_dom,
    double** sig_merge_for,
    int*     sig_merge_n);

/*	Merge rates, fx and corr term structures */
Err merge_rates_fx_corr_ts(
    double*  sig_dom_mat,
    double*  sig_dom,
    int      sig_n_dom,
    double*  sig_for_mat,
    double*  sig_for,
    int      sig_n_for,
    double*  sig_fx_mat,
    double*  sig_fx,
    int      sig_n_fx,
    double*  corr_mat,
    double*  corr_dom_for,
    double*  corr_dom_fx,
    double*  corr_for_fx,
    int      corr_n_mat,
    double** sig_merge_mat,
    double** sig_merge_dom,
    double** sig_merge_for,
    double** sig_merge_fx,
    double** corr_merge_dom_for,
    double** corr_merge_dom_fx,
    double** corr_merge_for_fx,
    int*     sig_merge_n);

/*	IR Phi */
Err Fx3DtsPhi(
    double maturity, double* maturity_sig, long nbSig, double* sig_curve, double lda, double* phi);

/*	Fx Implied vol */
Err Fx3DtsImpliedVol(
    double  opt_maturity,
    double  start_date,
    double  end_date,
    double* maturity_rates,
    long    nbMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* maturity_fx,
    double* sig_curve_fx,
    long    nbrMat_fx,
    double  correl_dom_for,
    double  correl_dom_fx,
    double  correl_for_fx,
    double* fx_vol);

Err Fx3DtsImpliedVolExtrapol(
    double  opt_maturity,
    double  start_date,
    double  end_date,
    double* maturity_rates,
    long    nbMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* maturity_fx,
    double* sig_curve_fx,
    long    nbrMat_fx,
    double  sig_inf,
    double  lda_vol,
    double  correl_dom_for,
    double  correl_dom_fx,
    double  correl_for_fx,
    double* fx_vol);

Err Fx3DtsImpliedVolExtrapol_corr(
    double  opt_maturity,
    double  start_date,
    double  end_date,
    double* maturity_rates,
    long    nbMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* maturity_fx,
    double* sig_curve_fx,
    long    nbrMat_fx,
    double  sig_inf,
    double  lda_vol,
    double* correl_times,
    double* correl_dom_for_ts,
    double* correl_dom_fx_ts,
    double* correl_for_fx_ts,
    long    nb_correl,
    double* fx_vol);

/*	Calibration */
Err Fx3DtsCalibration(
    double*  exercise_opt,
    double*  maturity_opt,
    double*  vol_opt,
    long     nbrOpt,
    double*  maturity_rates,
    long     nbrMat,
    double*  sig_curve_dom,
    double   lda_dom,
    double*  sig_curve_for,
    double   lda_for,
    double   correl_dom_for,
    double   correl_dom_fx,
    double   correl_for_fx,
    double** fx_vol_curve);

/*	Forward */
Err Fx3DtsFwdFx(
    SrtUndPtr dom_und, SrtUndPtr for_und, SrtUndPtr fx_und, long maturity_date, double* fwd_fx);

/* Fx3D functions that can be used with the name of the underlying Fx. Used in exportation to Xl */

/*	Calibration */
Err Fx3DCalibration(
    char*    dom_underlying,
    char*    for_underlying,
    double   correl_dom_for,
    double   correl_dom_fx,
    double   correl_for_fx,
    double*  exercise_opt,
    double*  maturity_opt,
    double*  vol_opt,
    long     nbropt,
    double** fx_vol_curve);

/*	Forward */
Err Fx3DFwdFx(char* fx_underlying, long maturity_date, double* fwd_fx);

/*	Implied vol */
Err Fx3DImpliedVol(
    char* fx_underlying, double val_time, double start_time, double end_time, double* vol);

/*	Option price */
Err Fx3DFxOption(
    char*   fx_underlying,
    double  strike,
    long    fix_date,
    long    val_date,
    long    pay_date,
    char*   pay_reveceive,
    double* price);

/*	Forward vol of spot fx */
Err Fx3DSpotVol(char* fx_underlying, double maturity, double* fx_vol);

/*	Compute the expectation of log ( S (T) / S (0) )
                under a variety of measures */

/*	Calculate adjustment between Q-Tpay and Q-beta, such that
        Q-beta expect log ( S ( Tfix ) / S ( 0 ) ) = Q-Tpay expect log ( S ( Tfix ) / S ( 0 ) ) +
   adj */
Err Fx3DtsFwdBetaAdjustment(
    double T0,   /*	Forward start date */
    double Tval, /*	Value date of the forward */
    double Tpay, /*	Pay date of the forward */
    double Tfix, /*	Fix date of the forward */
    /*	Model data */
    double* maturity_rates,
    long    nbrMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* maturity_fx,
    double* sig_curve_fx,
    long    nbrMat_fx,
    double  correl_dom_for,
    double  correl_dom_fx,
    double  correl_for_fx,
    /*	Result */
    double* adjust);

/*	Calculate adjustment between Q-Tpay1 and Q-Tpay2, such that
        Q-Tpay2 expect log ( S ( Tfix ) / S ( 0 ) ) = Q-Tpay1 expect log ( S ( Tfix ) / S ( 0 ) ) +
   adj */
Err Fx3DtsFwdPayAdjustment(
    double T0,    /*	Forward start date */
    double Tval,  /*	Value date of the forward */
    double Tpay1, /*	Original pay date */
    double Tpay2, /*	Adjusted pay date */
    double Tfix,  /*	Fix date of the forward */
    /*	Model data */
    double* maturity_rates,
    long    nbrMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* maturity_fx,
    double* sig_curve_fx,
    long    nbrMat_fx,
    double  correl_dom_for,
    double  correl_dom_fx,
    double  correl_for_fx,
    /*	Result */
    double* adjust);

/*	Calculate cumulative covariance of F ( mat1 ) and F ( mat2 ) between 0 and T */
Err Fx3DtsFwdCumCovar(
    double T0,    /*	Forward start date */
    double Tval1, /*	Value date of the 1st forward */
    double Tval2, /*	Value date of the 2nd forward */
    double Tfix,  /*	Fix date of both forwards */
    /*	Model data */
    double* maturity_rates,
    long    nbrMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* maturity_fx,
    double* sig_curve_fx,
    long    nbrMat_fx,
    double  correl_dom_for,
    double  correl_dom_fx,
    double  correl_for_fx,
    /*	Result */
    double* adjust);

/*	Q beta */
Err Fx3DFwdBeta(
    char* fx_underlying,
    /*	Forward fixing time */
    double Tfix,
    /*	Forward value time */
    double Tval,
    /*	0: expect [log Fx], 1: expect [Fx] */
    int log_flag,
    /*	Answer */
    double* expect);

/*	Q Tpay */
Err Fx3DFwdTpay(
    char* fx_underlying,
    /*	Forward fixing time */
    double Tfix,
    /*	Forward value time */
    double Tval,
    /*	Payment time */
    double Tpay,
    /*	0: expect [log Fx], 1: expect [Fx] */
    int log_flag,
    /*	Answer */
    double* expect);

/*	Q F */
Err Fx3DFwdQf(
    char* fx_underlying,
    /*	F1 fixing time */
    double Tbarfix,
    /*	F1 value time */
    double Tbarval,
    /*	F2 fixing time */
    double Tfix,
    /*	F2 value time */
    double Tval,
    /*	F2 Payment time */
    double Tpay,
    /*	0: expect [log Fx], 1: expect [Fx] */
    int log_flag,
    /*	Answer */
    double* expect);

/*	Compute covar ( log ( S ( T1 ) ) , log ( S ( T2 ) ) ) */
Err Fx3DCovar(
    char* fx_underlying, double Tfix1, double Tval1, double Tfix2, double Tval2, double* covar);

/* covariance between two fx underlyings		*/
/* all the term structures have to be merged... */
Err Fx3DtsFxFwdCov(
    double   option_maturity,
    double   start_date,
    double   end_date,
    double*  maturity,
    long     nbMat,
    double*  sig_curve_dom1,
    double   lda_dom1,
    double*  sig_curve_for1,
    double   lda_for1,
    double*  sig_curve_fx1,
    double*  sig_curve_dom2,
    double   lda_dom2,
    double*  sig_curve_for2,
    double   lda_for2,
    double*  sig_curve_fx2,
    double** correlations,
    double*  covar);

/* calculates any correlation within the 3F framework */
/* can be fx/fx fx/lgm lgm/lgm */
/* put NULL vector when not needed */
/* all the term structures have to be merged... */

Err Fx3DtsFxFwdCov_corr(
    double  option_maturity1,
    double  option_maturity2,
    double  start_date,
    double  end_date,
    double* maturity,
    long    nbMat,
    /* FX1 or LGM1 */
    double* sig_curve_dom1,
    double  lda_dom1,
    double* sig_curve_for1,
    double  lda_for1,
    double* sig_curve_fx1,
    /* FX2 or LGM2 */
    double*   sig_curve_dom2,
    double    lda_dom2,
    double*   sig_curve_for2,
    double    lda_for2,
    double*   sig_curve_fx2,
    double*** correlations,
    double*   var1,
    double*   var2,
    double*   covar);

/* Reading term structures from underlyings */
Err Get_LGM_TermStructure(
    char*    underlying,
    double** sigma_time,
    double** sigma,
    long*    sigma_n,
    double** tau_time,
    double** tau,
    long*    tau_n);

Err Get_LGM_TermStructure2(
    char* underlying, double** sigma_time, double** sigma, long* sigma_n, double* fixed_tau);

Err Get_LGM2F_TermStructure(
    char*    und_name,
    double** sigma_time,
    double** sigma,
    long*    sigma_n,
    double*  fixed_tau,
    double*  fixed_alpha,
    double*  fixed_gamma,
    double*  fixed_rho);

Err Get_LGM2F_TermStructure_ts(
    char*    underlying,
    double** sigma_time,
    double** sigma,
    long*    sigma_n,
    double*  fixed_tau,
    double*  fixed_alpha,
    double*  fixed_gamma,
    double*  fixed_rho,
    int*     tau_ts);

Err Get_LGM2F_TermStructure2(
    char*    underlying,
    double** sigma,
    double** sigma_time,
    long*    nb_sigma,
    double** lambda,
    double** lambda_time,
    long*    nb_lambda,
    double*  fixed_alpha,
    double*  fixed_gamma,
    double*  fixed_rho);

Err LGM2FDtsPhi(
    double  maturity,
    double* maturity_sig,
    long    nbSig,
    double* sig_curve,
    double* maturity_lam,
    long    nbLam,
    double* lam_curve,
    double  alpha,
    double  gamma,
    double  rho,
    double* phi1,
    double* phi2,
    double* phi12);

Err Get_FX_StochRate_TermStructures(
    char*    underlying,
    double** sigma_date_dom,
    double** sigma_dom,
    long*    sigma_n_dom,
    double** tau_date_dom,
    double** tau_dom,
    long*    tau_n_dom,
    double** sigma_date_for,
    double** sigma_for,
    long*    sigma_n_for,
    double** tau_date_for,
    double** tau_for,
    long*    tau_n_for,
    double** sigma_date_fx,
    double** sigma_fx,
    long*    sigma_n_fx,
    double*  rho_df,
    double*  rho_dfx,
    double*  rho_ffx);

Err display_FX1F_TermStruct(char* szFxUndName, long* sigma_n, double** sigma_date, double** sigma);

Err Fx3DDefImpVolBeta(/*	Total number of options */
                      int num_opt,

                      /*	Number of options to calibrate to */
                      int num_calib,

                      /*	Implied volatility */
                      double* opt_exe,
                      double* opt_mat,
                      double* ivol,

                      /*	Market parameters */
                      long   today,
                      long   spot_date_dom,
                      long   spot_date_for,
                      long   maturity_date,
                      double spot_fx,
                      char*  dom_yc,
                      char*  vol_dom,
                      char*  ref_rate_dom,
                      char*  basis_dom,
                      char*  cpd_dom,
                      char*  ccy_dom,
                      char*  for_yc,
                      char*  vol_for,
                      char*  ref_rate_for,
                      char*  basis_for,
                      char*  cpd_for,
                      char*  ccy_for,

                      /* Model parameters */
                      double lam_dom,
                      double lam_for,
                      double rho_dom_for,
                      double rho_dom_fx,
                      double rho_for_fx,
                      double beta,

                      /* Extrapolation parameters */
                      double vol_lim,
                      double vol_lam,

                      /* Calibration parameters */
                      char*  call_freq,
                      double dom_vol_shift,
                      double for_vol_shift,

                      long   nbSteps,
                      double disc_dt,
                      double fx_dt,
                      long   nbIterMax,

                      /* Vol Function */
                      Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                                          char*   vol_curve_name,
                                          double  start_date,
                                          double  end_date,
                                          double  cash_strike,
                                          int     zero,
                                          char*   ref_rate_name,
                                          double* vol,
                                          double* power),

                      /* Answer */
                      double* vol);

/*	Term correlation between a domestic IR index and the forward Fx */
Err quanto_correl(
    long   fwd_start_date, /*	Fwd start date */
    long   fix_date,       /*	Fixing time */
    long   start_date,     /*	Index start date */
    long   end_date,       /*	Index end date */
    long   pay_date,       /*	Payment time */
    double fx_ivol_pay,    /*	Implied vol of Fx for payment time */
    char*  dom_yc_name,    /*	Name of the domestic yield curve */
    char*  dom_vc_name,    /*	Name of the domestic market vol curve */
    char*  dom_ref_name,   /*	Name of the domestic reference rate */
    char*  dom_swap_freq,  /*	Frequency and basis of underlying swaptions */
    char*  dom_swap_basis,
    char*  for_yc_name,   /*	Name of the foreign yield curve */
    char*  for_vc_name,   /*	Name of the foreign market vol curve */
    char*  for_ref_name,  /*	Name of the foreign reference rate */
    char*  for_swap_freq, /*	Frequency and basis of underlying swaptions */
    char*  for_swap_basis,
    Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    /*	Lambdas */
    double dom_lam,
    double for_lam,
    /*	The 3 standard correlation inputs */
    double corr_dom_for,
    double corr_spot_fx_dom,
    double corr_spot_fx_for,
    /*	The "exotic" correl term structure, between domestic yield and domestic index */
    int     num_corr,
    double* corr_mats,
    double* corr_dom_yld_idx,
    /*	The output */
    double* the_corr,
    double* quanto_adj);

/*	Term correlation between a domestic IR index and the forward Fx */
Err quanto_correl_corr(
    long   fwd_start_date, /*	Fwd start date */
    long   fix_date,       /*	Fixing time */
    long   start_date,     /*	Index start date */
    long   end_date,       /*	Index end date */
    long   pay_date,       /*	Payment time */
    double fx_ivol_pay,    /*	Implied vol of Fx for payment time */
    char*  dom_yc_name,    /*	Name of the domestic yield curve */
    char*  dom_vc_name,    /*	Name of the domestic market vol curve */
    char*  dom_ref_name,   /*	Name of the domestic reference rate */
    char*  dom_swap_freq,  /*	Frequency and basis of underlying swaptions */
    char*  dom_swap_basis,
    char*  for_yc_name,   /*	Name of the foreign yield curve */
    char*  for_vc_name,   /*	Name of the foreign market vol curve */
    char*  for_ref_name,  /*	Name of the foreign reference rate */
    char*  for_swap_freq, /*	Frequency and basis of underlying swaptions */
    char*  for_swap_basis,
    Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    /*	Lambdas */
    double dom_lam,
    double for_lam,
    /*	The 3 standard correlation inputs */
    double* corr_times,
    double* corr_dom_for_ts,
    double* corr_dom_fx_ts,
    double* corr_for_fx_ts,
    long    corr_n_times,
    /*	The "exotic" correl term structure, between domestic yield and domestic index */
    int     num_corr,
    double* corr_mats,
    double* corr_dom_yld_idx,
    /*	The output */
    double* the_corr,
    double* quanto_adj);

/*	**********************************************************	*/
/*	All this part is a second version of the functions above	*/
/*	with a term-structure of coorelation						*/
/*	**********************************************************	*/
/*	**********************************************************	*/

Err Get_FX_StochRate_TermStructures_corr(
    char*    underlying,
    double** sigma_date_dom,
    double** sigma_dom,
    long*    sigma_n_dom,
    double** tau_date_dom,
    double** tau_dom,
    long*    tau_n_dom,
    double** sigma_date_for,
    double** sigma_for,
    long*    sigma_n_for,
    double** tau_date_for,
    double** tau_for,
    long*    tau_n_for,
    double** sigma_date_fx,
    double** sigma_fx,
    long*    sigma_n_fx,
    double** correl_date,
    double** correl_dom_for,
    double** correl_dom_fx,
    double** correl_for_fx,
    long*    correl_n);

/*	Fx Implied volatility */
/*  The rates maturity dates have to be merged ! */
Err Fx3DtsImpliedVol_corr(
    double  opt_maturity,
    double  start_date,
    double  end_date,
    double* maturity_rates,
    long    nbMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* maturity_fx,
    double* sig_curve_fx,
    long    nbrMat_fx,
    double* maturity_corr,
    double* correl_dom_for,
    double* correl_dom_fx,
    double* correl_for_fx,
    long    nbrMat_corr,
    double* fx_vol);

/*	Calibration of a fx term structure to a set of fx options  */
Err Fx3DtsCalibration_corr(
    double*  exercise_opt,
    double*  maturity_opt,
    double*  vol_opt,
    long     nbrOpt,
    double*  maturity_rates,
    long     nbrMat,
    double*  sig_curve_dom,
    double   lda_dom,
    double*  sig_curve_for,
    double   lda_for,
    double*  maturity_corr,
    double*  correl_dom_for,
    double*  correl_dom_fx,
    double*  correl_for_fx,
    long     nbrCorr,
    double** fx_vol_curve);

/*	Implied vol direct from underlying */
Err Fx3DImpliedVol_corr(
    char* fx_underlying, double val_time, double start_time, double end_time, double* vol);

/*	Fx calibration direct from underlying */
Err Fx3DCalibration_corr(
    char*    dom_underlying,
    char*    for_underlying,
    double*  maturity_correl,
    double*  correl_dom_for,
    double*  correl_dom_fx,
    double*  correl_for_fx,
    long     nb_correl,
    double*  exercise_opt,
    double*  maturity_opt,
    double*  vol_opt,
    long     nbropt,
    double** fx_vol_curve);

Err Fx3DtsFwdBetaAdjustment_corr(
    double T0,   /*	Forward start date */
    double Tval, /*	Value date of the forward */
    double Tpay, /*	Pay date of the forward */
    double Tfix, /*	Fix date of the forward */
    /*	Model data */
    double* maturity_rates,
    long    nbrMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* maturity_fx,
    double* sig_curve_fx,
    long    nbrMat_fx,
    double* maturity_corr,
    double* correl_dom_for,
    double* correl_dom_fx,
    double* correl_for_fx,
    long    nbrMat_corr,
    /*	Result */
    double* adjust);

/*	Q beta */
Err Fx3DFwdBeta_corr(
    char* fx_underlying,
    /*	Forward fixing time */
    double Tfix,
    /*	Forward value time */
    double Tval,
    /*	0: expect [log Fx], 1: expect [Fx] */
    int log_flag,
    /*	Answer */
    double* expect);

/*	Calculate adjustment between Q-Tpay1 and Q-Tpay2, such that
        Q-Tpay2 expect log ( S ( Tfix ) / S ( 0 ) ) = Q-Tpay1 expect log ( S ( Tfix ) / S ( 0 ) ) +
   adj */
Err Fx3DtsFwdPayAdjustment_corr(
    double T0,    /*	Forward start date */
    double Tval,  /*	Value date of the forward */
    double Tpay1, /*	Original pay date */
    double Tpay2, /*	Adjusted pay date */
    double Tfix,  /*	Fix date of the forward */
    /*	Model data */
    double* maturity_rates,
    long    nbrMat,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* maturity_fx,
    double* sig_curve_fx,
    long    nbrMat_fx,
    double* maturity_corr,
    double* correl_dom_for,
    double* correl_dom_fx,
    double* correl_for_fx,
    long    nbrMat_corr,
    /*	Result */
    double* adjust);

/*	Implied correlations */
Err Fx3DtsCorrel(
    double  start_mat,
    double  end_mat,
    double  val_mat,
    double* sig_dates,
    long    nb_sig_dates,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* sig_curve_fx,
    double* correl_dom_for_ts,
    double* correl_dom_fx_ts,
    double* correl_for_fx_ts,
    double* dom_vol,
    double* for_vol,
    double* fx_vol,
    double* dom_for_cor,
    double* dom_fx_cor,
    double* for_fx_cor);

/*	Implied correlations */
Err Fx3DCorrel(
    char*   fx_underlying,
    double  val_time,
    double  start_time,
    double  end_time,
    double* dom_vol,
    double* for_vol,
    double* fx_vol,
    double* dom_for_cor,
    double* dom_fx_cor,
    double* for_fx_cor);

Err srt_f_optrainbow3F(
    double  T,
    double  K,
    double  nX,
    double  nY,
    char*   fx_und1,
    char*   fx_und2,
    double* res,
    double* correl);

Err Fx3DFwdTpay_corr(
    char* fx_underlying, double Tfix, double Tval, double Tpay, int log_flag, double* expect);

Err Fx3DGenCorrelation(
    char*   und1_name,
    char*   und2_name,
    double  val_time,
    double  start_time,
    double  end_time,
    double* vol1,
    double* vol2,
    double* correl);

#endif