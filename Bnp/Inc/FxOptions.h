
#ifndef FxOPTIONSH
#define FxOPTIONSH

#include "CPDProdstruct.h"

Err ConvexFx1F(char *underlying, long maturity, double barrier, double strike,
               double typeBar, /*	0: minimum < Bar   1: maximum > bar */
               long npaths, long nstp, double *res1, double *error1,
               double *res2, double *error2);

Err SABR_to_Strikes(double Fwd, double T, double DfDom, double DfFor,
                    double Sigma, double Alpha, double Beta, double Rho,
                    SrtDiffusionType VolType, int DeltaType, int IsForeign,
                    double Precision, double *Strike1, double *Strike2,
                    double *Strike3, double *Vol1, double *Vol2, double *Vol3,
                    double *VolATM);

Err SABR_SolveDeltaNewton(double F, double T, double DfDom, double DfFor,
                          double Delta, int IsCall, double Sigma, double Alpha,
                          double Beta, double Rho, SrtDiffusionType VolType,
                          double (*GetDelta)(double F, double T, double DfDom,
                                             double DfFor, double K, double Vol,
                                             int IsCall),
                          int IsForeign, double Precision, double *K,
                          double *Sig);

Err SABR_SolveStraddleNewton(
    double F, double T, double DfDom, double DfFor, double Sigma, double Alpha,
    double Beta, double Rho, SrtDiffusionType VolType,
    double (*GetDelta)(double F, double T, double DfDom, double DfFor, double K,
                       double Vol, int IsCall),
    int IsForeign, double Precision, double *K, double *Sig);

Err SolveDeltaNewton(double F, double T, double DfDom, double DfFor,
                     double Delta, int IsCall, double Sigma,
                     double (*GetDelta)(double F, double T, double DfDom,
                                        double DfFor, double K, double Vol,
                                        int IsCall),
                     int IsForeign, double Precision, double *K);

Err SolveStraddleNewton(double F, double T, double DfDom, double DfFor,
                        double Sigma,
                        double (*GetDelta)(double F, double T, double DfDom,
                                           double DfFor, double K, double Vol,
                                           int IsCall),
                        int IsForeign, double Precision, double *K);

double DeltaSoho(double Fwd, double T, double DfDom, double DfFor, double K,
                 double Vol, int IsCall, int DeltaType);

Err Strikes_to_SABR(double Fwd, double T, double DfDom, double DfFor,
                    double Vol25Put, double VolStraddle, double Vol25Call,
                    double *Strike25Put, double *StrikeStraddle,
                    double *Strike25Call, double *volATM, double *Alpha,
                    double *Beta, double *Rho, int FreezeAlpha, int FreezeBeta,
                    int FreezeRho, int DeltaType, int IsForeign,
                    double Precision, double *fitting_error);

/*	Barrier Autocal model params */
typedef struct {
  long today;
  double spot_fx;
  double *sig;
  double *sig_time;
  long nb_sig;

  double alpha;
  double beta;
  double rho;
  double lambda;
  double cvxty;
  double floormu;

} barautocal_und, *BARAUTOCAL_UND;

void barautocal_free_und(BARAUTOCAL_UND und);

Err BarrierOption_autocal(
    /*	Market parameters */
    long today, char *dom_yc, char *for_yc,
    double spot_fx, /* as quoted in the market */

    /*	Product parameters */
    double notional, long exercise_date, long settlmt_date, double strike,
    double barrier_up, double barrier_down, double rebate_up,
    double rebate_down, int is_call, /* 1: Call      , 0: Put */
    int is_ko,                       /* 1: KO      , 0: KI */
    int is_cvx,                      /* 1: use 1/Fx      , 0: use Fx */
    int is_digital,  /* 1: digital payoff      , 0      , regular option payoff
                      */
    int is_american, /* 1: American      , 0: European */

    /*	Model parameters */
    double mdl_alpha, double mdl_beta, double mdl_rho, double mdl_lambda,
    double mdl_gamma, double mdl_floor,

    /*	Calib parameters */
    long *opt_matdates,           /* mkt options maturity */
    double *opt_atmvols,          /* mkt options ATM vols */
    int nb_opt, int calib_all_ts, /* 1: calib all the TS      , 0: calib only
                                     the relevant */

    int do_smile_calib, /* 1: calibrates alpha and rho      , 0: uses model
                           input */
    int cal_stochvol,   /* 1: calib alpha and rho      , 0: calib beta and cvxty
                         */
    long cal_maturity,  /* option maturity      , exercise is 2bd before */
    double cal_alpha,   /* SABR alpha */
    double cal_beta,    /* SABR beta */
    double cal_rho,     /* SABR rho */

    /*	Numerical parameters */
    int pr_nstept,   /* time steps for pricing */
    int pr_nstepfx,  /* Fx steps for pricing */
    int pr_nstepvol, /* Vol steps for pricing */

    int cal_nstept,       /* time steps for calibration */
    int cal_nstepfx,      /* Fx steps for calibration */
    int cal_nstepvol,     /* Vol steps for calibration */
    int cal_atmmaxiter,   /* Maximum number of iterations for ATM calibration */
    double cal_atmprec,   /* Precision on ATM vols */
    int cal_smilemaxiter, /* Maximum number of iterations for ATM calibration */
    double cal_smileprec, /* Precision on ATM vols */
    int cal_usetotalts, /* 1: when calibrating smile      , consider total TS ,
                       0: consider only the corresponding ATM vol */

    /* IOD / EOD flags */
    int eod_fix_flag, /*	EOD Fixing Flag 0: I      , 1: E */
    int eod_pay_flag, /*	EOD Payment Flag 0: I      , 1: E */
    int eod_ex_flag,  /*	EOD Exercise Flag 0: I      , 1: E */

    /*	Exercised flag */
    int exercised,         /*	Is exercised Flag */
    int knocked,           /*	Has knocked flag */
    int has_knocked_up,    /*	Knocked on the up barrier:1      , on the down
                          one:
                          0 */
    long ex_date_ex,       /*	Date when exercised */
    long ex_date_set,      /*	Corresponding settlement date */
    double spot_fx_fixing, /*	Fixing of the Fx */

    /*	Outputs */
    double *price, double *greeks, int export_ts, int use_old_ts,
    BARAUTOCAL_UND calib_und);

/* Price in the simple lognormal model (1 / FX(T) - K)+ * 1{Hd < FX(T) < Hu} */
Err ConvexFwdSimple(
    /*	Market parameters */
    long today, char *dom_yc, char *for_yc,
    double spot_fx, /* as quoted in the market */

    /*	Product parameters */
    double notional, long exercise_date, long settlmt_date, double strike,
    int is_call, double barrier_up, double barrier_down,

    /*	Calib parameters */
    long *opt_matdates,  /* mkt options maturity */
    double *opt_atmvols, /* mkt options ATM vols */
    int nb_opt,

    /* IOD / EOD flags */
    int eod_fix_flag,      /*	EOD Fixing Flag 0: I      , 1: E */
    int eod_pay_flag,      /*	EOD Payment Flag 0: I      , 1: E */
    int eod_ex_flag,       /*	EOD Exercise Flag 0: I      , 1: E */
    double spot_fx_fixing, /*	fixing of the spot fx if today > exercise */
    int exercised,

    double *price, double *vol);

Err ConvexFwdSABR_static_replication(
    /*	Market parameters */
    long today, char *dom_yc, char *for_yc,
    double spot_fx, /* as quoted in the market */

    /*	Product parameters */
    double notional, long exercise_date, long settlmt_date, double barrier_up,
    double barrier_down,

    /*	Calib parameters */
    long *opt_matdates,  /* mkt options maturity */
    double *opt_atmvols, /* mkt options ATM vols */
    int nb_opt,

    double alpha, double beta, double rho,

    /*	Replication parameters */
    double *strikes, long nb_strike, double floorvol_strike, long is_std,
    double dig_spread, int calib_strikes,

    /* IOD / EOD flags */
    int eod_fix_flag,      /*	EOD Fixing Flag 0: I      , 1: E */
    int eod_pay_flag,      /*	EOD Payment Flag 0: I      , 1: E */
    int eod_ex_flag,       /*	EOD Exercise Flag 0: I      , 1: E */
    double spot_fx_fixing, /*	fixing of the spot fx if today > exercise */
    int exercised,

    double *price, double *coef, double *vols, double **bar_adjust);

Err FxSabrQuadSmilets(long today, double *sigma_time_fx, double *sigma_fx,
                      long sigma_n_fx, char *dom_yc, char *for_yc,
                      double spot_fx, double alpha, double gamma, double beta,
                      double rho, double lambda, double floormu,
                      /*	Product data */
                      long mat_date, double *strike, int nb_strike, int nstp,
                      int nstpfx, int nstpvol, double **res);

/*	log cvx param for levenberg optimisation */
typedef struct {
  double mat;
  double df;
  double fwd;
  double barrier_down;
  double barrier_up;
  double log_vol;
  double log_std;
  long nb_strike;
  double *strikes;
  double *coef;
  double *coefJ;
  double dig_spread;

} logconvex_param, *LOGCONVEX_PARAM;

Err fwd_fxoption(
    long today, long spot_date, double spot_fx, /*	2bd fwd */
    char *dom_yc, char *dom_vc, char *dom_ref, char *dom_swap_freq,
    char *dom_swap_basis, double dom_lam, char *for_yc, char *for_vc,
    char *for_ref, char *for_swap_freq, char *for_swap_basis, double for_lam,
    double *corr_times, double *correl_dom_for, /*	Correlations */
    double *correl_dom_fx, double *correl_for_fx, long corr_n_times,
    Err (*get_cash_vol)(/*	Function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),

    /*	Fx vol from the market */
    long *fx_mkt_vol_date, double *fx_mkt_vol, int num_fx_mkt_vol,

    /*	The structure */
    double notional, long exercise_date, long fx_fixing_date, long settlmt_date,
    double strike, int call_put,

    /* IOD / EOD flags */
    int eod_fix_flag, /*	EOD Fixing Flag 0: I      , 1: E */
    int eod_pay_flag, /*	EOD Payment Flag 0: I      , 1: E */
    int eod_ex_flag,  /*	EOD Exercise Flag 0: I      , 1: E */

    int exercised,         /*	Is exercised Flag */
    double spot_fx_fixing, /*	Fixing of the Fx */

    /* Extra params */
    char *calib_freq, long calib_dom_tau, long calib_for_tau,
    double dom_vol_shift, double for_vol_shift, double fx_vol_shift,

    /* Outputs */
    double *price, double *vol,
    int export_ts, /*	1: Export TS      , 0: don't */
    CPD_UND und_exp);

#endif
