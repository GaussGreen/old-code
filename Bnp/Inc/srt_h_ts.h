/* SRT_H_TS.H */
#ifndef SRT_H_TS_H
#define SRT_H_TS_H

#include "utDates.h"
#include "utListStruct.h"

/*
#include "srt_h_all.h"
*/
/* -------------------------------------------------------------------------------
 */

/* ----------------------- The TERM STRUCT is ANOTHER SORT LIST
 * ------------------ */

typedef SrtList TermStruct;

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

enum VAL_ORIGIN {
  SIGMA_DATE,   /* value comes from a sigma date */
  TAU_DATE,     /* value comes from a tau date */
  BOTH_DATE,    /* value comes from both sigma and tau date */
  INTERPOLATED, /* value comes from interpolation */
  VASICEK_MEAN_REV_DATE,
  OTHER /* other */

};

/* -------------------------------------------------------------------------- */

typedef enum TermStruct_type {
  PIECEWISE_CONST_VOL_LAMBDA_TS,
  CONST_VOL_LAMBDA_TS,
  PIECEWISE_CONST_VOL_TAU_TS,
  TWOFAC_PIECEWISE_CONST_VOL_TAU_TS,
  BS_TS
} TermStructType;

/* -------------------------------------------------------------------------- */

/* -------------- INTEREST RATE MODELS TERM STRUCTURES CONTENTS ---------------
 */

typedef struct irmtermstructval {
  Ddate date;
  int val_origin; /* wether value comes from sigma or tau date */
  double time;
  double sig;
  double tau;

  double mean_rev_level;
  double vasicek_init_cond;
  double vasicek_mean_sr;
  double vasicek_mean_int_sr;

  double F;
  double Psi;
  double G;
  double H;

  /*For Jumping Numeraire LGM*/
  double I;
  double J;
  double K;
  double L;
  double N;
  double O;
  double Q;
  double R;
  double S;
  double T;
  double Phi;

  double int_phi;

  /* For LGM Jumping Quanto adjustment */
  double M;
  double dFxVolTimesCor;

  /* Stoch Vol  model */
  double vovol;
  double rho;
  double meanvol;

  double beta;
  double eta;
  double omega; /* == Beta 2 - BEta 1 */
  double Zeta; /* aka tau in my notes aka int^t_0{.5*alpha(t')*alpha(t')}dt' psh
                  paper */
  double **M_beta_eta;
  SRT_Boolean is_M_allocated_here;
  double Lambda;

} IrmTermStructVal;

/* -------------------------------------------------------------------------- */

typedef struct ExpoVolStr_ {
  double sig;
  double sig2; /* sigma square*/

  double tau; /* 1/lambda*/
  double lambda;

  double beta; /* The Power of r */

  double F; /* exp -[ int(u=0      ,u=s) lam(u) du  */
  double
      Psi; /* int(s=0      ,s=t) [ exp -[ int(u=0      ,u=s) lam(u) du ] ds ]*/
  double Kappa; /* int(s=0      ,s=t) [ exp [ int(u=0      ,u=s) lam(u) du ] ds
                   ] */

  double J; /* for quanto only*/
  double M; /* for quanto only*/
} ExpoVolStr;

/* -------------------------------------------------------------------------- */

/* These functions have been removed from the TermStruct because their number
   depends on the number of factor used...*/
typedef struct lgmrebldfuncts {
  double G;
  double H;
  double stdev_x; /* In a jumping numeraire from Ti to Ti+1:
                     this is the sqrt of a variance matrix...
                     The off-diagonal elements correspond to
                     the correlation value
                  */
} LGMRebFnc;

/* -------------------------------------------------------------------------- */

typedef struct TF_TSVal_ {
  Ddate date;
  double time;
  int val_origin; /* wether value comes from sigma or tau date */

  /* Correlation linked parameters */
  double alpha;
  double gamma;
  double rho;

  /* Smile parameter */
  double omega;

  /* Volatility strucutre for each factor */
  ExpoVolStr exp_ts[2];

  /* Reconstruction formula linked functions */
  LGMRebFnc reb_ts[2][2];

} TwoFacIrmTermStructVal;

/* ----------------------------------------------------------------------------------
 */
typedef struct StochRatesTsVal {
  double sig;
  double df;
  double fwd;
  double tau;
  double lambda;
  double F;
  double G;
  double H;
  double Psi;
  double I;
  double J;
  double K;
  double L;
  double N;
  double O;
  double Q;
  double Phi;

} StochRatesTsVal;

/* ---------------- EQUITY MODELS TERM STRUCTURES CONTENTS -------------------
 */

typedef struct {
  Ddate date;
  double time;
  double sig;

  /* SRVGS Model Parameters */
  double omega;
  double beta;
  double gamma;

  double basevol;
  double voldrift;
  double vovol;
  double rho_spot_vol;

  double int_sig_dt;
  double int_sig2_dt;

  StochRatesTsVal StochRate_Ts;
  double rho_e_ir;
  double V_ir_e;
  double W_ir_e;

} EquityTermStructVal;

/* -------------------------------------------------------------------------- */

typedef struct FxTermStructVal {

  Ddate date;
  double time;
  double sigx;
  double int_sig_dt;
  double int_sig2_dt;

  /* For FX Stochastic Rates Only with Jumping  */
  double rhofd;
  double rhofx;
  double rhodx;
  StochRatesTsVal StochRates_Ts[2]; /* One for each IR und ; 1 = domestic und */
  double M_fx;
  double N_fx;
  double V_dx;
  double W_dx;
  double V_fx;
  double W_fx;
  double O_fd;
  double P_fd;
  double R_fd;
  double M_fd;
  double Phi_fd;
  double H_fd;
  double Q_fd;
  double T_fd;
  double U_fd;
  double V_fd;
  double S_fd;
  double W_fd;
  double X_dx;
  double Y_dx;
  double X_fx;
  double Y_fx;

} FxTermStructVal;

/* -------------------------------------------------------------------------- */

#endif
