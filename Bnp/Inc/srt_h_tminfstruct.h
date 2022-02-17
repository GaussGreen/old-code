#ifndef SRT_H_TMINFSTRUCT_H
#define SRT_H_TMINFSTRUCT_H

/* -----------------------------------------------------------------
  TYPE            :SrtBasicTmInf
  AUTHOR          :Ka Lok Chau
  DESCRIPTION     :information needed for deter mdl at a particular time
  DEFINITION      :

  ---------------------------------------------------------------------*/
typedef struct {
  double df; /* B(0      ,t) */
} SrtBasicTmInf;

/* -----------------------------------------------------------------
  TYPE            :SrtIRTmInf
  DESCRIPTION     :information needed for the discretisation of an
                        INTEREST RATE MODEL at a particular time
  DEFINITION      :

  ---------------------------------------------------------------------*/
typedef struct {
  double df;         /* B(0      ,t) */
  SrtMdlDim mdl_dim; /* Number of factors		*/
  union {
    ExpoVolStr onef;
    ExpoVolStr twof[2];
  } ev;
  union {
    LGMRebFnc onef;
    LGMRebFnc twof[2][2];
  } rf;
  double correl_x;   /* For two factor only(for the moment) can be extended to
  real matrix when   many factors */
  SrtSample fwd_sam; /* forward sample */
  YTatt_param yp;    /* FOR DF FROM NODE TO NEXT NODE */

  double quanto_adjustment;
  double fxquanto_adjustment; /*if it is a foreign underlying*/
  double ffdquanto_adjustment;
  double sfdquanto_adjustment;

  /* For stochastic vol only */
  double vovol;
  double rho;
  double vovol_sqr;

  /* For the Eta-Beta model */
  double the_beta;
  double eta;
  double lambda_t;
  double zeta_t;
  double A_t;

  double mean_rev_level;
  double vasicek_init_cond;

} SrtIRMTmInf;

/* -----------------------------------------------------------------
  TYPE            :SrtLogTmInf
  AUTHOR          :Ka Lok Chau
  DESCRIPTION     :information needed for BS mdl at a particular time
  DEFINITION      :

  ---------------------------------------------------------------------*/
typedef struct {
  double df;

  double int_sig2_dt;        /* cum vol from node to next node */
  double sqrt_int_sig2_dt;   /* sqrt of cum vol from node to next node */
  double int_sig_dt;         /* integral of sig*dt from node to next node */
  double quanto_adjustment;  /*quanto_adjustment if it is a foreign under*/
  double init_fwd_val;       /* initial (expected) fwd */
  double vol2_cum;           /* Integral of sigma^2 * dt from 0 till 't */
  double ajust_init_fwd_val; /* Forward adjusted by exp(-.5*sig^2*T) */
  double inv_exp_int_spread_dt;

  double drift; /* The r for equities      , or rd-rf for FX */

  /* stochastic rate & volatility gamma smile model parameters */
  double sig;
  double omega;
  double beta;
  double gamma;
  double basevol;
  double voldrift;
  double vovol;
  double rho;

  double time;

} SrtLogTmInf;

/* -----------------------------------------------------------------
  TYPE            :SrtFXTmInf
  AUTHOR          :R. Benichou
  DESCRIPTION     :information needed for the discretisation of an
                        foreign underlying with stochastic rates
  DEFINITION      :

  ---------------------------------------------------------------------*/
typedef struct {

  double df;
  double int_sigx_dt;
  double int_sigx2_dt;
  double sqrt_int_sigx2_dt;
  double init_spot;

  double mean_ln_fx;
  double var_ln_fx;
  double l[3][3];
  StochRatesTsVal StochRatesVal[2];

  double dom_lambda_dt;
  double for_lambda_dt;

  double mean_int_dom_sr_dt;
  double mean_int_for_sr_dt;
  double var_int_dom_sr;
  double var_int_for_sr;

  double cov_int_dom_for;
  double cov_int_dom_fx;
  double cov_int_for_fx;

  double M_fx;
  double N_fx;
  double O_fx;
  double O_fd;
  double P_fd;
  double Phi_fd;
  double Q_fd;
  double S_fd;
  double T_fd;
  double U_fd;
  double V_fd;
  double W_fd;
  double V_dx;

  double quanto_adjustment; /*quanto_adjustment if it is a foreign under*/

} SrtFXTmInf;

/*<%%STA-----------------------------------------------------------------
  TYPE            :SrtOneFacTmInf
  AUTHOR          :E.Auld
  DESCRIPTION     :information needed for an interest rate model
                                at a particular time
  DEFINITION      :

typedef struct{
  double df;
  double sig;
  double sig2;
  double lambda;
  double tau;
  double F      ,G      ,H      ,J      ,Psi      ,stdev_r;
  SrtSample fwd_sam;
  YTatt_param  yp;
}SrtOneFacTmInf;
<%%END---------------------------------------------------------------------*/

/*<%%STA-----------------------------------------------------------------
  TYPE            :SrtTwoFacTmInf
  AUTHOR          :O. Van Eyseren
  DESCRIPTION     :information needed for the discretisation of a
                        two factor interest rate model at a particular time
  DEFINITION      :

typedef struct{
  double df;
  double sig;
  double int_sig2_dt;
  double lambda;
  double tau;
  double F[2]      ,Psi[2]      ,J[2];
  double stdev_x[2]      ,correl_x;
  double G[2][2]      ,H[2][2];
  SrtSample fwd_sam;
  YTatt_param  yp;
}SrtTwoFacTmInf;

<%%END---------------------------------------------------------------------*/

#endif
