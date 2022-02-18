#ifndef OPSABR_H
#define OPSABR_H

#include "utError.h"
#include "utTypes.h"

/*	New SABR pricing functions	*/
/*	Toto Dec99	*/

/*	Specific type for volatility	*/
typedef enum
{
    SABR_ATM_BETA, /*	ATM BETA VOL	*/
    SABR_ATM_LOG,  /*	ATM LOGNORMAL VOL	*/
    SABR_ATM_NORM, /*	ATM NORMAL VOL	*/
    SABR_STR_LOG,  /*	LOGNORMAL VOL AT THE STRIKE	*/
    SABR_STR_NORM  /*	NORMAL VOL AT THE STRIKE	*/
} SABR_VOL_TYPE;

/*	Risk parameters	*/
typedef struct SABR_RISK_PARAM
{
    /*	Shift sizes	*/
    double delta_shift;
    double gamma_shift;
    double vega_shift;
    double vanna_shift;
    double volga_shift;
    double theta_shift;
    double alpha_shift;
    double beta_shift;
    double rho_shift;

    /*	Multiplicative shifts on vega?	*/
    int vega_mult;
    int volga_mult;

    /*	Factors that multiply greeks	*/
    double delta_fac;
    double gamma_fac;
    double vega_fac;
    double volga_fac;
    double vanna_fac;
    double theta_fac;
    double alpha_fac;
    double beta_fac;
    double rho_fac;

    /*	Which vol to freeze/bump	*/
    SABR_VOL_TYPE delta_freeze_vol;
    SABR_VOL_TYPE theta_freeze_vol;
    SABR_VOL_TYPE alpha_beta_rho_freeze_vol;
    SABR_VOL_TYPE vega_bump_vol;
} SABR_RISK_PARAM;

/*	Interpret volatility type	*/
Err interp_sabr_vol_type(char* str, SABR_VOL_TYPE* val);

/*	Setup parameters	*/
Err op_sabr_set_param(int num_param, char** param_str, char** value_str, SABR_RISK_PARAM* param);

/*	Conversion from one vol type to another	*/
Err vol_conv(
    double        input_vol,
    SABR_VOL_TYPE input_vol_type,
    double*       output_vol,
    SABR_VOL_TYPE output_vol_type,
    /*	Forward	*/
    double forward,
    /*	Option specs	*/
    double strike,
    double maturity,
    /*	Model specs	*/
    double alpha,
    double beta,
    double rho);

/*	Pricing	*/
/*	All parameters checks are supposed to have been done before call	*/
Err op_sabr_pricing(
    /*	Black-Scholes parameters	*/
    double         forward,
    double         strike,
    double         maturity,
    double         disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha,
    double beta,
    double rho,
    /*	Vol input	*/
    double        input_vol,
    SABR_VOL_TYPE input_vol_type,
    /*	Premium	*/
    double* prem);

/*	Delta/Gamma	*/
/*	All parameters checks are supposed to have been done before call	*/
Err op_sabr_delta_gamma(
    /*	Black-Scholes parameters	*/
    double         forward,
    double         strike,
    double         maturity,
    double         disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha,
    double beta,
    double rho,
    /*	Vol input	*/
    double        input_vol,
    SABR_VOL_TYPE input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double delta_shift,
    double gamma_shift,
    /*	Which vol to freeze when bumping	*/
    SABR_VOL_TYPE freeze_vol_type,
    /*	Calculate delta or gamma	*/
    /*	0: Delta, 1: Gamma	*/
    int delta_gamma,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answer	*/
    double* sens);

/*	Vega/Volga	*/
/*	All parameters checks are supposed to have been done before call	*/
Err op_sabr_vega_volga(
    /*	Black-Scholes parameters	*/
    double         forward,
    double         strike,
    double         maturity,
    double         disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha,
    double beta,
    double rho,
    /*	Vol input	*/
    double        input_vol,
    SABR_VOL_TYPE input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double vega_shift,
    double volga_shift,
    /*	Multiplicative or additive	*/
    /*	0: Additive, 1: Multiplicative	*/
    int vega_mult,
    int volga_mult,
    /*	Which vol to bump	*/
    SABR_VOL_TYPE bump_vol_type,
    /*	Calculate vega or volga	*/
    /*	0: Vega, 1: Volga	*/
    int vega_volga,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answer	*/
    double* sens);

/*	Vanna	*/
/*	All parameters checks are supposed to have been done before call	*/
Err op_sabr_vanna(
    /*	Black-Scholes parameters	*/
    double         forward,
    double         strike,
    double         maturity,
    double         disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha,
    double beta,
    double rho,
    /*	Vol input	*/
    double        input_vol,
    SABR_VOL_TYPE input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double delta_shift,
    double vega_shift,
    /*	Multiplicative or additive vega shift */
    /*	0: Additive, 1: Multiplicative	*/
    int vega_mult,
    /*	Which vol to bump	*/
    SABR_VOL_TYPE bump_vol_type,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answer	*/
    double* vanna);

/*	Theta	*/
/*	All parameters checks are supposed to have been done before call	*/
Err op_sabr_theta(
    /*	Black-Scholes parameters	*/
    double         forward,
    double         strike,
    double         maturity,
    double         disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha,
    double beta,
    double rho,
    /*	Vol input	*/
    double        input_vol,
    SABR_VOL_TYPE input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double time_shift,
    /*	Which vol to freeze when bumping	*/
    SABR_VOL_TYPE freeze_vol_type,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answers	*/
    double* theta,
    double* thetagamma,
    double* thetavolga,
    double* thetavanna,
    double* thetacarry);

/*	Sensitivity to model parameters	*/
/*	All parameters checks are supposed to have been done before call	*/
Err op_sabr_model_sens(
    /*	Black-Scholes parameters	*/
    double         forward,
    double         strike,
    double         maturity,
    double         disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha,
    double beta,
    double rho,
    /*	Vol input	*/
    double        input_vol,
    SABR_VOL_TYPE input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double alpha_shift,
    double beta_shift,
    double rho_shift,
    /*	Which vol to freeze when bumping	*/
    SABR_VOL_TYPE freeze_vol_type,
    /*	Factor that multiplies results	*/
    double alpha_fac,
    double beta_fac,
    double rho_fac,
    /*	Answers	*/
    double* alpha_sens,
    double* beta_sens,
    double* rho_sens);

Err op_sabr_parabola(
    double forward,
    double maturity,

    /*	Model parameters	*/
    double alpha,
    double beta,
    double rho,

    /*	Vol input	*/
    double        input_vol,
    SABR_VOL_TYPE input_vol_type,

    /*	Options	*/
    double alpha_shift,
    double beta_shift,
    double rho_shift,
    /*	Which vol to freeze when bumping	*/
    SABR_VOL_TYPE freeze_vol_type,
    /*	Factor that multiplies results	*/
    double alpha_fac,
    double beta_fac,
    double rho_fac,
    /*	Answers	*/
    double* alpha_sens,
    double* beta_sens,
    double* rho_sens);

#endif