#ifndef OPSABRAF_H
#define	OPSABRAF_H

#include "utError.h"
#include "utTypes.h"

/*	New SABRAF pricing functions	*/
/*	Aug04	*/
/*	J.B / P.T */

/*	Risk parameters	*/
typedef struct SABRAF_RISK_PARAM
{
	/*	Shift sizes	*/
	double				delta_shift;
	double				gamma_shift;
	double				vega_shift;
	double				volga_shift;
	double				theta_shift;
	double				alpha_shift;
	double				beta_shift;
	double				rho_shift;
	double				zeta_shift;

	/*	Multiplicative shifts on vega?	*/
	int					vega_mult;
	int					volga_mult;

	int					alpha_mult; 
	double				beta_mult;
	double				rho_mult;	
	double				zeta_mult;	

	/*	Factors that multiply greeks	*/
	double				delta_fac;
	double				gamma_fac;
	double				vega_fac;
	double				volga_fac;
	double				vanna_fac;
	double				theta_fac;
	double				alpha_fac;
	double				beta_fac;
	double				rho_fac;
	double				zeta_fac;

	/*	Which vol to freeze/bump	*/
	SrtDiffusionType		delta_freeze_vol;
	SrtDiffusionType		theta_freeze_vol;
	SrtDiffusionType		alpha_beta_rho_zeta_freeze_vol;
	SrtDiffusionType		vega_bump_vol;
} SABRAF_RISK_PARAM;

/*	Setup parameters	*/
Err op_sabraf_set_param(
							int						num_param,
							char					**param_str,
							char					**value_str,
							SABRAF_RISK_PARAM			*param);

/*	Pricing	*/
/*	All parameters checks are supposed to have been done before call	*/
Err srt_f_op_sabraf_price(	double				Forward,
							double				Strike,
							double				Maturity,
							double				Disc,
							SrtCallPutType		call_put, 
							double				Sigma,
							double				Alpha,
							double				Beta,
							double				Rho,
							double				Zeta,
							SrtDiffusionType	input,
							double				*price);

/*	Delta/Gamma	*/
/*	All parameters checks are supposed to have been done before call	*/
Err op_sabraf_delta_gamma(
							/*	Black-Scholes parameters	*/
							double					forward,							
							double					strike,
							double					maturity,
							double					disc,
							SrtCallPutType			call_put, 
							/*	Model parameters	*/
							double					alpha,
							double					beta,
							double					rho,
							double					zeta,
							/*	Vol input	*/
							double					input_vol,
							SrtDiffusionType			input_vol_type,
							/*	Options	*/
							/*	Size of bumps	*/
							double					delta_shift,	
							double					gamma_shift,
							/*	Which vol to freeze when bumping	*/
							SrtDiffusionType		freeze_vol_type,
							/*	Calculate delta or gamma	*/
							/*	0: Delta, 1: Gamma	*/
							int						delta_gamma,		
							/*	Factor that multiplies result	*/
							double					fac,			
							/*	Answer	*/
							double					*sens);

/*	Vega/Volga	*/
/*	All parameters checks are supposed to have been done before call	*/
Err op_sabraf_vega_volga(
							/*	Black-Scholes parameters	*/
							double					forward,							
							double					strike,
							double					maturity,
							double					disc,
							SrtCallPutType			call_put, 
							/*	Model parameters	*/
							double					alpha,
							double					beta,
							double					rho,
							double					zeta,
							/*	Vol input	*/
							double					input_vol,
							SrtDiffusionType		input_vol_type,
							/*	Options	*/
							/*	Size of bumps	*/
							double					vega_shift,	
							double					volga_shift,
							/*	Multiplicative or additive	*/
							/*	0: Additive, 1: Multiplicative	*/
							int						vega_mult,
							int						volga_mult,
							/*	Which vol to bump	*/
							SrtDiffusionType			bump_vol_type,
							/*	Calculate vega or volga	*/
							/*	0: Vega, 1: Volga	*/
							int						vega_volga,		
							/*	Factor that multiplies result	*/
							double					fac,			
							/*	Answer	*/
							double					*sens);

/*	Vanna	*/
/*	All parameters checks are supposed to have been done before call	*/
Err op_SABRAF_vanna(
						/*	Black-Scholes parameters	*/
						double					forward,							
						double					strike,
						double					maturity,
						double					disc,
						SrtCallPutType			call_put, 
						/*	Model parameters	*/
						double					alpha,
						double					beta,
						double					rho,
						double					zeta,
						/*	Vol input	*/
						double					input_vol,
						SrtDiffusionType		input_vol_type,
						/*	Options	*/
						/*	Size of bumps	*/
						double					delta_shift,	
						double					vega_shift,
						/*	Multiplicative or additive vega shift */
						/*	0: Additive, 1: Multiplicative	*/
						int						vega_mult,
						/*	Which vol to bump	*/
						SrtDiffusionType		bump_vol_type,
						/*	Factor that multiplies result	*/
						double					fac,			
						/*	Answer	*/
						double					*vanna);

/*	Theta	*/
/*	All parameters checks are supposed to have been done before call	*/
Err op_SABRAF_theta(
						/*	Black-Scholes parameters	*/
						double					forward,							
						double					strike,
						double					maturity,
						double					disc,
						SrtCallPutType			call_put, 
						/*	Model parameters	*/
						double					alpha,
						double					beta,
						double					rho,
						double					zeta,
						/*	Vol input	*/
						double					input_vol,
						SrtDiffusionType		input_vol_type,
						/*	Options	*/
						/*	Size of bumps	*/
						double					time_shift,	
						/*	Which vol to freeze when bumping	*/
						SrtDiffusionType		freeze_vol_type,
						/*	Factor that multiplies result	*/
						double					fac,			
						/*	Answers	*/
						double					*theta,
						double					*thetacarry);

/*	Sensitivity to model parameters	*/
/*	All parameters checks are supposed to have been done before call	*/
Err op_SABRAF_model_sens(
						/*	Black-Scholes parameters	*/
						double					forward,							
						double					strike,
						double					maturity,
						double					disc,
						SrtCallPutType			call_put, 
						/*	Model parameters	*/
						double					alpha,
						double					beta,
						double					rho,
						double					zeta,
						/*	Vol input	*/
						double					input_vol,
						SrtDiffusionType		input_vol_type,
						/*	Options	*/
						/*	Size of bumps	*/
						double					alpha_shift,	
						double					beta_shift,	
						double					rho_shift,	
						double					zeta_shift,	
						/* mult or add */
						double					alpha_mult,	
						double					beta_mult,	
						double					rho_mult,	
						double					zeta_mult,	
						/*	Which vol to freeze when bumping	*/
						SrtDiffusionType		freeze_vol_type,
						/*	Factor that multiplies results	*/
						double					alpha_fac,
						double					beta_fac,
						double					rho_fac,
						double					zeta_fac,
						/*	Answers	*/
						double					*alpha_sens,
						double					*beta_sens,
						double					*rho_sens,
						double					*zeta_sens);
#endif