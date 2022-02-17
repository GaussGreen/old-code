/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_closedforms_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou,ocroissant
 * Initial version
 *
 */

/*! \file ARM_xl_gp_closedforms_local.h
 *
 *  \brief file for the closed forms library
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date February 2004
 */

#ifndef ARM_XL_GP_CLOSEDFORMS_H
#define ARM_XL_GP_CLOSEDFORMS_H

#include <ARM\libarm_local\firstToBeIncluded.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"


//////////////////////////////////
/// 
/// In order to factorise code
/// many function are calling
/// common function
/// these functions are not
/// declared here 
/// 
/// only exported functions are 
/// included here to facilitate 
/// the reading of the header
/// 
//////////////////////////////////

///////////////////////////////////
/// Creates a deal description object
/// Handles the case of 
///	previous creation of object
///////////////////////////////////

////////////////////////////////////////////
/// Vanilla Option Normal
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_VanillaOption_Normal(
	LPXLOPER XL_Undelying,
	LPXLOPER XL_Volatility,
	LPXLOPER XL_Strike,
	LPXLOPER XL_Maturity,
	LPXLOPER XL_CallOrPut);

__declspec(dllexport) LPXLOPER WINAPI Local_VanillaOption_Normal_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_Undelying,
	LPXLOPER XL_Volatility,
	LPXLOPER XL_Strike,
	LPXLOPER XL_Maturity,
	LPXLOPER XL_CallOrPut);

__declspec(dllexport) LPXLOPER WINAPI Local_VanillaOption_Normal_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_Undelying,
	LPXLOPER XL_Volatility,
	LPXLOPER XL_Strike,
	LPXLOPER XL_Maturity,
	LPXLOPER XL_CallOrPut);

__declspec(dllexport) LPXLOPER WINAPI Local_DoubleDigital_Normal(
	LPXLOPER XL_Fwd1,
	LPXLOPER XL_Fwd2,
	LPXLOPER XL_Maturity,
	LPXLOPER XL_Strike1,
	LPXLOPER XL_Spread1,
	LPXLOPER XL_Strike2,
	LPXLOPER XL_Spread2,
	LPXLOPER XL_Vol1plus,
	LPXLOPER XL_Vol1minus,
	LPXLOPER XL_Vol2plus,
	LPXLOPER XL_Vol2minus,
	LPXLOPER XL_Correl,
	LPXLOPER XL_CallOrPut1,
	LPXLOPER XL_CallOrPut2);

////////////////////////////////////////////
/// Spread Option LogNormal
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n );

__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_LogNormal_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n );

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_SpreadOption_Calibrate_Correlation(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n );

__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_LogNormal_SpreadOption_Calibrate_Correlation(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n );

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n );
__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_LogNormal_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n );

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n );

__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_LogNormal_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype,
	LPXLOPER XL_n );

////////////////////////////////////////////
/// Spread Option Normal
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Normal_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype);

__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_Normal_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype);

__declspec(dllexport) LPXLOPER WINAPI Local_Normal_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype);

__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_Normal_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype);

__declspec(dllexport) LPXLOPER WINAPI Local_Normal_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype);

__declspec(dllexport) LPXLOPER WINAPI Local_Smiled_Normal_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_rho,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_slope1,
	LPXLOPER XL_slope2,
	LPXLOPER XL_callput,
	LPXLOPER XL_optiontype);

////////////////////////////////////////////
/// Spread Option SABR Gaussian
////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER flag,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec	);

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_Digital_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER flag,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec	);


__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER flag,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec	);

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_Digital_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER flag,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec	);

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER flag,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec	);

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_Digital_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER flag,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Parameters_vec	);


__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_SABR_Power_SpreadOption_Certitude(
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_n
	);

////////////////////////////////////////////////////////////////////////////////////////
///
///			Black and Sholes Formula
///
////////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholes(
	LPXLOPER XL_forward,
	LPXLOPER XL_totalvolatility,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut
	);
__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholes_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_forward,
	LPXLOPER XL_totalvolatility,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut
	);

__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholes_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_forward,
	LPXLOPER XL_totalvolatility,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut
	);

__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholes_ImplicitVol(
	LPXLOPER XL_forward,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_optprice
	);

__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholes_ImplicitVolatility(
	LPXLOPER XL_forward,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_maturity,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_optprice
	);

/////////////////////////////////////////////////////////////////////////////////////
///
/// Spread Option ShiftedLognormal Gaussian
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_ShiftedLN_Power_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sigma1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_sigma2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	);


__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_ShiftedLN_Power_SpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sigma1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_sigma2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_ShiftedLN_Power_SpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_sigma1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_sigma2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Gaussian_ShiftedLN_Power_SpreadOption_Certitude(
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_n
	);

/////////////////////////////////////////////////////////////////////////////////////
///
/// Merton Jump Diffusion Formula
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_Merton_JumpDiffusion(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_t,
	LPXLOPER XL_sigma,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_callput,
	LPXLOPER XL_n
	);


__declspec(dllexport) LPXLOPER WINAPI Local_Merton_JumpDiffusion_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_t,
	LPXLOPER XL_sigma,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_callput,
	LPXLOPER XL_n
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Merton_JumpDiffusion_Der2(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_t,
	LPXLOPER XL_sigma,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_callput,
	LPXLOPER XL_n
	);

/////////////////////////////////////////////////////////////////////////////////////
///
///		SABR Implicit Vol Formula
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_SABR_ImplicitVol(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_ImplicitVol_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_ImplicitVol_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	);
/////////////////////////////////////////////////////////////////////////////////////
///
///		SABR Vanilla Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_SABR_VanillaOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_callput,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_FromGaussianToDistribution(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_aexp,
	LPXLOPER XL_atanh,
	LPXLOPER XL_ktanh
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_FromDistributionToGaussian(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_aexp,
	LPXLOPER XL_atanh,
	LPXLOPER XL_ktanh
	);


__declspec(dllexport) LPXLOPER WINAPI Local_SABR_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_callput,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_callput,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_From_Sigma_To_Alpha(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_sigma,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps
	);

/////////////////////////////////////////////////////////////////////////////////////
///
///		Single Barriere Option in Black and Sholes Model Formula
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroBarriere(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_b,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_d,
	LPXLOPER XL_callput,
	LPXLOPER XL_inout,
	LPXLOPER XL_updown
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroBarriere_ImpliedVol(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_b,
	LPXLOPER XL_r,
	LPXLOPER XL_op,
	LPXLOPER XL_t,
	LPXLOPER XL_d,
	LPXLOPER XL_callput,
	LPXLOPER XL_inout,
	LPXLOPER XL_updown
	);


__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroBarriere_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_b,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_d,
	LPXLOPER XL_callput,
	LPXLOPER XL_inout,
	LPXLOPER XL_updown
	);
__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroBarriere_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_b,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_d,
	LPXLOPER XL_callput,
	LPXLOPER XL_inout,
	LPXLOPER XL_updown
	);
////////////////////////////////////////////////////////////////////////////////////////////
///
///				Single partial barrier finishing at a date
///
////////////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_Start_SingleBarrier(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_bendtime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	);
__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_Start_SingleBarrier_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_bendtime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	);
__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_Start_SingleBarrier_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_bendtime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	);

////////////////////////////////////////////////////////////////////////////////////////////
///
///				Single partial barrier starting at a date
///
////////////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_End_SingleBarrier(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_Starttime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_End_SingleBarrier_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_bstarttime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BS_PartialTime_End_SingleBarrier_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_barrier,
	LPXLOPER XL_rebate,
	LPXLOPER XL_v,
	LPXLOPER XL_bstarttime,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	);

////////////////////////////////////////////////////////////////////////////////////////////
///
///			Single Barrier on two assets 
///					
////////////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_BS_SingleBarrier_2Asset(
	LPXLOPER XL_f1,
	LPXLOPER XL_k1,
	LPXLOPER XL_f2,
	LPXLOPER XL_k2, 
	LPXLOPER XL_v1, 
	LPXLOPER XL_v2,
	LPXLOPER XL_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BS_SingleBarrier_2Asset_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f1,
	LPXLOPER XL_k1,
	LPXLOPER XL_f2,
	LPXLOPER XL_k2, 
	LPXLOPER XL_v1, 
	LPXLOPER XL_v2,
	LPXLOPER XL_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BS_SingleBarrier_2Asset_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f1,
	LPXLOPER XL_k1,
	LPXLOPER XL_f2,
	LPXLOPER XL_k2, 
	LPXLOPER XL_v1, 
	LPXLOPER XL_v2,
	LPXLOPER XL_corr,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_optype
	);



/////////////////////////////////////////////////////////////////////////////////////
///
///		Double Barriere Option in Black and Sholes Model Formula
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroDoubleBarriere(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_bup,
	LPXLOPER XL_bdown,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_r,
	LPXLOPER XL_b,
	LPXLOPER XL_callput
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroDoubleBarriere_ImpliedVol(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_bup,
	LPXLOPER XL_bdown,
	LPXLOPER XL_opt,
	LPXLOPER XL_t,
	LPXLOPER XL_r,
	LPXLOPER XL_b,
	LPXLOPER XL_callput
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroDoubleBarriere_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_bup,
	LPXLOPER XL_bdown,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_callput
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BS_EuroDoubleBarriere_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_bup,
	LPXLOPER XL_bdown,
	LPXLOPER XL_v,
	LPXLOPER XL_t,
	LPXLOPER XL_callput
	);
////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the bivariate function 
///					
////////////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Bivariate(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_rho,
	LPXLOPER XL_p,
	LPXLOPER XL_q
	);
////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the gamma function 
///					
////////////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Gamma(
	LPXLOPER XL_x
	);
////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the imcomplete beta function 
///					
////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ImcompleteBeta(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_z
	);

__declspec(dllexport) LPXLOPER WINAPI Local_InverseImcompleteBeta(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_z0,
	LPXLOPER XL_z
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Hypergeometric2F1(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_z0,
	LPXLOPER XL_z
	);
////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of complex error function erf(z) 
///					
////////////////////////////////////////////////////////////////////////////////////////////




__declspec(dllexport) LPXLOPER WINAPI Local_RealPart_ComplexErf(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_n
	);

__declspec(dllexport) LPXLOPER WINAPI Local_ImaginaryPart_ComplexErf(
	LPXLOPER XL_x,
	LPXLOPER XL_y,
	LPXLOPER XL_n
	);

__declspec(dllexport) LPXLOPER WINAPI Local_NormalCummulative_Inverse(
	LPXLOPER XL_x
	);

/////////////////////////////////////////////////////////////////////////////////////
///
///		Generalized Heston  Vanilla option Formula
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_VanillaOption(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_sig,
	LPXLOPER XL_t,
	LPXLOPER XL_longtermV,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi,
	LPXLOPER XL_rho,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_sig,
	LPXLOPER XL_t,
	LPXLOPER XL_longtermV,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi,
	LPXLOPER XL_rho,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_sig,
	LPXLOPER XL_t,
	LPXLOPER XL_longtermV,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi,
	LPXLOPER XL_rho,
	LPXLOPER XL_lambda,
	LPXLOPER XL_muJ,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       CEV VanillaOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CEV_VanillaOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);
	

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       CEV DoubleBarrierOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CEV_DoubleBarrierOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_Bdown,
	LPXLOPER XL_Bup,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_DoubleBarrierOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_Bdown,
	LPXLOPER XL_Bup,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);
__declspec(dllexport) LPXLOPER WINAPI Local_CEV_DoubleBarrierOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_Bdown,
	LPXLOPER XL_Bup,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       CEV BarrierOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CEV_BarrierOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_B,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,
	LPXLOPER XL_optype,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_BarrierOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_B,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,
	LPXLOPER XL_optype,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_CEV_BarrierOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_B,
	LPXLOPER XL_drift,
	LPXLOPER XL_sig,
	LPXLOPER XL_beta,	
	LPXLOPER XL_optype,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	);
		

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	);
	
__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	);
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread Digital Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);
		

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption_Der(
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);
	
__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a1,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread Digital Option2 Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption2(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_b0,	
	LPXLOPER XL_b1,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption2_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_b0,	
	LPXLOPER XL_b1,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_LogNormal_TriSpreadDigitalOption2_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_S3,
	LPXLOPER XL_sig1,
	LPXLOPER XL_sig2,
	LPXLOPER XL_sig3,	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23,
	LPXLOPER XL_mu1,
	LPXLOPER XL_mu2,
	LPXLOPER XL_mu3,	
	LPXLOPER XL_a0,
	LPXLOPER XL_a2,
	LPXLOPER XL_a3,	
	LPXLOPER XL_b0,	
	LPXLOPER XL_b1,	
	LPXLOPER XL_t,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Lognormal_Asian_VanillaOption(
	LPXLOPER XL_S,
	LPXLOPER XL_k,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_alpha,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Lognormal_Asian_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_S,
	LPXLOPER XL_k,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_alpha,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Lognormal_Asian_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_S,
	LPXLOPER XL_k,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_v,
	LPXLOPER XL_alpha,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_GaussianIntegrals_Legendre_Coeffs(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_n
	);

__declspec(dllexport) LPXLOPER WINAPI Local_GaussianIntegrals_Hermite_Coeffs(
	LPXLOPER XL_n
	);

__declspec(dllexport) LPXLOPER WINAPI Local_GaussianIntegrals_Laguerre_Coeffs(
	LPXLOPER XL_a,
	LPXLOPER XL_n
	);


__declspec(dllexport) LPXLOPER Local_ConvertFXOptionToStrike(
    LPXLOPER XL_TargetDelta,
	LPXLOPER XL_Fwd,
	LPXLOPER XL_TotalVol,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_InitValue );

/////////////////////////////////////////////////////////////////////////////////////
///
///		Generalized Heston  Vanilla option Formula With Vector of Model
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_VanillaOption_ModelVector(
	LPXLOPER XL_K,
	LPXLOPER XL_t,
	LPXLOPER XL_callput,
	LPXLOPER XL_Interpolation_Method,
	LPXLOPER XL_Maturity_Vec,
	LPXLOPER XL_F_Vec,
	LPXLOPER XL_InitialVol_Vec,
	LPXLOPER XL_longtermV_Vec,
	LPXLOPER XL_theta_Vec,
	LPXLOPER XL_ksi_Vec,
	LPXLOPER XL_rho_Vec,
	LPXLOPER XL_lambda_Vec,
	LPXLOPER XL_muJ_Vec,
	LPXLOPER XL_sigmaJ_Vec,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_Implicit_Volatility_ModelVector(
	LPXLOPER XL_K,
	LPXLOPER XL_t,
	LPXLOPER XL_Interpolation_Method,
	LPXLOPER XL_Maturity_Vec,
	LPXLOPER XL_F_Vec,
	LPXLOPER XL_InitialVol_Vec,
	LPXLOPER XL_longtermV_Vec,
	LPXLOPER XL_theta_Vec,
	LPXLOPER XL_ksi_Vec,
	LPXLOPER XL_rho_Vec,
	LPXLOPER XL_lambda_Vec,
	LPXLOPER XL_muJ_Vec,
	LPXLOPER XL_sigmaJ_Vec,
	LPXLOPER XL_nb
	);

/////////////////////////////////////////////////////////////////////////////////////
///
///		Calibration of a SABR Model to a curve of points (smile)
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Model_Calibrate(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_f,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbstep,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_alpha0,
	LPXLOPER XL_beta0,
	LPXLOPER XL_rho0,
	LPXLOPER XL_nu0,
	LPXLOPER XL_Weight_Vec
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Model_Calibrate_BetaFixed(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_f,
	LPXLOPER XL_beta,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbstep,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_alpha0,
	LPXLOPER XL_rho0,
	LPXLOPER XL_nu0,
	LPXLOPER XL_Weight_Vec,
	LPXLOPER XL_Alpha_Exp,
	LPXLOPER XL_Alpha_Tanh,
	LPXLOPER XL_Kb_Tanh
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Model_Calibrate_BetaFixed_Linked(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_f,
	LPXLOPER XL_beta,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbstep,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_alpha0,
	LPXLOPER XL_rho0,
	LPXLOPER XL_nu0,
	LPXLOPER XL_alphap,
	LPXLOPER XL_rhop,
	LPXLOPER XL_nup,
	LPXLOPER XL_rweight_alpha,
	LPXLOPER XL_rweight__rho,
	LPXLOPER XL_rweight__nu,
	LPXLOPER XL_Weight_Vec
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Model_Calibrate_WForward(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_t,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbstep,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_alpha0,
	LPXLOPER XL_beta0,
	LPXLOPER XL_rho0,
	LPXLOPER XL_nu0,
	LPXLOPER XL_f0,
	LPXLOPER XL_Weight_Vec
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Calibrate_BetaFixedToOne(
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_fwd,
	LPXLOPER XL_mat,
	LPXLOPER XL_atmvol,
	LPXLOPER XL_Weight_Vec
	);
	

/////////////////////////////////////////////////////////////////////////////////////
///
///		Calibration of GHeston
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_Model_Calibrate_Total(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_f,
	LPXLOPER XL_t,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_V0	,
	LPXLOPER XL_omega,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi	,
	LPXLOPER XL_rho	,
	LPXLOPER XL_muJ	,
	LPXLOPER XL_sigmaJ,
	LPXLOPER XL_lambda
	);

__declspec(dllexport) LPXLOPER WINAPI Local_GHeston_Model_Calibrate_NoJump(
	
	LPXLOPER XL_K_Vec,
	LPXLOPER XL_ImpVol_Vec,
	LPXLOPER XL_f,
	LPXLOPER XL_t,
	LPXLOPER XL_muJ0	,
	LPXLOPER XL_sigmaJ0,
	LPXLOPER XL_lambda0,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_algorithm,
	LPXLOPER XL_V0	,
	LPXLOPER XL_omega,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi	,
	LPXLOPER XL_rho	
	);
/////////////////////////////////////////////////////////////////////////////////////
///
///		SABR + Student  Spreadoptions
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Student_SABR_Power_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_vec,
	LPXLOPER XL_t,
	LPXLOPER flag,
	LPXLOPER XL_a10,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_Params
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Student_SABR_Power_Digital_SpreadOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_S2,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2, 
	LPXLOPER XL_rho2, 
	LPXLOPER XL_nu2,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_t,
	LPXLOPER flag,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n,
	LPXLOPER XL_alpha_exp,
	LPXLOPER XL_alpha_tanh,
	LPXLOPER XL_kb_tanh
	);
//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Mepi Vanilla Option 
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Mepi_VanillaOption_STOBS(
	LPXLOPER XL_date,
	LPXLOPER XL_P0,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_ZeroCurve,
	LPXLOPER XL_b,
	LPXLOPER XL_YearlyFees,
	LPXLOPER XL_Cashspread,
	LPXLOPER XL_SABR,
	LPXLOPER XL_minExp,
	LPXLOPER XL_maxExp,
	LPXLOPER XL_riskFac,
	LPXLOPER XL_g0,
	LPXLOPER XL_g,
	LPXLOPER XL_AveragingPeriodNb,
	LPXLOPER XL_AsianReset,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_n,
	LPXLOPER XL_nz,
	LPXLOPER XL_nh);

__declspec(dllexport) LPXLOPER WINAPI Local_Mepi_VanillaOption_STOBS_delta(
	LPXLOPER XL_date,
	LPXLOPER XL_P0,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_ZeroCurve,
	LPXLOPER XL_b,
	LPXLOPER XL_YearlyFees,
	LPXLOPER XL_Cashspread,
	LPXLOPER XL_SABR,
	LPXLOPER XL_minExp,
	LPXLOPER XL_maxExp,
	LPXLOPER XL_riskFac,
	LPXLOPER XL_g0,
	LPXLOPER XL_g,
	LPXLOPER XL_AveragingPeriodNb,
	LPXLOPER XL_AsianReset,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_n,
	LPXLOPER XL_nz,
	LPXLOPER XL_nh);


__declspec(dllexport) LPXLOPER WINAPI Local_Mepi_VanillaOption_SABR(
	LPXLOPER XL_date,
	LPXLOPER XL_P0,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_ZeroCurve,
	LPXLOPER XL_b,
	LPXLOPER XL_YearlyFees,
	LPXLOPER XL_Cashspread,
	LPXLOPER XL_SABR,
	LPXLOPER XL_SpotName,
	LPXLOPER XL_minExp,
	LPXLOPER XL_maxExp,
	LPXLOPER XL_riskFac,
	LPXLOPER XL_g0,
	LPXLOPER XL_g,
	LPXLOPER XL_AveragingPeriodNb,
	LPXLOPER XL_AsianReset,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_n,
	LPXLOPER XL_nbMC);

__declspec(dllexport) LPXLOPER WINAPI Local_Mepi_VanillaOption_SABR_delta(
	LPXLOPER XL_date,
	LPXLOPER XL_P0,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_ZeroCurve,
	LPXLOPER XL_b,
	LPXLOPER XL_YearlyFees,
	LPXLOPER XL_Cashspread,
	LPXLOPER XL_SABR,
	LPXLOPER XL_SpotName,
	LPXLOPER XL_minExp,
	LPXLOPER XL_maxExp,
	LPXLOPER XL_riskFac,
	LPXLOPER XL_g0,
	LPXLOPER XL_g,
	LPXLOPER XL_AveragingPeriodNb,
	LPXLOPER XL_AsianReset,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_n,
	LPXLOPER XL_nbMC);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Fund VanillaOption Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_VanillaOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);


__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_Ari_VanillaOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);


__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_Ari_VanillaOption_Der(
	LPXLOPER XL_i,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_LN_Ari_VanillaOption_Der2(
	LPXLOPER XL_i,
	LPXLOPER XL_j,
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Stochastic Normal volatility Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_StochasticVol_N_VanillaOption(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_r,
	LPXLOPER XL_sig,
	LPXLOPER XL_VolDrift,
	LPXLOPER XL_VolVol,	
	LPXLOPER XL_callput,
	LPXLOPER XL_nb
	);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       SKEW optimization
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER Local_OptimizeSkewVector(
		LPXLOPER XL_C_tagetVect,
		LPXLOPER XL_C_weights,
        LPXLOPER XL_C_presicions,
		LPXLOPER XL_C_InitVector,
		LPXLOPER XL_C_LBoundVector,
		LPXLOPER XL_C_UBoundVector,
		LPXLOPER XL_C_algo);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_OptimizeSkewVector(
		LPXLOPER XL_C_tagetVect,
		LPXLOPER XL_C_weights,
        LPXLOPER XL_C_presicions,
		LPXLOPER XL_C_InitVector,
		LPXLOPER XL_C_LBoundVector,
		LPXLOPER XL_C_UBoundVector,
		LPXLOPER XL_C_algo);


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       GLambda Distribution
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_GLambda_From_SABR_Calibrate(
	
	LPXLOPER XL_forward,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_T,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_Scope,
	LPXLOPER XL_IniL1,
	LPXLOPER XL_IniL2,
	LPXLOPER XL_IniL3,
	LPXLOPER XL_IniL4,
	LPXLOPER XL_IniL5,
	LPXLOPER XL_IniL6,
	LPXLOPER XL_Algo
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Student_GLambda_Power_SpreadOption(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_l1b,
	LPXLOPER XL_l2b,
	LPXLOPER XL_l3b,
	LPXLOPER XL_l4b,
	LPXLOPER XL_l5b,
	LPXLOPER XL_l6b,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_b10,
	LPXLOPER XL_k10,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	);
__declspec(dllexport) LPXLOPER WINAPI Local_Student_GLambda_Power_Digital_SpreadOption(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_l1b,
	LPXLOPER XL_l2b,
	LPXLOPER XL_l3b,
	LPXLOPER XL_l4b,
	LPXLOPER XL_l5b,
	LPXLOPER XL_l6b,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Student_GLambda_Power_Index2Digital_SpreadOption(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_l1b,
	LPXLOPER XL_l2b,
	LPXLOPER XL_l3b,
	LPXLOPER XL_l4b,
	LPXLOPER XL_l5b,
	LPXLOPER XL_l6b,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	);
__declspec(dllexport) LPXLOPER WINAPI Local_Student_GLambda_Power_Index1Digital_SpreadOption(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_l1b,
	LPXLOPER XL_l2b,
	LPXLOPER XL_l3b,
	LPXLOPER XL_l4b,
	LPXLOPER XL_l5b,
	LPXLOPER XL_l6b,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_a20,
	LPXLOPER XL_b20,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	);
__declspec(dllexport) LPXLOPER WINAPI Local_Student_GLambda_SpreadOption(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_l1b,
	LPXLOPER XL_l2b,
	LPXLOPER XL_l3b,
	LPXLOPER XL_l4b,
	LPXLOPER XL_l5b,
	LPXLOPER XL_l6b,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_k20,
	LPXLOPER XL_n
	);

////////////////////////////////////////////////////////////////////////////////////////
///
///			Lambert Function
///
////////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Lambert_Function(
	LPXLOPER XL_x
	);

////////////////////////////////////////////////////////////////////////////////////////
///
///			Black and Sholes Asymptotic Time Value Formula
///
////////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_CF_LN_VanillaOption(
	LPXLOPER XL_forward,
	LPXLOPER XL_volatility,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_maturity,
	LPXLOPER XL_CallPut
	);


__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholesTimeValue(
	LPXLOPER XL_forward,
	LPXLOPER XL_totalvolatility,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut
	);

__declspec(dllexport) LPXLOPER WINAPI Local_CF_BlackSholesTimeValue_ImplicitVol(
	LPXLOPER XL_forward,
	LPXLOPER XL_bondprice,
	LPXLOPER XL_strike,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_optprice
	);

/////////////////////////////////////////////////////////////////////////////////
/// 
///		Quantile functions
///
/////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_ShiftedLogNormal_Quantile(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_sigma,
	LPXLOPER XL_shift
	);

__declspec(dllexport) LPXLOPER WINAPI Local_ShiftedLogNormal_Distribution(
	LPXLOPER XL_f,
	LPXLOPER XL_k,
	LPXLOPER XL_t,
	LPXLOPER XL_sigma,
	LPXLOPER XL_shift
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Quantile(
	LPXLOPER XL_f,
	LPXLOPER XL_strike,
	LPXLOPER XL_T,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_SABRType,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_alpha_exp,
	LPXLOPER XL_alpha_tanh,
	LPXLOPER XL_kb_tanh


	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Distribution(
	LPXLOPER XL_f,
	LPXLOPER XL_strike,
	LPXLOPER XL_T,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_alpha_exp,
	LPXLOPER XL_alpha_tanh,
	LPXLOPER XL_kb_tanh
	);

__declspec(dllexport) LPXLOPER WINAPI Local_NormalCummulative(
	LPXLOPER XL_x
	);


__declspec(dllexport) LPXLOPER WINAPI Local_GLambda_Distribution(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_x);

__declspec(dllexport) LPXLOPER WINAPI Local_GLambda_Quantile(
	LPXLOPER XL_l1a,
	LPXLOPER XL_l2a,
	LPXLOPER XL_l3a,
	LPXLOPER XL_l4a,
	LPXLOPER XL_l5a,
	LPXLOPER XL_l6a,
	LPXLOPER XL_x
	);
__declspec(dllexport) LPXLOPER WINAPI Local_Student_Quantile(
	LPXLOPER XL_degre,
	LPXLOPER XL_x
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Student_Distribution(
	LPXLOPER XL_degre,
	LPXLOPER XL_x
	);

__declspec(dllexport) LPXLOPER WINAPI Local_ImcompleteBeta_Inverse(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Student_QIntegral(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Normal_ImpliedVol(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x,
	LPXLOPER CallPut
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Normal_Digital_ImpliedVol(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x,
	LPXLOPER CallPut
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Hypergeometric_Whittaker_W(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Hypergeometric_Whittaker_M(
	LPXLOPER XL_a,
	LPXLOPER XL_b,
	LPXLOPER XL_x
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Bessel_Y(
	LPXLOPER XL_a,
	LPXLOPER XL_x
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Bessel_I(
	LPXLOPER XL_a,
	LPXLOPER XL_x
	);


_declspec(dllexport) LPXLOPER WINAPI Local_Bessel_J(
	LPXLOPER XL_a,
	LPXLOPER XL_x
	);


__declspec(dllexport) LPXLOPER WINAPI Local_Bessel_K(
	LPXLOPER XL_a,
	LPXLOPER XL_x
	);
/////////////////////////////////////////////////////////////////////////////////////
///
///		Shifted Heston  Vanilla option Formula
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_Shifted_Heston_VanillaOption(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_sig,
	LPXLOPER XL_t,
	LPXLOPER XL_longtermV,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi,
	LPXLOPER XL_rho,
	LPXLOPER XL_Shift,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb1,
	LPXLOPER XL_nb,
	LPXLOPER XL_nbS,
	LPXLOPER XL_nbO,
	LPXLOPER XL_prec
	);
/////////////////////////////////////////////////////////////////////////////////////
///
///		SABR Heston  Vanilla option Formula
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Heston_VanillaOption(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_sig,
	LPXLOPER XL_t,
	LPXLOPER XL_longtermV,
	LPXLOPER XL_theta,
	LPXLOPER XL_ksi,
	LPXLOPER XL_rho,
	LPXLOPER XL_beta,
	LPXLOPER XL_callput,
	LPXLOPER XL_nb1,
	LPXLOPER XL_nb,
	LPXLOPER XL_nbS,
	LPXLOPER XL_nbO,
	LPXLOPER XL_prec
	);

/////////////////////////////////////////////////////////////////////////////////////
///
///		Glambda spreadoption Calibration
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_GLambda_CompleteSpreadoption(
	
	LPXLOPER XL_l1a_Vec,
	LPXLOPER XL_l2a_Vec,
	LPXLOPER XL_l3a_Vec,
	LPXLOPER XL_l4a_Vec,
	LPXLOPER XL_l5a_Vec,
	LPXLOPER XL_l6a_Vec,
	LPXLOPER XL_l1b_Vec,
	LPXLOPER XL_l2b_Vec,
	LPXLOPER XL_l3b_Vec,
	LPXLOPER XL_l4b_Vec,
	LPXLOPER XL_l5b_Vec,
	LPXLOPER XL_l6b_Vec,
	LPXLOPER XL_Discount_Vec,
	LPXLOPER XL_copula_corr,
	LPXLOPER XL_copula_degre,
	LPXLOPER XL_k,
	LPXLOPER XL_n
	);


/////////////////////////////////////////////////////////////////////////////////////
///
///		JumpDiffusion Mepi  Vanilla option Formula
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_JumpDiffusion_Mepi_Call(
	LPXLOPER XL_P0,
	LPXLOPER XL_f0,
	LPXLOPER XL_T,
	LPXLOPER XL_K,
	LPXLOPER XL_R,
	LPXLOPER XL_Emin,
	LPXLOPER XL_Lmax,
	LPXLOPER XL_gamma0,
	LPXLOPER XL_gamma1,
	LPXLOPER XL_sig,
	LPXLOPER XL_lambda,
	LPXLOPER XL_sigJ,
	LPXLOPER XL_r,
	LPXLOPER XL_s,
	LPXLOPER XL_mu,
	LPXLOPER XL_fees,
	LPXLOPER XL_volDrift,
	LPXLOPER XL_volVol,
	LPXLOPER XL_callput,
	LPXLOPER XL_params
	);


__declspec(dllexport) LPXLOPER WINAPI Local_Util_TrigonalSolve(
	
	LPXLOPER XL_A_Vec,
	LPXLOPER XL_B_Vec,
	LPXLOPER XL_C_Vec,
	LPXLOPER XL_R_Vec
	);


/////////////////////////////////////////////////////////////////////////////////////
///
///	BiSABR Vanilla Spreadoption Formula
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_SpreadOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_flag
	);




__declspec(dllexport) LPXLOPER WINAPI Local_Hypergeometric_Appell(
	LPXLOPER XL_a,
	LPXLOPER XL_b1,
	LPXLOPER XL_b2,
	LPXLOPER XL_c,
	LPXLOPER XL_x,
	LPXLOPER XL_xim,
	LPXLOPER XL_y,
	LPXLOPER XL_yim,
	LPXLOPER XL_nb);

__declspec(dllexport) LPXLOPER WINAPI Local_Util_Eigenvalues4(
	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho14,
	LPXLOPER XL_rho23,
	LPXLOPER XL_rho24,
	LPXLOPER XL_rho34
	);
__declspec(dllexport) LPXLOPER WINAPI Local_Util_Eigenvalues3(
	
	LPXLOPER XL_rho12,
	LPXLOPER XL_rho13,
	LPXLOPER XL_rho23
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Util_BiSABR_CorrelationEvolution(
	
		LPXLOPER XL_rho1,
	LPXLOPER XL_rho2,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_newrho1,
	LPXLOPER XL_newrho2
	);


__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_Calibrate(
	
	LPXLOPER	XL_F1_vec, 
	LPXLOPER	XL_alpha1_vec,
	LPXLOPER	XL_beta1_vec ,
	LPXLOPER	XL_rho1_vec,
	LPXLOPER	XL_nu1_vec,
	LPXLOPER	XL_F2_vec,
	LPXLOPER	XL_alpha2_vec,
	LPXLOPER	XL_beta2_vec ,
	LPXLOPER	XL_rho2_vec,
	LPXLOPER	XL_nu2_vec,
	LPXLOPER	XL_strike_vec,
	LPXLOPER	XL_maturity_vec,
	LPXLOPER	XL_callput_vec,
	LPXLOPER	XL_price_vec,
	LPXLOPER	XL_weight_vec,
	LPXLOPER	XL_initialparams

	);
/////////////////////////////////////////////////////////////////////////////////////
///
///	Lognormal Vanilla digital option
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_LN_DigitalOption(
	LPXLOPER XL_forward,
	LPXLOPER XL_strike,
	LPXLOPER XL_maturity,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_volatility
	);

/////////////////////////////////////////////////////////////////////////////////////
///
///	Lognormal Ratio option
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_LN_RatioOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_Mu1,
	LPXLOPER XL_Sigma1,
	LPXLOPER XL_S2,
	LPXLOPER XL_Mu2,
	LPXLOPER XL_Sigma2,
	LPXLOPER XL_Rho,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut
	);

/////////////////////////////////////////////////////////////////////////////////////
///
///	Lognormal Product option
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_LN_ProductOption(
	LPXLOPER XL_S1,
	LPXLOPER XL_Mu1,
	LPXLOPER XL_Sigma1,
	LPXLOPER XL_S2,
	LPXLOPER XL_Mu2,
	LPXLOPER XL_Sigma2,
	LPXLOPER XL_Rho,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut
	);

/////////////////////////////////////////////////////////////////////////////////////
///
///	copula base SABR calls
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_GaussianSABRDigitalCall(
	LPXLOPER XL_f1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_flag1,
	LPXLOPER XL_f2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_flag2,
	LPXLOPER XL_rho,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_LegendreNb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_GaussianSABRDigitalCallPayingS1(
	LPXLOPER XL_f1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_flag1,
	LPXLOPER XL_f2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_flag2,
	LPXLOPER XL_rho,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_LegendreNb
	);

__declspec(dllexport) LPXLOPER WINAPI Local_SABR_GaussianSABRDigitalCallPayingS2(
	LPXLOPER XL_f1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_flag1,
	LPXLOPER XL_f2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_flag2,
	LPXLOPER XL_rho,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_LegendreNb
	);


__declspec(dllexport) LPXLOPER WINAPI Local_SABR_GaussianSABRDigitalCallPayingS3(
	LPXLOPER XL_f1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_flag1,
	LPXLOPER XL_f2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_flag2,
	LPXLOPER XL_f3,
	LPXLOPER XL_sigma3,
	LPXLOPER XL_Correlations_Vec,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_Params_Vec
	);

/////////////////////////////////////////////////////////////////////////////////////
///
///	Proxy
///
/////////////////////////////////////////////////////////////////////////////////////

_declspec(dllexport) LPXLOPER WINAPI Local_TarnProxy_Create(
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_Fwds,
	LPXLOPER XL_DensityFunctors,
	LPXLOPER XL_DiscountFactors,
	LPXLOPER XL_LevPrec,	
	LPXLOPER XL_Lev,		
	LPXLOPER XL_Fix,
	LPXLOPER XL_Cap,
	LPXLOPER XL_Floor,
	LPXLOPER XL_Fees,
	LPXLOPER XL_Dcf,
	LPXLOPER XL_Target,
	LPXLOPER XL_GlobalCap,
	LPXLOPER XL_GlobalFloor,
	LPXLOPER XL_CorrelInput,
	LPXLOPER XL_NbSimul
	);

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_TarnProxy_Create(
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_Fwds,
	LPXLOPER XL_DensityFunctors,
	LPXLOPER XL_DiscountFactors,
	LPXLOPER XL_LevPrec,	
	LPXLOPER XL_Lev,		
	LPXLOPER XL_Fix,
	LPXLOPER XL_Cap,
	LPXLOPER XL_Floor,
	LPXLOPER XL_Fees,
	LPXLOPER XL_Dcf,
	LPXLOPER XL_Target,
	LPXLOPER XL_GlobalCap,
	LPXLOPER XL_GlobalFloor,
	LPXLOPER XL_CorrelInput,
	LPXLOPER XL_NbSimul
	);

__declspec(dllexport) LPXLOPER WINAPI Local_TarnProxy_GetPrice(
	LPXLOPER XL_TarnProxyId,
    LPXLOPER XL_NbPayoff);

_declspec(dllexport) LPXLOPER WINAPI Local_VBMinMaxProxy_Create(
	LPXLOPER XL_AsOf,
	LPXLOPER XL_resetDates,
	LPXLOPER XL_fwdRates,
	LPXLOPER XL_totalVol,
	LPXLOPER XL_leftVol,
	LPXLOPER XL_rightVol,
	LPXLOPER XL_nu,
	LPXLOPER XL_rho,
	LPXLOPER XL_nbSimul,
	LPXLOPER XL_sabrDiff,
	LPXLOPER XL_typeprice,
	LPXLOPER XL_type1sens,
	LPXLOPER XL_type2sens,
	LPXLOPER XL_rate1Lev,
	LPXLOPER XL_rate1Add,
	LPXLOPER XL_capRateLev,
	LPXLOPER XL_capRateAdd,
	LPXLOPER XL_vbLev,
	LPXLOPER XL_maxChoice,
	LPXLOPER XL_minChoice,
	LPXLOPER XL_minmaxFreq
	);

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_VBMinMaxProxy_Create(
	LPXLOPER XL_AsOf,
	LPXLOPER XL_resetDates,
	LPXLOPER XL_fwdRates,
	LPXLOPER XL_totalVol,
	LPXLOPER XL_leftVol,
	LPXLOPER XL_rightVol,
	LPXLOPER XL_nu,
	LPXLOPER XL_rho,
	LPXLOPER XL_nbSimul,
	LPXLOPER XL_sabrDiff,
	LPXLOPER XL_typeprice,
	LPXLOPER XL_type1sens,
	LPXLOPER XL_type2sens,
	LPXLOPER XL_rate1Lev,
	LPXLOPER XL_rate1Add,
	LPXLOPER XL_capRateLev,
	LPXLOPER XL_capRateAdd,
	LPXLOPER XL_vbLev,
	LPXLOPER XL_maxChoice,
	LPXLOPER XL_minChoice,
	LPXLOPER XL_minmaxFreq
	);

_declspec(dllexport) LPXLOPER WINAPI Local_VBMinMaxProxy_GetInfo(
	LPXLOPER XL_ProxyId,
	LPXLOPER XL_Info
	);

_declspec(dllexport) LPXLOPER WINAPI Local_SpreadVBProxy_Create(
	LPXLOPER XL_AsOf,
	LPXLOPER XL_resetDates,
	LPXLOPER XL_fwdRates1,
	LPXLOPER XL_fwdRatesVol1,
	LPXLOPER XL_fwdRatesPartVol1,
	LPXLOPER XL_vols1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_fwdRates2,
	LPXLOPER XL_fwdRatesVol2,
	LPXLOPER XL_fwdRatesPartVol2,
	LPXLOPER XL_vols2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_Rate1Rate2Correl,
	LPXLOPER XL_typeprice,
	LPXLOPER XL_nbSimul,
	LPXLOPER XL_fwdLev,
	LPXLOPER XL_fwdStrikes,
	LPXLOPER XL_type1sens,
	LPXLOPER XL_type2sens,
	LPXLOPER XL_sabrDiff,
	LPXLOPER XL_CorrMeanRev,
	LPXLOPER XL_CorrVol,
	LPXLOPER XL_FixedCpn,
	LPXLOPER XL_Levier);

_declspec(dllexport) LPXLOPER WINAPI Local_SpreadVBProxy_GetInfo(
	LPXLOPER XL_SpreadVBProxyId,
	LPXLOPER XL_info
	);

_declspec(dllexport) LPXLOPER WINAPI Local_Berm2DatesProxy_Create(
	LPXLOPER XL_AsOf,
	LPXLOPER XL_resetDates,
	LPXLOPER XL_fwdRates,
	LPXLOPER XL_strike,
	LPXLOPER XL_DFs,
	LPXLOPER XL_fwdRatesVol,
	LPXLOPER XL_fwdRatesPartVol,
	LPXLOPER XL_vols,
	LPXLOPER XL_volvols,
	LPXLOPER XL_rho,
	LPXLOPER XL_Rho12,
	LPXLOPER XL_nbSimul,
	LPXLOPER XL_TypeDiff);

_declspec(dllexport) LPXLOPER WINAPI Local_Berm2DatesProxy_GetPrice(
	LPXLOPER XL_ProxyId,
	LPXLOPER XL_Info
	);


__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_Digital_SpreadOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_flag
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_DigitalPaysS1_SpreadOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_flag
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_DigitalPaysS2_SpreadOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_flag
	);


__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_DigitalPaysS3_SpreadOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_S3params,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_flag
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_Quantile(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_flag
	);
__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_Distribution(
	LPXLOPER XL_F1,
	LPXLOPER XL_alpha1,
	LPXLOPER XL_beta1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_F2,
	LPXLOPER XL_alpha2,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_flag
	);

__declspec(dllexport) LPXLOPER WINAPI Local_BetaEqualZeroSABR(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_mu,
	LPXLOPER XL_alpha,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_CallPut
	);


__declspec(dllexport) LPXLOPER WINAPI Local_Shifted2LogNormal_Distribution(
	LPXLOPER XL_f1,
	LPXLOPER XL_sigma1,
	LPXLOPER XL_f2,
	LPXLOPER XL_sigma2,
	LPXLOPER XL_alpha,
	LPXLOPER XL_rho,
	LPXLOPER XL_T,
	LPXLOPER XL_n,
	LPXLOPER XL_K
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Shifted2LogNormal_Quantile(
	LPXLOPER XL_f1,
	LPXLOPER XL_sigma1,
	LPXLOPER XL_f2,
	LPXLOPER XL_sigma2,
	LPXLOPER XL_alpha,
	LPXLOPER XL_rho,
	LPXLOPER XL_T,
	LPXLOPER XL_n,
	LPXLOPER XL_K
	);

/////////////////////////////////////////////////////////////////////////////////////
///
///	BiSABR : Local_BiSABR_S3_SpreadOption  pays A1*(S1-S2) +B1*S3 -K1 if A2*(S1-S2) +B2*S3 -K2 >0
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABR_S3_SpreadOption(
	LPXLOPER XL_S1params,
	LPXLOPER XL_S2params,
	LPXLOPER XL_S3params,
	LPXLOPER XL_rhos,
	LPXLOPER XL_rhov,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_Correlation,
	LPXLOPER XL_T,
	LPXLOPER XL_A1,
	LPXLOPER XL_B1,
	LPXLOPER XL_K1,
	LPXLOPER XL_A2,
	LPXLOPER XL_B2,
	LPXLOPER XL_K2,
	LPXLOPER XL_flag,
	LPXLOPER XL_nbsteps
	);
__declspec(dllexport) LPXLOPER WINAPI Local_Heston_OptionPrice(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_ResetTime,
	LPXLOPER XL_Forward,
	LPXLOPER XL_Strike,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_InitialVar,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Theta,
	LPXLOPER XL_VVol,
	LPXLOPER XL_Shift,
	LPXLOPER XL_Times,
	LPXLOPER XL_Level
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Heston2B_OptionPrice(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_ResetTime,
	LPXLOPER XL_Forward,
	LPXLOPER XL_Strike,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_InitialVar1,
	LPXLOPER XL_Kappa1,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Theta1,
	LPXLOPER XL_VVol1,
	LPXLOPER XL_InitialVar2,
	LPXLOPER XL_Kappa2,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_Theta2,
	LPXLOPER XL_VVol2,
	LPXLOPER XL_Shift,
	LPXLOPER XL_Times,
	LPXLOPER XL_Level
	);

__declspec(dllexport) LPXLOPER WINAPI Local_MixteHeston_OptionPrice(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_ResetTime,
	LPXLOPER XL_Forward,
	LPXLOPER XL_Strike,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_Sigma,
	LPXLOPER XL_Weight,
	LPXLOPER XL_InitialVar,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Theta,
	LPXLOPER XL_VVol,
	LPXLOPER XL_Shift,
	LPXLOPER XL_Times,
	LPXLOPER XL_Level
	);

/////////////////////////////////////////////////////////////////////////////////////
///
///	
///			Normal Heston Vanilla option
///
///
/////////////////////////////////////////////////////////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_Normal_Heston_VanillaCall(
	LPXLOPER XL_rho,
	LPXLOPER XL_lambdaV,
	LPXLOPER XL_thetaV,
	LPXLOPER XL_kappaV,
	LPXLOPER XL_V0,
	LPXLOPER XL_S0,
	LPXLOPER XL_k,
	LPXLOPER XL_T,
	LPXLOPER XL_lambdaB,
	LPXLOPER XL_callput,
	LPXLOPER XL_nbfirst,
	LPXLOPER XL_nb,
	LPXLOPER XL_NbStage,
	LPXLOPER XL_NbOscill,
	LPXLOPER XL_prec
	);



_declspec(dllexport) LPXLOPER WINAPI Local_SABR_To_Heston_SmileCalibration_Create(
	LPXLOPER XL_resetTime,
	LPXLOPER XL_fwdRate,
	LPXLOPER XL_ATMVol,
	LPXLOPER XL_Alpha,
	LPXLOPER XL_Beta,
	LPXLOPER XL_RhoSABR,
	LPXLOPER XL_Nu,
	LPXLOPER XL_Sabr_Type,
	LPXLOPER XL_InitialVar,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Theta,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Shift);

_declspec(dllexport) LPXLOPER WINAPI Local_SABR_To_Heston_SmileCalibration_GetValue(
	LPXLOPER XL_CalibrationId,
	LPXLOPER XL_Info
	);

_declspec(dllexport) LPXLOPER WINAPI Local_SmileCalib_Spread2Heston(
	LPXLOPER XL_ResetTime,
	LPXLOPER XL_Fwd1,
	LPXLOPER XL_MktVols1,
	LPXLOPER XL_Strikes1,
	LPXLOPER XL_ConstrVol1,
	LPXLOPER XL_ConstrK1,
	LPXLOPER XL_Fwd2,
	LPXLOPER XL_MktVols2,
	LPXLOPER XL_Strikes2,
	LPXLOPER XL_ConstrVol2,
	LPXLOPER XL_ConstrK2,
	LPXLOPER XL_MktVolsSpread,
	LPXLOPER XL_StrikesSpread,
	LPXLOPER XL_ConstrVolSpread,
	LPXLOPER XL_ConstrKSpread,
	LPXLOPER XL_v0,
	LPXLOPER XL_kappa,
	LPXLOPER XL_theta
);

_declspec(dllexport) LPXLOPER WINAPI Local_Spread2HestonVanilla(
	LPXLOPER XL_ResetTime,
	LPXLOPER XL_Fwd1,
	LPXLOPER XL_Fwd2,
	LPXLOPER XL_Strike,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_V0,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Theta,
	LPXLOPER XL_Nu,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_Shift1,
	LPXLOPER XL_Shift2,
	LPXLOPER XL_Level1,
	LPXLOPER XL_Level2,
	LPXLOPER XL_Correl,
	LPXLOPER XL_Index1Lev,
	LPXLOPER XL_Index2Lev
);

_declspec(dllexport) LPXLOPER WINAPI Local_Spread2Heston_TOTEMCalibration(
	LPXLOPER XL_TOTEMMat,
	LPXLOPER XL_TOTEMStrikes,
	LPXLOPER XL_TOTEMPrices,
	LPXLOPER XL_FullScheduleReset,
	LPXLOPER XL_FullScheduleAnnuity,
	LPXLOPER XL_FullScheduleFwd1,
	LPXLOPER XL_FullScheduleFwd2,
	LPXLOPER XL_FwdCalibReset,
	LPXLOPER XL_LongFwds,
	LPXLOPER XL_LongVolsId,
	LPXLOPER XL_LongStrikesId,
	LPXLOPER XL_ShortFwds,
	LPXLOPER XL_ShortVolsId,
	LPXLOPER XL_ShortStrikesId,
	LPXLOPER XL_v0,
	LPXLOPER XL_kappa,
	LPXLOPER XL_theta,
	LPXLOPER XL_ConstrCorrel
);
///----------------------------------------------
///----------------------------------------------
///             CIR
/// Inputs :
///
///		ModelParams
///
///----------------------------------------------
///---------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI Local_CIR_ModelParamsCreate(
	LPXLOPER XL_ModelParamsId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CIR_ModelParamsCreate(
	LPXLOPER XL_ModelParamsId);

__declspec(dllexport) LPXLOPER WINAPI Local_CIR_BondPrice(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_phi);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CIR_BondPrice(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_phi);

__declspec(dllexport) LPXLOPER WINAPI Local_CIR_BondDensity(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_r,
	LPXLOPER XL_frequency);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CIR_BondDensity(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_r,
	LPXLOPER XL_frequency);

__declspec(dllexport) LPXLOPER WINAPI Local_CIR_BondDistribution(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_nbDiscr,
	LPXLOPER XL_factor);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CIR_BondDistribution(
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_time,
	LPXLOPER XL_r0,
	LPXLOPER XL_nbDiscr,
	LPXLOPER XL_factor);

__declspec(dllexport) LPXLOPER Local_CIR_Bond_Weight(
	LPXLOPER XL_ModelParamsId);

__declspec(dllexport) LPXLOPER Local_CIR_Bond_Expectation(
	LPXLOPER XL_ModelParamsId);

__declspec(dllexport) LPXLOPER Local_CIR_Bond_Variance(
	LPXLOPER XL_ModelParamsId);


///----------------------------------------------
///----------------------------------------------
///             Riccati
/// time dependant
///
///----------------------------------------------
///---------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI Local_EXP_RICCATI_Create(
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta,
	LPXLOPER XL_delta,
	LPXLOPER XL_lambda,
	LPXLOPER XL_x0,
	LPXLOPER XL_x1,
	LPXLOPER XL_x2,
	LPXLOPER XL_y0,
	LPXLOPER XL_t0);

__declspec(dllexport) LPXLOPER WINAPI Local_EXP_RICCATI_Price(
	LPXLOPER XL_RiccatiEq,
	LPXLOPER XL_t);

__declspec(dllexport) LPXLOPER WINAPI Local_EXP_RICCATI_Int(
	LPXLOPER XL_RiccatiEq,
	LPXLOPER XL_t,
	LPXLOPER XL_T);


_declspec(dllexport) LPXLOPER WINAPI Local_SmileCalibration(
	LPXLOPER	XL_AsOfDate,
	LPXLOPER	XL_CalibTimes,
	LPXLOPER	XL_Forwards,
	LPXLOPER	XL_MktVols,
	LPXLOPER	XL_Strikes,
	LPXLOPER	XL_CalibParamId,
	LPXLOPER	XL_ConstraintStrikes,
	LPXLOPER	XL_ConstraintVols,
	LPXLOPER	XL_Weights);

__declspec(dllexport) LPXLOPER WINAPI Local_TRiSABR_VanillaOption(
	LPXLOPER XL_S1params,
	LPXLOPER XL_S2params,
	LPXLOPER XL_S3params,
	LPXLOPER XL_rhos12,
	LPXLOPER XL_rhos23,
	LPXLOPER XL_rhos13,
	LPXLOPER XL_rhov12,
	LPXLOPER XL_rhov23,
	LPXLOPER XL_rhov13,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21,
	LPXLOPER XL_rhoc23,
	LPXLOPER XL_rhoc32,
	LPXLOPER XL_rhoc13,
	LPXLOPER XL_rhoc31,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallOrPut,
	LPXLOPER XL_flag
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Util_BiSABR_Eigenvalues(
	
	LPXLOPER XL_rho1,
	LPXLOPER XL_rho2,
	LPXLOPER XL_rhoS,
	LPXLOPER XL_rhoV,
	LPXLOPER XL_rhoc12,
	LPXLOPER XL_rhoc21
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Util_TriSABR_Eigenvalues(
			LPXLOPER XL_rho1,   LPXLOPER XL_rho2,   LPXLOPER XL_rho3,
			LPXLOPER XL_rhos12, LPXLOPER XL_rhos23, LPXLOPER XL_rhos13,
			LPXLOPER XL_rhov12, LPXLOPER XL_rhov23, LPXLOPER XL_rhov13,
			LPXLOPER XL_rhoc12, LPXLOPER XL_rhoc13,
			LPXLOPER XL_rhoc21, LPXLOPER XL_rhoc23,
			LPXLOPER XL_rhoc31, LPXLOPER XL_rhoc32
	);

__declspec(dllexport) LPXLOPER WINAPI Local_sabr2b_ImplicitVol(
	LPXLOPER XL_f,
	LPXLOPER XL_K,
	LPXLOPER XL_tex,
	LPXLOPER XL_alpha,
	LPXLOPER XL_beta1,
	LPXLOPER XL_beta2,
	LPXLOPER XL_rho,
	LPXLOPER XL_nu,
	LPXLOPER XL_zero,
	LPXLOPER XL_lambda
	);


/////////////////////////////////////////////////////////////////////////////////////
///
///	NonParametric Complete Option
///
/////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_CompleteOption(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_Strike2_vec,
	LPXLOPER XL_Vol2_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_S2params_vec,
	LPXLOPER XL_correlation,
	LPXLOPER XL_maturity,
	LPXLOPER XL_a1,
	LPXLOPER XL_b1,
	LPXLOPER XL_k1,
	LPXLOPER XL_a2,
	LPXLOPER XL_b2,
	LPXLOPER XL_k2,
	LPXLOPER XL_nbsteps,
	LPXLOPER XL_algo,
	LPXLOPER XL_smiletype
	);


__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_LogVolatility(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_k
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_NormalVolatility(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_k
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_NormalDistribution(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_S,
	LPXLOPER XL_T,
	LPXLOPER XL_k
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_NormalQuantile(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_S,
	LPXLOPER XL_T,
	LPXLOPER XL_k
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_LogNormalDistribution(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_S,
	LPXLOPER XL_T,
	LPXLOPER XL_k
	);

__declspec(dllexport) LPXLOPER WINAPI Local_Nonparametric_LogNormalQuantile(
	LPXLOPER XL_Strike1_vec,
	LPXLOPER XL_Vol1_vec,
	LPXLOPER XL_S1params_vec,
	LPXLOPER XL_S,
	LPXLOPER XL_T,
	LPXLOPER XL_k
	);


/////////////////////////////////////////////////////////////////////////////////////
///
///	   Bi Heston
///
/////////////////////////////////////////////////////////////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_BiShiftedHeston_VanillaOption(
	LPXLOPER XL_F1,
	LPXLOPER XL_V1,
	LPXLOPER XL_Vinfini1,
	LPXLOPER XL_lambda1,
	LPXLOPER XL_nu1,
	LPXLOPER XL_rho1,
	LPXLOPER XL_gamma1,
	LPXLOPER XL_F2,
	LPXLOPER XL_V2,
	LPXLOPER XL_Vinfini2,
	LPXLOPER XL_lambda2,
	LPXLOPER XL_nu2,
	LPXLOPER XL_rho2,
	LPXLOPER XL_gamma2,
	LPXLOPER XL_Correlations,
	LPXLOPER XL_k,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_LambdaB,
	LPXLOPER XL_Flag
	);

__declspec(dllexport)LPXLOPER WINAPI Local_Merton_VanillaOption(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_Sigma,
	LPXLOPER XL_Lambda1,
	LPXLOPER XL_U1,
	LPXLOPER XL_Lambda2,
	LPXLOPER XL_U2,
	LPXLOPER XL_N);


__declspec(dllexport)LPXLOPER WINAPI Local_BSImpliedVol(
	LPXLOPER XL_F,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_Target);

__declspec(dllexport)LPXLOPER WINAPI Local_SuperNormal_Heston_VanillaOption(
	LPXLOPER XL_S0,
	LPXLOPER XL_K,
	LPXLOPER XL_T,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Theta1,
	LPXLOPER XL_Kappa1,
	LPXLOPER XL_Nu1,
	LPXLOPER XL_V01,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_Theta2,
	LPXLOPER XL_Kappa2,
	LPXLOPER XL_Nu2,
	LPXLOPER XL_V02,
	LPXLOPER XL_Nb);

#endif
