/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file sabrimpliedvol.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_EXTENDED_SABR_H
#define _GP_CF_EXTENDED_SABR_H

#include "firsttoinc.h"
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)

double SABR_implicit_vol_direct(double f,double K,double tex, double alpha, double beta, double rho, double nu,int flag);
double SABR_Direct_FromATMsigma_to_alpha(double f,double K,double tex, double sigma, double beta, double rho, double nu,int flag);
double SABR_implicit_vol_normal(double f,double K,double tex, double alpha, double beta, double rho, double nu,int flag);
double SABR_normal_FromATMsigma_to_alpha(double f,double K,double tex, double sigma, double beta, double rho, double nu,int flag);

double SABR_ImpVol_DerNu(double f, double K, double tex, double alpha, double beta, double rho, double nu);
double SABR_ImpVol_DerRho(double f, double K, double tex, double alpha, double beta, double rho, double nu);
double SABR_ImpVol_DerBeta(double f, double K, double tex, double alpha, double beta, double rho, double nu);
double SABR_ImpVol_DerAlpha(double f, double K, double tex, double alpha, double beta, double rho, double nu);
double SABR_ImpVol_DerMaturity(double f, double K, double tex, double alpha, double beta, double rho, double nu);
double SABR_ImpVol_DerForward(double f, double K, double tex, double alpha, double beta, double rho, double nu);
double SABR_ImpVol_DerStrike(double f, double K, double tex, double alpha, double beta, double rho, double nu);

double SABR_NumericalIntegration_Call(double f, double K, double tex, double alpha, double beta, double rho, double nu, int nbsteps);
double SABR_NumericalIntegration_ImplicitVol(double f, double K, double tex, double alpha, double beta, double rho, double nu, int nbsteps);
double SABR_NumericalIntegration_Call2(double f, double K, double tex, double alpha, double beta, double rho, double nu, int nbsteps);
double SABR_NumericalIntegration_ImplicitVol2(double f, double K, double tex, double alpha, double beta, double rho, double nu, int nbsteps);

double SABR_NumericalIntegration_FromATMsigma_to_alpha(double f,double K,double tex, double sigma, double beta, double rho, double nu,int nbsteps);
double SABR_NumericalIntegration_FromATMsigma_to_alpha2(double f,double K,double tex, double sigma, double beta, double rho, double nu,int nbsteps);

double SABR_implicit_vol_direct_DerAlpha(double f, double K, double tex, double alpha, double beta, double rho, double nu,int flag);
double SABR_implicit_vol_normal_DerAlpha(double f, double K, double tex, double alpha, double beta, double rho, double nu, int flag);

double SABR_BetEqOne_CompatibleKernel_ImplVol(double f,double K,double tex, double alpha,double rho, double nu);
double SABR_BetaEqOne_CompatibleKernel_FromATMsigma_to_alpha(double f,double K,double tex, double sigma, double rho, double nu);

double SABR_Direct_FromATMsigma_to_alpha(double f,
		double K,
		double tex,
		double sigma, 
		double beta,
		double rho, 
		double nu,
		int flag);

double SABR_Normal_FromATMsigma_to_alpha(double f,
       double K,
       double tex, 
       double sigma,
       double beta,
       double rho, 
       double nu,
       int flag);

////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of sum{0,t} of exp(-a/x+b*x)/sqrt(x)
///
////////////////////////////////////////////////////////////////////////////////////////////


double SABR_expintegral(double a,double b, double t, int nbsteps);

double sabr2b_h1(double x, double beta1, double beta2, double lambda);

double sabr2b_h2(double x, double beta1, double beta2, double lambda);

double sabr2b_dh1(double x, double beta1, double beta2, double lambda);

double sabr2b_dh2(double x, double beta1, double beta2, double lambda);

double sabr2b_f(double x, double beta1, double beta2, double lambda);

double sabr2b_df(double x, double beta1, double beta2, double lambda);

double sabr2b_d2f(double x, double beta1, double beta2, double lambda);

double sabr2b_C(double x, double beta1, double beta2, double lambda);

double sabr2b_dC(double x, double beta1, double beta2, double lambda);

double sabr2b_d2C(double x, double beta1, double beta2, double lambda);

double sabr2b_I(double f,double k, double beta1, double beta2, double lambda);

double sabr2b_atmvol(	double fwd, double maturity,
						double alpha, double beta1, double beta2, double rho,
						double nu, double zero, double lambda);

double sabr2b_implicit_vol(	double fwd, double k, double maturity,
							double alpha, double beta1, double beta2, double rho,
							double nu, double zero, double lambda);

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

