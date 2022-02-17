/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_closedforms.h,v $
 * Revision 1.1  2004/05/02 15:08:43  ocroissant
 * Initial version
 *
 */

 /*! \file ARM_local_gp_closedforms.h
 *
 *  \brief file for the various closed forms in the generic pricer
 *
 *	\author  O Croissant
 *	\version 1.0
 *	\date February 2004
 */

#ifndef ARMLOCAL_GP_CLOSEDFORMS_H
#define ARMLOCAL_GP_CLOSEDFORMS_H

#include "firstToBeIncluded.h"
#include "ARM_local_gp_genericaddin.h"
#include <vector>
using std::vector;
#include <string>
#include <deque>

/// forward declaration
class ARM_result;
class ARM_ZeroCurve;
//class ARM_SABR_Eq;
////////////////////////////////////////////
/// Vanilla Option Normal
////////////////////////////////////////////
extern long ARMLOCAL_VanillaOption_Normal(
	double underling,
	double volatility,
	double maturity,
	double strike,
	int callput,
	ARM_result& result
	);

extern long ARMLOCAL_VanillaOption_Normal_Der(
	double i,
	double underling,
	double volatility,
	double maturity,
	double strike,
	int callput,
	ARM_result& result
	);

extern long ARMLOCAL_VanillaOption_Normal_Der2(
	double i,
	double j,
	double underling,
	double volatility,
	double maturity,
	double strike,
	int callput,
	ARM_result& result
	);

extern long ARMLOCAL_DoubleDigital_Normal(
	   double fwd1, 
	   double fwd2,
	   double maturity,
	   double K1, double spread1,
	   double K2, double spread2,
	   double vol1plus, double vol1minus,
	   double vol2plus, double vol2minus,
	   double correl,
	   int callorput1,
	   int callorput2,
	   ARM_result& result);

////////////////////////////////////////////
/// Spread Option Lognormal
////////////////////////////////////////////
extern long ARMLOCAL_LogNormal_SpreadOption(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	);

extern long ARMLOCAL_LogNormal_SpreadOption_Calibrate_Correlation(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	);

extern long ARMLOCAL_Smiled_LogNormal_SpreadOption_Calibrate_Correlation(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	);

extern long ARMLOCAL_Smiled_LogNormal_SpreadOption(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	);


extern long ARMLOCAL_LogNormal_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	);

extern long ARMLOCAL_Smiled_LogNormal_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	);

extern long ARMLOCAL_LogNormal_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	);

extern long ARMLOCAL_Smiled_LogNormal_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	);

////////////////////////////////////////////
/// Spread Option Normal
////////////////////////////////////////////
extern long ARMLOCAL_Normal_SpreadOption(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	ARM_result& result
	);

extern long ARMLOCAL_Smiled_Normal_SpreadOption(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	ARM_result& result
	);


extern long ARMLOCAL_Normal_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	ARM_result& result
	);

extern long ARMLOCAL_Smiled_Normal_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	ARM_result& result
	);

extern long ARMLOCAL_Normal_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	ARM_result& result
	);

extern long ARMLOCAL_Smiled_Normal_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	ARM_result& result
	);

////////////////////////////////////////////
/// Spread Option SABR Gaussian
////////////////////////////////////////////

extern long ARMLOCAL_Gaussian_SABR_Power_SpreadOption(
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	);

extern long ARMLOCAL_Gaussian_SABR_Power_Digital_SpreadOption(
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	);

extern long ARMLOCAL_Gaussian_SABR_Power_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	);

extern long ARMLOCAL_Gaussian_SABR_Power_Digital_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	);

extern long ARMLOCAL_Gaussian_SABR_Power_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	);

extern long ARMLOCAL_Gaussian_SABR_Power_Digital_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	);

extern long ARMLOCAL_Gaussian_SABR_Power_SpreadOption_Certitude(
	double copula_corr,
	double t,
	int n,
	ARM_result& result
	);

////////////////////////////////////////////////////////////////////////////////////////
///
///			Black and Sholes Formula
///
////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_CF_BlackSholes(
								 double forward,
								 double totalvolatility,
								 double bondprice,
								 double strike,
								 double CallPut,
								 ARM_result& result
	);


extern long ARMLOCAL_CF_BlackSholes_Der(
									 int i,
									 double forward,
									 double totalvolatility,
									 double bondprice,
									 double strike,
									 double CallPut,
									 ARM_result& result
									 );
extern long ARMLOCAL_CF_BlackSholes_Der2(
									  int i,
									  int j,
									  double forward,
									  double totalvolatility,
									  double bondprice,
									  double strike,
									  double CallPut,
									  ARM_result& result
									  );

extern long ARMLOCAL_CF_BlackSholes_ImplicitVol(
								 double forward,
								 double bondprice,
								 double strike,
								 double CallPut,
								 double optprice,
								 int algo,
								 ARM_result& result
	);

extern long ARMLOCAL_CF_BlackSholes_ImplicitVolatility(
												double forward,
												double bondprice,
												double strike,
												double maturity,
												double CallPut,
												double optprice,
												int algo,
												ARM_result& result
												);
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        ShiftedLognormal Power Spreadoption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption(
	double S1,
	double S2,
	double sigma1,
	double alpha1,
	double sigma2,
	double alpha2,
	double copula_corr,
	double t,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	);




extern long ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double sigma1,
	double alpha1,
	double sigma2,
	double alpha2,
	double copula_corr,
	double t,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	);


extern long ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double sigma1,
	double alpha1,
	double sigma2,
	double alpha2,
	double copula_corr,
	double t,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	);

extern long ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption_Certitude(
	double copula_corr,
	double t,
	int n,
	ARM_result& result
	);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Merton Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_Merton_JumpDiffusion(
		double F,
		double K,
		double t,
		double sigma,
		double lambda,
		double muJ, 
		double sigmaJ,
		int callorput,int nb,
		ARM_result& result
		);




extern long ARMLOCAL_Merton_JumpDiffusion_Der(
		int i,
		double F,
		double K,
		double t,
		double sigma,
		double lambda,
		double muJ, 
		double sigmaJ,
		int callorput,
		int nb,
		ARM_result& result
	);


extern long ARMLOCAL_Merton_JumpDiffusion_Der2(
		int i,
		int j,
		double F,
		double K,
		double t,
		double sigma,
		double lambda,
		double muJ, 
		double sigmaJ,
		int callorput,
		int nb,
		ARM_result& result
	);
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        SABR Implicit Vol Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_SABR_ImplicitVol(
									  double f, 
									  double K, 
									  double tex,
									  double alpha, 
									  double beta, 
									  double rho, 
									  double nu,
									  int flag,
									  int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
									  ARM_result& result
									  );




extern long ARMLOCAL_SABR_ImplicitVol_Der(
										  int i,
										  double f, 
										  double K, 
										  double tex,
										  double alpha, 
										  double beta, 
										  double rho, 
										  double nu,
										  int flag,
										  int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
										  ARM_result& result
										  );


extern long ARMLOCAL_SABR_ImplicitVol_Der2(
										   int i,
										   int j,
										   double f, 
										   double K, 
										   double tex,
										   double alpha, 
										   double beta, 
										   double rho, 
										   double nu,
										   int flag,
										   int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
										   ARM_result& result
										   );

extern long ARMLOCAL_SABR2B_ImplicitVol(
									  double f, 
									  double K, 
									  double tex,
									  double alpha, 
									  double beta1, 
									  double beta2, 
									  double rho, 
									  double nu,
									  double zero,
									  double lambda,
									  ARM_result& result
									  );

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        SABR  vanilla options Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_SABR_VanillaOption(
									  double f, 
									  double K, 
									  double tex,
									  double alpha, 
									  double beta, 
									  double rho, 
									  double nu,
									  int callput,
									  int flag,
									  int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
									  ARM_result& result
									  );

extern long ARMLOCAL_SABR_FromGaussianToDistribution(
									  double f, 
									  double K, 
									  double tex,
									  double alpha, 
									  double beta, 
									  double rho, 
									  double nu,
									  int flag,
									  int nbsteps,double st1,double st2,double st3,
									  ARM_result& result
									  );

extern long ARMLOCAL_SABR_FromDistributionToGaussian(
									  double f, 
									  double K, 
									  double tex,
									  double alpha, 
									  double beta, 
									  double rho, 
									  double nu,
									  int flag,
									  int nbsteps,double st1,double st2,double st3,
									  ARM_result& result
									  );

extern long ARMLOCAL_SABR_VanillaOption_Der(
										  int i,
										  double f, 
										  double K, 
										  double tex,
										  double alpha, 
										  double beta, 
										  double rho, 
										  double nu,
										  int callput,
										  int flag,
										  int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
										  ARM_result& result
										  );


extern long ARMLOCAL_SABR_VanillaOption_Der2(
										   int i,
										   int j,
										   double f, 
										   double K, 
										   double tex,
										   double alpha, 
										   double beta, 
										   double rho, 
										   double nu,
										   int callput,
										   int flag,
										   int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
										   ARM_result& result
										   );

extern long ARMLOCAL_SABR_From_Sigma_To_Alpha (
							   double f, 
							   double K, 
							   double tex,
							   double sigma, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,
							   ARM_result& result);



/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Barriere Option in BS Model Formula (Single and double barriere)
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_BS_EuroBarriere(
									  double f, 
									  double k, 
									  double b,
									  double r, 
									  double v, 
									  double t, 
									  double d,
									  int callput,
									  int inout,
									  int updown,
									  ARM_result& result
									  );

extern long ARMLOCAL_BS_EuroBarriere_ImpliedVol(
									  double f, 
									  double k, 
									  double b,
									  double r, 
									  double op, 
									  double t, 
									  double d,
									  int callput,
									  int inout,
									  int updown,
									  ARM_result& result
									  );

extern long ARMLOCAL_BS_EuroDoubleBarriere(
									  double f, 
									  double k, 
									  double bup,
									  double bdown, 
									  double v, 
									  double t, 
									  double r,
									  double b,
									  int callput,
									  ARM_result& result
									  );

extern long ARMLOCAL_BS_EuroDoubleBarriere_ImpliedVol(
									  double f, 
									  double k, 
									  double bup,
									  double bdown, 
									  double opt, 
									  double t, 
									  double r,
									  double b,
									  int callput,
									  ARM_result& result
									  );


extern long ARMLOCAL_BS_PartialTime_Start_SingleBarrier(
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bendtime,
									  double t,
									  int    callput,
									  int    optype,
									  ARM_result& result
									  );

extern long ARMLOCAL_BS_PartialTime_End_SingleBarrier(
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bstarttime,
									  double t,
									  int callput,
									  int optype,
									  ARM_result& result
									  );

extern long ARMLOCAL_BS_SingleBarrier_2Asset(
									  double f1,
									  double k1,
									  double f2,
									  double k2, 
									  double v1, 
									  double v2,
									  double corr,
									  double t,
									  int    callput, 
									  int    optype,
									  ARM_result& result
									  );


extern long ARMLOCAL_BS_EuroBarriere_Der(
										 int i,
										 double f, 
										 double k, 
										 double b,
										 double r, 
										 double v, 
										 double t, 
										 double d,
										 int callput,
										 int inout,
										 int updown,
										 ARM_result& result
										 );

extern long ARMLOCAL_BS_EuroDoubleBarriere_Der(
										 int i,
										 double f, 
										 double k, 
										 double bup,
										 double bdown, 
										 double v, 
										 double t, 
										 int callput,
										 ARM_result& result
										 );

extern long ARMLOCAL_BS_PartialTime_Start_SingleBarrier_Der(
									  int i,
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bendtime,
									  double t,
									  int callput,
									  int optype,
									  ARM_result& result
									  );

extern long ARMLOCAL_BS_PartialTime_End_SingleBarrier_Der(
									  int i,
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bstarttime,
									  double t,
									  int callput,
									  int optype,
									  ARM_result& result
									  );

extern long ARMLOCAL_BS_SingleBarrier_2Asset_Der(
									  int i,
									  double f1,
									  double k1,
									  double f2,
									  double k2, 
									  double v1, 
									  double v2,
									  double corr,
									  double t,
									  int    callput, 
									  int    optype,
									  ARM_result& result
									  );




extern long ARMLOCAL_BS_EuroBarriere_Der2(
										 int i,
										 int j,
										 double f, 
										 double k, 
										 double b,
										 double r, 
										 double v, 
										 double t, 
										 double d,
										 int callput,
										 int inout,
										 int updown,
										 ARM_result& result
										 );


extern long ARMLOCAL_BS_EuroDoubleBarriere_Der2(
										 int i,
										 int j,
										 double f, 
										 double k, 
										 double bup,
										 double bdown, 
										 double v, 
										 double t, 
										 int callput,
										 ARM_result& result
										 );

extern long ARMLOCAL_BS_PartialTime_Start_SingleBarrier_Der2(
									  int i,
									  int j,
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bendtime,
									  double t,
									  int callput,
									  int optype,
									  ARM_result& result
									  );

extern long ARMLOCAL_BS_PartialTime_End_SingleBarrier_Der2(
									  int i,
									  int j,
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bstarttime,
									  double t,
									  int callput,
									  int optype,
									  ARM_result& result
									  );

extern long ARMLOCAL_BS_SingleBarrier_2Asset_Der2(
									  int i,
									  int j,
									  double f1,
									  double k1,
									  double f2,
									  double k2, 
									  double v1, 
									  double v2,
									  double corr,
									  double t,
									  int    callput, 
									  int    optype,
									  ARM_result& result
									  );

////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the bivariate function 
///					
////////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_Bivariate(
										 double x,
										 double y,
										 double rho,
										 int p,
										 int q,
										 ARM_result& result
										 );
////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the gamma function 
///					
////////////////////////////////////////////////////////////////////////////////////////////

extern long ARMLOCAL_Gamma(
										 double x,
										 ARM_result& result
										 );
////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the Imcomplete beta function 
///					
////////////////////////////////////////////////////////////////////////////////////////////

extern long ARMLOCAL_ImcompleteBeta(
										 double x,double y,double z,
										 ARM_result& result
										 );

extern long ARMLOCAL_InverseImcompleteBeta(	 double x,double y,double z0,double z,
										 ARM_result& result
										 );

extern long ARMLOCAL_Hypergeometric2F1(	 double x,double y,double z0,double z,
										 ARM_result& result
										 );
////////////////////////////////////////////////////////////////////////////////////////////
///
///			Exportation of	computation of complex error function erf(z) 
///					
////////////////////////////////////////////////////////////////////////////////////////////

extern long ARMLOCAL_RealPart_ComplexErf(
										 double x,
										 double y,
										 int n,
										 ARM_result& result
										 );

extern long ARMLOCAL_ImaginaryPart_ComplexErf(
										 double x,
										 double y,
										 int n,
										 ARM_result& result
										 );

extern long ARMLOCAL_cdfNormal_Inv(
										 double x,
										 ARM_result& result
										 );
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Generalized Heston  Vanilla Options Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_GHeston_VanillaOption(
										double F,
										double K,
										double sig,
										double t,
										double longtermV,
										double theta,
										double ksi,
										double rho,
										double lambda,
										double muJ, 
										double sigmaJ,
										int callorput,int nb,
										ARM_result& result
										);




extern long ARMLOCAL_GHeston_VanillaOption_Der(
											int i,
											double F,
											double K,
											double sig,
											double t,
											double longtermV,
											double theta,
											double ksi,
											double rho,
											double lambda,
											double muJ, 
											double sigmaJ,
											int callorput,
											int nb,
											ARM_result& result
											);


extern long ARMLOCAL_GHeston_VanillaOption_Der2(
											 int i,
											 int j,
											 double F,
											 double K,
											 double sig,
											 double t,
											 double longtermV,
											 double theta,
											 double ksi,
											 double rho,
											 double lambda,
											 double muJ, 
											 double sigmaJ,
											 int callorput,
											 int nb,
											 ARM_result& result
											 );


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread  Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_LogNormal_TriSpreadOption(
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int n,
										   ARM_result& result
										   );

extern long ARMLOCAL_LogNormal_TriSpreadOption_Der(int i,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int n,
										   ARM_result& result
										   );



extern long ARMLOCAL_LogNormal_TriSpreadOption_Der2(int i,int j,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int n,
										   ARM_result& result
										   );



/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread Digital Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption(
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int callput,int n,
										   ARM_result& result
										   );


extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption_Der(int i,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int callput,int n,
										   ARM_result& result
										   );



extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption_Der2(int i,int j,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int callput,int n,
										   ARM_result& result
										   );

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread Digital Option2 Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption2(
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a2,double a3,double b0,double b1,double t,int n,
										   ARM_result& result
										   );

extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption2_Der(int i,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a2,double a3,double b0,double b1,double t,int n,
										   ARM_result& result
										   );


extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption2_Der2(int i,int j,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a2,double a3,double b0,double b1,double t,int n,
										   ARM_result& result
										   );
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Forward value  CEV VanillaOption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_CEV_VanillaOption(double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   );

extern long ARMLOCAL_CEV_VanillaOption_Der(int i,
										double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   );


extern long ARMLOCAL_CEV_VanillaOption_Der2(int i,int j,
										double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   );

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Forward value  CEV DoubleBarrierOption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_CEV_DoubleBarrierOption(double f, double K, double T,double barrdown,double barrup,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   );

extern long ARMLOCAL_CEV_DoubleBarrierOption_Der(int i,
										double f, double K, double T,double barrdown,double barrup,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   );


extern long ARMLOCAL_CEV_DoubleBarrierOption_Der2(int i,int j,
										double f, double K, double T,double barrdown,double barrup,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   );
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Forward value  CEV BarrierOption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_CEV_BarrierOption(double f, double K, double T,double barrier,double drift, double sig, double beta,int optype,int callput, int nbsteps,
										   ARM_result& result
										   );

extern long ARMLOCAL_CEV_BarrierOption_Der(int i,
										double f, double K, double T,double barrier,double drift, double sig, double beta,int optype,int callput, int nbsteps,
										   ARM_result& result
										   );


extern long ARMLOCAL_CEV_BarrierOption_Der2(int i,int j,
										double f, double K, double T,double barrier,double drift, double sig, double beta,int optype,int callput, int nbsteps,
										   ARM_result& result
										   );
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Asian Lognormal vanilla option (geman yor formula)
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_Lognormal_Asian_VanillaOption(double f, double K, double T, double r, double sig,int alpha,int callput, int nbsteps,
										   ARM_result& result
										   );
extern long ARMLOCAL_Lognormal_Asian_VanillaOption_Der(int i,
										double f, double K, double T, double r, double sig,int alpha,int callput, int nbsteps,
										   ARM_result& result
										   );

extern long ARMLOCAL_Lognormal_Asian_VanillaOption_Der2(int i,int j,
										double f, double K, double T, double r, double sig,int alpha,int callput, int nbsteps,
										   ARM_result& result
										   );

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Numerics : Gauss LAgendre Coefficients
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

extern long ARMLOCAL_GaussianIntegrals_Legendre_Coeffs(double a,double b,int n,ARM_result& result
										   );

extern long ARMLOCAL_GaussianIntegrals_Hermite_Coeffs(int n,ARM_result& result
										   );
extern long ARMLOCAL_GaussianIntegrals_Laguerre_Coeffs(double a,int n,ARM_result& result
										   );


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       FX Smile Convertor from delta to strike
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
extern long ARMLOCAL_ConvertFXOptionToStrike( double C_TargetDelta, double C_Fwd, double C_TotalVol, int C_callPut, double C_initValue, ARM_result& result );

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Generalized Heston Pricer with a vector of models
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
extern long ARMLOCAL_GHeston_VanillaOption_ModelVector(
		double C_K,
		double C_t,
		int C_callput,
		int C_Interpolation_Method,
		const vector<double>& C_Maturity_Vec,
		const vector<double>& C_F_Vec,
		const vector<double>& C_InitialVol_Vec,
		const vector<double>& C_longtermV_Vec,
		const vector<double>& C_theta_Vec,
		const vector<double>& C_ksi_Vec,
		const vector<double>& C_rho_Vec,
		const vector<double>& C_lambda_Vec,
		const vector<double>& C_muJ_Vec,
		const vector<double>& C_sigmaJ_Vec,
		int C_nb,
		ARM_result&  C_result);

extern long ARMLOCAL_GHeston_Implicit_Volatility_ModelVector(
		double C_K,
		double C_t,
		int C_Interpolation_Method,
		const vector<double>& C_Maturity_Vec,
		const vector<double>& C_F_Vec,
		const vector<double>& C_InitialVol_Vec,
		const vector<double>& C_longtermV_Vec,
		const vector<double>& C_theta_Vec,
		const vector<double>& C_ksi_Vec,
		const vector<double>& C_rho_Vec,
		const vector<double>& C_lambda_Vec,
		const vector<double>& C_muJ_Vec,
		const vector<double>& C_sigmaJ_Vec,
		int C_nb,
		ARM_result&  C_result);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Calibration of a SABR Model to a curve of implied vol
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

extern long ARMLOCAL_SABR_Model_Calibrate(
										  const vector<double>& C_K_Vec,
										  const vector<double>& C_ImpVol_Vec,
										  const vector<double>& C_Weigth_Vec,
										  double f,
										  double t,
										  int flag,
										  int nb,int algorithm,
										  double alpha0,double beta0,double rho0,double nu0,
										  ARM_result&  result);

extern long ARMLOCAL_SABR_Model_Calibrate_WForward(
												   const vector<double>& C_K_Vec,
												   const vector<double>& C_ImpVol_Vec,
												   const vector<double>& C_Weigth_Vec,
												   double t,
												   int flag,
												   int nb,int algorithm,
												   double alpha0,double beta0,double rho0,double nu0,double f0,
												   ARM_result&  result);

long ARMLOCAL_SABR_Model_Calibrate_FixedBeta(
											 const vector<double>& C_K_Vec,
											 const vector<double>& C_ImpVol_Vec,
											 const vector<double>& C_Weigth_Vec,
											 double f,
											 double beta,
											 double t,
											 int flag ,double Alpha_Exp,double Alpha_Tanh,double Kb_Tanh,
											 int nbsteps,int algorithm,
											 double alpha0,double rho0,double nu0,
											 ARM_result&  result);

long ARMLOCAL_SABR_Model_Calibrate_FixedBeta_Linked(
													const vector<double>& C_K_Vec,
													const vector<double>& C_ImpVol_Vec,
													const vector<double>& C_Weigth_Vec,
													double f,
													double beta,
													double t,
													int flag ,int nbsteps,int algorithm,
													double alpha0,double rho0,double nu0,
													double alphap,double rhop, double nup,double rweight_alpha,double rweight_rho,double rweight_nu,
													ARM_result&  result);

long ARMLOCAL_SABR_Calibrate_BetaFixedToOne(
											 const vector<double>& C_K_Vec,
											 const vector<double>& C_ImpVol_Vec,
											 const vector<double>& C_Weigth_Vec,
											 double f,
											 double t,
											 double atmvol,
											 ARM_result&  result);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Calibration of a GHeston Model to a curve of implied vol
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

long ARMLOCAL_GHeston_Model_Calibrate_Total(
									   const vector<double>& C_K_Vec,
									   const vector<double>& C_ImpVol_Vec,
									   double C_f,
									   double C_t,
									   int C_nb,
									   int C_algorithm,
									   double C_V0,
									   double C_omega,
									   double C_theta,
									   double C_ksi,
									   double C_rho,
									   double C_muJ,
									   double C_sigmaJ,
									   double C_lambda,
									   ARM_result&  result);

long ARMLOCAL_GHeston_Model_Calibrate_NoJump(
			const vector<double>& C_K_Vec,
			const vector<double>& C_ImpVol_Vec,
			double C_f,
			double C_t,
			int C_nb,
			int C_algorithm,
			double C_V0,
			double C_omega,
			double C_theta,
			double C_ksi,
			double C_rho,
			double C_muJ0,
			double C_sigmaJ0,
			double C_lambda0,
			ARM_result&  result);

////////////////////////////////////////////
/// Spread Option SABR Student
////////////////////////////////////////////

extern long ARMLOCAL_Student_SABR_Power_SpreadOption(
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	const vector<double>&  C_Copula_vec,
	double t,
	int flag,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	);

extern long ARMLOCAL_Student_SABR_Power_Digital_SpreadOption(
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double copula_degre,
	double t,
	int flag,
	double k10,
	double a20,
	double b20,
	double k20,
	int n,
	double alpha_exp,
	double alpha_tanh,
	double kb_tanh,
	ARM_result& result
	);


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Mepi Vanilla Option 
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_Mepi_VanillaOption_STOBS(
	double date,
	double P0,
	double K,
	double T,
	int zero_curve_ptr,
	double b,
	double YearlyFees,
	double cashspread,
	int sabr_model_ptr,
	double minExp,
	double maxExp,
	double riskFac,
	double g0,
	double g,
	int AveragingPeriodNb,
	double reset,
	int	   callOrPut,
	int	   n,
	int nz,
	int nh,
	ARM_result& result
	);

extern long ARMLOCAL_Mepi_VanillaOption_STOBS_delta(
	double date,
	double P0,
	double K,
	double T,
	int zero_curve_ptr,
	double b,
	double YearlyFees,
	double cashspread,
	int sabr_model_ptr,
	double minExp,
	double maxExp,
	double riskFac,
	double g0,
	double g,
	int AveragingPeriodNb,
	double reset,
	int	   callOrPut,
	int	   n,
	int nz,
	int nh,
	ARM_result& result
	);

extern long ARMLOCAL_Mepi_VanillaOption_SABR(
	double date,
	double P0,
	double K,
	double T,
	int zero_curve_ptr,
	double b,
	double YearlyFees,
	double cashspread,
	int sabr_model_ptr,
	CCString C_SpotName,
	double minExp,
	double maxExp,
	double riskFac,
	double g0,
	double g,
	int AveragingPeriodNb,
	double reset,
	int	   callOrPut,
	int	   n,
	int nz,
	int nh,
	ARM_result& result
	);

extern long ARMLOCAL_Mepi_VanillaOption_SABR_delta(
	double date,
	double P0,
	double K,
	double T,
	int zero_curve_ptr,
	double b,
	double YearlyFees,
	double cashspread,
	int sabr_model_ptr,
	CCString C_SpotName,
	double minExp,
	double maxExp,
	double riskFac,
	double g0,
	double g,
	int AveragingPeriodNb,
	double reset,
	int	   callOrPut,
	int	   n,
	int nz,
	int nh,
	ARM_result& result
	);



/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Present value  of a VanillaOption on an underlying with stochastic volatility
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_StochasticVol_LN_VanillaOption(
													double f,
													double K, 
													double T,
													double drift, 
													double sig, 
													double VolDrift,
													double VolVol,double averaging,double reset,
													int callput,
													int nbsteps,
													ARM_result& result
													);

extern long ARMLOCAL_StochasticVol_LN_VanillaOption_Der(
													int i,
													double f,
													double K, 
													double T,
													double drift, 
													double sig, 
													double VolDrift,
													double VolVol,double averaging,double reset,
													int callput,
													int nbsteps,
													ARM_result& result
													);

extern long ARMLOCAL_StochasticVol_LN_VanillaOption_Der2(
													int i,
													int j,
													double f,
													double K, 
													double T,
													double drift, 
													double sig, 
													double VolDrift,
													double VolVol,double averaging,double reset,
													int callput,
													int nbsteps,
													ARM_result& result
													);

extern long ARMLOCAL_StochasticVol_LN_Ari_VanillaOption(
													double f,
													double K, 
													double T,
													double drift, 
													double sig, 
													double VolDrift,
													double VolVol,double averaging,double reset,
													int callput,
													int nbsteps,
													ARM_result& result
													);

extern long ARMLOCAL_StochasticVol_LN_Ari_VanillaOption_Der(
													int i,
													double f,
													double K, 
													double T,
													double drift, 
													double sig, 
													double VolDrift,
													double VolVol,double averaging,double reset,
													int callput,
													int nbsteps,
													ARM_result& result
													);

extern long ARMLOCAL_StochasticVol_LN_Ari_VanillaOption_Der2(
													int i,
													int j,
													double f,
													double K, 
													double T,
													double drift, 
													double sig, 
													double VolDrift,
													double VolVol,double averaging,double reset,
													int callput,
													int nbsteps,
													ARM_result& result
													);




extern long ARMLOCAL_OptimizeSkewVector(
	const vector<double>&	C_tagetVect,
	const vector<double>&   C_weights,
    const vector<double>&   C_presicions,
	const vector<double>&   C_InitVector,
	const vector<double>&   C_LBoundVector,
	const vector<double>&   C_UBoundVector,
	vector<double>&   C_DataResult,
	const long&				C_algo,
	ARM_result&				result);


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Present value  of a Normal VanillaOption on an underlying with stochastic volatility
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
extern long ARMLOCAL_StochasticVol_N_VanillaOption(double f, double K, double T,double drift, double sig, double VolDrift,double VolVol,int callput, int nbsteps,
										   ARM_result& result
										   );

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       GLambda Distribution
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

long ARMLOCAL_GLambda_From_SABR_Calibrate(

	double  C_forward,
	double	C_alpha,
	double	C_beta,
	double	C_rho,
	double	C_nu,
	double	C_T,
	double	C_SabrType,
	double	C_Scope,
	double	C_IniL1,
	double	C_IniL2,
	double	C_IniL3,
	double	C_IniL4,
	double	C_IniL5,
	double	C_IniL6,
	double	C_Algo,
	ARM_result&  result);



extern long ARMLOCAL_Student_GLambda_Power_SpreadOption(
	double l1a,
	double l2a,
	double l3a,
	double l4a,
	double l5a,
	double l6a,

	double l1b,
	double l2b,
	double l3b,
	double l4b,
	double l5b,
	double l6b,

	double copula_corr,
	double copula_degre,

	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	);

extern long ARMLOCAL_Student_GLambda_Power_Digital_SpreadOption(
	double l1a,
	double l2a,
	double l3a,
	double l4a,
	double l5a,
	double l6a,

	double l1b,
	double l2b,
	double l3b,
	double l4b,
	double l5b,
	double l6b,

	double copula_corr,
	double copula_degre,

	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	);

extern long ARMLOCAL_Student_GLambda_Power_Index1Digital_SpreadOption(
	double l1a,
	double l2a,
	double l3a,
	double l4a,
	double l5a,
	double l6a,

	double l1b,
	double l2b,
	double l3b,
	double l4b,
	double l5b,
	double l6b,

	double copula_corr,
	double copula_degre,

	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	);
extern long ARMLOCAL_Student_GLambda_Power_Index2Digital_SpreadOption(
	double l1a,
	double l2a,
	double l3a,
	double l4a,
	double l5a,
	double l6a,

	double l1b,
	double l2b,
	double l3b,
	double l4b,
	double l5b,
	double l6b,

	double copula_corr,
	double copula_degre,

	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	);
extern long ARMLOCAL_Student_GLambda_SpreadOption(
	double l1a,
	double l2a,
	double l3a,
	double l4a,
	double l5a,
	double l6a,

	double l1b,
	double l2b,
	double l3b,
	double l4b,
	double l5b,
	double l6b,

	double copula_corr,
	double copula_degre,

	double k20,
	int n,
	ARM_result& result
	);

////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the Lambert function 
///					
////////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_Lambert_Function(
										 double x,
										 ARM_result& result
										 );

////////////////////////////////////////////////////////////////////////////////////////
///
///			Black and Sholes Asymptotic time value Formula
///
////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_CF_BlackSholesTimeValue(
								 double forward,
								 double totalvolatility,
								 double bondprice,
								 double strike,
								 double CallPut,
								 ARM_result& result
	);

extern long ARMLOCAL_CF_BlackSholesTimeValue_ImplicitVol(
												double forward,
												double bondprice,
												double strike,
												double CallPut,
												double optprice,
												ARM_result& result
												);

extern long ARMLOCAL_CF_ShiftedLogNormal_Quantile(
												double f,
												double k,
												double t,
												double sigma,
												double alpha,
												ARM_result& result
												);

extern long ARMLOCAL_CF_ShiftedLogNormal_Distribution(
												double f,
												double k,
												double t,
												double sigma,
												double alpha,
												ARM_result& result
												);

extern long ARMLOCAL_CF_SABR_Quantile(
												double f,
												double k,
												double t,
												double alpha,
												double beta,
												double rho,
												double nu,
												double Sabr_Type,
												double n,
												double alpha_exp,
												double alpha_tanh,
												double kb_tanh,
												ARM_result& result
												);

extern long ARMLOCAL_CF_SABR_Distribution(
												double f,
												double k,
												double t,
												double alpha,
												double beta,
												double rho,
												double nu,
												double Sabr_Type,
												double n,
												double alpha_exp,
												double alpha_tanh,
												double kb_tanh,
												ARM_result& result
												);
extern long ARMLOCAL_cdfNormal(
										 double x,
										 ARM_result& result
										 );

extern long ARMLOCAL_GLambda_Distribution(
										  double l1,
										  double l2,
										  double l3,
										  double l4,
										  double l5,
										  double l6,
										  double x,
										  ARM_result& result
										  );


extern long ARMLOCAL_GLambda_Quantile(
										  double l1,
										  double l2,
										  double l3,
										  double l4,
										  double l5,
										  double l6,
										  double x,
										  ARM_result& result
										  );

extern long ARMLOCAL_Student_Quantile(
										  double deg,
										  double x,
										  ARM_result& result
										  );

extern long ARMLOCAL_Student_Distribution(
										  double deg,
										  double x,
										  ARM_result& result
										  );

extern long ARMLOCAL_ImcompleteBeta_Inverse(
										  double a,
										  double b,
										  double x,
										  ARM_result& result
										  );


extern long ARMLOCAL_Student_QIntegral(
										  double a,
										  double b,
										  double x,
										  ARM_result& result
										  );

extern long ARMLOCAL_Normal_ImpliedVol(
										  double a,
										  double b,
										  double x,
										  double cp,
										  ARM_result& result
										  );

extern long ARMLOCAL_Normal_Digital_ImpliedVol(
										  double a,
										  double b,
										  double x,
										  double callput,
										  ARM_result& result
										  );

extern long ARMLOCAL_Hypergeometric_Whittaker_W(
										  double a,
										  double b,
										  double x,
										  ARM_result& result
										  );
extern long ARMLOCAL_Hypergeometric_Whittaker_M(
										  double a,
										  double b,
										  double x,
										  ARM_result& result
										  );

extern long ARMLOCAL_Bessel_Y(
										  double a,
										  double x,
										  ARM_result& result
										  );



extern long ARMLOCAL_Bessel_I(
										  double a,
										  double x,
										  ARM_result& result
										  );

extern long ARMLOCAL_Bessel_J(
										  double a,
										  double x,
										  ARM_result& result
										  );

extern long ARMLOCAL_Bessel_K(
										  double a,
										  double x,
										  ARM_result& result
										  );

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Shifted Heston Vanilla Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_Shifted_Heston_VanillaOption(
										   double F,
										   double K,
										   double sig,
										   double t,
										   double longtermV,
										   double theta,
										   double ksi,
										   double rho,
										   double shift,
										   int callorput,
										   int nb1,int nb,int NbS, int NbO,double prec,
										   ARM_result& result
										   );

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        SABR Heston Vanilla Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_SABR_Heston_VanillaOption(
										   double F,
										   double K,
										   double sig,
										   double t,
										   double longtermV,
										   double theta,
										   double ksi,
										   double rho,
										   double beta,
										   int callorput,
										   int nb1,int nb,int NbS, int NbO,double prec,
										   ARM_result& result
										   );

/////////////////////////////////////////////////////////////////////////////////////
///
///		Glambda spreadoption Calibration
///
/////////////////////////////////////////////////////////////////////////////////////

long ARMLOCAL_GLambda_CompleteSpreadoption(
										   const vector<double>& C_l1a_Vec,
										   const vector<double>& C_l2a_Vec,
										   const vector<double>& C_l3a_Vec,
										   const vector<double>& C_l4a_Vec,
										   const vector<double>& C_l5a_Vec,
										   const vector<double>& C_l6a_Vec,
										   const vector<double>& C_l1b_Vec,
										   const vector<double>& C_l2b_Vec,
										   const vector<double>& C_l3b_Vec,
										   const vector<double>& C_l4b_Vec,
										   const vector<double>& C_l5b_Vec,
										   const vector<double>& C_l6b_Vec,
										   const vector<double>& C_Discount_Vec,
										   double C_copula_corr,
										   double C_copula_degre,
										   double C_k,
										   double C_n,
										   ARM_result&  result);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Jump Diffusion Mepi Vanilla Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_JumpDiffusion_Mepi_Call(
											   double   C_P0,
											   double	C_f0,
											   double	C_T,
											   double	C_K,
											   double	C_R,
											   double	C_Emin,
											   double	C_Lmax,
											   double	C_gamma0,
											   double	C_gamma1,
											   double	C_sig,
											   double	C_lambda,
											   double	C_sigJ,
											   double	C_r,
											   double	C_s,
											   double	C_mu,
											   double	C_fees,
											   double	C_volDrift,
											   double	C_volVol,
											   int	C_allPut,
											   const vector<double>& C_params,
											   ARM_result& C_result);


extern long ARMLOCAL_Util_TrigonalSolve(
										const vector<double>& C_A_Vec,
										const vector<double>& C_B_Vec,
										const vector<double>& C_C_Vec,
										const vector<double>& C_R_Vec,
										ARM_result& result);


/////////////////////////////////////////////////////////////////////////////////////
///
///	BiSABR Vanilla Spreadoption Formula
///
/////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_BiSABR_SpreadOption(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double   C_flag,
											   ARM_result& result);


extern long ARMLOCAL_Hypergeometric_Appell(
										   double   C_a,
										   double	C_b1,
										   double	C_b2,
										   double	C_c,
										   double	C_x,
										   double	C_xim,
										   double	C_y,
										   double	C_yim,
										   double	C_nb,
										   ARM_result& result);

extern long ARMLOCAL_Util_Eigenvalues4(
										const double rho12,
										const double rho13,
										const double rho14,
										const double rho23,
										const double rho24,
										const double rho34,
										ARM_result& result);

extern long ARMLOCAL_Util_Eigenvalues3(
										const double rho12,
										const double rho13,
										const double rho23,
										ARM_result& result);


extern long ARMLOCAL_Util_BiSABR_CorrelationEvolution(
										const double rho1,
										const double rho2,
										const double rhos,
										const double rhov,
										const double rhoc12,
										const double rhoc21,
										const double newrho1,
										const double newrho2,
										ARM_result& result);



extern long ARMLOCAL_BiSABR_Calibrate(
										const vector<double>& C_F1_vec		, 
										const vector<double>& C_alpha1_vec	,
										const vector<double>& C_beta1_vec	,
										const vector<double>& C_rho1_vec	,
										const vector<double>& C_nu1_vec		,
										const vector<double>& C_F2_vec		,
										const vector<double>& C_alpha2_vec	,
										const vector<double>& C_beta2_vec	,
										const vector<double>& C_rho2_vec	,
										const vector<double>& C_nu2_vec		,
										const vector<double>& C_strike_vec	,
										const vector<double>& C_maturity_vec,
										const vector<double>& C_callput_vec	,
										const vector<double>& C_price_vec	,
										const vector<double>& C_weight_vec	,
										const vector<double>& C_initialparams,
										ARM_result& result);

extern long ARMLOCAL_LN_DigitalOption(
											   double   C_forward,
											   double	C_strike,
											   double	C_maturity,
											   double	C_CallPut,
											   double	C_volatility,
											   ARM_result& result);




extern long ARMLOCAL_LN_RatioOption(
											   double   C_S1,
											   double	C_Mu1,
											   double	C_Sigma1,
											   double   C_S2,
											   double	C_Mu2,
											   double	C_Sigma2,
											   double	C_Rho,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   ARM_result& result);

extern long ARMLOCAL_LN_ProductOption(
											   double   C_S1,
											   double	C_Mu1,
											   double	C_Sigma1,
											   double   C_S2,
											   double	C_Mu2,
											   double	C_Sigma2,
											   double	C_Rho,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   ARM_result& result);




extern long ARMLOCAL_SABR_GaussianSABRDigitalCall(
		double C_f1,
		double C_alpha1,
		double C_beta1,
		double C_rho1,
		double C_nu1,
		double C_flag1,
		double C_f2,
		double C_alpha2,
		double C_beta2,
		double C_rho2,
		double C_nu2,
		double C_flag2,
		double C_rho,
		double C_K,
		double C_T,
		double C_LegendreNb,double C_alpha_exp,double  C_alpha_tanh,double  C_kb_tanh,
		ARM_result& C_result);

extern long ARMLOCAL_SABR_GaussianSABRDigitalCallPayingS1(
		double C_f1,
		double C_alpha1,
		double C_beta1,
		double C_rho1,
		double C_nu1,
		double C_flag1,
		double C_f2,
		double C_alpha2,
		double C_beta2,
		double C_rho2,
		double C_nu2,
		double C_flag2,
		double C_rho,
		double C_K,
		double C_T,
		double C_LegendreNb,double C_alpha_exp,double  C_alpha_tanh,double  C_kb_tanh,
		ARM_result& C_result);


extern long ARMLOCAL_SABR_GaussianSABRDigitalCallPayingS2(
		double C_f1,
		double C_alpha1,
		double C_beta1,
		double C_rho1,
		double C_nu1,
		double C_flag1,
		double C_f2,
		double C_alpha2,
		double C_beta2,
		double C_rho2,
		double C_nu2,
		double C_flag2,
		double C_rho,
		double C_K,
		double C_T,
		double C_LegendreNb,double C_alpha_exp,double  C_alpha_tanh,double  C_kb_tanh,
		ARM_result& C_result);

extern long ARMLOCAL_SABR_GaussianSABRDigitalCallPayingS3(
		double C_f1,
		double C_alpha1,
		double C_beta1,
		double C_rho1,
		double C_nu1,
		double C_flag1,
		double C_f2,
		double C_alpha2,
		double C_beta2,
		double C_rho2,
		double C_nu2,
		double C_flag2,
		double C_f3,
		double C_sigma3,
		const vector<double>& C_Correlations_Vec,
		double C_K,
		double C_T,
		const vector<double>& C_Params_Vec,
		ARM_result& C_result);

//////////////////////////////////////////////////
//// Function to create a ARM_SmiledFRM
//////////////////////////////////////////////////
extern long ARMLOCAL_TarnProxy_Create(
	const vector<double>&	resetDates,
	const vector<double>&	fwds,
	const VECTOR<long>&		densityFunctorId,
	const vector<double>&	df,
	const vector<double>&	levPrec,
	const vector<double>&	lev,
	const vector<double>&	fix,
	const vector<double>&	cap,
	const vector<double>&	floor,
	const vector<double>&	fees,
	const vector<double>&	dcf,
	const double&			target,
	const bool&				globalcap,
	const bool&				globalfloor,
	const double&			C_CorrelInput,
	const int&				C_NbSimul,
	ARM_result&		result, 
	long			objId );

extern long ARMLOCAL_TarnProxy_GetPrice(
	const long&	C_TarnProxyId, 
	const double& C_NbPayoff, 
	ARM_result&	result);

extern long ARMLOCAL_VBMinMaxProxy_Create(
	const double&			C_AsOf,
	const vector<double>&	C_resetDates,
	const vector<double>&	C_fwdRates,
	const vector<double>&	C_totalVol,
	const vector<double>&	C_leftVol,
	const vector<double>&	C_rightVol,
	const vector<double>&	C_nu,
	const vector<double>&	C_rho,
	const int&				C_nbSimul,
	const bool&				C_sabrDiff,
	const int&				C_typeprice,
	const int&				C_type1sens,
	const int&				C_type2sens,
	const vector<double>&	C_rate1Lev,
	const vector<double>&	C_rate1Add,
	const vector<double>&	C_capRateLev,
	const vector<double>&	C_capRateAdd,
	const vector<double>&	C_vbLev,
	const int&				C_maxChoice,
	const int&				C_minChoice,
	const double&			C_minmaxFreq,
	ARM_result&				result,
	long					objId);

extern long ARMLOCAL_VBMinMaxProxy_GetInfo(
	const long& VBMinMaxProxyId,
	const int& info,
	vector<double>&	vresult,
	ARM_result& result
	);

extern long ARMLOCAL_Berm2DatesProxy_Create(
	const double&			C_AsOf,
	const vector<double>&	C_resetDates,
	const vector<double>&	C_fwdRates,
	const double&			C_strike,
	const vector<double>&	C_DFs,
	const vector<double>&	C_fwdRatesVol,
	const vector<double>&	C_fwdRatesPartVol,
	const vector<double>&	C_vols,
	const vector<double>&	C_volvols,
	const vector<double>&	C_rho,
	const double&			C_Rho12,
	const int&				nbSimul,
	const int&				C_TypeDiff,
	ARM_result&				result,
	long					objId);

extern long ARMLOCAL_Berm2DatesProxy_GetPrice(
	const long& Berm2DatesProxyId,
	const int& info,
	ARM_result& result
	);

extern long ARMLOCAL_SpreadVBProxy_Create(
	const double&			C_AsOf,
	const vector<double>&	C_resetDates,
	const vector<double>&	C_fwdRates1,
	const vector<double>&	C_fwdRatesVol1,
	const vector<double>&	C_fwdRatesPartVol1,
	const vector<double>&	C_vols1,
	const vector<double>&	C_nu1,
	const vector<double>&	C_rho1,
	const vector<double>&	C_fwdRates2,
	const vector<double>&	C_fwdRatesVol2,
	const vector<double>&	C_fwdRatesPartVol2,
	const vector<double>&	C_vols2,
	const vector<double>&	C_nu2,
	const vector<double>&	C_rho2,
	const vector<double>&	C_Rate1Rate2Correl,
	const int&				nbSimul,
	const bool&				sabrDiff,
	const vector<double>&	C_fwdLev,
	const vector<double>&	C_fwdStrikes,
	const int&				C_typeprice,
	const int&				C_type1sens,
	const int&				C_type2sens,
	const double&			C_CorrMeanRev,
	const double&			C_CorrVol,
	const vector<double>&	C_Levier,
	const vector<double>&	C_Fixed,
	ARM_result&				result,
	long					objId);

extern long ARMLOCAL_SpreadVBProxy_GetInfo(
	const long& SpreadVBProxyId,
	const int& info,
	vector<double>&	vresult,
	ARM_result& result
	);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        BiSABR Digital Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_BiSABR_Digital_SpreadOption(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double	C_flag,
											  
											   ARM_result& result);


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        BiSABR Digital Pays S1 Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_BiSABR_DigitalPaysS1_SpreadOption(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double	C_flag,
											  
											   ARM_result& result);




/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        BiSABR Digital Pays S2 Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_BiSABR_DigitalPaysS2_SpreadOption(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double	C_flag,
											  
											   ARM_result& result);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        BiSABR Digital Pays S3 Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

extern long ARMLOCAL_BiSABR_DigitalPaysS3_SpreadOption(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   const vector<double>& C_S3params,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   double	C_flag,
											  
											   ARM_result& result);

extern long ARMLOCAL_BiSABR_Distribution(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double	C_K,
											   double	C_T,
											   double	C_flag,
											  
											   ARM_result& result);

extern long ARMLOCAL_BiSABR_Quantile(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double	C_K,
											   double	C_T,
											   double	C_flag,
											  
											   ARM_result& result);


extern long ARMLOCAL_BetaEqualZeroSABR(
											   double   C_F,
											   double	C_K,
											   double	C_T,
											   double	C_mu,
											   double	C_alpha,
											   double	C_rho,
											   double	C_nu,
											   double	C_CallPut,
											   ARM_result& result);


extern long ARMLOCAL_Shifted2LogNormal_Distribution(
											   double   C_F1,
											   double	C_sigma1,
											   double	C_F2,
											   double	C_sigma2,
											   double	C_alpha,
											   double	C_rho,
											   double	C_K,
											   double	C_T,
											   double	C_n,
											   ARM_result& result);

extern long ARMLOCAL_Shifted2LogNormal_Quantile(
											   double   C_F1,
											   double	C_sigma1,
											   double	C_F2,
											   double	C_sigma2,
											   double	C_alpha,
											   double	C_rho,
											   double	C_K,
											   double	C_T,
											   double	C_n,
											   ARM_result& result);

extern long ARMLOCAL_BiSABR_S3_SpreadOption(
											const vector<double>& C_S1params_vec,
											const vector<double>& C_S2params_vec,
											const vector<double>& C_S3params_vec,
											double	C_rhos,
											double	C_rhov,
											double	C_rhoc12,
											double	C_rhoc21,
											double	C_Correlation,
											double	C_T,
											double	C_A1,
											double	C_B1,
											double	C_K1,
											double	C_A2,
											double	C_B2,
											double	C_K2,
											double	C_flag,
											double  C_nbsteps,
											ARM_result& result);



extern long ARMLOCAL_Heston_OptionPrice(
												double	C_AsOfDate,
												double	C_ResetTime,
												double	C_Forward,
												double	C_Strike,
												int		C_CallOrPut,
												double	C_V0,
												double	C_Kappa,
												double  C_Rho,
												double  C_Theta,
												double	C_VVol,
												double	C_Shift,
												const vector<double>& C_Times,
												const vector<double>& C_Level,
												ARM_result& result );

extern long ARMLOCAL_Heston2b_OptionPrice(
												double	C_AsOfDate,
												double	C_ResetTime,
												double	C_Forward,
												double	C_Strike,
												int		C_CallOrPut,
												double	C_V01,
												double	C_Kappa1,
												double  C_Rho1,
												double  C_Theta1,
												double	C_VVol1,
												double	C_V02,
												double	C_Kappa2,
												double  C_Rho2,
												double  C_Theta2,
												double	C_VVol2,
												double	C_Shift,
												const vector<double>& C_Times,
												const vector<double>& C_Level,
												ARM_result& result );

extern long ARMLOCAL_MixteHeston_OptionPrice(
	const double&	asOf,
	const double&	resetDate,
	const double&	forward,
	const double&	strike,
	const int&		callPut,
	const double&	sigma,
	const double&	v0,
	const double&	kappa,
	const double&	rho,
	const double&	theta,
	const double&	vvol,
	const double&	shift,
	const vector<double>& times,
	const vector<double>& levels,
	ARM_result& result );

extern long ARMLOCAL_Normal_Heston_VanillaCall(
											   double   C_rho,
											   double	C_lambdaV,
											   double	C_thetaV,
											   double	C_kappaV,
											   double	C_V0,
											   double	C_S0,
											   double	C_k,
											   double	C_T,
											   double	C_lambdaB,
											   double	C_callput,
											   double	C_nbfirst,
											   double	C_nb,
											   double	C_NbStage,
											   double	C_NbOscill,
											   double	C_prec,
											   ARM_result& result);

extern long ARMLOCAL_SuperNormal_Heston_VanillaCall(
	const double& C_rho1,
	const double& C_lambda1,
	const double& C_theta1,
	const double& C_kappa1,
	const double& C_V01,
	const double& C_rho2,
	const double& C_lambda2,
	const double& C_theta2,
	const double& C_kappa2,
	const double& C_V02,
	const double& C_S0,
	const double& C_k,
	const double& C_T,
	const double& C_CallPut,
	const double& C_nb,
	ARM_result& result);


extern long ARMLOCAL_SABR_To_Heston_SmileCalibration_Create(
	const double&  C_resetTime,
	const double&  C_fwdRate,
	const double&  C_ATMVol,
	const double&  C_Alpha,
	const double&  C_Beta,
	const double&  C_RhoSABR,
	const double&  C_Nu,
	const double&  C_Sabr_Type,
	const double&  C_InitialVar,
	const double&  C_Kappa,
	const double&  C_Theta,
	const double&  C_Rho,
	const double&  C_Shift,
	const bool& C_CalibTheta,
	const bool& C_CalibRho,
	const bool& C_CalibShift,
	ARM_result&				result,
	long					objId);

extern long ARMLOCAL_SABR_To_Heston_SmileCalibration_GetValue(
	const long& CalibrationId,
	const int& info,
	ARM_result& result
	);
	
class ARM_Heston_CalibrateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_Heston_CalibrateFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_SABR_SmileCalibration_ParamFunctor : public ARM_GenericAddinFunctor 
{
public:
	ARM_SABR_SmileCalibration_ParamFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_SABR2B_SmileCalibration_ParamFunctor : public ARM_GenericAddinFunctor 
{
public:
	ARM_SABR2B_SmileCalibration_ParamFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_Heston_SmileCalibration_ParamFunctor : public ARM_GenericAddinFunctor 
{
public:
	ARM_Heston_SmileCalibration_ParamFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_Heston2b_SmileCalibration_ParamFunctor : public ARM_GenericAddinFunctor 
{
public:
	ARM_Heston2b_SmileCalibration_ParamFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_BiSABR_SmileCalibration_ParamFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_BiSABR_SmileCalibration_ParamFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_Merton_SmileCalibration_ParamFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_Merton_SmileCalibration_ParamFunctor() {}

	virtual long operator()( ARM_result& result, long objId );
};

extern long ARMLOCAL_SmileCalibration(
	double C_AsOfDate,
	const vector<double>& C_CalibTimes,
	const vector<double>& C_Forwards,
	const vector< vector<double> >& C_MktVols,
	const vector< vector<double> >& C_Strikes,
	const long& C_CalibParamId,
	const vector<double>& C_ConstraintStrikes,
	const vector<double>& C_ConstraintVols,
	const bool& C_CalibConstraint,
	const vector<double>& C_Weights,
	int& rowsResult,
	int& colsResult,
	vector<double>& outResult,
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& outResultBool,
	ARM_result& result );

extern long ARMLOCAL_SpreadSmileCalibration2Heston(
	double C_ResetTime,
	double C_Fwd1,
	const vector<double>& C_MktVols1,
	const vector<double>& C_Strikes1,
	double C_ConstrVol1,
	double C_ConstrK1,
	double C_Fwd2,
	const vector<double>& C_MktVols2,
	const vector<double>& C_Strikes2,
	double C_ConstrVol2,
	double C_ConstrK2,
	const vector<double>& C_MktVolsSpread,
	const vector<double>& C_StrikesSpread,
	double C_ConstrVolSpread,
	double C_ConstrKSpread,
	const double& C_v0,
	const double& C_kappa,
	const double& C_theta,
	const bool& C_calibTheta,
	int& rowsResult,
	int& colsResult,
	vector<double>& outResult,
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& outResultBool,
	ARM_result& result);

extern long ARMLOCAL_Spread2HestonVanilla(
	double C_ResetTime,
	double C_Fwd1,
	double C_Fwd2,
	double C_Strike,
	double C_CallPut,
	double C_V0,
	double C_Kappa,
	double C_Theta,
	double C_Nu,
	double C_Rho1,
	double C_Rho2,
	double C_Shift1,
	double C_Shift2,
	double C_Level1,
	double C_Level2,
	double C_Correl,
	double C_Index1Lev,
	double C_Index2Lev,
	ARM_result& result);

extern long ARMLOCAL_Spread2Heston_TOTEMCalibration(
	const vector<double>& TOTEMMat,
	const vector<double>& TOTEMStrikes,
	const vector<double>& TOTEMPrices,
	const vector<double>& FullScheduleReset,
	const vector<double>& FullScheduleAnnuity,
	const vector<double>& FullScheduleFwd1,
	const vector<double>& FullScheduleFwd2,
	const vector<double>& FwdCalibReset,
	const vector<double>& LongFwds,
	const long& LongVolsId,
	const long& LongStrikesId,
	const vector<double>& ShortFwds,
	const long& ShortVolsId,
	const long& ShortStrikesId,
	const double& v0,
	const double& kappa,
	const double& theta,
	const bool& constrCorrel,
	int& nbRowsRes,
	int& nbColsRes,
	vector<double>& Res,
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& outResBool,
	ARM_result& result
	);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        TriSABR Vanilla Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
extern long ARMLOCAL_TRiSABR_VanillaOption(
											const vector<double>& C_S1params_vec,
											const vector<double>& C_S2params_vec,
											const vector<double>& C_S3params_vec,
											double	C_rhos12,
											double	C_rhos23,
											double	C_rhos13,
											double	C_rhov12,
											double	C_rhov23,
											double	C_rhov13,
											double	C_rhoc12,
											double	C_rhoc21,
											double	C_rhoc23,
											double	C_rhoc32,
											double	C_rhoc13,
											double	C_rhoc31,
											double	C_K,
											double	C_T,	
											double	C_callput,
											double	C_flag,
											ARM_result& result);

extern long ARMLOCAL_Util_TriSABR_Eigenvalues(
			double rho1,   double rho2,   double rho3,
			double rhos12, double rhos23, double rhos13,
			double rhov12, double rhov23, double rhov13,
			double rhoc12, double rhoc13,
			double rhoc21, double rhoc23,
			double rhoc31, double rhoc32,
			ARM_result& result);


extern long ARMLOCAL_Nonparametric_CompleteOption(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_Strike2_vec,
											const vector<double>& C_Vol2_vec,
											const vector<double>& C_S1params_vec,
											const vector<double>& C_S2params_vec,
											double	C_correlation,
											double	C_maturity,
											double	C_a1,
											double	C_b1,
											double	C_k1,
											double	C_a2,
											double	C_b2,
											double	C_k2,
											double	C_nbsteps,
											double	C_algo,
											double  C_smiletype,
								
											ARM_result& result);
	

extern long ARMLOCAL_Nonparametric_LogVolatility(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
											double	C_k,
											ARM_result& result);

extern long ARMLOCAL_Nonparametric_NormalVolatility(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
											double	C_k,
											ARM_result& result);

	
extern long ARMLOCAL_Nonparametric_NormalDistribution(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
											double C_S,double C_T,double C_k,
											ARM_result& result);

extern long ARMLOCAL_Nonparametric_NormalQuantile(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
											double C_S,double C_T,double C_k,
											ARM_result& result);


extern long ARMLOCAL_Nonparametric_LogNormalDistribution(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
											double C_S,double C_T,double C_k,
											ARM_result& result);

extern long ARMLOCAL_Nonparametric_LogNormalQuantile(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
											double C_S,double C_T,double C_k,
											ARM_result& result);


////////////////////////////////////////////
//// Function to create an CIR_ModelParamCreate
////////////////////////////////////////////
extern long ARMLOCAL_CIR_ModelParamsCreate(
	const vector<long>&		modelParamVec,
	ARM_result&				result, 
	long					objId);

////////////////////////////////////////////
//// Function to create an CIR_BondPrice
////////////////////////////////////////////

extern long ARMLOCAL_CIRBondPrice(
	const long&				modelParamsId,
	const double&			time,
	const double&			r0,
	const double&			phi,
	ARM_result&				result);

////////////////////////////////////////////
//// Function to create an CIR_BondDensity
////////////////////////////////////////////

extern long ARMLOCAL_CIRBondDensity(
	const long&				modelParamsId,
	const double&			time,
	const double&			r0,
	const double&			r,
	const double&			frequency,
	ARM_result&				result);


////////////////////////////////////////////
//// Function to create an CIR_BondDistribution
////////////////////////////////////////////

extern long ARMLOCAL_CIRBondDistribution(
	const long&				modelParamsId,
	const double&			time,
	const double&			r0,
	const double&			nbDiscr,
	const double&			frequency,
	ARM_result&				result);

////////////////////////////////////////////
//// Function to create an CIRBondWeight
////////////////////////////////////////////
extern long ARMLOCAL_CIRBondWeight(
	const long&				modelParamsId,
	ARM_result&				result );

////////////////////////////////////////////
//// Function to create an CIRBondExpectation
////////////////////////////////////////////
extern long ARMLOCAL_CIRBondExpectation(
	const long&				modelParamsId,
	ARM_result&				result );

////////////////////////////////////////////
//// Function to create an CIRBondVariance
////////////////////////////////////////////
extern long ARMLOCAL_CIRBondVariance(
	const long&				modelParamsId,
	ARM_result&				result );

////////////////////////////////////////////
//// Function to evaluate exponential riccati  equation
////////////////////////////////////////////

extern long ARMLOCAL_EXP_RICCATI_Create(
	const double&			alpha,
	const double&			beta,
	const double&			delta,
	const double&			lambda,
	const double&			x0,
	const double&			x1,
	const double&			x2,
	const double&			y0,
	const double&			t0,
	ARM_result&				result,
	long					objId=-1);

extern long ARMLOCAL_EXP_RICCATI_Price(
	const long&				RiccatiEq,
	const double&			t,
	ARM_result&				result );

extern long ARMLOCAL_EXP_RICCATI_Int(
	const long&				RiccatiEq,
	const double&			t,
	const double&			T,
	ARM_result&				result );

////////////////////////////////////////////
//// BiShifted Heston
////////////////////////////////////////////
extern long ARMLOCAL_BiShiftedHeston_VanillaOption(
										const double C_F1,
										const double C_V1,
										const double C_Vinfini1,
										const double C_lambda1,
										const double C_nu1,
										const double C_rho1,
										const double C_gamma1,
										const double C_F2,
										const double C_V2,
										const double C_Vinfini2,
										const double C_lambda2,
										const double C_nu2,
										const double C_rho2,
										const double C_gamma2,
										const vector<double>&  C_Correlations_vec,
										const double C_k,
										const double C_T,
										const double C_CallPut,
										const double C_LambdaB,
										const double C_Flag,
										ARM_result& result);



extern long ARMLOCAL_MertonOption(
	const double& Fwd,
	const double& Strike,
	const double& t,
	const double& CallPut,
	const double& sigma,
	const double& lambda1,
	const double& U1,
	const double& lambda2,
	const double& U2,
	const int& N,
	ARM_result& result);

extern long ARMLOCAL_BSImpliedVol(
	const double& F,
	const double& K,
	const double& T,
	const double& CallPut,
	const double& Target,
	ARM_result& result);

#endif


