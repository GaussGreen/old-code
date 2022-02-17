/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2007
 *
 *  basic functions for the closed form framework 
 *
 *	\file distribution_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2007
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/numericconstant.h"

#include "gpclosedforms/bivariate_normal.h"
#include "gpclosedforms/gamma.h"
#include "gpclosedforms/incompletebeta.h"
#include "gpclosedforms/hypergeometric.h"
#include "gpclosedforms/lambert_function.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/sabr_calibration.h"
#include "gpclosedforms/heston_calibration.h"
#include "gpclosedforms/smile_shiftedlognormal.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/smile_bisabr.h"
#include "gpclosedforms/student_copula.h"
#include "gpclosedforms/whittaker.h"
#include "gpclosedforms/bessel.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/tridiagonalsolve.h"
#include "gpclosedforms/eigenvalues.h"
#include "gpclosedforms/smile_2logsmiled.h"
#include "gpbase/env.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/typedef.h"
#include "gpclosedforms/sabr_calibration.h"
#include "gpclosedforms/heston_calibration.h"
#include "gpclosedforms/nonparametric_quantile.h"
#include "gpclosedforms/nonparametric_spline.h"


using namespace std;

CC_BEGIN_NAMESPACE(ARM)



double Export_ShiftedLogNormal_Quantile(double f,double K,double tex,double sigma, double alpha)
{
	return ShiftedLogNormal_Smile::inverse_distribution(f,K,tex,sigma,alpha);
}


double Export_ShiftedLogNormal_Distribution(double f,double x,double t,double v, double m)
{
	ArgumentList a(f,v,m);
	return ShiftedLogNormal_Smile::probability_distribution(a,x,t);
}

double Export_SABR_Quantile(double f,double k,double T,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps)

{
	return SABR_smile::inverse_distribution( f, k, T,alpha,beta,rho,nu,Sabr_Type,nbsteps,f/4.,1.5,f/2.);
}

double Export_SABR_Quantile(double f,double k,double T,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps,
							double alpha_exp,double alpha_tanh,double kb_tanh)

{
	return SABR_smile::inverse_distribution( f, k, T,alpha,beta,rho,nu,Sabr_Type,nbsteps,alpha_exp,alpha_tanh,kb_tanh);
}

double Export_SABR_Distribution(double f,double x,double t,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps)

{
	ArgumentList a(f,alpha,beta,rho,nu,Sabr_Type,nbsteps,f/4.,1.5,f/2.);
	return SABR_smile::probability_distribution( a,x,t);
}

double Export_SABR_Distribution(double f,double x,double t,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps,
								double alpha_exp,double alpha_tanh,double kb_tanh)

{
	ArgumentList a(f,alpha,beta,rho,nu,Sabr_Type,nbsteps,alpha_exp,alpha_tanh,kb_tanh);
	return SABR_smile::probability_distribution( a,x,t);
}

double Export_SABR_Density(double f,double x,double t,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps)

{
	ArgumentList a(f,alpha,beta,rho,nu,Sabr_Type,nbsteps,f/4.,1.5,f/2.);
	return SABR_smile::probability_density( a,x,t);
}

double Export_SABR_Density(double f,double x,double t,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps,
						   double alpha_exp,double alpha_tanh,double kb_tanh)

{
	ArgumentList a(f,alpha,beta,rho,nu,Sabr_Type,nbsteps,alpha_exp,alpha_tanh,kb_tanh);
	return SABR_smile::probability_density( a,x,t);
}

double Export_BiSABR_Quantile(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag,
		double alpha_exp,double alpha_tanh,double kb_tanh)

{
	return BiSABR_smile::inverse_distribution(f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,x,T,flag,
		alpha_exp,alpha_tanh,kb_tanh);
}

double Export_BiSABR_Quantile(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag)

{
	return BiSABR_smile::inverse_distribution(f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,x,T,flag,(f1+f2)/8.,1.5,(f1+f2)/4.);
}

double Export_BiSABR_Distribution(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag,
		double alpha_exp,double alpha_tanh,double kb_tanh)

{
		ArgumentList a(f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,flag,alpha_exp,alpha_tanh,kb_tanh);
	return BiSABR_smile::probability_distribution( a,x,T);
}

double Export_BiSABR_Distribution(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag)

{
		ArgumentList a(f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,flag,(f1+f2)/8.,1.5,(f1+f2)/4.);
	return BiSABR_smile::probability_distribution( a,x,T);
}


double Export_BiSABR_Density(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag,
		double alpha_exp,double alpha_tanh,double kb_tanh)

{
	ArgumentList a(f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,flag,alpha_exp,alpha_tanh,kb_tanh);
	return BiSABR_smile::probability_density( a,x,T);
}

double Export_BiSABR_Density(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag)

{
	ArgumentList a(f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,flag,(f1+f2)/8.,1.5,(f1+f2)/2.);
	return BiSABR_smile::probability_density( a,x,T);
}

double Export_Shifted2LogNormal_Quantile(double f1,double sigma1,double f2,double sigma2,double alpha,double rho,double x,double T,int n)

{
	return Spread_Shifted2LogNormal_Smile::inverse_distribution( f1, f2, x, T, sigma1, sigma2, alpha,  rho,  n);
}

double Export_Shifted2LogNormal_Distribution(double f1,double sigma1,double f2,double sigma2,double alpha,double rho,double x, double T,int n)

{
		ArgumentList a(f1,sigma1,f2,sigma2,alpha,rho,n);
	return Spread_Shifted2LogNormal_Smile::probability_distribution( a,x,T);
}

double Export_Shifted2LogNormal_Density(double f1,double sigma1,double f2,double sigma2,double alpha,double rho,double x,double T,int n)

{
	ArgumentList a(f1,sigma1,f2,sigma2,alpha,rho,n);
	return Spread_Shifted2LogNormal_Smile::probability_density( a,x,T);
}


double Export_Student_Quantile(double rank, double x)

{
	ArgumentList a(0.5,rank);
	return StudentCopula::marginal_quantile(a, x,0.);
}

double Export_Student_Distribution(double rank, double x)

{
	ArgumentList a(0.5,rank);
	return StudentCopula::marginal_distribution(a, x,0.);
}


double Export_NonParametric_LN_Distribution(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag, double S,double T, double strike)
{
	int nb=x->size();
	ARM_GP_Vector y2(nb);
	double beginderivative, endderivative;
	NonParametric_LogVolatility_TailSlope_Compute(x,y, index_begin, index_end,
							 beginflag, endflag, &beginderivative,&endderivative );
	cubicspline_precompute(x,y, beginderivative, endderivative, y2);
	return NonParametric_LN_Distribution_Interpolate( x,y,&y2,index_begin, index_end,beginflag, endflag,  S, T, strike);
}


double Export_NonParametric_LN_Quantile(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag,double S,double T, double proba)
{
	int nb=x->size();
	ARM_GP_Vector y2(nb);
	double beginderivative, endderivative;
	NonParametric_LogVolatility_TailSlope_Compute(x,y, index_begin, index_end,
							 beginflag, endflag, &beginderivative,&endderivative );
	cubicspline_precompute(x,y, beginderivative, endderivative, y2);
	return NonParametric_LN_Quantile_Interpolate(x,y,&y2,index_begin, index_end,beginflag, endflag, S, T,  proba);
}

double Export_NonParametric_LogVolatility(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag, double strike)
{
	int nb=x->size();
	ARM_GP_Vector y2(nb);
	double beginderivative, endderivative;
	NonParametric_LogVolatility_TailSlope_Compute(x,y, index_begin, index_end,
							 beginflag, endflag, &beginderivative,&endderivative );
	cubicspline_precompute(x,y, beginderivative, endderivative, y2);
	return NonParametric_LogVolatility_Interpolate(x,y,&y2,index_begin, index_end,beginflag, endflag,strike);
}

double Export_NonParametric_N_Distribution(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag, double S,double T, double strike)
{
	int nb=x->size();
	ARM_GP_Vector y2(nb);
	double beginderivative, endderivative;
	NonParametric_NormalVolatility_TailSlope_Compute(x,y, index_begin, index_end,
							 beginflag, endflag, &beginderivative,&endderivative );
	cubicspline_precompute(x,y, beginderivative, endderivative, y2);
	return NonParametric_N_Distribution_Interpolate( x,y,&y2,index_begin, index_end,beginflag, endflag,  S, T, strike);
}


double Export_NonParametric_N_Quantile(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag,double S,double T, double proba)
{
	int nb=x->size();
	ARM_GP_Vector y2(nb);
	double beginderivative, endderivative;
	NonParametric_NormalVolatility_TailSlope_Compute(x,y, index_begin, index_end,
							 beginflag, endflag, &beginderivative,&endderivative );
	cubicspline_precompute(x,y, beginderivative, endderivative, y2);
	return NonParametric_N_Quantile_Interpolate(x,y,&y2,index_begin, index_end,beginflag, endflag, S, T,  proba);
}


double Export_NonParametric_NormalVolatility(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag, double strike)
{
	int nb=x->size();
	ARM_GP_Vector y2(nb);
	double beginderivative, endderivative;
	NonParametric_NormalVolatility_TailSlope_Compute(x,y, index_begin, index_end,
							 beginflag, endflag, &beginderivative,&endderivative );
	cubicspline_precompute(x,y, beginderivative, endderivative, y2);
	return NonParametric_NormalVolatility_Interpolate(x,y,&y2,index_begin, index_end,beginflag, endflag,strike);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
