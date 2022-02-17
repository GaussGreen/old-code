/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_shiftedlognormal.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"

#include <cmath>

#include <vector>
#include <iostream>
#include <iomanip>
#include "expt.h"
#include "gpnumlib/gaussiananalytics.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/spreadoption_lognormal.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/shifted_lognormal_formula.h"
//#include "gpcalib/densityfunctors.h"

/// gpnumlib
#include "gpnumlib/numfunction.h"
#include "gpnumlib/dichotomy.h"
#include "gpnumlib/newtonraphson.h"


#include "gpclosedforms/smile_2logsmiled.h"


CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_K_SHIFT_FOR_DERIVATION 0.0000001

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///    Utilitary functions
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double Spread_Shifted2LogNormal_Smile::call_option(double f1,double f2,double K,double tex,double sigma1, 
		double sigma2, double alpha,double rho, int n)
{
	return f1*CorrelationRobust_SpreadDigitalCall( f1*exp(sigma1*sigma1*tex/2.), f2*exp(rho*sigma1*sigma2*tex), sigma1, sigma2, rho, K+alpha, tex, n)-
		f2*CorrelationRobust_SpreadDigitalCall( f1*exp(rho*sigma1*sigma2*tex), f2*exp(sigma2*sigma2*tex/2.), sigma1, sigma2, rho, K+alpha, tex, n)-
		K*CorrelationRobust_SpreadDigitalCall( f1, f2, sigma1, sigma2, rho, K+alpha, tex, n);

}


double Spread_Shifted2LogNormal_Smile::digital_call_option(double f1,double f2,double K,double tex,double sigma1, 
		double sigma2, double alpha,double rho, int n)
{
 ///   suppose that K+alpha>0
	return  SpreadDigitalCall( f1, f2, sigma1, sigma2, rho, K+alpha, tex, n);
}


double Spread_Shifted2LogNormal_Smile::probability_distribution(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<7)||(argsize>7))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Spread_Shifted2LogNormal_Smile::probability_density  : bad argsize");
	}
	double f1=Underlying[0];
	double sigma1=Underlying[1];
	double f2=Underlying[2];
	double sigma2=Underlying[3];
	double alpha =Underlying[4];
	double rho=Underlying[5];
	int n=Underlying[6];
	
	if(x+alpha>0) return 1.-digital_call_option( f1, f2, x, t, sigma1, sigma2,  alpha, rho,  n);
	else return digital_call_option( f2, f1, -x, t, sigma2, sigma1,  -alpha, rho,  n);
}


double Spread_Shifted2LogNormal_Smile::inverse_distribution(double f1,double f2,double proba,double tex,double sigma1, 
		double sigma2,double alpha, double rho, int n)
{
	class DistributionToInverse : public DoubleToDoubleFunc 
	{
	public: 
		double f01;
		double sigma01;
		double f02;
		double sigma02;
		double alpha0;
		double rho0;
		double tex0;
		int n0;
		mutable ArgumentList arg;
		DistributionToInverse(double f1a,double sigma1a,double f2a,double sigma2a,double alphaa,int rhoa,
			int na, double texa):
		f01(f1a),tex0(texa),alpha0(alphaa),sigma01(sigma1a),sigma02(sigma2a),rho0(rhoa),n0(na),
			arg(f1a,sigma1a,f2a,sigma2a,alphaa,rhoa,na) {}
		
		virtual double operator() (double K0)  const
		{
			return Spread_Shifted2LogNormal_Smile::probability_distribution( arg,  K0,  tex0);
		}
	};
	
	DistributionToInverse x(f1,sigma1,f2,sigma2,alpha,rho,n,tex);

	/// initial guess
	double f=sqrt(f1*f1*sigma1*sigma1+f2*f2*sigma2*sigma2-2.*rho*f1*f2*sigma1*sigma2);
	double guess = alpha;
	{
		return Inverse(x,Inverse::REAL)(proba,guess,f/5.,1e-12);  // suppose that  Y is always postive
	}	
}


double Spread_Shifted2LogNormal_Smile::gaussian_to_distribution(double f1,double f2,double K,double tex,double sigma1,
														double sigma2,double alpha,double rho, int n)
{
	return inverse_distribution(f1, f2, ARM_GaussianAnalytics::cdfNormal(K), tex, sigma1, sigma2, alpha,  rho,  n);
}





////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a generic implementation that uses the following conventions :
/// 
///				For the Shifted Lognormal Distributions the  
///		 Underlying[0],	//INDEX1
///      Underlying[1],	//SIGMA1
///		 Underlying[2],	//INDEX2
///      Underlying[3],	//SIGMA2
///		 Underlying[4], //ALPHA
///		 Underlying[5],	//RHO
///		 Underlying[6],	//N
///		 
///
///
////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////

double Spread_Shifted2LogNormal_Smile::gaussian_to_distribution(const ArgumentList& Underlying, double x, double t)
{

	int argsize=Underlying.size();
	if ((argsize<7)||(argsize>7))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Spread_Shifted2LogNormal_Smile::gaussian_to_distribution  : bad argsize");
	}
	return gaussian_to_distribution(Underlying[0],Underlying[2],x,t,Underlying[1],Underlying[3],Underlying[4],Underlying[5],Underlying[6]);
}



double Spread_Shifted2LogNormal_Smile::distribution_to_gaussian(const ArgumentList& Underlying, double x, double t)
	{
	 int argsize=Underlying.size();
	 if ((argsize<7)||(argsize>7))
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Spread_Shifted2LogNormal_Smile::gaussian_to_distribution  : bad argsize");
	 }
	 return NormalCDFInverse(probability_distribution(
				 Underlying,	//INDEX1
				 x,
				 t
				 ));
	}

double Spread_Shifted2LogNormal_Smile::quantile(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<7)||(argsize>7))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Spread_Shifted2LogNormal_Smile::gaussian_to_distribution  : bad argsize");
	}
	double f1=Underlying[0];
	double sigma1=Underlying[1];
	double f2=Underlying[2];
	double sigma2=Underlying[3];
	double alpha =Underlying[4];
	double rho=Underlying[5];
	int n=Underlying[6];
	return  inverse_distribution( f1, f2, x, t, sigma1, sigma2, alpha,  rho,  n);
}

double Spread_Shifted2LogNormal_Smile::probability_density(const ArgumentList& Underlying, double x, double t)
	{
	 int argsize=Underlying.size();
	 if ((argsize<7)||(argsize>7))
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Spread_Shifted2LogNormal_Smile::probability_density  : bad argsize");
	 }
	double f1=Underlying[0];
	double sigma1=Underlying[1];
	double f2=Underlying[2];
	double sigma2=Underlying[3];
	double alpha =Underlying[4];
	double rho=Underlying[5];
	int n=Underlying[6];

	double s=ARM_CF_Shifted_Lognormal_Formula::specific_shift(1);
	double v1=digital_call_option( f1, f2, x+s/2., t, sigma1, sigma2,  alpha, rho,  n);
	double v2=digital_call_option( f1, f2, x-s/2., t, sigma1, sigma2,  alpha, rho,  n);
	return (v1-v2)/s;
}


double Spread_Shifted2LogNormal_Smile::probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x, double t)
{
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Spread_Shifted2LogNormal_Smile::probability_density  : Should not reach here !");
	}
	
	return 0;
}

ArgumentList*  Spread_Shifted2LogNormal_Smile:: HomotheticTransformation(const ArgumentList* Underlying, double positivenumber)
{
	 int argsize=Underlying->size();
	 if ((argsize<7)||(argsize>7))
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Spread_Shifted2LogNormal_Smile::HomotheticTransformation  : bad argsize");
	 }
	return new ArgumentList(
		positivenumber*(*Underlying)[0],
		(*Underlying)[1],
		positivenumber*(*Underlying)[2],
		(*Underlying)[3],
		positivenumber*(*Underlying)[4],
		(*Underlying)[5],
		(*Underlying)[6]
		);
}


CC_END_NAMESPACE()
 

#undef  ARM_CF_K_SHIFT_FOR_DERIVATION
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
