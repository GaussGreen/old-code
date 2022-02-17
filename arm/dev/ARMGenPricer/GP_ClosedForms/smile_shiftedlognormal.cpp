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
#include <glob/firsttoinc.h>

#include <cmath>

#include <vector>
#include <iostream>
#include <iomanip>
#include <glob/expt.h>
#include "gpnumlib/gaussiananalytics.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/shifted_lognormal_formula.h"
#include "gpclosedforms/normal.h"
#include "gpcalib/densityfunctors.h"


/// gpnumlib
#include "gpnumlib/numfunction.h"
#include "gpnumlib/dichotomy.h"
#include "gpnumlib/newtonraphson.h"


#include "gpclosedforms/smile_shiftedlognormal.h"


CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_K_SHIFT_FOR_DERIVATION 0.0000001

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///    Utilitary functions
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double ShiftedLogNormal_Smile::call_option(double f,double K,double tex, double sigma, double alpha)
{
	return BS(f+alpha,K+alpha,tex,sigma);
}


double ShiftedLogNormal_Smile::digital_call_option(double f,double K,double tex,double sigma,  double alpha)
{
	double  d2,totalvol=sigma*sqrt(tex);

    d2 = (log((f+alpha)/(K+alpha)))/totalvol-0.5*totalvol ;

    return ARM_GaussianAnalytics::cdfNormal(d2);

}

/*
double ShiftedLogNormal_Smile::probability_distribution(const ArgumentList& Underlying, double x, double t)
	{
	 int argsize=Underlying.size();
	 if ((argsize<3)||(argsize>3))
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ShiftedLogNormal_Smile::probability_density  : bad argsize");
	 }
	 ArgumentList arg(
		 Underlying[0],	//INDEX1
		 x,
		 t,
		 Underlying[1],	//ALPHA1
		 Underlying[2],	//BETA1
		 K_CALL
		 );	
	 double s=ARM_CF_Shifted_Lognormal_Formula::specific_shift(1);
	 double digitale=-ARM_CF_Shifted_Lognormal_Formula::value(1,arg,s);
	 if(digitale>1.0) return 0.0;
	 return 1.0-digitale;
}
*/
double ShiftedLogNormal_Smile::probability_distribution(const ArgumentList& Underlying, double x, double t)
	{
	 int argsize=Underlying.size();
	 if ((argsize<3)||(argsize>3))
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ShiftedLogNormal_Smile::probability_density  : bad argsize");
	 }
	double f=Underlying[0];
	double v=Underlying[1];
	double m=Underlying[2];
	if(x==-m)
	{
		return 0;
	}
	else
	{
		double arg=(f+m)/(m+x);	
		if(arg<=0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ShiftedLogNormal_Smile::probability_density  : bad arg<!");
		}
		return 1.0-NormalCDF(log(arg)-t*v*v/2.);
	}
}


double ShiftedLogNormal_Smile::inverse_distribution(double f,double k,double t,double v, double m)
{
	if(k==0)
	{
		return -m;
	}
	else
	{
		double q=exp(t*v*v/2.+sqrt(t)*v*NormalCDFInverse(1.-k));
		return (f+m*(1.-q))/q;
	}
	
}


double ShiftedLogNormal_Smile::gaussian_to_distribution(double f,double x,double tex,double sigma,  double alpha)
{
	return (f+alpha)*exp(-0.5*sigma*sigma*tex+sigma*sqrt(tex)*x)-alpha;
}


/// -----------------------------------------------
/// given two call prices @ different strikes, 
///	compute vol & shift parameters
/// -----------------------------------------------
void ShiftedLogNormal_Smile::volshift_calibration(  double f,  double tex, 
													double K1, double price1,
													double K2, double price2,
													/// result
													double& sigma,
													double& alpha)
{
	/// -- NOT TESTED YET -- ///

	SLNDensity_Approx func(f, tex, K1, price1, K2, price2);

	UnaryFuncWithNumDerivative<double> funcWithDev(func);

	const double DEFAULT_PRECISION = 1.e-6;

	T_DichotomySolver< UnaryFuncWithNumDerivative<double> > solver(funcWithDev,0,DEFAULT_PRECISION,DEFAULT_PRECISION);
	solver.setInitialGuess(0.);/// solver.setInitialGuess(0., shiftinf,shiftsup);

	alpha = solver.Solve();
	sigma = func.GetVol(alpha);
}



////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a generic implementation that uses the following conventions :
/// 
///				For the shifted lognormal  Distributions  
///      Underlying[0],	//INDEX1
///		 Underlying[1],	//SIGMA1
///		 Underlying[2],	//ALPHA1
///		 
///
///
////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////

double ShiftedLogNormal_Smile::gaussian_to_distribution(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<3)||(argsize>3))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ShiftedLogNormal_Smile::gaussian_to_distribution  : bad argsize");
	}
	return gaussian_to_distribution(Underlying[0],x,t,Underlying[1],Underlying[2]);
}



double ShiftedLogNormal_Smile::distribution_to_gaussian(const ArgumentList& Underlying, double x, double t)
	{
	 int argsize=Underlying.size();
	 if ((argsize<7)||(argsize>7))
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_smile::gaussian_to_distribution  : bad argsize");
	 }
	 double f= Underlying[0];	//INDEX1
	 double sigma=Underlying[1]; //ALPHA1
	 double beta=Underlying[2];	//BETA1
	 double ratio=exp(sqrt(t)*sigma);
	 double liminf=f/pow(ratio,4);
	 double limsup=f*pow(ratio,4);
	 double logaveragingperiod=log(ratio);
	 double f1,f2;
	
	 if ((x>liminf)&&(x<limsup))
	 {
			 return NormalCDFInverse(probability_distribution(
				 Underlying,	//INDEX1
				 x,
				 t
				 ));
		 }
		 else if (x>=limsup)
	 {
		  f1=NormalCDFInverse(probability_distribution(
			 Underlying,	//INDEX1
			 limsup*exp(-logaveragingperiod),
			 t
			 ));
		  f2=NormalCDFInverse(probability_distribution(
			 Underlying,	//INDEX1
			 limsup,
			 t
			 ));
		 return (log(x)-log(limsup))*(f2-f1)/logaveragingperiod+f2;
	 }
	 else		///  if (x<=liminf)
	 {
		  f1=NormalCDFInverse(probability_distribution(
			 Underlying,	//INDEX1
			 liminf*exp(logaveragingperiod),
			 t
			 ));
		  f2=NormalCDFInverse(probability_distribution(
			 Underlying,	//INDEX1
			 liminf,
			 t
			 ));
		 return (log(x)-log(liminf))*(f1-f2)/logaveragingperiod+f2;
	 }

	}

double ShiftedLogNormal_Smile::quantile(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<3)||(argsize>3))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ShiftedLogNormal_Smile::gaussian_to_distribution  : bad argsize");
	}
	return  inverse_distribution(Underlying[0],x,t,Underlying[1],Underlying[2]);
}

double ShiftedLogNormal_Smile::probability_density(const ArgumentList& Underlying, double x, double t)
	{
	 int argsize=Underlying.size();
	 if ((argsize<3)||(argsize>3))
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ShiftedLogNormal_Smile::probability_density  : bad argsize");
	 }
	 ArgumentList arg(
		 Underlying[0],	//INDEX1
		 x,
		 t,
		 Underlying[1],	//ALPHA1
		 Underlying[2],	//BETA1
		 K_CALL
		 );	
	 double s=ARM_CF_Shifted_Lognormal_Formula::specific_shift(1);
	 return -ARM_CF_Shifted_Lognormal_Formula::value(1,1,arg,s,s);
}

/*
double ShiftedLogNormal_Smile::probability_distribution(const ArgumentList& Underlying, double x, double t)
	{
	 int argsize=Underlying.size();
	 if ((argsize<3)||(argsize>3))
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ShiftedLogNormal_Smile::probability_density  : bad argsize");
	 }
	 ArgumentList arg(
		 Underlying[0],	//INDEX1
		 x,
		 t,
		 Underlying[1],	//ALPHA1
		 Underlying[2],	//BETA1
		 K_CALL
		 );	
	 double s=ARM_CF_Shifted_Lognormal_Formula::specific_shift(1);
	 double digitale=-ARM_CF_Shifted_Lognormal_Formula::value(1,arg,s);
	 if(digitale>1.0) return 0.0;
	 return 1.0-digitale;
}
*/

double ShiftedLogNormal_Smile::probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x, double t)
	{
	 int argsize=Underlying.size();
	 if ((argsize<3)||(argsize>3))
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ShiftedLogNormal_Smile::probability_density  : bad argsize");
	 }
	 ArgumentList arg(
		 Underlying[0],	//INDEX1
		 x,
		 t,
		 Underlying[1],	//ALPHA1
		 Underlying[2],	//BETA1
		 K_CALL
		 );	
	 double s=ARM_CF_Shifted_Lognormal_Formula::specific_shift(1);
	 double s2=ARM_CF_Shifted_Lognormal_Formula::specific_shift(i);
	 double digitale=-ARM_CF_Shifted_Lognormal_Formula::value(1,i,arg,s,s2);
	 return -digitale;
}

ArgumentList*  ShiftedLogNormal_Smile:: HomotheticTransformation(const ArgumentList* Underlying, double positivenumber)
{
	return new ArgumentList(
		positivenumber*(*Underlying)[0],
		(*Underlying)[1],
		positivenumber*(*Underlying)[2]
		);
}


CC_END_NAMESPACE()
 

#undef  ARM_CF_K_SHIFT_FOR_DERIVATION
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
