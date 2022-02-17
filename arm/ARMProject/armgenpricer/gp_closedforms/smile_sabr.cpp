/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_sabr.cpp
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
#include "gpbase/numericconstant.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/sabrimpliedvol.h"


CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_K_SHIFT_FOR_DERIVATION 0.0000001




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///    Class : SABR_smile
///
///    Approach following the work of Jean David Aube on Direct equivalence of the SABR approach with the LogNormal vol
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double SABR_smile::call_option(double f,double K,double tex, double alpha, double beta, double rho, double nu,int flag,int nbsteps,
							   double alpha_exp,double alpha_tanh,double kb_tanh)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,K_CALL,flag,nbsteps,alpha_exp,alpha_tanh,kb_tanh);
	Power_Expression<ARM_CF_SABR_VanillaOption_Formula> y;
	return y(a);
}


double SABR_smile::digital_call_option(double f,double K,double tex, double alpha, double beta, double rho, double nu,int flag,int nbsteps,
							   double alpha_exp,double alpha_tanh,double kb_tanh)
{
	if(K<=f*1e-9)
	{
		return 1.;
	}
	double shift=(f<K) ? f/10000.:K/10000.;
	double result= (-SABR_smile::call_option(f,K+shift,tex,alpha,beta,rho,nu,flag,nbsteps,alpha_exp,alpha_tanh,kb_tanh)+
		SABR_smile::call_option(f,K-shift,tex,alpha,beta,rho,nu,flag,nbsteps,alpha_exp,alpha_tanh,kb_tanh))/(2.*shift);
	if (result<=0) return 1e-12;
	else return result;
}



inline double cmin(double x, double y) {return ((x<=y) ? x : y);}

double SABR_smile::inverse_distribution(double f,double proba,double tex, double alpha, double beta, double rho, double nu,int flag,int nbsteps,
							   double alpha_exp,double alpha_tanh,double kb_tanh)
{
	class DistributionToInverse : public DoubleToDoubleFunc 
	{
	public: 
		double f0;
		double tex0;
		double alpha0;
		double beta0;
		double rho0;
		double nu0;
		int flag0;
		int nbsteps0;
		double alpha_exp0;
		double alpha_tanh0;
		double kb_tanh0;
		mutable ArgumentList arg;
		DistributionToInverse(double f1,double tex1, double alpha1, double beta1, double rho1, double nu1,int flag1,int nbsteps1,
							   double alpha_exp1,double alpha_tanh1,double kb_tanh1):
		f0(f1),tex0(tex1),alpha0(alpha1),beta0(beta1),rho0(rho1),nu0(nu1),flag0(flag1),nbsteps0(nbsteps1),alpha_exp0(alpha_exp1),
			alpha_tanh0(alpha_tanh1),kb_tanh0(kb_tanh1),
			arg(f1,0.,tex1,alpha1,beta1,rho1,nu1,K_CALL,flag1,nbsteps1,alpha_exp1,alpha_tanh1,kb_tanh0) {}
		
		virtual double operator() (double K0)  const
		{
			if(K0<=f0*1e-9)
			{
				return 0.;
			}
			arg.set_nth(1,K0);
			// ArgumentList a(f0,K0,tex0,alpha0,beta0,rho0,nu0,K_CALL,flag0,nbsteps0);
			double s=ARM_CF_SABR_VanillaOption_Formula::specific_shift(1);
			return 1.+ARM_CF_SABR_VanillaOption_Formula::value(1,arg,s);
		}
	};
	
	DistributionToInverse x(f,tex,alpha,beta,rho,nu,flag,nbsteps,alpha_exp,alpha_tanh,kb_tanh);

	const double EPS=0.0001;
	/// initial guess
	double vol = alpha * pow(f, beta-1.0);/// approx of atm lognormal vol (à la grosse louche)
	double guess = f/exp(0.5*tex*vol*vol+sqrt(tex)*vol*NormalCDFInverse(1.-proba));
	if(proba>EPS)
	{
		return Inverse(x,Inverse::ALWAYSPOSITIVE)(proba,guess,f/5.,1e-12);  // suppose that  Y is always postive
	}
	else
	{	
		double p0=Inverse(x,Inverse::ALWAYSPOSITIVE)(EPS,guess,f/5.,1e-12);  
		return p0*proba/EPS;

	}
}


double SABR_smile::gaussian_to_distribution(double f,double x,double tex, double alpha, double beta, double rho, double nu,int flag,int nbsteps,
											double alpha_exp,double alpha_tanh,double kb_tanh)
{
	return inverse_distribution(f,ARM_GaussianAnalytics::cdfNormal(x),tex,alpha,beta,rho,nu,flag,nbsteps,alpha_exp,alpha_tanh,kb_tanh);
}


////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a generic implementation that uses the following conventions :
/// 
///				For the SABR Distributions the  
///      Underlying[0],	//INDEX1
///		 Underlying[1],	//ALPHA1
///		 Underlying[2],	//BETA1
///		 Underlying[3],	//RHO1
///		 Underlying[4],	//NU1
///		 Underlying[5]	//FLAG
///		 Underlying[6]	//NBSTEPS
///		 Underlying[7]	//ALPHA_EXP
///		 Underlying[8]	//ALPHA_TANH
///		 Underlying[9]	//KB_TANH
///
///
////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////

double SABR_smile::gaussian_to_distribution(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<10)||(argsize>10))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_smile::gaussian_to_distribution  : bad argsize");
	}
	double liminf=-3.5;
	double limsup=5;
	double averagingperiod=1;
	if ((x>liminf)&&(x<limsup))
	{
		return inverse_distribution(
			Underlying[0],	//INDEX1
			ARM_GaussianAnalytics::cdfNormal(x),
			t,
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5], //FLAG
			Underlying[6],	//NBSTEPS
			Underlying[7],	//ALPHA_EXP
			Underlying[8],	//ALPHA_TANH
			Underlying[9]	//KB_TANH
			);
	}
	else if (x>=limsup)
	{
		double f1=log(inverse_distribution(
			Underlying[0],	//INDEX1
			ARM_GaussianAnalytics::cdfNormal(limsup-averagingperiod),
			t,
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5], //FLAG
			Underlying[6],	//NBSTEPS
			Underlying[7],	//ALPHA_EXP
			Underlying[8],	//ALPHA_TANH
			Underlying[9]	//KB_TANH
			));
		double f2=log(inverse_distribution(
			Underlying[0],	//INDEX1
			ARM_GaussianAnalytics::cdfNormal(limsup),
			t,
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5], //FLAG
			Underlying[6],	//NBSTEPS
			Underlying[7],	//ALPHA_EXP
			Underlying[8],	//ALPHA_TANH
			Underlying[9]	//KB_TANH
			));
		return exp((x-limsup)*(f2-f1)/averagingperiod+f2);
	}
	else // if (x<=liminf)
	{
		double f1=log(inverse_distribution(
			Underlying[0],	//INDEX1
			ARM_GaussianAnalytics::cdfNormal(liminf),
			t,
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5], //FLAG
			Underlying[6],	//NBSTEPS
			Underlying[7],	//ALPHA_EXP
			Underlying[8],	//ALPHA_TANH
			Underlying[9]	//KB_TANH
			));
		double f2=log(inverse_distribution(
			Underlying[0],	//INDEX1
			ARM_GaussianAnalytics::cdfNormal(liminf+averagingperiod),
			t,
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5], //FLAG;
			Underlying[6],	//NBSTEPS
			Underlying[7],	//ALPHA_EXP
			Underlying[8],	//ALPHA_TANH
			Underlying[9]	//KB_TANH
			));
		return exp((x-liminf)*(f2-f1)/averagingperiod+f1);
	}
	
}
double SABR_smile::distribution_to_gaussian(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<10)||(argsize>10))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_smile::gaussian_to_distribution  : bad argsize");
	}
	double f= Underlying[0];	//INDEX1
	double alpha=Underlying[1]; //ALPHA1
	double beta=Underlying[2];	//BETA1
	double rho=Underlying[3];	//RHO1
	double nu=Underlying[4];	//NU1
	double sigma=alpha*pow(f,beta-1.);
	double ratio=exp(sqrt(t)*sigma);
	double liminf=f/pow(ratio,4);
	double limsup=f*pow(ratio,4);
	double logaveragingperiod=log(ratio);
	double f1,f2;
	
	if ((x>liminf)&&(x<limsup))
	{
		return ARM_GaussianAnalytics::cdfNormal_Inv(probability_distribution(
			Underlying,	//INDEX1
			x,
			t
			));
	}
	else if (x>=limsup)
	{
		f1=ARM_GaussianAnalytics::cdfNormal_Inv(probability_distribution(
			Underlying,	//INDEX1
			limsup*exp(-logaveragingperiod),
			t
			));
		f2=ARM_GaussianAnalytics::cdfNormal_Inv(probability_distribution(
			Underlying,	//INDEX1
			limsup,
			t
			));
		return (log(x)-log(limsup))*(f2-f1)/logaveragingperiod+f2;
	}
	else /// last case (x<=liminf)
	{
		f1=ARM_GaussianAnalytics::cdfNormal_Inv(probability_distribution(
			Underlying,	//INDEX1
			liminf*exp(logaveragingperiod),
			t
			));
		f2=ARM_GaussianAnalytics::cdfNormal_Inv(probability_distribution(
			Underlying,	//INDEX1
			liminf,
			t
			));
		return (log(x)-log(liminf))*(f1-f2)/logaveragingperiod+f2;
	}
}

double SABR_smile::distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying, double x, double t)
{
	double s=ARM_CF_SABR_VanillaOption_Formula::specific_shift(i);
	ArgumentList* Underlying1;
	ArgumentList* Underlying2;
	switch(i)
	{	
	case ARM_CF_SABR_VanillaOption_Formula::FORWARD :										/// FORWARD
		{
			Underlying1=new ArgumentList(Underlying[0]-s,Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],
				Underlying[7],Underlying[8],Underlying[9]);
			Underlying2=new ArgumentList(Underlying[0]+s,Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],
				Underlying[7],Underlying[8],Underlying[9]);
			break;
		}
	case ARM_CF_SABR_VanillaOption_Formula::ALPHA :											/// ALPHA 
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1]-s,Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],
				Underlying[7],Underlying[8],Underlying[9]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1]+s,Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],
				Underlying[7],Underlying[8],Underlying[9]);
			break;
		}
	case ARM_CF_SABR_VanillaOption_Formula::BETA :											/// BETA
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2]-s,Underlying[3],Underlying[4],Underlying[5],Underlying[6],
				Underlying[7],Underlying[8],Underlying[9]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2]+s,Underlying[3],Underlying[4],Underlying[5],Underlying[6],
				Underlying[7],Underlying[8],Underlying[9]);
			break;
		}
	case ARM_CF_SABR_VanillaOption_Formula::RHO :											/// RHO
		{		
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3]-s,Underlying[4],Underlying[5],Underlying[6],
				Underlying[7],Underlying[8],Underlying[9]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3]+s,Underlying[4],Underlying[5],Underlying[6],
				Underlying[7],Underlying[8],Underlying[9]);
			break;
		}
	case ARM_CF_SABR_VanillaOption_Formula::NU :											/// NU
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4]-s,Underlying[5],Underlying[6],
				Underlying[7],Underlying[8],Underlying[9]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4]+s,Underlying[5],Underlying[6],
				Underlying[7],Underlying[8],Underlying[9]);
			break;
		}
		
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_smile::distribution_to_gaussian_first_derivative : incorrect input");
	};
	double v1=distribution_to_gaussian( *Underlying1,  x,  t);
	double v2=distribution_to_gaussian( *Underlying2,  x,  t);
	delete Underlying1;
	delete Underlying2;
	return (v2-v1)/(2.*s);
				
}



double SABR_smile::quantile(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<10)||(argsize>10))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_smile::quantile  : bad argsize");
	}
	return inverse_distribution(
		Underlying[0],	//INDEX1
		x,
		t,
		Underlying[1],	//ALPHA1
		Underlying[2],	//BETA1
		Underlying[3],	//RHO1
		Underlying[4],	//NU1
		Underlying[5], //FLAG
		Underlying[6],	//NBSTEPS
		Underlying[7],	//ALPHA_EXP
		Underlying[8],	//ALPHA_TANH
		Underlying[9]	//KB_TANH
		);		
}

double SABR_smile::probability_density(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<10)||(argsize>10))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_smile::quantile  : bad argsize");
	}
	ArgumentList arg(
		Underlying[0],	//INDEX1
		x,
		t,
		Underlying[1],	//ALPHA1
		Underlying[2],	//BETA1
		Underlying[3],	//RHO1
		Underlying[4],	//NU1
		K_CALL,
		Underlying[5], //FLAG
		Underlying[6],	//NBSTEPS
		Underlying[7],	//ALPHA_EXP
		Underlying[8],	//ALPHA_TANH
		Underlying[9]	//KB_TANH
		);
	if(x<=0) return 0.0;
	double s=ARM_CF_SABR_VanillaOption_Formula::specific_shift(1);
	return fabs(ARM_CF_SABR_VanillaOption_Formula::value(1,1,arg,s,s));
}

double SABR_smile::probability_distribution(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<10)||(argsize>10))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_smile::quantile  : bad argsize");
	}
	ArgumentList arg(
		Underlying[0],	//INDEX1
		x,
		t,
		Underlying[1],	//ALPHA1
		Underlying[2],	//BETA1
		Underlying[3],	//RHO1
		Underlying[4],	//NU1
		K_CALL,
		Underlying[5], //FLAG
		Underlying[6],	//NBSTEPS
		Underlying[7],	//ALPHA_EXP
		Underlying[8],	//ALPHA_TANH
		Underlying[9]	//KB_TANH
		);	
	if(x<=0) return 0.0;
	double s=ARM_CF_SABR_VanillaOption_Formula::specific_shift(1);
	double digitale=-ARM_CF_SABR_VanillaOption_Formula::value(1,arg,s);
	if(digitale>1.0) return 0.0;
	return 1.0-digitale;
}

double SABR_smile::probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x, double t)
{
	double v= distribution_to_gaussian( Underlying,x,t);
	double derive=distribution_to_gaussian_first_derivative(i,Underlying,x,t);
	return exp(-v*v/2.)*derive*ARM_NumericConstants::ARM_INVSQRT2PI;
}

ArgumentList*  SABR_smile:: HomotheticTransformation(const ArgumentList* Underlying, double positivenumber)
{
	return new ArgumentList(
		positivenumber*(*Underlying)[0],
		pow(positivenumber,1.-(*Underlying)[2])*(*Underlying)[1],
		(*Underlying)[2],
		(*Underlying)[3],
		(*Underlying)[4],
		(*Underlying)[5],
		(*Underlying)[6],
		positivenumber*(*Underlying)[7],
		(*Underlying)[8],
		positivenumber*(*Underlying)[9]
		);
}



CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/