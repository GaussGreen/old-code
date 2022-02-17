/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_bisabr.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date Octobre 2006
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
#include "gpclosedforms/bisabr_spreadoption_formula.h"
#include "gpclosedforms/smile_bisabr.h"


CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_K_SHIFT_FOR_DERIVATION 0.0000001


inline double cmin(double x, double y) {return ((x<=y) ? x : y);}
inline double cmax(double x, double y) {return ((x<=y) ? y : x);}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///    Class : BiSABR_smile
///
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double BiSABR_smile::call_option(double f1,double alpha1,double beta1, double rho1, double nu1,
								 double f2,double alpha2,double beta2, double rho2, double nu2,
								 double rhos,double rhov,double rhoc12, double rhoc21, double K,double T,int flag,
								 double alpha_exp,double alpha_tanh,double kb_tanh)
{
	ArgumentList a(f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,K,T,K_CALL,flag,alpha_exp,alpha_tanh,kb_tanh);
	Power_Expression<ARM_CF_BiSABR_SpreadOption_Formula> y;
	return y(a);
}


double BiSABR_smile::digital_call_option(double f1,double alpha1,double beta1, double rho1, double nu1,
								 double f2,double alpha2,double beta2, double rho2, double nu2,
								 double rhos,double rhov,double rhoc12, double rhoc21, double K,double T,int flag,
								 double alpha_exp,double alpha_tanh,double kb_tanh)
{
	double shift=cmax(f1/10000.,f2/10000.);
	double result= (-BiSABR_smile::call_option(f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,K+shift,T,flag,alpha_exp,alpha_tanh,kb_tanh)+
		BiSABR_smile::call_option(f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,K-shift,T,flag,alpha_exp,alpha_tanh,kb_tanh))/(2.*shift);
	if (result<=0) return 1e-12;
	else return result;
}


double BiSABR_smile::inverse_distribution(double f1,double alpha1,double beta1, double rho1, double nu1,
								 double f2,double alpha2,double beta2, double rho2, double nu2,
								 double rhos,double rhov,double rhoc12, double rhoc21, double proba,double T,int flag,
								 double alpha_exp,double alpha_tanh,double kb_tanh)
{
	class DistributionToInverse : public DoubleToDoubleFunc 
	{
	public: 
		double f10;
		double alpha10;
		double beta10;
		double rho10;
		double nu10;
		double f20;
		double alpha20;
		double beta20;
		double rho20;
		double nu20;
		double rhos0;
		double rhov0;
		double rhoc120;
		double rhoc210;
		double T0;
		int flag0;
		double alpha_exp0;
		double alpha_tanh0;
		double kb_tanh0;

		mutable ArgumentList arg;
		DistributionToInverse(double f1a,double alpha1a,double beta1a, double rho1a, double nu1a,
								 double f2a,double alpha2a,double beta2a, double rho2a, double nu2a,
								 double rhosa,double rhova,double rhoc12a, double rhoc21a, double Ta,int flaga,
								 double alpha_expa,double alpha_tanha,double kb_tanha):
			f10(f1a),alpha10(alpha1a),beta10(beta1a),rho10(rho1a),nu10(nu1a),
			f20(f2a),alpha20(alpha1a),beta20(beta1a),rho20(rho1a),nu20(nu1a),
			T0(Ta),flag0(flaga),alpha_exp0(alpha_expa),alpha_tanh0(alpha_tanha),kb_tanh0(kb_tanha),
			arg(f1a,alpha1a,beta1a,rho1a,nu1a,f2a,alpha2a,beta2a,rho2a,nu2a,rhosa,rhova,rhoc12a,rhoc21a,0,Ta,K_CALL,flaga,alpha_expa,alpha_tanha,kb_tanha) {}
		
		virtual double operator() (double K0)  const
		{
			arg.set_nth(14,K0);
			double s=ARM_CF_BiSABR_SpreadOption_Formula::specific_shift(14);
			return 1.+ARM_CF_BiSABR_SpreadOption_Formula::value(14,arg,s);
		}
	};
	
	DistributionToInverse x(f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,T,flag,alpha_exp,alpha_tanh,kb_tanh);

	/// initial guess
	double guess = f1-f2;
	double marge=alpha1*pow(f1,beta1)+alpha2*pow(f2,beta2);
	if(proba>0.0001)
	{
		return Inverse(x,Inverse::REAL)(proba,guess,marge,1e-12);  // suppose that  Y is always postive
	}
	else
	{	double eps=0.0001;
		double p0=Inverse(x,Inverse::REAL)(eps,guess,marge,1e-12);  
		return p0*proba/eps;

	}
}


double BiSABR_smile::gaussian_to_distribution(double f1,double alpha1,double beta1, double rho1, double nu1,
								 double f2,double alpha2,double beta2, double rho2, double nu2,
								 double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag,
								 double alpha_exp,double alpha_tanh,double kb_tanh)
{
	return inverse_distribution(f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,ARM_GaussianAnalytics::cdfNormal(x),T,flag,
		alpha_exp,alpha_tanh,kb_tanh);
}


////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////

double BiSABR_smile::gaussian_to_distribution(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<18)||(argsize>18))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_smile::gaussian_to_distribution  : bad argsize");
	}
	double liminf=-5;
	double limsup=5;
	double averagingperiod=1;
	if ((x>liminf)&&(x<limsup))
	{
		return inverse_distribution(
			Underlying[0],	//INDEX1
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5],	//INDEX2
			Underlying[6],	//ALPHA2
			Underlying[7],	//BETA2
			Underlying[8],	//RHO2
			Underlying[9],	//NU2
			Underlying[10],	//RHOS
			Underlying[11],	//RHOV
			Underlying[12],	//RHOC12
			Underlying[13],	//RHOC21
			ARM_GaussianAnalytics::cdfNormal(x),	//proba
			t,				//T
			Underlying[14],	//FLAG
			Underlying[15],	//ALPHA_EXP
			Underlying[16],	//ALPHA_TANH
			Underlying[17]	//KB_TANH
			);
	}
	else if (x>=limsup)
	{
		double f1=log(inverse_distribution(
			Underlying[0],	//INDEX1
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5],	//INDEX2
			Underlying[6],	//ALPHA2
			Underlying[7],	//BETA2
			Underlying[8],	//RHO2
			Underlying[9],	//NU2
			Underlying[10],	//RHOS
			Underlying[11],	//RHOV
			Underlying[12],	//RHOC12
			Underlying[13],	//RHOC21
			ARM_GaussianAnalytics::cdfNormal(limsup-averagingperiod),	//proba
			t,				//T
			Underlying[14],	//FLAG
			Underlying[15],	//ALPHA_EXP
			Underlying[16],	//ALPHA_TANH
			Underlying[17]	//KB_TANH
			));
		double f2=log(inverse_distribution(
			Underlying[0],	//INDEX1
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5],	//INDEX2
			Underlying[6],	//ALPHA2
			Underlying[7],	//BETA2
			Underlying[8],	//RHO2
			Underlying[9],	//NU2
			Underlying[10],	//RHOS
			Underlying[11],	//RHOV
			Underlying[12],	//RHOC12
			Underlying[13],	//RHOC21
			ARM_GaussianAnalytics::cdfNormal(limsup),	//proba
			t,				//T
			Underlying[14],	//FLAG
			Underlying[15],	//ALPHA_EXP
			Underlying[16],	//ALPHA_TANH
			Underlying[17]	//KB_TANH
			));
		return exp((x-limsup)*(f2-f1)/averagingperiod+f2);
	}
	else // if (x<=liminf)
	{
		double f1=log(-inverse_distribution(
			Underlying[0],	//INDEX1
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5],	//INDEX2
			Underlying[6],	//ALPHA2
			Underlying[7],	//BETA2
			Underlying[8],	//RHO2
			Underlying[9],	//NU2
			Underlying[10],	//RHOS
			Underlying[11],	//RHOV
			Underlying[12],	//RHOC12
			Underlying[13],	//RHOC21
			ARM_GaussianAnalytics::cdfNormal(liminf),	//proba
			t,				//T
			Underlying[14],	//FLAG
			Underlying[15],	//ALPHA_EXP
			Underlying[16],	//ALPHA_TANH
			Underlying[17]	//KB_TANH
			));
		double f2=log(-inverse_distribution(
			Underlying[0],	//INDEX1
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5],	//INDEX2
			Underlying[6],	//ALPHA2
			Underlying[7],	//BETA2
			Underlying[8],	//RHO2
			Underlying[9],	//NU2
			Underlying[10],	//RHOS
			Underlying[11],	//RHOV
			Underlying[12],	//RHOC12
			Underlying[13],	//RHOC21
			ARM_GaussianAnalytics::cdfNormal(liminf+averagingperiod),	//proba
			t,				//T
			Underlying[14],	//FLAG
			Underlying[15],	//ALPHA_EXP
			Underlying[16],	//ALPHA_TANH
			Underlying[17]	//KB_TANH
			));
		return -exp((x-liminf)*(f2-f1)/averagingperiod+f1);
	}
	
}
double BiSABR_smile::distribution_to_gaussian(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<15)||(argsize>15))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_smile::gaussian_to_distribution  : bad argsize");
	}
	double f1= Underlying[0];	//INDEX1
	double alpha1=Underlying[1]; //ALPHA1
	double beta1=Underlying[2];	//BETA1
	double rho1=Underlying[3];	//RHO1
	double nu1=Underlying[4];	//NU1
	double f2= Underlying[5];	//INDEX2
	double alpha2=Underlying[6]; //ALPHA2
	double beta2=Underlying[7];	//BETA2
	double rho2=Underlying[8];	//RHO2
	double nu2=Underlying[9];	//NU2
	double rhos=Underlying[10];	//RHOS
	double rhov=Underlying[11];	//RHOV
	double rhoc12=Underlying[12];	//RHOC12
	double rhoc21=Underlying[13];	//RHOC21
	int flag=Underlying[14];	//T
	
	double X1=pow(f1,beta1);double X2=pow(f2,beta2);
	double sigma=sqrt(X1*X1*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*X1*X2*alpha1*alpha2*rhos);
	double averagingperiod=sqrt(t)*sigma;
	double liminf=f1-f2-4.*averagingperiod;
	double limsup=f1-f2+4.*averagingperiod;

	double F1,F2;
	
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
		F1=ARM_GaussianAnalytics::cdfNormal_Inv(probability_distribution(
			Underlying,	//INDEX1
			limsup-averagingperiod,
			t
			));
		F2=ARM_GaussianAnalytics::cdfNormal_Inv(probability_distribution(
			Underlying,	//INDEX1
			limsup,
			t
			));
		return (x-limsup)*(F2-F1)/averagingperiod+F2;
	}
	else /// last case (x<=liminf)
	{
		F1=ARM_GaussianAnalytics::cdfNormal_Inv(probability_distribution(
			Underlying,	//INDEX1
			liminf-averagingperiod,
			t
			));
		F2=ARM_GaussianAnalytics::cdfNormal_Inv(probability_distribution(
			Underlying,	//INDEX1
			liminf,
			t
			));
		return (x-liminf)*(F1-F2)/averagingperiod+F2;
	}
}

double BiSABR_smile::distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying, double x, double t)
{
	double s=ARM_CF_SABR_VanillaOption_Formula::specific_shift(i);
	ArgumentList* Underlying1;
	ArgumentList* Underlying2;
	switch(i)
	{	
	case ARM_CF_BiSABR_SpreadOption_Formula::INDEX1 :										/// FORWARD1
		{
			Underlying1=new ArgumentList(Underlying[0]-s,Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0]+s,Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::ALPHA1 :											/// ALPHA1 
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1]-s,Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1]+s,Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::BETA1 :											/// BETA1
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2]-s,Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2]+s,Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::RHO1 :											/// RHO1
		{		
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3]-s,Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3]+s,Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::NU1 :											/// NU1
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4]-s,Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4]+s,Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::INDEX2 :										/// FORWARD1
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5]-s,Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5]+s,Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::ALPHA2 :											/// ALPHA1 
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6]-s,Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6]+s,Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::BETA2 :											/// BETA1
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7]-s,Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7]+s,Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::RHO2 :											/// RHO1
		{		
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8]-s,Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8]+s,Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::NU2 :											/// NU1
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9]-s,Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9]+s,Underlying[10],Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::RHOS :											/// RHOS
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10]-s,Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10]+s,Underlying[11],Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::RHOV :											/// RHOV
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11]-s,Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11]+s,Underlying[12],Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::RHOC12 :											/// RHOC12
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12]-s,Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12]+s,Underlying[13],t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::RHOC21 :											/// RHOC21
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13]-s,t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13]+s,t,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
	case ARM_CF_BiSABR_SpreadOption_Formula::T :											/// RHOC21
		{
			Underlying1=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t-s,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			Underlying2=new ArgumentList(Underlying[0],Underlying[1],Underlying[2],Underlying[3],Underlying[4],Underlying[5],Underlying[6],Underlying[7],Underlying[8],Underlying[9],Underlying[10],Underlying[11],Underlying[12],Underlying[13],t+s,Underlying[14],
				Underlying[15],Underlying[16],Underlying[17]);
			break;
		}
		
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_smile::distribution_to_gaussian_first_derivative : incorrect input");
	};
	double v1=distribution_to_gaussian( *Underlying1,  x,  t);
	double v2=distribution_to_gaussian( *Underlying2,  x,  t);
	delete Underlying1;
	delete Underlying2;
	return (v2-v1)/(2.*s);
				
}



double BiSABR_smile::quantile(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<18)||(argsize>18))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_smile::quantile  : bad argsize");
	}
	return inverse_distribution(
			Underlying[0],	//INDEX1
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5],	//INDEX2
			Underlying[6],	//ALPHA2
			Underlying[7],	//BETA2
			Underlying[8],	//RHO2
			Underlying[9],	//NU2
			Underlying[10],	//RHOS
			Underlying[11],	//RHOV
			Underlying[12],	//RHOC12
			Underlying[13],	//RHOC21
			x,	//proba
			t,				//T
			Underlying[14],	//FLAG
			Underlying[15],	//ALPHA_EXP
			Underlying[16],	//ALPHA_TANH
			Underlying[17]	//KB_TANH
		);		
}

double BiSABR_smile::probability_density(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<18)||(argsize>18))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_smile::density  : bad argsize");
	}
	ArgumentList arg(
			Underlying[0],	//INDEX1
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5],	//INDEX2
			Underlying[6],	//ALPHA2
			Underlying[7],	//BETA2
			Underlying[8],	//RHO2
			Underlying[9],	//NU2
			Underlying[10],	//RHOS
			Underlying[11],	//RHOV
			Underlying[12],	//RHOC12
			Underlying[13],	//RHOC21
			x,				//proba
			t,				// T
			K_CALL,			//CallorPut
			Underlying[14],	//FLAG
			Underlying[15],	//ALPHA_EXP
			Underlying[16],	//ALPHA_TANH
			Underlying[17]	//KB_TANH
		);
	if(x<=0) return 0.0;
	double s=ARM_CF_BiSABR_SpreadOption_Formula::specific_shift(1);
	return fabs(ARM_CF_BiSABR_SpreadOption_Formula::value(14,14,arg,s,s));
}

double BiSABR_smile::probability_distribution(const ArgumentList& Underlying, double x, double t)
{
	int argsize=Underlying.size();
	if ((argsize<18)||(argsize>18))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_smile::distribution  : bad argsize");
	}
	ArgumentList arg(
			Underlying[0],	//INDEX1
			Underlying[1],	//ALPHA1
			Underlying[2],	//BETA1
			Underlying[3],	//RHO1
			Underlying[4],	//NU1
			Underlying[5],	//INDEX2
			Underlying[6],	//ALPHA2
			Underlying[7],	//BETA2
			Underlying[8],	//RHO2
			Underlying[9],	//NU2
			Underlying[10],	//RHOS
			Underlying[11],	//RHOV
			Underlying[12],	//RHOC12
			Underlying[13],	//RHOC21
			x,				//K
			t,				//T
			K_CALL,			//Call or Put
			Underlying[14],	//FLAG
			Underlying[15],	//ALPHA_EXP
			Underlying[16],	//ALPHA_TANH
			Underlying[17]	//KB_TANH
		);	
	double s=ARM_CF_BiSABR_SpreadOption_Formula::specific_shift(14);
	double digitale=-ARM_CF_BiSABR_SpreadOption_Formula::value(14,arg,s);
	if(digitale>1.0) return 1.0;
	return digitale;
}

double BiSABR_smile::probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x, double t)
{
	double v= distribution_to_gaussian( Underlying,x,t);
	double derive=distribution_to_gaussian_first_derivative(i,Underlying,x,t);
	return exp(-v*v/2.)*derive*ARM_NumericConstants::ARM_INVSQRT2PI;
}

ArgumentList*  BiSABR_smile:: HomotheticTransformation(const ArgumentList* Underlying, double positivenumber)
{
	return new ArgumentList(
		positivenumber*(*Underlying)[0],
		pow(positivenumber,1.-(*Underlying)[2])*(*Underlying)[1],
		(*Underlying)[2],
		(*Underlying)[3],
		(*Underlying)[4],
		(*Underlying)[5],
		(*Underlying)[6],
		(*Underlying)[7],
		(*Underlying)[8],
		(*Underlying)[9],
		(*Underlying)[10],
		(*Underlying)[11],
		(*Underlying)[12],
		(*Underlying)[13],
		(*Underlying)[14],
		positivenumber*(*Underlying)[15],
		(*Underlying)[16],
		positivenumber*(*Underlying)[17]
	
		);
}

ArgumentList*  BiSABR_smile:: HomotheticTransformation2(const ArgumentList* Underlying, double positivenumber)
{
	return new ArgumentList(
		(*Underlying)[0],
		(*Underlying)[1],
		(*Underlying)[2],
		(*Underlying)[3],
		(*Underlying)[4],
		positivenumber*(*Underlying)[5],
		pow(positivenumber,1.-(*Underlying)[7])*(*Underlying)[6],
		(*Underlying)[7],
		(*Underlying)[8],
		(*Underlying)[9],
		(*Underlying)[10],
		(*Underlying)[11],
		(*Underlying)[12],
		(*Underlying)[13],
		(*Underlying)[14],
		positivenumber*(*Underlying)[15],
		(*Underlying)[16],
		positivenumber*(*Underlying)[17]
		);
}



CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/