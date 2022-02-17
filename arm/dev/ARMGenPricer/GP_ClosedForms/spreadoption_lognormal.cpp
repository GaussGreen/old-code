/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_lognormal.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"


#include <cmath>

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"

#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/spreadoption_lognormal.h"
#include "gpclosedforms/spreadoption_lognormal_formula.h"			/// for the flags
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/vanilla_bs.h"

#include <glob/expt.h>


inline double cmin(double x, double y) {return ((x<=y) ? x : y);}
inline double cmax(double x, double y) {return ((x<=y) ? y : x);}


CC_BEGIN_NAMESPACE(ARM)
//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions for rho different of 1
///
//////////////////////////////////////////////////////////////////////////////////////////


double SpreadDigitalCall_Integrand(double S1,double S2,double sig1,double sig2,double rho,double k,double t,double x)
{
	double denominator=(sqrt(1-rho*rho)*sig1*sqrt(t));
	double arglog=(k+S2*exp(-0.5*sig2*sig2*t+sig2*x*sqrt(t)))/(S1*exp(-0.5*sig1*sig1*t));
	double s1=log(arglog);
	double s2=-rho*sig1*sqrt(t)*x;
	double val=NormalCDF(-(log(arglog)-rho*sig1*sqrt(t)*x)/denominator);
	return val;
}


double SpreadDigitalCall_Discriminant(double S1,double S2,double sig1,double sig2,double rho,double k,double t)
{
	return SpreadDigitalCall_Integrand(S1,S2,sig1,sig2,rho,k,t,0);
}

			
/// k is assumed to be positive here in SpreadDigitalCall_Integrand
double SpreadDigitalCall(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int n)
{
	if(t<=0)
	{
		if((S1-S2-k)>=0)
		{
			return 1.;
		}
		else
		{
			return 0.;
		}
	}

	if(S2<K_NEW_DOUBLE_TOL)
	{
		double totalvolatility = sig1*sqrt(t);
		double bondprice = 1.0;
		double strike = k;
		double CallPut = 1;
		return DigitalBlackSholes_Formula(S1,
							totalvolatility,
							bondprice,
							strike,
							CallPut);	
	}

	else if (S1<K_NEW_DOUBLE_TOL)
	{
		double totalvolatility = sig2*sqrt(t);
		double bondprice = 1.0;
		double strike = -k;
		double CallPut = -1;
		return DigitalBlackSholes_Formula(S2,
							totalvolatility,
							bondprice,
							strike,
							CallPut);

	}

	else
	{
		if (rho>0.95)
		{
			/// la convention du signe de K est differente pour  SpreadDigitalCall_Integrand et CorrelationRobust_SpreadDigitalCall
			/// a fixer
			return 1.-CorrelationRobust_SpreadDigitalCall(S2,S1,sig2,sig1,rho,-k,t,n);

		}
		else
		{
			GaussHermite_Coefficients c(n);
			double Sum=0;
			double x;
			for(int i=0;i<n;i++){
				x=c.get_point(i)*ARM_NumericConstants::ARM_SQRT_2;
				Sum+= SpreadDigitalCall_Integrand(S1,S2,sig1,sig2,rho,k,t,x)*c.get_weight(i);
			}
			return Sum*ARM_NumericConstants::ARM_INVSQRTPI;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions for rho different as close of  1 as we want
///
//////////////////////////////////////////////////////////////////////////////////////////
double SpreadDigitalCall_Integrand_1(double S1,double S2,double sig1,double sig2,double rho,double k,double t,double x)
{
	double sqrtt=sqrt(t);
	double denominator=(sqrt(1-rho*rho)*sig2*sqrtt);
	double arglog=(S1*exp(0.5*(sig2*sig2-sig1*sig1)*t+sig1*x*sqrtt)-k*exp(sig2*sig2*t/2.))/S2;
	double argnorm=(log(arglog)-rho*sig2*sqrtt*x);
	if (rho>=1.)
	{
		if(argnorm>0) return 1.;
		if(argnorm<0) return 0.;
		return 0.5;
	}
	else
	{
		return NormalCDF((log(arglog)-rho*sig2*sqrtt*x)/denominator);
	}
	
}

double CorrelationRobust_SpreadDigitalCall(double S10,double S20,double sig10,double sig20,double rho0,double k0,double t0,int n)
{
	///     determination of x1 and x2 
	struct EquationToSolve : public DoubleToDoubleFunc 
	{
		
		virtual double operator() (double x)  const
		{
			double sqrtt=sqrt(t);
			double arglog=(S1*exp(0.5*(sig2*sig2-sig1*sig1)*t+sig1*x*sqrtt)-k*exp(sig2*sig2*t/2.))/S2;
			return log(arglog)-sig2*sqrtt*x;
		}
		double S1;
		double S2;
		double sig1;
		double sig2;
		double rho;
		double t;
		double k;
		
		EquationToSolve(double S1a,double S2a,double sig1a,double sig2a,double rhoa,double ta,double ka):
		S1(S1a),S2(S2a),sig1(sig1a),sig2(sig2a),rho(rhoa),t(ta),k(ka) {}
		
	};
	
	EquationToSolve eq(S10,S20,sig10,sig20,rho0,t0,k0);
	
	double x_Limit=12., ILimit,xminimum,yminimum,x1,x2,xini_first,xini_last;
	double I1=0.,I2=0.,I3=0.,x2down,x2up,Sum=0.,x;
	int i;
	GaussLegendre_Coefficients c_root(n);
	/// estimation of the limit of the roots
	if(fabs(sig10-sig20)<0.001*sig10)
	{
		if (fabs(k0)<0.0001*(S10+S20)) /// on reprends le cas si dessous en floorant le k 
		{
			double k0_a=-0.0001*(S10+S20);
			if(S10<S20)
			{
				x1=(sig10*sig10*t0/2.-log((S10-S20)/k0_a))/(sig10*sqrt(t0)); /// k should e different from 0 !
				ILimit=NormalCDF(x1);
				Sum=0.;
				GaussLegendre_Coefficients c1(&c_root,-x_Limit,x1);
				for(i=0;i<n;i++){
					x=c1.get_point(i);
					Sum+= (1.-SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0_a,t0,x))*exp(-x*x/2.)*c1.get_weight(i);
				}
				I1=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
				Sum=0.;
				GaussLegendre_Coefficients c2(&c_root,x1,x_Limit);
				for(i=0;i<n;i++){
					x=c2.get_point(i);
					Sum+= SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0_a,t0,x)*exp(-x*x/2.)*c2.get_weight(i);
				}
				I2=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
				return cmax(ILimit-I1+I2,0);
			}
			else
			{
				ILimit=1.;
				Sum=0.;
				GaussLegendre_Coefficients c1(&c_root,-x_Limit,x_Limit);
				for(i=0;i<n;i++){
					x=c1.get_point(i);
					Sum+= (1.-SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0_a,t0,x))*exp(-x*x/2.)*c1.get_weight(i);
				}
				I1=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
				I2=0.;
				return cmax(ILimit-I1+I2,0);
			}
		}
				
		
		else
		{
			
			if(S10<S20)
			{
				x1=(sig10*sig10*t0/2.-log((S10-S20)/k0))/(sig10*sqrt(t0)); /// k should e different from 0 !
				ILimit=NormalCDF(x1);
				Sum=0.;
				GaussLegendre_Coefficients c1(&c_root,-x_Limit,x1);
				for(i=0;i<n;i++){
					x=c1.get_point(i);
					Sum+= (1.-SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x))*exp(-x*x/2.)*c1.get_weight(i);
				}
				I1=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
				Sum=0.;
				GaussLegendre_Coefficients c2(&c_root,x1,x_Limit);
				for(i=0;i<n;i++){
					x=c2.get_point(i);
					Sum+= SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x)*exp(-x*x/2.)*c2.get_weight(i);
				}
				I2=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
				return cmax(ILimit-I1+I2,0);
			}
			else
			{
				ILimit=1.;
				Sum=0.;
				GaussLegendre_Coefficients c1(&c_root,-x_Limit,x_Limit);
				for(i=0;i<n;i++){
					x=c1.get_point(i);
					Sum+= (1.-SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x))*exp(-x*x/2.)*c1.get_weight(i);
				}
				I1=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
				I2=0.;
				return cmax(ILimit-I1+I2,0);
			}
		}
	}
	else
	{
		
		if (fabs(k0)<0.0001*(S10+S20))
		{
			/// computation of the only root
			if (sig10-sig20>0)
			{
				x1= -(log(S10/S20)-0.5*t0*(sig10*sig10-sig20*sig20))/(sqrt(t0)*fabs(sig10-sig20));
				if(x1>x_Limit) x1=x_Limit;
				if(x1<-x_Limit) x1=-x_Limit;
				/// corr =1 limit
				ILimit=1-NormalCDF(x1);
				Sum=0.;
				GaussLegendre_Coefficients c1(&c_root,x1,x_Limit);
				for(i=0;i<n;i++){
					x=c1.get_point(i);
					Sum+= (1-SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x))*exp(-x*x/2.)*c1.get_weight(i);
				}
				I1=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
				Sum=0.;
				GaussLegendre_Coefficients c2(&c_root,-x_Limit,x1);
				for(i=0;i<n;i++){
					x=c2.get_point(i);
					Sum+= SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x)*exp(-x*x/2.)*c2.get_weight(i);
				}
				I2=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
			}
			else
			{
				x1= (log(S10/S20)-0.5*t0*(sig10*sig10-sig20*sig20))/(sqrt(t0)*fabs(sig10-sig20));
				if(x1>x_Limit) x1=x_Limit;
				if(x1<-x_Limit) x1=-x_Limit;
				/// corr =1 limit
				ILimit=NormalCDF(x1);
				Sum=0.;
				GaussLegendre_Coefficients c1(&c_root,-x_Limit,x1);
				for(i=0;i<n;i++){
					x=c1.get_point(i);
					Sum+= (1-SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x))*exp(-x*x/2.)*c1.get_weight(i);
				}
				I1=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
				Sum=0.;
				GaussLegendre_Coefficients c2(&c_root,x1,x_Limit);
				for(i=0;i<n;i++){
					x=c2.get_point(i);
					Sum+= SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x)*exp(-x*x/2.)*c2.get_weight(i);
				}
				I2=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
			}
			return cmax(ILimit-I1+I2,0);
			
		}
		else
		{
			
			xini_first=-(log(-S20/k0)-sig20*sig20*t0/2.)/(sig20*sqrt(t0));
			
			xini_last=(log(S20/S10)+(sig10*sig10-sig20*sig20)*t0/2.)/((sig10-sig20)*sqrt(t0));
			
			
			/// estimation ofthe x atthe minimum for the double root case
			
			
			
			if ((sig10-sig20)>0)	/// 2 roots 
			{
				/// estimation ofthe x atthe minimum for the double root case
				
				xminimum=(t0*sig10*sig10-2.*log(-S10/(k0*sig20)*(sig10-sig20)))/(2.*sqrt(t0)*sig10);
				
				/// estimation of the value orf the arg of the N[] function at xminimum
				
				yminimum= 0.5*t0*(sig20-sig10)*sig20+sig20*log(S10*(-sig10+sig20)/(k0*sig20))/sig10+log(-k0*sig10/(S20*(sig10-sig20)));
				
				if (yminimum>=0)		/// no real roots the rho=1 limit is therefore 0
				{
					
					Sum=0.;
					GaussLegendre_Coefficients c1(&c_root,-x_Limit,x_Limit);
					for(i=0;i<n;i++){
						x=c1.get_point(i);
						Sum+= SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x)*exp(-x*x/2.)*c1.get_weight(i);
					}
					return cmax(Sum,0)*ARM_NumericConstants::ARM_INVSQRT2PI;
					
				}
				else		/// two real roots
				{
					/// so we can formulate the mean as a start 
					
					x1=Inverse(eq,Inverse::REAL)(0.,(xini_first+xminimum)/2.,(xminimum-xini_first)/2.,1e-12);
					
					x2=Inverse(eq,Inverse::REAL)(0.,(xini_last+xminimum)/2.,(xini_last-xminimum)/2.,1e-12);
					
					if(x1<-x_Limit)
					{
						I1=0.;
						x2down=-x_Limit;
						
					}
					else
					{
						Sum=0.;
						GaussLegendre_Coefficients c1(&c_root,-x_Limit,x1);
						for(i=0;i<n;i++){
							x=c1.get_point(i);
							Sum+= (1.-SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x))*exp(-x*x/2.)*c1.get_weight(i);
						}
						I1=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
						x2down=x1;
					}
					
					if(x2>x_Limit)
					{
						I3=0.;
						x2up=x_Limit;
						
					}
					else
					{
						Sum=0.;
						GaussLegendre_Coefficients c3(&c_root,x2,x_Limit);
						for(i=0;i<n;i++){
							x=c3.get_point(i);
							Sum+= (1.-SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x))*exp(-x*x/2.)*c3.get_weight(i);
						}
						I3=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
						x2up=x2;
					}
					if(x2down<x2up)
					{	
						Sum=0.;
						GaussLegendre_Coefficients c2(&c_root,x2down,x2up);
						for(i=0;i<n;i++){
							x=c2.get_point(i);
							Sum+=SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x)*exp(-x*x/2.)*c2.get_weight(i);
						}
						I2=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
					}
					else
					{
						I2=0.;
					}
					ILimit=NormalCDF(x1)+1.-NormalCDF(x2);
					
					return cmax(ILimit + I2 - I1 - I3,0);
				}
			}
			else		// 1 root
			{
				x1=Inverse(eq,Inverse::REAL)(0.,(xini_first+xini_last)/2.,0.2,1e-12);
				
				if(x1<-x_Limit)
				{
					ILimit=0;
					I1=0;
					Sum=0.;
					GaussLegendre_Coefficients c2(&c_root,-x_Limit,x_Limit);
					for(i=0;i<n;i++){
						x=c2.get_point(i);
						Sum+= SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x)*exp(-x*x/2.)*c2.get_weight(i);
					}
					I2=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
					
				}
				else if (x1>x_Limit)
				{
					ILimit=1.;
					Sum=0.;
					
					GaussLegendre_Coefficients c1(&c_root,-x_Limit,x_Limit);
					for(i=0;i<n;i++){
						x=c1.get_point(i);
						Sum+= (1.-SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x))*exp(-x*x/2.)*c1.get_weight(i);
					}
					I1=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
					I2=0;
				}
				else
				{	
					ILimit=NormalCDF(x1);
					Sum=0.;
					GaussLegendre_Coefficients c1(&c_root,-x_Limit,x1);
					for(i=0;i<n;i++){
						x=c1.get_point(i);
						Sum+= (1.-SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x))*exp(-x*x/2.)*c1.get_weight(i);
					}
					I1=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
					Sum=0.;
					GaussLegendre_Coefficients c2(&c_root,x1,x_Limit);
					for(i=0;i<n;i++){
						x=c2.get_point(i);
						Sum+= SpreadDigitalCall_Integrand_1(S10,S20,sig10,sig20,rho0,k0,t0,x)*exp(-x*x/2.)*c2.get_weight(i);
					}
					I2=Sum*ARM_NumericConstants::ARM_INVSQRT2PI;
				}
				return cmax(ILimit - I1+I2,0);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////


double Vega1SpreadOption_Integrand(double S1,double S2,double sig1,double sig2,double rho,double k,double t,double x)
{
	double h1=exp(rho*sig1*sig2*t+sig2*sqrt(t)*x);
	double q1=sig1*sig1*t/2.+log((k+exp(-sig2*sig2*t/2.+sig2*sqrt(t)*x)*S2)/S1);
	double q2=-sig1*sig1*t/2.+log((k+exp(rho*sig1*sig2*t-sig2*sig2*t/2.+sig2*sqrt(t)*x)*S2)/S1);
	double q3=sig1*sig1*t/2.-rho*sig1*sig2*t+log((k+exp(sig2*sig2*t/2.+sig2*sqrt(t)*x)*S2)/S1);
	double abase=(q1-rho*sig1*sqrt(t)*x);	
	double a1=exp(abase*abase/(2*(-1+rho*rho)*sig1*sig1*t));	//  sauf si rho =+ou- 1
	abase=(q2-rho*sig1*sqrt(t)*x);
	double a2=exp(abase*abase/(2*(-1+rho*rho)*sig1*sig1*t));
	abase=(q3-rho*sig1*sqrt(t)*x);
	double a3=exp(abase*abase/(2*(-1+rho*rho)*sig1*sig1*t));
	double b1=exp(sig2*sig2*t/2.);

	return (a1*k*(-q1+sig1*sig1*t)-a3*S2*(q3+sig1*(-sig1+rho*sig2)*t)+
		(a2*S1*(q2*(b1*k+h1*S2)+sig1*(b1*k*sig1+h1*S2*(sig1-rho*sig2))*t))/
		(b1*k+h1*S2))/
		(sqrt((1-rho*rho)*t)*sig1*sig1*ARM_NumericConstants::ARM_SQRT_2_PI);
}



double Vega1SpreadOptionCall_aux(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int n)
{
	GaussHermite_Coefficients c(n);
	double Sum=0;
	double x;
	for(int i=0;i<n;i++){
		x=c.get_point(i)*ARM_NumericConstants::ARM_SQRT_2;
		Sum+= Vega1SpreadOption_Integrand(S1,S2,sig1,sig2,rho,k,t,x)*c.get_weight(i);
	}
	return Sum*ARM_NumericConstants::ARM_INVSQRTPI;
}

double Vega1SpreadOptionCall(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int n)
{
	double d00=SpreadDigitalCall_Discriminant(S1,S2,sig1,sig2,rho,k,t);
	double d12=SpreadDigitalCall_Discriminant(S2,S1,sig2,sig1,rho,k,t);

	if(d00>d12) return Vega1SpreadOptionCall_aux(S1,S2,sig1,sig2,rho,k,t,n);
	else return Vega1SpreadOptionCall_aux(S2,S1,sig2,sig1,rho,k,t,n);
}



double Vega2SpreadOption_Integrand(double S1,double S2,double sig1,double sig2,double rho,double k,double t,double x)
{
	
	double q1=sig1*sig1*t/2.+log((k+exp(-sig2*sig2*t/2.+sig2*sqrt(t)*x)*S2)/S1);
	double q2=-sig1*sig1*t/2.+log((k+exp(rho*sig1*sig2*t-sig2*sig2*t/2.+sig2*sqrt(t)*x)*S2)/S1);
	double q3=sig1*sig1*t/2.-rho*sig1*sig2*t+log((k+exp(sig2*sig2*t/2.+sig2*sqrt(t)*x)*S2)/S1);
	double abase=(q1-rho*sig1*sqrt(t)*x);	
	double a1=abase*abase/(2*(-1+rho*rho)*sig1*sig1*t);
	abase=(q2-rho*sig1*sqrt(t)*x);
	double a2=abase*abase/(2*(-1+rho*rho)*sig1*sig1*t);
	abase=(q3-rho*sig1*sqrt(t)*x);
	double a3=abase*abase/(2*(-1+rho*rho)*sig1*sig1*t);
	double d1=exp(sig2*sig2*t/2.)*k+exp(sig2*sqrt(t)*x)*S2;
	double d2=exp(sig2*sig2*t/2.)*k+exp(rho*sig1*sig2*t+sig2*sqrt(t)*x)*S2;
	double d3=k+exp(sig2*sig2*t/2.+sig2*sqrt(t)*x)*S2;

	return S2*(
		-(exp(a1+sig2*sqrt(t)*x)*k*sqrt(t)*(sig2*sqrt(t)-x))/d1
		-exp(a2+rho*sig1*sig2*t+sig2*sqrt(t)*x)*S1*sqrt(t)*((rho*sig1-sig2)*sqrt(t)+x)/d2
		+exp(a3)*(-k*rho*sig1*t+exp(sig2*sig2*t/2.+sig2*sqrt(t)*x)*S2*sqrt(t)*((-rho*sig1+sig2)*sqrt(t)+x))/d3
		)/(sqrt((1-rho*rho)*t)*sig1*ARM_NumericConstants::ARM_SQRT_2_PI);
}

double Vega2SpreadOptionCall_aux(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int n)
{
	GaussHermite_Coefficients c(n);
	double Sum=0;
	double x;
	for(int i=0;i<n;i++){
		x=c.get_point(i)*ARM_NumericConstants::ARM_SQRT_2;
		Sum+= Vega2SpreadOption_Integrand(S1,S2,sig1,sig2,rho,k,t,x)*c.get_weight(i);
	}
	return Sum*ARM_NumericConstants::ARM_INVSQRTPI;
}

double Vega2SpreadOptionCall(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int n)
{
	double d00=SpreadDigitalCall_Discriminant(S1,S2,sig1,sig2,rho,k,t);
	double d12=SpreadDigitalCall_Discriminant(S2,S1,sig2,sig1,rho,k,t);

	if(d00>d12) return Vega2SpreadOptionCall_aux(S1,S2,sig1,sig2,rho,k,t,n);
	else return Vega2SpreadOptionCall_aux(S2,S1,sig2,sig1,rho,k,t,n);
}



//////////////////////////////////////////////////////////////////////////////////////////
///
///   Call/Put  k<0, k>0 Nuance Introductors
///
//////////////////////////////////////////////////////////////////////////////////////////

double SpreadDigitalOption(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int n)
{
	switch (callput)
	{
	case K_CALL :
		{
			if (k>0) 
			{
				return SpreadDigitalCall(S1,S2,sig1,sig2,rho,k,t,n );
			}
			else
			{
				return 1.- SpreadDigitalCall(S2,S1,sig2,sig1,rho,-k,t,n );
			}
			
			break;
		}
	case K_PUT :
		{
			if (k>0) 
			{
				return 1.-SpreadDigitalCall(S1,S2,sig1,sig2,rho,k,t,n );
			}
			else
			{
				return SpreadDigitalCall(S2,S1,sig2,sig1,rho,-k,t,n );
			}
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SpreadDigitalOption : callput , bad input :");
			break;
		}
	}
	
}



double Vega1SpreadOption(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int n)
{
	switch (callput)
	{
	case K_CALL :
		{
			if (k>=0) 
			{
				return Vega1SpreadOptionCall(S1,S2,sig1,sig2,rho,k,t,n );
			}
			else
			{
				return - Vega1SpreadOptionCall(S2,S1,sig2,sig1,rho,-k,t,n );
			}
			
			break;
		}
	case K_PUT :
		{
			if (k>=0) 
			{
				return -Vega1SpreadOptionCall(S1,S2,sig1,sig2,rho,k,t,n );
			}
			else
			{
				return Vega1SpreadOptionCall(S2,S1,sig2,sig1,rho,-k,t,n );
			}
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Vega1SpreadOption : callput , bad input :");
			break;
		}
	}
	
}



double Vega2SpreadOption(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int n)
{
	switch (callput)
	{
	case K_CALL :
		{
			if (k>=0) 
			{
				return Vega2SpreadOptionCall(S1,S2,sig1,sig2,rho,k,t,n );
			}
			else
			{
				return -Vega2SpreadOptionCall(S2,S1,sig2,sig1,rho,-k,t,n );
			}
			
			break;
		}
	case K_PUT :
		{
			if (k>=0) 
			{
				return -Vega2SpreadOptionCall(S1,S2,sig1,sig2,rho,k,t,n );
			}
			else
			{
				return Vega2SpreadOptionCall(S2,S1,sig2,sig1,rho,-k,t,n );
			}
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Vega2SpreadOption : callput , bad input :");
			break;
		}
	}
	
}




///////////////////////////////////////////////////////////////////////
///  
///			Calibration Functions Functions 
///
///////////////////////////////////////////////////////////////////////

double  LogNormal_SpreadOption_Calibrate_Correlation(double S1,double S2,double sig1,double sig2,double optionprice,
											 double k,double t,int callput,int optiontype,int n)
{
struct PricingFunctionToInverse : public DoubleToDoubleFunc 
	{
		double S10;
		double S20;
		double sig10;
		double sig20;
		double k0;
		double t0;
		int callput0;
		int optiontype0;
		int n0;
		PricingFunctionToInverse(double S1a,double S2a, double sig1a, double sig2a, 
			double ka, double ta, int callputa,int optiontypea,int na):
		S10(S1a),S20(S2a),sig10(sig1a),sig20(sig2a),k0(ka),t0(ta),callput0(callputa),optiontype0(optiontypea),n0(na)
		{}

		virtual double operator() (double rho)  const
		{
			ArgumentList a(t0,S10,S20,sig10,sig20,rho,k0,callput0,n0);
			switch (optiontype0) 
			{
			case ARM_CF_SpreadDigitalOption_Formula::DIGITALOPTION :
				{
					Power_Expression<ARM_CF_SpreadDigitalOption_Formula> y;
					return y(a);
					break;
				}
			case ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION :
				{
					Power_Expression<ARM_CF_SpreadOption_Formula> y;
					return y(a);
					break;
				}
			case ARM_CF_SpreadDigitalOption_Formula::PAYFIRST :
				{
					Power_Expression<ARM_CF_Index1Paying_SpreadDigitalOption_Formula> y;
					return y(a);
					break;
				}
			case ARM_CF_SpreadDigitalOption_Formula::PAYSECOND :
				{
					Power_Expression<ARM_CF_Index2Paying_SpreadDigitalOption_Formula> y;
					return y(a);
					break;
				}
			default :
				{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_LogNormal_SpreadOption : incorrect optiontype toggle ");
				}
			}
		}
	};

	PricingFunctionToInverse x(S1,S2,sig1,sig2,k,t,callput,optiontype,n);
	return Inverse(x,Inverse::CORRELATION)(optionprice,0.0,0.0001,1e-12);  // we assume the the correlation is always between -1 and 1 !
}

double  Smiled_LogNormal_SpreadOption_Calibrate_Correlation(double S1,double S2,double sig1,double sig2,double optionprice,
											 double k,double t,double slope1,double slope2,int callput,int optiontype,int n)
{
class PricingFunctionToInverse : public DoubleToDoubleFunc 
	{
	public: 
		double S10;
		double S20;
		double sig10;
		double sig20;
		double k0;
		double t0;
		double slope10;
		double slope20;
		int callput0;
		int optiontype0;
		int n0;
		PricingFunctionToInverse(double S1a,double S2a, double sig1a, double sig2a, 
			double ka, double ta, double slope1a,double slope2a,int callputa,int optiontypea,int na):
		S10(S1a),S20(S2a),sig10(sig1a),sig20(sig2a),k0(ka),t0(ta),slope10(slope1a),slope20(slope2a),callput0(callputa),optiontype0(optiontypea),n0(na)
		{}

		virtual double operator() (double rho) const 
		{
			ArgumentList a(t0,S10,S20,sig10,sig20,rho,k0,callput0,n0,slope10,slope20);  // on met bien les variables supplementaires (les slopes) a la fin pour traitement par les templates !
			ArgumentList a0(t0,S10,S20,sig10,sig20,rho,k0,callput0,n0);
			switch (optiontype0) 
			{
			case ARM_CF_SpreadDigitalOption_Formula::DIGITALOPTION :
		{
			Power_Expression<ARM_CF_Smiled_SpreadDigitalOption_Formula> y;
			return y(a);
			break;
		}
			case ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION :
				{
					Power_Expression<ARM_CF_SpreadOption_Formula> y;
					return y(a0);
					break;
				}
			case ARM_CF_SpreadDigitalOption_Formula::PAYFIRST :
				{
					Power_Expression<ARM_CF_Index1Paying_Smiled_SpreadDigitalOption_Formula> y;
					return y(a);
					break;
				}
			case ARM_CF_SpreadDigitalOption_Formula::PAYSECOND :
				{
					Power_Expression<ARM_CF_Index2Paying_Smiled_SpreadDigitalOption_Formula> y;
					return y(a);
					break;
				}
			default :
				{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_LogNormal_SpreadOption : incorrect optiontype toggle ");
				}
			}
		}
	};

	PricingFunctionToInverse x(S1,S2,sig1,sig2,k,t,slope1,slope2,callput,optiontype,n);
	return Inverse(x,Inverse::CORRELATION)(optionprice,0.0,0.0001,1e-12);  // we assume the the correlation is always between -1 and 1 !
}



CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/