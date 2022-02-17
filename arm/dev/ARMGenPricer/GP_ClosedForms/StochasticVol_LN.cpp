/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file merton.cpp
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

#include "gpclosedforms/StochasticVol_LN.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpnumlib/gaussiananalytics.h"


CC_BEGIN_NAMESPACE(ARM)

//////////////////////////////////////////////////////////////////////////////////////////////////////
///
///				Arithmetic Asian Option
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

/*
double AsianCallArithmetic(double S,double K,double T,double r, double sig,double averaging,int callput)
{
	if(K<=0)
	{
		if(callput==K_CALL)
		{
			return (S-K>0)?exp(-r*T)*(S-K):0;
		}
		else
		{
			return (K-S>0)?exp(-r*T)*(K-S):0;
		}
	}
	else
	{
		double r_plus_sig2over2=r+sig*sig/2.0;
		double M1,M2a;
		if(averaging>0.00001)
		{
			if(r>0.00001)
			{
				M1=S*exp(r*(T-averaging))*(exp(r*averaging)-1.0)/(r*averaging);
				M2a=S*exp(r_plus_sig2over2*(T-averaging))*(exp(r_plus_sig2over2*averaging)-1.0)/(r_plus_sig2over2*averaging);
			}
			else
			{
				M1=S*exp(r*(T-averaging));
				M2a=S*exp(r_plus_sig2over2*(T-averaging))*(exp(r_plus_sig2over2*averaging)-1.0)/(r_plus_sig2over2*averaging);
				
			}
			
		}
		else
		{
			M1=S*exp(r*(T-averaging));
			M2a=S*exp(r_plus_sig2over2*(T-averaging));
		}
		double M2=M2a*M2a-M1*M1;
		double voltotal=sqrt(M2)/M1;
		double disc=S/M1;
		return BlackSholes_Formula(M1,voltotal,disc,K,callput);
	}
}
*/




/// Exact computation of the first moment, but continuized between the discret values for T and avg
double ArithmeticAsianFirstMomentExact(double mu,double sig,double T,double avg,double dt)
{
	if((fabs(mu)<1e-10))
	{
		return 1;
	}
	else
	{
		double N=T/dt;
		double L=avg/dt;
		return exp((T-avg)*mu)*(exp((L+1)*mu*dt)-1.0)/(exp(mu*dt)-1.0)/(L+1.0);
		
	}
}

/// Exact computation of the second moment, but continuized between the discret values for T and avg
double ArithmeticAsianSecondMomentExact(double mu,double sig,double T,double avg,double dt)
{
	double sig2=sig*sig;
	double N=T/dt;
	double L=avg/dt;
	if((fabs(mu)<1e-10))
	{
		return exp((N-L)*sig2*dt)*
			(1.0	+exp((1.0+L)*dt*sig2)+	exp((2.0+L)*dt*sig2)+	2.0*L-exp(dt*sig2)*(3.0+2.0*L))
			/((L+1.0)*(L+1.0)*(exp(sig2*dt)-1.0)*(exp(sig2*dt)-1.0));
	}
	else
	{
		return exp((N-L)*(2.0*mu+sig2)*dt)*	
			(
			-1.0	-exp(dt*mu)	+2.0*exp((1.0+L)*mu*dt)	+exp(dt*(sig2+mu))	+exp(dt*(sig2+2.0*mu))	-exp((1.0+L)*dt*(sig2+2.0*mu))	+exp((2.0+L)*dt*(sig2+2.0*mu))
			-2.0*exp(((3.0+L)*mu+sig2)*dt)	+exp(((3.0+2.0*L)*mu+(1.0+L)*sig2)*dt)	-exp(((3.0+2.0*L)*mu+(2.0+L)*sig2)*dt)
			)
			/((L+1.0)*(L+1.0)*(exp(mu*dt)-1.0)*(exp((mu+sig2)*dt)-1.0)*(exp((2.0*mu+sig2)*dt)-1.0));
		
	}
}

/// approximated computation of the first moment, but continuized between the discret values for T and avg
// and valid even for dt=0 ro so
double ArithmeticAsianFirstMomentShorted(double mu,double sig,double T,double avg,double dt)
{
	if((fabs(mu)<1e-10))
	{
		return 1;
	}
	else
	{
		return exp((T-avg)*mu)*(exp((avg+dt)*mu)-1.0)/(avg*mu);
		
	}
}

/// approximated computation of the second moment, but continuized between the discret values for T and avg
// and valid even for dt=0 ro so
double ArithmeticAsianSecondMomentShorted(double mu,double sig,double T,double avg,double dt)
{
	double sig2=sig*sig;
	if((fabs(mu)<1e-10))
	{
		return -exp((T-avg)*sig2)*(1.0+exp(dt*sig2)-2.0*exp((avg+dt)*sig2)+(2.0*avg+dt)*sig2)/(avg*avg*sig2*sig2);
	}
	else
	{
		return exp((T-avg)*(sig2+2.0*mu))*	(

			(1.0	+exp(dt*sig2)	-4.0*exp((avg+dt)*mu)	-exp(dt*(sig2+2.0*mu))	+2.0*exp((avg+dt)*(2.0*mu+sig2))	)*mu
			+(1.0	+exp(dt*mu)	- 2.0*exp((avg+dt)*mu)	)*sig2

												)/(avg*avg*mu*(2.0*mu*mu+3.0*mu*sig2+sig2*sig2));
		
	}
}




/// uses the true distribution of the average  but discount with the riskfree rate
/// the discount rate may be different from the assumed "risk neutral" drift of the underlying

double AsianCallArithmetic(double S,double K,double T,double mu, double sig,double r,double avg,double dt,int callput)
{
	double M1,M2;
	if(fabs(dt)<1e-8)
	{
		M1=S*ArithmeticAsianFirstMomentShorted( mu, sig, T, avg, dt);
		M2=S*S*ArithmeticAsianSecondMomentShorted( mu, sig, T, avg, dt)-M1*M1;
	}
	else
	{
		M1=S*ArithmeticAsianFirstMomentExact( mu, sig, T, avg, dt);
		M2=S*S*ArithmeticAsianSecondMomentExact( mu, sig, T, avg, dt)-M1*M1;
	}

	double voltotal=sqrt(M2)/M1;
	double discount=exp(-r*T);

	if(K<=0)
	{
		if(callput==K_CALL)
		{
			return (M1-K>0)?discount*(M1-K):0;
		}
		else
		{
			return (K-M1>0)?discount*(K-M1):0;
		}
	}
	else
	{
		
		return BlackSholes_Formula(M1,voltotal,discount,K,callput);
	}
}


double StochasticVol_LN_Arithmetic_VanillaOption(double f,double K,double T,double mu,double sig,double r,double VolDrift,double VolOfVol,double averaging,double dt,double callput,int LegendreNb)
{
	double Vol;
	double mu2=(VolDrift-VolOfVol*VolOfVol/2.0)*T/2.0;
	double sigma2=sqrt(2.0/3.0)*VolOfVol*sqrt(T);
	if(LegendreNb>1)
	{
		ReducedGaussHermite_Coefficients c(LegendreNb);
		double Sum=0;
		int i;
		for(i=0;i<LegendreNb;i++){
			Vol=sig*exp(sigma2*c.get_point(i)+mu2);
			Sum+= AsianCallArithmetic(f,K,T,mu,Vol,r,averaging,dt,(int)callput) *c.get_weight(i);
		}
		return Sum;
	}
	else
	{	
		Vol=sig;
		return AsianCallArithmetic(f,K,T,mu,Vol,r,averaging,dt,(int)callput) ;
	}
}


/// we give up here the genrality associated with the distinction bewenn mu and r

double StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(double f,double K,double T,double r,double sig,double VolDrift,
												 double VolOfVol,double averaging,double reset,double callput,int LegendreNb)
{
	double mu=r;
	double dt=1./12.;
	if(T<0.0001)
	{
		if(callput==K_CALL)
		{
			return (reset-K>0)?exp(-r*T)*(reset-K):0;
		}
		else
		{
			return (K-reset>0)?exp(-r*T)*(K-reset):0;
		}
	}
	else
	{
		double Vol;
		double mu2=(VolDrift-VolOfVol*VolOfVol/2.0)*T/2.0;
		double sigma2=sqrt(2.0/3.0)*VolOfVol*sqrt(T);
		double K1;
		if(averaging>0.0001) K1=K-reset*(averaging-T)/averaging;;
		double f1;
		double sigmatotal=sig*sqrt(T);
		if(averaging>T)
		{
			if(LegendreNb>1)
			{
				ReducedGaussHermite_Coefficients c(LegendreNb);
				double Sum=0;
				int i;
				for(i=0;i<LegendreNb;i++){
					Vol=sig*exp(sigma2*c.get_point(i)+mu2);
					if(averaging>0.0001)
					{
						f1=T/averaging*f;
						Sum+= AsianCallArithmetic(f1,K1,T,mu,Vol,r,T,dt,(int)callput) *c.get_weight(i);
					}
					else
					{
						Sum+= exp(-r*T)*BlackSholes_Formula(f,sig*sqrt(T),1.0,K,callput) *c.get_weight(i);
					}
				}
				return Sum;
			}
			else
			{	
				if(averaging>0.0001)
				{
					f1=T/averaging*f;
					Vol=sig;
					return  AsianCallArithmetic(f1,K1,T,mu,Vol,r,T,dt,(int)callput) ;
				}
				else
				{
					return exp(-r*T)*BlackSholes_Formula(f,sig*sqrt(T),1.0,K,callput);
				}
			}
		}
		else
		{
			return StochasticVol_LN_Arithmetic_VanillaOption(f,K,T,mu,sig,r,VolDrift,VolOfVol,averaging,dt,callput,LegendreNb);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
///
///				Geometric Asian Option
///
/////////////////////////////////////////////////////////////////////////////////////////////////////:


double StochasticVol_LN_Geometric_VanillaOption(double f,double K,double T,double r,double sig,double VolDrift,double VolOfVol,double averaging,double callput,int LegendreNb)
{
	double f1,Voltotal,sigma1total;
	double mu2=(VolDrift-VolOfVol*VolOfVol/2.0)*T/2.0;
	double sigma2=sqrt(2.0/3.0)*VolOfVol*sqrt(T);
	double sigmatotal=sig*sqrt(T);
	double r1;
	if(LegendreNb>1)
	{
		ReducedGaussHermite_Coefficients c(LegendreNb);
		double Sum=0;
		int i;
		for(i=0;i<LegendreNb;i++){
			Voltotal=sigmatotal*exp(sigma2*c.get_point(i)+mu2);
			f1=f*exp(r*(T-averaging/2.0)-1/12.0*Voltotal*Voltotal*averaging/T);
			sigma1total=Voltotal*sqrt(1.0-2.0/3.0*averaging/T);
			r1=r*(1-averaging/T/2.0)-1/12.0*averaging*Voltotal*Voltotal/T/T;
			Sum+= exp(-r1*T)*BlackSholes_Formula(f1,sigma1total,1.0,K,callput)*c.get_weight(i);
		}
		return Sum;
	}
	else
	{	
		Voltotal=sigmatotal;
		f1=f*exp(r*(T-averaging/2.0)-1/12.0*Voltotal*Voltotal*averaging/T);
		r1=r*(1-averaging/T/2.0)-1/12.0*averaging/T*sig*sig;
		sigma1total=Voltotal*sqrt(1.0-2.0/3.0*averaging/T);
		return  BlackSholes_Formula(f1,sigma1total,1.0,K,callput)*exp(-r1*T);
	}
}

double StochasticVol_LN_Geometric_VanillaOption_with_Reset(double f,double K,double T,double r,double sig,double VolDrift,
												 double VolOfVol,double averaging,double reset,double callput,int LegendreNb)
{
	if(T<0.0001)
	{
		if(callput==K_CALL)
		{
			return (reset-K>0)?exp(-r*T)*(reset-K):0;
		}
		else
		{
			return (K-reset>0)?exp(-r*T)*(K-reset):0;
		}
	}
	else
	{
		double f1,Voltotal,sigma1total;
		double mu2=(VolDrift-VolOfVol*VolOfVol/2.0)*T/2.0;
		double sigma2=sqrt(2.0/3.0)*VolOfVol*sqrt(T);
		double K1=K-reset*(averaging-T)/averaging;
		double sigmatotal=sig*sqrt(T);
		if(averaging>T)
		{
			if(LegendreNb>1)
			{
				ReducedGaussHermite_Coefficients c(LegendreNb);
				double Sum=0;
				int i;
				for(i=0;i<LegendreNb;i++){
					Voltotal=sigmatotal*exp(sigma2*c.get_point(i)+mu2);
					f1=f*T/averaging*exp(0.5*(r-Voltotal*Voltotal/(6.0*T))*T);
					sigma1total=sqrt(exp(Voltotal*Voltotal/3.0) - 1.0);
					Sum+= BlackSholes_Formula(f1,sigma1total,1.0,K1,callput)*c.get_weight(i);
				}
				return Sum*exp(-r*T);
			}
			else
			{	
				Voltotal=sigmatotal;
				f1=f*T/averaging*exp(0.5*(r-Voltotal*Voltotal/(6.0*T))*T);
				sigma1total=sqrt(exp(Voltotal*Voltotal/3.0) - 1.0);
				return  BlackSholes_Formula(f1,sigma1total,1.0,K1,callput)*exp(-r*T);
			}
		}
		else
		{
			return StochasticVol_LN_Geometric_VanillaOption(f,K,T,r,sig,VolDrift,VolOfVol,averaging,callput,LegendreNb);
		}
	}
	
}



///////////////////////////////////////////////////////////////////////////////////
///
///
///					Alternative Formulas
///
///
////////////////////////////////////////////////////////////////////////////////////
double AsianCallLevy(double S,double X, double r,double b,double T,double T2,double sig, double SA,int callput)
{
	double SE=S/(T*b)*(exp((b-r)*T2)-exp(-r*T2));
	double Xs=X-(T-T2)/T*SA;
	double M=2.0*S*S/(sig*sig+b)*((exp((2.0*b+sig*sig)*T2)-1.0)/(sig*sig+2.0*b)-(exp(b*T2)-1.0)/b);
	double Dd=M/(T*T);
	double V=log(Dd)-2.0*(r*T2+log(SE));
	double d1=1.0/sqrt(V)*(log(Dd)/2.0-log(Xs));
	double d2=d1-sqrt(V);
	if(callput==K_CALL)
	{
		return SE*ARM_GaussianAnalytics::cdfNormal(d1)-Xs*exp(-r*T2)*ARM_GaussianAnalytics::cdfNormal(d2);
	}
	else
	{
		return SE*ARM_GaussianAnalytics::cdfNormal(d1)-Xs*exp(-r*T2)*ARM_GaussianAnalytics::cdfNormal(d2)-SE+Xs*exp(-r*T2);
	}
	
}

double AsianCurranApprox(int callput, double S, double SA, double X, double t1,
						 double T, double n, double m, double r, double b, double v)
{	
	// By Espen Gaarder Haug, The Collector, 2000, based on Curran's method
	// See The Complete Guide To Option Pricing Formulas for more Asian option formulas
	double dt, my, myi;
	double vxi, vi, vx;
	double Km, sum1, sum2;
	double ti, EA;
	int z, i;
	z=callput;
	if (m==n-1.0) // Only one fix left use black Scholes weighted with time
	{
		X=n*X-(n-1.0)*SA;
		return exp(-r*T)*BlackSholes_Formula(S*exp((r-b)*T),v*sqrt(T),1.0,X,callput) * 1.0 / n;
	};
	dt=(T-t1)/(n-1.0);
	
	if (b==0.0)
	{
		EA=S;
	}
	else
	{
		EA=S/n*exp(b*t1)*(1.0-exp(b*dt*n))/(1.0-exp(b*dt));
	};
	if (m>0.0)
	{
		if (SA>(n/m*X) ) //Exercise is certain for call, put must be out-of-the-money
		{
			if (callput==K_CALL)
			{
				return 0.0;
			}
			else if (callput==K_PUT)
			{
				SA=SA*m/n+EA*(n-m)/n;
				return (SA-X)*exp(-r*T);
			};};};
	if (m>0.0) X=n/(n-m)*X-m/(n-m)*SA;
	vx = v*sqrt(t1+dt*(n-1.0)*(2*n-1.0)/(6.0*n));
	my =log(S)+(b-v*v*0.5)*(t1+(n-1.0)*dt/2.0);
	sum1 = 0.0;
	for (i=1;i<=n;i++)
	{
		ti=dt*i+t1-dt;
		vi=v*sqrt(t1+(i-1.0)*dt);
		vxi=v*v*(t1+dt*((i-1.0)-i*(i-1.0)/(2.0*n)));
		myi=log(S)+(b-v*v*0.5)*ti;
		sum1=sum1+exp(myi+vxi/(vx*vx)*
			(log(X)-my)+(vi*vi-vxi*vxi/(vx*vx))*0.5);
	};
	Km=2.0*X-1.0/n*sum1;
	sum2=0.0;
	for (i=1;i<=n;i++)
	{
		ti=dt*i+t1-dt;
		vi=v*sqrt(t1+(i-1.0)*dt);
		vxi=v*v*(t1+dt*((i-1.0)-i*(i-1.0)/(2.0*n)));
		myi=log(S)+(b-v*v*0.5)*ti;
		sum2=sum2+exp(myi+vi*vi*0.5)*ARM_GaussianAnalytics::cdfNormal(z*((my-log(Km))/vx+vxi/vx));
	};
	return exp(-r*T)*z*(1.0/n*sum2-X*ARM_GaussianAnalytics::cdfNormal(z*(my-log(Km))/vx))*(n-m)/n;
};


CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
