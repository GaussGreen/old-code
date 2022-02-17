/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file SABR_Analytics.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/cev.h"
#include "gpclosedforms/bessel.h"
#include "gpclosedforms/hypergeometric.h"
#include "gpclosedforms/whittaker.h"
#include "gpclosedforms/barriere_bs.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/gamma.h"

#include "gpbase/numericconstant.h"




CC_BEGIN_NAMESPACE(ARM)



//////////////////////////////////////////////////////////////////////////////////////////////
///
///			In the following, the process is assumed to be
///
///				dS(t)=r*S(t) dt+delta*S(t)^(1-1/(2*nu))
///					and the delta input is comparable to  : 
///				sigmaLN = delta*S(0)^(-1/(2*nu))
///

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///			CEV Vanilla option through the integration of bessel functions (not the inifinite sum of gammas)
///
///			
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////


double CEV_Integral1_SmallZ(double x, double a, double nu, double c,int nbsteps)
{
	GaussLaguerre_Coefficients p(0,nbsteps);
	double Sum=0;
	double z;
	for(int i=0;i<nbsteps;i++){
		z=p.get_point(i)+a;
		Sum+=bessel_fractional_I(nu,z)*exp(z-a-z*z/(4*x)-x)*z/(2*x) *pow(z/(2*c),nu)*p.get_weight(i);		// Exp(z-a) supplementaire a cause de Laguerre
	}
	return Sum;
}


double CEV_Integral2_SmallZ(double x, double a, double nu, int nbsteps)
{
	GaussLaguerre_Coefficients p(0,nbsteps);
	double Sum=0;
	double z;
	for(int i=0;i<nbsteps;i++){
		z=p.get_point(i)+a;
		Sum+=bessel_fractional_I(nu,z)*exp(z-a-z*z/(4*x)-x) *pow(2.*x/z,(nu-1.))*p.get_weight(i);		// Exp(z-a) supplementaire a cause de Laguerre
	}
	return Sum;
}


double CEV_Integral1_LargeZ(double x, double a, double nu, double c,int nbsteps)
{
	GaussLaguerre_Coefficients p(0,nbsteps);
	double Sum=0;
	double z;
	for(int i=0;i<nbsteps;i++){
		z=p.get_point(i)+a;
		Sum+=exp(bessel_fractional_I_LargeZ_Log(nu,z,20)+z-a-z*z/(4*x)-x+log(z/(2*x)) +log(z/(2*c))*nu)*p.get_weight(i);		// Exp(z-a) supplementaire a cause de Laguerre
	}
	return Sum;
}

double CEV_Integral2_LargeZ(double x, double a, double nu, int nbsteps)
{
	GaussLaguerre_Coefficients p(0,nbsteps);
	double Sum=0;
	double z;
	for(int i=0;i<nbsteps;i++){
		z=p.get_point(i)+a;
		Sum+=exp(bessel_fractional_I_LargeZ_Log(nu,z,20)+z-a-z*z/(4*x)-x+log(2.*x/z)*(nu-1.))*p.get_weight(i);		// Exp(z-a) supplementaire a cause de Laguerre
	}
	return Sum;
}

double CEV_Integral1(double x, double a, double nu, double c,int nbsteps)
{
	if(a<100.) return CEV_Integral1_SmallZ(x,a,nu,c,nbsteps);
	else return CEV_Integral1_LargeZ(x,a,nu,c,nbsteps);
	
}

double CEV_Integral2(double x, double a, double nu, int nbsteps)
{
	if(a<100.) return CEV_Integral2_SmallZ(x,a,nu,nbsteps);
	else return CEV_Integral2_LargeZ(x,a,nu,nbsteps);
}




double CEV_Call(double f, double K, double T,double drift, double vol, double nu,int nbsteps)
{
	double c,x,sig;
	sig=vol*pow(f,1./(2.*nu));
	c=2.*nu*drift/(sig*sig*(exp(drift*T/nu)-1));
	x=pow(f,1./nu)*c*exp(drift*T/nu);
	double a=2.*sqrt(c*x*pow(K,1./nu));
	return exp(-drift)*(CEV_Integral1(x,a,nu,c,nbsteps)-K*CEV_Integral2(x,a,nu,nbsteps));
}
//////////////////////////////////////////////////////////////////////////////////////////////
///
///			In the following, the process is assumed to be
///
///				dS(t)=mu*S(t) dt+delta*S(t)^(beta+1)
///
///
////////////////////////////////////////////////////////////////////////////////////////////////
///
///					1)	Double barrier call
///			
////////////////////////////////////////////////////////////////////////////////////////////////


double CEV_DoubleBarrierCall(double S,double K,double T,double r,double beta,double mu,double delta,double L,double U,int  n)
{ 
	if(K>=U) return 0;
	else
	{
		int i;
		long double kn, Nn,cn,ln,fn, sum=0.;
		for(i=1;i<=n;i++)
		{
			kn = CEV_DoubleBarrierCallComputekn(beta, mu, delta, L,U, i);
			Nn = CEV_DoubleBarrierCallNormalisation(K, beta, mu, delta, L,U, kn);
			cn=CEV_DoubleBarrierCallComputecn(K, beta, mu, delta,L, U, kn, Nn);
			ln=CEV_DoubleBarrierCallComputeEigenValue(beta, mu, delta,L, U, kn);
			fn=CEV_DoubleBarrierCallComputeEigenFunction(S, beta, mu, delta,L, U, kn, Nn);
			sum+= cn*exp(-(r + ln)*T)*fn;
		}
		return sum;
	}
}



///
///				 the roots of the wronskian are the eigen values of the elliptic operator
///				 

/// wronskian : deltakm function

double CEV_DoubleBarrierCalldeltakm2(double k,double m,double a,double b)
{
	/*
	
	cout<<"\n\n k="<<k<<" m="<<m<<" b="<<b<<" wM ="<<Hypergeometric_Whittaker_M(k,m,b);
	cout<< "\n\n inside deltakn ************";
	cout<<"\n\n w1="<<Hypergeometric_Whittaker_W(k,m,a);
	cout<<"\n\n m1="<<Hypergeometric_Whittaker_M(k,m,b);
	cout<<"\n\n w2="<<Hypergeometric_Whittaker_M(k,m,a);
	cout<<"\n\n m2="<<Hypergeometric_Whittaker_W(k,m,b);
	cout<<"\n gam1="<<gammaLD(0.5+m-k).todouble();
	cout<<"\n gam2="<<gammaLD(1.+2.*m).todouble();
	*/
	
	double resu=Hypergeometric_Whittaker_W(k,m,a)*Hypergeometric_Whittaker_M(k,m,b)-
			Hypergeometric_Whittaker_M(k,m,a)*Hypergeometric_Whittaker_W(k,m,b);
	long_double r=gammaLD(0.5+m-k)/gammaLD(1.+2.*m);
	return r.todouble()*resu;
}


double CEV_DoubleBarrierCalldeltakm(double k,double m,double a,double b)
{	

//	cout<<"\n\n k="<<k<<" m="<<m<<"  a="<<a<<" b="<<b;

	long_double resu1=Hypergeometric_Whittaker_W_L(k,m,a)*Hypergeometric_Whittaker_M_L(k,m,b);
	long_double resu2=Hypergeometric_Whittaker_M_L(k,m,a)*Hypergeometric_Whittaker_W_L(k,m,b);
	long_double factor=gammaLD(0.5+m-k)/gammaLD(1.+2.*m);
	long_double resu1a=resu1*factor;
	long_double resu2a=resu2*factor;
	return resu1a.todouble()-resu2a.todouble();	
}

double CEV_DoubleBarrierCallComputekn(double beta,double mu, double delta, double L, double U, int n)
{
	double b=fabs(beta);
	double B=sqrt(2.)/(delta*b)*(pow(U,b)-pow(L,b));
	double a0=-(pow(delta,-2)*pow(L,-b)*pow(U,-b)*
      (3*(-1 + pow(b,2))*pow(delta,4) + 
        12*(1 + 2*beta)*mu*pow(delta,2)*pow(L,b)*
         pow(U,b) - 4*pow(L,b)*pow(mu,2)*pow(U,b)*
         (pow(L,2*b) + pow(L,b)*pow(U,b) + pow(U,2*b)))
      )/24.;
	double a2=(pow(b,-6)*pow(beta,-4)*pow(delta,-6)*pow(L,-3*b)*
     pow(U,-3*b)*pow(pow(L,b) - pow(U,b),2)*
     (360*mu*pow(b,8)*(-1 - 2*beta + pow(b,2) + 
          2*pow(beta,3))*pow(delta,6)*pow(L,2*b)*
        pow(U,2*b) + 15*
        (1 - 26*pow(b,2) + 25*pow(b,4))*pow(b,8)*
        pow(delta,8)*(pow(L,2*b) + pow(L,b)*pow(U,b) + 
          pow(U,2*b)) - 
       480*(1 + 2*beta)*pow(b,8)*pow(delta,2)*
        pow(L,3*b)*pow(mu,3)*
        (pow(L,2*b) + pow(L,b)*pow(U,b) + pow(U,2*b))*
        pow(U,3*b) + 144*pow(b,8)*pow(L,3*b)*pow(mu,4)*
        pow(U,3*b)*(pow(L,4*b) + pow(L,3*b)*pow(U,b) + 
          pow(L,2*b)*pow(U,2*b) + 
          pow(L,b)*pow(U,3*b) + pow(U,4*b)) - 
       5*pow(b,8)*pow(L,b)*pow(U,b)*
        (72*mu*(-1 - 2*beta + pow(b,2) + 
             2*pow(beta,3))*pow(delta,6)*pow(L,b)*
           pow(U,b) - 96*(1 + 2*beta)*pow(delta,2)*
           pow(L,2*b)*pow(mu,3)*pow(U,2*b)*
           (pow(L,2*b) + pow(L,b)*pow(U,b) + 
             pow(U,2*b)) - 
          24*pow(delta,4)*pow(L,b)*pow(mu,2)*pow(U,b)*
           ((-1 + pow(b,2))*pow(L,2*b) - 
             2*(-1 + 7*pow(b,2))*pow(L,b)*pow(U,b) + 
             (-1 + pow(b,2))*pow(U,2*b)) + 
          9*pow(delta,8)*pow(-1 + pow(b,2),2) + 
          16*pow(L,2*b)*pow(mu,4)*pow(U,2*b)*
           pow(pow(L,2*b) + pow(L,b)*pow(U,b) + 
             pow(U,2*b),2))))/5760.;

// FIXMEFRED: mig.vc8 (22/05/2007 18:24:39):pow(int, ..) doesnt exist
	return (-2 + pow(fabs(beta),-1) + 
     2*pow(mu,-1)*(a0 + 
        pow(static_cast<double>(n),-4)*pow(ARM_NumericConstants::ARM_PI,-4)*
		( a2*pow(static_cast<double>(n),2)*pow(ARM_NumericConstants::ARM_PI,2) + 
		pow(B,-2)*pow(static_cast<double>(n),6)*pow(ARM_NumericConstants::ARM_PI,6)))*
      pow(fabs(beta),-1))/4.;
}

double CEV_DoubleBarrierCallComputeEigenValue (double beta,double mu, double delta, double L, double U, double kn)
{
	return 2.*mu*fabs(beta)*(0.5 + kn - 1./(fabs(beta)*4.));
}


double CEV_DoubleBarrierCallComputeEigenFunction (double S, double beta,double mu, double delta, double L, double U, double kn, double Nn)
{
	double m,l,u,xS;
	m=1/(4*fabs(beta));
	l=mu/(delta*delta*fabs(beta))*pow(L,(-2.*beta));
	u=mu/(delta*delta*fabs(beta))*pow(U,(-2.*beta));
	xS=mu/(delta*delta*fabs(beta))*pow(S,(-2.*beta));
	double f=CEV_DoubleBarrierCalldeltakm(kn,m,l,xS);
	double g=Nn*pow(S,0.5+beta)*exp(-xS/2.);
	return f*g;	
}

double CEV_DoubleBarrierCallNormalisation (double K, double beta,double mu, double delta, double L, double U, double kn)
{
	double m,l,u,kappa,Dnm;
	m=1/(4*fabs(beta));
	l=mu/(delta*delta*fabs(beta))*pow(L,(-2.*beta));
	u=mu/(delta*delta*fabs(beta))*pow(U,(-2.*beta));
	kappa=mu/(delta*delta*fabs(beta))*pow(K,(-2.*beta));
	double Dnm1=CEV_DoubleBarrierCalldeltakm(kn+0.001,m,l,u);
	double Dnm2=CEV_DoubleBarrierCalldeltakm(kn-0.001,m,l,u);
	Dnm=(Dnm1-Dnm2)/0.002;
	long_double resu=Hypergeometric_Whittaker_W_L(kn,m,u)/Hypergeometric_Whittaker_W_L(kn,m,l);
	return sqrt(delta*delta*fabs(beta)*resu.todouble()/Dnm);	
}


double CEV_DoubleBarrierCallComputecn (double K, double beta,double mu, double delta, double L, double U, double kn, double Nn)
{
	double m,l,u,kappa,In,Jn;
	m=1/(4*fabs(beta));
	l=mu/(delta*delta*fabs(beta))*pow(L,(-2.*beta));
	u=mu/(delta*delta*fabs(beta))*pow(U,(-2.*beta));
	kappa=mu/(delta*delta*fabs(beta))*pow(K,(-2.*beta));
	long_double factor_global=gammaLD(0.5-kn+m)/gammaLD(1.+2.*m);
	long_double factor_Ina=Hypergeometric_Whittaker_W_L(kn,m,l);
	long_double factor_Jna=Hypergeometric_Whittaker_M_L(kn,m,l);
	long_double factor_In=factor_global*factor_Ina;
	long_double factor_Jn=factor_global*factor_Jna;

	long_double In_sub1=Hypergeometric_Whittaker_M_L(0.5 + kn,-0.5 + m,kappa)*factor_In;
	long_double In_sub2=Hypergeometric_Whittaker_M_L(0.5 + kn,-0.5 + m,u)*factor_In;
	long_double In_sub3=Hypergeometric_Whittaker_M_L(0.5 + kn,0.5 + m,kappa)*factor_In;
	long_double In_sub4=Hypergeometric_Whittaker_M_L(0.5 + kn,0.5 + m,u)*factor_In;

	long_double Jn_sub1=Hypergeometric_Whittaker_W_L(0.5 + kn,-0.5 + m,kappa)*factor_Jn;
	long_double Jn_sub2=Hypergeometric_Whittaker_W_L(0.5 + kn,-0.5 + m,u)*factor_Jn;
	long_double Jn_sub3=Hypergeometric_Whittaker_W_L(0.5 + kn,0.5 + m,kappa)*factor_Jn;
	long_double Jn_sub4=Hypergeometric_Whittaker_W_L(0.5 + kn,0.5 + m,u)*factor_Jn;

	In=1./(delta*sqrt(mu*fabs(beta)))*
		(2*m*exp(kappa/2.)*sqrt(K)/(-0.5 - kn + m)*
		In_sub1.todouble()- 
		2*K*m*exp(u/2.)/(-0.5 - kn + m)/sqrt(U)*
		In_sub2.todouble()- 
		exp(kappa/2.)*sqrt(K)*pow(1 + 2*m,-1)*
		In_sub3.todouble()+ 
		exp(u/2.)/(1. + 2*m)*sqrt(U)*
		In_sub4.todouble());
	Jn=1./(delta*sqrt(mu*fabs(beta)))*
		(exp(kappa/2.)*sqrt(K)/(0.5 + kn - m)*
		Jn_sub1.todouble()- 
		K*exp(u/2.)/(0.5 + kn - m)/sqrt(U)*
		Jn_sub2.todouble()- 
		exp(kappa/2.)*sqrt(K)/(0.5 + kn + m)*
		Jn_sub3.todouble()+ 
		exp(u/2.)/(0.5 + kn + m)*sqrt(U)*
		Jn_sub4.todouble());
	return Nn*(In-Jn);	
}

////////////////////////////////////////////////////////////////////////////////////////////////
///
///				2)		Single barrier Up and Out Call
///			
////////////////////////////////////////////////////////////////////////////////////////////////

double CEV_SingleBarrierUpandOutCall(double S,double K,double T,double r,double beta,double mu,double delta,double U,int  n)
{ 
	if(K>=U) return 0;
	else
	{
		int i;
		double kn, Nn,cn,ln,fn, sum=0.;
		for(i=1;i<=n;i++)
		{
			kn = CEV_SingleBarrierUpandOutCallComputekn(beta, mu, delta, U, i);
			Nn = CEV_SingleBarrierUpandOutCallNormalisation(K, beta, mu, delta, U, kn);
			cn=CEV_SingleBarrierUpandOutCallComputecn(K, beta, mu, delta, U, kn, Nn);
			ln=CEV_SingleBarrierUpandOutCallComputeEigenValue(beta, mu, delta, U, kn);
			fn=CEV_SingleBarrierUpandOutCallComputeEigenFunction(S, beta, mu, delta, U, kn, Nn);
			sum+= cn*exp(-(r + ln)*T)*fn;
		}
		return sum;
	}
}


double CEV_SingleBarrierUpandOutCallComputekn(double beta,double mu,double delta,double U,double n)
{
	double m=1./(4.*fabs(beta));
	double u=mu/(delta*delta*fabs(beta))*pow(U,-2.*beta);
	double kn0=(n+m-0.25)*(n+m-0.25)*ARM_NumericConstants::ARM_PI*ARM_NumericConstants::ARM_PI/(4.*u);
	class whittakerM: public DoubleToDoubleFunc
	{
	public:
		double m0;
		double u0;
		whittakerM(double m1, double u1): m0(m1), u0(u1) {}
		virtual double operator() (double knx0) const
		{
			return Hypergeometric_Whittaker_M(knx0,m0,u0);
		}
	};
	whittakerM w(m,u);
	double knx=Inverse(w,Inverse::ALWAYSPOSITIVE)(0,kn0,kn0/20.,1e-12);
	return knx;
}

double CEV_SingleBarrierUpandOutCallNormalisation(double K,double beta,double mu,double delta,double U,double kn)
{
	double m=1./(4.*fabs(beta));
	double u=mu/(delta*delta*fabs(beta))*pow(U,-2.*beta);
	double Mnm=(Hypergeometric_Whittaker_M(kn+0.001,m,u)-Hypergeometric_Whittaker_M(kn-0.001,m,u))/0.002;
	long_double Unm=(Hypergeometric_Whittaker_W_L(kn,m,u)*gammaLD(0.5-kn+m))/gammaLD(1.+2.*m);
	return sqrt(fabs(beta)*delta*delta*Unm.todouble()/Mnm);
}

double CEV_SingleBarrierUpandOutCallComputecn(double K,double beta,double mu,double delta,double U,double kn,double Nn)
{
	double m,u,kappa,In;
	m=1/(4*fabs(beta));
	u=mu/(delta*delta*fabs(beta))*pow(U,(-2.*beta));
	kappa=mu/(delta*delta*fabs(beta))*pow(K,(-2.*beta));

	long_double In_sub1=Hypergeometric_Whittaker_M_L(0.5 + kn,-0.5 + m,kappa);
	long_double In_sub2=Hypergeometric_Whittaker_M_L(0.5 + kn,-0.5 + m,u);
	long_double In_sub3=Hypergeometric_Whittaker_M_L(0.5 + kn,0.5 + m,kappa);
	long_double In_sub4=Hypergeometric_Whittaker_M_L(0.5 + kn,0.5 + m,u);

	In=1./(delta*sqrt(mu*fabs(beta)))*
		(2*m*exp(kappa/2.)*sqrt(K)/(-0.5 - kn + m)*
		In_sub1.todouble()- 
		2*K*m*exp(u/2.)/(-0.5 - kn + m)/sqrt(U)*
		In_sub2.todouble()- 
		exp(kappa/2.)*sqrt(K)/(1.+ 2*m)*
		In_sub3.todouble()+ 
		exp(u/2.)/(1. + 2*m)*sqrt(U)*
		In_sub4.todouble());
	return Nn*In;	
}

double CEV_SingleBarrierUpandOutCallComputeEigenValue(double beta,double mu,double delta,double U,double kn)
{
	return 2.*mu*fabs(beta)*(0.5 + kn - 1./(fabs(beta)*4.));
}

double CEV_SingleBarrierUpandOutCallComputeEigenFunction(double S,double beta,double mu,double delta,double U,double kn,double Nn)
{
	double m,u,xS;
	m=1/(4*fabs(beta));
	u=mu/(delta*delta*fabs(beta))*pow(U,(-2.*beta));
	xS=mu/(delta*delta*fabs(beta))*pow(S,(-2.*beta));
	double f=Hypergeometric_Whittaker_M(kn,m,xS);
	double g=Nn*pow(S,0.5+beta)*exp(-xS/2.);
	return f*g;	
}


////////////////////////////////////////////////////////////////////////////////////////////////
///
///				3)		Single barrier Down and Out Put
///			
////////////////////////////////////////////////////////////////////////////////////////////////

double CEV_SingleBarrierDownandOutPut(double S,double K,double T,double r,double beta,double mu,double delta,double L,int  n)
{ 
	if(K<=L) return 0;
	else
	{
		int i;
		double kn, Nn,cn,ln,fn, sum=0.;
		for(i=1;i<=n;i++)
		{
			kn = CEV_SingleBarrierDownandOutPutComputekn(beta, mu, delta, L, i);
			Nn = CEV_SingleBarrierDownandOutPutNormalisation(K, beta, mu, delta, L, kn);
			cn=CEV_SingleBarrierDownandOutPutComputecn(K, beta, mu, delta, L, kn, Nn);
			ln=CEV_SingleBarrierDownandOutPutComputeEigenValue(beta, mu, delta, L, kn);
			fn=CEV_SingleBarrierDownandOutPutComputeEigenFunction(S, beta, mu, delta, L, kn, Nn);
			sum+= cn*exp(-(r + ln)*T)*fn;
		}
		return sum;
	}
}


double CEV_SingleBarrierDownandOutPutComputekn(double beta,double mu,double delta,double L,double n)
{
	double m=1./(4.*fabs(beta));
	double l=mu/(delta*delta*fabs(beta))*pow(L,-2.*beta);
	double pi=ARM_NumericConstants::ARM_PI;
	double kn0=n-0.25+2.*l/(pi*pi)+2./pi*sqrt((n-0.25)*l+l*l/(pi*pi));
	class whittakerM: public DoubleToDoubleFunc
	{
	public:
		double m0;
		double l0;
		whittakerM(double m1, double l1): m0(m1), l0(l1) {}
		virtual double operator() (double knx0) const
		{
			return Hypergeometric_Whittaker_W(knx0,m0,l0);
		}
	};
	whittakerM w(m,l);
	double knx=Inverse(w,Inverse::ALWAYSPOSITIVE)(0,kn0,kn0/20.,1e-12);
	return knx;
}

double CEV_SingleBarrierDownandOutPutNormalisation(double K,double beta,double mu,double delta,double L,double kn)
{
	double m=1./(4.*fabs(beta));
	double l=mu/(delta*delta*fabs(beta))*pow(L,-2.*beta);
	double Mnm=(Hypergeometric_Whittaker_W(kn+0.001,m,l)-Hypergeometric_Whittaker_W(kn-0.001,m,l))/0.002;
	long_double Unm=(Hypergeometric_Whittaker_M_L(kn,m,l)*gammaLD(0.5-kn+m))/gammaLD(1.+2.*m);
	if ((log(Unm.todouble())-log(Mnm))<-700.)
	{
		return 0.;
	}
	else
	{
	return sqrt(fabs(beta)*delta*delta*Unm.todouble()/Mnm);
	}
}

double CEV_SingleBarrierDownandOutPutComputecn(double K,double beta,double mu,double delta,double L,double kn,double Nn)
{
	double m,l,kappa,In;
	m=1/(4*fabs(beta));
	l=mu/(delta*delta*fabs(beta))*pow(L,(-2.*beta));
	kappa=mu/(delta*delta*fabs(beta))*pow(K,(-2.*beta));

	long_double In_sub1=Hypergeometric_Whittaker_W_L(0.5 + kn,0.5 + m,l);
	long_double In_sub2=Hypergeometric_Whittaker_W_L(0.5 + kn,-0.5 + m,l);
	long_double In_sub3=Hypergeometric_Whittaker_W_L(0.5 + kn,0.5 + m,kappa);
	long_double In_sub4=Hypergeometric_Whittaker_W_L(0.5 + kn,-0.5 + m,kappa);

	In=1./(delta*sqrt(mu*fabs(beta)))*
		(
		exp(l/2.)*sqrt(L)/(0.5 + kn + m)*In_sub1.todouble()
		-K/sqrt(L)*exp(l/2.)/(0.5 + kn - m)*In_sub2.todouble()
		-sqrt(K)*exp(kappa/2.)/(0.5 + kn + m)*In_sub3.todouble()
		+sqrt(K)*exp(kappa/2.)/(0.5 + kn - m)*In_sub4.todouble()
		);
	return Nn*In;	
}

double CEV_SingleBarrierDownandOutPutComputeEigenValue(double beta,double mu,double delta,double L,double kn)
{
	return 2.*mu*fabs(beta)*(0.5 + kn - 1./(fabs(beta)*4.));
}

double CEV_SingleBarrierDownandOutPutComputeEigenFunction(double S,double beta,double mu,double delta,double L,double kn,double Nn)
{
	double m,l,xS;
	m=1/(4*fabs(beta));
	l=mu/(delta*delta*fabs(beta))*pow(L,(-2.*beta));
	xS=mu/(delta*delta*fabs(beta))*pow(S,(-2.*beta));
	double f=Hypergeometric_Whittaker_W(kn,m,xS);
	double g=Nn*pow(S,0.5+beta)*exp(-xS/2.);
	return f*g;	
}


CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/