/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_spreadoption.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2006
 */

#include <glob/firsttoinc.h>

#include <cmath>
#include <complex>

#include <vector>
#include <iostream>
#include <iomanip>
#include <glob/expt.h>
#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/erf.h"
#include "gpclosedforms/hypergeometric.h"
#include "gpclosedforms/bisabr_spreadoption.h"
#include "gpclosedforms/spreadoption_lognormal_formula.h"
#include "gpclosedforms/spreadoption_lognormal.h"
#include "gpclosedforms/spreadoption_lognormal_interface.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/normal_heston.h"

using namespace std; 


CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_K_SHIFT_FOR_DERIVATION 0.0000001


////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of sum{0,x} of exp(-a^2*s^2-b^2/s^2
///
////////////////////////////////////////////////////////////////////////////////////////////

double GIntegralAnalytical(double a1, double b1,double x1)
{
	double x=fabs(x1);
	double b=fabs(b1);
	double a=fabs(a1);
	double signx1;
	if(x1>0) signx1=1.0; else signx1=-1.0;
	if(a>0.0001)
	{
		return signx1*(exp(-2.*a*b)*ARM_NumericConstants::ARM_SQRT_PI*
			(1.-exp(4.*a*b)-NormalCDF2(ARM_NumericConstants::ARM_SQRT_2*(b/x-a*x))+
			exp(4.*a*b)*NormalCDF2(ARM_NumericConstants::ARM_SQRT_2*(b/x+a*x))))/(2.*a);
	}
	else
	{
		return signx1*(exp(-b*b/(x*x))*x*(30.+10*a*a*(2.*b*b-x*x)+a*a*a*a*(4.*b*b*b*b-
			2.*b*b*x*x+3.*x*x*x*x))+4.*b*(15.+10.*a*a*b*b+2.*a*a*a*a*b*b*b*b)*
			ARM_NumericConstants::ARM_SQRT_PI*(-1.+NormalCDF2(ARM_NumericConstants::ARM_SQRT_2*b/x)))/30.;
	}
}

///  GIntegralAnalytical2(double a2, double b1,double x1) == GIntegralAnalytical(double sqrt(a2), double b1,double x1)
/// also works when a2 <0 ! it is a taylor expansion valid for abs(a2)<1
////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of sum{0,x} of exp(-a2*s^2-b1^2/s^2
///
////////////////////////////////////////////////////////////////////////////////////////////

double GIntegralAnalytical2(double a2, double b1,double x1)
{
	/*
	if(a2>0)
	{
		double a=sqrt(a2);
		double b=fabs(b1);
		if(((b/x1+a*x1)>3.)&&((b/x1-a*x1)>3.))
		{
			
			double x2=x1*x1;
			double x4=x2*x2;
			double A=b-a*x2;
			double A2=A*A;
			double A4=A2*A2;
			double B=b+a*x2;
			double B2=B*B;
			double B4=B2*B2;
			double res= (16.+105.*x4*x4/(A4*A4)-30.*x2*x4/(A2*A4)+12.*x4/A4-8.*x2/A2)/A+
				(-105.*x4*x4+30.*x2*x4*B2-12.*x4*B4+8.*x2*B2*B4-16.*B4*B4)/(B4*B4*B);
			return exp(-a2*x2-b1*b1/x2)*x1/16.*res;		
		}
		else
		{
			return GIntegralAnalytical( sqrt(a2),  b1, x1);
		}
	}
	else
	*/
	{
		double x=x1;
		double b=fabs(b1);
		double ERF=2.*NormalCDF2(ARM_NumericConstants::ARM_SQRT_2*b/x)-2.;
		double expb2=exp(b*b/(x*x));
		double b2=b*b;
		double b3=b2*b;
		double b5=b3*b2;
		double x2=x*x;
		double x3=x2*x;
		double x5=x3*x*x;
		double a4=a2*a2;
		double a0=ERF*b*ARM_NumericConstants::ARM_SQRT_PI+x/expb2;
		double b0=2.*b3*ERF*ARM_NumericConstants::ARM_SQRT_PI+(2.*b2*x-x3)/expb2;
		double c0=4.*b3*b2*ERF*ARM_NumericConstants::ARM_SQRT_PI+(4.*b3*b*x-2.*b2*x3+3.*x3*x2)/expb2;
		double d0=8.*b5*b2*ERF*ARM_NumericConstants::ARM_SQRT_PI+(8.*b3*b3*x-4.*b3*b*x3+6.*b2*x3*x2-15.*x5*x2)/expb2;
		double e0=16.*b5*b2*b2*ERF*ARM_NumericConstants::ARM_SQRT_PI+(16.*b5*b3*x-8.*b5*b*x3+12.*b3*b*x5-30.*b2*x5*x2+105.*x5*x3*x)/expb2;
		return a0+b0*a2/3.+c0*a4/30.+d0*a4*a2/630.+e0*a4*a4/22680.;
	}
}





double BiSABRIntegral(double A,double B,double G,double F,double gamma)
{
	double xpoint,xweight;
	int nbl=120;
	int i;
	GaussLegendre_Coefficients d(nbl, A, B);
	double Sum=0;
	for(i=0;i<nbl;i++)
	{
		xpoint=d.get_point(i);
		xweight=d.get_weight(i);
		Sum+=pow(xpoint,gamma)/sqrt(xpoint*xpoint-xpoint*G+F)*xweight;
	}
	return Sum;
}

double BiSABRIntegral_Analytical(double A,double B,double G,double F,double gamma)
{	
	int Nb=150;
	complex<double> Ac(A,0);
	complex<double> Bc(B,0);
	complex<double> Gc(G,0);
	complex<double> Fc(F,0);
	complex<double> gammac(gamma,0);
	complex<double> un(1.,0);
	complex<double> undemi(0.5,0);
	complex<double> deux(2.,0);
	complex<double> quatre(4.,0);
	
	complex<double> sq=sqrt(Gc*Gc-quatre*Fc);
	
	double x1=real(HypergeometricAppellF1(un+gammac,undemi,undemi,deux+gammac,deux*Ac/(G-sq),deux*Ac/(G+sq),Nb));
	double x2=real(HypergeometricAppellF1(un+gammac,undemi,undemi,deux+gammac,deux*Bc/(G-sq),deux*Bc/(G+sq),Nb));
	
	return (pow(B,1.+gamma)*x2-pow(A,1.+gamma)*x1)/(sqrt(F)*(1.+gamma));
}



double BiSABRz(double K,double F1,double F2,double alpha1,double alpha2,double rhos,double beta1,double beta2)
{
	double F2beta2=pow(F2,beta2);
	double alpha2overalpha1=alpha2/alpha1;
	if(K+F2>0)
	{
	return BiSABRIntegral(pow(K+F2,beta1),pow(F1,beta1),2.*F2beta2*alpha2overalpha1*rhos,
		F2beta2*F2beta2*alpha2overalpha1*alpha2overalpha1,(1.-beta1)/beta1)/(alpha1*beta1);
	}
	else
	{
	return BiSABRIntegral(0,pow(F1,beta1),2.*F2beta2*alpha2overalpha1*rhos,
		F2beta2*F2beta2*alpha2overalpha1*alpha2overalpha1,(1.-beta1)/beta1)/(alpha1*beta1);

	}
}


double SABRgeneric5(double intrinsec,double z,double b1,double b2,double zavg,double sqB0Balphaz,double T,double mu,double alpha,double rho,double nu)
{
	double x,I0,I1,I2,kappa,f1overf0;
	double Iz2=1.+z*nu*(z*nu-2.*rho);
	double sqrho=sqrt(1.-rho*rho);
	if(nu>0.00000001)
	{
		x=log((-rho+nu*z+sqrt(Iz2))/(1.-rho))/nu;
		f1overf0=exp(((2.*mu+b1*alpha*nu*rho)*(2.*rho*atan((z*nu-rho)/sqrho)+2.*rho*atan(rho/sqrho)+sqrho*log(Iz2)))/(4.*nu*nu*sqrho));
		
	}
	else
	{
		x=z;
		f1overf0=exp(mu*z*z/2.);
		
	}
	I0=sqrt(1.+zavg*nu*(zavg*nu-2.*rho));
	I1=(zavg*nu-rho)/I0;
	I2=(1.-rho*rho)/(I0*I0*I0);
	kappa=(12.*mu+(2*b2-3.*b1*b1)*alpha*alpha
		-4.*z*z*mu*mu+(2.*I0*I2-I1*I1)*nu*nu
		+(6.*b1*alpha-16.*zavg*mu)*nu*rho)
		/(8.*(1.0+zavg*nu*(z*nu-2.*rho)));
	double gintegral=GIntegralAnalytical2(-kappa,x/ARM_NumericConstants::ARM_SQRT_2,sqrt(T));
	return intrinsec+alpha*ARM_NumericConstants::ARM_INVSQRT2PI*sqB0Balphaz*f1overf0*sqrt(I0)*gintegral;
	
}


///  Modification of the kappa that garantie a long term behaviour more stable and more compatible
///  with the Hagan formula
/// by adding a factor (3.*I0*I0-1.) instead of 2 inside the kappa, making the kappa match the hagan one for vanilla sabr
double SABRgeneric6(double intrinsec,double z,double b1,double b2,double zavg,double sqB0Balphaz,double T,double mu,double alpha,double rho,double nu)
{
	double x,I0,I1,I2,kappa,f1overf0;
	double Iz2=1.+z*nu*(z*nu-2.*rho);
	double sqrho=sqrt(1.-rho*rho);
	if(nu>0.00000001)
	{
		x=log((-rho+nu*z+sqrt(Iz2))/(1.-rho))/nu;

		f1overf0=exp(((2.*mu+b1*alpha*nu*rho)*(2.*rho*atan((z*nu-rho)/sqrho)+2.*rho*atan(rho/sqrho)+sqrho*log(Iz2)))/(4.*nu*nu*sqrho));
		
	}
	else
	{
		x=z;
		f1overf0=exp(mu*z*z/2.);
		
	}
	I0=sqrt(1.+z*nu*(z*nu-2.*rho));
	I1=(z*nu-rho)/I0;
	I2=(1.-rho*rho)/(I0*I0*I0);
	kappa=(12.*mu+(2*b2-3.*b1*b1)*alpha*alpha
		-4.*z*z*mu*mu+((3.*I0*I0-1.)*I0*I2-I1*I1)*nu*nu
		+(6.*b1*alpha-16.*z*mu)*nu*rho)
		/8.;
	double gintegral=GIntegralAnalytical2(-kappa,x/ARM_NumericConstants::ARM_SQRT_2,sqrt(T));
	return intrinsec+alpha*ARM_NumericConstants::ARM_INVSQRT2PI*sqB0Balphaz*f1overf0*sqrt(I0)*gintegral;
	
}

double SABRgeneric7(double intrinsec,double z,double b1,double b2,double zavg,double sqB0Balphaz,double T,double mu,double alpha,double rho,double nu)
{
	double x,I0,I1,I2,kappa,f1overf0;
	double Iz2=1.+z*nu*(z*nu-2.*rho);
	double sqrho=sqrt(1.-rho*rho);
	if(nu>0.00000001)
	{
		x=log((-rho+nu*z+sqrt(Iz2))/(1.-rho))/nu;

		f1overf0=exp(((2.*mu+b1*alpha*nu*rho)*(2.*rho*atan((z*nu-rho)/sqrho)+2.*rho*atan(rho/sqrho)+sqrho*log(Iz2)))/(4.*nu*nu*sqrho));
		
	}
	else
	{
		x=z;
		f1overf0=exp(mu*z*z/2.);
		
	}
	I0=sqrt(1.+z*nu*(z*nu-2.*rho));
	I1=(z*nu-rho)/I0;
	I2=(1.-rho*rho)/(I0*I0*I0);
	kappa=(12.*mu+(2*b2-3.*b1*b1)*alpha*alpha
		-4.*z*z*mu*mu+(sqrt(5.*I0*I0-1.)*I0*I2-I1*I1)*nu*nu
		+(6.*b1*alpha-16.*z*mu)*nu*rho)
		/8.;
	double gintegral=GIntegralAnalytical2(-kappa,x/ARM_NumericConstants::ARM_SQRT_2,sqrt(T));
	return intrinsec+alpha*ARM_NumericConstants::ARM_INVSQRT2PI*sqB0Balphaz*f1overf0*sqrt(I0)*gintegral;
	
}


/////////////////////////////////////////////////////////
///
///
///							 MonoSABR (Sans Flag) beta < 1, mu is the vol drift (in usual sabr, mu=0)
///
///
////////////////////////////////////////////////////////


double InfinitlyDerivableInterpolator(double x)  
	// monotone C-infini function  from {0,1}  to {0,1} with all derivatives =0 at boundaries
{
	return (tanh(tan((x-0.5)*ARM_NumericConstants::ARM_PI))+1.)/2.;
}





double SimplifiedSABR(double f,double K,double T,double mu,double alpha,double beta,double rho,double nu,int callput)
{
	if (K<=0)
	{
		switch (callput)
		{
		case K_CALL:
			{
				return f-K;
				break;
			}
		case K_PUT:
			{
				return 0.;
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SimplifiedSABR : callput  bad input :");
				break;
			}
		}
	}
	double intrinsec,intrinsec0;
	
	switch (callput)
	{
	case K_CALL:
		{
			intrinsec=(f>=K)?f-K:0;
			break;
		}
	case K_PUT:
		{
		intrinsec=(K>=f)?K-f:0;
			break;
		}
		
	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SimplifiedSABR : callput  bad input :");
			break;
		}
	}

	double LimK= T*f/100; /// limit for the interpolation schema at k=0

	double z;
	if (K> LimK)
	{
		if (beta<0.99999)
		{
			z=(pow(f,1.-beta)-pow(K,1.-beta))/(alpha*(1.-beta));
		}
		else
		{
			double z=log(f/K)/alpha;
		}
		double favg=(f+K)/2.;
		double b1=beta*pow(favg,beta-1.);
		double b2=beta*(beta-1.)*pow(favg,2.*beta-2.)+b1*b1;
		double sqB0Balphaz=pow(K*f,beta/2.);
		return SABRgeneric5( intrinsec, z, b1, b2, z, sqB0Balphaz, T, mu, alpha, rho, nu);
	}
	else
	{
		if (beta<0.99999)
		{
			z=(pow(f,1.-beta)-pow(LimK,1.-beta))/(alpha*(1.-beta));
		}
		else
		{
			z=log(f/K)/alpha;
		}
		double favg=(f+LimK)/2.;
		double b1=beta*pow(favg,beta-1.);
		double b2=beta*(beta-1.)*pow(favg,2.*beta-2.)+b1*b1;
		double sqB0Balphaz=pow(LimK*f,beta/2.);
		switch (callput)
		{
		case K_CALL:
			{
				intrinsec0=f-LimK;
				break;
			}
		case K_PUT:
			{
				intrinsec0=LimK-f;
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SimplifiedSABR : callput  bad input :");
				break;
			}
		}
		double V0= SABRgeneric5( intrinsec0, z, b1, b2, z, sqB0Balphaz, T, mu, alpha, rho, nu)-intrinsec0;
		if (beta<0.99999)
		{
			z=(pow(f,1.-beta)-pow(K,1.-beta))/(alpha*(1.-beta));
		}
		else
		{
			double z=log(f/K)/alpha;
		}
		favg=(f+K)/2.;
		b1=beta*pow(favg,beta-1.);
		b2=beta*(beta-1.)*pow(favg,2.*beta-2.)+b1*b1;
		sqB0Balphaz=pow(K*f,beta/2.);
		double Vt= SABRgeneric5( intrinsec, z, b1, b2, z, sqB0Balphaz, T, mu, alpha, rho, nu)-intrinsec;
		double epsilon=InfinitlyDerivableInterpolator(K/LimK);
		double alpha0 = V0/(LimK*LimK);
		return intrinsec+(Vt*epsilon+(1.-epsilon)*(K*K)*alpha0);

	}
	
}

///////////////////////////////////////////////////////////////////////////////////////////
///
///
///     beta =0    SABR variable 
///
///
///
///////////////////////////////////////////////////////////////////////////////////////////////




double BetaEqualZeroSABR(double f,double K,double T,double mu,double alpha,double rho,double nu,int callput)
{
	double intrinsec;
	
	switch (callput)
	{
	case K_CALL:
		{
			intrinsec=(f>=K)?f-K:0;
			break;
		}
	case K_PUT:
		{
		intrinsec=(K>=f)?K-f:0;
			break;
		}
		
	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BetaEqualZeroSABR : callput  bad input :");
			break;
		}
	}

	double LimK= T*f/100; /// limit for the interpolation schema at k=0
	double z;
	z=fabs(f-K)/alpha;
	double b1=0;
	double b2=0;
	double sqB0Balphaz=1.;
	if (K> LimK)
	{
		double favg=(f+K)/2.;
		return SABRgeneric5( intrinsec, z, b1, b2, z, sqB0Balphaz, T, mu, alpha, rho, nu);
	}
	else
	{
		double favg=(f+LimK)/2.;
		return SABRgeneric6( intrinsec, z, b1, b2, z, sqB0Balphaz, T, mu, alpha, rho, nu);	
	}
	
}
/////////////////////////////////////////////////////////
///
///
///							 biSABR
///
///
////////////////////////////////////////////////////////

double nuBiSABR(double F1,double alpha1,double beta1,double rho1,double nu1,
				
				double F2,double alpha2,double beta2,double rho2,double nu2,
				double rhos,double rhov,double rhoc12,double rhoc21,double X1,double X2,double alphas)
{
	double term1,term2,term3,term4,term5;
	double alpha1SQ=alpha1*alpha1;double alpha2SQ=alpha2*alpha2;
	double beta1SQ=beta1*beta1;double beta2SQ=beta2*beta2;
	double nu1SQ=nu1*nu1;double nu2SQ=nu2*nu2;
	double X1SQ=X1*X1;double X2SQ=X2*X2;
	double rhosSQ=rhos*rhos;
	double X1overF1=X1/F1;double X2overF2=X2/F2;
	double alphasSQ=alphas*alphas;

	term1=X1SQ*X1SQ*alpha1SQ*alpha1SQ*(X1overF1*X1overF1*alpha1SQ*beta1SQ+nu1SQ+2.*X1overF1*alpha1*beta1*nu1*rho1);
	term2=X2SQ*X2SQ*alpha2SQ*alpha2SQ*(X2overF2*X2overF2*alpha2SQ*beta2SQ+nu2SQ+2.*X2overF2*alpha2*beta2*nu2*rho2);
	
	term3=alpha1SQ*alpha2SQ*X1SQ*X2SQ*(nu1SQ*rhosSQ + nu2SQ*rhosSQ + 2*nu1*nu2*rhov + 2*nu1*nu2*rhosSQ*rhov +
		(alpha1SQ*beta1SQ*rhosSQ*X1SQ)/(F1*F1)+(2*alpha2*beta2*(nu2*rho2*rhosSQ+nu1*rhoc21*(1+rhosSQ))*X2)/F2 + 
		(2*alpha1*beta1*X1*(nu1*rho1*rhosSQ + nu2*rhoc12*(1 + rhosSQ) +
		(alpha2*beta2*rhos*(1 + rhosSQ)*X2)/F2))/F1 + (alpha2SQ*beta2SQ*rhosSQ*X2SQ)/(F2*F2));
	
	term4=2.*alpha1*alpha2*alpha2SQ*rhos*X1*(alpha2SQ*beta2SQ*pow(F2,-2. + 2.*beta2) + 
		(alpha2*beta2*pow(F2,-1. + beta2)*(2*F1*nu2*rho2 + F1*nu1*rhoc21 + alpha1*beta1*pow(F1,beta1)*rhos))/
		F1 + nu2*(nu2 + nu1*rhov + (alpha1*beta1*rhoc12*X1)/F1))*X2*X2SQ;

	term5=2.*alpha1*alpha1SQ*alpha2*rhos*X1*X1SQ*X2*
		((alpha1SQ*beta1SQ*X1SQ)/(F1*F1) + nu1*(nu1 + nu2*rhov + (alpha2*beta2*rhoc21*X2)/F2) + 
		(alpha1*beta1*X1*(2.*nu1*rho1 + nu2*rhoc12 + (alpha2*beta2*rhos*X2)/F2))/F1);
	return sqrt(term1+term2+term3-term4-term5)/alphasSQ;
}

double rhoBiSABR(double F1,double alpha1,double beta1,double rho1,double nu1,
				
				double F2,double alpha2,double beta2,double rho2,double nu2,
				double rhos,double rhov,double rhoc12,double rhoc21,double X1,double X2,double alphas,double nus)
{
	
	double alpha1SQ=alpha1*alpha1;double alpha2SQ=alpha2*alpha2;
	double X1SQ=X1*X1;double X2SQ=X2*X2;double alphasSQ=alphas*alphas;
	
	double result= (alpha1*alpha1SQ*X1*(nu1*rho1 + (alpha1*beta1*X1)/F1)*X1SQ - 
		alpha1SQ*alpha2*X1SQ*X2*(nu1*(rhoc21 + rho1*rhos) + 
		rhos*(nu2*rhoc12 + (2.*alpha1*beta1*X1)/F1 + (alpha2*beta2*rhos*X2)/F2)) - 
		alpha2*alpha2SQ*X2*(nu2*rho2 + (alpha2*beta2*X2)/F2)*X2SQ + 
		alpha1*alpha2SQ*X1*(nu2*(rhoc12 + rho2*rhos) + 
		rhos*(nu1*rhoc21 + (alpha1*beta1*rhos*X1)/F1 + 
		(2*alpha2*beta2*X2)/F2))*X2SQ)/(alphas*alphasSQ*nus);
	if (result>0.9999) return 0.9999;
	if (result<-0.9999) return -0.9999;
	return result;
}


double muBiSABR(double F1,double alpha1,double beta1,double rho1,double nu1,
				
				double F2,double alpha2,double beta2,double rho2,double nu2,
				double rhos,double rhov,double rhoc12,double rhoc21,double X1,double X2,double alphas)
{
	double alpha1SQ=alpha1*alpha1;double alpha2SQ=alpha2*alpha2;
	double X1SQ=X1*X1;double X2SQ=X2*X2;
	double F1SQ=F1*F1;double F2SQ=F2*F2;
	double nu1SQ=nu1*nu1;double nu2SQ=nu2*nu2;
	double rhosSQ=rhos*rhos;
	
	double term1=alpha1*alpha1SQ*alpha1SQ*beta1*F2SQ*X1*(2.*F1*nu1*rho1 + alpha1*(-1. + beta1)*X1)*X1SQ*X1SQ + 
		alpha2SQ*alpha2SQ*alpha2SQ*(-1. + beta2)*beta2*F1SQ*X2SQ*X2SQ*X2SQ;

	double term2=-3.*alpha1SQ*alpha1SQ*alpha2*beta1*F2SQ*rhos*(2.*F1*nu1*rho1 + alpha1*(-1. + beta1)*X1)*X1SQ*X1SQ*X2 + 
		alpha1*alpha2SQ*alpha2SQ*beta2*F1SQ*X1*(-6.*F2*nu2*rho2*rhos + 
		alpha1*(-1. + 2*beta2 + (-2. + beta2)*rhosSQ)*X1)*X2SQ*X2SQ + 
		alpha2*alpha2SQ*alpha2SQ*beta2*F1SQ*(2.*F2*nu2*rho2 - 3.*alpha1*(-1. + beta2)*rhos*X1)*X2*X2SQ*X2SQ;
	
	double term3=alpha1SQ*alpha2*alpha2SQ*(2*F1*F2*(-(beta1*F2*nu1*rho1*rhos) + 
        beta2*F1*(nu1*rhoc21*(-1. + rhosSQ) + nu2*rho2*(2. + rhosSQ))) + 
		alpha1*rhos*(-((-1. + beta2)*beta2*F1SQ) - (-1 + beta1)*beta1*F2SQ + 2*beta1*beta2*F1*F2*(-1. + rhosSQ))*X1)*
		X1SQ*X2*X2SQ;
	
	double term4=alpha1SQ*alpha2SQ*F2*X1SQ*(-(F1SQ*F2*(-1. + rhosSQ)*(nu1SQ + nu2SQ - 2*nu1*nu2*rhov)) + 
		2.*alpha1*F1*(-(beta2*F1*nu2*rho2*rhos) + beta1*F2*(nu2*rhoc12*(-1. + rhosSQ) + nu1*rho1*(2. + rhosSQ)))*X1 + 
		alpha1SQ*beta1*F2*(-1. + 2.*beta1 + (-2. + beta1)*rhosSQ)*X1SQ)*X2SQ;

	return (term1+term2+term3+term4)/(2.*F1SQ*F2SQ*alphas*alphas*alphas);
}

double BiSABRSpreadC(double F1,double F2,double alpha1,double alpha2,double beta1,double beta2,double rhos,double X1,double X2,double alphas,double s)
{
	double sF1beta;
	if(s+F2>0)
	{
		sF1beta=pow(s+F2,beta1);
	}
	else
	{
		sF1beta=0.0;
	}
	return sqrt(sF1beta*sF1beta*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*sF1beta*X2*alpha1*alpha2*rhos)/alphas;
}

double BiSABRSpreadC_derivative(double F1,double F2,double alpha1,double alpha2,double beta1,double beta2,double rhos,double X1,double X2,double alphas,double s)
{
	double sF1beta;
	if(s+F2>0)
	{
		sF1beta=pow(s+F2,beta1);
	}
	else
	{
		sF1beta=0.0;
	}
	double alphacurrent= sqrt(sF1beta*sF1beta*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*sF1beta*X2*alpha1*alpha2*rhos);
	return (sF1beta*sF1beta/(s+F2)*alpha1*alpha1*beta1-sF1beta/(s+F2)*X2*alpha1*alpha2*beta1*rhos)/(alphas*alphacurrent);
}

double BiSABRSpreadC_derivative2(double F1,double F2,double alpha1,double alpha2,double beta1,double beta2,double rhos,double X1,double X2,double alphas,double s)
{
	double sF1beta;
	if(s+F2>0)
	{
		sF1beta=pow(s+F2,beta1);
	}
	else
	{
		sF1beta=0.0;
	}
	double alphacurrent= sqrt(sF1beta*sF1beta*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*sF1beta*X2*alpha1*alpha2*rhos);
	double element1=(2.*sF1beta*sF1beta/(s+F2)*alpha1*alpha1*beta1-2.*sF1beta/(s+F2)*X2*alpha1*alpha2*beta1*rhos);
	double element2=(2.*sF1beta*sF1beta/(s+F2)/(s+F2)*alpha1*alpha1*beta1*(-1.+2.*beta1)-2.*sF1beta/(s+F2)/(s+F2)*X2*alpha1*alpha2*beta1*rhos*(-1.+beta1));
	return (element1*element1)/(4.*alphacurrent*alphacurrent*alphacurrent)+element2/(2.*alphacurrent);
}

double BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,double rhos,double rhov,double rhoc12,double rhoc21)
{
	/*
	if (F2<1e-10) return SimplifiedSABR( F1, K, T, 0.0, alpha1, beta1, rho1, nu1, K_CALL);
	if (F1<1e-10) return SimplifiedSABR( F2, -K, T, 0.0, alpha2, beta2, rho2, nu2, K_PUT);
	*/


	int SABRflag=ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL;
	if (F2<1e-10) 
	{	
		if(K<=0) return F1-K;
		double efficient_K=(K<F1/100 )?F1/100 : K;
		return BlackSholes_Formula(F1,
							SABR_ComputeImpliedVol(F1, efficient_K,T,alpha1,beta1,rho1, nu1,SABRflag )*sqrt(T),
							1.0,K,K_CALL);
	}
	if (F1<1e-10) 
	{
		if(K>=0) return 0;
		double efficient_K=(-K<F2/100 )?F2/100 : -K;
		return BlackSholes_Formula(F2,
							SABR_ComputeImpliedVol(F2,efficient_K,T,alpha2,beta2,rho2, nu2,SABRflag )*sqrt(T),
							1.0,efficient_K,K_PUT);
	}


	double intrinsec,z,b1,b2,sqB0Balphaz,favg,zavg,Ck,alphasp,nusp,rhosp,musp,X1,X2,Csecond,Cfavg;
	X1=pow(F1,beta1);X2=pow(F2,beta2);alphasp=sqrt(X1*X1*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*X1*X2*alpha1*alpha2*rhos);
	nusp= nuBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	rhosp= rhoBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp,nusp);
	musp= muBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	z=BiSABRz( K, F1, F2, alpha1, alpha2, rhos, beta1, beta2);
	favg=(F1-F2+K)/2.;
	zavg=BiSABRz( favg, F1, F2, alpha1, alpha2, rhos, beta1, beta2);
	b1= BiSABRSpreadC_derivative( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Csecond=BiSABRSpreadC_derivative2( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Cfavg=BiSABRSpreadC( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Ck=BiSABRSpreadC( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, K);
	b2=Csecond*Cfavg+b1*b1;
	sqB0Balphaz=sqrt(Ck);
	if (F1-F2-K>=0) 
	{
		intrinsec=F1-F2-K;
	}
	else
	{
		intrinsec=0;
	}

	return SABRgeneric5(intrinsec, z, b1, b2, zavg, sqB0Balphaz, T, musp, alphasp, rhosp, nusp);

}


/// correct the instability problem arising for high maturities and far from the money options
/// by using SABRgeneric6
///  inside  SABRgeneric6 we add a factor (3.*I0*I0-1.) instead of 2 inside the kappa, making the kappa 
/// match the hagan one for vanilla sabr
double BiSABR_SpreadOption2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,double rhos,double rhov,double rhoc12,double rhoc21)
{
	/*
	if (F2<1e-10) return SimplifiedSABR( F1, K, T, 0.0, alpha1, beta1, rho1, nu1, K_CALL);
	if (F1<1e-10) return SimplifiedSABR( F2, -K, T, 0.0, alpha2, beta2, rho2, nu2, K_PUT);
	*/


	int SABRflag=ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL;
	if (F2<1e-10) 
	{	
		if(K<=0) return F1-K;
		double efficient_K=(K<F1/100 )?F1/100 : K;
		return BlackSholes_Formula(F1,
							SABR_ComputeImpliedVol(F1, efficient_K,T,alpha1,beta1,rho1, nu1,SABRflag )*sqrt(T),
							1.0,K,K_CALL);
	}
	if (F1<1e-10) 
	{
		if(K>=0) return 0;
		double efficient_K=(-K<F2/100 )?F2/100 : -K;
		return BlackSholes_Formula(F2,
							SABR_ComputeImpliedVol(F2,efficient_K,T,alpha2,beta2,rho2, nu2,SABRflag )*sqrt(T),
							1.0,efficient_K,K_PUT);
	}


	double intrinsec,z,b1,b2,sqB0Balphaz,favg,zavg,Ck,alphasp,nusp,rhosp,musp,X1,X2,Csecond,Cfavg;
	X1=pow(F1,beta1);X2=pow(F2,beta2);alphasp=sqrt(X1*X1*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*X1*X2*alpha1*alpha2*rhos);
	nusp= nuBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	rhosp= rhoBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp,nusp);
	musp= muBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	z=BiSABRz( K, F1, F2, alpha1, alpha2, rhos, beta1, beta2);
	favg=(F1-F2+K)/2.;
	zavg=BiSABRz( favg, F1, F2, alpha1, alpha2, rhos, beta1, beta2);
	b1= BiSABRSpreadC_derivative( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Csecond=BiSABRSpreadC_derivative2( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Cfavg=BiSABRSpreadC( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Ck=BiSABRSpreadC( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, K);
	b2=Csecond*Cfavg+b1*b1;
	sqB0Balphaz=sqrt(Ck);
	if (F1-F2-K>=0) 
	{
		intrinsec=F1-F2-K;
	}
	else
	{
		intrinsec=0;
	}

	return SABRgeneric6(intrinsec, z, b1, b2, zavg, sqB0Balphaz, T, musp, alphasp, rhosp, nusp);

}


/// correct the instability problem arising for high maturities and far from the money options
/// by using SABRgeneric7
double BiSABR_SpreadOption3(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,double rhos,double rhov,double rhoc12,double rhoc21)
{
	/*
	if (F2<1e-10) return SimplifiedSABR( F1, K, T, 0.0, alpha1, beta1, rho1, nu1, K_CALL);
	if (F1<1e-10) return SimplifiedSABR( F2, -K, T, 0.0, alpha2, beta2, rho2, nu2, K_PUT);
	*/


	int SABRflag=ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL;
	if (F2<1e-10) 
	{	
		if(K<=0) return F1-K;
		double efficient_K=(K<F1/100 )?F1/100 : K;
		return BlackSholes_Formula(F1,
							SABR_ComputeImpliedVol(F1, efficient_K,T,alpha1,beta1,rho1, nu1,SABRflag )*sqrt(T),
							1.0,K,K_CALL);
	}
	if (F1<1e-10) 
	{
		if(K>=0) return 0;
		double efficient_K=(-K<F2/100 )?F2/100 : -K;
		return BlackSholes_Formula(F2,
							SABR_ComputeImpliedVol(F2,efficient_K,T,alpha2,beta2,rho2, nu2,SABRflag )*sqrt(T),
							1.0,efficient_K,K_PUT);
	}


	double intrinsec,z,b1,b2,sqB0Balphaz,favg,zavg,Ck,alphasp,nusp,rhosp,musp,X1,X2,Csecond,Cfavg;
	X1=pow(F1,beta1);X2=pow(F2,beta2);alphasp=sqrt(X1*X1*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*X1*X2*alpha1*alpha2*rhos);
	nusp= nuBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	rhosp= rhoBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp,nusp);
	musp= muBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	z=BiSABRz( K, F1, F2, alpha1, alpha2, rhos, beta1, beta2);
	favg=(F1-F2+K)/2.;
	zavg=BiSABRz( favg, F1, F2, alpha1, alpha2, rhos, beta1, beta2);
	b1= BiSABRSpreadC_derivative( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Csecond=BiSABRSpreadC_derivative2( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Cfavg=BiSABRSpreadC( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Ck=BiSABRSpreadC( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, K);
	b2=Csecond*Cfavg+b1*b1;
	sqB0Balphaz=sqrt(Ck);
	if (F1-F2-K>=0) 
	{
		intrinsec=F1-F2-K;
	}
	else
	{
		intrinsec=0;
	}

	return SABRgeneric7(intrinsec, z, b1, b2, zavg, sqB0Balphaz, T, musp, alphasp, rhosp, nusp);

}

/// correct the instability problem arising for high maturities and far from the money options
/// by using a projection on normal heston process instead of using hagan approximation
/// formula for transforming local vol model into an option price (a SABRgeneric** function)
double BiSABR_SpreadOption4(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,double rhos,double rhov,double rhoc12,double rhoc21)
{
	double kappaV,lambdaV,thetaV,lambdaShift,S0;
	int SABRflag=ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL;
	if (F2<1e-10) 
	{	
		if(K<=0) return F1-K;
		double efficient_K=(K<F1/100 )?F1/100 : K;
		return BlackSholes_Formula(F1,
							SABR_ComputeImpliedVol(F1, efficient_K,T,alpha1,beta1,rho1, nu1,SABRflag )*sqrt(T),
							1.0,K,K_CALL);
	}
	if (F1<1e-10) 
	{
		if(K>=0) return 0;
		double efficient_K=(-K<F2/100 )?F2/100 : -K;
		return BlackSholes_Formula(F2,
							SABR_ComputeImpliedVol(F2,efficient_K,T,alpha2,beta2,rho2, nu2,SABRflag )*sqrt(T),
							1.0,efficient_K,K_PUT);
	}


	double alphasp,nusp,rhosp,musp,X1,X2;
	X1=pow(F1,beta1);X2=pow(F2,beta2);
	alphasp=sqrt(X1*X1*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*X1*X2*alpha1*alpha2*rhos);
	nusp= nuBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	rhosp= rhoBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp,nusp);
	musp= muBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	kappaV=2.*nusp*alphasp;
	lambdaV=(kappaV*kappaV/(2.*alphasp*alphasp)-(2.*musp+nusp*nusp));
	thetaV=kappaV*kappaV/(2.*lambdaV);
	lambdaShift=0.1;
	S0=F1-F2;
	double prec=1.e-8;
	double nbfirst=60,nb=60,NbStage=-1,NbOscill=0;
	return NormalHeston(rhosp,lambdaV,thetaV,kappaV,alphasp*alphasp,S0,K,T,K_CALL,lambdaShift,nbfirst,nb,NbStage,NbOscill,prec);

}

///
///  Correct the slight bias when comparing with bi log normal spreadoption
///
double BiSABR_SpreadOption_bilog_corrected(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,double rhos,double rhov,double rhoc12,double rhoc21)
{
	int nb_gauss=120;
	double option=BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, rhos, rhov, rhoc12, rhoc21);
	double option_bisabr_limit=BiSABR_SpreadOption(F1,alpha1,1.0,0,0,F2,alpha2,1.0,0,0,K,T,rhos,0,0,0);
	double option_log_normal_limit=Export_LogNormal_SpreadOption( F1, F2, alpha1, alpha2, rhos, K, T, K_CALL, 
		ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION, nb_gauss);
	
	if((option_bisabr_limit>0)&&(option_log_normal_limit>1e-10))
	{
		return option*pow((option_log_normal_limit/option_bisabr_limit),pow((beta1+beta2)/2.,8.));
	}
	else
	{
		return option;
	}
}

double Complete_BiSABR_SpreadOption_bilog_corrected(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
		switch (CallPut)
		{
		case K_CALL:
			{
///  Corrects the put call parity bias if it happens !				
								 
				return 0.5*(BiSABR_SpreadOption_bilog_corrected( F1, alpha1, beta1, rho1, nu1, 
												 F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
							BiSABR_SpreadOption_bilog_corrected( F2, alpha2, beta2, rho2, nu2, 
												 F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
							F1-F2-K);

				break;
			}
		case K_PUT:
			{
				return 0.5*(BiSABR_SpreadOption_bilog_corrected( F1, alpha1, beta1, rho1, nu1, 
												 F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
							BiSABR_SpreadOption_bilog_corrected( F2, alpha2, beta2, rho2, nu2, 
												 F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
							K-F1+F2);
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Complete_BiSABR_SpreadOption_bilog_corrected : callput  bad input :");
				break;
			}
		}
}

double Complete_BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
		switch (CallPut)
		{
		case K_CALL:
			{				 

				return 0.5*(BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, 
												 F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
							BiSABR_SpreadOption( F2, alpha2, beta2, rho2, nu2, 
												 F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
							F1-F2-K);

				break;
			}
		case K_PUT:
			{

				return 0.5*(BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, 
												 F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
							BiSABR_SpreadOption( F2, alpha2, beta2, rho2, nu2, 
												 F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
							K-F1+F2);
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Complete_BiSABR_SpreadOption : callput  bad input :");
				break;
			}
		}
}

double Complete_BiSABR_SpreadOption_2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
		switch (CallPut)
		{
		case K_CALL:
			{				 

				return 0.5*(BiSABR_SpreadOption2( F1, alpha1, beta1, rho1, nu1, 
												 F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
							BiSABR_SpreadOption2( F2, alpha2, beta2, rho2, nu2, 
												 F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
							F1-F2-K);

				break;
			}
		case K_PUT:
			{

				return 0.5*(BiSABR_SpreadOption2( F1, alpha1, beta1, rho1, nu1, 
												 F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
							BiSABR_SpreadOption2( F2, alpha2, beta2, rho2, nu2, 
												 F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
							K-F1+F2);
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Complete_BiSABR_SpreadOption : callput  bad input :");
				break;
			}
		}
}



////////////////////////////////////////////////////////////////////////////
///
///
///               Packaged form
///
///
///////////////////////////////////////////////////////////////////////////

complex<double> cplx(double x)
{
	complex<double> y(x,0);
	return y;
}

double Packaged_BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
								  double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21,int flag)
{
	switch (CallPut)
	{
	case K_CALL:
		{
			switch (flag)
			{
				case 0:
				{
					return 0.5*(BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, 
						F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
						BiSABR_SpreadOption( F2, alpha2, beta2, rho2, nu2, 
						F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
						F1-F2-K);
					break;
				}
				case 2:
				{
					return 0.5*(BiSABR_SpreadOption2( F1, alpha1, beta1, rho1, nu1, 
						F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
						BiSABR_SpreadOption2( F2, alpha2, beta2, rho2, nu2, 
						F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
						F1-F2-K);
					break;
				}
				case 4:
				{
					return 0.5*(BiSABR_SpreadOption4( F1, alpha1, beta1, rho1, nu1, 
						F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
						BiSABR_SpreadOption4( F2, alpha2, beta2, rho2, nu2, 
						F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
						F1-F2-K);
					break;
				}
				case 10:
				{
					return 0.5*(real(Complexified_BiSABR_SpreadOption( cplx(F1), cplx(alpha1), cplx(beta1), cplx(rho1), cplx(nu1), 
						cplx(F2), cplx(alpha2), cplx(beta2), cplx(rho2), cplx(nu2),cplx(K), cplx(T), cplx(rhos), cplx(rhov), cplx(rhoc12), cplx(rhoc21)))+
						real(Complexified_BiSABR_SpreadOption( cplx(F2), cplx(alpha2), cplx(beta2), cplx(rho2), cplx(nu2), 
						cplx(F1), cplx(alpha1), cplx(beta1), cplx(rho1), cplx(nu1),-cplx(K), cplx(T), cplx(rhos), cplx(rhov), cplx(rhoc21), cplx(rhoc12)))+
						F1-F2-K);
					break;
				}
				case 1:
				{
					return 0.5*(BiSABR_SpreadOption_bilog_corrected( F1, alpha1, beta1, rho1, nu1, 
						F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
						BiSABR_SpreadOption_bilog_corrected( F2, alpha2, beta2, rho2, nu2, 
						F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
						F1-F2-K);
					break;
				}
				default:
					{
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_BiSABR_SpreadOption : flag  bad input :");
						break;
					}
					
			}
			break;
		}
	case K_PUT:
		{
			
			switch (flag)
			{
			case 0:
				{
					return 0.5*(BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, 
						F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
						BiSABR_SpreadOption( F2, alpha2, beta2, rho2, nu2, 
						F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
						K-F1+F2);
				}
			case 2:
				{
					return 0.5*(BiSABR_SpreadOption2( F1, alpha1, beta1, rho1, nu1, 
						F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
						BiSABR_SpreadOption2( F2, alpha2, beta2, rho2, nu2, 
						F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
						K-F1+F2);
				}
			case 4:
				{
					return 0.5*(BiSABR_SpreadOption4( F1, alpha1, beta1, rho1, nu1, 
						F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
						BiSABR_SpreadOption4( F2, alpha2, beta2, rho2, nu2, 
						F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
						K-F1+F2);
				}
			case 10:
				{
					return 0.5*(real(Complexified_BiSABR_SpreadOption( cplx(F1), cplx(alpha1), cplx(beta1), cplx(rho1), cplx(nu1), 
						cplx(F2), cplx(alpha2), cplx(beta2), cplx(rho2), cplx(nu2),cplx(K), cplx(T), cplx(rhos), cplx(rhov), cplx(rhoc12), cplx(rhoc21)))+
						real(Complexified_BiSABR_SpreadOption( cplx(F2), cplx(alpha2), cplx(beta2), cplx(rho2), cplx(nu2), 
						cplx(F1), cplx(alpha1), cplx(beta1), cplx(rho1), cplx(nu1),-cplx(K), cplx(T), cplx(rhos), cplx(rhov), cplx(rhoc21), cplx(rhoc12)))+
						K-F1+F2);
				}
			case 1:
				{
					return 0.5*(BiSABR_SpreadOption_bilog_corrected( F1, alpha1, beta1, rho1, nu1, 
						F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21)+
						BiSABR_SpreadOption_bilog_corrected( F2, alpha2, beta2, rho2, nu2, 
						F1, alpha1, beta1, rho1, nu1,-K, T, rhos, rhov, rhoc21, rhoc12)+
						K-F1+F2);
				}
			default:
				{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_BiSABR_SpreadOption : flag  bad input :");
					break;
				}
				
			}
			break;
		}
		
	default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_BiSABR_SpreadOption : callput  bad input :");
				break;
			}
		}
}

/// the following extension is mainly used for using negative F1 and/or  F2



double Packaged_Complexified_BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,
								  double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int callput,double rhos,double rhov,double rhoc12,double rhoc21,int flag)
{
		switch (callput)
	{
	case K_CALL:
		{
			return 0.5*(real(Complexified_BiSABR_SpreadOption( cplx(F1), cplx(alpha1), cplx(beta1), cplx(rho1), cplx(nu1), 
					cplx(F2), cplx(alpha2), cplx(beta2), cplx(rho2), cplx(nu2),cplx(K), cplx(T), cplx(rhos), cplx(rhov), cplx(rhoc12), cplx(rhoc21)))+
					real(Complexified_BiSABR_SpreadOption( cplx(F2), cplx(alpha2), cplx(beta2), cplx(rho2), cplx(nu2), 
					cplx(F1), cplx(alpha1), cplx(beta1), cplx(rho1), cplx(nu1),-cplx(K), cplx(T), cplx(rhos), cplx(rhov), cplx(rhoc21), cplx(rhoc12)))+
					F1-F2-K);
			break;
		}
	case K_PUT:
		{
			return 0.5*(real(Complexified_BiSABR_SpreadOption( cplx(F1), cplx(alpha1), cplx(beta1), cplx(rho1), cplx(nu1), 
					cplx(F2), cplx(alpha2), cplx(beta2), cplx(rho2), cplx(nu2),cplx(K), cplx(T), cplx(rhos), cplx(rhov), cplx(rhoc12), cplx(rhoc21)))+
					real(Complexified_BiSABR_SpreadOption( cplx(F2), cplx(alpha2), cplx(beta2), cplx(rho2), cplx(nu2), 
					cplx(F1), cplx(alpha1), cplx(beta1), cplx(rho1), cplx(nu1),-cplx(K), cplx(T), cplx(rhos), cplx(rhov), cplx(rhoc21), cplx(rhoc12)))+
					K-F1+F2);
			break;
		}
		
	default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_Complexified_BiSABR_SpreadOption : callput  bad input :");
				break;
			}
		}
}


CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/