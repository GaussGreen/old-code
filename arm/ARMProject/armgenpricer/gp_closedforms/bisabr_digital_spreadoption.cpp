/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_digital_spreadoption.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date August 2006
 */

#include "firsttoinc.h"

#include <cmath>
#include <complex>

#include <vector>
#include <iostream>
#include <iomanip>
#include "expt.h"
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
#include "gpclosedforms/extendedsabrformula.h"

using namespace std; 


CC_BEGIN_NAMESPACE(ARM)


#define ARM_CF_NBSTEPS_FOR_CERF 10
#define ARM_CF_K_SHIFT_FOR_DERIVATION 0.0001
#define ARM_CF_K_S1_S2_PARTITION_FACTOR 0.5

////////////////////////////////////////////////////////////////////////////////////////////
///
///			complex Arc Tangent
///
////////////////////////////////////////////////////////////////////////////////////////////
complex<double> atan(complex<double> x)
{
	complex<double> I(0,1.);
	complex<double> deux(2.,0);
	return -I/deux * log((I-x)/(I+x));
}


////////////////////////////////////////////////////////////////////////////////////////////
///
///			complex power redefined to correct a bug and introduce some optimization
///
////////////////////////////////////////////////////////////////////////////////////////////
complex<double> pow(complex<double> x,complex<double> y)
{
	complex<double> moinsun(-1.0,0);
	complex<double> deux(2.0,0);
	complex<double> trois(3.0,0);
	if(y==moinsun) return 1./x;
	if(y==deux) return x*x;
	if(y==trois) return x*x*x;
	complex<double> loga=log(x);
	return exp(loga*y);
}

double pow(double x,double y)
{
	if(y==-1.0) return 1./x;
	if(y==2.) return x*x;
	if(y==3.) return x*x*x;
	return exp(log(x)*y);
}


/// All the following computation ar done in complex to 
/// do the analytic continuation of th bisabr pricing function 

////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of sum{0,x} of exp(-a^2*s^2-b^2/s^2
///
////////////////////////////////////////////////////////////////////////////////////////////

complex<double> GIntegralAnalytical(complex<double> a1, complex<double> b1,complex<double> x1)
{
	complex<double> x=abs(x1);
	complex<double> b=abs(b1);
	complex<double> a=abs(a1);
	complex<double> signx1;
	complex<double> un(1.,0);
	complex<double> deux(2.,0);
	complex<double> quatre(4.,0);
	complex<double> trois(3.,0);
	complex<double> dix(10.,0);
	complex<double> trente(30.,0);
	complex<double> quinze(15.,0);

	complex<double> SQRT_PI(ARM_NumericConstants::ARM_SQRT_PI,0);
	complex<double> SQRT_2(ARM_NumericConstants::ARM_SQRT_2,0);
	if(real(x1)>0) signx1=1.0; else signx1=-1.0;
	if(abs(a)>0.0001)
	{
		return signx1*(exp(-deux*a*b)*SQRT_PI*
			(un-exp(quatre*a*b)-cnormal(SQRT_2*(b/x-a*x),ARM_CF_NBSTEPS_FOR_CERF)+
			exp(quatre*a*b)*cnormal(SQRT_2*(b/x+a*x),ARM_CF_NBSTEPS_FOR_CERF)))/(deux*a);
	}
	else
	{
		return signx1*(exp(-b*b/(x*x))*x*(trente+dix*a*a*(deux*b*b-x*x)+a*a*a*a*(quatre*b*b*b*b-
			deux*b*b*x*x+trois*x*x*x*x))+quatre*b*(quinze+dix*a*a*b*b+deux*a*a*a*a*b*b*b*b)*
			SQRT_PI*(-un+cnormal(SQRT_2*b/x,ARM_CF_NBSTEPS_FOR_CERF)))/trente;
	}
}

///  GIntegralAnalytical2(complex<double> a2, complex<double> b1,complex<double> x1) == GIntegralAnalytical(complex<double> sqrt(a2), complex<double> b1,complex<double> x1)
/// also works when a2 <0 ! it is a taylor expansion valid for abs(a2)<1
////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of sum{0,x} of exp(-a2*s^2-b1^2/s^2
///
////////////////////////////////////////////////////////////////////////////////////////////

complex<double> GIntegralAnalytical2(complex<double> a2, complex<double> b1,complex<double> x1)
{

	{
		
		complex<double> x=abs(x1);
		complex<double> b=abs(b1);
		complex<double> argERF=ARM_NumericConstants::ARM_SQRT_2*b/x;
		complex<double> ERF=2.*cnormal(argERF,ARM_CF_NBSTEPS_FOR_CERF)-2.;
		complex<double> expb2=exp(b*b/(x*x));
		complex<double> b2=b*b;
		complex<double> b3=b2*b;
		complex<double> b5=b3*b2;
		complex<double> x2=x*x;
		complex<double> x3=x2*x;
		complex<double> x5=x3*x*x;
		complex<double> a4=a2*a2;
		complex<double> a0=ERF*b*ARM_NumericConstants::ARM_SQRT_PI+x/expb2;
		complex<double> b0=2.*b3*ERF*ARM_NumericConstants::ARM_SQRT_PI+(2.*b2*x-x3)/expb2;
		complex<double> c0=4.*b3*b2*ERF*ARM_NumericConstants::ARM_SQRT_PI+(4.*b3*b*x-2.*b2*x3+3.*x3*x2)/expb2;
		complex<double> d0=8.*b5*b2*ERF*ARM_NumericConstants::ARM_SQRT_PI+(8.*b3*b3*x-4.*b3*b*x3+6.*b2*x3*x2-15.*x5*x2)/expb2;
		complex<double> e0=16.*b5*b2*b2*ERF*ARM_NumericConstants::ARM_SQRT_PI+(16.*b5*b3*x-8.*b5*b*x3+12.*b3*b*x5-30.*b2*x5*x2+105.*x5*x3*x)/expb2;
		complex<double> res;
		if(real(x1)>0)
		{
			res= a0+b0*a2/3.+c0*a4/30.+d0*a4*a2/630.+e0*a4*a4/22680.;
		}
		else
		{
			res= -(a0+b0*a2/3.+c0*a4/30.+d0*a4*a2/630.+e0*a4*a4/22680.);
		}
		return res;
	}
}


complex<double> BiSABRIntegral(complex<double> A,complex<double> B,complex<double> G,complex<double> F,complex<double> gamma)
{
	complex<double> xpoint,xweight;
	int nbl=180;
	int i;
	GaussLegendre_Coefficients d(nbl, real(A), real(B));
	complex<double> Sum=0;
	for(i=0;i<nbl;i++)
	{
		xpoint=d.get_point(i);
		xweight=d.get_weight(i);
		Sum+=pow(xpoint,gamma)/sqrt(xpoint*xpoint-xpoint*G+F)*xweight;
	}
	return Sum;
}

complex<double> BiSABRz(complex<double> K,complex<double> F1,complex<double> F2,complex<double> alpha1,complex<double> alpha2,complex<double> rhos,complex<double> beta1,complex<double> beta2)
{
	complex<double> F2beta2=pow(F2,beta2);
	complex<double> alpha2overalpha1=alpha2/alpha1;
	complex<double> un(1.,0);
	complex<double> deux(2.,0);
	if(real(K+F2)>0)
	{
	return BiSABRIntegral(pow(K+F2,beta1),pow(F1,beta1),deux*F2beta2*alpha2overalpha1*rhos,
		F2beta2*F2beta2*alpha2overalpha1*alpha2overalpha1,(un-beta1)/beta1)/(alpha1*beta1);
	}
	else
	{
	return BiSABRIntegral(0,pow(F1,beta1),deux*F2beta2*alpha2overalpha1*rhos,
		F2beta2*F2beta2*alpha2overalpha1*alpha2overalpha1,(un-beta1)/beta1)/(alpha1*beta1);

	}
}


complex<double> SABRgeneric5(complex<double> intrinsec,complex<double> z,complex<double> b1,complex<double> b2,complex<double> zavg,complex<double> sqB0Balphaz,complex<double> T,complex<double> mu,complex<double> alpha,complex<double> rho,complex<double> nu)
{
	complex<double> x,I0,I1,I2,kappa,f1overf0;
	complex<double> un(1.,0);
	complex<double> deux(2.,0);
	complex<double> trois(3.,0);
	complex<double> quatre(4.,0);
	complex<double> six(6.,0);
	complex<double> seize(16.,0);
	complex<double> douze(12.,0);
	complex<double> huit(8.,0);
	complex<double> INVSQRT2PI(ARM_NumericConstants::ARM_INVSQRT2PI,0);
	complex<double> SQRT_2(ARM_NumericConstants::ARM_SQRT_2,0);
	complex<double> Iz2;
	complex<double> sqrho;
	
	Iz2=un+z*nu*(z*nu-deux*rho);
	sqrho=sqrt(un-rho*rho);
	if(real(nu)>0.00000001)
	{
		x=log((-rho+nu*z+sqrt(Iz2))/(un-rho))/nu;
		f1overf0=exp(((deux*mu+b1*alpha*nu*rho)*(deux*rho*atan((z*nu-rho)/sqrho)+deux*rho*atan(rho/sqrho)+sqrho*log(Iz2)))/(deux*deux*nu*nu*sqrho));
		
	}
	else
	{
		x=z;
		f1overf0=exp(mu*z*z/2.);
		
	}
	I0=sqrt(un+zavg*nu*(zavg*nu-deux*rho));
	I1=(zavg*nu-rho)/I0;
	I2=(un-rho*rho)/(I0*I0*I0);
	kappa=(douze*mu+(deux*b2-trois*b1*b1)*alpha*alpha
		-quatre*z*z*mu*mu+(deux*I0*I2-I1*I1)*nu*nu
		+(six*b1*alpha-seize*zavg*mu)*nu*rho)
		/(huit*(un+zavg*nu*(z*nu-deux*rho)));
	complex<double> gintegral=GIntegralAnalytical2(-kappa,x/SQRT_2,sqrt(T)); /// l'argument fournit en premiere position de GIntegralAnalytical2 est le carre de celui similaire dans mathematica
	return intrinsec+alpha*INVSQRT2PI*sqB0Balphaz*f1overf0*sqrt(I0)*gintegral;
	
}

/// correct for the long term instability
complex<double> SABRgeneric6(complex<double> intrinsec,complex<double> z,complex<double> b1,complex<double> b2,complex<double> zavg,complex<double> sqB0Balphaz,complex<double> T,complex<double> mu,complex<double> alpha,complex<double> rho,complex<double> nu)
{
	complex<double> x,I0,I1,I2,kappa,f1overf0;
	complex<double> un(1.,0);
	complex<double> deux(2.,0);
	complex<double> trois(3.,0);
	complex<double> quatre(4.,0);
	complex<double> six(6.,0);
	complex<double> seize(16.,0);
	complex<double> douze(12.,0);
	complex<double> huit(8.,0);
	complex<double> INVSQRT2PI(ARM_NumericConstants::ARM_INVSQRT2PI,0);
	complex<double> SQRT_2(ARM_NumericConstants::ARM_SQRT_2,0);
	complex<double> Iz2;
	complex<double> sqrho;
	
	Iz2=un+z*nu*(z*nu-deux*rho);
	sqrho=sqrt(un-rho*rho);
	if(real(nu)>0.00000001)
	{
		x=log((-rho+nu*z+sqrt(Iz2))/(un-rho))/nu;
		f1overf0=exp(((deux*mu+b1*alpha*nu*rho)*(deux*rho*atan((z*nu-rho)/sqrho)+deux*rho*atan(rho/sqrho)+sqrho*log(Iz2)))/(deux*deux*nu*nu*sqrho));
		
	}
	else
	{
		x=z;
		f1overf0=exp(mu*z*z/2.);
		
	}
	I0=sqrt(un+zavg*nu*(zavg*nu-deux*rho));
	I1=(zavg*nu-rho)/I0;
	I2=(un-rho*rho)/(I0*I0*I0);
	kappa=(douze*mu+(deux*b2-trois*b1*b1)*alpha*alpha
		-quatre*z*z*mu*mu+((trois*I0*I0-un)*I0*I2-I1*I1)*nu*nu
		+(six*b1*alpha-seize*zavg*mu)*nu*rho)
		/huit;
	complex<double> gintegral=GIntegralAnalytical2(-kappa,x/SQRT_2,sqrt(T)); /// l'argument fournit en premiere position de GIntegralAnalytical2 est le carre de celui similaire dans mathematica
	return intrinsec+alpha*INVSQRT2PI*sqB0Balphaz*f1overf0*sqrt(I0)*gintegral;
	
}


complex<double> nuBiSABR(complex<double> F1,complex<double> alpha1,complex<double> beta1,complex<double> rho1,complex<double> nu1,
				
				complex<double> F2,complex<double> alpha2,complex<double> beta2,complex<double> rho2,complex<double> nu2,
				complex<double> rhos,complex<double> rhov,complex<double> rhoc12,complex<double> rhoc21,complex<double> X1,complex<double> X2,complex<double> alphas)
{
	complex<double> un(1.,0);
	complex<double> deux(2.,0);
	complex<double> term1,term2,term3,term4,term5;
	complex<double> alpha1SQ=alpha1*alpha1;complex<double> alpha2SQ=alpha2*alpha2;
	complex<double> beta1SQ=beta1*beta1;complex<double> beta2SQ=beta2*beta2;
	complex<double> nu1SQ=nu1*nu1;complex<double> nu2SQ=nu2*nu2;
	complex<double> X1SQ=X1*X1;complex<double> X2SQ=X2*X2;
	complex<double> rhosSQ=rhos*rhos;
	complex<double> X1overF1=X1/F1;complex<double> X2overF2=X2/F2;
	complex<double> alphasSQ=alphas*alphas;

	term1=X1SQ*X1SQ*alpha1SQ*alpha1SQ*(X1overF1*X1overF1*alpha1SQ*beta1SQ+nu1SQ+deux*X1overF1*alpha1*beta1*nu1*rho1);
	term2=X2SQ*X2SQ*alpha2SQ*alpha2SQ*(X2overF2*X2overF2*alpha2SQ*beta2SQ+nu2SQ+deux*X2overF2*alpha2*beta2*nu2*rho2);
	
	term3=alpha1SQ*alpha2SQ*X1SQ*X2SQ*(nu1SQ*rhosSQ + nu2SQ*rhosSQ + deux*nu1*nu2*rhov + deux*nu1*nu2*rhosSQ*rhov +
		(alpha1SQ*beta1SQ*rhosSQ*X1SQ)/(F1*F1)+(deux*alpha2*beta2*(nu2*rho2*rhosSQ+nu1*rhoc21*(un+rhosSQ))*X2)/F2 + 
		(deux*alpha1*beta1*X1*(nu1*rho1*rhosSQ + nu2*rhoc12*(un+ rhosSQ) +
		(alpha2*beta2*rhos*(un + rhosSQ)*X2)/F2))/F1 + (alpha2SQ*beta2SQ*rhosSQ*X2SQ)/(F2*F2));
	
	term4=deux*alpha1*alpha2*alpha2SQ*rhos*X1*(alpha2SQ*beta2SQ*pow(F2,-deux + deux*beta2) + 
		(alpha2*beta2*pow(F2,-un + beta2)*(deux*F1*nu2*rho2 + F1*nu1*rhoc21 + alpha1*beta1*pow(F1,beta1)*rhos))/
		F1 + nu2*(nu2 + nu1*rhov + (alpha1*beta1*rhoc12*X1)/F1))*X2*X2SQ;

	term5=deux*alpha1*alpha1SQ*alpha2*rhos*X1*X1SQ*X2*
		((alpha1SQ*beta1SQ*X1SQ)/(F1*F1) + nu1*(nu1 + nu2*rhov + (alpha2*beta2*rhoc21*X2)/F2) + 
		(alpha1*beta1*X1*(2.*nu1*rho1 + nu2*rhoc12 + (alpha2*beta2*rhos*X2)/F2))/F1);
	return sqrt(term1+term2+term3-term4-term5)/alphasSQ;
}

complex<double> rhoBiSABR(complex<double> F1,complex<double> alpha1,complex<double> beta1,complex<double> rho1,complex<double> nu1,
				
				complex<double> F2,complex<double> alpha2,complex<double> beta2,complex<double> rho2,complex<double> nu2,
				complex<double> rhos,complex<double> rhov,complex<double> rhoc12,complex<double> rhoc21,complex<double> X1,complex<double> X2,complex<double> alphas,complex<double> nus)
{
	complex<double> un(1.,0);
	complex<double> deux(2.,0);
	complex<double> alpha1SQ=alpha1*alpha1;complex<double> alpha2SQ=alpha2*alpha2;
	complex<double> X1SQ=X1*X1;complex<double> X2SQ=X2*X2;complex<double> alphasSQ=alphas*alphas;
	
	return (alpha1*alpha1SQ*X1*(nu1*rho1 + (alpha1*beta1*X1)/F1)*X1SQ - 
		alpha1SQ*alpha2*X1SQ*X2*(nu1*(rhoc21 + rho1*rhos) + 
		rhos*(nu2*rhoc12 + (deux*alpha1*beta1*X1)/F1 + (alpha2*beta2*rhos*X2)/F2)) - 
		alpha2*alpha2SQ*X2*(nu2*rho2 + (alpha2*beta2*X2)/F2)*X2SQ + 
		alpha1*alpha2SQ*X1*(nu2*(rhoc12 + rho2*rhos) + 
		rhos*(nu1*rhoc21 + (alpha1*beta1*rhos*X1)/F1 + 
		(deux*alpha2*beta2*X2)/F2))*X2SQ)/(alphas*alphasSQ*nus);
}


complex<double> muBiSABR(complex<double> F1,complex<double> alpha1,complex<double> beta1,complex<double> rho1,complex<double> nu1,
				
				complex<double> F2,complex<double> alpha2,complex<double> beta2,complex<double> rho2,complex<double> nu2,
				complex<double> rhos,complex<double> rhov,complex<double> rhoc12,complex<double> rhoc21,complex<double> X1,complex<double> X2,complex<double> alphas)
{
	complex<double> alpha1SQ=alpha1*alpha1;complex<double> alpha2SQ=alpha2*alpha2;
	complex<double> X1SQ=X1*X1;complex<double> X2SQ=X2*X2;
	complex<double> F1SQ=F1*F1;complex<double> F2SQ=F2*F2;
	complex<double> nu1SQ=nu1*nu1;complex<double> nu2SQ=nu2*nu2;
	complex<double> rhosSQ=rhos*rhos;

	complex<double> un(1.,0);
	complex<double> deux(2.,0);
	complex<double> trois(3.,0);
	
	complex<double> term1=alpha1*alpha1SQ*alpha1SQ*beta1*F2SQ*X1*(2.*F1*nu1*rho1 + alpha1*(-un + beta1)*X1)*X1SQ*X1SQ + 
		alpha2SQ*alpha2SQ*alpha2SQ*(-un + beta2)*beta2*F1SQ*X2SQ*X2SQ*X2SQ;

	complex<double> term2=-trois*alpha1SQ*alpha1SQ*alpha2*beta1*F2SQ*rhos*(deux*F1*nu1*rho1 + alpha1*(-un + beta1)*X1)*X1SQ*X1SQ*X2 + 
		alpha1*alpha2SQ*alpha2SQ*beta2*F1SQ*X1*(-6.*F2*nu2*rho2*rhos + 
		alpha1*(-un + deux*beta2 + (-deux + beta2)*rhosSQ)*X1)*X2SQ*X2SQ + 
		alpha2*alpha2SQ*alpha2SQ*beta2*F1SQ*(deux*F2*nu2*rho2 - trois*alpha1*(-un + beta2)*rhos*X1)*X2*X2SQ*X2SQ;
	
	complex<double> term3=alpha1SQ*alpha2*alpha2SQ*(deux*F1*F2*(-(beta1*F2*nu1*rho1*rhos) + 
        beta2*F1*(nu1*rhoc21*(-un + rhosSQ) + nu2*rho2*(deux+ rhosSQ))) + 
		alpha1*rhos*(-((-un + beta2)*beta2*F1SQ) - (-un + beta1)*beta1*F2SQ + deux*beta1*beta2*F1*F2*(-un + rhosSQ))*X1)*
		X1SQ*X2*X2SQ;
	
	complex<double> term4=alpha1SQ*alpha2SQ*F2*X1SQ*(-(F1SQ*F2*(-un + rhosSQ)*(nu1SQ + nu2SQ - deux*nu1*nu2*rhov)) + 
		deux*alpha1*F1*(-(beta2*F1*nu2*rho2*rhos) + beta1*F2*(nu2*rhoc12*(-un + rhosSQ) + nu1*rho1*(deux+ rhosSQ)))*X1 + 
		alpha1SQ*beta1*F2*(-un + deux*beta1 + (-deux + beta1)*rhosSQ)*X1SQ)*X2SQ;

	return (term1+term2+term3+term4)/(deux*F1SQ*F2SQ*alphas*alphas*alphas);
}

complex<double> BiSABRSpreadC(complex<double> F1,complex<double> F2,complex<double> alpha1,complex<double> alpha2,complex<double> beta1,complex<double> beta2,complex<double> rhos,complex<double> X1,complex<double> X2,complex<double> alphas,complex<double> s)
{
	complex<double> sF1beta;
	if(real(s+F2)>0)
	{
		sF1beta=pow(s+F2,beta1);
	}
	else
	{
		sF1beta=0.0;
	}
	return sqrt(sF1beta*sF1beta*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*sF1beta*X2*alpha1*alpha2*rhos)/alphas;
}

complex<double> BiSABRSpreadC_derivative(complex<double> F1,complex<double> F2,complex<double> alpha1,complex<double> alpha2,complex<double> beta1,complex<double> beta2,complex<double> rhos,complex<double> X1,complex<double> X2,complex<double> alphas,complex<double> s)
{
	complex<double> sF1beta;
	if(real(s+F2)>0)
	{
		sF1beta=pow(s+F2,beta1);
	}
	else
	{
		sF1beta=0.0;
	}
	complex<double> alphacurrent= sqrt(sF1beta*sF1beta*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*sF1beta*X2*alpha1*alpha2*rhos);
	return (sF1beta*sF1beta/(s+F2)*alpha1*alpha1*beta1-sF1beta/(s+F2)*X2*alpha1*alpha2*beta1*rhos)/(alphas*alphacurrent);
}

complex<double> BiSABRSpreadC_derivative2(complex<double> F1,complex<double> F2,complex<double> alpha1,complex<double> alpha2,complex<double> beta1,complex<double> beta2,complex<double> rhos,complex<double> X1,complex<double> X2,complex<double> alphas,complex<double> s)
{
	complex<double> sF1beta;
	if(real(s+F2)>0)
	{
		sF1beta=pow(s+F2,beta1);
	}
	else
	{
		sF1beta=0.0;
	}
	complex<double> alphacurrent= sqrt(sF1beta*sF1beta*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*sF1beta*X2*alpha1*alpha2*rhos);
	complex<double> element1=(2.*sF1beta*sF1beta/(s+F2)*alpha1*alpha1*beta1-2.*sF1beta/(s+F2)*X2*alpha1*alpha2*beta1*rhos);
	complex<double> element2=(2.*sF1beta*sF1beta/(s+F2)/(s+F2)*alpha1*alpha1*beta1*(-1.+2.*beta1)-2.*sF1beta/(s+F2)/(s+F2)*X2*alpha1*alpha2*beta1*rhos*(-1.+beta1));
	return (element1*element1)/(4.*alphacurrent*alphacurrent*alphacurrent)+element2/(2.*alphacurrent);
}

///  complex spreadoptions

complex<double> Complexified_BiSABR_SpreadOption(complex<double> F1,complex<double> alpha1,complex<double> beta1,complex<double> rho1,complex<double> nu1,complex<double> F2,complex<double> alpha2,complex<double> beta2,complex<double> rho2,complex<double> nu2,
					complex<double> K,complex<double> T,complex<double> rhos,complex<double> rhov,complex<double> rhoc12,complex<double> rhoc21)
{

	int SABRflag=ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL;
	if (abs(F2)<1e-10) return BlackSholes_Formula(real(F1),
							SABR_ComputeImpliedVol(real(F1), real(K),real(T),real(alpha1),real(beta1),real(rho1), real(nu1),SABRflag )*sqrt(real(T)),
							1.0,real(K),K_CALL);
	if (abs(F1)<1e-10) return BlackSholes_Formula(real(F2),
							SABR_ComputeImpliedVol(real(F2), -real(K),real(T),real(alpha2),real(beta2),real(rho2), real(nu2),SABRflag )*sqrt(real(T)),
							1.0,-real(K),K_PUT);


	complex<double> intrinsec,z,b1,b2,sqB0Balphaz,favg,zavg,Ck,alphasp,nusp,rhosp,musp,X1,X2,Csecond,Cfavg;
	X1=pow(F1,beta1);X2=pow(F2,beta2);alphasp=sqrt(X1*X1*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*X1*X2*alpha1*alpha2*rhos);
	nusp= nuBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	rhosp= rhoBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp,nusp);
	musp= muBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	z=BiSABRz( K, F1, F2, alpha1, alpha2, rhos, beta1, beta2);
	favg=(F1-F2+K)/2.;
	zavg=BiSABRz(favg, F1, F2, alpha1, alpha2, rhos, beta1, beta2);;
	b1= BiSABRSpreadC_derivative( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Csecond=BiSABRSpreadC_derivative2( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Cfavg=BiSABRSpreadC( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, favg);
	Ck=BiSABRSpreadC( F1, F2, alpha1, alpha2, beta1, beta2, rhos, X1, X2, alphasp, K);
	b2=Csecond*Cfavg+b1*b1;
	sqB0Balphaz=sqrt(Ck);
	if (real(F1-F2-K)>=0) 
	{
		intrinsec=F1-F2-K;
	}
	else
	{
		intrinsec=0;
	}

	return SABRgeneric6(intrinsec, z, b1, b2, zavg, sqB0Balphaz, T, musp, alphasp, rhosp, nusp);

}

complex<double> Complexified_BiSABR_SpreadOption_4(complex<double> F1,complex<double> alpha1,complex<double> beta1,complex<double> rho1,complex<double> nu1,complex<double> F2,complex<double> alpha2,complex<double> beta2,complex<double> rho2,complex<double> nu2,
					complex<double> K,complex<double> T,complex<double> rhos,complex<double> rhov,complex<double> rhoc12,complex<double> rhoc21)
{

	int SABRflag=ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL;
	if (abs(F2)<1e-10) return BlackSholes_Formula(real(F1),
							SABR_ComputeImpliedVol(real(F1), real(K),real(T),real(alpha1),real(beta1),real(rho1), real(nu1),SABRflag )*sqrt(real(T)),
							1.0,real(K),K_CALL);
	if (abs(F1)<1e-10) return BlackSholes_Formula(real(F2),
							SABR_ComputeImpliedVol(real(F2), -real(K),real(T),real(alpha2),real(beta2),real(rho2), real(nu2),SABRflag )*sqrt(real(T)),
							1.0,-real(K),K_PUT);


	complex<double> intrinsec,z,b1,b2,sqB0Balphaz,favg,zavg,Ck,alphasp,nusp,rhosp,musp,X1,X2,Csecond,Cfavg;
	X1=pow(F1,beta1);X2=pow(F2,beta2);alphasp=sqrt(X1*X1*alpha1*alpha1+X2*X2*alpha2*alpha2-2.*X1*X2*alpha1*alpha2*rhos);
	nusp= nuBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	rhosp= rhoBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp,nusp);
	musp= muBiSABR( F1, alpha1, beta1, rho1, nu1,F2, alpha2, beta2, rho2, nu2,rhos, rhov, rhoc12, rhoc21, X1, X2, alphasp);
	complex<double> deux(2.0,0);
	complex<double> kappaV=deux*nusp*alphasp;
	complex<double> lambdaV=(kappaV*kappaV/(deux*alphasp*alphasp)-(deux*musp+nusp*nusp));
	complex<double> thetaV=kappaV*kappaV/(lambdaV+lambdaV);
	double lambdaShift=0.1;
	complex<double> S0=F1-F2;
	double prec=1.e-8;
	double nbfirst=60,nb=60,NbStage=-1,NbOscill=0;
	return Complexified_NormalHeston(rhosp,lambdaV,thetaV,kappaV,alphasp*alphasp,S0,K,T,
		K_CALL,lambdaShift,nbfirst,nb,NbStage,NbOscill,prec);

}



double Complete_BiSABR_SpreadOption_ThroughComplex(complex<double> F1,complex<double> alpha1,complex<double> beta1,complex<double> rho1,complex<double> nu1,complex<double> F2,complex<double> alpha2,complex<double> beta2,complex<double> rho2,complex<double> nu2,
					complex<double> K,complex<double> T,int CallPut,complex<double> rhos,complex<double> rhov,complex<double> rhoc12,complex<double> rhoc21)
{
	complex<double> m1,m2;
	switch (CallPut)
		{
		case K_CALL:
			{				 
				m1=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, 
												 F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21);
				return real(m1);

				break;
			}
		case K_PUT:
			{
				m1=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, 
												 F2, alpha2, beta2, rho2, nu2,K, T, rhos, rhov, rhoc12, rhoc21);
				return real(K-F1+F2+m1);
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Complete_BiSABR_SpreadOption_ThroughComplex : callput  bad input :");
				break;
			}
		}
}



/// compute the SABR parameters associated with  S1/S2 and K/S2 for the change of numeraires

void BiSABREqivalents(double F1,double alpha1,double beta1,double rho1,double nu1,		
				double F2,double alpha2,double beta2,double rho2,double nu2,
				double rhos,double rhov,double rhoc12,double rhoc21	,double K,
				/* output */
				double& alphaZ,double& betaZ,double& rhoZ, double& nuZ,
				complex<double>& alphaY,double& betaY,double& rhoY, double& nuY,
				double& rhoYZ,double& rhoalphaYalphaZ,double& rhoalphaYZ,double& rhoYalphaZ
				)
{
	/// K should never be 0 !
	betaY=2.-beta1;
	complex<double> cK(K,0.0);
	complex<double> cbeta1(beta1-1.0,0.0);

		alphaY=pow(cK,cbeta1)*alpha1;

	nuY=nu1;
	rhoY=-rho1;
	betaZ=((2*pow(alpha1,2)*pow(F1,-4 + 2*beta1) - (-4 + 2*beta1)*pow(alpha1,2)*pow(F1,-6 + 2*beta1)*pow(F2,2) - 
		2*alpha1*alpha2*(1 + beta2)*rhos*pow(F1,-3 + beta1)*pow(F2,-1 + beta2) + 2*pow(alpha2,2)*pow(F1,-4)*pow(F2,2*beta2) + 
		2*alpha1*alpha2*(-3 + beta1)*rhos*pow(F1,-5 + beta1)*pow(F2,1 + beta2) + 
		2*beta2*pow(alpha2,2)*pow(F1,-2)*pow(F2,-2 + 2*beta2))*pow(pow(F1,-2) + pow(F2,-2),-1)*
		pow(pow(alpha1,2)*pow(F1,-4 + 2*beta1)*pow(F2,2) + pow(alpha2,2)*pow(F1,-2)*pow(F2,2*beta2) - 
		2*alpha1*alpha2*rhos*pow(F1,-3 + beta1)*pow(F2,1 + beta2),-1))/2.;
	alphaZ=pow(F2*pow(F1,-1),-betaZ)*pow(pow(alpha1,2)*pow(F1,-4 + 2*beta1)*pow(F2,2) + pow(alpha2,2)*pow(F1,-2)*pow(F2,2*beta2) - 
		2*alpha1*alpha2*rhos*pow(F1,-3 + beta1)*pow(F2,1 + beta2),0.5);
	rhoYZ=alpha1*F2*pow(F1,-2 + beta1)*pow(pow(alpha1,2)*pow(F1,-4 + 2*beta1)*pow(F2,2) + pow(alpha2,2)*pow(F1,-2)*pow(F2,2*beta2) - 
		2*alpha1*alpha2*rhos*pow(F1,-3 + beta1)*pow(F2,1 + beta2),-0.5) - 
		alpha2*rhos*pow(F1,-1)*pow(F2,beta2)*pow(pow(alpha1,2)*pow(F1,-4 + 2*beta1)*pow(F2,2) + 
		pow(alpha2,2)*pow(F1,-2)*pow(F2,2*beta2) - 2*alpha1*alpha2*rhos*pow(F1,-3 + beta1)*pow(F2,1 + beta2),-0.5);
	rhoalphaYZ=-(alpha1*F2*rhos*pow(F1,-2 + beta1)*pow(pow(alpha1,2)*pow(F1,-4 + 2*beta1)*pow(F2,2) + 
        pow(alpha2,2)*pow(F1,-2)*pow(F2,2*beta2) - 2*alpha1*alpha2*rhos*pow(F1,-3 + beta1)*pow(F2,1 + beta2),-0.5)) + 
		alpha2*pow(F1,-1)*pow(F2,beta2)*pow(pow(alpha1,2)*pow(F1,-4 + 2*beta1)*pow(F2,2) + 
		pow(alpha2,2)*pow(F1,-2)*pow(F2,2*beta2) - 2*alpha1*alpha2*rhos*pow(F1,-3 + beta1)*pow(F2,1 + beta2),-0.5);
	nuZ=pow(F2*pow(F1,-1),betaZ)*pow(pow(F1,-4)*(pow(alpha1,2)*pow(F1,2*beta1)*pow(F2,2) + pow(alpha2,2)*pow(F1,2)*pow(F2,2*beta2) - 
		2*alpha1*alpha2*rhos*pow(F1,1 + beta1)*pow(F2,1 + beta2)),-0.5)*
		pow(pow(F1,-6)*(pow(alpha1,6)*pow(-2 + beta1 + betaZ,2)*pow(F1,6*beta1)*pow(F2,4) + 
		2*(-2 + beta1 + betaZ)*pow(alpha1,5)*pow(F1,1 + 5*beta1)*pow(F2,3)*
        (F2*nu1*rho1 - alpha2*(-4 + beta1 + 3*betaZ)*rhos*pow(F2,beta2)) + 
		2*alpha1*pow(alpha2,3)*pow(F1,5 + beta1)*(-(nu2*rhos*(nu2 + nu1*rhov)*pow(F2,2)) - 
		(2 + beta2 - 3*betaZ)*(beta2 - betaZ)*rhos*pow(alpha2,2)*pow(F2,2*beta2) + 
		alpha2*((-1 + betaZ)*nu2*rhoc12 + (-1 - 2*beta2 + 3*betaZ)*nu2*rho2*rhos + (-beta2 + betaZ)*nu1*rhoc21*rhos)*
		pow(F2,1 + beta2))*pow(F2,-1 + 3*beta2) + 
		pow(alpha2,4)*pow(F1,6)*pow(F2,-2 + 4*beta2)*(pow(alpha2,2)*pow(beta2 - betaZ,2)*pow(F2,2*beta2) + 
		2*alpha2*(beta2 - betaZ)*nu2*rho2*pow(F2,1 + beta2) + pow(F2,2)*pow(nu2,2)) + 
		pow(alpha1,4)*pow(F1,2 + 4*beta1)*pow(F2,2)*(-2*alpha2*
		((-1 + betaZ)*nu1*rhoc21 + (-5 + 2*beta1 + 3*betaZ)*nu1*rho1*rhos + (-2 + beta1 + betaZ)*nu2*rhoc12*rhos)*
		pow(F2,1 + beta2) + pow(F2,2)*pow(nu1,2) + 
		pow(alpha2,2)*pow(F2,2*beta2)*((-1 + betaZ)*(-5 + 2*beta1 + 3*betaZ) + 
		(19 + 4*beta2 - 2*beta1*(5 + beta2 - 5*betaZ) - 2*(16 + beta2)*betaZ + pow(beta1,2) + 12*pow(betaZ,2))*pow(rhos,2)))\
        + pow(alpha1,2)*pow(alpha2,2)*pow(F1,4 + 2*beta1)*pow(F2,2*beta2)*
        (pow(alpha2,2)*pow(F2,2*beta2)*((-1 + betaZ)*(-1 - 2*beta2 + 3*betaZ) + 
		(3 - 2*beta2*(-5 + beta1 + 5*betaZ) + 2*betaZ*(-8 + beta1 + 6*betaZ) + pow(beta2,2))*pow(rhos,2)) + 
		pow(F2,2)*((pow(nu1,2) + pow(nu2,2))*pow(rhos,2) + 2*nu1*nu2*rhov*(1 + pow(rhos,2))) - 
		2*alpha2*pow(F2,1 + beta2)*((-1 + betaZ)*nu1*rho1*rhos + (-4 + beta1 + 3*betaZ)*nu2*rhoc12*rhos + 
		nu2*rho2*(-1 + betaZ - (1 + beta2 - 2*betaZ)*pow(rhos,2)) - 
		nu1*rhoc21*(beta2 - betaZ + (1 + beta2 - 2*betaZ)*pow(rhos,2)))) + 
		2*alpha2*pow(alpha1,3)*pow(F1,3 + 3*beta1)*pow(F2,1 + beta2)*
        (-(nu1*rhos*(nu1 + nu2*rhov)*pow(F2,2)) + rhos*pow(alpha2,2)*pow(F2,2*beta2)*
		(-5 + beta1 - 3*beta2 + beta1*beta2 + 12*betaZ - 2*beta1*betaZ + 2*beta2*betaZ - 6*pow(betaZ,2) + 
		(1 + beta2 - 2*betaZ)*(-3 + beta1 + 2*betaZ)*pow(rhos,2)) + 
		alpha2*pow(F2,1 + beta2)*((-1 + betaZ)*nu2*rho2*rhos - (2 + beta2 - 3*betaZ)*nu1*rhoc21*rhos + 
		nu1*rho1*(-1 + betaZ + (-3 + beta1 + 2*betaZ)*pow(rhos,2)) + 
		nu2*rhoc12*(-2 + beta1 + betaZ + (-3 + beta1 + 2*betaZ)*pow(rhos,2)))))*pow(F2*pow(F1,-1),-2*betaZ)*
		pow(pow(alpha1,2)*pow(F1,2*beta1)*pow(F2,2) + pow(alpha2,2)*pow(F1,2)*pow(F2,2*beta2) - 
		2*alpha1*alpha2*rhos*pow(F1,1 + beta1)*pow(F2,1 + beta2),-1),0.5);
	rhoYalphaZ=-(pow(alphaZ,-1)*pow(nuZ,-1)*(alpha2*nu2*rhoc12*pow(F1,-3)*pow(F2,beta2)*
        (-(alpha1*F2*rhos*pow(F1,beta1)) + alpha2*F1*pow(F2,beta2))*pow(F2*pow(F1,-1),-betaZ)*
        pow(pow(F1,-4)*(pow(alpha1,2)*pow(F1,2*beta1)*pow(F2,2) + pow(alpha2,2)*pow(F1,2)*pow(F2,2*beta2) - 
		2*alpha1*alpha2*rhos*pow(F1,1 + beta1)*pow(F2,1 + beta2)),-0.5) + 
		alpha1*F2*nu1*rho1*pow(F1,-4 + beta1)*(alpha1*F2*pow(F1,beta1) - alpha2*F1*rhos*pow(F2,beta2))*pow(F2*pow(F1,-1),-betaZ)*
        pow(pow(F1,-4)*(pow(alpha1,2)*pow(F1,2*beta1)*pow(F2,2) + pow(alpha2,2)*pow(F1,2)*pow(F2,2*beta2) - 
		2*alpha1*alpha2*rhos*pow(F1,1 + beta1)*pow(F2,1 + beta2)),-0.5) + 
		alpha2*rhos*pow(F1,-4)*pow(F2,-1 + beta2)*(-((-1 + betaZ)*pow(alpha1,2)*pow(F1,2*beta1)*pow(F2,2)) + 
		(beta2 - betaZ)*pow(alpha2,2)*pow(F1,2)*pow(F2,2*beta2) - 
		alpha1*alpha2*(1 + beta2 - 2*betaZ)*rhos*pow(F1,1 + beta1)*pow(F2,1 + beta2))*pow(F2*pow(F1,-1),-betaZ)*
        pow(pow(F1,-4)*(pow(alpha1,2)*pow(F1,2*beta1)*pow(F2,2) + pow(alpha2,2)*pow(F1,2)*pow(F2,2*beta2) - 
		2*alpha1*alpha2*rhos*pow(F1,1 + beta1)*pow(F2,1 + beta2)),-0.5) + 
		alpha1*pow(F1,-5 + beta1)*((-2 + beta1 + betaZ)*pow(alpha1,2)*pow(F1,2*beta1)*pow(F2,2) + 
		(-1 + betaZ)*pow(alpha2,2)*pow(F1,2)*pow(F2,2*beta2) - 
		alpha1*alpha2*(-3 + beta1 + 2*betaZ)*rhos*pow(F1,1 + beta1)*pow(F2,1 + beta2))*pow(F2*pow(F1,-1),-betaZ)*
        pow(pow(F1,-4)*(pow(alpha1,2)*pow(F1,2*beta1)*pow(F2,2) + pow(alpha2,2)*pow(F1,2)*pow(F2,2*beta2) - 
		2*alpha1*alpha2*rhos*pow(F1,1 + beta1)*pow(F2,1 + beta2)),-0.5)));
	rhoalphaYalphaZ=pow(alphaZ,-1)*pow(F1,-6)*pow(nuZ,-1)*((-2 + beta1 + betaZ)*rhos*pow(alpha1,3)*pow(F1,3*beta1)*pow(F2,3) + 
		pow(alpha2,2)*pow(F1,3)*(F2*nu2*rho2 + alpha2*(beta2 - betaZ)*pow(F2,beta2))*pow(F2,2*beta2) - 
		alpha1*alpha2*rhos*pow(F1,2 + beta1)*(F2*(nu2*rho2 + nu1*rhoc21) + alpha2*(2 + beta2 - 3*betaZ)*pow(F2,beta2))*
		pow(F2,1 + beta2) + pow(alpha1,2)*pow(F1,1 + 2*beta1)*pow(F2,2)*
		(F2*nu1*rhoc21 - alpha2*pow(F2,beta2)*(-1 + betaZ + (-3 + beta1)*pow(rhos,2) + 2*betaZ*pow(rhos,2))))*
		pow(F2*pow(F1,-1),-1 - betaZ)*pow(pow(F1,-4)*(pow(alpha1,2)*pow(F1,2*beta1)*pow(F2,2) + 
		pow(alpha2,2)*pow(F1,2)*pow(F2,2*beta2) - 2*alpha1*alpha2*rhos*pow(F1,1 + beta1)*pow(F2,1 + beta2)),-0.5);
	rhoZ=pow(alphaZ,-1)*pow(F1,-4)*pow(nuZ,-1)*(-((-2 + beta1 + betaZ)*pow(alpha1,4)*pow(F1,4*beta1)*pow(F2,4)) + 
		pow(alpha1,3)*pow(F1,1 + 3*beta1)*pow(F2,3)*(-(F2*nu1*rho1) + 2*alpha2*(-3 + beta1 + 2*betaZ)*rhos*pow(F2,beta2)) + 
		pow(alpha2,3)*pow(F1,4)*(F2*nu2*rho2 + alpha2*(beta2 - betaZ)*pow(F2,beta2))*pow(F2,3*beta2) - 
		alpha1*pow(alpha2,2)*pow(F1,3 + beta1)*(F2*(nu1*rhoc21*rhos + nu2*(rhoc12 + rho2*rhos)) + 
        2*alpha2*(1 + beta2 - 2*betaZ)*rhos*pow(F2,beta2))*pow(F2,1 + 2*beta2) - 
		alpha2*pow(alpha1,2)*pow(F1,2 + 2*beta1)*pow(F2,2 + beta2)*
		(-(F2*(nu2*rhoc12*rhos + nu1*(rhoc21 + rho1*rhos))) + 
        alpha2*pow(F2,beta2)*(-2 + (-4 + beta1 - beta2)*pow(rhos,2) + betaZ*(2 + 4*pow(rhos,2)))))*pow(F2*pow(F1,-1),-1 - betaZ)*
		pow(pow(alpha1,2)*pow(F1,2*beta1)*pow(F2,2) + pow(alpha2,2)*pow(F1,2)*pow(F2,2*beta2) - 
		2*alpha1*alpha2*rhos*pow(F1,1 + beta1)*pow(F2,1 + beta2),-1);
		
}

/////////////////////////////////////////////////////////////////////////////////////////
///
/// compute the necessary digital spread options  when F1 and F2 can be >0 or <0
///
/////////////////////////////////////////////////////////////////////////////////////////


double BiSABR_Digital_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
	double resultat;
	if ((F1>0)&&(F2>0))
	{
		double vplus,vmoins;
		switch (CallPut)
		{
		case K_CALL:
			{		
				double vplus,vmoins;
				vplus=BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= -(vplus-vmoins)/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
		case K_PUT:
			{
				vplus=BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= 1.0+(vplus-vmoins)/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : callput  bad input :");
				break;
			}
		}
	}
	else if (((F1<0)&&(F2>0))||((F1>0)&&(F2<0)))
	{
		complex<double> vplus,vmoins;
		switch (CallPut)
		{
		case K_CALL:
			{				 
				vplus=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= -real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
		case K_PUT:
			{
				vplus=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= 1.0+real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : callput  bad input :");
				break;
			}
		}
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : F1 or F2 =0 :");
	}
	
	if (resultat>1.0) return 1.0;
	if (resultat <0.0) return 0.0;
	return resultat;

}

/// use the corrected  BiSABR_SpreadOption2 for the long term smoothing 
/// by using SABRgeneric6
///  inside  SABRgeneric6 we add a factor (3.*I0*I0-1.) instead of 2 inside the kappa, making the kappa 
/// match the hagan one for vanilla sabr

double BiSABR_Digital_SpreadOption_2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
	double resultat;
	if ((F1>0)&&(F2>0))
	{
		double vplus,vmoins;
		switch (CallPut)
		{
		case K_CALL:
			{		
				double vplus,vmoins;
				vplus=BiSABR_SpreadOption2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=BiSABR_SpreadOption2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= -(vplus-vmoins)/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
		case K_PUT:
			{
				vplus=BiSABR_SpreadOption2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=BiSABR_SpreadOption2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= 1.0+(vplus-vmoins)/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : callput  bad input :");
				break;
			}
		}
	}
	else if (((F1<0)&&(F2>0))||((F1>0)&&(F2<0)))
	{
		complex<double> vplus,vmoins;
		switch (CallPut)
		{
		case K_CALL:
			{				 
				vplus=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= -real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
		case K_PUT:
			{
				vplus=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= 1.0+real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : callput  bad input :");
				break;
			}
		}
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : F1 or F2 =0 :");
	}
	if (resultat>1.0) return 1.0;
	if (resultat <0.0) return 0.0;
	return resultat;
	

}

/// use the corrected  BiSABR_SpreadOption3 for the long term smoothing 
/// by using SABRgeneric7
double BiSABR_Digital_SpreadOption_3(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
	double resultat;
	if ((F1>0)&&(F2>0))
	{
		double vplus,vmoins;
		switch (CallPut)
		{
		case K_CALL:
			{		
				double vplus,vmoins;
				vplus=BiSABR_SpreadOption3( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=BiSABR_SpreadOption3( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= -(vplus-vmoins)/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
		case K_PUT:
			{
				vplus=BiSABR_SpreadOption3( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=BiSABR_SpreadOption3( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= 1.0+(vplus-vmoins)/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : callput  bad input :");
				break;
			}
		}
	}
	else if (((F1<0)&&(F2>0))||((F1>0)&&(F2<0)))
	{
		complex<double> vplus,vmoins;
		switch (CallPut)
		{
		case K_CALL:
			{				 
				vplus=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= -real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
		case K_PUT:
			{
				vplus=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= 1.0+real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : callput  bad input :");
				break;
			}
		}
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : F1 or F2 =0 :");
	}
	if (resultat>1.0) return 1.0;
	if (resultat <0.0) return 0.0;
	return resultat;	
}


/// use the corrected  BiSABR_SpreadOption4 for the long term smoothing 
/// that uses the projection on a normal_heston process
///  
/// 
double BiSABR_Digital_SpreadOption_4(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
	double resultat;
	if ((F1>0)&&(F2>0))
	{
		double vplus,vmoins;
		switch (CallPut)
		{
		case K_CALL:
			{		
				double vplus,vmoins;
				vplus=BiSABR_SpreadOption4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=BiSABR_SpreadOption4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= -(vplus-vmoins)/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
		case K_PUT:
			{
				vplus=BiSABR_SpreadOption4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=BiSABR_SpreadOption4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= 1.0+(vplus-vmoins)/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : callput  bad input :");
				break;
			}
		}
	}
	else if (((F1<0)&&(F2>0))||((F1>0)&&(F2<0)))
	{
		complex<double> vplus,vmoins;
		switch (CallPut)
		{
		case K_CALL:
			{				 
				vplus=Complexified_BiSABR_SpreadOption_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=Complexified_BiSABR_SpreadOption_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= -real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
		case K_PUT:
			{
				vplus=Complexified_BiSABR_SpreadOption_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				vmoins=Complexified_BiSABR_SpreadOption_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
				resultat= 1.0+real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
				break;
			}
			
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption4 : callput  bad input :");
				break;
			}
		}
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : F1 or F2 =0 :");
	}
	if (resultat>1.0) return 1.0;
	if (resultat <0.0) return 0.0;
	return resultat;
	

}

/////////////////////////////////////////////////////////////////////////////////////////
///
/// compute the necessary digital spread options when alpha2 is complex
///
/////////////////////////////////////////////////////////////////////////////////////////
double BiSABR_Digital_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,complex<double> alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
	double resultat;
	complex<double> vplus,vmoins;
	switch (CallPut)
	{
	case K_CALL:
		{				 
			vplus=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
			vmoins=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
			resultat=  -real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
			break;
		}
	case K_PUT:
		{
			vplus=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
			vmoins=Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
			resultat=  1.0+real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
			break;
		}
		
	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : callput  bad input :");
			break;
		}
	}
	if (resultat>1.0) return 1.0;
	if (resultat <0.0) return 0.0;
	return resultat;

}

double BiSABR_Digital_SpreadOption_4(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,complex<double> alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
	double resultat;
	complex<double> vplus,vmoins;
	switch (CallPut)
	{
	case K_CALL:
		{				 
			vplus=Complexified_BiSABR_SpreadOption_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
			vmoins=Complexified_BiSABR_SpreadOption_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
			resultat=  -real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
			break;
		}
	case K_PUT:
		{
			vplus=Complexified_BiSABR_SpreadOption_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K+ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
			vmoins=Complexified_BiSABR_SpreadOption_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
				K-ARM_CF_K_SHIFT_FOR_DERIVATION/2.0, T, rhos, rhov, rhoc12, rhoc21);
			resultat=  1.0+real((vplus-vmoins))/ARM_CF_K_SHIFT_FOR_DERIVATION;
			break;
		}
		
	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiSABR_Digital_SpreadOption : callput  bad input :");
			break;
		}
	}
	if (resultat>1.0) return 1.0;
	if (resultat <0.0) return 0.0;
	return resultat;

}



/////////////////////////////////////////////////////////////////////////////////////////
///
/// compute the necessary digital spread option That Pays S1 :
///     use necessarly the complexified spreadoption pricing
///
/////////////////////////////////////////////////////////////////////////////////////////
double BiSABR_Digital_SpreadOption_PayS1(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
	double alphaZ,betaZ,rhoZ,nuZ,betaY,rhoY,nuY,rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ;
	complex<double> alphaY;
	double Keff;
	if(K==0) Keff=0.0000001; else Keff=K;

	BiSABREqivalents( F1, alpha1, beta1, rho1, nu1,		
				 F2, alpha2, beta2, rho2, nu2,
				 rhos, rhov, rhoc12, rhoc21	, -Keff,
				/* output */
				 alphaZ,betaZ,rhoZ,nuZ,alphaY,betaY,rhoY,nuY,rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ);

	return F1*(1.0- BiSABR_Digital_SpreadOption( F2/F1,  alphaZ,betaZ,rhoZ,nuZ, -Keff/F1, alphaY,betaY,rhoY,nuY,
					 1.0, T, CallPut, rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ));
}

double BiSABR_Digital_SpreadOption_PayS1_4(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
	double alphaZ,betaZ,rhoZ,nuZ,betaY,rhoY,nuY,rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ;
	complex<double> alphaY;
	double Keff;
	if(K==0) Keff=0.0000001; else Keff=K;

	BiSABREqivalents( F1, alpha1, beta1, rho1, nu1,		
				 F2, alpha2, beta2, rho2, nu2,
				 rhos, rhov, rhoc12, rhoc21	, -Keff,
				/* output */
				 alphaZ,betaZ,rhoZ,nuZ,alphaY,betaY,rhoY,nuY,rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ);

	return F1*(1.0- BiSABR_Digital_SpreadOption_4( F2/F1,  alphaZ,betaZ,rhoZ,nuZ, -Keff/F1, alphaY,betaY,rhoY,nuY,
					 1.0, T, CallPut, rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ));
} 
/// forward declaration !
double BiSABR_Digital_SpreadOption_PayS2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

/// Corrects the bias constated in the S1-S2-K reconstruction of the spreadoption

double Corrected_BiSABR_Digital_SpreadOption_PayS1(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{

		double paysS1,paysS2,spreadoption,newspreadoption,payK;

		paysS1=BiSABR_Digital_SpreadOption_PayS1( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
		paysS2=BiSABR_Digital_SpreadOption_PayS2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
		payK=K*BiSABR_Digital_SpreadOption_2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
		spreadoption=Complete_BiSABR_SpreadOption_2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
			switch (CallPut)
		{
		case K_CALL:
			{				 
				newspreadoption=paysS1-paysS2-payK;
				break;
			}
		case K_PUT:
			{
				newspreadoption=payK-paysS1+paysS2;
				break;
			}
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Complete_BiSABR_SpreadOption_ThroughComplex : callput  bad input :");
				break;
			}
		}
			return paysS1+ARM_CF_K_S1_S2_PARTITION_FACTOR*(spreadoption-newspreadoption);
}



/////////////////////////////////////////////////////////////////////////////////////////
///
/// compute the necessary digital spread option That Pays S2
///  use necessarly the complexified spreadoption pricing
/// 
/////////////////////////////////////////////////////////////////////////////////////////

double BiSABR_Digital_SpreadOption_PayS2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
		double alphaZ,betaZ,rhoZ,nuZ,betaY,rhoY,nuY,rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ;
			complex<double> alphaY;
		double Keff;
	if(K==0) Keff=0.0000001; else Keff=K;

	BiSABREqivalents( F2, alpha2, beta2, rho2, nu2,		
				 F1, alpha1, beta1, rho1, nu1,
				 rhos, rhov, rhoc21, rhoc12	, Keff,
				/* output */
				 alphaZ,betaZ,rhoZ,nuZ,alphaY,betaY,rhoY,nuY,rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ);

	return F2* BiSABR_Digital_SpreadOption( F1/F2,  alphaZ,betaZ,rhoZ,nuZ, Keff/F2, alphaY,betaY,rhoY,nuY,
					 1.0, T, CallPut, rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ);
}

double BiSABR_Digital_SpreadOption_PayS2_4(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{
		double alphaZ,betaZ,rhoZ,nuZ,betaY,rhoY,nuY,rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ;
			complex<double> alphaY;
		double Keff;
	if(K==0) Keff=0.0000001; else Keff=K;

	BiSABREqivalents( F2, alpha2, beta2, rho2, nu2,		
				 F1, alpha1, beta1, rho1, nu1,
				 rhos, rhov, rhoc21, rhoc12	, Keff,
				/* output */
				 alphaZ,betaZ,rhoZ,nuZ,alphaY,betaY,rhoY,nuY,rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ);

	return F2* BiSABR_Digital_SpreadOption_4( F1/F2,  alphaZ,betaZ,rhoZ,nuZ, Keff/F2, alphaY,betaY,rhoY,nuY,
					 1.0, T, CallPut, rhoYZ,rhoalphaYalphaZ,rhoalphaYZ,rhoYalphaZ);
}

double Corrected_BiSABR_Digital_SpreadOption_PayS2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21)
{

		double paysS1,paysS2,spreadoption,newspreadoption,payK;
		paysS1=BiSABR_Digital_SpreadOption_PayS1( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
		paysS2=BiSABR_Digital_SpreadOption_PayS2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
		payK=K*BiSABR_Digital_SpreadOption_2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
		spreadoption=Complete_BiSABR_SpreadOption_2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
					 K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
			switch (CallPut)
		{
		case K_CALL:
			{				 
				newspreadoption=paysS1-paysS2-payK;
				break;
			}
		case K_PUT:
			{
				newspreadoption=payK-paysS1+paysS2;
				break;
			}
		default:
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Complete_BiSABR_SpreadOption_ThroughComplex : callput  bad input :");
				break;
			}
		}
			return paysS2-(1.0-ARM_CF_K_S1_S2_PARTITION_FACTOR)*(spreadoption-newspreadoption);
}


/////////////////////////////////////////////////////////////////////////////////////////
///
/// compute the necessary digital spread option That Pays S3
///   use necessarly the complexified spreadoption pricing
///
/////////////////////////////////////////////////////////////////////////////////////////

/// we assume that S3 ~ a S1 +b S2 + c W + d
double BiSABR_Digital_SpreadOption_PayS3(
									   double F1,double alpha1,double beta1,double rho1,double nu1,
									   double F2,double alpha2,double beta2,double rho2,double nu2,
									   double rhos,double rhov,double rhoc12,double rhoc21,
									   double F3,
									   double sigma30,double rho13,double rho23,
									   double K,double T,int CallPut)
{
	double sigma3=sigma30*F3;
	double sigma1=pow(F1,beta1)*alpha1;
	double sigma2=pow(F2,beta2)*alpha2;
	double b= (rho23*sigma1-rho13*rhos*sigma2)*sigma3/(sigma2*(sigma2-rhos/sigma1));
	double a= (rho13*sigma3-b*rhos*sigma2)/sigma1;
	double payS1=BiSABR_Digital_SpreadOption_PayS1( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
		K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
	double payS2=BiSABR_Digital_SpreadOption_PayS2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
		K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
	double payD=       BiSABR_Digital_SpreadOption_2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
		K, T, CallPut, rhos, rhov, rhoc12, rhoc21);

	return a*payS1+b*payS2+(F3-a*F1-b*F2)*payD;



}


double BiSABR_Digital_SpreadOption_PayS3_4(
									   double F1,double alpha1,double beta1,double rho1,double nu1,
									   double F2,double alpha2,double beta2,double rho2,double nu2,
									   double rhos,double rhov,double rhoc12,double rhoc21,
									   double F3,
									   double sigma30,double rho13,double rho23,
									   double K,double T,int CallPut)
{
	double sigma3=sigma30*F3;
	double sigma1=SABR_ComputeImpliedVolAroundATM(F1,F1,T,alpha1,beta1,rho1, nu1,ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL );
	double sigma2=SABR_ComputeImpliedVolAroundATM(F2,F2,T,alpha2,beta2,rho2, nu2,ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL );
	/*
	double sigma1=pow(F1,beta1)*alpha1;
	double sigma2=pow(F2,beta2)*alpha2;
	double b= (rho23*sigma1-rho13*rhos*sigma2)*sigma3/(sigma2*(sigma2-rhos/sigma1));
	double a= (rho13*sigma3-b*rhos*sigma2)/sigma1;
	*/
	double unmoinsrhos2=(1.-rhos*rhos);
	double a=(rho13 - rhos*rho23)*sigma3/(sigma1*unmoinsrhos2);
	double b=(rhos*rho13-rho23)*sigma3/(sigma2*unmoinsrhos2);
	double payS1=BiSABR_Digital_SpreadOption_PayS1_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
		K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
	double payS2=BiSABR_Digital_SpreadOption_PayS2_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
		K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
	double payD=       BiSABR_Digital_SpreadOption_4( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
		K, T, CallPut, rhos, rhov, rhoc12, rhoc21);

	return a*payS1+b*payS2+(F3-a*F1-b*F2)*payD;



}


double Corrected_BiSABR_Digital_SpreadOption_PayS3(double F1,double alpha1,double beta1,double rho1,double nu1,
									   double F2,double alpha2,double beta2,double rho2,double nu2,
									   double rhos,double rhov,double rhoc12,double rhoc21,
									   double F3,
									   double sigma30,double rho13,double rho23,
									   double K,double T,int CallPut)									   
{
	double sigma3=sigma30*F3;
	double sigma1=pow(F1,beta1)*alpha1;
	double sigma2=pow(F2,beta2)*alpha2;
	double b= (rho23*sigma1-rho13*rhos*sigma2)*sigma3/(sigma2*(sigma2-rhos/sigma1));
	double a= (rho13*sigma3-b*rhos*sigma2)/sigma1;
	double payS1=Corrected_BiSABR_Digital_SpreadOption_PayS1( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
		K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
	double payS2=Corrected_BiSABR_Digital_SpreadOption_PayS2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
		K, T, CallPut, rhos, rhov, rhoc12, rhoc21);
	double payD=       BiSABR_Digital_SpreadOption_2( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
		K, T, CallPut, rhos, rhov, rhoc12, rhoc21);

	return a*payS1+b*payS2+(F3-a*F1-b*F2)*payD;

}



CC_END_NAMESPACE()

#undef ARM_CF_K_S1_S2_PARTITION_FACTOR
#undef  ARM_CF_K_SHIFT_FOR_DERIVATION
#undef  ARM_CF_NBSTEPS_FOR_CERF
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/