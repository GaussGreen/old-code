/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file incompletebeta.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/port.h"

#include <cmath>
#include <complex>
#include <limits>

#include "gpclosedforms/incompletebeta.h"
#include "gpclosedforms/gamma.h"
#include "gpclosedforms/inverse.h"

#include "gpbase/numericconstant.h"

/// ARM Kernel
#include "glob/expt.h"   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-13
#define ARM_CF_MAXIT 2000



/////////////////////////////////////////////////////////////////////////
///
///    imcomplete beta function(x,a,b)
///     = 1/Beta(a,b)  Integral({t,0,x}, t^(a-1) (1-t)^(b-1) )
///
/////////////////////////////////////////////////////////////////////////

double BetaContinuousFraction(const double a, const double b, const double x)
{
		const double EPS=numeric_limits<double>::epsilon();
		const double FPMIN=numeric_limits<double>::min()/EPS;
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=ARM_CF_MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) <= EPS) break;
	}
	if (m > ARM_CF_MAXIT) 
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"a or b too big, or MAXIT too small in BetaContinuousFraction" );
	}
	return h;
}


double IncompleteBeta(const double& a, const double& b, const double& x)
{
	double bt;

	if (x < 0.0 || x > 1.0) 
	{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"Bad x in routine IncompleteBeta");
	}
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammalog(a+b)-gammalog(a)-gammalog(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*BetaContinuousFraction(a,b,x)/a;
	else
		return 1.0-bt*BetaContinuousFraction(b,a,1.0-x)/b;
}


double IncompleteBeta(const double& a, const double& b, const double& x0,const double& x1)
{
	return IncompleteBeta(a,b,x1)-IncompleteBeta(a, b,  x0);
}
/////////////////////////////////////////////////////////////////////////
///
///    incomplete beta inverse function (x,a,b) 
///     solution in y of IncompleteBeta(a,b,y)==x
///
/////////////////////////////////////////////////////////////////////////

double ApproxIncompleteBetaInverse(const double& a, const double& b, const double& x)
{
	double azbeta=a*x*IncompleteBeta(a,b,x);
	double f1= pow(azbeta,1.0/a)+(-1.+b)*pow(azbeta,2./a)/(1.+a)+(-1.+b)*(-4.-a+a*a+5.*b+3.*a*b)*pow(azbeta,3./a)/((1.+a)*(1.+a)*(2.+a));
	azbeta=a*x*IncompleteBeta(a,b,1.-x);
	double f2= pow(azbeta,1.0/a)+(-1.+b)*pow(azbeta,2./a)/(1.+a)+(-1.+b)*(-4.-a+a*a+5.*b+3.*a*b)*pow(azbeta,3./a)/((1.+a)*(1.+a)*(2.+a));
	f2=1.-f2;
	return (1.-x)*f1+x*f2;

}

double IncompleteBetaInverse_aux(const double& a, const double& b, const double& K)
{
		class DistributionToInverse : public DoubleToDoubleFunc 
	{
	public: 
		double a0;
		double b0;
		DistributionToInverse(double a1,double b1):
		a0(a1),b0(b1) {}
		
		virtual double operator() (double K0)  const
		{
			if(K0<=1e-9)
			{
				return 0.;
			}
				
			return IncompleteBeta(a0,b0,K0);
		}
	};
	double s=ApproxIncompleteBetaInverse(a,b,K);
	DistributionToInverse x(a,b);
	return Inverse(x,Inverse::BOUNDEDBY0AND1)(K,s,s/10.,1e-12);  // suppose that  Y is always postive
}

double IncompleteBetaInverse(const double& a, const double& b, const double& x)
{
	if (x==1) return 1;
	if (x==0) return 0;
	if(x>0.6)
	{
		return 1.-IncompleteBetaInverse_aux(b,a,1.-x);
	}
	else
	{
		return IncompleteBetaInverse_aux(a,b,x);
	}
}


double IncompleteBetaInverse(const double& a, const double& b,const double& y0, const double& K)
{
		return IncompleteBetaInverse( a, b, K+IncompleteBeta(a,b,y0));
}



CC_END_NAMESPACE()


 
#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/