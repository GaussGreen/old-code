/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file erf.cpp
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
#include <complex>

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"


#include "gpclosedforms/long_double.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/gamma.h"

#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-13
#define ARM_CF_MAXIT 2000


inline double SIGN(const double &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;} 


////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of real error function erf(x) 
///					erf(t) =   2/Sqrt[Pi]  Integrate[Exp[-x^2], {x, 0, t}]
///
////////////////////////////////////////////////////////////////////////////////////////////

double erf(double x)
 
{
	return 2.*ARM_GaussianAnalytics::cdfNormal(ARM_NumericConstants::ARM_SQRT_2*x)-1.;
}


////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of imaginary error function erfi(x) 
///					erfi(t) = erf(i t) / i  =   2/Sqrt[Pi]  Integrate[Exp[x^2], {x, 0, t}]
///
////////////////////////////////////////////////////////////////////////////////////////////


double erfi(const double x)
{
	const int NMAX=6;
	const double H=0.4, A1=2.0/3.0, A2=0.4, A3=2.0/7.0;
	int i,n0;
	static bool init = true;
	double d1,d2,e1,e2,sum,x2,xp,xx,ans,q;
	static vector<double> c(NMAX);

	if (init) {
		init=false;
		for (i=0;i<NMAX;i++) 
		{
			q=(2.0*i+1.0)*H;
			c[i]=exp(-q*q);
		}
	}
	if (fabs(x) < 0.2) {
		x2=x*x;
		ans=x*(1.0-A1*x2*(1.0-A2*x2*(1.0-A3*x2)));
	} else {
		xx=fabs(x);
		n0=2*int(0.5*xx/H+0.5);
		xp=xx-n0*H;
		e1=exp(2.0*xp*H);
		e2=e1*e1;
		d1=n0+1;
		d2=d1-2.0;
		sum=0.0;
		for (i=0;i<NMAX;i++,d1+=2.0,d2-=2.0,e1*=e2)
			sum += c[i]*(e1/d1+1.0/(d2*e1));

		ans=0.56418958354775629*SIGN(exp(-xp*xp),x)*sum;
	}
	return 2*ans*exp(x*x)/ARM_NumericConstants::ARM_SQRT_PI;
	//return ans
		;
}

////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of complex error function erf(z) 
///					
////////////////////////////////////////////////////////////////////////////////////////////

complex<double> cerf(complex<double> z,int nbsteps)
{
	double x=z.real();
	double y=z.imag();
	if(fabs(x)<1e-12)
	{
		return complex<double>(0,erfi(y));
	}
	else
	{
		double facn,fn,gn;
		double h=exp(-x*x)/ARM_NumericConstants::ARM_PI;
		double x1=erf(x)+h*(1-cos(2.*x*y))/(2.*x);
		double y1=h*sin(2.*x*y)/(2.*x);
		for(int i=1;i<nbsteps;i++)
		{
			facn=2.*exp(-i*i/4.)*h/(i*i+4.*x*x);
			fn=2.*x*(1.-cosh(y*i)*cos(2.*x*y))+sinh(y*i)*i*sin(2.*x*y);
			gn=2.*x*cosh(y*i)*sin(2.*x*y)+sinh(y*i)*i*cos(2.*x*y);
			x1+=facn*fn;
			y1+=facn*gn;
		}
		return complex<double>(x1,y1);
	}
}


complex<double> cnormal(complex<double> z,int nbsteps)
{
	complex<double> SQRT_2(ARM_NumericConstants::ARM_SQRT_2,0);
	complex<double> un(1.,0);
	complex<double> deux(2.,0);
	return (cerf(z/SQRT_2,nbsteps)+un)/deux;
}


////////////////////////////////////////////////////////////////////////////////////////////
///
///			Exportation of	computation of complex error function erf(z) 
///					
////////////////////////////////////////////////////////////////////////////////////////////

double Export_RealPart_ComplexErf(double realpart, double imaginarypart, int nbterm)
{
	complex<double> z(realpart,imaginarypart);
	complex<double> result=cerf(z,nbterm);
	return real(result);
}

double Export_ImaginaryPart_ComplexErf(double realpart, double imaginarypart, int nbterm)
{
	complex<double> z(realpart,imaginarypart);
	complex<double> result=cerf(z,nbterm);
	return imag(result);
}



CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
