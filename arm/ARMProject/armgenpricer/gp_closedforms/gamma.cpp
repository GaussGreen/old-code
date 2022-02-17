/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file gamma.cpp
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
#include <complex>

#include "gpclosedforms/gamma.h"

#include "gpbase/numericconstant.h"

#include "expt.h"   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-13
#define ARM_CF_MAXIT 2000

#define G_ITMAX 100
#define G_FPMIN 1.0e-30
#define G_EPS 1.2e-7



/////////////////////////////////////////////////////////////////////////
///
///    Log of Gamma(x)  ( Gamma(x) = Factorial(x-1) for x integer )
///
/////////////////////////////////////////////////////////////////////////

long double factorial(int n)
{
	
	if (n>18) 
	{
		long_double r=gammaLD(1.+n);
		return gammaLD(1.+n).todouble();
	}
	switch (n)
	{
	case 0:
		return 1;
	case 1:
		return 1;
	case 2:
		return 2;
	case 3:
		return 6;
	case 4:
		return 24;
	case 5:
		return 120;
	case 6:
		return 720;
	case 7:
		return 5040;
	case 8:
		return 40320;
	case 9:
		return 362880;
	case 10:
		return 3628800;
	case 11:
		return 39916800;
	case 12:
		return 479001600;
	case 13:
		return 6227020800;
	case 14:
		return 87178291200;
	case 15:
		return 1307674368000;
	case 16:
		return 20922789888000;
	case 17:
		return 355687428096000;
	case 18:
		return 6402373705728000;
	default:
		return 0;
		
	}
}


/// assume that xx is >1
double  gammalog(const double& xx)
{
	int j;
	long double x,y,tmp,ser;
	static const double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
		-0.5395239384953e-5};

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<6;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

long_double  gammaLD(const double& xx)
{
	int sign=1;
	double loga=0;
	if(xx<1.)
	{
	double z=ARM_NumericConstants::ARM_PI/sin(ARM_NumericConstants::ARM_PI*(1.-xx));
	loga=log(fabs(z))-gammalog(1.-xx);
	if (z<0) sign=-1;	
	}
	else
	{
		loga=gammalog(xx);
	}
	return long_double(loga,sign);
	
}

double  gamma(const double& xx)
{
	int sign=1;
	double loga=0;
	if(xx<1.)
	{
		double z=ARM_NumericConstants::ARM_PI/sin(ARM_NumericConstants::ARM_PI*(1.-xx));
		loga=log(fabs(z))-gammalog(1.-xx);
		if (z<0) sign=-1;	
	}
	else
	{
		loga=gammalog(xx);
	}
	return sign*exp(loga);
	
}

/////////////////////////////// Complex Gamma /////////////////////////////////////

/// Lanczos approximation  (NR p218)

complex<double> GammaLog(complex<double> z)
{
	int j;
	complex<double> un(1.,0);
	complex<double> x,y,tmp,ser;
	static const complex<double> cof[6]={
		76.18009172947146,
		-86.50532032941677,
		24.01409824083091,
		-1.231739572450155,
		0.1208650973866179e-2,
		-0.5395239384953e-5
	};

	y=x=z;
	tmp=x+5.5;
	tmp -= (x+0.5)*std::log(tmp);
	ser=1.000000000190015;
	
	for (j=0;j<6;j++) {
		y+=un;
		ser += cof[j]/y;
	}
	return -tmp+std::log(2.5066282746310005*ser/x);
}




complex<double> Gamma(complex<double> z)
{	
	if(real(z)<1)
	{
		complex<double> un(1.,0);
		complex<double> u=un-z;
		complex<double> pi(ARM_NumericConstants::ARM_PI,0);
		complex<double> x=GammaLog(u);
		complex<double> y=sin(pi*u);
		return std::exp(-x)*pi/y;
	}
	else
	{
		complex<double> x=GammaLog(z);
		return std::exp(x);
	}
}


//Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser.
void gser(double *gamser, double a, double x, double *gln)
{
	int n;
	double sum,del,ap;
	*gln=gammalog(a);
	if (x <= 0.0){
		if (x < 0.0) return;
		*gamser=0.0;
		return;
	} 
	else{
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=G_ITMAX;n++){
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*G_EPS){
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		return;
	}
}

//Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation
void gcf(double *gammcf, double a, double x, double *gln)
{
	int i;
	double an,b,c,d,del,h;
	*gln=gammalog(a);
	b=x+1.0-a; 
	c=1.0/G_FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=G_ITMAX;i++){ 
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < G_FPMIN) d=G_FPMIN;
		c=b+an/c;
		if (fabs(c) < G_FPMIN) c=G_FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < G_EPS) break;
	}
	if (i > G_ITMAX) return;
	*gammcf=exp(-x+a*log(x)-(*gln))*h; 
}

// Returns the incomplete gamma function P(a, x).
// = (1/Gamm(a))*Int[0,x]{e-t x t^a-1 dt}
double gammp(double a, double x)
{
	double gamser,gammcf,gln;
	if (x < 0.0 || a <= 0.0) return 0.;
	if (x < (a+1.0)){ 
		gser(&gamser,a,x,&gln);
		return gamser;
	} 
	else{ 
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf; 
	}
}


CC_END_NAMESPACE()


 
#undef ARM_CF_EPS
#undef ARM_CF_MAXIT
#undef G_ITMAX
#undef G_EPS
#undef G_FPMIN


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/