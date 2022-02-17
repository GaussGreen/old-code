/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bivariate_normal.cpp
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
#include "gpnumlib/normalinvcum.h"
#include "gpclosedforms/lambert_function.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpbase/utilityport.h"

#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)


double AsymptoticNegativeNormalCDF(double x)
{
	if(x<-37.)
	{
		return 0.0;
	}
	else
	{
		return exp(-x*x/2.)/(ARM_NumericConstants::ARM_SQRT_PI*(-x/ARM_NumericConstants::ARM_SQRT_2+sqrt(x*x/2.+2.)));
	}
}




double NormalCDF(double x)
{	
	if(x>5.)
	{
		return 1.-AsymptoticNegativeNormalCDF(-x);
	}
	else if (x<-5.)
	{
		return AsymptoticNegativeNormalCDF(x);
	}
	else
	{
		return ARM_GaussianAnalytics::cdfNormal(x);
	}		
}

double NormalPDF(double x)
{
	return exp(-x*x/2.)/ARM_NumericConstants::ARM_SQRT_2_PI;
}


double AsymptoticNegativeNormalCDFInverse(double x)
{
	if(x<1e-10) return -37.;
	double z1=-sqrt(lambertfunc(1./(2.*ARM_NumericConstants::ARM_PI*x*x)));
	double z2=-sqrt(2.*fabs(log(-ARM_NumericConstants::ARM_SQRT_PI*x*(z1/ARM_NumericConstants::ARM_SQRT_2-sqrt(2.+z1*z1/2.)))));
	double z3=-sqrt(2.*fabs(log(-ARM_NumericConstants::ARM_SQRT_PI*x*(z2/ARM_NumericConstants::ARM_SQRT_2-sqrt(2.+z2*z2/2.)))));
	return -sqrt(2.*fabs(log(-ARM_NumericConstants::ARM_SQRT_PI*x*(z3/ARM_NumericConstants::ARM_SQRT_2-sqrt(2.+z3*z3/2.)))));
}


double NormalCDFInverse(double x)
{
	if(x>0.999999)
	{
		return -AsymptoticNegativeNormalCDFInverse(1.-x);
	}
	else if (x<0.000001)
	{
		return AsymptoticNegativeNormalCDFInverse(x);
	}
	else
	{
		return ARM_NormalInvCum::Inverse_erf_Moro(x);
	}

}


double AsymptoticNegativeNormalCDF2(double x)
{
	if(x<-37.)
	{
		return 0.0;
	}
	else
	{
		double x2=x*x;
		double x4=x2*x2;
		double x5=x4*x;
		double x9=x4*x4*x;
		return exp(-x*x/2.)/ARM_NumericConstants::ARM_SQRT_2_PI*(-10395.0/(x9*x4)+945.0/(x9*x2)-105.0/(x9)+15.0/(x5*x2)-3.0/x5+1.0/(x2*x)-1.0/x);
	}
}


////////////////////////////////////////////////////
///	
///		     Normal distribution probabilities accurate to 1.e-15. 
///		     Based upon algorithm 5666 for the error function, from: 
///		     Hart, J.F. et al, 'Computer Approximations', Wiley 1968 
///		     Programmer: Alan Miller 
///		     Latest revision - 30 March 1986 
///////////////////////////////////////////////////
double RegularNormalCDF2(double d)
{
   	/// System generated locals
    double ret_val, d__1;
	
    /// Local variables
    double zabs, p, expntl;
    zabs = fabs(d);
	
	///     |Z| > 37
	
    if (zabs > 37.) {
		p = 0.;
    } else {
		
		///     |Z| <= 37
		
		/// Computing 2nd power
		d__1 = zabs;
		expntl = exp(-(d__1 * d__1) / 2);
		
		///     |Z| < CUTOFF = 10/SQRT(2)
		if (zabs < 7.071067811865475) 
		{
			p = expntl * ((((((zabs * .03526249659989109 + .7003830644436881) 
				* zabs + 6.37396220353165) * zabs + 33.912866078383) * 
				zabs + 112.0792914978709) * zabs + 221.2135961699311) * 
				zabs + 220.2068679123761) / (((((((zabs * 
				.08838834764831844 + 1.755667163182642) * zabs + 
				16.06417757920695) * zabs + 86.78073220294608) * zabs + 
				296.5642487796737) * zabs + 637.3336333788311) * zabs + 
				793.8265125199484) * zabs + 440.4137358247522);
			
			///     |Z| >= CUTOFF
		} else {
			p = expntl / (zabs + 1 / (zabs + 2 / (zabs + 3 / (zabs + 4 / (
				zabs + .65))))) / 2.506628274631001;
		}
    }
    if (d > 0.) {
		p = 1 - p;
    }
    ret_val = p;
    return (ret_val);
}


double NormalCDF2(double x)
{	
	if(x>4.)
	{
		return 1.-AsymptoticNegativeNormalCDF2(-x);
	}
	else if (x<-4.)
	{
		return AsymptoticNegativeNormalCDF2(x);
	}
	else
	{
		return RegularNormalCDF2(x);
	}		
}

double NormalCDF3(double x)
{	
	if(x>4.)
	{
		return 1.-AsymptoticNegativeNormalCDF2(-x);
	}
	else if (x<-4.)
	{
		return AsymptoticNegativeNormalCDF2(x);
	}
	else if (x> 1)
	{
		return ARM_GaussianAnalytics::cdfNormal(x);
	}
	else if (x<-1)
	{
		return ARM_GaussianAnalytics::cdfNormal(x);
	}
	else 
	{
		return RegularNormalCDF2(x);
	}		
}

double NormalCDF(double x, double y, double rho,
				 double nGaussLegendre_1, double nGaussLegendre_2, double lowerBound, double upperBound)
{	
	/// switch x & y to avoid pbs if rho == 1
	/// si rho = 1, l'input de NormalCDF2 est plus ou moins l'infini
	/// quand on passe de l'un à l'autre, la normal cdf passe de 1 à 0
	/// le pb, c'est que l'endroit où ça change n'est pas dans les abscisses 
	/// d'intégration. Dans ce cas, on splitte l'integrale en deux pour supprimer
	/// ce biais de discrétisation de l'intégrale.
	///


	y = CC_Min(y, upperBound);
	if(x<lowerBound) return 0.;
	if(y<lowerBound) return 0.;

	if ( rho != 0 && (x / rho < y) )
	{
		int i; 
		double sum = 0.0;
		double sqr1mrho = sqrt(1.0-rho*rho);

		GaussLegendre_Coefficients c0(nGaussLegendre_1);

		double y0 = x / rho;
		if(y0 > lowerBound)
		{
			GaussStratifiedHermiteLegendre_Coefficients c1(&c0,nGaussLegendre_2, lowerBound, y0);
			int n1 = c1.get_order();
			
			for(i=0;i<n1;i++)
			{
				sum+=NormalCDF2((x-c1.get_point(i)*rho)/sqr1mrho)*c1.get_weight(i);
			}
		}
		else
			y0 = lowerBound;

		GaussStratifiedHermiteLegendre_Coefficients c2(&c0,nGaussLegendre_2, y0, y);
		int n2 = c2.get_order();
					
		for(i=0;i<n2;i++)
		{
			sum+=NormalCDF2((x-c2.get_point(i)*rho)/sqr1mrho)*c2.get_weight(i);
		}
		return sum;
	}
	else
	{
		GaussLegendre_Coefficients c0(nGaussLegendre_1);
		GaussStratifiedHermiteLegendre_Coefficients c(&c0,nGaussLegendre_2,lowerBound,y);
		int i;int n = c.get_order();
		double sum=0.0;double sqr1mrho=sqrt(1.0-rho*rho);
		for(i=0;i<n;i++)
		{
			sum+=NormalCDF2((x-c.get_point(i)*rho)/sqr1mrho)*c.get_weight(i);
		}
		return sum;
	}
}

/// Integral of x * Bivariatedensity(x,y,rho)   up to X,Y
double Normal_X_Expectation(double X, double y, double rho,
				 double nGaussLegendre_1, double nGaussLegendre_2, double lowerBound, double upperBound)

{
	X = CC_Min(X, upperBound);
	if(X<lowerBound) return 0.;
	if(y<lowerBound) return 0.;
	GaussLegendre_Coefficients c0(nGaussLegendre_1);
	GaussStratifiedHermiteLegendre_Coefficients c(&c0,nGaussLegendre_2,lowerBound,X);
	int i;int n=c.get_order();
	double sum=0.0;double sqr1mrho=sqrt(1.0-rho*rho);double x;
	for(i=0;i<n;i++)
	{
		x=c.get_point(i);
		sum+=x*NormalCDF2((y-x*rho)/sqr1mrho)*c.get_weight(i);
	}
	return sum;
}

double Normal_XX_Expectation(double X, double y, double rho,
				 double nGaussLegendre_1, double nGaussLegendre_2, double lowerBound, double upperBound)

{
	X = CC_Min(X, upperBound);
	if(X<lowerBound) return 0.;
	if(y<lowerBound) return 0.;
	GaussLegendre_Coefficients c0(nGaussLegendre_1);
	GaussStratifiedHermiteLegendre_Coefficients c(&c0,nGaussLegendre_2,lowerBound,X);
	int i;int n=c.get_order();
	double sum=0.0;double sqr1mrho=sqrt(1.0-rho*rho);double x;
	for(i=0;i<n;i++)
	{
		x=c.get_point(i);
		sum+=x*x*NormalCDF2((y-x*rho)/sqr1mrho)*c.get_weight(i);
	}
	return sum;
}

/// Integral of y * Bivariatedensity(x,y,rho)   up to X,Y
double Normal_Y_Expectation(double x, double Y, double rho,
				 double nGaussLegendre_1, double nGaussLegendre_2, double lowerBound, double upperBound)

{
	Y = CC_Min(Y, upperBound);
	if(x<lowerBound) return 0.;
	if(Y<lowerBound) return 0.;
	GaussLegendre_Coefficients c0(nGaussLegendre_1);
	GaussStratifiedHermiteLegendre_Coefficients c(&c0,nGaussLegendre_2,lowerBound,Y);
	int i;int n=c.get_order();
	double sum=0.0;double sqr1mrho=sqrt(1.0-rho*rho);
	for(i=0;i<n;i++)
	{
		double y=c.get_point(i);
		sum+=y*NormalCDF2((x-y*rho)/sqr1mrho)*c.get_weight(i);
	}
	return sum;
}

double Normal_YY_Expectation(double x, double Y, double rho,
				 double nGaussLegendre_1, double nGaussLegendre_2, double lowerBound, double upperBound)

{
	Y = CC_Min(Y, upperBound);
	if(x<lowerBound) return 0.;
	if(Y<lowerBound) return 0.;
	GaussLegendre_Coefficients c0(nGaussLegendre_1);
	GaussStratifiedHermiteLegendre_Coefficients c(&c0,nGaussLegendre_2,lowerBound,Y);
	int i;int n=c.get_order();
	double sum=0.0;double sqr1mrho=sqrt(1.0-rho*rho);
	for(i=0;i<n;i++)
	{
		double y=c.get_point(i);
		sum+=y*y*NormalCDF2((x-y*rho)/sqr1mrho)*c.get_weight(i);
	}
	return sum;
}



/// Integral of x*y * Bivariatedensity(x,y,rho)   up to X,Y
double Normal_XY_Expectation(double x, double y, double rho,
				 double nGaussLegendre_1, double nGaussLegendre_2, double lowerBound, double upperBound)

{
	if(x<lowerBound) return 0.;
	if(y<lowerBound) return 0.;
	double unmrho=1-rho*rho;
	double sunmrho=sqrt(unmrho);
	double u=(y-rho*x)/sunmrho;
	return exp(-x*x/2.)*unmrho/(ARM_NumericConstants::ARM_SQRT_2_PI)*
		(exp(-u*u/2.)*sunmrho/ARM_NumericConstants::ARM_SQRT_2_PI-x*rho*NormalCDF2(u))+
		rho*Normal_YY_Expectation(x,y,rho,nGaussLegendre_1,nGaussLegendre_2,lowerBound,upperBound);
}

/// Integral of the bivariate with respect to y (derivative of the biavriate cummulative w.r. to x)
double DxBivariateCummulative2D(double x , double y, double rho)
{
	return exp(-x*x/2.)*NormalCDF2((y-x*rho)/sqrt(1.0-rho*rho))/(ARM_NumericConstants::ARM_SQRT_2_PI);
}

/// Integral of the bivariate with respect to x (derivative of the biavriate cummulative w.r. to y)
double DyBivariateCummulative2D(double x , double y, double rho)
{
	return exp(-y*y/2.)*NormalCDF2((x-y*rho)/sqrt(1.0-rho*rho))/(ARM_NumericConstants::ARM_SQRT_2_PI);
}		
		
		/// Derivative of the gaussian 2-copula with respect to x

double DxBivariateCopula(double x,double y, double rho)
{
	if(x<0.00001)
	{
		return 0.0;
	}
	else
	{
		double invx=NormalCDFInverse(x);
		double invy=NormalCDFInverse(y);
		return DxBivariateCummulative2D(invx,invy,rho)/NormalPDF(invx);
	}
	
	
}
		/// Derivative of the gaussian 2-copula with respect to y

double DyBivariateCopula(double x,double y, double rho)
{
	if(y<0.00001)
	{
		return 0.0;
	}
	else
	{
		double invx=NormalCDFInverse(x);
		double invy=NormalCDFInverse(y);
		return DyBivariateCummulative2D(invx,invy,rho)/NormalPDF(invy);
	}
	
	
}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/