/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file gaussiananalytics.h
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2003
 */


#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"
#include "gpnumlib/normalinvcum.h"
#include <cmath>


CC_BEGIN_NAMESPACE( ARM )




////////////////////////////////////////////////////
///	Class  : ARM_GaussianAnalytics
///	Routine: dNormal
///	Returns: double
///	Action : returns the density of Gaussian distribution (mean 0, variance 1)
///				at input value d
////////////////////////////////////////////////////

double ARM_GaussianAnalytics::dNormal(double d)
{
	double result = exp(-0.5*d*d)*ARM_NumericConstants::ARM_INVSQRT2PI;
    
	return result;
}



////////////////////////////////////////////////////
///	Class  : ARM_GaussianAnalytics
///	Routine: cdfNormal
///	Action : inverse the cumulative function of a normal distribution
///				using the Moro algorithm(more precise than the standard one)
////////////////////////////////////////////////////

double ARM_GaussianAnalytics::cdfNormal(double d)
{
    /// truncation at 6 standard deviation
    if ( d < -6 )
    {
       return 0;
    }
    else if ( d > 6 )
    {
        return 1;
    }
    else
    {
       double x,y,z;
       double b1=0.2316419,b2=0.319381530,b3=-0.356563782;
       double b4=1.781477937,b5=-1.821255978,b6=1.330274429;
       double sqrtdeuxpi = 2.5066282746;
 
       if ( d >= 0 )
       {
          x=1/(1+b1*d);
          y=exp(-d*d/2)/sqrtdeuxpi;
          z=x*(b2+x*(b3+x*(b4+x*(b5+x*b6))));
          z=1-y*z;
       }
       else
       {
          x=1/(1-b1*d);
          y=exp(-d*d/2)/sqrtdeuxpi;
          z=x*(b2+x*(b3+x*(b4+x*(b5+x*b6))));
          z=y*z;
       }
 
       return(z);
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_GaussianAnalytics
///	Routine: cdfNormal2
///	Action : Calculation of  the Normale Cumulative
///		     Normal distribution probabilities accurate to 1.e-15. 
///		     Z = no. of standard deviations from the mean. 
///	
///		     Based upon algorithm 5666 for the error function, from: 
///		     Hart, J.F. et al, 'Computer Approximations', Wiley 1968 
///		     Programmer: Alan Miller 
///		     Latest revision - 30 March 1986 
///////////////////////////////////////////////////
double ARM_GaussianAnalytics::cdfNormal2(double d)
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


////////////////////////////////////////////////////
///	Class  : ARM_GaussianAnalytics
///	Routine: cdfNormal_Inv
///	Action : Calculation of  the Normale Cumulative inverse
///		     
///		precise algorithm (1e-12)
///////////////////////////////////////////////////
double ARM_GaussianAnalytics::cdfNormal_Inv(double d)
{
   /// return ARM_NumericConstants::ARM_SQRT_2*ARM_NormalInvCum::Inverse_erf_Moro(2.0*d-1.0);
	 return ARM_NormalInvCum::Inverse_erf_Moro(d);
}

CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----