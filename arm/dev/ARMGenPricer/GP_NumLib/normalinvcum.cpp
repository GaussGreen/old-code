/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file normalinvcum.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2003
 */


#include "gpnumlib/normalinvcum.h"
#include <cmath>


CC_BEGIN_NAMESPACE( ARM )

double ARM_NormalInvCum::HUGE_NB = 6.7;


////////////////////////////////////////////////////
///	Class  : ARM_NormalInvCum
///	Routine: Inverse_erf
///	Returns: double
///	Action : inverse the cumulative function of a normal distribution
////////////////////////////////////////////////////

double ARM_NormalInvCum::Inverse_erf(double x)
{
	/// handle huge number!
	if( x>1)
		return ARM_NormalInvCum::HUGE_NB;
	if( x<0)
		return -ARM_NormalInvCum::HUGE_NB;

	static double c0= 2.515517;
	static double c1= 0.802853;
	static double c2= 0.010328;
	static double d1= 1.432788;
	static double d2= 0.189269;
	static double d3= 0.001308;

	double t,sign;

	if(x>0.5) 
	{
	  sign= +1.0;
	  x= 1.0-x;
	} 
	else 
	{
	  sign= -1.0;
	}
	t=sqrt(-2.0*log(x));
	return sign*(t -((c2*t+c1)*t+c0)/(1.0+t*(d1+t*(d2+d3*t))));
}



////////////////////////////////////////////////////
///	Class  : ARM_NormalInvCum
///	Routine: Inverse_erf_Moro
///	Returns: double
///	Action : inverse the cumulative function of a normal distribution
///				using the Moro algorithm(more precise than the standard one)
/// Ref    : The Full Monte, by Boris Moro, Union Bank of Switzerland RISK 1995(2)
///				same algorithm as in the kernel!
////////////////////////////////////////////////////

double ARM_NormalInvCum::Inverse_erf_Moro(double x)    
{
    static double a[4] = {2.50662823884, -18.61500062529, 41.39119773534, 
                   -25.44106049637};

    static double b[4] = {-8.47351093090, 23.08336743743, -21.06224101826, 
                   3.13082909833};

    static double c[9] = {0.3374754822726147, 0.9761690190917186, 0.1607979714918209, 
                   0.0276438810333863, 0.0038405729373609, 0.0003951896511919,
                   0.0000321767881768, 0.0000002888167364, 0.0000003960315187};

	double	u = x-0.5,
			r;

	if (fabs(u)<0.42) 
	{
		r=u*u;
		r=u*(((a[3]*r+a[2])*r+a[1])*r+a[0])/((((b[3]*r+b[2])*r+b[1])*r+b[0])*r+1.0);
		return(r);
	}
	
	r=x;
	if(u>0.0) 
		r=1.0-x;
	r=log(-log(r));
	r=c[0]+r*(c[1]+r*(c[2]+r*(c[3]+r*(c[4]+r*(c[5]+r*c[6]+r*(c[7]+r*c[8]))))));
	if (u<0.0) 
		r=-r;
	return r;
}

CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----