/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Complex Numbers
 * Revision 1.1  2003/10/08 16:45:06  atriki
 * Initial revision
 *
 *
 *
 */

#ifndef _CPLX_H_
#define _CPLX_H_

#include "gpbase/port.h"

#include "gpbase/numericconstant.h"
#include <math.h>

CC_BEGIN_NAMESPACE( ARM )

class ARM_Cplx
{
// -- static direct functions (much faster....)
public:
	inline static void inverse(double re, double im, /* result-> */ double& re_inv, double& im_inv);
	inline static void exp(double re, double im, /* result-> */ double& re_exp, double& im_exp);
	inline static void sqrt(double re, double im, /* result-> */ double& re_sqrt, double& im_sqrt);

	inline static void log(double re, double im, /* result-> */ double& re_log, double& im_log);

	inline static void sum(double re1, double im1, double re2, double im2, /* result-> */ double& re_sum,  double& im_sum);
	inline static void difference(double re1, double im1, double re2, double im2, /* result-> */ double& re_diff, double& im_diff);
	inline static void product(double re1, double im1, double re2, double im2, /* result-> */ double& re_prod, double& im_prod);
	inline static void quotient(double re1, double im1, double re2, double im2, /* result-> */ double& re_quot, double& im_quot);
	
	// for robustness / uses STL. Much slower...
	static void robust_quotient(double re1, double im1, double re2, double im2, /* result-> */ double& re_quot, double& im_quot);



};

////////////////////////////////////////////////////
///	Class  : ARM_ARM_Cplx
///	Routine: inverse
///	Returns: void
///	Action : compute 1/u
////////////////////////////////////////////////////
inline void ARM_Cplx::inverse(double re, double im, /* result-> */ double& re_inv, double& im_inv)
{
	double sqrmod = re*re + im*im;

#if defined(__GP_STRICT_VALIDATION)
	if(sqrmod<ARM_NumericConstants::ARM_LOWEST_POSITIVE_NUMBER)
		/// Non representable number in double format
		ARM_THROW( ERR_INVALID_ARGUMENT, "Trying to inverse a null number !" );
#endif

	re_inv =  re/sqrmod;
	im_inv = -im/sqrmod;
}

////////////////////////////////////////////////////
///	Class  : ARM_ARM_Cplx
///	Routine: exp
///	Returns: void
///	Action : compute exp(u)
////////////////////////////////////////////////////
inline void ARM_Cplx::exp (double re, double im, /* result-> */ double& re_exp, double& im_exp)
{
	double exp_re = ::exp(re);
	re_exp = exp_re * cos (im);
	im_exp = exp_re * sin (im);
}

////////////////////////////////////////////////////
///	Class  : ARM_ARM_Cplx
///	Routine: sqrt
///	Returns: void
///	Action : compute sqrt(u)
////////////////////////////////////////////////////
inline void ARM_Cplx::sqrt (double re, double im, /* result-> */ double& re_sqrt, double& im_sqrt)
{
	if (im!=0)
	{
		re_sqrt = ::sqrt ( 0.5 * (re + ::sqrt(re * re + im * im) ) ) ;
		im_sqrt = 0.5 * im / re_sqrt ;
	}
	else	
	{	if (re>0)	
		{	im_sqrt = 0;
			re_sqrt = ::sqrt (re);
		}
		else
		{	re_sqrt = 0;
			im_sqrt = ::sqrt(-re);
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ARM_Cplx
///	Routine: sum
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
inline void ARM_Cplx::sum (double re1, double im1, double re2, double im2, /* result-> */ double& re_sum,  double& im_sum)
{
	re_sum = re1 + re2;
	im_sum = im1 + im2;
}

////////////////////////////////////////////////////
///	Class  : ARM_ARM_Cplx
///	Routine: difference
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
inline void ARM_Cplx::difference (double re1, double im1, double re2, double im2, /* result-> */ double& re_sum,  double& im_sum)
{
	re_sum = re1 - re2;
	im_sum = im1 - im2;
}

////////////////////////////////////////////////////
///	Class  : ARM_ARM_Cplx
///	Routine: product
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
inline void ARM_Cplx::product (double re1, double im1, double re2, double im2, /* result-> */ double& re_prod,  double& im_prod)
{
	re_prod = re1 * re2 - im1 * im2;
	im_prod = re1 * im2 + im1 * re2;
}

////////////////////////////////////////////////////
///	Class  : ARM_ARM_Cplx
///	Routine: quotient
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
inline void ARM_Cplx::quotient (double re1, double im1, double re2, double im2, /* result-> */ double& re_quot,  double& im_quot)
{
    double sqrmod2,absRe2=fabs(re2),absIm2=fabs(im2);
    double x=1.0,y=1.0,maxAbsReIm2 = (absRe2 > absIm2 ? absRe2 : absIm2);
    if(maxAbsReIm2<ARM_NumericConstants::ARM_SQRT_LOWEST_POSITIVE_NUMBER)
    {
        if(absIm2>=absRe2*1.0e+8)
		{
			x = re2/im2;
            sqrmod2 = im2; // mantisse will loose re^2 then drop it (re could not be represented then = 0.0)
		}
        else
		{
			y = im2/re2;
			if(absRe2>=absIm2*1.0e+8)
				sqrmod2 = re2; // mantisse will loose im^2 then drop it (im could not be represented then = 0.0)
			else
				sqrmod2 = re2*(1.0+y*y);
        }

#if defined(__GP_STRICT_VALIDATION)
	if(sqrmod2<ARM_NumericConstants::ARM_LOWEST_POSITIVE_NUMBER)
		/// Non representable number in double format
		ARM_THROW( ERR_INVALID_ARGUMENT, "Trying to inverse a null number !" );
#endif

		re_quot = (re1*x + im1*y) / sqrmod2;
		im_quot = (im1*x - re1*y) / sqrmod2;
	}
	else
	{
		sqrmod2	= re2*re2 + im2*im2;
		re_quot =  (re1 * re2 + im1 * im2) / sqrmod2;
		im_quot =  (im1 * re2 - re1 * im2) / sqrmod2;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ARM_Cplx
///	Routine: log
///	Returns: void
///	Action : ln(u)
////////////////////////////////////////////////////
inline void ARM_Cplx::log (double re, double im, /* result-> */ double& re_log, double& im_log)
{
    double mod,absRe=fabs(re),absIm=fabs(im);
    double maxAbsReIm = (absRe > absIm ? absRe : absIm);
    if(maxAbsReIm<ARM_NumericConstants::ARM_SQRT_LOWEST_POSITIVE_NUMBER)
    {
        if(absIm>=absRe*1.0e+8)
            mod = absIm; // mantisse will loose re^2 then drop it (re could not be represented then = 0.0)
        else if(absRe>=absIm*1.0e+8)
            mod = absRe; // mantisse will loose im^2 then drop it (im could not be represented then = 0.0)
        else
        {
            double x=im/re;
            mod = absRe * ::sqrt(1.0+x*x);
        }
        if(mod == 0.0)
        {
            /// re & im are not representable
            re_log = ::log(ARM_NumericConstants::ARM_LOWEST_POSITIVE_NUMBER); 
	        im_log = 0.0 ;
            return;
        }
    }
    else
	    mod = ::sqrt(re*re + im*im);

	double asinIm = asin(im/mod);
	double arg = (re>0) ? asinIm : (im>0 ? ARM_NumericConstants::ARM_PI : -ARM_NumericConstants::ARM_PI) - asinIm;

	re_log = ::log(mod); 
	im_log = arg;
}

CC_END_NAMESPACE()

#endif //_CPLX_H_
