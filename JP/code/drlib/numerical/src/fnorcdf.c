#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  ANORDF/DNORDF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 17, 1985

    Purpose:    Evaluate the standard normal (Gaussian) distribution
                function.

    Usage:      ANORDF(X)

    Arguments:
       X      - Argument for which the normal distribution function is to
                be evaluated.  (Input)
       ANORDF - Function value, the probability that a normal random
                variable takes a value less than or equal to X.  (Output)

    Keywords:   P-value; Probability integral; Error function;
                Probability distribution; Continuous random variables;
                Cumulative distribution function; CDF

    GAMS:       L5a1n; C8a

    Chapters:   STAT/LIBRARY Probability Distribution Functions and
                             Inverses
                SFUN/LIBRARY Probability Distribution Functions and
                             Inverses

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	SQR1D2	0.707106781186547524400844362105e0
#ifdef ANSI
Mfloat imsl_f_normal_cdf(Mfloat x)
#else
Mfloat imsl_f_normal_cdf(x)
	Mfloat           x;
#endif
{
	Mfloat           anordf_v;


	E1PSH("imsl_f_normal_cdf", "imsl_d_normal_cdf");
	/* Constant SQR1D2 = SQRT(1.0/2.0) */
	anordf_v = F_HALF * imsl_f_erfc(-x * SQR1D2);
	/*
	 * Eliminate underflow message from ERFC if it occurred.
	 */
	if (imsl_n1rty(0) == 2) {
		imsl_e1mes(0, 0, " ");
	}
	E1POP("imsl_f_normal_cdf", "imsl_d_normal_cdf");
	return (anordf_v);
}				/* end of function */
