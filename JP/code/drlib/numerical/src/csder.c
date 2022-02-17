


























#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  CSDER/DCSDER (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 15, 1988

    Purpose:    Evaluate the derivative of a cubic spline.

    Usage:      CSDER(IDERIV, X, NINTV, BREAK, CSCOEF)

    Arguments:
       IDERIV - Order of the derivative to be evaluated.  (Input)
                In particular, IDERIV = 0 returns the value of the
                polynomial.
       X      - Point at which the polynomial is to be evaluated.
                (Input)
       NINTV  - Number of polynomial pieces.  (Input)
       BREAK  - Array of length NINTV+1 containing the breakpoints for
                the piecewise cubic representation.  (Input)
                BREAK must be strictly increasing.
       CSCOEF - Matrix of size 4 by NINTV+1 containing the local
                coefficients of the cubic pieces.  (Input)
       CSDER  - Value of the IDERIV-th derivative of the polynomial at X.
                (Output)

    Keyword:    Differentiate

    GAMS:       E3; K6

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_csder(Mint *ideriv, Mfloat *x, Mint *nintv, Mfloat break_[], Mfloat *cscoef)
#else
Mfloat imsl_csder(ideriv, x, nintv, break_, cscoef)
	Mint            *ideriv;
	Mfloat          *x;
	Mint            *nintv;
	Mfloat           break_[], *cscoef;
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
	Mfloat           imsl_csder_v;


	imsl_e1psh("IMSL_CSDER ");

	imsl_csder_v = imsl_ppder(*ideriv, *x, 4, *nintv, break_, cscoef);

L_9000:
	;
	imsl_e1pop("IMSL_CSDER ");
	return (imsl_csder_v);
}				/* end of function */
