#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  CSITG/DCSITG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 15, 1988

    Purpose:    Evaluate the integral of a cubic spline.

    Usage:      CSITG(A, B, NINTV, BREAK, CSCOEF)

    Arguments:
       A      - Lower limit of integration.  (Input)
       B      - Upper limit of integration.  (Input)
       NINTV  - Number of polynomial pieces.  (Input)
       BREAK  - Array of length NINTV+1 containing the breakpoints
                for the piecewise cubic representation.  (Input)
                BREAK must be strictly increasing.
       CSCOEF - Matrix of size 4 by NINTV+1 containing the local
                coefficients of the cubic pieces.  (Input)
       CSITG  - Value of the integral of the spline from A to B.
                (Output)

    Keyword:    Quadrature

    GAMS:       E3; K6

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_csitg(Mfloat *a, Mfloat *b, Mint *nintv, Mfloat break_[], Mfloat *cscoef)
#else
Mfloat imsl_csitg(a, b, nintv, break_, cscoef)
	Mfloat          *a, *b;
	Mint            *nintv;
	Mfloat           break_[], *cscoef;
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
	Mfloat           imsl_csitg_v;


	imsl_e1psh("IMSL_CSITG ");

	imsl_csitg_v = imsl_ppitg(*a, *b, 4, *nintv, break_, cscoef);

L_9000:
	;
	imsl_e1pop("IMSL_CSITG ");
	return (imsl_csitg_v);
}				/* end of function */
