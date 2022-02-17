#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/25/90 at 15:25:43
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CADD (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Add a scalar to each component of a vector, x = x + a,
                all d_complex.

    Usage:      CALL CADD (N, CA, CX, INCX)

    Arguments:
       N      - Length of vector X.  (Input)
       CA     - Complex scalar added to each element of X.  (Input)
       CX     - Complex vector of length MAX(N*IABS(INCX),1).
                (Input/Output)
                CADD replaces X(I) with X(I) + CA for I = 1,...N.
                X(I) refers to a specific element of CX.
       INCX   - Displacement between elements of CX.  (Input)
                X(I) is defined to be
                   CX(1+(I-1)*INCX) if INCX.GE.0  or
                   CX(1+(I-N)*INCX) if INCX.LT.0.

    Keyword:    Level 1 BLAS

    GAMS:       D1a

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_cadd(Mint *n, Mf_complex *ca, Mf_complex *cx, Mint *incx)
#else
void imsl_cadd(n, ca, cx, incx)
	Mint            *n;
	Mf_complex      *ca, cx[];
	Mint            *incx;
#endif
{
	Mint             i, ix;
	Mfloat           ccabs;


	ccabs = fabs(imsl_fc_convert(*ca)) + fabs(imsl_c_aimag(*ca));
	if (ccabs == F_ZERO)
		return;
	if (*n > 0) {
		if (*incx != 1) {
			/* CODE FOR INCREMENT NOT EQUAL TO 1 */
			ix = 1;
			if (*incx < 0)
				ix = (-*n + 1) ** incx + 1;
			for (i = 1; i <= *n; i++) {
				cx[ix - 1] = imsl_c_add(*ca, cx[ix - 1]);
				ix += *incx;
			}
		} else {
			/* CODE FOR INCREMENT EQUAL TO 1 */
			for (i = 1; i <= *n; i++) {
#if 0
				cx[i - 1] = imsl_c_add(*ca, cx[i - 1]);
#endif
				cx[i-1].re += ca->re;
				cx[i-1].im += ca->im;
			}
		}
	}
	return;
}				/* end of function */
