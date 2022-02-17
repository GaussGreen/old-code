#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/24/90 at 17:35:17
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CSET (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Set the components of a vector to a scalar, all d_complex.

    Usage:      CALL CSET (N, CA, CX, INCX)

    Arguments:
       N      - Length of vector X.  (Input)
       CA     - Complex scalar.  (Input)
       CX     - Complex vector of length N*INCX.  (Input/Output)
                CSET replaces X(I) with CA for I=1,...,N. X(I) refers to
                a specific element of CX. See INCX argument description.
       INCX   - Displacement between elements of CX.  (Input)
                X(I) is defined to be CX(1+(I-1)*INCX). INCX must be
                greater than zero.

    Keyword:    Level 1 BLAS

    GAMS:       D1a1

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_cset(Mint *n, Mf_complex *ca, Mf_complex *cx, Mint *incx)
#else
void imsl_cset(n, ca, cx, incx)
	Mint            *n;
	Mf_complex      *ca, cx[];
	Mint            *incx;
#endif
{
	Mint             _d_l, _d_m, _do0, _do1, i, nincx;


	if (*n > 0) {
		if (*incx != 1) {
			/* CODE FOR INCREMENT NOT EQUAL TO 1 */
			nincx = *n ** incx;
			for (i = 1, _do0 = DOCNT(1, nincx, _do1 = *incx); _do0 > 0; i += _do1, _do0--) {
				cx[i - 1] = *ca;
			}
		} else {
			/*
			 * CODE FOR INCREMENT EQUAL TO 1 CLEAN-UP LOOP
			 */
			for (i = 1; i <= *n; i++) {
				cx[i - 1] = *ca;
			}
		}
	}
	return;
}				/* end of function */
