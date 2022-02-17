#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/17/90 at 15:59:54
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CSCAL (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Multiply a vector by a scalar, y = ay, both d_complex.

    Usage:      CALL CSCAL (N, CA, CX, INCX)

    Arguments:
       N      - Length of vectors X.  (Input)
       CA     - Complex scalar.  (Input)
       CX     - Complex vector of length MAX(N*IABS(INCX),1).
                   (Input/Output)
                CSCAL replaces X(I) with CA*X(I) for I = 1,...,N.
                X(I) refers to a specific element of CX.
       INCX   - Displacement between elements of CX.  (Input)
                X(I) is defined to be CX(1+(I-1)*INCX). INCX must be
                greater than 0.

    Keywords:   Level 1 BLAS; CSCAL

    GAMS:       D1a6

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_cscal(Mint *n, Mf_complex *ca, Mf_complex *cx, Mint *incx)
#else
void imsl_cscal(n, ca, cx, incx)
	Mint            *n;
	Mf_complex      *ca, cx[];
	Mint            *incx;
#endif
{
    Mint	i, nincx;
    Mf_complex  ctemp;

    if (*n <= 0) return;
    if (*incx != 1) {
	/* CODE FOR INCREMENT NOT EQUAL TO 1 */
	nincx = *n * (*incx);
	for (i = 1;  i <= nincx;  i += (*incx)) {
	    cx[i - 1] = imsl_c_mul(*ca, cx[i - 1]);
	}
    } else {
	/* CODE FOR INCREMENT EQUAL TO 1 */
	for (i = 1; i <= *n; i++) {
#if 0
	    cx[i - 1] = imsl_c_mul(*ca, cx[i - 1]);
#endif
	    ctemp.re = ca->re*cx[i-1].re - ca->im*cx[i-1].im;
	    ctemp.im = ca->re*cx[i-1].im + ca->im*cx[i-1].re;
	    cx[i-1] = ctemp;
	}
    }
}
