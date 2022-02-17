#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/24/90 at 17:40:12
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CAXPY (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Compute the scalar times a vector plus a vector,
                y = ax + y, all d_complex.

    Usage:      CALL CAXPY (N, CA, CX, INCX, CY, INCY)

    Arguments:
       N      - Length of vectors X and Y.  (Input)
       CA     - Complex scalar.  (Input)
       CX     - Complex vector of length MAX(N*IABS(INCX),1).  (Input)
       INCX   - Displacement between elements of CX.  (Input)
                X(I) is defined to be
                   CX(1+(I-1)*INCX) if INCX.GE.0  or
                   CX(1+(I-N)*INCX) if INCX.LT.0.
       CY     - Complex vector of length MAX(N*IABS(INCY),1).
                   (Input/Output)
                CAXPY replaces Y(I) with CA*X(I) + Y(I) for I = 1,...N.
                X(I) and Y(I) refer to specific elements of CX and CY.
       INCY   - Displacement between elements of CY.  (Input)
                Y(I) is defined to be
                   CY(1+(I-1)*INCY) if INCY.GE.0  or
                   CY(1+(I-N)*INCY) if INCY.LT.0.

    Keywords:   Level 1 BLAS; CAXPY

    GAMS:       D1a7

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_caxpy(Mint *n, Mf_complex *ca, Mf_complex *cx, Mint *incx, 
		Mf_complex *cy, Mint *incy)
#else
void /* FUNCTION */ 
imsl_caxpy(n, ca, cx, incx, cy, incy)
	Mint            *n;
	Mf_complex      *ca, cx[];
	Mint            *incx;
	Mf_complex       cy[];
	Mint            *incy;
#endif
{
	Mint             i, ix, iy;
	Mf_complex	ctemp;


	if (*n > 0) {
		if (fabs(imsl_fc_convert(*ca)) + fabs(imsl_c_aimag(*ca)) != F_ZERO) {
			if (*incx != 1 || *incy != 1) {
				/*
				 * CODE FOR UNEQUAL INCREMENTS OR EQUAL
				 * INCREMENTS NOT EQUAL TO 1
				 */
				ix = 1;
				iy = 1;
				if (*incx < 0)
					ix = (-*n + 1) ** incx + 1;
				if (*incy < 0)
					iy = (-*n + 1) ** incy + 1;
				for (i = 1; i <= *n; i++) {
					cy[iy - 1] = imsl_c_add(cy[iy - 1], imsl_c_mul(*ca, cx[ix - 1]));
					ix += *incx;
					iy += *incy;
				}
				/* CODE FOR BOTH INCREMENTS EQUAL TO 1 */
			} else {
				for (i = 1; i <= *n; i++) {
#if 0
					cy[i - 1] = imsl_c_add(cy[i - 1], imsl_c_mul(*ca, cx[i - 1]));
#endif
					ctemp.re = ca->re*cx[i-1].re - ca->im*cx[i-1].im;
					ctemp.im = ca->re*cx[i-1].im + ca->im*cx[i-1].re;
					cy[i-1].re += ctemp.re;
					cy[i-1].im += ctemp.im;
				}
			}
		}
	}
	return;
}				/* end of function */
