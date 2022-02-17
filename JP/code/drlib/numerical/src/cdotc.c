#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/24/90 at 17:40:40
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CDOTC (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Compute the d_complex conjugate dot product, imsl_c_conjg(x)*y.

    Usage:      CDOTC(N, CX, INCX, CY, INCY)

    Arguments:
       N      - Length of vectors X and Y.  (Input)
       CX     - Complex vector of length MAX(N*IABS(INCX),1).  (Input)
       INCX   - Displacement between elements of CX.  (Input)
                X(I) is defined to be
                   CX(1+(I-1)*INCX) if INCX.GE.0  or
                   CX(1+(I-N)*INCX) if INCX.LT.0.
       CY     - Complex vector of length MAX(N*IABS(INCY),1).  (Input)
       INCY   - Displacement between elements of CY.  (Input)
                Y(I) is defined to be
                   CY(1+(I-1)*INCY) if INCY.GE.0  or
                   CY(1+(I-N)*INCY) if INCY.LT.0.
       CDOTC  - Complex sum from I=1 to N of CONJ(X(I))*Y(I).  (Output)
                X(I) and Y(I) refer to specific elements of CX and CY
                respectively.

    Keywords:   Level 1 BLAS; CDOTC; Inner product; Scalar product

    GAMS:       D1a4

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mf_complex imsl_cdotc(Mint *n, Mf_complex *cx, Mint *incx,
		Mf_complex *cy, Mint *incy)
#else
Mf_complex imsl_cdotc(n, cx, incx, cy, incy)
	Mint            *n;
	Mf_complex       cx[];
	Mint            *incx;
	Mf_complex       cy[];
	Mint            *incy;
#endif
{
	Mint             i, ix, iy;
	Mf_complex       cdotc_v;
	Mf_complex	ctemp;


	cdotc_v = imsl_cf_convert(F_ZERO, F_ZERO);
	if (*n > 0) {
		if (*incx != 1 || *incy != 1) {
			/*
			 * CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
			 * NOT EQUAL TO 1
			 */
			ix = 1;
			iy = 1;
			if (*incx < 0)
				ix = (-*n + 1) ** incx + 1;
			if (*incy < 0)
				iy = (-*n + 1) ** incy + 1;
			for (i = 1; i <= *n; i++) {
				cdotc_v = imsl_c_add(cdotc_v, imsl_c_mul(imsl_c_conjg(cx[ix - 1]),
							       cy[iy - 1]));
				ix += *incx;
				iy += *incy;
			}
		} else {
			/* CODE FOR BOTH INCREMENTS EQUAL TO 1 */
			for (i = 1; i <= *n; i++) {
#if 0
				cdotc_v = imsl_c_add(cdotc_v, imsl_c_mul(imsl_c_conjg(cx[i - 1]), cy[i - 1]));
#endif
				ctemp.re = cx[i-1].re*cy[i-1].re + cx[i-1].im*cy[i-1].im;
				ctemp.im = cx[i-1].re*cy[i-1].im - cx[i-1].im*cy[i-1].re;
				cdotc_v.re += ctemp.re;
				cdotc_v.im += ctemp.im;
			}
		}
	}
	return (cdotc_v);
}				/* end of function */
