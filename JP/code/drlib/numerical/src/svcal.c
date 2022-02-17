#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 07/30/90 at 09:54:58
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SVCAL (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    May 7, 1987

    Purpose:    Multiply a vector by a scalar and store the result in
                another vector, y = ax, all single precision.

    Usage:      CALL SVCAL (N, SA, SX, INCX, SY, INCY)

    Arguments:
       N      - Length of vectors X.  (Input)
       SA     - Real scalar.  (Input)
       SX     - Real vector of length MAX(N*IABS(INCX),1).  (Input)
                SVCAL computes SA*X(I) for I = 1,...,N. X(I) refers
                to a specific element of SX.
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
                greater than 0.
       SY     - Real vector of length MAX(N*IABS(INCY),1).  (Output)
                SVCAL sets Y(I) equal to SA*X(I) for I = 1,...,N.
                Y(I) refers to a specific element of SY.
       INCY   - Displacement between elements of SY.  (Input)
                Y(I) is defined to be SY(1+(I-1)*INCY). INCY must be
                greater than 0.

    Keyword:    Level 1 BLAS

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
void imsl_svcal(Mint n, Mfloat sa, Mfloat sx[], Mint incx, Mfloat sy[], Mint incy)
#else
void imsl_svcal(n, sa, sx, incx, sy, incy)
	Mint            n;
	Mfloat          sa, sx[];
	Mint            incx;
	Mfloat           sy[];
	Mint            incy;
#endif
{
	Mint             i, ix, iy;


	if (n > 0) {
		if (incx != 1 || incy != 1) {
			/*
			 * CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
			 * NOT EQUAL TO 1
			 */
			ix = 1;
			iy = 1;
			for (i = 1; i <= n; i++) {
				sy[iy - 1] = sa * sx[ix - 1];
				ix += incx;
				iy += incy;
			}
		} else {
			for (i = 1; i <= n; i++) {
				sy[i - 1] = sa * sx[i - 1];
			}
		}
	}
	return;
}				/* end of function */
