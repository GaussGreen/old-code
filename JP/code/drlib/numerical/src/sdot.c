#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 01/08/90 at 17:18:57
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SDOT (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Compute the single-precision dot product x*y.

    Usage:      SDOT(N, SX, INCX, SY, INCY)

    Arguments:
       N      - Length of vectors X and Y.  (Input)
       SX     - Real vector of length MAX(N*IABS(INCX),1).  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be.. SX(1+(I-1)*INCX) if INCX .GE. 0
                or SX(1+(I-N)*INCX) if INCX .LT. 0.
       SY     - Real vector of length MAX(N*IABS(INCY),1).  (Input)
       INCY   - Displacement between elements of SY.  (Input)
                Y(I) is defined to be.. SY(1+(I-1)*INCY) if INCY .GE. 0
                or SY(1+(I-N)*INCY) if INCY .LT. 0.
       SDOT   - Sum from I=1 to N of X(I)*Y(I).  (Output)
                X(I) and Y(I) refer to specific elements of SX and SY,
                respectively.  See INCX and INCY argument descriptions.

    Keywords:   Level 1 BLAS; SDOT; Inner product; Scalar product

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
Mfloat imsl_sdot(Mint n, Mfloat *sx, Mint incx, Mfloat *sy, Mint incy)
#else
Mfloat imsl_sdot(n, sx, incx, sy, incy)
    Mint             n;
    Mfloat           sx[];
    Mint             incx;
    Mfloat           sy[];
    Mint             incy;
#endif
{
	Mint             i, ix, iy;
	Mfloat           sdot_v;


	sdot_v = F_ZERO;
	if (n > 0) {
		if (incx != 1 || incy != 1) {
			/* CODE FOR UNEQUAL INCREMENTS */
			ix = 1;
			iy = 1;
			if (incx < 0)
				ix = (-n + 1) * incx + 1;
			if (incy < 0)
				iy = (-n + 1) * incy + 1;
			for (i = 1; i <= n; i++) {
				sdot_v += sx[ix - 1] * sy[iy - 1];
				ix += incx;
				iy += incy;
			}
		} else {
			for (i = 1; i <= n; i++) {
				sdot_v += sx[i - 1] * sy[i - 1];
			}
		}
	}
	return (sdot_v);
}				/* end of function */
