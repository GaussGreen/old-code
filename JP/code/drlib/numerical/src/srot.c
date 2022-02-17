/*Translated by FOR_C++, v0.1, on 08/16/90 at 11:30:30 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
#include <stdio.h>
#include <math.h>
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/16/90 at 11:30:27
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SROT (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Apply a Givens plane rotation in single precision.

    Usage:      CALL SROT (N, SX, INCX, SY, INCY, C, S)

    Arguments:
       N      - Length of vectors X and Y.  (Input)
       SX     - Real vector of length MAX(N*IABS(INCX),1).
                (Input/Output)
                SROT replaces X(I) with SC*X(I) + SS*Y(I) for I=1,...,N.
                X(I) and Y(I) refer to specific elements of SX and SY.
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be
                   SX(1+(I-1)*INCX) if INCX.GE.0  or
                   SX(1+(I-N)*INCX) if INCX.LT.0.
       SY     - Real vector of length MAX(N*IABS(INCY),1).
                (Input/Output)
                SROT replaces Y(I) with -SS*X(I) + SC*Y(I) for I=1,...,N.
                X(I) and Y(I) refer to specific elements of SX and SY.
       INCY   - Displacement between elements of SY.  (Input)
                Y(I) is defined to be
                   SY(1+(I-1)*INCY) if INCY.GE.0  or
                   SY(1+(I-N)*INCY) if INCY.LT.0.
       C      - Real scalar containing elements of the rotation matrix.
                (Input)
       S      - Real scalar containing elements of the rotation matrix.
                (Input)

    Remark:
                    ( SC SS )    (X(1) ... X(N))
       SROT applies (       ) to (             )
                    (-SS SC )    (Y(1) ... Y(N))

    Keywords:   Level 1 BLAS; SROT

    GAMS:       D1a8

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_srot(Mint n, Mfloat sx[], Mint incx, Mfloat sy[], Mint incy,
               Mfloat c, Mfloat s)
#else
void imsl_srot(n, sx, incx, sy, incy, c, s)
	Mint            n;
	Mfloat          sx[];
	Mint            incx;
	Mfloat          sy[];
	Mint            incy;
	Mfloat          c, s;
#endif
{
	Mint             i, ix, iy;
	Mfloat           stemp;


	if (n > 0) {
		if (incx != 1 || incy != 1) {
			/*
			 * CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
			 * NOT EQUAL TO 1
			 */
			ix = 1;
			iy = 1;
			if (incx < 0)
				ix = (-n + 1) * incx + 1;
			if (incy < 0)
				iy = (-n + 1) * incy + 1;
			for (i = 1; i <= n; i++) {
				stemp = c * sx[ix - 1] + s * sy[iy - 1];
				sy[iy - 1] = c * sy[iy - 1] - s * sx[ix - 1];
				sx[ix - 1] = stemp;
				ix += incx;
				iy += incy;
			}
		} else {
			/* CODE FOR BOTH INCREMENTS EQUAL TO 1 */
			for (i = 1; i <= n; i++) {
				stemp = c * sx[i - 1] + s * sy[i - 1];
				sy[i - 1] = c * sy[i - 1] - s * sx[i - 1];
				sx[i - 1] = stemp;
			}
		}
	}
	return;
}				/* end of function */
