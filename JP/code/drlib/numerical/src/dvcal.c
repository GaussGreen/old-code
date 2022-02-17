#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 08/02/90 at 15:36:24
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  DVCAL (Double precision version)

    Computer:   FORC/DOUBLE

    Revised:    May 7, 1987

    Purpose:    Multiply a vector by a scalar and store the result in
                another vector, y = ax, all double precision.

    Usage:      CALL DVCAL (N, DA, DX, INCX, DY, INCY)

    Arguments:
       N      - Length of vectors X.  (Input)
       DA     - Double precision scalar.  (Input)
       DX     - Double precision vector of length MAX(N*IABS(INCX),1).
                   (Input)
                DVCAL computes DA*X(I) for I = 1,...,N. X(I) refers
                to a specific element of DX.
       INCX   - Displacement between elements of DX.  (Input)
                X(I) is defined to be DX(1+(I-1)*INCX). INCX must be
                greater than 0.
       DY     - Double precision vector of length MAX(N*IABS(INCY),1).
                   (Output)
                DVCAL sets Y(I) equal to DA*X(I) for I = 1,...,N.
                Y(I) refers to a specific element of DY.
       INCY   - Displacement between elements of DY.  (Input)
                Y(I) is defined to be DY(1+(I-1)*INCY). INCY must be
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
void imsl_dvcal(Mint n, Mdouble da, Mdouble dx[], Mint incx, Mdouble dy[],
                Mint incy)
#else
void imsl_dvcal(n, da, dx, incx, dy, incy)
	Mint            n;
	Mdouble         da, dx[];
	Mint            incx;
	Mdouble          dy[];
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
				dy[iy - 1] = da * dx[ix - 1];
				ix += incx;
				iy += incy;
			}
		} else {
			for (i = 1; i <= n; i++) {
				dy[i - 1] = da * dx[i - 1];
			}
		}
	}
	return;
}				/* end of function */
