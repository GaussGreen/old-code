#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 01/08/90 at 17:17:59
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  ISAMAX (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Find the smallest index of the component of a
                single-precision vector having maximum absolute value.

    Usage:      ISAMAX(N, SX, INCX)

    Arguments:
       N      - Length of vector X.  (Input)
       SX     - Real vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
                greater than zero.
       ISAMAX - The smallest index I such that ABS(X(I)) is the maximum
                of ABS(X(J)) for J=1 to N.  (Output)
                X(I) refers to a specific element of SX. see INCX
                argument description.

    Keywords:   Level 1 BLAS; ISAMAX

    GAMS:       D1a2

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
*/
#ifdef ANSI
Mint imsl_isamax(Mint n, Mfloat *sx, Mint incx)
#else
Mint imsl_isamax(n, sx, incx)
    Mint             n;
    Mfloat           sx[];
    Mint             incx;
#endif
{
    Mint             i, ii, isamax_v, ns;
    Mfloat           smax, xmag;


	isamax_v = 0;
	if (n >= 1) {
		isamax_v = 1;
		if (n > 1) {
			if (incx != 1) {
				/* CODE  FOR INCREMENTS NOT EQUAL TO 1. */
				smax = fabs(sx[0]);
				ns = n * incx;
				ii = 1;
/*				for (i = 1, _do0 = DOCNT(1, ns, _do1 = incx); _do0 > 0; i += _do1, _do0--) { */
				for (i = 1;  i <= ns;  i += incx) {
					xmag = fabs(sx[i - 1]);
					if (xmag > smax) {
						isamax_v = ii;
						smax = xmag;
					}
					ii += 1;
				}
			} else {
				/* CODE FOR INCREMENTS EQUAL TO 1. */
				smax = fabs(sx[0]);
				for (i = 2; i <= n; i++) {
					xmag = fabs(sx[i - 1]);
					if (xmag > smax) {
						isamax_v = i;
						smax = xmag;
					}
				}
			}
		}
	}
	return (isamax_v);
}				/* end of function */
