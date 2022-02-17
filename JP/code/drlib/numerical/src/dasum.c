#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 08/01/90 at 15:46:54
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  DASUM (Double precision version)

    Computer:   FORC/DOUBLE

    Revised:    August 9, 1986

    Purpose:    Sum the absolute values of the components of a
                double precision vector.

    Usage:      DASUM(N, DX, INCX)

    Arguments:
       N      - Length of vectors X.  (Input)
       DX     - Double precision vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of DX.  (Input)
                X(I) is defined to be DX(1+(I-1)*INCX).  INCX must be
                greater than 0.
       DASUM  - Double precision sum from I=1 to N of DABS(X(I)).
                (Output)
                X(I) refers to a specific element of DX.

    Keywords:   Level 1 BLAS; DASUM

    GAMS:       D1a3a

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mdouble imsl_dasum(Mint n, Mdouble *dx, Mint incx)
#else
Mdouble imsl_dasum(n, dx, incx)
	Mint            n;
	Mdouble         dx[];
	Mint            incx;
#endif
{
	Mint             _d_l, _d_m, _do0, _do1, i, nincx;
	Mdouble          dasum_v;


	dasum_v = F_ZERO;
	if (n > 0) {
		if (incx != 1) {
			/* CODE FOR INCREMENT NOT EQUAL TO 1 */
			nincx = n * incx;
			for (i = 1, _do0 = DOCNT(1, nincx, _do1 = incx); _do0 > 0; i += _do1, _do0--) {
				dasum_v += fabs(dx[i - 1]);
			}
		} else {
			for (i = 1; i <= n; i++) {
				dasum_v += fabs(dx[i - 1]);
			}
		}
	}
	return (dasum_v);
}				/* end of function */
