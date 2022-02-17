#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 07/26/90 at 14:39:22
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SASUM (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Sum the absolute values of the components of a
                single precision vector.

    Usage:      SASUM(N, SX, INCX)

    Arguments:
       N      - Length of vectors X.  (Input)
       SX     - Real vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
                greater than 0.
       SASUM  - Single precision sum from I=1 to N of ABS(X(I)).
                (Output)
                X(I) refers to a specific element of SX.

    Keywords:   Level 1 BLAS; SASUM

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
Mfloat imsl_sasum(Mint n, Mfloat *sx, Mint incx)
#else
Mfloat imsl_sasum(n, sx, incx)
	Mint            n;
	Mfloat          sx[];
	Mint            incx;
#endif
{
	Mint             _d_l, _d_m, _do0, _do1, i, nincx;
	Mfloat           sasum_v;


	sasum_v = F_ZERO;
	if (n > 0) {
		if (incx != 1) {
			/* CODE FOR INCREMENT NOT EQUAL TO 1 */
			nincx = n * incx;
			for (i = 1, _do0 = DOCNT(1, nincx, _do1 = incx); _do0 > 0; i += _do1, _do0--) {
				sasum_v += fabs(sx[i - 1]);
			}
		} else {
			for (i = 1; i <= n; i++) {
				sasum_v += fabs(sx[i - 1]);
			}
		}
	}
	return (sasum_v);
}				/* end of function */
