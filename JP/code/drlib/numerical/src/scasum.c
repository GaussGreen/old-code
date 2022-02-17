#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/24/90 at 17:36:07
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SCASUM (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Sum the absolute values of the real part together with
                the absolute values of the imaginary part of the
                components of a d_complex vector.

    Usage:      SCASUM(N, CX, INCX)

    Arguments:
       N      - Length of vectors X.  (Input)
       CX     - Complex vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of CX.  (Input)
                X(I) is defined to be CX(1+(I-1)*INCX).  INCX must be
                greater than 0.
       SCASUM - Sum from I=1 to N of ABS(REAL(X(I)))+ABS(AIMAG(X(I)))).
                (Output)
                X(I) refers to a specific element of CX.

    Keywords:   Level 1 BLAS; SCASUM

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
Mfloat imsl_scasum(Mint *n, Mf_complex *cx, Mint *incx)
#else
Mfloat imsl_scasum(n, cx, incx)
	Mint            *n;
	Mf_complex       cx[];
	Mint            *incx;
#endif
{
	Mint             _d_l, _d_m, _do0, _do1, i, nincx;
	Mfloat           scasum_v;


	scasum_v = F_ZERO;
	if (*n > 0) {
		if (*incx != 1) {
			/* CODE FOR INCREMENT NOT EQUAL TO 1 */
			nincx = *n ** incx;
			for (i = 1, _do0 = DOCNT(1, nincx, _do1 = *incx); _do0 > 0; i += _do1, _do0--) {
				scasum_v += fabs(imsl_fc_convert(cx[i - 1])) + fabs(imsl_c_aimag(cx[i - 1]));
			}
		} else {
			/* CODE FOR INCREMENT EQUAL TO 1 */
			for (i = 1; i <= *n; i++) {
#if 0
				scasum_v += fabs(imsl_fc_convert(cx[i - 1])) + fabs(imsl_c_aimag(cx[i - 1]));
#endif
				scasum_v += fabs(cx[i - 1].re) + fabs(cx[i - 1].im);
			}
		}
	}
	return (scasum_v);
}				/* end of function */
