#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/24/90 at 17:35:41
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  ICAMAX (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 8, 1989

    Purpose:    Find the smallest index of the component of a d_complex
                vector having maximum magnitude.

    Usage:      ICAMAX(N, CX, INCX)

    Arguments:
       N      - Length of vector X.  (Input)
       CX     - Complex vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of CX.  (Input)
                X(I) is defined to be CX(1+(I-1)*INCX). INCX must be
                greater than zero.
       ICAMAX - The smallest index I such that CABS(X(I)) is the maximum
                of CABS(X(J)) for J=1 to N.  (Output)
                X(I) refers to a specific element of CX.

    Keywords:   Level 1 BLAS; ICAMAX

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
Mint imsl_icamax(Mint *n, Mf_complex *cx, Mint *incx)
#else
Mint imsl_icamax(n, cx, incx)
	Mint            *n;
	Mf_complex       cx[];
	Mint            *incx;
#endif
{
	Mint             i, icamax_v, ix;
	Mfloat           smax;
	Mf_complex       cdum;


#define CABS1(cdum)	(Mfloat)(fabs( imsl_fc_convert( (cdum) ) ) + fabs( imsl_c_aimag( (cdum) ) ))

	icamax_v = 0;
	if (*n >= 1) {
		icamax_v = 1;
		if (*n != 1) {
			if (*incx != 1) {
				/* CODE FOR INCREMENT NOT EQUAL TO 1 */
				ix = 1;
				smax = CABS1(cx[ix - 1]);
				ix += *incx;
				for (i = 2; i <= *n; i++) {
					if (CABS1(cx[ix - 1]) > smax) {
						icamax_v = i;
						smax = CABS1(cx[ix - 1]);
					}
					ix += *incx;
				}
			} else {
				/* CODE FOR INCREMENT EQUAL TO 1 */
				smax = CABS1(cx[0]);
				for (i = 2; i <= *n; i++) {
					if (CABS1(cx[i - 1]) > smax) {
						icamax_v = i;
						smax = CABS1(cx[i - 1]);
					}
				}
			}
		}
	}
	return (icamax_v);
#undef	CABS1
}				/* end of function */
