#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/18/90 at 12:46:26
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CSWAP (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Interchange vectors X and Y, both d_complex.

    Usage:      CALL CSWAP (N, CX, INCX, CY, INCY)

    Arguments:
       N      - Length of vectors X and Y.  (Input)
       CX     - Complex vector of length MAX(N*IABS(INCX),1).
                (Input/Output)
       INCX   - Displacement between elements of CX.  (Input)
                X(I) is defined to be
                   CX(1+(I-1)*INCX) if INCX.GE.0  or
                   CX(1+(I-N)*INCX) if INCX.LT.0.
       CY     - Complex vector of length MAX(N*IABS(INCY),1).
                (Input/Output)
       INCY   - Displacement between elements of CY.  (Input)
                Y(I) is defined to be
                   CY(1+(I-1)*INCY) if INCY.GE.0  or
                   CY(1+(I-N)*INCY) if INCY.LT.0.

    Keywords:   Level 1 BLAS; CSWAP; Swap; Exchange

    GAMS:       D1a5

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_cswap(Mint *n, Mf_complex *cx, Mint *incx, Mf_complex *cy,
			Mint *incy)
#else
void imsl_cswap(n, cx, incx, cy, incy)
	Mint            *n;
	Mf_complex       cx[];
	Mint            *incx;
	Mf_complex       cy[];
	Mint            *incy;
#endif
{
	Mint             i, ix, iy;
	Mf_complex       ctemp;


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
				ctemp = cx[ix - 1];
				cx[ix - 1] = cy[iy - 1];
				cy[iy - 1] = ctemp;
				ix += *incx;
				iy += *incy;
			}
		} else {
			/* CODE FOR BOTH INCREMENTS EQUAL TO 1 */
			for (i = 1; i <= *n; i++) {
				ctemp = cx[i - 1];
				cx[i - 1] = cy[i - 1];
				cy[i - 1] = ctemp;
			}
		}
	}
	return;
}				/* end of function */
