#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 07/13/90 at 11:15:49
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  PPITG/DPPITG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 15, 1988

    Purpose:    Evaluate the integral of a piecewise polynomial.

    Usage:      PPITG(A, B, KORDER, NINTV, BREAK, PPCOEF)

    Arguments:
       A      - Lower limit of integration.  (Input)
       B      - Upper limit of integration.  (Input)
       KORDER - Order of the polynomial.  (Input)
       NINTV  - Number of piecewise polynomial pieces.  (Input)
       BREAK  - Array of length NINTV+1 containing the breakpoints
                for the piecewise polynomial.  (Input)
                BREAK must be strictly increasing.
       PPCOEF - Array of size KORDER*NINTV containing the
                local coefficients of the piecewise polynomial pieces.
                (Input)
                PPCOEF is treated internally as a matrix of size
                KORDER by NINTV.
       PPITG  - Value of the integral from A to B of the piecewise
                polynomial.  (Output)

    Keyword:    Quadrature

    GAMS:       E3; K6

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_ppitg(Mfloat temp_a, Mfloat temp_b, Mint korder, Mint nintv,
                  Mfloat break_[], Mfloat *ppcoef)
#else
Mfloat imsl_ppitg(temp_a, temp_b, korder, nintv, break_, ppcoef)
	Mfloat           temp_a, temp_b;
	Mint             korder, nintv;
	Mfloat           break_[], *ppcoef;
#endif
{
#define PPCOEF(I_,J_)	(ppcoef+(I_)*(korder)+(J_))
	Mint             i, ia, ib, isign = 0, j;
	Mfloat           fmm, h, ppitg_v, ra, rb, sum, value,a,b;


	imsl_e1psh("IMSL_PPITG");
	/* In case of errors */
        a = (Mfloat)temp_a;
        b = (Mfloat)temp_b;
	value = F_ZERO;
	/* Check NINTV */
	if (nintv < 1) {
		imsl_e1sti(1, nintv);

/*		(5, 1, "The number of intervals must be at least 1 while NINTV = %(i1) is given. ");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NINTV_NOT_POSITIVE);
	}
	/* Check argument KORDER */
	if (korder <= 0) {
		imsl_e1sti(1, korder);

/*	        (5, 2, "The order of the interpolating polynomial must be positive while KORDER = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_KORDER_NOT_POSITIVE);
	}
	/* Check for errors */
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/*
	 * Assign RA to MIN(A,B) Assign RB to MAX(A,B)
	 */
	if (a < b) {
		ra = a;
		rb = b;
		isign = 1;
	} else if (a > b) {
		ra = b;
		rb = a;
		isign = -1;
	} else {
		isign = 0;
		goto L_9000;
	}
	/*
	 * Find indices IA,IB of largest breakpoints to the left of A,B
	 */
	imsl_p3der(korder, nintv, break_, ra, &ia);
	imsl_p3der(korder, nintv, break_, rb, &ib);
	/* Integrate from RA to BREAK(IA) */
	fmm = korder + F_ONE;
	h = ra - break_[ia - 1];
	sum = F_ZERO;
	for (j = korder; j >= 1; j--) {
		sum = (sum / fmm) * h + *PPCOEF(ia - 1, j - 1);
		fmm -= F_ONE;
	}
	value = -sum * h;
	/*
	 * Integrate from BREAK(IA) to BREAK(IB)
	 */
	for (i = ia; i <= (ib - 1); i++) {
		h = break_[ia] - break_[ia - 1];
		sum = F_ZERO;
		fmm = korder + 1;
		for (j = korder; j >= 1; j--) {
			sum = (sum / fmm) * h + *PPCOEF(ia - 1, j - 1);
			fmm -= F_ONE;
		}
		value += sum * h;
		ia += 1;
	}
	/* Integrate from BREAK(IB) to B */
	h = rb - break_[ib - 1];
	sum = F_ZERO;
	fmm = korder + 1;
	for (j = korder; j >= 1; j--) {
		sum = (sum / fmm) * h + *PPCOEF(ib - 1, j - 1);
		fmm -= F_ONE;
	}
	value += sum * h;

L_9000:
	;
	ppitg_v = value * (Mfloat) (isign);
	imsl_e1pop("IMSL_PPITG");
	return (ppitg_v);
}				/* end of function */
