#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  PPDER/DPPDER (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 11, 1986

    Purpose:    Evaluate the derivative of a piecewise polynomial.

    Usage:      PPDER(IDERIV, X, KORDER, NINTV, BREAK, PPCOEF)

    Arguments:
       IDERIV - Order of the derivative to be evaluated.  (Input)
                In particular, IDERIV = 0 returns the value of the
                polynomial.
       X      - Point at which the polynomial is to be evaluated.
                (Input)
       KORDER - Order of the polynomial.  (Input)
       NINTV  - Number of polynomial pieces.  (Input)
       BREAK  - Array of length NINTV+1 containing the breakpoints of
                the piecewise polynomial representation.  (Input)
                BREAK must be strictly increasing.
       PPCOEF - Array of size KORDER*NINTV containing the
                local coefficients of the piecewise polynomial pieces.
                (Input)
                PPCOEF is treated internally as a matrix of size
                KORDER by NINTV.
       PPDER  - Value of the IDERIV-th derivative of the piecewise
                polynomial at X.  (Output)

    Keyword:    Differentiate

    GAMS:       E3; K6

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_ppder(Mint ideriv, Mfloat temp_x, Mint korder, Mint nintv,
                  Mfloat break_[], Mfloat *ppcoef)
#else
Mfloat imsl_ppder(ideriv, temp_x, korder, nintv, break_, ppcoef)
	Mint             ideriv;
	Mfloat           temp_x;
	Mint             korder, nintv;
	Mfloat           break_[], *ppcoef;
#endif
{
#define PPCOEF(I_,J_)	(ppcoef+(I_)*(korder)+(J_))
	Mint             j, left;
	Mfloat           fmm, h, ppder_v, value,x;


	imsl_e1psh("IMSL_PPDER");
	/* IN CASE OF ERRORS */
        x = (Mfloat)temp_x;
	value = F_ZERO;
	/* Check argument NINTV */
	if (nintv < 1) {
		imsl_e1sti(1, nintv);

/*		(5, 1, "The number of intervals must be at least 1 while NINTV = %(i1) is given. ");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NINTV_NOT_POSITIVE);
	}
	/* Check argument IDERIV */
	if (ideriv < 0) {
		imsl_e1sti(1, ideriv);

/*		(5, 2, "The order of the derivative must be positive while IDERIV = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_IDERIV_NOT_POSITIVE);
	}
	/* Check argument KORDER */
	if (korder <= 0) {
		imsl_e1sti(1, korder);

/*		(5, 3, "The order of the interpolating polynomial must be positive while KORDER = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_KORDER_NOT_POSITIVE);
	}
	/* Check for errors */
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/*
	 * Derivatives of order KORDER or higher are identically zero
	 */
	if (ideriv >= korder)
		goto L_9000;
	/*
	 * Find index I of largest breakpoint to the left of X
	 */
	imsl_p3der(korder, nintv, break_, x, &left);
	/*
	 * Evaluate jderiv-th derivative of I-th polynomial piece at X
	 */
	fmm = korder - ideriv;
	h = x - break_[left - 1];
	for (j = korder; j >= (ideriv + 1); j--) {
		value = (value / fmm) * h + *PPCOEF(left - 1, j - 1);
		fmm -= F_ONE;
	}

L_9000:
	;
	ppder_v = value;
	imsl_e1pop("IMSL_PPDER");
	return (ppder_v);
}				/* end of function */

/*  -----------------------------------------------------------------------
    IMSL Name:  P3DER/DP3DER (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 11, 1986

    Purpose:    Compute MAX (I, 1 .LE. NINTV. .AND. BREAK(I) .LE. X)

    Usage:      CALL P3DER (KORD, NINTV, BREAK, X, LEFT)

    Arguments:
       KORD   - Order of the polynomial.  (Input)
       NINTV  - Number of polynomial pieces.  (Input)
       BREAK  - Vector of length NINTV+1 containing the breakpoints
                of the piecewise polynomial representation.  (Input)
       X      - The point whose location in BREAK is to be found.
       LEFT   - Integer whose value is
                  LEFT
                    1      IF                       X .LT.  BREAK(1)
                    I      IF         BREAK(I).LE. X .LT. BREAK(I+1)
                   NINTV   IF                    BREAK(NINTV) .LE. X
                The asymmetric treatment of the interval is due to the
                decision to make all PP functions continuous from the
                right.  (Output)

    Remark:
       This routine is based in INTERV in de Boor, p92-93.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* KORD is not used here, but leave the calling sequence intact. */
#ifdef ANSI
void imsl_p3der(Mint kord, Mint nintv, Mfloat *break_, Mfloat temp_x, Mint *left)
#else
void imsl_p3der(kord, nintv, break_, temp_x, left)
	Mint             kord, nintv;
	Mfloat           break_[], temp_x;
	Mint            *left;
#endif
{
	Mint             ihi, istep, middle;
	static Mint      ilo = 1;
        Mfloat           x;

        x = (Mfloat)temp_x;
	ihi = ilo + 1;
	if (ihi >= nintv) {
		if (x >= break_[nintv - 1]) {
			*left = nintv;
			goto L_9000;
		} else if (nintv <= 1) {
			*left = 1;
			goto L_9000;
		}
		ilo = nintv - 1;
		ihi = nintv;
	}
	if (x < break_[ihi - 1]) {
		if (x >= break_[ilo - 1]) {
			*left = ilo;
			goto L_9000;
		}
		/*
		 * Now X .LT. BREAK(ILO) - decrease ILO to capture X
		 */
		istep = 1;
L_10:
		;
		ihi = ilo;
		ilo = ihi - istep;
		if (ilo > 1) {
			if (x >= break_[ilo - 1])
				goto L_30;
			istep *= 2;
			goto L_10;
		}
		ilo = 1;
		if (x < break_[0]) {
			*left = 1;
			goto L_9000;
		}
		goto L_30;
	}
	/*
	 * Now X .GE. BREAK(IHI) - increase IHI to capture X
	 */
	istep = 1;
L_20:
	;
	ilo = ihi;
	ihi = ilo + istep;
	if (ihi < nintv) {
		if (x < break_[ihi - 1])
			goto L_30;
		istep *= 2;
		goto L_20;
	}
	if (x >= break_[nintv - 1]) {
		*left = nintv;
		goto L_9000;
	}
	ihi = nintv;
	/*
	 * Now BREAK(ILO) .LE. X .LT. BREAK(IHI) - narrow the inteval
	 */
L_30:
	;
	middle = (ilo + ihi) / 2;
	if (middle == ilo) {
		*left = ilo;
		goto L_9000;
	}
	/*
	 * It is assumed that MIDDLE = ILO in case IHI = ILO+1
	 */
	if (x < break_[middle - 1]) {
		ihi = middle;
	} else {
		ilo = middle;
	}
	goto L_30;
L_9000:
	;
	return;
}				/* end of function */
