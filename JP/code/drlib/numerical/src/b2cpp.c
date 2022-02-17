#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  B2CPP/DB2CPP (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    August 7, 1986

    Purpose:    Convert a spline in B-spline representation to a
                piecewise polynomial representation.

    Usage:      CALL B2CPP (KORDER, XKNOT, NCOEF, BSCOEF, NPPCF, BREAK,
                            PPCOEF, WK)

    Arguments:
       KORDER - Order of the spline.  (Input)
       XKNOT  - Array of length KORDER+NCOEF containing the knot
                sequence.  (Input)
                XKNOT must be nondecreasing.
       NCOEF  - Number of B-spline coefficients.  (Input)
       BSCOEF - Array of length NCOEF containing the B-spline
                coefficients.  (Input)
       NPPCF  - Number of piecewise polynomial pieces.  (Output)
                NPPCF is always less than or equal to NCOEF-KORDER+1.
       BREAK  - Array of length (NPPCF+1) containing the breakpoints
                of the piecewise polynomial representation.  (Output)
                BREAK must be dimensioned at least NCOEF-KORDER+2.
       PPCOEF - Array of length KORDER*NPPCF containing the local
                coefficients of the polynomial pieces.  (Output)
                PPCOEF is treated internally as a matrix of size
                KORDER by NPPCF.
       WK     - Work array of length (KORDER+3)*KORDER.

    Remark:
       Informational errors
       Type Code
         4   4  Multiplicity of the knots cannot exceed the order
                of the spline.
         4   5  The knots must be nondecreasing.

    GAMS:       E3; K6

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b2cpp(Mint *korder, Mfloat xknot[], Mint *ncoef, Mfloat bscoef[], Mint *nppcf, Mfloat break_[],
	   Mfloat ppcoef[], Mfloat wk[])
#else
void imsl_b2cpp(korder, xknot, ncoef, bscoef, nppcf, break_,
	   ppcoef, wk)
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           bscoef[];
	Mint            *nppcf;
	Mfloat           break_[], ppcoef[], wk[];
#endif
{

	imsl_e1psh("IMSL_B2CPP");
	/* CHECK KORDER */
	if (*korder < 1) {
		imsl_e1sti(1, *korder);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);
		goto L_9000;
	}
	/* CHECK NCOEF */
	if (*ncoef < *korder) {
		imsl_e1sti(1, *ncoef);
		imsl_e1sti(2, *korder);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_COEFFS);
		goto L_9000;
	}
	/* CHECK ARGUMENTS FOR ERRORS */
	imsl_b3int(korder, xknot, ncoef);
	if (imsl_n1rty(0) != 0)
		goto L_9000;

	imsl_b3cpp(korder, xknot, ncoef, bscoef, nppcf, break_, ppcoef, &wk[0],
		   &wk[*korder], &wk[*korder * 2], &wk[*korder * 3]);

L_9000:
	;
	imsl_e1pop("IMSL_B2CPP");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B3CPP/DB3CPP (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    August 7, 1986

    Purpose:    Convert a spline in B-spline representation to
                piecewise polynomial representation.

    Usage:      CALL B3CPP (KORDER, XKNOT, NCOEF, BSCOEF, NPPCF, BREAK,
                            PPCOEF, BIATX, DELTAL, DELTAR, SCRTCH)

    Arguments:
       KORDER - Order of the spline.  (Input)
       XKNOT  - Array of length KORDER+NCOEF containing the knot
                sequence.  (Input)
                XKNOT must be nondecreasing.
       NCOEF  - Length of BSCOEF.  (Input)
       BSCOEF - Array of length NCOEF containing the B-spline
                coefficients.  (Input)
       NPPCF  - Number of piecewise polynomial pieces.  (Output)
                It satisfies NPPCF .LE. NCOEF-KORDER+1.
       BREAK  - Array of length (NPPCF+1) containing the breakpoints
                of the piecewise polynomial representation.  (Output)
       PPCOEF - KORDER by NPPCF array containing the piecewise
                polynomial coefficients.  (Output)
       BIATX  - Work array of length KORDER for the nonzero splines at
                a point.
       DELTAL - Work array of length KORDER used by B3INT.
       DELTAR - Work array of length KORDER used by B3INT.
       SCRTCH - Work array of length KORDER**2.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b3cpp(Mint *korder, Mfloat xknot[], Mint *ncoef, Mfloat bscoef[], Mint *nppcf, Mfloat break_[],
	   Mfloat *ppcoef, Mfloat biatx[], Mfloat deltal[], Mfloat deltar[], Mfloat *scrtch)
#else
void imsl_b3cpp(korder, xknot, ncoef, bscoef, nppcf, break_,
	   ppcoef, biatx, deltal, deltar, scrtch)
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           bscoef[];
	Mint            *nppcf;
	Mfloat           break_[], *ppcoef, biatx[], deltal[], deltar[],
	               *scrtch;
#endif
{
#define PPCOEF(I_,J_)	(ppcoef+(I_)*(*korder)+(J_))
#define SCRTCH(I_,J_)	(scrtch+(I_)*(*korder)+(J_))
	Mint             i, j, left, lsofar;
	Mfloat           imsl_diff, saved, term, xnext;


	imsl_e1psh("IMSL_B3CPP");

	lsofar = 0;
	break_[0] = xknot[*korder - 1];
	for (left = *korder; left <= *ncoef; left++) {
		/*
		 * Find the next nontrivial knot interval.
		 */
		if (xknot[left] == xknot[left - 1])
			goto L_50;
		lsofar += 1;
		xnext = xknot[left];
		break_[lsofar] = xnext;
		if (*korder <= 1) {
			*PPCOEF(lsofar - 1, 0) = bscoef[left - 1];
			goto L_50;
		}
		/*
		 * Store the KORDER B-spline coeff.S relevant to current knot
		 * interval in SCRTCH(.,1) .
		 */
		scopy(*korder, &bscoef[left - *korder], 1, scrtch, 1);

		/*
		 * FOR J=1,...,KORDER-1, compute the KORDER-J B-spline
		 * coeff.s relevant to current knot interval for the J-th
		 * derivative by differencing those for the (J-1)st
		 * derivative, and store in SCRTCH(.,J+1) .
		 */
		for (j = 1; j <= (*korder - 1); j++) {
			for (i = 1; i <= (*korder - j); i++) {
				imsl_diff = xknot[left + i - 1] - xknot[left + i - *korder + j - 1];
				if (imsl_diff > F_ZERO) {
					*SCRTCH(j, i - 1) = ((*SCRTCH(j - 1, i) - *SCRTCH(j - 1, i - 1)) /
					 imsl_diff) * (Mfloat) (*korder - j);
				}
			}
		}
		/*
		 * FOR J = 0, ..., KORDER-1, find the values at XKNOT(LEFT)
		 * of the J+1 B-splines of order J+1 whose support contains
		 * the current knot interval from those of order J (in BIATX
		 * ), then combine with the B-spline coeff.s (in
		 * SCRTCH(.,KORDER-J) ) found earlier to compute the
		 * (KORDER-J-1)st derivative at XKNOT(LEFT) of the given
		 * spline.
		 */
		biatx[0] = F_ONE;
		*PPCOEF(lsofar - 1, *korder - 1) = *SCRTCH(*korder - 1, 0);
		for (j = 1; j <= (*korder - 1); j++) {
			deltar[j - 1] = xknot[left + j - 1] - xknot[left - 1];
			deltal[j - 1] = xknot[left - 1] - xknot[left - j];
			saved = F_ZERO;
			for (i = 1; i <= j; i++) {
				term = biatx[i - 1] / (deltar[i - 1] + deltal[j - i]);
				biatx[i - 1] = saved + deltar[i - 1] * term;
				saved = deltal[j - i] * term;
			}
			biatx[j] = saved;
			*PPCOEF(lsofar - 1, *korder - j - 1) = imsl_sdot(j + 1, biatx,
					  1, SCRTCH(*korder - j - 1, 0), 1);
		}
L_50:
		;
	}
	*nppcf = lsofar;

	imsl_e1pop("IMSL_B3CPP");
	return;
}				/* end of function */
