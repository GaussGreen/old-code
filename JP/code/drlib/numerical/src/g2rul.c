#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  G2RUL/DG2RUL (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 1, 1985

    Purpose:    Compute a Gauss, Gauss-Radau or Gauss-Lobatto quadrature
                rule with various classical weight functions.

    Usage:      CALL G2RUL (N, IWEIGH, ALPHA, BETA, NFIX, QXFIX,
                            QX, QW, WK)

    Arguments:  (See GQRUL)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_g2rul(Mint *n, Mint *iweigh, Mfloat *alpha, Mfloat *imsl_beta, 
	        Mint *nfix, Mfloat qxfix[], Mfloat qx[], Mfloat qw[], Mfloat wk[])
#else
void imsl_g2rul(n, iweigh, alpha, imsl_beta, 
	        nfix, qxfix, qx, qw, wk)
	Mint            *n, *iweigh;
	Mfloat          *alpha, *imsl_beta;
	Mint            *nfix;
	Mfloat           qxfix[], qx[], qw[], wk[];
#endif
{

	imsl_e1psh("G2RUL ");
	/* CHECK N */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		(5, 1, "The number of quadrature points is n = %(i1).  N must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NUM_QUADRATURE_POINTS);
		goto L_10;
	}
	/* CHECK NFIX */
        /*Note that there should be no way to execute this error
message since UNKNOWN ARGUMENT would be triggered first in the supercode*/
	if (*nfix < 0 || *nfix > 2) {
		imsl_e1sti(1, *nfix);
/*	        (5, 3, "The number of fixed points NFIX = %(i1).  It must be 0, 1 or 2.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_NFIX_VALUE);
		goto L_10;
	}
	/* N GREATER THAN NFIX */
	if (*n <= *nfix) {
		imsl_e1sti(1, *n);
		imsl_e1sti(2, *nfix);
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_MORE_QUAD_POINTS);
		goto L_10;
	}
	/* CHECK NFIX AND QXFIX */
	if (*iweigh == 6) {
		if (*nfix == 2) {

/*			(5, 5, "The number of fixed quadrature points is 2.  Only 0 or 1 fixed points may be specified when using Generalized Laguerre weighting.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NFIX_CANNOT_EQUAL_TWO);
			goto L_10;
		} else if (*nfix == 1) {
			if (qxfix[0] > F_ZERO) {
				imsl_e1str(1, qxfix[0]);

/*				(5, 6, "The fixed quadrature point a = %(r1).  The fixed point must be chosen less than or equal to zero when using Generalized Laguerre weighting.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_FIXED_POINT_VALUE);
				goto L_10;
			}
		}
	} else if (*iweigh == 4 || *iweigh == 7) {
		if (*nfix > 0) {

/*			(5, 7, "Fixed quadrature points have been specified.  No fixed points may be specified for Hermite or COSH weighting.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_NO_FIXED_POINTS_ALLOWED);
			goto L_9000;
		}
	} else {
		if (*nfix >= 1) {
			if (qxfix[0] > -F_ONE && qxfix[0] < F_ONE) {
				imsl_e1str(1, qxfix[0]);

/*				(5, 8, "The fixed point, a = %(r1), must be outside the interval (-1,+1).");
*/
                                imsl_ermes(IMSL_TERMINAL,
				IMSL_A_MUST_BE_OUTSIDE_INTERVAL);
				goto L_9000;
			}
		}
		if (*nfix == 2) {
			if (qxfix[1] > -F_ONE && qxfix[1] < F_ONE) {
				imsl_e1str(1, qxfix[1]);

/*				(5, 9, "The fixed point, b = %(r1), must be outside the interval (-1,+1).");
*/
				imsl_ermes(IMSL_TERMINAL,
				IMSL_B_MUST_BE_OUTSIDE_INTERVAL);
				goto L_9000;
			}
			if ((qxfix[0] <= -F_ONE && qxfix[1] < F_ONE) || (qxfix[0] >=
						  F_ONE && qxfix[1] > -F_ONE)) {
				imsl_e1str(1, qxfix[0]);
				imsl_e1str(2, qxfix[1]);

/*				(5, 10, "The fixed points, a = %(r1) and b = %(r2), must be on opposite sides of the interval (-1, 1).");
*/
				imsl_ermes(IMSL_TERMINAL,
				IMSL_A_AND_B_ON_OPPOSITE_SIDES);
				goto L_9000;
			}
		}
	}
	/* COMPUTE RECURRENCE COEFFICIENTS */
	imsl_reccf(n, iweigh, alpha, imsl_beta, qx, qw);
	/* COMPUTE QUADRATURE RULE */
	imsl_g2rcf(n, qx, qw, nfix, qxfix, qx, qw, wk);

L_10:
	;
L_9000:
	;
	imsl_e1pop("G2RUL ");
	return;
}				/* end of function */
