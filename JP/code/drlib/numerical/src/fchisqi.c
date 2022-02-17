#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  CHIIN/DCHIIN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 25, 1985

    Purpose:    Evaluate the inverse of the imsl_chi-squared distribution
                function.

    Usage:      CHIIN(P, DF)

    Arguments:
       P      - Probability for which the inverse of the imsl_chi-squared
                distribution function is to be evaluated.  (Input)
                P must be in the open interval (0.0, 1.0).
       DF     - Number of degrees of freedom of the imsl_chi-squared
                distribution.  (Input)
                DF must be greater than or equal to 0.5.
       CHIIN  - Function value.  (Output)
                The probability that a imsl_chi-squared random variable
                takes a value less than or equal to CHIIN is P.

    Remark:
       Informational errors
       Type Code
         3   3  The bounds which enclose P could not be found. An
                approximation for CHIIN is returned.
         3   4  The value of the inverse imsl_chi-squared could not be found
                within a specified number of iterations. An approximation
                for CHIIN is returned.

    Keywords:   Percentage point; Probability distribution; Continuous
                random variables; Percentile; Fractile

    GAMS:       L5a2c; C7e

    Chapters:   STAT/LIBRARY Probability Distribution Functions and
                             Inverses
                SFUN/LIBRARY Probability Distribution Functions and
                             Inverses

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_chi_squared_inverse_cdf(Mfloat p, Mfloat df)
#else
Mfloat imsl_f_chi_squared_inverse_cdf(p, df)
	Mfloat          p, df;
#endif
{
	short int       ibisec;
	int             iter;
	Mfloat          alam=0.0, chiin_v, delta, eps, f1, f2, f3, fd, slope,
	                x0, x1, x2, x3, xd, xint, xm;


	E1PSH("imsl_f_chi_squared_inverse_cdf","imsl_d_chi_squared_inverse_cdf");
	chiin_v = imsl_amach(6);
	/* Check P */
	if (p <= F_ZERO || p >= F_ONE) {
		imsl_e1str(1, p);

/*		imsl_ermes(5, 1, "The input probability, P = %(r1), must be in the open interval (0.0, 1.0).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_P_OUTSIDE_OPEN_INTERVAL);
	}
	/* Check DF */
	if (df < F_HALF) {
		imsl_e1str(1, df);

/*		imsl_ermes(5, 2, "The degrees of freedom, DF = %(r1), must be greater than or equal to 0.5.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_DF_MUST_BE_GE_POINT_5);
	}
	if (imsl_n1rty(0) == 5)
		goto L_9000;

	eps = F_TEN * imsl_amach(4);
	/* Use adjusted normal as starting value */
	xint = imsl_f_normal_inverse_cdf(p);
	x0 = F_TWO / (F_NINE * df);
	x1 = df * pow(F_ONE - x0 + xint * sqrt(x0), F_THREE);

	if (x1 < F_ZERO)
		x1 = F_ZERO;
	f1 = imsl_f_chi_squared_cdf(x1, df) - p;
	if (f1 == F_ZERO) {
		x2 = x1;
		goto L_30;
	}
	if (fabs(xint) >= F_ONE) {
		x2 = xint * 1.05;
		xd = x2 - x1;
	} else {
		x2 = xint + 0.05;
		xd = 0.05;
	}
	if (x2 < F_ZERO)
		x2 = F_ZERO;
	f2 = imsl_f_chi_squared_cdf(x2, df) - p;
	/* Bracket estimate */
	slope = imsl_f_max(.01, (f2 - f1) / xd);
	delta = -f1 / slope;
	iter = 0;
L_10:
	delta *= F_TWO;
	iter += 1;
	if (iter > 100) {
		imsl_e1str(1, p);
		imsl_e1str(2, df);
		imsl_e1str(3, alam);

/*		imsl_ermes(5, 6, "Unable to bracket the value of the inverse noncentral imsl_chi-squared at P = %(r1), with parameters DF = %(r2) and ALAM = %(r3).");
*/
                imsl_ermes(IMSL_WARNING, IMSL_UNABLE_TO_BRACKET_VALUE);
		goto L_9000;
	}
	x2 = x1 + delta;
	if (x2 < F_ZERO)
		x2 = F_ZERO;
	f2 = imsl_f_chi_squared_cdf(x2, df) - p;
	if (f1 * f2 >= F_ZERO) {
		x1 = x2;
		goto L_10;
	}
	/* Use regula falsi to get the estimate */
	ibisec = 0;
	for (iter = 1; iter <= 100; iter++) {
		xm = (x1 + x2) / F_TWO;
		fd = f2 - f1;
		xd = x2 - x1;
		if (fabs(xm) != F_ZERO) {
			if (fabs(xd / xm) < eps)
				goto L_30;
		} else {
			if (fabs(xd) < eps)
				goto L_30;
		}
		/*
		 * X3 is the intersection of the secant line through(X1,F1),
		 * (X2,F2) and the X axis
		 */
		if (ibisec) {
			x3 = xm;
		} else {
			x3 = x2 - f2 * xd / fd;
		}
		ibisec = 0;
		if (x3 < F_ZERO)
			x3 = F_ZERO;
		f3 = imsl_f_chi_squared_cdf(x3, df) - p;
		if (f3 * f2 <= F_ZERO) {
			/* Root was trapped, so use regula falsi */
			x1 = x2;
			f1 = f2;
			x2 = x3;
			f2 = f3;
		} else {
			/*
			 * Root was not trapped use Illinois modification
			 */
			x2 = x3;
			f2 = f3;
			f1 /= F_TWO;
			if (fabs(f2) > fabs(f1)) {
				/* Use bisection */
				f1 *= F_TWO;
				ibisec = 1;
			}
		}
	}

/*	imsl_ermes(4, 1, "Over 100 iterations have occurred without convergence.  Convergence is assumed.");
*/
        imsl_ermes(IMSL_WARNING, IMSL_CHI_2_INV_CDF_CONVERGENCE);
L_30:
	chiin_v = (x1 + x2) / F_TWO;

L_9000:
	E1POP("imsl_f_chi_squared_inverse_cdf","imsl_d_chi_squared_inverse_cdf");
	return (chiin_v);
}				/* end of function */
