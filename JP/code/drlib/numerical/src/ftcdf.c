#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  TDF/DTDF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1984

    Purpose:    Evaluate the Student's t distribution function.

    Usage:      TDF(T, DF)

    Arguments:
       T      - Argument for which the Student's t distribution function
                is to be evaluated.  (Input)
       DF     - Degrees of freedom.  (Input)
                DF must be greater than or equal to 1.0.
       TDF    - Function value, the probability that a Student's t random
                variable takes a value less than or equal to the input T.
                (Output)

    Keywords:   P-value; Probability integral; Probability distribution;
                Continuous random variables; Cumulative distribution
                function; CDF

    GAMS:       L5a1t

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
Mfloat imsl_f_t_cdf(Mfloat t, Mfloat df)
#else
Mfloat imsl_f_t_cdf(t, df)
	Mfloat           t, df;
#endif
{
	Mint            n;
	Mfloat          a, an, b, q, tdf_v, temp, w, xj, y, z;
	static Mfloat   con1 = 0.63661977236758;



	E1PSH("imsl_f_t_cdf","imsl_d_t_cdf");
	tdf_v = imsl_amach(6);
	/* Check DF */
	if (df < F_ONE) {
		imsl_e1str(1, df);

/*		imsl_ermes(5, 1, "The input number of degrees of freedom, DF = %(r1), must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_DF_MUST_BE_AT_LEAST_ONE);  
		goto L_9000;
	}
	temp = t * t;
	if (df > temp) {
		temp = t;
		an = df;
		n = an;
		temp *= temp;
		y = temp / an;
		b = F_ONE + y;
		if ((an != n || an >= 20.0) || an > 200.0) {
			/* Asymptotic series for large AN */
			w = b - F_ONE;
			if (w != F_ZERO)
				y *= log(b) / w;
			a = an - F_HALF;
			b = 48.0 * a * a;
			y *= a;
			y = (((((-0.4 * y - 3.3) * y - 24.0) * y - 85.5) / (0.8 * (y * y) +
			    100.0 + b) + y + F_THREE) / b + F_ONE) * sqrt(y);
			if (y < 18.8125) {
				/*
				 * Overflow (or underflow?) could occur on
				 * some machines on call to ERFC
				 */
				q = imsl_f_erfc(y * sqrt(F_HALF));
			} else {
				q = F_ZERO;
			}
		} else {
			if (an < 20.0 && temp < F_FOUR) {
				/* Nested summation of *COSINE* series */
				y = sqrt(y);
				a = y;
				if (an == F_ONE)
					a = F_ZERO;
		L_10:
				an -= F_TWO;
				if (an > F_ONE) {
					a = (an - F_ONE) / (b * an) * a + y;
					goto L_10;
				}
				if (an == F_ZERO)
					a /= sqrt(b);
				if (an != F_ZERO)
					a = (atan(y) + a / b) * con1;
				q = F_ONE - a;
			} else {
				/*
				 * *TAIL* series expansion for large T-values
				 */
				a = F_ONE;
				y = an;
				xj = F_ZERO;
				z = F_ZERO;
		L_20:
				if (a != z) {
					xj += F_TWO;
					z = a;
					y = y * (xj - F_ONE) / (b * xj);
					a += y / (an + xj);
					goto L_20;
				}
		L_30:
				if (an > F_ONE && a >= 1.0e-30) {
					/*
					 * NOTE:  The conditional,
					 * A.GE.1.0E-30, included above is
					 * needed for the division by A below
					 * (overflow). According to the
					 * original basis deck, some machines
					 * require a less restriction on A or
					 * none at all.
					 */
					a *= (an - F_ONE) / (b * an);
					an -= F_TWO;
					goto L_30;
				}
				if (an != F_ZERO)
					a = sqrt(b) * con1 * a / b;
				q = a;
			}
		}
	} else {
		temp = df / (df + temp);
		a = F_HALF * df;
		b = F_HALF;
		q = imsl_f_beta_incomplete(temp, a, b);
	}

	if (t > F_ZERO) {
		tdf_v = F_ONE - F_HALF * q;
	} else {
		tdf_v = F_HALF * q;
	}

L_9000:
	E1POP("imsl_f_t_cdf","imsl_d_t_cdf");

	return (tdf_v);
}				/* end of function */
