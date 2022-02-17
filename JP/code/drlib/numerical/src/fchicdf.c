#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  CHIDF/DCHIDF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Evaluate the imsl_chi-squared distribution function.

    Usage:      CHIDF(CHSQ, DF)

    Arguments:
       CHSQ   - Argument for which the imsl_chi-squared distribution function
                is to be evaluated.  (Input)
       DF     - Number of degrees of freedom of the imsl_chi-squared
                distribution.  (Input)
                DF must be greater than or equal to 0.5.
       CHIDF  - Function value, the probability that a imsl_chi-squared random
                variable takes a value less than or equal to CHSQ.
                (Output)

    Remark:
       Informational errors
       Type Code
         1   1  The input argument, CHSQ, is less than zero.
         2   3  Using the normal distribution for large degrees of
                freedom, underflow would have occured.

    Keywords:   P-value; Probability integral; Probability distribution;
                Continuous random variables; Cumulative distribution
                function; CDF

    GAMS:       L5a1c; C7e

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
Mfloat imsl_f_chi_squared_cdf(Mfloat chsq, Mfloat df)
#else
Mfloat imsl_f_chi_squared_cdf(chsq, df)
	Mfloat          chsq, df;
#endif
{
	Mint            i;
	Mfloat          a, b, chidf_v, eps, gam, p, pt2, w, w1, x, z, z1;


#define FUNC(w,a,z)	(Mfloat) ((w)*exp( (a)*log( (z) ) - (z) ))

	E1PSH("imsl_f_chi_squared_cdf", "imsl_d_chi_squared_cdf");
	chidf_v = imsl_amach(6);
	/* Check CHSQ */
	if (chsq < F_ZERO) {
		imsl_e1str(1, chsq);

                imsl_e1stl(1, "chi_squared");
/*		imsl_ermes(1, 1, "Since CHSQ = %(r1) is less than zero, the distribution function is zero at CHSQ.");
*/
                imsl_ermes(IMSL_NOTE, IMSL_ARG_LESS_THAN_ZERO);
		chidf_v = F_ZERO;
		goto L_9000;
	}
	/*
	 * Check DF NOTE: In the documentation for the Library (9.2) version
	 * of this code (MDCH), DF had to be less than 66 for H36 machines.
	 */
	if (df < F_HALF) {
		imsl_e1str(1, df);

/*		imsl_ermes(5, 2, "The number of degrees of freedom of the distribution, DF = %(r1), must be at least equal to 0.5.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_DF_MUST_BE_GE_POINT_5);
		goto L_9000;
	}
	eps = imsl_amach(4);
	/*
	 * Set P=0. If CHSQ is less than or equal to 10.**(-12).
	 */
	if (chsq <= 1.0e-12) {
		p = F_ZERO;
	} else {
		if (df > 65.0) {
			/*
			 * NOTE:  If the upper bound of DF does not change,
			 * the following block of code is unnecessary.
			 * 
			 * Use normal distribution approximation for large
			 * degrees of freedom.
			 */
			if (chsq < F_TWO) {
				p = F_ZERO;
			} else {
				pt2 = F_TWO / F_NINE;
				x = (pow(chsq / df, F_ONE / F_THREE) - (F_ONE - pt2 / df)) /
					sqrt(pt2 / df);
				if (x > F_FIVE) {
					p = F_ONE;
				} else {
					if (x < -18.8055) {
						/*
						 * Warning error - underflow
						 * would have occurred.
						 */

/*						imsl_ermes(2, 3, "The normal distribution is used for large degrees of freedom.  However, it has produced underflow.  Therefore, the probability is set to 0.");
*/
                                                imsl_ermes(IMSL_ALERT,
						IMSL_NORMAL_UNDERFLOW);
						chidf_v = F_ZERO;
						goto L_9000;
					} else {
						p = imsl_f_normal_cdf(x);
					}
				}
			}
		} else {
			/*
			 * Initialization for calculation using incomplete
			 * GAMMA function.
			 */
			if (chsq > 200.0) {
				p = F_ONE;
			} else {
				a = F_HALF * df;
				z = F_HALF * chsq;
				gam = imsl_f_gamma(a);
				w = imsl_f_max(F_HALF*a, 13.0);
				if (z < w) {
					if (df > 25.0 && chsq < F_TWO) {
						p = F_ZERO;
					} else {
						/*
						 * Calculate using equation
						 * no. 6.5.29.
						 */
						w = F_ONE / (gam * a);
						w1 = w;
						for (i = 1; i <= 50; i++) {
							b = i;
							w1 = w1 * z / (a + b);
							if (w1 <= eps * w)
								goto L_20;
							w += w1;
						}
				L_20:
						p = FUNC(w, a, z);
					}
				} else {
					/*
					 * Calculate using equation no.
					 * 6.5.32.
					 */
					z1 = F_ONE / z;
					b = a - F_ONE;
					w1 = b * z1;
					w = F_ONE + w1;
					for (i = 2; i <= 50; i++) {
						b -= F_ONE;
						w1 *= b * z1;
						if (w1 <= eps * w)
							goto L_40;
						w += w1;
					}
			L_40:
					w = z1 * FUNC(w, a, z);
					p = F_ONE - w / gam;
				}
			}
		}
	}

	chidf_v = p;

L_9000:
	E1POP("imsl_f_chi_squared_cdf", "imsl_d_chi_squared_cdf");
	return (chidf_v);
#undef	FUNC
}				/* end of function */
