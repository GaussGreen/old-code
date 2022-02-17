#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  ANORIN/DNORIN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 8, 1985

    Purpose:    Evaluate the inverse of the standard normal (Gaussian)
                distribution function.

    Usage:      ANORIN(P)

    Arguments:
       P      - Probability for which the inverse of the normal
                distribution function is to be evaluated.  (Input)
                P must be in the open interval (0.0, 1.0).
       ANORIN - Function value.  (Output)
                The probability that a standard normal random variable
                takes a value less than or equal to ANORIN is P.

    Keywords:   Percentage point; Error function; Probability
                distribution; Continuous random variables; Percentile;
                Fractile

    GAMS:       L5a2n; C8a

    Chapters:   STAT/LIBRARY Probability Distribution Functions and
                             Inverses
                SFUN/LIBRARY Probability Distribution Functions and
                             Inverses

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_normal_inverse_cdf(Mfloat p)
#else
Mfloat imsl_f_normal_inverse_cdf(p)
	Mfloat          p;
#endif
{
	Mfloat          a, anorin_v, b, f, sd, sn, w, wi, z, z2;
	static Mfloat   eps;
	static Mfloat   sqrt2 = 1.4142135623730950488;
	static Mfloat   e0 = -.5668422e-1;
	static Mfloat   e1 = .3937021;
	static Mfloat   e2 = -.3166501;
	static Mfloat   e3 = .6208963e-1;
	static Mfloat   g0 = .1851159e-3;
	static Mfloat   g1 = -.2028152e-2;
	static Mfloat   g2 = -.1498384;
	static Mfloat   g3 = .1078639e-1;
	static Mfloat   h0 = .9952975e-1;
	static Mfloat   h1 = .5211733;
	static Mfloat   h2 = -.6888301e-1;
	static Mfloat   a1 = -.5751703;
	static Mfloat   a2 = -1.896513;
	static Mfloat   a3 = -.5496261e-1;
	static Mfloat   b0 = -.113773;
	static Mfloat   b1 = -3.293474;
	static Mfloat   b2 = -2.374996;
	static Mfloat   b3 = -1.187515;
	static Mfloat   c0 = -.1146666;
	static Mfloat   c1 = -.1314774;
	static Mfloat   c2 = -.2368201;
	static Mfloat   c3 = .5073975e-1;
	static Mfloat   d0 = -44.27977;
	static Mfloat   d1 = 21.98546;
	static Mfloat   d2 = -7.586103;
	static Mfloat   f0 = -6.266786;
	static Mfloat   f1 = 4.666263;
	static Mfloat   f2 = -2.962883;


	anorin_v = imsl_amach(6);

	eps = imsl_amach(4);

	if (p <= F_ZERO || p >= F_ONE) {
#ifdef DOUBLE
		imsl_e1psh("imsl_d_normal_inverse_cdf");
#else
		imsl_e1psh("imsl_f_normal_inverse_cdf");
#endif
		imsl_e1str(1, p);

/*		imsl_ermes(5, 1, "The argument to the function must be greater than 0.0 and less than 1.0.  P = %(r1).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_ZERO_LT_P_LT_ONE);
#ifdef DOUBLE
		imsl_e1pop("imsl_d_normal_inverse_cdf");
#else
		imsl_e1pop("imsl_f_normal_inverse_cdf");
#endif
		goto L_9000;
	}
	if (p > eps && p < (F_ONE - F_THREE * eps)) {
		z = fabs(F_ONE - (p + p));
		if (z > .85) {
			a = F_ONE - z;
			b = z;
			/*
			 * Reduced argument is in (.85,1.0), obtain the
			 * transformed variable
			 */
			w = sqrt(-log(a + a * b));
			if (w >= 2.5) {
				if (w >= F_FOUR) {
					/*
					 * W greater than 4.0, approximate F
					 * by a rational function in 1.0/W
					 */
					wi = F_ONE / w;
					sn = ((g3 * wi + g2) * wi + g1) * wi;
					sd = ((wi + h2) * wi + h1) * wi + h0;
					f = w + w * (g0 + sn / sd);
				} else {
					/*
					 * W between 2.5 and 4.0, approximate
					 * F by a rational function in W
					 */
					sn = ((e3 * w + e2) * w + e1) * w;
					sd = ((w + f2) * w + f1) * w + f0;
					f = w + w * (e0 + sn / sd);
				}
			} else {
				/*
				 * W between 1.1322 and 2.5, approximate F by
				 * a rational function in W
				 */
				sn = ((c3 * w + c2) * w + c1) * w;
				sd = ((w + d2) * w + d1) * w + d0;
				f = w + w * (c0 + sn / sd);
			}
		} else {
			/*
			 * Z between 0.0 and .85, approximate F by a rational
			 * function in Z
			 */
			z2 = z * z;
			f = z + z * (b0 + a1 * z2 / (b1 + z2 + a2 / (b2 + z2 + a3 / (b3 +
								     z2))));
		}
		/*
		 * Determine the proper sign and form the solution
		 */
		if (p >= F_HALF) {
			anorin_v = sqrt2 * f;
		} else {
			anorin_v = -sqrt2 * f;
		}
	} else {
		/*
		 * Probability too small or too close to 1, so compute
		 * function value by a rational function in 1.0/W
		 */
		if (p >= F_HALF) {
			a = F_TWO * (F_ONE - p);
		} else {
			a = p + p;
		}
		w = sqrt(-log(a + (a - a * a)));
		wi = F_ONE / w;
		sn = ((g3 * wi + g2) * wi + g1) * wi;
		sd = ((wi + h2) * wi + h1) * wi + h0;
		f = w + w * (g0 + sn / sd);
		if (p >= F_HALF) {
			anorin_v = sqrt2 * f;
		} else {
			anorin_v = -sqrt2 * f;
		}
	}

L_9000:
	return (anorin_v);
}				/* end of function */
