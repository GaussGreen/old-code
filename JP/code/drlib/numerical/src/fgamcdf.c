#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* 
  -----------------------------------------------------------------------
    IMSL Name:  GAMDF/DGAMDF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1984

    Purpose:    Evaluate the imsl_gamma distribution function.

    Usage:      GAMDF(X, A)

    Arguments:
       X      - Argument for which the imsl_gamma distribution function is to
                be evaluated.  (Input)
       A      - The shape parameter of the imsl_gamma distribution.  (Input)
                This parameter must be positive.
       GAMDF  - Function value, the probability that a imsl_gamma random
                variable takes a value less than or equal to X.  (Output)

    Remark:
       Informational error
       Type Code
         1   1  The input argument, X, is less than zero.

    Keywords:   P-value; Probability integral; Erlang distribution;
                Probability distribution; Continuous random variables;
                Cumulative distribution function; CDF

    GAMS:       L5a1g; C7e

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
Mfloat imsl_f_gamma_cdf(Mfloat x, Mfloat a)
#else
Mfloat imsl_f_gamma_cdf(x, a)
	Mfloat          x, a;
#endif
{
	Mint            i;
	Mfloat          ax, big, cnt, cut, gamdf_v, p, pnlg, ratio, reduc,
	                xmin, y, ycnt, z;
	Mfloat		*v, *v1;
        char           _e0[8*sizeof(Mfloat)];
        v = (Mfloat *) _e0;
        v1 = (Mfloat *) (_e0 + 2*sizeof(Mfloat));

	E1PSH("imsl_f_gamma_cdf","imsl_d_gamma_cdf");

	gamdf_v = imsl_amach(6);
	/* Check X */
	if (x < F_ZERO) {
		imsl_e1str(1, x);
                imsl_e1stl(1,"x");
/*		imsl_ermes(1, 1, "Since X = %(r1) is less than zero, the distribution function is zero at X. ");
*/
                imsl_ermes(IMSL_NOTE, IMSL_ARG_LESS_THAN_ZERO);
		gamdf_v = F_ZERO;
		goto L_9000;
	} else if (x == F_ZERO) {
		gamdf_v = F_ZERO;
	} else {
		/* Check A */
		if (a <= F_ZERO) {
			imsl_e1str(1, a);
/*			imsl_ermes(5, 2, "The shape parameter of the imsl_gamma distribution must be positive.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_SHAPE_PARAMETER_NEGATIVE);
			goto L_9000;
		}
		if (log(x) >= log(imsl_amach(2)) / F_THREE) {
			if (a >= F_HALF * x) {
				imsl_e1str(1, x);
				imsl_e1str(2, a);

/*				imsl_ermes(5, 3, "Since X = %(r1) and A = %(r2) are so large, the algorithm would overflow. ");
*/
                                imsl_ermes(IMSL_FATAL,
				IMSL_X_AND_A_ARE_TOO_LARGE);
				goto L_9000;
			} else {
				gamdf_v = F_ONE;
				goto L_9000;
			}
		}
		xmin = log(imsl_amach(1)) + imsl_amach(4);
		if (x == F_ZERO) {
			p = F_ZERO;
		} else {
			/* Define LOG-GAMMA and initialize */
			pnlg = imsl_f_log_gamma(a);
			cnt = a * log(x);
			ycnt = x + pnlg;
			if ((cnt - ycnt) <= xmin) {
				ax = F_ZERO;
			} else {
				ax = exp(cnt - ycnt);
			}
			big = 1.0e30;
			cut = 1.0e-8;
			/* Choose algorithmic method */
			if ((x > F_ONE) && (x >= a)) {
				/* Continued fraction expansion */
				y = F_ONE - a;
				z = x + y + F_ONE;
				cnt = F_ZERO;
				v[0] = F_ONE;
				v[1] = x;
				v[2] = x + F_ONE;
				v[3] = z * x;
				p = v[2] / v[3];
		L_10:
				cnt += F_ONE;
				y += F_ONE;
				z += F_TWO;
				ycnt = y * cnt;
				v[4] = v1[0] * z - v[0] * ycnt;
				v[5] = v1[1] * z - v[1] * ycnt;
				if (v[5] == F_ZERO) {
					for (i = 1; i <= 4; i++) {
						v[i - 1] = v1[i - 1];
					}
					if (fabs(v[4]) < big && fabs(v[5]) < big)
						goto L_10;
					/*
					 * Scale terms down to prevent
					 * overflow
					 */
					for (i = 1; i <= 4; i++) {
						v[i - 1] /= big;
					}
					goto L_10;
				}
				ratio = v[4] / v[5];
				reduc = fabs(p - ratio);
				if (reduc <= ratio * cut * F_TEN && reduc <= cut) {
					p = F_ONE - p * ax;
				} else {
					p = ratio;
					for (i = 1; i <= 4; i++) {
						v[i - 1] = v1[i - 1];
					}
					if (fabs(v[4]) < big && fabs(v[5]) < big)
						goto L_10;
					/*
					 * Scale terms down to prevent
					 * overflow
					 */
					for (i = 1; i <= 4; i++) {
						v[i - 1] /= big;
					}
					goto L_10;
				}
			} else {
				/* Series expansion */
				ratio = a;
				cnt = F_ONE;
				p = F_ONE;
		L_60:
				ratio += F_ONE;
				cnt = cnt * x / ratio;
				p += cnt;
				if (cnt > cut)
					goto L_60;
				p = p * ax / a;
			}
		}

		gamdf_v = p;
	}
L_9000:
	E1POP("imsl_f_gamma_cdf","imsl_d_gamma_cdf");

	return (gamdf_v);
}				/* end of function */
