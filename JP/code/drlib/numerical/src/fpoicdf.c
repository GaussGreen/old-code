#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* 
  -----------------------------------------------------------------------
    IMSL Name:  POIDF/DPOIDF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 15, 1986

    Purpose:    Evaluate the Poisson distribution function.

    Usage:      POIDF(K, THETA)

    Arguments:
       K      - Argument for which the Poisson distribution function
                is to be evaluated.  (Input)
       THETA  - Mean of the Poisson distribution.  (Input)
                THETA must be positive.
       POIDF  - Function value, the probability that a Poisson random
                variable takes a value less than or equal to K.  (Output)

    Remark:
       Informational error
       Type Code
         1   1  The input argument, K, is less than zero.

    Keywords:   P-value; Cumulative distribution function; CDF;
                Probability distribution; Discrete random variables

    GAMS:       L5a1p

    Chapters:   STAT/LIBRARY Probability Distribution Functions and
                             Inverses
                SFUN/LIBRARY Probability Distribution Functions and
                             Inverses

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_poisson_cdf(Mint k, Mfloat theta)
#else
Mfloat imsl_f_poisson_cdf(k, theta)
	Mint            k;
	Mfloat          theta;
#endif
{
	Mint            icnt, j, jj, k1, kcnt;
	Mfloat          alnsml, eps, g, h, p1, p2, pe, poidf_v, sml, temp,
	                x, x2, y, y2;


	E1PSH("imsl_f_poisson_cdf","imsl_d_poisson_cdf");
	poidf_v = imsl_amach(6);
	/* Check THETA */
	if (theta <= F_ZERO) {
		imsl_e1str(1, theta);

/*		imsl_ermes(5, 2, "The mean of the Poisson distribution, THETA = %(r1), must be positive.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_THETA_MUST_BE_POSITIVE);
		goto L_9000;
	}
	/* Check K */
	if (k < 0) {
		imsl_e1sti(1, k);

/*		imsl_ermes(1, 1, "Since the input value, K = %(i1), is less than 0, the distribution function is set to 0. ");
*/
                imsl_ermes(IMSL_NOTE, IMSL_LESS_THAN_ZERO);
		poidf_v = F_ZERO;
		goto L_9000;
	}
	k1 = k + 1;
	eps = imsl_amach(4);
	sml = 2 * imsl_amach(1);
	alnsml = log(sml);
	/* Special case - LAMBDA = 0. */
	if (theta <= eps) {
		pe = F_ONE;
	} else {
		/*
		 * Initialization for forward calculation.
		 */
		x = theta;
		y = F_ONE;
		jj = 1;
		p1 = -theta;
		icnt = p1 / alnsml;
		p1 += -icnt * alnsml;
		p1 = exp(p1);
		/*
		 * Initialization for backward calculation.
		 */
		x2 = k;
		y2 = theta;
		g = x2 * log(y2);
		h = k1;
		h = imsl_f_log_gamma(h);
		p2 = -y2 + g - h;
		kcnt = p2 / alnsml;
		p2 += -kcnt * alnsml;
		p2 = exp(p2);
		g = F_ONE;
		h = F_ONE;
		if (icnt == 0)
			g = F_ONE - p1;
		if (kcnt == 0)
			h = F_ONE - p2;
		pe = F_ZERO;
		/*
		 * Determine at which end caluclation is occurring.
		 */
L_10:
		j = icnt - kcnt;
		/* Forward calculation. */
		if (j > F_ZERO || (j == F_ZERO && p1 <= p2)) {
			/*
			 * Scaling is unnecessary. Store individual term.
			 */
			if (icnt == 0)
				pe += p1;
			/* All terms are accounted for. */
			if (jj != k1) {
				/* Calculate next term recursively. */
				p1 = p1 * x / y;
				if (p1 >= h) {
					/* Rescale. */
					temp = p1 * sml;
					if (temp != F_ZERO) {
						p1 = temp;
						icnt -= 1;
					}
				}
				jj += 1;
				y += F_ONE;
				goto L_10;
			}
			/* Backward calculation. */
		} else {
			/*
			 * Scaling is unnecessary--(J.LT.0.0 .OR. (J.EQ.0.0
			 * .AND. P1.GT.P2)). Store individual term.
			 */
			if (kcnt == 0)
				pe += p2;
			/* All terms are accounted for. */
			if (jj != k1) {
				/* Calculate next term recursively. */
				p2 = p2 * x2 / y2;
				if (p2 >= g) {
					/* Rescale. */
					temp = p2 * sml;
					if (temp != F_ZERO) {
						p2 = temp;
						kcnt -= 1;
					}
				}
				k1 -= 1;
				x2 -= F_ONE;
				goto L_10;
			}
		}
	}

	if (pe > F_ONE)
		pe = F_ONE;
	poidf_v = pe;

L_9000:
	E1POP("imsl_f_poisson_cdf","imsl_d_poisson_cdf");
	return (poidf_v);
}				/* end of function */
