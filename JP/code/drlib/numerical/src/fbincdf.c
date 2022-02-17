#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  BINDF/DBINDF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 14, 1986

    Purpose:    Evaluate the binomial distribution function.

    Usage:      BINDF(K, N, P)

    Arguments:
       K      - Argument for which the binomial distribution function
                is to be evaluated.  (Input)
       N      - Number of Bernoulli trials.  (Input)
       P      - Probability of success on each trial.  (Input)
       BINDF  - Function value, the probability that a binomial random
                variable takes a value less than or equal to K.  (Output)
                BINDF is the probability that K or fewer successes occur
                in N independent Bernoulli trials each of which has a
                P probability of success.

    Remark:
       Informational errors
       Type Code
         1   3  The input argument, K, is less than zero.
         1   4  The input argument, K, is greater than the number of
                Bernoulli trials, N.

    Keywords:   P-value; Bernoulli distribution; Cumulative distribution
                function; CDF; Probability distribution; Discrete random
                variable

    GAMS:       L5a1b

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
Mfloat imsl_f_binomial_cdf(Mint k, Mint n, Mfloat p)
#else
Mfloat imsl_f_binomial_cdf(k, n, p)
	Mint            k, n;
	Mfloat          p;
#endif
{
	short int       scale;
	Mint            icnt, j, k1;
	Mfloat          alqn, alsml, bindf_v, p1, q1, qp, sml, xj, xn,
	                xx;
	Mfloat          prkcm;
	static Mint     ksave = 0;
	static Mint     nsave = 0;
	static Mfloat   psave = 0.0;
	static Mfloat   cdf = 0.0;
	static Mfloat   prk = 0.0;


        E1PSH("imsl_f_binomial_cdf","imsl_d_binomial_cdf");

	bindf_v = imsl_amach(6);

	if (n <= 0) {
		imsl_e1sti(1, n);

/*		imsl_ermes(5, 1, "N must be greater than zero, N = %(i1).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GREATER_THAN_ZERO);
	}
	if (p > F_ONE || p < F_ZERO) {
		imsl_e1str(1, p);

/*		imsl_ermes(5, 2, "P = %(r1), but P must be nonnegative and no greater than 1.0, since P is a probability.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_PROBABILITY_VALUE);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;

	if (k < 0) {
		imsl_e1sti(1, k);

/*		imsl_ermes(1, 3, "Since K = %(i1) is less than zero, the distribution function is set to zero.");
*/
                imsl_ermes(IMSL_NOTE, IMSL_LESS_THAN_ZERO);
		prkcm = F_ZERO;
		bindf_v = F_ZERO;
		goto L_9000;
	}
	if (k > n) {
		imsl_e1sti(1, k);
		imsl_e1sti(2, n);

/*		imsl_ermes(1, 4, "Since K = %(i1) is greater than N = %(i2), the distribution function is set to 1.0.");
*/
                imsl_ermes(IMSL_NOTE, IMSL_GREATER_THAN_N);
		prkcm = F_ZERO;
		bindf_v = F_ONE;
		goto L_9000;
	}
	/*
	 * Check to see if any computations are necessary, else initialize
	 * save variables.
	 */
	if ((k == ksave && n == nsave) && p == psave)
		goto L_20;
	ksave = k;
	nsave = n;
	psave = p;

	sml = F_TEN * imsl_amach(1);
	alsml = log(sml);
	prk = F_ZERO;
	cdf = F_ZERO;
	if (p == F_ZERO) {
		/* SPECIAL CASE FOR P = 0.0 */
		cdf = F_ONE;
		if (k != 0)
			goto L_20;
		prk = F_ONE;
		goto L_20;
	} else if (p == F_ONE) {
		/* SPECIAL CASE FOR P = 1.0 */
		if (k != n)
			goto L_20;
		prk = F_ONE;
		cdf = F_ONE;
		goto L_20;
	}
	p1 = p;
	q1 = F_ONE - p;
	k1 = k;
	xn = n;
	xx = xn * p;
	if (k > xx) {
		p1 = q1;
		q1 = p;
		k1 = n - k;
	}
	alqn = xn * log(q1);
	icnt = alqn / alsml;
	alqn += -icnt * alsml;
	prk = exp(alqn);
	if (k1 != 0) {
		qp = p1 / q1;
		xj = F_ZERO;
		xn += F_ONE;
		scale = 1;
		for (j = 1; j <= k1; j++) {
			if (scale)
				cdf += prk;
			scale = 1;
			xj += F_ONE;
			prk = prk * (qp * (xn - xj)) / xj;
			if (icnt > 0 && prk >= F_ONE) {
				cdf += prk;
				prk *= sml;
				cdf *= sml;
				icnt -= 1;
				scale = 0;
			}
		}
		if (scale)
			cdf += prk;
	} else if (icnt == 0) {
		cdf = prk;
	}
	if (icnt != 0)
		prk = F_ZERO;
	if (k > xx) {
		cdf = F_ONE - cdf + prk;
	}
	if (cdf > F_ONE)
		cdf = F_ONE;
	if (cdf < F_ZERO)
		cdf = F_ZERO;
	/* Set value in COMMON for BINPR. */
L_20:
	prkcm = prk;
	bindf_v = cdf;

L_9000:
        E1POP("imsl_f_binomial_cdf","imsl_d_binomial_cdf");
	return (bindf_v);
}				/* end of function */
