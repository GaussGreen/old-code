#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  HYPDF/DHYPDF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 14, 1986

    Purpose:    Evaluate the hypergeometric distribution function.

    Usage:      HYPDF(K, N, M, L)

    Arguments:
       K      - Argument for which the hypergeometric distribution
                function is to be evaluated.  (Input)
       N      - Sample size.  (Input)
                N must be greater than zero and greater than or equal to
                K.
       M      - Number of defectives in the lot.  (Input)
       L      - Lot size.  (Input)
                L must be greater than or equal to N and M.
       HYPDF  - Function value, the probability that a hypergeometric
                random variable takes a value less than or equal to K.
                (Output)
                HYPDF is the probability that K or fewer defectives occur
                in a sample of size N drawn from a lot of size L that
                contains M defectives.

    Remark:
       Informational errors
       Type Code
         1   3  The input argument, K, is less than zero.
         1   4  The input argument, K, is greater than the sample size.

    Keywords:   P-value; Acceptance sampling; Probability distribution;
                Continuous random variables; Cumulative distribution
                function; CDF

    GAMS:       L5a1h

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
Mfloat imsl_f_hypergeometric_cdf(Mint k, Mint n, Mint m, Mint l)
#else
Mfloat imsl_f_hypergeometric_cdf(k, n, m, l)
	Mint            k, n, m, l;
#endif
{
	Mint            icnt, iflag, j, k0, mm, mm1, nmmax, nmmin;
	Mfloat          a, a1, aa, aj, al, anmmax, anmmin, b1, bb, hypdf_v,
	                sml, u, xval1, xval2;
	static Mint     ksave = 0;
	static Mint     nsave = 0;
	static Mint     msave = 0;
	static Mint     lsave = 0;
	static Mfloat   p = 0.0;
	static Mfloat   sp = 0.0;

/* note: the routine name is shortened to get around a error.c bug  (8/23/90)*/
	E1PSH("imsl_f_hyper_cdf","imsl_d_hyper_cdf");
	hypdf_v = imsl_amach(6);

	if (n <= 0 || m <= 0) {
		imsl_e1sti(1, n);
		imsl_e1sti(2, m);

/*		imsl_ermes(5, 1, "The sample size and the number of defectives must be greater than zero. Sample size = %(i1).  Number of defectives = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_ARGUMENT_GT_ZERO);                
	}
	if (l < n || l < m) {
		imsl_e1sti(1, l);
		imsl_e1sti(2, n);
		imsl_e1sti(3, m);

/*		imsl_ermes(5, 2, "The lot size must be greater than the sample size and the number of defectives in the lot.  Lot size = %(i1).  Sample size = %(i2).  Number of defectives in the lot = %(i3).");
*/
                imsl_ermes(IMSL_FATAL, IMSL_LOT_SIZE_TOO_SMALL);
	}
	if (imsl_n1rcd(0) > 0)
		goto L_9000;

	if (k < 0) {
		imsl_e1sti(1, k);
                imsl_e1stl(1, "k");

/*		imsl_ermes(1, 3, "Since the argument to the function is less than 0, the probability is set to 0. Argument = %(i1).");
*/
                imsl_ermes(IMSL_NOTE, IMSL_LESS_THAN_ZERO);
		hypdf_v = F_ZERO;
		goto L_9000;
	}
	if (k > n) {
		imsl_e1sti(1, k);
		imsl_e1sti(2, n);
/*		imsl_ermes(1, 4, "Since the argument to the function is greater than the sample size, the probability is set to 1. Argument = %(i1). Sample size = %(i2).");
*/
                imsl_ermes(IMSL_NOTE, IMSL_K_GREATER_THAN_N);
		hypdf_v = F_ONE;
		goto L_9000;
	}
	/*
	 * Check to see if any calculations are necessary, else initialize
	 * save variables.
	 */
	if (((k == ksave && n == nsave) && m == msave) && l == lsave)
		goto L_30;
	ksave = k;
	nsave = n;
	msave = m;
	lsave = l;
	/* Handle special cases */
	if (k > m) {
		hypdf_v = F_ONE;
		sp = F_ONE;
		goto L_9000;
	}
	if (n - k > l - m) {
		hypdf_v = F_ZERO;
		sp = F_ZERO;
		goto L_9000;
	}
	/* Initialization for recursive formula */
	p = F_ONE;
	al = l;
	nmmin = imsl_i_min(n, m);
	anmmin = nmmin;
	nmmax = imsl_i_max(n, m);
	anmmax = nmmax;
	/*
	 * Choose more efficient direction for computation of probabilities
	 */
	xval1 = k * (l + 2);
	xval2 = (m + 1) * (n + 1);

	if (xval1 <= xval2) {
		/*
		 * Initialization for calculation of P(X .LE. K)
		 */
		iflag = 0;
		aa = anmmax - anmmin;
		bb = al - anmmax - anmmin;
		k0 = n + (m - l);
		if (k0 < 0)
			k0 = 0;
		aj = k0;
		a1 = anmmin - aj;
		b1 = aj + F_ONE;
		mm1 = k - k0;
		a = al - anmmax + aj;
		mm = nmmin - k0;
	} else {
		/*
		 * Initialization for calculation of P(X .GT. K)
		 */
		iflag = 1;
		aa = al - anmmax - anmmin;
		bb = anmmax - anmmin;
		a1 = anmmin;
		b1 = F_ONE;
		mm1 = nmmin - k;
		mm = l - nmmax;
		a = al - anmmin;
		if (nmmin < mm) {
			mm = nmmin;
			a = anmmax;
		}
	}

	icnt = 0;
	sml = imsl_amach(1) * F_TEN;

	if (mm != 0) {
		/* GET PRELIMINARY PROBABILITIES */
		for (j = 1; j <= mm; j++) {
			u = a / al;
			/*
			 * DETECT AND CORRECT FOR POTENTIAL UNDERFLOW
			 */
			if (u < sml / p) {
				p /= sml;
				icnt += 1;
			}
			p *= u;
			a -= F_ONE;
			al -= F_ONE;
		}
	}
	sp = F_ZERO;
	if (mm1 != 0) {
		for (j = 1; j <= mm1; j++) {
			if (icnt == 0)
				sp += p;
			u = a1 * (a1 + aa) / (b1 * (b1 + bb));
			p *= u;
			if (p >= F_ONE) {
				/* INVERSE UNDERFLOW ADJUSTMENT */
				p *= sml;
				icnt -= 1;
			}
			a1 -= F_ONE;
			b1 += F_ONE;
		}
	}
	if (icnt != 0)
		p = F_ZERO;

	if (iflag == 0) {
		/* P(X .LE. K)=P(X .EQ. K)+P(X .LT. K) */
		sp += p;
	} else {
		sp = F_ONE - sp;
	}

L_30:
	hypdf_v = sp;

L_9000:
	E1POP("imsl_f_hyper_cdf","imsl_d_hyper_cdf");
	return (hypdf_v);
}				/* end of function */
