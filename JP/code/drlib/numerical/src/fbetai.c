#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  BETAI/DBETAI (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 23, 1988

    Purpose:    Evaluate the incomplete imsl_beta function ratio.

    Usage:      BETAI(X, PIN, QIN)

    Arguments:
       X      - Upper limit of integration.  (Input)
                X must be in the interval (0.0,1.0) inclusive.
       PIN    - First imsl_beta distribution parameter.  (Input)
       QIN    - Second imsl_beta distribution parameter.  (Input)
       BETAI  - Probability that a random variable from a imsl_beta
                distribution having parameters PIN and QIN will be
                less than or equal to X.  (Output)

    Remarks:
    1. PIN and QIN must both be positive.

    2. Based on Bosten and Battiste, Remark on Algorithm 179, Comm. ACM,
       V 17, P 153, (1974).

    GAMS:       C7b

    Chapter:    SFUN/LIBRARY Gamma Function and Related Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_beta_incomplete(Mfloat x, Mfloat pin, Mfloat qin)
#else
Mfloat imsl_f_beta_incomplete(x, pin, qin)
	Mfloat          x, pin, qin;
#endif
{
	Mint            i, ib, n;
	Mfloat          betai_v, c, finsum, p, p1, ps, q, term, xb, y;
	static Mfloat   eps = 0.0;
	static Mfloat   alneps = 0.0;
	static Mfloat   sml = 0.0;
	static Mfloat   alnsml = 0.0;

	E1PSH("imsl_f_beta_incomplete", "imsl_d_beta_incomplete");
	betai_v = imsl_amach(6);

	if (eps == F_ZERO) {
		eps = imsl_amach(3);
		alneps = log(eps);
		sml = imsl_amach(1);
		alnsml = 0.95 * log(sml);
	}
	if (x < F_ZERO || x > F_ONE) {
	           /* %(L1) must be between %(R1) and %(R2),
                      but %(L1) = %(R3). */
	        imsl_e1stl(1, "x");
		imsl_e1str(1, F_ZERO);
		imsl_e1str(2, F_ONE);
		imsl_e1str(3, x);
		imsl_ermes(IMSL_TERMINAL, IMSL_REAL_OUT_OF_RANGE);
	}
	if (pin <= F_ZERO || qin <= F_ZERO) {
	            /* Both A = %(r1) and B = %(r2) must be
                       greater than zero. */
		imsl_e1str(1, pin);
		imsl_e1str(2, qin);
		imsl_ermes(IMSL_TERMINAL, IMSL_BETA_NEG_ARG);
	}
	if (imsl_n1rty(0) == 5)
		goto L_9000;

	y = x;
	p = pin;
	q = qin;
	if (q <= p && x < 0.8)
		goto L_10;
	if (x < 0.2)
		goto L_10;
	y = F_ONE - y;
	p = qin;
	q = pin;

L_10:
	if ((p + q) * y / (p + F_ONE) < eps)
		goto L_70;
	/*
	 * EVALUATE THE INFINITE SUM FIRST. TERM WILL EQUAL Y**P/BETA(PS,P) *
	 * (1.0-PS)I * Y**I / FAC(I)
	 */
	ps = q - (Mint) q;
	if (ps == F_ZERO)
		ps = F_ONE;
	xb = p * log(y) - imsl_f_log_beta(ps, p) - log(p);
	betai_v = F_ZERO;
	if (xb < alnsml)
		goto L_30;

	betai_v = exp(xb);
	term = betai_v * p;
	if (ps == F_ONE)
		goto L_30;

	n = imsl_f_max(alneps / log(y), F_FOUR);

	for (i = 1; i <= n; i++) {
		term = term * ((float) (i) - ps) * y / (float) (i);
		betai_v += term / (p + (float) (i));
	}
	/* NOW EVALUATE THE FINITE SUM, MAYBE. */
L_30:
	if (q <= F_ONE)
		goto L_60;

	xb = p * log(y) + q * log(F_ONE - y) - imsl_f_log_beta(p, q) - log(q);
	ib = imsl_f_max(xb / alnsml, F_ZERO);
	term = exp(xb - (float) (ib) * alnsml);
	c = F_ONE / (F_ONE - y);
	p1 = q * c / (p + q - F_ONE);

	finsum = F_ZERO;
	n = q;
	if (q == (float) (n))
		n -= 1;
	for (i = 1; i <= n; i++) {
		if (p1 <= F_ONE && term / eps <= finsum)
			goto L_50;
		term = (q - (float) (i - 1)) * c * term / (p + q - (float) (i));

		if (term > F_ONE) {
			ib -= 1;
			term *= sml;
		}
		if (ib == 0)
			finsum += term;
	}

L_50:
	betai_v += finsum;
L_60:
	if (y != x || p != pin)
		betai_v = F_ONE - betai_v;
	betai_v = imsl_f_max(imsl_f_min(betai_v, F_ONE), F_ZERO);
	goto L_9000;

L_70:
	betai_v = F_ZERO;
	xb = p * log(imsl_f_max(y, sml)) - log(p) - imsl_f_log_beta(p, q);
	if (xb > alnsml && y != F_ZERO)
		betai_v = exp(xb);
	if (y != x || p != pin)
		betai_v = F_ONE - betai_v;

L_9000:
	E1POP("imsl_f_beta_incomplete", "imsl_d_beta_incomplete");
	return (betai_v);
}				/* end of function */
