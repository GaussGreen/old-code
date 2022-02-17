#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static Mfloat	l_gamit(Mfloat , Mfloat);
static Mfloat	l_r9gmit(Mfloat, Mfloat, Mfloat, Mfloat, Mfloat);
static void	l_algams(Mfloat, Mfloat*, Mfloat*);
static Mfloat	l_r9lgic(Mfloat, Mfloat, Mfloat);
static Mfloat	l_r9lgit(Mfloat, Mfloat, Mfloat);
static Mfloat	l_gamr(Mfloat);
#else
static Mfloat	l_gamit();
static Mfloat	l_r9gmit();
static void 	l_algams();
static Mfloat 	l_r9lgic();
static Mfloat 	l_r9lgit();
static Mfloat 	l_gamr();
#endif


/* Structured by FOR_STRUCT, v0.2, on 09/12/90 at 16:33:16
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  GAMI/DGAMI (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 23, 1988

    Purpose:    Evaluate the incomplete imsl_gamma function.

    Usage:      GAMI(A, X)

    Arguments:
       A      - The integrand exponent parameter as per the remarks.
                (Input)
       X      - The upper limit of the integral definition of GAMI.
                (Input)
       GAMI   - Function value.  (Output)

    Remarks:
    1. GAMI equals the integral from T = 0 to X of
       EXP(-T) * T**(A-1.0) dT.

    2. GAMI is evaluated for positive values of A and nonnegative values
       of X. Because logarithmic variables are used, a slight
       deterioration of 2 or 3 digits of accuracy will occur when GAMI
       is very large or very small.

    GAMS:       C7

    Chapter:    SFUN/LIBRARY Gamma Function and Related Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_gamma_incomplete(Mfloat a, Mfloat x)
#else
Mfloat imsl_f_gamma_incomplete(a, x)
	Mfloat           a, x;
#endif
{
	Mfloat           factor, gami_v;


#ifdef DOUBLE
	imsl_e1psh("imsl_d_gamma_incomplete");
#else
	imsl_e1psh("imsl_f_gamma_incomplete");
#endif
	gami_v = imsl_amach(6);

	if (a <= F_ZERO) {
		imsl_e1str(1, a);

		imsl_ermes(IMSL_TERMINAL, IMSL_FIRST_ARG_LT_ZERO);
	}
	if (x < F_ZERO) {
		imsl_e1str(1, x);

		imsl_ermes(IMSL_TERMINAL, IMSL_SECOND_ARG_LT_ZERO);
	}
	if (imsl_n1rty(0) == 5)
		goto L_9000;

	if (x == F_ZERO) {
		gami_v = F_ZERO;

	} else {
		/*
		 * THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL
		 * OVERFLOW.
		 */
		factor = exp(imsl_f_log_gamma(a) + a * log(x));
		gami_v = factor * l_gamit(a, x);
	}

L_9000:
#ifdef DOUBLE
	imsl_e1pop("imsl_d_gamma_incomplete");
#else
	imsl_e1pop("imsl_f_gamma_incomplete");
#endif
	return (gami_v);
}				/* end of function */

/* Structured by FOR_STRUCT, v0.2, on 09/12/90 at 16:35:02
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  GAMIT/DGAMIT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 23, 1988

    Purpose:    Evaluate the Tricomi form of the incomplete imsl_gamma
                function.

    Usage:      GAMIT(A, X)

    Arguments:
       A      - The integrand exponent parameter as per the remarks.
                (Input)
       X      - The upper limit of the integral definition of GAMIT.
                (Input)
       GAMIT  - Function value.  (Output)

    Remarks:
    1. Informational error,
       Type Code
         3   2  Result of GAMIT(A,X) is accurate to less than one half
                precision because A is too close to a negative integer.

    2. The Tricomi incomplete imsl_gamma function GAMIT equals
       X**(-A)/GAMMA(A) * integral T = 0 to X of EXP(-T) * T**(A-1.0) dT
       and analytic continuation for A .LE. 0.0.

    3. GAMIT is evaluated for arbitrary REAL values of A and for
       nonnegative values of X (even though GAMIT is defined for
       X .LT. 0.0).

    4. A slight deterioration of 2 or 3 digits of accuracy will occur
       when GAMIT is very large or very small in absolute value because
       logarithmic variables are used. Also, if the parameter A is very
       close to a negative integer (but not quite a negative integer),
       there is a loss of accuracy which is reported if the result is
       less than half machine precision.

    5. References - Gautschi, W., (1979), A Computational Procedure for
                    Incomplete Gamma Functions, ACM TRANS. MATH.
                    SOFTWARE, 5, 466-481.
                  - Gautschi, W., (1979), Algorithm 542:  Incomplete
                    Gamma Functions, ACM TRANS. MATH. SOFTWARE, 5,
                    482-489.

    GAMS:       C7

    Chapter:    SFUN/LIBRARY Gamma Function and Related Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_gamit(Mfloat a, Mfloat x)
#else
static Mfloat l_gamit(a, x)
	Mfloat           a, x;
#endif
{
	Mfloat           aeps, ainta, algap1, alng,
	                alx, gamit_v, h, sga, sgngam, t;
	static Mfloat    alneps = 0.0;
	static Mfloat    sqeps = 0.0;
	static Mfloat    bot = 0.0;



	imsl_e1psh("l_gamit");
	gamit_v = imsl_amach(6);

	if (alneps == F_ZERO) {
		alneps = -log(imsl_amach(3));
		sqeps = sqrt((double) imsl_amach(4));
		bot = log(imsl_amach(1));
	}
	if (x < F_ZERO) {
		imsl_e1str(1, x);

		imsl_ermes(IMSL_TERMINAL, IMSL_SECOND_ARG_LT_ZERO);
		goto L_9000;
	}
	if (x != F_ZERO)
		alx = log(x);
	sga = F_ONE;
	if (a != F_ZERO)
		sga = sign(F_ONE, a);
/*	ainta = imsl_aint(ADR(_f0, a + 0.5 * sga)); */
	ainta = (Mint)(a + F_HALF * sga);
	aeps = a - ainta;

	if (x <= F_ZERO) {
		gamit_v = F_ZERO;
		if (ainta > F_ZERO || aeps != F_ZERO)
			gamit_v = l_gamr(a + F_ONE);

	} else if (x <= F_ONE) {
		if (a >= -F_HALF || aeps != F_ZERO) {
			l_algams(a + F_ONE, &algap1, &sgngam);
		}
		gamit_v = l_r9gmit(a, x, algap1, sgngam, alx);

	} else if (x <= a) {
		t = l_r9lgit(a, x, imsl_f_log_gamma(a + F_ONE));
		if (t < bot) {
			if (imsl_n1rty(1) < 5)
				imsl_e1mes(0, 0, " ");
		}
		gamit_v = exp(t);

	} else {
		alng = l_r9lgic(a, x, alx);
		/*
		 * EVALUATE GAMIT IN TERMS OF aLOG(GAMIC(A,X))
		 */
		h = F_ONE;
		if (aeps == F_ZERO && ainta <= F_ZERO)
			goto L_10;
		l_algams(a + F_ONE, &algap1, &sgngam);
		t = log(fabs(a)) + alng - algap1;
		if (t > alneps)
			goto L_20;
		if (t > -alneps)
			h = F_ONE - sga * sgngam * exp(t);
		if (fabs(h) <= sqeps) {
			imsl_e1str(1, a);

/*			imsl_ermes(3, 2, "The result is accurate to less than one half precision because X = %(r1) is too close to a negative
integer.");
*/
                        imsl_ermes(IMSL_WARNING, IMSL_NEAR_NEG_INT_WARN);
		}
L_10:
		t = -a * alx + log(fabs(h));
		if (t < bot) {
			if (imsl_n1rty(1) < 5)
				imsl_e1mes(0, 0, " ");
		}
		gamit_v = sign(exp(t), h);
		goto L_9000;

L_20:
		t += -a * alx;
		if (t < bot) {
			if (imsl_n1rty(1) < 5)
				imsl_e1mes(0, 0, " ");
		}
		gamit_v = -sga * sgngam * exp(t);
	}

L_9000:
	imsl_e1pop("l_gamit");
	return (gamit_v);
}				/* end of function */



/* Structured by FOR_STRUCT, v0.2, on 09/13/90 at 08:58:58
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  R9GMIT

    Computer:   FORC/SINGLE

    Revised:    January 1, 1984

    Purpose:    Evaluate the Tricomi incomplete imsl_gamma function for small
                arguments.

    Usage:      R9GMIT (A, X, ALGAP1, SGNGAM, ALX)

    Arguments:
       A      - The integrand exponent parameter as per the remarks.
                (Input)
       X      - The upper limit of the integral definition of R9GMIT.
                (Input)
       ALGAP1 - Log of GAMMA(A+1). (Input)
       SGNGAM - Sign of GAMMA(A+1). (Input)
       ALX    - Log of X. (Input)
       R9GMIT - Function value. (Output)
    Remarks:
    1. R9GMIT returns X**(-A)/GAMMA(A) * integral from 0 to X of
       EXP(-T) * T**(A-1.0) dT.

    Copyright:  1984 by IMSL, Inc. All rights reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* ALX is not used here, but leave in the calling sequence. */
#ifdef ANSI
static Mfloat l_r9gmit(Mfloat a, Mfloat x, Mfloat algap1, Mfloat sgngam,
			Mfloat alx)
#else
static Mfloat l_r9gmit(a, x, algap1, sgngam, alx)
	Mfloat           a, x, algap1, sgngam, alx;
#endif
{
	Mint             k, m, ma;
	Mfloat           ae, aeps, alg2, algs, fk, r9gmit_v, s, sgng2, t,
	                te;
	static Mfloat    eps = 0.0;
	static Mfloat    bot = 0.0;



	imsl_e1psh("R9GMIT");
	r9gmit_v = imsl_amach(6);

	if (eps == F_ZERO) {
		eps = F_HALF * imsl_amach(3);
		bot = log(imsl_amach(1));
	}
	if (x <= F_ZERO) {
		imsl_e1str(1, x);

/*		imsl_ermes(5, 5, "X = %(r1) must be greater than zero.");*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NON_POSITIVE_ARGUMENT);
		goto L_9000;
	}
	if (a >= F_ZERO) {
		ma = a + F_HALF;
	} else {
		ma = a - F_HALF;
	}
	aeps = a - (Mfloat) (ma);

	if (a >= -F_HALF) {
		ae = a;
	} else {
		ae = aeps;
	}

	t = F_ONE;
	te = ae;
	s = t;
	for (k = 1; k <= 200; k++) {
		fk = k;
		te = -x * te / fk;
		t = te / (ae + fk);
		s += t;
		if (fabs(t) < eps * fabs(s))
			goto L_20;
	}

/*	imsl_ermes(5, 6, "The function did not converge in 200 terms of Taylor series.");
*/
        imsl_ermes(IMSL_FATAL, IMSL_NO_CONV_200_TS_TERMS);     
 	goto L_9000;

L_20:
	if (a >= -F_HALF) {
		algs = -algap1 + log(s);
		goto L_50;
	}
	algs = -imsl_f_log_gamma(F_ONE + aeps) + log(s);
	s = F_ONE;
	m = -ma - 1;
	if (m == 0)
		goto L_40;
	t = F_ONE;
	for (k = 1; k <= m; k++) {
		t = x * t / (aeps - (Mfloat) (m + 1 - k));
		s += t;
		if (fabs(t) < eps * fabs(s))
			goto L_40;
	}

L_40:
	r9gmit_v = F_ZERO;
	algs += -(Mfloat) (ma) * log(x);
	if (s == F_ZERO || aeps == F_ZERO)
		goto L_50;

	sgng2 = sgngam * sign(F_ONE, s);
	alg2 = -x - algap1 + log(fabs(s));

	if (alg2 > bot)
		r9gmit_v = sgng2 * exp(alg2);
	if (algs > bot)
		r9gmit_v += exp(algs);
	goto L_9000;

L_50:
	r9gmit_v = exp(algs);

L_9000:
	imsl_e1pop("R9GMIT");
	return (r9gmit_v);
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  ALGAMS/DLGAMS (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    November 23, 1988

    Purpose:    Return the logarithm of the absolute value of the imsl_gamma
                function and the sign of imsl_gamma.

    Usage:      CALL ALGAMS (X, ALGM, S)

    Arguments:
       X      - Argument for which the logarithm of the absolute value
                of the imsl_gamma function is desired.  (Input)
       ALGM   - Result of the calculation.  (Output)
       S      - Sign of imsl_gamma(X).  (Output)
                If imsl_gamma(X) is greater than or equal to zero, S = 1.0.
                If imsl_gamma(X) is less than zero, S = -1.0.

    Remark:
       Informational error,
       Type Code
         3   2  Result of ALGAMS is accurate to less than one half
                precision because X is too near a negative integer.

    GAMS:       C7

    Chapter:    SFUN/LIBRARY Gamma Function and Related Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_algams(Mfloat x, Mfloat *algm, Mfloat *s)
#else
static void l_algams(x, algm, s)
	Mfloat           x, *algm, *s;
#endif
{
	Mint             int_;


	imsl_e1psh("aLGAMS");
	*algm = imsl_amach(6);
	*s = imsl_amach(6);

	*algm = imsl_f_log_gamma(x);

	if (x > F_ZERO) {
		*s = F_ONE;
	} else {
		int_ = fmod(-(Mint)(x), F_TWO) + 0.1;
		if (int_ == 0) {
			*s = -F_ONE;
		} else {
			*s = F_ONE;
		}
	}

	imsl_e1pop("aLGAMS");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  R9LGIC

    Computer:   FORC/SINGLE

    Revised:    January 1, 1984

    Purpose:    Evaluate the log complementary incomplete imsl_gamma function
                for large X and for A less than or equal to X.

    Usage:      R9LGIC (A, X, ALX)

    Arguments:
       A      - The integrand exponent parameter as per the remarks.
                (Input)
       X      - The lower limit of the integral definition of R9LGIC.
                (Input)
       ALX    - The log of X. (Input)
       R9LGIC - Function value. (Output)
    Remarks:
    1. R9LGIC returns the log of the integral from X to infinity of
       EXP(-T) * T**(A-1.0) dT.

    Copyright:  1984 by IMSL, Inc. All rights reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_r9lgic(Mfloat a, Mfloat x, Mfloat alx)
#else
static Mfloat l_r9lgic(a, x, alx)
	Mfloat           a, x, alx;
#endif
{
	Mint             k;
	Mfloat           fk, p, r, r9lgic_v, s, t, xma, xpa;
	static Mfloat    eps = 0.0;



	imsl_e1psh("R9LGIC");
	r9lgic_v = imsl_amach(6);

	if (eps == F_ZERO)
		eps = F_HALF * imsl_amach(3);

	xpa = x + F_ONE - a;
	xma = x - F_ONE - a;

	r = F_ZERO;
	p = F_ONE;
	s = p;
	for (k = 1; k <= 200; k++) {
		fk = k;
		t = fk * (a - fk) * (F_ONE + r);
		r = -t / ((xma + F_TWO * fk) * (xpa + F_TWO * fk) + t);
		p *= r;
		s += p;
		if (fabs(p) < eps * s)
			goto L_20;
	}

/*	imsl_ermes(5, 5, "The function did not converge in 200 terms of continued fraction.");
*/
        imsl_ermes(IMSL_FATAL, IMSL_NO_CONV_200_CF_TERMS);
	goto L_9000;

L_20:
	r9lgic_v = a * alx - x + log(s / xpa);

L_9000:
	imsl_e1pop("R9LGIC");
	return (r9lgic_v);
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  R9LGIT

    Computer:   FORC/SINGLE

    Revised:    January 1, 1984

    Purpose:    Evaluate the log of the Tricomi incomplete imsl_gamma function
                with the Perron continued fraction for large X and for
                A greater than or equal to X.

    Usage:      R9LGIT (A, X, ALGAP1)

    Arguments:
       A      - The integrand exponent parameter as per the remarks.
                (Input)
       X      - The upper limit of the integral definition of R9LGIT.
                (Input)
       ALGAP1 - Log of GAMMA(A+1). (Input)
       R9LGIT - Function value. (Output)
    Remarks:
    1. R9LGIT returns the log of X**(-A)/GAMMA(A) * integral from 0 to X
       of EXP(-T) * T**(A-1.0) dT.
    2. Informational errors,
       Type Code Meaning
         3   2   Result of R9LGIT is accurate to less than one half
                 precision.
       Users wishing to override the print/stop attributes associated
       with error messages issued by this routine are referred to the
       Error Handling section of the introduction to the manual.

    Copyright:  1984 by IMSL, Inc. All rights reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_r9lgit(Mfloat a, Mfloat x, Mfloat algap1)
#else
static Mfloat l_r9lgit(a, x, algap1)
	Mfloat           a, x, algap1;
#endif
{
	Mint             k;
	Mfloat           a1x, ax, fk, hstar, p, r, r9lgit_v, s, t;
	static Mfloat    eps = 0.0;
	static Mfloat    sqeps = 0.0;



	imsl_e1psh("R9LGIT");
	r9lgit_v = imsl_amach(6);

	if (eps == F_ZERO) {
		eps = F_HALF * imsl_amach(3);
		sqeps = sqrt(imsl_amach(4));
	}
	if (x <= F_ZERO || a < x) {
		imsl_e1str(1, x);
		imsl_e1str(2, a);

/*		imsl_ermes(5, 5, "X = %(r1) must be greater than zero and less than or equal to A = %(r2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_ZERO_LT_X_LE_A);
		goto L_9000;
	}
	ax = a + x;
	a1x = ax + F_ONE;
	r = F_ZERO;
	p = F_ONE;
	s = p;
	for (k = 1; k <= 200; k++) {
		fk = k;
		t = (a + fk) * x * (F_ONE + r);
		r = t / ((ax + fk) * (a1x + fk) - t);
		p *= r;
		s += p;
		if (fabs(p) < eps * s)
			goto L_20;
	}

/*	imsl_ermes(5, 6, "The function did not converge in 200 terms of continued fraction.");
*/
        imsl_ermes(IMSL_FATAL, IMSL_NO_CONV_200_CF_TERMS);
	goto L_9000;

L_20:
	hstar = F_ONE - x * s / a1x;
	if (hstar < sqeps) {

/*		imsl_ermes(3, 2, "The result is accurate to less than one half precision.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_LT_HALF_ACCURATE);
	}
	r9lgit_v = -x - algap1 - log(hstar);

L_9000:
	imsl_e1pop("R9LGIT");
	return (r9lgit_v);
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  GAMR/DGAMR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 23, 1988

    Purpose:    Evaluate the reciprocal imsl_gamma function.

    Usage:      GAMR(X)

    Arguments:
       X      - Argument for which the reciprocal imsl_gamma function is
                desired.  (Input)
       GAMR   - Function value.  (Output)

    GAMS:       C7

    Chapter:    SFUN/LIBRARY Gamma Function and Related Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_gamr(Mfloat x)
#else
static Mfloat l_gamr(x)
	Mfloat           x;
#endif
{
	Mfloat           alngx, gamr_v, sgngx;


	imsl_e1psh("GAMR  ");
	gamr_v = imsl_amach(6);

	if (x <= F_ZERO && (Mfloat)((Mint)(x)) == x) {
		gamr_v = F_ZERO;

	} else if (fabs(x) <= F_TEN) {
		gamr_v = F_ONE / imsl_f_gamma(x);
		if (imsl_n1rty(1) == 3)
			imsl_e1mes(0, 0, " ");

	} else {
		l_algams(x, &alngx, &sgngx);
		if (imsl_n1rty(1) == 3)
			imsl_e1mes(0, 0, " ");
		gamr_v = sgngx * exp(-alngx);
	}

	imsl_e1pop("GAMR  ");
	return (gamr_v);
}				/* end of function */
