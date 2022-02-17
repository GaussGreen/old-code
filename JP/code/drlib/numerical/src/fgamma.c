#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef DOUBLE
#define   imsl_f_beta       imsl_d_beta
#endif

static void	PROTO(l_r9gaml,(Mfloat *xmin, Mfloat *xmax));

static Mfloat lv_xmin = 0.0;
static Mfloat lv_xmax = 0.0;

#define COUNT(A)    (sizeof(A)/sizeof(Mfloat))
	/*
	 * Series for GAMMA on the interval 0.0 to 1.00000D+00 
	 * with weighted error        5.79D-32
	 * log weighted error         31.24
	 * significant figures reqd.  30.00
	 * decimal places required    32.05
	 */
static Mfloat lv_gamcs[] = {
     .857119559098933142192006239994e-2,
     .441538132484100675719131577165e-2,
     .568504368159936337863266458879e-1,
     -.421983539641856050101250018662e-2,
     .132680818121246022058400679635e-2,
     -.189302452979888043252394702389e-3,
     .360692532744124525657808221723e-4,
     -.605676190446086421848554829037e-5,
     .105582954630228334473182350909e-5,
     -.181196736554238404829185589117e-6,
     .311772496471532227779025459317e-7,
     -.535421963901968714087408102435e-8,
     .919327551985958894688778682594e-9,
     -.157794128028833976176742327395e-9,
     .270798062293495454326654043309e-10,
     -.464681865382573014408166105893e-11,
     .797335019200741965646076717536e-12,
     -.136807820983091602579949917231e-12,
     .234731948656380065723347177169e-13,
     -.40274326149490669327665705347e-14,
     .691005174737210091213833697526e-15,
     -.118558450022199290705238712619e-15,
     .203414854249637395520102605193e-16,
     -.349005434171740584927401294911e-17,
     .598799385648530556713505106603e-18,
     -.102737805787222807449006977843e-18,
     .176270281606052982494275966075e-19,
     -.302432065373530626095877211204e-20,
     .518891466021839783971783355051e-21,
     -.890277084245657669244925160107e-22,
     .152747406849334260227459689131e-22,
     -.26207312561873629002573283328e-23,
     .449646404783053867033104657067e-24,
     -.771471273133687791170390152533e-25,
     .132363545312604403648657271467e-25,
     -.227099941294292881670231381333e-26,
     .389641899800399144932081664e-27,
     -.6685198115125953327792128e-28,
     .114699866314002438434761386667e-28,
     -.1967938586345134677295104e-29,
     .337644881658533809033489066667e-30,
     -.579307033578213578462549333333e-31
};

	/*
	 * SERIES FOR ALGM ON THE INTERVAL 0.0 TO  1.00000D-02 WITH WEIGHTED
	 * ERROR        1.28D-31 LOG WEIGHTED ERROR        30.89 SIGNIFICANT
	 * FIGURES REQD. 29.81 DECIMAL PLACES REQUIRED   31.48
	 */
static Mfloat lv_algmcs[] = {
     .166638948045186324720572965082e0,
     -.138494817606756384073298605914e-4,
     .981082564692472942615717154749e-8,
     -.180912947557249419426330626672e-10,
     .622109804189260522712601554342e-13,
     -.339961500541772194430333059967e-15,
     .268318199848269874895753884667e-17,
     -.28680424353346432841446224e-19,
     .396283706104643480367930666667e-21,
     -.6831888753985766870112e-23,
     .142922735594249814757333333333e-24,
     -.35475981581010705472e-26,
     .1025680058010470912e-27,
     -.34011022543167488e-29,
     .127664219563006293333333333333e-30
};

/* Structured by FOR_STRUCT, v0.2, on 01/08/90 at 17:50:08
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  GAMMA/DGAMMA (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 23, 1988

    Purpose:    Evaluate the complete gamma function.

    Usage:      GAMMA(X)

    Arguments:
       X      - Argument for which the complete gamma function is
                desired.  (Input)
       GAMMA  - Function value.  (Output)

    Remark:
       Informational errors
       Type Code
         2   1  The function underflows because X is too small.
         3   2  Result is accurate to less than one-half precision
                because X is too near a negative integer.

    Keyword:    Utilities

    GAMS:       C7a

    Chapters:   SFUN/LIBRARY Gamma Function and Related Functions
                STAT/LIBRARY Mathematical Support
                MATH/LIBRARY Utilities

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_gamma(Mfloat x)
#else
Mfloat imsl_f_gamma(x)
	Mfloat           x;
#endif
{
	Mint             i, n;
	Mfloat           gamma_v, sinpiy, xi, xn, y;
	static Mfloat    pi = 3.14159265358979323846264338328e0;
		      /* sq2pil is 0.5*log(2*PI) = log(sqrt(2*PI)) */
	static Mfloat    sq2pil = 0.918938533204672741780329736406e0;
	static Mint      n_gamcs = 0;
	static Mfloat    dxrel = 0.0;
	static Mfloat    xsml = 0.0;

	E1PSH("imsl_f_gamma", "imsl_d_gamma");

	gamma_v = imsl_amach(6);

	if (lv_xmax == F_ZERO) l_r9gaml(&lv_xmin, &lv_xmax);
	if (n_gamcs == 0) {
		/*
		 * INITIALIZE. FIND LEGAL BOUNDS FOR X, AND DETERMINE THE
		 * NUMBER OF TERMS IN THE SERIES REQUIRED TO ATTAIN AN
		 * ACCURACY TEN TIMES BETTER THAN MACHINE PRECISION.
		 */
		n_gamcs = imsl_inits(lv_gamcs, COUNT(lv_gamcs), 0.1 * imsl_amach(3));
		dxrel = sqrt(imsl_amach(4));
		xsml = exp(imsl_f_max(log(imsl_amach(1)), -log(imsl_amach(2)))+0.01);
	}
	/* START EVALUATING GAMMA(X) */
	y = fabs(x);

	if (y > F_TEN)
		goto L_40;
	/*
	 * COMPUTE GAMMA(X) FOR ABS(X).LE.10.0. REDUCE INTERVAL AND FIND
	 * GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
	 */
	n = x;
	if (x < F_ZERO)
		n -= 1;
	xn = n;
	y = x - xn;
	n -= 1;
	gamma_v = 0.9375 + imsl_csevl(F_TWO * y - F_ONE, lv_gamcs, n_gamcs);
	if (n == 0)
		goto L_9000;

	if (n > 0)
		goto L_20;
	/* COMPUTE GAMMA(X) FOR X .LT. 1. */
	n = -n;

	if (x == F_ZERO) {
		    /* The argument for the function can not be zero. */
		imsl_ermes(IMSL_TERMINAL, IMSL_ARG_ZERO);
		gamma_v = imsl_amach(6);
		goto L_9000;
	}
	/*
	 * PREVENT UNDERFLOW IN COMPARISON ON CDC
	 */
	if (1.0e20 * y < 1.0e20 * xsml) {
		    /* The function overflows because X = %(r1) is  */
		    /* too close to zero.			    */
		imsl_e1str(1, x);
		imsl_ermes(IMSL_FATAL, IMSL_ZERO_ARG_OVERFLOW);
		gamma_v = imsl_amach(6);
		goto L_9000;
	}
	xn = n - 2;
	if (x < F_ZERO && x + xn == F_ZERO) {
		imsl_e1str(1, x);
		    /* The argument for the function can not be a   */
		    /* negative integer. Argument X = %(r1).	    */
		imsl_ermes(IMSL_TERMINAL, IMSL_NEGATIVE_INTEGER);
		gamma_v = imsl_amach(6);
		goto L_9000;
	}
	if (x < -F_HALF && fabs((x - ((Mint)(x - F_HALF))) / x) < dxrel) {
		    /* The result is accurate to less than one half	*/
		    /* precision because X = %(r1) is too close to a	*/
		    /* negative integer.				*/
		imsl_e1str(1, x);
		imsl_ermes(IMSL_WARNING, IMSL_NEAR_NEG_INT_WARN);
	}
	xi = F_ZERO;
	for (i = 1; i <= n; i++) {
		gamma_v /= x + xi;
		xi += F_ONE;
	}
	goto L_9000;
	/* GAMMA(X) FOR X .GE. 2.0 */
L_20:
	xi = F_ONE;
	for (i = 1; i <= n; i++) {
		gamma_v *= y + xi;
		xi += F_ONE;
	}
	goto L_9000;
	/*
	 * GAMMA(X) FOR ABS(X) .GT. 10.0. RECALL Y = ABS(X).
	 */
L_40:
	if (x > lv_xmax) {
		    /* The function overflows because X = %(r1) is  */
		    /* greater than %(R2). */
		imsl_e1str(1, x);
		imsl_e1str(2, lv_xmax);
		imsl_ermes(IMSL_FATAL, IMSL_LARGE_ARG_OVERFLOW);
		gamma_v = imsl_amach(6);
		goto L_9000;
	}
	gamma_v = F_ZERO;
	if (x < lv_xmin) {
		    /* The function underflows because X = %(r1) is */
		    /* less than %(r2). The result is set to zero.  */
		imsl_e1str(1, x);
		imsl_e1str(2, lv_xmin);
		imsl_ermes(IMSL_ALERT, IMSL_SMALL_ARG_UNDERFLOW);
		goto L_9000;
	}
	gamma_v = exp((y - F_HALF) * log(y) - y + sq2pil + imsl_r9lgmc(y));
	if (x > F_ZERO)
		goto L_9000;

	if (fabs((x - ((Mint)(x - F_HALF))) / x) < dxrel) {
		imsl_e1str(1, x);
		    /* The result is accurate to less than one half	*/
		    /* precision because X = %(r1) is too close to a	*/
		    /* negative integer.				*/
		imsl_ermes(IMSL_WARNING, IMSL_NEAR_NEG_INT_WARN);
	}
	sinpiy = sin(pi * y);
	if (sinpiy == F_ZERO) {
		    /* The argument for the function can not be a	*/
		    /* negative integer. Argument X = %(r1).		*/
		imsl_e1str(1, x);
		imsl_ermes(IMSL_FATAL, IMSL_NEAR_NEG_INT_FATAL);
		gamma_v = imsl_amach(6);
		goto L_9000;
	}
	gamma_v = -pi / (y * sinpiy * gamma_v);

L_9000:
	E1POP("imsl_f_gamma", "imsl_d_gamma");
	return (gamma_v);
}				/* end of function */


/* Structured by FOR_STRUCT, v0.2, on 01/08/90 at 17:51:55
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  R9GAML (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1984

    Purpose:    Compute the underflow and overflow limits for the gamma
                function and several related functions.

    Usage:      CALL R9GAML (XMIN, XMAX)

    Arguments:
       XMIN   - Underflow limit.  (Output)
       XMAX   - Overflow limit.  (Output)

    Remark:
       XMIN and XMAX are not the only bounds for the function, but they
       are the only nontrivial ones to calculate.

    GAMS:       C7a

    Chapter:    SFUN/LIBRARY Gamma Function and Related Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r9gaml(Mfloat *xmin, Mfloat *xmax)
#else
static void l_r9gaml(xmin, xmax)
	Mfloat          *xmin, *xmax;
#endif
{
	Mint             i;
	Mfloat           alnbig, alnsml, xln, xold;

	imsl_e1psh("l_r9gaml");
	*xmin = imsl_amach(6);
	*xmax = imsl_amach(6);

	alnsml = log(imsl_amach(1));
	*xmin = -alnsml;
	for (i = 1; i <= 10; i++) {
		xold = *xmin;
		xln = log(*xmin);
		*xmin += -*xmin * ((*xmin + F_HALF) * xln - *xmin - 0.2258 + alnsml) /
			(*xmin * xln + F_HALF);
		if (fabs(*xmin - xold) < 0.005)
			goto L_20;
	}
	imsl_ermes(IMSL_FATAL, IMSL_CANNOT_FIND_XMIN);
	*xmin = imsl_amach(6);
	goto L_9000;

L_20:
	*xmin = -*xmin + 0.01;

	alnbig = log(imsl_amach(2));
	*xmax = alnbig;
	for (i = 1; i <= 10; i++) {
		xold = *xmax;
		xln = log(*xmax);
		*xmax += -*xmax * ((*xmax - F_HALF) * xln - *xmax + 0.9189 - alnbig) /
			(*xmax * xln - F_HALF);
		if (fabs(*xmax - xold) < 0.005)
			goto L_40;
	}
	imsl_ermes(IMSL_FATAL, IMSL_CANNOT_FIND_XMAX);
	*xmax = imsl_amach(6);
	goto L_9000;

L_40:
	*xmax -= 0.01;
	*xmin = imsl_f_max(*xmin, -*xmax + F_ONE);

L_9000:
	imsl_e1pop("l_r9gaml");
	return;
}				/* end of function */


/*
  -----------------------------------------------------------------------
    IMSL Name:  R9LGMC (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1984

    Purpose:    Evaluate the log gamma_c correction term for argument
                values greater than or equal to 10.0.

    Usage:      R9LGMC(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       R9LGMC - Function value.  (Output)

    Remarks:
    1. Informational error:
       Type Code
         2   1  The function underflows because X is too large.

    2. R9LGMC calculates the log gamma_c correction factor such that
       ALOG (GAMMA(X)) = ALOG(SQRT(2*PI))+(X-.5)*ALOG(X)-X+R9LGMC(X).

    Chapter:    SFUN/LIBRARY Gamma Function and Related Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_r9lgmc(Mfloat x)
#else
Mfloat imsl_r9lgmc(x)
    Mfloat           x;
#endif
{
	Mfloat           r9lgmc_v;
	static Mint      nalgm = 0;
	static Mfloat    xbig = 0.0;
	static Mfloat    xmax = 0.0;


	imsl_e1psh("imsl_r9lgmc");
	r9lgmc_v = imsl_amach(6);

	if (nalgm == 0) {
		nalgm = imsl_inits(lv_algmcs, COUNT(lv_algmcs), imsl_amach(3));
		xbig = F_ONE / sqrt(imsl_amach(3));
		xmax = exp(imsl_f_min(log(imsl_amach(2) / 12.0), -log(12.0 * imsl_amach(1))));
	}
	if (x < F_TEN) {
		    /* The argument X = %(r1) must be greater than  */
		    /* or equal to %(R2).			    */
		imsl_e1str(1, x);
		imsl_e1str(2, F_TEN);
		imsl_ermes(IMSL_FATAL, IMSL_SMALL_ARGUMENT);
	} else if (x < xmax) {
		if (x >= xbig) {
			r9lgmc_v = F_ONE / (12.0 * x);
		} else {
			r9lgmc_v = imsl_csevl(F_TWO * pow(F_TEN / x, F_TWO) - F_ONE, lv_algmcs, nalgm) /
				x;
		}

	} else {
		r9lgmc_v = F_ZERO;
		    /* The function underflows because X = %(r1) is */
		    /* greater than %(r2).  The result is set to zero */
		imsl_e1str(1, x);
		imsl_e1str(2, xmax);
		imsl_ermes(IMSL_ALERT, IMSL_LARGE_ARG_UNDERFLOW);
	}
	imsl_e1pop("imsl_r9lgmc");
	return (r9lgmc_v);
}				/* end of function */



/* Structured by FOR_STRUCT, v0.2, on 09/11/90 at 10:57:46
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BETA/DBETA (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 23, 1988

    Purpose:    Evaluate the complete beta function.

    Usage:      BETA(A, B)

    Arguments:
       A      - First beta parameter.  (Input)
                A must be positive.
       B      - Second beta parameter.  (Input)
                B must be positive.
       BETA   - Function value.  (Output)

    Remark:
       Informational error
       Type Code
         2   1  The function underflows because A and/or B is too large.

    GAMS:       C7b

    Chapters:   SFUN/LIBRARY Gamma Function and Related Functions
                STAT/LIBRARY Mathematical Support

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_beta(Mfloat a, Mfloat b)
#else
Mfloat imsl_f_beta(a, b)
	Mfloat           a, b;
#endif
{
	Mfloat           beta_v;
	static Mfloat    alnsml = 0.0;

	E1PSH("imsl_f_beta", "imsl_d_beta");
	beta_v = imsl_amach(6);

	if (lv_xmax == F_ZERO) l_r9gaml(&lv_xmin, &lv_xmax);
	if (alnsml == F_ZERO)  alnsml = log(imsl_amach(1));

	if (a <= F_ZERO || b <= F_ZERO) {
                    /* Both A = %(r1) and B = %(r2) must be greater than zero. */
		imsl_e1str(1, a);
		imsl_e1str(2, b);
		imsl_ermes(IMSL_TERMINAL, IMSL_BETA_NEG_ARG);
	} else if (a + b < lv_xmax) {
		beta_v = imsl_f_gamma(a) * imsl_f_gamma(b) / imsl_f_gamma(a + b);
	} else {
		beta_v = imsl_f_log_beta(a, b);
		if (beta_v >= alnsml) {
			beta_v = exp(beta_v);
		} else {
                    /* The function underflows because A = %(r1) and/or 
                        B = %(r2) is too large. */
			beta_v = F_ZERO;
			imsl_e1str(1, a);
			imsl_e1str(2, b);
			imsl_ermes(IMSL_ALERT, IMSL_BETA_UNDERFLOW);
		}
	}
	E1POP("imsl_f_beta", "imsl_d_beta");
	return (beta_v);
}				/* end of function */

