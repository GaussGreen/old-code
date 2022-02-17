#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static Mfloat l_alnrel(Mfloat);
#else
static Mfloat l_alnrel();
#endif
/* 
  -----------------------------------------------------------------------
    IMSL Name:  ALBETA/DLBETA (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 23, 1988

    Purpose:    Evaluate the natural logarithm of the complete imsl_beta
                function for positive arguments.

    Usage:      ALBETA(A, B)

    Arguments:
       A      - The first argument of the BETA function.  (Input)
       B      - The second argument of the BETA function.  (Input)
       ALBETA - Function value.  (Output)

    Remarks:
    1. ALBETA returns the natural log of
       BETA(A,B) = ALOG (GAMMA(A)*GAMMA(B))/GAMMA(A+B)

    2. Note than the natural log of BETA(A,B) equals the natural log
       of BETA(B,A).

    3. ALBETA (A,B) returns accurate results even when A or B is very
       small.

    4. The arguments, A and B, must be greater than 0.

    GAMS:       C7b

    Chapter:    SFUN/LIBRARY Gamma Function and Related Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_log_beta(Mfloat a, Mfloat b)
#else
Mfloat imsl_f_log_beta(a, b)
	Mfloat          a, b;
#endif
{
	Mfloat          albeta_v, corr, p, q, temp;
	static Mfloat   sq2pil = .918938533204672741780329736406e0;

	E1PSH("imsl_f_log_beta","imsl_d_log_beta");
	albeta_v = imsl_amach(6);

	p = imsl_f_min(a, b);
	q = imsl_f_max(a, b);

	if (p <= F_ZERO) {

/*		imsl_ermes(5, 5, "Both arguments for the function must be greater than zero.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BOTH_ARGS_ARE_LT_ZERO);
	} else if (p >= F_TEN) {
		/* P AND Q ARE BIG */
		corr = imsl_r9lgmc(p) + imsl_r9lgmc(q) - imsl_r9lgmc(p + q);
		/* CHECK FOR UNDERFLOW FROM R9LGMC */
		if (imsl_n1rcd(1) == 1)
			imsl_e1mes(0, 0, " ");
		temp = l_alnrel(-p / (p + q));
		albeta_v = -F_HALF * log(q) + sq2pil + corr + (p - F_HALF) * log(p /
							(p + q)) + q * temp;

	} else if (q >= F_TEN) {

		/* P IS SMALL, BUT Q IS BIG */
		corr = imsl_r9lgmc(q) - imsl_r9lgmc(p + q);
		/* CHECK FOR UNDERFLOW FROM R9LGMC */
		if (imsl_n1rcd(1) == 1)
			imsl_e1mes(0, 0, " ");
		albeta_v = imsl_f_log_gamma(p) + corr + p - p * log(p + q) + (q -
					   F_HALF) * l_alnrel(-p / (p + q));

	} else {
		/* P AND Q ARE SMALL */
		albeta_v = log(imsl_f_gamma(p) * (imsl_f_gamma(q) / imsl_f_gamma(p + q)));
	}

	E1POP("imsl_f_log_beta","imsl_d_log_beta");
	return (albeta_v);
}				/* end of function */
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  ALNREL/DLNREL (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 31, 1988

    Purpose:    Evaluate the natural logarithm of one plus the argument.

    Usage:      ALNREL(X)

    Arguments:
       X      - Argument for the function.  (Input)
       ALNREL - Function value.  (Output)

    Remarks:
    1. Informational error
       Type Code
         3   2  Result of ALNREL(X) is accurate to less than one half
                precision because X is too near -1.0.

    2. ALNREL evaluates the natural logarithm of (1+X) accurate in
       the sense of relative error even when X is very small. This
       routine (as opposed to ALOG) should be used to maintain relative
       accuracy whenever X is small and accurately known.

    Keyword:    Elementary functions

    GAMS:       C4b

    Chapter:    SFUN/LIBRARY Elementary Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_alnrel(Mfloat x)
#else
static Mfloat l_alnrel(x)
	Mfloat          x;
#endif
{
	Mfloat          alnrel_v;
	static Mint     nlnrel = 0;
	static Mfloat   xmin = 0.0;
	/*
	 * SERIES FOR ALNR ON THE INTERVAL -3.75000E-01 TO  3.75000E-01 WITH
	 * WEIGHTED ERROR        1.93E-17 LOG WEIGHTED ERROR        16.72
	 * SIGNIFICANT FIGURES REQD. 16.44 DECIMAL PLACES REQUIRED   17.40
	 */
	static Mfloat alnrcs[] = {
		.103786935627437698006862677191e1,
		-.133643015049089180987660415531e0,
		.194082491355205633579261993748e-1,
		-.301075511275357776903765377766e-2,
		.486946147971548500904563665091e-3,
		-.810548818931753560668099430086e-4,
		.137788477995595247829382514961e-4,
		-.238022108943589702513699929149e-5,
		.41640416213865183476391859902e-6,
		-.73595828378075994984266837032e-7,
		.13117611876241674949152294345e-7,
		-.235467093177424251366960923302e-8,
		.425227732760349977756380529626e-9,
		-.771908941348407968261081074933e-10,
		.140757464813590699092153564722e-10,
		-.257690720580246806275370786276e-11,
		.473424066662944218491543950059e-12,
		-.872490126747426417453012632927e-13,
		.161246149027405514657398331191e-13,
		-.298756520156657730067107924168e-14,
		.554807012090828879830413216973e-15,
		-.103246191582715695951413339619e-15,
		.192502392030498511778785032449e-16,
		-.359550734652651500111897078443e-17,
		.672645425378768578921945742268e-18,
		-.126026241687352192520824256376e-18,
		.236448844086062100449161589555e-19,
		-.444193770508079368988783891797e-20,
		.835465944640342590162412939947e-21,
		-.157315594164795625748992535211e-21,
		.296531287402474226861543697067e-22,
		-.559495834818159472921560132267e-23,
		.105663542688356810481872841387e-23,
		-.199724836806702045483149994667e-24,
		.37782977818839361421049856e-25,
		-.715315868890817403450381653333e-26,
		.135524884636742136465020245333e-26,
		-.256946730484875674300798293333e-27,
		.4874775606621694907645952e-28,
		-.925421125308497153211323733333e-29,
		.1757859784176023923326976e-29,
		-.334100266777310103513770666667e-30,
		.635339361802361873541802666667e-31};

	imsl_e1psh("l_alnrel");
	alnrel_v = imsl_amach(6);

	if (nlnrel == 0) {
#ifdef DOUBLE
		nlnrel = imsl_inits(alnrcs, 43, 0.1 * imsl_amach(3));
#else
		nlnrel = imsl_inits(alnrcs, 23, 0.1 * imsl_amach(3));
#endif
		xmin = -F_ONE + sqrt(imsl_amach(4));
	}
	if (x <= -F_ONE) {
		imsl_e1str(1, x);

/*		imsl_ermes(5, 5, "The argument X = %(r1) must be greater than -1.0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_X_IS_LESS_THAN_MINUS_1);

	} else {
		if (fabs(x) <= 0.375) {
			alnrel_v = x * (F_ONE - x * imsl_csevl(x / .375, alnrcs, nlnrel));
		} else {
			alnrel_v = log(F_ONE + x);
		}

		if (x < xmin) {
			imsl_e1str(1, x);
			imsl_e1str(2, xmin);

/*			imsl_ermes(3, 2, "The result is accurate to less than one half precision because X = %(r1) is too close to -1.0. X must be greater than %(r2).");
*/
                        imsl_ermes(IMSL_WARNING,
			IMSL_X_IS_TOO_CLOSE_TO_NEG_1);
		}
	}

	imsl_e1pop("l_alnrel");
	return (alnrel_v);
}				/* end of function */
