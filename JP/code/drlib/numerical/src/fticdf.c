#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  TIN/DTIN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1984

    Purpose:    Evaluate the inverse of the Student's t distribution
                function.

    Usage:      TIN(P, DF)

    Arguments:
       P      - Probability for which the inverse of the Student's t
                distribution function is to be evaluated.  (Input)
                P must be in the open interval (0.0, 1.0).
       DF     - Degrees of freedom.  (Input)
                DF must be greater than or equal to 1.0.
       TIN    - Function value.  (Output)
                The probability that a Student's t random variable takes
                a value less than or equal to TIN is P.

    Remark:
       Informational error
       Type Code
         4   3  TIN is set to machine infinity since overflow would occur
                upon modifying the inverse value for the F distribution
                with the result obtained from the inverse BETA
                distribution.

    Keywords:   Percentage point; Probability distribution; Continuous
                random variables; Percentile; Fractile

    GAMS:       L5a2t

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
Mfloat imsl_f_t_inverse_cdf(Mfloat temp_p, Mfloat temp_df)
#else
Mfloat imsl_f_t_inverse_cdf(temp_p, temp_df)
	Mfloat          temp_p, temp_df;
#endif
{
	Mfloat          a, b, c, d, dfd, dfn, eps, f, pp, ptmp, q, 
	                tin_v, x, xx, y, p, df;
	static Mfloat   sdelp, smexe, xinf;
	static Mfloat   hpi = 1.570796;



	E1PSH("imsl_f_t_inverse","imsl_d_t_inverse");
	tin_v = imsl_amach(6);
        p = (Mfloat)temp_p;
        df = (Mfloat)temp_df;
	/* Check DF */
	if (df < F_ONE) {
		imsl_e1str(1, df);

/*		imsl_ermes(5, 1, "The input number of degrees of freedom DF = %(r1), must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_DF_MUST_BE_AT_LEAST_ONE);
		goto L_9000;
	}
	/* Check P */
	if (p <= F_ZERO || p >= F_ONE) {
		imsl_e1str(1, p);

/*		imsl_ermes(5, 2, "The input probability value, P = %(r1), must be in the open interval (0, 1).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_P_OUTSIDE_OPEN_INTERVAL);
		goto L_9000;
	}
	/*
	 * If P is close to 0.5 then set TIN = 0.0E0
	 */
	eps = imsl_amach(4);
	if (p > (F_HALF - eps) && p < (F_HALF + eps)) {
		tin_v = F_ZERO;
		goto L_9000;
	}
	/* Initializations */
	if (p < F_HALF) {
		q = F_TWO * p;
	} else {
		q = F_TWO * (F_ONE - p);
	}
	sdelp = F_ONE - imsl_amach(4);
	xinf = imsl_amach(7);
	smexe = log(imsl_amach(1));
	/*
	 * Exact integral for 2 degrees of freedom
	 */
	if (fabs(df - F_TWO) <= imsl_amach(4)) {
		x = sqrt(F_TWO / (q * (F_TWO - q)) - F_TWO);
	} else {
		/*
		 * Exact integral for 1 degree of freedom
		 */
		if (fabs(df - F_ONE) <= imsl_amach(4)) {
			a = q * hpi;
			x = cos(a) / sin(a);
		} else {
			if (df <= F_TWO) {
				/*
				 * DF is less than 2.0 Adjust DF for use in
				 * BETA distribution.
				 */
				dfn = F_ONE;
				dfd = df;
				/*
				 * Intialize a number slightly smaller than
				 * 1.0
				 */
				/* Reparameterize for BETA usage */
				a = F_HALF * dfn;
				b = F_HALF * dfd;
				/*
				 * Adjust P for use in the F distribution
				 */
				if (p < F_HALF) {
					ptmp = 1 - 2 * p;
				} else {
					ptmp = -(1 - 2 * p);
				}

				pp = F_ONE - ptmp;
				x = imsl_betin(pp, b, a);
				if (x == F_ZERO) {
					tin_v = xinf;

/*					imsl_ermes(4, 3, "Upon modifying the inverse value for the F distribution from the value obtained from IMSL inverse BETA distribution, BETIN, overflow would occur.  Therefore, the value is set to machine infinity.");
*/
                                        imsl_ermes(IMSL_FATAL,IMSL_F_INVERSE_OVERFLOW);
					goto L_9000;
				} else {
					f = (F_ONE / x - F_ONE) * dfd / dfn;
					x = sqrt(f);
				}
			} else {
				/* Expansion for DF greater than 2 */

				a = F_ONE / (df - F_HALF);
				b = 48.0e0 / (a * a);
				c = ((20700.0e0 * a / b - 98.0e0) * a - 16.0e0) * a + 96.36e0;
				d = ((94.5e0 / (b + c) - F_THREE) / b + F_ONE) * sqrt(a * hpi) *
					df;
				xx = d * q;
				y = pow(xx, F_TWO / df);
				if (y <= (a + 0.05e0)) {
					y = ((F_ONE / (((df + F_SIX) / (df * y) - 0.089e0 * d -
							0.822e0) * (df + F_TWO) * F_THREE) + F_HALF / (df + F_FOUR)) *
					     y - F_ONE) * (df + F_ONE) / (df + F_TWO) + F_ONE /
						y;
					x = sqrt(df * y);
				} else {
					/*
					 * Asymptotic inverse expansion about
					 * normal
					 */
					x = F_HALF * q;
					xx = imsl_f_normal_inverse_cdf(x);
					y = xx * xx;
					if (df < F_FIVE)
						c += 0.3e0 * (df - 4.5e0) * (xx + 0.6e0);
					c += (((0.05e0 * d * xx - F_FIVE) * xx - F_SEVEN) * xx -
					      F_TWO) * xx + b;
					y = (((((0.4e0 * y + 6.3e0) * y + 36.0e0) * y + 94.5e0) /
					   c - y - F_THREE) / b + F_ONE) * xx;
					y *= a * y;
					d = y;
					if (d <= 0.002e0)
						y += F_HALF * y * y;
					if (d > 0.002e0)
						y = exp(y) - F_ONE;
					x = sqrt(df * y);
				}
			}
		}
	}

	if (p < F_HALF)
		x = -x;
	tin_v = x;

L_9000:
	E1POP("imsl_f_t_inverse","imsl_d_t_inverse");

	return (tin_v);
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  BETIN/DBETIN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    May 30, 1990

    Purpose:    Evaluate the inverse of the imsl_beta distribution function.

    Usage:      BETIN(P, PIN, QIN)

    Arguments:
       P      - Probability for which the inverse of the imsl_beta
                distribution function is to be evaluated.  (Input)
                P must be in the open interval (0.0, 1.0).
       PIN    - First imsl_beta distribution parameter.  (Input)
                PIN must be positive.
       QIN    - Second imsl_beta distribution parameter.  (Input)
                QIN must be positive.
       BETIN  - Function value.  (Output)
                The probability that a imsl_beta random variable takes a value
                less than or equal to BETIN is P.

    Remark:
       Informational error
       Type Code
         3   1  The value for the inverse Beta distribution could not be
                found in 100 iterations. The best approximation is used.

    Keywords:   Percentage point; Probability distribution; Continuous
                random variables; Percentile; Fractile

    GAMS:       L5a2b; C7b

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
Mfloat imsl_betin(Mfloat p, Mfloat pin, Mfloat qin)
#else
Mfloat imsl_betin(p, pin, qin)
	Mfloat          p, pin, qin;
#endif
{
	Mint            ic, nc;
	Mfloat          a, afn, alnsml, b, betin_v, c, dtemp, eps, fcs,
	                fn, fxl, p1, q0, qx, sig, temp, x, xc, xl,
	                xr, xrmxl, xt, zi, zz;

	/*
	 * NOTE:  Upper bounds and lower bounds of P, PIN, QIN in relation to
	 * one another may have to be investigated when testing is done.  (To
	 * avoid computational overflow)
	 */
	E1PSH("imsl_betin","imsl_dbetin");

	betin_v = imsl_amach(6);
	/* Check PIN. */
	if (pin <= F_ZERO) {
		imsl_e1str(1, pin);

/*		imsl_ermes(5, 1, "The first parameter of the imsl_beta distribution, PIN = %(r1), must be positive.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_PIN_MUST_BE_POSITIVE);
		goto L_9000;
	}
	/* Check QIN. */
	if (qin <= F_ZERO) {
		imsl_e1str(1, qin);

/*		imsl_ermes(5, 2, "The second parameter of the imsl_beta distribution, QIN = %(r1), must be positive.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_QIN_MUST_BE_POSITIVE);
		goto L_9000;
	}
	/* Check P. */
	if (p <= F_ZERO || p >= F_ONE) {
		imsl_e1str(1, p);

/*		imsl_ermes(5, 3, "The probability which the inverse of the imsl_beta distribution is based on, P = %(r1), must be in the exclusive interval (0.0,
1.0).");*/
                imsl_ermes(IMSL_TERMINAL, IMSL_P_OUTSIDE_EXCLUSIVE_INT);
		goto L_9000;
	}
	eps = imsl_amach(4);
	sig = imsl_amach(4);
	alnsml = log(imsl_amach(1)) + F_TEN;
	xl = imsl_f_min(pin, qin);
	if (xl > F_ONE) {
		xr = imsl_f_max(pin, qin);
		if (F_TEN * xl > xr) {
			ic = 0;
			xl = F_ZERO;
			xr = F_ONE;
			fxl = -p;
			/* Bisection method. */
	L_10:
			x = (xl + xr) * F_HALF;
			p1 = imsl_f_beta_incomplete(x, pin, qin);
			fcs = p1 - p;
			if (fcs * fxl <= F_ZERO) {
				xr = x;
			} else {
				xl = x;
				fxl = fcs;
			}
			xrmxl = xr - xl;
			if (xrmxl <= sig && fabs(fcs) <= eps)
				goto L_40;
			ic += 1;
			if (ic <= 30)
				goto L_10;
		}
	}
	/* Use Newtons method for skewed cases. */
	if (p <= F_HALF) {
		a = pin;
		b = qin;
		q0 = log(p);
	} else {
		q0 = log(F_ONE - p);
		a = qin;
		b = pin;
	}
	xt = a / (a + b);
	dtemp = imsl_f_log_gamma(a + b) - imsl_f_log_gamma(a) - imsl_f_log_gamma(b);
	dtemp += -(a + b) * log(a + b) + (a - F_HALF) * log(a) + (b - F_HALF) *
		log(b);
	dtemp += F_HALF * log(b / a) + a * log(F_ONE + b / a) + b * log(F_ONE + a / b);
	for (nc = 1; nc <= 100; nc++) {
		temp = log(15.0e0 + a + b);
		fn = 0.7e0 * temp * temp + imsl_f_max(xt * (a + b) - a, F_ZERO);
		temp = a + fn + fn;
		afn = (Mint) fn + F_ONE;
		c = F_ONE - (a + b) * xt / temp;
		zi = F_TWO / (c + sqrt(c * c - F_FOUR * fn * (fn - b) * xt / (temp * temp)));
L_20:
		afn -= F_ONE;
		if (afn >= F_HALF) {
			temp = a + afn + afn;
			zi = (temp - F_TWO) * (temp - F_ONE - afn * (afn - b) * xt * zi /
					       temp);
			temp = a + afn - F_ONE;
			zi = F_ONE / (F_ONE - temp * (temp + b) * xt / zi);
			goto L_20;
		}
		zz = zi;
		temp = log(xt);
		if (temp <= alnsml) {
			if (p <= F_HALF) {
				x = F_ZERO;
			} else {
				x = F_ONE;
			}
			goto L_40;
		}
		qx = dtemp + a * temp + b * log(F_ONE - xt) + log(zz);
		xc = (q0 - qx) * (F_ONE - xt) * zz / a;
		xc = imsl_f_max(xc, -0.99e0);
		temp = F_HALF / xt - F_HALF;
		xc = imsl_f_min(xc, temp);
		xt *= F_ONE + xc;
		if (fabs(xc) < sig) {
			if (p <= F_HALF) {
				x = xt;
			} else {
				x = F_ONE - xt;
			}
			goto L_40;
		}
	}

	if (p <= F_HALF) {
		x = xt;
	} else {
		x = F_ONE - xt;
	}
	imsl_e1str(1, p);
	imsl_e1str(2, x);

/*	imsl_ermes(3, 1, "BETIN for the value P = %(r1) could not be found in 100 iterations of Newtons method.  The best approximation calculated is BETIN = %(r2).");
*/
        imsl_ermes(IMSL_WARNING, IMSL_BEST_BETIN_APPROXIMATION);
L_40:
	betin_v = x;
L_9000:
	E1POP("imsl_betin","imsl_dbetin");

	return (betin_v);
}				/* end of function */
