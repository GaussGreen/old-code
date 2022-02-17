#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static Mfloat   *lv_beta;
#ifdef ANSI
static void l_rnbet(Mint, Mfloat, Mfloat, Mfloat[]);
static VA_LIST_HACK  l_random_beta(Mint, Mfloat, Mfloat, va_list);
#else
static void l_rnbet();
static VA_LIST_HACK  l_random_beta();
#endif
 
#ifdef ANSI
Mfloat *imsl_f_random_beta(Mint n_numbers, Mfloat a, Mfloat b, ...)
#else
Mfloat *imsl_f_random_beta(n_numbers, a, b, va_alist)
    Mint        n_numbers;
    Mfloat      a;
    Mfloat      b;
    va_dcl
#endif
{
    va_list     argptr;
    VA_START(argptr, b);
 
    E1PSH("imsl_f_random_beta", "imsl_d_random_beta");

    lv_beta = NULL;
    IMSL_CALL(l_random_beta(n_numbers, a, b, argptr));
    va_end(argptr);

    E1POP("imsl_f_random_beta", "imsl_d_random_beta");

    return(lv_beta);
}

#ifdef ANSI
static VA_LIST_HACK l_random_beta(Mint n_numbers, Mfloat a, Mfloat b, va_list argptr)
#else
static VA_LIST_HACK l_random_beta(n_numbers, a, b, argptr)
    Mint        n_numbers;
    Mfloat      a;
    Mfloat      b;
    va_list     argptr;
#endif
{
    Mint        code, user_beta = 0,ner, arg_number = 4;

    code = va_arg(argptr, Mint);
    if (code == (Mint)IMSL_RETURN_USER) {
        lv_beta = va_arg(argptr, Mfloat*);
        user_beta = 1;
    }
    else if (code != 0){
        imsl_e1sti (1, code);
        imsl_e1sti (2, arg_number);
        imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
        goto RETURN;
    }

    ner = 1;
    if (n_numbers <1) {
         imsl_c1iarg(n_numbers, "n_random", 1, 0, &ner);
         goto RETURN;
     }

    if (!user_beta) {
        lv_beta = (Mfloat *) imsl_malloc (n_numbers*sizeof(Mfloat));
        if (!lv_beta){
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
            goto RETURN;
        }
    }
    l_rnbet(n_numbers, a, b, lv_beta);
    if (imsl_n1rty(0)>3 && !user_beta) {
        imsl_free(lv_beta);
        lv_beta = NULL;
    }

RETURN:

    return(argptr);
}
/*
  -----------------------------------------------------------------------
    IMSL Name:  RNBET/DRNBET (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 10, 1985

    Purpose:    Generate pseudorandom numbers from a imsl_beta distribution.

    Usage:      CALL RNBET (NR, PIN, QIN, R)

    Arguments:
       NR     - Number of random numbers to generate.  (Input)
       PIN    - First imsl_beta distribution parameter.  (Input)
                PIN must be positive.
       QIN    - Second imsl_beta distribution parameter.  (Input)
                QIN must be positive.
       R      - Vector of length NR containing the random standard imsl_beta
                deviates.  (Output)

    Remark:
       The IMSL routine RNSET can be used to initialize the seed of the
       random number generator.  The routine RNOPT can be used to select
       the form of the generator.

    Keywords:   Monte Carlo; Simulation; Arcsine distribution;
                Power function distribution; Univariate continuous;
                Random numbers

    GAMS:       L6a2

    Chapter:    STAT/LIBRARY Random Number Generation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_rnbet(Mint nr, Mfloat pin, Mfloat qin, Mfloat r[])
#else
static void l_rnbet(nr, pin, qin, r)
	Mint            nr;
	Mfloat          pin, qin, r[];
#endif
{
	short int       check, pmin;
	Mint             i, imsl_i1x, icount, ii, ner;
	Mfloat           a, a1, a5, ain, alpha, alv, b, imsl_beta, bin,
	                c, cons3, cons4, d, f1, f2, f4, f5, gama, omg,
	                omt, p1, p10, p2, p3, p4, p5, p6, p7, p8, p9, pma,
	                pp, q1, qmb, qq, rr, s, s1, s2, t, temp,
	                temp1, temp2, tempa, tempb, trr, u, ua[2], v, w,
	                w2, x, x1, x2, x3, x4, x5, y, z;
	static Mfloat    eps = 0.9537e-06;
	static Mfloat    cons1 = 1.386294;
	static Mfloat    cons2 = 2.609438;



	cons3 = F_ONE - imsl_amach(4);
	cons4 = imsl_amach(3);
	if ((nr <= 0 || pin <= F_ZERO) || qin <= F_ZERO) {
		E1PSH("imsl_f_random_beta","imsl_d_random_beta");
		ner = 1;
		imsl_c1iarg(nr, "nr", 1, 0, &ner);
		if (pin <= F_ZERO) {
			imsl_e1str(1, pin);

/*			imsl_ermes(5, 2, "PIN = %(r1).  But, it must be positive.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_NEED_POSITIVE_PIN_VALUE);
		}
		if (qin <= F_ZERO) {
			imsl_e1str(1, qin);

/*			imsl_ermes(5, 3, "QIN = %(r1).  But, it must be positive.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_NEED_POSITIVE_QIN_VALUE);
		}
		E1POP("imsl_f_random_beta","imsl_d_random_beta");
		goto L_9000;
	}
	if (pin < F_ONE && qin < F_ONE) {
		/*
		 * JOHNKS method - both PIN and QIN less than 1.0
		 */
		temp = log(imsl_amach(2));
		p1 = F_ONE / pin;
		q1 = F_ONE / qin;
		/*
		 * If PIN and QIN are both less than .5, check for overflow
		 * when dividing
		 */
		if (p1 > F_TWO && q1 > F_TWO) {
			check = 1;
		} else {
			check = 0;
		}
		for (i = 1; i <= nr; i++) {
			/* Obtain uniform deviates */
			icount = 1;
	L_10:
			imsl_f_random_uniform(2, IMSL_RETURN_USER, ua, 0);
			x = pow(ua[0], p1);
			y = pow(ua[1], q1);
			y += x;
			/* Reject or accept */
			if (y > F_ONE)
				goto L_10;
			if (check) {
				if (x < cons4 || y < cons4) {
					icount += 1;
				} else if ((log(x) - log(y)) > temp) {
					icount += 1;
				} else {
					goto L_30;
				}
				if (icount > 1000) {
					imsl_e1psh("RNBET ");
					imsl_e1sti(1, i);

/*					imsl_ermes(5, 4, "1000 unsuccessful attempts were made to assign a imsl_beta random deviate to R(%(i1)), since PIN and QIN are so small.  All elements of R beyond %(i1) are set to not-a-number.");
*/
                                        imsl_ermes(IMSL_TERMINAL,
					IMSL_ELMNTS_SET_TO_NAN);
					for (ii = i; ii <= nr; ii++) {
						r[ii - 1] = imsl_amach(6);
					}
					imsl_e1pop("RNBET ");
					goto L_9000;
				}
				goto L_10;
			}
	L_30:
			r[i - 1] = x / y;
			if (r[i - 1] >= F_ONE)
				r[i - 1] = cons3;
			if (r[i - 1] <= cons4)
				r[i - 1] = cons4;
		}
	} else if (pin == F_ONE || qin == F_ONE) {
		/*
		 * Use the inverse of the cumulative distribution function
		 * for PIN or QIN equal to 1.0
		 */
		imsl_f_random_uniform(nr, IMSL_RETURN_USER, r, 0);
		if (qin == F_ONE) {
			/* QIN is 1.0 */
			p1 = F_ONE / pin;
			for (i = 1; i <= nr; i++) {
				r[i - 1] = pow(r[i - 1], p1);
				if (r[i - 1] >= F_ONE)
					r[i - 1] = cons3;
				if (r[i - 1] <= cons4)
					r[i - 1] = cons4;
			}
		} else {
			/* PIN is 1.0 */
			q1 = F_ONE / qin;
			for (i = 1; i <= nr; i++) {
				r[i - 1] = F_ONE - pow(r[i - 1], q1);
				if (r[i - 1] >= F_ONE)
					r[i - 1] = cons3;
				if (r[i - 1] <= cons4)
					r[i - 1] = cons4;
			}
		}
	} else if ((pin < F_ONE && qin > F_ONE) || (qin < F_ONE && pin > F_ONE)) {
		/*
		 * ATKINSONS method - one of PIN AND QIN less than 1.0 and
		 * the other greater than 1.0
		 */
		if (pin < F_ONE) {
			a = pin;
			b = qin;
			imsl_i1x = 0;
		} else {
			a = qin;
			b = pin;
			imsl_i1x = 1;
		}
		tempa = a - F_ONE;
		tempb = b - F_ONE;
		t = -tempa / (b + F_ONE - a);
		s1 = F_ONE;
		s2 = pow(t, tempa);
		temp = b * t;
		gama = temp / (temp + a * (pow(F_ONE - t, b)));
		pma = pin - a;
		qmb = qin - b;
		omt = F_ONE - t;
		omg = F_ONE - gama;
		ain = F_ONE / a;
		bin = F_ONE / b;
		for (i = 1; i <= nr; i++) {
	L_70:
			imsl_f_random_uniform(2, IMSL_RETURN_USER, ua, 0);
			if (ua[0] <= gama) {
				x = t * pow(ua[0] / gama, ain);
				f1 = pow(F_ONE - x, tempb);
				if (imsl_i1x > 1) {
					f1 *= pow(x, pma);
				}
				if (s1 * ua[1] > f1)
					goto L_70;
			} else {
				x = F_ONE - omt * pow((F_ONE - ua[0]) / omg, bin);
				f2 = pow(x, tempa);
				if (imsl_i1x > 1) {
					f2 *= pow(F_ONE - x, qmb);
				}
				if (s2 * ua[1] > f2)
					goto L_70;
			}
			r[i - 1] = x;
			if (imsl_i1x == 1)
				r[i - 1] = F_ONE - x;
			if (r[i - 1] >= F_ONE)
				r[i - 1] = cons3;
			if (r[i - 1] <= cons4)
				r[i - 1] = cons4;
		}
	} else {
		if (nr > 4) {
			/*
			 * SCHMEISERS method - for both PIN and QIN greater
			 * than 1.0 and NR greater than or equal to 4
			 */
			pp = pin - F_ONE;
			qq = qin - F_ONE;
			rr = pp + qq;
			trr = -rr - rr;
			c = rr * log(rr) - pp * log(pp) - qq * log(qq);
			x1 = F_ZERO;
			x2 = F_ZERO;
			f1 = F_ZERO;
			f2 = F_ZERO;
			f4 = F_ZERO;
			f5 = F_ZERO;
			x3 = pp / rr;
			x4 = F_ONE;
			x5 = F_ONE;

			if (rr > (F_ONE + eps)) {
				d = sqrt(pp * qq / (rr - F_ONE)) / rr;
				if (d * (F_ONE + eps) < x3) {
					x2 = x3 - d;
					x1 = x2 - (x2 * (F_ONE - x2)) / (pp - rr * x2);
					a1 = pp / x1 - qq / (F_ONE - x1);
					f1 = exp(c + pp * log(x1) + qq * log(F_ONE - x1));
					f2 = exp(c + pp * log(x2) + qq * log(F_ONE - x2));
				}
				if (d < (F_ONE - x3 - eps)) {
					x4 = x3 + d;
					x5 = x4 - ((x4 * (F_ONE - x4)) / (pp - rr * x4));
					a5 = qq / (F_ONE - x5) - pp / x5;
					f4 = exp(c + pp * log(x4) + qq * log(F_ONE - x4));
					f5 = exp(c + pp * log(x5) + qq * log(F_ONE - x5));
				}
			}
			/* Calculate areas of the ten regions */
			p1 = F_ZERO;
			if (f2 > F_ZERO)
				p1 = f2 * (x3 - x2);
			p2 = p1;
			if (f4 > F_ZERO)
				p2 = f4 * (x4 - x3) + p1;
			p3 = p2;
			if (f1 > F_ZERO)
				p3 = f1 * (x2 - x1) + p2;
			p4 = p3;
			if (f5 > F_ZERO)
				p4 = f5 * (x5 - x4) + p3;
			p5 = (F_ONE - f2) * (x3 - x2) + p4;
			p6 = (F_ONE - f4) * (x4 - x3) + p5;
			p7 = p6;
			if (f2 > f1)
				p7 = (f2 - f1) * (x2 - x1) * F_HALF + p6;
			p8 = p7;
			if (f4 > f5)
				p8 = (f4 - f5) * (x5 - x4) * F_HALF + p7;
			p9 = p8;
			if (f1 > F_ZERO)
				p9 = f1 / a1 + p8;
			p10 = p9;
			if (f5 > F_ZERO)
				p10 = f5 / a5 + p9;
			/* Rejection procedure begins here */
			for (i = 1; i <= nr; i++) {
		L_90:
				imsl_f_random_uniform(1, IMSL_RETURN_USER, &u, 0);
				u *= p10;
				/*
				 * The four regions with 0.0 probability of
				 * rejection
				 */
				if (u <= p4) {
					if (u <= p1) {
						x = x2 + u / f2;
					} else if (u <= p2) {
						x = x3 + (u - p1) / f4;
					} else if (u <= p3) {
						x = x1 + (u - p2) / f1;
					} else {
						x = x4 + (u - p3) / f5;
					}
				} else {
					/*
					 * The two regions using rectangular
					 * rejection
					 */
					imsl_f_random_uniform(1, IMSL_RETURN_USER, &w, 0);
					if (u <= p5) {
						x = x2 + (x3 - x2) * w;
						if ((u - p4) / (p5 - p4) <= w)
							goto L_100;
						v = f2 + (u - p4) / (x3 - x2);
					} else if (u <= p6) {
						x = x3 + (x4 - x3) * w;
						if ((p6 - u) / (p6 - p5) >= w)
							goto L_100;
						v = f4 + (u - p5) / (x4 - x3);
						/* The two triangular regions */
					} else if (u <= p8) {
						imsl_f_random_uniform(1, IMSL_RETURN_USER, &w2, 0);
						if (w2 > w)
							w = w2;
						if (u <= p7) {
							x = x1 + (x2 - x1) * w;
							v = f1 + F_TWO * w * (u - p6) / (x2 - x1);
							if (v < f2 * w)
								goto L_100;
						} else {
							x = x5 - w * (x5 - x4);
							v = f5 + F_TWO * w * (u - p7) / (x5 - x4);
							if (v <= f4 * w)
								goto L_100;
						}
						/*
						 * The two exponential
						 * regions
						 */
					} else if (u <= p9) {
						u = (p9 - u) / (p9 - p8);
						x = x1 + log(u) / a1;
						if (x <= F_ZERO)
							goto L_90;
						if (w <= (a1 * (x - x1) + F_ONE) / u)
							goto L_100;
						v = w * f1 * u;
					} else {
						u = (p10 - u) / (p10 - p9);
						x = x5 - log(u) / a5;
						if (x >= F_ONE)
							goto L_90;
						if (w <= (a5 * (x5 - x) + F_ONE) / u)
							goto L_100;
						v = w * f5 * u;
					}
					/*
					 * Check easy rejection via
					 * comparison with normal density
					 * function
					 */
					alv = log(v);
					if (alv > (x - x3) * (x - x3) * trr)
						goto L_90;
					/* Perform the standard rejection */
					if (alv > (c + pp * log(x) + qq * log(F_ONE - x)))
						goto L_90;
				}
		L_100:
				r[i - 1] = x;
				if (r[i - 1] >= F_ONE)
					r[i - 1] = cons3;
				if (r[i - 1] <= cons4)
					r[i - 1] = cons4;
			}
		} else {
			/*
			 * BB method - both PIN and QIN greater than 1.0 and
			 * NR less than 4
			 */
			pmin = 0;
			a = qin;
			b = pin;
			if (qin >= pin) {
				pmin = 1;
				a = pin;
				b = qin;
			}
			alpha = a + b;
			imsl_beta = sqrt((alpha - F_TWO) / (F_TWO * a * b - alpha));
			gama = a + F_ONE / imsl_beta;
			for (i = 1; i <= nr; i++) {
		L_120:
				imsl_f_random_uniform(2, IMSL_RETURN_USER, ua, 0);
				v = imsl_beta * log(ua[0] / (F_ONE - ua[0]));
				w = a * exp(v);
				z = ua[0] * ua[0] * ua[1];
				x = gama * v - cons1;
				s = a + x - w;
				temp2 = b + w;
				if ((s + cons2) < F_FIVE * z) {
					t = log(z);
					if (s < t) {
						temp1 = x + alpha * log(alpha / temp2);
						if (temp1 < t)
							goto L_120;
					}
				}
				if (pmin) {
					r[i - 1] = w / temp2;
					if (r[i - 1] >= F_ONE)
						r[i - 1] = cons3;
					if (r[i - 1] <= cons4)
						r[i - 1] = cons4;
				} else {
					r[i - 1] = b / temp2;
					if (r[i - 1] >= F_ONE)
						r[i - 1] = cons3;
					if (r[i - 1] <= cons4)
						r[i - 1] = cons4;
				}
			}
		}
	}
L_9000:
	return;
}				/* end of function */
