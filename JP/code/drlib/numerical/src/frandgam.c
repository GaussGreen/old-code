#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static Mfloat	*lv_gamma;
#ifdef ANSI
static void l_rngam(Mint, Mfloat, Mfloat[]);
static VA_LIST_HACK  l_random_gamma(Mint, Mfloat, va_list);
#else
static void l_rngam();
static VA_LIST_HACK  l_random_gamma();
#endif

#ifdef ANSI
Mfloat *imsl_f_random_gamma(Mint n_numbers, Mfloat a, ...)
#else
Mfloat *imsl_f_random_gamma(n_numbers, a, va_alist)
    Mint        n_numbers;
    Mfloat	a;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr, a);

    E1PSH("imsl_f_random_gamma", "imsl_d_random_gamma");

    lv_gamma = NULL;
    IMSL_CALL(l_random_gamma(n_numbers, a, argptr));
    va_end(argptr);

    E1POP("imsl_f_random_gamma", "imsl_d_random_gamma");

    return(lv_gamma);
}

#ifdef ANSI
static VA_LIST_HACK l_random_gamma(Mint n_numbers, Mfloat a, va_list argptr)
#else
static VA_LIST_HACK l_random_gamma(n_numbers, a, argptr)
    Mint        n_numbers;
    Mfloat	a;
    va_list     argptr;
#endif
{
    Mint        code, user_gamma = 0, ner, arg_number = 3;

    code = va_arg(argptr, Mint);
    if (code == (Mint)IMSL_RETURN_USER) {
        lv_gamma = va_arg(argptr, Mfloat*);
        user_gamma = 1;
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


    if (!user_gamma) {
        lv_gamma = (Mfloat *) imsl_malloc (n_numbers*sizeof(Mfloat));
        if (!lv_gamma){
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
            goto RETURN;
        }
    }
    l_rngam(n_numbers, a, lv_gamma);
    if (imsl_n1rty(0)>3 && !user_gamma) {
        imsl_free(lv_gamma);
        lv_gamma = NULL;
    }

RETURN:

    return(argptr);
}


/*
  -----------------------------------------------------------------------
    IMSL Name:  RNGAM/DRNGAM (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 20, 1990

    Purpose:    Generate pseudorandom numbers from a standard imsl_gamma
                distribution.

    Usage:      CALL RNGAM (NR, A, R)

    Arguments:
       NR     - Number of random numbers to generate.  (Input)
       A      - The shape parameter of the imsl_gamma distribution.  (Input)
                This parameter must be positive.
       R      - Vector of length NR containing the random standard imsl_gamma
                deviates.  (Output)

    Remark:
       The IMSL routine RNSET can be used to initialize the seed of the
       random number generator.  The routine RNOPT can be used to select
       the form of the generator.

    Keywords:   Monte Carlo; Simulation; Chi-squared distribution;
                Erlang distribution; Compartmental analysis; Univariate
                continuous; Random numbers

    GAMS:       L6a7

    Chapter:    STAT/LIBRARY Random Number Generation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_rngam(Mint nr, Mfloat a, Mfloat r[])
#else
static void l_rngam(nr, a, r)
	Mint           nr;
	Mfloat        a, r[];
#endif
{
	Mint            i, ner;
	Mfloat          aa, aleta, ex, p, temp, u, v, w1, w2, x, xarg, xmax, xx;
	static Mfloat   e = 2.718281828459045235360287;
	Mfloat   a1 = F_ZERO;
	Mfloat   xll = F_ZERO;
	Mfloat   xlr = F_ZERO;
	Mfloat   x1 = F_ZERO;
	Mfloat   x2 = F_ZERO;
	Mfloat   x3 = F_ZERO;
	Mfloat   x4 = F_ZERO;
	Mfloat   x5 = F_ZERO;
	Mfloat   f1 = F_ZERO;
	Mfloat   f2 = F_ZERO;
	Mfloat   f4 = F_ZERO;
	Mfloat   f5 = F_ZERO;
	Mfloat   p1 = F_ZERO;
	Mfloat   p10 = F_ZERO;
	Mfloat   p2 = F_ZERO;
	Mfloat   p3 = F_ZERO;
	Mfloat   p4 = F_ZERO;
	Mfloat   p5 = F_ZERO;
	Mfloat   p6 = F_ZERO;
	Mfloat   p7 = F_ZERO;
	Mfloat   p8 = F_ZERO;
	Mfloat   p9 = F_ZERO;



	if (a <= F_ZERO || nr <= 0) {
		E1PSH("imsl_f_random_gamma","imsl_d_random_gamma");
		ner = 1;
		imsl_c1iarg(nr, "nr", 1, 0, &ner);
		if (a <= F_ZERO) {
			imsl_e1str(1, a);

			imsl_ermes(IMSL_TERMINAL, IMSL_SHAPE_PARAMETER_A);
		}
		E1POP("imsl_f_random_gamma","imsl_d_random_gamma");
		goto L_9000;
	}
	if (a <= F_ONE) {
		/* First, try ad hoc procedures */
		if (a == F_HALF) {
			/*
			 * For A = 0.5, use normal**2 deviate divided by 2
			 */
			imsl_f_random_normal(nr, IMSL_RETURN_USER, r, 0);
			for (i = 1; i <= nr; i++) {
				r[i - 1] *= r[i - 1] * F_HALF;
			}
		} else if (a == F_ONE) {
			/* For A = 1.0, generate exponentials */
			imsl_f_random_uniform(nr, IMSL_RETURN_USER, r, 0);
			for (i = 0; i < nr; i++) {
			    r[i] = -log(r[i]);
		}
		} else {
			/* A less than 1.0 */
			w1 = (e + a) / e;
			w2 = F_ONE / a;
			aa = imsl_amach(1);
			aleta = log(aa);
			i = 1;
			/*
			 * Generate a uniform (0,1) deviate and an
			 * exponential deviate.
			 */
	L_20:
			imsl_f_random_uniform(1, IMSL_RETURN_USER, &xx, 0);
			ex = -log(xx);
			/* Rejection test */
			imsl_f_random_uniform(1, IMSL_RETURN_USER, &p, 0);
			p = w1 * p;

			if (p <= F_ONE) {
				/* Small X, check for underflow */
				if (w2 * log(p) < aleta) {
					temp = imsl_amach(1);
				} else {
					temp = pow(p, w2);
				}
				x = temp;
			} else {
				/* large X */
				x = -log((w1 - p) * w2);
				temp = -(a - F_ONE) * log(x);
			}
			if (temp > ex)
				goto L_20;
			r[i - 1] = x;
			i += 1;
			if (i <= nr)
				goto L_20;
		}
		goto L_9000;
	}
	/* A greater than 1.0 */
	i = 1;
	if (a != a1) {
		/* Initialization */
		xmax = log(imsl_amach(2)) - 0.001;
		x3 = a - F_ONE;
		temp = sqrt(x3);
		x2 = F_ZERO;
		x1 = F_ZERO;
		f1 = F_ZERO;
		f2 = F_ZERO;
		xll = F_ONE;
		if (temp * (F_ONE + imsl_amach(4)) < x3) {
			x2 = x3 - temp;
			x1 = x2 * (F_ONE - F_ONE / temp);
			xll = F_ONE - x3 / x1;
			xarg = x3 * log(x1 / x3) + x3 - x1;
			if (xarg >= xmax) {
				imsl_e1psh("RNGAM ");
                                imsl_ermes(IMSL_TERMINAL,
				IMSL_A_IS_TOO_LARGE);
				imsl_e1pop("RNGAM ");
				goto L_9000;
			}
			f1 = exp(xarg);
			xarg = x3 * log(x2 / x3) + x3 - x2;
			if (xarg >= xmax) {
				imsl_e1psh("RNGAM ");
                                imsl_ermes(IMSL_TERMINAL,
				IMSL_A_IS_TOO_LARGE);
				imsl_e1pop("RNGAM ");
				goto L_9000;
			}
			f2 = exp(xarg);
		}
		x4 = x3 + temp;
		x5 = x4 * (F_ONE + F_ONE / temp);
		xlr = F_ONE - x3 / x5;
		xarg = x3 * log(x4 / x3) + x3 - x4;
		if (xarg >= xmax) {
			imsl_e1psh("RNGAM ");

                        imsl_ermes(IMSL_TERMINAL, IMSL_A_IS_TOO_LARGE);
			imsl_e1pop("RNGAM ");
			goto L_9000;
		}
		f4 = exp(xarg);
		xarg = x3 * log(x5 / x3) + x3 - x5;
		if (xarg >= xmax) {
			imsl_e1psh("RNGAM ");
                        imsl_ermes(IMSL_TERMINAL, IMSL_A_IS_TOO_LARGE);
			imsl_e1pop("RNGAM ");
			goto L_9000;
		}
		f5 = exp(xarg);
		/*
		 * Calculate probability for each of the ten regions
		 */
		p1 = f2 * (x3 - x2);
		p2 = f4 * (x4 - x3) + p1;
		p3 = f1 * (x2 - x1) + p2;
		p4 = f5 * (x5 - x4) + p3;
		p5 = (F_ONE - f2) * (x3 - x2) + p4;
		p6 = (F_ONE - f4) * (x4 - x3) + p5;
		p7 = (f2 - f1) * (x2 - x1) * F_HALF + p6;
		p8 = (f4 - f5) * (x5 - x4) * F_HALF + p7;
		p9 = -f1 / xll + p8;
		p10 = f5 / xlr + p9;
		/*
		 * Set the save variable to the value of the shape parameter
		 */
		a1 = a;
	}
	/*
	 * Begin the generation process, first, generate one uniform(0,1)
	 * deviate
	 */
L_30:
	imsl_f_random_uniform(1, IMSL_RETURN_USER, &u, 0);
	u = u * p10;
	/*
	 * The four regions with zero probability of rejection
	 */
	if (u <= p4) {
		if (u <= p1) {
			x = x2 + u / f2;
			goto L_60;
		}
		if (u <= p2) {
			x = x3 + (u - p1) / f4;
			goto L_60;
		} else if (u <= p3) {
			x = x1 + (u - p2) / f1;
			goto L_60;
		}
		x = x4 + (u - p3) / f5;
		goto L_60;
	}
	/*
	 * The two regions using rectangular rejection
	 */
	imsl_f_random_uniform(1, IMSL_RETURN_USER, &w1, 0);
	if (u <= p5) {
		x = x2 + (x3 - x2) * w1;
		if ((u - p4) / (p5 - p4) <= w1)
			goto L_60;
		v = f2 + (u - p4) / (x3 - x2);
		goto L_50;
	} else if (u <= p6) {
		x = x3 + (x4 - x3) * w1;
		if ((p6 - u) / (p6 - p5) >= w1)
			goto L_60;
		v = f4 + (u - p5) / (x4 - x3);
		goto L_50;
		/* The two triangular regions */
	} else if (u <= p8) {
		imsl_f_random_uniform(1, IMSL_RETURN_USER, &w2, 0);
		if (w2 > w1)
			w1 = w2;
		if (u > p7)
			goto L_40;
		x = x1 + (x2 - x1) * w1;
		v = f1 + F_TWO * w1 * (u - p6) / (x2 - x1);
		if (v <= f2 * w1)
			goto L_60;
		goto L_50;
L_40:
		x = x5 - w1 * (x5 - x4);
		v = f5 + F_TWO * w1 * (u - p7) / (x5 - x4);
		if (v <= f4 * w1)
			goto L_60;
		goto L_50;
		/* The two exponential regions */
	} else if (u < p9) {
		u = (p9 - u) / (p9 - p8);
		x = x1 - log(u) / xll;
		if (x <= F_ZERO)
			goto L_30;
		if (w1 < (xll * (x1 - x) + F_ONE) / u)
			goto L_60;
		v = w1 * f1 * u;
		goto L_50;
	}
	if (p10 == u)
		goto L_30;
	u = (p10 - u) / (p10 - p9);
	x = x5 - log(u) / xlr;
	if (w1 < (xlr * (x5 - x) + F_ONE) / u)
		goto L_60;
	v = w1 * f5 * u;
	/* preform the standard rejection */
L_50:
	if (log(v) > x3 * log(x / x3) + x3 - x)
		goto L_30;
	/* Accept the variate and assign to R */
L_60:
	r[i - 1] = x;
	i += 1;
	if (i <= nr)
		goto L_30;
L_9000:
	return;
}				/* end of function */
