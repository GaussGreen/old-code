#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static void	    l_rnpoi(Mint, Mfloat, Mint[]);
static VA_LIST_HACK	    l_random_poisson(Mint, Mfloat, va_list);
#else
static void	    l_rnpoi();
static VA_LIST_HACK	    l_random_poisson();
#endif

static Mint	*lv_poisson;

#ifdef ANSI
Mint *imsl_random_poisson(Mint n_numbers, Mfloat theta, ...)
#else
Mint *imsl_random_poisson(n_numbers, theta, va_alist)
    Mint	n_numbers;
    Mfloat	theta;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr, theta);

    imsl_e1psh("imsl_random_poisson");

    lv_poisson = NULL;
    IMSL_CALL(l_random_poisson(n_numbers, theta, argptr));
    va_end(argptr);

    imsl_e1pop("imsl_random_poisson");

    return(lv_poisson);
}

#ifdef ANSI
static VA_LIST_HACK l_random_poisson(Mint n_numbers, Mfloat theta, va_list argptr)
#else
static VA_LIST_HACK l_random_poisson(n_numbers, theta, argptr)
    Mint	n_numbers;
    Mfloat	theta;
    va_list	argptr;
#endif
{
    Mint	code, user_poisson = 0, ner, arg_number = 3;

    code = va_arg(argptr, Mint);
    if (code == IMSL_RETURN_USER) {
	lv_poisson = va_arg(argptr, Mint*);
	user_poisson = 1;
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

    if (!user_poisson) {
	lv_poisson = (Mint *) imsl_malloc (n_numbers*sizeof(Mint));
	if (!lv_poisson){
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
	    goto RETURN;

	}
    }

    l_rnpoi(n_numbers, theta, lv_poisson);
    if (imsl_n1rty(0)>3 && !user_poisson){
      imsl_free(lv_poisson);
      lv_poisson = NULL;
    }

RETURN:

    return(argptr);
}
/*
  -----------------------------------------------------------------------
    IMSL Name:  RNPOI (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    June 21, 1990

    Purpose:    Generate pseudorandom numbers from a Poisson
                distribution.

    Usage:      CALL RNPOI (NR, THETA, IR)

    Arguments:
       NR     - Number of random numbers to generate.  (Input)
       THETA  - Mean of the Poisson distribution.  (Input)
                THETA must be positive.
       IR     - Vector of length NR containing the random Poisson
                deviates.  (Output)

    Remark:
       The IMSL routine RNSET can be used to initialize the seed of the
       random number generator.  The routine RNOPT can be used to select
       the form of the generator.

    Keywords:   Monte Carlo; Simulation; Univariate discrete; Random
                numbers

    GAMS:       L6a16

    Chapter:    STAT/LIBRARY Random Number Generation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_rnpoi(Mint nr, Mfloat theta, Mint ir[])
#else
static void l_rnpoi(nr, theta, ir)
	Mint            nr;
	Mfloat          theta;
	Mint            ir[];
#endif
{
	Mint            i, ix, j, ner;
	static Mint     m;
	Mfloat          alv, d, f, fp, q, u, ub, v, x;
	static Mfloat   al, c, fm, p1, p2, p3, p4, xl, xll, xlr, xm, xr;
	static Mfloat   ymu = -1.0;



	if (theta <= F_ZERO || nr <= 0) {
		imsl_e1psh("l_rnpoi");
		ner = 1;
		imsl_c1iarg(nr, "nr", 1, 0, &ner);
		if (theta <= F_ZERO) {
			imsl_e1str(1, theta);

			imsl_ermes(IMSL_TERMINAL, IMSL_NEED_POSITIVE_THETA);
		}
		imsl_e1pop("l_rnpoi");
		goto L_9000;
	}
	if (theta >= 15.0) {
		/*
		 * Use Schmeiser-Kachitvichyanukul four region
		 * acceptance/rejection method.
		 */
		if (theta != ymu) {
			/* Setup */
			ymu = theta;
			m = ymu;
			fm = m;
			p1 = (int) (2.195 * sqrt(fm) - 2.2) + F_HALF;
			c = 0.133 + 8.56 / (6.83 + ymu);
			xm = m + F_HALF;
			xl = xm - p1;
			xr = xm + p1;
			al = (ymu - xl) / ymu;
			xll = al * (F_ONE + F_HALF * al);
			al = (xr - ymu) / xr;
			xlr = al * (F_ONE + F_HALF * al);
			p2 = p1 * (F_ONE + c + c);
			p3 = p2 + (0.109 + 8.25 / (10.86 + ymu)) / xll;
			p4 = p3 + c / xlr;
		}
		for (i = 1; i <= nr; i++) {
			/* Generate uniforms to start. */
	L_10:
			imsl_f_random_uniform(1, IMSL_RETURN_USER, &u, 0);
			u *= p4;
			imsl_f_random_uniform(1, IMSL_RETURN_USER, &v, 0);
			if (u <= p1) {
				/* Triangular region. */
				ix = xm - p1 * v + u;
			} else {
				/*
				 * Use acceptance/rejection over other
				 * regions.
				 */
				if (u <= p2) {
					/* Parallelogram region. */
					x = xl + (u - p1) / c;
					v = v * c + F_ONE - fabs(fm - x + F_HALF) / p1;
					/* If V .gt. 1, reject. */
					if (v > F_ONE)
						goto L_10;
					ix = x;
				} else if (u <= p3) {
					/* Left tail. */
					ix = xl + log(v) / xll;
					/* If IX .Lt. 0, reject. */
					if (ix < 0)
						goto L_10;
					v *= (u - p2) * xll;
				} else {
					/* Right tail. */
					ix = xr - log(v) / xlr;
					v *= (u - p3) * xlr;
				}
				/*
				 * Now do acceptance/rejection by comparing V
				 * to the scaled Poisson mass function.
				 */
				if (m < 100 || ix <= 50) {
					/* Evaluate explicitly. */
					f = F_ONE;
					if (m < ix) {
						for (j = m + 1; j <= ix; j++) {
							f = f * ymu / j;
						}
					} else if (m > ix) {
						for (j = ix + 1; j <= m; j++) {
							f = f * j / ymu;
						}
					}
					/* If V .gt. F, reject. */
					if (v > f)
						goto L_10;
				} else {
					/* Squeeze using ALOG(F(X)) */
					x = ix;
					q = (ymu - x) / x;
					ub = x - ymu + (x + F_HALF) * q * (F_ONE + q * (-F_HALF + q /
							    F_THREE)) + 0.00084;
					alv = log(v);
					/* If ALV .gt. UB, reject. */
					if (alv > ub)
						goto L_10;
					d = (x + F_HALF) * 0.25 * powl(q * q, 2);
					if (q < F_ZERO)
						d /= F_ONE + q;
					if (alv >= ub - d - 0.004) {
						/*
						 * Use Stirling's formula for
						 * the final
						 * acceptance/rejection test.
						 */
						if (alv > (fm + F_HALF) * log(fm / ymu) + (x +
											F_HALF) * log(ymu / x) - fm + x + (F_ONE / fm - F_ONE /
															x) / 12.0 + F_ONE / (360.0 * x * x * x) - F_ONE / (360.0 *
							      fm * fm * fm))
							goto L_10;
					}
				}
			}
			ir[i - 1] = ix;
		}
	} else {
		/* Use inverse CDF. */
		if (theta != ymu) {
			/* Setup */
			ymu = theta;
			fm = exp(-ymu);
		}
		for (i = 1; i <= nr; i++) {
			x = F_ZERO;
			fp = fm;
			imsl_f_random_uniform(1, IMSL_RETURN_USER, &u, 0);
	L_50:
			if (u > fp) {
				x += F_ONE;
				u -= fp;
				fp = fp * ymu / x;
				goto L_50;
			}
			ir[i - 1] = x;
		}
	}
L_9000:
	return;
}				/* end of function */
