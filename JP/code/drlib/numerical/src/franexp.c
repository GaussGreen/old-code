#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif


#ifdef ANSI
static VA_LIST_HACK  l_random_exponential(Mint, va_list);
static void     l_rnexp (Mint, Mfloat[]);
#else
static void      l_rnexp ();
static VA_LIST_HACK   l_random_exponential();
#endif

static Mfloat   *lv_exponential;

#ifdef ANSI
Mfloat *imsl_f_random_exponential(Mint n_numbers, ...)
#else
Mfloat *imsl_f_random_exponential(n_numbers, va_alist)
    Mint        n_numbers;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr, n_numbers);

    E1PSH("imsl_f_random_exponential", "imsl_d_random_exponential");

    lv_exponential = NULL;
    IMSL_CALL(l_random_exponential(n_numbers, argptr));
    va_end(argptr);

    E1POP("imsl_f_random_exponential", "imsl_d_random_exponential");

    return(lv_exponential);
}

#ifdef ANSI
static VA_LIST_HACK l_random_exponential(Mint n_numbers, va_list argptr)
#else
static VA_LIST_HACK l_random_exponential(n_numbers, argptr)
    Mint        n_numbers;
    va_list     argptr;
#endif
{
    Mint        code, user_exponential = 0, arg_number = 2, ner = 1;

    code = va_arg(argptr, Mint);
    if (code == (Mint)IMSL_RETURN_USER) {
        lv_exponential = va_arg(argptr, Mfloat*);
        user_exponential = 1;
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


    if (!user_exponential) {
        lv_exponential = (Mfloat *) imsl_malloc (n_numbers*sizeof(Mfloat));
        if (!lv_exponential){
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
            goto RETURN;

        }
    }
    l_rnexp(n_numbers, lv_exponential);
    if (imsl_n1rty(0)>3 && !user_exponential) {
	imsl_free(lv_exponential);
	lv_exponential = NULL;
    }

RETURN:

    return(argptr);
}
/*Translated by FOR_C++, v0.1, on 09/10/91 at 16:56:36 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 09/10/91 at 16:56:35
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  RNEXP/DRNEXP (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    May 30, 1991

    Purpose:    Generate pseudorandom numbers from a standard
                exponential distribution.

    Usage:      CALL RNEXP (NR, R)

    Arguments:
       NR     - Number of random numbers to generate.  (Input)
       R      - Vector of length NR containing the random standard
                exponential deviates.  (Output)

    Remark:
       The routine RNSET  can be used to initialize the seed
       of the random number generator.  The routine RNOPT  can
       be used to select the form of the generator.

    Keywords:   Monte Carlo; Simulation; Poisson process; Univariate
                continuous; Random numbers

    GAMS:       L6a5

    Chapter:    STAT/LIBRARY Random Number Generation

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_rnexp(Mint nr, Mfloat r[])
#else
static void l_rnexp(nr, r)
	Mint            nr;
	Mfloat           r[];
#endif
{
	Mint             _l0, _l1, i, ner;


	if (nr <= 0) {
		imsl_e1psh("l_rnexp  ");
		ner = 1;
		imsl_c1iarg(nr, "nr", 1, 0, &ner);
		imsl_e1pop("l_rnexp");
	} else {
		/*
		 * Get uniform deviates and transform to exponential using
		 * 1.0 - CDF**-1.
		 */
                imsl_f_random_uniform(nr, IMSL_RETURN_USER, r, 0);
		for (i = 1; i <= nr; i++) {
			r[i - 1] = -log(r[i - 1]);
		}
	}
	return;
}				/* end of function */
