#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif


#ifdef ANSI
static VA_LIST_HACK  l_random_normal(Mint, va_list);
static void     l_rnnor (Mint, Mfloat[]);
#else
static void      l_rnnor ();
static VA_LIST_HACK   l_random_normal();
#endif

static Mfloat   *lv_normal;

#ifdef ANSI
Mfloat *imsl_f_random_normal(Mint n_numbers, ...)
#else
Mfloat *imsl_f_random_normal(n_numbers, va_alist)
    Mint        n_numbers;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr, n_numbers);

    E1PSH("imsl_f_random_normal", "imsl_d_random_normal");

    lv_normal = NULL;
    IMSL_CALL(l_random_normal(n_numbers, argptr));
    va_end(argptr);

    E1POP("imsl_f_random_normal", "imsl_d_random_normal");

    return(lv_normal);
}

#ifdef ANSI
static VA_LIST_HACK l_random_normal(Mint n_numbers, va_list argptr)
#else
static VA_LIST_HACK l_random_normal(n_numbers, argptr)
    Mint        n_numbers;
    va_list     argptr;
#endif
{
    Mint        code, user_normal = 0, arg_number = 2, ner = 1;

    code = va_arg(argptr, Mint);
    if (code == (Mint)IMSL_RETURN_USER) {
        lv_normal = va_arg(argptr, Mfloat*);
        user_normal = 1;
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


    if (!user_normal) {
        lv_normal = (Mfloat *) imsl_malloc (n_numbers*sizeof(Mfloat));
        if (!lv_normal){
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
            goto RETURN;

        }
    }
    l_rnnor(n_numbers, lv_normal);
    if (imsl_n1rty(0)>3 && !user_normal) {
	imsl_free(lv_normal);
	lv_normal = NULL;
    }

RETURN:

    return(argptr);
}

/* 
  -----------------------------------------------------------------------
    IMSL Name:  RNNOR/DRNNOR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 21, 1990

    Purpose:    Generate pseudorandom numbers from a standard normal
                distribution using an inverse CDF method.

    Usage:      CALL RNNOR (NR, R)

    Arguments:
       NR     - Number of random numbers to generate.  (Input)
       R      - Vector of length NR containing the random standard normal
                deviates.  (Output)

    Remark:
       The IMSL routine RNSET can be used to initialize the seed of the
       random number generator.  The routine RNOPT can be used to select
       the form of the generator.

    Keywords:   Monte Carlo; Simulation; Gaussian distribution;
                Univariate continuous; Random numbers

    GAMS:       L6a14

    Chapter:    STAT/LIBRARY Random Number Generation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_rnnor(Mint nr, Mfloat r[])
#else
static void l_rnnor(nr, r)
	Mint            nr;
	Mfloat          r[];
#endif
{
	Mint            i, ner;


	if (nr <= 0) {
		imsl_e1psh("l_rnnor");
		ner = 1;
		imsl_c1iarg(nr, "nr", 1, 0, &ner);
		imsl_e1pop("l_rnnor");
	} else {
		/* Get NR random uniforms */
		imsl_f_random_uniform(nr, IMSL_RETURN_USER, r, 0);
		/* Transform each uniform deviate */
		for (i = 0; i < nr; i++) {
			r[i] = imsl_f_normal_inverse_cdf(r[i]);
		}
	}
	return;
}				/* end of function */
