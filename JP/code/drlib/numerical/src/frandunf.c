#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#include <time.h>
#ifdef COMPUTER_PMXUX
#include <sys/types.h>
#include <sys/timeb.h>
#endif
#ifdef COMPUTER_MIPS
#include <bsd/sys/types.h>
#include <bsd/sys/timeb.h>
#endif
static Mfloat   *lv_uniform;

static VA_LIST_HACK PROTO(l_random_uniform,(Mint, va_list));
static void PROTO(l_rnun,(Mint, Mfloat[]));
static void l_r1ins();

static Mfloat lv_wk[128];

#ifdef DOUBLE
extern
#endif 
struct imsl_random_common_structure {
	double          d2p31a, dseed;/*, dwk[128];*/
	Mint             dinttb;
	Mint             indctr;
	Mint             inttb;
/*	Mfloat           wk[128];*/
} imsl_random_common;

#ifdef ANSI
Mfloat *imsl_f_random_uniform(Mint n_numbers, ...)
#else
Mfloat *imsl_f_random_uniform(n_numbers, va_alist)
    Mint        n_numbers;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr, n_numbers);

    E1PSH("imsl_f_random_uniform", "imsl_d_random_uniform");

    lv_uniform = NULL;
    IMSL_CALL(l_random_uniform(n_numbers, argptr));
    va_end(argptr);

    E1POP("imsl_f_random_uniform", "imsl_d_random_uniform");

    return(lv_uniform);
}

#ifdef ANSI
static VA_LIST_HACK l_random_uniform(Mint n_numbers, va_list argptr)
#else
static VA_LIST_HACK l_random_uniform(n_numbers, argptr)
    Mint        n_numbers;
    va_list     argptr;
#endif
{
    Mint        code, user_uniform = 0, ner, arg_number = 2;

    code = va_arg(argptr, Mint);
    if (code == (Mint)IMSL_RETURN_USER) {
        lv_uniform = va_arg(argptr, Mfloat*);
        user_uniform = 1;
    }
    else if (code != 0){
        imsl_e1sti (1, code);
        imsl_e1sti (2, arg_number);
        imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
        goto RETURN;
    }
    if (n_numbers <= 0) {
        ner = 1;
        imsl_c1iarg(n_numbers, "n_random", 1, 0, &ner);
        goto RETURN;
      }

    if (!user_uniform) {
        lv_uniform = (Mfloat *) imsl_malloc (n_numbers*sizeof(Mfloat));
        if (!lv_uniform){
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
            goto RETURN;

        }
    }
    l_rnun (n_numbers, lv_uniform);
    if (imsl_n1rty(0)>3 && !user_uniform) {
        imsl_free(lv_uniform);
        lv_uniform = NULL;
    }

RETURN:

    return(argptr);
}

/*
  -----------------------------------------------------------------------
    IMSL Name:  RNUN/DRNUN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 13, 1989

    Purpose:    Generate pseudorandom numbers from a uniform (0,1)
                distribution.

    Usage:      CALL RNUN (NR, R)

    Arguments:
       NR     - Number of random numbers to generate.  (Input)
       R      - Vector of length NR containing the random uniform (0,1)
                deviates.  (Output)

    Remark:
       The IMSL routine RNSET can be used to initialize the seed of the
       random number generator.  The routine RNOPT can be used to select
       the form of the generator.

    Keywords:   Utilities; Monte Carlo; Simulation; Rectangular
                distribution

    GAMS:       L6a21

    Chapters:   STAT/LIBRARY Random Number Generation
                MATH/LIBRARY Utilities

    Copyright:  1987 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_rnun(Mint nr, Mfloat r[])
#else
static void l_rnun(nr, r)
	Mint            nr;
	Mfloat          r[];
#endif
{
	Mint            i, j, ner;
	Mfloat          x;
	double          dseed1, dseed2;
	static int first = 1;
	static double   d2p31m = 2147483647.0e0;


        ner = 0;
	if (nr <= 0) {
		imsl_e1psh("l_rnun");
		ner = 1;
		imsl_c1iarg(nr, "nr", 1, 0, &ner);
		imsl_e1pop("l_rnun");
	} else {
		if (first) {
			imsl_r1int(0);
			first = 0;
		}
		if (imsl_random_common.indctr == 1) {
			/* Use the multiplier 16807 (=7**5) */
			for (i = 1; i <= nr; i++) {
				imsl_random_common.dseed = fmod(16807.0e0 * imsl_random_common.dseed, d2p31m);
				r[i - 1] = imsl_random_common.dseed / imsl_random_common.d2p31a;
			}
		} else if (imsl_random_common.indctr == 2) {
			/*
			 * Use the multiplier 16807 with shuffling
			 */
			if (!imsl_random_common.inttb) {
				/* Initialize table */
				l_r1ins();
				imsl_random_common.inttb = 1;
			}
			for (i = 1; i <= nr; i++) {
				/* Generate a new Uniform (0,1) deviate */
				imsl_random_common.dseed = fmod(16807.0e0 * imsl_random_common.dseed, d2p31m);
				/* Now shuffle */
				j = fmod(imsl_random_common.dseed, 128.0e0) + F_ONE;
				x = imsl_random_common.dseed / imsl_random_common.d2p31a;
				r[i - 1] = lv_wk[j - 1];
				lv_wk[j - 1] = x;
			}
		} else if (imsl_random_common.indctr == 3) {
			/* Use the multiplier 397204094 */
			for (i = 1; i <= nr; i++) {
				dseed1 = fmod(32768.0e0 * imsl_random_common.dseed, d2p31m);
				dseed2 = fmod(23166.0e0 * imsl_random_common.dseed, d2p31m);
				imsl_random_common.dseed = fmod(12121.0e0 * dseed1 + dseed2, d2p31m);
				r[i - 1] = imsl_random_common.dseed / imsl_random_common.d2p31a;
			}
		} else if (imsl_random_common.indctr == 4) {
			/*
			 * Use the multiplier 397204094 with shuffling
			 */
			if (!imsl_random_common.inttb) {
				/* Initialize table */
				l_r1ins();
				imsl_random_common.inttb = 1;
			}
			for (i = 1; i <= nr; i++) {
				/* Generate a new Uniform (0,1) deviate */
				dseed1 = fmod(32768.0e0 * imsl_random_common.dseed, d2p31m);
				dseed2 = fmod(23166.0e0 * imsl_random_common.dseed, d2p31m);
				imsl_random_common.dseed = fmod(12121.0e0 * dseed1 + dseed2, d2p31m);
				/* Now shuffle */
				j = fmod(imsl_random_common.dseed, 128.0e0) + F_ONE;
				x = imsl_random_common.dseed / imsl_random_common.d2p31a;
				r[i - 1] = lv_wk[j - 1];
				lv_wk[j - 1] = x;
			}
		} else if (imsl_random_common.indctr == 5) {
			/* Use the multiplier 950706376 */
			for (i = 1; i <= nr; i++) {
				dseed1 = fmod(32768.0e0 * imsl_random_common.dseed, d2p31m);
				dseed2 = fmod(8392.0e0 * imsl_random_common.dseed, d2p31m);
				imsl_random_common.dseed = fmod(29013.0e0 * dseed1 + dseed2, d2p31m);
				r[i - 1] = imsl_random_common.dseed / imsl_random_common.d2p31a;
			}
		} else if (imsl_random_common.indctr == 6) {
			/*
			 * Use the multiplier 950706376 with shuffling
			 */
			if (!imsl_random_common.inttb) {
				/* Initialize table */
				l_r1ins();
				imsl_random_common.inttb = 1;
			}
			for (i = 1; i <= nr; i++) {
				/* Generate a new Uniform (0,1) deviate */
				dseed1 = fmod(32768.0e0 * imsl_random_common.dseed, d2p31m);
				dseed2 = fmod(8392.0e0 * imsl_random_common.dseed, d2p31m);
				imsl_random_common.dseed = fmod(29013.0e0 * dseed1 + dseed2, d2p31m);
				/* Now shuffle */
				j = fmod(imsl_random_common.dseed, 128.0e0) + F_ONE;
				x = imsl_random_common.dseed / imsl_random_common.d2p31a;
				r[i - 1] = lv_wk[j - 1];
				lv_wk[j - 1] = x;
			}
		}
	}
	return;
}				/* end of function */


/*
  -----------------------------------------------------------------------
    IMSL Name:  R1INS/DR1INS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 13, 1989

    Purpose:    Randomly initialize the table for the IMSL shuffled
                random number generators.

    Usage:      CALL R1INS

    Arguments:  (none)

    GAMS:       L6c

    Chapter:    STAT/LIBRARY Random Number Generation (not documented)

    Copyright:  1987 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
static void l_r1ins()
{
	Mint            i;
	Mdouble          dseed1, dseed2;
	static Mdouble   d2p31m = 2147483647.0e0;
	static Mint first = 1;



	if (first) {
		imsl_r1int(0);
		first = 0;
	}
	if (imsl_random_common.indctr == 2) {
		/* Use the multiplier 16807 */
		for (i = 1; i <= 128; i++) {
			imsl_random_common.dseed = fmod(16807.0e0 * imsl_random_common.dseed, d2p31m);
			lv_wk[i - 1] = imsl_random_common.dseed / imsl_random_common.d2p31a;
		}
	} else if (imsl_random_common.indctr == 4) {
		/* Use the multiplier 397204094 */
		for (i = 1; i <= 128; i++) {
			dseed1 = fmod(32768.0e0 * imsl_random_common.dseed, d2p31m);
			dseed2 = fmod(23166.0e0 * imsl_random_common.dseed, d2p31m);
			imsl_random_common.dseed = fmod(12121.0e0 * dseed1 + dseed2, d2p31m);
			lv_wk[i - 1] = imsl_random_common.dseed / imsl_random_common.d2p31a;
		}
	} else if (imsl_random_common.indctr == 6) {
		/* Use the multiplier 950706376 */
		for (i = 1; i <= 128; i++) {
			dseed1 = fmod(32768.0e0 * imsl_random_common.dseed, d2p31m);
			dseed2 = fmod(8392.0e0 * imsl_random_common.dseed, d2p31m);
			imsl_random_common.dseed = fmod(29013.0e0 * dseed1 + dseed2, d2p31m);
			lv_wk[i - 1] = imsl_random_common.dseed / imsl_random_common.d2p31a;
		}
	}
	return;
}				/* end of function */


/*
  -----------------------------------------------------------------------
    IMSL Name:  RNSES/DRNSES (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 19, 1988

    Purpose:    Initialize the table in the IMSL random number generators
                that use shuffling.

    Usage:      CALL RNSES (TABLE)

    Argument:
       TABLE  - Vector of length 128 used in the IMSL random number
                generators that use shuffling.  (Input)
                The values in TABLE should be uniform (0,1) random
                deviates.  Except for the element which may be
                nonpositive, all elements of TABLE must be positive and
                strictly less than 1.0.  If TABLE(1) .LE. 0.0, then
                the values in TABLE are not stored in the shuffling
                table; instead, the first 128 numbers resulting from the
                next call to a shuffled generator are used to initialize
                the table (prior to any shuffling).

    Remark:
       This routine is generally used when it is desired to restart a
       simulation that uses a shuffled generator.  At the end of one
       simulation run, RNGES is used to obtain the current value of TABLE
       and then RNSES is used at the beginning of the next simulation run
       to restore TABLE.

    Keyword:    Utility

    GAMS:       L6c

    Chapter:    STAT/LIBRARY Random Number Generation

    Copyright:  1987 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_rnses(Mfloat table[])
#else
void imsl_rnses(table)
	Mfloat          table[];
#endif
{
	Mint            i;


	if (table[0] <= F_ZERO) {
		imsl_random_common.inttb = 0;
	} else {
		for (i = 1; i <= 128; i++) {
			if (table[i - 1] <= F_ZERO || table[i - 1] >= F_ONE) {
				E1PSH("imsl_rnses","imsl_drnses");
#ifdef DOUBLE
				imsl_e1std(1, table[i - 1]);
#else
				imsl_e1str(1, table[i - 1]);
#endif
				imsl_e1sti(1, i);

/*				imsl_ermes(5, 1, "TABLE(%(i1)) = %(r1).  All elements of TABLE must be positive and strictly less than 1.0.");
*/
                                imsl_ermes(IMSL_TERMINAL,
				IMSL_STRICTLY_POS_TABLE_ELMNTS);
				E1POP("imsl_rnses","imsl_drnses");
			} else {
				lv_wk[i - 1] = table[i - 1];
			}
		}
		imsl_random_common.inttb = 1;
	}
	return;
}				/* end of function */

#ifndef DOUBLE
/*
  -----------------------------------------------------------------------
    IMSL Name:  RNGET (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    April 19, 1988

    Purpose:    Retrieve the current value of the seed used in the IMSL
                random number generators.

    Usage:      CALL RNGET (ISEED)

    Argument:
       ISEED  - The seed of the random number generator.  (Output)
                ISEED is in the range (1, 2147483646).

    Keyword:    Utilities

    GAMS:       L6c

    Chapters:   STAT/LIBRARY Random Number Generation
                MATH/LIBRARY Utilities

    Copyright:  1987 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
Mint imsl_random_seed_get()
{
	imsl_r1int(0);
	return((Mint)(imsl_random_common.dseed+F_HALF));
}				/* end of function */

/*
  -----------------------------------------------------------------------
    IMSL Name:  RNSET (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    April 19, 1988

    Purpose:    Initialize a random seed for use in the IMSL random
                number generators.

    Usage:      CALL RNSET (ISEED)

    Argument:
       ISEED  - The seed of the random number generator.  (Input)
                ISEED must be in the range (0, 2147483646).  If ISEED is
                zero, a value is computed using the system clock; and,
                hence, the results of programs using the IMSL random
                number generators will be different at different times.

    Keyword:    Utilities

    GAMS:       L6c

    Chapters:   STAT/LIBRARY Random Number Generation
                MATH/LIBRARY Utilities

    Copyright:  1987 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_random_seed_set(Mint iseed)
#else
void imsl_random_seed_set(iseed)
	Mint            iseed;
#endif
{
	Mint            iiseed;


	if (iseed < 0 || iseed >= 2147483647) {
		imsl_r1int(0);
	} else if (iseed == 0) {
		imsl_r1clk(&iiseed);
		imsl_random_common.dseed = (double) iiseed;
	} else {
		imsl_r1int(-1);
		imsl_random_common.dseed = (double) iseed;
	}
	return;
}				/* end of function */


/*
  -----------------------------------------------------------------------
    IMSL Name:  R1CLK (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    July 11, 1985

    Purpose:    Set a value in the range (1, 2147483646) using the
                system clock.

    Usage:      CALL R1CLK (IISEED)

    Argument:
       IISEED - A value in the range (1, 2147483646) computed from a
                call to the system clock.  (Output)

    GAMS:       L6c

    Chapter:    STAT/LIBRARY Not-user-callable

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_r1clk(Mint *iiseed)
#else
void imsl_r1clk(iiseed)
	Mint            *iiseed;
#endif
{
        struct tm	*tim;
        static time_t    clock;
#ifndef COMPUTER_SGRUXS
	time_t           time();
#endif
	/* Set seed using the system clock */

	clock = time(&clock);
        tim = gmtime(&clock);
        if (tim == NULL) tim = localtime(&clock);
	*iiseed = tim->tm_hour * 10000 + tim->tm_min * 100 + tim->tm_sec + 1;
	if (*iiseed <= 1) {
		/*
		 * Seed is nonrandomly set.  Use of the system clock has not
		 * been devloped.
		 */
		*iiseed = 123457;
	}
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  R1INT (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    April 19, 1988

    Purpose:    Initialize variables in common used by the IMSL
                multiplicative congruential random number generators.

    Usage:      CALL R1INT (IOPT)

    Argument:
       IOPT   - Option indicator.  (Input)
                IOPT       ACTION
                Positive   INDCTR in COMMON is set to IOPT.
                 0         A check is made to see if the seed has been
                           set.  If a seed has not been set, it is
                           initialized using the system clock.
                -1         The seed is being set in the calling routine.
                -2         No action, other than the initialization (on
                           the first call only) of the variables in
                           COMMON that the user has no control over.

    Remark:
       On the first call to R1INT, the default indicator is set and the
       scaling factor for (1,2**31-1) into (0,1) is determined.

    GAMS:       L6c

    Chapter:    STAT/LIBRARY Random Number Generation (not documented)

    Copyright:  1987 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_r1int(Mint iopt)
#else
void imsl_r1int(iopt)
	Mint            iopt;
#endif
{
	Mint            iiseed;
	Mfloat          temp;
	static Mint first = 1;
	static Mint setsed = 0;



	if (first) {
		/* Determine scaling factor for U(0,1) */
		imsl_random_common.d2p31a = 2147483647.0e0;
		temp = 2147483646.0e0 / imsl_random_common.d2p31a;
		if (temp >= F_ONE) {
			imsl_random_common.d2p31a = 2147483655.0e0;
			temp = 2147483646.0e0 / imsl_random_common.d2p31a;
			if (temp >= F_ONE) {
				imsl_random_common.d2p31a = 2147483711.0e0;
		L_10:
				;
				temp = 2147483646.0e0 / imsl_random_common.d2p31a;
				if (temp >= F_ONE) {
					imsl_random_common.d2p31a += F_EIGHT;
					goto L_10;
				}
			}
		}
		imsl_random_common.indctr = 1;
		imsl_random_common.inttb = 0;
		imsl_random_common.dinttb = 0;
		first = 0;
	}
	if (iopt == -1) {
		/* Seed is being set by calling routine */
		setsed = 1;
	} else if (iopt == 0) {
		if (!setsed) {
			/* Set seed using the system clock */
			imsl_r1clk(&iiseed);
			imsl_random_common.dseed = (double) iiseed;
			setsed = 1;
		}
	} else if (iopt > 0) {
		/*
		 * Select the Uniform random number generator
		 */
		imsl_random_common.indctr = iopt;
	}
	return;
}				/* end of function */

/*
  -----------------------------------------------------------------------
    IMSL Name:  RNOPT (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    July 15, 1985

    Purpose:    Select the uniform (0,1) multiplicative congruential
                pseudorandom number generator.

    Usage:      CALL RNOPT (IOPT)

    Argument:
       IOPT   - Indicator of the generator.  (Input)
                The random number generator is a multiplicative
                congruential generator with modulus 2**31 - 1.  IOPT
                is used to choose the multiplier and whether or not
                shuffling is done.
                IOPT   Generator
                  1    The multiplier 16807 is used.
                  2    The multiplier 16807 is used with shuffling.
                  3    The multiplier 397204094 is used.
                  4    The multiplier 397204094 is used with shuffling.
                  5    The multiplier 950706376 is used.
                  6    The multiplier 950706376 is used with shuffling.

    Keyword:    Utilities

    GAMS:       L6c

    Chapters:   STAT/LIBRARY Random Number Generation
                MATH/LIBRARY Utilities

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void  imsl_random_option(Mint iopt)
#else
void  imsl_random_option(iopt)
	Mint            iopt;
#endif
{
	Mint            ner = 1;

	imsl_e1psh("imsl_random_option");
	imsl_c1iarg(iopt, "generator_option", 1, 6, &ner);
	if (imsl_n1rcd(0) == 0)
		imsl_r1int(iopt);

	imsl_e1pop("imsl_random_option");

	return;
}				/* end of function */
#endif  /* end ifndef DOUBLE */
