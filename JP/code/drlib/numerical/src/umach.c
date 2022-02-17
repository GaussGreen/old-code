#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 01/08/90 at 17:17:34
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  UMACH (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    March 21, 1984

    Purpose:    Set or retrieve input or output device unit numbers.

    Usage:      CALL UMACH (N, NUNIT)

    Arguments:
       N      - Index of desired unit.  (Input)
                The values of N are defined as follows:
                N = 1, corresponds to the standard input unit.
                N = 2, corresponds to the standard output unit.
       NUNIT  - I/O unit.  (Input or Output)
                If the value of N is negative, the unit corresponding
                to the index is reset to the value given in NUNIT.
                Otherwise, the value corresponding to the index is
                returned in NUNIT.

    GAMS:       R1

    Chapters:   MATH/LIBRARY Reference Material
                STAT/LIBRARY Reference Material
                SFUN/LIBRARY Reference Material

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

static Mint	lv_init = 0;
static FILE	*lv_unit[] = {NULL, NULL, NULL};

#ifdef ANSI
void imsl_umach(Mint n, FILE **nunit)
#else
void imsl_umach(n, nunit)
        FILE           **nunit;
	Mint             n;
#endif
{
	Mint             nn;

	if (!lv_init) {
	    lv_unit[0] = stdin;
	    lv_unit[1] = stdout;
	    lv_unit[2] = stderr;
	    lv_init = 1;
	}

	nn = abs(n);
	if (nn < 1 || nn > 3) {
		/* ERROR.  INVALID RANGE FOR N. */
	    imsl_e1sti(1,n);
/*	    (IMSL_TERMINAL, 10000, "The absolute value of the index variable must be 1 or 2.  N = %(I1).");
*/
            imsl_ermes(IMSL_TERMINAL, IMSL_INDEX_VARIABLE_VALUE);
	} else if (n < 0) {	/*  RESET */
		lv_unit[nn-1] = *nunit;
	} else {	/*  RETRIEVE */
		*nunit = lv_unit[n-1];
	}
	return;
}
