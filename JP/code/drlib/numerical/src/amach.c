#include "imsl_inc.h"


#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* nada */

/*-----------------------------------------------------------------------
    IMSL Name:  AMACH (Single precision version)

    Computer:   SUN/SINGLE

    Revised:    January 8, 1990

    Purpose:    Retrieve single-precision machine constants.

    Usage:      AMACH(N)

    Arguments:
       N      - Index of desired constant.  (Input)
       AMACH  - Machine constant.  (Output)
                AMACH(1) = B**(EMIN-1), the smallest positive magnitude.
                AMACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
                AMACH(3) = B**(-T), the smallest relative spacing.
                AMACH(4) = B**(1-T), the largest relative spacing.
                AMACH(5) = LOG10(B), the log, base 10, of the radix.
                AMACH(6) = not-a-number.
                AMACH(7) = positive machine infinity.
                AMACH(8) = negative machine infinity.

    GAMS:       R1

    Chapters:   MATH/LIBRARY Reference Material
                STAT/LIBRARY Reference Material
                SFUN/LIBRARY Reference Material

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#undef imsl_amach

#ifdef ANSI
Mfloat imsl_amach(Mint n)
#else
Mfloat imsl_amach(n)
    Mint n;
#endif
{
    if (n<1 || n>8) {
	imsl_e1psh("imsl_amach");
	imsl_e1sti(1, n);
	imsl_ermes(IMSL_TERMINAL, IMSL_ARG_OUT_OF_RANGE);
	imsl_e1pop("imsl_amach");
	return F_ZERO;
    }

    return imsl_machine.f[n-1];
}
