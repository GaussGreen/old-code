/*-----------------------------------------------------------------------
    IMSL Name:  DMACH (Single precision version)

    Computer:   FORC/DOUBLE

    Revised:    January 8, 1990

    Purpose:    Retrieve single-precision machine constants.

    Usage:      DMACH(N)

    Arguments:
       N      - Index of desired constant.  (Input)
       DMACH  - Machine constant.  (Output)
                DMACH(1) = B**(EMIN-1), the smallest positive magnitude.
                DMACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
                DMACH(3) = B**(-T), the smallest relative spacing.
                DMACH(4) = B**(1-T), the largest relative spacing.
                DMACH(5) = LOG10(B), the log, base 10, of the radix.
                DMACH(6) = not-a-number.
                DMACH(7) = positive machine infinity.
                DMACH(8) = negative machine infinity.

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
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#undef imsl_dmach

#ifdef ANSI
Mdouble imsl_dmach(Mint n)
#else
Mdouble imsl_dmach(n)
    Mint n;
#endif
{
    if (n<1 || n>8) {
	imsl_e1psh("imsl_dmach");
	imsl_e1sti(1, n);
/*	imsl_ermes(5, 5, "The argument must be between 1 and 8 inclusive. N = %(I1)");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_ARG_OUT_OF_RANGE);
	imsl_e1pop("imsl_dmach");
	return F_ZERO;
    }
    return imsl_machine.d[n-1];
}
