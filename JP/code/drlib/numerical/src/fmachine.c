#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

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

#ifdef ANSI
Mfloat imsl_f_machine(Mint n)
#else
Mfloat imsl_f_machine(n)
    Mint n;
#endif
{
    if (n<1 || n>8) {
            /* %(L1) must be between %(I1) and %(I2), but %(L1) = %(I3). */
	E1PSH("imsl_f_machine", "imsl_d_machine");
        imsl_e1stl(1, "n");
	imsl_e1sti(1, 1);
	imsl_e1sti(2, 8);
	imsl_e1sti(3, n);
        imsl_ermes(IMSL_TERMINAL, IMSL_INTEGER_OUT_OF_RANGE);
	E1POP("imsl_f_machine", "imsl_d_machine");
	return 0;
    }
#ifdef DOUBLE
    return imsl_machine.d[n-1];
#else
    return imsl_machine.f[n-1];
#endif
}
