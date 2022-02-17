#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 07/26/90 at 14:45:21
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  IMACH (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    March 26, 1984

    Purpose:    Retrieve integer machine constants.

    Usage:      IMACH(N)

    Arguments:
       N      - Index of desired constant.  (Input)
       IMACH  - Machine constant.  (Output)

    Remark:
       Following is a description of the assorted integer machine
       constants.

       Words

          IMACH( 1) = Number of bits per integer storage unit.
          IMACH( 2) = Number of characters per integer storage unit.

       Integers

          Assume integers are represented in the S-DIGIT, BASE-A form
          SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
          where 0 .LE. X(I) .LT. A for I=0,...,S-1.  Then

          IMACH( 3) = A, the base.
          IMACH( 4) = S, number of BASE-A digits.
          IMACH( 5) = A**S - 1, largest magnitude.

       Floating-point numbers

          Assume floating-point numbers are represented in the T-DIGIT,
          BASE-B form SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
          where 0 .LE. X(I) .LT. B for I=1,...,T,
          0 .LT. X(1), and EMIN .LE. E .LE. EMAX.  Then

          IMACH( 6) = B, the base.

          Single precision

             IMACH( 7) = T, number of BASE-B digits.
             IMACH( 8) = EMIN, smallest exponent E.
             IMACH( 9) = EMAX, largest exponent E.

          Double precision

             IMACH(10) = T, number of BASE-B digits.
             IMACH(11) = EMIN, smallest exponent E.
             IMACH(12) = EMAX, largest exponent E.

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

#ifdef ANSI
Mint imsl_imach(Mint n)
#else
Mint imsl_imach(n)
	Mint            n;
#endif
{
	Mint             imach_v;
	static Mint      imachv[12];
	static Mint      _aini = 1;
	if (_aini) {
		imachv[0] = 32;
		imachv[1] = 4;
		imachv[2] = 2;
		imachv[3] = 31;
		imachv[4] = 2147483647;
		imachv[5] = 2;
		imachv[6] = 24;
		imachv[7] = -125;
		imachv[8] = 128;
		imachv[9] = 53;
		imachv[10] = -1021;
		imachv[11] = 1024;
		_aini = 0;
	}
	/* DEFINE CONSTANTS */

	if (n < 1 || n > 12) {
		/* ERROR.  INVALID RANGE FOR N. */
                imsl_e1sti(1, 1);
                imsl_e1sti(2, 12);
                imsl_e1sti(3, n);
                imsl_e1stl(1,"n");
                imsl_ermes(IMSL_TERMINAL, IMSL_INTEGER_OUT_OF_RANGE);
	} else {
		imach_v = imachv[n - 1];
	}

	return (imach_v);
}				/* end of function */
