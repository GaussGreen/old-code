#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 05/23/90 at 14:55:13
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CSEVL (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1984

    Purpose:    Evaluate the N-term Chebyshev series.

    Usage:      CSEVL(X, CS, N)

    Arguments:
       X      - Argument at which the series is to be evaluated.  (Input)
       CS     - Vector of length N containing the terms of a Chebyshev
                series.  (Input)
                In evaluating CS, only half of the first coefficient is
                summed.
       N      - Number of terms in the vector CS.  (Input)
       CSEVL  - Function value.  (Output)

    Remark:
       Informational error:
       Type Code
         3   7  X is outside the interval (-1.1,+1.1)

    GAMS:       C19

    Chapter:    SFUN/LIBRARY Miscellaneous Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_csevl(Mfloat x, Mfloat *cs, Mint n)
#else
Mfloat imsl_csevl(x, cs, n)
	Mfloat           x, cs[];
	Mint             n;
#endif
{
	Mint             i, ni;
	Mfloat           b0, b1, b2, csevl_v, twox;


	imsl_e1psh("imsl_csevl");
	csevl_v = imsl_amach(6);

	if (n < 1) {
		imsl_e1sti(1, n);
                imsl_ermes(IMSL_TERMINAL, IMSL_CHEBY_NEG_TERMS);
	} else if (n > 1000) {
		imsl_e1sti(1, n);
                imsl_ermes(IMSL_TERMINAL, IMSL_CHEBY_TOO_MANY_TERMS);
	} else {
		if (x < -1.1 || x > 1.1) {
			imsl_e1str(1, x);
                        imsl_ermes(IMSL_WARNING, IMSL_CHEBY_RANGE);
		}
		b1 = F_ZERO;
		b0 = F_ZERO;
		twox = x + x;
		for (i = 1; i <= n; i++) {
			b2 = b1;
			b1 = b0;
			ni = n + 1 - i;
			b0 = twox * b1 - b2 + cs[ni - 1];
		}

		csevl_v = F_HALF * (b0 - b2);
	}

	imsl_e1pop("imsl_csevl");
	return (csevl_v);
}
