#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* Structured by FOR_STRUCT, v0.2, on 05/23/90 at 14:56:08
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  INITS (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1984

    Purpose:    Initialize the orthogonal series so the function value
                is the number of terms needed to insure the error is no
                larger than the requested accuracy.

    Usage:      INITS(OS, NOS, ETA)

    Arguments:
       OS     - Vector of length NOS containing coefficients in an
                orthogonal series.  (Input)
       NOS    - Number of coefficients in OS.  (Input)
       ETA    - Requested accuracy of the series.  (Input)
       INITS  - Number of terms needed to insure the error is no larger
                than ETA.  (Output)

    Remark:
       ETA will usually be chosen to be one tenth of machine precision.

    GAMS:       C7b

    Chapter:    SFUN/LIBRARY Miscellaneous Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mint imsl_inits(Mfloat *os, Mint nos, Mfloat eta)
#else
Mint imsl_inits(os, nos, eta)
	Mfloat           os[];
	Mint             nos;
	Mfloat           eta;
#endif
{
	Mint             i, ii, inits_v;
	Mfloat           err;


	imsl_e1psh("imsl_inits");
	inits_v = 0;

	if (nos < 1) {
		imsl_e1sti(1, nos);
                imsl_ermes(IMSL_TERMINAL, IMSL_CHEBY_EVAL_TERMS);
	} else {
		err = F_ZERO;
		for (ii = 1; ii <= nos; ii++) {
			i = nos + 1 - ii;
			err += fabs(os[i - 1]);
			if (err > eta) break;
		}
		if (i == nos) {
                    imsl_ermes(IMSL_TERMINAL, IMSL_CHEBY_EVAL_TOL);
		} else {
			inits_v = i;
		}
	}

	imsl_e1pop("imsl_inits");
	return (inits_v);
}				/* end of function */
