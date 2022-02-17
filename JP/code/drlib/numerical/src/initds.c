#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/*
  -----------------------------------------------------------------------
    IMSL Name:  INITDS (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1984

    Purpose:    Initialize the orthogonal series so the function value
                is the number of terms needed to insure the error is no
                larger than the requested accuracy.

    Usage:      INITDS(DOS, NOS, ETA)

    Arguments:
       DOS    - Double precision vector of length NOS containing
                coefficients in an orthogonal series.  (Input)
       NOS    - Number of coefficients in OS.  (Input)
       ETA    - Requested accuracy of the series.  (Input)
       INITDS - Number of terms needed to insure the error is no larger
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
Mint imsl_initds(Mdouble dos[], Mint nos, Mdouble eta)
#else
Mint imsl_initds(dos, nos, eta)
	Mdouble         dos[];
	Mint            nos;
	Mdouble         eta;
#endif
{
	Mint            i, ii, initds_v;
	float           err;


	imsl_e1psh("imsl_initds");
	initds_v = 0;

	if (nos < 1) {
		imsl_e1sti(1, nos);
/*		(5, 5, "The number of coefficients is less than 1. NOS = %(i1).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_CHEBY_EVAL_TERMS);
	} else {
		err = F_ZERO;
		for (ii = 1; ii <= nos; ii++) {
			i = nos + 1 - ii;
			err += fabs((float) (dos[i - 1]));
			if (err > eta)
				goto L_20;
		}

L_20:
		if (i == nos) {

/*			(5, 6, "Too much accuracy may be requested. ETA should be increased.");
*/
                        imsl_ermes(IMSL_TERMINAL, IMSL_CHEBY_EVAL_TOL);
		} else {
			initds_v = i;
		}
	}

	imsl_e1pop("imsl_initds");
	return (initds_v);
}				/* end of function */
