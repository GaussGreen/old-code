#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  FIN/DFIN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Evaluate the inverse of the F distribution function.

    Usage:      FIN(P, DFN, DFD)

    Arguments:
       P      - Probability for which the inverse of the F
                distribution function is to be evaluated.  (Input)
                P must be in the open interval (0.0, 1.0).
       DFN    - Numerator degrees of freedom.  (Input)
                DFN must be positive.
       DFD    - Denominator degrees of freedom.  (Input)
                DFD must be positive.
       FIN    - Function value.  (Output)
                The probability that an F random variable takes a value
                less than or equal to FIN is P.

    Remark:
       Informational error
       Type Code
         4   3  FIN is set to machine infinity since overflow would occur
                upon modifying the inverse value for the F distribution
                with the result obtained from the inverse BETA
                distribution.

    Keywords:   Percentage point; Snedecor's F; Variance ratio;
                Probability distribution; Continuous random variables;
                Percentile; Fractile

    GAMS:       L5a2f

    Chapters:   STAT/LIBRARY Probability Distribution Functions and
                             Inverses
                SFUN/LIBRARY Probability Distribution Functions and
                             Inverses

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_F_inverse_cdf(Mfloat p, Mfloat dfn, Mfloat dfd)
#else
Mfloat imsl_f_F_inverse_cdf(p, dfn, dfd)
	Mfloat          p, dfn, dfd;
#endif
{
	Mfloat          a, b, fin_v, pp, rdelp, x;


	/*
	 * NOTE:  This code calls BETIN. Please look at the NOTE in BETIN.
	 */
	E1PSH ("imsl_f_F_inverse_cdf","imsl_d_F_inverse_cdf");
	fin_v = imsl_amach(6);
	/* Check P */
	if (p <= F_ZERO || p >= F_ONE) {
		imsl_e1str(1, p);

/*		imsl_ermes(5, 1, "The input probability, P = %(r1), is not in the OPEN interval (0.0, 1.0).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_P_OUTSIDE_OPEN_INTERVAL);
		goto L_9000;
	}
	/* Check DFN and DFD */
	if (dfn <= F_ZERO || dfd <= F_ZERO) {
		imsl_e1str(1, dfn);
		imsl_e1str(2, dfd);

/*		imsl_ermes(5, 2, "The input value for the numerator or denominator degrees of freedom for the F distribution, DFN = %(r1), DFD = %(r2), must be positive.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_DFN_OR_DFD_IS_NEGATIVE);
		goto L_9000;
	}
	/*
	 * Intialize a number slightly smaller than 1.0
	 */
	rdelp = F_ONE - imsl_amach(4);
	/* Reparameterize for BETA usage */
	a = F_HALF * dfn;
	b = F_HALF * dfd;
	if (p <= F_HALF) {
		x = imsl_betin(p, a, b);
		if (x >= rdelp) {
			fin_v = imsl_amach(7);

/*			imsl_ermes(4, 3, "Upon modifying the inverse value for the F distribution from the value obtained from IMSL inverse BETA distribution, BETIN, overflow would occur.  Therefore, the value is set to machine infinity.");
*/
                        imsl_ermes(IMSL_FATAL, IMSL_F_INVERSE_OVERFLOW);
			goto L_9000;
		} else {
			/*
			 * Modify BETA output for X from F distribution
			 */
			x = dfd * x / (dfn * (F_ONE - x));
		}
	} else {
		pp = F_ONE - p;
		x = imsl_betin(pp, b, a);
		if (x == F_ZERO) {
			fin_v = imsl_amach(7);

/*			imsl_ermes(4, 3, "Upon modifying the inverse value for the F distribution from the value obtained from IMSL inverse BETA distribution, BETIN, overflow would occur.  Therefore, the value is set to machine infinity.");
*/
                        imsl_ermes(IMSL_FATAL, IMSL_F_INVERSE_OVERFLOW);
			goto L_9000;
		} else {
			x = (F_ONE / x - F_ONE) * dfd / dfn;
		}
	}

	fin_v = x;

L_9000:
	E1POP ("imsl_f_F_inverse_cdf","imsl_d_F_inverse_cdf");

	return (fin_v);
}				/* end of function */
