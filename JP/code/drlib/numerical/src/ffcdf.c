#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  FDF/DFDF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 20, 1990

    Purpose:    Evaluate the F distribution function.

    Usage:      FDF(F, DFN, DFD)

    Arguments:
       F      - Argument for which the F distribution function is to
                be evaluated.  (Input)
       DFN    - Numerator degrees of freedom.  (Input)
                DFN must be positive.
       DFD    - Denominator degrees of freedom.  (Input)
                DFD must be positive.
       FDF    - Function value, the probability that an F random
                variable takes a value less than or equal to the input F.
                (Output)

    Remark:
      Information error
       Type Code
         1   2  Since the input argument F is not positive, the
                distribution function is zero at F.

    Keywords:   P-value; Probability integral; Snedecor's F; Variance
                ratio; Probability distribution; Continuous random
                variables; Cumulative distribution function; CDF

    GAMS:       L5a1f

    Chapters:   STAT/LIBRARY Probability Distribution Functions and
                             Inverses
                SFUN/LIBRARY Probability Distribution Functions and
                             Inverses

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_F_cdf(Mfloat f, Mfloat dfn, Mfloat dfd)
#else
Mfloat imsl_f_F_cdf(f, dfn, dfd)
	Mfloat          f, dfn, dfd;
#endif
{
	Mfloat          fdf_v, fx, pcomp, pin, qin;


	E1PSH("imsl_f_F_cdf","imsl_d_F_cdf");

	if (dfn <= F_ZERO || dfd <= F_ZERO) {
		imsl_e1str(1, dfn);
		imsl_e1str(2, dfd);

/*		imsl_ermes(5, 1, "Both DFN = %(r1) and DFD = %(r2) must be greater than
0.0.");*/
                imsl_ermes(IMSL_TERMINAL, IMSL_DFN_OR_DFD_IS_NEGATIVE);
		fdf_v = imsl_amach(6);

	} else if (f <= F_ZERO) {
		imsl_e1str(1, f);

/*		imsl_ermes(1, 2, "Since F = %(r1) is not positive, the distribution function is zero at F.");
*/
                imsl_ermes(IMSL_NOTE, IMSL_DIST_FCN_SET_TO_ZERO);
		fdf_v = F_ZERO;

	} else {
		/* TRANSFORM F TO A BETA. */
		fx = dfd / (dfd + dfn * f);
		pin = F_HALF * dfd;
		qin = F_HALF * dfn;
		pcomp = imsl_f_beta_incomplete(fx, pin, qin);
		fdf_v = F_ONE - pcomp;
	}

	E1POP("imsl_f_F_cdf","imsl_d_F_cdf");
	return (fdf_v);
}				/* end of function */
