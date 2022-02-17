#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif



/*Translated by FOR_C++, v0.1, on 09/09/91 at 14:56:58 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 09/09/91 at 14:56:56
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BETDF/DBETDF (Single/Double precision version)

    Computer:   SUN/SINGLE

    Revised:    May 29, 1991

    Purpose:    Evaluate the imsl_beta probability distribution function.

    Usage:      BETDF(X, PIN, QIN)

    Arguments:
       X      - Argument for which the imsl_beta distribution function is to
                be evaluated.  (Input)
       PIN    - First imsl_beta distribution parameter.  (Input)
                PIN must be positive.
       QIN    - Second imsl_beta distribution parameter.  (Input)
                QIN must be positive.
       BETDF  - Probability that a random variable from a imsl_beta
                distribution having parameters PIN and QIN will be less
                than or equal to X.  (Output)

    Remark:
       Informational errors
       Type Code
         1   1  Since the input argument X is less than or equal to
                zero, the distribution function is equal to zero at X.
         1   2  Since the input argument X is greater than or equal to
                one, the distribution function is equal to one at X.

    Keywords:   P-value; Probability integral; Probability distribution;
                Continuous random variables; Cumulative distribution
                function; CDF

    GAMS:       L5a1b; C7f

    Chapters:   STAT/LIBRARY Probability Distribution Functions and
                Inverses
                MATH/LIBRARY SPECIAL FUNCTIONS Probability Distribution
                Functions and Inverses

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_beta_cdf(Mfloat x, Mfloat pin, Mfloat qin)
#else
Mfloat imsl_f_beta_cdf(x, pin, qin)
	Mfloat           x, pin, qin;
#endif
{
	Mint             ner;
	Mfloat           betdf_v;


	E1PSH("imsl_f_beta_cdf", "imsl_d_beta_cdf");
	ner = 1;

	if (pin <= 0.0e0) {
		imsl_e1str(1, pin);

		imsl_e1mes(5, ner + 2, "PIN = %(r1) must be greater than 0.0.");
		betdf_v = imsl_amach(6);
	} else if (qin <= 0.0e0) {
		imsl_e1str(1, qin);

		imsl_e1mes(5, ner + 3, "QIN = %(r1) must be greater than 0.0.");
		betdf_v = imsl_amach(6);
	} else if (x <= 0.0e0) {
		imsl_e1str(1, x);

		imsl_e1mes(1, ner, "Since X = %(r1) is less than or equal to zero, the distribution function is zero at X.");
		betdf_v = 0.0e0;
	} else if (x >= 1.0e0) {
		imsl_e1str(1, x);

		imsl_e1mes(1, ner + 1, "Since X = %(r1) is greater than or equal to one, the distribution function is one at X.");
		betdf_v = 1.0e0;
	} else {
		betdf_v = imsl_f_beta_incomplete(x, pin, qin);
	}

	E1POP("imsl_f_beta_cdf", "imsl_d_beta_cdf");

	return (betdf_v);
}				/* end of function */
