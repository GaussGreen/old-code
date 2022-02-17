#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#define	    TRUNC(A) ((Mfloat)((Mint)(A)))

/*
  -----------------------------------------------------------------------
    IMSL Name:  ALNGAM/DLNGAM (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 23, 1988

    Purpose:    Evaluate the logarithm of the absolute value of the
                gamma_c function.

    Usage:      ALNGAM(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       ALNGAM - Function value.  (Output)

    Remark:
       Informational error:
       Type Code
         3   2  Result of ALNGAM(X) is accurate to less than one half
                precision because X is too near a negative integer.

    GAMS:       C7a

    Chapter:    SFUN/LIBRARY Gamma Function and Related Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_log_gamma(Mfloat x)
#else
Mfloat imsl_f_log_gamma(x)
	Mfloat          x;
#endif
{
	Mfloat          alngam_v, sinpiy, y;
	static Mfloat   pi = 3.14159265358979323846264338328e0;
	/*
	 * sq2pil = alog (sqrt(2.*pi)) sqpi2l = alog(SQRT(PI/2.))
	 */
	static Mfloat   sq2pil = 0.918938533204672741780329736406e0;
	static Mfloat   sqpi2l = .225791352644727432363097614947e0;
	static Mfloat   xmax = 0.0;
	static Mfloat   dxrel = 0.0;

	E1PSH("imsl_f_log_gamma", "imsl_d_log_gamma");

	alngam_v = imsl_amach(6);

	if (xmax == F_ZERO) {
		xmax = imsl_amach(2) / log(imsl_amach(2));
		dxrel = sqrt(imsl_amach(4));
	}
	y = fabs(x);
	/*
	 * ALOG (ABS (GAMMA(X))) FOR ABS(X) .LE. 10.0
	 */
	if (y <= F_TEN) {
		alngam_v = log(fabs(imsl_f_gamma(x)));
		/*
		 * ALOG (ABS (GAMMA(X))) FOR ABS(X) .GT. 10.0
		 */
	} else if (y > xmax) {
                    /* The function overflows because ABS(%(r1)) is greater
                       than %(r2). */
		imsl_e1str(1, x);
		imsl_e1str(2, xmax);
		imsl_ermes(IMSL_FATAL, IMSL_LARGE_ABS_ARG_OVERFLOW);
	} else if (x > F_ZERO) {
		alngam_v = sq2pil + (x - F_HALF) * log(x) - x + imsl_r9lgmc(y);
	} else {
		sinpiy = fabs(sin(pi * y));
		if (sinpiy == F_ZERO) {
                            /* The argument for the function can not be a
                               negative integer. Argument X = %(r1). */
			imsl_e1str(1, x);
			imsl_ermes(IMSL_FATAL, IMSL_NEGATIVE_INTEGER);
		} else {
			alngam_v = sqpi2l + (x-F_HALF) * log(y) - x - 
                                log(sinpiy) - imsl_r9lgmc(y);
			if (fabs((x - TRUNC(x-F_HALF)) * alngam_v/x) < dxrel) {
                                /* The result is accurate to less than one half
                                   precision because X = %(r1) is too close to a
                                   negative integer. */
				imsl_e1str(1, x);
				imsl_ermes(IMSL_WARNING, IMSL_NEAR_NEG_INT_WARN);
			}
		}
	}
	E1POP("imsl_f_log_gamma", "imsl_d_log_gamma");
	return (alngam_v);
}
