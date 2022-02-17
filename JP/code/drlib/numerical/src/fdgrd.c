/*Translated by FOR_C++, v0.1, on 08/16/90 at 13:11:29 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
#include <stdio.h>
#include <math.h>
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef INCL_PROTOTYPES
#include "/home/usr2/imsl1/clib/newclib/include/imsl_int.h"
#endif
/* Structured by FOR_STRUCT, v0.2, on 08/16/90 at 13:11:27
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  FDGRD/DFDGRD  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Approximate the gradient using forward differences.

    Usage:      CALL FDGRD (FCN, N, XC, XSCALE, FC, EPSFCN, GC)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate the function to be
                minimized.  The usage is
                CALL FCN (N, X, F), where
                N      - Length of X.  (Input)
                X      - The point at which the function is evaluated.
                         (Input)
                         X should not be changed by FCN.
                F      - The computed function value at the point X.
                         (Output)
                FCN must be declared EXTERNAL in the calling program.
       N      - Dimension of the problem.  (Input)
       XC     - Vector of length N containing the point at which the
                gradient is to be estimated.  (Input)
       XSCALE - Vector of length N containing the diagonal scaling matrix
                for the variables.  (Input)
                In the absence of other information, set all entries
                to 1.0.
       FC     - Scalar containing the value of the function at XC.
                (Input)
       EPSFCN - Estimate of the relative noise in the function.  (Input)
                EPSFCN must be less than or equal to 0.1.  In the absence
                of other information, set EPSFCN to 0.0.
       GC     - Vector of length N containing the estimated gradient
                at XC.  (Output)

    Remark:
       This is Algorithm A5.6.3, Dennis and Schnabel, 1983, imsl_page 322.

    Keywords:   Forward difference; Gradient; Service routine

    GAMS:       G4f

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
void imsl_fdgrd(Mfloat (*fcn) (Mint, Mfloat*, Mfloat*), Mint n, Mfloat xc[],
                Mfloat xscale[], Mfloat *fc, Mfloat *epsfcn, Mfloat gc[])
#else
void imsl_fdgrd(Mfloat (*fcn) (Mint, Mfloat[], Mfloat*), Mint n, Mfloat xc[],
                Mfloat xscale[], Mfloat *fc, Mfloat *epsfcn, Mfloat gc[])
#endif
#else
void imsl_fdgrd(fcn, n, xc, xscale, fc, epsfcn, gc)
	Mfloat           (*fcn) ();
	Mint             n;
	Mfloat           xc[], xscale[], *fc, *epsfcn, gc[];
#endif
{
	Mint             i, j;
	Mfloat           eps, fnew, stepsz, xtempj;


	imsl_e1psh("FDGRD ");

	if (n <= 0) {
		imsl_e1sti(1, n);

/*		imsl_ermes(5, 1, "The number of variables must be positive while N = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_N_MUST_BE_POSITIVE);
	} else if (*epsfcn > 0.1e0 || *epsfcn < F_ZERO) {
		imsl_e1str(1, *epsfcn);

/*		imsl_ermes(5, 2, "The estimate for the relative noise in the function must be between 0.0 and 0.1 while EPSFCN = %(r1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_EPSFCN_VALUE);
	} else {
		for (i = 1; i <= n; i++) {
			if (xscale[i - 1] <= F_ZERO) {
				imsl_e1sti(1, i);
				imsl_e1str(1, xscale[i - 1]);

/*				imsl_ermes(5, 3, "The values for the diagonal scaling matrix must be positive while XSCALE(%(i1)) = %(r1) is given.");
*/
                                imsl_ermes(IMSL_TERMINAL,
				IMSL_POS_XSCALE_ELMNTS_NEEDED);
				goto L_9000;
			}
		}
	}

	if (imsl_n1rcd(0) == 0) {
		eps = sqrt(imsl_f_max(*epsfcn, imsl_amach(4)));
		for (j = 1; j <= n; j++) {
			stepsz = (eps) * imsl_f_max(fabs(xc[j - 1]), 
                                                    F_ONE / xscale[j - 1]);
			if (xc[j - 1] < F_ZERO)
				stepsz = -stepsz;
			xtempj = xc[j - 1];
			xc[j - 1] = xtempj + stepsz;
			(*fcn) (n, xc, &fnew);
			xc[j - 1] = xtempj;
			gc[j - 1] = (fnew - *fc) / stepsz;
		}
	}
L_9000:
	imsl_e1pop("FDGRD ");
	return;
}				/* end of function */
