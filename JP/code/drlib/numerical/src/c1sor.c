#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif


/* Structured by FOR_STRUCT, v0.2, on 07/09/90 at 13:13:24
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  C1SOR/DC1SOR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 5, 1986

    Purpose:    Copy and sort XDATA into BREAK and FDATA into the first
                row of CSCOEF.

    Usage:      CALL C1SOR (NDATA, XDATA, FDATA, BREAK, CSCOEF, LDC,
                            IWORK)

    Arguments:
       NDATA  - Number of data points (must be at least 2).  (Input)
       XDATA  - Array of length NDATA containing the data point
                abscissas.  (Input)
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       BREAK  - Array of length NDATA containing the breakpoints
                for the piecewise cubic representation.  (Output)
       CSCOEF - Matrix of size 4 by NDATA containing the cubic spline
                (piecewise polynomial) representation in the first
                NDATA-1 columns.  (Output)
       LDC    - Leading dimension of CSCOEF exactly as specified in the
                dimension statement of the calling program.  (Input)
       IWORK  - Work array containing the permutation vector.  (Output)

    GAMS:       E1a

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c1sor(Mint ndata, Mfloat xdata[], Mfloat fdata[], Mfloat break_[],
                Mfloat *cscoef, Mint ldc, Mint iwork[])
#else
void imsl_c1sor(ndata, xdata, fdata, break_, cscoef, ldc, iwork)
	Mint             ndata;
	Mfloat           xdata[], fdata[], break_[], *cscoef;
	Mint             ldc, iwork[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(ldc)+(J_))
	Mint             i, j;

	/* CHECK FOR SORTED ARRAY */
	imsl_e1psh("IMSL_C1SOR");
	for (i = 2; i <= ndata; i++) {
		if (xdata[i - 2] >= xdata[i - 1]) {
			/* CHECK THAT XDATA VALUES ARE DISTINCT */
			if (xdata[i - 2] == xdata[i - 1]) {
				j = i - 1;
				imsl_e1sti(1, j-1);
				imsl_e1sti(2, i-1);
				imsl_e1str(1, xdata[i - 1]);

/*				"Points in the data point abscissas array, XDATA, must be distinct, but XDATA(%(i1)) = XDATA(%(i2)) = %(r1).");
*/
                                imsl_ermes(IMSL_FATAL,
				IMSL_DUPLICATE_XDATA_VALUES);

				goto L_9000;
			} else {
				goto L_20;
			}
		}
	}
	/*
	 * DATA IS ALREADY SORTED.  M1VE TO BREAK AND CSCOEF.
	 */
    scopy(ndata, xdata, 1, break_, 1);
	scopy(ndata, fdata, 1, cscoef, ldc);
	goto L_9000;
	/* SET INITIAL PERMUTATION */
L_20:
	for (i = 1; i <= ndata; i++) {
		iwork[i - 1] = i;
	}
	/* FIND SORTING PERMUTATION */
	imsl_svrgp(ndata, xdata, break_, iwork);
	/* APPLY PERMUTATION */
	for (i = 1; i <= ndata; i++) {
		*CSCOEF(i - 1, 0) = fdata[iwork[i - 1] - 1];
	}
	/* CHECK THE XDATA VALUES ARE DISTINCT */
	for (i = 2; i <= ndata; i++) {
		if (break_[i - 2] == break_[i - 1]) {
			imsl_e1sti(1, iwork[i - 2] -1);
			imsl_e1sti(2, iwork[i - 1] -1);
			imsl_e1str(1, break_[i - 1]);

/*			"Points in the data point abscissas array, XDATA, must be distinct, but XDATA(%(i1)) = XDATA(%(i2)) = %(r1).");
*/
                                imsl_ermes(IMSL_FATAL,
				IMSL_DUPLICATE_XDATA_VALUES);
			goto L_9000;
		}
	}

L_9000:
	;
	imsl_e1pop("IMSL_C1SOR");
	return;
}				/* end of function */


