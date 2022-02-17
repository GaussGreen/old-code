#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  B2NAK/DB2NAK (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 1, 1986

    Purpose:    Compute the 'not-a-knot' spline knot sequence.

    Usage:      CALL B2NAK (NDATA, XDATA, KORDER, XKNOT, XSRT, IWK)

    Arguments:
       NDATA  - Number of data points.  (Input)
       XDATA  - Array of length NDATA containing the location of the
                data points.  (Input)
       KORDER - Order of the spline.  (Input)
       XKNOT  - Array of length NDATA+KORDER containing the knot
                sequence.  (Output)
       XSRT   - Work array of length NDATA to hold the sorted XDATA
                values.
                If XDATA is not needed, XSRT may be the same as XDATA.
       IWK    - Work array of length NDATA to hold the permutation
                of XDATA.

    Remarks:
    1. Informational error
       Type Code
         4   4  The XDATA values must be distinct.

    2. The first knot is at the left endpoint and the last knot is
       slightly beyond the last endpoint.  Both endpoints have
       multiplicity KORDER.

    3. Interior knots have multiplicity one.

    Keyword:    B-spline interpolation

    GAMS:       E3

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b2nak(Mint *ndata, Mfloat xdata[], Mint *korder,
                Mfloat xknot[], Mfloat xsrt[],Mint iwk[])
#else
void imsl_b2nak(ndata, xdata, korder, xknot, xsrt, iwk)
	Mint            *ndata;
	Mfloat           xdata[];
	Mint            *korder;
	Mfloat           xknot[], xsrt[];
	Mint             iwk[];
#endif
{
	Mint             i, j;


	imsl_e1psh("IMSL_B2NAK ");
	/* Check argument KORDER */
	if (*korder <= 1) {
		imsl_e1sti(1, *korder);

/*		imsl_ermes(5, 1, "The order of the spline must be at least 2 while KORDER = %(i1) is given.");
*/
		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER_1);
	}
	/* Check NDATA */
	if (*ndata < *korder) {
		imsl_e1sti(1, *ndata);
		imsl_e1sti(2, *korder);

/*		imsl_ermes(5, 2, "The number of data points must be at least as large as the order of the spline while NDATA = %(i1) and KORDER = %(i2) are given.");
*/
		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_NEED_DATA_PTS);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Check argument XDATA */
	for (i = 2; i <= *ndata; i++) {
		if (xdata[i - 2] >= xdata[i - 1]) {
			/* Check that XDATA values are distinct */
			if (xdata[i - 2] == xdata[i - 1]) {
				j = i - 1;
				imsl_e1sti(1, j);
				imsl_e1sti(2, i-2);
				imsl_e1str(1, xdata[i - 1]);

/*				imsl_ermes(4, 4, "Points in the data point abscissas array, XDATA, must be distinct, but XDATA(%(i1)) = XDATA(%(i2)) = %(r1).");
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
	 * Data is already sorted.  Move XDATA to XSRT.
	 */
	scopy(*ndata, xdata, 1, xsrt, 1);
	goto L_50;
	/* Set initial permutation */
L_20:
	for (i = 1; i <= *ndata; i++) {
		iwk[i - 1] = i;
	}
	/* Find sorting permutation */
	imsl_svrgp(*ndata, xdata, xsrt, iwk);
	/* Check the XDATA values are distinct */
	for (i = 2; i <= *ndata; i++) {
		if (xsrt[i - 2] == xsrt[i - 1]) {
			imsl_e1sti(1, iwk[i - 2]-1);
			imsl_e1sti(2, iwk[i - 1]-1);
			imsl_e1str(1, xsrt[i - 1]);

/*			imsl_ermes(5, 5, "Points in the data point abscissas array, XDATA, must be distinct, but XDATA(%(i1)) = XDATA(%(i2)) = %(r1).");
*/
/*
       It is intenionally set to IMSL_FATAL not IMSL_TERMINAL.
*/
                        imsl_ermes(IMSL_FATAL,
			IMSL_DUPLICATE_XDATA_VALUES);
			goto L_9000;
		}
	}
L_50:
	;
	/*
	 * The first and last knots have multiplicity KORDER.  The last knot
	 * is extended to the right.
	 */
	for (i = 1; i <= *korder; i++) {
		xknot[i - 1] = xsrt[0];
		xknot[*ndata + i - 1] = xsrt[*ndata - 1] * (F_ONE + sign(F_ONE,
				   xsrt[*ndata - 1]) * imsl_amach(4) * 100);
		if (xsrt[*ndata - 1] == F_ZERO)
			xknot[*ndata + i - 1] = xsrt[*ndata - 1] + 100 * imsl_amach(4);
	}
	/* The middle knots have multiplicty one */
	if (mod(*korder, 2) == 0) {
		/* KORDER even - knots at XSRT points. */
		scopy(*ndata - *korder, &xsrt[*korder / 2], 1, &xknot[*korder],
			   1);
	} else {
		/*
		 * KORDER odd - knots between XSRT points.
		 */
		for (i = *korder + 1; i <= *ndata; i++) {
			xknot[i - 1] = F_HALF * (xsrt[i - *korder / 2 - 2] + xsrt[i - *korder / 2 - 1]);
		}
	}

L_9000:
	;
	imsl_e1pop("IMSL_B2NAK ");
	return;
}				/* end of function */
