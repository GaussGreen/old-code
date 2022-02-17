#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  C1NOT/DC1NOT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 11, 1986

    Purpose:    Check knot and data point sequences for B-spline
                interpolation.

    Usage:      CALL C1NOT (VNAME, KNAME, NDATA, XDATA, KORDER, XKNOT)

    Arguments:
       VNAME  - Character string containing name of dimension being
                checked.  (Input)
       KNAME  - Character string containing name the order being
                checked.  (Input)
       NDATA  - Number of data points.  (Input)
       XDATA  - Data point arrary of length NDATA.  (Input)
       KORDER - Order of the spline.  (Input)
       XKNOT  - Array of length NDATA+KORDER containing the knot
                sequence.  (Input)
                It must be nondecreasing.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c1not(Mchar *vname, Mchar *kname, Mint *ndata, Mfloat xdata[],
	   Mint *korder, Mfloat xknot[])
#else
void imsl_c1not(vname, kname, ndata, xdata,
	   korder, xknot)
	Mchar           *vname;
	Mchar           *kname;
	Mint            *ndata;
	Mfloat           xdata[];
	Mint            *korder;
	Mfloat           xknot[];
#endif
{
	Mchar            nname[7];
	Mint             i, mult;
	Mfloat           var1, var2, var3;


	if (strcmp(kname, "KORDER") == 0) {
		strcpy(nname, "NDATA");
	} else {
       strcpy(nname,"N");
       strncat(nname,vname,1);
       strncat(nname,"DATA",4);
	}

	/* Check KORDER */
	if (*korder <= 0) {
		imsl_e1sti(1, *korder);
		imsl_e1stl(1, vname);
		imsl_e1stl(2, kname);

/*		imsl_ermes(5, 10, "The order of the spline in the %(l1) direction must be at least 1 while %(l2) = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_ARB_2);
		goto L_9000;
	}
	/* Check argument NDATA */
	if (*ndata < *korder) {
		imsl_e1sti(1, *ndata);
		imsl_e1sti(2, *korder);
		imsl_e1stl(1, vname);
		imsl_e1stl(2, kname);
		imsl_e1stl(3, nname);

/*		imsl_ermes(5, 11, "The number of data points in the %(l1) direction must be at least as large as the order of the spline in the %(l1) direction while %(l3) = %(i1) and %(l2) = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_DATA_PTS_XY);
		goto L_9000;
	}
	/* Check knot sequence */
	mult = 1;
	for (i = 2; i <= (*ndata + *korder); i++) {
		if (xknot[i - 1] == xknot[i - 2]) {
			mult += 1;
			if (mult > *korder) {
				imsl_e1sti(1, (i-1) - mult + 1);
				imsl_e1sti(2, i-1);
				imsl_e1str(1, xknot[i - 1]);
				imsl_e1sti(3, *korder);
				imsl_e1stl(1, vname);
				imsl_e1stl(2, kname);

/*				imsl_ermes(4, 13, "The knots %(l1)KNOT(%(i1)) through %(l1)KNOT(%(i2)) are all equal to %(r1).  The multiplicity of the knots must not exceed %(l2) = %(i3).");
*/
                                imsl_ermes(IMSL_FATAL,
				IMSL_KNOT_MULTIPLICITY);
				goto L_9000;
			}
		} else if (xknot[i - 1] < xknot[i - 2]) {
			imsl_e1sti(1, i - 2);
			imsl_e1sti(2, i - 1);
			imsl_e1str(1, xknot[i - 2]);
			imsl_e1str(2, xknot[i - 1]);
			imsl_e1stl(1, vname);

/*			imsl_ermes(4, 14, "The knot %(l1)KNOT(%(i1)) = %(r1) and %(l1)KNOT(%(i2)) = %(r2).  The knots must be nondecreasing.");
*/
                        imsl_ermes(IMSL_FATAL, IMSL_KNOT_NOT_INCREASING);
			goto L_9000;
		} else {
			mult = 1;
		}
	}
	/* Check data points */
	if (xdata[0] < xknot[*korder] && xdata[0] >= xknot[0]) {
		for (i = 2; i <= (*ndata - 1); i++) {
			if (!(xdata[i - 1] < xknot[i + *korder - 1] && xdata[i - 1] >=
			      xknot[i - 1])) {
				var1 = xdata[i - 1];
				var2 = xknot[i - 1];
				var3 = xknot[i + *korder - 1];
				imsl_e1stl(1, vname);
				imsl_e1stl(2, kname);
				imsl_e1sti(1, i -1);
				imsl_e1sti(2, *korder + i -1);
				imsl_e1str(1, var1);
				imsl_e1str(2, var2);
				imsl_e1str(3, var3);
/*             "The Ith smallest element (I=%(i1)) in */
/* %(l1)DATA must be in the interval (%(l1)KNOT(I),%(l1)KNOT(I+%(l2))) */
/* while %(l1)DATA(%(i1)) = %(r1), %(l1)KNOT(%(i1)) = %(r2), %(l1)KNOT(%(i2))*/
/* = %(r3) are given.  This causes the interpolation matrix to be */
/*singular.");*/
/*				imsl_ermes(4, 15, "The Ith s.el. (I=%(i1)) in %(l1)DATA m.b.i.t.i (%(l1)KNOT(I),%(l1)KNOT(I+%(l2))) w. %(l1)DATA(%(i1)) = %(r1), %(l1)KNOT(%(i1)) = %(r2), %(l1)KNOT(%(i2)) = %(r3) a.g.sing.SEE CODE..");
*/
                                imsl_ermes(IMSL_FATAL,
				IMSL_KNOT_DATA_INTERLACING);
				goto L_9000;
			}
		}
		/* Check the largest data point. */
		if ((xdata[*ndata - 1] > xknot[*ndata + *korder - 1]) || (xdata[*ndata - 1] <=
						       xknot[*ndata - 1])) {
			var1 = xdata[*ndata - 1];
			var2 = xknot[*ndata - 1];
			var3 = xknot[*ndata + *korder - 1];
			imsl_e1stl(1, vname);
			imsl_e1stl(2, kname);
			imsl_e1str(1, var1);
			imsl_e1str(2, var2);
			imsl_e1str(3, var3);

/*			imsl_ermes(4, 16, "The largest element in %(l1)DATA must be greater than %(l1)KNOT(NDATA) and less than or equal to %(l1)KNOT(NDATA+%(l2)) while %(l1)DATA(NDATA) = %(r1), %(l1)KNOT(NDATA) = %(r2), %(l1)KNOT(%(i2)) = %(r3) are given.  ");
*/
                        imsl_ermes(IMSL_FATAL, IMSL_DATA_TOO_LARGE);
			goto L_9000;
		}
	} else {
		var1 = xdata[0];
		var2 = xknot[0];
		var3 = xknot[*korder];
		imsl_e1stl(1, vname);
		imsl_e1stl(2, kname);
		imsl_e1sti(2, *korder + 1);
		imsl_e1str(1, var1);
		imsl_e1str(2, var2);
		imsl_e1str(3, var3);
/*		"The smallest element in %(l1)DATA must be 
greater than or equal to %(l1)KNOT(1) and less than %(l1)KNOT(1+%(l2)) while
 %(l1)DATA(1) = %(r1), %(l1)KNOT(1) = %(r2), %(l1)KNOT(%(i2)) = %(r3) are 
 given.  This causes the interpolation matrix to be singular.");
*/
/*                imsl_ermes(4, 17, "ERROR - see code!!!.  This causes the interpolation matrix to be singular.");*/
                imsl_ermes(IMSL_FATAL, IMSL_DATA_TOO_SMALL);

	}

L_9000:
	;
	return;
}				/* end of function */
