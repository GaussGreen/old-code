#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  B2INT/DB2INT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 25, 1984

    Purpose:    Compute the spline interpolant, returning the B-spline
                representation.

    Usage:      CALL B2INT (NDATA, XDATA, FDATA, KORDER, XKNOT,
                            BSCOEF, WORK1, WORK2, WORK3, IWORK)

    Arguments:
       NDATA  - Number of data points.  (Input)
       XDATA  - Array of length NDATA containing the data point
                abscissas.  (Input)
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       KORDER - Order of the spline.  (Input)
                KORDER must be less than or equal to NDATA.
       XKNOT  - Array of length NDATA+KORDER containing the
                knot sequence.  (Input)
                XKNOT must be nondecreasing.
       BSCOEF - Array of length NDATA containing the B-spline
                coefficients.  (Output)
       WORK1  - Work array of length (5*KORDER-2)*NDATA.
       WORK2  - Work array of length NDATA.
       WORK3  - Work array of length NDATA.
       IWORK  - Work array of length NDATA.

    Remark:
       Informational errors
       Type Code
         3   1  The interpolation matrix is ill-conditioned.
         4   3  The XDATA values must be distinct.
         4   4  Multiplicity of the knots cannot exceed the order
                of the spline.
         4   5  The knots must be nondecreasing.
         4   15 The Ith smallest element of the data point array
                must be greater than the Ith knot and less than
                the (I+KORDER)th knot.
         4   16 The largest element of the data point array
                must be greater than the (NDATA)th knot and less
                than or equal to the (NDATA+KORDER)th knot.
         4   17 The smallest element of the data point array
                must be greater than or equal to the first knot
                and less than the (KORDER+1)st knot.

    GAMS:       E1a

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b2int(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mint *korder, Mfloat xknot[], Mfloat bscoef[],
	   Mfloat work1[], Mfloat work2[], Mfloat work3[], Mint iwork[])
#else
void /* FUNCTION */ 
imsl_b2int(ndata, xdata, fdata, korder, xknot, bscoef,
	   work1, work2, work3, iwork)
	Mint            *ndata;
	Mfloat           xdata[], fdata[];
	Mint            *korder;
	Mfloat           xknot[], bscoef[], work1[], work2[], work3[];
	Mint             iwork[];
#endif
{
	Mint             _l0, i, iaband, ifactr;


	imsl_e1psh("IMSL_B2INT");
	/* CHECK FOR ARGUMENTS FOR ERRORS */
	imsl_b3int(korder, xknot, ndata);
	if (imsl_n1rty(1) != 0)
		goto L_9000;
	/*
	 * SORT XDATA INTO WORK2 AND FDATA INTO WORK3. SET INITIAL
	 * PERMUTATION
	 */
	for (i = 1; i <= *ndata; i++) {
		iwork[i - 1] = i;
	}
	imsl_svrgp(*ndata, xdata, work2, iwork);
	/*
	 * CHECK XDATA FOR UNIQUENESS AND REARRANGE FDATA.
	 */
	work3[0] = fdata[iwork[0] - 1];
	for (i = 2; i <= *ndata; i++) {
		work3[i - 1] = fdata[iwork[i - 1] - 1];
		if (work2[i - 2] == work2[i - 1]) {
			imsl_e1sti(1, iwork[i - 2]-1);
			imsl_e1sti(2, iwork[i - 1]-1);
			imsl_e1str(1, work2[i - 1]);

			imsl_ermes(IMSL_FATAL, IMSL_DUPLICATE_XDATA_VALUES);
			goto L_9000;
		}
	}
	/*
	 * Make sure data points are in the correct interval.
	 */
	imsl_c1not("X", "KORDER", ndata, work2,
		   korder, xknot);
	if (imsl_n1rty(0) != 0)
		goto L_9000;

	iaband = 1;
	ifactr = iaband + (2 ** korder - 1) ** ndata;
        _l0 = 2 ** korder -1;
	imsl_b5int(ndata, work2, work3, korder, xknot, bscoef, &work1[iaband - 1],
		   &_l0, &work1[ifactr - 1], iwork);

L_9000:
	;
	imsl_e1pop("IMSL_B2INT");
	return;
}				/* end of function */
/*  -----------------------------------------------------------------------
    IMSL Name:  B3INT/DB3INT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 25, 1984

    Purpose:    Check B-spline parameters.

    Usage:      CALL B3INT (KORDER, XKNOT, NDATA)

    Arguments:
       KORDER - Order of the spline.  (Input)
       XKNOT  - Array of length NDATA+KORDER containing the knot
                sequence.  (Input)
                It must be nondecreasing.
       NDATA  - Length of BSCOEF.  (Input)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
void imsl_b3int(Mint *korder, Mfloat xknot[], Mint *ndata)
#else 
void imsl_b3int(korder, xknot, ndata)
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ndata;
#endif
{
	Mint             i, mult;

	/* CHECK KORDER */
	imsl_e1psh("IMSL_B3INT ");

	if (*korder < 1) {
		imsl_e1sti(1, *korder);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);
		goto L_9000;
	}
	/* CHECK NDATA */
	if (*ndata < *korder) {
		imsl_e1sti(1, *ndata);
		imsl_e1sti(2, *korder);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_NEED_DATA_PTS);
		goto L_9000;
	}
	/* CHECK KNOT SEQUENCE */
	mult = 1;
	for (i = 2; i <= (*ndata + *korder); i++) {
		if (xknot[i - 1] == xknot[i - 2]) {
			mult += 1;
			if (mult > *korder) {
				imsl_e1sti(1, (i-1) - mult + 1);
				imsl_e1sti(2, i - 1);
				imsl_e1str(1, xknot[i - 1]);
				imsl_e1sti(3, *korder);
                                imsl_e1stl(1,"X");

                                imsl_ermes(IMSL_FATAL,
				IMSL_KNOT_MULTIPLICITY);
				goto L_9000;
			}
		} else if (xknot[i - 1] < xknot[i - 2]) {
			imsl_e1sti(1, i - 2);
			imsl_e1sti(2, i - 1);
			imsl_e1str(1, xknot[i - 2]);
			imsl_e1str(2, xknot[i - 1]);
                        imsl_e1stl(1,"X");

			imsl_ermes(IMSL_FATAL, IMSL_KNOT_NOT_INCREASING);
			goto L_9000;
		} else {
			mult = 1;
		}
	}

L_9000:
	;
	imsl_e1pop("IMSL_B3INT ");
	return;
}				/* end of function */
/*  -----------------------------------------------------------------------
    IMSL Name:  B4INT/DB4INT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 6, 1986

    Purpose:    Generate the value of all B-splines with support at
                the given point.

    Usage:      CALL B4INT (XKNOT, JHIGH, X, LEFT, BIATX, DELTAL,
                            DELTAR)

    Arguments:
       XKNOT  - The knot sequence.  (Input)
       JHIGH  - The highest order B-spline to be computed.  (Input)
       X      - The point at which the B-splines are computed.  (Input)
       LEFT   - An integer such that (usually)
                XKNOT(LEFT) .LE. X .LE. XKNOT(LEFT+1)  (Input)
       BIATX  - Values of the B-splines at X.  (Output)
       DELTAL - Work array of length JHIGH.  (Input/Output)
       DELTAR - Work array of length JHIGH.  (Input/Output)

    Remark:
       Based on BSPLVB in deBoor, p 134.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
void imsl_b4int(Mfloat xknot[], Mint *jhigh, Mfloat *x, Mint *left, Mfloat biatx[], Mfloat deltal[], Mfloat deltar[])
#else
void imsl_b4int(xknot, jhigh, x, left, biatx, deltal, deltar)
	Mfloat           xknot[];
	Mint            *jhigh;
	Mfloat          *x;
	Mint            *left;
	Mfloat           biatx[], deltal[], deltar[];
#endif
{
	Mint             i, j;
	Mfloat           saved, term;


	biatx[0] = F_ONE;
	for (j = 1; j <= (*jhigh - 1); j++) {
		deltar[j - 1] = xknot[*left + j - 1] - *x;
		deltal[j - 1] = *x - xknot[*left - j];
		saved = F_ZERO;
		for (i = 1; i <= j; i++) {
			term = biatx[i - 1] / (deltar[i - 1] + deltal[j - i]);
			biatx[i - 1] = saved + deltar[i - 1] * term;
			saved = deltal[j - i] * term;
		}
		biatx[j] = saved;
	}

	return;
}				/* end of function */
/*  -----------------------------------------------------------------------
    IMSL Name:  B5INT/DB5INT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 6, 1986

    Purpose:    Compute the spline interpolant, returning the B-spline
                representation.

    Usage:      CALL B5INT (NDATA, XDATA, FDATA, KORDER, XKNOT, BSCOEF,
                            ABAND, LBAND, WORK, IPVT)

    Arguments:
       NDATA  - See BSINT.
       XDATA  - See BSINT.
       FDATA  - See BSINT.
       KORDER - See BSINT.
       XKNOT  - See BSINT.
       BSCOEF - See BSINT.
       ABAND  - 2*KORDER-1 by NDATA work array containing the
                interpolation matrix in band format.
       LBAND  - Bandwidth = 2*KORDER-1.  (Input)
       WORK   - 3*KORDER-1 by NDATA work array used by the linear
                equation solver as workspace.
       IPVT   - Work array of length NDATA containing pivot information.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b5int(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mint *korder, Mfloat xknot[], Mfloat bscoef[],
	   Mfloat *aband, Mint *lband, Mfloat work[], Mint ipvt[])
#else
void /* FUNCTION */ 
imsl_b5int(ndata, xdata, fdata, korder, xknot, bscoef,
	   aband, lband, work, ipvt)
	Mint            *ndata;
	Mfloat           xdata[], fdata[];
	Mint            *korder;
	Mfloat           xknot[], bscoef[], *aband;
	Mint            *lband;
	Mfloat           work[];
	Mint             ipvt[];
#endif
{
#define ABAND(I_,J_)	(aband+(I_)*(*lband)+(J_))
	Mint             _l0, _l1, _l2, i,left;
	Mfloat           xknoti;


	imsl_e1psh("IMSL_B5INT");



	/* Zero out all entries of ABAND */
	sset(*lband ** ndata, F_ZERO, aband, 1);
	/*
	 * Loop over I to construct the NDATA interpolation equations
	 */

	left = *korder;
	for (i = 1; i <= *ndata; i++) {
		xknoti = xdata[i - 1];
		/*
		 * Find LEFT in the closed interval (I,I+KORDER-1) such that
		 * XKNOT(LEFT) .LE. XDATA(I) .LT. XKNOT(LEFT+1). Matrix is
		 * singular if this is not possible.
		 */
		left = imsl_i_max(left, i);
L_10:
		if (xknoti >= xknot[left]) {
			left += 1;
			if (left < imsl_i_min(i + *korder, *ndata + 1))
				goto L_10;
			left -= 1;
		}
		/*
		 * The I-th equation enforces interpolation at XKNOTI, hence
		 * A(I,J) = B(J,KORDER,XKNOT)(XKNOTI), , all J. Only the
		 * KORDER entries with J = LEFT-KORDER+1, ..., LEFT actually
		 * might be nonzero. These KORDER numbers are returned by the
		 * following (in BSCOEF, for temp storage)
		 */
		imsl_b4int(xknot, korder, &xknoti, &left, bscoef, &work[0], &work[*korder]);
		/*
		 * We therefore want BSCOEF(J) = B(LEFT-K+J)(XKNOTI) to go
		 * into A(I,LEFT-KORDER+J) =
		 * ABAND(I-LEFT+2*KORDER-J,LEFT-KORDER +J) since the number
		 * of lower codiagonals is KORDER-1.
		 */
		scopy(*korder, bscoef, 1, ABAND( left-*korder,i - left + *korder * 2 - 2),
			   (*lband - 1));
	}
	/*
	 * Solve the linear system
	 */

        _l0 = *korder -1;
        _l1 = *korder -1;
        _l2 = 1;
	imsl_l2lrb(ndata, aband, lband, &_l0, &_l1,
		   fdata, &_l2, bscoef, &work[*ndata], ipvt, work);

	imsl_e1pop("IMSL_B5INT");
	return;
}				/* end of function */
