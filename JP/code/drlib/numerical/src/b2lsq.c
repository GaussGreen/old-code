#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  B2LSQ/DB2LSQ (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 14, 1986

    Purpose:    Compute the least squares spline approximation.

    Usage:      CALL B2LSQ (NDATA, XDATA, FDATA, WEIGHT, KORDER, XKNOT,
                            NCOEF, BSCOEF, WK, XSORT, FSORT, WSORT, IPVT)

    Arguments:
       NDATA  - Number of data points.  (Input)
       XDATA  - Array of length NDATA containing the data point
                abscissas.  (Input)
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       WEIGHT - Array of length NDATA containing the weights.  (Input)
       KORDER - Order of the spline.  (Input)
                KORDER must be less than or equal to NDATA.
       XKNOT  - Array of length NCOEF+KORDER containing the knot
                sequence.  (Input)
                XKNOT must be nondecreasing.
       NCOEF  - Number of B-spline coefficients.  (Input)
                NCOEF cannot be greater than NDATA.
       BSCOEF - Array of length NCOEF containing the B-spline
                coefficients.  (Output)
       WK     - Work array of length (3+NCOEF)*KORDER.
       XSORT  - Work array of length NDATA.
       FSORT  - Work array of length NDATA.
       WSORT  - Work array of length NDATA.
       IPVT   - Work array of length NDATA.

    Remark:
       Informational errors
       Type Code
         4   5  Multiplicity of the knots cannot exceed the order
                of the spline.
         4   6  The knots must be nondecreasing.
         4   7  All weights must be greater than zero.
         4   8  The smallest element of the data point array
                must be greater than or equal to the KORDth knot.
         4   9  The largest element of the data point array
                must be less than or equal to the (NCOEF+1)st knot.

    GAMS:       K1a1a1

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b2lsq(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat weight[], Mint *korder, Mfloat xknot[],
	   Mint *ncoef, Mfloat bscoef[], Mfloat wk[], Mfloat xsort[], Mfloat fsort[], Mfloat wsort[], Mint ipvt[])
#else
void imsl_b2lsq(ndata, xdata, fdata, weight, korder, xknot,
	   ncoef, bscoef, wk, xsort, fsort, wsort, ipvt)
	Mint            *ndata;
	Mfloat           xdata[], fdata[], weight[];
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           bscoef[], wk[], xsort[], fsort[], wsort[];
	Mint             ipvt[];
#endif
{
	Mint             i;


	imsl_e1psh("IMSL_B2LSQ ");
	/* CHECK FOR ARGUMENTS FOR ERRORS */
	imsl_b3lsq(ndata, korder, xknot, weight, ncoef);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/*
	 * SORT XDATA INTO XSORT AND FDATA INTO FSORT. SET INITIAL
	 * PERMUTATION
	 */
	for (i = 1; i <= *ndata; i++) {
		ipvt[i - 1] = i;
	}
	imsl_svrgp(*ndata, xdata, xsort, ipvt);
	/*
	 * CHECK XDATA FOR UNIQUENESS AND REARRANGE FDATA. THIS CHECK WAS
	 * REMOVED 11/89. UNIQUE XDATA IS NOT NECESSARY FOR THE ALGORITHM TO
	 * WORK.
	 */
	fsort[0] = fdata[ipvt[0] - 1];
	wsort[0] = weight[ipvt[0] - 1];
	for (i = 2; i <= *ndata; i++) {
		fsort[i - 1] = fdata[ipvt[i - 1] - 1];
		wsort[i - 1] = weight[ipvt[i - 1] - 1];
		/*
		 * ..... IF ($COMPUTER.EQ.CRAY2U .OR. $COMPUTER.EQ.CRXMPC
		 * .OR. &    $COMPUTER.EQ.CRXPMU .OR. $COMPUTER.EQ.CRAY1C)
		 * THEN TMP = XSORT(I-1) - XSORT(I) IF (TMP .EQ. 0.0) THEN
		 * .....    ELSE IF (XSORT(I-1) .EQ. XSORT(I)) THEN .....
		 * END IF CALL E1STI (1, IPVT(I-1)) CALL E1STI (2, IPVT(I))
		 * CALL E1STR (1, XSORT(I)) CALL E1MES (5, 4, 'Points in the
		 * interpolation vector '// &                 ', XDATA, must
		 * be distinct, but XDATA(%(I1)) '// &                 '=
		 * XDATA(%(I2)) = %(R1).') GO TO 9000 END IF
		 */
	}
	imsl_b4lsq(ndata, xsort, fsort, wsort, korder, xknot, ncoef, bscoef,
		   &wk[0], &wk[*korder ** ncoef]);

L_9000:
	;
	imsl_e1pop("IMSL_B2LSQ ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B3LSQ/DB3LSQ (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 14, 1986

    Purpose:    Check B-spline parameters.

    Usage:      CALL B3LSQ (NDATA, KORDER, XKNOT, WEIGHT, NCOEF)

    Arguments:
       NDATA  - Number of data points.  (Input)
       KORDER - Order of the spline.  (Input)
       XKNOT  - Array of length NCOEF+KORDER containing the knot
                sequence.  (Input)
                XKNOT must be nondecreasing.
       WEIGHT - Array of length NDATA containing the weights.  (Input)
       NCOEF  - Length of BSCOEF.  (Input)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b3lsq(Mint *ndata, Mint *korder, Mfloat xknot[], Mfloat weight[], Mint *ncoef)
#else
void imsl_b3lsq(ndata, korder, xknot, weight, ncoef)
	Mint            *ndata, *korder;
	Mfloat           xknot[], weight[];
	Mint            *ncoef;
#endif
{
	Mint             i, mult;
        Mint             num_weights_zero = 0;

	/* CHECK KORDER */
	if (*korder < 1) {
		imsl_e1sti(1, *korder);

/*		imsl_ermes(5, 1, "The order of the spline must be at least 1 while KORDER = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);
		goto L_9000;
	}
	/* CHECK NCOEF */
	if (*ncoef < *korder) {
		imsl_e1sti(1, *ncoef);
		imsl_e1sti(2, *korder);

/*		imsl_ermes(5, 2, "The number of coefficients must be at least as large as the order of the spline while NCOEF = %(i1) and KORDER = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_COEFFS);
		goto L_9000;
	}
	/* CHECK ARGUMENT NCOEF */
	if (*ncoef > *ndata) {
		imsl_e1sti(1, *ncoef);
		imsl_e1sti(2, *ndata);

/*		imsl_ermes(5, 3, "The number of coefficients must be less than or equal to the number of data points while NCOEF = %(i1) and NDATA = %(i2) are given.");
*/
		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_COEFFS_2);
		goto L_9000;
	}
	/* CHECK KNOT SEQUENCE */
	mult = 1;
	for (i = 2; i <= (*ncoef + *korder); i++) {
		if (xknot[i - 1] == xknot[i - 2]) {
			mult += 1;
			if (mult > *korder) {
				imsl_e1sti(1, (i-1) - mult + 1);
				imsl_e1sti(2, i-1);
				imsl_e1str(1, xknot[i - 1]);
				imsl_e1sti(3, *korder);
                                imsl_e1stl(1,"X");

/*				imsl_ermes(4, 5, "The knots XKNOT(%(i1)) through XKNOT(%(i2)) are all equal to %(r1).  The multiplicity of the knots must not exceed KORDER = %(i3).");
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
                        imsl_e1stl(1,"X");


/*			imsl_ermes(4, 6, "The knot XKNOT(%(i1)) = %(r1) and XKNOT(%(i2)) = %(r2).  The knots must be nondecreasing.");
*/
			imsl_ermes(IMSL_FATAL, IMSL_KNOT_NOT_INCREASING);
			goto L_9000;
		} else {
			mult = 1;
		}
	}
        /* CHECK WEIGHTS */
        for (i = 1; i <= *ndata; i++) {
                if (weight[i-1] == F_ZERO) num_weights_zero++;
                if (weight[i - 1] < F_ZERO) {
                        imsl_e1sti(1, i - 1);
                        imsl_e1str(1, weight[i - 1]);
                        imsl_e1stl(1,"X");

/*                        imsl_ermes(4, 7, "All elements of the argument WEIGHT must be greater than or equal to zero, but WEIGHT(%(i1)) = %(r1).");
*/
                        imsl_e1stl(1, "X");
			imsl_ermes(IMSL_FATAL, IMSL_NEGATIVE_WEIGHTS);
                        goto L_9000;
                }
        }
        if (num_weights_zero == *ndata){
                          imsl_e1stl(1,"X");


/*                        imsl_ermes(5, 10, "At least one element of WEIGHT must be greater than zero, but WEIGHT(I) = 0.0 for I =0,...,NDATA-1. ");
*/
			imsl_ermes(IMSL_FATAL, IMSL_SPLINE_NO_POS_ELMNT);                              
                        goto L_9000;
                }
            

L_9000:
	;
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B4LSQ/DB4LSQ (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 14, 1986

    Purpose:    Compute the least squares spline approximation.

    Usage:      CALL B4LSQ (NDATA, XDATA, FDATA, WEIGHT, KORDER, XKNOT,
                            NCOEF, BSCOEF, Q, WK)

    Arguments:  (See BSLSQ)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b4lsq(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat weight[], Mint *korder, Mfloat xknot[],
	   Mint *ncoef, Mfloat bscoef[], Mfloat *q, Mfloat wk[])
#else
void imsl_b4lsq(ndata, xdata, fdata, weight, korder, xknot,
	   ncoef, bscoef, q, wk)
	Mint            *ndata;
	Mfloat           xdata[], fdata[], weight[];
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           bscoef[], *q, wk[];
#endif
{
#define Q(I_,J_)	(q+(I_)*(*korder)+(J_))
	Mint             i, left, mm;
	Mfloat           dw;


	imsl_e1psh("IMSL_B4LSQ ");
	/* Test XDATA(1) .GE. XKNOT(KORDER) */
	if (xdata[0] < xknot[*korder - 1]) {
		imsl_e1str(1, xdata[0]);
		imsl_e1str(2, xknot[*korder - 1]);

/*		imsl_ermes(4, 8, "The smallest element in XDATA must be greater than or equal to XKNOT(KORDER) while the smallest element in XDATA = %(r1) and XKNOT(KORDER) = %(r2).");
*/
		imsl_ermes(IMSL_FATAL, IMSL_XDATA_TOO_SMALL);
		goto L_9000;
	}
	/* Test XDATA(NDATA) .LE. XKNOT(NDATA+1) */
	if (xdata[*ndata - 1] > xknot[*ncoef]) {
		imsl_e1str(1, xdata[*ndata - 1]);
		imsl_e1str(2, xknot[*ncoef]);

/*		imsl_ermes(4, 9, "The largest element in XDATA must be less than or equal to XKNOT(NCOEF+1) while the largest element in XDATA = %(r1) and XKNOT(NCOEF+1) = %(r2).");
*/
		imsl_ermes(IMSL_FATAL, IMSL_XDATA_TOO_LARGE);
		goto L_9000;
	}
	/* ZERO OUT ALL ENTRIES OF Q AND BSCOEF */
	sset(*korder ** ncoef, F_ZERO, q, 1);
	sset(*ncoef, F_ZERO, bscoef, 1);

	left = *korder;
	for (i = 1; i <= *ndata; i++) {
		/*
		 * FIND LEFT IN THE CLOSED INTERVAL (I,I+KORDER-1) SUCH THAT
		 * XKNOT(LEFT) .LE. XDATA(I) .LT. XKNOT(LEFT+1)
		 */
L_10:
		if (left < *ncoef && xdata[i - 1] >= xknot[left]) {
			left += 1;
			goto L_10;
		}
		/*
		 * FIND WK(MM) = B(LEFT-KORDER+MM) (XDATA(I))
		 */
		imsl_b4int(xknot, korder, &xdata[i - 1], &left, &wk[0], &wk[*korder],
			   &wk[*korder * 2]);
		for (mm = 1; mm <= *korder; mm++) {
			dw = wk[mm - 1] * weight[i - 1];
			bscoef[left - *korder + mm - 1] += dw * fdata[i - 1];
			saxpy(*korder - mm + 1, dw, &wk[mm - 1], 1, Q(left - *korder + mm - 1, 0),
				   1);
		}
	}
	/*
	 * CONSTRUCT THE CHOLESKY FACTORIZATION FOR C IN Q AND USE IT TO
	 * SOLVE THE NORMAL EQUATIONS C*X = BSCOEF FOR X AND STORE X IN
	 * BSCOEF
	 * 
	 * 1. FACTOR
	 */
	imsl_b5lsq(q, korder, ncoef);
	/* 2. SOLVE */
	imsl_b6lsq(q, korder, ncoef, bscoef);

L_9000:
	;
	imsl_e1pop("IMSL_B4LSQ ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B5LSQ/DB5LSQ (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 14, 1986

    Purpose:    Compute the least squares spline approximation.

    Usage:      CALL B5LSQ (W, NBANDS, NROW)

    Arguments:
       W      - Normal equation matrix in band symmetric form.
       NBANDS - Number of bands in W.  (Input)
       NROW   - Number of rows of W.  (Input)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b5lsq(Mfloat *w, Mint *nbands, Mint *nrow)
#else
void imsl_b5lsq(w, nbands, nrow)
	Mfloat          *w;
	Mint            *nbands, *nrow;
#endif
{
#define W(I_,J_)	(w+(I_)*(*nbands)+(J_))
	Mint             i, imax, j, jmax, n;
	Mfloat           ratio;


	if (*nrow > 1) {
		/*
		 * STORE DIAGONAL OF C IN DIAG CALL SCOPY (NROW, W, NBANDS,
		 * DIAG, 1) FACTORIZATION
		 */
		for (n = 1; n <= *nrow; n++) {
			if (*W(n - 1, 0) <= F_ONE / imsl_amach(2)) {
				/* SET TO ZERO IF SMALL */
				sset(*nbands, F_ZERO, W(n - 1, 0), 1);
			} else {
				*W(n - 1, 0) = F_ONE / *W(n - 1, 0);
				imax = imsl_i_min(*nbands - 1, *nrow - n);
				if (imax >= 1) {
					jmax = imax;
					for (i = 1; i <= imax; i++) {
						ratio = *W(n - 1, i) ** W(n - 1, 0);
						for (j = 1; j <= jmax; j++) {
							*W(n + i - 1, j - 1) += -*W(n - 1, j + i - 1) *
								ratio;
						}
						jmax -= 1;
						*W(n - 1, i) = ratio;
					}
				}
			}
		}
	} else {
		if (*W(0, 0) > F_ZERO)
			*W(0, 0) = F_ONE / *W(0, 0);
	}

	return;
}				/* end of function */
#undef W
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B6LSQ/DB6LSQ (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 14, 1986

    Purpose:    Compute the least squares spline approximation.

    Usage:      CALL B6LSQ (W, NBANDS, NROW, B)

    Arguments:
       W      - Normal equation matrix in band symmetric form.
       NBANDS - Number of bands in W.  (Input)
       NROW   - Number of rows of W.  (Input)
       B      - Array containing the RHS on input and the solution on
                output.  (Input/Output)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b6lsq(Mfloat *w, Mint *nbands, Mint *nrow, Mfloat b[])
#else
void imsl_b6lsq(w, nbands, nrow, b)
	Mfloat          *w;
	Mint            *nbands, *nrow;
	Mfloat           b[];
#endif
{
#define W(I_,J_)	(w+(I_)*(*nbands)+(J_))
	Mint             j, jmax, n, nbndm1;


	if (*nrow > 1) {
		/*
		 * FORWARD SUBSTITUTION. SOLVE L*Y=B FOR Y, STORE IN B.
		 */
		nbndm1 = *nbands - 1;
		for (n = 1; n <= *nrow; n++) {
			jmax = imsl_i_min(nbndm1, *nrow - n);
			if (jmax >= 1) {
				for (j = 1; j <= jmax; j++) {
					b[j + n - 1] += -*W(n - 1, j) * b[n - 1];
				}
			}
		}
		/*
		 * BACK SUBSTITUTION. SOLVE SOLVE TRANS(L)*X = D**(-1)*Y FOR
		 * X, STORE IN B
		 */
		for (n = *nrow; n >= 1; n--) {
			b[n - 1] *= *W(n - 1, 0);
			jmax = imsl_i_min(nbndm1, *nrow - n);
			if (jmax >= 1) {
				for (j = 1; j <= jmax; j++) {
					b[n - 1] += -*W(n - 1, j) * b[j + n - 1];
				}
			}
		}
	} else {
		b[0] *= *W(0, 0);
	}

	return;
}				/* end of function */
#undef W
