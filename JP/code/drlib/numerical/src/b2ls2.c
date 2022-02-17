#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef TRUE
#undef TRUE
#define TRUE 1
#else
#define TRUE 1
#endif
#ifdef FALSE
#undef FALSE
#define FALSE 0
#else
#define FALSE 0
#endif

/* Structured by FOR_STRUCT, v0.2, on 08/23/90 at 11:53:24
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  B2LS2/DB2LS2 (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 17, 1986

    Purpose:    Compute the B-spline coefficients of a tensor product
                spline by solving a least squares problem with tensor
                product data.

    Usage:      CALL B2LS2 (NXDATA, XDATA, NYDATA, YDATA, FDATA, LDF,
                            KXORD, KYORD, XKNOT, YKNOT, NXCOEF, NYCOEF,
                            XWEIGH, YWEIGH, BSCOEF, WK)

    Arguments:
       NXDATA - Number of data points in the X-direction.  (Input)
       XDATA  - Array of length NXDATA containing the data points in the
                X-direction.  (Input)
                XDATA must be nondecreasing.
       NYDATA - Number of data points in the Y-direction.  (Input)
       YDATA  - Array of length NYDATA containing the data points in the
                Y-direction.  (Input)
                YDATA must be nondecreasing.
       FDATA  - Array of size NXDATA by NYDATA containing the values on
                the X-Y grid to be interpolated.  (Input)
                FDATA(I,J) contains the value at (XDATA(I),YDATA(I)).
       LDF    - Leading dimension of FDATA exactly as specified in
                the dimension statement of calling program.  (Input)
       KXORD  - Order of the spline in the X-direction.  (Input)
       KYORD  - Order of the spline in the Y-direction.  (Input)
       XKNOT  - Array of length KXORD+NXCOEF containing the knots in the
                X-direction.  (Input)
                XKNOT must be nondecreasing.
       YKNOT  - Array of length KYORD+NYCOEF containing the knots in the
                Y-direction.  (Input)
                YKNOT must be nondecreasing.
       NXCOEF - Number of B-spline coefficients in the X-direction.
                (Input)
       NYCOEF - Number of B-spline coefficients in the Y-direction.
                (Input)
       XWEIGH - Array of length NXDATA containing the positive weights
                of XDATA.  (Input)
       YWEIGH - Array of length NYDATA containing the positive weights
                of YDATA.  (Input)
       BSCOEF - Array of length NXCOEF*NYCOEF that contains the
                tensor product B-spline coefficients.  (Output)
                BSCOEF is treated internally as an array of size
                NXCOEF by NYCOEF.
       WK     - Work array of length (NXCOEF+1)*NYDATA+KXORD*NXCOEF+
                KYORD*NYCOEF+3*MAX(KXORD,KYORD).

    Remark:
       Informational errors
       Type Code
         3   14 There may be less than one digit of accuracy in the
                least squares fit.  Try using higher precision if
                possible.
         4   5  Multiplicity of the knots cannot exceed the order
                of the spline.
         4   6  The knots must be nondecreasing.
         4   7  All weights must be greater than zero.
         4   9  The data point values must be strictly increasing.
         4   10 The smallest element of the data point array
                must be greater than or equal to the K_ORDth knot.
         4   11 The largest element of the data point array
                must be less than or equal to the (N_COEF+1)st knot.

    GAMS:       K1a1b

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
 void imsl_b2ls2(Mint *nxdata, Mfloat xdata[], Mint *nydata, Mfloat ydata[],
                    Mfloat fdata[], Mint *ldf, Mint *kxord, Mint *kyord, 
                    Mfloat xknot[], Mfloat yknot[], Mint *nxcoef, Mint *nycoef,
                    Mfloat xweigh[], Mfloat yweigh[], Mfloat bscoef[], Mfloat wk[])
#else
 void imsl_b2ls2(nxdata, xdata, nydata, ydata, fdata, ldf,
	 kxord, kyord, xknot, yknot, nxcoef, nycoef, xweigh, yweigh, bscoef,
	   wk)
	Mint            *nxdata;
	Mfloat           xdata[];
	Mint            *nydata;
	Mfloat           ydata[], fdata[];
	Mint            *ldf, *kxord, *kyord;
	Mfloat           xknot[], yknot[];
	Mint            *nxcoef, *nycoef;
	Mfloat           xweigh[], yweigh[], bscoef[], wk[];
#endif
{
 	Mint             _l0, ia, ic, iwbcof, iwfdta, iwk;

	imsl_e1psh("imsl_b2ls2");
	/* Check LDF */
	if (*ldf < *nxdata) {
		imsl_e1sti(1, *nxdata);
		imsl_e1sti(2, *ldf);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_LD_FDATA_1);
		goto L_9000;
	}
	/* Check X data */
        _l0 = FALSE;
	imsl_b3ls2(nxdata, xdata, kxord, xknot, nxcoef, xweigh, "X", sizeof("X"
							), (Mint *)&_l0);

	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Check Y data (scalars only) */
        _l0 = FALSE;
	imsl_b3ls2(nydata, ydata, kyord, yknot, nycoef, yweigh, "Y", sizeof("Y"
							), (Mint *)&_l0);

	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/*
	 * Partition workspace WBCOEF(NXCOEF,NYDATA) WFDATA(NYDATA)
	 * A(KXORD,NXCOEF) C(KYORD,NYCOEF) WK(3*KORD) where
	 * KORD=MAX(KXORD,KYORD)
	 */
	iwbcof = 1;
	iwfdta = iwbcof + *nxcoef ** nydata;
	ia = iwfdta + *nydata;
	ic = ia + *kxord ** nxcoef;
	iwk = ic + *kyord ** nycoef;
	/* Perform least-squares approximation */
	imsl_b4ls2(nxdata, xdata, nydata, ydata, fdata, ldf, kxord, kyord,
		   xknot, yknot, nxcoef, nycoef, xweigh, yweigh, bscoef, &wk[iwbcof - 1],
		   &wk[iwfdta - 1], &wk[ia - 1], &wk[ic - 1], &wk[iwk - 1]);

L_9000:
	imsl_e1pop("imsl_b2ls2");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B3LS2/DB3LS2 (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 17, 1986

    Purpose:    Perform error checking for 2-dimensional tensor product
                B-spline.

    Usage:      CALL B3LS2 (NXDATA, XDATA, KXORD, XKNOT, NXCOEF, XWEIGH,
                            XY, SCALAR)

    Arguments:
       NXDATA - Number of data points in the X-direction.  (Input)
       XDATA  - Data points in the X-direction.  (Input)
                The points must be distinct and in ascending order.
       KXORD  - Order of the spline in the X-direction.  (Input)
       XKNOT  - Vector of length KXORD+NXCOEF containing the knots in the
                X-direction.  (Input)
                The knots must be increasing.
       NXCOEF - Number of B-spline coefficients in the X-direction.
                (Input)
       XWEIGH - Vector of length NXDATA containing the positive weights
                of XDATA.  (Input)
       XY     - Character describing the data.  (Input)
                XY should be either 'X' or 'Y'.
       SCALAR - Logical variable telling if only scalars are to be
                checked.  (Input)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* XY_S is not used, but leave the calling sequence intact.*/
#ifdef ANSI
 void imsl_b3ls2(Mint *nxdata, Mfloat xdata[], Mint *kxord, Mfloat xknot[],
                    Mint *nxcoef, Mfloat xweigh[], Mchar *xy, Mint xy_s,
                    Mint *scalar)
#else
 void imsl_b3ls2(nxdata, xdata, kxord, xknot, nxcoef, xweigh,
       	            xy, xy_s, scalar)
	Mint            *nxdata;
	Mfloat           xdata[];
	Mint            *kxord;
	Mfloat           xknot[];
	Mint            *nxcoef;
	Mfloat           xweigh[];
	Mchar           *xy;
	Mint             xy_s;
	Mint            *scalar;
#endif
{
	Mint             _l0, i, j, mult;
        Mint             num_xweigh_zero = 0;


	imsl_e1psh("imsl_b3ls2");
	/* Check NXDATA */
	if (*nxdata < 1) {
		imsl_e1stl(1, xy);
		imsl_e1sti(1, *nxdata);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_POS_DATA_PTS);
	}
	/* Check KXORD */
	if (*kxord < 1) {
		imsl_e1stl(1, xy);
		imsl_e1sti(1, *kxord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_POSI);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Check NXCOEF */
	if (*nxcoef < *kxord) {
		imsl_e1stl(1, xy);
		imsl_e1sti(1, *nxcoef);
		imsl_e1sti(2, *kxord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_COEFF_XY);
	}
	/* Check NXCOEF */
	if (*nxcoef > *nxdata) {
		imsl_e1stl(1, xy);
		imsl_e1sti(1, *nxcoef);
		imsl_e1sti(2, *nxdata);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_MORE_COEF_REQ);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;

	if (!*scalar) {
		/* Check knot sequence */
		mult = 1;
		for (i = 2; i <= (*nxcoef + *kxord); i++) {
			if (xknot[i - 1] == xknot[i - 2]) {
				mult += 1;
				if (mult > *kxord) {
					imsl_e1stl(1, xy);
					imsl_e1sti(1, (i-1) - mult + 1);
					imsl_e1sti(2, i - 1);
					imsl_e1str(1, xknot[i - 1]);
					imsl_e1sti(3, *kxord);

					imsl_ermes(IMSL_FATAL,
					IMSL_KNOT_MULTIPLICITY);
					goto L_9000;
				}
			} else if (xknot[i - 1] < xknot[i - 2]) {
				imsl_e1stl(1, xy);
				imsl_e1sti(1, i - 2);
				imsl_e1sti(2, i - 1);
				imsl_e1str(1, xknot[i - 2]);
				imsl_e1str(2, xknot[i - 1]);

				imsl_ermes(IMSL_FATAL,
				IMSL_KNOT_NOT_INCREASING);
			} else {
				mult = 1;
			}
		}
		/* CHECK WEIGHTS */
		for (i = 1; i <= *nxdata; i++) {
                        if (xweigh[i-1] == F_ZERO) num_xweigh_zero++;
			if (xweigh[i - 1] < F_ZERO) {
				imsl_e1stl(1, xy);
				imsl_e1sti(1, i - 1);
				imsl_e1str(1, xweigh[i - 1]);
				/* Print error message */
                                _l0 = 7;
				imsl_b7ls2(&_l0);
				goto L_9000;
			}
		}
                if (num_xweigh_zero == *nxdata){
		                imsl_e1stl(1, xy);
                                _l0 = 12;
				imsl_b7ls2(&_l0);
				goto L_9000;
		}
		for (i = 2; i <= *nxdata; i++) {
			if (xdata[i - 2] >= xdata[i - 1]) {
				/*
				 * Check that XDATA values are distinct This
				 * check was eliminated 11/89.  The routine
				 * should work with nondistinct data.  Note
				 * that XDATA must be nondecreasing.
				 */
				if (xdata[i - 2] == xdata[i - 1]) {
					/*
					 * J = I - 1 CALL E1STL (1, XY) CALL
					 * E1STI (1, J) CALL E1STI (2, I)
					 * CALL E1STR (1, XDATA(I)) Print
					 * error message CALL B7LS2 (8) C
					 * CALL E1MES (5, 8, 'Points in the
					 * %(L1) data point'// C     &
					 * ' array, %(L1)DATA, must be '// C
					 * &
					 * 'distinct, but %(L1)DATA(%(I1)) =
					 * '// C     &
					 * '%(L1)DATA(%(I2)) = %(R1).') GO TO
					 * 9000
					 */
				} else {
					j = i - 1;
					imsl_e1stl(1, xy);
					imsl_e1sti(1, j - 1);
					imsl_e1sti(2, i - 1);
					imsl_e1str(1, xdata[j - 1]);
					imsl_e1str(2, xdata[i - 1]);
					/* Print error message */
					_l0 = 9;
					imsl_b7ls2(&_l0);
					/*
					 * C                  CALL E1MES (5,
					 * 9, 'Points in the %(L1) data
					 * point'// C     &
					 * ' array, %(L1)DATA, must be '// C
					 * &
					 * 'nondecreasing, but
					 * %(L1)DATA(%(I1)) = ' C     &
					 * '%(R1) and
					 * %(L1)DATA(%(I2))=%(R2).')
					 */
					goto L_9000;
				}
			}
		}
		/* Test XDATA(1) .GE. XKNOT(KXORD) */
		if (xdata[0] < xknot[*kxord - 1]) {
			imsl_e1stl(1, xy);
			imsl_e1str(1, xdata[0]);
			imsl_e1str(2, xknot[*kxord - 1]);
			/* Print error message */
			_l0 = 10;
			imsl_b7ls2(&_l0);
			/*
			 * C            CALL E1MES (5, 10, 'The smallest
			 * element in %(L1)DATA '// C     &
			 * 'must greater than or equal to '// C     &
			 * '%(L1)KNOT(K%(L1)ORD) while the smallest '// C
			 * &                  'element in %(L1)DATA = %(R1)
			 * and '// C     &
			 * '%(L1)KNOT(K%(L1)ORD) = %(R2).')
			 */
			goto L_9000;
		}
		/* Test XDATA(NXDATA).LE.XKNOT(NCOEF+1) */
		if (xdata[*nxdata - 1] > xknot[*nxcoef]) {
			imsl_e1stl(1, xy);
			imsl_e1str(1, xdata[*nxdata - 1]);
			imsl_e1str(2, xknot[*nxcoef]);
			/* Print error message */
                        _l0 = 11;
			imsl_b7ls2(&_l0);
			goto L_9000;
		}
	}
L_9000:
	;
	imsl_e1pop("imsl_b3ls2");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B4LS2/DB4LS2 (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 17, 1986

    Purpose:    Calculate the B-spline coefficients of a tensor product
                spline by solving a least squares problem with tensor
                product data.

    Usage:      CALL B4LS2 (NXDATA, XDATA, NYDATA, YDATA, FDATA, LDF,
                            KXORD, KYORD, XKNOT, YKNOT, NXCOEF, NYCOEF,
                            XWEIGH, YWEIGH, BSCOEF, WBCOEF, WFDATA,
                            A, C, WK)

    Arguments:
       NXDATA - Number of data points in the X-direction.  (Input)
       XDATA  - Data points in the X-direction.  (Input)
                The points must be nondecreasing.
       NYDATA - Number of data points in the Y-direction.  (Input)
       YDATA  - Data points in the Y-direction.  (Input)
                The points must be nondecreasing.
       FDATA  - Values on the X-Y grid to be fit.  (Input)
                FDATA(I,J) contains the function value at
                (XDATA(I),YDATA(I)).
       LDF    - Leading dimension of FDATA exactly as specified in
                the calling program.  (Input)
       KXORD  - Order of the spline in the X-direction.  (Input)
       KYORD  - Order of the spline in the Y-direction.  (Input)
       XKNOT  - Vector of length KXORD+NXCOEF containing the knots in the
                X-direction.  (Input)
                The knots must be non-decreasing.
       YKNOT  - Vector of length KYORD+NYCOEF containing the knots in the
                Y-direction.  (Input)
                The knots must be non-decreasing.
       NXCOEF - Number of B-spline coefficients in the X-direction.
                (Input)
       NYCOEF - Number of B-spline coefficients in the Y-direction.
                (Input)
       XWEIGH - Vector of length NXDATA containing the positive weights
                of XDATA.  (Input)
       YWEIGH - Vector of length NYDATA containing the positive weights
                of YDATA.  (Input)
       BSCOEF - Matrix of size NXCOEF by NYCOEF which contains
                the coefficients of the tensor product spline.  (Output)
                BSCOEF is dimension NYCOEF by NXCOEF in this routine
                to allow BSCOEF to be used as workspace.
       WBCOEF - Work vector of size NXCOEF*NYDATA.
       WFDATA - Work vector of size NYDATA.
       A      - Work vector of size KXORD*NXCOEF.
       C      - Work vector of size KYORD*NYCOEF.
       WK     - Work vector of size 3*MAX(KXORD,KYORD).

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
 void imsl_b4ls2(Mint *nxdata, Mfloat xdata[], Mint *nydata, Mfloat ydata[],
                    Mfloat *fdata, Mint *ldf, Mint *kxord, Mint *kyord,
                    Mfloat xknot[], Mfloat yknot[], Mint *nxcoef, Mint *nycoef,
                    Mfloat xweigh[], Mfloat yweigh[], Mfloat *bscoef, Mfloat *wbcoef,
                    Mfloat wfdata[], Mfloat *a, Mfloat *c, Mfloat wk[])
#else
 void imsl_b4ls2(nxdata, xdata, nydata, ydata, fdata, ldf,
	 kxord, kyord, xknot, yknot, nxcoef, nycoef, xweigh, yweigh, bscoef,
	   wbcoef, wfdata, a, c, wk)
	Mint            *nxdata;
	Mfloat           xdata[];
	Mint            *nydata;
	Mfloat           ydata[], *fdata;
	Mint            *ldf, *kxord, *kyord;
	Mfloat           xknot[], yknot[];
	Mint            *nxcoef, *nycoef;
	Mfloat           xweigh[], yweigh[], *bscoef, *wbcoef, wfdata[],
	               *a, *c, wk[];
#endif
{
#define FDATA(I_,J_)	(fdata+(I_)*(*ldf)+(J_))
#define BSCOEF(I_,J_)	(bscoef+(I_)*(*nycoef)+(J_))
#define WBCOEF(I_,J_)	(wbcoef+(I_)*(*nxcoef)+(J_))
#define A(I_,J_)	(a+(I_)*(*kxord)+(J_))
#define C(I_,J_)	(c+(I_)*(*kyord)+(J_))
	Mint             first;
	Mint             _l0, _l1, i, nminus, nplus;
	Mfloat           ainnrm, anr1m, cinnrm, cnr1m, cnum;


	imsl_e1psh("imsl_b4ls2");
	/*
	 * Fill least squares matrix A on the first call to B5LS2 and then
	 * fill up WBCOEF with the answer to the RHS (RHS - FDATA) of least
	 * squares system.
	 */
	for (i = 1; i <= *nydata; i++) {
		first = FALSE;
		if (i == 1)
			first = TRUE;
		imsl_b5ls2(nxdata, xdata, FDATA(i - 1, 0), xweigh, kxord, xknot,
			   nxcoef, WBCOEF(i - 1, 0), a, wk, &first);
	}

	/*
	 * Fill least squares matrix C on the first call to B5LS2 and then
	 * fill up BSCOEF with RHS (RHS - BSCOEF of previous system) of least
	 * squares system.
	 */
	for (i = 1; i <= *nxcoef; i++) {
		/* Use a row of WBCOEF for FDATA */
		scopy(*nydata, WBCOEF(0, i - 1), *nxcoef, wfdata, 1);
		first = FALSE;
		if (i == 1)
			first = TRUE;
		imsl_b5ls2(nydata, ydata, wfdata, yweigh, kyord, yknot, nycoef,
			   BSCOEF(i - 1, 0), c, wk, &first);
	}
        _l0 = *kxord - 1;
        _l1 = 0;
	imsl_nr1rb(nxcoef, a, kxord, &_l0, &_l1, &anr1m);
        _l0 = *kyord - 1;
        _l1 = 0;
	imsl_nr1rb(nycoef, c, kyord, &_l0, &_l1, &cnr1m);
	/*
	 * Now solve CINV*TRAN(R) where TRAN(R) is in BSCOEF. Factor C
	 */
	imsl_b5lsq(c, kyord, nycoef);
	for (i = 1; i <= *nxcoef; i++) {
		/* Solve with multiple RHS */
		imsl_b6lsq(c, kyord, nycoef, BSCOEF(i - 1, 0));
		/* Place answer in WBCOEF */
		scopy(*nycoef, BSCOEF(i - 1, 0), 1, WBCOEF(0, i - 1), *nxcoef);
	}
	/*
	 * Now solve AINV*WBCOEF We must call B6LS2 to change the dimension
	 * of BSCOEF The results are in BSCOEF
	 */
	imsl_b6ls2(nxcoef, nycoef, kxord, a, wbcoef, bscoef);
	/* Calculate condition number */
	if (*nxcoef / 2 == (Mfloat) (*nxcoef) / F_TWO) {
		nplus = *nxcoef / 2;
		nminus = *nxcoef / 2;
	} else {
		nplus = *nxcoef / 2 + 1;
		nminus = *nxcoef / 2;
	}
	sset(nplus, F_ONE, wbcoef, 2);
	sset(nminus, -F_ONE, WBCOEF(0, 1), 2);
	imsl_b6lsq(a, kxord, nxcoef, wbcoef);
	i = imsl_isamax(*nxcoef, wbcoef, 1);
	ainnrm = fabs(*WBCOEF(0, i - 1));
	if (*nycoef / 2 == (Mfloat) (*nycoef) / F_TWO) {
		nplus = *nycoef / 2;
		nminus = *nycoef / 2;
	} else {
		nplus = *nycoef / 2 + 1;
		nminus = *nycoef / 2;
	}
	sset(nplus, F_ONE, wfdata, 2);
	sset(nminus, -F_ONE, &wfdata[1], 2);
	imsl_b6lsq(c, kyord, nycoef, wfdata);
	i = imsl_isamax(*nycoef, wfdata, 1);
	cinnrm = fabs(wfdata[i - 1]);
	cnum = ainnrm * anr1m * cnr1m * cinnrm;
	if (cnum > F_ONE / (F_TEN * imsl_amach(4))) {

		imsl_ermes(IMSL_WARNING, IMSL_SPLINE_LOW_ACCURACY);
	}
	imsl_e1pop("imsl_b4ls2");
	return;
}				/* end of function */
#undef FDATA
#undef BSCOEF
#undef WBCOEF
#undef A
#undef C
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B5LS2/DB5LS2 (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 17, 1986

    Purpose:    Compute the least squares spline approximation.

    Usage:      CALL B5LS2 (NDATA, XDATA, FDATA, WEIGHT, KORDER, XKNOT,
                            NCOEF, BSCOEF, Q, WK, FIRST)

    Arguments:  (See BSLS2)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
 void imsl_b5ls2(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat weight[],
                    Mint *korder, Mfloat xknot[], Mint *ncoef, Mfloat bscoef[],
                    Mfloat *q, Mfloat wk[], Mint *first)
#else
 void imsl_b5ls2(ndata, xdata, fdata, weight, korder, xknot,
        	   ncoef, bscoef, q, wk, first)
	Mint            *ndata;
	Mfloat           xdata[], fdata[], weight[];
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ncoef;
	Mfloat           bscoef[], *q, wk[];
	Mint            *first;
#endif
{
#define Q(I_,J_)	(q+(I_)*(*korder)+(J_))
	Mint             i, left, mm;
	Mfloat           dw;

	/*
	 * Zero out all entries of Q and BSCOEF. If FIRST is true then
	 * compute the least squares matrix; otherwise compte BSCOEF only.
	 */
	if (*first)
		sset(*korder ** ncoef, F_ZERO, q, 1);
	sset(*ncoef, F_ZERO, bscoef, 1);

	left = *korder;
	for (i = 1; i <= *ndata; i++) {
		/*
		 * Find LEFT in the closed interval (I,I+KORDER-1) such that
		 * XKNOT(LEFT) .LE. XDATA(I) .LT. XKNOT(LEFT+1)
		 */
L_10:
		if (left < *ncoef && xdata[i - 1] >= xknot[left]) {
			left += 1;
			goto L_10;
		}
		/*
		 * Find WK(MM) = B(LEFT-KORDER+MM) (XDATA(I))
		 */
		imsl_b4int(xknot, korder, &xdata[i - 1], &left, &wk[0], &wk[*korder],
			   &wk[*korder * 2]);
		for (mm = 1; mm <= *korder; mm++) {
			dw = wk[mm - 1] * weight[i - 1];
			bscoef[left - *korder + mm - 1] += dw * fdata[i - 1];
			if (*first)
				saxpy(*korder - mm + 1, dw, &wk[mm - 1], 1, Q(left - *korder + mm - 1, 0),
					   1);
		}
	}
	return;
}				/* end of function */
#undef  Q
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B6LS2/DB6LS2 (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 17, 1986

    Purpose:    Solve system AINV*TRANS(X) giving BSCOEF.

    Usage:      CALL B6LS2 (NXCOEF, NYCOEF, KXORD, A, WBCOEF, BSCOEF)

    Arguments:  (See B4LS2)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
 void imsl_b6ls2(Mint *nxcoef, Mint *nycoef, Mint *kxord, Mfloat *a,
                    Mfloat *wbcoef, Mfloat *bscoef)
#else
 void imsl_b6ls2(nxcoef, nycoef, kxord, a, wbcoef, bscoef)
	Mint            *nxcoef, *nycoef, *kxord;
	Mfloat          *a, *wbcoef, *bscoef;
#endif
{
#define A(I_,J_)	(a+(I_)*(*kxord)+(J_))
#define WBCOEF(I_,J_)	(wbcoef+(I_)*(*nxcoef)+(J_))
#define BSCOEF(I_,J_)	(bscoef+(I_)*(*nxcoef)+(J_))
	Mint             i;

	/*
	 * Solve system with NYCOEF right hand sides. The RHS are in WBCOEF,
	 * the answer is put back in WBCOEF and copied into BSCOEF.
	 * 
	 * Factor
	 */
	imsl_b5lsq(a, kxord, nxcoef);
	for (i = 1; i <= *nycoef; i++) {
		/* Solve */
		imsl_b6lsq(a, kxord, nxcoef, WBCOEF(i - 1, 0));
		scopy(*nxcoef, WBCOEF(i - 1, 0), 1, BSCOEF(i - 1, 0), 1);
	}
	return;
}				/* end of function */
#undef A
#undef WBCOEF
#undef BSCOEF
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B7LS2/DB7LS2 (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 17, 1986

    Purpose:    Perform error checking for 2-dimensional tensor product
                B-spline.

    Usage:      CALL B7LS2 (ICODE)

    Arguments:
       ICODE  - Integer flag containing an error code.  (Input)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
 void imsl_b7ls2(Mint *icode)
#else
 void imsl_b7ls2(icode)
	Mint            *icode;
#endif
{


	if (*icode == 7) {

		imsl_ermes(IMSL_FATAL, IMSL_NEGATIVE_WEIGHTS);
	} else if (*icode == 8) {
		/*
		 * THIS CHECK WAS ELIMINATED FOR BSLS2 AS OF 11/89. CALL
		 * E1MES (5, 8, 'Points in the %(L1) data point '// &
		 * 'array, %(L1)DATA, must be distinct, but '// &
		 * '%(L1)DATA(%(I1)) = %(L1)DATA(%(I2)) = %(R1).')
		 */
	} else if (*icode == 9) {

		imsl_ermes(IMSL_FATAL, IMSL_SPLINE_NONDEC_DATA);

	} else if (*icode == 10) {

		imsl_ermes(IMSL_FATAL, IMSL_SPLINE_SMLST_ELEMNT);

	} else if (*icode == 11) {

		imsl_ermes(IMSL_FATAL, IMSL_SPLINE_LRGST_ELEMNT);

	} else if (*icode == 12) {
	        imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_NO_POS_ELMNT);

	}
	return;
}				/* end of function */
