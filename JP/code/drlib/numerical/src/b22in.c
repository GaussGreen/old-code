#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  B22IN/DB22IN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 6, 1986

    Purpose:    Compute a two-dimensional tensor product spline
                interpolant, returning the tensor product B-spline
                representation.

    Usage:      CALL B22IN (NXDATA, XDATA, NYDATA, YDATA, FDATA, LDF,
                            KXORD, KYORD, XKNOT, YKNOT, BSCOEF, WK, IWK)

    Arguments:
       NXDATA - Number of data points in the X-direction.  (Input)
       XDATA  - Array of length NXDATA containing the data points
                in the X-direction.  (Input)
                XDATA must be strictly increasing.
       NYDATA - Number of data points in the Y-direction.  (Input)
       YDATA  - Array of length NYDATA containing the data points
                in the Y-direction.  (Input)
                YDATA must be strictly increasing.
       FDATA  - Array of size NXDATA by NYDATA containing the values
                to be interpolated.  (Input)
                FDATA(I,J) is the value at (XDATA(I),YDATA(J)).
       LDF    - The leading dimension of FDATA exactly as specified in
                the dimension statement of the calling program.  (Input)
       KXORD  - Order of the spline in the X-direction.  (Input)
                KXORD must be less than or equal to NXDATA.
       KYORD  - Order of the spline in the Y-direction.  (Input)
                KYORD must be less than or equal to NYDATA.
       XKNOT  - Array of length NXDATA+KXORD containing the knot sequence
                in the X-direction.  (Input)
                XKNOT must be nondecreasing.
       YKNOT  - Array of length NYDATA+KYORD containing the knot sequence
                in the Y-direction.  (Input)
                YKNOT must be nondecreasing.
       BSCOEF - Array of length NXDATA*NYDATA containing the
                tensor-product B-spline coefficients.  (Output)
                BSCOEF is treated internally as a matrix of size
                NXDATA by NYDATA.
       WK     - Work array of length NXDATA*NYDATA +
                MAX((2*KXORD-1)*NXDATA,(2*KYORD-1)*NYDATA) +
                MAX((3*KXORD-2)*NXDATA,(3*KYORD-2)*NYDATA) +
                2*MAX(NXDATA,NYDATA).
       IWK    - Work array of length MAX(NXDATA,NYDATA).

    Remark:
       Informational errors
       Type Code
         3   1  Interpolation matrix is nearly singular.  LU
                factorization failed.
         3   2  Interpolation matrix is nearly singular.  Iterative
                refinement failed.
         4   6  The XDATA values must be strictly increasing.
         4   7  The YATA values must be strictly increasing.
         4   13 Multiplicity of the knots cannot exceed the order
                of the spline.
         4   14 The knots must be nondecreasing.
         4   15 The Ith smallest element of the data point array
                must be greater than the Ith knot and less than
                the (I+K_ORD)th knot.
         4   16 The largest element of the data point array
                must be greater than the (N_DATA)th knot and less
                than or equal to the (N_DATA+K_ORD)th knot.
         4   17 The smallest element of the data point array
                must be greater than or equal to the first knot
                and less than the (K_ORD+1)st knot.

    GAMS:       E2a

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b22in(Mint *nxdata, Mfloat xdata[], Mint *nydata, Mfloat ydata[], Mfloat *fdata, Mint *ldf,
	   Mint *kxord, Mint *kyord, Mfloat xknot[], Mfloat yknot[], Mfloat bscoef[], Mfloat wk[], Mint iwk[])
#else
void imsl_b22in(nxdata, xdata, nydata, ydata, fdata, ldf,
	   kxord, kyord, xknot, yknot, bscoef, wk, iwk)
	Mint            *nxdata;
	Mfloat           xdata[];
	Mint            *nydata;
	Mfloat           ydata[], *fdata;
	Mint            *ldf, *kxord, *kyord;
	Mfloat           xknot[], yknot[], bscoef[], wk[];
	Mint             iwk[];
#endif
{
#define FDATA(I_,J_)	(fdata+(I_)*(*ldf)+(J_))
	Mint             _l0, i, iaband, ifac, iwork1, iwork2, iwork3, laband,
	                lfac, lwork1, lwork2;


	imsl_e1psh("IMSL_B22IN ");
	/* Check X arguments */
	if (*kxord < 1) {
		imsl_e1sti(1, *kxord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_X);
	}
	/* Check KYORD */
	if (*kyord < 1) {
		imsl_e1sti(1, *kyord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_Y);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Check NXDATA */
	if (*nxdata < *kxord) {
		imsl_e1sti(1, *nxdata);
		imsl_e1sti(2, *kxord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_DATA_X);
	}
	/* Check NYDATA */
	if (*nydata < *kyord) {
		imsl_e1sti(1, *nydata);
		imsl_e1sti(2, *kyord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_DATA_Y);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Check LDF */
	if (*ldf < *nxdata) {
		imsl_e1sti(1, *ldf);
		imsl_e1sti(2, *nxdata);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_LD_FDATA);
		goto L_9000;
	}
	/* Check that XDATA is increasing */
	for (i = 2; i <= *nxdata; i++) {
		if (xdata[i - 2] >= xdata[i - 1]) {
			imsl_e1sti(1, i - 2);
			imsl_e1sti(2, i - 1);
			imsl_e1str(1, xdata[i - 2]);
			imsl_e1str(2, xdata[i - 1]);

			imsl_ermes(IMSL_FATAL, IMSL_XDATA_NOT_INCREASING);
			goto L_9000;
		}
	}
	/* Check that YDATA is increasing */
	for (i = 2; i <= *nydata; i++) {
		if (ydata[i - 2] >= ydata[i - 1]) {
			imsl_e1sti(1, i - 2);
			imsl_e1sti(2, i - 1);
			imsl_e1str(1, ydata[i - 2]);
			imsl_e1str(2, ydata[i - 1]);

			imsl_ermes(IMSL_FATAL, IMSL_YDATA_NOT_INCREASING);
			goto L_9000;
		}
	}
	imsl_c1not("X", "KXORD", nxdata, xdata,
		   kxord, xknot);
	imsl_c1not("Y", "KYORD", nydata, ydata,
		   kyord, yknot);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Partition workspace */
	laband = imsl_i_max((2 ** kxord - 1) ** nxdata, (2 ** kyord - 1) ** nydata);
	lfac = imsl_i_max((3 ** kxord - 2) ** nxdata, (3 ** kyord - 2) ** nydata);
	lwork1 = imsl_i_max(*nxdata, *nydata);
	lwork2 = lwork1;
	iaband = 1;
	ifac = iaband + laband;
	iwork1 = ifac + lfac;
	iwork2 = iwork1 + lwork1;
	iwork3 = iwork2 + lwork2;
	/*
	 * B-spline interpolation in the X-direction results put into WORK3
	 */
        _l0 = 2 ** kxord -1;
	imsl_b42in("X", nxdata, nydata, xdata, fdata, ldf, kxord,
	  xknot, &wk[iwork3 - 1], &wk[iaband - 1], &_l0,
		   &wk[ifac - 1], &wk[iwork1 - 1], &wk[iwork2 - 1], iwk);
	/* Check for errors */
	if (imsl_n1rty(0) != 0 && imsl_n1rty(0) != 3)
		goto L_9000;
	/*
	 * B-spline interpolation in the Y-direction.
	 */
        _l0 = 2 ** kyord -1;
	imsl_b42in("Y",nydata, nxdata, ydata, &wk[iwork3 - 1],
	nydata, kyord, yknot, bscoef, &wk[iaband - 1], &_l0,
        &wk[ifac - 1], &wk[iwork1 - 1], &wk[iwork2 - 1], iwk);

L_9000:
	;
	imsl_e1pop("IMSL_B22IN ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B32IN/DB32IN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 6, 1986

    Purpose:    Check parameters for tensor product B-spline
                interpolation routines.

    Usage:      CALL B32IN (VNAME, KORDER, XKNOT, NDATA)

    Arguments:
       VNAME  - Character string containing name of dimension being
                checked.  (Input)
       KORDER - Order of the spline.  (Input)
       XKNOT  - Array of length NDATA+KORDER containing the knot
                sequence.  (Input)
                It must be nondecreasing.
       NDATA  - Number of data points.  (Input)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b32in(Mchar *vname, Mint *korder, Mfloat xknot[], Mint *ndata)
#else
void imsl_b32in(vname, korder, xknot, ndata)
	Mchar           *vname;
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *ndata;
#endif
{
	Mint             i, mult;

	/* Check KORDER */
	if (*korder <= 0) {
		imsl_e1sti(1, *korder);
		imsl_e1stl(1, vname);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_X_OR_Y);
		goto L_9000;
	}
	/* Check argument NDATA */
	if (*ndata < *korder) {
		imsl_e1sti(1, *ndata);
		imsl_e1sti(2, *korder);
		imsl_e1stl(1, vname);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_XYDATA);
		goto L_9000;
	}
	/* Check knot sequence */
	mult = 1;
	for (i = 2; i <= (*ndata + *korder); i++) {
		if (xknot[i - 1] == xknot[i - 2]) {
			mult += 1;
			if (mult > *korder) {
				imsl_e1sti(1, (i-1) - mult + 1);
				imsl_e1sti(2, i - 1);
				imsl_e1str(1, xknot[i - 1]);
				imsl_e1sti(3, *korder);
				imsl_e1stl(1, vname);

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

			imsl_ermes(IMSL_FATAL, IMSL_KNOT_NOT_INCREASING);
			goto L_9000;
		} else {
			mult = 1;
		}
	}

L_9000:
	;
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B42IN/DB42IN (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 6, 1986

    Purpose:    Compute a two-dimensional tensor product spline
                interpolant, returning the tensor product B-spline
                representation.

    Usage:      CALL B42IN (CHR, NXDATA, NYDATA, XDATA, FDATA, LDF,
                            KORDER, XKNOT, BSCOEF, ABAND, LBAND, FAC,
                            WORK1, WORK2, IPVT)

    Arguments:  (See BS2IN)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b42in(Mchar *chr, Mint *nxdata, Mint *nydata, Mfloat xdata[], Mfloat *fdata,
	   Mint *ldf, Mint *korder, Mfloat xknot[], Mfloat *bscoef, Mfloat *aband,
           Mint *lband, Mfloat imsl_fac[], Mfloat work1[], Mfloat work2[],
	   Mint ipvt[])
#else
void imsl_b42in(chr, nxdata, nydata, xdata, fdata,
	   ldf, korder, xknot, bscoef, aband, lband, imsl_fac, work1, work2,
	   ipvt)
 	Mchar           *chr;
	Mint            *nxdata, *nydata;
	Mfloat           xdata[], *fdata;
	Mint            *ldf, *korder;
	Mfloat           xknot[], *bscoef, *aband;
	Mint            *lband;
	Mfloat           imsl_fac[], work1[], work2[];
	Mint             ipvt[];
#endif
{
#define FDATA(I_,J_)	(fdata+(I_)*(*ldf)+(J_))
#define BSCOEF(I_,J_)	(bscoef+(I_)*(*nydata)+(J_))
#define ABAND(I_,J_)	(aband+(I_)*(*lband)+(J_))
	Mint             _l0, _l1, _l2, _l3, i, j, left;
	Mfloat           rcond, taui;


	imsl_e1psh("IMSL_B42IN ");
	/* Zero out all entries of ABAND */
	sset(*nxdata ** lband, F_ZERO, aband, 1);
	/*
	 * Loop over I to construct the NXDATA interpolation equations
	 */
	left = *korder;
	for (i = 1; i <= *nxdata; i++) {
		taui = xdata[i - 1];
		/*
		 * Find LEFT in the closed interval (I,I+KORDER-1) such that
		 * XKNOT(LEFT) .LE. XDATA(I) .LT. XKNOT(LEFT+1) Matrix is
		 * singular if this is not possible.
		 */
		left = imsl_i_max(left, i);
L_10:
		if (taui >= xknot[left]) {
			left += 1;
			if (left < imsl_i_min(i + *korder, *nxdata + 1))
				goto L_10;
			left -= 1;
		}
		/*
		 * The I-TH equation enforces interpolation at TAUI, hence
		 * A(I,J) = B(J,KORDER,XKNOT)(TAUI), ALL J. Only the KORDER
		 * entries with J = LEFT-KORDER+1, ..., LEFT actually might
		 * be nonzero. These KORDER numbers are returned by the
		 * following in WORK2. FAC is used as workspace.
		 */
		imsl_b4int(xknot, korder, &taui, &left, work2, &imsl_fac[0], &imsl_fac[*korder]);
		/*
		 * We therefore want WORK2(J) = B(LEFT-K+J)(TAUI) to go into
		 * A(I,LEFT-KORDER+J) = ABAND(I-LEFT+2*KORDER-J,
		 * LEFT-KORDER+J) since the number of lower codiagonals is
		 * KORDER-1.
		 */
		scopy(*korder, work2, 1, ABAND(left - *korder, i - left + *korder * 2 - 2),
			   *lband - 1);
	}
	/* Factor the linear system */
        _l0 = *korder -1;
        _l1 = *korder -1;
        _l2 = 3 ** korder -2;
	imsl_l2crb(nxdata, aband, lband, &_l0, &_l1, imsl_fac, &_l2, ipvt, &rcond, work1);
	if (imsl_n1rty(0) != 0) {
		imsl_e1stl(1, chr);

		imsl_ermes(IMSL_WARNING, IMSL_ILL_COND_INTERP_PROB);
	}
	/*
	 * Solve the linear system with iterative refinement.
	 */
        _l0 = 3 ** korder -2;
        _l1 = *korder -1;
        _l2 = *korder -1;
        _l3 = 1;
	for (j = 1; j <= *nydata; j++) {
		imsl_lfsrb(nxdata, imsl_fac, &_l0, &_l1, &_l2, ipvt, FDATA(j - 1, 0), &_l3,work2);
							
/* --------------------------------------------------------------------------------------*		    
*       ITERATIVE REFINEMENT WAS ELIMINATED FOR CLIB 			                 *
*		imsl_lfirb(nxdata, aband, lband, ADR(_l0, *korder - 1),                  *
*                           ADR(_l1, *korder - 1), imsl_fac, ADR(_l2, 3 ** korder - 2),  *
*                           ipvt, FDATA(j - 1, 0), ADR(_l3, 1), work2, work1);           *
*		if (imsl_n1rcd(1) != 0) {                                                *
*			imsl_e1stl(1, chr);                                              *
*			imsl_ermes(IMSL_WARNING, IMSL_ILL_COND_INTERP_PROB);             *
*		}                                                                        *
*----------------------------------------------------------------------------------------*		    
*/
		scopy(*nxdata, work2, 1, BSCOEF(0, j - 1), *nydata);
	}
	imsl_e1pop("IMSL_B42IN ");
	return;
}				/* end of function */
