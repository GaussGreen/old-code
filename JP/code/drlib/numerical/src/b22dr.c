#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  B22DR/DB22DR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 17, 1990

    Purpose:    Evaluate the derivatives of a two-dimensional tensor
                product spline, given its tensor product B-spline
                representation.

    Usage:      B22DR(IXDER, IYDER, X, Y, KXORD, KYORD, XKNOT, YKNOT,
                      NXCOEF, NYCOEF, BSCOEF, WK)

    Arguments:
       IXDER  - Order of the derivative in the X-direction.  (Input)
       IYDER  - Order of the derivative in the Y-direction.  (Input)
       X      - X-coordinate of the point at which the spline is to be
                evaluated.  (Input)
       Y      - Y-coordinate of the point at which the spline is to be
                evaluated.  (Input)
       KXORD  - Order of the spline in the X-direction.  (Input)
       KYORD  - Order of the spline in the Y-direction.  (Input)
       XKNOT  - Array of length NXCOEF+KXORD containing the knot
                sequence in the X-direction.  (Input)
                XKNOT must be nondecreasing.
       YKNOT  - Array of length NYCOEF+KYORD containing the knot
                sequence in the Y-direction.  (Input)
                YKNOT must be nondecreasing.
       NXCOEF - Number of B-spline coefficients in the X-direction.
                (Input)
       NYCOEF - Number of B-spline coefficients in the Y-direction.
                (Input)
       BSCOEF - Array of length NXCOEF*NYCOEF containing the
                tensor-product B-spline coefficients.  (Input)
                BSCOEF is treated internally as a matrix of size
                NXCOEF by NYCOEF.
       WK     - Work array of length 3*MAX(KXORD,KYORD) + KYORD.
       B22DR  - Value of the (IXDER,IYDER) derivative of the spline
                at (X,Y).  (Output)

    Remark:
       Informational errors
       Type Code
         3   1  The point X is not in the closed interval XKNOT(KXORD)
                to XKNOT(NXCOEF+1).
         3   2  The point Y is not in the closed interval
                YKNOT(KYORD) to YKNOT(NYCOEF+1).

    Keyword:    Differentiate

    GAMS:       E3

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_b22dr(Mint *ixder, Mint *iyder, Mfloat *x, Mfloat *y,
           Mint *kxord, Mint *kyord, Mfloat xknot[],
	   Mfloat yknot[], Mint *nxcoef, Mint *nycoef, Mfloat *bscoef, Mfloat wk[])
#else
Mfloat imsl_b22dr(ixder, iyder, x, y, kxord, kyord, xknot,
	   yknot, nxcoef, nycoef, bscoef, wk)
	Mint            *ixder, *iyder;
	Mfloat          *x, *y;
	Mint            *kxord, *kyord;
	Mfloat           xknot[], yknot[];
	Mint            *nxcoef, *nycoef;
	Mfloat          *bscoef, wk[];
#endif
{
#define BSCOEF(I_,J_)	(bscoef+(I_)*(*nxcoef)+(J_))
	Mint             _l0, iaj, idl, idr, iybscf, j, lefty, mflag, mxkord;
	Mfloat           b22dr_v, value;


	imsl_e1psh("IMSL_B22DR ");
	/* In case of errors */
	value = F_ZERO;
	/* Check KXORD */
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
	/* Check IXDER */
	if (*ixder < 0) {
		imsl_e1sti(1, *ixder);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_DERIV_X);
		goto L_9000;
	}
	/* Check IYDER */
	if (*iyder < 0) {
		imsl_e1sti(1, *iyder);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_DERIV_Y); 
		goto L_9000;
	}
	/* Check NXCOEF */
	if (*nxcoef < *kxord) {
		imsl_e1sti(1, *nxcoef);
		imsl_e1sti(2, *kxord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_COEFF_X);
	}
	/* Check NYCOEF */
	if (*nycoef < *kyord) {
		imsl_e1sti(1, *nycoef);
		imsl_e1sti(2, *kyord);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_COEFF_Y);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/*
	 * Check X and Y in proper interval Note at a later date this will
	 * probably come out and we will make B3DER a nuclei and not call
	 * B3DER the last time.
	 */
	if (*x < xknot[*kxord - 1] || *x > xknot[*nxcoef]) {
		imsl_e1str(1, *x);

		imsl_ermes(IMSL_WARNING, IMSL_X_NOT_WITHIN_KNOTS);
		value = F_ZERO;
		goto L_9000;
	}
	if (*y < yknot[*kyord - 1] || *y > yknot[*nycoef]) {
		imsl_e1str(1, *y);

		imsl_ermes(IMSL_WARNING, IMSL_Y_NOT_WITHIN_KNOTS);
		value = F_ZERO;
		goto L_9000;
	}
	/* Partition workspace */
	mxkord = imsl_i_max(*kxord, *kyord);
	iaj = 1;
	idl = iaj + mxkord;
	idr = idl + mxkord;
	iybscf = idr + mxkord;
	/* Find the correct Y interval */
        _l0 =  *nycoef + *kyord;
	imsl_b4der(yknot,&_l0, y, &lefty, &mflag);
	if (mflag != 0)
		goto L_9000;
	/*
	 * Interpolate to compute B-spline coefs.
	 */
	for (j = 1; j <= *kyord; j++) {
		wk[iybscf + j - 2] = imsl_b3der(ixder, x, kxord, xknot, nxcoef,
						BSCOEF(lefty - *kyord + j - 1, 0), &wk[iaj - 1], &wk[idl - 1],
						&wk[idr - 1]);
		if (imsl_n1rty(0) != 0)
			goto L_9000;
	}
	/*
	 * Interpolate to find the value at (X,Y).
	 */
	value = imsl_b3der(iyder, y, kyord, &yknot[lefty - *kyord], kyord,
		 &wk[iybscf - 1], &wk[iaj - 1], &wk[idl - 1], &wk[idr - 1]);

L_9000:
	;
	b22dr_v = value;
	imsl_e1pop("IMSL_B22DR ");
	return (b22dr_v);
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  B32DR/DB32DR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 7, 1986

    Purpose:    Check parameters for tensor product B-spline
                evaluation routines.

    Usage:      CALL B32DR (VNAME, KORDER, NCOEF)

    Arguments:
       VNAME  - Character string containing name of dimension being
                checked.  (Input)
       KORDER - Order of the B-spline.  (Input)
       NCOEF  - Number of B-spline coefficients.  (Input)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b32dr(Mchar *vname, Mint *korder, Mint *ncoef)
#else
void imsl_b32dr(vname, korder, ncoef)
	Mchar           *vname;
	Mint            *korder, *ncoef;
#endif
{

	/* Check KORDER */
	if (*korder <= 0) {
		imsl_e1sti(1, *korder);
		imsl_e1stl(1, vname);

                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_ARB);
		goto L_9000;
	}
	/* Check NCOEF */
	if (*ncoef < *korder) {
		imsl_e1sti(1, *ncoef);
		imsl_e1sti(2, *korder);
		imsl_e1stl(1, vname);

		imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_COEFF_XY);
	}
L_9000:
	return;
}				/* end of function */
