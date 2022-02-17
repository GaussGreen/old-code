#define ADR(t,x)    ( t = x, &t )

#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static void l_b32gd (Mfloat t[], Mint *k, Mfloat *x, Mint *left,
                Mfloat *a, Mfloat *dbiatx, Mint *nderiv);
static void l_b42gd (Mfloat t[], Mint *jhigh, Mint *index, Mfloat *x,
                Mint *left, Mfloat biatx[]);
#else
static void l_b32gd();
static void l_b42gd();
#endif


/* Structured by FOR_STRUCT, v0.2, on 11/04/91 at 14:35:24
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  B22GD/DB22GD (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 17, 1990

    Purpose:    Evaluate the derivative of a two-dimensional
                tensor-product spline, given its tensor-product
                B-spline representation on a grid.

    Usage:      CALL B22GD (IXDER, IYDER, NX, XVEC, NY, YVEC, KXORD,
                            KYORD, XKNOT, YKNOT, NXCOEF, NYCOEF,
                            BSCOEF, VALUE, LDVALU, LEFTX, LEFTY, A, B,
                            DBIATX, DBIATY, BX, BY)

    Arguments:
       IXDER  - Order of the derivative in the X-direction.  (Input)
       IYDER  - Order of the derivative in the Y-direction.  (Input)
       NX     - Number of grid points in the x-direction.  (Input)
       XVEC   - Array of length NX containing the x-coordinates at
                which the spline is to be evaluated.  (Input)
                The points in XVEC should be strictly increasing.
       NY     - Number of grid points in the y-direction.  (Input)
       YVEC   - Array of length NY containing the y-coordinates at
                which the spline is to be evaluated.  (Input)
                The points in YVEC should be strictly increasing.
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
       VALUE  - Value of the (IXDER,IYDER) derivative of the spline
                on the NX by NY grid.  (Output)
                VALUE(I,J)  contains the derivative of the spline at
                the point (XVEC(I), YVEC(J)).
       LDVALU - Leading dimension of VALUE exactly as specified in
                the dimension statement of the calling program.  (Input)
       LEFTX  - Integer work array of length NX
       LEFTY  - Integer work array of length NY
       A      - Work array of length KXORD*KXORD
       B      - Work array of length KYORD*KYORD
       DBIATX - Work array of length KXORD*(IXDER+1)
       DBIATY - Work array of length KYORD*(IYDER+1)
       BX     - Work array of length KXORD*NX
       BY     - Work array of length KYORD*NY

    Remark:
       Informational errors
       Type Code
         3   1  XVEC(I) does not satisfy XKNOT(KXORD) .LE. XVEC(I) .LE.
                XKNOT(NXCOEF+1)
         3   2  YVEC(I) does not satisfy YKNOT(KXORD) .LE. YVEC(I) .LE.
                YKNOT(NYCOEF+1)
         4   3  XVEC is not strictly increasing.
         4   4  YVEC is not strictly increasing.

    GAMS:       E3

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b22gd (Mint *ixder, Mint *iyder, Mint *nx, Mfloat xvec[],
		Mint *ny, Mfloat yvec[], Mint *kxord, Mint *kyord, 
		Mfloat xknot[], Mfloat yknot[],  Mint *nxcoef,
		Mint *nycoef, Mfloat *bscoef, Mfloat *value, Mint *ldvalu,
		Mint leftx[], Mint lefty[], Mfloat *a, Mfloat *b,
		Mfloat *dbiatx, Mfloat *dbiaty, Mfloat *bx, Mfloat *by)
#else
void imsl_b22gd (ixder, iyder, nx, xvec, ny, yvec, kxord, kyord,
                xknot, yknot, nxcoef, nycoef, bscoef, value, ldvalu, leftx, lefty,
                a, b, dbiatx, dbiaty, bx, by)
    Mint        *ixder, *iyder, *nx;
    Mfloat       xvec[];
    Mint        *ny;
    Mfloat       yvec[];
    Mint        *kxord, *kyord;
    Mfloat       xknot[], yknot[];
    Mint        *nxcoef, *nycoef;
    Mfloat      *bscoef, *value;
    Mint        *ldvalu, leftx[], lefty[];
    Mfloat      *a, *b, *dbiatx, *dbiaty, *bx, *by;
#endif
{
#define BSCOEF(I_,J_)	(bscoef+(I_)*(*nxcoef)+(J_))
#define VALUE(I_,J_)	(value+(I_)*(*ldvalu)+(J_))
#define A(I_,J_)	(a+(I_)*(*kxord)+(J_))
#define B(I_,J_)	(b+(I_)*(*kyord)+(J_))
#define DBIATX(I_,J_)	(dbiatx+(I_)*(*kxord)+(J_))
#define DBIATY(I_,J_)	(dbiaty+(I_)*(*kyord)+(J_))
#define BX(I_,J_)	(bx+(I_)*(*kxord)+(J_))
#define BY(I_,J_)	(by+(I_)*(*kyord)+(J_))
    Mint         _l0, i, ikx, iky, ix, iy, j, mflag;


    imsl_e1psh ("B22GD");
    /* In case of errors */
    sset (*nx ** ny, 0.0, value, 1);
    /* Check KXORD */
    if (*kxord < 1) {
	imsl_e1sti (1, *kxord);
        imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_X);
    }
    /* Check KYORD */
    if (*kyord < 1) {
	imsl_e1sti (1, *kyord);
        imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_Y);
    }
    if (imsl_n1rty (0) != 0)
	goto L_9000;
    /* Check IXDER */
    if (*ixder < 0) {
	imsl_e1sti (1, *ixder);
        imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_DERIV_X);
	goto L_9000;
    }
    /* Check IYDER */
    if (*iyder < 0) {
	imsl_e1sti (1, *iyder);
        imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_DERIV_Y);
	goto L_9000;
    }
    /* Check NXCOEF */
    if (*nxcoef < *kxord) {
	imsl_e1sti (1, *nxcoef);
	imsl_e1sti (2, *kxord);
        imsl_ermes(IMSL_TERMINAL,IMSL_SPLINE_COEFF_X);
    }
    /* Check NYCOEF */
    if (*nycoef < *kyord) {
	imsl_e1sti (1, *nycoef);
	imsl_e1sti (2, *kyord);
        imsl_ermes(IMSL_TERMINAL,IMSL_SPLINE_COEFF_Y);
    }
    if (imsl_n1rty (0) != 0)
	goto L_9000;
    /*
     * Check X and Y in proper interval Note at a later date this will
     * probably come out and we will make B3DER a nuclei and not call B3DER
     * the last time.
     */
    for (i = 1; i <= *nx; i++) {
	if (xvec[i - 1] < xknot[*kxord - 1] || xvec[i - 1] > xknot[*nxcoef]) {
	    imsl_e1str (1, xvec[i - 1]);
	    imsl_e1sti (2, i-1);
            imsl_ermes(IMSL_WARNING, IMSL_X_NOT_WITHIN_KNOTS);
	    goto L_9000;
	}
    }

    for (i = 1; i <= *ny; i++) {
	if (yvec[i - 1] < yknot[*kyord - 1] || yvec[i - 1] > yknot[*nycoef]) {
	    imsl_e1str (1, yvec[i - 1]);
	    imsl_e1sti (2, i);
            imsl_ermes(IMSL_WARNING, IMSL_Y_NOT_WITHIN_KNOTS);
	    goto L_9000;
	}
    }

    for (ix = 1; ix <= *nx; ix++) {
	imsl_b4der (xknot, ADR (_l0, *kxord + *nxcoef), &xvec[ix - 1], &leftx[ix - 1],
	    &mflag);
    }
    for (iy = 1; iy <= *ny; iy++) {
	imsl_b4der (yknot, ADR (_l0, *kyord + *nycoef), &yvec[iy - 1], &lefty[iy - 1],
	    &mflag);
    }

    for (ix = 1; ix <= *nx; ix++) {
	l_b32gd (xknot, kxord, &xvec[ix - 1], &leftx[ix - 1], a, dbiatx,
	    ADR (_l0, *ixder + 1));
	for (ikx = 1; ikx <= *kxord; ikx++) {
	    *BX (ix - 1, ikx - 1) = *DBIATX (*ixder, ikx - 1);
	}
    }

    for (iy = 1; iy <= *ny; iy++) {
	l_b32gd (yknot, kyord, &yvec[iy - 1], &lefty[iy - 1], b, dbiaty,
	    ADR (_l0, *iyder + 1));
	for (iky = 1; iky <= *kyord; iky++) {
	    *BY (iy - 1, iky - 1) = *DBIATY (*iyder, iky - 1);
	}
    }

    for (ix = 1; ix <= *nx; ix++) {
	for (iy = 1; iy <= *ny; iy++) {
	    *VALUE (iy - 1, ix - 1) = 0.0;
	    for (i = 1; i <= *kxord; i++) {
		for (j = 1; j <= *kyord; j++) {
		    *VALUE (iy - 1, ix - 1) += *BSCOEF (lefty[iy - 1] -
			*kyord + j - 1, leftx[ix - 1] - *kxord + i - 1) *
			*BX (ix - 1, i - 1) ** BY (iy - 1, j - 1);
		}
	    }
	}
    }

L_9000:
    ;
    imsl_e1pop ("B22GD");
    return;
}				/* end of function */

#undef BSCOEF
#undef VALUE
#undef A
#undef B
#undef DBIATX
#undef DBIATY
#undef BX
#undef BY

/*----------------------------------------------------------------------- */

/*  IMSL Name:  B32GD/DB32GD (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 31, 1990

    Purpose:    Calculate the value and derivatives of all B-splines
                which do not vanish at X.

    Usage:      CALL B32GD (T, K, X, LEFT, A, DBIATX, NDERIV)

    Arguments:
       T      - Array of length LEFT+K containing the knot sequence.
                (Input)
       K      - Order of the spline.  (Input)
       X      - Point at which derivative is to be evaluated.  (Input)
       LEFT   - The left endpoint of the interval of interest.  (Input)
                The K B-splines whose support contains the interval
                (T(LEFT,T(LEFT+1)) are to be considered.
       A      - Work array of order (K,K).  (Input)
       DBIATX - Array of order (K,NDERIV).  (Output)
                Entry (I,M) contains the value of (M-1)st derivative
                of (LEFT-K+I)th B-spline of order K for knot sequence
                T, I = 1, ... , NDERIV.
       NDERIV - Integer indicating that values of B-splines and their
                derivatives up to but not including the NDERIV-th are
                requested.  (Input)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_b32gd (Mfloat t[], Mint *k, Mfloat *x, Mint *left,
		Mfloat *a, Mfloat *dbiatx, Mint *nderiv)
#else
static void l_b32gd (t, k, x, left, a, dbiatx, nderiv)
    Mfloat       t[];
    Mint        *k;
    Mfloat      *x;
    Mint        *left;
    Mfloat      *a, *dbiatx;
    Mint        *nderiv;
#endif
{
#define A(I_,J_)	(a+(I_)*(*k)+(J_))
#define DBIATX(I_,J_)	(dbiatx+(I_)*(*k)+(J_))
    Mint         _l0, _l1, i, ideriv, il, j, jlow, jp1mid, kp1, kp1mm, ldummy, m,
                mhigh;
    Mfloat       factor, fkp1mm, sum;


    mhigh = imsl_i_max (imsl_i_min (*nderiv, *k), 1);

    kp1 = *k + 1;
    l_b42gd (t, ADR (_l0, kp1 - mhigh), ADR (_l1, 1), x, left, dbiatx);
    if (mhigh == 1)
	goto L_90;

    ideriv = mhigh;
    for (m = 2; m <= mhigh; m++) {
	jp1mid = 1;
	for (j = ideriv; j <= *k; j++) {
	    *DBIATX (ideriv - 1, j - 1) = *DBIATX (0, jp1mid - 1);
	    jp1mid += 1;
	}
	ideriv -= 1;
	l_b42gd (t, ADR (_l0, kp1 - ideriv), ADR (_l1, 2), x, left, dbiatx);
    }

    jlow = 1;
    for (i = 1; i <= *k; i++) {
	for (j = jlow; j <= *k; j++) {
	    *A (i - 1, j - 1) = 0.0;
	}
	jlow = i;
	*A (i - 1, i - 1) = 1.0;
    }

    for (m = 2; m <= mhigh; m++) {
	kp1mm = kp1 - m;
	fkp1mm = (Mfloat) (kp1mm);
	il = *left;
	i = *k;
	for (ldummy = 1; ldummy <= kp1mm; ldummy++) {
	    factor = fkp1mm / (t[il + kp1mm - 1] - t[il - 1]);
	    for (j = 1; j <= i; j++) {
		*A (j - 1, i - 1) = (*A (j - 1, i - 1) - *A (j - 1, i - 2)) *
		    factor;
	    }
	    il -= 1;
	    i -= 1;
	}

	for (i = 1; i <= *k; i++) {
	    sum = 0.0;
	    jlow = imsl_i_max (i, m);
	    for (j = jlow; j <= *k; j++) {
		sum += *A (i - 1, j - 1) ** DBIATX (m - 1, j - 1);
	    }
	    *DBIATX (m - 1, i - 1) = sum;
	}
    }
L_90:
    return;
}				/* end of function */

#undef A
#undef DBIATX

/*----------------------------------------------------------------------- */

/*  IMSL Name:  B42GD/DB42GD (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 31, 1990

    Purpose:    Calculate the value of all possibly nonzero B-splines
                at X of order JOUT = MAX(JHIGH,(J+1*(INDEX-1)) with
                knot sequence T.

    Usage:      CALL B42GD (T, JHIGH, INDEX, X, LEFT, BIATX)

    Arguments:
       T      - Array of length JOUT+LEFT containing the knot sequence.
                (Input)
       JHIGH  - Integer used to determine the order.  (Input)
       INDEX  - Integer used to determine the order.  (Input)
       X      - Point at which the B-splines are to be evaluated.
                (Input)
       LEFT   - Integer chosen so that T(LEFT).LE.X.LE.T(LEFT+1).
                (Input)
       BIATX  - Array of length JOUT, with BIATX(I) containing the
                value at X of the polynomial of order JOUT which agrees
                with the B-spline B(LEFT-JOUT+I,JOUT,T) on the interval
                (T(LEFT),T(LEFT+1)).  (Output)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	JMAX	20

#ifdef ANSI
static void l_b42gd (Mfloat t[], Mint *jhigh, Mint *index, Mfloat *x,
		Mint *left, Mfloat biatx[])
#else
static void l_b42gd (t, jhigh, index, x, left, biatx)
    Mfloat       t[];
    Mint        *jhigh, *index;
    Mfloat      *x;
    Mint        *left;
    Mfloat       biatx[];
#endif
{
    Mint         i, jp1;
    Mfloat       saved, term;
    static Mfloat deltal[JMAX];
    static Mfloat deltar[JMAX];
    static Mint  j = 1;



    if (*index == 1) {
	j = 1;
	biatx[0] = 1.0;
	if (j >= *jhigh)
	    goto L_30;
    }
L_10:
    jp1 = j + 1;
    deltar[j - 1] = t[*left + j - 1] - *x;
    deltal[j - 1] = *x - t[*left - j];
    saved = 0.0;
    for (i = 1; i <= j; i++) {
	term = biatx[i - 1] / (deltar[i - 1] + deltal[jp1 - i - 1]);
	biatx[i - 1] = saved + deltar[i - 1] * term;
	saved = deltal[jp1 - i - 1] * term;
    }
    biatx[jp1 - 1] = saved;
    j = jp1;
    if (j < *jhigh)
	goto L_10;
L_30:
    return;
}				/* end of function */
