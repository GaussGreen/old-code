#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static Mfloat PROTO(l_ssum,(Mint *n, Mfloat sx[], Mint *incx));
/*  -----------------------------------------------------------------------
    IMSL Name:  B2OPK/DB2OPK (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 1, 1986

    Purpose:    Compute the 'optimal' spline knot sequence.

    Usage:      CALL B2OPK (NDATA, XDATA, KORDER, XKNOT, MAXIT, WK, IWK)

    Arguments:  (See BSOPK)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b2opk(Mint *ndata, Mfloat xdata[], Mint *korder, Mfloat xknot[],
                Mint *maxit, Mfloat wk[], Mint iwk[])
#else
void imsl_b2opk(ndata, xdata, korder, xknot, maxit, wk, iwk)
	Mint            *ndata;
	Mfloat           xdata[];
	Mint            *korder;
	Mfloat           xknot[];
	Mint            *maxit;
	Mfloat           wk[];
	Mint             iwk[];
#endif
{
	Mint             i, ia, ibiatx, id, idl, idr, iexdat, iw, iwork,
	                ixi, ixsrt, j;


	imsl_e1psh("IMSL_B2OPK ");
	/* Check KORDER */
	if (*korder < 3) {
		imsl_e1sti(1, *korder);

/*		imsl_ermes(5, 1, "The order of the spline must be at least 3 while KORDER = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER_2);
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
				imsl_e1sti(1, j-1);
				imsl_e1sti(2, i-1);
				imsl_e1str(1, xdata[i - 1]);

/*				imsl_ermes(4, 3, "Points in the data point abscissas array, XDATA, must be distinct, but XDATA(%(i1)) = XDATA(%(i2)) = %(r1).");
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
	scopy(*ndata, xdata, 1, wk, 1);
	goto L_50;
	/* Set initial permutation */
L_20:
	for (i = 1; i <= *ndata; i++) {
		iwk[i - 1] = i;
	}
	/* Find sorting permutation */
	imsl_svrgp(*ndata, xdata, wk, iwk);
	/* Check the XDATA values are distinct */
	for (i = 2; i <= *ndata; i++) {
		if (wk[i - 2] == wk[i - 1]) {
			imsl_e1sti(1, iwk[i - 2]-1);
			imsl_e1sti(2, iwk[i - 1]-1);
			imsl_e1str(1, wk[i - 1]);

/*			imsl_ermes(4, 3, "Points in the data point abscissas array, XDATA, must be distinct, but XDATA(%(i1)) = XDATA(%(i2)) = %(r1).");
*/
                        imsl_ermes(IMSL_FATAL, IMSL_DUPLICATE_XDATA_VALUES);
			goto L_9000;
		}
	}
L_50:
	;

	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Partition workspace */
	ixsrt = 1;
	iexdat = ixsrt + *ndata;
	ixi = iexdat + *ndata + 2 ** korder;
	ia = ixi + *ndata - *korder + 2;
	id = ia + *ndata;
	ibiatx = id + *ndata - *korder;
	idl = ibiatx + *korder + 1;
	idr = idl + *korder + 1;
	iwork = idr + *korder + 1;
	iw = iwork + (*ndata - *korder);

	imsl_b3opk(ndata, &wk[ixsrt - 1], korder, xknot, &wk[iexdat - 1],
		   &wk[ixi - 1], &wk[ia - 1], &wk[ibiatx - 1], &wk[id - 1], &wk[idl - 1],
		   &wk[idr - 1], &wk[iwork - 1], maxit, &wk[iw - 1], iwk);

L_9000:
	;
	imsl_e1pop("IMSL_B2OPK ");
	return;
}				/* end of function */
/*  -----------------------------------------------------------------------
    IMSL Name:  B3OPK/DB3OPK (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 7, 1986

    Purpose:    Compute the 'optimal' spline knot sequence.

    Usage:      CALL B3OPK (NDATA, XDATA, KORDER, XKNOT, EXDATA, XI, A,
                            BIATX, D, DL, DR, WORK, MAXIT, W, IPVT)

    Arguments:
       NDATA  - Number of data points.  (Input)
       XDATA  - Array of length NDATA containing the sorted data points.
                (Input)
       KORDER - Order of the spline.  (Input)
       XKNOT  - Array of length NDATA+KORDER containing the knot
                sequence.  (Output)
       EXDATA - Array of length NDATA+2*KORDER for the extended knot
                sequence.  (Workspace)
       XI     - Array of length NDATA-KORDER+2 for the intermediate
                knot locations.  (Workspace)
       A      - Array of length NDATA.  (Workspace)
       BIATX  - Array of length KORDER+1 to hold spline values.
                (Workspace)
       D      - Array of length NDATA-KORDER to hold solution to
                Newton's method.  (Workspace)
       DL     - Array of length KORDER+1 used by B4INT.  (Workspace)
       DR     - Array of length KORDER+1 used by B4INT.  (Workspace)
       WORK   - Array of length NDATA-KORDER used by linear equation
                solver.  (Workspace)
       MAXIT  - Maximum number of iterations of Newton's method.  (Input)
       W      - Array of length (NDATA-KORDER)*(3*KORDER-2) to hold the
                band matrix of order NDATA-KORDER with KORDER-1 upper
                and lower codiagonals.  (Workspace)
                The extra NLCA (KORDER-1)  rows are workspace.
       IPVT   - Array of length NDATA.  (Workspace)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b3opk(Mint *ndata, Mfloat xdata[], Mint *korder, Mfloat xknot[],
                Mfloat exdata[], Mfloat xi[], Mfloat a[],
	        Mfloat biatx[], Mfloat d[], Mfloat dl[], Mfloat dr[], Mfloat work[],
                Mint *maxit, Mfloat w[], Mint ipvt[])
#else
void imsl_b3opk(ndata, xdata, korder, xknot, exdata, xi, a,
	   biatx, d, dl, dr, work, maxit, w, ipvt)
	Mint            *ndata;
	Mfloat           xdata[];
	Mint            *korder;
	Mfloat           xknot[], exdata[], xi[], a[], biatx[], d[], dl[],
	                dr[], work[];
	Mint            *maxit;
	Mfloat           w[];
	Mint             ipvt[];
#endif
{
	Mint             _l0, _l1, _l2, i, j, nbands, newtmx, newton;
	Mfloat           del, delmax, eps, rcond, siign, siignst, tol;


	imsl_e1psh("IMSL_B3OPK ");
	newtmx = *maxit;
	/* Stack endpoints if N=K */
	if (*ndata == *korder)
		goto L_70;
	/*
	 * Extend XDATA to a knot sequence and store in EXDATA.
	 */
	sset(*korder, xdata[0], exdata, 1);
	sset(*korder, xdata[*ndata - 1], &exdata[*korder + *ndata], 1);
	scopy(*ndata, xdata, 1, &exdata[*korder], 1);
	/* First guess for XI */
	xi[0] = xdata[0];
	xi[*ndata - *korder + 1] = xdata[*ndata - 1];
	for (j = 1; j <= (*ndata - *korder); j++) {
                _l0 = *korder - 1;
                _l1 = 1;
		xi[j] = l_ssum(&_l0, &xdata[j], &_l1) /
			(Mfloat) (*korder - 1);
	}
	/* Last entry of - A is always */
	a[*ndata - 1] = F_HALF;
	/* Set initial sign */
	if (mod(*ndata - *korder, 2) == 0) {
		siignst = -F_ONE;
	} else {
		siignst = F_ONE;
	}
	/* Start newton iteration. */
	newton = 1;
	tol = sqrt(imsl_amach(3)) * (xdata[*ndata - 1] - xdata[0]) / (Mfloat) (*ndata -
								   *korder);
	/*
	 * Start Newton step compute the 2K-1 bands of the matrix C and store
	 * in W and compute the vector A
	 */
L_20:
	;
        _l0 = 3 ** korder - 2;
	imsl_b4opk(ndata, korder, &siignst, exdata, xi, a, biatx, dl, dr, w,
		   &_l0);
	/* Compute D from A */
	for (i = *ndata; i >= 2; i--) {
		a[i - 2] += a[i - 1];
	}
	for (i = 1; i <= (*ndata - *korder); i++) {
		d[i - 1] = a[i - 1] * (xdata[i + *korder - 1] - xdata[i - 1]) /
			(Mfloat) (*korder);
	}
	/* Compute X by solving linear system */
	nbands = imsl_i_min(*ndata - *korder - 1, *korder - 1);
        _l0 = *ndata - *korder;
        _l1 = 3 ** korder - 2;
        _l2 = 3 ** korder - 2;
	imsl_l2crb(&_l0, w, &_l1, &nbands,
		 &nbands, w, &_l2, ipvt, &rcond, work);
	if (imsl_n1rty(1) != 0) {

/*		imsl_ermes(4, 4, "The algorithm has detected ill-conditioning in a linear system and can not proceed.  Higher precision may help.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_ILL_COND_LIN_SYS);
		goto L_9000;
	}
        _l0 = *ndata - *korder;
        _l1 = 3 ** korder - 2;
        _l2 = 1;
	imsl_lfsrb(&_l0, w, &_l1, &nbands,
		   &nbands, ipvt, d, &_l2, d);
	/*
	 * Compute D = Change in XI.  Modify, IF necessary, to prevent new XI
	 * from moving more than 1/3 of the way to its neighbors. Then add to
	 * XI to obtain new XI.
	 */
	delmax = F_ZERO;
	siign = siignst;
	for (i = 1; i <= (*ndata - *korder); i++) {
		del = siign * d[i - 1];
		delmax = imsl_f_max(delmax, fabs(del));
		if (del > F_ZERO) {
			del = imsl_f_min(del, (xi[i + 1] - xi[i]) / F_THREE);
		} else {
			del = imsl_f_max(del, (xi[i - 1] - xi[i]) / F_THREE);
		}
		siign = -siign;
		xi[i] += del;
	}
	/*
	 * Call it a day in case change in XI was small enough or too many
	 * steps were taken.
	 */
	if (delmax < tol)
		goto L_60;
	newton += 1;
	if (newton <= newtmx) {
		goto L_20;
	} else {
		imsl_e1sti(1, newtmx);

/*		imsl_ermes(3, 6, "No convergence after %(i1) Newton iterations.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_NO_CONV_NEWTON);
	}
L_60:
	;
	/* Set output knot sequence XKNOT */
L_70:
	sset(*korder, xdata[0], xknot, 1);
	scopy(*ndata - *korder, &xi[1], 1, &xknot[*korder], 1);
	/*
	 * Move the last endpoint slightly to the right of the last data
	 * point.
	 */
	eps = 100.0e0 * imsl_amach(4);
L_80:
	if (xdata[*ndata - 1] + eps <= xdata[*ndata - 1]) {
		eps *= F_TEN;
		goto L_80;
	}
	sset(*korder, xdata[*ndata - 1] + eps, &xknot[*ndata], 1);

L_9000:
	;
	imsl_e1pop("IMSL_B3OPK ");
	return;
}				/* end of function */
/*  -----------------------------------------------------------------------
    IMSL Name:  B4OPK/DB4OPK (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 7, 1986

    Purpose:    Compute the 'optimal' B-spline knot sequence.

    Usage:      CALL B4OPK (NDATA, KORDER, SIGNST, EXDATA, XI, A,
                            BIATX, DL, DR, W, LDW)

    Arguments:
       NDATA  - Number of data points.  (Input)
       KORDER - Order of the B-spline.  (Input)
       SIGNST - Flag that is -1 if NDATA-KORDER is even and 1 otherwise.
                (Input)
       EXDATA - Array of length NDATA+2*KORDER for the extended knot
                sequence.  (Workspace)
       XI     - Array of length NDATA-KORDER+2 for the intermediate
                knot locations.  (Workspace)
       A      - Array of length NDATA.  (Workspace)
       BIATX  - Array of length NDATA-KORDER to hold B-spline values.
                (Workspace)
       DL     - Array of length KORDER+1 used by B4INT.  (Workspace)
       DR     - Array of length KORDER+1 used by B4INT.  (Workspace)
       W      - Array of length (NDATA-KORDER)*(3*KORDER-2) to hold the
                band matrix of order NDATA-KORDER with KORDER-1 upper
                and lower codiagonals.  (Workspace)
       LDW    - Leading dimension of W = 3*KORDER-2.  (Input)
                Note that the last KORDER-1 columns of W are workspace.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_b4opk(Mint *ndata, Mint *korder, Mfloat *siignst, Mfloat exdata[],
                Mfloat xi[], Mfloat a[], Mfloat biatx[],
	        Mfloat dl[], Mfloat dr[], Mfloat *w, Mint *ldw)
#else
void imsl_b4opk(ndata, korder, siignst, exdata, xi, a, biatx,
	   dl, dr, w, ldw)
	Mint            *ndata, *korder;
	Mfloat          *siignst, exdata[], xi[], a[], biatx[], dl[], dr[],
	               *w;
	Mint            *ldw;
#endif
{
#define W(I_,J_)	(w+(I_)*(*ldw)+(J_))
	Mint             _l0, id, inc, j, jj, left, ll, nbands;
	Mfloat           siign;

	/* Zero out W and A(1:NDATA-1) */
	sset(*ldw * (*ndata - *korder), F_ZERO, w, 1);
	sset(*ndata - 1, F_ZERO, a, 1);

	siign = *siignst;
	left = *korder + 1;
	for (j = 1; j <= (*ndata - *korder); j++) {
L_10:
		if (xi[j] < exdata[left])
			goto L_20;
		left += 1;
		if (left < *korder + *ndata)
			goto L_10;
		left -= 1;
		/* Compute B-splines at XI(J) */
L_20:
		imsl_b4int(exdata, korder, &xi[j], &left, biatx, dl, dr);
		/*
		 * The knot sequence in EXDATA is preceded by KORDER
		 * additional knots therefore, BIATX(LL) now contains
		 * B(LEFT-2K+LL)(XI(J)) which is destined for
		 * C(LEFT-2K+LL,J). Also, C being of order NDATA-KORDER, we
		 * would want 1 .LE. LEFT-2K+LL .LE. N-K OR 2K-LEFT+1 .LE. LL
		 * .LE. N-LEFT+K
		 */
		for (ll = imsl_i_max(1, 2 ** korder - left); ll <= imsl_i_min(*korder,
					   *ndata - left + *korder); ll++) {
			*W(j - 1, left - *korder + ll - j - 1) = biatx[ll - 1];
		}
		/*
		 * Compute B-splines of order KORDER+1 for A.
		 */
                _l0 = *korder + 1;
		imsl_b4int(exdata, &_l0, &xi[j], &left, biatx,
			   dl, dr);
		id = imsl_i_max(0, left - 2 ** korder - 1);
		for (ll = 1 - imsl_i_min(0, left - 2 ** korder - 1); ll <= (*korder +
								 1); ll++) {
			id += 1;
			a[id - 1] += -siign * biatx[ll - 1];
		}
		siign = -siign;
	}

	/*
	 * If there are not K-1 bands then we must move the rows of W up
	 * accordingly.
	 */
	nbands = imsl_i_min(*korder - 1, *ndata - *korder - 1);
	if (nbands != *korder - 1) {
		inc = *korder - 1 - nbands;
		for (jj = 1; jj <= (2 * nbands + 1); jj++) {
			scopy(*ndata - *korder, W(0, jj + inc - 1), *ldw, W(0, jj - 1),
				   *ldw);
		}
	}
	return;
}				/* end of function */
/*  -----------------------------------------------------------------------
    IMSL Name:  SSUM (Single precision version)
 
    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Sum the values of a single precision vector.

    Usage:      SSUM(N, SX, INCX)

    Arguments:
       N      - Length of vectors X.  (Input)
       SX     - Real vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
                greater than 0.
       SSUM   - Single precision sum from I=1 to N of X(I).  (Output)
                X(I) refers to a specific element of SX.

    Keyword:    Level 1 BLAS

    GAMS:       D1a

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_ssum(Mint *n, Mfloat sx[], Mint *incx)
#else
static Mfloat l_ssum(n, sx, incx)
        Mint            *n;
        Mfloat          *sx;
        Mint            *incx;
#endif
{
        Mint            _d_l, _d_m, _do0, _do1, i, nincx;
        Mfloat          ssum_v;


        ssum_v = F_ZERO;
        if (*n > 0) {
                if (*incx != 1) {
                        /* CODE FOR INCREMENT NOT EQUAL TO 1 */
                        nincx = *n * *incx;
                        for (i = 1, _do0 = DOCNT(1, nincx, _do1 = *incx); _do0 > 0; i += _do1, _do0--) {
                                ssum_v += sx[i - 1];
                        }
                } else {
                        for (i = 1; i <= *n; i++) {
                                ssum_v += sx[i - 1];
                        }
                }
        }
        return (ssum_v);
}  
