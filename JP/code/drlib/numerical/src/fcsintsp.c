#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static Mf_ppoly *pp = NULL;

static void PROTO (l_c2akm, (Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat break_[],
	            Mfloat *cscoef, Mint ipvt[]));
    static VA_LIST_HACK PROTO (l_cub_spline_interp_shape, (Mint ndata, Mfloat *xdata,
	            Mfloat *fdata, va_list argptr));
    static void PROTO (l_c2con, (Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mint *ibreak,
	            Mfloat break_[], Mfloat *cscoef,
	            Mint *itmax, Mfloat xsrt[], Mfloat fsrt[], Mfloat a[], Mfloat y[],
	            Mfloat divd[], Mint id[], Mfloat wk[]));
    static void PROTO (l_c4con, (Mint *nm2, Mfloat a[], Mfloat x[], Mint *itmax, Mfloat *eps,
	            Mint *iflag, Mfloat *tl, Mfloat *tr, Mfloat divd[],
	            Mint id[], Mfloat sup[], Mfloat diag[], Mfloat sub[],
	            Mfloat h[], Mfloat xold[]));
    static void PROTO (l_c5con, (Mfloat a[], Mfloat xdata[], Mfloat *cscoef, Mint *ndata,
	            Mfloat fdata[], Mint *ibreak, Mfloat break_[],
	            Mfloat divd[], Mint id[]));
    static void PROTO (l_c6con, (Mfloat sub[], Mfloat diag[], Mfloat sup[], Mfloat b[], Mint *n));


#ifdef ANSI
    Mf_ppoly   *imsl_f_cub_spline_interp_shape (Mint ndata, Mfloat xdata[],
                Mfloat fdata[],...)
#else
    Mf_ppoly   *imsl_f_cub_spline_interp_shape (ndata, xdata, fdata, va_alist)
    Mint        ndata;
    Mfloat      xdata[];
    Mfloat      fdata[];
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, fdata);
#ifdef DOUBLE
    imsl_e1psh ("imsl_d_cub_spline_interp_shape");
#else
    imsl_e1psh ("imsl_f_cub_spline_interp_shape");
#endif
    pp = NULL;
    IMSL_CALL (l_cub_spline_interp_shape (ndata, xdata, fdata, argptr));
    va_end (argptr);
#ifdef DOUBLE
    imsl_e1pop ("imsl_d_cub_spline_interp_shape");
#else
    imsl_e1pop ("imsl_f_cub_spline_interp_shape");
#endif
    return pp;
}


#ifdef ANSI
static VA_LIST_HACK l_cub_spline_interp_shape (Mint ndata, Mfloat *xdata,
                Mfloat *fdata, va_list argptr)
#else
static VA_LIST_HACK l_cub_spline_interp_shape (ndata, xdata, fdata, argptr)
    Mint        ndata;
    Mfloat      xdata[];
    Mfloat      fdata[];
    va_list     argptr;
#endif
{
    Mint        arg_number = 3;
    Mint        spline_type = 1;/* DEFAULT TO AKIMA SPLINE */
    Mint       *orders = NULL;
    Mint       *num_breakpoints = NULL;
    Mint        four = 4;
    Mint        tmp;
    Mint        code;
    Mint       *iwork = NULL;
    Mfloat     *coef_work = NULL;
    Mint        i;
    Mint        itmax = 25;
    Mint        domain_dim;
    Mint        target_dim;
    Mfloat     *xsrt = NULL;
    Mfloat     *fsrt = NULL;
    Mfloat     *a = NULL;
    Mfloat     *y = NULL;
    Mfloat     *divd = NULL;
    Mint       *id = NULL;
    Mfloat     *wk = NULL;
    Mint        free_the_structure = 0;
    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_CONCAVE:
	    spline_type = 2;
	    break;
	case IMSL_CONCAVE_ITMAX:
	    itmax = va_arg (argptr, Mint);
	    arg_number++;
	case 0:
	    break;
	default:
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, arg_number);
	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    goto RETURN;
	}
    }
    /****************** AKIMA SPLINE *********************/
    if (spline_type == 1) {
	/* CHECK ARGUMENT NDATA */
	if (ndata < 4) {
	    imsl_e1sti (1, ndata);

/*              imsl_ermes(5, 1, "The number of data points must be 4 or more while NDATA = %(i1) is given.");
*/
	    imsl_ermes (IMSL_TERMINAL, IMSL_NEED_AT_LEAST_4_PTS);
	    goto RETURN;
	}
	domain_dim = 1;
	target_dim = 1;
	orders = &four;
	tmp = ndata;
	num_breakpoints = &tmp;
	/* CREATE THE STRUCTURE */
	pp = imsl_f_ppoly_create (domain_dim, target_dim, orders, num_breakpoints, 0);
	if (imsl_n1rty (1) == 4) {
	    imsl_e1mes (0, 0, " ");
	    imsl_e1stl (1, "ndata");
	    imsl_e1sti (1, ndata);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto RETURN;
	}
	/* GET THE WORKSPACE    */
	/*
	 * NEED TO USE WORKSPACE FOR THE COEFFICIENTS SINCE THE SPACE NEEDED
	 * BY C2INT IS LARGER THAN THE SPACE PROVIDED IN THE STRUCTURE
	 * IMSL_PPOLY. THIS IS CAUSED BY C2INT USING THE LAST 4 MEMORY SPACES
	 * OF CSCOEF FOR WORKSPACE
	 */
	coef_work = (Mfloat *) imsl_malloc (4 * ndata * sizeof (*coef_work));
	iwork = (Mint *) imsl_malloc (ndata * sizeof (*iwork));
	if ((coef_work == NULL) || (iwork == NULL)) {
	    free_the_structure = 1;
	    imsl_e1stl (1, "ndata");
	    imsl_e1sti (1, ndata);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto AKIMA_FREE_SPACE;
	}
	/* CALL THE SPLINE ROUTINE */
	l_c2akm (&ndata, xdata, fdata, pp->breakpoints[0], coef_work, iwork);
	if (imsl_n1rty (1) > 3) {
	    free_the_structure = 1;
	    goto AKIMA_FREE_SPACE;
	}

	/* COPY THE COEFFICIENTS INTO THE STRUCTURE */
	for (i = 0; i < pp->num_coef[0]; i++)
	    pp->coef[0][i] = coef_work[i];
	goto AKIMA_FREE_SPACE;
    }
    /****************** CONCAVE SPLINE *********************/
    else if (spline_type == 2) {
	if (ndata <= 2) {
	    imsl_e1sti (1, ndata);

/*              imsl_ermes(5, 1, "The number of data points must be at least 3 while NDATA = %(i1) is given.");
*/
	    imsl_ermes (IMSL_TERMINAL, IMSL_NEED_AT_LEAST_3_PTS);
	    goto RETURN;
	}
	domain_dim = 1;
	target_dim = 1;
	orders = &four;
	tmp = 2 * ndata;
	num_breakpoints = &tmp;
	/* CREATE THE STRUCTURE */
	pp = imsl_f_ppoly_create (domain_dim, target_dim, orders, num_breakpoints, 0);
	if (imsl_n1rty (1) == 4) {
	    imsl_e1mes (0, 0, " ");
	    imsl_e1stl (1, "ndata");
	    imsl_e1sti (1, ndata);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto RETURN;
	}
	/* GET THE WORKSPACE    */
	/*
	 * NEED TO USE WORKSPACE FOR THE COEFFICIENTS SINCE THE SPACE NEEDED
	 * BY C2INT IS LARGER THAN THE SPACE PROVIDED IN THE STRUCTURE
	 * IMSL_PPOLY. THIS IS CAUSED BY C2INT USING THE LAST 4 MEMORY SPACES
	 * OF CSCOEF FOR WORKSPACE
	 */
	coef_work = (Mfloat *) imsl_malloc (2 * 4 * ndata * sizeof (*coef_work));
	xsrt = (Mfloat *) imsl_malloc (ndata * sizeof (*xsrt));
	fsrt = (Mfloat *) imsl_malloc (ndata * sizeof (*fsrt));
	a = (Mfloat *) imsl_malloc (ndata * sizeof (*a));
	y = (Mfloat *) imsl_malloc ((ndata - 2) * sizeof (*y));
	divd = (Mfloat *) imsl_malloc ((ndata - 2) * sizeof (*divd));
	id = (Mint *) imsl_malloc (ndata * sizeof (*id));
	wk = (Mfloat *) imsl_malloc ((5 * (ndata - 2)) * sizeof (*wk));
	if ((coef_work == NULL) || (xsrt == NULL) || (fsrt == NULL) || (a == NULL) || (y == NULL) || (divd == NULL) || (id == NULL) || (wk == NULL)) {
	    free_the_structure = 1;
	    imsl_e1stl (1, "ndata");
	    imsl_e1sti (1, ndata);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto CONCAVE_FREE_SPACE;
	}

	/* CALL THE SPLINE ROUTINE */
	l_c2con (&ndata, xdata, fdata, pp->num_breakpoints, pp->breakpoints[0], coef_work, &itmax, xsrt, fsrt, a, y, divd, id, wk);
	if (imsl_n1rty (1) > 3) {
	    free_the_structure = 1;
	    goto CONCAVE_FREE_SPACE;
	}
	pp->num_coef[0] = ((pp->num_breakpoints[0]) - 1) * (pp->order[0]);
	/* COPY THE COEFFICIENTS INTO THE STRUCTURE */
	for (i = 0; i < pp->num_coef[0]; i++)
	    pp->coef[0][i] = coef_work[i];
	goto CONCAVE_FREE_SPACE;
    }
    /* FREE THE WORKSPACE USED */
AKIMA_FREE_SPACE:
    if (free_the_structure == 1) {
	if (pp != NULL)
	    imsl_free (pp);
	pp = NULL;
    }
    if (iwork != NULL)
	imsl_free (iwork);
    if (coef_work != NULL)
	imsl_free (coef_work);
    goto RETURN;
CONCAVE_FREE_SPACE:
    if (free_the_structure == 1) {
	if (pp != NULL)
	    imsl_free (pp);
	pp = NULL;
    }
    if (coef_work != NULL)
	imsl_free (coef_work);
    if (xsrt != NULL)
	imsl_free (xsrt);
    if (fsrt != NULL)
	imsl_free (fsrt);
    if (a != NULL)
	imsl_free (a);
    if (y != NULL)
	imsl_free (y);
    if (divd != NULL)
	imsl_free (divd);
    if (id != NULL)
	imsl_free (id);
    if (wk != NULL)
	imsl_free (wk);
    goto RETURN;
RETURN:
    return argptr;
}






/*  -----------------------------------------------------------------------
    IMSL Name:  C2CON/DC2CON (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 26, 1985

    Purpose:    Compute a cubic spline interpolant which is consistent
                with the concavity of the data.

    Usage:      CALL C2CON (NDATA, XDATA, FDATA, IBREAK, BREAK, CSCOEF,
                            ITMAX, XSRT, FSRT, A, Y, DIVD, ID, WK)

    Arguments:  (See CSCON)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c2con (Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mint *ibreak, Mfloat break_[], Mfloat *cscoef,
                Mint *itmax, Mfloat xsrt[], Mfloat fsrt[], Mfloat a[], Mfloat y[], Mfloat divd[], Mint id[], Mfloat wk[])
#else
static void l_c2con (ndata, xdata, fdata, ibreak, break_, cscoef,
                itmax, xsrt, fsrt, a, y, divd, id, wk)
    Mint       *ndata;
    Mfloat      xdata[], fdata[];
    Mint       *ibreak;
    Mfloat      break_[], *cscoef;
    Mint       *itmax;
    Mfloat      xsrt[], fsrt[], a[], y[], divd[];
    Mint        id[];
    Mfloat      wk[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
    Mint        i, iflag = 0, ittmp, k, n;
    Mfloat      eps, tl, tr;


    imsl_e1psh ("C2CON");

    if (*ndata <= 2) {
	imsl_e1sti (1, *ndata);

/*		imsl_ermes(5, 1, "The number of data points must be at least 3 while NDATA = %(i1) is given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_AT_LEAST_3_PTS);
	goto L_9000;
    }
    if (*itmax <= 1) {
	imsl_e1sti (1, *itmax);

/*		imsl_ermes(5, 2, "The maximum number of iterations must be at least 2 while ITMAX = %(i1) is given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_LARGER_ITMAX);
	goto L_9000;
    }
    /*
     * Sort XDATA,FDATA and check XDATA for uniqueness.
     */
    for (i = 1; i <= *ndata; i++) {
	id[i - 1] = i;
    }
    imsl_svrgp (*ndata, xdata, xsrt, id);
    fsrt[0] = fdata[id[0] - 1];
    for (i = 2; i <= *ndata; i++) {
	fsrt[i - 1] = fdata[id[i - 1] - 1];
	if (xsrt[i - 2] == xsrt[i - 1]) {
	    imsl_e1sti (1, id[i - 2] - 1);
	    imsl_e1sti (2, id[i - 1] - 1);
	    imsl_e1str (1, xsrt[i - 1]);

/*			imsl_ermes(4, 3, "Points in the data point abscissas array, XDATA, must be distinct, but XDATA(%(i1)) = XDATA(%(i2)) = %(r1).");*/
	    imsl_ermes (IMSL_FATAL, IMSL_DUPLICATE_XDATA_VALUES);
	    goto L_9000;
	}
    }

    n = *ndata - 2;
    /*
     * (EPS) is a small positive number used to test for convergence in
     * Newton'S method - SUBROUTINE C4CON. EPS = MACHINE EPSILON*10
     */
    eps = imsl_amach (4) * F_TEN;
    /*
     * X is the knot sequence, XSRT, with the endpoints TL and TR deleted.
     */
    tl = xsrt[0];
    tr = xsrt[*ndata - 1];
    /*
     * DIVD consists of the scaled second divided differences.
     */
    for (k = 1; k <= n; k++) {
	divd[k - 1] = (fsrt[k + 1] - fsrt[k]) / (xsrt[k + 1] - xsrt[k]) -
	    (fsrt[k] - fsrt[k - 1]) / (xsrt[k] - xsrt[k - 1]);
	if (fabs (divd[k - 1]) <= eps)
	    divd[k - 1] = F_ZERO;
    }
    /*
     * The initial guess (Y) for Newton'S method will yield the second
     * derivative of the natural spline solution, except possibly when
     * DIVD(K)= 0.0 for some K.
     */
    for (k = 1; k <= n; k++) {
	y[k - 1] = -sign (F_ONE, -divd[k - 1]);
    }
    /*
     * ID(K)= 1 indicates that the interpolating function is constrained to
     * be convex on the closed interval (XSRT(K),XSRT(K+1)) and, hence, its
     * second derivative is constrained to be nonpositive on this interval.
     * ID(K)= -1 indicates that the interpolating function is constrained to
     * be concave on the closed interval (XSRT(K),XSRT(K+1)) and, hence, its
     * second derivative is constrained to be nonpositive on this interval.
     * ID(K)= 0 indicates that the interpolating function is unconstrained on
     * the closed interval (XSRT(K),XSRT(K+1)).
     */
    for (i = 1; i <= (n - 1); i++) {
	id[i] = 0;
	if (divd[i - 1] >= F_ZERO && divd[i] >= F_ZERO)
	    id[i] = 1;
	if (divd[i - 1] <= F_ZERO && divd[i] <= F_ZERO)
	    id[i] = -1;
    }
    if (divd[0] >= F_ZERO) {
	id[0] = 1;
    }
    else {
	id[0] = -1;
    }

    if (divd[n - 1] >= F_ZERO) {
	id[n] = 1;
    }
    else {
	id[n] = -1;
    }
    /*
     * If a nonzero data value D(I) lies between two zero data values
     * DIVD(I-1) and DIVD(I+1), then DIVD(I) is taken to be zero for
     * computational purposes.
     */
    for (i = 2; i <= (n - 1); i++) {
	if (divd[i - 2] == F_ZERO && divd[i] == F_ZERO)
	    divd[i - 1] = F_ZERO;
    }
    /*
     * SUBROUTINE C4CON calculates the piecewise linear second derivative of
     * the shape- preserving interpolant.
     */
    ittmp = *itmax;
    l_c4con (&n, y, &xsrt[1], &ittmp, &eps, &iflag, &tl, &tr, divd,
	id, &wk[0], &wk[n], &wk[n * 2], &wk[n * 3], &wk[n * 4]);
    if (iflag == 2) {
	imsl_e1sti (1, *itmax);

/*		imsl_ermes(3, 16, "The maximum number of iterations was reached in Newtons Method.  The best answer will be returned.  ITMAX = %(i1) was used, a larger value may help.");
*/
	imsl_ermes (IMSL_WARNING, IMSL_MAX_ITERATIONS_REACHED);
    }
    a[0] = F_ZERO;
    a[*ndata - 1] = F_ZERO;
    scopy (n, y, 1, &a[1], 1);
    /*
     * SUBROUTINE C5CON integrates the result from SUBROUTINE C4CON.
     */
    l_c5con (a, xsrt, cscoef, ndata, fsrt, ibreak, break_, divd, id);

L_9000:
    imsl_e1pop ("C2CON");
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  C4CON/DC4CON (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 26, 1985

    Purpose:    Compute a cubic spline interpolant which is
                consistent with the concavity of the data.

    Usage:      CALL C4CON (NM2, A, X, ITMAX, EPS, IFLAG, TL, TR, DIVD,
                            ID, SUP, DIAG, SUB, H, XOLD)

    Arguments:
       NM2    - The size of the A and X.  (Input)
                NM2 is equal to NDATA-2.
       A      - Vector of length NM2 containing the initial estimate
                for Newton's method.  (Input/Output)
                On output A contains the calculated zero from Newton's
                Method.
       X      - Vector of length NDATA-2 containing the knot sequence
                with the endpoints deleted.  (Input)
       ITMAX  - Number of iterations for Newton's method.  (Input/Output)
                On input, ITMAX contains the maximum number of iterations
                for Newton's Method.  On output, ITMAX is set to the
                actual number of iterations that occured.
       EPS    - Value used to test for convergence of Newton's method.
                (Input)
       IFLAG  - Flag to indicate convergence/non-convergence.  (Ouput)
                IFLAG = 1: Convergence indicated by comparing the L1
                norms of the iterates.
                IFLAG = 2: Number of iterations exceeded ITMAX.
       TL     - Left endpoint of the interval.  (Input)
       TR     - Right endpoint of the interval.  (Input)
       DIVD   - Vector of divided differences.  (Input)
       ID     - Integer vector of length NM2.  (Workspace)
       SUP    - Holds the super diagonal of the matrix.  (Workspace)
       DIAG   - Holds the main diagonal of the matrix.  (Workspace)
       SUB    - Holds the sub-diagonal of the matrix.  (Workspace)
       H      - Temporary solution vector.  (Workspace)
       XOLD   - Vector that holds the previous solution vector.
                (Workspace)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c4con (Mint *nm2, Mfloat a[], Mfloat x[], Mint *itmax, Mfloat *eps,
                Mint *iflag, Mfloat *tl, Mfloat *tr, Mfloat divd[],
                Mint id[], Mfloat sup[], Mfloat diag[], Mfloat sub[], Mfloat h[], Mfloat xold[])
#else
static void l_c4con (nm2, a, x, itmax, eps, iflag, tl, tr, divd,
                id, sup, diag, sub, h, xold)
    Mint       *nm2;
    Mfloat      a[], x[];
    Mint       *itmax;
    Mfloat     *eps;
    Mint       *iflag;
    Mfloat     *tl, *tr, divd[];
    Mint        id[];
    Mfloat      sup[], diag[], sub[], h[], xold[];
#endif
{
    Mint        i, ii, j1, j2, k, lj, n;
    Mfloat      al, ar, big, da, dt, gleft, grigh, r1nrm, small, t, w,
                xl, xr;


    imsl_e1psh ("C4CON");
    small = imsl_amach (1);
    big = imsl_amach (2);
    if (small * big < F_ONE)
	small = F_ONE / big;
    n = *nm2;
    small = sqrt (small);
    for (lj = 1; lj <= *itmax; lj++) {
	/*
	 * The arrays (SUB), (DIAG), and (SUP) contain the elements of the
	 * tridiagonal positive-definite Jacobian matrix (J), evaluated at
	 * the vector (A). It should be noted that the matrix equation
	 * solver, the SUBROUTINE (TRID), does not take advantage of the
	 * symmetry of (J). Hence (SUB) and (SUP) are both necessary,
	 * although SUB(K)=SUP(K-1), equations for both arrays are written
	 * out in full. IF DIVD(K)=0.0 for some K, then the number of
	 * unknowns (and equations) reduce. In order to permit the
	 * computation of one Jacobian matrix the program sets
	 * SUB(K)=SUP(K-1)=0.0 and DIAG(K)=1.0.
	 */
	for (k = 1; k <= n; k++) {

	    if (k == 1) {
		al = F_ZERO;
		xl = *tl;
	    }
	    else {
		al = a[k - 2];
		xl = x[k - 2];
	    }

	    if (k == n) {
		ar = F_ZERO;
		xr = *tr;
	    }
	    else {
		ar = a[k];
		xr = x[k];
	    }

	    if (al >= F_ZERO && a[k - 1] >= F_ZERO)
		j1 = 1;
	    if (al < F_ZERO && a[k - 1] >= F_ZERO)
		j1 = 2;
	    if (al >= F_ZERO && a[k - 1] < F_ZERO)
		j1 = 3;
	    if (al <= F_ZERO && a[k - 1] <= F_ZERO)
		j1 = 4;

	    if (a[k - 1] >= F_ZERO && ar >= F_ZERO)
		j2 = 1;
	    if (a[k - 1] < F_ZERO && ar >= F_ZERO)
		j2 = 2;
	    if (a[k - 1] >= F_ZERO && ar < F_ZERO)
		j2 = 3;
	    if (a[k - 1] <= F_ZERO && ar <= F_ZERO)
		j2 = 4;

	    dt = x[k - 1] - xl;
	    da = a[k - 1] - al;

	    if (id[k - 1] == 1) {

		if (k != 1) {

		    if (j1 == 1) {

			sub[k - 1] = dt / F_SIX;
			gleft = dt / F_THREE;

		    }
		    else if (j1 == 2) {

			t = xl - (dt / da) * al;
			w = F_HALF * (x[k - 1] + t);
			sub[k - 1] = (x[k - 1] - t) / F_SIX * (((t - xl) /
				dt) * ((x[k - 1] - t) / dt) + F_FOUR * ((w - xl) / dt) *
			    ((x[k - 1] - w) / dt));
			gleft = (x[k - 1] - t) / F_SIX * (imsl_fi_power ((t - xl) /
				dt, 2) + F_FOUR * (imsl_fi_power ((w - xl) / dt, 2)) + F_ONE);

		    }
		    else if (j1 == 3) {

			t = xl - (dt / da) * al;
			w = F_HALF * (t + xl);
			sub[k - 1] = (t - xl) / F_SIX * (F_FOUR * ((w - xl) / dt) *
			    ((x[k - 1] - w) / dt) + ((t - xl) / dt) * ((x[k - 1] -
				    t) / dt));
			gleft = (t - xl) / F_SIX * (F_FOUR * (imsl_fi_power ((w - xl) / dt, 2)) +
			    imsl_fi_power ((t - xl) / dt, 2));

		    }
		    else if (j1 == 4) {

			sub[k - 1] = F_ZERO;
			gleft = F_ZERO;

		    }
		}
		else if (k == 1) {

		    sub[0] = F_ZERO;
		    gleft = F_ZERO;
		    if (j1 == 1)
			gleft = dt / F_THREE;

		}
	    }
	    else if (id[k - 1] == 0) {

		sub[k - 1] = dt / F_SIX;
		gleft = dt / F_THREE;

	    }
	    else if (id[k - 1] == -1) {

		if (k != 1) {

		    if (j1 == 4) {

			sub[k - 1] = dt / F_SIX;
			gleft = dt / F_THREE;

		    }
		    else if (j1 == 3) {

			t = xl - (dt / da) * al;
			w = F_HALF * (x[k - 1] + t);
			sub[k - 1] = (x[k - 1] - t) / F_SIX * (((t - xl) /
				dt) * ((x[k - 1] - t) / dt) + F_FOUR * ((w - xl) / dt) *
			    ((x[k - 1] - w) / dt));
			gleft = (x[k - 1] - t) / F_SIX * (imsl_fi_power ((t - xl) /
				dt, 2) + F_FOUR * (imsl_fi_power ((w - xl) / dt, 2)) + F_ONE);

		    }
		    else if (j1 == 2) {

			t = xl - (dt / da) * al;
			w = F_HALF * (t + xl);
			sub[k - 1] = (t - xl) / F_SIX * (F_FOUR * ((w - xl) / dt) *
			    ((x[k - 1] - w) / dt) + ((t - xl) / dt) * ((x[k - 1] -
				    t) / dt));
			gleft = (t - xl) / F_SIX * (F_FOUR * (imsl_fi_power ((w - xl) / dt, 2)) +
			    imsl_fi_power ((t - xl) / dt, 2));

		    }
		    else if (j1 == 1) {

			sub[k - 1] = F_ZERO;
			gleft = F_ZERO;

		    }
		}
		else if (k == 1) {

		    sub[0] = F_ZERO;
		    gleft = F_ZERO;
		    if (j1 == 4)
			gleft = dt / F_THREE;

		}
	    }
	    if (k != 1) {
		if (divd[k - 2] == F_ZERO) {
		    sub[k - 1] = F_ZERO;
		    gleft = F_ZERO;
		}
	    }
	    dt = xr - x[k - 1];
	    da = ar - a[k - 1];

	    if (id[k] == 1) {

		if (k != n) {

		    if (j2 == 1) {

			sup[k - 1] = dt / F_SIX;
			grigh = dt / F_THREE;

		    }
		    else if (j2 == 2) {

			t = x[k - 1] - (dt / da) * a[k - 1];
			w = F_HALF * (xr + t);
			sup[k - 1] = (xr - t) / F_SIX * (((t - x[k - 1]) /
				dt) * ((xr - t) / dt) + F_FOUR * ((w - x[k - 1]) / dt) *
			    ((xr - w) / dt));
			grigh = (xr - t) / F_SIX * (imsl_fi_power ((xr - t) / dt, 2) +
			    F_FOUR * (imsl_fi_power ((xr - w) / dt, 2)));

		    }
		    else if (j2 == 3) {

			t = x[k - 1] - (dt / da) * a[k - 1];
			w = F_HALF * (t + x[k - 1]);
			sup[k - 1] = (t - x[k - 1]) / F_SIX * (F_FOUR * ((w -
				    x[k - 1]) / dt) * ((xr - w) / dt) + ((t - x[k - 1]) /
				dt) * ((xr - t) / dt));
			grigh = (t - x[k - 1]) / F_SIX * (F_ONE + F_FOUR * (imsl_fi_power ((xr -
					w) / dt, 2)) + imsl_fi_power ((xr - t) / dt, 2));

		    }
		    else if (j2 == 4) {

			sup[k - 1] = F_ZERO;
			grigh = F_ZERO;

		    }
		}
		else if (k == n) {

		    sup[n - 1] = F_ZERO;
		    grigh = F_ZERO;
		    if (j2 == 1)
			grigh = dt / F_THREE;

		}
	    }
	    else if (id[k] == 0) {

		sup[k - 1] = dt / F_SIX;
		grigh = dt / F_THREE;

	    }
	    else if (id[k] == -1) {

		if (k != n) {

		    if (j2 == 4) {

			sup[k - 1] = dt / F_SIX;
			grigh = dt / F_THREE;

		    }
		    else if (j2 == 3) {

			t = x[k - 1] - (dt / da) * a[k - 1];
			w = F_HALF * (xr + t);
			sup[k - 1] = (xr - t) / F_SIX * (((t - x[k - 1]) /
				dt) * ((xr - t) / dt) + F_FOUR * ((w - x[k - 1]) / dt) *
			    ((xr - w) / dt));
			grigh = (xr - t) / F_SIX * (imsl_fi_power ((xr - t) / dt, 2) +
			    F_FOUR * (imsl_fi_power ((xr - w) / dt, 2)));

		    }
		    else if (j2 == 2) {

			t = x[k - 1] - (dt / da) * a[k - 1];
			w = F_HALF * (t + x[k - 1]);
			sup[k - 1] = (t - x[k - 1]) / F_SIX * (F_FOUR * ((w -
				    x[k - 1]) / dt) * ((xr - w) / dt) + ((t - x[k - 1]) /
				dt) * ((xr - t) / dt));
			grigh = (t - x[k - 1]) / F_SIX * (F_ONE + F_FOUR * (imsl_fi_power ((xr -
					w) / dt, 2)) + imsl_fi_power ((xr - t) / dt, 2));

		    }
		    else if (j2 == 1) {

			sup[k - 1] = F_ZERO;
			grigh = F_ZERO;

		    }
		}
		else if (k == n) {

		    sup[n - 1] = F_ZERO;
		    grigh = F_ZERO;
		    if (j2 == 4)
			grigh = dt / F_THREE;

		}
	    }
	    if (k != n) {
		if (divd[k] == F_ZERO) {
		    sup[k - 1] = F_ZERO;
		    grigh = F_ZERO;
		}
	    }
	    diag[k - 1] = gleft + grigh;

	    if (divd[k - 1] == F_ZERO) {
		diag[k - 1] = F_ONE;
		sub[k - 1] = F_ZERO;
		sup[k - 1] = F_ZERO;
	    }
	}

	scopy (n, divd, 1, h, 1);
	/*
	 * We solve the matrix equation JX=H, The array (H) being identical
	 * to the array (DIVD). The solution is returned in the array (H).
	 * 
	 * First we check DIAG for close to values.
	 */
	for (ii = 1; ii <= n; ii++) {
	    if (diag[ii - 1] <= small) {

/*				imsl_ermes(5, 4, "The diagonal of the matrix has a near zero element.  The matrix is ill-conditioned.");
*/
		imsl_ermes (IMSL_TERMINAL,
		    IMSL_DIAG_ELMNT_NEAR_ZERO);
		goto L_9000;
	    }
	}

	l_c6con (sub, diag, sup, h, &n);

	scopy (n, h, 1, a, 1);

	r1nrm = F_ZERO;
	if (lj != 1) {
	    for (i = 1; i <= n; i++) {
		r1nrm += fabs (a[i - 1] - xold[i - 1]);
	    }
	    *iflag = 1;
	    if (r1nrm <= (Mfloat) (n) ** eps)
		goto L_50;
	}
	scopy (n, a, 1, xold, 1);
    }
    *iflag = 2;
L_50:
    ;
    *itmax = lj;
L_9000:
    imsl_e1pop ("C4CON");
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  C5CON/DC5CON (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 26, 1985

    Purpose:    Compute a cubic spline interpolant which is
                consistent with the concavity of the data.

    Usage:      CALL C5CON (A, XDATA, CSCOEF, NDATA, FDATA, IBREAK,
                            BREAK, DIVD, ID)

    Arguments:
       A      - Array of length NDATA-2.  (Workspace)
       XDATA  - Array of length NDATA containing the data point
                abscissas.  (Input)
                The data point abscissas must be distinct.
       CSCOEF - Matrix of size 4 by N where N is the dimension of BREAK.
                (Output)
                The first IBREAK-1 columns of CSCOEF contain the cubic
                spline coefficients.
       NDATA  - Number of data points (must be at least 3).  (Input)
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       IBREAK - The number of breakpoints.  (Output)
                It will be less than 3*NDATA/2.
       BREAK  - Array of length greater than or equal IBREAK containing
                the breakpoints of the piecewise cubic representation
                in its first IBREAK positions.  (Output)
                It must be dimensioned greater than 3*NDATA/2 .
       DIVD   - Array of length NDATA-2 containing the divided
                differences.  (Input)
       ID     - Array of length NDATA-2.  (Workspace)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c5con (Mfloat a[], Mfloat xdata[], Mfloat *cscoef, Mint *ndata, Mfloat fdata[],
                Mint *ibreak, Mfloat break_[], Mfloat divd[], Mint id[])
#else
static void l_c5con (a, xdata, cscoef, ndata, fdata, ibreak, break_,
                divd, id)
    Mfloat      a[], xdata[], *cscoef;
    Mint       *ndata;
    Mfloat      fdata[];
    Mint       *ibreak;
    Mfloat      break_[], divd[];
    Mint        id[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
    Mint        j, jp, l, li, m, mn1;
    Mfloat      c, da, df, dt, e, tau;
    /*
     * SUBROUTINE C5CON integrates back twice the positive part of the
     * piecewise linear second derivative where the data suggests that the
     * interpolating curve should be convex, the negative part of the
     * piecewise linear second derivative where the data suggests that the
     * interpolating curve should be concave, and the remaining portion of
     * the piecewise linear second derivative on the transition intervals.
     */
    li = 1;
    mn1 = *ndata - 1;
    m = *ndata;
    for (l = 1; l <= mn1; l++) {
	df = fdata[l] - fdata[l - 1];
	dt = xdata[l] - xdata[l - 1];
	da = a[l] - a[l - 1];
	jp = 0;
	if (l == 1) {
	    if (divd[0] == F_ZERO)
		jp = 1;
	}
	else if (l == mn1) {
	    if (divd[m - 3] == F_ZERO)
		jp = 1;
	}
	else {
	    c = divd[l - 2] * divd[l - 1];
	    if (c == F_ZERO)
		jp = 1;
	}

	if (jp == 1) {

	    *CSCOEF (li - 1, 3) = F_ZERO;
	    *CSCOEF (li - 1, 2) = F_ZERO;
	    *CSCOEF (li - 1, 1) = df / dt;
	    *CSCOEF (li - 1, 0) = fdata[l - 1];
	    break_[li - 1] = xdata[l - 1];
	    li += 1;

	}
	else if (jp == 0) {

	    if (a[l - 1] >= F_ZERO && a[l] >= F_ZERO)
		j = 1;
	    if (a[l - 1] < F_ZERO && a[l] > F_ZERO)
		j = 2;
	    if (a[l - 1] > F_ZERO && a[l] < F_ZERO)
		j = 3;
	    if (a[l - 1] <= F_ZERO && a[l] <= F_ZERO)
		j = 4;

	    if (id[l - 1] == 1) {

		if (j == 1) {

		    c = df / dt - (da / F_SIX + a[l - 1] / F_TWO) * dt;
		    *CSCOEF (li - 1, 3) = da / dt;
		    *CSCOEF (li - 1, 2) = a[l - 1];
		    *CSCOEF (li - 1, 1) = c;
		    *CSCOEF (li - 1, 0) = fdata[l - 1];
		    break_[li - 1] = xdata[l - 1];
		    li += 1;

		}
		else if (j == 2) {

		    tau = xdata[l - 1] - a[l - 1] * dt / da;
		    c = df / dt - (imsl_fi_power (a[l], 3)) * dt / (F_SIX * da * da);
		    *CSCOEF (li - 1, 3) = F_ZERO;
		    *CSCOEF (li - 1, 2) = F_ZERO;
		    *CSCOEF (li - 1, 1) = c;
		    *CSCOEF (li - 1, 0) = fdata[l - 1];
		    *CSCOEF (li, 3) = da / dt;
		    *CSCOEF (li, 2) = F_ZERO;
		    *CSCOEF (li, 1) = c;
		    *CSCOEF (li, 0) = c * (tau - xdata[l - 1]) + fdata[l - 1];
		    break_[li - 1] = xdata[l - 1];
		    break_[li] = tau;
		    li += 2;

		}
		else if (j == 3) {

		    tau = xdata[l - 1] - a[l - 1] * dt / da;
		    e = fdata[l - 1] - (imsl_fi_power (a[l - 1], 3)) * dt * dt / (F_SIX *
			da * da);
		    c = df / dt + (imsl_fi_power (a[l - 1], 3)) * dt / (F_SIX * da * da);
		    *CSCOEF (li - 1, 3) = da / dt;
		    *CSCOEF (li - 1, 2) = a[l - 1];
		    *CSCOEF (li - 1, 1) = c + a[l - 1] * a[l - 1] * dt * F_HALF /
			da;
		    *CSCOEF (li - 1, 0) = fdata[l - 1];
		    *CSCOEF (li, 3) = F_ZERO;
		    *CSCOEF (li, 2) = F_ZERO;
		    *CSCOEF (li, 1) = c;
		    *CSCOEF (li, 0) = c * (tau - xdata[l - 1]) + e;
		    break_[li - 1] = xdata[l - 1];
		    break_[li] = tau;
		    li += 2;

		}
		else if (j == 4) {

		    *CSCOEF (li - 1, 3) = F_ZERO;
		    *CSCOEF (li - 1, 2) = F_ZERO;
		    *CSCOEF (li - 1, 1) = df / dt;
		    *CSCOEF (li - 1, 0) = fdata[l - 1];
		    break_[li - 1] = xdata[l - 1];
		    li += 1;

		}
	    }
	    else if (id[l - 1] == 0) {

		c = df / dt - (da / F_SIX + a[l - 1] / F_TWO) * dt;
		*CSCOEF (li - 1, 3) = da / dt;
		*CSCOEF (li - 1, 2) = a[l - 1];
		*CSCOEF (li - 1, 1) = c;
		*CSCOEF (li - 1, 0) = fdata[l - 1];
		break_[li - 1] = xdata[l - 1];
		li += 1;

	    }
	    else if (id[l - 1] == -1) {

		if (j == 4) {

		    c = df / dt - (da / F_SIX + a[l - 1] / F_TWO) * dt;
		    *CSCOEF (li - 1, 3) = da / dt;
		    *CSCOEF (li - 1, 2) = a[l - 1];
		    *CSCOEF (li - 1, 1) = c;
		    *CSCOEF (li - 1, 0) = fdata[l - 1];
		    break_[li - 1] = xdata[l - 1];
		    li += 1;

		}
		else if (j == 3) {

		    tau = xdata[l - 1] - a[l - 1] * dt / da;
		    c = df / dt - (imsl_fi_power (a[l], 3)) * dt / (F_SIX * da * da);
		    *CSCOEF (li - 1, 3) = F_ZERO;
		    *CSCOEF (li - 1, 2) = F_ZERO;
		    *CSCOEF (li - 1, 1) = c;
		    *CSCOEF (li - 1, 0) = fdata[l - 1];
		    *CSCOEF (li, 3) = da / dt;
		    *CSCOEF (li, 2) = F_ZERO;
		    *CSCOEF (li, 1) = c;
		    *CSCOEF (li, 0) = c * (tau - xdata[l - 1]) + fdata[l - 1];
		    break_[li - 1] = xdata[l - 1];
		    break_[li] = tau;
		    li += 2;

		}
		else if (j == 2) {

		    tau = xdata[l - 1] - a[l - 1] * dt / da;
		    e = fdata[l - 1] - (imsl_fi_power (a[l - 1], 3)) * dt * dt / (F_SIX *
			da * da);
		    c = df / dt + (imsl_fi_power (a[l - 1], 3)) * dt / (F_SIX * da * da);
		    *CSCOEF (li - 1, 3) = da / dt;
		    *CSCOEF (li - 1, 2) = a[l - 1];
		    *CSCOEF (li - 1, 1) = c + a[l - 1] * a[l - 1] * dt * F_HALF /
			da;
		    *CSCOEF (li - 1, 0) = fdata[l - 1];
		    *CSCOEF (li, 3) = F_ZERO;
		    *CSCOEF (li, 2) = F_ZERO;
		    *CSCOEF (li, 1) = c;
		    *CSCOEF (li, 0) = c * (tau - xdata[l - 1]) + e;
		    break_[li - 1] = xdata[l - 1];
		    break_[li] = tau;
		    li += 2;

		}
		else if (j == 1) {

		    *CSCOEF (li - 1, 3) = F_ZERO;
		    *CSCOEF (li - 1, 2) = F_ZERO;
		    *CSCOEF (li - 1, 1) = df / dt;
		    *CSCOEF (li - 1, 0) = fdata[l - 1];
		    break_[li - 1] = xdata[l - 1];
		    li += 1;

		}
	    }
	}
    }
    *CSCOEF (li - 1, 3) = F_ZERO;
    *CSCOEF (li - 1, 2) = F_ZERO;
    *CSCOEF (li - 1, 1) = F_ZERO;
    *CSCOEF (li - 1, 0) = fdata[m - 1];
    break_[li - 1] = xdata[m - 1];
    *ibreak = li;
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  C6CON/DC6CON (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 26, 1985

    Purpose:    Compute a cubic spline interpolant which is
                consistent with the concavity of the data.

    Usage:      CALL C6CON (SUB, DIAG, SUP, B, N)

    Arguments:
       SUB    - Array containing the sub-diagonal of the matrix.  (Input)
       DIAG   - Array containing the main diagonal of the matrix.
                (Input)
       SUP    - Array containing the super diagonal of the matrix.
                (Input)
       B      - Array that hols the RHS of the matrix.  (Input)
       N      - Number of equations.   (Input)
                N is equal to NDATA-2.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c6con (Mfloat sub[], Mfloat diag[], Mfloat sup[], Mfloat b[], Mint *n)
#else
static void l_c6con (sub, diag, sup, b, n)
    Mfloat      sub[], diag[], sup[], b[];
    Mint       *n;
#endif
{
    Mint        i;


    if (*n <= 1) {
	b[0] /= diag[0];
    }
    else {
	for (i = 2; i <= *n; i++) {
	    sub[i - 1] /= diag[i - 2];
	    diag[i - 1] += -sub[i - 1] * sup[i - 2];
	    b[i - 1] += -sub[i - 1] * b[i - 2];
	}
	b[*n - 1] /= diag[*n - 1];
	for (i = *n - 1; i >= 1; i--) {
	    b[i - 1] = (b[i - 1] - sup[i - 1] * b[i]) / diag[i - 1];
	}
    }
    return;
}				/* end of function */



























/*  -----------------------------------------------------------------------
    IMSL Name:  C2AKM/DC2AKM (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1984

    Purpose:    Compute the Akima cubic spline interpolant.

    Usage:      CALL C2AKM (NDATA, XDATA, FDATA, BREAK, CSCOEF, IPVT)

    Arguments:  (See CSAKM)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c2akm (Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat break_[],
                Mfloat *cscoef, Mint ipvt[])
#else
static void l_c2akm (ndata, xdata, fdata, break_, cscoef, ipvt)
    Mint       *ndata;
    Mfloat      xdata[], fdata[], break_[], *cscoef;
    Mint        ipvt[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
    Mint        i;
    Mfloat      b, rm1, rm2, rm3, rm4, t1, t2;


    imsl_e1psh ("C2AKM ");
    /* CHECK ARGUMENT NDATA */
    if (*ndata < 4) {
	imsl_e1sti (1, *ndata);

/*		imsl_ermes(5, 1, "The number of data points must be 4 or more while NDATA = %(i1) is given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_AT_LEAST_4_PTS);
	goto L_9000;
    }
    /* COPY AND SORT INPUT DATA */
    imsl_c1sor (*ndata, xdata, fdata, break_, cscoef, 4, ipvt);
    if (imsl_n1rty (0) != 0)
	goto L_9000;

    rm3 = (*CSCOEF (1, 0) - *CSCOEF (0, 0)) / (break_[1] - break_[0]);
    t1 = rm3 - (*CSCOEF (1, 0) - *CSCOEF (2, 0)) / (break_[1] - break_[2]);
    rm2 = rm3 + t1;
    rm1 = rm2 + t1;
    /* NOW GET THE SLOPES */
    for (i = 1; i <= *ndata; i++) {
	if (i < *ndata - 1) {
	    rm4 = (*CSCOEF (i + 1, 0) - *CSCOEF (i, 0)) / (break_[i + 1] -
		break_[i]);
	}
	else {
	    rm4 = rm3 - rm2 + rm3;
	}
	t1 = fabs (rm4 - rm3);
	t2 = fabs (rm2 - rm1);
	b = t1 + t2;
	if (b == F_ZERO) {
	    /* IF DENOMINATOR IS ZERO, GET AVERAGE */
	    *CSCOEF (i - 1, 1) = F_HALF * (rm2 + rm3);
	}
	else {
	    *CSCOEF (i - 1, 1) = (t1 * rm2 + t2 * rm3) / b;
	}
	rm1 = rm2;
	rm2 = rm3;
	rm3 = rm4;
    }
    /*
     * COMPUTE THE COEFFICIENTS FOR THE NDATA-1 INTERVALS
     */
    for (i = 1; i <= (*ndata - 1); i++) {
	t1 = F_ONE / (break_[i] - break_[i - 1]);
	t2 = (*CSCOEF (i, 0) - *CSCOEF (i - 1, 0)) * t1;
	b = (*CSCOEF (i - 1, 1) + *CSCOEF (i, 1) - t2 - t2) * t1;
	*CSCOEF (i - 1, 3) = F_SIX * b * t1;
	*CSCOEF (i - 1, 2) = F_TWO * (-b + (t2 - *CSCOEF (i - 1, 1)) * t1);
    }

L_9000:
    ;
    imsl_e1pop ("C2AKM ");
    return;
}				/* end of function */
