#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static VA_LIST_HACK PROTO (l_lin_sol_gen, (Mint n, Mf_complex *a,
                                 Mf_complex *b, va_list argptr));
    static void PROTO (l_lfscg, (Mint *n, Mf_complex *imsl_fac,
                                 Mint *ldfac, Mint *ipvt,
                                 Mf_complex *b, Mint *ipath,
                                 Mf_complex *x));
    static void PROTO (l_l2tcg, (Mint *n, Mf_complex *a, Mint *lda,
                                 Mf_complex *imsl_fac, Mint *ldfac,
                                 Mint *ipvt, Mf_complex *scale));
    static void PROTO (l_l4tcg, (Mint *n, Mf_complex *x,
                                 Mf_complex *scale, Mint *k));
    static void PROTO (l_l2ccg, (Mint *n, Mf_complex *a, Mint *lda,
                                 Mf_complex *imsl_fac, Mint *ldfac,
                                 Mint *ipvt, Mfloat *rcond,
                                 Mf_complex *z));
    static void PROTO (l_l2ncg, (Mint *n, Mf_complex *a, Mint *lda,
                                 Mf_complex *ainv, Mint *ldainv,
                                 Mf_complex *wk, Mint *iwk));
    static void PROTO (l_linct, (Mint *n, Mf_complex *a, Mint *lda,
                                 Mint *ipath, Mf_complex *ainv,
                                 Mint *ldainv));

    static Mf_complex *lv_x;
    static Mfloat lv_rcond;
#ifdef ANSI
    Mf_complex *imsl_c_lin_sol_gen (Mint n, Mf_complex *a,
                                    Mf_complex *b,...)
#else
    Mf_complex *imsl_c_lin_sol_gen (n, a, b, va_alist)
    Mint        n;
    Mf_complex *a;
    Mf_complex *b;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, b);

    E1PSH ("imsl_c_lin_sol_gen", "imsl_z_lin_sol_gen");
    lv_x = NULL;

    IMSL_CALL (l_lin_sol_gen (n, a, b, argptr));
    va_end (argptr);

    E1POP ("imsl_c_lin_sol_gen", "imsl_z_lin_sol_gen");
    return lv_x;
}


#ifdef ANSI
static VA_LIST_HACK l_lin_sol_gen (Mint n, Mf_complex *a, Mf_complex *b,
                              va_list argptr)
#else
static VA_LIST_HACK l_lin_sol_gen (n, a, b, argptr)
    Mint        n;
    Mf_complex *a;
    Mf_complex *b;
    va_list     argptr;
#endif
{
    Mint        code = 1;
    Mint        arg_number = 3;
    Mint        path = 1;
    Mf_complex **factor_ptr = NULL;
    Mf_complex *factor = NULL;
    Mf_complex **inva_ptr = NULL;
    Mf_complex *inva = NULL;
    Mint        a_col_dim = n;
    Mint        fac_col_dim = n;
    Mint        inv_col_dim = n;
    Mfloat     *cond = NULL;
    Mint      **ipvt_ptr = NULL;
    Mint       *ipvt = NULL;

    Mf_complex *fwork = NULL;
    Mint       *iwork = NULL;

    Mint        return_factor = 0;
    Mint        factor_user = 0;
    Mint        return_inverse = 0;
    Mint        inverse_user = 0;
    Mint        inverse_only = 0;
    Mint        factor_only = 0;
    Mint        solve_only = 0;
    Mint        return_cond = 0;
    Mint        result_user = 0;
    Mint        error = 0;
    Mint        lda = n;
    Mint        i;
    Mfloat      rcond;
    static Mf_complex c_zero = {0.0, 0.0};
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_RETURN_USER:
	    result_user = 1;
	    lv_x = va_arg (argptr, Mf_complex *);
	    if (!lv_x) {
		imsl_e1stl (1, "x");
		imsl_e1stl (2, "IMSL_RETURN_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    arg_number++;
	    break;
	case IMSL_A_COL_DIM:
	    a_col_dim = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_TRANSPOSE:
	    path = 2;
	    break;
	case IMSL_FACTOR:
	    return_factor = 1;
	    factor_user = 0;
	    ipvt_ptr = va_arg (argptr, Mint **);
	    factor_ptr = va_arg (argptr, Mf_complex **);
	    arg_number += 2;
	    break;
	case IMSL_FACTOR_USER:
	    return_factor = 1;
	    factor_user = 1;
	    ipvt = va_arg (argptr, Mint *);
	    if (!ipvt) {
		imsl_e1stl (1, "ipvt");
		imsl_e1stl (2, "IMSL_FACTOR");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    factor = va_arg (argptr, Mf_complex *);
	    if (!factor) {
		imsl_e1stl (1, "factor");
		imsl_e1stl (2, "IMSL_FACTOR");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_2);
		++error;
	    }
	    arg_number += 2;
	    break;
	case IMSL_FAC_COL_DIM:
	    fac_col_dim = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_INVERSE:
	    return_inverse = 1;
	    inverse_user = 0;
	    inva_ptr = va_arg (argptr, Mf_complex **);
	    arg_number++;
	    break;
	case IMSL_INVERSE_USER:
	    return_inverse = 1;
	    inverse_user = 1;
	    inva = va_arg (argptr, Mf_complex *);
	    if (!inva) {
		imsl_e1stl (1, "inva");
		imsl_e1stl (2, "IMSL_INVERSE_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    arg_number++;
	    break;
	case IMSL_INV_COL_DIM:
	    inv_col_dim = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_CONDITION:
	    return_cond = 1;
	    cond = va_arg (argptr, Mfloat *);
            arg_number++;
	    if (!cond) {
		imsl_e1stl (1, "cond");
		imsl_e1stl (2, "IMSL_CONDITION");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    break;
	case IMSL_FACTOR_ONLY:
	    factor_only = 1;
	    break;
	case IMSL_SOLVE_ONLY:
	    solve_only = 1;
	    break;
	case IMSL_INVERSE_ONLY:
	    inverse_only = 1;
	    break;
	case 0:
	    break;
	default:
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, arg_number);
	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    ++error;
	    break;
	}
    }

    if (!a && !solve_only) {
	imsl_e1stl (1, "a");
	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
	++error;
    }
    if (!b && !(factor_only || inverse_only)) {
	imsl_e1stl (1, "b");
	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
	++error;
    }

    if (error)
	return argptr;

    if (n <= 0) {
	imsl_e1sti (1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
	++error;
    }
    else if (n > a_col_dim) {
	imsl_e1sti (1, n);
	imsl_e1sti (2, a_col_dim);
	imsl_e1stl (1, "a");
	imsl_ermes (IMSL_TERMINAL, IMSL_COL_DIM_LESS_ORDER);
	++error;
    }
    if (n > fac_col_dim) {
	imsl_e1sti (1, n);
	imsl_e1sti (2, fac_col_dim);
	imsl_e1stl (1, "factor");
	imsl_ermes (IMSL_TERMINAL, IMSL_COL_DIM_LESS_ORDER);
	++error;
    }
    if (solve_only + factor_only + inverse_only > 1) {
	imsl_ermes (IMSL_TERMINAL, IMSL_BAD_SOLVE_FACTOR_INVERSE);
	++error;
    }

    if (error)
	return argptr;

    if (solve_only) {
	if (!factor_user) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_SPECIFY_SOLVE_ONLY);
	    ++error;
	}
	if (return_cond) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_CONDITION_ONLY_SPECIFIER);
	    ++error;
	}
	if (error)
	    return argptr;
    }
    if (inverse_only && !return_inverse) {
	imsl_ermes (IMSL_TERMINAL, IMSL_INVERSE_ONLY_SPECIFIER);
	++error;
    }
    if (factor_only && !return_factor) {
	imsl_ermes (IMSL_TERMINAL, IMSL_FACTOR_ONLY_SPECIFIER);
	++error;
    }
    if (error)
	return argptr;

    if (!inverse_only) {
	if (!solve_only) {
	    Mint        lda = n;
	    Mint        ldfac = n;
	    if (!factor_user) {
		factor = (Mf_complex *) imsl_malloc (n * fac_col_dim * sizeof (Mf_complex));
		ipvt = (Mint *) imsl_malloc (n * sizeof (Mint));
		if (!factor || !ipvt) {
		    imsl_e1stl (1, "n");
		    imsl_e1sti (1, n);
		    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
		    ++error;
		    goto FREE_SPACE;
		}
	    }
	    else if (fac_col_dim > n) {
		imsl_c_m1ran (n, fac_col_dim, factor, factor);
		if (error = (imsl_n1rty (1) > 3) ? 1 : 0)
		    goto FREE_SPACE;
	    }
	    fwork = (Mf_complex *) imsl_malloc (n * sizeof (Mf_complex));
	    if (!fwork) {
		imsl_e1stl (1, "n");
		imsl_e1sti (1, n);
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
		++error;
		goto FREE_SPACE;
	    }
	    imsl_c_m1ran (n, a_col_dim, a, a);
	    if (!(error = imsl_n1rty (1) > 3)) {
		if (factor_only && !return_cond) {
		    l_l2tcg (&n, a, &lda, factor, &ldfac, ipvt, fwork);
		    error = imsl_n1rty (1) > 3 ? 1 : 0;
		}
		else {
		    l_l2ccg (&n, a, &lda, factor, &ldfac, ipvt, &rcond, fwork);
		    error = imsl_n1rty (1) > 3 ? 1 : 0;
		    if (!error && return_cond) {
			if (rcond > imsl_amach (1)) {
			    *cond = 1.0 / rcond;
			}
			else {
			    *cond = imsl_amach (7);
			}
		    }
		}
		imsl_c_m1ran (a_col_dim, n, a, a);
	    }
	    imsl_free (fwork);
	    fwork = NULL;
	}
	if (!factor_only && !error) {
	    if (!result_user) {
		lv_x = (Mf_complex *) imsl_malloc (n * sizeof (Mf_complex));
		if (lv_x == NULL) {
		    imsl_e1stl (1, "n");
		    imsl_e1sti (1, n);
		    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
		    ++error;
		    goto FREE_SPACE;
		}
	    }
	    if (solve_only) {
		imsl_c_m1ran (n, fac_col_dim, factor, factor);
		if (error = imsl_n1rty (1) > 3 ? 1 : 0)
		    goto FREE_SPACE;
	    }
	    l_lfscg (&n, factor, &n, ipvt, b, &path, lv_x);
	    error = imsl_n1rty (1) > 3 ? 1 : 0;
	}
    }
    if (return_inverse && !error) {
	Mint        ldinv;
	if (!inverse_user) {
	    inva = (Mf_complex *) imsl_malloc (n * inv_col_dim * sizeof (Mf_complex));
	    if (!inva) {
		imsl_e1stl (1, "n");
		imsl_e1sti (1, n);
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
		++error;
		goto FREE_SPACE;
	    }
	}

	fwork = (Mf_complex *) imsl_malloc ((n + n * (n - 1) / 2) * sizeof (Mf_complex));
	iwork = (Mint *) imsl_malloc (n * sizeof (Mint));
	if (!iwork || !fwork) {
	    imsl_e1stl (1, "n");
	    imsl_e1sti (1, n);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    ++error;
	    goto FREE_SPACE;
	}
	if (inverse_only) {
	    imsl_c_m1ran (n, a_col_dim, a, a);
	    if (!(error = (imsl_n1rty (1) > 3) ? 1 : 0)) {
		ldinv = n;
		if ((inv_col_dim > n) && inverse_user) {
		    imsl_c_m1ran (n, inv_col_dim, inva, inva);
		    error = (imsl_n1rty (1) > 3) ? 1 : 0;
		}
	    }
	}
	else {
	    ldinv = inv_col_dim;
	    lda = a_col_dim;
	}
	if (!error) {
	    l_l2ncg (&n, a, &lda, inva, &ldinv, fwork, iwork);
	    error = imsl_n1rty (1) > 3 ? 1 : 0;
	    if (inverse_only) {
		if (!error && return_cond) {
		    if (lv_rcond > imsl_amach (1)) {
			*cond = 1.0 / lv_rcond;
		    }
		    else {
			*cond = imsl_amach (7);
		    }
		}
		imsl_c_m1ran (a_col_dim, n, a, a);
		imsl_c_m1ran (inv_col_dim, n, inva, inva);
	    }
	}
    }
FREE_SPACE:

    if (!factor_user) {
	if (error || !return_factor) {
	    if (factor)
		imsl_free (factor);
	    if (ipvt)
		imsl_free (ipvt);
	    factor = NULL;
	}
	else {
	    imsl_c_m1ran (fac_col_dim, n, factor, factor);
	    for (i = n; i < fac_col_dim; i++) {
		imsl_cset (&n, &c_zero, (factor + i), &fac_col_dim);
	    }
	    *factor_ptr = factor;
	    *ipvt_ptr = ipvt;
	}
	factor = NULL;
	ipvt = NULL;
    }
    else {
	imsl_c_m1ran (fac_col_dim, n, factor, factor);
    }
    if (!inverse_user) {
	if (error || !return_inverse) {
	    if (inva)
		imsl_free (inva);
	}
	else {
	    for (i = n; i < inv_col_dim; i++) {
		imsl_cset (&n, &c_zero, (inva + i), &inv_col_dim);
	    }

	    *inva_ptr = inva;
	}
	inva = NULL;
    }
    if (!result_user) {
	if (error && lv_x) {
	    imsl_free (lv_x);
	    lv_x = NULL;
	}
    }
    if (fwork)
	imsl_free (fwork);
    if (iwork)
	imsl_free (iwork);

RETURN:
    return (argptr);
}


/* Structured by FOR_STRUCT, v0.2, on 08/28/90 at 10:49:18
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  LFSCG/DLFSCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 25, 1987

    Purpose:    Solve a d_complex general system of linear equations given
                the LU factorization of the coefficient matrix.

    Usage:      CALL LFSCG (N, FAC, LDFAC, IPVT, B, IPATH, X)

    Arguments:
       N      - Number of equations.  (Input)
       FAC    - Complex N by N matrix containing the LU factorization of
                the coefficient matrix A as output from subroutine
                LFCCG/DLFCCG or LFTCG/DLFTCG.  (Input)
       LDFAC  - Leading dimension of FAC exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPVT   - Vector of length N containing the pivoting information
                for the LU factorization of A as output from subroutine
                LFCCG/DLFCCG or LFTCG/DLFTCG.  (Input)
       B      - Complex vector of length N containing the right-hand side
                of the linear system.  (Input)
       IPATH  - Path indicator.  (Input)
                IPATH = 1 means the system A*X = B is solved.
                IPATH = 2 means the system ctrans(A)*X = B is solved,
                        where ctrans(A) is the conjugate transpose of A.
       X      - Complex vector of length N containing the solution to the
                linear system.  (Output)
                If B is not needed, B and X can share the same storage
                locations.

    GAMS:       D2c1

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_lfscg (Mint *n, Mf_complex *imsl_fac, Mint *ldfac, Mint *ipvt,
                Mf_complex *b, Mint *ipath, Mf_complex *x)
#else
static void l_lfscg (n, imsl_fac, ldfac, ipvt, b, ipath, x)
    Mint       *n;
    Mf_complex *imsl_fac;
    Mint       *ldfac, ipvt[];
    Mf_complex  b[];
    Mint       *ipath;
    Mf_complex  x[];
#endif
{
#define FAC(I_,J_)	(imsl_fac+(I_)*(aldfac)+(J_))
    Mint        aldfac = *ldfac;

    Mint        _l0, _l1, _l2, k, l;
    Mfloat      big, small;
    Mf_complex  t;


#define L3CCT(zdum)	(Mfloat)(fabs( imsl_fc_convert( (zdum) ) ) + fabs( imsl_c_aimag( (zdum) ) ))

    imsl_e1psh ("l_lfscg");

    if (*n <= 0) {
	imsl_e1sti (1, *n);

	/*
	 * (5, 1, "The number of equations must be positive while N = %(i1)
	 * is given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NUM_OF_EQUATIONS);
	goto L_9000;
    }
    if (*n > *ldfac) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *ldfac);

	/*
	 * (5, 2, "The order of the matrix must be less than or equal to its
	 * leading dimension while N = %(i1) and LDFAC = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_LDFAC_LESS_ORDER);
	goto L_9000;
    }
    /*
     * COPY B INTO X AND USE X TO PRESERVE INPUT
     */
    _l0 = 1;
    _l1 = 1;
    imsl_ccopy (n, b, &_l0, x, &_l1);

    small = imsl_amach (1);
    big = imsl_amach (2);
    if (small * big < F_ONE)
	small = F_ONE / big;
    if (*ipath == 1) {
	/*
	 * IPATH = 1 , SOLVE A * X = B FIRST SOLVE L*Y = B
	 */
	for (k = 1; k <= (*n - 1); k++) {
	    l = ipvt[k - 1];
	    t = x[l - 1];
	    if (l != k) {
		x[l - 1] = x[k - 1];
		x[k - 1] = t;
	    }
	    _l0 = *n - k;
	    _l1 = 1;
	    _l2 = 1;
	    imsl_caxpy (&_l0, &t, FAC (k - 1, k), &_l1,
		&x[k], &_l2);
	}
	/* NOW SOLVE U*X = Y */
	for (k = *n; k >= 1; k--) {
	    if (L3CCT (*FAC (k - 1, k - 1)) <= small) {

		/*
		 * (5, 3, "The input matrix is singular.  Some of the
		 * diagonal elements of the upper triangular matrix U of the
		 * LU factorization are close to zero.");
		 */
		imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_MATRIX);
		goto L_9000;
	    }
	}
	_l0 = 1;
	imsl_ctrsv ("U", sizeof ("U"), "N", sizeof ("N"), "N", sizeof ("N"), n, FAC (0, 0), ldfac, x, &_l0);

    }
    else if (*ipath == 2) {
	/*
	 * IPATH = 2, SOLVE CTRANS(A) * X = B FIRST SOLVE CTRANS(U)*Y = B
	 */
	for (k = 1; k <= *n; k++) {
	    if (L3CCT (*FAC (k - 1, k - 1)) <= small) {

		/*
		 * (5, 4, "The input matrix is singular.  Some of the
		 * diagonal elements of the upper triangular matrix U of the
		 * LU factorization are close to zero.");
		 */
		imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_MATRIX);
		goto L_9000;
	    }
	}
	_l0 = 1;
	imsl_ctrsv ("U", sizeof ("U"), "C", sizeof ("C"), "N", sizeof ("N"), n, FAC (0, 0), ldfac, x, &_l0);
	/* NOW SOLVE TRANS(L)*X = Y */
	for (k = *n - 1; k >= 1; k--) {
	    _l0 = *n - k;
	    _l1 = 1;
	    _l2 = 1;
	    x[k - 1] = imsl_c_add (x[k - 1], imsl_cdotc (&_l0, FAC (k - 1, k),
		    &_l1, &x[k], &_l2));
	    l = ipvt[k - 1];
	    if (l != k) {
		t = x[l - 1];
		x[l - 1] = x[k - 1];
		x[k - 1] = t;
	    }
	}

    }
    else {
	imsl_e1sti (1, *ipath);

	/*
	 * (5, 5, "IPATH must be either 1 or 2 while a value of %(i1) is
	 * given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_IPATH_RANGE_3);
    }

L_9000:
    imsl_e1pop ("l_lfscg");
    return;
#undef	L3CCT
#undef	FAC
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  L2TCG/DL2TCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 25, 1990

    Purpose:    Compute the LU factorization of a d_complex general matrix.

    Usage:      CALL L2TCG (N, A, LDA, FAC, LDFAC, IPVT, SCALE)

    Arguments:  See LFTCG/DLFTCG.

    Remarks:    See LFTCG/DLFTCG.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------

 */
#ifdef ANSI
static void l_l2tcg (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *imsl_fac,
                Mint *ldfac, Mint *ipvt, Mf_complex *scale)
#else
static void l_l2tcg (n, a, lda, imsl_fac, ldfac, ipvt, scale)
    Mint       *n;
    Mf_complex *a;
    Mint       *lda;
    Mf_complex *imsl_fac;
    Mint       *ldfac, ipvt[];
    Mf_complex  scale[];
#endif
{
#define A(I_,J_)	(a+(I_)*(alda)+(J_))
#define FAC(I_,J_)	(imsl_fac+(I_)*(aldfac)+(J_))
    Mint        aldfac = *ldfac;
    Mint        _l0, i, indj, info, j, k, ktemp, m0, m1, m2, m3, m4, m5, m6,
                m7;
    Mfloat      big, small;
    Mf_complex  t, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    static Mf_complex fpm1 = {-1.0e0, 0.0e0};

    Mf_complex ctemp1, ctemp2, ctemp3, ctemp4;
    Mf_complex ctemp5, ctemp6, ctemp7, ctemp8;



#define L3CCT(z1)	(Mfloat)(fabs( imsl_fc_convert( (z1) ) ) + fabs( imsl_c_aimag( (z1) ) ))

    /*
     * this code is for computer types: aliant and necv
     */
    imsl_e1psh ("l_l2tcg");

    if (*n <= 0) {
	imsl_e1sti (1, *n);

	/*
	 * (5, 1, "The order of the matrix must be positive while N = %(i1)
	 * is given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
	goto L_9900;
    }
    if (*n > *lda) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *lda);

	/*
	 * (5, 2, "The order of the matrix must be less than or equal to its
	 * leading dimension while N = %(i1) and LDA = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
	goto L_9900;
    }
    if (*n > *ldfac) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *ldfac);

	/*
	 * (5, 3, "The order of the matrix must be less than or equal to its
	 * leading dimension while N = %(i1) and LDFAC = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_LDFAC_LESS_ORDER);
	goto L_9900;
    }
    /* Preserve a copy of the input matrix */
    imsl_ccgcg (n, a, lda, imsl_fac, ldfac);

    /*
     * Gaussian elimination with scaled partial pivoting using method lu***
     * 
     * A brief description of the algorithm follows:
     * 
     * For A = | 11  12 |   block 11 is 8 x 8 | 21  22 |
     * 
     * step 1. Factor n by 8 matrix trans(11 21)
     * 
     * step 2. Update 12 and 22 from back to front along the columns
     * 
     * step 3. Repeat procedure on the 22 block
     */
    info = 0;
    j = 1;
    ktemp = 1;
    small = imsl_amach (1);
    big = imsl_amach (2);
    if (small * big < F_ONE)
	small = F_ONE / big;
    j = 1;
    ktemp = 1;
    /*
     * determine scale factors for scaled partial pivoting
     */
    for (i = 1; i <= *n; i++) {
	indj = imsl_icamax (n, FAC (0, i - 1), ldfac);
	scale[i - 1] = imsl_cf_convert (imsl_c_abs (*FAC (indj - 1, i - 1)), F_ZERO);
	if (imsl_fc_convert (scale[i - 1]) < small) {
	    scale[i - 1] = imsl_cf_convert (F_ONE, F_ZERO);
	}
	else {
	    scale[i - 1] = imsl_c_div (imsl_cf_convert (F_ONE, F_ZERO), scale[i - 1]);
	}
    }
    /* begin factorization step */
L_20:
    if (j >= *n)
	goto L_170;
    /*
     * determine the pivot element for column j of L
     */
    _l0 = *n - j + 1;
    l_l4tcg (&_l0, FAC (j - 1, j - 1), &scale[j - 1], &m0);
    m0 += j - 1;
    ipvt[j - 1] = m0;

    if (L3CCT (*FAC (j - 1, m0 - 1)) > small) {
	if (m0 != j)
	    ktemp = -ktemp;
	/*
	 * swap element j,j with the pivot element and the scaling vector
	 */
	t = scale[j - 1];
	scale[j - 1] = scale[m0 - 1];
	scale[m0 - 1] = t;
	t = *FAC (j - 1, m0 - 1);
	*FAC (j - 1, m0 - 1) = *FAC (j - 1, j - 1);
	*FAC (j - 1, j - 1) = t;
	t = imsl_c_div (fpm1, t);
	/* update U(j,j+1) */
	t1 = *FAC (j, m0 - 1);
	*FAC (j, m0 - 1) = *FAC (j, j - 1);
	*FAC (j, j - 1) = t1;
	/* update columns j of L and j+1 of A */
	for (i = j + 1; i <= *n; i++) {
	    *FAC (j - 1, i - 1) = imsl_c_mul (*FAC (j - 1, i - 1), t);
	    *FAC (j, i - 1) = imsl_c_add (*FAC (j, i - 1), imsl_c_mul (t1, *FAC (j - 1, i - 1)));
	}
    }
    else {
	ktemp = 0;
	info = m0;
    }

    if (j + 1 >= *n)
	goto L_170;
    /*
     * determine the pivot element for column j+1
     */
    _l0 = *n - j;
    l_l4tcg (&_l0, FAC (j, j), &scale[j], &m1);
    m1 += j;
    ipvt[j] = m1;
    t = scale[j];
    scale[j] = scale[m1 - 1];
    scale[m1 - 1] = t;
    /*
     * swap element j+1,j+1 with the pivot element
     */
    t = *FAC (j, m1 - 1);
    *FAC (j, m1 - 1) = *FAC (j, j);
    *FAC (j, j) = t;

    if (L3CCT (t) > small) {
	t = imsl_c_div (fpm1, t);
	if (m1 != j + 1)
	    ktemp = -ktemp;
    }
    else {
	ktemp = 0;
	info = m1;
    }
    /*
     * swap subdiagonal elements of row j+1 with the pivot row in L
     */
    t1 = *FAC (j - 1, m1 - 1);
    *FAC (j - 1, m1 - 1) = *FAC (j - 1, j);
    *FAC (j - 1, j) = t1;
    /* update U(j,j+2) to U(j+1,j+2) */
    t1 = *FAC (j + 1, m0 - 1);
    *FAC (j + 1, m0 - 1) = *FAC (j + 1, j - 1);
    *FAC (j + 1, j - 1) = t1;
    t2 = imsl_c_add (*FAC (j + 1, m1 - 1), imsl_c_mul (t1, *FAC (j - 1, j)));
    *FAC (j + 1, m1 - 1) = *FAC (j + 1, j);
    *FAC (j + 1, j) = t2;
    /* update column j+1 of L and j+2 of A */
    for (i = j + 2; i <= *n; i++) {
	*FAC (j, i - 1) = imsl_c_mul (*FAC (j, i - 1), t);
	*FAC (j + 1, i - 1) = imsl_c_add (imsl_c_add (*FAC (j + 1, i - 1), imsl_c_mul (t1,
		    *FAC (j - 1, i - 1))), imsl_c_mul (t2, *FAC (j, i - 1)));
    }

    if (j + 2 >= *n)
	goto L_160;
    /*
     * determine the pivot element for column j+2
     */
    _l0 = *n - j - 1;
    l_l4tcg (&_l0, FAC (j + 1, j + 1), &scale[j + 1], &m2);
    m2 += j + 1;
    ipvt[j + 1] = m2;
    t = scale[j + 1];
    scale[j + 1] = scale[m2 - 1];
    scale[m2 - 1] = t;
    /*
     * swap element j+2,j+2 with the pivot element
     */
    t = *FAC (j + 1, m2 - 1);
    *FAC (j + 1, m2 - 1) = *FAC (j + 1, j + 1);
    *FAC (j + 1, j + 1) = t;

    if (L3CCT (t) > small) {
	t = imsl_c_div (fpm1, t);
	if (m2 != j + 2)
	    ktemp = -ktemp;
    }
    else {
	ktemp = 0;
	info = m2;
    }
    /*
     * swap subdiagonal elements of row j+2 with the pivot row in L
     */
    t1 = *FAC (j - 1, m2 - 1);
    *FAC (j - 1, m2 - 1) = *FAC (j - 1, j + 1);
    *FAC (j - 1, j + 1) = t1;
    t1 = *FAC (j, m2 - 1);
    *FAC (j, m2 - 1) = *FAC (j, j + 1);
    *FAC (j, j + 1) = t1;
    /* update U(j,j+3) to U(j+2,j+3) */
    t1 = *FAC (j + 2, m0 - 1);
    *FAC (j + 2, m0 - 1) = *FAC (j + 2, j - 1);
    *FAC (j + 2, j - 1) = t1;
    t2 = imsl_c_add (*FAC (j + 2, m1 - 1), imsl_c_mul (t1, *FAC (j - 1, j)));
    *FAC (j + 2, m1 - 1) = *FAC (j + 2, j);
    *FAC (j + 2, j) = t2;
    t3 = imsl_c_add (imsl_c_add (*FAC (j + 2, m2 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 1))),
	imsl_c_mul (t2, *FAC (j, j + 1)));
    *FAC (j + 2, m2 - 1) = *FAC (j + 2, j + 1);
    *FAC (j + 2, j + 1) = t3;
    /* update column j+2 of L and j+3 of A */
    for (i = j + 3; i <= *n; i++) {
	*FAC (j + 1, i - 1) = imsl_c_mul (*FAC (j + 1, i - 1), t);
	*FAC (j + 2, i - 1) = imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 2, i - 1), imsl_c_mul (t1,
			*FAC (j - 1, i - 1))), imsl_c_mul (t2, *FAC (j, i - 1))), imsl_c_mul (t3, *FAC (j + 1, i - 1)));
    }

    if (j + 3 >= *n)
	goto L_150;
    /*
     * determine the pivot element for column j+3
     */
    _l0 = *n - j - 2;
    l_l4tcg (&_l0, FAC (j + 2, j + 2), &scale[j + 2], &m3);
    m3 += j + 2;
    ipvt[j + 2] = m3;
    t = scale[j + 2];
    scale[j + 2] = scale[m3 - 1];
    scale[m3 - 1] = t;
    /*
     * swap element j+3,j+3 with the pivot element
     */
    t = *FAC (j + 2, m3 - 1);
    *FAC (j + 2, m3 - 1) = *FAC (j + 2, j + 2);
    *FAC (j + 2, j + 2) = t;

    if (L3CCT (t) > small) {
	t = imsl_c_div (fpm1, t);
	if (m3 != j + 3)
	    ktemp = -ktemp;
    }
    else {
	ktemp = 0;
	info = m3;
    }
    /*
     * swap the subdiagonal elements of row j+3 with the pivot row in L
     */
    t1 = *FAC (j - 1, m3 - 1);
    *FAC (j - 1, m3 - 1) = *FAC (j - 1, j + 2);
    *FAC (j - 1, j + 2) = t1;
    t1 = *FAC (j, m3 - 1);
    *FAC (j, m3 - 1) = *FAC (j, j + 2);
    *FAC (j, j + 2) = t1;
    t1 = *FAC (j + 1, m3 - 1);
    *FAC (j + 1, m3 - 1) = *FAC (j + 1, j + 2);
    *FAC (j + 1, j + 2) = t1;
    /* update U(j,j+4) to U(j+3,j+4) */
    t1 = *FAC (j + 3, m0 - 1);
    *FAC (j + 3, m0 - 1) = *FAC (j + 3, j - 1);
    *FAC (j + 3, j - 1) = t1;
    t2 = imsl_c_add (*FAC (j + 3, m1 - 1), imsl_c_mul (t1, *FAC (j - 1, j)));
    *FAC (j + 3, m1 - 1) = *FAC (j + 3, j);
    *FAC (j + 3, j) = t2;
    t3 = imsl_c_add (imsl_c_add (*FAC (j + 3, m2 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 1))),
	imsl_c_mul (t2, *FAC (j, j + 1)));
    *FAC (j + 3, m2 - 1) = *FAC (j + 3, j + 1);
    *FAC (j + 3, j + 1) = t3;
    t4 = imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 3, m3 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 2))),
	    imsl_c_mul (t2, *FAC (j, j + 2))), imsl_c_mul (t3, *FAC (j + 1, j + 2)));
    *FAC (j + 3, m3 - 1) = *FAC (j + 3, j + 2);
    *FAC (j + 3, j + 2) = t4;
    /* update column j+3 of L and j+4 of A */
    for (i = j + 4; i <= *n; i++) {
	*FAC (j + 2, i - 1) = imsl_c_mul (*FAC (j + 2, i - 1), t);
	*FAC (j + 3, i - 1) = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 3, i - 1),
			imsl_c_mul (t1, *FAC (j - 1, i - 1))), imsl_c_mul (t2, *FAC (j, i - 1))), imsl_c_mul (t3,
		    *FAC (j + 1, i - 1))), imsl_c_mul (t4, *FAC (j + 2, i - 1)));
    }

    if (j + 4 >= *n)
	goto L_140;
    /*
     * determine the pivot element for column j+4
     */
    _l0 = *n - j - 3;
    l_l4tcg (&_l0, FAC (j + 3, j + 3), &scale[j + 3], &m4);
    m4 += j + 3;
    ipvt[j + 3] = m4;
    t = scale[j + 3];
    scale[j + 3] = scale[m4 - 1];
    scale[m4 - 1] = t;
    /*
     * swap element j+4,j+4 with the pivot element
     */
    t = *FAC (j + 3, m4 - 1);
    *FAC (j + 3, m4 - 1) = *FAC (j + 3, j + 3);
    *FAC (j + 3, j + 3) = t;

    if (L3CCT (t) > small) {
	t = imsl_c_div (fpm1, t);
	if (m4 != j + 4)
	    ktemp = -ktemp;
    }
    else {
	ktemp = 0;
	info = m4;
    }
    /*
     * swap the subdiagonal elements of row j+4 with the pivot row in L
     */
    t1 = *FAC (j - 1, m4 - 1);
    *FAC (j - 1, m4 - 1) = *FAC (j - 1, j + 3);
    *FAC (j - 1, j + 3) = t1;
    t1 = *FAC (j, m4 - 1);
    *FAC (j, m4 - 1) = *FAC (j, j + 3);
    *FAC (j, j + 3) = t1;
    t1 = *FAC (j + 1, m4 - 1);
    *FAC (j + 1, m4 - 1) = *FAC (j + 1, j + 3);
    *FAC (j + 1, j + 3) = t1;
    t1 = *FAC (j + 2, m4 - 1);
    *FAC (j + 2, m4 - 1) = *FAC (j + 2, j + 3);
    *FAC (j + 2, j + 3) = t1;
    /* update U(j,j+5) to U(j+4,j+5) */
    t1 = *FAC (j + 4, m0 - 1);
    *FAC (j + 4, m0 - 1) = *FAC (j + 4, j - 1);
    *FAC (j + 4, j - 1) = t1;
    t2 = imsl_c_add (*FAC (j + 4, m1 - 1), imsl_c_mul (t1, *FAC (j - 1, j)));
    *FAC (j + 4, m1 - 1) = *FAC (j + 4, j);
    *FAC (j + 4, j) = t2;
    t3 = imsl_c_add (imsl_c_add (*FAC (j + 4, m2 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 1))),
	imsl_c_mul (t2, *FAC (j, j + 1)));
    *FAC (j + 4, m2 - 1) = *FAC (j + 4, j + 1);
    *FAC (j + 4, j + 1) = t3;
    t4 = imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 4, m3 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 2))),
	    imsl_c_mul (t2, *FAC (j, j + 2))), imsl_c_mul (t3, *FAC (j + 1, j + 2)));
    *FAC (j + 4, m3 - 1) = *FAC (j + 4, j + 2);
    *FAC (j + 4, j + 2) = t4;
    t5 = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 4, m4 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 3))),
		imsl_c_mul (t2, *FAC (j, j + 3))), imsl_c_mul (t3, *FAC (j + 1, j + 3))), imsl_c_mul (t4,
	    *FAC (j + 2, j + 3)));
    *FAC (j + 4, m4 - 1) = *FAC (j + 4, j + 3);
    *FAC (j + 4, j + 3) = t5;
    /* update column j+4 of L and j+5 of A */
    for (i = j + 5; i <= *n; i++) {
	*FAC (j + 3, i - 1) = imsl_c_mul (*FAC (j + 3, i - 1), t);
	*FAC (j + 4, i - 1) = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 4, i - 1),
			    imsl_c_mul (t1, *FAC (j - 1, i - 1))), imsl_c_mul (t2, *FAC (j, i - 1))), imsl_c_mul (t3,
			*FAC (j + 1, i - 1))), imsl_c_mul (t4, *FAC (j + 2, i - 1))), imsl_c_mul (t5,
		*FAC (j + 3, i - 1)));
    }

    if (j + 5 >= *n)
	goto L_130;
    /*
     * determine the pivot element for column j+5
     */
    _l0 = *n - j - 4;
    l_l4tcg (&_l0, FAC (j + 4, j + 4), &scale[j + 4], &m5);
    m5 += j + 4;
    ipvt[j + 4] = m5;
    t = scale[j + 4];
    scale[j + 4] = scale[m5 - 1];
    scale[m5 - 1] = t;
    /*
     * swap element j+5,j+5 with the pivot element
     */
    t = *FAC (j + 4, m5 - 1);
    *FAC (j + 4, m5 - 1) = *FAC (j + 4, j + 4);
    *FAC (j + 4, j + 4) = t;

    if (L3CCT (t) > small) {
	t = imsl_c_div (fpm1, t);
	if (m5 != j + 5)
	    ktemp = -ktemp;
    }
    else {
	ktemp = 0;
	info = m5;
    }
    /*
     * swap the subdiagonal elements of row j+5 with the pivot row in L
     */
    t1 = *FAC (j - 1, m5 - 1);
    *FAC (j - 1, m5 - 1) = *FAC (j - 1, j + 4);
    *FAC (j - 1, j + 4) = t1;
    t1 = *FAC (j, m5 - 1);
    *FAC (j, m5 - 1) = *FAC (j, j + 4);
    *FAC (j, j + 4) = t1;
    t1 = *FAC (j + 1, m5 - 1);
    *FAC (j + 1, m5 - 1) = *FAC (j + 1, j + 4);
    *FAC (j + 1, j + 4) = t1;
    t1 = *FAC (j + 2, m5 - 1);
    *FAC (j + 2, m5 - 1) = *FAC (j + 2, j + 4);
    *FAC (j + 2, j + 4) = t1;
    t1 = *FAC (j + 3, m5 - 1);
    *FAC (j + 3, m5 - 1) = *FAC (j + 3, j + 4);
    *FAC (j + 3, j + 4) = t1;
    /* update U(j,j+6) to U(j+5,j+6) */
    t1 = *FAC (j + 5, m0 - 1);
    *FAC (j + 5, m0 - 1) = *FAC (j + 5, j - 1);
    *FAC (j + 5, j - 1) = t1;
    t2 = imsl_c_add (*FAC (j + 5, m1 - 1), imsl_c_mul (t1, *FAC (j - 1, j)));
    *FAC (j + 5, m1 - 1) = *FAC (j + 5, j);
    *FAC (j + 5, j) = t2;
    t3 = imsl_c_add (imsl_c_add (*FAC (j + 5, m2 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 1))),
	imsl_c_mul (t2, *FAC (j, j + 1)));
    *FAC (j + 5, m2 - 1) = *FAC (j + 5, j + 1);
    *FAC (j + 5, j + 1) = t3;
    t4 = imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 5, m3 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 2))),
	    imsl_c_mul (t2, *FAC (j, j + 2))), imsl_c_mul (t3, *FAC (j + 1, j + 2)));
    *FAC (j + 5, m3 - 1) = *FAC (j + 5, j + 2);
    *FAC (j + 5, j + 2) = t4;
    t5 = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 5, m4 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 3))),
		imsl_c_mul (t2, *FAC (j, j + 3))), imsl_c_mul (t3, *FAC (j + 1, j + 3))), imsl_c_mul (t4,
	    *FAC (j + 2, j + 3)));
    *FAC (j + 5, m4 - 1) = *FAC (j + 5, j + 3);
    *FAC (j + 5, j + 3) = t5;
    t6 = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 5, m5 - 1), imsl_c_mul (t1,
			    *FAC (j - 1, j + 4))), imsl_c_mul (t2, *FAC (j, j + 4))), imsl_c_mul (t3, *FAC (j + 1, j + 4))),
	    imsl_c_mul (t4, *FAC (j + 2, j + 4))), imsl_c_mul (t5, *FAC (j + 3, j + 4)));
    *FAC (j + 5, m5 - 1) = *FAC (j + 5, j + 4);
    *FAC (j + 5, j + 4) = t6;
    /* update column j+5 of L and j+6 of A */
    for (i = j + 6; i <= *n; i++) {
	*FAC (j + 4, i - 1) = imsl_c_mul (*FAC (j + 4, i - 1), t);
	*FAC (j + 5, i - 1) = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 5, i - 1),
				imsl_c_mul (t1, *FAC (j - 1, i - 1))), imsl_c_mul (t2, *FAC (j, i - 1))), imsl_c_mul (t3,
			    *FAC (j + 1, i - 1))), imsl_c_mul (t4, *FAC (j + 2, i - 1))), imsl_c_mul (t5,
		    *FAC (j + 3, i - 1))), imsl_c_mul (t6, *FAC (j + 4, i - 1)));
    }

    if (j + 6 >= *n)
	goto L_120;
    /*
     * determine the pivot element for column j+6
     */
    _l0 = *n - j - 5;
    l_l4tcg (&_l0, FAC (j + 5, j + 5), &scale[j + 5], &m6);
    m6 += j + 5;
    ipvt[j + 5] = m6;
    t = scale[j + 5];
    scale[j + 5] = scale[m6 - 1];
    scale[m6 - 1] = t;
    /*
     * swap element j+6,j+6 with the pivot element
     */
    t = *FAC (j + 5, m6 - 1);
    *FAC (j + 5, m6 - 1) = *FAC (j + 5, j + 5);
    *FAC (j + 5, j + 5) = t;

    if (L3CCT (t) > small) {
	t = imsl_c_div (fpm1, t);
	if (m6 != j + 6)
	    ktemp = -ktemp;
    }
    else {
	ktemp = 0;
	info = m6;
    }
    /*
     * swap the subdiagonal elements of row j+6 with the pivot row in L
     */
    t1 = *FAC (j - 1, m6 - 1);
    *FAC (j - 1, m6 - 1) = *FAC (j - 1, j + 5);
    *FAC (j - 1, j + 5) = t1;
    t1 = *FAC (j, m6 - 1);
    *FAC (j, m6 - 1) = *FAC (j, j + 5);
    *FAC (j, j + 5) = t1;
    t1 = *FAC (j + 1, m6 - 1);
    *FAC (j + 1, m6 - 1) = *FAC (j + 1, j + 5);
    *FAC (j + 1, j + 5) = t1;
    t1 = *FAC (j + 2, m6 - 1);
    *FAC (j + 2, m6 - 1) = *FAC (j + 2, j + 5);
    *FAC (j + 2, j + 5) = t1;
    t1 = *FAC (j + 3, m6 - 1);
    *FAC (j + 3, m6 - 1) = *FAC (j + 3, j + 5);
    *FAC (j + 3, j + 5) = t1;
    t1 = *FAC (j + 4, m6 - 1);
    *FAC (j + 4, m6 - 1) = *FAC (j + 4, j + 5);
    *FAC (j + 4, j + 5) = t1;
    /* update U(j,j+7) to U(j+6,j+7) */
    t1 = *FAC (j + 6, m0 - 1);
    *FAC (j + 6, m0 - 1) = *FAC (j + 6, j - 1);
    *FAC (j + 6, j - 1) = t1;
    t2 = imsl_c_add (*FAC (j + 6, m1 - 1), imsl_c_mul (t1, *FAC (j - 1, j)));
    *FAC (j + 6, m1 - 1) = *FAC (j + 6, j);
    *FAC (j + 6, j) = t2;
    t3 = imsl_c_add (imsl_c_add (*FAC (j + 6, m2 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 1))),
	imsl_c_mul (t2, *FAC (j, j + 1)));
    *FAC (j + 6, m2 - 1) = *FAC (j + 6, j + 1);
    *FAC (j + 6, j + 1) = t3;
    t4 = imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 6, m3 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 2))),
	    imsl_c_mul (t2, *FAC (j, j + 2))), imsl_c_mul (t3, *FAC (j + 1, j + 2)));
    *FAC (j + 6, m3 - 1) = *FAC (j + 6, j + 2);
    *FAC (j + 6, j + 2) = t4;
    t5 = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 6, m4 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 3))),
		imsl_c_mul (t2, *FAC (j, j + 3))), imsl_c_mul (t3, *FAC (j + 1, j + 3))), imsl_c_mul (t4,
	    *FAC (j + 2, j + 3)));
    *FAC (j + 6, m4 - 1) = *FAC (j + 6, j + 3);
    *FAC (j + 6, j + 3) = t5;
    t6 = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 6, m5 - 1), imsl_c_mul (t1,
			    *FAC (j - 1, j + 4))), imsl_c_mul (t2, *FAC (j, j + 4))), imsl_c_mul (t3, *FAC (j + 1, j + 4))),
	    imsl_c_mul (t4, *FAC (j + 2, j + 4))), imsl_c_mul (t5, *FAC (j + 3, j + 4)));
    *FAC (j + 6, m5 - 1) = *FAC (j + 6, j + 4);
    *FAC (j + 6, j + 4) = t6;
    t7 = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 6, m6 - 1), imsl_c_mul (t1,
				*FAC (j - 1, j + 5))), imsl_c_mul (t2, *FAC (j, j + 5))), imsl_c_mul (t3, *FAC (j + 1, j + 5))),
		imsl_c_mul (t4, *FAC (j + 2, j + 5))), imsl_c_mul (t5, *FAC (j + 3, j + 5))), imsl_c_mul (t6,
	    *FAC (j + 4, j + 5)));
    *FAC (j + 6, m6 - 1) = *FAC (j + 6, j + 5);
    *FAC (j + 6, j + 5) = t7;
    /* update column j+6 of L and j+7 of A */
    for (i = j + 7; i <= *n; i++) {
	*FAC (j + 5, i - 1) = imsl_c_mul (*FAC (j + 5, i - 1), t);
	*FAC (j + 6, i - 1) = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + 6, i - 1),
				    imsl_c_mul (t1, *FAC (j - 1, i - 1))), imsl_c_mul (t2, *FAC (j, i - 1))), imsl_c_mul (t3,
				*FAC (j + 1, i - 1))), imsl_c_mul (t4, *FAC (j + 2, i - 1))), imsl_c_mul (t5,
			*FAC (j + 3, i - 1))), imsl_c_mul (t6, *FAC (j + 4, i - 1))), imsl_c_mul (t7,
		*FAC (j + 5, i - 1)));
    }

    if (j + 7 >= *n)
	goto L_110;
    /*
     * determine the pivot element for column j+7
     */
    _l0 = *n - j - 6;
    l_l4tcg (&_l0, FAC (j + 6, j + 6), &scale[j + 6], &m7);
    m7 += j + 6;
    ipvt[j + 6] = m7;
    t = scale[j + 6];
    scale[j + 6] = scale[m7 - 1];
    scale[m7 - 1] = t;
    /*
     * swap element j+7,j+7 with the pivot element
     */
    t = *FAC (j - 1, m7 - 1);
    *FAC (j - 1, m7 - 1) = *FAC (j - 1, j + 6);
    *FAC (j - 1, j + 6) = t;
    t = *FAC (j, m7 - 1);
    *FAC (j, m7 - 1) = *FAC (j, j + 6);
    *FAC (j, j + 6) = t;
    t = *FAC (j + 1, m7 - 1);
    *FAC (j + 1, m7 - 1) = *FAC (j + 1, j + 6);
    *FAC (j + 1, j + 6) = t;
    t = *FAC (j + 2, m7 - 1);
    *FAC (j + 2, m7 - 1) = *FAC (j + 2, j + 6);
    *FAC (j + 2, j + 6) = t;
    t = *FAC (j + 3, m7 - 1);
    *FAC (j + 3, m7 - 1) = *FAC (j + 3, j + 6);
    *FAC (j + 3, j + 6) = t;
    t = *FAC (j + 4, m7 - 1);
    *FAC (j + 4, m7 - 1) = *FAC (j + 4, j + 6);
    *FAC (j + 4, j + 6) = t;
    t = *FAC (j + 5, m7 - 1);
    *FAC (j + 5, m7 - 1) = *FAC (j + 5, j + 6);
    *FAC (j + 5, j + 6) = t;
    t = *FAC (j + 6, m7 - 1);
    *FAC (j + 6, m7 - 1) = *FAC (j + 6, j + 6);
    *FAC (j + 6, j + 6) = t;

    if (L3CCT (t) > small) {
	t = imsl_c_div (fpm1, t);
	for (i = j + 8; i <= *n; i++) {
	    *FAC (j + 6, i - 1) = imsl_c_mul (*FAC (j + 6, i - 1), t);
	}
	if (m7 != j + 7)
	    ktemp = -ktemp;
    }
    else {
	ktemp = 0;
	info = m7;
    }

    /*
     * VD$L NODEPCHK VD$L CNCALL VD$L NOSYNC
     */
    for (k = *n; k >= (j + 8); k--) {
	t1 = *FAC (k - 1, m0 - 1);
	*FAC (k - 1, m0 - 1) = *FAC (k - 1, j - 1);
	*FAC (k - 1, j - 1) = t1;
	t2 = imsl_c_add (*FAC (k - 1, m1 - 1), imsl_c_mul (t1, *FAC (j - 1, j)));
	*FAC (k - 1, m1 - 1) = *FAC (k - 1, j);
	*FAC (k - 1, j) = t2;
	t3 = imsl_c_add (imsl_c_add (*FAC (k - 1, m2 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 1))),
	    imsl_c_mul (t2, *FAC (j, j + 1)));
	*FAC (k - 1, m2 - 1) = *FAC (k - 1, j + 1);
	*FAC (k - 1, j + 1) = t3;
	t4 = imsl_c_add (imsl_c_add (imsl_c_add (*FAC (k - 1, m3 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 2))),
		imsl_c_mul (t2, *FAC (j, j + 2))), imsl_c_mul (t3, *FAC (j + 1, j + 2)));
	*FAC (k - 1, m3 - 1) = *FAC (k - 1, j + 2);
	*FAC (k - 1, j + 2) = t4;
	t5 = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (k - 1, m4 - 1), imsl_c_mul (t1, *FAC (j - 1, j + 3))),
		    imsl_c_mul (t2, *FAC (j, j + 3))), imsl_c_mul (t3, *FAC (j + 1, j + 3))), imsl_c_mul (t4,
		*FAC (j + 2, j + 3)));
	*FAC (k - 1, m4 - 1) = *FAC (k - 1, j + 3);
	*FAC (k - 1, j + 3) = t5;
	t6 = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (k - 1, m5 - 1), imsl_c_mul (t1,
				*FAC (j - 1, j + 4))), imsl_c_mul (t2, *FAC (j, j + 4))), imsl_c_mul (t3, *FAC (j + 1, j + 4))),
		imsl_c_mul (t4, *FAC (j + 2, j + 4))), imsl_c_mul (t5, *FAC (j + 3, j + 4)));
	*FAC (k - 1, m5 - 1) = *FAC (k - 1, j + 4);
	*FAC (k - 1, j + 4) = t6;
	t7 = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (k - 1, m6 - 1),
				imsl_c_mul (t1, *FAC (j - 1, j + 5))), imsl_c_mul (t2, *FAC (j, j + 5))), imsl_c_mul (t3,
			    *FAC (j + 1, j + 5))), imsl_c_mul (t4, *FAC (j + 2, j + 5))), imsl_c_mul (t5,
		    *FAC (j + 3, j + 5))), imsl_c_mul (t6, *FAC (j + 4, j + 5)));
	*FAC (k - 1, m6 - 1) = *FAC (k - 1, j + 5);
	*FAC (k - 1, j + 5) = t7;
	t8 = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (k - 1, m7 - 1),
				    imsl_c_mul (t1, *FAC (j - 1, j + 6))), imsl_c_mul (t2, *FAC (j, j + 6))), imsl_c_mul (t3,
				*FAC (j + 1, j + 6))), imsl_c_mul (t4, *FAC (j + 2, j + 6))), imsl_c_mul (t5,
			*FAC (j + 3, j + 6))), imsl_c_mul (t6, *FAC (j + 4, j + 6))), imsl_c_mul (t7,
		*FAC (j + 5, j + 6)));
	*FAC (k - 1, m7 - 1) = *FAC (k - 1, j + 6);
	*FAC (k - 1, j + 6) = t8;
	/*
	 * rank 8 update of the lower right block from rows j+8 to n and
	 * columns j+8 to n
	 */

	for (i = j + 8; i <= *n; i++) {
#if 0
	    *FAC (k - 1, i - 1) = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (k - 1, i - 1),
					    imsl_c_mul (t1, *FAC (j - 1, i - 1))), imsl_c_mul (t2, *FAC (j, i - 1))),
				    imsl_c_mul (t3, *FAC (j + 1, i - 1))), imsl_c_mul (t4, *FAC (j + 2, i - 1))),
			    imsl_c_mul (t5, *FAC (j + 3, i - 1))), imsl_c_mul (t6, *FAC (j + 4, i - 1))),
		    imsl_c_mul (t7, *FAC (j + 5, i - 1))), imsl_c_mul (t8, *FAC (j + 6, i - 1)));
#endif

	    ctemp1.re = t1.re*(FAC (j-1, i-1))->re - t1.im*(FAC (j-1, i-1))->im;
	    ctemp1.im = t1.re*(FAC (j-1, i-1))->im + t1.im*(FAC (j-1, i-1))->re;
		
	    ctemp2.re = t2.re*(FAC (j, i-1))->re - t2.im*(FAC (j, i-1))->im;
	    ctemp2.im = t2.re*(FAC (j, i-1))->im + t2.im*(FAC (j, i-1))->re;
		
	    ctemp3.re = t3.re*(FAC (j+1, i-1))->re - t3.im*(FAC (j+1, i-1))->im;
	    ctemp3.im = t3.re*(FAC (j+1, i-1))->im + t3.im*(FAC (j+1, i-1))->re;
		
	    ctemp4.re = t4.re*(FAC (j+2, i-1))->re - t4.im*(FAC (j+2, i-1))->im;
	    ctemp4.im = t4.re*(FAC (j+2, i-1))->im + t4.im*(FAC (j+2, i-1))->re;
		
	    ctemp5.re = t5.re*(FAC (j+3, i-1))->re - t5.im*(FAC (j+3, i-1))->im;
	    ctemp5.im = t5.re*(FAC (j+3, i-1))->im + t5.im*(FAC (j+3, i-1))->re;
		
	    ctemp6.re = t6.re*(FAC (j+4, i-1))->re - t6.im*(FAC (j+4, i-1))->im;
	    ctemp6.im = t6.re*(FAC (j+4, i-1))->im + t6.im*(FAC (j+4, i-1))->re;
		
	    ctemp7.re = t7.re*(FAC (j+5, i-1))->re - t7.im*(FAC (j+5, i-1))->im;
	    ctemp7.im = t7.re*(FAC (j+5, i-1))->im + t7.im*(FAC (j+5, i-1))->re;
		
	    ctemp8.re = t8.re*(FAC (j+6, i-1))->re - t8.im*(FAC (j+6, i-1))->im;
	    ctemp8.im = t8.re*(FAC (j+6, i-1))->im + t8.im*(FAC (j+6, i-1))->re;
		
	    (FAC (k-1, i-1))->re += ctemp1.re + ctemp2.re + ctemp3.re +
	    		ctemp4.re + ctemp5.re + ctemp6.re + ctemp7.re +
			ctemp8.re;
	    (FAC (k-1, i-1))->im += ctemp1.im + ctemp2.im + ctemp3.im +
	    		ctemp4.im + ctemp5.im + ctemp6.im + ctemp7.im +
			ctemp8.im;
	}
    }
    /*
     * swap elements in L to get FAC into IMSL factored form
     */

    t8 = *FAC (j - 1, m7 - 1);
    *FAC (j - 1, m7 - 1) = *FAC (j - 1, j + 6);
    *FAC (j - 1, j + 6) = t8;
    t8 = *FAC (j, m7 - 1);
    *FAC (j, m7 - 1) = *FAC (j, j + 6);
    *FAC (j, j + 6) = t8;
    t8 = *FAC (j + 1, m7 - 1);
    *FAC (j + 1, m7 - 1) = *FAC (j + 1, j + 6);
    *FAC (j + 1, j + 6) = t8;
    t8 = *FAC (j + 2, m7 - 1);
    *FAC (j + 2, m7 - 1) = *FAC (j + 2, j + 6);
    *FAC (j + 2, j + 6) = t8;
    t8 = *FAC (j + 3, m7 - 1);
    *FAC (j + 3, m7 - 1) = *FAC (j + 3, j + 6);
    *FAC (j + 3, j + 6) = t8;
    t8 = *FAC (j + 4, m7 - 1);
    *FAC (j + 4, m7 - 1) = *FAC (j + 4, j + 6);
    *FAC (j + 4, j + 6) = t8;
    t8 = *FAC (j + 5, m7 - 1);
    *FAC (j + 5, m7 - 1) = *FAC (j + 5, j + 6);
    *FAC (j + 5, j + 6) = t8;

L_110:
    t9 = *FAC (j - 1, m6 - 1);
    *FAC (j - 1, m6 - 1) = *FAC (j - 1, j + 5);
    *FAC (j - 1, j + 5) = t9;
    t9 = *FAC (j, m6 - 1);
    *FAC (j, m6 - 1) = *FAC (j, j + 5);
    *FAC (j, j + 5) = t9;
    t9 = *FAC (j + 1, m6 - 1);
    *FAC (j + 1, m6 - 1) = *FAC (j + 1, j + 5);
    *FAC (j + 1, j + 5) = t9;
    t9 = *FAC (j + 2, m6 - 1);
    *FAC (j + 2, m6 - 1) = *FAC (j + 2, j + 5);
    *FAC (j + 2, j + 5) = t9;
    t9 = *FAC (j + 3, m6 - 1);
    *FAC (j + 3, m6 - 1) = *FAC (j + 3, j + 5);
    *FAC (j + 3, j + 5) = t9;
    t9 = *FAC (j + 4, m6 - 1);
    *FAC (j + 4, m6 - 1) = *FAC (j + 4, j + 5);
    *FAC (j + 4, j + 5) = t9;

L_120:
    t0 = *FAC (j - 1, m5 - 1);
    *FAC (j - 1, m5 - 1) = *FAC (j - 1, j + 4);
    *FAC (j - 1, j + 4) = t0;
    t0 = *FAC (j, m5 - 1);
    *FAC (j, m5 - 1) = *FAC (j, j + 4);
    *FAC (j, j + 4) = t0;
    t0 = *FAC (j + 1, m5 - 1);
    *FAC (j + 1, m5 - 1) = *FAC (j + 1, j + 4);
    *FAC (j + 1, j + 4) = t0;
    t0 = *FAC (j + 2, m5 - 1);
    *FAC (j + 2, m5 - 1) = *FAC (j + 2, j + 4);
    *FAC (j + 2, j + 4) = t0;
    t0 = *FAC (j + 3, m5 - 1);
    *FAC (j + 3, m5 - 1) = *FAC (j + 3, j + 4);
    *FAC (j + 3, j + 4) = t0;

L_130:
    t = *FAC (j - 1, m4 - 1);
    *FAC (j - 1, m4 - 1) = *FAC (j - 1, j + 3);
    *FAC (j - 1, j + 3) = t;
    t = *FAC (j, m4 - 1);
    *FAC (j, m4 - 1) = *FAC (j, j + 3);
    *FAC (j, j + 3) = t;
    t = *FAC (j + 1, m4 - 1);
    *FAC (j + 1, m4 - 1) = *FAC (j + 1, j + 3);
    *FAC (j + 1, j + 3) = t;
    t = *FAC (j + 2, m4 - 1);
    *FAC (j + 2, m4 - 1) = *FAC (j + 2, j + 3);
    *FAC (j + 2, j + 3) = t;

L_140:
    t1 = *FAC (j - 1, m3 - 1);
    *FAC (j - 1, m3 - 1) = *FAC (j - 1, j + 2);
    *FAC (j - 1, j + 2) = t1;
    t1 = *FAC (j, m3 - 1);
    *FAC (j, m3 - 1) = *FAC (j, j + 2);
    *FAC (j, j + 2) = t1;
    t1 = *FAC (j + 1, m3 - 1);
    *FAC (j + 1, m3 - 1) = *FAC (j + 1, j + 2);
    *FAC (j + 1, j + 2) = t1;

L_150:
    t2 = *FAC (j - 1, m2 - 1);
    *FAC (j - 1, m2 - 1) = *FAC (j - 1, j + 1);
    *FAC (j - 1, j + 1) = t2;
    t2 = *FAC (j, m2 - 1);
    *FAC (j, m2 - 1) = *FAC (j, j + 1);
    *FAC (j, j + 1) = t2;

L_160:
    t3 = *FAC (j - 1, m1 - 1);
    *FAC (j - 1, m1 - 1) = *FAC (j - 1, j);
    *FAC (j - 1, j) = t3;
    j += 8;
    goto L_20;

L_170:
    ipvt[*n - 1] = *n;
    if (L3CCT (*FAC (*n - 1, *n - 1)) <= small)
	info = *n;

    if (info != 0) {

	/*
	 * (4, 2, "The input matrix is singular.  Some of the diagonal
	 * elements of the upper triangular matrix U of the LU factorization
	 * are close to zero.");
	 */
	imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_MATRIX);
    }
L_9900:
    imsl_e1pop ("l_l2tcg");
    return;
#undef	L3CCT
#undef  A
#undef  FAC
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  L4TCG/DL4TCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 9, 1989

    Purpose:    Find the pivot element using scaling.

    Usage:      CALL L4TCG (N, X, SCALE, K)

    Arguments:
       N      - Number of elements in the row.  (Input)
       X      - Row of the matrix.  (Input)
       SCALE  - Norms of the rows.  (Input)
       K      - Value of I for which ABS(X(I))*SCALE(I) is largest.
               (Output)

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	NTEMP	400
#ifdef ANSI
static void l_l4tcg (Mint *n, Mf_complex *x, Mf_complex *scale, Mint *k)
#else
static void l_l4tcg (n, x, scale, k)
    Mint       *n;
    Mf_complex  x[], scale[];
    Mint       *k;
#endif
{
    Mint        _l0, _l1, j, j1, j2, jmax;
    Mfloat      curmax, value;
    Mf_complex  temp[NTEMP];


#define L3CCT(zdum)	(Mfloat)(fabs( imsl_fc_convert( (zdum) ) ) + fabs( imsl_c_aimag( (zdum) ) ))

    *k = 1;
    curmax = F_ZERO;
    for (j1 = 1; j1 <= *n; j1 += NTEMP) {
	j2 = imsl_i_min (j1 + NTEMP - 1, *n);
	for (j = j1; j <= j2; j++) {
#if 0
	    temp[j - j1] = imsl_c_mul (x[j - 1], scale[j - 1]);
#endif
	    temp[j-j1].re = x[j-1].re*scale[j-1].re - x[j-1].im*scale[j-1].im;
	    temp[j-j1].im = x[j-1].re*scale[j-1].im + x[j-1].im*scale[j-1].re;
	}
	_l0 = j2 - j1 + 1;
	_l1 = 1;
	jmax = imsl_icamax (&_l0, temp, &_l1);
	value = L3CCT (temp[jmax - 1]);
	if (value > curmax) {
	    curmax = value;
	    *k = jmax + j1 - 1;
	}
    }

    return;
#undef	L3CCT
#undef  NTEMP
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  L2CCG/DL2CCG  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 13, 1987

    Purpose:    Compute the LU factorization of a d_complex general matrix
                and estimate its L1 condition number.

    Usage:      CALL L2CCG (N, A, LDA, FAC, LDFAC, IPVT, RCOND, Z)

    Arguments:  See LFCCG/DLFCCG.

    Remarks:    See LFCCG/DLFCCG.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_l2ccg (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *imsl_fac,
                Mint *ldfac, Mint *ipvt, Mfloat *rcond, Mf_complex *z)
#else
static void l_l2ccg (n, a, lda, imsl_fac, ldfac, ipvt, rcond, z)
    Mint       *n;
    Mf_complex *a;
    Mint       *lda;
    Mf_complex *imsl_fac;
    Mint       *ldfac, ipvt[];
    Mfloat     *rcond;
    Mf_complex  z[];
#endif
{
#define A(I_,J_)	(a+(I_)*(alda)+(J_))
#define FAC(I_,J_)	(imsl_fac+(I_)*(aldfac)+(J_))
    Mint        alda = *lda;
    Mint        aldfac = *ldfac;

    Mint        _l0, _l1, _l2, j, k, kp1, l;
    Mfloat      anorm, s, sm, ynorm;
    Mf_complex  ek, t, wk, wkm, zero;


#define L3CCT(zdum)	(Mfloat)(fabs( imsl_fc_convert( (zdum) ) ) + fabs( imsl_c_aimag( (zdum) ) ))
#define L4CCT(zdum1,zdum2)	(imsl_c_mul(imsl_cf_convert(L3CCT( (zdum1) ),0.),\
	 (imsl_c_div((zdum2),imsl_cf_convert(L3CCT( (zdum2) ),F_ZERO)))))

    imsl_e1psh ("l_l2ccg");

    if (*n <= 0) {
	imsl_e1sti (1, *n);

	/*
	 * (5, 1, "The order of the matrix must be positive while N = %(i1)
	 * is given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
	goto L_9000;
    }
    if (*n > *lda) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *lda);

	/*
	 * (5, 2, "The order of the matrix must be less than or equal to its
	 * leading dimension while N = %(i1) and LDA = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
	goto L_9000;
    }
    if (*n > *ldfac) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *ldfac);

	/*
	 * (5, 3, "The order of the matrix must be less than or equal to its
	 * leading dimension while N = %(i1) and LDFAC = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_LDFAC_LESS_ORDER);
	goto L_9000;
    }
    /* COMPUTE 1-NORM OF A */
    *rcond = F_ZERO;
    anorm = F_ZERO;
    for (j = 1; j <= *n; j++) {
	_l0 = 1;
	anorm = imsl_f_max (anorm, imsl_scasum (n, A (j - 1, 0), &_l0));
    }
    /*
     * FACTORIZATION STEP
     */
    l_l2tcg (n, a, lda, imsl_fac, ldfac, ipvt, z);
    if (imsl_n1rty (1) == 4)
	goto L_9000;
    /*
     * RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))).  ESTIMATE =
     * NORM(Z)/NORM(Y) WHERE A*Z = Y AND CTRANS(A)*Y = E.  CTRANS(A) IS THE
     * TRANSPOSE OF A.  THE COMPONENTS OF E ARE CHOSEN TO CAUSE MAXIMUM LOCAL
     * GROWTH IN THE ELEMENTS OF W WHERE CTRANS(U)*W = E. THE VECTORS ARE
     * FREQUENTLY RESCALED TO AVOID OVERFLOW. SOLVE CTRANS(U)*W = E
     */
    zero = imsl_cf_convert (F_ZERO, F_ZERO);
    ek = imsl_cf_convert (F_ONE, F_ZERO);
    _l0 = 1;
    imsl_cset (n, &zero, z, &_l0);
    for (k = 1; k <= *n; k++) {
	wk = imsl_c_neg (z[k - 1]);
	if (L3CCT (z[k - 1]) != F_ZERO)
	    ek = L4CCT (ek, wk);
	if (L3CCT (imsl_c_sub (ek, z[k - 1])) > L3CCT (*FAC (k - 1, k - 1))) {
	    s = L3CCT (*FAC (k - 1, k - 1)) / L3CCT (imsl_c_sub (ek, z[k - 1]));
	    _l0 = 1;
	    imsl_csscal (n, &s, z, &_l0);
	    ek = imsl_c_mul (imsl_cf_convert (s, F_ZERO), ek);
	}
	wk = imsl_c_sub (ek, z[k - 1]);
	wkm = imsl_c_sub (imsl_c_neg (ek), z[k - 1]);
	s = L3CCT (wk);
	sm = L3CCT (wkm);
	if (L3CCT (*FAC (k - 1, k - 1)) != F_ZERO) {
	    wk = imsl_c_div (wk, imsl_c_conjg (*FAC (k - 1, k - 1)));
	    wkm = imsl_c_div (wkm, imsl_c_conjg (*FAC (k - 1, k - 1)));
	}
	else {
	    wk = imsl_cf_convert (F_ONE, F_ZERO);
	    wkm = imsl_cf_convert (F_ONE, F_ZERO);
	}
	kp1 = k + 1;
	if (kp1 <= *n) {
	    for (j = kp1; j <= *n; j++) {
		sm += L3CCT (imsl_c_add (z[j - 1], imsl_c_mul (wkm, imsl_c_conjg (*FAC (j - 1, k - 1)))));
		z[j - 1] = imsl_c_add (z[j - 1], imsl_c_mul (wk, imsl_c_conjg (*FAC (j - 1, k - 1))));
		s += L3CCT (z[j - 1]);
	    }
	    if (s < sm) {
		t = imsl_c_sub (wkm, wk);
		wk = wkm;
		for (j = kp1; j <= *n; j++) {
		    z[j - 1] = imsl_c_add (z[j - 1], imsl_c_mul (t, imsl_c_conjg (*FAC (j - 1, k - 1))));
		}
	    }
	}
	z[k - 1] = wk;
    }
    _l0 = 1;
    s = F_ONE / imsl_scasum (n, z, &_l0);
    _l0 = 1;
    imsl_csscal (n, &s, z, &_l0);
    /* SOLVE CTRANS(L)*Y = W */
    for (k = *n; k >= 1; k--) {
	_l0 = *n - k;
	_l1 = 1;
	_l2 = 1;
	if (k < *n)
	    z[k - 1] = imsl_c_add (z[k - 1], imsl_cdotc (&_l0, FAC (k - 1, k),
		    &_l1, &z[k], &_l2));
	if (L3CCT (z[k - 1]) > F_ONE) {
	    s = F_ONE / L3CCT (z[k - 1]);
	    _l0 = 1;
	    imsl_csscal (n, &s, z, &_l0);
	}
	l = ipvt[k - 1];
	t = z[l - 1];
	z[l - 1] = z[k - 1];
	z[k - 1] = t;
    }
    _l0 = 1;
    s = F_ONE / imsl_scasum (n, z, &_l0);
    _l0 = 1;
    imsl_csscal (n, &s, z, &_l0);

    ynorm = F_ONE;
    /* SOLVE L*V = Y */
    for (k = 1; k <= *n; k++) {
	l = ipvt[k - 1];
	t = z[l - 1];
	z[l - 1] = z[k - 1];
	z[k - 1] = t;
	_l0 = *n - k;
	_l1 = 1;
	_l2 = 1;
	if (k < *n)
	    imsl_caxpy (&_l0, &t, FAC (k - 1, k), &_l1,
		&z[k], &_l2);
	if (L3CCT (z[k - 1]) > F_ONE) {
	    s = F_ONE / L3CCT (z[k - 1]);
	    _l0 = 1;
	    imsl_csscal (n, &s, z, &_l0);
	    ynorm *= s;
	}
    }
    _l0 = 1;
    s = F_ONE / imsl_scasum (n, z, &_l0);
    _l0 = 1;
    imsl_csscal (n, &s, z, &_l0);
    ynorm *= s;
    /* SOLVE U*Z = V */
    for (k = *n; k >= 1; k--) {
	if (L3CCT (z[k - 1]) > L3CCT (*FAC (k - 1, k - 1))) {
	    s = L3CCT (*FAC (k - 1, k - 1)) / L3CCT (z[k - 1]);
	    _l0 = 1;
	    imsl_csscal (n, &s, z, &_l0);
	    ynorm *= s;
	}
	if (L3CCT (*FAC (k - 1, k - 1)) != F_ZERO) {
	    z[k - 1] = imsl_c_div (z[k - 1], *FAC (k - 1, k - 1));
	}
	else {
	    z[k - 1] = imsl_cf_convert (F_ONE, F_ZERO);
	}
	t = imsl_c_neg (z[k - 1]);
	_l0 = k - 1;
	_l1 = 1;
	_l2 = 1;
	imsl_caxpy (&_l0, &t, FAC (k - 1, 0), &_l1, &z[0],
	    &_l2);
    }
    /* MAKE ZNORM = 1.0 */
    _l0 = 1;
    s = F_ONE / imsl_scasum (n, z, &_l0);
    _l0 = 1;
    imsl_csscal (n, &s, z, &_l0);
    ynorm *= s;

    if (anorm != F_ZERO)
	*rcond = ynorm / anorm;
    if (*rcond <= imsl_amach (4)) {
	imsl_e1str (1, *rcond);

	/*
	 * (3, 1, "The matrix is algorithmically singular.  An estimate of
	 * the reciprocal of its L1 condition number is RCOND = %(r1).");
	 */
	imsl_ermes (IMSL_WARNING, IMSL_ILL_CONDITIONED);
    }
L_9000:
    imsl_e1pop ("l_l2ccg");
    return;
#undef	L4CCT
#undef	L3CCT
#undef  A
#undef  FAC
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  L2NCG/DL2NCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 27, 1985

    Purpose:    Compute the inverse of a d_complex general matrix.

    Usage:      CALL L2NCG (N, A, LDA, AINV, LDAINV, WK, IWK)

    Arguments:  See LINCG/DLINCG.

    Remarks:    See LINCG/DLINCG.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_l2ncg (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *ainv,
                Mint *ldainv, Mf_complex *wk, Mint *iwk)
#else
static void l_l2ncg (n, a, lda, ainv, ldainv, wk, iwk)
    Mint       *n;
    Mf_complex *a;
    Mint       *lda;
    Mf_complex *ainv;
    Mint       *ldainv;
    Mf_complex  wk[];
    Mint        iwk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(alda)+(J_))
#define AINV(I_,J_)	(ainv+(I_)*(aldainv)+(J_))
    Mint        aldainv = *ldainv;
    Mint        _l0, _l1, _l2, i, inc, j, k, l;
    Mf_complex  _cx0, _cx1;


    imsl_e1psh ("l_l2ncg");

    if (*n <= 0) {
	imsl_e1sti (1, *n);

	/*
	 * (5, 1, "The order of the matrix must be positive while N = %(i1)
	 * is given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
	goto L_9000;
    }
    if (*n > *lda) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *lda);

	/*
	 * (5, 2, "The order of the matrix must be less than or equal to its
	 * leading dimension while N = %(i1) and LDA = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
	goto L_9000;
    }
    if (*n > *ldainv) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *ldainv);

	/*
	 * (5, 3, "The order of the matrix must be less than or equal to its
	 * leading dimension while N = %(i1) and LDAINV = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_LDAINV_LESS_ORDER);
	goto L_9000;
    }
    /*
     * COMPUTE THE LU FACTORIZATION OF A AND ESTIMATE ITS CONDITION NUMBER
     */
    inc = *n * (*n - 1) / 2;
    l_l2ccg (n, a, lda, ainv, ldainv, iwk, &lv_rcond, &wk[inc]);
    if (imsl_n1rty (1) == 4)
	goto L_9000;

    j = inc;
    k = 0;
    for (i = 1; i <= (*n - 1); i++) {
	j -= k;
	_l0 = 1;
	_l1 = 1;
	imsl_ccopy (&i, AINV (*n - i - 1, *n - i), &_l0, &wk[j - 1],
	    &_l1);
	k = i + 1;
    }
    /* COMPUTE INVERSE(U) */
    _l0 = 2;
    l_linct (n, ainv, ldainv, &_l0, ainv, ldainv);
    j = inc;
    k = 0;
    for (i = 1; i <= (*n - 1); i++) {
	j -= k;
	_l0 = 1;
	_l1 = 1;
	imsl_ccopy (&i, &wk[j - 1], &_l0, AINV (*n - i - 1, *n - i),
	    &_l1);
	k = i + 1;
    }
    /* FORM INVERSE(U)*INVERSE(L) */
    for (k = *n - 1; k >= 1; k--) {
	_l0 = *n - k;
	_l1 = 1;
	_l2 = 1;
	imsl_ccopy (&_l0, AINV (k - 1, k), &_l1, &wk[k + inc], &_l2);
	_l0 = *n - k;
	_cx0 = imsl_cf_convert (F_ZERO, F_ZERO);
	_l1 = 1;
	imsl_cset (&_l0, &_cx0, AINV (k - 1, k), &_l1);
	_l0 = *n - k;
	_cx0 = imsl_cf_convert (F_ONE, F_ZERO);
	_l1 = 1;
	_cx1 = imsl_cf_convert (F_ONE, F_ZERO);
	_l2 = 1;
	imsl_cgemv ("N", sizeof ("N"), n, &_l0, &_cx0,
	    AINV (k, 0), ldainv, &wk[k + inc], &_l1, &_cx1,
	    AINV (k - 1, 0), &_l2);
	l = iwk[k - 1];
	_l0 = 1;
	_l1 = 1;
	if (l != k)
	    imsl_cswap (n, AINV (k - 1, 0), &_l0, AINV (l - 1, 0), &_l1);
    }

    if (lv_rcond <= imsl_amach (4)) {
	imsl_e1str (1, lv_rcond);

	/*
	 * (3, 1, "The matrix is too ill-conditioned. An estimate of the
	 * reciprocal of its L1 condition number is RCOND = %(r1).  The
	 * inverse might not be accurate.");
	 */
	imsl_ermes (IMSL_WARNING, IMSL_ILL_CONDITIONED);
    }
L_9000:
    imsl_e1pop ("l_l2ncg");
    return;
}				/* end of function */
#undef  A
#undef  AINV
/*----------------------------------------------------------------------- */

/*  IMSL Name:  LINCT/DLINCT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 25, 1987

    Purpose:    Compute the inverse of a d_complex triangular matrix.

    Usage:      CALL LINCT (N, A, LDA, IPATH, AINV, LDAINV)

    Arguments:
       N      - Order of the matrix.  (Input)
       A      - Complex N by N matrix containing the triangular matrix
                to be inverted.  (Input)
                For a lower triangular matrix, only the lower triangle
                of A is referenced.  For an upper triangular matrix, only
                the upper triangle of A is referenced.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPATH  - Path indicator.  (Input)
                IPATH = 1 means A is lower triangular,
                IPATH = 2 means A is upper triangular.
       AINV   - Complex N by N matrix containing the inverse of A.
                (Output)
                If A is lower triangular, AINV is also lower triangular.
                If A is upper triangular, AINV is also upper triangular.
                If A is not needed, A and AINV can share the same storage
                locations.
       LDAINV - Leading dimension of AINV exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remark:
       Informational error
       Type Code
         4   1  The input triangular matrix is singular.  Some of its
                diagonal elements are close to zero.

    GAMS:       D2c3

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_linct (Mint *n, Mf_complex *a, Mint *lda, Mint *ipath,
                Mf_complex *ainv, Mint *ldainv)
#else
static void l_linct (n, a, lda, ipath, ainv, ldainv)
    Mint       *n;
    Mf_complex *a;
    Mint       *lda, *ipath;
    Mf_complex *ainv;
    Mint       *ldainv;
#endif
{
#define A(I_,J_)	(a+(I_)*(alda)+(J_))
#define AINV(I_,J_)	(ainv+(I_)*(aldainv)+(J_))
    Mint        alda = *lda;
    Mint        aldainv = *ldainv;
    Mint        _l0, _l1, _l2, info, j, k;
    Mfloat      big, small;
    Mf_complex  _cx0, cone, czero, temp;


#define L3CCT(zdum)	(Mfloat)(fabs( imsl_fc_convert( (zdum) ) ) + fabs( imsl_c_aimag( (zdum) ) ))

    imsl_e1psh ("l_linct");

    if (*n <= 0) {
	imsl_e1sti (1, *n);

	/*
	 * (5, 1, "The order of the matrix must be positive while N = %(i1)
	 * is given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
    }
    else {
	if (*n > *lda) {
	    imsl_e1sti (1, *n);
	    imsl_e1sti (2, *lda);

	    /*
	     * (5, 2, "The order of the matrix must be less than or equal to
	     * its leading dimension while N = %(i1) and LDA = %(i2) are
	     * given.");
	     */
	    imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
	}
	if (*n > *ldainv) {
	    imsl_e1sti (1, *n);
	    imsl_e1sti (2, *ldainv);

	    /*
	     * (5, 3, "The order of the matrix must be less than or equal to
	     * its leading dimension while N = %(i1) and LDAINV = %(i2) are
	     * given.");
	     */
	    imsl_ermes (IMSL_TERMINAL, IMSL_LDAINV_LESS_ORDER);
	}
	if (*ipath != 1 && *ipath != 2) {
	    imsl_e1sti (1, *ipath);

	    /*
	     * (5, 4, "IPATH must be either 1 or 2 while a value of %(i1) is
	     * given.");
	     */
	    imsl_ermes (IMSL_TERMINAL, IMSL_IPATH_RANGE_3);
	}
	if (imsl_n1rcd (0) == 0) {
	    cone = imsl_cf_convert (F_ONE, F_ZERO);
	    czero = imsl_cf_convert (F_ZERO, F_ZERO);
	    small = imsl_amach (1);
	    big = imsl_amach (2);
	    if (small * big < F_ONE)
		small = F_ONE / big;

	    if (*ipath == 1) {
		/*
		 * MAKE A COPY OF A IN AINV AND WORK WITH AINV
		 */
		for (j = 1; j <= *n; j++) {
		    _l0 = *n - j + 1;
		    _l1 = 1;
		    _l2 = 1;
		    imsl_ccopy (&_l0, A (j - 1, j - 1), &_l1,
			AINV (j - 1, j - 1), &_l2);
		    _l0 = j - 1;
		    _l1 = 1;
		    imsl_cset (&_l0, &czero, AINV (j - 1, 0), &_l1);
		}
		/*
		 * COMPUTE INVERSE OF LOWER TRIANGULAR MATRIX
		 */
		for (k = *n; k >= 1; k--) {
		    info = k;
		    if (L3CCT (*AINV (k - 1, k - 1)) <= small)
			goto L_50;
		    *AINV (k - 1, k - 1) = imsl_c_div (cone, *AINV (k - 1, k - 1));
		    temp = imsl_c_neg (*AINV (k - 1, k - 1));
		    if (k < *n) {
			_l0 = *n - k;
			_l1 = 1;
			imsl_cscal (&_l0, &temp, AINV (k - 1, k),
			    &_l1);
			_l0 = *n - k;
			_l1 = k - 1;
			_cx0 = imsl_cf_convert (F_ONE, F_ZERO);
			_l2 = 1;
			imsl_cgeru (&_l0, &_l1, &_cx0, AINV (k - 1, k), &_l2, AINV (0, k - 1),
			    ldainv, AINV (0, k), ldainv);
		    }
		    _l0 = k - 1;
		    imsl_cscal (&_l0, AINV (k - 1, k - 1), AINV (0, k - 1),
			ldainv);
		}
		info = 0;
	    }
	    else {
		/*
		 * MAKE A COPY OF A IN AINV AND WORK WITH AINV
		 */
		for (j = 1; j <= *n; j++) {
		    _l0 = 1;
		    _l1 = 1;
		    imsl_ccopy (&j, A (j - 1, 0), &_l0, AINV (j - 1, 0),
			&_l1);
		    _l0 = *n - j;
		    _l1 = 1;
		    if (j < *n)
			imsl_cset (&_l0, &czero, AINV (j - 1, j),
			    &_l1);
		}
		/*
		 * COMPUTE INVERSE OF AN UPPER TRIANGULAR MATRIX
		 */
		for (k = 1; k <= *n; k++) {
		    info = k;
		    if (L3CCT (*AINV (k - 1, k - 1)) <= small)
			goto L_50;
		    *AINV (k - 1, k - 1) = imsl_c_div (cone, *AINV (k - 1, k - 1));
		    temp = imsl_c_neg (*AINV (k - 1, k - 1));
		    _l0 = k - 1;
		    _l1 = 1;
		    imsl_cscal (&_l0, &temp, AINV (k - 1, 0), &_l1);
		    _l0 = k - 1;
		    _l1 = *n - k;
		    _cx0 = imsl_cf_convert (F_ONE, F_ZERO);
		    _l2 = 1;
		    if (k < *n) {
			imsl_cgeru (&_l0, &_l1, &_cx0, AINV (k - 1, 0), &_l2, AINV (k, k - 1),
			    ldainv, AINV (k, 0), ldainv);
			_l0 = *n - k;
			imsl_cscal (&_l0, AINV (k - 1, k - 1),
			    AINV (k, k - 1), ldainv);
		    }
		}
		info = 0;
	    }

    L_50:
	    if (info != 0) {
/*		imsl_e1sti (1, info);  */

		/*
		 * (4, 1, "The matrix to be inverted is singular.  The index
		 * of the first zero diagonal element of A is %(i1).");
		 */
		imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_MATRIX);
	    }
	}
    }

    imsl_e1pop ("l_linct");
    return;
#undef	L3CCT
}				/* end of function */
#undef  A
#undef  AINV
/* Structured by FOR_STRUCT, v0.2, on 08/28/90 at 14:52:32
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CTRSV  (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 26, 1989

    Purpose:    Solve one of the triangular systems,
                    x = inv(A)*x,
                    x = inv(trans(A))*x,
                 or
                    x = inv(ctrans(A))*x,
                where A is a triangular matrix, trans(A) is the transpose
                of the matrix, and ctrans(A) is the conjugate transpose
                of the matrix.

    Usage:      CALL CTRSV (UPLO, TRANS, DIAG, N, A, LDA, X, INCX)

    Arguments:
       UPLO   - Character specifing the storage structure.
                (Input)
                   UPLO              Structure
                'U' or 'u'      Matrix is upper triangular
                'L' or 'l'      Matrix is lower triangular
       TRANS  - Character specifing if the transpose solution is to be
                computed.  (Input)
                   TRANS              Meaning
                'N' or 'n'      Compute x = inv(A)*x
                'T' or 't'      Compute x = inv(trans(A))*x
                'C' or 'c'      Compute x = inv(ctrans(A))*x
       DIAG   - Character specifing if A is unit triangular.
                (Input)
                   DIAG               Meaning
                'U' or 'u'      A is assumed to be unit triangular.
                'N' or 'n'      A is not assumed to be unit triangular.
                If DIAG is 'U' or 'u' then the diagonal elements
                of A are assumed to be one and are not referenced.  If
                DIAG is 'N' or 'n' then the actual diagonal elements of
                array A are used.  (Input)
       N      - Order of the matrix A.  (Input)
       A      - Complex triangular matrix of order N.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                calling routine.  (Input)
       X      - Complex vector of length (N-1)*IABS(INCX)+1.
                (Input/Output)
       INCX   - Displacement between elements of X.  (Input)

    GAMS:       D1b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_ctrsv (Mchar *uplo, unsigned uplo_s, Mchar *trans, unsigned trans_s,
                Mchar *diag, unsigned diag_s, Mint *n, Mf_complex *a,
                Mint *lda, Mf_complex *x, Mint *incx)
#else
void imsl_ctrsv (uplo, uplo_s, trans, trans_s, diag, diag_s,
                n, a, lda, x, incx)
    Mchar      *uplo;
    unsigned    uplo_s;
    Mchar      *trans;
    unsigned    trans_s;
    Mchar      *diag;
    unsigned    diag_s;
    Mint       *n;
    Mf_complex *a;
    Mint       *lda;
    Mf_complex  x[];
    Mint       *incx;
#endif
{
#define A(I_,J_)	(a+(I_)*(alda)+(J_))
    Mint        alda = *lda;
    Mlong        imsl_ctran, lower, ndiag, ntran, tran, udiag, upper;
    Mint        _l0, _l1, i, ix;


    upper = imsl_l1ame (uplo, uplo_s, "U", sizeof ("U"));
    lower = imsl_l1ame (uplo, uplo_s, "L", sizeof ("L"));
    udiag = imsl_l1ame (diag, diag_s, "U", sizeof ("U"));
    ndiag = imsl_l1ame (diag, diag_s, "N", sizeof ("N"));
    ntran = imsl_l1ame (trans, trans_s, "N", sizeof ("N"));
    tran = imsl_l1ame (trans, trans_s, "T", sizeof ("T"));
    imsl_ctran = imsl_l1ame (trans, trans_s, "C", sizeof ("C"));
    /*
     * Test the input parameters.
     */
    if (*n < 0) {
	imsl_e1psh ("imsl_ctrsv");
	imsl_e1sti (1, *n);

	/*
	 * (5, 1, "N must be greater than or equal to zero while %(i1) is
	 * given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
	imsl_e1pop ("imsl_ctrsv");
	goto L_9000;
    }
    else if ((*lda < *n) || (*lda == 0)) {
	imsl_e1psh ("imsl_ctrsv");
	imsl_e1sti (1, *lda);
	imsl_e1sti (2, *n);

	/*
	 * (5, 2, "LDA must be greater than or equal to N and greater than
	 * zero while LDA = %(i1) and N = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_INVALID_LDA_VALUE_GIVEN);
	imsl_e1pop ("imsl_ctrsv");
	goto L_9000;
    }
    else if (*incx == 0) {
	imsl_e1psh ("imsl_ctrsv");
	imsl_e1sti (1, *incx);

	/* (5, 3, "INCX must not be equal to zero while %(i1) is given."); */
	imsl_ermes (IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
	imsl_e1pop ("imsl_ctrsv");
	goto L_9000;
    }
    else if (((!ntran) && (!tran)) && (!imsl_ctran)) {
	imsl_e1psh ("imsl_ctrsv");
	imsl_e1stl (1, trans);

	/*
	 * (5, 4, "TRANS must be set equal to N or T or C while %(l1) is
	 * given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_TRANS_MUST_EQUAL_N_T_OR_C);
	imsl_e1pop ("imsl_ctrsv");
	goto L_9000;
    }
    else if ((!upper) && (!lower)) {
	imsl_e1psh ("imsl_ctrsv");
	imsl_e1stl (1, uplo);

	/* (5, 5, "UPLO must be set equal to U or L while %(l1) is given."); */
	imsl_ermes (IMSL_TERMINAL, IMSL_INVALID_UPLO_VALUE);
	imsl_e1pop ("imsl_ctrsv");
	goto L_9000;
    }
    else if ((!udiag) && (!ndiag)) {
	imsl_e1psh ("imsl_ctrsv");
	imsl_e1stl (1, diag);

	/* (5, 6, "DIAG must be set equal to U or N while %(l1) is given."); */
	imsl_ermes (IMSL_TERMINAL, IMSL_DIAG_MUST_EQUAL_U_OR_N);
	imsl_e1pop ("imsl_ctrsv");
	goto L_9000;
    }
    /*
     * Quick return if possible.
     */
    if (*n == 0)
	goto L_9000;

    if (upper) {
	if (tran) {
	    if (*incx > 0) {
		ix = 1;
		for (i = 1; i <= *n; i++) {
		    _l0 = i - 1;
		    _l1 = 1;
		    x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotu (&_l0, A (i - 1, 0), &_l1, x, incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], *A (i - 1, i - 1));
		    ix += *incx;
		}
	    }
	    else {
		ix = (-*n + 1) ** incx + 1;
		for (i = 1; i <= *n; i++) {
		    _l0 = i - 1;
		    _l1 = 1;
		    x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotu (&_l0, A (i - 1, 0), &_l1, &x[ix - *incx - 1],
			    incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], *A (i - 1, i - 1));
		    ix += *incx;
		}
	    }
	}
	else if (imsl_ctran) {
	    if (*incx > 0) {
		ix = 1;
		for (i = 1; i <= *n; i++) {
		    _l0 = i - 1;
		    _l1 = 1;
		    x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotc (&_l0, A (i - 1, 0), &_l1, x, incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], imsl_c_conjg (*A (i - 1, i - 1)));
		    ix += *incx;
		}
	    }
	    else {
		ix = (-*n + 1) ** incx + 1;
		for (i = 1; i <= *n; i++) {
		    _l0 = i - 1;
		    _l1 = 1;
		    x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotc (&_l0, A (i - 1, 0), &_l1, &x[ix - *incx - 1],
			    incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], imsl_c_conjg (*A (i - 1, i - 1)));
		    ix += *incx;
		}
	    }
	}
	else {
	    if (*incx > 0) {
		ix = (*n - 1) ** incx + 1;
		for (i = *n; i >= 1; i--) {
		    _l0 = *n - i;
		    if (i < *n)
			x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotu (&_l0, A (i, i - 1), lda, &x[ix + *incx - 1],
				incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], *A (i - 1, i - 1));
		    ix -= *incx;
		}
	    }
	    else {
		ix = 1;
		for (i = *n; i >= 1; i--) {
		    _l0 = *n - i;
		    if (i < *n)
			x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotu (&_l0, A (i, i - 1), lda, x, incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], *A (i - 1, i - 1));
		    ix -= *incx;
		}
	    }
	}
    }
    else {
	if (tran) {
	    if (*incx > 0) {
		ix = (*n - 1) ** incx + 1;
		for (i = *n; i >= 1; i--) {
		    _l0 = *n - i;
		    _l1 = 1;
		    if (i < *n)
			x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotu (&_l0, A (i - 1, i), &_l1, &x[ix + *incx - 1],
				incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], *A (i - 1, i - 1));
		    ix -= *incx;
		}
	    }
	    else {
		ix = 1;
		for (i = *n; i >= 1; i--) {
		    _l0 = *n - i;
		    _l1 = 1;
		    if (i < *n)
			x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotu (&_l0, A (i - 1, i), &_l1, x, incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], *A (i - 1, i - 1));
		    ix -= *incx;
		}
	    }
	}
	else if (imsl_ctran) {
	    if (*incx > 0) {
		ix = (*n - 1) ** incx + 1;
		for (i = *n; i >= 1; i--) {
		    _l0 = *n - i;
		    _l1 = 1;
		    if (i < *n)
			x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotc (&_l0, A (i - 1, i), &_l1, &x[ix + *incx - 1],
				incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], imsl_c_conjg (*A (i - 1, i - 1)));
		    ix -= *incx;
		}
	    }
	    else {
		ix = 1;
		for (i = *n; i >= 1; i--) {
		    _l0 = *n - i;
		    _l1 = 1;
		    if (i < *n)
			x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotc (&_l0, A (i - 1, i), &_l1, x, incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], imsl_c_conjg (*A (i - 1, i - 1)));
		    ix -= *incx;
		}
	    }
	}
	else {
	    if (*incx > 0) {
		ix = 1;
		for (i = 1; i <= *n; i++) {
		    _l0 = i - 1;
		    x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotu (&_l0, A (0, i - 1), lda, x, incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], *A (i - 1, i - 1));
		    ix += *incx;
		}
	    }
	    else {
		ix = (-*n + 1) ** incx + 1;
		for (i = 1; i <= *n; i++) {
		    _l0 = i - 1;
		    x[ix - 1] = imsl_c_sub (x[ix - 1], imsl_cdotu (&_l0, A (0, i - 1), lda, &x[ix - *incx - 1], incx));
		    if (ndiag)
			x[ix - 1] = imsl_c_div (x[ix - 1], *A (i - 1, i - 1));
		    ix += *incx;
		}
	    }
	}
    }

L_9000:
    return;
}				/* end of function */
#undef  A
/*----------------------------------------------------------------------- */

/*  IMSL Name:  CGERU  (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    July 12, 1989

    Purpose:    Perform the rank-one matrix update:
                    A = A + alpha*x*trans(y),
                where trans(y) is the transpose of the vector.

    Usage:      CALL CGERU (M, N, ALPHA, X, INCX, Y, INCY, A, LDA)

    Arguments:
       M      - Number of rows in A.  (Input)
       N      - Number of columns in A.  (Input)
       ALPHA  - Complex scalar.  (Input)
       X      - Complex vector of length (M-1)*IABS(INCX)+1.  (Input)
       INCX   - Displacement between elements of X.  (Input)
       Y      - Complex vector of length (N-1)*IABS(INCY)+1.  (Input)
       INCY   - Displacement between elements of Y.  (Input)
       A      - Complex array of size M by N.  (Input/Output)
                On input, A contains the matrix to be updated.
                On output, A contains the updated matrix.
       LDA    - Leading dimension of A exactly as specified in the
                calling routine.  (Input)

    GAMS:       D1b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_cgeru (Mint *m, Mint *n, Mf_complex *alpha, Mf_complex *x,
                Mint *incx, Mf_complex *y, Mint *incy, Mf_complex *a,
                Mint *lda)
#else
void imsl_cgeru (m, n, alpha, x, incx, y, incy, a, lda)
    Mint       *m, *n;
    Mf_complex *alpha, x[];
    Mint       *incx;
    Mf_complex  y[];
    Mint       *incy;
    Mf_complex  a[];
    Mint       *lda;
#endif
{
    Mint        _l0, imsl_i1x, iy, j;
    Mf_complex  _cx0;
    /*
     * Test the input parameters.
     */
    if (*m < 0) {
	imsl_e1psh ("imsl_cgeru");
	imsl_e1sti (1, *m);

	/*
	 * (5, 1, "M must be greater than or equal to zero while %(i1) is
	 * given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_M_GE_ZERO);
	imsl_e1pop ("imsl_cgeru");
	goto L_9000;
    }
    else if (*n < 0) {
	imsl_e1psh ("imsl_cgeru");
	imsl_e1sti (1, *n);

	/*
	 * (5, 2, "N must be greater than or equal to zero while %(i1) is
	 * given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
	imsl_e1pop ("imsl_cgeru");
	goto L_9000;
    }
    else if ((*lda < *m) || (*lda == 0)) {
	imsl_e1psh ("imsl_cgeru");
	imsl_e1sti (1, *lda);
	imsl_e1sti (2, *m);

	/*
	 * (5, 3, "LDA must be greater than or equal to M and greater than
	 * zero while LDA = %(i1) and M = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_LDA_GE_M);
	imsl_e1pop ("imsl_cgeru");
	goto L_9000;
    }
    else if (*incx == 0) {
	imsl_e1psh ("imsl_cgeru");
	imsl_e1sti (1, *incx);

	/* (5, 4, "INCX must not be equal to zero while %(i1) is given."); */
	imsl_ermes (IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
	imsl_e1pop ("imsl_cgeru");
	goto L_9000;
    }
    else if (*incy == 0) {
	imsl_e1psh ("imsl_cgeru");
	imsl_e1sti (1, *incy);

	/* (5, 5, "INCY must not be equal to zero while %(i1) is given."); */
	imsl_ermes (IMSL_TERMINAL, IMSL_INCY_EQUALS_ZERO);
	imsl_e1pop ("imsl_cgeru");
	goto L_9000;
    }
    /* Quick return if possible */
    if ((*m == 0 || *n == 0) || imsl_c_eq (*alpha, imsl_cf_convert (F_ZERO, F_ZERO)))
	goto L_9000;
    iy = 1;
    if (*incy < 0)
	iy = (-*n + 1) ** incy + 1;

    imsl_i1x = 1;
    for (j = 1; j <= *n; j++) {
	_cx0 = imsl_c_mul (*alpha, y[iy - 1]);
	_l0 = 1;
	imsl_caxpy (m, &_cx0, x, incx, &a[imsl_i1x - 1],
	    &_l0);
	iy += *incy;
	imsl_i1x += *lda;
    }

L_9000:
    return;
}				/* end of function */
