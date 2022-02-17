#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static VA_LIST_HACK PROTO (l_lin_sol_posdef, (Mint n, Mf_complex a[], Mf_complex b[],
	            va_list argptr));
    static void PROTO (l_lftdh, (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *imsl_fac,
	            Mint *ldfac));
    static void PROTO (l_l2cdh, (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *imsl_fac,
	            Mint *ldfac, Mfloat *rcond, Mf_complex *z));
    static void PROTO (l_lslct, (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *b,
	            Mint *ipath, Mf_complex *x));
    static void PROTO (l_chfcg, (Mint *n, Mf_complex *a, Mint *lda));
    static void PROTO (l_ctrsv, (Mchar *uplo, unsigned uplo_s, Mchar *trans, unsigned
	            trans_s, Mchar *diag, unsigned diag_s, Mint *n, Mf_complex *a,
	            Mint *lda, Mf_complex *x, Mint *incx));

    static Mf_complex *lv_x;
#ifdef ANSI
    Mf_complex *imsl_c_lin_sol_posdef (Mint n, Mf_complex *a, Mf_complex *b,...)
#else
    Mf_complex *imsl_c_lin_sol_posdef (n, a, b, va_alist)
    Mint        n;
    Mf_complex *a;
    Mf_complex *b;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, b);

    E1PSH ("imsl_c_lin_sol_posdef", "imsl_z_lin_sol_posdef");
    lv_x = NULL;
    IMSL_CALL (l_lin_sol_posdef (n, a, b, argptr));
    va_end (argptr);

    E1POP ("imsl_c_lin_sol_posdef", "imsl_z_lin_sol_posdef");
    return lv_x;
}


#ifdef ANSI
static VA_LIST_HACK l_lin_sol_posdef (Mint n, Mf_complex a[], Mf_complex b[], va_list argptr)
#else
static VA_LIST_HACK l_lin_sol_posdef (n, a, b, argptr)
    Mint        n;
    Mf_complex *a;
    Mf_complex *b;
    va_list     argptr;
#endif
{
    Mint        code = 1;
    Mint        arg_number = 3;
    Mint        a_col_dim = n;
    Mf_complex **factor_ptr = NULL;
    Mf_complex *factor = NULL;
    Mfloat     *condition = NULL;
    Mint        fac_col_dim = n;
    Mint        return_factor = 0;
    Mint        user_factor = 0;
    Mint        factor_only = 0;
    Mint        solve_only = 0;
    Mint        user_solution = 0;
    Mint        error = 0;
    Mint        factor_transpose = 1;
    Mint        return_condition = 0;
    Mint        ldfac = n;
    Mf_complex *work = NULL;
    Mint        i;
    static Mf_complex c_zero = {0.0, 0.0};


    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_RETURN_USER:
	    lv_x = va_arg (argptr, Mf_complex *);
            arg_number++;
	    if (!lv_x) {
		imsl_e1stl (1, "x");
		imsl_e1stl (2, "IMSL_RETURN_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    user_solution = 1;
	    break;
	case IMSL_A_COL_DIM:
	    a_col_dim = va_arg (argptr, Mint);
            arg_number++; 
	    break;
	case IMSL_FACTOR:
	    return_factor = 1;
	    factor_ptr = va_arg (argptr, Mf_complex **);
            arg_number++; 
	    break;
	case IMSL_FACTOR_USER:
	    return_factor = 1;
	    user_factor = 1;
	    factor = va_arg (argptr, Mf_complex *);
            arg_number++; 
	    if (!factor) {
		imsl_e1stl (1, "factor");
		imsl_e1stl (2, "IMSL_FACTOR_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    break;
	case IMSL_FAC_COL_DIM:
	    fac_col_dim = va_arg (argptr, Mint);
            arg_number++; 
	    break;
	case IMSL_FACTOR_ONLY:
	    factor_only = 1;
	    break;
	case IMSL_SOLVE_ONLY:
	    solve_only = 1;
	    break;
	case IMSL_CONDITION:
	    condition = (Mfloat *) va_arg (argptr, Mfloat *);
            arg_number++; 
	    if (!condition) {
		imsl_e1stl (1, "cond");
		imsl_e1stl (2, "IMSL_CONDITION");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    return_condition = 1;
	    break;
	case 0:
	    break;
	default:
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, arg_number);
	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    ++error;
	}
    }
    if (!solve_only && a == NULL) {
	imsl_e1stl (1, "a");
	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
	++error;
    }
    if (!factor_only && b == NULL) {
	imsl_e1stl (1, "b");
	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
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
    if (solve_only + factor_only > 1) {
	imsl_ermes (IMSL_TERMINAL, IMSL_BAD_SOLVE_FACTOR);
	error = 1;
    }
    if (error)
	return argptr;

    if (solve_only) {
	if (!user_factor) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_SPECIFY_SOLVE_ONLY);
	    ++error;
	}
	if (return_condition) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_CONDITION_ONLY_SPECIFIER);
	    ++error;
	}
    }
    if (!solve_only && !error) {
	if (!user_factor) {
	    factor = (Mf_complex *) imsl_malloc (fac_col_dim * n * sizeof (Mf_complex));
	    if (!factor) {
		imsl_ermes (IMSL_TERMINAL, IMSL_NO_MEM_FOR_FAC);
		++error;
	    }
	}
	else if (fac_col_dim > n) {
	    imsl_c_m1ran (n, fac_col_dim, factor, factor);
	    if (error = imsl_n1rty (1) > 3 ? 1 : 0)
		factor_transpose = 0;
	}
	if (!error) {
	    imsl_c_m1ran (n, a_col_dim, a, a);
	    error = imsl_n1rty (1) > 3 ? 1 : 0;

	    if (!error) {
		Mint        lda = n;
		if (factor_only && !return_condition) {
		    l_lftdh (&n, a, &lda, factor, &ldfac);
		    error = imsl_n1rty (1) > 3 ? 1 : 0;
		}
		else {
		    if (!(work = (Mf_complex *) imsl_malloc (n * sizeof (Mf_complex)))) {
			++error;
		    }
		    else {
			Mfloat      cond;
			l_l2cdh (&n, a, &lda, factor, &ldfac, &cond, work);
			if (!(error = imsl_n1rty (1) > 3 ? 1 : 0) && return_condition) {
			    if (cond > imsl_amach (1)) {
				*condition = 1.0 / cond;
			    }
			    else {
				*condition = imsl_amach (7);
			    }
			}
			if (work != NULL) imsl_free (work);
		    }
		}
		imsl_c_m1ran (a_col_dim, n, a, a);
	    }
	}
    }
    if (!factor_only && !error) {
	if (!user_solution) {
	    lv_x = (Mf_complex *) imsl_malloc (n * sizeof (Mf_complex));
	    if (!lv_x) {
		imsl_ermes (IMSL_TERMINAL, IMSL_NO_MEM_FOR_SYS);
		++error;
	    }
	}
	if (!error && solve_only) {
	    imsl_c_m1ran (n, fac_col_dim, factor, factor);
	    if (error = imsl_n1rty (1) > 3 ? 1 : 0)
		factor_transpose = 0;
	}
	if (!error) {
	    Mint        l = 4;
	    /* solve ctrans(r)*y = b */
	    l_lslct (&n, factor, &ldfac, b, &l, lv_x);
	    error = imsl_n1rty (1) > 3 ? 1 : 0;

	    if (!error) {
		/* solve r*x = y */
		l = 2;
		l_lslct (&n, factor, &ldfac, lv_x, &l, lv_x);
		error = imsl_n1rty (1) > 3 ? 1 : 0;
	    }
	}
    }

    if (!user_factor) {
	if (!return_factor || error) {
	    if (factor != NULL) imsl_free (factor);
	}
	else {
	    imsl_c_m1ran (fac_col_dim, n, factor, factor);
	    for (i = n; i < fac_col_dim; i++) {
		imsl_cset (&n, &c_zero, (factor + i), &fac_col_dim);
	    }
	    *factor_ptr = factor;
	}
	factor = NULL;
    }
    else if (factor_transpose) {
	imsl_c_m1ran (fac_col_dim, n, factor, factor);
    }
    if (error && !user_solution) {
	if (lv_x != NULL) imsl_free (lv_x);
	lv_x = NULL;
    }
    return argptr;
}

#if 0 /* old lftdh */

/*
  -----------------------------------------------------------------------
    IMSL Name:  LFTDH/DLFTDH (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 27, 1990

    Purpose:    Compute the hermite(R)*R factorization of a d_complex
                Hermitian positive definite matrix.

    Usage:      CALL LFTDH (N, A, LDA, FAC, LDFAC)

    Arguments:
       N      - Order of the matrix.  (Input)
       A      - Complex N by N Hermitian positive definite matrix to be
                factored.  (Input)
                Only the upper triangle of A is referenced.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       FAC    - Complex N by N matrix containing the upper triangular
                matrix R of the factorization of A in the upper triangle.
                (Output)
                Only the upper triangle of FAC will be used.  If A is not
                needed, A and FAC can share the same storage locations.
       LDFAC  - Leading dimension of FAC exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remark:
       Informational errors
       Type Code
         3   4  The input matrix is not Hermitian.  It has a diagonal
                entry with a small imaginary part.
         4   2  The input matrix is not positive definite.
         4   4  The input matrix is not Hermitian.  It has a diagonal
                entry with an imaginary part.

    Keyword:    Cholesky factorization

    GAMS:       D2d1b

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------

 */
#define	NB	129
#ifdef ANSI
static void l_lftdh (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *imsl_fac,
                Mint *ldfac)
#else
static void l_lftdh (n, a, lda, imsl_fac, ldfac)
    Mint       *n;
    Mf_complex *a;
    Mint       *lda;
    Mf_complex *imsl_fac;
    Mint       *ldfac;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define FAC(I_,J_)	(imsl_fac+(I_)*(*ldfac)+(J_))
    Mint        _l0, _l1, i, info, j, k, l, ll, nlenb, nrcfac;
    Mfloat      big, eps, r0, r1, r2, r3, r4, r5, r6, r7, small;
    Mf_complex  rtemp, t[8][NB];


    /*
     * this code is for computer types: fosivv and rtxlxs
     */
    imsl_e1psh ("l_lftdh");

    /* PRESERVE A COPY OF THE INPUT MATRIX */
    for (i = 1; i <= *n; i++) {
	_l0 = 1;
	_l1 = 1;
	imsl_ccopy (&i, A (i - 1, 0), &_l0, FAC (i - 1, 0), &_l1);
    }
    /* Check that A is Hermitian */
    eps = F_TEN * imsl_amach (4);
    for (i = 1; i <= *n; i++) {
	if (fabs (imsl_c_aimag (*FAC (i - 1, i - 1))) != F_ZERO) {
	    if (fabs (imsl_c_aimag (*FAC (i - 1, i - 1))) > eps * fabs (imsl_fc_convert (*FAC (i - 1, i - 1)))) {
		imsl_e1sti (1, i-1);
		imsl_e1stc (1, *FAC (i - 1, i - 1));
		/*
		 * (4, 4, "The matrix element A(%(i1),%(i1)) = %(c1).  The
		 * diagonal of a Hermitian must be real.");
		 */
		imsl_ermes (IMSL_FATAL, IMSL_HERMITIAN_DIAG_REAL);
		goto L_9000;
	    }
	    else {
		imsl_e1sti (1, i-1);
		imsl_e1stc (1, *FAC (i - 1, i - 1));
		/*
		 * (3, 4, "The matrix element A(%(i1),%(i1)) = %(c1).  The
		 * diagonal of a Hermitian matrix must be real.  The
		 * imaginary part will be used as zero in the algorithm.");
		 */
		imsl_ermes (IMSL_WARNING, IMSL_HERMITIAN_DIAG_REAL_2);
		*FAC (i - 1, i - 1) = imsl_cf_convert (imsl_fc_convert (*FAC (i - 1, i - 1)), F_ZERO);
	    }
	}
    }

    /*
     * Cholesky decomposition using method LLT**
     * 
     * A brief description of the algorithm follows: For a symmetric positive
     * definite matrix at the k-th step :
     * 
     * k   |  11  12  13 | A  = |  21  22  23 | , |  31  32  33 |
     * 
     * where trans(A11 A21 A31) is n x (k*8) and factored. assume trans(A22 A32)
     * is the active block of 8 columns. The factorization is accomplished
     * by:
     * 
     * Step 1. Factor trans(A22 A32) and store -trans(A22 A32) into a local work
     * array
     * 
     * Step 2. Update A32 with matrix multiplication between trans(A22 A32) and
     * the local work array.
     * 
     * Step 3. Repeat
     * 
     * Fill in the lower triangle
     */
    l_chfcg (n, imsl_fac, ldfac);

    small = imsl_amach (1);
    big = imsl_amach (2);
    if (small * big < F_ONE)
	small = F_ONE / big;
    info = 0;

    nrcfac = mod (*n, 8);
    for (j = 1; j <= (*n - nrcfac); j += 8) {
	nlenb = imsl_i_min (NB, *n - j);
	/* prepare j-th column */
	if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= F_ZERO) {
	    info = j;
	    goto L_460;
	}
	*FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
	    F_ZERO);
	r0 = imsl_fc_convert (imsl_c_div (C_ONE, *FAC (j - 1, j - 1)));
	/*
	 * Form the j-th mutiplier and load t(1:nlenb,1) with
	 * -imsl_fac(j+1:nlenb,j)
	 */
	for (i = 1; i <= nlenb; i++) {
	    *FAC (j - 1, j + i - 1) = imsl_c_mul (imsl_cf_convert (r0, F_ZERO), *FAC (j - 1, j + i - 1));
	    t[0][i - 1] = imsl_c_neg (imsl_c_conjg (*FAC (j - 1, j + i - 1)));
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j - 1, i - 1) = imsl_c_mul (imsl_cf_convert (r0, F_ZERO), *FAC (j - 1, i - 1));
	}
	/* update columns j+1 thru j+7 */
	for (k = 1; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) = imsl_c_add (*FAC (j + k - 1, j + i - 1),
		    imsl_c_mul (*FAC (j - 1, j + i - 1), t[0][k - 1]));
	    }
	}
	/* prepare (j+1)-th column */
	if (imsl_fc_convert (*FAC (j, j)) <= F_ZERO) {
	    info = j + 1;
	    goto L_460;
	}
	*FAC (j, j) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j, j))), F_ZERO);
	r1 = imsl_fc_convert (imsl_c_div (C_ONE, *FAC (j, j)));
	/*
	 * Form the (j+1)-th mutiplier and load t(2:nlenb,2) with
	 * -imsl_fac(j+2:nlenb,j+1)
	 */
	for (i = 2; i <= nlenb; i++) {
	    *FAC (j, j + i - 1) = imsl_c_mul (imsl_cf_convert (r1, F_ZERO), *FAC (j, j + i - 1));
	    t[1][i - 1] = imsl_c_neg (imsl_c_conjg (*FAC (j, j + i - 1)));
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j, i - 1) = imsl_c_mul (imsl_cf_convert (r1, F_ZERO), *FAC (j, i - 1));
	}
	/* update columns j+2 thru j+7 */
	for (k = 2; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) = imsl_c_add (*FAC (j + k - 1, j + i - 1),
		    imsl_c_mul (*FAC (j, j + i - 1), t[1][k - 1]));
	    }
	}
	/* prepare (j+2)-th column */
	if (imsl_fc_convert (*FAC (j + 1, j + 1)) <= F_ZERO) {
	    info = j + 2;
	    goto L_460;
	}
	*FAC (j + 1, j + 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j + 1, j + 1))),
	    F_ZERO);
	r2 = imsl_fc_convert (imsl_c_div (C_ONE, *FAC (j + 1, j + 1)));
	/*
	 * Form the (j+2)-th mutiplier and load t(3:nlenb,3) with
	 * -imsl_fac(j+3:nlenb,j+2)
	 */
	for (i = 3; i <= nlenb; i++) {
	    *FAC (j + 1, j + i - 1) = imsl_c_mul (imsl_cf_convert (r2, F_ZERO), *FAC (j + 1, j + i - 1));
	    t[2][i - 1] = imsl_c_neg (imsl_c_conjg (*FAC (j + 1, j + i - 1)));
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 1, i - 1) = imsl_c_mul (imsl_cf_convert (r2, F_ZERO), *FAC (j + 1, i - 1));
	}
	/* update columns j+3 thru j+7 */
	for (k = 3; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) = imsl_c_add (*FAC (j + k - 1, j + i - 1),
		    imsl_c_mul (*FAC (j + 1, j + i - 1), t[2][k - 1]));
	    }
	}
	/* prepare (j+3)-th column */
	if (imsl_fc_convert (*FAC (j + 2, j + 2)) <= F_ZERO) {
	    info = j + 3;
	    goto L_460;
	}
	*FAC (j + 2, j + 2) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j + 2, j + 2))),
	    F_ZERO);
	r3 = imsl_fc_convert (imsl_c_div (C_ONE, *FAC (j + 2, j + 2)));
	/*
	 * Form the (j+3)-th mutiplier and load t(4:nlenb,4) with
	 * -imsl_fac(j+3:nlenb,j+3)
	 */
	for (i = 4; i <= nlenb; i++) {
	    *FAC (j + 2, j + i - 1) = imsl_c_mul (imsl_cf_convert (r3, F_ZERO), *FAC (j + 2, j + i - 1));
	    t[3][i - 1] = imsl_c_neg (imsl_c_conjg (*FAC (j + 2, j + i - 1)));
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 2, i - 1) = imsl_c_mul (imsl_cf_convert (r3, F_ZERO), *FAC (j + 2, i - 1));
	}
	/* update columns j+4 thru j+7 */
	for (k = 4; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) = imsl_c_add (*FAC (j + k - 1, j + i - 1),
		    imsl_c_mul (*FAC (j + 2, j + i - 1), t[3][k - 1]));
	    }
	}
	/* prepare (j+4)-th column */
	if (imsl_fc_convert (*FAC (j + 3, j + 3)) <= F_ZERO) {
	    info = j + 4;
	    goto L_460;
	}
	*FAC (j + 3, j + 3) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j + 3, j + 3))),
	    F_ZERO);
	r4 = imsl_fc_convert (imsl_c_div (C_ONE, *FAC (j + 3, j + 3)));
	/*
	 * Form the (j+4)-th mutiplier and load t(5:nlenb,5) with
	 * -imsl_fac(j+4:nlenb,j+4)
	 */
	for (i = 5; i <= nlenb; i++) {
	    *FAC (j + 3, j + i - 1) = imsl_c_mul (imsl_cf_convert (r4, F_ZERO), *FAC (j + 3, j + i - 1));
	    t[4][i - 1] = imsl_c_neg (imsl_c_conjg (*FAC (j + 3, j + i - 1)));
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 3, i - 1) = imsl_c_mul (imsl_cf_convert (r4, F_ZERO), *FAC (j + 3, i - 1));
	}
	/* update columns j+5 thru j+7 */
	for (k = 5; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) = imsl_c_add (*FAC (j + k - 1, j + i - 1),
		    imsl_c_mul (*FAC (j + 3, j + i - 1), t[4][k - 1]));
	    }
	}
	/* prepare (j+5)-th column */
	if (imsl_fc_convert (*FAC (j + 4, j + 4)) <= F_ZERO) {
	    info = j + 5;
	    goto L_460;
	}
	*FAC (j + 4, j + 4) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j + 4, j + 4))),
	    F_ZERO);
	r5 = imsl_fc_convert (imsl_c_div (C_ONE, *FAC (j + 4, j + 4)));
	/*
	 * Form the (j+5)-th mutiplier and load t(5:nlenb,5) with
	 * -imsl_fac(j+3:nlenb,j+3)
	 */
	for (i = 6; i <= nlenb; i++) {
	    *FAC (j + 4, j + i - 1) = imsl_c_mul (imsl_cf_convert (r5, F_ZERO), *FAC (j + 4, j + i - 1));
	    t[5][i - 1] = imsl_c_neg (imsl_c_conjg (*FAC (j + 4, j + i - 1)));
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 4, i - 1) = imsl_c_mul (imsl_cf_convert (r5, F_ZERO), *FAC (j + 4, i - 1));
	}
	/* update columns j+6 thru j+7 */
	for (k = 6; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) = imsl_c_add (*FAC (j + k - 1, j + i - 1),
		    imsl_c_mul (*FAC (j + 4, j + i - 1), t[5][k - 1]));
	    }
	}
	/* prepare (j+6)-th column */
	if (imsl_fc_convert (*FAC (j + 5, j + 5)) <= F_ZERO) {
	    info = j + 6;
	    goto L_460;
	}
	*FAC (j + 5, j + 5) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j + 5, j + 5))),
	    F_ZERO);
	r6 = imsl_fc_convert (imsl_c_div (C_ONE, *FAC (j + 5, j + 5)));
	/*
	 * Form the (j+6)-th mutiplier and load t(7:nlenb,7) with
	 * -imsl_fac(j+7:nlenb,j+6)
	 */
	for (i = 7; i <= nlenb; i++) {
	    *FAC (j + 5, j + i - 1) = imsl_c_mul (imsl_cf_convert (r6, F_ZERO), *FAC (j + 5, j + i - 1));
	    t[6][i - 1] = imsl_c_neg (imsl_c_conjg (*FAC (j + 5, j + i - 1)));
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 5, i - 1) = imsl_c_mul (imsl_cf_convert (r6, F_ZERO), *FAC (j + 5, i - 1));
	}
	/* update column j+7 */
	for (i = 7; i <= (*n - j); i++) {
	    *FAC (j + 6, j + i - 1) = imsl_c_add (*FAC (j + 6, j + i - 1), imsl_c_mul (*FAC (j + 5, j + i - 1),
		    t[6][6]));
	}
	/* prepare (j+7)-th column */
	if (imsl_fc_convert (*FAC (j + 6, j + 6)) <= F_ZERO) {
	    info = j + 7;
	    goto L_460;
	}
	*FAC (j + 6, j + 6) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j + 6, j + 6))),
	    F_ZERO);
	r7 = imsl_fc_convert (imsl_c_div (C_ONE, *FAC (j + 6, j + 6)));
	/*
	 * Form the (j+7)-th mutiplier and load t(8:nlenb,8) with
	 * -imsl_fac(j+8:nlenb,j+7)
	 */
	for (i = 8; i <= nlenb; i++) {
	    *FAC (j + 6, j + i - 1) = imsl_c_mul (imsl_cf_convert (r7, F_ZERO), *FAC (j + 6, j + i - 1));
	    t[7][i - 1] = imsl_c_neg (imsl_c_conjg (*FAC (j + 6, j + i - 1)));
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 6, i - 1) = imsl_c_mul (imsl_cf_convert (r7, F_ZERO), *FAC (j + 6, i - 1));
	}
	/*
	 * Perform update on the lower triangle on columns j+7 thru j + nlenb
	 * rows j+7 thru n
	 */
	for (k = 8; k <= nlenb; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (j + k - 1, j + i - 1),
						imsl_c_mul (*FAC (j - 1, j + i - 1), t[0][k - 1])), imsl_c_mul (*FAC (j, j + i - 1),
						t[1][k - 1])), imsl_c_mul (*FAC (j + 1, j + i - 1), t[2][k - 1])),
				    imsl_c_mul (*FAC (j + 2, j + i - 1), t[3][k - 1])), imsl_c_mul (*FAC (j + 3, j + i - 1),
				    t[4][k - 1])), imsl_c_mul (*FAC (j + 4, j + i - 1), t[5][k - 1])),
			imsl_c_mul (*FAC (j + 5, j + i - 1), t[6][k - 1])), imsl_c_mul (*FAC (j + 6, j + i - 1),
			t[7][k - 1]));
	    }
	}

	for (ll = j + NB + 1; ll <= *n; ll += NB) {
	    l = ll - 1;
	    nlenb = imsl_i_min (*n - l, NB);
	    /*
	     * form mutipliers j,j+1,j+2,j+3 rows ll thru ll+nlenb
	     */
	    for (k = 0; k <= 7; k++) {
		for (i = 1; i <= nlenb; i++) {
		    t[k][i - 1] = imsl_c_neg (imsl_c_conjg (*FAC (j + k - 1, l + i - 1)));
		}
	    }
	    /*
	     * Update lower triangle from rows ll thru ll+n columns ll thru
	     * ll + nlenb
	     */
	    for (k = 1; k <= nlenb; k++) {
		for (i = k; i <= (*n - l); i++) {
		    *FAC (l + k - 1, l + i - 1) = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (l + k - 1, l + i - 1),
						    imsl_c_mul (*FAC (j - 1, l + i - 1), t[0][k - 1])), imsl_c_mul (*FAC (j, l + i - 1),
						    t[1][k - 1])), imsl_c_mul (*FAC (j + 1, l + i - 1), t[2][k - 1])),
					imsl_c_mul (*FAC (j + 2, l + i - 1), t[3][k - 1])), imsl_c_mul (*FAC (j + 3, l + i - 1),
					t[4][k - 1])), imsl_c_mul (*FAC (j + 4, l + i - 1), t[5][k - 1])),
			    imsl_c_mul (*FAC (j + 5, l + i - 1), t[6][k - 1])), imsl_c_mul (*FAC (j + 6, l + i - 1),
			    t[7][k - 1]));
		}
	    }
	}

    }
    /*
     * Take care of remaining nrcfac columns of imsl_fac
     */
    j = *n - nrcfac + 1;
    if (nrcfac < 7)
	goto L_400;

    if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= F_ZERO) {
	info = j;
	goto L_460;
    }
    *FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
	F_ZERO);
    rtemp = imsl_c_div (C_ONE, *FAC (j - 1, j - 1));
    *FAC (j - 1, j) = imsl_c_mul (rtemp, *FAC (j - 1, j));
    *FAC (j - 1, j + 1) = imsl_c_mul (rtemp, *FAC (j - 1, j + 1));
    *FAC (j - 1, j + 2) = imsl_c_mul (rtemp, *FAC (j - 1, j + 2));
    *FAC (j - 1, j + 3) = imsl_c_mul (rtemp, *FAC (j - 1, j + 3));
    *FAC (j - 1, j + 4) = imsl_c_mul (rtemp, *FAC (j - 1, j + 4));
    *FAC (j - 1, j + 5) = imsl_c_mul (rtemp, *FAC (j - 1, j + 5));
    *FAC (j, j) = imsl_c_sub (*FAC (j, j), imsl_c_mul (*FAC (j - 1, j), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 1) = imsl_c_sub (*FAC (j, j + 1), imsl_c_mul (*FAC (j - 1, j + 1), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 2) = imsl_c_sub (*FAC (j, j + 2), imsl_c_mul (*FAC (j - 1, j + 2), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 3) = imsl_c_sub (*FAC (j, j + 3), imsl_c_mul (*FAC (j - 1, j + 3), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 4) = imsl_c_sub (*FAC (j, j + 4), imsl_c_mul (*FAC (j - 1, j + 4), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 5) = imsl_c_sub (*FAC (j, j + 5), imsl_c_mul (*FAC (j - 1, j + 5), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j + 1, j + 1) = imsl_c_sub (*FAC (j + 1, j + 1), imsl_c_mul (*FAC (j - 1, j + 1),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 1, j + 2) = imsl_c_sub (*FAC (j + 1, j + 2), imsl_c_mul (*FAC (j - 1, j + 2),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 1, j + 3) = imsl_c_sub (*FAC (j + 1, j + 3), imsl_c_mul (*FAC (j - 1, j + 3),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 1, j + 4) = imsl_c_sub (*FAC (j + 1, j + 4), imsl_c_mul (*FAC (j - 1, j + 4),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 1, j + 5) = imsl_c_sub (*FAC (j + 1, j + 5), imsl_c_mul (*FAC (j - 1, j + 5),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 2, j + 2) = imsl_c_sub (*FAC (j + 2, j + 2), imsl_c_mul (*FAC (j - 1, j + 2),
	    imsl_c_conjg (*FAC (j - 1, j + 2))));
    *FAC (j + 2, j + 3) = imsl_c_sub (*FAC (j + 2, j + 3), imsl_c_mul (*FAC (j - 1, j + 3),
	    imsl_c_conjg (*FAC (j - 1, j + 2))));
    *FAC (j + 2, j + 4) = imsl_c_sub (*FAC (j + 2, j + 4), imsl_c_mul (*FAC (j - 1, j + 4),
	    imsl_c_conjg (*FAC (j - 1, j + 2))));
    *FAC (j + 2, j + 5) = imsl_c_sub (*FAC (j + 2, j + 5), imsl_c_mul (*FAC (j - 1, j + 5),
	    imsl_c_conjg (*FAC (j - 1, j + 2))));
    *FAC (j + 3, j + 3) = imsl_c_sub (*FAC (j + 3, j + 3), imsl_c_mul (*FAC (j - 1, j + 3),
	    imsl_c_conjg (*FAC (j - 1, j + 3))));
    *FAC (j + 3, j + 4) = imsl_c_sub (*FAC (j + 3, j + 4), imsl_c_mul (*FAC (j - 1, j + 4),
	    imsl_c_conjg (*FAC (j - 1, j + 3))));
    *FAC (j + 3, j + 5) = imsl_c_sub (*FAC (j + 3, j + 5), imsl_c_mul (*FAC (j - 1, j + 5),
	    imsl_c_conjg (*FAC (j - 1, j + 3))));
    *FAC (j + 4, j + 4) = imsl_c_sub (*FAC (j + 4, j + 4), imsl_c_mul (*FAC (j - 1, j + 4),
	    imsl_c_conjg (*FAC (j - 1, j + 4))));
    *FAC (j + 4, j + 5) = imsl_c_sub (*FAC (j + 4, j + 5), imsl_c_mul (*FAC (j - 1, j + 5),
	    imsl_c_conjg (*FAC (j - 1, j + 4))));
    *FAC (j + 5, j + 5) = imsl_c_sub (*FAC (j + 5, j + 5), imsl_c_mul (*FAC (j - 1, j + 5),
	    imsl_c_conjg (*FAC (j - 1, j + 5))));

    j += 1;
    nrcfac -= 1;
L_400:
    if (nrcfac < 6)
	goto L_410;

    if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= F_ZERO) {
	info = j;
	goto L_460;
    }
    *FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
	F_ZERO);
    rtemp = imsl_c_div (C_ONE, *FAC (j - 1, j - 1));
    *FAC (j - 1, j) = imsl_c_mul (rtemp, *FAC (j - 1, j));
    *FAC (j - 1, j + 1) = imsl_c_mul (rtemp, *FAC (j - 1, j + 1));
    *FAC (j - 1, j + 2) = imsl_c_mul (rtemp, *FAC (j - 1, j + 2));
    *FAC (j - 1, j + 3) = imsl_c_mul (rtemp, *FAC (j - 1, j + 3));
    *FAC (j - 1, j + 4) = imsl_c_mul (rtemp, *FAC (j - 1, j + 4));
    *FAC (j, j) = imsl_c_sub (*FAC (j, j), imsl_c_mul (*FAC (j - 1, j), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 1) = imsl_c_sub (*FAC (j, j + 1), imsl_c_mul (*FAC (j - 1, j + 1), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 2) = imsl_c_sub (*FAC (j, j + 2), imsl_c_mul (*FAC (j - 1, j + 2), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 3) = imsl_c_sub (*FAC (j, j + 3), imsl_c_mul (*FAC (j - 1, j + 3), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 4) = imsl_c_sub (*FAC (j, j + 4), imsl_c_mul (*FAC (j - 1, j + 4), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j + 1, j + 1) = imsl_c_sub (*FAC (j + 1, j + 1), imsl_c_mul (*FAC (j - 1, j + 1),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 1, j + 2) = imsl_c_sub (*FAC (j + 1, j + 2), imsl_c_mul (*FAC (j - 1, j + 2),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 1, j + 3) = imsl_c_sub (*FAC (j + 1, j + 3), imsl_c_mul (*FAC (j - 1, j + 3),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 1, j + 4) = imsl_c_sub (*FAC (j + 1, j + 4), imsl_c_mul (*FAC (j - 1, j + 4),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 2, j + 2) = imsl_c_sub (*FAC (j + 2, j + 2), imsl_c_mul (*FAC (j - 1, j + 2),
	    imsl_c_conjg (*FAC (j - 1, j + 2))));
    *FAC (j + 2, j + 3) = imsl_c_sub (*FAC (j + 2, j + 3), imsl_c_mul (*FAC (j - 1, j + 3),
	    imsl_c_conjg (*FAC (j - 1, j + 2))));
    *FAC (j + 2, j + 4) = imsl_c_sub (*FAC (j + 2, j + 4), imsl_c_mul (*FAC (j - 1, j + 4),
	    imsl_c_conjg (*FAC (j - 1, j + 2))));
    *FAC (j + 3, j + 3) = imsl_c_sub (*FAC (j + 3, j + 3), imsl_c_mul (*FAC (j - 1, j + 3),
	    imsl_c_conjg (*FAC (j - 1, j + 3))));
    *FAC (j + 3, j + 4) = imsl_c_sub (*FAC (j + 3, j + 4), imsl_c_mul (*FAC (j - 1, j + 4),
	    imsl_c_conjg (*FAC (j - 1, j + 3))));
    *FAC (j + 4, j + 4) = imsl_c_sub (*FAC (j + 4, j + 4), imsl_c_mul (*FAC (j - 1, j + 4),
	    imsl_c_conjg (*FAC (j - 1, j + 4))));

    j += 1;
    nrcfac -= 1;
L_410:
    if (nrcfac < 5)
	goto L_420;

    if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= F_ZERO) {
	info = j;
	goto L_460;
    }
    *FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
	F_ZERO);
    rtemp = imsl_c_div (C_ONE, *FAC (j - 1, j - 1));
    *FAC (j - 1, j) = imsl_c_mul (rtemp, *FAC (j - 1, j));
    *FAC (j - 1, j + 1) = imsl_c_mul (rtemp, *FAC (j - 1, j + 1));
    *FAC (j - 1, j + 2) = imsl_c_mul (rtemp, *FAC (j - 1, j + 2));
    *FAC (j - 1, j + 3) = imsl_c_mul (rtemp, *FAC (j - 1, j + 3));
    *FAC (j, j) = imsl_c_sub (*FAC (j, j), imsl_c_mul (*FAC (j - 1, j), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 1) = imsl_c_sub (*FAC (j, j + 1), imsl_c_mul (*FAC (j - 1, j + 1), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 2) = imsl_c_sub (*FAC (j, j + 2), imsl_c_mul (*FAC (j - 1, j + 2), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 3) = imsl_c_sub (*FAC (j, j + 3), imsl_c_mul (*FAC (j - 1, j + 3), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j + 1, j + 1) = imsl_c_sub (*FAC (j + 1, j + 1), imsl_c_mul (*FAC (j - 1, j + 1),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 1, j + 2) = imsl_c_sub (*FAC (j + 1, j + 2), imsl_c_mul (*FAC (j - 1, j + 2),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 1, j + 3) = imsl_c_sub (*FAC (j + 1, j + 3), imsl_c_mul (*FAC (j - 1, j + 3),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 2, j + 2) = imsl_c_sub (*FAC (j + 2, j + 2), imsl_c_mul (*FAC (j - 1, j + 2),
	    imsl_c_conjg (*FAC (j - 1, j + 2))));
    *FAC (j + 2, j + 3) = imsl_c_sub (*FAC (j + 2, j + 3), imsl_c_mul (*FAC (j - 1, j + 3),
	    imsl_c_conjg (*FAC (j - 1, j + 2))));
    *FAC (j + 3, j + 3) = imsl_c_sub (*FAC (j + 3, j + 3), imsl_c_mul (*FAC (j - 1, j + 3),
	    imsl_c_conjg (*FAC (j - 1, j + 3))));

    j += 1;
    nrcfac -= 1;
L_420:
    if (nrcfac < 4)
	goto L_430;

    if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= F_ZERO) {
	info = j;
	goto L_460;
    }
    *FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
	F_ZERO);
    rtemp = imsl_c_div (C_ONE, *FAC (j - 1, j - 1));
    *FAC (j - 1, j) = imsl_c_mul (rtemp, *FAC (j - 1, j));
    *FAC (j - 1, j + 1) = imsl_c_mul (rtemp, *FAC (j - 1, j + 1));
    *FAC (j - 1, j + 2) = imsl_c_mul (rtemp, *FAC (j - 1, j + 2));
    *FAC (j, j) = imsl_c_sub (*FAC (j, j), imsl_c_mul (*FAC (j - 1, j), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 1) = imsl_c_sub (*FAC (j, j + 1), imsl_c_mul (*FAC (j - 1, j + 1), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 2) = imsl_c_sub (*FAC (j, j + 2), imsl_c_mul (*FAC (j - 1, j + 2), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j + 1, j + 1) = imsl_c_sub (*FAC (j + 1, j + 1), imsl_c_mul (*FAC (j - 1, j + 1),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 1, j + 2) = imsl_c_sub (*FAC (j + 1, j + 2), imsl_c_mul (*FAC (j - 1, j + 2),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));
    *FAC (j + 2, j + 2) = imsl_c_sub (*FAC (j + 2, j + 2), imsl_c_mul (*FAC (j - 1, j + 2),
	    imsl_c_conjg (*FAC (j - 1, j + 2))));

    j += 1;
    nrcfac -= 1;
L_430:
    if (nrcfac < 3)
	goto L_440;

    if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= F_ZERO) {
	info = j;
	goto L_460;
    }
    *FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
	F_ZERO);
    rtemp = imsl_c_div (C_ONE, *FAC (j - 1, j - 1));
    *FAC (j - 1, j) = imsl_c_mul (rtemp, *FAC (j - 1, j));
    *FAC (j - 1, j + 1) = imsl_c_mul (rtemp, *FAC (j - 1, j + 1));
    *FAC (j, j) = imsl_c_sub (*FAC (j, j), imsl_c_mul (*FAC (j - 1, j), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 1) = imsl_c_sub (*FAC (j, j + 1), imsl_c_mul (*FAC (j - 1, j + 1), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j + 1, j + 1) = imsl_c_sub (*FAC (j + 1, j + 1), imsl_c_mul (*FAC (j - 1, j + 1),
	    imsl_c_conjg (*FAC (j - 1, j + 1))));

    j += 1;
    nrcfac -= 1;
L_440:
    if (nrcfac < 2)
	goto L_450;

    if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= F_ZERO) {
	info = j;
	goto L_460;
    }
    *FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
	F_ZERO);
    rtemp = imsl_c_div (C_ONE, *FAC (j - 1, j - 1));
    *FAC (j - 1, j) = imsl_c_mul (rtemp, *FAC (j - 1, j));
    *FAC (j, j) = imsl_c_sub (*FAC (j, j), imsl_c_mul (*FAC (j - 1, j), imsl_c_conjg (*FAC (j - 1, j))));

    j += 1;
    nrcfac -= 1;
L_450:
    if (nrcfac < 1)
	goto L_460;

    if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= F_ZERO) {
	info = j;
	goto L_460;
    }
    *FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
	F_ZERO);
L_460:
    if (info != 0) {
	imsl_e1sti (1, info);

	/*
	 * (4, 2, "The leading %(i1) by %(i1) submatrix of the input matrix
	 * is not positive definite.");
	 */
	imsl_ermes (IMSL_FATAL, IMSL_NONPOSITIVE_MATRIX);
    }
    /* Fill in the upper triangle */
    for (i = 1; i <= (*n - 1); i++) {
	for (j = 1; j <= i; j++) {
	    *FAC (i, j - 1) = imsl_c_conjg (*FAC (j - 1, i));
	}
	/*
	 * call imsl_ccopy (n-i, imsl_fac(i+1,i), 1, imsl_fac(i,i+1), ldfac)
	 */
    }

L_9000:
    imsl_e1pop ("l_lftdh");

    return;
}				/* end of function */
#undef  NB
#undef  FAC
#undef  A

#endif /* old lftdh */

/* new lftdh */

/*Translated by FOR_C++, v0.1, on 11/29/91 at 12:53:00 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 11/29/91 at 12:52:53
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  LFTDH/DLFTDH (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 26, 1991

    Purpose:    Compute the hermite(R)*R factorization of a d_complex
                Hermitian positive definite matrix.

    Usage:      CALL LFTDH (N, A, LDA, FAC, LDFAC)

    Arguments:
       N      - Order of the matrix.  (Input)
       A      - Complex N by N Hermitian positive definite matrix to be
                factored.  (Input)
                Only the upper triangle of A is referenced.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       FAC    - Complex N by N matrix containing the upper triangular
                matrix R of the factorization of A in the upper
                triangle.  (Output)
                Only the upper triangle of FAC will be used.  If A is
                not needed, A and FAC can share the same storage
                locations.
       LDFAC  - Leading dimension of FAC exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remark:
       Informational errors
       Type Code
         3   4  The input matrix is not Hermitian.  It has a diagonal
                entry with a small imaginary part.
         4   2  The input matrix is not positive definite.
         4   4  The input matrix is not Hermitian.  It has a diagonal
                entry with an imaginary part.

    Keyword:    Cholesky factorization

    GAMS:       D2d1b

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
static void l_lftdh (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *imsl_fac,
                Mint *ldfac)
#else
static void l_lftdh (n, a, lda, imsl_fac, ldfac)
    Mint        *n;
    Mf_complex  *a;
    Mint        *lda;
    Mf_complex  *imsl_fac;
    Mint        *ldfac;
#endif
{
#define A(I_,J_)        (a+(I_)*(*lda)+(J_))
#define FAC(I_,J_)      (imsl_fac+(I_)*(*ldfac)+(J_))
    Mint         _l0, _l1, i, info, j, k, kmod2, nrcfac;
    Mfloat       big, eps, one, r0, r1, r2, r3, small, zero;
    Mf_complex   rtemp, s1, s2, s3, s4, t1, t2, t3, t4;
    Mf_complex   ctemp1, ctemp2, ctemp3, ctemp4;


    /*
     * this code is for computer types: fosivv and fujxv
     */
    imsl_e1psh ("LFTDH");
    zero = 0.0e0;
    one = 1.0e0;
    small = imsl_amach (1);
    big = imsl_amach (2);
    if (small * big < one)
        small = one / big;

    if (*n <= 0) {
        imsl_e1sti (1, *n);
/*
        (5, 1, "The order of the matrix must be positive
	 while N = %(i1) is given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
    }

    if (*n > *lda) {
        imsl_e1sti (1, *n);
        imsl_e1sti (2, *lda);
/*
        (5, 2, "The order of the matrix must be less than or equal to 
	its leading dimension while N = %(i1) and LDA = %(i2) are given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
    }
 
    if (*n > *ldfac) {
        imsl_e1sti (1, *n);
        imsl_e1sti (2, *ldfac);
/*
        (5, 3, "The order of the matrix must be less than or 
	equal to its leading dimension while N = %(i1) and LDFAC = %(i2) 
	are given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_LDFAC_LESS_ORDER);
    }
    if (imsl_n1rcd (0) != 0)
        goto L_9000;
    /* PRESERVE A COPY OF THE INPUT MATRIX */
    _l0 = 1;
    _l1 = 1;
    for (i = 1; i <= *n; i++) {
        imsl_ccopy (&i, A (i - 1, 0), &_l0, FAC (i - 1, 0), &_l1);
    }
    /*
     * Check that A is Hermitian This is ad hoc code that has little value.
     * A better idea is to just use the real part and ignore the imaginary
     * part.
     */
    eps = 10.0 * imsl_amach (4);
    for (i = 1; i <= *n; i++) {
        if (fabs (imsl_c_aimag (*FAC (i - 1, i - 1))) != zero) {
            if (fabs (imsl_c_aimag (*FAC (i - 1, i - 1))) > eps * fabs (imsl_fc_convert (*FAC (i - 1, i - 1)))) {
                imsl_e1sti (1, i-1);
                imsl_e1stc (1, *FAC (i - 1, i - 1));
 /*
                (4, 4, "The matrix element A(%(i1),%(i1)) = %(c1).  The 
		diagonal of a Hermitian must be real.");
*/
		imsl_ermes (IMSL_FATAL, IMSL_HERMITIAN_DIAG_REAL);
                goto L_9000;
            }
            else {
                imsl_e1sti (1, i-1);
                imsl_e1stc (1, *FAC (i - 1, i - 1));
 /*
                (3, 4, "The matrix element A(%(i1),%(i1)) = %(c1).  The 
		diagonal of a Hermitian matrix must be real.  The imaginary 
		part will be used as zero in the algorithm.");
*/
		imsl_ermes (IMSL_WARNING, IMSL_HERMITIAN_DIAG_REAL_2);
                *FAC (i - 1, i - 1) = imsl_cf_convert (imsl_fc_convert (*FAC (i - 1, i - 1)),
                    zero);
            }
        }
    }

    /*
     * Cholesky decomposition using method LLT**
     * 
     * A brief description of the algorithm follows: For a symmetric positive
     * definite matrix at the k-th step :
     * 
     * k   |  11  12  13 | A  = |  21  22  23 | , |  31  32  33 |
     * 
     * where trans(A11 A21 A31) is n x (k*4) and factored. assume trans(A22 A32)     * is the active block of 4 columns. The factorization is accomplished
     * by:
     * 
     * Step 1. Factor trans(A22 A32)
     * 
     * Step 2. Update the lower half of A33 with matrix multiplication using
     * trans(A22 A32).
     * 
     * Step 3. repeat
     * 
     * Fill in the lower triangle
     */
    for (j = 1; j <= (*n - 1); j++) {
        for (i = j; i <= (*n - 1); i++) {
            *FAC (j - 1, i) = imsl_c_conjg (*FAC (i, j - 1));
        }
    }

    info = 0;

    nrcfac = mod (*n, 4);
    for (j = 1; j <= (*n - nrcfac); j += 4) {
        /* prepare j-th column */
        if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= small) {
            info = j;
            goto L_140;
        }
        *FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
            zero);
        r0 = imsl_fc_convert (imsl_c_div (imsl_cf_convert (one, 0.), *FAC (j - 1, j - 1)));
        /* Form the j-th mutiplier */
        *FAC (j - 1, j) = imsl_c_mul (imsl_cf_convert (r0, 0.), *FAC (j - 1, j));
        *FAC (j - 1, j + 1) = imsl_c_mul (imsl_cf_convert (r0, 0.), *FAC (j - 1, j + 1));
        /* update columns j+1 thru j+3 */
        t1 = imsl_c_neg (imsl_c_conjg (*FAC (j - 1, j)));
        t2 = imsl_c_neg (imsl_c_conjg (*FAC (j - 1, j + 1)));
        t3 = imsl_c_neg (imsl_c_conjg (imsl_c_mul (imsl_cf_convert (r0, 0.), *FAC (j - 1, j + 2))));
        *FAC (j, j) = imsl_c_add (*FAC (j, j), imsl_c_mul (*FAC (j - 1, j), t1));
        *FAC (j, j + 1) = imsl_c_add (*FAC (j, j + 1), imsl_c_mul (*FAC (j - 1, j + 1),
                t1));
        *FAC (j + 1, j + 1) = imsl_c_add (*FAC (j + 1, j + 1), imsl_c_mul (*FAC (j - 1, j + 1),
                t2));
        for (i = 3; i <= (*n - j); i++) {
            *FAC (j - 1, j + i - 1) = imsl_c_mul (imsl_cf_convert (r0, 0.), *FAC (j - 1, j + i - 1));
            *FAC (j, j + i - 1) = imsl_c_add (*FAC (j, j + i - 1), imsl_c_mul (*FAC (j - 1, j + i - 1),
                    t1));
            *FAC (j + 1, j + i - 1) = imsl_c_add (*FAC (j + 1, j + i - 1), imsl_c_mul (*FAC (j - 1, j + i - 1),
                    t2));
            *FAC (j + 2, j + i - 1) = imsl_c_add (*FAC (j + 2, j + i - 1), imsl_c_mul (*FAC (j - 1, j + i - 1),
                    t3));
        }
        /* prepare (j+1)-th column */
        if (imsl_fc_convert (*FAC (j, j)) <= small) {
            info = j + 1;
            goto L_140;
        }
        *FAC (j, j) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j, j))), zero);
        r1 = imsl_fc_convert (imsl_c_div (imsl_cf_convert (one, 0.), *FAC (j, j)));
        *FAC (j, j + 1) = imsl_c_mul (imsl_cf_convert (r1, 0.), *FAC (j, j + 1));
        /*
         * form j+1 multiplier, and update columns j+2 thru j+3
         */
        t2 = imsl_c_neg (imsl_c_conjg (*FAC (j, j + 1)));
        t3 = imsl_c_neg (imsl_c_conjg (imsl_c_mul (imsl_cf_convert (r1, 0.), *FAC (j, j + 2))));
        *FAC (j + 1, j + 1) = imsl_c_add (*FAC (j + 1, j + 1), imsl_c_mul (*FAC (j, j + 1),
                t2));
        for (i = 3; i <= (*n - j); i++) {
            *FAC (j, j + i - 1) = imsl_c_mul (imsl_cf_convert (r1, 0.), *FAC (j, j + i - 1));
            *FAC (j + 1, j + i - 1) = imsl_c_add (*FAC (j + 1, j + i - 1), imsl_c_mul (*FAC (j, j + i - 1),
                    t2));
            *FAC (j + 2, j + i - 1) = imsl_c_add (*FAC (j + 2, j + i - 1), imsl_c_mul (*FAC (j, j + i - 1),
                    t3));
        }
        /* prepare (j+2)-th column */
        if (imsl_fc_convert (*FAC (j + 1, j + 1)) <= small) {
            info = j + 2;
            goto L_140;
        }
        *FAC (j + 1, j + 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j + 1, j + 1))),
            zero);
        r2 = imsl_fc_convert (imsl_c_div (imsl_cf_convert (one, 0.), *FAC (j + 1, j + 1)));
        /*
         * form j+2 multiplier, and update column j+3
         */
        t3 = imsl_c_neg (imsl_c_conjg (imsl_c_mul (imsl_cf_convert (r2, 0.), *FAC (j + 1, j + 2))));
        for (i = 3; i <= (*n - j); i++) {
            *FAC (j + 1, j + i - 1) = imsl_c_mul (*FAC (j + 1, j + i - 1), imsl_cf_convert (r2, 0.));
            *FAC (j + 2, j + i - 1) = imsl_c_add (*FAC (j + 2, j + i - 1), imsl_c_mul (*FAC (j + 1, j + i - 1),
                    t3));
        }
        /* prepare (j+3)-th column */
        if (imsl_fc_convert (*FAC (j + 2, j + 2)) <= small) {
            info = j + 3;
            goto L_140;
        }
        *FAC (j + 2, j + 2) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j + 2, j + 2))),
            zero);
        r3 = imsl_fc_convert (imsl_c_div (imsl_cf_convert (one, 0.), *FAC (j + 2, j + 2)));
        /*
         * Form the (j+3)-th mutiplier and load t(4:nlenb,4) with
         * -imsl_fac(j+3:nlenb,j+3)
         */
        for (i = 4; i <= (*n - j); i++) {
            *FAC (j + 2, j + i - 1) = imsl_c_mul (imsl_cf_convert (r3, 0.), *FAC (j + 2, j + i - 1));
        }
        /*
         * Perform update on the lower triangle on columns j+4 thru n rows
         * j+4 thru n
         */
        kmod2 = mod (*n - j - 3, 2);
        if (kmod2 > 0) {
            *FAC (*n - 1, *n - 1) = imsl_c_sub (imsl_c_sub (imsl_c_sub (imsl_c_sub (*FAC (*n - 1, *n - 1),
                            imsl_c_mul (*FAC (j - 1, *n - 1), imsl_c_conjg (*FAC (j - 1, *n - 1)))),
                        imsl_c_mul (*FAC (j, *n - 1), imsl_c_conjg (*FAC (j, *n - 1)))), imsl_c_mul (*FAC (j + 1, *n - 1),
                        imsl_c_conjg (*FAC (j + 1, *n - 1)))), imsl_c_mul (*FAC (j + 2, *n - 1),
                    imsl_c_conjg (*FAC (j + 2, *n - 1))));
        }
        for (k = *n - kmod2 - 1; k >= (j + 4); k += -2) {
#if 0
            t1 = imsl_c_neg (imsl_c_conjg (*FAC (j - 1, k - 1)));
            t2 = imsl_c_neg (imsl_c_conjg (*FAC (j, k - 1)));
            t3 = imsl_c_neg (imsl_c_conjg (*FAC (j + 1, k - 1)));
            t4 = imsl_c_neg (imsl_c_conjg (*FAC (j + 2, k - 1)));
            s1 = imsl_c_neg (imsl_c_conjg (*FAC (j - 1, k)));
            s2 = imsl_c_neg (imsl_c_conjg (*FAC (j, k)));
            s3 = imsl_c_neg (imsl_c_conjg (*FAC (j + 1, k)));
            s4 = imsl_c_neg (imsl_c_conjg (*FAC (j + 2, k)));
#endif
            t1 = *FAC (j - 1, k - 1); t1.re = -t1.re;
            t2 = *FAC (j, k - 1); t2.re = -t2.re;
            t3 = *FAC (j + 1, k - 1); t3.re = -t3.re;
            t4 = *FAC (j + 2, k - 1); t4.re = -t4.re;
            s1 = *FAC (j - 1, k); s1.re = -s1.re;
            s2 = *FAC (j, k); s2.re = -s2.re;
            s3 = *FAC (j + 1, k); s3.re = -s3.re;
            s4 = *FAC (j + 2, k); s4.re = -s4.re;
            for (i = k; i <= *n; i++) {
#if 0
                *FAC (k - 1, i - 1) = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (k - 1, i - 1),
                                imsl_c_mul (*FAC (j - 1, i - 1), t1)), imsl_c_mul (*FAC (j, i - 1),
                                t2)), imsl_c_mul (*FAC (j + 1, i - 1), t3)), imsl_c_mul (*FAC (j + 2, i - 1),
                        t4));
                *FAC (k, i - 1) = imsl_c_add (imsl_c_add (imsl_c_add (imsl_c_add (*FAC (k, i - 1),
                                imsl_c_mul (*FAC (j - 1, i - 1), s1)), imsl_c_mul (*FAC (j, i - 1),
                                s2)), imsl_c_mul (*FAC (j + 1, i - 1), s3)), imsl_c_mul (*FAC (j + 2, i - 1),
                        s4));
#endif
		ctemp1.re = (FAC(j-1, i-1))->re * t1.re - (FAC(j-1, i-1))->im * t1.im;
		ctemp1.im = (FAC(j-1, i-1))->re * t1.im + (FAC(j-1, i-1))->im * t1.re;
		ctemp2.re = (FAC(j, i-1))->re * t2.re - (FAC(j, i-1))->im * t2.im;
		ctemp2.im = (FAC(j, i-1))->re * t2.im + (FAC(j, i-1))->im * t2.re;
		ctemp3.re = (FAC(j+1, i-1))->re * t3.re - (FAC(j+1, i-1))->im * t3.im;
		ctemp3.im = (FAC(j+1, i-1))->re * t3.im + (FAC(j+1, i-1))->im * t3.re;
		ctemp4.re = (FAC(j+2, i-1))->re * t4.re - (FAC(j+2, i-1))->im * t4.im;
		ctemp4.im = (FAC(j+2, i-1))->re * t4.im + (FAC(j+2, i-1))->im * t4.re;
		(FAC(k-1, i-1))->re += ctemp1.re + ctemp2.re + ctemp3.re + ctemp4.re;
		(FAC(k-1, i-1))->im += ctemp1.im + ctemp2.im + ctemp3.im + ctemp4.im;

                ctemp1.re = (FAC(j-1, i-1))->re * s1.re - (FAC(j-1, i-1))->im * s1.im;
                ctemp1.im = (FAC(j-1, i-1))->re * s1.im + (FAC(j-1, i-1))->im * s1.re;
                ctemp2.re = (FAC(j, i-1))->re * s2.re - (FAC(j, i-1))->im * s2.im;
                ctemp2.im = (FAC(j, i-1))->re * s2.im + (FAC(j, i-1))->im * s2.re;
                ctemp3.re = (FAC(j+1, i-1))->re * s3.re - (FAC(j+1, i-1))->im * s3.im;
                ctemp3.im = (FAC(j+1, i-1))->re * s3.im + (FAC(j+1, i-1))->im * s3.re;
                ctemp4.re = (FAC(j+2, i-1))->re * s4.re - (FAC(j+2, i-1))->im * s4.im;
                ctemp4.im = (FAC(j+2, i-1))->re * s4.im + (FAC(j+2, i-1))->im * s4.re;
                (FAC(k, i-1))->re += ctemp1.re + ctemp2.re + ctemp3.re + ctemp4.re;
                (FAC(k, i-1))->im += ctemp1.im + ctemp2.im + ctemp3.im + ctemp4.im;
            }
        }

    }
    /*
     * Take care of remaining nrcfac columns of imsl_fac
     */
    j = *n - nrcfac + 1;
    if (nrcfac < 3)
        goto L_120;
    if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= small) {
        info = j;
        goto L_140;
    }
    *FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
        zero);
    rtemp = imsl_c_div (imsl_cf_convert (one, 0.), *FAC (j - 1, j - 1));
    *FAC (j - 1, j) = imsl_c_mul (rtemp, *FAC (j - 1, j));
    *FAC (j - 1, j + 1) = imsl_c_mul (rtemp, *FAC (j - 1, j + 1));
    *FAC (j, j) = imsl_c_sub (*FAC (j, j), imsl_c_mul (*FAC (j - 1, j), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j, j + 1) = imsl_c_sub (*FAC (j, j + 1), imsl_c_mul (*FAC (j - 1, j + 1), imsl_c_conjg (*FAC (j - 1, j))));
    *FAC (j + 1, j + 1) = imsl_c_sub (*FAC (j + 1, j + 1), imsl_c_mul (*FAC (j - 1, j + 1),
            imsl_c_conjg (*FAC (j - 1, j + 1))));
    j += 1;
    nrcfac -= 1;
L_120:
    if (nrcfac < 2)
        goto L_130;
    if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= small) {
        info = j;
        goto L_140;
    }
    *FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
        zero);
    rtemp = imsl_c_div (imsl_cf_convert (one, 0.), *FAC (j - 1, j - 1));
    *FAC (j - 1, j) = imsl_c_mul (rtemp, *FAC (j - 1, j));
    *FAC (j, j) = imsl_c_sub (*FAC (j, j), imsl_c_mul (*FAC (j - 1, j), imsl_c_conjg (*FAC (j - 1, j))));
    j += 1;
    nrcfac -= 1;
L_130:
    if (nrcfac < 1)
        goto L_140;
    if (imsl_fc_convert (*FAC (j - 1, j - 1)) <= small) {
        info = j;
        goto L_140;
    }
    *FAC (j - 1, j - 1) = imsl_cf_convert (sqrt (imsl_fc_convert (*FAC (j - 1, j - 1))),
        zero);
L_140:
    if (info != 0) {
        imsl_e1sti (1, info);
/*
        (4, 2, "The leading %(i1) by %(i1) submatrix of the 
	input matrix is not positive definite.");
*/
	imsl_ermes (IMSL_FATAL, IMSL_NONPOSITIVE_MATRIX);
    }
    /* Fill in the upper triangle */
    for (i = 1; i <= (*n - 1); i++) {
        for (j = 1; j <= i; j++) {
#if 0
            *FAC (i, j - 1) = imsl_c_conjg (*FAC (j - 1, i));
#endif
	    (FAC (i, j-1))->re =  (FAC (j-1, i))->re;
	    (FAC (i, j-1))->im = -(FAC (j-1, i))->im;
        }
    }

L_9000:
    imsl_e1pop ("LFTDH");

    return;
}                               /* end of function */

#undef A
#undef FAC

#ifdef ANSI
static void l_l2cdh (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *imsl_fac,
                Mint *ldfac, Mfloat *rcond, Mf_complex *z)
#else
static void l_l2cdh (n, a, lda, imsl_fac, ldfac, rcond, z)
    Mint       *n;
    Mf_complex *a;
    Mint       *lda;
    Mf_complex *imsl_fac;
    Mint       *ldfac;
    Mfloat     *rcond;
    Mf_complex  z[];
#endif
{
#define A(I_,J_)        (a+(I_)*(*lda)+(J_))
#define FAC(I_,J_)      (imsl_fac+(I_)*(*ldfac)+(J_))
    Mint        _l0, _l1, _l2, j, k;
    Mfloat      anorm, s, sm, suma, ynorm;
    Mf_complex  _cx0, ek, t, wk, wkm;


#define L3CCT(zdum)     (Mfloat)(fabs( imsl_fc_convert( (zdum) ) ) + fabs( imsl_c_aimag( (zdum))))
#define L4CCT(zdum1,zdum2)      (imsl_c_mul(imsl_cf_convert(L3CCT( (zdum1) ),0.),(imsl_c_div((zdum2),imsl_cf_convert(L3CCT( (zdum2) ),0.)))))

    imsl_e1psh ("l_l2cdh");

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
	else if (*n > *ldfac) {
	    imsl_e1sti (1, *n);
	    imsl_e1sti (2, *ldfac);

	    /*
	     * (5, 3, "The order of the matrix must be less than or equal to
	     * its leading dimension while N = %(i1) and LDFAC = %(i2) are
	     * given.");
	     */
	    imsl_ermes (IMSL_TERMINAL, IMSL_LDFAC_LESS_ORDER);
	}
	else {
	    /*
	     * COMPUTE 1-NORM OF A USING ONLY THE UPPER TRIANGLE OF A
	     */
	    *rcond = F_ZERO;
	    anorm = F_ZERO;
	    for (j = 1; j <= *n; j++) {
		_l0 = 1;
		suma = imsl_scasum (&j, A (j - 1, 0), &_l0);
		_l0 = *n - j;
		if (j < *n)
		    suma += imsl_scasum (&_l0, A (j, j - 1), lda);
		anorm = imsl_f_max (anorm, suma);
	    }
	    /* FACTOR */
	    l_lftdh (n, a, lda, imsl_fac, ldfac);
	    if (imsl_n1rty (1) < 4) {
		/*
		 * RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
		 * ESTIMATE = NORM(Z)/NORM(Y) WHERE A*Z = Y AND A*Y = E . THE
		 * COMPONENTS OF E ARE CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH
		 * IN THE ELEMENTS OF W WHERE CTRANS(R)*W = E. THE VECTORS
		 * ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. SOLVE
		 * CTRANS(R)*W = E
		 */
		ek = C_ONE;
		_cx0 = imsl_cf_convert (F_ZERO, F_ZERO);
		_l0 = 1;
		imsl_cset (n, &_cx0, z, &_l0);
		for (k = 1; k <= *n; k++) {
		    wk = imsl_c_neg (z[k - 1]);
		    if (L3CCT (z[k - 1]) != F_ZERO)
			ek = L4CCT (ek, wk);

		    if (L3CCT (imsl_c_sub (ek, z[k - 1])) > imsl_fc_convert (*FAC (k - 1, k - 1))) {
			s = imsl_fc_convert (*FAC (k - 1, k - 1)) / L3CCT (imsl_c_sub (ek,
				z[k - 1]));
			_l0 = 1;
			imsl_csscal (n, &s, z, &_l0);
			ek = imsl_c_mul (imsl_cf_convert (s, F_ZERO), ek);
		    }
		    wk = imsl_c_sub (ek, z[k - 1]);
		    wkm = imsl_c_sub (imsl_c_neg (ek), z[k - 1]);
		    s = L3CCT (wk);
		    sm = L3CCT (wkm);
		    wk = imsl_c_div (wk, *FAC (k - 1, k - 1));
		    wkm = imsl_c_div (wkm, *FAC (k - 1, k - 1));

		    if (k + 1 <= *n) {
			for (j = k + 1; j <= *n; j++) {
			    sm += L3CCT (imsl_c_add (z[j - 1], imsl_c_mul (wkm, imsl_c_conjg (*FAC (j - 1, k - 1)))));
			    z[j - 1] = imsl_c_add (z[j - 1], imsl_c_mul (wk, imsl_c_conjg (*FAC (j - 1, k - 1))));
			    s += L3CCT (z[j - 1]);
			}

			if (s < sm) {
			    t = imsl_c_sub (wkm, wk);
			    wk = wkm;
			    for (j = k + 1; j <= *n; j++) {
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
		/* SOLVE R*Y = W */
		for (k = *n; k >= 1; k--) {
		    if (L3CCT (z[k - 1]) > imsl_fc_convert (*FAC (k - 1, k - 1))) {
			s = imsl_fc_convert (*FAC (k - 1, k - 1)) / L3CCT (z[k - 1]);
			_l0 = 1;
			imsl_csscal (n, &s, z, &_l0);
		    }
		    z[k - 1] = imsl_c_div (z[k - 1], *FAC (k - 1, k - 1));
		    t = imsl_c_neg (z[k - 1]);
		    _l0 = k - 1;
		    _l1 = 1;
		    _l2 = 1;
		    imsl_caxpy (&_l0, &t, FAC (k - 1, 0), &_l1,
			&z[0], &_l2);
		}
		_l0 = 1;
		s = F_ONE / imsl_scasum (n, z, &_l0);
		_l0 = 1;
		imsl_csscal (n, &s, z, &_l0);

		ynorm = F_ONE;
		/* SOLVE CTRANS(R)*V = Y */
		for (k = 1; k <= *n; k++) {
		    _l0 = k - 1;
		    _l1 = 1;
		    _l2 = 1;
		    z[k - 1] = imsl_c_sub (z[k - 1], imsl_cdotc (&_l0, FAC (k - 1, 0), &_l1, &z[0], &_l2));
		    if (L3CCT (z[k - 1]) > imsl_fc_convert (*FAC (k - 1, k - 1))) {
			s = imsl_fc_convert (*FAC (k - 1, k - 1)) / L3CCT (z[k - 1]);
			_l0 = 1;
			imsl_csscal (n, &s, z, &_l0);
			ynorm *= s;
		    }
		    z[k - 1] = imsl_c_div (z[k - 1], *FAC (k - 1, k - 1));
		}
		_l0 = 1;
		s = F_ONE / imsl_scasum (n, z, &_l0);
		_l0 = 1;
		imsl_csscal (n, &s, z, &_l0);
		ynorm *= s;
		/* SOLVE R*Z = V */
		for (k = *n; k >= 1; k--) {
		    if (L3CCT (z[k - 1]) > imsl_fc_convert (*FAC (k - 1, k - 1))) {
			s = imsl_fc_convert (*FAC (k - 1, k - 1)) / L3CCT (z[k - 1]);
			_l0 = 1;
			imsl_csscal (n, &s, z, &_l0);
			ynorm *= s;
		    }
		    z[k - 1] = imsl_c_div (z[k - 1], *FAC (k - 1, k - 1));
		    t = imsl_c_neg (z[k - 1]);
		    _l0 = k - 1;
		    _l1 = 1;
		    _l2 = 1;
		    imsl_caxpy (&_l0, &t, FAC (k - 1, 0), &_l1,
			&z[0], &_l2);
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
		     * (3, 1, "The matrix is algorithmically singular.  An
		     * estimate of the reciprocal of its L1 condition number
		     * is RCOND = %(r1).");
		     */
		    imsl_ermes (IMSL_WARNING, IMSL_ILL_CONDITIONED);
		}
	    }
	}
    }

    imsl_e1pop ("l_l2cdh");
    return;
}				/* end of function */
#undef  L4CCT
#undef  L3CCT
#undef  A
#undef  FAC





/*
-----------------------------------------------------------------------

  IMSL Name:  LSLCT/DLSLCT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 25, 1987

    Purpose:    Solve a d_complex triangular system of linear equations.

    Usage:      CALL LSLCT (N, A, LDA, B, IPATH, X)

    Arguments:
       N      - Number of equations.  (Input)
       A      - Complex N by N matrix containing the coefficient matrix
                of the triangular linear system.  (Input)
                For a lower triangular system, only the lower triangle
                of A is referenced.  For an upper triangular system, only
                the upper triangle of A is referenced.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       B      - Complex vector of length N containing the right-hand side
                of the linear system.  (Input)
       IPATH  - Path indicator.  (Input)
                IPATH = 1 means solve A*X = B, A lower triangular
                IPATH = 2 means solve A*X = B, A upper triangular
                IPATH = 3 means solve ctrans(A)*X = B, A lower triangular
                IPATH = 4 means solve ctrans(A)*X = B, A upper triangular
                        where ctrans(A) is the conjugate transpose of A.
       X      - Complex vector of length N containing the solution to the
                linear system.  (Output)
                If B is not needed, B and X can share the same storage
                locations.

    Remark:
       Informational error
       Type Code
         4   1  The input triangular matrix is singular.  Some of its
                diagonal elements are near zero.

    GAMS:       D2c3

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_lslct (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *b,
                Mint *ipath, Mf_complex *x)
#else
static void l_lslct (n, a, lda, b, ipath, x)
    Mint       *n;
    Mf_complex *a;
    Mint       *lda;
    Mf_complex  b[];
    Mint       *ipath;
    Mf_complex  x[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
    Mint        _l0, _l1, i;
    Mfloat      big, small;


#define L3CCT(zdum)	(Mfloat)(fabs( imsl_fc_convert( (zdum) ) ) + fabs( imsl_c_aimag( (zdum) ) ))

    imsl_e1psh ("l_lslct");

    /* CHECK FOR ZERO DIAGONAL ELEMENTS. */
    small = imsl_amach (1);
    big = imsl_amach (2);
    if (small * big < F_ONE)
	small = F_ONE / big;
    for (i = 1; i <= *n; i++) {
	if (L3CCT (*A (i - 1, i - 1)) < small) {
	    imsl_e1sti (1, i-1);
	    /*
	     * (4, 1, "The input triangular matrix is singular.  The index of
	     * the first zero diagonal element is equal to %(i1).");
	     */
	    imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_TRI_MATRIX);
	    goto L_9000;
	}
    }
    /* MAKE A COPY OF B IN X AND WORK WITH X */
    _l0 = 1;
    _l1 = 1;
    imsl_ccopy (n, b, &_l0, x, &_l1);
    /* Solve A*X=B for a lower triangular */
    if (*ipath == 1) {
	_l0 = 1;
	l_ctrsv ("L", sizeof ("L"), "N", sizeof ("N"), "N", sizeof ("N"
	    ), n, a, lda, x, &_l0);
	/* Solve A*X=B for an upper triangular. */
    }
    else if (*ipath == 2) {
	_l0 = 1;
	l_ctrsv ("U", sizeof ("U"), "N", sizeof ("N"), "N", sizeof ("N"
	    ), n, a, lda, x, &_l0);
	/*
	 * Solve CTRANS(A)*X=B for a lower triangular.
	 */
    }
    else if (*ipath == 3) {
	_l0 = 1;
	l_ctrsv ("L", sizeof ("L"), "C", sizeof ("C"), "N", sizeof ("N"
	    ), n, a, lda, x, &_l0);
	/*
	 * Solve CTRANS(A)*X=B for an upper triangular.
	 */
    }
    else if (*ipath == 4) {
	_l0 = 1;
	l_ctrsv ("U", sizeof ("U"), "C", sizeof ("C"), "N", sizeof ("N"
	    ), n, a, lda, x, &_l0);
    }
    else {
	imsl_e1sti (1, *ipath);

	/*
	 * (5, 3, "IPATH must be either 1, 2, 3 or 4 while a value of %(i1)
	 * is given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_IPATH_RANGE_4);
    }

L_9000:
    imsl_e1pop ("l_lslct");
    return;
}				/* end of function */
#undef	L3CCT
#undef	A
/*
-----------------------------------------------------------------------

  IMSL Name:  CHFCG/DCHFCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 23, 1985

    Purpose:    Extend a d_complex Hermitian matrix defined in its upper
                triangle to its lower triangle.

    Usage:      CALL CHFCG (N, A, LDA)

    Arguments:
       N      - Order of the matrix.  (Input)
       A      - Complex Hermitian matrix of order N.  (Input/Output)
                On input, the upper triangle of A defines a Hermitian
                matrix.  On output, the lower triangle of A is defined
                so that A is Hermitian.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remark:
       Informational errors
       Type Code
         3   1  The matrix is not Hermitian.  It has a diagonal entry
                with a small imaginary part.
         4   2  The matrix is not Hermitian.  It has a diagonal entry
                with an imaginary part.

    Keywords:   Basic matrix operation; Matrix conversion

    GAMS:       D1b9

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_chfcg (Mint *n, Mf_complex *a, Mint *lda)
#else
static void l_chfcg (n, a, lda)
    Mint       *n;
    Mf_complex *a;
    Mint       *lda;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
    Mint        i, j;
    Mfloat      eps;
    imsl_e1psh ("l_chfcg");
    /* Check that A is Hermitian */
    eps = F_TEN * imsl_amach (4);
    for (i = 1; i <= *n; i++) {
	if (fabs (imsl_c_aimag (*A (i - 1, i - 1))) != F_ZERO) {
	    if (fabs (imsl_c_aimag (*A (i - 1, i - 1))) > eps * fabs (imsl_fc_convert (*A (i - 1, i - 1)))) {
		imsl_e1sti (1, i-1);
		imsl_e1stc (1, *A (i - 1, i - 1));
		/*
		 * (4, 2, "The matrix element A(%(i1),%(i1)) = %(c1).  The
		 * diagonal of a Hermitian matrix must be real.");
		 */
		imsl_ermes (IMSL_FATAL, IMSL_HERMITIAN_DIAG_REAL);
		goto L_9000;
	    }
	    else {
		imsl_e1sti (1, i-1);
		imsl_e1stc (1, *A (i - 1, i - 1));
		/*
		 * (3, 1, "The matrix element A(%(i1),%(i1)) = %(c1).  The
		 * diagonal of a Hermitian matrix must be real; its imaginary
		 * part is set to zero.");
		 */
		imsl_ermes (IMSL_WARNING, IMSL_HERMITIAN_DIAG_REAL_1);
		*A (i - 1, i - 1) = imsl_cf_convert (imsl_fc_convert (*A (i - 1, i - 1)), F_ZERO);
	    }
	}
    }

    for (j = 1; j <= (*n - 1); j++) {
	for (i = j + 1; i <= *n; i++) {
	    *A (j - 1, i - 1) = imsl_c_conjg (*A (i - 1, j - 1));
	}
    }

L_9000:
    ;
    imsl_e1pop ("l_chfcg");
    return;
}				/* end of function */
#undef  A
/*
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
static void l_ctrsv (Mchar *uplo, unsigned uplo_s, Mchar *trans, unsigned trans_s,
                Mchar *diag, unsigned diag_s, Mint *n, Mf_complex *a,
                Mint *lda, Mf_complex *x, Mint *incx)
#else
static void l_ctrsv (uplo, uplo_s, trans, trans_s, diag, diag_s,
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
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
    Mlong        imsl_ctran, ndiag, tran, upper;
    Mint        _l0, _l1, i, ix;


    upper = imsl_l1ame (uplo, uplo_s, "U", sizeof ("U"));
    ndiag = imsl_l1ame (diag, diag_s, "N", sizeof ("N"));
    tran = imsl_l1ame (trans, trans_s, "T", sizeof ("T"));
    imsl_ctran = imsl_l1ame (trans, trans_s, "C", sizeof ("C"));
    /*
     * Test the input parameters.
     */
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

    return;
}				/* end of function */
#undef  A
