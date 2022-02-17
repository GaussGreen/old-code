#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static VA_LIST_HACK PROTO (l_lin_sol_nonnegdef, (Mint n, Mfloat *a, Mfloat *b,
	            va_list argptr));
    static void PROTO (l_chfac, (Mint *n, Mfloat *a, Mint *lda, Mfloat *tol,
	            Mint *irank, Mfloat *r, Mint *ldr));
    static void PROTO (l_girts, (Mint *n, Mfloat *r, Mint *ldr, Mint *nb, Mfloat *b,
	            Mint *ldb, Mint *ipath, Mint *irank, Mfloat *x,
	            Mint *ldx, Mfloat *rinv, Mint *ldrinv));
    static void PROTO (l_c1r, (Mint *n, Mfloat *r, Mint *ldr, Mint *ner));
    static void PROTO (l_c1trg, (Mint *n, Mfloat *a, Mint *lda));
    static Mfloat PROTO (l_a1ot, (Mint *n, Mfloat *sx, Mint *incx, Mfloat *sy,
	            Mint *incy));

    static Mfloat *lv_x;
#ifdef ANSI
    Mfloat     *imsl_f_lin_sol_nonnegdef (Mint n, Mfloat *a, Mfloat *b,...)
#else
    Mfloat     *imsl_f_lin_sol_nonnegdef (n, a, b, va_alist)
    Mint        n;
    Mfloat     *a, *b;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, b);

    E1PSH ("imsl_f_lin_sol_nonnegdef", "imsl_d_lin_sol_nonnegdef");

    lv_x = NULL;
    IMSL_CALL (l_lin_sol_nonnegdef (n, a, b, argptr));
    va_end (argptr);

    E1POP ("imsl_f_lin_sol_nonnegdef", "imsl_d_lin_sol_nonnegdef");

    return lv_x;
}
#ifdef ANSI
static VA_LIST_HACK l_lin_sol_nonnegdef (Mint n, Mfloat *a, Mfloat *b, va_list
                argptr)
#else
static VA_LIST_HACK l_lin_sol_nonnegdef (n, a, b, argptr)
    Mint        n;
    Mfloat     *a, *b;
    va_list     argptr;
#endif
{
    Mint        a_col_dim = n;
    Mfloat    **p_factor;
    Mfloat     *factor = NULL;
    Mint        fac_col_dim = n;
    Mfloat    **p_inva;
    Mfloat     *inva = NULL;
    Mint        inva_col_dim = n;
    Mfloat      tol;

    Mint        factor_only = 0;
    Mint        inverse_only = 0;
    Mint        solve_only = 0;
    Mint        user_result = 0;
    Mint        user_factor = 0;
    Mint        return_factor = 0;
    Mint        user_inverse = 0;
    Mint        return_inverse = 0;
    Mint        shared_memory = 0;

    Mint         code = 1;
    Mint         arg_number = 3;
    Mint         error = 0;

    Mint        irank;
    Mint        ldr;
    Mint        i;

    tol = 100.0 * imsl_amach (4);
    lv_x = NULL;

    while (code > 0) {
	code = va_arg (argptr, Mint);
	++arg_number;
	switch (code) {
	case IMSL_RETURN_USER:
	    lv_x = va_arg (argptr, Mfloat *);
	    if (!lv_x) {
		imsl_e1stl (1, "x");
		imsl_e1stl (2, "IMSL_RETURN_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    user_result = 1;
	    arg_number++;
	    break;
	case IMSL_A_COL_DIM:
	    a_col_dim = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_FACTOR:
	    p_factor = va_arg (argptr, Mfloat **);
	    *p_factor = NULL;
	    return_factor = 1;
	    arg_number++;
	    break;
	case IMSL_FACTOR_USER:
	    factor = va_arg (argptr, Mfloat *);
	    if (!factor) {
		imsl_e1stl (1, "factor");
		imsl_e1stl (2, "IMSL_FACTOR_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    return_factor = 1;
	    user_factor = 1;
	    arg_number++;
	    break;
	case IMSL_FAC_COL_DIM:
	    fac_col_dim = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_INVERSE:
	    p_inva = va_arg (argptr, Mfloat **);
	    *p_inva = NULL;
	    return_inverse = 1;
	    arg_number++;
	    break;
	case IMSL_INVERSE_USER:
	    inva = va_arg (argptr, Mfloat *);
	    if (!inva) {
		imsl_e1stl (1, "inva");
		imsl_e1stl (2, "IMSL_INVERSE_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    return_inverse = 1;
	    user_inverse = 1;
	    arg_number++;
	    break;
	case IMSL_INV_COL_DIM:
	    inva_col_dim = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_TOLERANCE:
	    tol = (Mfloat) va_arg (argptr, Mdouble);
	    arg_number++;
	    break;
	case IMSL_TOLERANCE_ADR:
	    tol = *( va_arg (argptr, Mfloat *));
	    arg_number++;
	    break;
	case IMSL_FACTOR_ONLY:
	    factor_only = 1;
	    arg_number++;
	    break;
	case IMSL_SOLVE_ONLY:
	    solve_only = 1;
	    arg_number++;
	    break;
	case IMSL_INVERSE_ONLY:
	    inverse_only = 1;
	    arg_number++;
	    break;
	case 0:
	    break;
	default:
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, arg_number);
	    imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
	    return argptr;
	}
    }
    if (!a && !solve_only) {
	imsl_e1stl (1, "a");
	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
	++error;
    }
    if (!b && (!inverse_only && !factor_only)) {
	imsl_e1stl (1, "b");
	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
	++error;
    }

    if (error)
	return argptr;
    if (n < 1) {
	imsl_e1sti (1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
	++error;
    }
    if (solve_only + factor_only + inverse_only > 1) {
	imsl_ermes (IMSL_TERMINAL, IMSL_BAD_SOLVE_FACTOR_INVERSE);
	++error;
    }
    if (error)
	return argptr;

    if (solve_only) {
	if (!user_factor) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_SPECIFY_SOLVE_ONLY);
	    ++error;
	}
	else if (return_inverse) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_BAD_SOLVE_ONLY);
	    ++error;
	}
    }
    else if (factor_only) {
	if (!return_factor) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_FACTOR_ONLY_SPECIFIER);
	    ++error;
	}
	if (return_inverse) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_BAD_FACTOR_ONLY);
	    ++error;
	}
    }
    else if (inverse_only) {
	if (!return_inverse) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_INVERSE_ONLY_SPECIFIER);
	    ++error;
	}
	if (return_factor) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_BAD_INVERSE_ONLY);
	    ++error;
	}
    }
    if (error)
	return argptr;

    if (return_inverse && inva_col_dim < n) {
	imsl_e1sti (1, inva_col_dim);
	imsl_e1stl (1, "inva");
	imsl_e1sti (2, n);
	imsl_e1stl (2, "n");
	imsl_ermes (IMSL_TERMINAL, IMSL_COLUMN_DIM_ERROR);
	++error;
    }
    else if (!return_inverse) {
	inva_col_dim = n;
    }
    if (return_factor && fac_col_dim < n) {
	imsl_e1sti (1, fac_col_dim);
	imsl_e1stl (1, "factor");
	imsl_e1sti (2, n);
	imsl_e1stl (2, "n");
	imsl_ermes (IMSL_TERMINAL, IMSL_COLUMN_DIM_ERROR);
	++error;
    }
    else if (!return_factor) {
	if (return_inverse) {
	    fac_col_dim = inva_col_dim;
	}
	else {
	    fac_col_dim = n;
	}
    }
    ldr = fac_col_dim;
    if (error)
	return argptr;

    if (!solve_only && !error) {
	if (!user_factor) {
	    if (!(factor = (Mfloat *) imsl_malloc (n * fac_col_dim * sizeof (Mfloat)))) {
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_SPACE);
		++error;
	    }
	}
	if (!error) {
	    Mint        lda = a_col_dim;
	    l_chfac (&n, a, &lda, &tol, &irank, factor, &ldr);
	    error = (imsl_n1rty (1) > 3) ? 1 : 0;
	}
    }

    if (!factor_only && !error) {
	if (!user_result && !inverse_only) {
	    if (!(lv_x = (Mfloat *) imsl_malloc (n * sizeof (Mfloat)))) {
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_SPACE);
		++error;
	    }
	}
	if (!error) {
	    Mint        nb = 1;
	    Mint        ldb = n;
	    Mint        ldx = n;
	    Mint        ldrinv = inva_col_dim;
	    Mint        ipath, jpath;
	    Mint        j;
	    Mfloat     *rinv, *rhs;
	    if (return_inverse) {
		ipath = 4;
		if (!inverse_only) {
		    jpath = 2;
		    rhs = b;
		}
		else {
		    rhs = NULL;
		    jpath = 3;
		    nb  = 0;
		}
		if (!user_inverse) {
		    if (!return_factor) {
			inva = factor;
			ldrinv = fac_col_dim;
			shared_memory = 1;
		    }
		    else {
			inva = (Mfloat *) imsl_malloc (n * inva_col_dim * sizeof (Mfloat));
			if (!inva)
			    ++error;
		    }
		}
		rinv = inva;
	    }
	    else {
		rhs = b;
		ipath = 2;
		jpath = 0;
		rinv = factor;
		ldrinv = fac_col_dim;
	    }
	    if (!error) {
		while (ipath > jpath && !error) {
		    l_girts (&n, factor, &ldr, &nb, rhs, &ldb, &ipath, &irank, lv_x, &ldx, rinv, &ldrinv);
		    error = (imsl_n1rty (1) > 3) ? 1 : 0;
		    --ipath;
		    rhs = lv_x;
		}
		if (error) {
		    if (return_inverse && !user_inverse) {
			if (inva != NULL) imsl_free (inva);
			inva = NULL;
		    }
		}
		else if (return_inverse) {
		    for (i = 0; i < n; i++) {
			for (j = i; j < n; j++) {
			    inva[i * inva_col_dim + j] =
				imsl_sdot (n - j, (inva + j * inva_col_dim + i), inva_col_dim,
				(inva + j * inva_col_dim + j), inva_col_dim);
			}
			scopy (n - i - 1, (inva + i * inva_col_dim + i + 1), 1, (inva + (i + 1) * inva_col_dim + i), inva_col_dim);
		    }
		    if (!user_inverse) {
			for (i = n; i < inva_col_dim; i++) {
			    sset (n, 0.0, (inva + i), inva_col_dim);
			}
			*p_inva = inva;
			inva = NULL;
		    }
		}
	    }
	}
    }

    if (error) { 
	if (!user_factor && !shared_memory) {
	    if (factor)
		imsl_free (factor);
	    factor = NULL;
	}
    } 
    else if (return_factor && !solve_only) {
	for (i = 1; i < n; i++) {
	    scopy (i, (factor + i * fac_col_dim), 1, (factor + i), fac_col_dim);
	}
	if (!user_factor) {
	    for (i = n; i < fac_col_dim; i++) {
		sset (n, 0.0, (factor + i), fac_col_dim);
	    }
	    *p_factor = factor;
	    factor = NULL;
	}
    }

    if (!user_result) {
	if (error || inverse_only) {
	    if (lv_x != NULL) imsl_free (lv_x);
	    lv_x = NULL;
	}
    }
    if (!user_factor && !shared_memory) {
        if (factor)
       	    imsl_free (factor);
        factor = NULL;
    }

    return argptr;
}

/* Structured by FOR_STRUCT, v0.2, on 08/28/90 at 16:28:49
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CHFAC/DCHFAC (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 9, 1985

    Purpose:    Compute an upper triangular factorization of a real
                symmetric nonnegative definite matrix.

    Usage:      CALL CHFAC (N, A, LDA, TOL, IRANK, R, LRD)

    Arguments:
       N      - Order of the matrix.  (Input)
       A      - N by N symmetric nonnegative definite matrix for which an
                upper triangular factorization is desired.  (Input)
                Only elements in the upper triangle of A are referenced.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement in the calling program.  (Input)
       TOL    - Tolerance used in determining linear dependence.  (Input)
                For CHFAC, TOL = 100.*AMACH(4) is a common choice.   For
                DCHFAC, TOL = 100.D0*DMACH(4) is a common choice.  See
                documentation for IMSL routine AMACH/DMACH.
       IRANK  - Rank of A.  (Output)
                N - IRANK is the number of effective zero pivots.
       R      - N by N upper triangular matrix containing the R matrix
                from a Cholesky decomposition trans(R)*R of A.  (Output)
                The elements of the appropriate rows of R are set to 0.0
                if linear dependence of the columns of A is declared.
                (There are N-IRANK rows of R whose elements are set
                to 0.0.)  If A is not needed, then R and A can share
                the same storage locations.
       LDR    - Leading dimension of R exactly as specified in the
                dimension statement in the calling program.  (Input)

    Remarks:
    1. Informational error
       Type Code
         3   1  The input matrix is not nonnegative definite within the
                tolerance defined by TOL.

    2. Elements of row I of R are set to 0.0 if a linear dependence is
       declared.  Linear dependence is declared if
                      I-1
         ABS(A(I,I) - SUM R(J,I)**2) .LE. TOL*ABS(A(I,I))
                      J=1

    Keywords:   Cholesky decomposition; Square root method

    GAMS:       D9; D2b1b

    Chapter:    STAT/LIBRARY Mathematical Support

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_chfac (Mint *n, Mfloat *a, Mint *lda, Mfloat *tol, Mint *irank,
                Mfloat *r, Mint *ldr)
#else
static void l_chfac (n, a, lda, tol, irank, r, ldr)
    Mint       *n;
    Mfloat      a[];
    Mint       *lda;
    Mfloat     *tol;
    Mint       *irank;
    Mfloat      r[];
    Mint       *ldr;
#endif
{
    Mint        info, j, k, ner;
    Mfloat      s, t, x;


    imsl_e1psh ("l_chfac");
    /* CHECK FOR TERMINAL ERRORS */
    ner = 1;

    imsl_c1dim (1, *n, "n", *lda, "lda", &ner);

    imsl_c1dim (1, *n, "*n", *ldr, "ldr", &ner);
    if (*tol < 0.0 || *tol > 1.0) {
	imsl_e1str (1, *tol);
	imsl_ermes (IMSL_TERMINAL, IMSL_TOL_OUT_OF_RANGE);
    }
    if (imsl_n1rty (0) != 0)
	goto L_9000;
    info = 0;
    /* PRESERVE A COPY OF THE INPUT MATRIX */
    for (j = 1; j <= *n; j++) {
	scopy (j, &a[*lda * (j - 1)], 1, &r[*ldr * (j - 1)], 1);
    }

    *irank = 0;
    for (j = 1; j <= *n; j++) {
	s = 0.0e0;
	x = *tol * sqrt (fabs (r[j + *ldr * (j - 1) - 1]));
	for (k = 1; k <= (j - 1); k++) {
	    t = r[k + *ldr * (j - 1) - 1] - imsl_sdot (k - 1, &r[*ldr * (k - 1)],
		1, &r[*ldr * (j - 1)], 1);
	    if (r[k + *ldr * (k - 1) - 1] != 0.0) {
		t /= r[k + *ldr * (k - 1) - 1];
		r[k + *ldr * (j - 1) - 1] = t;
		s += t * t;
	    }
	    else {
		if (info == 0) {
		    if (fabs (t) > x * imsl_snrm2 (k - 1, &r[*ldr * (k - 1)],
			    1)) {
			info = j;
		    }
		}
		r[k + *ldr * (j - 1) - 1] = 0.0;
	    }
	}
	s = r[j + *ldr * (j - 1) - 1] - s;
	if (fabs (s) <= *tol * fabs (r[j + *ldr * (j - 1) - 1])) {
	    s = 0.0e0;
	}
	else if (s < 0.0e0) {
	    s = 0.0e0;
	    if (info == 0)
		info = j;
	}
	else {
	    *irank += 1;
	}
	r[j + *ldr * (j - 1) - 1] = sqrt (s);
    }

    if (info != 0) {
	imsl_e1sti (1, info);
	imsl_e1str (1, *tol);
	imsl_ermes (IMSL_WARNING, IMSL_NOT_NONNEG_DEFINITE);
    }
    /* FILL LOWER TRIANGLE OF R WITH ZEROS */
    l_c1trg (n, r, ldr);
    /* EXIT SECTION */
L_9000:
    imsl_e1pop ("l_chfac");
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  GIRTS/DGIRTS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 10, 1986

    Purpose:    Solve a triangular (possibly singular) set of linear
                systems and/or compute a generalized inverse of an upper
                triangular matrix.

    Usage:      CALL GIRTS (N, R, LDR, NB, B, LDB, IPATH, IRANK, X, LDX,
                            RINV, LDRINV)

    Arguments:
       N      - Order of the upper triangular matrix R.  (Input)
       R      - N by N upper triangular matrix.  (Input)
                If R contains a zero along the diagonal, the remaining
                elements of the row must also be zero.  Only
                the upper triangle of R is referenced.
       LDR    - Leading dimension of R exactly as specified in the
                dimension statement of the calling program.  (Input)
       NB     - Number of columns in B.  (Input)
                NB must be nonnegative.  If NB is zero, no linear systems
                are solved.
       B      - N by NB matrix containing the right hand sides of the
                linear system.  (Input, if NB .GT. 0)
                If NB = 0, B is not referenced and can be a vector of
                length one.
       LDB    - Leading dimension of B exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPATH  - Path option.  (Input)
                IPATH  Action
                  1    Solve R*X = B.
                  2    Solve trans(R)*X = B.
                  3    Solve R*X = B and compute RINV.
                  4    Solve trans(R)*X = B and compute RINV.
       IRANK  - Rank of R.  (Output)
       X      - N by NB matrix containing the solution matrix
                corresponding to the right hand side B.  (Output, if
                NB .GT.0)
                If B is not needed, then X and B can share the same
                storage locations.  If NB = 0, X is not referenced
                and can be a vector of length one.
       LDX    - Leading dimension of X exactly as specified in the
                dimension statement of the calling program.  (Input)
       RINV   - N by N upper triangular matrix that is the inverse of R
                when R is nonsingular.  (Output, if IPATH equals 3 or 4)
                (When R is singular, RINV is a g3 inverse.  See Remark 3
                for an explanation of g3 inverses.)   If IPATH = 1 or 2,
                RINV is not referenced and can be a vector of length one.
                If IPATH = 3 or 4 and R is not needed, then R and RINV
                can share the same storage locations.
       LDRINV - Leading dimension of RINV exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remarks:
    1. Informational error
       Type Code
         3   1  The linear system of equations is inconsistent.

    2. GIRTS assumes that a singular R is represented by zero rows
       in R.  No other forms of singularity in R are allowed.

    3. RINV is a g3 inverse means that it satisfies conditions 1, 2, and
       3 for the Moore-Penrose inverse but generally fails condition 4.
       The four conditions for AINV to be a Moore-Penrose inverse of A
       are as follows:
          1.  A*AINV*A = A
          2.  AINV*A*AINV = AINV
          3.  A*AINV is symmetric
          4.  AINV*A is symmetric.

    GAMS:       D9

    Chapter:    STAT/LIBRARY Mathematical Support

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_girts (Mint *n, Mfloat *r, Mint *ldr, Mint *nb, Mfloat *b,
                Mint *ldb, Mint *ipath, Mint *irank, Mfloat *x,
                Mint *ldx, Mfloat *rinv, Mint *ldrinv)
#else
static void l_girts (n, r, ldr, nb, b, ldb, ipath, irank, x, ldx,
                rinv, ldrinv)
    Mint       *n;
    Mfloat      r[];
    Mint       *ldr, *nb;
    Mfloat      b[];
    Mint       *ldb, *ipath, *irank;
    Mfloat      x[];
    Mint       *ldx;
    Mfloat      rinv[];
    Mint       *ldrinv;
#endif
{
    Mint        _l0, _l1, _l2, i, j, k, ner;
    Mfloat      temp1, temp2;


    imsl_e1psh ("l_girts");
    /* Check for terminal errors */
    ner = 1;

    _l0 = 1;
    imsl_c1dim (_l0, *n, "n", *ldr, "ldr", &ner);
    _l0 = 0;
    _l1 = -1;
    imsl_c1iarg (*nb, "nb", _l0, _l1, &ner);

    _l0 = 1;
    _l1 = 4;
    imsl_c1iarg (*ipath, "ipath", _l0, _l1, &ner);
    if (*nb > 0) {

	_l0 = 1;
	imsl_c1dim (_l0, *n, "*n", *ldb, "ldb", &ner);

	_l0 = 1;
	imsl_c1dim (_l0, *n, "*n", *ldx, "ldx", &ner);
    }
    else {
	ner += 6;
    }
    if (*ipath == 3 || *ipath == 4) {

	_l0 = 1;
	imsl_c1dim (_l0, *n, "*n", *ldrinv, "ldrinv", &ner);
    }
    if (imsl_n1rty (0) != 0)
	goto L_9000;
    /*
     * Linear dependent rows of R must be represented only by rows whose
     * elements are zero.
     */
    l_c1r (n, r, ldr, &ner);
    if (imsl_n1rty (0) != 0)
	goto L_9000;
    /* Get rank */
    *irank = 0;
    for (i = 1; i <= *n; i++) {
	if (r[i + *ldr * (i - 1) - 1] != 0.0)
	    *irank += 1;
    }
    /*
     * Make A copy of B in X and work with X
     */
    for (j = 1; j <= *nb; j++) {
	scopy (*n, &b[*ldb * (j - 1)], 1, &x[*ldx * (j - 1)], 1);
    }

    if (*ipath == 1 || *ipath == 3) {
	/* Solve R*X = B */
	if (*irank < *n) {
	    for (i = 1; i <= *nb; i++) {
		for (j = *n; j >= 1; j--) {
		    if (r[j + *ldr * (j - 1) - 1] == 0.0) {
			if (x[j + *ldx * (i - 1) - 1] != 0.0) {
			    imsl_e1sti (1, j-1);
			    imsl_e1sti (2, i-1);
			    imsl_e1str (1, x[j + *ldx * (i - 1) - 1]);
			    imsl_ermes (IMSL_WARNING, IMSL_INCONSISTENT_EQUATIONS);
			}
			x[j + *ldx * (i - 1) - 1] = 0.0;
		    }
		    else {
			x[j + *ldx * (i - 1) - 1] /= r[j + *ldr * (j - 1) - 1];
			temp1 = -x[j + *ldx * (i - 1) - 1];
			saxpy (j - 1, temp1, &r[*ldr * (j - 1)], 1,
			    &x[*ldx * (i - 1)], 1);
		    }
		}
	    }
	}
	else {
	    for (i = 1; i <= *nb; i++) {

		imsl_strsv ("UPPER", "NOT-TRANS", "NOT-DIAG", *n, r, *ldr, &x[*ldx * (i - 1)], 1);
	    }
	}

    }
    else if (*ipath == 2 || *ipath == 4) {
	/* Solve R'*X=B */
	if (*irank < *n) {
	    /* Case of singular R */
	    for (i = 1; i <= *nb; i++) {
		for (j = 1; j <= *n; j++) {
		    temp1 = x[j + *ldx * (i - 1) - 1] - imsl_sdot (j - 1,
			&r[*ldr * (j - 1)], 1, &x[*ldx * (i - 1)], 1);
		    if (r[j + *ldr * (j - 1) - 1] == 0.0) {
			_l0 = j - 1;
			_l1 = 1;
			_l2 = 1;
			temp2 = fabs (x[j + *ldx * (i - 1) - 1]) +
			    l_a1ot (&_l0, &r[*ldr * (j - 1)], &_l1,
			    &x[*ldx * (i - 1)], &_l2);
			temp2 *= 200.0 * imsl_amach (4);
			if (fabs (temp1) > temp2) {
			    imsl_e1sti (1, j-1);
			    imsl_ermes (IMSL_WARNING, IMSL_INCONSISTENT_EQUATIONS_2);
			}
			x[j + *ldx * (i - 1) - 1] = 0.0;
		    }
		    else {
			x[j + *ldx * (i - 1) - 1] = temp1 / r[j + *ldr * (j - 1) - 1];
		    }
		}
	    }
	}
	else {
	    /* Case of nonsingular R */
	    for (i = 1; i <= *nb; i++) {
		imsl_strsv ("UPPER", "TRANSPOSE", "NOT-UNIT", *n, r, *ldr, &x[*ldx * (i - 1)], 1);
	    }
	}
    }
    if (*ipath == 3 || *ipath == 4) {
	/* Invert R */
	for (j = 1; j <= *n; j++) {
	    scopy (j, &r[*ldr * (j - 1)], 1, &rinv[*ldrinv * (j - 1)],
		1);
	}
	for (k = 1; k <= *n; k++) {
	    if (rinv[k + *ldrinv * (k - 1) - 1] == 0.0) {
		sset (k, 0.0, &rinv[*ldrinv * (k - 1)], 1);
		sset (*n - k, 0.0, &rinv[k + *ldrinv * k - 1], *ldrinv);
	    }
	    else {
		rinv[k + *ldrinv * (k - 1) - 1] = 1.0 / rinv[k + *ldrinv * (k - 1) - 1];
		temp1 = -rinv[k + *ldrinv * (k - 1) - 1];
		sscal (k - 1, temp1, &rinv[*ldrinv * (k - 1)], 1);
		if (k < *n) {
		    imsl_sger (k - 1, *n - k, 1.0, &rinv[*ldrinv * (k - 1)],
			1, &rinv[k + *ldrinv * k - 1], *ldrinv, &rinv[*ldrinv * k],
			*ldrinv);
		    sscal (*n - k, rinv[k + *ldrinv * (k - 1) - 1],
			&rinv[k + *ldrinv * k - 1], *ldrinv);
		}
	    }
	}
	/*
	 * Fill lower triangle of RINV with zeros
	 */
	l_c1trg (n, rinv, ldrinv);
    }
    /* Exit section */
L_9000:
    imsl_e1pop ("l_girts");
    return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 08/24/90 at 16:32:47
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  C1R/DC1R (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 17, 1985

    Purpose:    Call E1MES if a diagonal element of the Cholesky factor R
                is zero but a remaining element in that row of R is not.

    Usage:      CALL C1R (N, R, LDR, NER)

    Arguments:
       N      - Order of R.  (Input)
       R      - Upper triagular matrix, N by N, with the Cholesky factor.
       LDR    - Leading dimension of R exactly as specified in the
                dimension statement of the calling program.  (Input)
       NER    - Error code used in call to E1MES.  (Input/Output)
                On output NER is incremented by one.

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c1r (Mint *n, Mfloat *r, Mint *ldr, Mint *ner)
#else
static void l_c1r (n, r, ldr, ner)
    Mint       *n;
    Mfloat      r[];
    Mint       *ldr, *ner;
#endif
{
    Mint        i, j;


    for (i = 1; i <= *n; i++) {
	if (r[i + *ldr * (i - 1) - 1] == 0.0) {
	    for (j = i + 1; j <= *n; j++) {
		if (r[i + *ldr * (j - 1) - 1] != 0.0) {
		    imsl_e1sti (1, i-1);
		    imsl_e1sti (2, j-1);
		    imsl_e1str (1, r[i + *ldr * (i - 1) - 1]);
		    imsl_e1str (2, r[i + *ldr * (j - 1) - 1]);
		    imsl_ermes (IMSL_TERMINAL,
			IMSL_ROW_ELMNTS_MUST_BE_ZERO);
		    goto L_9000;
		}
	    }
	}
    }
    *ner += 1;
L_9000:
    return;
}				/* end of function */
/*---------------------------------------------------------------------- */

/*  IMSL Name:  C1TRG/DC1TRG  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 7, 1985

    Purpose:    Fill out the lower triangular portion of an upper
                triangular matrix with zeros.

    Usage:      C1TRG (N, A, LDA)

    Arguments:
       N      - Order of the matrix A.  (Input)
       A      - N x N upper triangular matrix to be filled out.
                (Input/Output)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)

    GAMS:       D1b9;

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Utilities

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c1trg (Mint *n, Mfloat *a, Mint *lda)
#else
static void l_c1trg (n, a, lda)
    Mint       *n;
    Mfloat      a[];
    Mint       *lda;
#endif
{
    Mint        i;


    imsl_e1psh ("l_c1trg");

    if (*n <= 0) {
	imsl_e1sti (1, *n);

/*		imsl_ermes(5, 1, "N = %(i1).  The order of A, N, must be greater than 0.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_A_AND_N_GT_ZERO);
    }
    if (*n > *lda) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *lda);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_N_LE_LDA);
    }
    if (imsl_n1rcd (0) != 0)
	goto L_9000;
    /* Fill lower triangle with zeros */
    for (i = 1; i <= (*n - 1); i++) {
	sset (*n - i, 0.0, &a[(i + 1) + *lda * (i - 1) - 1], 1);
    }
    /* Exit section */
L_9000:
    imsl_e1pop ("l_c1trg");
    return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  A1OT/DA1OT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Compute sum of absolute values of products.

    Usage:      A1OT(N, SX, INCX, SY, INCY)

    Arguments:
       N      - Length of vectors X and Y.  (Input)
       SX     - Real vector of length MAX(N*IABS(INCX),1).  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be.. SX(1+(I-1)*INCX) if INCX .GE. 0
                or SX(1+(I-N)*INCX) if INCX .LT. 0.
       SY     - Real vector of length MAX(N*IABS(INCY),1).  (Input)
       INCY   - Displacement between elements of SY.  (Input)
                Y(I) is defined to be.. SY(1+(I-1)*INCY) if INCY .GE. 0
                or SY(1+(I-N)*INCY) if INCY .LT. 0.
       A1OT   - Sum from I=1 to N of ABS(X(I)*Y(I)).  (Output)
                X(I) and Y(I) refer to specific elements of SX and SY,
                respectively.  See INCX and INCY argument descriptions.

    GAMS:       D1a4

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_a1ot (Mint *n, Mfloat *sx, Mint *incx, Mfloat *sy, Mint *incy)
#else
static Mfloat l_a1ot (n, sx, incx, sy, incy)
    Mint       *n;
    Mfloat      sx[];
    Mint       *incx;
    Mfloat      sy[];
    Mint       *incy;
#endif
{
    Mint        i, ix, iy;
    Mfloat      a1ot_v;


    a1ot_v = 0.0;
    if (*n > 0) {
	if (*incx != 1 || *incy != 1) {
	    /* CODE FOR UNEQUAL INCREMENTS */
	    ix = 1;
	    iy = 1;
	    if (*incx < 0)
		ix = (-*n + 1) ** incx + 1;
	    if (*incy < 0)
		iy = (-*n + 1) ** incy + 1;
	    for (i = 1; i <= *n; i++) {
		a1ot_v += fabs (sx[ix - 1] * sy[iy - 1]);
		ix += *incx;
		iy += *incy;
	    }
	}
	else {
	    for (i = 1; i <= *n; i++) {
		a1ot_v += fabs (sx[i - 1] * sy[i - 1]);
	    }
	}
    }
    return (a1ot_v);
}				/* end of function */
