#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK 	l_eig_symgen(Mint n, Mfloat *a, Mfloat *b,
				 va_list argptr);
static void l_csfrg(Mint *n, Mfloat *a, Mint *lda);
static void l_e4lsf(Mint *n, Mfloat *a, Mint *lda, Mfloat eval[],
                        Mfloat *wk, Mint iwk[]);
static void l_e5csf(Mint *n, Mfloat *a, Mint *lda, Mfloat eval[],
                        Mfloat *evec, Mint *ldevec, Mfloat *wk,
                        Mint iwk[]);
static void l_e5lsf(Mint *n, Mfloat d[], Mfloat e2[], Mint iwk[]);
static void l_e6csf(Mint *n, Mfloat *a, Mint *lda, Mfloat eval[],
                        Mfloat e[], Mfloat e2[], Mfloat *evec,
                        Mint *ldevec, Mint *eigvec, Mfloat *scale);
static void l_e7csf(Mint *lb, Mint *nb, Mint *n, Mfloat eval[],
                        Mfloat dwk[], Mfloat e[], Mfloat *fnorm,
                        Mfloat *evec, Mint *ldevec, Mfloat wk[]);
static void l_g3csp(Mint *n, Mfloat *a, Mint *lda, Mfloat *b, 
                        Mint *ldb, Mfloat eval[], Mfloat *evec,
                        Mint *ldevec, Mint iwk[], Mfloat *wk,
                        Mfloat *s);
static void l_g3lsp(Mint *n, Mfloat *a, Mint *lda, Mfloat *b,
                        Mint *ldb, Mfloat eval[], Mint iwk[],
                        Mfloat *wk, Mfloat *s);
static void l_lftds(Mint *n, Mfloat *a, Mint *lda, Mfloat *imsl_fac,
                        Mint *ldfac);
static void l_strsm(Mchar *side, unsigned side_s, Mchar *uplo,
                        unsigned uplo_s, Mchar *transa, unsigned
                        transa_s, Mchar *diag, unsigned diag_s,
                        Mint *m, Mint *n, Mfloat *alpha, Mfloat *a,
                        Mint *lda, Mfloat *b, Mint *ldb);
#else
static VA_LIST_HACK	l_eig_symgen();
static void 	l_csfrg();
static void 	l_e4lsf();
static void 	l_e5csf();
static void 	l_e5lsf();
static void	l_e6csf();
static void	l_e7csf();
static void	l_g3csp();
static void	l_g3lsp();
static void	l_lftds();
static void	l_strsm();
#endif

static Mfloat	*lv_eval;

#ifdef ANSI
Mfloat *imsl_f_eig_symgen(Mint n, Mfloat *a, Mfloat *b, ...)
#else
Mfloat *imsl_f_eig_symgen(n, a, b, va_alist)
    Mint	n;
    Mfloat	*a;
    Mfloat	*b;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr,b);

#ifdef DOUBLE
    imsl_e1psh("imsl_d_eig_symgen");
#else
    imsl_e1psh("imsl_f_eig_symgen");
#endif
    lv_eval = NULL;
    IMSL_CALL(l_eig_symgen(n, a, b, argptr));
    va_end(argptr);
#ifdef DOUBLE
    imsl_e1pop("imsl_d_eig_symgen");
#else
    imsl_e1pop("imsl_f_eig_symgen");
#endif
    return lv_eval;
}


#ifdef ANSI
static VA_LIST_HACK l_eig_symgen(Mint n, Mfloat *a, Mfloat *b, va_list argptr)
#else
static VA_LIST_HACK l_eig_symgen(n, a, b, argptr)
    Mint	n;
    Mfloat	*a;
    Mfloat	*b;
    va_list	argptr;
#endif
{
    Mint	    code;
    Mint	    arg_number  = 2;
    Mint	    a_col_dim   = n;
    Mint	    b_col_dim   = n;
    Mfloat	    **evec	= NULL;
    Mfloat	    *evecu	= NULL;
    Mint	    evecu_col_dim = n;
    Mint	    vectors	= 0;
    Mint	    vectors_user = 0;
    Mint	    evals_user	= 0;
    Mint	    *iwork	= NULL;
    Mfloat	    *acopy	= NULL;
    Mfloat	    *work	= NULL;
    Mfloat	    *s		= NULL;

    code = 1;
    while (code > 0) {
	code = va_arg(argptr, Mint);
	arg_number++;
	switch (code) {
	    case IMSL_VECTORS:
		evec = va_arg(argptr, Mfloat**);
		arg_number++;
		vectors = 1;
		break;
	    case IMSL_VECTORS_USER:
		evecu = va_arg(argptr, Mfloat*);
		arg_number++;
		vectors_user = 1;
		break;
	    case IMSL_RETURN_USER:
		lv_eval = va_arg(argptr, Mfloat*);
		arg_number++;
		evals_user = 1;
		break;
	    case IMSL_A_COL_DIM:
		a_col_dim = va_arg(argptr, Mint);
		arg_number++;
		break;
	    case IMSL_B_COL_DIM:
		b_col_dim = va_arg(argptr, Mint);
		arg_number++;
		break;
	    case IMSL_EVECU_COL_DIM:
		evecu_col_dim = va_arg(argptr, Mint);
		arg_number++;
		break;
	    case 0:
		break;
	    default:
		imsl_e1sti (1, code);
		imsl_e1sti (2, arg_number);
		imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
		break;
	}
    } 

    if (imsl_n1rty(0)) goto RETURN;

    if (n <= 0) {
	imsl_e1sti(1, n);
	imsl_ermes(IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
    } else if (n > a_col_dim) {
	imsl_e1sti(1, n);
	imsl_e1sti(2, a_col_dim);
	imsl_e1stl(1, "a");
	imsl_ermes(IMSL_TERMINAL, IMSL_COL_DIM_LESS_ORDER);
    } else if (n > b_col_dim) {
	imsl_e1sti(1, n);
	imsl_e1sti(2, b_col_dim);
	imsl_e1stl(1, "b");
	imsl_ermes(IMSL_TERMINAL, IMSL_COL_DIM_LESS_ORDER);
    }
    if (imsl_n1rty(0)) goto RETURN;

    if (!vectors && !vectors_user) {
	work    = (Mfloat *) imsl_malloc(2*n*sizeof(*work));
	iwork	= (Mint *) imsl_malloc(n*sizeof(*iwork));
	s	= (Mfloat *) imsl_malloc((n*n+n)*sizeof(*s));
    }
    if (vectors || vectors_user) {
	work    = (Mfloat *) imsl_malloc(3*n*sizeof(*work));
	iwork	= (Mint *) imsl_malloc(n*sizeof(*iwork));
	s	= (Mfloat *) imsl_malloc((n*n+n)*sizeof(*s));
    }
    if (vectors)
	*evec	= (Mfloat *) imsl_malloc(n*n*sizeof(**evec));
    if ( work==NULL || iwork==NULL || s==NULL ||
		(vectors && (*evec==NULL)) ) {
	imsl_e1stl(1, "n");
	imsl_e1sti(1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto FREE_SPACE;
    }
    if (lv_eval == NULL) {
	lv_eval = (Mfloat *)imsl_malloc (n*sizeof(*lv_eval));
	if (lv_eval == NULL) {
	    imsl_e1stl(1, "n");
	    imsl_e1sti(1, n);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE;
	}
    }

    if (!vectors && !vectors_user)
	l_g3lsp(&n, a, &a_col_dim, b, &b_col_dim, lv_eval, iwork, work, s);
    if (vectors) { 
        l_g3csp(&n, a, &a_col_dim, b, &b_col_dim, lv_eval, *evec, &n,
			 iwork, work, s);
	imsl_trnrr(n, n, *evec, n, n, n, *evec, n);
	}
    if (vectors_user) {
        l_g3csp(&n, a, &a_col_dim, b, &b_col_dim, lv_eval, evecu, 
			&evecu_col_dim, iwork, work, s);
	imsl_trnrr(n, n, evecu, evecu_col_dim, n, n, evecu,
			 evecu_col_dim);
	}

FREE_SPACE:
    if (work != NULL) imsl_free(work);
    if (acopy != NULL) imsl_free(acopy);
    if (iwork != NULL) imsl_free(iwork);
    if (s != NULL) imsl_free(s);
RETURN:
    return (argptr);
}




/* Structured by FOR_STRUCT, v0.2, on 09/24/90 at 11:23:59
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CSFRG/DCSFRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 23, 1985

    Purpose:    Extend a real symmetric matrix defined in its upper
                triangle to its lower triangle.

    Usage:      CALL CSFRG (N, A, LDA)

    Arguments:
       N      - Order of the matrix A.  (Input)
       A      - N by N symmetric matrix of order N to be filled out.
                (Input/Output)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)

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
static void l_csfrg(Mint *n, Mfloat *a, Mint *lda)
#else
static void l_csfrg(n, a, lda)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda;
#endif
{
#define A(I_,J_)	(a+(I_)*(l_lda)+(J_))
	Mint             i;
	Mint             l_lda = *lda;


	imsl_e1psh("CSFRG ");

	if (*n <= 0) {
		imsl_e1sti(1, *n);
/*		imsl_ermes(5, 1, "N = %(i1).  The order of A, N, must be greater than 0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_A_AND_N_GT_ZERO);
	}
	if (*n > *lda) {
		imsl_e1sti(1, *n);
		imsl_e1sti(2, *lda);

/*		imsl_ermes(5, 2, "N = %(i1) and LDA = %(i2).  The order of A, N, must be less than or equal to the leading dimension of A, LDA.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_LE_LDA);
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;
	/*
	 * Copy upper triangular values to lower triangular values
	 */
	for (i = 1; i <= (*n - 1); i++) {
		scopy(*n - i, A(i, i - 1), *lda, A(i - 1, i), 1);
	}
	/* Exit section */
L_9000:
	imsl_e1pop("CSFRG ");
	return;
}				/* end of function */


#undef A

/* Structured by FOR_STRUCT, v0.2, on 09/24/90 at 13:11:26
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E4LSF/DE4LSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Purpose:    Compute all of the eigenvalues of a real symmetric
                matrix.

    Usage:      CALL E4LSF (N, A, LDA, EVAL, WK, IWK)

    Arguments:  (See EVLSF)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e4lsf(Mint *n, Mfloat *a, Mint *lda, Mfloat eval[],
			Mfloat *wk, Mint iwk[])
#else
static void l_e4lsf(n, a, lda, eval, wk, iwk)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda;
	Mfloat           eval[], *wk;
	Mint             iwk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(l_lda)+(J_))
#define WK(I_,J_)	(wk+(I_)*(l_n)+(J_))
	Mint             l_lda = *lda;
	Mint             l_n = *n;
	Mint             _l0, i;
	Mfloat           scale;


	imsl_e1psh("E4LSF ");
	/* Check N */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "The argument N = %(i1).  The order of the matrix must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Check LDA */
	if (*lda < *n) {
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 2, "The argument LDA = %(i1).  The leading dimension of the matrix must be at least equal to the order of the matrix, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/*
	 * The lower triangle is used. The upper part of the matrix is
	 * preserved. A copy of 'A' is not needed. Reduce to tridiagonal Do
	 * not accumulate transformations
	 */
	_l0 = 0;
	l_e6csf(n, a, lda, eval, WK(0, 0), WK(1, 0), a, lda, &_l0, &scale);
	/* Matrix is zero; exit */
	if (scale == F_ZERO)
		goto L_9000;
	/* Find eigenvalues */
	l_e5lsf(n, eval, WK(1, 0), iwk);
	/* Scale back the eigenvalues. */
	if (scale != F_ONE)
		sscal(*n, scale, eval, 1);
	/*
	 * Sort eigenvalues into ascending order of magnitude.
	 */
	imsl_svrbn(n, eval, WK(1, 0));
	/*
	 * Resort into descending order of magnitude
	 */
	for (i = 1; i <= *n; i++) {
		eval[i - 1] = *WK(1, *n - i);
	}
	/*
	 * Retrieve the original symmetric matrix. Copy the upper triangle to
	 * lower so that a is restored.
	 */
	for (i = 1; i <= (*n - 1); i++) {
		scopy(*n - i, A(i, i - 1), *lda, A(i - 1, i), 1);
	}

L_9000:
	imsl_e1pop("E4LSF ");
	return;
}				/* end of function */


#undef A
#undef WK



/* Structured by FOR_STRUCT, v0.2, on 09/24/90 at 11:29:45
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E5CSF/DE5CSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Purpose:    Compute all of the eigenvalues and eigenvectors of a real
                symmetric matrix.

    Usage:      CALL E5CSF (N, A, LDA, EVAL, EVEC, LDEVEC, WK,IWK)

    Arguments:  (See EVCSF)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  ----------------------------------------------------------------------- */
#ifdef ANSI
static void l_e5csf(Mint *n, Mfloat *a, Mint *lda, Mfloat eval[],
			Mfloat *evec, Mint *ldevec, Mfloat *wk,
			Mint iwk[])
#else
static void l_e5csf(n, a, lda, eval, evec, ldevec, wk, iwk)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda;
	Mfloat           eval[], *evec;
	Mint            *ldevec;
	Mfloat          *wk;
	Mint             iwk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(l_lda)+(J_))
#define EVEC(I_,J_)	(evec+(I_)*(l_ldevec)+(J_))
#define WK(I_,J_)	(wk+(I_)*(l_n)+(J_))
	Mint             l_lda = *lda;
	Mint             l_ldevec = *ldevec;
	Mint             l_n = *n;
	Mint       first;
	Mint             _l0, i, j, k, lb, nb;
	Mfloat           fnorm, h, scale, test;


	imsl_e1psh("l_e5csf ");
	/*
	 * wk(1,1) is used in imsl_e5lsf and imsl_e7csf wk(1,2) and wk(1,3)
	 * are the subdiagonal and the squares of the subdiagonal elements
	 * respectively Check N
	 */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "The argument N = %(i1).  The order of the matrix must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);
		goto L_9000;
	}
	/* Check LDA */
	if (*lda < *n) {
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 2, "The argument LDA = %(i1).  The leading dimension of the matrix must be at least equal to the order of the matrix, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);
	}
	/* Check LDEVEC */
	if (*ldevec < *n) {
		imsl_e1sti(1, *ldevec);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 3, "The argument LDEVEC = %(i1).  The leading dimension of the eigenvector matrix must be at least the order of the matrix, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDEVEC_VALUE_TOO_SMALL);
	}
	if (imsl_n1rty(0) > 0)
		goto L_9000;

	/* Reduction is trivial. */
	if (*n == 1) {
		*EVEC(0, 0) = F_ONE;
		eval[0] = *A(0, 0);
		goto L_9000;
	}
	/*
	 * The logical variable first is is used to prevent us from forming
	 * the square of the off-diagonal elements since they are returned in
	 * wk(1,3) from imsl_e6csf.
	 * 
	 * The lower half of 'a' is used. The upper half will remain intact.
	 * Thus, a copy of 'A' is not needed.
	 */
/*	first = TRUE; */
	first = 1;
	/*
	 * Reduce matrix to tridiagonal form and accumulate the
	 * transformation returned in evec(*,*).
	 */
	_l0 = 1;
	l_e6csf(n, a, lda, eval, WK(1, 0), WK(2, 0), evec, ldevec, &_l0,
		   &scale);
	/*
	 * If scale is equal to zero, then evec is the identity matrix and
	 * eval is set to the zero vector. Exit.
	 */
	if (scale == F_ZERO)
		goto L_9000;
	/*
	 * Form the Frobenius norm of the matrix used for a scale factor
	 */
	fnorm = F_ZERO;
	for (i = 1; i <= *n; i++) {
		fnorm += imsl_fi_power(eval[i - 1], 2) + *WK(2, i - 1);
	}
	fnorm = sqrt(fnorm * scale);
	/*
	 * Test = 1.0e0/imsl_amach(2) is the bound on smallest number that
	 * can be reciprocated.
	 */
	test = imsl_amach(1);
	if (test * imsl_amach(2) < F_ONE)
		test = F_ONE / imsl_amach(2);
	/*
	 * lb and nb are the start and end of current block respectively.
	 */
	lb = 1;
	nb = *n;
	/* Find the start of the block. */
L_20:
	for (i = lb; i <= (nb - 1); i++) {
		h = imsl_amach(4) * (fabs(eval[i - 1]) + fabs(*WK(1, i)));
		if (test < h)
			test = h;
		if (fabs(*WK(1, i)) > test)
			goto L_40;
	}
	/* All eigenvalues and vectors found. */
	goto L_80;
	/* lb is the start of the block. */
L_40:
	lb = i;
	/*
	 * Find the end of the block and store into nb
	 */
	for (i = lb + 1; i <= (nb - 1); i++) {
		if (fabs(*WK(1, i)) <= test) {
			nb = i;
			goto L_60;
		}
	}
	/*
	 * Make a copy of the diagonal elements of the tridiagonal matrix
	 * eval
	 */
L_60:
	scopy(nb - lb + 1, &eval[lb - 1], 1, WK(0, lb - 1), 1);
	/*
	 * Square off-diagonal elements if not the first time.
	 */
	if (!first) {
		for (k = 1; k <= (nb - lb); k++) {
			/*
			 * lb + k is used instead of lb because wk(lb,3) is a
			 * dummy element
			 */
			*WK(2, lb + k - 1) = imsl_fi_power(*WK(1, lb + k - 1), 2);
		}
	}
	/*
	 * Find the eigenvalues of the tridiagonal block.
	 */
	_l0 = nb - lb + 1;
	l_e5lsf(&_l0, WK(0, lb - 1), WK(2, lb - 1), iwk);
	/*
	 * Sort eigenvalues into ascending order of magnitude.
	 */
	_l0 = nb - lb + 1;
	imsl_svrbn(&_l0, WK(0, lb - 1), WK(0, lb - 1));
	/*
	 * Find the eigenvectors of the tri- diagonal block and update evec.
	 * Employing the explicit perfect shift, accumulate the remaining
	 * similarity transforms to find the eigenvectors.
	 */
	l_e7csf(&lb, &nb, n, eval, WK(0, 0), WK(1, 0), &fnorm, evec, ldevec,
		   WK(2, 0));
	/*
	 * If lb less than n implies the matrix has split or that a Wilkinson
	 * shift was used.
	 */
	if (lb < nb) {
		/*
		 * This is to ensure that the off-diagon elements are squared
		 * for use in imsl_e5lsf
		 */
/*		first = FALSE; */
		first = 0;
		goto L_20;
		/*
		 * lb equal to nb implies that we are done with the current
		 * block.
		 */
	} else if (lb == nb) {
		lb = nb + 1;
		nb = *n;
		goto L_20;
	}
	/* Scale back the eigenvalues. */
L_80:
	if (scale != F_ONE)
		sscal(*n, scale, eval, 1);
	/*
	 * Final sort of eigenvalues and eigenvectors.
	 */
	for (i = 1; i <= *n; i++) {
		iwk[i - 1] = i;
		*WK(2, i - 1) = -fabs(eval[i - 1]);
	}
	/*
	 * Eigenvalues are sorted into ascending order by magnitude
	 */
	imsl_svrgp(*n, WK(2, 0), WK(2, 0), iwk);
	/* Resort the record of the pivots. */
	for (i = 1; i <= *n; i++) {
		for (j = i; j <= *n; j++) {
			if (iwk[j - 1] == i) {
				k = iwk[i - 1];
				iwk[i - 1] = j;
				iwk[j - 1] = k;
				goto L_110;
			}
		}
L_110:
		;
	}
	/*
	 * Finally leave the eigenvalues and eigenvectors sorted in
	 * descending order of magnitude.
	 */
	for (i = *n - 1; i >= 1; i--) {
		if (i == iwk[i - 1])
			goto L_120;
		sswap(*n, EVEC(i - 1, 0), 1, EVEC(iwk[i - 1] - 1, 0), 1);
		sswap(1, &eval[i - 1], 1, &eval[iwk[i - 1] - 1], 1);
L_120:
		;
	}
	/*
	 * Normalize each eigenvector so that its biggest component is
	 * positive. The eigenvectors then form a right-hand system.
	 */
	for (j = 1; j <= *n; j++) {
		i = imsl_isamax(*n, EVEC(j - 1, 0), 1);
		if (*EVEC(j - 1, i - 1) < F_ZERO) {
			for (k = 1; k <= *n; k++) {
				*EVEC(j - 1, k - 1) = -*EVEC(j - 1, k - 1);
			}
		}
	}
	/*
	 * Retrieve the original symmetric matrix. Copy the upper triangle to
	 * lower.
	 */
	for (i = 1; i <= (*n - 1); i++) {
		scopy(*n - i, A(i, i - 1), *lda, A(i - 1, i), 1);
	}

L_9000:
	imsl_e1pop("l_e5csf ");
	return;
}				/* end of function */


#undef A
#undef EVEC
#undef WK


/* Structured by FOR_STRUCT, v0.2, on 09/24/90 at 11:34:54
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E5LSF/DE5LSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Purpose:    Compute eigenvalues of a symmetric triadiagonal matrix.

    Usage:      CALL E5LSF (N, D, E2,  IWK)

    Arguments:
       N      - Order of the matrix.  (Input)
       D      - Real vector of length N.  On input, the diagonal elements
                of the matrix. On output, the eigenvalues of the
                matrix.  (Input/Output)
       E2     - Real vector of length N.  On input, the squares of
                the elements of the off diagonal.  E2(1) is arbitrary.
                On output, E2 is destroyed.  (Input/Output)
       IWK    - Integer work array.  (Input)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  ----------------------------------------------------------------------- */
#ifdef ANSI
static void l_e5lsf(Mint *n, Mfloat d[], Mfloat e2[], Mint iwk[])
#else
static void l_e5lsf(n, d, e2, iwk)
	Mint            *n;
	Mfloat           d[], e2[];
	Mint             iwk[];
#endif
{
	Mint             i, iter, j, k, l, m;
	Mfloat           cossq, e, epslon, gama, h, p, prevc, prevga, recipe,
	                sigma, sinsq;


	imsl_e1psh("l_e5lsf ");

	if (*n == 1)
		goto L_9000;
	/*
	 * The following code is based on the PWK algorithm, pages 164-169,
	 * "Symmetric Eigenvalue Problem" B. Parlett.
	 */
	epslon = imsl_fi_power(imsl_amach(4), 2);
	e2[0] = F_ZERO;
	iter = 0;
	/*
	 * Bound on smallest number that can be reciprocated.
	 */
	h = imsl_amach(1);
	if (h * imsl_amach(2) < F_ONE)
		h = F_ONE / imsl_amach(2);
	/*
	 * iwk() is used for the tracking when the current block decouples j
	 * is the pointer.
	 */
	j = 1;
	iwk[j - 1] = 1;
	/*
	 * k is the start of the current block and m is the end of the
	 * current block.
	 */
	k = 1;
	m = *n;
	/*
	 * Deflate the current block i times if possible.
	 */
L_10:
	e = epslon * (F_HALF * e2[m - 1] + d[m - 1] * d[m - 1]);
	if (e > h)
		h = e;
	l = m;
	for (i = l; i >= (k + 1); i--) {
		if (e2[i - 1] > h)
			goto L_30;
		m -= 1;
		iter = 0;
	}
	/*
	 * Check to see if current block can be decoupled further due to
	 * appearance of new zeros at the top of the current block.
	 */
L_30:
	for (i = m - 2; i >= k; i--) {
		if (e2[i] > h)
			goto L_40;
		j += 1;
		iwk[j - 1] = i + 1;
		goto L_50;
L_40:
		;
	}
	/*
	 * Make sure, in the current block, that k < m.
	 */
L_50:
	l = j;
	for (i = l; i >= 1; i--) {
		if (iwk[i - 1] < m)
			goto L_70;
		j -= 1;
		m -= 1;
	}
	/* Exit for m <= 1. */
L_70:
	if (m <= 1)
		goto L_9000;
	k = iwk[j - 1];
	iter += 1;
	/*
	 * Initialize the implicit shift sigma The Wilkinson shift is used.
	 * This is the eigenvalue of the SE 2 by 2 nearest the SE entry.
	 */
	e = (d[m - 2] - d[m - 1]) * F_HALF;
	sigma = d[m - 1] - e2[m - 1] / sign(sqrt(imsl_fi_power(e, 2) + e2[m - 1]) +
					    fabs(e), e);
	/* Initialize variables for the QR */
	cossq = F_ONE;
	sinsq = F_ZERO;
	gama = d[k - 1] - sigma;
	p = gama * gama;
	/*
	 * The next 13 lines are executed outside the main loop to avoid a
	 * conditional statement. Start the QR factorization. Pre and Post
	 * multiply the current block by the first plane rotation that
	 * attempts to eliminate the off diagonal element b(k) of the block.
	 */
	e = p + e2[k];
	recipe = F_ONE / e;
	prevc = cossq;
	cossq = p * recipe;
	sinsq = e2[k] * recipe;
	prevga = gama;
	gama = cossq * (d[k] - sigma) - sinsq * prevga;
	d[k - 1] = prevga + (d[k] - gama);
	p = e2[k];
	if (cossq != F_ZERO)
		p = imsl_fi_power(gama, 2) / cossq;
	/*
	 * Finish the QR factorization of the current block. The following
	 * loop "chases the bulge".
	 */
	for (i = k + 1; i <= (m - 1); i++) {
		e = p + e2[i];
		recipe = F_ONE / e;
		e2[i - 1] = sinsq * e;
		prevc = cossq;
		/*
		 * cossq and sinsq repesent the i-th plane rotation, where
		 * cossq + sinsq = 1.
		 */
		cossq = p * recipe;
		sinsq = e2[i] * recipe;
		prevga = gama;
		gama = cossq * (d[i] - sigma) - sinsq * prevga;
		/*
		 * update i-th diagonal element that represents the shifted
		 * i-th eigenvalue.
		 */
		d[i - 1] = prevga + (d[i] - gama);
		if (cossq == F_ZERO) {
			p = prevc * e2[i];
		} else {
			p = gama * gama / cossq;
		}
	}
	/*
	 * Update the last off-diagonal element of the current block.
	 */
	e2[m - 1] = sinsq * p;
	/*
	 * Add shift back in, a(m) is the current m-th eigenvalue.
	 */
	d[m - 1] = gama + sigma;
	/*
	 * Iterate back to the top The iteration counter(.le. 100) is
	 * arbitrary. Deflation is done if 100 iterations were taken. Usually
	 * only 2 - 6 iterations are required.
	 */
	if (iter <= 100) {
		goto L_10;
	} else {

/*		imsl_ermes(3, 1, "The iteration for the eigenvalue failed  to converge in 100 iterations before               deflating ");
*/
                imsl_ermes(IMSL_WARNING, IMSL_SLOW_CONVERGENCE_SYM);
	}
	/*
	 * Use current value of a(m) as the eigenvalue. Deflate problem and
	 * iterate back to the top.
	 */
	m -= 1;
	goto L_10;

L_9000:
	imsl_e1pop("l_e5lsf ");

	return;
}				/* end of function */



/* Structured by FOR_STRUCT, v0.2, on 09/24/90 at 11:37:11
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E6CSF/DE6CSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Purpose:    Reduce a real symmetric matrix to a tridiagonal
                matrix by orthogonal similarity transformations.

    Usage:      CALL E6CSF (N, A, LDA, EVAL, E, E2, EVEC, EIGVEC, SCALE)

    Arguments:
       N      - Order of the matrix.  (Input)
       A      - Real symmetric matrix (Only lower half used).
                (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       EVAL   - Real vector of length N.  (Input/Output)
       E      - Real vector of length N containing the subdiagonal
                elements of the tridiagonal matrix (E(1) = 0).  (Output)
       E2     - Real vector of length N containing the squares of the
                subdiagonal elements (E2(1) = 0)).  (Output)
       EVEC   - Accumulated product of transformations. (Output)

       EIGVEC - Logical variable (Input)

       SCALE  - Real variable used for scaling.  (Output)

    Remark:
       THIS SUBROUTINE IS BASED ON THE ALGOL PROCEDURE TRED1,
       NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
       HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).

       This code is based on code developed by A. Dubrulle for the
       IBM 3090-VF.

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

       ------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e6csf(Mint *n, Mfloat *a, Mint *lda, Mfloat eval[],
			Mfloat e[], Mfloat e2[], Mfloat *evec,
			Mint *ldevec, Mint *eigvec, Mfloat *scale)
#else
static void l_e6csf(n, a, lda, eval, e, e2, evec, ldevec, eigvec, scale)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda;
	Mfloat           eval[], e[], e2[], *evec;
	Mint		*ldevec;
	Mint      *eigvec;
	Mfloat          *scale;
#endif
{
#define A(I_,J_)	(a+(I_)*(l_lda)+(J_))
#define EVEC(I_,J_)	(evec+(I_)*(l_ldevec)+(J_))
        Mint            l_lda = *lda;
        Mint            l_ldevec = *ldevec;
	Mint            i, j, k;
	Mfloat           c, p, s, t, tol;


	imsl_e1psh("l_e6csf ");

	/*
	 * ------------------------------------ The reduction of a symmetric
	 * matrix A to symmetric tridiagonal is performed by the recurrence k
	 * k-1 A = H  *  A  *  H   , k=1,..,n-2 k+1         k+1
	 * 
	 * 0        n-2 where  A  = A.  A   is then the required symmetric
	 * tridiagonal matrix unitarily similar to A. H   is the Householder
	 * matrix k+1 associated with the k-th column of k
	 * t A,  H  =  I + u(1) * u  *  u . k+1              k+1   k+1
	 * ------------------------------------
	 * 
	 * Tol = imsl_amach(1)/imsl_amach(4) is used to ensure that a
	 * transformation is skipped whenever there is a danger that it might
	 * not be accurately orthogonal.
	 */
	tol = imsl_amach(1) / imsl_amach(4);
	/*
	 * Preserve the diagonal elements and find the element of maximum
	 * magnitude
	 */
	*scale = F_ZERO;
	for (j = 1; j <= *n; j++) {
		eval[j - 1] = *A(j - 1, j - 1);
		for (i = j; i <= *n; i++) {
			if (fabs(*A(j - 1, i - 1)) > *scale)
				*scale = fabs(*A(j - 1, i - 1));
		}
	}
	/* Matrix is all zero. */
	if (*scale == F_ZERO) {
		/* Set the vector eval to zero */
		sset(*n, F_ZERO, eval, 1);
		/* Set evec to the identity matrix */
		if (*eigvec) {
			for (i = 1; i <= *n; i++) {
				sset(*n, F_ZERO, EVEC(i - 1, 0), 1);
				*EVEC(i - 1, i - 1) = F_ONE;
			}
		}
		/* Exit */
		goto L_220;
	}
	/* Reduction is trival */
	if (*n == 1) {
		*scale = F_ONE;
		if (*eigvec)
			*EVEC(0, 0) = F_ONE;
		goto L_220;
	}
	/*
	 * Form the ratio imsl_alog10(scale) to imsl_alog(base of machine)
	 */
/*	k = trunc(log10(*scale) / imsl_amach(5)); */
	k = (Mint)(log10(*scale) / imsl_amach(5));
	/* imsl_imach(3) is the base of the machine */
	s = F_ONE / (Mfloat) (imsl_i_machine(1));
	/* Adjust for scale < 1.0 */
	if (*scale < F_ONE) {
		s = F_ONE / s;
		k = -k;
	}
	/*
	 * Compute scale. This product is typically exact.
	 */
	*scale = F_ONE;
	for (i = 1; i <= k; i++) {
		*scale *= s;
	}
	/*
	 * Scale the matrix (lower triangle) This product is typically exact.
	 */
	if (*scale != F_ONE) {
		for (j = 1; j <= *n; j++) {
			for (i = j; i <= *n; i++) {
				*A(j - 1, i - 1) *= *scale;
			}
		}
	}
	for (k = 1; k <= (*n - 2); k++) {
		/*
		 * swap of the original and transformed diagonal elements
		 */
		s = eval[k - 1];
		eval[k - 1] = *A(k - 1, k - 1);
		*A(k - 1, k - 1) = s;
		/*
		 * computation of the square of the elements of a(k+2:n,k)
		 */
		e2[k] = F_ZERO;
		e[k] = *A(k - 1, k);
		for (i = k + 2; i <= *n; i++) {
			e[i - 1] = *A(k - 1, i - 1);
			e2[k] += e[i - 1] * e[i - 1];
		}
		/*
		 * Skip the k-th transformation if the sum of squares is too
		 * small for orthogonality to be guaranteed.
		 */
		if ((e2[k] + imsl_fi_power(e[k], 2)) <= tol) {
			e[k] = F_ZERO;
			*A(k - 1, k) = F_ZERO;
			e2[k] = F_ZERO;
			/*
			 * If e2(k+1) is zero then column k is in tridiagonal
			 * form, no reduction necesssary, skip to column k+1.
			 */
		} else if (e2[k] <= tol) {
			e[k] = *A(k - 1, k);
			*A(k - 1, k) = F_ZERO;
			e2[k] = e[k] * e[k];
			/*
			 * Reduce column k to upper hessenberg form using
			 * Householder reflectors.
			 */
		} else {
			/*
			 * compute the l2 norm of a(k+1:n, k) and the
			 * codiagonal element.
			 */
			e2[k] += imsl_fi_power(e[k], 2);
			p = -sign(sqrt(e2[k]), e[k]);
			t = F_ONE / p;
			/*
			 * form u   ; store in a(k+1:n, k) k+1 Note that
			 * t*e(k+1) <= 0.0.
			 */
			*A(k - 1, k) = t * e[k] - F_ONE;
			for (i = k + 2; i <= *n; i++) {
				*A(k - 1, i - 1) = t * e[i - 1];
			}
			/* Normalizing factor of the reflection */
			t = F_ONE / *A(k - 1, k);
			/*
			 * ----------------------------------- The active
			 * submatrix is now updated. k An equivalent form for
			 * A can be employed: k    k-1           t          t
			 * A  = A    +   u  * f   +  f * u k+1
			 * k+1
			 * 
			 * k-1 f = imsl_gamma*( A * u   +  sigma*u ) k+1
			 * k+1
			 * 
			 * t    k-1 sigma = 0.5*u(1)*( u  * A  * u ) k+1
			 * k+1 ----------------------------------- k-1
			 * Compute A * u ; store in e(k+1,...,) k+1
			 */
			e[*n - 1] = *A(k - 1, *n - 1) ** A(*n - 1, *n - 1);
			for (j = *n - 1; j >= (k + 1); j--) {
				s = *A(k - 1, j - 1) ** A(j - 1, j - 1);
				c = *A(k - 1, j - 1);
				for (i = *n; i >= (j + 1); i--) {
					s += *A(k - 1, i - 1) ** A(j - 1, i - 1);
					e[i - 1] += c ** A(j - 1, i - 1);
				}
				e[j - 1] = s;
			}
			/*
			 * t    k-1 Compute  u  * A  * u k+1       k+1
			 */
			s = F_ZERO;
			for (i = k + 1; i <= *n; i++) {
				s += e[i - 1] ** A(k - 1, i - 1);
			}
			/* Compute sigma */
			s *= F_HALF * t;
			/* Compute f */
			for (i = k + 1; i <= *n; i++) {
				e[i - 1] = t * (e[i - 1] + s ** A(k - 1, i - 1));
			}
			/*
			 * k Compute A ;transformation of the unreduced
			 * submatrix .
			 */
			for (j = k + 1; j <= *n; j++) {
				/*
				 * ******* IBM uses directive ( IGNORE
				 * RECDEPS ******
				 */
				s = e[j - 1];
				c = *A(k - 1, j - 1);
				for (i = j; i <= *n; i++) {
					*A(j - 1, i - 1) += s ** A(k - 1, i - 1) + c * e[i - 1];
				}
			}
			e[k] = p;
		}
	}
	/*
	 * Take care of the last elements of the matrix.
	 */
	s = eval[*n - 2];
	eval[*n - 2] = *A(*n - 2, *n - 2);
	*A(*n - 2, *n - 2) = s;
	e2[*n - 1] = imsl_fi_power(*A(*n - 2, *n - 1), 2);
	e[*n - 1] = *A(*n - 2, *n - 1);
	s = eval[*n - 1];
	eval[*n - 1] = *A(*n - 1, *n - 1);
	*A(*n - 1, *n - 1) = s;
	/* Reciprocate scale for output */
	*scale = F_ONE / *scale;
	/*
	 * Construct right operator of the similarity transformation by
	 * accumulating the Householder matricies.
	 */
	if (*eigvec) {
		/* Set the matrix Z to zero */
		for (i = 1; i <= *n; i++) {
			sset(*n, F_ZERO, EVEC(i - 1, 0), 1);
		}
		if (*n > 1)
			*EVEC(*n - 2, *n - 2) = F_ONE;
		*EVEC(*n - 1, *n - 1) = F_ONE;
		/*--------------------------------------------------------------------
		                                    Form the product H * H * .. * H,
		                                                      1   2        n-2
		
		                                    where H   is the Householder matrix
		                                           k+1
		                                    associated with the k-th column of
		                                     k                          t
		                                    A,  H  =  I + u(1) * u  *  u .
		                                         k+1              k+1   k+1
		
		                                    Store into the matrix Z.
		                                           n-2     t
		                                    Then  A    =  Z * A * Z.
		  -------------------------------------------------------------------- */
		for (k = *n - 2; k >= 1; k--) {
			*EVEC(k - 1, k - 1) = F_ONE;
			/*
			 * u(1) associated with u   is the k+1 reciprocal of
			 * a(k+1,k).
			 */
			if (*A(k - 1, k) != F_ZERO) {
				t = F_ONE / *A(k - 1, k);
				for (j = k + 1; j <= *n; j++) {
					s = F_ZERO;
					/*
					 * Form the dot product of a(k+1:n,k)
					 * with evec(k+1:n,j)
					 */
					for (i = k + 1; i <= *n; i++) {
						s += *A(k - 1, i - 1) ** EVEC(j - 1, i - 1);
					}
					/*
					 * Normalize the dot products by
					 * a(k+1,k)
					 */
					s *= t;
					/*
					 * Form the imsl_saxpy of
					 * evec(k+1:n,j) <- evec(k+1:n,j) + s *
					 * a(k+1:n,k)
					 */
					for (i = k + 1; i <= *n; i++) {
						*EVEC(j - 1, i - 1) += s ** A(k - 1, i - 1);
					}
				}
			}
		}
	}
L_220:
	e2[0] = F_ZERO;
	e[0] = F_ZERO;
	imsl_e1pop("l_e6csf ");
	return;
}				/* end of function */


#undef A
#undef EVEC


/* Structured by FOR_STRUCT, v0.2, on 09/24/90 at 11:35:59
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E7CSF/DE7CSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Purpose:    Compute all of the eigenvalues and eigenvectors of
                a real symmetric tridiagonal matrix.

    Usage:      CALL E7CSF (LB,NB,N,EVAL,DWK,E,FNORM,EVEC,LDEVEC,WK)
    Arguments:
       LB     - Start of the current block. Initially set to 1.
                (Input/Output)
       NB     - End of the current block. Initially set to N.
                (Input/Output)
       N      -  Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       EVAL   - Real vector of length N. Contains the diagonal
                terms of the tridiagonal matrix.  (Input)
       DWK    - The eigenvalues as returned by imsl_e5lsf. (Input)
       E      - Real vector of length N containing the n-1 off-diagonal
                elements of the tridiagonal matrix (E(1) = 0).  (Output)
       FNORM  - Frobenius norm used for a scale factor.(Input)
       EVEC   - Accumulated householder transformations. (Input/Output)
                Used to reduce original matrix to tridiagonal form.
                If original matrix is tridiagonal then evec must be
                set to the identity matrix.  Upon output evec is the
                matrix of eigenvectors.
       LDEVEC - Leading dimension of evec.  (Input)
       WK     - Work array. For use on vector machines.

    Remark:
       THIS SUBROUTINE IS BASED ON THE ALGOL PROCEDURE TQL2,
       NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
       WILKINSON.
       HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).

       VERSION FOR THE IBM 3090VF DATED NOVEMBER 1987.

    Keyword:    QL algorithm with explicit shift

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

       ------------------------------------------------------------------ */
/* WK is not used here but leave the calling sequence intact. */
#ifdef ANSI
static void l_e7csf(Mint *lb, Mint *nb, Mint *n, Mfloat eval[],
			Mfloat dwk[], Mfloat e[], Mfloat *fnorm,
			Mfloat *evec, Mint *ldevec, Mfloat wk[])
#else
static void l_e7csf(lb, nb, n, eval, dwk, e, fnorm, evec, ldevec,
	   wk)
	Mint            *lb, *nb, *n;
	Mfloat           eval[], dwk[], e[], *fnorm, *evec;
	Mint            *ldevec;
	Mfloat           wk[];
#endif
{
#define EVEC(I_,J_)	(evec+(I_)*(l_ldevec)+(J_))
        Mint             l_ldevec = *ldevec;
	Mint             i, ics, j, k, l, m;
	Mfloat           c, c2, g, h, p, r, s, s2, sigma, temp, tst1;


	imsl_e1psh("l_e7csf ");

	if (*n == 1)
		goto L_9000;
	/*
	 * Tst1 = 1.0e0/imsl_amach(2) is the bound on smallness. This is used
	 * to prevent soft underflow that is possible on IEEE machines.
	 */
	tst1 = imsl_amach(1);
	if (tst1 * imsl_amach(2) < F_ONE)
		tst1 = F_ONE / imsl_amach(2);
	/*
	 * Simulated DO LOOP. This is done in case there is a finicky
	 * complier that does not like the upper bound on the DO LOOP changed
	 * during iteration. l - DO LOOP variable, initialized to lb and
	 * terminates at nb.
	 */
	l = *lb;
L_10:
	if (l == *nb)
		goto L_9000;
	/*
	 * j represents how many iterations the QR was applied in order to
	 * obtain the l-th eigenvalue.
	 */
	j = 0;
	/*
	 * First test for convergence. dwk(1:n) are the perfect shifts as
	 * returned by imsl_e5lsf. l-th eigenvalue found ?
	 */
L_20:
	if (fabs(e[l]) <= fabs(dwk[l - 1]) * imsl_amach(4)) {
		/* Increase counter l. */
		l += 1;
		/* Check if a Wilkinson shift was used. */
		if (j > 2) {
			/* Set lb to the start of the block. */
			*lb = l;
			/*
			 * Wilkinson's shift was used for the l-th
			 * eigenvalue, exit.
			 */
			goto L_9000;
		}
		/* Start on next eigenvalue. */
		goto L_10;
	}
	/*
	 * Set tst1 to a threshold based on the l-th row. Threshold varies as
	 * l increases.
	 */
	tst1 = imsl_f_max(tst1, imsl_amach(4) * (fabs(eval[l - 1]) + fabs(e[l])));
	/*
	 * Second test for convergence. Same reasoning as in first test.
	 */
	if (fabs(e[l]) <= tst1) {
		l += 1;
		if (j > 2) {
			*lb = l;
			goto L_9000;
		}
		goto L_10;
	}
	/*
	 * Check if current block has decoupled.
	 */
	for (m = l + 1; m <= (*nb - 1); m++) {
		if (fabs(e[m]) <= tst1) {
			*lb = l;
			/*
			 * The current block has decoupled. Exit this
			 * routine. Return to E5LSF.
			 */
			goto L_9000;
		}
	}
	m = *nb;

	/* Increase iteration counter. */
	j += 1;
	/*
	 * Initialize the shift sigma. Use the perfect shift for the first
	 * two iterations.
	 */
	sigma = dwk[l - 1];
	if (j > 2) {
		/* Extra check */
		if ((j >= 5) && (fabs(e[l]) <= *fnorm * imsl_amach(4))) {
			*lb = l + 1;
			goto L_9000;
		}
		/*
		 * Initialize the shift sigma. The Wilkinson shift is used if
		 * more than two shifts are needed.
		 */
		p = (eval[l] - eval[l - 1]) * F_HALF;
		sigma = eval[l - 1] - e[l] / sign(sqrt(imsl_fi_power(p, 2) + imsl_fi_power(e[l], 2)) +
						  fabs(p), p);
	}
	/*
	 * Shift the diagonal elements of the current block.
	 */
	for (i = l; i <= *nb; i++) {
		eval[i - 1] -= sigma;
	}
	/* ql imsl_sweep */
	p = eval[*nb - 1];
	c = F_ONE;
	s = F_ZERO;
	ics = 0;
	for (i = *nb - 1; i >= l; i--) {
		ics += 1;
		c2 = c;
		s2 = s;
		g = c * e[i];
		h = c * p;
		if (fabs(p) < fabs(e[i])) {
			c = p / e[i];
			r = sqrt(F_ONE + c * c);
			e[i + 1] = s2 * e[i] * r;
			s = F_ONE / r;
			c *= s;
		} else {
			s = e[i] / p;
			r = sqrt(F_ONE + s * s);
			e[i + 1] = s2 * p * r;
			c = F_ONE / r;
			s *= c;
		}
		p = c * eval[i - 1] - s * g;
		eval[i] = h + s * (c * g + s * eval[i - 1]);
		/*
		 * eigenvectors ********substitute the input vector wk
		 * instead of scalar temp****** ********on vector machines                                  *******
		 * apply two rotations at once
		 */
		if (ics == 2) {
			for (k = 1; k <= *n; k++) {
				temp = c2 ** EVEC(i, k - 1) - s2 ** EVEC(i + 1, k - 1);
				*EVEC(i + 1, k - 1) = s2 ** EVEC(i, k - 1) + c2 ** EVEC(i + 1, k - 1);
				*EVEC(i, k - 1) = s ** EVEC(i - 1, k - 1) + c * temp;
				*EVEC(i - 1, k - 1) = c ** EVEC(i - 1, k - 1) - s * temp;
			}
			ics = 0;
		}
	}
	/* one rotation may remain */
	if (ics == 1) {
		for (k = 1; k <= *n; k++) {
			temp = s ** EVEC(l - 1, k - 1) + c ** EVEC(l, k - 1);
			*EVEC(l - 1, k - 1) = c ** EVEC(l - 1, k - 1) - s ** EVEC(l, k - 1);
			*EVEC(l, k - 1) = temp;
		}
	}
	e[l] = s * p;
	eval[l - 1] = c * p;
	/* Add shift back. */
	for (i = l; i <= *nb; i++) {
		eval[i - 1] += sigma;
	}
	/*
	 * Iterate back to the first test for convergence.
	 */
	if (j < 100)
		goto L_20;
	/*
	 * If j .ge. 100 then accept eval(l) as the l-th eigenvalue and
	 * deflate.
	 */
	l += 1;
	goto L_10;

L_9000:
	imsl_e1pop("l_e7csf ");
	return;
}				/* end of function */


#undef EVEC





/* Structured by FOR_STRUCT, v0.2, on 09/24/90 at 13:08:40
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  G3LSP/DG3LSP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 31, 1990

    Purpose:    Compute all of the eigenvalues of the generalized real
                symmetric eigenvalue problem A*z = w*B*z, with B
                symmetric positive definite.

    Usage:      CALL G3LSP (N, A, LDA, B, LDB, EVAL, IWK, WK, S)

    Arguments:  (See GVCSP).

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_g3lsp(Mint *n, Mfloat *a, Mint *lda, Mfloat *b,
			Mint *ldb, Mfloat eval[], Mint iwk[],
			Mfloat *wk, Mfloat *s)
#else
static void l_g3lsp(n, a, lda, b, ldb, eval, iwk, wk, s)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda;
	Mfloat          *b;
	Mint            *ldb;
	Mfloat           eval[];
	Mint             iwk[];
	Mfloat          *wk, *s;
#endif
{
#define A(I_,J_)	(a+(I_)*(l_lda)+(J_))
#define B(I_,J_)	(b+(I_)*(l_ldb)+(J_))
#define WK(I_,J_)	(wk+(I_)*(l_n)+(J_))
#define S(I_,J_)	(s+(I_)*(l_n + 1)+(J_))
        Mint            l_lda = *lda;
        Mint            l_ldb = *ldb;
        Mint            l_n   = *n;
	Mint             _l0, i, nucols, nurows;
	Mfloat           _f0, _f1;

	/* First executable statement */
	imsl_e1psh("G3LSP ");
	/* Check N */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "The argument N = %(i1).  The order of the matrix must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);
		goto L_9000;
	}
	/* Check LDA */
	if (*lda < *n) {
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 2, "The argument LDA = %(i1).  The leading dimension of the matrix A must be at least equal to the order of the matrix, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);
	}
	/* Check LDB */
	if (*ldb < *n) {
		imsl_e1sti(1, *ldb);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 3, "The argument LDB = %(i1).  The leading dimension of the matrix B must be at least equal to the order of the matrix, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDB_GE_MATRIX_ORDER);
	}
	if (imsl_n1rty(0) > 0)
		goto L_9000;
	/* 1. Save A, B in n^2 + n locations */
	for (i = 1; i <= *n; i++) {
		scopy(i, A(i - 1, 0), 1, S(i - 1, 0), 1);
		scopy(*n - i + 1, B(i - 1, i - 1), 1, S(i - 1, i), 1);
	}
	/* 2. Cholesky decompose B */
	l_lftds(n, b, ldb, b, ldb);
	if (imsl_n1rty(0) == 4) {

/*		imsl_ermes(4, 2, "Matrix B is not positive definite.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_MATRIX_B_NOT_POS_DEFINITE);
		goto L_9000;
	}
	/*
	 * 3. Compute just the lower triangle of F = A*inv(R); we actually
	 * compute trans(F) = inv(trans(R))*A and store it in the upper
	 * triangular part of A
	 */
	for (i = 1; i <= *n; i++) {

		imsl_strsv("u", "t", "n", i, b, *ldb, A(i - 1, 0), 1);
	}
	/*
	 * 4. Compute just the lower part of G = inv(trans(R))*F; [please see
	 * Wilkinson's AEP, p. 339!]
	 */
	for (i = 1; i <= *n; i++) {
		nurows = i - 1;
		nucols = *n - nurows;
		/*
		 * a. Due to symmetry, we may use relevant, previously
		 * computed elements of G to reduce the size of the system to
		 * be solved
		 */
		_f0 = -F_ONE;
		_l0 = 1;
		_f1 = F_ONE;
		imsl_sgemv("t", sizeof("t"), &nurows, &nucols, &_f0,
		B(i - 1, 0), ldb, A(i - 1, 0), &_l0, &_f1,
			   A(i - 1, i - 1), lda);
		/* b. Now solve the system */

		imsl_strsv("u", "t", "n", nucols, B(i - 1, i - 1), *ldb, A(i - 1, i - 1), *lda);
	}
	/* c. Extend G to a symmetric matrix. */
	l_csfrg(n, a, lda);
	/* 5. Find the eigenvalues of the system */
	l_e4lsf(n, a, lda, eval, wk, iwk);
	/*
	 * Restore A and B from SAVEAB; copy into upper triangles...
	 */
	for (i = 1; i <= *n; i++) {
		scopy(i, S(i - 1, 0), 1, A(i - 1, 0), 1);
		scopy(*n - i + 1, S(i - 1, i), 1, B(i - 1, i - 1), *ldb);
	}
	/* ...and now extend to full matrices. */
	l_csfrg(n, a, lda);
	l_csfrg(n, b, ldb);

L_9000:
	imsl_e1pop("G3LSP ");
	return;
}				/* end of function */


#undef A
#undef B
#undef WK
#undef S



/* Structured by FOR_STRUCT, v0.2, on 09/24/90 at 11:25:53
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  LFTDS/DLFTDS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 13, 1989

    Purpose:    Compute the trans(R)*R Cholesky factorization of a real
                symmetric positive definite matrix.

    Usage:      CALL LFTDS (N, A, LDA, FAC, LDFAC)

    Arguments:
       N      - Order of the matrix.  (Input)
       A      - N by N symmetric positive definite matrix to be factored.
                (Input)
                Only the upper triangle of A is referenced.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       FAC    - N by N matrix containing the upper triangular matrix
                R of the factorization of A in the upper triangle.
                (Output)
                Only the upper triangle of FAC will be used.  If A is not
                needed, A and FAC can share the same storage locations.
       LDFAC  - Leading dimension of FAC exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remark:
       Informational error
       Type Code
         4   2  The input matrix is not positive definite.

    Keyword:    Cholesky factorization

    GAMS:       D2b1b

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------

 */
#define	NB	129

#ifdef ANSI
static void l_lftds(Mint *n, Mfloat *a, Mint *lda, Mfloat *imsl_fac,
			Mint *ldfac)
#else
static void l_lftds(n, a, lda, imsl_fac, ldfac)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda;
	Mfloat          *imsl_fac;
	Mint            *ldfac;
#endif
{
#define A(I_,J_)	(a+(I_)*(l_lda)+(J_))
#define FAC(I_,J_)	(imsl_fac+(I_)*(l_ldfac)+(J_))
        Mint            l_lda = *lda;
        Mint            l_ldfac = *ldfac;
	Mint             i, info, j, k, l, ll, nlenb, nrcfac;
	Mfloat           big, r0, r1, r2, r3, r4, r5, r6, r7, rtemp, small,
	                t[8][NB];


	/*
	 * this code is for computer types: fosivv, rtxlxs, and vxvmsv
	 */
	imsl_e1psh("LFTDS ");

	if (*n <= 0) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "The order of the matrix must be positive while N = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_NOT_POSITIVE);
	}
	if (*n > *lda) {
		imsl_e1sti(1, *n);
		imsl_e1sti(2, *lda);

/*		imsl_ermes(5, 2, "The order of the matrix must be less than or equal to its leading dimension while N = %(i1) and LDA = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_MATRIX_ORDER_LE_LDA);
	}
	if (*n > *ldfac) {
		imsl_e1sti(1, *n);
		imsl_e1sti(2, *ldfac);

/*		imsl_ermes(5, 3, "The order of the matrix must be less than or equal to its leading dimension while N = %(i1) and LDFAC = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_LE_LDFAC);
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;
	/* Preserve a copy of the input matrix */
	for (i = 1; i <= *n; i++) {
		scopy(i, A(i - 1, 0), 1, FAC(i - 1, 0), 1);
	}

	/*
	 * Cholesky decomposition using method LLT**
	 * 
	 * A brief description of the algorithm follows: For a symmetric
	 * positive definite matrix at the k-th step :
	 * 
	 * k   |  11  12  13 | A  = |  21  22  23 | , |  31  32  33 |
	 * 
	 * where trans(A11 A21 A31) is n x (k*8) and factored. assume trans(A22
	 * A32) is the active block of 8 columns. The factorization is
	 * accomplished by:
	 * 
	 * Step 1. Factor trans(A22 A32) and store -trans(A22 A32) into a local
	 * work array
	 * 
	 * Step 2. Update A32 with matrix multiplication between trans(A22 A32)
	 * and the local work array.
	 * 
	 * Step 3. repeat
	 * 
	 * Fill in the lower triangle
	 */
	l_csfrg(n, imsl_fac, ldfac);

	small = imsl_amach(1);
	big = imsl_amach(2);
	if (small * big < F_ONE)
		small = F_ONE / big;
	info = 0;

	nrcfac = mod(*n, 8);
	for (j = 1; j <= (*n - nrcfac); j += 8) {
		nlenb = imsl_i_min(NB, *n - j);
		/* prepare j-th column */
		if (*FAC(j - 1, j - 1) <= F_ZERO) {
			info = j;
			goto L_450;
		}
		*FAC(j - 1, j - 1) = sqrt(*FAC(j - 1, j - 1));
		r0 = F_ONE / *FAC(j - 1, j - 1);
		/*
		 * Form the j-th mutiplier and load t(1:nlenb,1) with
		 * -imsl_fac(j+1:nlenb,j)
		 */
		for (i = 1; i <= nlenb; i++) {
			*FAC(j - 1, j + i - 1) *= r0;
			t[0][i - 1] = -*FAC(j - 1, j + i - 1);
		}
		for (i = j + nlenb + 1; i <= *n; i++) {
			*FAC(j - 1, i - 1) *= r0;
		}
		/* update columns j+1 thru j+7 */
		for (k = 1; k <= 7; k++) {
			for (i = k; i <= (*n - j); i++) {
				*FAC(j + k - 1, j + i - 1) += *FAC(j - 1, j + i - 1) *
					t[0][k - 1];
			}
		}
		/* prepare (j+1)-th column */
		if (*FAC(j, j) <= F_ZERO) {
			info = j + 1;
			goto L_450;
		}
		*FAC(j, j) = sqrt(*FAC(j, j));
		r1 = F_ONE / *FAC(j, j);
		/*
		 * Form the (j+1)-th mutiplier and load t(2:nlenb,2) with
		 * -imsl_fac(j+2:nlenb,j+1)
		 */
		for (i = 2; i <= nlenb; i++) {
			*FAC(j, j + i - 1) *= r1;
			t[1][i - 1] = -*FAC(j, j + i - 1);
		}
		for (i = j + nlenb + 1; i <= *n; i++) {
			*FAC(j, i - 1) *= r1;
		}
		/* update columns j+2 thru j+7 */
		for (k = 2; k <= 7; k++) {
			for (i = k; i <= (*n - j); i++) {
				*FAC(j + k - 1, j + i - 1) += *FAC(j, j + i - 1) * t[1][k - 1];
			}
		}
		/* prepare (j+2)-th column */
		if (*FAC(j + 1, j + 1) <= F_ZERO) {
			info = j + 2;
			goto L_450;
		}
		*FAC(j + 1, j + 1) = sqrt(*FAC(j + 1, j + 1));
		r2 = F_ONE / *FAC(j + 1, j + 1);
		/*
		 * Form the (j+2)-th mutiplier and load t(3:nlenb,3) with
		 * -imsl_fac(j+3:nlenb,j+2)
		 */
		for (i = 3; i <= nlenb; i++) {
			*FAC(j + 1, j + i - 1) *= r2;
			t[2][i - 1] = -*FAC(j + 1, j + i - 1);
		}
		for (i = j + nlenb + 1; i <= *n; i++) {
			*FAC(j + 1, i - 1) *= r2;
		}
		/* update columns j+3 thru j+7 */
		for (k = 3; k <= 7; k++) {
			for (i = k; i <= (*n - j); i++) {
				*FAC(j + k - 1, j + i - 1) += *FAC(j + 1, j + i - 1) *
					t[2][k - 1];
			}
		}
		/* prepare (j+3)-th column */
		if (*FAC(j + 2, j + 2) <= F_ZERO) {
			info = j + 3;
			goto L_450;
		}
		*FAC(j + 2, j + 2) = sqrt(*FAC(j + 2, j + 2));
		r3 = F_ONE / *FAC(j + 2, j + 2);
		/*
		 * Form the (j+3)-th mutiplier and load t(4:nlenb,4) with
		 * -imsl_fac(j+3:nlenb,j+3)
		 */
		for (i = 4; i <= nlenb; i++) {
			*FAC(j + 2, j + i - 1) *= r3;
			t[3][i - 1] = -*FAC(j + 2, j + i - 1);
		}
		for (i = j + nlenb + 1; i <= *n; i++) {
			*FAC(j + 2, i - 1) *= r3;
		}
		/* update columns j+4 thru j+7 */
		for (k = 4; k <= 7; k++) {
			for (i = k; i <= (*n - j); i++) {
				*FAC(j + k - 1, j + i - 1) += *FAC(j + 2, j + i - 1) *
					t[3][k - 1];
			}
		}
		/* prepare (j+4)-th column */
		if (*FAC(j + 3, j + 3) <= F_ZERO) {
			info = j + 4;
			goto L_450;
		}
		*FAC(j + 3, j + 3) = sqrt(*FAC(j + 3, j + 3));
		r4 = F_ONE / *FAC(j + 3, j + 3);
		/*
		 * Form the (j+4)-th mutiplier and load t(5:nlenb,5) with
		 * -imsl_fac(j+4:nlenb,j+4)
		 */
		for (i = 5; i <= nlenb; i++) {
			*FAC(j + 3, j + i - 1) *= r4;
			t[4][i - 1] = -*FAC(j + 3, j + i - 1);
		}
		for (i = j + nlenb + 1; i <= *n; i++) {
			*FAC(j + 3, i - 1) *= r4;
		}
		/* update columns j+5 thru j+7 */
		for (k = 5; k <= 7; k++) {
			for (i = k; i <= (*n - j); i++) {
				*FAC(j + k - 1, j + i - 1) += *FAC(j + 3, j + i - 1) *
					t[4][k - 1];
			}
		}
		/* prepare (j+5)-th column */
		if (*FAC(j + 4, j + 4) <= F_ZERO) {
			info = j + 5;
			goto L_450;
		}
		*FAC(j + 4, j + 4) = sqrt(*FAC(j + 4, j + 4));
		r5 = F_ONE / *FAC(j + 4, j + 4);
		/*
		 * Form the (j+5)-th mutiplier and load t(5:nlenb,5) with
		 * -imsl_fac(j+3:nlenb,j+3)
		 */
		for (i = 6; i <= nlenb; i++) {
			*FAC(j + 4, j + i - 1) *= r5;
			t[5][i - 1] = -*FAC(j + 4, j + i - 1);
		}
		for (i = j + nlenb + 1; i <= *n; i++) {
			*FAC(j + 4, i - 1) *= r5;
		}
		/* update columns j+6 thru j+7 */
		for (k = 6; k <= 7; k++) {
			for (i = k; i <= (*n - j); i++) {
				*FAC(j + k - 1, j + i - 1) += *FAC(j + 4, j + i - 1) *
					t[5][k - 1];
			}
		}
		/* prepare (j+6)-th column */
		if (*FAC(j + 5, j + 5) <= F_ZERO) {
			info = j + 6;
			goto L_450;
		}
		*FAC(j + 5, j + 5) = sqrt(*FAC(j + 5, j + 5));
		r6 = F_ONE / *FAC(j + 5, j + 5);
		/*
		 * Form the (j+6)-th mutiplier and load t(7:nlenb,7) with
		 * -imsl_fac(j+7:nlenb,j+6)
		 */
		for (i = 7; i <= nlenb; i++) {
			*FAC(j + 5, j + i - 1) *= r6;
			t[6][i - 1] = -*FAC(j + 5, j + i - 1);
		}
		for (i = j + nlenb + 1; i <= *n; i++) {
			*FAC(j + 5, i - 1) *= r6;
		}
		/* update column j+7 */
		for (i = 7; i <= (*n - j); i++) {
			*FAC(j + 6, j + i - 1) += *FAC(j + 5, j + i - 1) * t[6][6];
		}
		/* prepare (j+7)-th column */
		if (*FAC(j + 6, j + 6) <= F_ZERO) {
			info = j + 7;
			goto L_450;
		}
		*FAC(j + 6, j + 6) = sqrt(*FAC(j + 6, j + 6));
		r7 = F_ONE / *FAC(j + 6, j + 6);
		/*
		 * Form the (j+7)-th mutiplier and load t(8:nlenb,8) with
		 * -imsl_fac(j+8:nlenb,j+7)
		 */
		for (i = 8; i <= nlenb; i++) {
			*FAC(j + 6, j + i - 1) *= r7;
			t[7][i - 1] = -*FAC(j + 6, j + i - 1);
		}
		for (i = j + nlenb + 1; i <= *n; i++) {
			*FAC(j + 6, i - 1) *= r7;
		}
		/*
		 * Perform update on the lower triangle on columns j+7 thru j
		 * + nlenb rows j+7 thru n
		 */
		for (k = 8; k <= nlenb; k++) {
			for (i = k; i <= (*n - j); i++) {
				*FAC(j + k - 1, j + i - 1) += *FAC(j - 1, j + i - 1) *
					t[0][k - 1] + *FAC(j, j + i - 1) * t[1][k - 1] + *FAC(j + 1, j + i - 1) *
					t[2][k - 1] + *FAC(j + 2, j + i - 1) * t[3][k - 1] +
					*FAC(j + 3, j + i - 1) * t[4][k - 1] + *FAC(j + 4, j + i - 1) *
					t[5][k - 1] + *FAC(j + 5, j + i - 1) * t[6][k - 1] +
					*FAC(j + 6, j + i - 1) * t[7][k - 1];
			}
		}

		for (ll = j + NB + 1; ll <= *n; ll += NB) {
			l = ll - 1;
			nlenb = imsl_i_min(*n - l, NB);
			/*
			 * form mutipliers j,j+1,j+2,j+3 rows ll thru
			 * ll+nlenb
			 */
			for (k = 0; k <= 7; k++) {
				for (i = 1; i <= nlenb; i++) {
					t[k][i - 1] = -*FAC(j + k - 1, l + i - 1);
				}
			}
			/*
			 * Update lower triangle from rows ll thru ll+n
			 * columns ll thru ll + nlenb
			 */
			for (k = 1; k <= nlenb; k++) {
				for (i = k; i <= (*n - l); i++) {
					*FAC(l + k - 1, l + i - 1) += *FAC(j - 1, l + i - 1) *
						t[0][k - 1] + *FAC(j, l + i - 1) * t[1][k - 1] +
						*FAC(j + 1, l + i - 1) * t[2][k - 1] + *FAC(j + 2, l + i - 1) *
						t[3][k - 1] + *FAC(j + 3, l + i - 1) * t[4][k - 1] +
						*FAC(j + 4, l + i - 1) * t[5][k - 1] + *FAC(j + 5, l + i - 1) *
						t[6][k - 1] + *FAC(j + 6, l + i - 1) * t[7][k - 1];
				}
			}
		}

	}
	/*
	 * Take care of remaining nrcfac columns of imsl_fac
	 */
	j = *n - nrcfac + 1;
	if (nrcfac < 7)
		goto L_390;
	if (*FAC(j - 1, j - 1) <= F_ZERO) {
		info = j;
		goto L_450;
	}
	*FAC(j - 1, j - 1) = sqrt(*FAC(j - 1, j - 1));
	rtemp = F_ONE / *FAC(j - 1, j - 1);
	*FAC(j - 1, j) *= rtemp;
	*FAC(j - 1, j + 1) *= rtemp;
	*FAC(j - 1, j + 2) *= rtemp;
	*FAC(j - 1, j + 3) *= rtemp;
	*FAC(j - 1, j + 4) *= rtemp;
	*FAC(j - 1, j + 5) *= rtemp;
	*FAC(j, j) += -imsl_fi_power(*FAC(j - 1, j), 2);
	*FAC(j, j + 1) += -*FAC(j - 1, j + 1) ** FAC(j - 1, j);
	*FAC(j, j + 2) += -*FAC(j - 1, j + 2) ** FAC(j - 1, j);
	*FAC(j, j + 3) += -*FAC(j - 1, j + 3) ** FAC(j - 1, j);
	*FAC(j, j + 4) += -*FAC(j - 1, j + 4) ** FAC(j - 1, j);
	*FAC(j, j + 5) += -*FAC(j - 1, j + 5) ** FAC(j - 1, j);
	*FAC(j + 1, j + 1) += -imsl_fi_power(*FAC(j - 1, j + 1), 2);
	*FAC(j + 1, j + 2) += -*FAC(j - 1, j + 2) ** FAC(j - 1, j + 1);
	*FAC(j + 1, j + 3) += -*FAC(j - 1, j + 3) ** FAC(j - 1, j + 1);
	*FAC(j + 1, j + 4) += -*FAC(j - 1, j + 4) ** FAC(j - 1, j + 1);
	*FAC(j + 1, j + 5) += -*FAC(j - 1, j + 5) ** FAC(j - 1, j + 1);
	*FAC(j + 2, j + 2) += -imsl_fi_power(*FAC(j - 1, j + 2), 2);
	*FAC(j + 2, j + 3) += -*FAC(j - 1, j + 3) ** FAC(j - 1, j + 2);
	*FAC(j + 2, j + 4) += -*FAC(j - 1, j + 4) ** FAC(j - 1, j + 2);
	*FAC(j + 2, j + 5) += -*FAC(j - 1, j + 5) ** FAC(j - 1, j + 2);
	*FAC(j + 3, j + 3) += -imsl_fi_power(*FAC(j - 1, j + 3), 2);
	*FAC(j + 3, j + 4) += -*FAC(j - 1, j + 4) ** FAC(j - 1, j + 3);
	*FAC(j + 3, j + 5) += -*FAC(j - 1, j + 5) ** FAC(j - 1, j + 3);
	*FAC(j + 4, j + 4) += -imsl_fi_power(*FAC(j - 1, j + 4), 2);
	*FAC(j + 4, j + 5) += -*FAC(j - 1, j + 5) ** FAC(j - 1, j + 4);
	*FAC(j + 5, j + 5) += -imsl_fi_power(*FAC(j - 1, j + 5), 2);
	j += 1;
	nrcfac -= 1;
L_390:
	if (nrcfac < 6)
		goto L_400;
	if (*FAC(j - 1, j - 1) <= F_ZERO) {
		info = j;
		goto L_450;
	}
	*FAC(j - 1, j - 1) = sqrt(*FAC(j - 1, j - 1));
	rtemp = F_ONE / *FAC(j - 1, j - 1);
	*FAC(j - 1, j) *= rtemp;
	*FAC(j - 1, j + 1) *= rtemp;
	*FAC(j - 1, j + 2) *= rtemp;
	*FAC(j - 1, j + 3) *= rtemp;
	*FAC(j - 1, j + 4) *= rtemp;
	*FAC(j, j) += -imsl_fi_power(*FAC(j - 1, j), 2);
	*FAC(j, j + 1) += -*FAC(j - 1, j + 1) ** FAC(j - 1, j);
	*FAC(j, j + 2) += -*FAC(j - 1, j + 2) ** FAC(j - 1, j);
	*FAC(j, j + 3) += -*FAC(j - 1, j + 3) ** FAC(j - 1, j);
	*FAC(j, j + 4) += -*FAC(j - 1, j + 4) ** FAC(j - 1, j);
	*FAC(j + 1, j + 1) += -imsl_fi_power(*FAC(j - 1, j + 1), 2);
	*FAC(j + 1, j + 2) += -*FAC(j - 1, j + 2) ** FAC(j - 1, j + 1);
	*FAC(j + 1, j + 3) += -*FAC(j - 1, j + 3) ** FAC(j - 1, j + 1);
	*FAC(j + 1, j + 4) += -*FAC(j - 1, j + 4) ** FAC(j - 1, j + 1);
	*FAC(j + 2, j + 2) += -imsl_fi_power(*FAC(j - 1, j + 2), 2);
	*FAC(j + 2, j + 3) += -*FAC(j - 1, j + 3) ** FAC(j - 1, j + 2);
	*FAC(j + 2, j + 4) += -*FAC(j - 1, j + 4) ** FAC(j - 1, j + 2);
	*FAC(j + 3, j + 3) += -imsl_fi_power(*FAC(j - 1, j + 3), 2);
	*FAC(j + 3, j + 4) += -*FAC(j - 1, j + 4) ** FAC(j - 1, j + 3);
	*FAC(j + 4, j + 4) += -imsl_fi_power(*FAC(j - 1, j + 4), 2);
	j += 1;
	nrcfac -= 1;
L_400:
	if (nrcfac < 5)
		goto L_410;
	if (*FAC(j - 1, j - 1) <= F_ZERO) {
		info = j;
		goto L_450;
	}
	*FAC(j - 1, j - 1) = sqrt(*FAC(j - 1, j - 1));
	rtemp = F_ONE / *FAC(j - 1, j - 1);
	*FAC(j - 1, j) *= rtemp;
	*FAC(j - 1, j + 1) *= rtemp;
	*FAC(j - 1, j + 2) *= rtemp;
	*FAC(j - 1, j + 3) *= rtemp;
	*FAC(j, j) += -imsl_fi_power(*FAC(j - 1, j), 2);
	*FAC(j, j + 1) += -*FAC(j - 1, j + 1) ** FAC(j - 1, j);
	*FAC(j, j + 2) += -*FAC(j - 1, j + 2) ** FAC(j - 1, j);
	*FAC(j, j + 3) += -*FAC(j - 1, j + 3) ** FAC(j - 1, j);
	*FAC(j + 1, j + 1) += -imsl_fi_power(*FAC(j - 1, j + 1), 2);
	*FAC(j + 1, j + 2) += -*FAC(j - 1, j + 2) ** FAC(j - 1, j + 1);
	*FAC(j + 1, j + 3) += -*FAC(j - 1, j + 3) ** FAC(j - 1, j + 1);
	*FAC(j + 2, j + 2) += -imsl_fi_power(*FAC(j - 1, j + 2), 2);
	*FAC(j + 2, j + 3) += -*FAC(j - 1, j + 3) ** FAC(j - 1, j + 2);
	*FAC(j + 3, j + 3) += -imsl_fi_power(*FAC(j - 1, j + 3), 2);
	j += 1;
	nrcfac -= 1;
L_410:
	if (nrcfac < 4)
		goto L_420;
	if (*FAC(j - 1, j - 1) <= F_ZERO) {
		info = j;
		goto L_450;
	}
	*FAC(j - 1, j - 1) = sqrt(*FAC(j - 1, j - 1));
	rtemp = F_ONE / *FAC(j - 1, j - 1);
	*FAC(j - 1, j) *= rtemp;
	*FAC(j - 1, j + 1) *= rtemp;
	*FAC(j - 1, j + 2) *= rtemp;
	*FAC(j, j) += -imsl_fi_power(*FAC(j - 1, j), 2);
	*FAC(j, j + 1) += -*FAC(j - 1, j + 1) ** FAC(j - 1, j);
	*FAC(j, j + 2) += -*FAC(j - 1, j + 2) ** FAC(j - 1, j);
	*FAC(j + 1, j + 1) += -imsl_fi_power(*FAC(j - 1, j + 1), 2);
	*FAC(j + 1, j + 2) += -*FAC(j - 1, j + 2) ** FAC(j - 1, j + 1);
	*FAC(j + 2, j + 2) += -imsl_fi_power(*FAC(j - 1, j + 2), 2);
	j += 1;
	nrcfac -= 1;
L_420:
	if (nrcfac < 3)
		goto L_430;
	if (*FAC(j - 1, j - 1) <= F_ZERO) {
		info = j;
		goto L_450;
	}
	*FAC(j - 1, j - 1) = sqrt(*FAC(j - 1, j - 1));
	rtemp = F_ONE / *FAC(j - 1, j - 1);
	*FAC(j - 1, j) *= rtemp;
	*FAC(j - 1, j + 1) *= rtemp;
	*FAC(j, j) += -imsl_fi_power(*FAC(j - 1, j), 2);
	*FAC(j, j + 1) += -*FAC(j - 1, j + 1) ** FAC(j - 1, j);
	*FAC(j + 1, j + 1) += -imsl_fi_power(*FAC(j - 1, j + 1), 2);
	j += 1;
	nrcfac -= 1;
L_430:
	if (nrcfac < 2)
		goto L_440;
	if (*FAC(j - 1, j - 1) <= F_ZERO) {
		info = j;
		goto L_450;
	}
	*FAC(j - 1, j - 1) = sqrt(*FAC(j - 1, j - 1));
	rtemp = F_ONE / *FAC(j - 1, j - 1);
	*FAC(j - 1, j) *= rtemp;
	*FAC(j, j) += -imsl_fi_power(*FAC(j - 1, j), 2);
	j += 1;
	nrcfac -= 1;
L_440:
	if (nrcfac < 1)
		goto L_450;
	if (*FAC(j - 1, j - 1) <= F_ZERO) {
		info = j;
		goto L_450;
	}
	*FAC(j - 1, j - 1) = sqrt(*FAC(j - 1, j - 1));
L_450:
	if (info != 0) {
		imsl_e1sti(1, info);

/*		imsl_ermes(4, 2, "The leading %(i1) by %(i1) submatrix of the input matrix is not positive definite.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_SUBMATRIX_NOT_POS_DEFINITE);
	}
	/* Fill in the upper triangle */
	for (i = 1; i <= (*n - 1); i++) {
		scopy(*n - i, FAC(i - 1, i), 1, FAC(i, i - 1), *ldfac);
	}

L_9000:
	imsl_e1pop("LFTDS ");

	return;
}				/* end of function */


#undef NB
#undef A
#undef FAC


/* Structured by FOR_STRUCT, v0.2, on 09/24/90 at 11:30:45
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  STRSM (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    July 27, 1989

    Purpose:    Solves one of the matrix equations:
                    op( A )*X = alpha*B,
                 or
                    X*op( A ) = alpha*B,
                where alpha is a scalar, X and B are m by n matrices, A
                is a unit, or non-unit, upper or lower triangular matrix
                and op( A ) is one of:
                    op( A ) = A
                 or
                    op( A ) = trans(A).
                Trans represents the transpose of the matrix.
                The matrix X is overwritten on B.

    Usage:      CALL STRSM (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A,
                            LDA, B, LDB)

    Arguments:
       SIDE   - Character string specifying whether op( A ) appears on
                the left or right of X.  (Input)
                The operations are as follows:
                   If SIDE = 'L' or 'l', then op( A )*X = alpha*B.
                   If SIDE = 'R' or 'r', then X*op( A ) = alpha*B.
       UPLO   - Character string specifying whether the matrix A is an
                upper or lower triangular matrix.  (Input)
                The operations are as follows:
                   If UPLO = 'U' or 'u', then A is an upper triangular
                   matrix.
                   If UPLO = 'L' or 'l', then A is a lower triangular
                   matrix.
       TRANSA - Character string specifying the form of op( A ) to be
                used in the matrix multiplication.  (Input)
                The operations are as follows:
                   If TRANSA = 'N' or 'n', then op( A ) = A.
                   If TRANSA = 'T' or 't', then op( A ) = trans(A).
                   If TRANSA = 'C' or 'c', then op( A ) = trans(A).
       DIAG   - Character string specifying whether or not A is unit
                triangular.  (Input)
                The operations are as follows:
                   If DIAG = 'U' or 'u', then A is assumed to be unit
                   triangular.
                   If DIAG = 'N' or 'n', then A is assumed to be non-unit
                   triangular.
       M      - M specifies the number of rows of B.  (Input)
       N      - N specifies the number of columns of B.  (Input)
       ALPHA  - Scalar multiplier for B.  (Input)
                When alpha is zero then A is not referenced and B need
                not be set on input.
       A      - Array of dimension ( LDA, K ), where K is M when
                SIDE = 'L' or 'l' and is N when  SIDE = 'R' or 'r'.
                (Input)
                If UPLO = 'U' or 'u', the strictly lower triangular part
                of A is not referenced.  If UPLO = 'L' or 'l', the
                strictly upper triangular part of A is not referenced.
                Note that when DIAG = 'U' or 'u', the diagonal elements
                of A are assumed to be unity.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       B      - Array of dimension ( LDB, N ).  (Input/Output)
                On input, B must contain the right-hand side matrix.
                On output, B contains the solution matrix.
       LDB    - Leading dimension of B exactly as specified in the
                dimension statement of the calling program.  (Input)

    GAMS:       D1b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	ONE	1.0e0
#define	ZERO	0.0e0

#ifdef ANSI
static void l_strsm(Mchar *side, unsigned side_s, Mchar *uplo,
			unsigned uplo_s, Mchar *transa, unsigned
			transa_s, Mchar *diag, unsigned diag_s,
			Mint *m, Mint *n, Mfloat *alpha, Mfloat *a,
			Mint *lda, Mfloat *b, Mint *ldb)
#else
static void l_strsm(side, side_s, uplo, uplo_s, transa, transa_s,
	   diag, diag_s, m, n, alpha, a, lda, b, ldb)
	Mchar           *side;
	unsigned        side_s;
	Mchar           *uplo;
	unsigned        uplo_s;
	Mchar           *transa;
	unsigned        transa_s;
	Mchar           *diag;
	unsigned        diag_s;
	Mint            *m, *n;
	Mfloat          *alpha, *a;
	Mint            *lda;
	Mfloat          *b;
	Mint            *ldb;
#endif
{
#define A(I_,J_)	(a+(I_)*(l_lda)+(J_))
#define B(I_,J_)	(b+(I_)*(l_ldb)+(J_))
        Mint            l_lda = *lda;
        Mint            l_ldb = *ldb;
	Mint            imsl_ctran, lower, lside, ndiag,
	                ntran, rside, tran, udiag, upper;
	Mint             i, j, k;
	Mfloat           temp;


	lside = imsl_l1ame(side, side_s, "L", sizeof("L"));
	rside = imsl_l1ame(side, side_s, "R", sizeof("R"));
	upper = imsl_l1ame(uplo, uplo_s, "U", sizeof("U"));
	lower = imsl_l1ame(uplo, uplo_s, "L", sizeof("L"));
	udiag = imsl_l1ame(diag, diag_s, "U", sizeof("U"));
	ndiag = imsl_l1ame(diag, diag_s, "N", sizeof("N"));
	ntran = imsl_l1ame(transa, transa_s, "N", sizeof("N"));
	tran = imsl_l1ame(transa, transa_s, "T", sizeof("T"));
	imsl_ctran = imsl_l1ame(transa, transa_s, "C", sizeof("C"));
	/*
	 * Test the input parameters.
	 */
	if (*m < 0) {
		imsl_e1psh("STRSM ");
		imsl_e1sti(1, *m);

/*		imsl_ermes(5, 1, "M must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_M_GE_ZERO);
		imsl_e1pop("STRSM ");
		goto L_9000;
	} else if (*n < 0) {
		imsl_e1psh("STRSM ");
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 2, "N must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
		imsl_e1pop("STRSM ");
		goto L_9000;
	} else if ((*ldb < *m) || (*ldb == 0)) {
		imsl_e1psh("STRSM ");
		imsl_e1sti(1, *ldb);
		imsl_e1sti(2, *m);

/*		imsl_ermes(5, 3, "LDB must be greater than or equal to M and greater than zero while LDB = %(i1) and M = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDB_GE_M_AND_GT_ZERO);
		imsl_e1pop("STRSM ");
		goto L_9000;
	} else if (lside) {
		if ((*lda < *m) || (*lda == 0)) {
			imsl_e1psh("STRSM ");
			imsl_e1sti(1, *lda);
			imsl_e1sti(2, *m);

/*			imsl_ermes(5, 4, "LDA must be greater than or equal to M and greater than zero when SIDE = L, but LDA = %(i1) and M = %(i2) are given.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_SIDE_EQ_L_LDA_TOO_SMALL);
			imsl_e1pop("STRSM ");
			goto L_9000;
		}
	} else if (rside) {
		if ((*lda < *n) || (*lda == 0)) {
			imsl_e1psh("STRSM ");
			imsl_e1sti(1, *lda);
			imsl_e1sti(2, *n);

/*			imsl_ermes(5, 5, "LDA must be greater than or equal to N and greater than zero when SIDE = R, but LDA = %(i1) and N = %(i2) are given.");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_NEED_NEW_LDA_SIDE_EQUALS_R);
			imsl_e1pop("STRSM ");
			goto L_9000;
		}
	} else {
		imsl_e1psh("STRSM ");
		imsl_e1stl(1, side);

/*		imsl_ermes(5, 6, "SIDE must be set equal to L or R while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_SIDE_MUST_EQUAL_L_OR_R);
		imsl_e1pop("STRSM ");
		goto L_9000;
	}

	if ((!upper) && (!lower)) {
		imsl_e1psh("STRSM ");
		imsl_e1stl(1, uplo);

/*		imsl_ermes(5, 7, "UPLO must be set equal to U or L while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_UPLO_MUST_EQUAL_U_OR_L);
		imsl_e1pop("STRSM ");
		goto L_9000;
	} else if (((!ntran) && (!tran)) && (!imsl_ctran)) {
		imsl_e1psh("STRSM");
		imsl_e1stl(1, transa);

/*		imsl_ermes(5, 8, "TRANSA must be set equal to N, T, or C while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_TRANSA_MUST_EQUAL_N_T_OR_C);
		imsl_e1pop("STRSM");
		goto L_9000;
	} else if ((!udiag) && (!ndiag)) {
		imsl_e1psh("STRSM");
		imsl_e1stl(1, diag);

/*		imsl_ermes(5, 9, "DIAG must be set equal to U or N while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_DIAG_MUST_EQUAL_U_OR_N);
		imsl_e1pop("STRSM");
		goto L_9000;
	}
	/*
	 * Quick return if possible.
	 */
	if (*n == 0)
		goto L_9000;

	if (*alpha == ZERO) {
		for (j = 1; j <= *n; j++) {
			for (i = 1; i <= *m; i++) {
				*B(j - 1, i - 1) = ZERO;
			}
		}
		goto L_9000;
	}
	/*
	 * Main computational section
	 */
	if (lside) {
		if (ntran) {
			/*
			 * Compute: B = alpha*inv( A )*B.
			 */
			if (upper) {
				/* Upper triangular case */
				for (j = 1; j <= *n; j++) {

					if (*alpha != ONE) {
						for (i = 1; i <= *m; i++) {
							*B(j - 1, i - 1) *= *alpha;
						}
					}
					for (k = *m; k >= 1; k--) {

						if (*B(j - 1, k - 1) != ZERO) {
							if (ndiag)
								*B(j - 1, k - 1) /= *A(k - 1, k - 1);
							for (i = 1; i <= (k - 1); i++) {
								*B(j - 1, i - 1) += -*B(j - 1, k - 1) *
									*A(k - 1, i - 1);
							}
						}
					}
				}
			} else {
				/* Lower triangular case */
				for (j = 1; j <= *n; j++) {

					if (*alpha != ONE) {
						for (i = 1; i <= *m; i++) {
							*B(j - 1, i - 1) *= *alpha;
						}
					}
					for (k = 1; k <= *m; k++) {

						if (*B(j - 1, k - 1) != ZERO) {
							if (ndiag)
								*B(j - 1, k - 1) /= *A(k - 1, k - 1);
							for (i = k + 1; i <= *m; i++) {
								*B(j - 1, i - 1) += -*B(j - 1, k - 1) *
									*A(k - 1, i - 1);
							}
						}
					}
				}
			}
		} else {
			/*
			 * Compute: B = alpha*inv( trans(A) )*B.
			 */
			if (upper) {
				/* Upper triangular case */
				for (j = 1; j <= *n; j++) {
					for (i = 1; i <= *m; i++) {
						temp = *alpha ** B(j - 1, i - 1);

						for (k = 1; k <= (i - 1); k++) {
							temp += -*A(i - 1, k - 1) ** B(j - 1, k - 1);
						}

						if (ndiag)
							temp /= *A(i - 1, i - 1);
						*B(j - 1, i - 1) = temp;
					}
				}
				/* Lower triangular case */
			} else {
				for (j = 1; j <= *n; j++) {
					for (i = *m; i >= 1; i--) {
						temp = *alpha ** B(j - 1, i - 1);

						for (k = i + 1; k <= *m; k++) {
							temp += -*A(i - 1, k - 1) ** B(j - 1, k - 1);
						}

						if (ndiag)
							temp /= *A(i - 1, i - 1);
						*B(j - 1, i - 1) = temp;
					}
				}
			}
		}
	} else {
		if (ntran) {
			/*
			 * Compute: B = alpha*B*inv( A ).
			 */
			if (upper) {
				/* Upper triangular case */
				for (j = 1; j <= *n; j++) {

					if (*alpha != ONE) {
						for (i = 1; i <= *m; i++) {
							*B(j - 1, i - 1) *= *alpha;
						}
					}
					for (k = 1; k <= (j - 1); k++) {

						if (*A(j - 1, k - 1) != ZERO) {
							for (i = 1; i <= *m; i++) {
								*B(j - 1, i - 1) += -*A(j - 1, k - 1) *
									*B(k - 1, i - 1);
							}
						}
					}

					if (ndiag) {
						temp = ONE / *A(j - 1, j - 1);
						for (i = 1; i <= *m; i++) {
							*B(j - 1, i - 1) *= temp;
						}
					}
				}
				/* Lower triangular case */
			} else {
				for (j = *n; j >= 1; j--) {

					if (*alpha != ONE) {
						for (i = 1; i <= *m; i++) {
							*B(j - 1, i - 1) *= *alpha;
						}
					}
					for (k = j + 1; k <= *n; k++) {
						if (*A(j - 1, k - 1) != ZERO) {

							for (i = 1; i <= *m; i++) {
								*B(j - 1, i - 1) += -*A(j - 1, k - 1) *
									*B(k - 1, i - 1);
							}

						}
					}

					if (ndiag) {
						temp = ONE / *A(j - 1, j - 1);
						for (i = 1; i <= *m; i++) {
							*B(j - 1, i - 1) *= temp;
						}
					}
				}
			}
		} else {
			/*
			 * Compute: B = alpha*B*inv( trans(A) ).
			 */
			if (upper) {
				/* Upper triangular case */
				for (k = *n; k >= 1; k--) {

					if (ndiag) {
						temp = ONE / *A(k - 1, k - 1);
						for (i = 1; i <= *m; i++) {
							*B(k - 1, i - 1) *= temp;
						}
					}
					for (j = 1; j <= (k - 1); j++) {
						if (*A(k - 1, j - 1) != ZERO) {
							temp = *A(k - 1, j - 1);
							for (i = 1; i <= *m; i++) {
								*B(j - 1, i - 1) += -temp ** B(k - 1, i - 1);
							}
						}
					}

					if (*alpha != ONE) {
						for (i = 1; i <= *m; i++) {
							*B(k - 1, i - 1) *= *alpha;
						}
					}
				}
				/* Lower triangular case */
			} else {
				for (k = 1; k <= *n; k++) {

					if (ndiag) {
						temp = ONE / *A(k - 1, k - 1);
						for (i = 1; i <= *m; i++) {
							*B(k - 1, i - 1) *= temp;
						}
					}
					for (j = k + 1; j <= *n; j++) {
						if (*A(k - 1, j - 1) != ZERO) {
							temp = *A(k - 1, j - 1);
							for (i = 1; i <= *m; i++) {
								*B(j - 1, i - 1) += -temp ** B(k - 1, i - 1);
							}
						}
					}

					if (*alpha != ONE) {
						for (i = 1; i <= *m; i++) {
							*B(k - 1, i - 1) *= *alpha;
						}
					}
				}
			}
		}
	}
L_9000:
	return;
}				/* end of function */


#undef ONE
#undef ZERO
#undef A
#undef B
/* Structured by FOR_STRUCT, v0.2, on 11/16/90 at 11:58:40
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  G3CSP/DG3CSP (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 31, 1990

    Purpose:    Compute all of the eigenvalues and eigenvectors of the
                generalized real symmetric eigenvalue problem
                A*z = w*B*z, with B symmetric positive definite.

    Usage:      CALL G3CSP (N, A, LDA, B, LDB, EVAL, EVEC, LDEVEC,
                            IWK, WK1, WK2)

    Arguments:  (See GVCSP).

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_g3csp(Mint *n, Mfloat *a, Mint *lda, Mfloat *b,
                        Mint *ldb, Mfloat eval[], Mfloat *evec,
                        Mint *ldevec, Mint iwk[], Mfloat *wk1,
                        Mfloat *wk2)
#else
static void l_g3csp (n, a, lda, b, ldb, eval, evec, ldevec, iwk,
                wk1, wk2)
    Mint        *n;
    Mfloat      *a;
    Mint        *lda;
    Mfloat      *b;
    Mint        *ldb;
    Mfloat       eval[], *evec;
    Mint        *ldevec, iwk[];
    Mfloat      *wk1, *wk2;
#endif
{
#define A(I_,J_)	(a+(I_)*(l_lda)+(J_))
#define B(I_,J_)	(b+(I_)*(l_ldb)+(J_))
#define EVEC(I_,J_)	(evec+(I_)*(l_ldevec)+(J_))
#define WK1(I_,J_)	(wk1+(I_)*(l_n)+(J_))
#define WK2(I_,J_)	(wk2+(I_)*(l_n + 1)+(J_))
        Mint            l_lda = *lda;
        Mint            l_ldb = *ldb;
        Mint            l_ldevec = *ldevec;
        Mint            l_n = *n;

    Mint         _l0, i, j, k, nucols, nurows;
    Mfloat       _f0, _f1, scale;
    /* First executable statement */
    imsl_e1psh ("G3CSP ");
    /* Check N */
    if (*n < 1) {
	imsl_e1sti (1, *n);

	imsl_ermes (IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);
	goto L_9000;
    }
    /* Check LDA */
    if (*lda < *n) {
	imsl_e1sti (1, *lda);
	imsl_e1sti (2, *n);

	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);
    }
    /* Check LDB */
    if (*ldb < *n) {
	imsl_e1sti (1, *ldb);
	imsl_e1sti (2, *n);

	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_LDB_GE_MATRIX_ORDER);
    }
    /* Check LDEVEC */
    if (*ldevec < *n) {
	imsl_e1sti (1, *ldevec);
	imsl_e1sti (2, *n);

	imsl_ermes (IMSL_TERMINAL, IMSL_LDEVEC_VALUE_TOO_SMALL);
    }
    if (imsl_n1rty (0) > 0)
	goto L_9000;
    /* 1. Save A, B in n^2 + n locations */
    for (i = 1; i <= *n; i++) {
	imsl_scopy (i, A (i - 1, 0), 1, WK2 (i - 1, 0), 1);
	imsl_scopy (*n - i + 1, B (i - 1, i - 1), 1, WK2 (i - 1, i), 1);
    }
    /* 2. Cholesky decompose B */
    l_lftds (n, b, ldb, b, ldb);
    if (imsl_n1rty (0) == 4) {

	imsl_ermes (IMSL_FATAL, IMSL_MATRIX_B_NOT_POS_DEFINITE);
	goto L_9000;
    }
    /*
     * 3. Compute just the lower triangle of F = A*inv(R); we actually
     * compute trans(F) = inv(trans(R))*A and store it in the upper
     * triangular part of A
     */
    for (i = 1; i <= *n; i++) {

	imsl_strsv ("u", "t", "n", i, b, *ldb, A (i - 1, 0), 1);
    }
    /*
     * 4. Compute just the lower part of G = inv(trans(R))*F; [please see
     * Wilkinson's AEP, p. 339!]
     */
    for (i = 1; i <= *n; i++) {
	nurows = i - 1;
	nucols = *n - nurows;
	/*
	 * a. Due to symmetry, we may use relevant, previously computed
	 * elements of G to reduce the size of the system to be solved
	 */
	_f0 = -1.0;
	_l0 = 1;
	_f1 = 1.0;
	imsl_sgemv ("t", sizeof ("t"), &nurows, &nucols, &_f0,
	    B (i - 1, 0), ldb, A (i - 1, 0), &_l0, &_f1,
	    A (i - 1, i - 1), lda);
	/* b. Now solve the system */

	imsl_strsv ("u", "t", "n", nucols, B (i - 1, i - 1), *ldb, A (i - 1, i - 1), *lda);
    }
    /* c. Extend G to a symmetric matrix. */
    imsl_csfrg (n, a, lda);
    /*
     * 5. Compute eigensystem expansion of symmetric matrix G = X*D*trans(X)
     */
    l_e5csf (n, a, lda, eval, evec, ldevec, wk1, iwk);
    /* 6. Restore coordinates, X = inv(R)*X */
	_f0 = 1.0;
    l_strsm ("l", sizeof ("l"), "u", sizeof ("u"), "n", sizeof ("n"), "n", sizeof ("n"
	), n, n, &_f0, b, ldb, evec, ldevec);
    /*
     * 7. Normalize eigenvectors using the Euclidean Norm
     */
    for (i = 1; i <= *n; i++) {
	scale = imsl_snrm2 (*n, EVEC (i - 1, 0), 1);
	if (scale > 0.0) {
	    sscal (*n, 1.0 / scale, EVEC (i - 1, 0), 1);
	}
    }
    /*
     * Normalize each eigenvector so that its biggest component is positive.
     * The eigenvectors then form a right-hand system.
     */
    for (j = 1; j <= *n; j++) {
	i = imsl_isamax (*n, EVEC (j - 1, 0), 1);
	if (*EVEC (j - 1, i - 1) < 0.0) {
	    for (k = 1; k <= *n; k++) {
		*EVEC (j - 1, k - 1) = -*EVEC (j - 1, k - 1);
	    }
	}
    }
    /*
     * Restore A and B from SAVEAB; copy into upper triangles...
     */
    for (i = 1; i <= *n; i++) {
	imsl_scopy (i, WK2 (i - 1, 0), 1, A (i - 1, 0), 1);
	imsl_scopy (*n - i + 1, WK2 (i - 1, i), 1, B (i - 1, i - 1), *ldb);
    }
    /* ...and now extend to full matrices. */
    imsl_csfrg (n, a, lda);
    imsl_csfrg (n, b, ldb);

L_9000:
    ;
    imsl_e1pop ("G3CSP ");
    return;
}				/* end of function */

#undef A
#undef B
#undef EVEC
#undef WK1
#undef WK2
