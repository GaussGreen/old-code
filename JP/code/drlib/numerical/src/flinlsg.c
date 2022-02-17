#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static VA_LIST_HACK PROTO (l_lin_least_squares_gen, (Mint nra, Mint nca, Mfloat *a, Mfloat *b, va_list argptr));
static void PROTO (l_l2err, (Mint *nrqr, Mint *ncqr, Mfloat *qr, Mint *ldqr,
	            Mfloat *qraux, Mfloat *q, Mint *ldq, Mfloat *work));
    static void PROTO (l_permu, (Mint *n, Mfloat *x, Mint *ipermu, Mint *ipath, Mfloat *xpermu));
    static Mint PROTO (l_ismax, (Mint *n, Mfloat *sx, Mint *incx));

    static Mfloat *lv_x;
#ifdef ANSI
    Mfloat     *imsl_f_lin_least_squares_gen (Mint nra, Mint nca, Mfloat *a, Mfloat *b,...)
#else
    Mfloat     *imsl_f_lin_least_squares_gen (nra, nca, a, b, va_alist)
    Mint        nra, nca;
    Mfloat     *a, *b;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, b);

    E1PSH ("imsl_f_lin_least_squares_gen", "imsl_d_lin_least_squares_gen");

    lv_x = NULL;
    IMSL_CALL (l_lin_least_squares_gen (nra, nca, a, b, argptr));
    va_end (argptr);

    E1POP ("imsl_f_lin_least_squares_gen", "imsl_d_lin_least_squares_gen");

    return lv_x;
}
#ifdef ANSI
static VA_LIST_HACK l_lin_least_squares_gen (Mint nra, Mint nca, Mfloat *a, Mfloat *b, va_list argptr)
#else
static VA_LIST_HACK l_lin_least_squares_gen (nra, nca, a, b, argptr)
    Mint        nra, nca;
    Mfloat     *a, *b;
    va_list     argptr;
#endif
{
    Mint        col_dim_a = nca;
    Mfloat      tol;
    Mint       *kbasis = NULL;
    Mint       *ipvt = NULL;
    Mfloat    **p_res = NULL;
    Mfloat     *res = NULL;
    Mfloat    **p_qraux;
    Mfloat     *qraux = NULL;
    Mfloat    **p_qr;
    Mfloat     *qr = NULL;
    Mint        col_dim_fac = nca;
    Mfloat    **p_q;
    Mfloat     *q = NULL;
    Mint        col_dim_q = nra;

    Mint        user_res = 0, return_res = 0;
    Mint        user_factor = 0, return_factor = 0;
    Mint        user_q = 0, return_q = 0;
    Mint        factor_only = 0, solve_only = 0;
    Mint        return_user = 0;
    Mint        error = 0;
    Mint        user_ipvt = 0;
    Mint        pivot = 0;
    Mint        i;

    Mint        code = 1, arg_number = 4;

    Mfloat     *work = NULL;
    tol = sqrt (imsl_amach (4));

    while (code > 0) {
	code = va_arg (argptr, Mint);
	++arg_number;
	switch (code) {
	case IMSL_A_COL_DIM:
	    col_dim_a = va_arg (argptr, Mint);
	    ++arg_number;
	    break;
	case IMSL_RETURN_USER:
	    lv_x = va_arg (argptr, Mfloat *);
	    if (!lv_x) {
		imsl_e1stl (1, "x");
		imsl_e1stl (2, "IMSL_RETURN_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    ++arg_number;
	    return_user = 1;
	    break;
	case IMSL_BASIS:
	    tol = (Mfloat) va_arg (argptr, Mdouble);
	    ++arg_number;
	    kbasis = va_arg (argptr, Mint *);
	    ++arg_number;
	    break;
	case IMSL_BASIS_ADR:
	    tol =  *(va_arg (argptr, Mfloat *));
	    ++arg_number;
	    kbasis = va_arg (argptr, Mint *);
	    ++arg_number;
	    break;
	case IMSL_PIVOT:
	    ipvt = va_arg (argptr, Mint *);
	    if (!ipvt) {
		imsl_e1stl (1, "ipvt");
		imsl_e1stl (2, "IMSL_PIVOT");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    ++arg_number;
	    user_ipvt = 1;
	    break;
	case IMSL_RESIDUAL:
	    p_res = va_arg (argptr, Mfloat **);
	    ++arg_number;
	    user_res = 0;
	    return_res = 1;
	    break;
	case IMSL_RESIDUAL_USER:
	    res = va_arg (argptr, Mfloat *);
	    if (!res) {
		imsl_e1stl (1, "res");
		imsl_e1stl (2, "IMSL_RESIDUAL_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    ++arg_number;
	    user_res = 1;
	    return_res = 1;
	    break;
	case IMSL_FACTOR:
	    p_qraux = va_arg (argptr, Mfloat **);
	    ++arg_number;
	    p_qr = va_arg (argptr, Mfloat **);
	    ++arg_number;
	    user_factor = 0;
	    return_factor = 1;
	    break;
	case IMSL_FACTOR_USER:
	    qraux = va_arg (argptr, Mfloat *);
	    if (!qraux) {
		imsl_e1stl (1, "qraux");
		imsl_e1stl (2, "IMSL_FACTOR_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    ++arg_number;
	    qr = va_arg (argptr, Mfloat *);
	    if (!qr) {
		imsl_e1stl (1, "qr");
		imsl_e1stl (2, "IMSL_FACTOR_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_2);
		++error;
	    }
	    ++arg_number;
	    user_factor = 1;
	    return_factor = 1;
	    break;
	case IMSL_FAC_COL_DIM:
	    col_dim_fac = va_arg (argptr, Mint);
	    ++arg_number;
	    break;
	case IMSL_Q:
	    p_q = va_arg (argptr, Mfloat **);
	    ++arg_number;
	    user_q = 0;
	    return_q = 1;
	    break;
	case IMSL_Q_USER:
	    q = va_arg (argptr, Mfloat *);
	    if (!q) {
		imsl_e1stl (1, "q");
		imsl_e1stl (2, "IMSL_Q_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_2);
		++error;
	    }
	    ++arg_number;
	    user_q = 1;
	    return_q = 1;
	    break;
	case IMSL_Q_COL_DIM:
	    col_dim_q = va_arg (argptr, Mint);
	    ++arg_number;
	    break;
	case IMSL_FACTOR_ONLY:
	    factor_only = 1;
	    break;
	case IMSL_SOLVE_ONLY:
	    solve_only = 1;
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
    if ((!solve_only)&&(!a)) {
	imsl_e1stl (1, "a");
	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
	++error;
    }
    if (!b && !factor_only) {
	imsl_e1stl (1, "b");
	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
	++error;
    }
    if (error)
	return argptr;

    if (nra <= 0 || nca <= 0) {
	imsl_e1sti (1, nra);
	imsl_e1sti (2, nca);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NRA_AND_NCA_GT_ZERO);
	++error;
    }
    if (nca > col_dim_a) {
	imsl_e1sti (1, nra);
	imsl_e1sti (2, col_dim_a);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_COL_A_LT_DIM);
	++error;
    }
    if (col_dim_fac < nca) {
	imsl_e1sti (1, nra);
	imsl_e1sti (2, col_dim_fac);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_QT_GT_NUM_COL);
	++error;
    }
    if (col_dim_q < nra) {
	imsl_e1sti (1, nra);
	imsl_e1sti (2, col_dim_q);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_Q_GT_NUM_COL);
	++error;
    }
    if (factor_only + solve_only > 1) {
	imsl_ermes (IMSL_TERMINAL, IMSL_BAD_SOLVE_FACTOR);
	++error;
    }
    if (error)
	return argptr;

    if (factor_only && !return_factor) {
	imsl_ermes (IMSL_TERMINAL, IMSL_SPECIFY_FACTOR_ONLY);
	return argptr;
    }
    if (solve_only) {
	if (!kbasis) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_BASIS_SOLVE_ONLY);
	    ++error;
	}
	if (!user_factor) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_SPECIFY_SOLVE_ONLY);
	    ++error;
	}
    }
    if (error)
	return argptr;

    if (solve_only) {
	Mint        ldqr = nra;
	Mint        path = 100;
	Mfloat      ab[1], qb[1], *qtb;
	qtb = (Mfloat *) imsl_malloc (nra * sizeof (Mfloat));
	if (!user_res)
	    res = (Mfloat *) imsl_malloc (nra * sizeof (Mfloat));
	if (!return_user)
	    lv_x = (Mfloat *) imsl_malloc (nca * sizeof (Mfloat));
	if (!lv_x || !qtb || !res) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
	    ++error;
	}
	else {
	    if (return_res)
		path += 10;

	    imsl_f_m1ran (nra, col_dim_fac, qr, qr);
	    if (!(error = (imsl_n1rty (1) > 3) ? 1 : 0)) {
		Mint        i;
		Mfloat      temp;
		imsl_lqrsl (&nra, kbasis, qr, &ldqr, qraux, b, &path, qb, qtb, lv_x, res, ab);
		error = (imsl_n1rty (1) > 3) ? 1 : 0;
		if (!error && user_ipvt) {
		    for (i = 0; i < *kbasis - 1; i++) {
			if (ipvt[i] > 0 && ipvt[i] <= *kbasis) {
			    temp = lv_x[ipvt[i] - 1];
			    lv_x[ipvt[i] - 1] = lv_x[i];
			    lv_x[i] = temp;
			}
		    }
		}
		imsl_free (qtb);
	    }
	}
    }
    else {
	imsl_f_m1ran (nra, col_dim_a, a, a);
	if (!(error = (imsl_n1rty (1) > 3) ? 1 : 0)) {
	    if (!user_factor) {
		qr = (Mfloat *) imsl_malloc (col_dim_fac * nra * sizeof (Mfloat));
		qraux = (Mfloat *) imsl_malloc (nca * sizeof (Mfloat));
	    }
	    else if (user_factor && col_dim_fac > nca) {
		imsl_f_m1ran (nra, col_dim_fac, qr, qr);
		error = imsl_n1rty (1) > 3 ? 1 : 0;
	    }
	    if (!user_ipvt) {
		ipvt = (Mint *) imsl_malloc (nca * sizeof (Mint));
		iset (nca, 0, ipvt, 1);
	    }
	    work = (Mfloat *) imsl_malloc ((2 * nca - 1) * sizeof (Mfloat));
	    if (!work || !ipvt || !qr || !qraux) {
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
		++error;
	    }
	    else if (factor_only) {
		Mint        lda = nra;
		Mint        ldqr = nra;
		Mint        IMSLTRUE = 1;
		if (!error) {
		    imsl_l2rrr (&nra, &nca, a, &lda, &IMSLTRUE, ipvt, qr, &ldqr, qraux, qraux, work);
		    error = (imsl_n1rty (1) > 3) ? 1 : 0;
		    if (kbasis && !error) {
			Mfloat      eps;
			Mint        i = 0;
			Mfloat     *r = qr;
			*kbasis = 0;
			eps = tol * fabs (qr[0]);
			while ((i++ < nca) && (fabs (*r) > eps)) {
			    ++(*kbasis);
			    r += nca + i;
			}
		    }
		}
		else if (kbasis) {
		    *kbasis = 0;
		}
	    }
	    else {
		if (!user_res)
		    res = (Mfloat *) imsl_malloc (nra * sizeof (Mfloat));
		if (!return_user)
		    lv_x = (Mfloat *) imsl_malloc (nca * sizeof (Mfloat));
		if (!lv_x || !res) {
		    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY);
		    ++error;
		}
		else {
		    Mint        lda = nra;
		    Mint        basis;
		    imsl_l2qrr (&nra, &nca, a, &lda, b, &tol, lv_x, res, &basis, qr, qraux, ipvt, work);
		    error = (imsl_n1rty (1) > 3) ? 1 : 0;
		    if (kbasis && !error)
			*kbasis = basis;
		}
	    }
	    imsl_f_m1ran (col_dim_a, nra, a, a);
	}
    }
    if (work)
	imsl_free (work);
    if (!user_ipvt && ipvt)
	imsl_free (ipvt);

    if (return_q && !error) {
	if (!user_q) {
	    q = (Mfloat *) imsl_malloc (nra * col_dim_q * sizeof (Mfloat));
	    if (!q) {
		++error;
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_FOR_Q);
	    }
	}
	else if (col_dim_q > nra){
	  imsl_f_m1ran (nra, col_dim_q, q, q);
	  error = (imsl_n1rty(1) > 3)? 1: 0;
	}
	if (!error) {
	    work = (Mfloat *) imsl_malloc (2 * nra * sizeof (Mfloat));
	    if (!work) {
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_FOR_Q);
	    }
	    else {
		Mint        ldqr = nra;
		Mint        ldq = nra;
		l_l2err (&nra, &nca, qr, &ldqr, qraux, q, &ldq, work);
		error = (imsl_n1rty (1) > 3) ? 1 : 0;
		if (!user_q) {
		    if (error) {
			imsl_free (q);
		    }
		    else {
			imsl_f_m1ran (col_dim_q, nra, q, q);
			for (i = nra; i < col_dim_q; i++) {
			  sset (nra, 0.0, (q+i), col_dim_q);
			}
			*p_q = q;
		    }
		    q = NULL;
		}
		else if (!error)
		    imsl_f_m1ran (col_dim_q, nra, q, q);

		imsl_free (work);
	    }
	}
    }
    if (!user_res) {
	if (error || !return_res) {
	    if (res)
		imsl_free (res);
	}
	else if (return_res)
	    *p_res = res;
	res = NULL;
    }

    if (!user_factor) {
	if (error || !return_factor) {
	    if (qraux)
		imsl_free (qraux);
	    if (qr)
		imsl_free (qr);
	}
	else if (return_factor) {
	    imsl_f_m1ran (col_dim_fac, nra, qr, qr);
	    for (i = nca; i < col_dim_fac; i++){
	      sset(nra, 0.0, (qr+i), col_dim_fac);
	    }
	    *p_qr = qr;
	    *p_qraux = qraux;
	}
	qraux = NULL;
	qr = NULL;
    }
    else if (return_factor)
	imsl_f_m1ran (col_dim_fac, nra, qr, qr);

    if (error) {
	if (!return_user) {
	    if (lv_x)
		imsl_free (lv_x);
	}
	lv_x = NULL;
    }

    return argptr;
}

/*
  -----------------------------------------------------------------------
    IMSL Name:  L2QRR/DL2QRR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Solve a linear least squares problem without iterative
                refinement.

    Usage:      CALL L2QRR (NRA, NCA, A, LDA, B, TOL, X, RES, KBASIS,
                            QR, QRAUX, IPVT, WORK)

    Arguments:  See LSQRR.

    Remarks:    See LSQRR.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_l2qrr (Mint *nra, Mint *nca, Mfloat *a, Mint *lda, Mfloat *b,
                Mfloat *tol, Mfloat *x, Mfloat *res, Mint *kbasis,
                Mfloat *qr, Mfloat *qraux, Mint *ipvt, Mfloat *work)
#else
void        imsl_l2qrr (nra, nca, a, lda, b, tol, x, res, kbasis,
                qr, qraux, ipvt, work)
    Mint       *nra, *nca;
    Mfloat     *a;
    Mint       *lda;
    Mfloat      b[], *tol, x[], res[];
    Mint       *kbasis;
    Mfloat      qr[], qraux[];
    Mint        ipvt[];
    Mfloat      work[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
    Mint        _l0, kk;
    Mfloat      dum[1];
    Mint        IMSLTRUE = 1;
    imsl_e1psh ("imsl_l2qrr");

    if (*nra <= 0 || *nca <= 0) {
	imsl_e1sti (1, *nra);
	imsl_e1sti (2, *nca);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NRA_AND_NCA_GT_ZERO);
	goto L_9000;
    }
    if (*nra > *lda) {
	imsl_e1sti (1, *nra);
	imsl_e1sti (2, *lda);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_SMALLER_NRA_VALUE);
	goto L_9000;
    }
    /*
     * Initialize IPVT so that all columns are free
     * 
     * iset(*nca, 0, ipvt, 1); (now initialized in lin_least_square_gen)
     * 
     * QR decomposition of A with column pivoting
     */
    imsl_l2rrr (nra, nca, a, lda, &IMSLTRUE, ipvt, qr, nra, qraux,
	qraux, work);
    /* Determine which columns to use */
    *kbasis = 0;
    for (kk = 1; kk <= imsl_i_min (*nra, *nca); kk++) {
	if (fabs (qr[(kk - 1) ** nra + kk - 1]) <= *tol * fabs (qr[0]))
	    goto L_20;
	*kbasis = kk;
    }
    /*
     * Solve the truncated least squares problem
     */
L_20:
    if (*kbasis != 0) {
	_l0 = 110;
	imsl_lqrsl (nra, kbasis, qr, nra, qraux, b, &_l0, dum,
	    res, x, res, dum);
	if (imsl_n1rcd (1) != 0)
	    goto L_9000;
    }
    /*
     * Set unused components of solution to zero
     */
    if (*kbasis < *nca)
	sset (*nca - *kbasis, F_ZERO, &x[*kbasis], 1);
    /*
     * Unscramble the solution, i.e. apply the permutation
     */
    _l0 = 2;
    l_permu (nca, x, ipvt, &_l0, x);

L_9000:
    imsl_e1pop ("imsl_l2qrr");
    return;
}				/* end of function */
#undef A
/*
  -----------------------------------------------------------------------
    IMSL Name:  L2ERR/DL2ERR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Accumulate the orthogonal matrix Q from its factored form
                given the QR factorization of a rectangular matrix A.

    Usage:      CALL L2ERR (NRQR, NCQR, QR, LDQR, QRAUX, Q, LDQ, WORK)

    Arguments:  See LQERR.

    Remarks:    See LQERR.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_l2err (Mint *nrqr, Mint *ncqr, Mfloat *qr, Mint *ldqr,
                Mfloat *qraux, Mfloat *q, Mint *ldq, Mfloat *work)
#else
static void l_l2err (nrqr, ncqr, qr, ldqr, qraux, q, ldq, work)
    Mint       *nrqr, *ncqr;
    Mfloat     *qr;
    Mint       *ldqr;
    Mfloat      qraux[], *q;
    Mint       *ldq;
    Mfloat      work[];
#endif
{
#define QR(I_,J_)	(qr+(I_)*(*ldqr)+(J_))
#define Q(I_,J_)	(q+(I_)*(*ldq)+(J_))
    Mint        _l0, _l1, j, k, minmn, mmk, mmkp1;
    Mfloat      _f0, _f1, big, small;


    imsl_e1psh ("l_l2err");

    if (*nrqr <= 0 || *ncqr <= 0) {
	imsl_e1sti (1, *nrqr);
	imsl_e1sti (2, *ncqr);

	/*
	 * (5, 1, "Both the number of rows and the number of columns of the
	 * input matrix have to be positive while NRQR = %(i1) and NCQR =
	 * %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_POS_COL_ROWS);
    }
    else if (*nrqr > *ldqr) {
	imsl_e1sti (1, *nrqr);
	imsl_e1sti (2, *ldqr);

	/*
	 * (5, 2, "The number of rows of QR must be less than or equal to its
	 * leading dimension while NRQR = %(i1) and LDQR = %(i2) are
	 * given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NRQR_GT_LDQR);
    }
    else if (*nrqr > *ldq) {
	imsl_e1sti (1, *nrqr);
	imsl_e1sti (2, *ldq);

	/*
	 * (5, 2, "The number of rows of Q must be less than or equal to its
	 * leading dimension while NRQR = %(i1) and LDQ = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_ROWS_LE_LEAD_DIM);
    }
    else {
	/* COPY RELEVANT PORTION OF QR INTO Q */
	minmn = imsl_i_min (*nrqr, *ncqr);
	for (j = 1; j <= minmn; j++) {
	    scopy (*nrqr, QR (j - 1, 0), 1, Q (j - 1, 0), 1);
	}
	/*
	 * ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(NRQR,NCQR) COLUMNS.
	 */
	for (j = 1; j <= minmn; j++) {
	    sset (j, F_ZERO, Q (j - 1, 0), 1);
	}
	/*
	 * INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.
	 */
	for (j = minmn + 1; j <= *nrqr; j++) {
	    sset (*nrqr, F_ZERO, Q (j - 1, 0), 1);
	    *Q (j - 1, j - 1) = F_ONE;
	}
	small = imsl_amach (1);
	big = imsl_amach (2);
	if (small * big < F_ONE)
	    small = F_ONE / big;
	/* ACCUMULATE Q FROM ITS FACTORED FORM. */
	for (k = minmn; k >= 1; k--) {
	    mmk = *nrqr - k;
	    work[k - 1] = qraux[k - 1];
	    if (mmk != 0)
		scopy (mmk, Q (k - 1, k), 1, &work[k], 1);
	    *Q (k - 1, k - 1) = F_ONE;
	    if (mmk != 0)
		sset (mmk, F_ZERO, Q (k - 1, k), 1);
	    if (fabs (work[k - 1]) >= small) {
		mmkp1 = *nrqr - k + 1;
		_f0 = F_ONE / work[k - 1];
		_l0 = 1;
		_f1 = F_ZERO;
		_l1 = 1;
		imsl_sgemv ("T", sizeof ("T"), &mmkp1, &mmkp1, &_f0, Q (k - 1, k - 1), ldq, &work[k - 1],
		    &_l0, &_f1, &work[*nrqr], &_l1);
		imsl_sger (mmkp1, mmkp1, -F_ONE, &work[k - 1], 1, &work[*nrqr],
		    1, Q (k - 1, k - 1), *ldq);
	    }
	}
    }

    imsl_e1pop ("l_l2err");
    return;
}				/* end of function */
#undef Q
#undef QR

#if 0 /* old l2rrr */
/*
  -----------------------------------------------------------------------
    IMSL Name:  L2RRR/DL2RRR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 31, 1989

    Purpose:    Compute the QR decomposition using Householder
                transformations.

    Usage:      CALL L2RRR (NRA, NCA, A, LDA, PIVOT, IPVT, QR, LDQR,
                            QRAUX, CONORM, WORK)

    Arguments:  See LQRRR/DLQRRR.

    Remarks:    See LQRRR/DLQRRR.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_l2rrr (Mint *nra, Mint *nca, Mfloat *a, Mint *lda, Mint *pivot,
                Mint *ipvt, Mfloat *qr, Mint *ldqr, Mfloat *qraux,
                Mfloat *conorm, Mfloat *work)
#else
void        imsl_l2rrr (nra, nca, a, lda, pivot, ipvt, qr, ldqr, qraux, conorm, work)
    Mint       *nra, *nca;
    Mfloat     *a;
    Mint       *lda;
    Mint       *pivot;
    Mint        ipvt[];
    Mfloat     *qr;
    Mint       *ldqr;
    Mfloat      qraux[], conorm[], work[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define QR(I_,J_)	(qr+(I_)*(*ldqr)+(J_))
    Mint        negj, swapj;
    Mint        _l0, _l1, _l2, _l3, ipl, ipu, j, jp, l, lup, maxj;
    Mfloat      _f0, _f1, anorm, t, tt, xlnrm;


    imsl_e1psh ("imsl_l2rrr");

    if (*nra <= 0 || *nca <= 0) {
	imsl_e1sti (1, *nra);
	imsl_e1sti (2, *nca);

	/*
	 * (5, 1, "Both the number of rows and the number of columns of the
	 * input matrix must be positive while NRA = %(i1) and NCA = %(i2)
	 * are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NRA_AND_NCA_GT_ZERO);

    }
    else if (*nra > *lda) {
	imsl_e1sti (1, *nra);
	imsl_e1sti (2, *lda);

	/*
	 * (5, 2, "The number of rows of A must be less than or equal to its
	 * leading dimension while NRA = %(i1) and LDA = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_SMALLER_NRA_VALUE);

    }
    else if (*nra > *ldqr) {
	imsl_e1sti (1, *nra);
	imsl_e1sti (2, *ldqr);

	/*
	 * (5, 3, "The number of rows of QR must be less than or equal to its
	 * leading dimension while NRA = %(i1) and LDQR = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NRQR_GREATER_THAN_LDQR);
    }
    if (imsl_n1rcd (0) != 0)
	goto L_9000;
    /*
     * MAKE A COPY OF A IN QR AND WORK WITH QR ONLY
     */
    for (j = 1; j <= *nca; j++) {
	scopy (*nra, A (j - 1, 0), 1, QR (j - 1, 0), 1);
    }

    ipl = 1;
    ipu = 0;
    /*
     * PIVOTING HAS BEEN REQUESTED. REARRANGE COLUMNS ACCORDING TO IPVT.
     */
    if (*pivot) {
	for (j = 1; j <= *nca; j++) {
	    swapj = ipvt[j - 1] > 0;
	    negj = ipvt[j - 1] < 0;
	    ipvt[j - 1] = j;
	    if (negj)
		ipvt[j - 1] = -j;

	    if (swapj) {
		if (j != ipl)
		    sswap (*nra, QR (ipl - 1, 0), 1, QR (j - 1, 0), 1);
		ipvt[j - 1] = ipvt[ipl - 1];
		ipvt[ipl - 1] = j;
		ipl += 1;
	    }
	}

	ipu = *nca;
	for (j = *nca; j >= 1; j--) {
	    if (ipvt[j - 1] < 0) {
		ipvt[j - 1] = -ipvt[j - 1];

		if (j != ipu) {
		    sswap (*nra, QR (ipu - 1, 0), 1, QR (j - 1, 0), 1);
		    jp = ipvt[ipu - 1];
		    ipvt[ipu - 1] = ipvt[j - 1];
		    ipvt[j - 1] = jp;
		}
		ipu -= 1;
	    }
	}
    }
    /*
     * COMPUTE NORMS OF THE COLUMNS OF A AND INITIALIZE THE WORK ARRAY.
     */
    for (j = 1; j <= *nca; j++) {
	conorm[j - 1] = imsl_snrm2 (*nra, QR (j - 1, 0), 1);
    }
    scopy (ipu - ipl + 1, &conorm[ipl - 1], 1, &qraux[ipl - 1], 1);
    if (*pivot)
	scopy (ipu - ipl + 1, &qraux[ipl - 1], 1, &work[ipl + *nca - 2],
	    1);
    /* PERFORM HOUSEHOLDER REDUCTION OF X. */
    lup = imsl_i_min (*nra, *nca);
    for (l = 1; l <= lup; l++) {
	/*
	 * LOCATE COLUMN OF LARGEST NORM AND BRING IT INTO PIVOT POSITION.
	 */
	if (l >= ipl && l < ipu) {
	    _l0 = ipu - l + 1;
	    _l1 = 1;
	    maxj = l_ismax (&_l0, &qraux[l - 1], &_l1) +
		l - 1;
	    anorm = qraux[maxj - 1];

	    if (maxj != l) {
		sswap (*nra, QR (l - 1, 0), 1, QR (maxj - 1, 0), 1);
		qraux[maxj - 1] = qraux[l - 1];
		if (*pivot)
		    work[maxj + *nca - 2] = work[l + *nca - 2];
		jp = ipvt[maxj - 1];
		ipvt[maxj - 1] = ipvt[l - 1];
		ipvt[l - 1] = jp;
	    }
	}
	qraux[l - 1] = F_ZERO;
	/*
	 * COMPUTE HOUSEHOLDER TRANSFORMATION COLUMN L.
	 */
	if (l != *nra) {
	    xlnrm = imsl_snrm2 (*nra - l + 1, QR (l - 1, l - 1), 1);
	    if (xlnrm != F_ZERO) {
		if (*QR (l - 1, l - 1) < F_ZERO)
		    xlnrm = -xlnrm;
		sscal (*nra - l + 1, F_ONE / xlnrm, QR (l - 1, l - 1),
		    1);
		*QR (l - 1, l - 1) += F_ONE;
		/*
		 * APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
		 * UPDATING THE NORMS
		 */
		if (l < *nca) {
		    _l0 = *nra - l + 1;
		    _l1 = *nca - l;
		    _f0 = -F_ONE / (*QR (l - 1, l - 1));
		    _l2 = 1;
		    _f1 = F_ZERO;
		    _l3 = 1;
		    imsl_sgemv ("T", sizeof ("T"), &_l0,
			&_l1, &_f0, QR (l, l - 1), ldqr, QR (l - 1, l - 1),
			&_l2, &_f1, work, &_l3);
		    imsl_sger (*nra - l + 1, *nca - l, F_ONE, QR (l - 1, l - 1),
			1, work, 1, QR (l, l - 1), *ldqr);
		}
		for (j = l + 1; j <= *nca; j++) {
		    if (j >= ipl && j <= ipu) {
			if (qraux[j - 1] != F_ZERO) {
			    tt = F_ONE - imsl_fi_power (fabs (*QR (j - 1, l - 1)) /
				qraux[j - 1], 2);
			    tt = imsl_f_max (tt, F_ZERO);
			    t = tt;
			    if (*pivot) {
				tt = F_ONE + 5.0e-2 * tt * imsl_fi_power (qraux[j - 1] /
				    work[j + *nca - 2], 2);
			    }
			    else {
				tt = F_ONE + 5.0e-2 * tt;
			    }
			    if (tt != F_ONE) {
				qraux[j - 1] *= sqrt (t);
			    }
			    else {
				qraux[j - 1] = imsl_snrm2 (*nra - l, QR (j - 1, l),
				    1);
				if (*pivot)
				    work[j + *nca - 2] = qraux[j - 1];
			    }
			}
		    }
		}
		/* SAVE THE TRANSFORMATION. */
		qraux[l - 1] = *QR (l - 1, l - 1);
		*QR (l - 1, l - 1) = -xlnrm;
	    }
	}
    }

L_9000:
    imsl_e1pop ("imsl_l2rrr");
    return;
}				/* end of function */
#undef  QR
#undef  A

#endif /*old l2rrr */

/*Translated by FOR_C++, v0.1, on 11/29/91 at 12:54:13 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 11/29/91 at 12:54:11
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  L2RRR/DL2RRR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 24, 1991

    Purpose:    Compute the QR decomposition using Householder
                transformations.

    Usage:      CALL L2RRR (NRA, NCA, A, LDA, PIVOT, IPVT, QR, LDQR,
                            QRAUX, CONORM, WORK)

    Arguments:  See LQRRR/DLQRRR.

    Remarks:    See LQRRR/DLQRRR.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.
 
  -----------------------------------------------------------------------
 */
void  imsl_l2rrr (nra, nca, a, lda, pivot, ipvt, qr, ldqr, qraux,
                conorm, work)
    Mint        *nra, *nca;
    Mfloat      *a;
    Mint        *lda;
    Mint        *pivot;
    Mint         ipvt[];
    Mfloat      *qr;
    Mint        *ldqr;
    Mfloat       qraux[], conorm[], work[];
{
#define A(I_,J_)        (a+(I_)*(*lda)+(J_))
#define QR(I_,J_)       (qr+(I_)*(*ldqr)+(J_))
    Mint   negj, swapj;
    Mint         _l0, _l1, _l2, _l3, ipl, ipu, j, jp, l,
                lup, maxj;
    Mfloat       _f0, anorm, big, one, small, t, tt, xlnrm, zero;


    imsl_e1psh ("L2RRR ");
    zero = 0.0e0;
    one = 1.0e0;
    small = imsl_amach (1);
    big = imsl_amach (2);
    if (small * big < one)
        small = one / big;

    if (*nra <= 0 || *nca <= 0) {
        imsl_e1sti (1, *nra);
        imsl_e1sti (2, *nca);
/*
        (5, 1, "Both the number of rows and the number of columns of 
	the input matrix must be positive while NRA = %(i1) and NCA = 
	%(i2) are given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NRA_AND_NCA_GT_ZERO);
    }
    else if (*nra > *lda) {
        imsl_e1sti (1, *nra);
        imsl_e1sti (2, *lda);
/*
        (5, 2, "The number of rows of A must be less than or equal to 
	its leading dimension while NRA = %(i1) and LDA = %(i2) are given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_SMALLER_NRA_VALUE);
    }
    else if (*nra > *ldqr) {
        imsl_e1sti (1, *nra);
        imsl_e1sti (2, *ldqr);
/*
        (5, 3, "The number of rows of QR must be less than or equal to 
	its leading dimension while NRA = %(i1) and LDQR = %(i2) are given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NRQR_GREATER_THAN_LDQR);
    }
    if (imsl_n1rcd (0) != 0)
        goto L_9000;
    /*
     * MAKE A COPY OF A IN QR AND WORK WITH QR ONLY
     */
    for (j = 1; j <= *nca; j++) {
        scopy (*nra, A (j - 1, 0), 1, QR (j - 1, 0), 1);
    }

    ipl = 1;
    ipu = 0;
    /*
     * PIVOTING HAS BEEN REQUESTED. REARRANGE COLUMNS ACCORDING TO IPVT.
     */
    if (*pivot) {
        for (j = 1; j <= *nca; j++) {
            swapj = ipvt[j - 1] > 0;
            negj = ipvt[j - 1] < 0;
            ipvt[j - 1] = j;
            if (negj)
                ipvt[j - 1] = -j;

            if (swapj) {
                if (j != ipl)
                    sswap (*nra, QR (ipl - 1, 0), 1, QR (j - 1, 0), 1);
                ipvt[j - 1] = ipvt[ipl - 1];
                ipvt[ipl - 1] = j;
                ipl += 1;
            }
        }

        ipu = *nca;
        for (j = *nca; j >= 1; j--) {
            if (ipvt[j - 1] < 0) {
                ipvt[j - 1] = -ipvt[j - 1];

                if (j != ipu) {
                    sswap (*nra, QR (ipu - 1, 0), 1, QR (j - 1, 0), 1);
                    jp = ipvt[ipu - 1];
                    ipvt[ipu - 1] = ipvt[j - 1];
                    ipvt[j - 1] = jp;
                }

                ipu -= 1;
            }
        }
    }
    /*
     * COMPUTE NORMS OF THE COLUMNS OF A AND INITIALIZE THE WORK ARRAY.
     */
    for (j = 1; j <= *nca; j++) {
        conorm[j - 1] = imsl_snrm2 (*nra, QR (j - 1, 0), 1);
    }
    scopy (ipu - ipl + 1, &conorm[ipl - 1], 1, &qraux[ipl - 1], 1);
    if (*pivot)
        scopy (ipu - ipl + 1, &qraux[ipl - 1], 1, &work[ipl + *nca - 2],
            1);
    /* PERFORM HOUSEHOLDER REDUCTION OF X. */
    lup = imsl_i_min (*nra, *nca);
    for (l = 1; l <= lup; l++) {
        /*
         * LOCATE COLUMN OF LARGEST NORM AND BRING IT INTO PIVOT POSITION.
         */
        if (l >= ipl && l < ipu) {
            _l0 = ipu - l + 1;
            _l1 = 1;
            maxj = l_ismax (&_l0, &qraux[l - 1], &_l1) +
                l - 1;
            anorm = qraux[maxj - 1];

            if (maxj != l) {
                sswap (*nra, QR (l - 1, 0), 1, QR (maxj - 1, 0), 1);
                qraux[maxj - 1] = qraux[l - 1];
                if (*pivot)
                    work[maxj + *nca - 2] = work[l + *nca - 2];
                jp = ipvt[maxj - 1];
                ipvt[maxj - 1] = ipvt[l - 1];
                ipvt[l - 1] = jp;
            }
        }
        qraux[l - 1] = zero;
        /*
         * COMPUTE HOUSEHOLDER TRANSFORMATION COLUMN L.
         */
        if (l != *nra) {
            xlnrm = imsl_snrm2 (*nra - l + 1, QR (l - 1, l - 1), 1);
            if (xlnrm >= small) {
                if (*QR (l - 1, l - 1) < zero)
                    xlnrm = -xlnrm;
                sscal (*nra - l + 1, one / xlnrm, QR (l - 1, l - 1), 1);
                *QR (l - 1, l - 1) += one;
                /*
                 * APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
                 * UPDATING THE NORMS
                 */
                if (l < *nca) {
                    _l0 = *nra - l + 1;
                    _l1 = *nca - l;
                    _f0 = -one / *QR (l - 1, l - 1);
                    _l2 = 1;
                    _l3 = 1;
                    imsl_sgemv ("T", sizeof ("T"), &_l0,
                        &_l1, &_f0,
                        QR (l, l - 1), ldqr, QR (l - 1, l - 1), &_l2,
                        &zero, work, &_l3);
                    imsl_sger (*nra - l + 1, *nca - l, one, QR (l - 1, l - 1),
                        1, work, 1, QR (l, l - 1), *ldqr);
                }

                for (j = l + 1; j <= *nca; j++) {
                    if (j >= ipl && j <= ipu) {
                        if (qraux[j - 1] != zero) {
                            tt = one - imsl_fi_power (fabs (*QR (j - 1, l - 1)) /
                                qraux[j - 1], 2);
                            tt = imsl_f_max (tt, zero);
                            t = tt;
                            if (*pivot) {
                                tt = one + 5.0e-2 * tt * imsl_fi_power (qraux[j - 1] /
                                    work[j + *nca - 2], 2);
                            }
                            else {
                                tt = one + 5.0e-2 * tt;
                            }
                            if (tt != one) {
                                qraux[j - 1] *= sqrt (t);
                            }
                            else {
                                qraux[j - 1] = imsl_snrm2 (*nra - l, QR (j - 1, l),
                                    1);
                                if (*pivot)
                                    work[j + *nca - 2] = qraux[j - 1];
                            }
                        }

                    }
                }
                /* SAVE THE TRANSFORMATION. */
                qraux[l - 1] = *QR (l - 1, l - 1);
                *QR (l - 1, l - 1) = -xlnrm;
            }
        }
    }

L_9000:
    ;
    imsl_e1pop ("L2RRR ");
    return;
}                               /* end of function */

#undef A
#undef QR


/*
  -----------------------------------------------------------------------
    IMSL Name:  LQRSL/DLQRSL (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the coordinate transformation, projection, and
                solution for the least squares problem.

    Usage:      CALL LQRSL (NRA, KBASIS, QR, LDQR, QRAUX, B, IPATH, QB,
                            QTB, X, RES, AB)

    Arguments:
       NRA    - Number of rows of A.  (Input)
       KBASIS - Number of columns of the submatrix of AK of A.  (Input)
                KBASIS must be less than or equal to MIN(NRA,NCA), where
                NCA is the number of columns in A.  See Remarks.
       QR     - NRA by NCA matrix containing information about the QR
                factorization of A as output from routine LQRRR/DLQRRR.
                (Input)
                NCA is the number of columns in A.
       LDQR   - Leading dimension of QR exactly as specified in the
                dimension statement of the calling program.  (Input)
       QRAUX  - Vector of length NCA containing information about the QR
                factorization of A as output from routine LQRRR/DLQRRR.
                (Input)
                NCA is the number of columns in A.
       B      - Vector of length NRA to be manipulated.  (Input)
       IPATH  - Option parameter specifying what is to be computed.
                (Input)
                IPATH has the decimal expansion IJKLM, such that
                   I .NE. 0 means compute Q*B,
                   J .NE. 0 means compute trans(Q)*B,
                   K .NE. 0 means compute trans(Q)*B and X,
                   L .NE. 0 means compute trans(Q)*B and RES,
                   M .NE. 0 means compute trans(Q)*B and AB.
                For example, if IPATH = 01101, then I = 0, J = 1, K = 1,
                L = 0 and M=1.
       QB     - Vector of length NRA containing Q*B if requested in
                IPATH.  (Output)
       QTB    - Vector of length NRA containing trans(Q)*B if requested
                in IPATH.  (Output)
       X      - Vector of length KBASIS containing the solution of the
                least squares problem AK*X = B if this is requested if
                IPATH.  (Output)
                If pivoting was requested in LQRRR/DLQRRR, then the J-th
                component of X will be associated with column IPVT(J) of
                the original matrix A.  See Remarks.
       RES    - Vector of length NRA containing the residuals of the
                least squares problem if requested in IPATH.  (Output)
                This is also the orthogonal projection of B onto the
                orthogonal complement of the column space of A.
       AB     - Vector of length NRA containing the least squares
                approximation A*X if requested in IPATH.  (Output)
                This is also the orthogonal projection of B onto the
                column space of A.

    Remarks:
    1. Informational error
       Type Code
         4   1  Computation of the least squares solution of AK*X = B
                is requested, but the upper triangular matrix R from the
                QR factorization is singular.

    2. This routine is designed to be used together with LQRRR.  It
       assumes that LQRRR/DLQRR has been called to get QR, QRAUX and
       IPVT.  The submatrix AK mentioned above is actually equal to

                AK = (A(IPVT(1)), A(IPVT(2)), ..., A(IPVT(KBASIS))),

       where A(IPVT(I)) is the IPVT(I)-th column of the original matrix.

    Keywords:   QR factorization; Overdetermined system; Underdetermined
                system

    GAMS:       D5; D9

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void        imsl_lqrsl (Mint *nra, Mint *kbasis, Mfloat *qr, Mint *ldqr, Mfloat *qraux,
                Mfloat *b, Mint *ipath, Mfloat *qb, Mfloat *qtb, Mfloat *x,
                Mfloat *res, Mfloat *ab)
#else
void        imsl_lqrsl (nra, kbasis, qr, ldqr, qraux, b, ipath, qb,
                qtb, x, res, ab)
    Mint       *nra, *kbasis;
    Mfloat     *qr;
    Mint       *ldqr;
    Mfloat      qraux[], b[];
    Mint       *ipath;
    Mfloat      qb[], qtb[], x[], res[], ab[];
#endif
{
#define QR(I_,J_)	(qr+(I_)*(*ldqr)+(J_))
    Mint        cax, cqb, cqtb, cres, csol;
    Mint        info, j, ju;
    Mfloat      t, temp;


    imsl_e1psh ("imsl_lqrsl");

    if (*nra <= 0) {
	imsl_e1sti (1, *nra);

	/*
	 * (5, 1, "The number of rows of the input matrix must be positive
	 * while NRA = %(i1) is given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NRA_MUST_BE_POSITIVE);
    }
    else if (*nra > *ldqr) {
	imsl_e1sti (1, *nra);
	imsl_e1sti (2, *ldqr);

	/*
	 * (5, 2, "The number of rows of QR must be less than or equal to its
	 * leading dimension while NRA = %(i1) and LDQR = %(i2) are given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_NRQR_GREATER_THAN_LDQR);
    }
    else if (*kbasis <= 0) {
	imsl_e1sti (1, *kbasis);

	/*
	 * (5, 3, "The number of columns of the submatrix AK of A must be
	 * positive while KBASIS = %(i1) is given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_KBASIS_IS_NEGATIVE);
    }
    if (*ipath <= 0 || *ipath > 99999) {
	imsl_e1sti (1, *ipath);

	/*
	 * (5, 4, "Illegal format for IPATH.  IPATH must be between 1 and
	 * 99999 while a value of %(i1) is given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_KBASIS_IS_NEGATIVE);
    }
    if (imsl_n1rcd (0) == 0) {
	/* DETERMINE WHAT IS TO BE COMPUTED */
	info = 0;
	cqb = *ipath / 10000 != 0;
	cqtb = mod (*ipath, 10000) != 0;
	csol = mod (*ipath, 1000) / 100 != 0;
	cres = mod (*ipath, 100) / 10 != 0;
	cax = mod (*ipath, 10) != 0;
	ju = imsl_i_min (*kbasis, *nra - 1);
	info = 0;
	/* SPECIAL CASE WHEN NRA = 1 */
	if (ju == 0) {
	    if (cqb)
		qb[0] = b[0];
	    if (cqtb)
		qtb[0] = b[0];
	    if (cres)
		res[0] = F_ZERO;
	    if (cax)
		ab[0] = b[0];
	    if (csol) {
		if (*QR (0, 0) == F_ZERO) {
		    info = 1;
		}
		else {
		    x[0] = b[0] / *QR (0, 0);
		}
	    }
	    goto L_50;
	}
	/*
	 * --- GENERAL CASE--- SET UP TO COMPUTE QB OR QTB
	 */
	if (cqb)
	    scopy (*nra, b, 1, qb, 1);
	if (cqtb)
	    scopy (*nra, b, 1, qtb, 1);
	/* COMPUTE QB */
	if (cqb) {
	    for (j = ju; j >= 1; j--) {
		if (qraux[j - 1] != F_ZERO) {
		    temp = *QR (j - 1, j - 1);
		    *QR (j - 1, j - 1) = qraux[j - 1];
		    t = imsl_sdot (*nra - j + 1, QR (j - 1, j - 1), 1, &qb[j - 1],
			1);
		    t = -F_ONE * t / *QR (j - 1, j - 1);
		    saxpy (*nra - j + 1, t, QR (j - 1, j - 1), 1, &qb[j - 1],
			1);
		    *QR (j - 1, j - 1) = temp;
		}
	    }
	}
	/* COMPUTE QTB */
	if (cqtb) {
	    for (j = 1; j <= ju; j++) {
		if (qraux[j - 1] != F_ZERO) {
		    temp = *QR (j - 1, j - 1);
		    *QR (j - 1, j - 1) = qraux[j - 1];
		    t = imsl_sdot (*nra - j + 1, QR (j - 1, j - 1), 1, &qtb[j - 1],
			1);
		    t = -F_ONE * t / *QR (j - 1, j - 1);
		    saxpy (*nra - j + 1, t, QR (j - 1, j - 1), 1, &qtb[j - 1],
			1);
		    *QR (j - 1, j - 1) = temp;
		}
	    }
	}
	/* SET UP TO COMPUTE X, RES OR AB */
	if (csol)
	    scopy (*kbasis, qtb, 1, x, 1);
	if (cax)
	    scopy (*kbasis, qtb, 1, ab, 1);
	if (cres && *kbasis < *nra)
	    scopy (*nra - *kbasis, &qtb[*kbasis], 1, &res[*kbasis],
		1);
	if (cax && (*kbasis + 1 <= *nra))
	    sset (*nra - *kbasis, F_ZERO, &ab[*kbasis], 1);
	if (cres)
	    sset (*kbasis, F_ZERO, res, 1);
	/* COMPUTE X */
	if (csol) {
	    for (j = *kbasis; j >= 1; j--) {
		if (*QR (j - 1, j - 1) == F_ZERO) {
		    info = j;
		    goto L_50;
		}
		else {
		    x[j - 1] /= *QR (j - 1, j - 1);
		    t = -x[j - 1];
		    saxpy (j - 1, t, QR (j - 1, 0), 1, x, 1);
		}
	    }
	}
	/* COMPUTE RES OR AB AS REQUIRED */
	if (cres || cax) {
	    for (j = ju; j >= 1; j--) {
		if (qraux[j - 1] != F_ZERO) {
		    temp = *QR (j - 1, j - 1);
		    *QR (j - 1, j - 1) = qraux[j - 1];
		    if (cres) {
			t = imsl_sdot (*nra - j + 1, QR (j - 1, j - 1), 1,
			    &res[j - 1], 1);
			t = -F_ONE * t / *QR (j - 1, j - 1);
			saxpy (*nra - j + 1, t, QR (j - 1, j - 1), 1,
			    &res[j - 1], 1);
		    }
		    if (cax) {
			t = imsl_sdot (*nra - j + 1, QR (j - 1, j - 1), 1,
			    &ab[j - 1], 1);
			t = -F_ONE * t / *QR (j - 1, j - 1);
			saxpy (*nra - j + 1, t, QR (j - 1, j - 1), 1,
			    &ab[j - 1], 1);
		    }
		    *QR (j - 1, j - 1) = temp;
		}
	    }
	}
L_50:
	if (info != 0) {
	    imsl_e1sti (1, info);

	    /*
	     * (4, 1, "Computation of the least squares solution has been
	     * requested but the upper triangular matrix R is exactly
	     * singular; the index of the first zero diagonal element of R is
	     * %(i1).");
	     */
	    imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_TRI_MATRIX);
	}
    }
    imsl_e1pop ("imsl_lqrsl");
    return;
}				/* end of function */
#undef  QR
/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 16:10:55
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  ISMAX (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Find the smallest index of the component of a
                single-precision vector having maximum value.

    Usage:      ISMAX(N, SX, INCX)

    Arguments:
       N      - Length of vector X.  (Input)
       SX     - Real vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
                greater than zero.
       ISMAX  - The smallest index I such that X(I)  is the maximum of
                X(J) for J=1 to N.  (Output)
                X(I) refers to a specific element of SX. See INCX
                argument description.

    Keyword:    Level 1 BLAS

    GAMS:       D1a2

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mint l_ismax (Mint *n, Mfloat *sx, Mint *incx)
#else
static Mint l_ismax (n, sx, incx)
    Mint       *n;
    Mfloat      sx[];
    Mint       *incx;
#endif
{
    Mint        i, ismax_v, ix;
    Mfloat      smax;


    ismax_v = 0;
    if (*n >= 1) {
	ismax_v = 1;
	if (*n != 1) {
	    if (*incx != 1) {
		/* CODE FOR INCREMENT NOT EQUAL TO 1 */
		ix = 1;
		smax = sx[0];
		ix += *incx;
		for (i = 2; i <= *n; i++) {
		    if (sx[ix - 1] > smax) {
			ismax_v = i;
			smax = sx[ix - 1];
		    }
		    ix += *incx;
		}
	    }
	    else {
		/* CODE FOR INCREMENT EQUAL TO 1 */
		smax = sx[0];
		for (i = 2; i <= *n; i++) {
		    if (sx[i - 1] > smax) {
			ismax_v = i;
			smax = sx[i - 1];
		    }
		}
	    }
	}
    }
    return (ismax_v);
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 16:08:50
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  PERMU/DPERMU (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 30, 1985

    Purpose:    Rearrange the elements of an array as specified by a
                permutation.

    Usage:      CALL PERMU (N, X, IPERMU, IPATH, XPERMU)

    Arguments:
       N      - Length of the arrays X and XPERMU.  (Input)
       X      - Real vector of length N containing the array to be
                permuted.  (Input)
       IPERMU - Integer vector of length N containing a permutation
                IPERMU(1), ..., IPERMU(N) of the integers 1, ..., N.
                (Input)
       IPATH  - Integer flag.  (Input)
                IPATH = 1 means IPERMU represents a forward permutation,
                          i.e., X(IPERMU(I)) is moved to XPERMU(I).
                IPATH = 2 means IPERMU represents a backward permutation,
                          i.e., X(I) is moved to XPERMU(IPERMU(I)).
       XPERMU - Real vector of length N containing the array X permuted.
                (Output)
                If X is not needed, X and XPERMU can share the same
                storage locations.

    Keywords:   Utilities; Forward permutation; Backward permutation

    GAMS:       N8

    Chapters:   MATH/LIBRARY Utilities
                STAT/LIBRARY Utilities

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_permu (Mint *n, Mfloat *x, Mint *ipermu, Mint *ipath, Mfloat *xpermu)
#else
static void l_permu (n, x, ipermu, ipath, xpermu)
    Mint       *n;
    Mfloat      x[];
    Mint        ipermu[], *ipath;
    Mfloat      xpermu[];
#endif
{
    Mint        i, j, k;
    Mfloat      temp;


    imsl_e1psh ("l_permu");

    if (*n <= 0) {
	imsl_e1sti (1, *n);

	/*
	 * (5, 1, "The length of the arrays X and XPERMU must be positive
	 * while N = %(i1) is given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_X_AND_XPERMU_LENGTH_LE_0);
    }
    if (*ipath != 1 && *ipath != 2) {
	imsl_e1sti (1, *ipath);

	/*
	 * (5, 2, "IPATH must be equal to 1 or 2 while a value of %(i1) is
	 * given.");
	 */
	imsl_ermes (IMSL_TERMINAL, IMSL_IPATH_RANGE_3);
    }
    if (imsl_n1rcd (0) != 0)
	goto L_9000;

    /*
     * MAKE A COPY OF X IN XPERMU AND WORK WITH XPERMU
     */
    scopy (*n, x, 1, xpermu, 1);
    if (*n == 1)
	goto L_9000;
    for (i = 1; i <= *n; i++) {
	if (ipermu[i - 1] < 1 || ipermu[i - 1] > *n) {
	    imsl_e1sti (1, i);
	    imsl_e1sti (2, *n);
	    imsl_e1sti (3, ipermu[i - 1]);

	    /*
	     * (5, 3, "IPERMU(%(i1)) = %(i3) is not allowed.  It must be
	     * between 1 and N = %(i2).");
	     */
	    imsl_ermes (IMSL_TERMINAL, IMSL_IPERMU_RANGE);
	}
	else {
	    ipermu[i - 1] = -ipermu[i - 1];
	}
    }
    if (imsl_n1rcd (0) != 0)
	goto L_9000;

    if (*ipath == 1) {

	/*
	 * FORWARD PERMUTATION
	 */
	for (i = 1; i <= *n; i++) {
	    if (ipermu[i - 1] <= 0) {
		j = i;
		ipermu[j - 1] = -ipermu[j - 1];
		k = ipermu[j - 1];
	L_20:
		;
		if (ipermu[k - 1] <= 0) {
		    temp = xpermu[j - 1];
		    xpermu[j - 1] = xpermu[k - 1];
		    xpermu[k - 1] = temp;
		    ipermu[k - 1] = -ipermu[k - 1];
		    j = k;
		    k = ipermu[k - 1];
		    goto L_20;
		}
	    }
	}
    }
    else {

	/*
	 * BACKWARD PERMUTATION
	 */
	for (i = 1; i <= *n; i++) {
	    if (ipermu[i - 1] <= 0) {
		ipermu[i - 1] = -ipermu[i - 1];
		j = ipermu[i - 1];
	L_40:
		;
		if (j != i) {
		    temp = xpermu[i - 1];
		    xpermu[i - 1] = xpermu[j - 1];
		    xpermu[j - 1] = temp;
		    ipermu[j - 1] = -ipermu[j - 1];
		    j = ipermu[j - 1];
		    goto L_40;
		}
	    }
	}
    }

L_9000:
    imsl_e1pop ("l_permu");
    return;
}				/* end of function */
