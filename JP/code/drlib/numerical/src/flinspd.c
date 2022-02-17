#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
static VA_LIST_HACK  l_lin_sol_posdef (Mint n, Mfloat a[], Mfloat b[],
                    va_list argptr);
static void     l_l2cds (Mint *n, Mfloat a[], Mint *lda, Mfloat imsl_fac[],
                    Mint *ldfac, Mfloat *rcond, Mfloat z[]);
static void     l_l2nds (Mint *n, Mfloat a[], Mint *lda, Mfloat ainv[], Mint *ldainv, Mfloat wk[]);
static void     l_lftds (Mint *n, Mfloat a[], Mint *lda, Mfloat imsl_fac[], Mint *ldfac);
static void     l_linrt (Mint *n, Mfloat a[], Mint *lda, Mint *ipath,
                    Mfloat ainv[], Mint *ldainv);
static void     l_lslrt (Mint *n, Mfloat a[], Mint *lda, Mfloat b[], Mint *ipath, Mfloat x[]);
static void     l_ssyr (Mchar *uplo, unsigned uplo_s, Mint *n, Mfloat *alpha,
                    Mfloat x[], Mint *incx, Mfloat a[], Mint *lda);
#else
static VA_LIST_HACK  l_lin_sol_posdef ();
static void     l_l2cds ();
static void     l_l2nds ();
static void     l_lftds ();
static void     l_linrt ();
static void     l_lslrt ();
static void     l_ssyr ();
#endif

static Mfloat  *lv_x;
static Mfloat   lv_cond;
#ifdef ANSI
Mfloat         *imsl_f_lin_sol_posdef (Mint n, Mfloat *a, Mfloat *b,...)
#else
Mfloat         *imsl_f_lin_sol_posdef (n, a, b, va_alist)
Mint            n;
Mfloat         *a;
Mfloat         *b;
va_dcl
#endif
{
    va_list         argptr;
    VA_START (argptr, b);

    E1PSH ("imsl_f_lin_sol_posdef", "imsl_d_lin_sol_posdef");
    lv_x = NULL;
    IMSL_CALL (l_lin_sol_posdef (n, a, b, argptr));
    va_end (argptr);
#ifdef DOUBLE
    imsl_e1pop ("imsl_d_lin_sol_posdef");
#else
    imsl_e1pop ("imsl_f_lin_sol_posdef");
#endif
    return lv_x;
}


#ifdef ANSI
static VA_LIST_HACK  l_lin_sol_posdef (Mint n, Mfloat a[], Mfloat b[], va_list argptr)
#else
static VA_LIST_HACK  l_lin_sol_posdef (n, a, b, argptr)
Mint            n;
Mfloat         *a;
Mfloat         *b;
va_list         argptr;
#endif
{
    Mint            code = 1;
    Mint            arg_number = 3;
    Mint            a_col_dim = n;
    Mfloat        **factor_ptr = NULL;
    Mfloat         *factor = NULL;
    Mfloat        **inva_ptr = NULL;
    Mfloat         *inva = NULL;
    Mfloat         *condition = NULL;
    Mint            fac_col_dim = n;
    Mint            inv_col_dim = n;
    short int       return_inverse = 0;
    short int       user_inverse = 0;
    short int       inverse_only = 0;
    short int       return_factor = 0;
    short int       user_factor = 0;
    short int       factor_only = 0;
    short int       solve_only = 0;
    short int       user_solution = 0;
    short int       error = 0;
    short int       factor_transpose = 1;
    short int       return_condition = 0;
    Mint            ldfac = n;
    Mfloat         *work = NULL;
    Mint            i;


    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	    case IMSL_RETURN_USER:
		lv_x = va_arg (argptr, Mfloat *);
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
		factor_ptr = va_arg (argptr, Mfloat **);
		arg_number++;
		break;

	    case IMSL_FACTOR_USER:
		return_factor = 1;
		user_factor = 1;
		factor = va_arg (argptr, Mfloat *);
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
		ldfac = fac_col_dim;
		break;
	    case IMSL_INVERSE:
		return_inverse = 1;
		inva_ptr = va_arg (argptr, Mfloat **);
		arg_number++;
		break;
	    case IMSL_INVERSE_USER:
		user_inverse = 1;
		return_inverse = 1;
		inva = va_arg (argptr, Mfloat *);
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
	    case IMSL_FACTOR_ONLY:
		factor_only = 1;
		break;
	    case IMSL_INVERSE_ONLY:
		inverse_only = 1;
		break;
	    case IMSL_SOLVE_ONLY:
		solve_only = 1;
		break;
	    case IMSL_CONDITION:
		condition = va_arg (argptr, Mfloat *);
		if (!condition) {
		    imsl_e1stl (1, "cond");
		    imsl_e1stl (2, "IMSL_CONDITION");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		    ++error;
		}
		arg_number++;
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
    if (!(factor_only || inverse_only) && b == NULL) {
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
    if (solve_only + factor_only + inverse_only > 1) {
        imsl_ermes (IMSL_TERMINAL, IMSL_BAD_SOLVE_FACTOR_INVERSE);
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
    if (!inverse_only) {
	if (!solve_only && !error) {
	    if (!user_factor) {
		factor = (Mfloat *) imsl_malloc (fac_col_dim * n * sizeof (Mfloat));
		if (!factor) {
		    imsl_ermes (IMSL_TERMINAL, IMSL_NO_MEM_FOR_FAC);
		    ++error;
		}
	    }
	    if (!error) {

		Mint            lda = a_col_dim;
		if (factor_only && !return_condition) {
		    l_lftds (&n, a, &lda, factor, &ldfac);
		    error = imsl_n1rty (1) > 3 ? 1 : 0;
		}
		else {
		    work = (Mfloat *) imsl_malloc (n * sizeof (Mfloat));
		    if (!work) {
			++error;
		    }
		    else {
			Mfloat          cond;
			l_l2cds (&n, a, &lda, factor, &ldfac, &cond, work);
			error = imsl_n1rty (1) > 3 ? 1 : 0;
			if (!error && return_condition) {
			    if (cond > imsl_amach (1)) {
				*condition = 1.0 / cond;
			    }
			    else {
				*condition = imsl_amach (7);
			    }
			}
			if (work) imsl_free (work);
			work = NULL;
		    }
		}
	    }
	}
	if (!factor_only && !error) {
	    if (!user_solution) {
		lv_x = (Mfloat *) imsl_malloc (n * sizeof (Mfloat));
		if (!lv_x) {
		    imsl_ermes (IMSL_TERMINAL, IMSL_NO_MEM_FOR_SYS);
		    ++error;
		}
	    }
	    if (!error) {
		Mint            l = 4;
	    /* solve ctrans(r)*y = b */
		l_lslrt (&n, factor, &ldfac, b, &l, lv_x);
		error = imsl_n1rty (1) > 3 ? 1 : 0;

		if (!error) {
		/* solve r*x = y */
		    l = 2;
		    l_lslrt (&n, factor, &ldfac, lv_x, &l, lv_x);
		    error = imsl_n1rty (1) > 3 ? 1 : 0;
		}
	    }
	}
    }
    if (return_inverse) {
	Mint            ldinva = n;
	work = (Mfloat *) imsl_malloc (n * sizeof (Mfloat));
	if (!work) {
	    ++error;
	}
	else {
	    if (!user_inverse) {
		inva = (Mfloat *) imsl_malloc (n * inv_col_dim * sizeof (Mfloat));
		if (!inva)
		    ++error;
	    }
	    else if (inv_col_dim > n) {
		imsl_f_m1ran (n, inv_col_dim, inva, inva);
		error = (imsl_n1rty (1) > 3) ? 1 : 0;
	    }
	    if (!error) {
		l_l2nds (&n, a, &a_col_dim, inva, &ldinva, work);
		if (!(error = (imsl_n1rty (1) > 3) ? 1 : 0)) {
		    imsl_f_m1ran (inv_col_dim, n, inva, inva);
		    if (!user_inverse) {
			for (i = n; i < inv_col_dim; i++) {
			    sset (n, 0.0, (inva + i), inv_col_dim);
			}
			*inva_ptr = inva;
			inva = NULL;
		    }

		    if (return_condition) {
			if (lv_cond > imsl_amach (1)) {
			    *condition = 1.0 / lv_cond;
			}
			else {
			    *condition = imsl_amach (7);
			}
		    }
		}
		else if (!user_inverse) {
		    imsl_free (inva);
		    inva = NULL;
		}
	    }
	    imsl_free (work);
	    work = NULL;
	}
    }


    if (!user_factor) {
	if (!return_factor || error) {
	    if (factor != NULL) imsl_free (factor);
	}
	else {
	    for (i = n; i < fac_col_dim; i++) {
		sset (n, 0.0, (factor + i), fac_col_dim);
	    }
	    *factor_ptr = factor;
	}
	factor = NULL;
    }

    if (error && !user_solution) {
	if (lv_x != NULL) imsl_free (lv_x);
	lv_x = NULL;
    }
    return argptr;
}



/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 19:11:16
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  L2CDS/DL2CDS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the trans(R)*R Cholesky factorization of a real
                symmetric positive definite matrix and estimate its L1
                condition number.

    Usage:      CALL L2CDS (N, A, LDA, FAC, LDFAC, RCOND, Z)

    Arguments:  See LFCDS/DLFCDS.

    Remarks:    See LFCDS/DLFCDS.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void     l_l2cds (Mint *n, Mfloat a[], Mint *lda, Mfloat imsl_fac[],
                    Mint *ldfac, Mfloat *rcond, Mfloat z[])
#else
static void     l_l2cds (n, a, lda, imsl_fac, ldfac, rcond, z)
Mint           *n;
Mfloat         *a;
Mint           *lda;
Mfloat         *imsl_fac;
Mint           *ldfac;
Mfloat         *rcond, z[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define FAC(I_,J_)	(imsl_fac+(I_)*(*ldfac)+(J_))
    Mint            _l0, _l1, j, k;
    Mfloat          anorm, big, ek, s, sm, small, suma, t, wk, wkm, ynorm;


    imsl_e1psh ("l_l2cds");

    if (*n <= 0) {
	imsl_e1sti (1, *n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
    }
    if (*n > *lda) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *lda);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
    }
    if (*n > *ldfac) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *ldfac);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDFAC_LESS_ORDER);
    }
    if (imsl_n1rcd (0) != 0)
	goto L_9000;
 /*
  * COMPUTE 1-NORM OF A USING ONLY THE UPPER TRIANGLE OF A
  */
    *rcond = F_ZERO;
    anorm = F_ZERO;
    _l0 = 1;
    for (j = 1; j <= *n; j++) {
	if (j < *n) {
	    _l1 = *n - j;
	    suma = imsl_sasum (j, A (j - 1, 0), _l0) + imsl_sasum (_l1,
		A (j, j - 1), *lda);
	}
	else {
	    suma = imsl_sasum (j, A (j - 1, 0), _l0);
	}
	anorm = imsl_f_max (anorm, suma);
    }
 /* FACTOR */
    l_lftds (n, a, lda, imsl_fac, ldfac);
    if (imsl_n1rcd (1) != 0)
	goto L_9000;
 /*
  * RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . ESTIMATE =
  * NORM(Z)/NORM(Y) WHERE A*Z = Y AND A*Y = E . THE COMPONENTS OF E ARE
  * CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W WHERE
  * TRANS(R)*W = E . THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
  * SOLVE TRANS(R)*W = E
  */
    ek = F_ONE;
    sset (*n, F_ZERO, z, 1);
    small = imsl_amach (1);
    big = imsl_amach (2);
    if (small * big <= F_ONE)
	small = F_ONE / big;

    for (k = 1; k <= *n; k++) {
	if (z[k - 1] != F_ZERO)
	    ek = sign (ek, -z[k - 1]);
	if (fabs (ek - z[k - 1]) > *FAC (k - 1, k - 1)) {
	    s = *FAC (k - 1, k - 1) / fabs (ek - z[k - 1]);
	    sscal (*n, s, z, 1);
	    ek *= s;
	}
	wk = ek - z[k - 1];
	wkm = -ek - z[k - 1];
	s = fabs (wk);
	sm = fabs (wkm);

	if (fabs (*FAC (k - 1, k - 1)) > small) {
	    wk /= *FAC (k - 1, k - 1);
	    wkm /= *FAC (k - 1, k - 1);
	}
	if (k + 1 <= *n) {
	    for (j = k + 1; j <= *n; j++) {
		sm += fabs (z[j - 1] + wkm ** FAC (j - 1, k - 1));
	    }
	    saxpy (*n - k, wkm, FAC (k, k - 1), *ldfac, &z[k], 1);
	    _l0 = *n - k;
	    _l1 = 1;
	    s = imsl_sasum (_l0, &z[k], _l1);

	    if (s < sm) {
		t = wkm - wk;
		wk = wkm;
		saxpy (*n - k, t, FAC (k, k - 1), *ldfac, &z[k], 1);
	    }
	}
	z[k - 1] = wk;
    }
    _l0 = 1;
    s = F_ONE / imsl_sasum (*n, z, _l0);
    sscal (*n, s, z, 1);
 /* SOLVE R*Y = W */
    for (k = *n; k >= 1; k--) {
	if (fabs (z[k - 1]) > *FAC (k - 1, k - 1)) {
	    s = *FAC (k - 1, k - 1) / fabs (z[k - 1]);
	    sscal (*n, s, z, 1);
	}
	if (fabs (*FAC (k - 1, k - 1)) > small)
	    z[k - 1] /= *FAC (k - 1, k - 1);
	t = -z[k - 1];
	saxpy (k - 1, t, FAC (k - 1, 0), 1, &z[0], 1);
    }
    _l0 = 1;
    s = F_ONE / imsl_sasum (*n, z, _l0);
    sscal (*n, s, z, 1);

    ynorm = F_ONE;
 /* SOLVE TRANS(R)*V = Y */
    for (k = 1; k <= *n; k++) {
	z[k - 1] -= imsl_sdot (k - 1, FAC (k - 1, 0), 1, &z[0], 1);
	if (fabs (z[k - 1]) > *FAC (k - 1, k - 1)) {
	    s = *FAC (k - 1, k - 1) / fabs (z[k - 1]);
	    sscal (*n, s, z, 1);
	    ynorm *= s;
	}
	if (fabs (*FAC (k - 1, k - 1)) > small)
	    z[k - 1] /= *FAC (k - 1, k - 1);
    }
    _l0 = 1;
    s = F_ONE / imsl_sasum (*n, z, _l0);
    sscal (*n, s, z, 1);
    ynorm *= s;
 /* SOLVE R*Z = V */
    for (k = *n; k >= 1; k--) {
	if (fabs (z[k - 1]) > *FAC (k - 1, k - 1)) {
	    s = *FAC (k - 1, k - 1) / fabs (z[k - 1]);
	    sscal (*n, s, z, 1);
	    ynorm *= s;
	}
	if (fabs (*FAC (k - 1, k - 1)) > small)
	    z[k - 1] /= *FAC (k - 1, k - 1);
	t = -z[k - 1];
	saxpy (k - 1, t, FAC (k - 1, 0), 1, &z[0], 1);
    }
 /* MAKE ZNORM = 1.0 */
    _l0 = 1;
    s = F_ONE / imsl_sasum (*n, z, _l0);
    sscal (*n, s, z, 1);
    ynorm *= s;

    if (anorm != F_ZERO)
	*rcond = ynorm / anorm;
    if (*rcond <= imsl_amach (4)) {
	imsl_e1str (1, *rcond);

    /*
     * (3, 1, "The matrix is algorithmically singular.  An estimate of the
     * reciprocal of its L1 condition number is RCOND = %(r1).");
     */
	imsl_ermes (IMSL_WARNING, IMSL_ILL_CONDITIONED);
    }
L_9000:
    imsl_e1pop ("l_l2cds");
    return;
}	       /* end of function */
#undef FAC
#undef A
/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 19:15:55
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  L2NDS/DL2NDS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the inverse of a real symmetric positive definite
                matrix.

    Usage:      CALL L2NDS (N, A, LDA, AINV, LDAINV, WK)

    Arguments:  See LINDS/DLINDS.

    Remarks:    See LINDS/DLINDS.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void     l_l2nds (Mint *n, Mfloat a[], Mint *lda, Mfloat ainv[], Mint *ldainv, Mfloat wk[])
#else
static void     l_l2nds (n, a, lda, ainv, ldainv, wk)
Mint           *n;
Mfloat         *a;
Mint           *lda;
Mfloat         *ainv;
Mint           *ldainv;
Mfloat          wk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define AINV(I_,J_)	(ainv+(I_)*(*ldainv)+(J_))
    Mint            _l0, _l1, j;
    Mfloat          _f0, temp;


    imsl_e1psh ("l_l2nds");

    if (*n <= 0) {
	imsl_e1sti (1, *n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
    }
    if (*n > *lda) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *lda);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
    }
    if (*n > *ldainv) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *ldainv);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDAINV_LESS_ORDER);
    }
    if (imsl_n1rcd (0) != 0)
	goto L_9000;
 /*
  * COMPUTE THE FACTORIZATION OF A AND ESTIMATE ITS CONDITION NUMBER
  */
    l_l2cds (n, a, lda, ainv, ldainv, &lv_cond, wk);
    if (imsl_n1rty (1) == 4)
	goto L_9000;
 /* COMPUTE INVERSE(R) */
    _l0 = 2;
    l_linrt (n, ainv, ldainv, &_l0, ainv, ldainv);
 /* FORM INVERSE(R) * TRANS(INVERSE(R)) */
    _f0 = F_ONE;
    _l1 = 1;
    for (j = 1; j <= *n; j++) {
	_l0 = j - 1;
	l_ssyr ("U", sizeof ("U"), &_l0, &_f0, AINV (j - 1, 0),
	    &_l1, ainv, ldainv);
	temp = *AINV (j - 1, j - 1);
	sscal (j, temp, AINV (j - 1, 0), 1);
    }

    if (lv_cond <= imsl_amach (4)) {
	imsl_e1str (1, lv_cond);
	imsl_ermes (IMSL_WARNING, IMSL_ILL_CONDITIONED);
    }
 /* FILL LOWER TRIANGLE OF AINV */
    imsl_csfrg (n, ainv, ldainv);

L_9000:
    imsl_e1pop ("l_l2nds");
    return;
}	       /* end of function */
#undef AINV
#undef A

/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 19:14:42
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
static void     l_lftds (Mint *n, Mfloat a[], Mint *lda, Mfloat imsl_fac[], Mint *ldfac)
#else
static void     l_lftds (n, a, lda, imsl_fac, ldfac)
Mint           *n;
Mfloat         *a;
Mint           *lda;
Mfloat         *imsl_fac;
Mint           *ldfac;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define FAC(I_,J_)	(imsl_fac+(I_)*(*ldfac)+(J_))
    Mint            i, info, j, k, l, ll, nlenb, nrcfac;
    Mfloat          big, r0, r1, r2, r3, r4, r5, r6, r7, rtemp, small,
                    t[8][NB];


 /*
  * this code is for computer types: fosivv, rtxlxs, and vxvmsv
  */
    imsl_e1psh ("l_lftds");

    if (*n <= 0) {
	imsl_e1sti (1, *n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
    }
    if (*n > *lda) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *lda);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
    }
    if (*n > *ldfac) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *ldfac);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDFAC_LESS_ORDER);
    }
    if (imsl_n1rcd (0) != 0)
	goto L_9000;
 /* Preserve a copy of the input matrix */
    for (i = 1; i <= *n; i++) {
	scopy (i, A (i - 1, 0), 1, FAC (i - 1, 0), 1);
    }

 /*
  * Cholesky decomposition using method LLT**
  * 
  * A brief description of the algorithm follows: For a symmetric positive
  * definite matrix at the k-th step :
  * 
  * k   |  11  12  13 | A  = |  21  22  23 | , |  31  32  33 |
  * 
  * where trans(A11 A21 A31) is n x (k*8) and factored. assume trans(A22 A32) is
  * the active block of 8 columns. The factorization is accomplished by:
  * 
  * Step 1. Factor trans(A22 A32) and store -trans(A22 A32) into a local work
  * array
  * 
  * Step 2. Update A32 with matrix multiplication between trans(A22 A32) and the
  * local work array.
  * 
  * Step 3. repeat
  * 
  * Fill in the lower triangle
  */
    imsl_csfrg (n, imsl_fac, ldfac);

    small = imsl_amach (1);
    big = imsl_amach (2);
    if (small * big < F_ONE)
	small = F_ONE / big;
    info = 0;

    nrcfac = mod (*n, 8);
    for (j = 1; j <= (*n - nrcfac); j += 8) {
	nlenb = imsl_i_min (NB, *n - j);
    /* prepare j-th column */
	if (*FAC (j - 1, j - 1) <= F_ZERO) {
	    info = j;
	    goto L_450;
	}
	*FAC (j - 1, j - 1) = sqrt (*FAC (j - 1, j - 1));
	r0 = F_ONE / *FAC (j - 1, j - 1);
    /*
     * Form the j-th mutiplier and load t(1:nlenb,1) with
     * -imsl_fac(j+1:nlenb,j)
     */
	for (i = 1; i <= nlenb; i++) {
	    *FAC (j - 1, j + i - 1) *= r0;
	    t[0][i - 1] = -*FAC (j - 1, j + i - 1);
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j - 1, i - 1) *= r0;
	}
    /* update columns j+1 thru j+7 */
	for (k = 1; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) += *FAC (j - 1, j + i - 1) *
		    t[0][k - 1];
	    }
	}
    /* prepare (j+1)-th column */
	if (*FAC (j, j) <= F_ZERO) {
	    info = j + 1;
	    goto L_450;
	}
	*FAC (j, j) = sqrt (*FAC (j, j));
	r1 = F_ONE / *FAC (j, j);
    /*
     * Form the (j+1)-th mutiplier and load t(2:nlenb,2) with
     * -imsl_fac(j+2:nlenb,j+1)
     */
	for (i = 2; i <= nlenb; i++) {
	    *FAC (j, j + i - 1) *= r1;
	    t[1][i - 1] = -*FAC (j, j + i - 1);
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j, i - 1) *= r1;
	}
    /* update columns j+2 thru j+7 */
	for (k = 2; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) += *FAC (j, j + i - 1) * t[1][k - 1];
	    }
	}
    /* prepare (j+2)-th column */
	if (*FAC (j + 1, j + 1) <= F_ZERO) {
	    info = j + 2;
	    goto L_450;
	}
	*FAC (j + 1, j + 1) = sqrt (*FAC (j + 1, j + 1));
	r2 = F_ONE / *FAC (j + 1, j + 1);
    /*
     * Form the (j+2)-th mutiplier and load t(3:nlenb,3) with
     * -imsl_fac(j+3:nlenb,j+2)
     */
	for (i = 3; i <= nlenb; i++) {
	    *FAC (j + 1, j + i - 1) *= r2;
	    t[2][i - 1] = -*FAC (j + 1, j + i - 1);
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 1, i - 1) *= r2;
	}
    /* update columns j+3 thru j+7 */
	for (k = 3; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) += *FAC (j + 1, j + i - 1) *
		    t[2][k - 1];
	    }
	}
    /* prepare (j+3)-th column */
	if (*FAC (j + 2, j + 2) <= F_ZERO) {
	    info = j + 3;
	    goto L_450;
	}
	*FAC (j + 2, j + 2) = sqrt (*FAC (j + 2, j + 2));
	r3 = F_ONE / *FAC (j + 2, j + 2);
    /*
     * Form the (j+3)-th mutiplier and load t(4:nlenb,4) with
     * -imsl_fac(j+3:nlenb,j+3)
     */
	for (i = 4; i <= nlenb; i++) {
	    *FAC (j + 2, j + i - 1) *= r3;
	    t[3][i - 1] = -*FAC (j + 2, j + i - 1);
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 2, i - 1) *= r3;
	}
    /* update columns j+4 thru j+7 */
	for (k = 4; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) += *FAC (j + 2, j + i - 1) *
		    t[3][k - 1];
	    }
	}
    /* prepare (j+4)-th column */
	if (*FAC (j + 3, j + 3) <= F_ZERO) {
	    info = j + 4;
	    goto L_450;
	}
	*FAC (j + 3, j + 3) = sqrt (*FAC (j + 3, j + 3));
	r4 = F_ONE / *FAC (j + 3, j + 3);
    /*
     * Form the (j+4)-th mutiplier and load t(5:nlenb,5) with
     * -imsl_fac(j+4:nlenb,j+4)
     */
	for (i = 5; i <= nlenb; i++) {
	    *FAC (j + 3, j + i - 1) *= r4;
	    t[4][i - 1] = -*FAC (j + 3, j + i - 1);
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 3, i - 1) *= r4;
	}
    /* update columns j+5 thru j+7 */
	for (k = 5; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) += *FAC (j + 3, j + i - 1) *
		    t[4][k - 1];
	    }
	}
    /* prepare (j+5)-th column */
	if (*FAC (j + 4, j + 4) <= F_ZERO) {
	    info = j + 5;
	    goto L_450;
	}
	*FAC (j + 4, j + 4) = sqrt (*FAC (j + 4, j + 4));
	r5 = F_ONE / *FAC (j + 4, j + 4);
    /*
     * Form the (j+5)-th mutiplier and load t(5:nlenb,5) with
     * -imsl_fac(j+3:nlenb,j+3)
     */
	for (i = 6; i <= nlenb; i++) {
	    *FAC (j + 4, j + i - 1) *= r5;
	    t[5][i - 1] = -*FAC (j + 4, j + i - 1);
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 4, i - 1) *= r5;
	}
    /* update columns j+6 thru j+7 */
	for (k = 6; k <= 7; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) += *FAC (j + 4, j + i - 1) *
		    t[5][k - 1];
	    }
	}
    /* prepare (j+6)-th column */
	if (*FAC (j + 5, j + 5) <= F_ZERO) {
	    info = j + 6;
	    goto L_450;
	}
	*FAC (j + 5, j + 5) = sqrt (*FAC (j + 5, j + 5));
	r6 = F_ONE / *FAC (j + 5, j + 5);
    /*
     * Form the (j+6)-th mutiplier and load t(7:nlenb,7) with
     * -imsl_fac(j+7:nlenb,j+6)
     */
	for (i = 7; i <= nlenb; i++) {
	    *FAC (j + 5, j + i - 1) *= r6;
	    t[6][i - 1] = -*FAC (j + 5, j + i - 1);
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 5, i - 1) *= r6;
	}
    /* update column j+7 */
	for (i = 7; i <= (*n - j); i++) {
	    *FAC (j + 6, j + i - 1) += *FAC (j + 5, j + i - 1) * t[6][6];
	}
    /* prepare (j+7)-th column */
	if (*FAC (j + 6, j + 6) <= F_ZERO) {
	    info = j + 7;
	    goto L_450;
	}
	*FAC (j + 6, j + 6) = sqrt (*FAC (j + 6, j + 6));
	r7 = F_ONE / *FAC (j + 6, j + 6);
    /*
     * Form the (j+7)-th mutiplier and load t(8:nlenb,8) with
     * -imsl_fac(j+8:nlenb,j+7)
     */
	for (i = 8; i <= nlenb; i++) {
	    *FAC (j + 6, j + i - 1) *= r7;
	    t[7][i - 1] = -*FAC (j + 6, j + i - 1);
	}
	for (i = j + nlenb + 1; i <= *n; i++) {
	    *FAC (j + 6, i - 1) *= r7;
	}
    /*
     * Perform update on the lower triangle on columns j+7 thru j + nlenb
     * rows j+7 thru n
     */
	for (k = 8; k <= nlenb; k++) {
	    for (i = k; i <= (*n - j); i++) {
		*FAC (j + k - 1, j + i - 1) += *FAC (j - 1, j + i - 1) *
		    t[0][k - 1] + *FAC (j, j + i - 1) * t[1][k - 1] + *FAC (j + 1, j + i - 1) *
		    t[2][k - 1] + *FAC (j + 2, j + i - 1) * t[3][k - 1] +
		    *FAC (j + 3, j + i - 1) * t[4][k - 1] + *FAC (j + 4, j + i - 1) *
		    t[5][k - 1] + *FAC (j + 5, j + i - 1) * t[6][k - 1] +
		    *FAC (j + 6, j + i - 1) * t[7][k - 1];
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
		    t[k][i - 1] = -*FAC (j + k - 1, l + i - 1);
		}
	    }
	/*
	 * Update lower triangle from rows ll thru ll+n columns ll thru ll +
	 * nlenb
	 */
	    for (k = 1; k <= nlenb; k++) {
		for (i = k; i <= (*n - l); i++) {
		    *FAC (l + k - 1, l + i - 1) += *FAC (j - 1, l + i - 1) *
			t[0][k - 1] + *FAC (j, l + i - 1) * t[1][k - 1] +
			*FAC (j + 1, l + i - 1) * t[2][k - 1] + *FAC (j + 2, l + i - 1) *
			t[3][k - 1] + *FAC (j + 3, l + i - 1) * t[4][k - 1] +
			*FAC (j + 4, l + i - 1) * t[5][k - 1] + *FAC (j + 5, l + i - 1) *
			t[6][k - 1] + *FAC (j + 6, l + i - 1) * t[7][k - 1];
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
    if (*FAC (j - 1, j - 1) <= F_ZERO) {
	info = j;
	goto L_450;
    }
    *FAC (j - 1, j - 1) = sqrt (*FAC (j - 1, j - 1));
    rtemp = F_ONE / *FAC (j - 1, j - 1);
    *FAC (j - 1, j) *= rtemp;
    *FAC (j - 1, j + 1) *= rtemp;
    *FAC (j - 1, j + 2) *= rtemp;
    *FAC (j - 1, j + 3) *= rtemp;
    *FAC (j - 1, j + 4) *= rtemp;
    *FAC (j - 1, j + 5) *= rtemp;
    *FAC (j, j) += -imsl_fi_power (*FAC (j - 1, j), 2);
    *FAC (j, j + 1) += -*FAC (j - 1, j + 1) ** FAC (j - 1, j);
    *FAC (j, j + 2) += -*FAC (j - 1, j + 2) ** FAC (j - 1, j);
    *FAC (j, j + 3) += -*FAC (j - 1, j + 3) ** FAC (j - 1, j);
    *FAC (j, j + 4) += -*FAC (j - 1, j + 4) ** FAC (j - 1, j);
    *FAC (j, j + 5) += -*FAC (j - 1, j + 5) ** FAC (j - 1, j);
    *FAC (j + 1, j + 1) += -imsl_fi_power (*FAC (j - 1, j + 1), 2);
    *FAC (j + 1, j + 2) += -*FAC (j - 1, j + 2) ** FAC (j - 1, j + 1);
    *FAC (j + 1, j + 3) += -*FAC (j - 1, j + 3) ** FAC (j - 1, j + 1);
    *FAC (j + 1, j + 4) += -*FAC (j - 1, j + 4) ** FAC (j - 1, j + 1);
    *FAC (j + 1, j + 5) += -*FAC (j - 1, j + 5) ** FAC (j - 1, j + 1);
    *FAC (j + 2, j + 2) += -imsl_fi_power (*FAC (j - 1, j + 2), 2);
    *FAC (j + 2, j + 3) += -*FAC (j - 1, j + 3) ** FAC (j - 1, j + 2);
    *FAC (j + 2, j + 4) += -*FAC (j - 1, j + 4) ** FAC (j - 1, j + 2);
    *FAC (j + 2, j + 5) += -*FAC (j - 1, j + 5) ** FAC (j - 1, j + 2);
    *FAC (j + 3, j + 3) += -imsl_fi_power (*FAC (j - 1, j + 3), 2);
    *FAC (j + 3, j + 4) += -*FAC (j - 1, j + 4) ** FAC (j - 1, j + 3);
    *FAC (j + 3, j + 5) += -*FAC (j - 1, j + 5) ** FAC (j - 1, j + 3);
    *FAC (j + 4, j + 4) += -imsl_fi_power (*FAC (j - 1, j + 4), 2);
    *FAC (j + 4, j + 5) += -*FAC (j - 1, j + 5) ** FAC (j - 1, j + 4);
    *FAC (j + 5, j + 5) += -imsl_fi_power (*FAC (j - 1, j + 5), 2);
    j += 1;
    nrcfac -= 1;
L_390:
    if (nrcfac < 6)
	goto L_400;
    if (*FAC (j - 1, j - 1) <= F_ZERO) {
	info = j;
	goto L_450;
    }
    *FAC (j - 1, j - 1) = sqrt (*FAC (j - 1, j - 1));
    rtemp = F_ONE / *FAC (j - 1, j - 1);
    *FAC (j - 1, j) *= rtemp;
    *FAC (j - 1, j + 1) *= rtemp;
    *FAC (j - 1, j + 2) *= rtemp;
    *FAC (j - 1, j + 3) *= rtemp;
    *FAC (j - 1, j + 4) *= rtemp;
    *FAC (j, j) += -imsl_fi_power (*FAC (j - 1, j), 2);
    *FAC (j, j + 1) += -*FAC (j - 1, j + 1) ** FAC (j - 1, j);
    *FAC (j, j + 2) += -*FAC (j - 1, j + 2) ** FAC (j - 1, j);
    *FAC (j, j + 3) += -*FAC (j - 1, j + 3) ** FAC (j - 1, j);
    *FAC (j, j + 4) += -*FAC (j - 1, j + 4) ** FAC (j - 1, j);
    *FAC (j + 1, j + 1) += -imsl_fi_power (*FAC (j - 1, j + 1), 2);
    *FAC (j + 1, j + 2) += -*FAC (j - 1, j + 2) ** FAC (j - 1, j + 1);
    *FAC (j + 1, j + 3) += -*FAC (j - 1, j + 3) ** FAC (j - 1, j + 1);
    *FAC (j + 1, j + 4) += -*FAC (j - 1, j + 4) ** FAC (j - 1, j + 1);
    *FAC (j + 2, j + 2) += -imsl_fi_power (*FAC (j - 1, j + 2), 2);
    *FAC (j + 2, j + 3) += -*FAC (j - 1, j + 3) ** FAC (j - 1, j + 2);
    *FAC (j + 2, j + 4) += -*FAC (j - 1, j + 4) ** FAC (j - 1, j + 2);
    *FAC (j + 3, j + 3) += -imsl_fi_power (*FAC (j - 1, j + 3), 2);
    *FAC (j + 3, j + 4) += -*FAC (j - 1, j + 4) ** FAC (j - 1, j + 3);
    *FAC (j + 4, j + 4) += -imsl_fi_power (*FAC (j - 1, j + 4), 2);
    j += 1;
    nrcfac -= 1;
L_400:
    if (nrcfac < 5)
	goto L_410;
    if (*FAC (j - 1, j - 1) <= F_ZERO) {
	info = j;
	goto L_450;
    }
    *FAC (j - 1, j - 1) = sqrt (*FAC (j - 1, j - 1));
    rtemp = F_ONE / *FAC (j - 1, j - 1);
    *FAC (j - 1, j) *= rtemp;
    *FAC (j - 1, j + 1) *= rtemp;
    *FAC (j - 1, j + 2) *= rtemp;
    *FAC (j - 1, j + 3) *= rtemp;
    *FAC (j, j) += -imsl_fi_power (*FAC (j - 1, j), 2);
    *FAC (j, j + 1) += -*FAC (j - 1, j + 1) ** FAC (j - 1, j);
    *FAC (j, j + 2) += -*FAC (j - 1, j + 2) ** FAC (j - 1, j);
    *FAC (j, j + 3) += -*FAC (j - 1, j + 3) ** FAC (j - 1, j);
    *FAC (j + 1, j + 1) += -imsl_fi_power (*FAC (j - 1, j + 1), 2);
    *FAC (j + 1, j + 2) += -*FAC (j - 1, j + 2) ** FAC (j - 1, j + 1);
    *FAC (j + 1, j + 3) += -*FAC (j - 1, j + 3) ** FAC (j - 1, j + 1);
    *FAC (j + 2, j + 2) += -imsl_fi_power (*FAC (j - 1, j + 2), 2);
    *FAC (j + 2, j + 3) += -*FAC (j - 1, j + 3) ** FAC (j - 1, j + 2);
    *FAC (j + 3, j + 3) += -imsl_fi_power (*FAC (j - 1, j + 3), 2);
    j += 1;
    nrcfac -= 1;
L_410:
    if (nrcfac < 4)
	goto L_420;
    if (*FAC (j - 1, j - 1) <= F_ZERO) {
	info = j;
	goto L_450;
    }
    *FAC (j - 1, j - 1) = sqrt (*FAC (j - 1, j - 1));
    rtemp = F_ONE / *FAC (j - 1, j - 1);
    *FAC (j - 1, j) *= rtemp;
    *FAC (j - 1, j + 1) *= rtemp;
    *FAC (j - 1, j + 2) *= rtemp;
    *FAC (j, j) += -imsl_fi_power (*FAC (j - 1, j), 2);
    *FAC (j, j + 1) += -*FAC (j - 1, j + 1) ** FAC (j - 1, j);
    *FAC (j, j + 2) += -*FAC (j - 1, j + 2) ** FAC (j - 1, j);
    *FAC (j + 1, j + 1) += -imsl_fi_power (*FAC (j - 1, j + 1), 2);
    *FAC (j + 1, j + 2) += -*FAC (j - 1, j + 2) ** FAC (j - 1, j + 1);
    *FAC (j + 2, j + 2) += -imsl_fi_power (*FAC (j - 1, j + 2), 2);
    j += 1;
    nrcfac -= 1;
L_420:
    if (nrcfac < 3)
	goto L_430;
    if (*FAC (j - 1, j - 1) <= F_ZERO) {
	info = j;
	goto L_450;
    }
    *FAC (j - 1, j - 1) = sqrt (*FAC (j - 1, j - 1));
    rtemp = F_ONE / *FAC (j - 1, j - 1);
    *FAC (j - 1, j) *= rtemp;
    *FAC (j - 1, j + 1) *= rtemp;
    *FAC (j, j) += -imsl_fi_power (*FAC (j - 1, j), 2);
    *FAC (j, j + 1) += -*FAC (j - 1, j + 1) ** FAC (j - 1, j);
    *FAC (j + 1, j + 1) += -imsl_fi_power (*FAC (j - 1, j + 1), 2);
    j += 1;
    nrcfac -= 1;
L_430:
    if (nrcfac < 2)
	goto L_440;
    if (*FAC (j - 1, j - 1) <= F_ZERO) {
	info = j;
	goto L_450;
    }
    *FAC (j - 1, j - 1) = sqrt (*FAC (j - 1, j - 1));
    rtemp = F_ONE / *FAC (j - 1, j - 1);
    *FAC (j - 1, j) *= rtemp;
    *FAC (j, j) += -imsl_fi_power (*FAC (j - 1, j), 2);
    j += 1;
    nrcfac -= 1;
L_440:
    if (nrcfac < 1)
	goto L_450;
    if (*FAC (j - 1, j - 1) <= F_ZERO) {
	info = j;
	goto L_450;
    }
    *FAC (j - 1, j - 1) = sqrt (*FAC (j - 1, j - 1));
L_450:
    if (info != 0) {
	imsl_e1sti (1, info);

    /*
     * (4, 2, "The leading %(i1) by %(i1) submatrix of the input matrix is
     * not positive definite.");
     */
	imsl_ermes (IMSL_FATAL, IMSL_NONPOSITIVE_MATRIX);
    }
 /* Fill in the upper triangle */
    for (i = 1; i <= (*n - 1); i++) {
	scopy (*n - i, FAC (i - 1, i), 1, FAC (i, i - 1), *ldfac);
    }

L_9000:
    imsl_e1pop ("l_lftds");

    return;
}	       /* end of function */
#undef FAC
#undef A
/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 19:17:49
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  LINRT/DLINRT  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the inverse of a real triangular matrix.

    Usage:      CALL LINRT (N, A, LDA, IPATH, AINV, LDAINV)

    Arguments:
       N      - Order of the matrix.  (Input)
       A      - N by N matrix containing the triangular matrix to be
                inverted.  (Input)
                For a lower triangular matrix, only the lower triangular
                part and diagonal of A are referenced.  For an upper
                triangular matrix, only the upper triangular part and
                diagonal of A are referenced.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPATH  - Path indicator.  (Input)
                IPATH = 1 means A is lower triangular,
                IPATH = 2 means A is upper triangular.
       AINV   - N by N matrix containing the inverse of A.  (Output)
                If A is lower triangular, AINV is also lower triangular.
                If A is upper triangular, AINV is also upper triangular.
                If A is not needed, A and AINV can share the same storage
                locations.
       LDAINV - Leading dimension of AINV exactly as specified in the
                dimension statement of the calling program.  (Input)

    GAMS:       D2a3

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void     l_linrt (Mint *n, Mfloat a[], Mint *lda, Mint *ipath,
                    Mfloat ainv[], Mint *ldainv)
#else
static void     l_linrt (n, a, lda, ipath, ainv, ldainv)
Mint           *n;
Mfloat         *a;
Mint           *lda, *ipath;
Mfloat         *ainv;
Mint           *ldainv;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define AINV(I_,J_)	(ainv+(I_)*(*ldainv)+(J_))
    Mint            info, j, k;
    Mfloat          big, small, temp;


    imsl_e1psh ("l_linrt");

    if (*n <= 0) {
	imsl_e1sti (1, *n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
    }
    else if (*n > *lda) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *lda);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
    }
    else if (*n > *ldainv) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *ldainv);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDAINV_LESS_ORDER);
    }
    else if (*ipath != 1 && *ipath != 2) {
	imsl_e1sti (1, *ipath);
	imsl_ermes (IMSL_TERMINAL, IMSL_IPATH_RANGE_3);
    }
    if (imsl_n1rcd (0) != 0)
	goto L_9000;

    small = imsl_amach (1);
    big = imsl_amach (2);
    if (small * big < F_ONE)
	small = F_ONE / big;
    if (*ipath == 1) {
    /*
     * MAKE A COPY OF A IN AINV AND ZERO THE STRICTLY UPPER TRIANGLE OF AINV
     */
	for (j = 1; j <= *n; j++) {
	    sset (j - 1, F_ZERO, AINV (j - 1, 0), 1);
	    scopy (*n - j + 1, A (j - 1, j - 1), 1, AINV (j - 1, j - 1),
		1);
	}
    /*
     * COMPUTE INVERSE OF LOWER TRIANGULAR MATRIX
     */
	for (k = *n; k >= 1; k--) {
	    info = k;
	    if (fabs (*AINV (k - 1, k - 1)) <= small)
		goto L_50;
	    *AINV (k - 1, k - 1) = F_ONE / *AINV (k - 1, k - 1);
	    temp = -*AINV (k - 1, k - 1);
	    if (k < *n) {
		sscal (*n - k, temp, AINV (k - 1, k), 1);
		imsl_sger (*n - k, k - 1, F_ONE, AINV (k - 1, k), 1, AINV (0, k - 1),
		    *ldainv, AINV (0, k), *ldainv);
	    }
	    sscal (k - 1, *AINV (k - 1, k - 1), AINV (0, k - 1), *ldainv);
	}
	info = 0;
    }
    else if (*ipath == 2) {
    /*
     * MAKE A COPY OF A IN AINV AND ZERO THE STRICTLY LOWER TRIANGLE OF AINV
     */
	for (j = 1; j <= *n; j++) {
	    if (j < *n)
		sset (*n - j, F_ZERO, AINV (j - 1, j), 1);
	    scopy (j, A (j - 1, 0), 1, AINV (j - 1, 0), 1);
	}
    /*
     * COMPUTE INVERSE OF AN UPPER TRIANGULAR MATRIX
     */
	for (k = 1; k <= *n; k++) {
	    info = k;
	    if (fabs (*AINV (k - 1, k - 1)) <= small)
		goto L_50;
	    *AINV (k - 1, k - 1) = F_ONE / *AINV (k - 1, k - 1);
	    temp = -*AINV (k - 1, k - 1);
	    sscal (k - 1, temp, AINV (k - 1, 0), 1);
	    if (k < *n) {
		imsl_sger (k - 1, *n - k, F_ONE, AINV (k - 1, 0), 1, AINV (k, k - 1),
		    *ldainv, AINV (k, 0), *ldainv);
		sscal (*n - k, *AINV (k - 1, k - 1), AINV (k, k - 1),
		    *ldainv);
	    }
	}
	info = 0;
    }
L_50:
    if (info != 0) {
/*	imsl_e1sti (1, info);  */

    /*
     * (5, 5, "The matrix to be inverted is singular.  The index of the first
     * zero diagonal element of A is %(i1).");
     */
	imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_MATRIX);
    }
L_9000:
    imsl_e1pop ("l_linrt");
    return;
}	       /* end of function */
#undef AINV
#undef A
/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 19:13:10
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  LSLRT/DLSLRT  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Solve a real triangular system of linear equations.

    Usage:      CALL LSLRT (N, A, LDA, B, IPATH, X)

    Arguments:
       N      - Number of equations.  (Input)
       A      - N by N matrix containing the coefficient matrix for
                the triangular linear system.  (Input)
                For a lower triangular system, only the lower triangular
                part and diagonal of A are referenced.  For an upper
                triangular system, only the upper triangular part and
                diagonal of A are referenced.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       B      - Vector of length N containing the right-hand side of
                the linear system.  (Input)
       IPATH  - Path indicator.  (Input)
                IPATH = 1 means solve A*X = B, A lower triangular,
                IPATH = 2 means solve A*X = B, A upper triangular,
                IPATH = 3 means solve trans(A)*X = B, A lower triangular,
                IPATH = 4 means solve trans(A)*X = B, A upper triangular.
                   trans(A) indicates the transpose of A.
       X      - Vector of length N containing the solution to the linear
                system.  (Output)
                If B is not needed, B and X can share the same storage
                locations.

    GAMS:       D2a3

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void     l_lslrt (Mint *n, Mfloat a[], Mint *lda, Mfloat b[], Mint *ipath, Mfloat x[])
#else
static void     l_lslrt (n, a, lda, b, ipath, x)
Mint           *n;
Mfloat         *a;
Mint           *lda;
Mfloat          b[];
Mint           *ipath;
Mfloat          x[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
    Mint            i;


    imsl_e1psh ("l_lslrt");

    if (*n <= 0) {
	imsl_e1sti (1, *n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NUM_OF_EQUATIONS);
	goto L_9000;
    }
    if (*n > *lda) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *lda);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
	goto L_9000;
    }
 /* CHECK FOR ZERO DIAGONAL ELEMENTS. */
    for (i = 1; i <= *n; i++) {
	if (*A (i - 1, i - 1) == F_ZERO) {
	    imsl_e1sti (1, i-1);
	    imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_TRI_MATRIX);
	    goto L_9000;
	}
    }
 /*
  * MAKE A COPY OF B IN X AND WORK WITH X
  */
    scopy (*n, b, 1, x, 1);

    if (*ipath == 1) {
    /* SOLVE A*X=B FOR A LOWER TRIANGULAR */

	imsl_strsv ("L", "N", "N", *n, a, *lda, x, 1);
    }
    else if (*ipath == 2) {
    /* SOLVE A*X=B FOR AN UPPER TRIANGULAR. */

	imsl_strsv ("U", "N", "N", *n, a, *lda, x, 1);
    }
    else if (*ipath == 3) {
    /*
     * SOLVE TRANS(A)*X=B FOR A LOWER TRIANGULAR.
     */

	imsl_strsv ("L", "T", "N", *n, a, *lda, x, 1);
    }
    else if (*ipath == 4) {
    /*
     * SOLVE TRANS(A)*X=B FOR AN UPPER TRIANGULAR.
     */

	imsl_strsv ("U", "T", "N", *n, a, *lda, x, 1);
    }
    else {
	imsl_e1sti (1, *ipath);

    /*
     * (5, 4, "IPATH must be either 1, 2, 3 or 4 while a value of %(i1) is
     * given.");
     */
	imsl_ermes (IMSL_TERMINAL, IMSL_IPATH_RANGE_4);
    }

L_9000:
    imsl_e1pop ("l_lslrt");
    return;
}	       /* end of function */
#undef A
/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 19:16:54
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SSYR  (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 26, 1989

    Purpose:    Perform the rank-one symmetric matrix update
                A = A + alpha*x*trans(x), where A is a symmetric matrix,
                and trans(x) represents the transpose of the vector.

    Usage:      CALL SSYR (UPLO, N, ALPHA, X, INCX, A, LDA)

    Arguments:
       UPLO  -  Character indicating storage form of the matrix.  (Input)
                If UPLO is 'U' or 'u' then the matrix is stored in the
                upper half of A.  If UPLO is 'L' or 'l' then the matrix
                is stored in the lower half of A.
       N      - Order of the symmetric matrix A.  (Input)
       ALPHA  - Scalar multiplier for the vector-vector product.  (Input)
       X      - Vector of length (N-1)*IABS(INCX)+1.  (Input)
       INCX   - Displacement between elements of X.  (Input)
       A      - Symmetric matrix of order N.  (Input/Output)
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
static void     l_ssyr (Mchar *uplo, unsigned uplo_s, Mint *n, Mfloat *alpha,
                    Mfloat x[], Mint *incx, Mfloat a[], Mint *lda)
#else
static void     l_ssyr (uplo, uplo_s, n, alpha, x, incx, a, lda)
Mchar          *uplo;
unsigned        uplo_s;
Mint           *n;
Mfloat         *alpha, x[];
Mint           *incx;
Mfloat          a[];
Mint           *lda;
#endif
{
/*      LOGICAL32        imsl_l1ame(), lower, upper; */
    Mint            lower, upper;
    Mint            ix, j;


    upper = imsl_l1ame (uplo, uplo_s, "U", sizeof ("U"));
    lower = imsl_l1ame (uplo, uplo_s, "L", sizeof ("L"));
 /*
  * Test the input parameters.
  */
    if (*n < 0) {
	imsl_e1psh ("l_ssyr");
	imsl_e1sti (1, *n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
	imsl_e1pop ("l_ssyr");
	goto L_9000;
    }
    else if ((*lda < *n) || (*lda == 0)) {
	imsl_e1psh ("l_ssyr");
	imsl_e1sti (1, *lda);
	imsl_e1sti (2, *n);
	imsl_ermes (IMSL_TERMINAL, IMSL_INVALID_LDA_VALUE_GIVEN);
	imsl_e1pop ("l_ssyr");
	goto L_9000;
    }
    else if (*incx == 0) {
	imsl_e1psh ("l_ssyr");
	imsl_e1sti (1, *incx);
	imsl_ermes (IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
	imsl_e1pop ("l_ssyr");
	goto L_9000;
    }
    else if ((!upper) && (!lower)) {
	imsl_e1psh ("l_ssyr");
	imsl_e1stl (1, uplo);
	imsl_ermes (IMSL_TERMINAL, IMSL_INVALID_UPLO_VALUE);
	imsl_e1pop ("l_ssyr");
	goto L_9000;
    }
 /* Quick return if possible */
    if (*n == 0 || *alpha == F_ZERO)
	goto L_9000;

    ix = 1;
    if (*incx < 0)
	ix = (-*n + 1) ** incx + 1;

    for (j = 1; j <= *n; j++) {
	if (upper) {
	    if (*incx >= 0) {
		saxpy (j, *alpha * x[ix - 1], x, *incx, &a[*lda * (j - 1)],
		    1);
	    }
	    else {
		saxpy (j, *alpha * x[ix - 1], &x[ix - 1], *incx, &a[*lda * (j - 1)],
		    1);
	    }
	}
	else {
	    if (*incx >= 0) {
		saxpy (j, *alpha * x[ix - 1], x, *incx, &a[j - 1], *lda);
	    }
	    else {
		saxpy (j, *alpha * x[ix - 1], &x[ix - 1], *incx, &a[j - 1],
		    *lda);
	    }
	}
	ix += *incx;
    }

L_9000:
    return;
}	       /* end of function */
