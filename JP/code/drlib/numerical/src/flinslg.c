
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static VA_LIST_HACK  PROTO (l_lin_sol_gen, (Mint, Mfloat *, Mfloat *,
                                       va_list argptr));
static void     PROTO (l_l4trg, (Mint, Mfloat *, Mfloat *, Mint *));
static void     PROTO (l_l2nrg, (Mint, Mfloat *, Mint, Mfloat *,
                                 Mint, Mfloat *, Mint *));
static void     PROTO (l_nr1rr, (Mint *, Mint *, Mfloat *, Mint *,
                                 Mfloat *));
static void     PROTO (l_linrt, (Mint *n, Mfloat *a, Mint *lda,
                                 Mint *ipath, Mfloat *ainv,
                                 Mint *ldainv));
static void     PROTO (l_l2crg, (Mint *n, Mfloat *a, Mint *lda,
                                 Mfloat *imsl_fac, Mint *ldfac,
                                 Mint ipvt[], Mfloat *rcond,
	                         Mfloat z[]));

static Mfloat  *lv_x;
static Mfloat   lv_rcond;
#ifdef ANSI
Mfloat         *imsl_f_lin_sol_gen (Mint n, Mfloat *a, Mfloat *b,...)
#else
Mfloat         *imsl_f_lin_sol_gen (n, a, b, va_alist)
Mint            n;
Mfloat         *a;
Mfloat         *b;
va_dcl
#endif
{
    va_list         argptr;
    VA_START (argptr, b);

    E1PSH ("imsl_f_lin_sol_gen", "imsl_d_lin_sol_gen");
    lv_x = NULL;

    IMSL_CALL (l_lin_sol_gen (n, a, b, argptr));
    va_end (argptr);

    E1POP ("imsl_f_lin_sol_gen", "imsl_d_lin_sol_gen");
    return lv_x;
}


#ifdef ANSI
static VA_LIST_HACK  l_lin_sol_gen (Mint n, Mfloat *a, Mfloat *b,
                               va_list argptr)
#else
static VA_LIST_HACK  l_lin_sol_gen (n, a, b, argptr)
Mint            n;
Mfloat         *a;
Mfloat         *b;
va_list         argptr;
#endif
{
    Mint            code =1;
    Mint            arg_number = 3;
    Mint            path = 1;
    Mint            a_col_dim = n;
    Mfloat        **factor_ptr = NULL;
    Mfloat         *factor = NULL;
    Mfloat        **inva_ptr = NULL;
    Mfloat         *inva = NULL;
    Mfloat         *cond = NULL;
    Mint          **ipvt_ptr = NULL;
    Mint           *ipvt = NULL;
    Mint            fac_col_dim = n;
    Mint            inv_col_dim = n;

    Mfloat         *fwork = NULL;
    Mint           *iwork = NULL;

    Mint            factor_user = 0;
    Mint            inverse_user = 0;
    Mint            result_user = 0;
    Mint            return_factor = 0;
    Mint            return_inverse = 0;
    Mint            return_cond = 0;
    Mint            inverse_only = 0;
    Mint            factor_only = 0;
    Mint            solve_only = 0;
    Mint            error = 0;
    Mint            lda = n;
    Mint            i;
    Mfloat          rcond;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	    case IMSL_RETURN_USER:
		lv_x = va_arg (argptr, Mfloat *);
		if (!lv_x) {
		    imsl_e1stl (1, "x");
		    imsl_e1stl (2, "IMSL_RETURN_USER");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		}
		result_user = 1;
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
		factor_user = 0;
		return_factor = 1;
		ipvt_ptr = va_arg (argptr, Mint **);
		factor_ptr = va_arg (argptr, Mfloat **);
		arg_number += 2;
		break;
	    case IMSL_FACTOR_USER:
		factor_user = 1;
		return_factor = 1;
		ipvt = va_arg (argptr, Mint *);
		if (!ipvt) {
		    imsl_e1stl (1, "ipvt");
		    imsl_e1stl (2, "IMSL_FACTOR");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		}
		factor = va_arg (argptr, Mfloat *);
		if (!factor) {
		    imsl_e1stl (1, "factor");
		    imsl_e1stl (2, "IMSL_FACTOR");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_2);
		}
		arg_number += 2;
		break;
	    case IMSL_FAC_COL_DIM:
		fac_col_dim = va_arg (argptr, Mint);
		arg_number++;
		break;
	    case IMSL_INVERSE:
		inverse_user = 0;
		return_inverse = 1;
		inva_ptr = va_arg (argptr, Mfloat **);
		arg_number++;
		break;
	    case IMSL_INVERSE_USER:
		inverse_user = 1;
		return_inverse = 1;
		inva = va_arg (argptr, Mfloat *);
		if (!inva) {
		    imsl_e1stl (1, "inva");
		    imsl_e1stl (2, "IMSL_INVERSE_USER");
		    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
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
		}
		break;
	    case IMSL_INVERSE_ONLY:
		inverse_only = 1;
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
		imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
		break;
	}
    }

    if (!a && !solve_only) {
	imsl_e1stl (1, "a");
	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
    }
    if (!b && !(factor_only || inverse_only)) {
	imsl_e1stl (1, "b");
	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
    }

    if (imsl_n1rty (0))
	goto RETURN;

    if (n <= 0) {
	imsl_e1sti (1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
    }
    else {
	if (n > a_col_dim) {
	    imsl_e1sti (1, n);
	    imsl_e1sti (2, a_col_dim);
	    imsl_e1stl (1, "a");
	    imsl_ermes (IMSL_TERMINAL, IMSL_COL_DIM_LESS_ORDER);
	}
	if (n > fac_col_dim) {
	    imsl_e1sti (1, n);
	    imsl_e1sti (2, fac_col_dim);
	    imsl_e1stl (1, "factor");
	    imsl_ermes (IMSL_TERMINAL, IMSL_COL_DIM_LESS_ORDER);
	}
    }
    if (solve_only + factor_only + inverse_only > 1) {
        imsl_ermes (IMSL_TERMINAL, IMSL_BAD_SOLVE_FACTOR_INVERSE);
        error = 1;
    }
    if (imsl_n1rty (0))
	goto RETURN;

    if (solve_only) {
	if (!factor_user) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_SPECIFY_SOLVE_ONLY);
	    error = 1;
	}
	if (return_cond) {
	    imsl_ermes (IMSL_TERMINAL, IMSL_CONDITION_ONLY_SPECIFIER);
	    error = 1;
	}
	if (error)
	    goto RETURN;
    }
    if (inverse_only && !return_inverse) {
	imsl_ermes (IMSL_TERMINAL, IMSL_INVERSE_ONLY_SPECIFIER);
	goto RETURN;
    }
    if (factor_only && !return_factor) {
	imsl_ermes (IMSL_TERMINAL, IMSL_FACTOR_ONLY_SPECIFIER);
	goto RETURN;
    }
    if (!inverse_only) {
	if (!solve_only) {
	    Mint            ldfac = n;
	    if (!factor_user) {
		factor = (Mfloat *) imsl_malloc (n * fac_col_dim * sizeof (Mfloat));
		ipvt = (Mint *) imsl_malloc (n * sizeof (Mint));
	    }
	    else if (fac_col_dim > n) {
		imsl_f_m1ran (n, fac_col_dim, factor, factor);
		if (error = (imsl_n1rty (1) > 3) ? 1 : 0)
		    goto FREE_SPACE;
	    }
	    fwork = (Mfloat *) imsl_malloc (n * sizeof (Mfloat));
	    if (!fwork) {
		imsl_e1stl (1, "n");
		imsl_e1sti (1, n);
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
		error = 1;
		goto FREE_SPACE;
	    }
	    imsl_f_m1ran (n, a_col_dim, a, a);
	    if (!(error = imsl_n1rty (1) > 3)) {
		if (factor_only && !return_cond) {
		    imsl_l2trg (n, a, lda, factor, ldfac, ipvt, fwork);
		    error = imsl_n1rty (1) > 3 ? 1 : 0;
		}
		else {
		    l_l2crg (&n, a, &lda, factor, &ldfac, ipvt, &rcond, fwork);
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
		imsl_f_m1ran (a_col_dim, n, a, a);
	    }
	    if (fwork) imsl_free (fwork);
	    fwork = NULL;
	}
	if (!factor_only && !error) {
	    if (!result_user) {
		lv_x = (Mfloat *) imsl_malloc (n * sizeof (Mfloat));
		if (lv_x == NULL) {
		    imsl_e1stl (1, "n");
		    imsl_e1sti (1, n);
		    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
		    error = 1;
		}
	    }
	    if (solve_only){
	      imsl_f_m1ran(n, fac_col_dim, factor, factor);
	      error  = imsl_n1rty(1) > 3 ? 1 : 0;
	    }
	    if (!error) {
		imsl_lfsrg (n, factor, n, ipvt, b, &path, lv_x);
		error = imsl_n1rty (1) > 3 ? 1 : 0;
	    }
	}
    }
    if (return_inverse && !error) {
	if (!inverse_user) {
	    inva = (Mfloat *) imsl_malloc (n * inv_col_dim * sizeof (Mfloat));
	    if (!inva) {
		imsl_e1stl (1, "n");
		imsl_e1sti (1, n);
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
		error = 1;
		goto FREE_SPACE;
	    }
	}
	fwork = (Mfloat *) imsl_malloc ((n + n * (n - 1) / 2) * sizeof (Mfloat));
	iwork = (Mint *) imsl_malloc (n * sizeof (Mint));
	if (!iwork || !fwork) {
	    imsl_e1stl (1, "n");
	    imsl_e1sti (1, n);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    error = 1;
	}
	if (!error) {
	    Mint            ldinv;
	    if (inverse_only) {
		imsl_f_m1ran (n, a_col_dim, a, a);
		if (!(error = (imsl_n1rty (1) > 3) ? 1 : 0)) {
		    ldinv = n;
		    if ((inv_col_dim > n) && inverse_user) {
			imsl_f_m1ran (n, inv_col_dim, inva, inva);
			error = (imsl_n1rty (1) > 3) ? 1 : 0;
		    }
		}
	    }
	    else {
		ldinv = inv_col_dim;
		lda = a_col_dim;
	    }
	    if (!error) {
		l_l2nrg (n, a, lda, inva, ldinv, fwork, iwork);
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
		    imsl_f_m1ran (inv_col_dim, n, inva, inva);
		    imsl_f_m1ran (a_col_dim, n, a, a);
		}
	    }
	}
    }

FREE_SPACE:

    if (!factor_user) {
	if (error || !return_factor) {
	    if (factor)	imsl_free (factor);
	    if (ipvt)	imsl_free (ipvt);
	    factor = NULL;
	}
	else {
	    imsl_f_m1ran (fac_col_dim, n, factor, factor);
	    for (i = n; i < fac_col_dim; i++) {
		sset (n, F_ZERO, (factor + i), fac_col_dim);
	    }
	    *factor_ptr = factor;
	    *ipvt_ptr = ipvt;
	}
	factor = NULL;
	ipvt = NULL;
    }
    else {
	imsl_f_m1ran (fac_col_dim, n, factor, factor);
    }
    if (!inverse_user) {
	if (error || !return_inverse) {
	    if (inva) imsl_free (inva);
	}
	else {
	    for (i = n; i < inv_col_dim; i++) {
		sset (n, F_ZERO, (inva + i), inv_col_dim);
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
    if (fwork) imsl_free (fwork);
    if (iwork) imsl_free (iwork);
RETURN:
    return (argptr);
}


/* Structured by FOR_STRUCT, v0.2, on 07/23/90 at 10:45:39
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  L2TRG/DL2TRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 25, 1989

    Purpose:    Compute the LU factorization of a real general matrix.

    Usage:      CALL L2TRG (N, A, LDA, FAC, LDFAC, IPVT, SCALE)

    Arguments:  See LFTRG/DLFTRG.

    Remarks:    See LFTRG/DLFTRG.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------

 */
#ifdef ANSI
void            imsl_l2trg (Mint n, Mfloat *a, Mint lda, Mfloat *fac, Mint ldfac,
                    Mint *ipvt, Mfloat *scale)
#else
void            imsl_l2trg (n, a, lda, fac, ldfac, ipvt, scale)
Mint            n;
Mfloat         *a;
Mint            lda;
Mfloat         *fac;
Mint            ldfac, ipvt[];
Mfloat          scale[];
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
#define FAC(I_,J_)	(fac+(I_)*(ldfac)+(J_))
    Mint            i, indj, info, j, k, ktemp, m0, m1, m2, m3, m4, m5, m6,
                    m7;
    Mfloat          big, small, t, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    static Mfloat   fpm1 = -1.0e0;
 /*
  * this code is for computer types: aliant and necv
  */
    E1PSH ("imsl_l2trg", "imsl_dl2trg");

    if (n <= 0) {
	imsl_e1sti (1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
	goto L_9000;
    }
    if (n > lda) {
	imsl_e1sti (1, n);
	imsl_e1sti (2, lda);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
	goto L_9000;
    }
    if (n > ldfac) {
	imsl_e1sti (1, n);
	imsl_e1sti (2, ldfac);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDFAC_LESS_ORDER);
	goto L_9000;
    }
 /* Preserve a copy of the input matrix */
    imsl_crgrg (n, a, lda, fac, ldfac);

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
    for (i = 1; i <= n; i++) {
	indj = imsl_isamax (n, FAC (0, i - 1), ldfac);
	scale[i - 1] = fabs (*FAC (indj - 1, i - 1));
	if (scale[i - 1] < small) {
	    scale[i - 1] = F_ONE;
	}
	else {
	    scale[i - 1] = F_ONE / scale[i - 1];
	}
    }
 /* begin factorization step */
L_20:
    if (j >= n)
	goto L_170;
 /*
  * determine the pivot element for column j of L
  */
    l_l4trg (n - j + 1, FAC (j - 1, j - 1), &scale[j - 1], &m0);
    m0 += j - 1;
    ipvt[j - 1] = m0;
    if (fabs (*FAC (j - 1, m0 - 1)) > small) {
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
	t = fpm1 / t;
    /* update U(j,j+1) */
	t1 = *FAC (j, m0 - 1);
	*FAC (j, m0 - 1) = *FAC (j, j - 1);
	*FAC (j, j - 1) = t1;
    /* update columns j of L and j+1 of A */
	for (i = j + 1; i <= n; i++) {
	    *FAC (j - 1, i - 1) *= t;
	    *FAC (j, i - 1) += t1 ** FAC (j - 1, i - 1);
	}
    }
    else {
	ktemp = 0;
	info = m0;
    }

    if (j + 1 >= n)
	goto L_170;
 /*
  * determine the pivot element for column j+1
  */
    l_l4trg (n - j, FAC (j, j), &scale[j], &m1);
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
    if (fabs (t) > small) {
	t = fpm1 / t;
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
    t2 = *FAC (j + 1, m1 - 1) + t1 ** FAC (j - 1, j);
    *FAC (j + 1, m1 - 1) = *FAC (j + 1, j);
    *FAC (j + 1, j) = t2;
 /* update column j+1 of L and j+2 of A */
    for (i = j + 2; i <= n; i++) {
	*FAC (j, i - 1) *= t;
	*FAC (j + 1, i - 1) += t1 ** FAC (j - 1, i - 1) + t2 ** FAC (j, i - 1);
    }

    if (j + 2 >= n)
	goto L_160;
 /*
  * determine the pivot element for column j+2
  */
    l_l4trg (n - j - 1, FAC (j + 1, j + 1), &scale[j + 1], &m2);
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
    if (fabs (t) > small) {
	t = fpm1 / t;
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
    t2 = *FAC (j + 2, m1 - 1) + t1 ** FAC (j - 1, j);
    *FAC (j + 2, m1 - 1) = *FAC (j + 2, j);
    *FAC (j + 2, j) = t2;
    t3 = *FAC (j + 2, m2 - 1) + t1 ** FAC (j - 1, j + 1) + t2 *
	*FAC (j, j + 1);
    *FAC (j + 2, m2 - 1) = *FAC (j + 2, j + 1);
    *FAC (j + 2, j + 1) = t3;
 /* update column j+2 of L and j+3 of A */
    for (i = j + 3; i <= n; i++) {
	*FAC (j + 1, i - 1) *= t;
	*FAC (j + 2, i - 1) += t1 ** FAC (j - 1, i - 1) + t2 ** FAC (j, i - 1) +
	    t3 ** FAC (j + 1, i - 1);
    }

    if (j + 3 >= n)
	goto L_150;
 /*
  * determine the pivot element for column j+3
  */
    l_l4trg (n - j - 2, FAC (j + 2, j + 2), &scale[j + 2], &m3);
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
    if (fabs (t) > small) {
	t = fpm1 / t;
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
    t2 = *FAC (j + 3, m1 - 1) + t1 ** FAC (j - 1, j);
    *FAC (j + 3, m1 - 1) = *FAC (j + 3, j);
    *FAC (j + 3, j) = t2;
    t3 = *FAC (j + 3, m2 - 1) + t1 ** FAC (j - 1, j + 1) + t2 *
	*FAC (j, j + 1);
    *FAC (j + 3, m2 - 1) = *FAC (j + 3, j + 1);
    *FAC (j + 3, j + 1) = t3;
    t4 = *FAC (j + 3, m3 - 1) + t1 ** FAC (j - 1, j + 2) + t2 *
	*FAC (j, j + 2) + t3 ** FAC (j + 1, j + 2);
    *FAC (j + 3, m3 - 1) = *FAC (j + 3, j + 2);
    *FAC (j + 3, j + 2) = t4;
 /* update column j+3 of L and j+4 of A */
    for (i = j + 4; i <= n; i++) {
	*FAC (j + 2, i - 1) *= t;
	*FAC (j + 3, i - 1) += t1 ** FAC (j - 1, i - 1) + t2 ** FAC (j, i - 1) +
	    t3 ** FAC (j + 1, i - 1) + t4 ** FAC (j + 2, i - 1);
    }

    if (j + 4 >= n)
	goto L_140;
 /*
  * determine the pivot element for column j+4
  */
    l_l4trg (n - j - 3, FAC (j + 3, j + 3), &scale[j + 3], &m4);
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
    if (fabs (t) > small) {
	t = fpm1 / t;
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
    t2 = *FAC (j + 4, m1 - 1) + t1 ** FAC (j - 1, j);
    *FAC (j + 4, m1 - 1) = *FAC (j + 4, j);
    *FAC (j + 4, j) = t2;
    t3 = *FAC (j + 4, m2 - 1) + t1 ** FAC (j - 1, j + 1) + t2 *
	*FAC (j, j + 1);
    *FAC (j + 4, m2 - 1) = *FAC (j + 4, j + 1);
    *FAC (j + 4, j + 1) = t3;
    t4 = *FAC (j + 4, m3 - 1) + t1 ** FAC (j - 1, j + 2) + t2 *
	*FAC (j, j + 2) + t3 ** FAC (j + 1, j + 2);
    *FAC (j + 4, m3 - 1) = *FAC (j + 4, j + 2);
    *FAC (j + 4, j + 2) = t4;
    t5 = *FAC (j + 4, m4 - 1) + t1 ** FAC (j - 1, j + 3) + t2 *
	*FAC (j, j + 3) + t3 ** FAC (j + 1, j + 3) + t4 ** FAC (j + 2, j + 3);
    *FAC (j + 4, m4 - 1) = *FAC (j + 4, j + 3);
    *FAC (j + 4, j + 3) = t5;
 /* update column j+4 of L and j+5 of A */
    for (i = j + 5; i <= n; i++) {
	*FAC (j + 3, i - 1) *= t;
	*FAC (j + 4, i - 1) += t1 ** FAC (j - 1, i - 1) + t2 ** FAC (j, i - 1) +
	    t3 ** FAC (j + 1, i - 1) + t4 ** FAC (j + 2, i - 1) + t5 *
	    *FAC (j + 3, i - 1);
    }

    if (j + 5 >= n)
	goto L_130;
 /*
  * determine the pivot element for column j+5
  */
    l_l4trg (n - j - 4, FAC (j + 4, j + 4), &scale[j + 4], &m5);
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
    if (fabs (t) > small) {
	t = fpm1 / t;
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
    t2 = *FAC (j + 5, m1 - 1) + t1 ** FAC (j - 1, j);
    *FAC (j + 5, m1 - 1) = *FAC (j + 5, j);
    *FAC (j + 5, j) = t2;
    t3 = *FAC (j + 5, m2 - 1) + t1 ** FAC (j - 1, j + 1) + t2 *
	*FAC (j, j + 1);
    *FAC (j + 5, m2 - 1) = *FAC (j + 5, j + 1);
    *FAC (j + 5, j + 1) = t3;
    t4 = *FAC (j + 5, m3 - 1) + t1 ** FAC (j - 1, j + 2) + t2 *
	*FAC (j, j + 2) + t3 ** FAC (j + 1, j + 2);
    *FAC (j + 5, m3 - 1) = *FAC (j + 5, j + 2);
    *FAC (j + 5, j + 2) = t4;
    t5 = *FAC (j + 5, m4 - 1) + t1 ** FAC (j - 1, j + 3) + t2 *
	*FAC (j, j + 3) + t3 ** FAC (j + 1, j + 3) + t4 ** FAC (j + 2, j + 3);
    *FAC (j + 5, m4 - 1) = *FAC (j + 5, j + 3);
    *FAC (j + 5, j + 3) = t5;
    t6 = *FAC (j + 5, m5 - 1) + t1 ** FAC (j - 1, j + 4) + t2 *
	*FAC (j, j + 4) + t3 ** FAC (j + 1, j + 4) + t4 ** FAC (j + 2, j + 4) +
	t5 ** FAC (j + 3, j + 4);
    *FAC (j + 5, m5 - 1) = *FAC (j + 5, j + 4);
    *FAC (j + 5, j + 4) = t6;
 /* update column j+5 of L and j+6 of A */
    for (i = j + 6; i <= n; i++) {
	*FAC (j + 4, i - 1) *= t;
	*FAC (j + 5, i - 1) += t1 ** FAC (j - 1, i - 1) + t2 ** FAC (j, i - 1) +
	    t3 ** FAC (j + 1, i - 1) + t4 ** FAC (j + 2, i - 1) + t5 *
	    *FAC (j + 3, i - 1) + t6 ** FAC (j + 4, i - 1);
    }

    if (j + 6 >= n)
	goto L_120;
 /*
  * determine the pivot element for column j+6
  */
    l_l4trg (n - j - 5, FAC (j + 5, j + 5), &scale[j + 5], &m6);
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
    if (fabs (t) > small) {
	t = fpm1 / t;
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
    t2 = *FAC (j + 6, m1 - 1) + t1 ** FAC (j - 1, j);
    *FAC (j + 6, m1 - 1) = *FAC (j + 6, j);
    *FAC (j + 6, j) = t2;
    t3 = *FAC (j + 6, m2 - 1) + t1 ** FAC (j - 1, j + 1) + t2 *
	*FAC (j, j + 1);
    *FAC (j + 6, m2 - 1) = *FAC (j + 6, j + 1);
    *FAC (j + 6, j + 1) = t3;
    t4 = *FAC (j + 6, m3 - 1) + t1 ** FAC (j - 1, j + 2) + t2 *
	*FAC (j, j + 2) + t3 ** FAC (j + 1, j + 2);
    *FAC (j + 6, m3 - 1) = *FAC (j + 6, j + 2);
    *FAC (j + 6, j + 2) = t4;
    t5 = *FAC (j + 6, m4 - 1) + t1 ** FAC (j - 1, j + 3) + t2 *
	*FAC (j, j + 3) + t3 ** FAC (j + 1, j + 3) + t4 ** FAC (j + 2, j + 3);
    *FAC (j + 6, m4 - 1) = *FAC (j + 6, j + 3);
    *FAC (j + 6, j + 3) = t5;
    t6 = *FAC (j + 6, m5 - 1) + t1 ** FAC (j - 1, j + 4) + t2 *
	*FAC (j, j + 4) + t3 ** FAC (j + 1, j + 4) + t4 ** FAC (j + 2, j + 4) +
	t5 ** FAC (j + 3, j + 4);
    *FAC (j + 6, m5 - 1) = *FAC (j + 6, j + 4);
    *FAC (j + 6, j + 4) = t6;
    t7 = *FAC (j + 6, m6 - 1) + t1 ** FAC (j - 1, j + 5) + t2 *
	*FAC (j, j + 5) + t3 ** FAC (j + 1, j + 5) + t4 ** FAC (j + 2, j + 5) +
	t5 ** FAC (j + 3, j + 5) + t6 ** FAC (j + 4, j + 5);
    *FAC (j + 6, m6 - 1) = *FAC (j + 6, j + 5);
    *FAC (j + 6, j + 5) = t7;
 /* update column j+6 of L and j+7 of A */
    for (i = j + 7; i <= n; i++) {
	*FAC (j + 5, i - 1) *= t;
	*FAC (j + 6, i - 1) += t1 ** FAC (j - 1, i - 1) + t2 ** FAC (j, i - 1) +
	    t3 ** FAC (j + 1, i - 1) + t4 ** FAC (j + 2, i - 1) + t5 *
	    *FAC (j + 3, i - 1) + t6 ** FAC (j + 4, i - 1) + t7 ** FAC (j + 5, i - 1);
    }

    if (j + 7 >= n)
	goto L_110;
 /*
  * determine the pivot element for column j+7
  */
    l_l4trg (n - j - 6, FAC (j + 6, j + 6), &scale[j + 6], &m7);
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
    if (fabs (t) > small) {
	t = fpm1 / t;
	for (i = j + 8; i <= n; i++) {
	    *FAC (j + 6, i - 1) *= t;
	}
	if (m7 != j + 7)
	    ktemp = -ktemp;
    }
    else {
	ktemp = 0;
	info = m7;
    }

 /*
  * update the remaining rectangular block of U, rows j to j+7 and columns
  * j+8 to n
  */
    for (k = n; k >= (j + 8); k--) {
	t1 = *FAC (k - 1, m0 - 1);
	*FAC (k - 1, m0 - 1) = *FAC (k - 1, j - 1);
	*FAC (k - 1, j - 1) = t1;
	t2 = *FAC (k - 1, m1 - 1) + t1 ** FAC (j - 1, j);
	*FAC (k - 1, m1 - 1) = *FAC (k - 1, j);
	*FAC (k - 1, j) = t2;
	t3 = *FAC (k - 1, m2 - 1) + t1 ** FAC (j - 1, j + 1) +
	    t2 ** FAC (j, j + 1);
	*FAC (k - 1, m2 - 1) = *FAC (k - 1, j + 1);
	*FAC (k - 1, j + 1) = t3;
	t4 = *FAC (k - 1, m3 - 1) + t1 ** FAC (j - 1, j + 2) +
	    t2 ** FAC (j, j + 2) + t3 ** FAC (j + 1, j + 2);
	*FAC (k - 1, m3 - 1) = *FAC (k - 1, j + 2);
	*FAC (k - 1, j + 2) = t4;
	t5 = *FAC (k - 1, m4 - 1) + t1 ** FAC (j - 1, j + 3) +
	    t2 ** FAC (j, j + 3) + t3 ** FAC (j + 1, j + 3) + t4 ** FAC (j + 2, j + 3);
	*FAC (k - 1, m4 - 1) = *FAC (k - 1, j + 3);
	*FAC (k - 1, j + 3) = t5;
	t6 = *FAC (k - 1, m5 - 1) + t1 ** FAC (j - 1, j + 4) +
	    t2 ** FAC (j, j + 4) + t3 ** FAC (j + 1, j + 4) + t4 ** FAC (j + 2, j + 4) +
	    t5 ** FAC (j + 3, j + 4);
	*FAC (k - 1, m5 - 1) = *FAC (k - 1, j + 4);
	*FAC (k - 1, j + 4) = t6;
	t7 = *FAC (k - 1, m6 - 1) + t1 ** FAC (j - 1, j + 5) +
	    t2 ** FAC (j, j + 5) + t3 ** FAC (j + 1, j + 5) + t4 ** FAC (j + 2, j + 5) +
	    t5 ** FAC (j + 3, j + 5) + t6 ** FAC (j + 4, j + 5);
	*FAC (k - 1, m6 - 1) = *FAC (k - 1, j + 5);
	*FAC (k - 1, j + 5) = t7;
	t8 = *FAC (k - 1, m7 - 1) + t1 ** FAC (j - 1, j + 6) +
	    t2 ** FAC (j, j + 6) + t3 ** FAC (j + 1, j + 6) + t4 ** FAC (j + 2, j + 6) +
	    t5 ** FAC (j + 3, j + 6) + t6 ** FAC (j + 4, j + 6) + t7 *
	    *FAC (j + 5, j + 6);
	*FAC (k - 1, m7 - 1) = *FAC (k - 1, j + 6);
	*FAC (k - 1, j + 6) = t8;
    /*
     * rank 8 update of the lower right block from rows j+8 to n and columns
     * j+8 to n
     */
	for (i = j + 8; i <= n; i++) {
	    *FAC (k - 1, i - 1) += t1 ** FAC (j - 1, i - 1) +
		t2 ** FAC (j, i - 1) + t3 ** FAC (j + 1, i - 1) + t4 *
		*FAC (j + 2, i - 1) + t5 ** FAC (j + 3, i - 1) +
		t6 ** FAC (j + 4, i - 1) + t7 ** FAC (j + 5, i - 1) +
		t8 ** FAC (j + 6, i - 1);
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
    ipvt[n - 1] = n;
    if (fabs (*FAC (n - 1, n - 1)) <= small)
	info = n;

    if (info != 0) {
	imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_MATRIX);
    }
L_9000:
    E1POP ("imsl_l2trg", "imsl_dl2trg");
    return;
}	       /* end of function */
#undef A
#undef FAC


/*----------------------------------------------------------------------- */
/*  IMSL Name:  L4TRG/DL4TRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 29, 1988

    Purpose:    Find the pivot element using scaling.

    Usage:      CALL L4TRG (N, X, SCALE, K)

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
static void     l_l4trg (Mint n, Mfloat *x, Mfloat *scale, Mint *k)
#else
static void     l_l4trg (n, x, scale, k)
Mint            n;
Mfloat          x[], scale[];
Mint           *k;
#endif
{
    Mint            j, j1, j2, jmax;
    Mfloat          curmax, temp[NTEMP], value;
 /*
  * Code for VECTOR machines with coded ISAMAX
  */
    *k = 1;
    curmax = F_ZERO;
    for (j1 = 1; j1 <= n; j1 += NTEMP) {
	j2 = imsl_i_min (j1 + NTEMP - 1, n);
	for (j = j1; j <= j2; j++) {
	    temp[j - j1] = x[j - 1] * scale[j - 1];
	}
	jmax = imsl_isamax (j2 - j1 + 1, temp, 1);
	value = fabs (temp[jmax - 1]);
	if (value > curmax) {
	    curmax = value;
	    *k = jmax + j1 - 1;
	}
    }

    return;
}
#undef NTEMP


/* Structured by FOR_STRUCT, v0.2, on 01/08/90 at 17:22:58
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  LFSRG/DLFSRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 27, 1985

    Purpose:    Solve a real general system of linear equations given the
                LU factorization of the coefficient matrix.

    Usage:      CALL LFSRG (N, FAC, LDFAC, IPVT, B, IPATH, X)

    Arguments:
       N      - Number of equations.  (Input)
       FAC    - N by N matrix containing the LU factorization of the
                coefficient matrix A as output from subroutine
                LFCRG/DLFCRG or LFTRG/DLFTRG.  (Input)
       LDFAC  - Leading dimension of FAC exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPVT   - Vector of length N containing the pivoting information
                for the LU factorization of A as output from subroutine
                LFCRG/DLFCRG or LFTRG/DLFTRG.  (Input)
       B      - Vector of length N containing the right-hand side of the
                linear system.  (Input)
       IPATH  - Path indicator.  (Input)
                IPATH = 1 means the system A*X = B is solved.
                IPATH = 2 means the system trans(A)*X = B is solved,
                          where trans(A) is the transpose of A.
       X      - Vector of length N containing the solution to the linear
                system.  (Output)
                If B is not needed, B and X can share the same storage
                locations.

    GAMS:       D2a1

    Chapters:   MATH/LIBRARY Linear Systems
                STAT/LIBRARY Mathematical Support

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void            imsl_lfsrg (Mint n, Mfloat *fac, Mint ldfac, Mint *ipvt, Mfloat *b,
                    Mint *ipath, Mfloat *x)
#else
void            imsl_lfsrg (n, fac, ldfac, ipvt, b, ipath, x)
Mint            n;
Mfloat         *fac;
Mint            ldfac, ipvt[];
Mfloat          b[];
Mint           *ipath;
Mfloat          x[];
#endif
{
#define FAC(I_,J_)	(fac+(I_)*(ldfac)+(J_))
    Mint            k, l;
    Mfloat          t;


    E1PSH ("imsl_lfsrg", "imsl_dlfsrg");

    if (n <= 0) {
	imsl_e1sti (1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
	goto L_9000;
    }
    if (n > ldfac) {
	imsl_e1sti (1, n);
	imsl_e1sti (2, ldfac);
	imsl_ermes (IMSL_TERMINAL, IMSL_LDFAC_LESS_ORDER);
	goto L_9000;
    }
 /*
  * COPY B INTO X AND USE X TO PRESERVE INPUT
  */
    scopy (n, b, 1, x, 1);

    if (*ipath == 1) {
    /*
     * IPATH = 1 , SOLVE  A * X = B FIRST SOLVE  L*Y = B
     */
	for (k = 1; k <= (n - 1); k++) {
	    l = ipvt[k - 1];
	    t = x[l - 1];
	    if (l != k) {
		x[l - 1] = x[k - 1];
		x[k - 1] = t;
	    }
	    saxpy (n - k, t, FAC (k - 1, k), 1, &x[k], 1);
	}
    /* NOW SOLVE  U*X = Y */
	for (k = n; k >= 1; k--) {
	    if (fabs (*FAC (k - 1, k - 1)) == F_ZERO) {
		imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_MATRIX);
		goto L_9000;
	    }
	}
	imsl_strsv ("U", "N", "N", n, fac, ldfac, x, 1);
    }
    else if (*ipath == 2) {
    /*
     * IPATH = 2, SOLVE  TRANS(A) * X = B FIRST SOLVE  TRANS(U)*Y = B
     */
	for (k = 1; k <= n; k++) {
	    if (fabs (*FAC (k - 1, k - 1)) == F_ZERO) {
		imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_MATRIX);
		goto L_9000;
	    }
	}

	imsl_strsv ("U", "T", "N", n, fac, ldfac, x, 1);
    /* NOW SOLVE TRANS(L)*X = Y */
	for (k = n - 1; k >= 1; k--) {
	    x[k - 1] += imsl_sdot (n - k, FAC (k - 1, k), 1, &x[k], 1);
	    l = ipvt[k - 1];
	    if (l != k) {
		t = x[l - 1];
		x[l - 1] = x[k - 1];
		x[k - 1] = t;
	    }
	}
    }
    else {
	imsl_e1sti (1, 1);
	imsl_e1sti (2, 2);
	imsl_e1sti (3, *ipath);
	imsl_ermes (IMSL_TERMINAL, IMSL_IPATH_RANGE);
    }

L_9000:
    E1POP ("imsl_lfsrg", "imsl_dlfsrg");
    return;
}	       /* end of function */
#undef FAC

/* Structured by FOR_STRUCT, v0.2, on 09/14/90 at 10:21:55
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  L2NRG/DL2NRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 27, 1985

    Purpose:    Compute the inverse of a real general matrix.

    Usage:      CALL L2NRG (N, A, LDA, AINV, LDAINV, WK, IWK)

    Arguments:  See LINRG/DLINRG.

    Remarks:    See LINRG/DLINRG.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void     l_l2nrg (Mint n, Mfloat *a, Mint lda, Mfloat *ainv,
                    Mint ldainv, Mfloat *wk, Mint *iwk)
#else
static void     l_l2nrg (n, a, lda, ainv, ldainv, wk, iwk)
Mint            n;
Mfloat         *a;
Mint            lda;
Mfloat         *ainv;
Mint            ldainv;
Mfloat          wk[];
Mint            iwk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
#define AINV(I_,J_)	(ainv+(I_)*(ldainv)+(J_))
    Mint            _l0, _l1, i, inc, j, k, l;
    Mfloat          _f0;
    imsl_e1psh ("L2NRG ");

    if (n <= 0) {
	imsl_e1sti (1, n);

    /*
     * (5, 1, "The order of the matrix must be positive while N = %(i1) is
     * given.");
     */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
	goto L_9000;
    }
    if (n > lda) {
	imsl_e1sti (1, n);
	imsl_e1sti (2, lda);

    /*
     * (5, 2, "The order of the matrix must be less than or equal to its
     * leading dimension while N = %(i1) and LDA = %(i2) are given.");
     */
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
	goto L_9000;
    }
    if (n > ldainv) {
	imsl_e1sti (1, n);
	imsl_e1sti (2, ldainv);

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
    inc = n * (n - 1) / 2;
    l_l2crg (&n, a, &lda, ainv, &ldainv, iwk, &lv_rcond, &wk[inc]);
    if (imsl_n1rty (1) == 4)
	goto L_9000;
 /* COMPUTE INVERSE(U) */
    j = inc;
    k = 0;
    for (i = 1; i <= (n - 1); i++) {
	j -= k;
	scopy (i, AINV (n - i - 1, n - i), 1, &wk[j - 1], 1);
	k = i + 1;
    }

    _l0 = 2;
    l_linrt (&n, ainv, &ldainv, &_l0, ainv, &ldainv);

    j = inc;
    k = 0;
    for (i = 1; i <= (n - 1); i++) {
	j -= k;
	scopy (i, &wk[j - 1], 1, AINV (n - i - 1, n - i), 1);
	k = i + 1;
    }
 /* FORM INVERSE(U)*INVERSE(L) */
    for (k = n - 1; k >= 1; k--) {
	scopy (n - k, AINV (k - 1, k), 1, &wk[k + inc], 1);
	sset (n - k, F_ZERO, AINV (k - 1, k), 1);
	_l0 = n - k;
	_f0 = F_ONE;
	_l1 = 1;
	imsl_sgemv ("N", sizeof ("N"), &n, &_l0, &_f0,
	    AINV (k, 0), &ldainv, &wk[k + inc], &_l1, &_f0,
	    AINV (k - 1, 0), &_l1);
	l = iwk[k - 1];
	if (l != k)
	    sswap (n, AINV (k - 1, 0), 1, AINV (l - 1, 0), 1);
    }

    if (lv_rcond <= imsl_amach (4)) {
	imsl_e1str (1, lv_rcond);

    /*
     * (3, 1, "The matrix is too ill-conditioned. An estimate of the
     * reciprocal of its L1 condition number is RCOND = %(r1).  The inverse
     * might not be accurate.");
     */
	imsl_ermes (IMSL_WARNING, IMSL_ILL_CONDITIONED);
    }
L_9000:
    imsl_e1pop ("L2NRG ");
    return;
}	       /* end of function */

#undef A
#undef AINV


/* Structured by FOR_STRUCT, v0.2, on 09/14/90 at 11:45:48
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  L2CRG/DL2CRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 27, 1985

    Purpose:    Compute the LU factorization of a real general matrix and
                estimate its L1 condition number.

    Usage:      CALL L2CRG (N, A, LDA, FAC, LDFAC, IPVT, RCOND, Z)

    Arguments:  See LFCRG/DLFCRG.

    Remarks:    See LFCRG/DLFCRG.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void     l_l2crg (Mint *n, Mfloat *a, Mint *lda, Mfloat *fac,
                    Mint *ldfac, Mint ipvt[], Mfloat *rcond,
                    Mfloat z[])
#else
static void     l_l2crg (n, a, lda, fac, ldfac, ipvt, rcond, z)
Mint           *n;
Mfloat         *a;
Mint           *lda;
Mfloat         *fac;
Mint           *ldfac, ipvt[];
Mfloat         *rcond, z[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define FAC(I_,J_)	(fac+(I_)*(*ldfac)+(J_))
    Mint            _l0, j, k, kp1, l;
    Mfloat          anorm, ek, s, sm, t, wk, wkm, ynorm;


    imsl_e1psh ("L2CRG ");

    if (*n <= 0) {
	imsl_e1sti (1, *n);

    /*
     * (5, 1, "The order of the matrix must be positive while N = %(i1) is
     * given.");
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
    *rcond = F_ZERO;
 /* COMPUTE 1-NORM OF A */
    l_nr1rr (n, n, a, lda, &anorm);
 /*
  * FACTORIZATION STEP
  */
    imsl_l2trg (*n, a, *lda, fac, *ldfac, ipvt, z);
    if (imsl_n1rty (1) == 4)
	goto L_9000;
 /*
  * RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))). ESTIMATE =
  * NORM(Z)/NORM(Y) WHERE A*Z = Y AND TRANS(A)*Y = E . TRANS(A) IS THE
  * TRANSPOSE OF A. THE COMPONENTS OF E ARE CHOSEN TO CAUSE MAXIMUM LO- CAL
  * GROWTH IN THE ELEMENTS OF W WHERE TRANS(U)*W = E. THE VECTORS ARE
  * FREQUENTLY RESCALED TO AVOID OVERFLOW. SOLVE TRANS(U)*W = E
  */
    ek = F_ONE;
    sset (*n, F_ZERO, z, 1);
    for (k = 1; k <= *n; k++) {
	if (z[k - 1] != F_ZERO)
	    ek = sign (ek, -z[k - 1]);
	if (fabs (ek - z[k - 1]) > fabs (*FAC (k - 1, k - 1))) {
	    s = fabs (*FAC (k - 1, k - 1)) / fabs (ek - z[k - 1]);
	    sscal (*n, s, z, 1);
	    ek *= s;
	}
	wk = ek - z[k - 1];
	wkm = -ek - z[k - 1];
	s = fabs (wk);
	sm = fabs (wkm);
	if (*FAC (k - 1, k - 1) != F_ZERO) {
	    wk /= *FAC (k - 1, k - 1);
	    wkm /= *FAC (k - 1, k - 1);
	}
	else {
	    wk = F_ONE;
	    wkm = F_ONE;
	}
	kp1 = k + 1;
	if (kp1 <= *n) {
	    for (j = kp1; j <= *n; j++) {
		sm += fabs (z[j - 1] + wkm ** FAC (j - 1, k - 1));
		z[j - 1] += wk ** FAC (j - 1, k - 1);
		s += fabs (z[j - 1]);
	    }
	    if (s < sm) {
		t = wkm - wk;
		wk = wkm;
		saxpy (*n - k, t, FAC (kp1 - 1, k - 1), *ldfac, &z[kp1 - 1],
		    1);
	    }
	}
	z[k - 1] = wk;
    }
    _l0 = 1;
    s = F_ONE / imsl_sasum (*n, z, _l0);
    sscal (*n, s, z, 1);
 /* SOLVE TRANS(L)*Y = W */
    for (k = *n; k >= 1; k--) {
	if (k < *n)
	    z[k - 1] += imsl_sdot (*n - k, FAC (k - 1, k), 1, &z[k], 1);
	if (fabs (z[k - 1]) > F_ONE) {
	    s = F_ONE / fabs (z[k - 1]);
	    sscal (*n, s, z, 1);
	}
	l = ipvt[k - 1];
	t = z[l - 1];
	z[l - 1] = z[k - 1];
	z[k - 1] = t;
    }
    _l0 = 1;
    s = F_ONE / imsl_sasum (*n, z, _l0);
    sscal (*n, s, z, 1);

    ynorm = F_ONE;
 /* SOLVE L*V = Y */
    for (k = 1; k <= *n; k++) {
	l = ipvt[k - 1];
	t = z[l - 1];
	z[l - 1] = z[k - 1];
	z[k - 1] = t;
	if (k < *n)
	    saxpy (*n - k, t, FAC (k - 1, k), 1, &z[k], 1);
	if (fabs (z[k - 1]) > F_ONE) {
	    s = F_ONE / fabs (z[k - 1]);
	    sscal (*n, s, z, 1);
	    ynorm *= s;
	}
    }
    _l0 = 1;
    s = F_ONE / imsl_sasum (*n, z, _l0);
    sscal (*n, s, z, 1);
    ynorm *= s;
 /* SOLVE U*Z = V */
    for (k = *n; k >= 1; k--) {
	if (fabs (z[k - 1]) > fabs (*FAC (k - 1, k - 1))) {
	    s = fabs (*FAC (k - 1, k - 1)) / fabs (z[k - 1]);
	    sscal (*n, s, z, 1);
	    ynorm *= s;
	}
	if (*FAC (k - 1, k - 1) != F_ZERO) {
	    z[k - 1] /= *FAC (k - 1, k - 1);
	}
	else {
	    z[k - 1] = F_ONE;
	}
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
    imsl_e1pop ("L2CRG ");
    return;
}	       /* end of function */


#undef A
#undef FAC

/* Structured by FOR_STRUCT, v0.2, on 09/14/90 at 11:46:51
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
static void     l_linrt (Mint *n, Mfloat *a, Mint *lda, Mint *ipath,
                    Mfloat *ainv, Mint *ldainv)
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


    imsl_e1psh ("LINRT ");

    if (*n <= 0) {
	imsl_e1sti (1, *n);

    /*
     * (5, 1, "The order of the matrix must be positive while N = %(i1) is
     * given.");
     */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
    }
    else if (*n > *lda) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *lda);

    /*
     * (5, 2, "The order of the matrix must be less than or equal to its
     * leading dimension while N = %(i1) and LDA = %(i2) are given.");
     */
	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_LESS_ORDER);
    }
    else if (*n > *ldainv) {
	imsl_e1sti (1, *n);
	imsl_e1sti (2, *ldainv);

    /*
     * (5, 3, "The order of the matrix must be less than or equal to its
     * leading dimension while N = %(i1) and LDAINV = %(i2) are given.");
     */
	imsl_ermes (IMSL_TERMINAL, IMSL_LDAINV_LESS_ORDER);
    }
    else if (*ipath != 1 && *ipath != 2) {
	imsl_e1sti (1, *ipath);

    /*
     * (5, 4, "IPATH must be either 1 or 2 while a value of %(i1) is
     * given.");
     */
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
/*	imsl_e1sti (1, info); */

    /*
     * (5, 5, "The matrix to be inverted is singular.  The index of the first
     * zero diagonal element of A is %(i1).");
     */
	imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_MATRIX);
    }
L_9000:
    imsl_e1pop ("LINRT ");
    return;
}	       /* end of function */

#undef A
#undef AINV


/* Structured by FOR_STRUCT, v0.2, on 09/14/90 at 11:59:13
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  NR1RR/DNR1RR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 17, 1985

    Purpose:    Compute the 1-norm of a real matrix.

    Usage:      CALL NR1RR (NRA, NCA, A, LDA, ANORM)

    Arguments:
       NRA    - Number of rows of A.  (Input)
       NCA    - Number of columns of A.  (Input)
       A      - Real NRA by NCA matrix whose 1-norm is to be computed.
                (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       ANORM  - Real scalar containing the 1-norm of A.  (Output)

    Keyword:    Basic matrix operation

    GAMS:       D1b2

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void     l_nr1rr (Mint *nra, Mint *nca, Mfloat *a, Mint *lda,
                    Mfloat *anorm)
#else
static void     l_nr1rr (nra, nca, a, lda, anorm)
Mint           *nra, *nca;
Mfloat         *a;
Mint           *lda;
Mfloat         *anorm;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
    Mint            _l0, j;
    Mfloat          anorm1;


    imsl_e1psh ("NR1RR ");
 /* CHECK FOR INPUT ERRORS */
    if (*nra > *lda) {
	imsl_e1sti (1, *nra);
	imsl_e1sti (2, *lda);

    /*
     * (5, 1, "The number of rows of the input matrix must be less than or
     * equal to the leading dimension while NRA = %(i1) and LDA = %(i2) are
     * given.");
     */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NUM_ROWS_LT_LEAD);
	goto L_9000;
    }
    if (*nra <= 0) {
	imsl_e1sti (1, *nra);

    /*
     * (5, 2, "The number of rows of the input matrix must be greater than
     * zero while NRA = %(i1) is given.");
     */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NUM_ROWS_GT_ZERO);
    }
    if (*nca <= 0) {
	imsl_e1sti (1, *nca);

    /*
     * (5, 3, "The number of columns of the input matrix must be greater than
     * zero while NCA = %(i1) is given.");
     */
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NUM_COLS_GT_ZERO);
    }
    if (imsl_n1rcd (0) != 0)
	goto L_9000;
 /* CALCULATE THE L1 NORM FOR A. */
    *anorm = F_ZERO;
    _l0 = 1;
    for (j = 1; j <= *nca; j++) {
	anorm1 = imsl_sasum (*nra, A (j - 1, 0), _l0);
	*anorm = imsl_f_max (anorm1, *anorm);
    }

L_9000:
    imsl_e1pop ("NR1RR ");
    return;
}	       /* end of function */

#undef A
