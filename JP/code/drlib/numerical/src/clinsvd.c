#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif


static VA_LIST_HACK PROTO (l_lin_svd_gen, (Mint nra, Mint nca, Mf_complex *a,
	            va_list argptr));
    static void PROTO (l_l2vcr, (Mint *nra, Mint *nca, Mf_complex *a, Mint *lda,
	            Mint *ipath, Mfloat *tol, Mint *irank, Mf_complex *s,
	            Mf_complex *u, Mint *ldu, Mf_complex *v, Mint *ldv,
	            Mf_complex *wka, Mf_complex *wk));
    static void PROTO (l_l3vcr, (Mint *nra, Mint *nca, Mf_complex *a, Mint *lda,
	            Mint *ipath, Mfloat *tol, Mint *irank, Mf_complex *s,
	            Mf_complex *u, Mint *ldu, Mf_complex *v, Mint *ldv,
	            Mf_complex *e, Mf_complex *work));
    static void PROTO (l_csrot, (Mint *n, Mf_complex *cx, Mint *incx, Mf_complex *cy,
	            Mint *incy, Mfloat *c, Mfloat *s));

    static Mf_complex *lv_s;


#ifdef ANSI
    Mf_complex *imsl_c_lin_svd_gen (Mint nra, Mint nca, Mf_complex *a,...)
#else
    Mf_complex *imsl_c_lin_svd_gen (nra, nca, a, va_alist)
    Mint        nra, nca;
    Mf_complex *a;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, a);

    E1PSH ("imsl_c_lin_svd_gen", "imsl_z_lin_svd_gen");

    lv_s = NULL;
    IMSL_CALL (l_lin_svd_gen (nra, nca, a, argptr));
    va_end (argptr);

    E1POP ("imsl_c_lin_svd_gen", "imsl_z_lin_svd_gen");

    return lv_s;
}
#ifdef ANSI
static VA_LIST_HACK l_lin_svd_gen (Mint nra, Mint nca, Mf_complex *a, va_list argptr)
#else
static VA_LIST_HACK l_lin_svd_gen (nra, nca, a, argptr)
    Mint        nra, nca;
    Mf_complex *a;
    va_list     argptr;
#endif
{
    Mint        a_col_dim = nca;
    Mfloat      tol;
    Mint       *rank = NULL;

    Mint        upath = 0;
    Mf_complex **p_u;
    Mf_complex *u_vector = NULL;
    Mint        u_col_dim = 0;
    Mint        vpath = 0;
    Mf_complex **p_v;
    Mf_complex *v_vector = NULL;
    Mint        v_col_dim = 0;
    Mf_complex **p_ginva;
    Mf_complex *ginva = NULL;
    Mint        ginva_col_dim = nra;

    Mint         code = 1;
    short int   arg_number = 3;

    short int   return_u = 0;
    short int   user_u = 0;
    short int   return_v = 0;
    short int   user_v = 0;
    short int   return_rank = 0;
    short int   user_s = 0;
    short int   return_ginva = 0;
    short int   user_ginva = 0;
    short int   error = 0;
    short int   v_user_memory = 0;
    short int   s_user_memory = 0;
    short int   free_s = 0;
    short int   free_v = 0;

    Mint        ncol_v, length_s, max_dim, min_dim, i;
    Mf_complex *wk = NULL;
    Mf_complex *wka = NULL;
    Mf_complex *u = NULL;
    Mf_complex *v = NULL;
    Mf_complex *s = NULL;
    Mint        one = 1;
    Mf_complex  c_zero;
    tol = 100.0 * imsl_amach (4);
    c_zero = imsl_cf_convert (F_ZERO, F_ZERO);
    lv_s = NULL;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	++arg_number;
	switch (code) {
	case IMSL_RETURN_USER:
	    lv_s = va_arg (argptr, Mf_complex *);
	    user_s = 1;
	    if (!lv_s) {
		imsl_e1stl (1, "S");
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
	case IMSL_RANK:
	    tol = (Mfloat) va_arg (argptr, Mdouble);
	    rank = va_arg (argptr, Mint *);
	    return_rank = 1;
	    arg_number += 2;
	    break;
	case IMSL_RANK_ADR:
	    tol = *( va_arg (argptr, Mfloat *));
	    rank = va_arg (argptr, Mint *);
	    return_rank = 1;
	    arg_number += 2;
	    break;
	case IMSL_U:
	    p_u = va_arg (argptr, Mf_complex **);
	    upath = 2;
	    user_u = 0;
	    return_u = 1;
	    arg_number += 2;
	    break;
	case IMSL_U_USER:
	    u_vector = va_arg (argptr, Mf_complex *);
	    if (!u_vector) {
		imsl_e1stl (1, "U");
		imsl_e1stl (2, "IMSL_U_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_2);
		++error;
	    }
	    upath = 2;
	    user_u = 1;
	    return_u = 1;
	    arg_number += 2;
	    break;
	case IMSL_U_COL_DIM:
	    u_col_dim = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_V:
	    p_v = va_arg (argptr, Mf_complex **);
	    vpath = 1;
	    user_v = 0;
	    return_v = 1;
	    ++arg_number;
	    break;
	case IMSL_V_USER:
	    v_vector = va_arg (argptr, Mf_complex *);
	    if (!v_vector) {
		imsl_e1stl (1, "V");
		imsl_e1stl (2, "IMSL_V_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    vpath = 1;
	    user_v = 1;
	    return_v = 1;
	    ++arg_number;
	    break;
	case IMSL_V_COL_DIM:
	    v_col_dim = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_INVERSE:
	    p_ginva = va_arg (argptr, Mf_complex **);
	    user_ginva = 0;
	    return_ginva = 1;
	    ++arg_number;
	    break;
	case IMSL_INVERSE_USER:
	    ginva = va_arg (argptr, Mf_complex *);
	    if (!ginva) {
		imsl_e1stl (1, "gen_inva");
		imsl_e1stl (2, "IMSL_INVERSE_USER");
		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
		++error;
	    }
	    user_ginva = 1;
	    return_ginva = 1;
	    ++arg_number;
	    break;
	case IMSL_INV_COL_DIM:
	    ginva_col_dim = va_arg (argptr, Mint);
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
    if (!a) {
	imsl_e1stl (1, "A");
	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);
	++error;
    }

    if (error)
	return argptr;

    if (nra < 1) {
	imsl_e1sti (1, nra);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NUM_ROWS_GT_ZERO);
	++error;
    }
    if (nca < 1) {
	imsl_e1sti (1, nca);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NUM_COLS_GT_ZERO);
	++error;
    }
    if (error)
	return argptr;

    if (a_col_dim < nca) {
	imsl_e1sti (1, a_col_dim);
	imsl_e1sti (2, nca);
	imsl_e1stl (1, "A");
	imsl_e1stl (2, "nca");
	imsl_ermes (IMSL_TERMINAL, IMSL_COLUMN_DIM_ERROR);
	++error;
    }
    min_dim = (nra < nca) ? nra : nca;
    if (!u_col_dim)
	u_col_dim = min_dim;
    if (!v_col_dim)
	v_col_dim = min_dim;
    if (return_u) {
	if (u_col_dim < min_dim) {
	    imsl_e1sti (1, u_col_dim);
	    imsl_e1sti (2, min_dim);
	    imsl_e1stl (1, "U");
	    imsl_e1stl (2, "min(nra,nca)");
	    imsl_ermes (IMSL_TERMINAL, IMSL_COLUMN_DIM_ERROR);
	    ++error;
	}
    }
    ncol_v = nca;
    if (return_v) {
	if (v_col_dim < min_dim) {
	    imsl_e1sti (1, v_col_dim);
	    imsl_e1sti (2, min_dim);
	    imsl_e1stl (1, "V");
	    imsl_e1stl (2, "min(nra,nca)");
	    imsl_ermes (IMSL_TERMINAL, IMSL_COLUMN_DIM_ERROR);
	    ++error;
	}
	else {
	    if (v_col_dim > ncol_v)
		ncol_v = v_col_dim;
	}
    }

    if (return_ginva && ginva_col_dim < nra) {
	imsl_e1sti (1, ginva_col_dim);
	imsl_e1sti (2, nra);
	imsl_e1stl (1, "gen_inva");
	imsl_e1stl (2, "nra");
	imsl_ermes (IMSL_TERMINAL, IMSL_COLUMN_DIM_ERROR);
	++error;
    }

    if (error)
	return argptr;

    length_s = (nra + 1 < nca) ? nra + 1 : nca;
    max_dim = (nra > nca) ? nra : nca;
    if (!user_s || (nra + 1 < nca)) {
	s = (Mf_complex *) imsl_malloc (length_s * sizeof (Mf_complex));
	if (!user_s) {
	    lv_s = s;
	    free_s = 0;
	}
	else {
	    free_s = 1;
	}
	s_user_memory = 0;
    }
    else {
	s = lv_s;
	s_user_memory = 1;
	free_s = 0;
    }
    if (!user_u && (return_u || return_ginva)) {
	u = (Mf_complex *) imsl_malloc (u_col_dim * nra * sizeof (Mf_complex));
	if (!u) {
	    imsl_e1sti (1, nra);
	    imsl_e1stl (1, "nra");
	    imsl_e1sti (2, nca);
	    imsl_e1stl (2, "nca");
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
	    ++error;
	}
    }
    else if (user_u) {
	if (u_col_dim > min_dim) {
	    imsl_c_m1ran (nra, u_col_dim, u_vector, u_vector);
	    error = imsl_n1rty (1) > 3 ? 1 : 0;
	}
	u = u_vector;
    }


    if (return_v || return_ginva) {
	if (nca == min_dim && v_col_dim >= nca && user_v) {
	    v = v_vector;
	    if (v_col_dim > nca) {
		imsl_c_m1ran (nca, v_col_dim, v, v);
		error = (imsl_n1rty (1) > 3) ? 1 : 0;
	    }
	    v_user_memory = 1;
	    free_v = 0;
	}
	else {
	    v_user_memory = 0;
	    if (ncol_v > v_col_dim) {
		free_v = 1;
	    }
	    else {
		free_v = 0;
	    }
	    v = (Mf_complex *) imsl_malloc (ncol_v * nca * sizeof (Mf_complex));
	    if (!v) {
		imsl_e1sti (1, nra);
		imsl_e1stl (1, "nra");
		imsl_e1sti (2, nca);
		imsl_e1stl (2, "nca");
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
		++error;
	    }
	}
    }
    if (!error) {
	wka = (Mf_complex *) imsl_malloc (nra * nca * sizeof (Mf_complex));
	wk = (Mf_complex *) imsl_malloc ((nra + nca + max_dim - 1) * sizeof (Mf_complex));
	if (!wk || !wka || !s) {
	    imsl_e1sti (1, nra);
	    imsl_e1stl (1, "nra");
	    imsl_e1sti (2, nca);
	    imsl_e1stl (2, "nca");
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
	    ++error;
	}
	else {
	    Mint        lda = nra;
	    Mint        ldu = nra;
	    Mint        ldv = nca;
	    Mint        ipath;
	    Mint        irank;
	    imsl_c_m1ran (nra, a_col_dim, a, a);
	    if (!(error = (imsl_n1rty (1) > 3) ? 1 : 0)) {
		if (return_ginva) {
		    upath = 2;
		    vpath = 1;
		}
		ipath = 10 * upath + vpath;
		l_l2vcr (&nra, &nca, a, &lda, &ipath, &tol, &irank, s, u, &ldu, v, &ldv, wka, wk);
		if (!(error = (imsl_n1rty (1) > 3) ? 1 : 0)) {
		    if (return_rank)
			*rank = irank;
		    if (return_ginva) {
			Mint        j;
			if (!user_ginva) {
			    imsl_free (wk);
			    wk = NULL;
			    imsl_free (wka);
			    wka = NULL;
			    ginva = (Mf_complex *) imsl_malloc (nca * ginva_col_dim * sizeof (Mf_complex));
			    if (!ginva) {
				imsl_e1sti (1, nca);
				imsl_e1sti (2, ginva_col_dim);
				imsl_e1stl (1, "nca");
				imsl_e1stl (2, "gen_inva_col_dim");
/*                            imsl_ermes(5, 2, "Not enough space to compute the generalized inverse.");
*/
				imsl_ermes (IMSL_TERMINAL, IMSL_NO_SPACE_FOR_INVERSE);
				++error;

			    }
			}
			if (!error) {
			    Mf_complex  c_div;
			    for (j = 0; j < nca; j++) {
				imsl_cset (&nra, &c_zero, (ginva + j * ginva_col_dim), &one);
			    }

#if 0
			    for (j = 0; j < irank; j++) {
				for (i = 0; i < nca; i++) {
				    c_div = imsl_c_div (imsl_c_conjg (v[j * nca + i]), lv_s[j]);
				    imsl_caxpy (&nra, &c_div, (u + j * nra), &one, (ginva + i * ginva_col_dim), &one);
				}
			    }
#endif
                            for (j = 0; j < irank; j++) {
				for (i=0; i<nra; i++)
					(u+j*nra+i)->im *= -1.0;
                                for (i = 0; i < nca; i++) {
                                    c_div = imsl_c_div (v[j*nca+i], lv_s[j]);
                                    imsl_caxpy (&nra, &c_div, (u + j * nra), &one, (ginva + i * ginva_col_dim), &one);
                                }
				for (i=0; i<nra; i++)
					(u+j*nra+i)->im *= -1.0;
                            }
			}
		    }
		}
		imsl_c_m1ran (a_col_dim, nra, a, a);
	    }
	}
    }
    if (wk)
	imsl_free (wk);
    if (wka)
	imsl_free (wka);
    if (!return_v || error) {
	if (v)
	    imsl_free (v);
    }
    else {
	if (user_v) {
	    if (v_user_memory) {
		imsl_c_m1ran (v_col_dim, nca, v, v);
	    }
	    else {
		for (i = 0; i < nca; i++) {
		    imsl_ccopy (&min_dim, (v + i), &nca, (v_vector + i * v_col_dim), &one);
		}
		imsl_free(v);
	    }
	}
	else {
	    if (free_v) {
		Mf_complex *v_new;
		v_new = (Mf_complex *) imsl_malloc (v_col_dim * nca * sizeof (Mf_complex));
		for (i = 0; i < nca; i++) {
		    imsl_ccopy (&min_dim, (v + i), &nca, (v_new + i * v_col_dim), &one);
		}
		imsl_free (v);
		v = v_new;
	    }
	    else {
		imsl_c_m1ran (v_col_dim, nca, v, v);
	    }
	    for (i = min_dim; i < v_col_dim; i++) {
		imsl_cset (&nca, &c_zero, (v + i), &v_col_dim);
	    }
	    *p_v = v;
	}
    }
    v = NULL;

    if (!user_u) {
	if (!return_u || error) {
	    if (u)
		imsl_free (u);
	}
	else {
	    imsl_c_m1ran (u_col_dim, nra, u, u);
	    for (i = min_dim; i < u_col_dim; i++) {
		imsl_cset (&nra, &c_zero, (u + i), &u_col_dim);
	    }
	    *p_u = u;
	}
    }
    else {
	imsl_c_m1ran (u_col_dim, nra, u, u);
    }
    u = NULL;

    if (!user_ginva && return_ginva) {
	if (error) {
	    if (ginva)
		imsl_free (ginva);
	}
	else {
	    for (i = nra; i < ginva_col_dim; i++) {
		imsl_cset (&nca, &c_zero, (ginva + i), &ginva_col_dim);
	    }
	    *p_ginva = ginva;
	}
	ginva = NULL;
    }
    if (error) {
	if (!user_s || !s_user_memory) {
	    if (s)
		imsl_free (s);
	    lv_s = NULL;
	}
    }
    else if (free_s) {
	imsl_ccopy (&min_dim, s, &one, lv_s, &one);
	imsl_free (s);
    }
    s = NULL;

    return argptr;
}


/*
  -----------------------------------------------------------------------
    IMSL Name:  L2VCR/DL2VCR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 1, 1985

    Purpose:    Compute the singular value decomposition of a d_complex
                matrix.

    Usage:      CALL L2VCR (NRA, NCA, A, LDA, IPATH, TOL, IRANK, S,
                            U, LDU, V, LDV, WKA, WK)

    Arguments:  See LSVCR.

    Remarks:    See LSVCR.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_l2vcr (Mint *nra, Mint *nca, Mf_complex *a, Mint *lda,
                Mint *ipath, Mfloat *tol, Mint *irank, Mf_complex *s,
                Mf_complex *u, Mint *ldu, Mf_complex *v, Mint *ldv,
                Mf_complex *wka, Mf_complex *wk)
#else
static void l_l2vcr (nra, nca, a, lda, ipath, tol, irank, s, u,
                ldu, v, ldv, wka, wk)
    Mint       *nra, *nca;
    Mf_complex *a;
    Mint       *lda, *ipath;
    Mfloat     *tol;
    Mint       *irank;
    Mf_complex  s[], u[];
    Mint       *ldu;
    Mf_complex  v[];
    Mint       *ldv;
    Mf_complex  wka[], wk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
    Mint        _l0, _l1, inde, indw, j, jobu, jobv;


    imsl_e1psh ("l_l2vcr");

    if (*nra <= 0 || *nca <= 0) {
	imsl_e1sti (1, *nra);
	imsl_e1sti (2, *nca);

/*		imsl_ermes(5, 1, "Both the number of rows and the number of columns of the input matrix have to be positive while NRA = %(i1) and NCA = %(i2) are given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NRA_AND_NCA_GT_ZERO);
	goto L_9000;
    }
    if (*nra > *lda) {
	imsl_e1sti (1, *nra);
	imsl_e1sti (2, *lda);

/*		imsl_ermes(5, 2, "The number of rows of A must be less than or equal to its leading dimension while NRA = %(i1) and LDA = %(i2) are given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_SMALLER_NRA_VALUE);
	goto L_9000;
    }
    if (*ldu <= 0) {
	imsl_e1sti (1, *ldu);

/*		imsl_ermes(5, 4, "The leading dimension of U, LDU, must be greater than zero while LDU = %(i1) is given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_LDU_GT_ZERO);
	goto L_9000;
    }
    if (*ldv <= 0) {
	imsl_e1sti (1, *ldv);

/*		imsl_ermes(5, 5, "The leading dimension of V, LDV, must be greater than zero while LDV = %(i1) is given.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_LDV_GT_ZERO);
	goto L_9000;
    }
    jobu = mod (*ipath, 100) / 10;
    jobv = mod (*ipath, 10);
    if (((jobu != 0 && jobu != 1) && jobu != 2) || (jobv != 0 && jobv !=
	    1)) {
	imsl_e1sti (1, jobu);
	imsl_e1sti (2, jobv);
/*		imsl_ermes(5, 7, "Error in computation control flag. The IJ decimal expansion of IPATH is I = %(i1) and J = %(i2).  I must be either 0, 1 or 2 and J must be either 0 or 1.");
*/
	imsl_ermes (IMSL_TERMINAL, IMSL_CONTROL_FLAG_ERROR);
	goto L_9000;
    }
    /*
     * MAKE A COPY OF A IN WKA AND WORK WITH WKA ONLY.
     */
    for (j = 1; j <= *nca; j++) {
	_l0 = 1;
	_l1 = 1;
	imsl_ccopy (nra, A (j - 1, 0), &_l0, &wka[(j - 1) ** nra], &_l1);
    }

    inde = 1;
    indw = *nca + 1;
    l_l3vcr (nra, nca, wka, nra, ipath, tol, irank, s, u, ldu, v, ldv,
	&wk[inde - 1], &wk[indw - 1]);

L_9000:
    imsl_e1pop ("l_l2vcr");
    return;
}				/* end of function */
#undef  A
/*-----------------------------------------------------------------------

    IMSL Name:  L3VCR/DL3VCR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 25, 1987

    Purpose:    Compute the singular value decomposition of a d_complex
                matrix.

    Usage:      CALL L3VCR (NRA, NCA, A, LDA, IPATH, TOL, IRANK, S,
                            U, LDU, V, LDV, E, WORK)

    Arguments:  See LSVCR.

    Remarks:    See LSVCR.

    Chapter:    MATH/LIBRARY Linear Systems

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_l3vcr (Mint *nra, Mint *nca, Mf_complex *a, Mint *lda,
                Mint *ipath, Mfloat *tol, Mint *irank, Mf_complex *s,
                Mf_complex *u, Mint *ldu, Mf_complex *v, Mint *ldv,
                Mf_complex *e, Mf_complex *work)
#else
static void l_l3vcr (nra, nca, a, lda, ipath, tol, irank, s, u,
                ldu, v, ldv, e, work)
    Mint       *nra, *nca;
    Mf_complex *a;
    Mint       *lda, *ipath;
    Mfloat     *tol;
    Mint       *irank;
    Mf_complex  s[], *u;
    Mint       *ldu;
    Mf_complex *v;
    Mint       *ldv;
    Mf_complex  e[], work[];
#endif
{
/*#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define U(I_,J_)	(u+(I_)*(*ldu)+(J_))
#define V(I_,J_)	(v+(I_)*(*ldv)+(J_))*/
#define A(I_,J_)        (a+(I_)*(l_lda)+(J_))
#define U(I_,J_)        (u+(I_)*(l_ldu)+(J_))
#define V(I_,J_)        (v+(I_)*(l_ldv)+(J_))
     Mint l_lda = *lda;
     Mint l_ldu = *ldu;
     Mint l_ldv = *ldv;
    short int   wantu, wantv;
    Mint        _l0, _l1, _l2, _l3, i, iend, info, istart, iter, j, jobu,
                k, kase, l, ll, lls, lp1, ls, lu, m, maxit, minmn, mm,
                mp1, nct, nctp1, ncu, nrt, nrtp1;
    Mfloat      anorm, b, c, cs, el, emm1, f, g, scale, shift, sl, sm,
                smm1, sn, stol, t1, test, ztest;
    Mf_complex  _cx0, _cx1, r, t;
    /* *SQRT */

#define L3CCT(zdum)	(Mfloat)(fabs( imsl_fc_convert( (zdum) ) ) + fabs( imsl_c_aimag( (zdum) ) ))

    imsl_e1psh ("l_l3vcr");
    /*
     * IF TOL.LT.0, COMPUTE INFINITY NORM OF A AND SINGULAR VALUE TOLERANCE
     * FOR RANK ESTIMATION.
     */
    if (*tol < F_ZERO) {
	anorm = F_ZERO;
	for (i = 1; i <= *nra; i++) {
	    anorm = imsl_f_max (anorm, imsl_scasum (nca, A (0, i - 1), lda));
	}
	stol = fabs (*tol) * anorm;
    }
    else {
	stol = *tol;
    }
    /* SET THE MAXIMUM NUMBER OF ITERATIONS */
    maxit = 30;
    /* DETERMINE WHAT IS TO BE COMPUTED. */
    wantu = 0;
    wantv = 0;
    jobu = mod (*ipath, 100) / 10;
    ncu = *nra;
    if (jobu == 2)
	ncu = imsl_i_min (*nra, *nca);
    if (jobu != 0)
	wantu = 1;
    if (mod (*ipath, 10) != 0)
	wantv = 1;
    /*
     * REDUCE A TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS IN S AND
     * THE SUPER-DIAGONAL ELEMENTS IN E.
     */
    info = 0;
    nct = imsl_i_min (*nra - 1, *nca);
    nrt = imsl_i_max (0, imsl_i_min (*nca - 2, *nra));
    lu = imsl_i_max (nct, nrt);
    for (l = 1; l <= lu; l++) {
	s[l - 1] = imsl_cf_convert (F_ZERO, F_ZERO);
	lp1 = l + 1;
	if (l <= nct) {
	    /*
	     * COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND PLACE THE
	     * L-TH DIAGONAL IN S(L).
	     */
	    _l0 = *nra - l + 1;
	    _l1 = 1;
	    s[l - 1] = imsl_cf_convert (imsl_scnrm2 (&_l0, A (l - 1, l - 1), &_l1), F_ZERO);
	    if (L3CCT (s[l - 1]) != F_ZERO) {
		if (L3CCT (*A (l - 1, l - 1)) != F_ZERO)
		    s[l - 1] = imsl_c_mul (imsl_cf_convert (imsl_c_abs (s[l - 1]), F_ZERO), (imsl_c_div (*A (l - 1, l - 1),
				imsl_cf_convert (imsl_c_abs (*A (l - 1, l - 1)), F_ZERO))));
		_l0 = *nra - l + 1;
		_l1 = 1;
		_cx0 = imsl_c_div (imsl_cf_convert (F_ONE, F_ZERO), s[l - 1]);
		imsl_cscal (&_l0, &_cx0, A (l - 1, l - 1), &_l1);
		*A (l - 1, l - 1) = imsl_c_add (imsl_cf_convert (F_ONE, F_ZERO), *A (l - 1, l - 1));
	    }
	    s[l - 1] = imsl_c_neg (s[l - 1]);
	}
	if (l <= nct && L3CCT (s[l - 1]) != F_ZERO) {
	    /* APPLY THE TRANSFORMATION. */
	    _l0 = *nra - l + 1;
	    _l1 = *nca - l;
	    _cx0 = imsl_c_neg (imsl_c_div (imsl_cf_convert (F_ONE,
			F_ZERO), imsl_c_conjg (*A (l - 1, l - 1))));
	    _l2 = 1;
	    _cx1 = imsl_cf_convert (F_ZERO, F_ZERO);
	    _l3 = 1;
	    imsl_cgemv ("C", sizeof ("C"), &_l0, &_l1, &_cx0,
		A (lp1 - 1, l - 1), lda, A (l - 1, l - 1), &_l2, &_cx1,
		&work[lp1 - 1], &_l3);
	    _l0 = *nra - l + 1;
	    _l1 = *nca - l;
	    _cx0 = imsl_cf_convert (F_ONE, F_ZERO);
	    _l2 = 1;
	    _l3 = 1;
	    imsl_cgerc (&_l0, &_l1, &_cx0,
		A (l - 1, l - 1), &_l2, &work[lp1 - 1], &_l3,
		A (lp1 - 1, l - 1), lda);
	    /*
	     * PLACE THE L-TH ROW OF X INTO E FOR THE SUBSEQUENT CALCULATION
	     * OF THE ROW TRANSFORMATION.
	     */
	}
	for (j = lp1; j <= *nca; j++) {
	    e[j - 1] = imsl_c_conjg (*A (j - 1, l - 1));
	}
	if (wantu && l <= nct) {
	    /*
	     * PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
	     * MULTIPLICATION.
	     */
	    _l0 = *nra - l + 1;
	    _l1 = 1;
	    _l2 = 1;
	    imsl_ccopy (&_l0, A (l - 1, l - 1), &_l1,
		U (l - 1, l - 1), &_l2);
	}
	if (l <= nrt) {
	    /*
	     * COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE L-TH
	     * SUPER-DIAGONAL IN E(L).
	     */
	    _l0 = *nca - l;
	    _l1 = 1;
	    e[l - 1] = imsl_cf_convert (imsl_scnrm2 (&_l0, &e[lp1 - 1],
		    &_l1), F_ZERO);
	    if (L3CCT (e[l - 1]) != F_ZERO) {
		if (L3CCT (e[lp1 - 1]) != F_ZERO)
		    e[l - 1] = imsl_c_mul (imsl_cf_convert (imsl_c_abs (e[l - 1]), F_ZERO), (imsl_c_div (e[lp1 - 1],
				imsl_cf_convert (imsl_c_abs (e[lp1 - 1]), F_ZERO))));
		_l0 = *nca - l;
		_l1 = 1;
		_cx0 = imsl_c_div (imsl_cf_convert (F_ONE, F_ZERO), e[l - 1]);
		imsl_cscal (&_l0, &_cx0,
		    &e[lp1 - 1], &_l1);
		e[lp1 - 1] = imsl_c_add (imsl_cf_convert (F_ONE, F_ZERO), e[lp1 - 1]);
	    }
	    e[l - 1] = imsl_c_neg (imsl_c_conjg (e[l - 1]));
	    if (lp1 <= *nra && L3CCT (e[l - 1]) != F_ZERO) {
		/* APPLY THE TRANSFORMATION. */
		_l0 = *nra - l;
		_l1 = *nca - l;
		_cx0 = imsl_cf_convert (F_ONE, F_ZERO);
		_l2 = 1;
		_cx1 = imsl_cf_convert (F_ZERO, F_ZERO);
		_l3 = 1;
		imsl_cgemv ("N", sizeof ("N"), &_l0, &_l1, &_cx0, A (lp1 - 1, lp1 - 1),
		    lda, &e[lp1 - 1], &_l2, &_cx1,
		    &work[lp1 - 1], &_l3);
		_l0 = *nra - l;
		_l1 = *nca - l;
		_cx0 =
		    imsl_c_neg (imsl_c_div (imsl_cf_convert (F_ONE,
			    F_ZERO), imsl_c_conjg (e[lp1 - 1])));
		_l2 = 1;
		_l3 = 1;
		imsl_cgerc (&_l0, &_l1, &_cx0, &work[lp1 - 1], &_l2, &e[lp1 - 1],
		    &_l3, A (lp1 - 1, lp1 - 1), lda);
	    }
	    if (wantv) {
		/*
		 * PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT BACK
		 * MULTIPLICATION.
		 */
		_l0 = *nca - l;
		_l1 = 1;
		_l2 = 1;
		imsl_ccopy (&_l0, &e[lp1 - 1], &_l1,
		    V (l - 1, lp1 - 1), &_l2);
	    }
	}
    }
    /*
     * SET UP THE FINAL BIDIAGONAL MATRIX OF ORDER M.
     */
    m = imsl_i_min (*nca, *nra + 1);
    nctp1 = nct + 1;
    nrtp1 = nrt + 1;
    if (nct < *nca)
	s[nctp1 - 1] = *A (nctp1 - 1, nctp1 - 1);
    if (*nra < m)
	s[m - 1] = imsl_cf_convert (F_ZERO, F_ZERO);
    if (nrtp1 < m)
	e[nrtp1 - 1] = *A (m - 1, nrtp1 - 1);
    e[m - 1] = imsl_cf_convert (F_ZERO, F_ZERO);
    /* IF REQUIRED, GENERATE U. */
    if (wantu) {
	for (j = nctp1; j <= ncu; j++) {
	    _cx0 = imsl_cf_convert (F_ZERO, F_ZERO);
	    _l0 = 1;
	    imsl_cset (nra, &_cx0, U (j - 1, 0), &_l0);
	    *U (j - 1, j - 1) = imsl_cf_convert (F_ONE, F_ZERO);
	}
	for (l = nct; l >= 1; l--) {
	    lp1 = l + 1;
	    if (L3CCT (s[l - 1]) != F_ZERO) {
		_l0 = *nra - l + 1;
		_l1 = ncu - l;
		_cx0 = imsl_c_neg (imsl_c_div (imsl_cf_convert (F_ONE, F_ZERO), imsl_c_conjg (*U (l - 1, l - 1))));
		_l2 = 1;
		_cx1 = imsl_cf_convert (F_ZERO, F_ZERO);
		_l2 = 1;
		_l3 = 1;
		imsl_cgemv ("C", sizeof ("C"), &_l0, &_l1, &_cx0,
		    U (lp1 - 1, l - 1), ldu, U (l - 1, l - 1), &_l2,
		    &_cx1, &work[*nra], &_l3);
		_l1 = ncu - l;
		_l2 = 1;
		_l3 = 1;
		_l0 = *nra - l + 1;
		_cx0 = imsl_cf_convert (F_ONE, F_ZERO);
		imsl_cgerc (&_l0, &_l1, &_cx0,
		    U (l - 1, l - 1), &_l2, &work[*nra], &_l3,
		    U (lp1 - 1, l - 1), ldu);
		_l0 = *nra - l + 1;
		_l1 = 1;
		_cx0 = imsl_cf_convert (-F_ONE, F_ZERO);
		imsl_cscal (&_l0, &_cx0,
		    U (l - 1, l - 1), &_l1);
		*U (l - 1, l - 1) = imsl_c_add (imsl_cf_convert (F_ONE, F_ZERO), *U (l - 1, l - 1));
		_l0 = l - 1;
		_l1 = 1;
		_cx0 = imsl_cf_convert (F_ZERO, F_ZERO);
		imsl_cset (&_l0, &_cx0,

		    U (l - 1, 0), &_l1);
	    }
	    else {
		_cx0 = imsl_cf_convert (F_ZERO, F_ZERO);
		_l0 = 1;
		imsl_cset (nra, &_cx0, U (l - 1, 0),
		    &_l0);
		*U (l - 1, l - 1) = imsl_cf_convert (F_ONE, F_ZERO);
	    }
	}
    }
    /* IF IT IS REQUIRED, GENERATE V. */
    if (wantv) {
	for (l = *nca; l >= 1; l--) {
	    lp1 = l + 1;
	    if (l <= nrt && L3CCT (e[l - 1]) != F_ZERO) {
		_l0 = *nca - l;
		_l1 = *nca - l;
		_cx0 =
		    imsl_c_neg (imsl_c_div (imsl_cf_convert (F_ONE,
			    F_ZERO), imsl_c_conjg (*V (l - 1, lp1 - 1))));
		_l2 = 1;
		_l3 = 1;
		_cx1 = imsl_cf_convert (F_ZERO, F_ZERO);
		imsl_cgemv ("C", sizeof ("C"), &_l0, &_l1, &_cx0,
		    V (lp1 - 1, lp1 - 1), ldv, V (l - 1, lp1 - 1), &_l2,
		    &_cx1, &work[*nra], &_l3);
		_l0 = *nca - l;
		_l1 = *nca - l;
		_cx0 = imsl_cf_convert (F_ONE, F_ZERO);
		_l2 = 1;
		_l3 = 1;
		imsl_cgerc (&_l0, &_l1, &_cx0,
		    V (l - 1, lp1 - 1), &_l2, &work[*nra], &_l3,
		    V (lp1 - 1, lp1 - 1), ldv);
	    }
	    _cx0 = imsl_cf_convert (F_ZERO, F_ZERO);
	    _l0 = 1;
	    imsl_cset (nca, &_cx0, V (l - 1, 0), &_l0);
	    *V (l - 1, l - 1) = imsl_cf_convert (F_ONE, F_ZERO);
	}
    }
    /*
     * TRANSFORM S AND E SO THAT THEY ARE REAL
     */
    for (i = 1; i <= m; i++) {
	if (L3CCT (s[i - 1]) != F_ZERO) {
	    t = imsl_cf_convert (imsl_c_abs (s[i - 1]), F_ZERO);
	    r = imsl_c_div (s[i - 1], t);
	    s[i - 1] = t;
	    if (i < m)
		e[i - 1] = imsl_c_div (e[i - 1], r);
	    _l0 = 1;
	    if (wantu)
		imsl_cscal (nra, &r, U (i - 1, 0), &_l0);
	}
	if (i != m) {
	    if (L3CCT (e[i - 1]) != F_ZERO) {
		t = imsl_cf_convert (imsl_c_abs (e[i - 1]), F_ZERO);
		r = imsl_c_div (t, e[i - 1]);
		e[i - 1] = t;
		s[i] = imsl_c_mul (s[i], r);
		_l0 = 1;
		if (wantv)
		    imsl_cscal (nca, &r, V (i, 0), &_l0);
	    }
	}
    }
    /*
     * MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
     */
    mm = m;
    iter = 0;
L_80:
    ;
    /*
     * QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND. ...EXIT
     */
    if (m != 0 && iter < maxit) {
	/*
	 * THIS SECTION OF THE PROGRAM INSPECTS FOR NEGLIGIBLE ELEMENTS IN
	 * THE S AND E ARRAYS. ON COMPLETION THE VARIABLES KASE AND L ARE SET
	 * AS FOLLOWS. KASE = 1 IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
	 * KASE = 2 IF S(L) IS NEGLIGIBLE AND L.LT.M KASE = 3 IF E(L-1) IS
	 * NEGLIGIBLE, L.LT.M, AND S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR
	 * STEP). KASE = 4 IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
	 */
	for (ll = 1; ll <= m; ll++) {
	    l = m - ll;
	    if (l == 0)
		goto L_100;
	    test = imsl_c_abs (s[l - 1]) + imsl_c_abs (s[l]);
	    ztest = test + imsl_c_abs (e[l - 1]);
	    if (ztest == test) {
		e[l - 1] = imsl_cf_convert (F_ZERO, F_ZERO);
		goto L_100;
	    }
	}

L_100:
	if (l == (m - 1)) {
	    kase = 4;
	}
	else {
	    lp1 = l + 1;
	    mp1 = m + 1;
	    for (lls = lp1; lls <= mp1; lls++) {
		ls = m - lls + lp1;
		if (ls == l)
		    goto L_120;
		test = F_ZERO;
		if (ls != m)
		    test += imsl_c_abs (e[ls - 1]);
		if (ls != (l + 1))
		    test += imsl_c_abs (e[ls - 2]);
		ztest = test + imsl_c_abs (s[ls - 1]);
		if (ztest == test) {
		    s[ls - 1] = imsl_cf_convert (F_ZERO, F_ZERO);
		    goto L_120;
		}
	    }

    L_120:
	    if (ls == l) {
		kase = 3;
	    }
	    else if (ls == m) {
		kase = 1;
	    }
	    else {
		kase = 2;
		l = ls;
	    }
	}
	l += 1;
	/* PERFORM THE TASK INDICATED BY KASE. */
	if (kase == 1) {
	    /* DEFLATE NEGLIGIBLE S(M). */
	    f = imsl_fc_convert (e[m - 2]);
	    e[m - 2] = imsl_cf_convert (F_ZERO, F_ZERO);
	    for (k = m - 1; k >= l; k--) {
		t1 = imsl_fc_convert (s[k - 1]);
		imsl_srotg (&t1, &f, &cs, &sn);
		s[k - 1] = imsl_cf_convert (t1, F_ZERO);
		if (k != l) {
		    f = -sn * imsl_fc_convert (e[k - 2]);
		    e[k - 2] = imsl_c_mul (imsl_cf_convert (cs, F_ZERO), e[k - 2]);
		}
		_l0 = 1;
		_l1 = 1;
		if (wantv)
		    l_csrot (nca, V (k - 1, 0), &_l0, V (m - 1, 0),
			&_l1, &cs, &sn);
	    }

	}
	else if (kase == 2) {
	    /* SPLIT AT NEGLIGIBLE S(L). */
	    f = imsl_fc_convert (e[l - 2]);
	    e[l - 2] = imsl_cf_convert (F_ZERO, F_ZERO);
	    for (k = l; k <= m; k++) {
		t1 = imsl_fc_convert (s[k - 1]);
		imsl_srotg (&t1, &f, &cs, &sn);
		s[k - 1] = imsl_cf_convert (t1, F_ZERO);
		f = -sn * imsl_fc_convert (e[k - 1]);
		e[k - 1] = imsl_c_mul (imsl_cf_convert (cs, F_ZERO), e[k - 1]);
		_l0 = 1;
		_l1 = 1;
		if (wantu)
		    l_csrot (nra, U (k - 1, 0), &_l0, U (l - 2, 0),
			&_l1, &cs, &sn);
	    }

	}
	else if (kase == 3) {
	    /*
	     * PERFORM ONE QR STEP. CALCULATE THE SHIFT.
	     */
	    scale = imsl_f_vmax (5, imsl_c_abs (s[m - 1]), imsl_c_abs (s[m - 2]), imsl_c_abs (e[m - 2]),
		imsl_c_abs (s[l - 1]), imsl_c_abs (e[l - 1]));
	    sm = imsl_fc_convert (s[m - 1]) / scale;
	    smm1 = imsl_fc_convert (s[m - 2]) / scale;
	    emm1 = imsl_fc_convert (e[m - 2]) / scale;
	    sl = imsl_fc_convert (s[l - 1]) / scale;
	    el = imsl_fc_convert (e[l - 1]) / scale;
	    b = ((smm1 + sm) * (smm1 - sm) + imsl_fi_power (emm1, 2)) / F_TWO;
	    c = imsl_fi_power (sm * emm1, 2);
	    shift = F_ZERO;
	    if (b != F_ZERO || c != F_ZERO) {
		shift = sqrt (b * b + c);
		if (b < F_ZERO)
		    shift = -shift;
		shift = c / (b + shift);
	    }
	    f = (sl + sm) * (sl - sm) - shift;
	    g = sl * el;
	    /* CHASE ZEROS. */
	    for (k = l; k <= (m - 1); k++) {
		imsl_srotg (&f, &g, &cs, &sn);
		if (k != l)
		    e[k - 2] = imsl_cf_convert (f, F_ZERO);
		f = cs * imsl_fc_convert (s[k - 1]) + sn * imsl_fc_convert (e[k - 1]);
		e[k - 1] = imsl_c_sub (imsl_c_mul (imsl_cf_convert (cs, F_ZERO), e[k - 1]), imsl_c_mul (imsl_cf_convert (sn, F_ZERO),
			s[k - 1]));
		g = sn * imsl_fc_convert (s[k]);
		s[k] = imsl_c_mul (imsl_cf_convert (cs, F_ZERO), s[k]);
		_l0 = 1;
		_l1 = 1;
		if (wantv)
		    l_csrot (nca, V (k - 1, 0), &_l0, V (k, 0), &_l1,
			&cs, &sn);
		imsl_srotg (&f, &g, &cs, &sn);
		s[k - 1] = imsl_cf_convert (f, F_ZERO);
		f = cs * imsl_fc_convert (e[k - 1]) + sn * imsl_fc_convert (s[k]);
		s[k] = imsl_c_add (imsl_c_neg (imsl_c_mul (imsl_cf_convert (sn, F_ZERO), e[k - 1])), imsl_c_mul (imsl_cf_convert (cs, F_ZERO),
			s[k]));
		g = sn * imsl_fc_convert (e[k]);
		e[k] = imsl_c_mul (imsl_cf_convert (cs, F_ZERO), e[k]);
		_l0 = 1;
		_l1 = 1;
		if (wantu && k < *nra)
		    l_csrot (nra, U (k - 1, 0), &_l0, U (k, 0), &_l1,
			&cs, &sn);
	    }
	    e[m - 2] = imsl_cf_convert (f, F_ZERO);
	    iter += 1;

	}
	else if (kase == 4) {
	    /*
	     * CONVERGENCE. MAKE THE SINGULAR VALUE POSITIVE.
	     */
	    if (imsl_fc_convert (s[l - 1]) < F_ZERO) {
		s[l - 1] = imsl_c_neg (s[l - 1]);
		_l0 = 1;
		_cx0 = imsl_cf_convert (-F_ONE, F_ZERO);
		if (wantv)
		    imsl_cscal (nca, &_cx0, V (l - 1, 0),
			&_l0);
	    }
	    /* ORDER THE SINGULAR VALUE. */
    L_160:
	    if (l == mm)
		goto L_170;
	    if (imsl_fc_convert (s[l - 1]) >= imsl_fc_convert (s[l]))
		goto L_170;
	    t = s[l - 1];
	    s[l - 1] = s[l];
	    s[l] = t;
	    _l0 = 1;
	    _l1 = 1;
	    if (wantv && l < *nca)
		imsl_cswap (nca, V (l - 1, 0), &_l0, V (l, 0), &_l1);
	    _l0 = 1;
	    _l1 = 1;
	    if (wantu && l < *nra)
		imsl_cswap (nra, U (l - 1, 0), &_l0, U (l, 0), &_l1);
	    l += 1;
	    goto L_160;
    L_170:
	    ;
	    iter = 0;
	    m -= 1;
	}
	goto L_80;
    }
    if (iter > maxit)
	info = m;
    if (info == 0) {
	/* ESTIMATE RANK OF A */
	ls = imsl_i_min (*nra + 1, *nca);
	minmn = imsl_i_min (*nra, *nca);
	*irank = minmn;
	for (i = 1; i <= minmn; i++) {
	    if (L3CCT (s[i - 1]) <= stol) {
		*irank = i - 1;
		goto L_190;
	    }
	}
L_190:
	;
    }
    else {
	istart = info + 1;
	iend = imsl_i_min (*nca, *nra);
	imsl_e1sti (1, istart);
	imsl_e1sti (2, iend);

/*		imsl_ermes(4, 1, "Convergence can only be obtained for the %(i1),..., %(i2) singular values and their corresponding singular vectors.");
*/
	imsl_ermes (IMSL_FATAL, IMSL_SLOWCONVERGENT_MATRIX);
    }

    imsl_e1pop ("l_l3vcr");
    return;
}				/* end of function */
#undef	L3CCT
#undef	A
#undef  U
#undef  V
/*
  -----------------------------------------------------------------------
    IMSL Name:  CSROT (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Apply a d_complex Givens plane rotation.

    Usage:      CALL CSROT (N, CX, INCX, CY, INCY, C, S)

    Arguments:
       N      - Length of vectors X and Y.  (Input)
       CX     - Complex vector of length MAX(N*IABS(INCX),1).
                (Input/Output)
                CSROT replaces X(I) with C*X(I) + S*Y(I) for I=1,...,N.
                X(I) and Y(I) refer to specific elements of CX and CY,
                respectively.  See INCX and INCY argument descriptions.
       INCX   - Displacement between elements of CX.  (Input)
                X(I) is defined to be... CX(1+(I-1)*INCX) if INCX .GE. 0
                or CX(1+(I-N)*INCX) if INCX .LT. 0.
       CY     - Complex vector of length MAX(N*IABS(INCY),1).
                (Input/Output)
                CSROT replaces Y(I) with -S*X(I) + C*Y(I) for I=1,...,N.
                X(I) and Y(I) refer to specific elements of CX and CY,
                respectively.  See INCX and INCY argument descriptions.
       INCY   - Displacement between elements of CY.  (Input)
                Y(I) is defined to be... CY(1+(I-1)*INCY) if INCY .GE. 0
                or CY(1+(I-N)*INCY) if INCY .LT. 0.
       C      - Real scalar containing elements of the rotation matrix.
                (Input)
       S      - Real scalar containing elements of the rotation matrix.
                (Input)

    Remark:
                       ( C  S )        (X(1) ... X(N))
       CSROT applies   (      )   to   (             )
                       (-S  C )        (Y(1) ... Y(N))

    Keywords:   Level 1 BLAS; CSROT

    GAMS:       D1b10

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_csrot (Mint *n, Mf_complex *cx, Mint *incx, Mf_complex *cy,
                Mint *incy, Mfloat *c, Mfloat *s)
#else
static void l_csrot (n, cx, incx, cy, incy, c, s)
    Mint       *n;
    Mf_complex  cx[];
    Mint       *incx;
    Mf_complex  cy[];
    Mint       *incy;
    Mfloat     *c, *s;
#endif
{
    Mint        i, ix, iy;
    Mf_complex  ctemp;


    if (*n > 0) {
	if (*incx != 1 || *incy != 1) {
	    /*
	     * CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1
	     */
	    ix = 1;
	    iy = 1;
	    if (*incx < 0)
		ix = (-*n + 1) ** incx + 1;
	    if (*incy < 0)
		iy = (-*n + 1) ** incy + 1;
	    for (i = 1; i <= *n; i++) {
		ctemp = imsl_c_add (imsl_c_mul (imsl_cf_convert (*c, F_ZERO), cx[ix - 1]), imsl_c_mul (imsl_cf_convert (*s, F_ZERO),
			cy[iy - 1]));
		cy[iy - 1] = imsl_c_sub (imsl_c_mul (imsl_cf_convert (*c, F_ZERO), cy[iy - 1]),
		    imsl_c_mul (imsl_cf_convert (*s, F_ZERO), cx[ix - 1]));
		cx[ix - 1] = ctemp;
		ix += *incx;
		iy += *incy;
	    }
	}
	else {
	    /* CODE FOR BOTH INCREMENTS EQUAL TO 1 */
	    for (i = 1; i <= *n; i++) {
#if 0
		ctemp = imsl_c_add (imsl_c_mul (imsl_cf_convert (*c, F_ZERO), cx[i - 1]), imsl_c_mul (imsl_cf_convert (*s, F_ZERO),
			cy[i - 1]));
		cy[i - 1] = imsl_c_sub (imsl_c_mul (imsl_cf_convert (*c, F_ZERO), cy[i - 1]), imsl_c_mul (imsl_cf_convert (*s, F_ZERO),
			cx[i - 1]));
#endif
		ctemp.re = *c*cx[i-1].re + *s*cy[i-1].re;
		ctemp.im = *c*cx[i-1].im + *s*cy[i-1].im;
		cy[i-1].re = *c*cy[i-1].re - *s*cx[i-1].re;
		cy[i-1].im = *c*cy[i-1].im - *s*cx[i-1].im;
		cx[i - 1] = ctemp;
	    }
	}
    }
    return;
}				/* end of function */
