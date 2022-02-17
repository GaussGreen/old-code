#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif



#ifdef TRUE

#undef TRUE

#define TRUE 1

#else

#define TRUE 1

#endif

#ifdef FALSE

#undef FALSE

#define FALSE 0

#else

#define FALSE 0

#endif



static VA_LIST_HACK PROTO (l_lin_svd_gen, (Mint nra, Mint nca, Mfloat *a, va_list

	            argptr));

    static void PROTO (l_l2vrr, (Mint *nra, Mint *nca, Mfloat *a, Mint *lda, Mint *ipath,

	            Mfloat *tol, Mint *irank, Mfloat *s, Mfloat *u, Mint *ldu,

	            Mfloat *v, Mint *ldv, Mfloat *wka, Mfloat *wk));

    static void PROTO (l_l3vrr, (Mint *nra, Mint *nca, Mfloat *a, Mint *lda, Mint *ipath,

	            Mfloat *tol, Mint *irank, Mfloat *s, Mfloat *u, Mint *ldu,

	            Mfloat *v, Mint *ldv, Mfloat *e, Mfloat *work));



    static Mfloat *lv_s;

#ifdef ANSI

    Mfloat     *imsl_f_lin_svd_gen (Mint nra, Mint nca, Mfloat *a,...)

#else

    Mfloat     *imsl_f_lin_svd_gen (nra, nca, a, va_alist)

    Mint        nra, nca;

    Mfloat     *a;

va_dcl

#endif

{

    va_list     argptr;

    VA_START (argptr, a);



    E1PSH ("imsl_f_lin_svd_gen", "imsl_d_lin_svd_gen");



    lv_s = NULL;

    IMSL_CALL (l_lin_svd_gen (nra, nca, a, argptr));

    va_end (argptr);



    E1POP ("imsl_f_lin_svd_gen", "imsl_d_lin_svd_gen");



    return lv_s;

}

#ifdef ANSI

static VA_LIST_HACK l_lin_svd_gen (Mint nra, Mint nca, Mfloat *a, va_list argptr)

#else

static VA_LIST_HACK l_lin_svd_gen (nra, nca, a, argptr)

    Mint        nra, nca;

    Mfloat     *a;

    va_list     argptr;

#endif

{

    Mint        a_col_dim = nca;

    Mfloat      tol;

    Mint       *rank = NULL;



    Mint        upath = 0;

    Mfloat    **p_u;

    Mfloat     *u_vector = NULL;

    Mint        u_col_dim = 0;

    Mint        vpath = 0;

    Mfloat    **p_v;

    Mfloat     *v_vector = NULL;

    Mint        v_col_dim = 0;

    Mfloat    **p_ginva;

    Mfloat     *ginva = NULL;

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



    Mint        min_dim, length_s, max_dim, i, ncol_v;

    Mfloat     *wk = NULL;

    Mfloat     *wka = NULL;

    Mfloat     *u = NULL;

    Mfloat     *v = NULL;

    Mfloat     *s = NULL;

    tol = 100.0 * imsl_amach (4);

    lv_s = NULL;

    while (code > 0) {

	code = va_arg (argptr, Mint);

	arg_number++;

	switch (code) {

	case IMSL_RETURN_USER:

	    lv_s = va_arg (argptr, Mfloat *);

	    if (!lv_s) {

		imsl_e1stl (1, "S");

		imsl_e1stl (2, "IMSL_RETURN_USER");

		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);

		++error;

	    }

	    user_s = 1;

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

	    p_u = va_arg (argptr, Mfloat **);

	    user_u = 0;

	    return_u = 1;

	    upath = 2;

	    ++arg_number;

	    break;

	case IMSL_U_USER:

	    u_vector = va_arg (argptr, Mfloat *);

	    if (!u_vector) {

		imsl_e1stl (1, "U");

		imsl_e1stl (2, "IMSL_U_USER");

		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);

		++error;

	    }

	    user_u = 1;

	    return_u = 1;

	    upath = 2;

	    ++arg_number;

	    break;

	case IMSL_U_COL_DIM:

	    u_col_dim = va_arg (argptr, Mint);

	    arg_number++;

	    break;

	case IMSL_V:

	    p_v = va_arg (argptr, Mfloat **);

	    user_v = 0;

	    return_v = 1;

	    vpath = 1;

	    ++arg_number;

	    break;

	case IMSL_V_USER:

	    v_vector = va_arg (argptr, Mfloat *);

	    if (!v_vector) {

		imsl_e1stl (1, "V");

		imsl_e1stl (2, "IMSL_V_USER");

		imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);

		++error;

	    }

	    user_v = 1;

	    return_v = 1;

	    vpath = 1;

	    ++arg_number;

	    break;

	case IMSL_V_COL_DIM:

	    v_col_dim = va_arg (argptr, Mint);

	    arg_number++;

	    break;

	case IMSL_INVERSE:

	    p_ginva = va_arg (argptr, Mfloat **);

	    user_ginva = 0;

	    return_ginva = 1;

	    arg_number += 2;

	    break;

	case IMSL_INVERSE_USER:

	    ginva = va_arg (argptr, Mfloat *);

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

	    if (v_col_dim >= ncol_v)

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

	s = (Mfloat *) imsl_malloc (length_s * sizeof (Mfloat));

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

	u = (Mfloat *) imsl_malloc (u_col_dim * nra * sizeof (Mfloat));

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

	    imsl_f_m1ran (nra, u_col_dim, u_vector, u_vector);

	    error = imsl_n1rty (1) > 3 ? 1 : 0;

	}

	if (!error)

	    u = u_vector;

    }

    if (return_v || return_ginva) {

	if (nca == min_dim && v_col_dim >= nca && user_v) {

	    v = v_vector;

	    if (v_col_dim > nca) {

		imsl_f_m1ran (nca, v_col_dim, v, v);

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

	    v = (Mfloat *) imsl_malloc (ncol_v * nca * sizeof (Mfloat));

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

	wka = (Mfloat *) imsl_malloc (nra * nca * sizeof (Mfloat));

	wk = (Mfloat *) imsl_malloc ((nra + nca + max_dim - 1) * sizeof (Mfloat));

	if (!wk || !wka || !s) {

	    imsl_e1sti (1, nra);

	    imsl_e1stl (1, "nra");

	    imsl_e1sti (2, nca);

	    imsl_e1stl (2, "nca");

	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);

	    ++error;

	}

	else if (!error) {

	    Mint        lda = nra;

	    Mint        ldu = nra;

	    Mint        ldv = nca;

	    Mint        ipath;

	    Mint        irank;

	    if (return_ginva) {

		if (upath == 0)

		    upath = 2;

		if (vpath == 0)

		    vpath = 1;

	    }

	    ipath = 10 * upath + vpath;

	    imsl_f_m1ran (nra, a_col_dim, a, a);

	    if (!(error = (imsl_n1rty (1) > 3) ? 1 : 0)) {

		l_l2vrr (&nra, &nca, a, &lda, &ipath, &tol, &irank, s, u, &ldu, v, &ldv, wka, wk);

		if (!(error = (imsl_n1rty (1) > 3) ? 1 : 0)) {

		    if (return_rank)

			*rank = irank;

		    if (return_ginva) {

			Mint        j;

			imsl_free (wk);

			wk = NULL;

			imsl_free (wka);

			wka = NULL;



			if (!user_ginva) {

			    ginva = (Mfloat *) imsl_malloc (nca * ginva_col_dim * sizeof (Mfloat));

			    if (!ginva) {

				imsl_e1sti (1, nca);

				imsl_e1sti (2, ginva_col_dim);

				imsl_e1stl (1, "nca");

				imsl_e1stl (2, "gen_inva_col_dim");

/*                              (5, 2, "Not enough space to compute the generalized inverse.");

*/

				imsl_ermes (IMSL_TERMINAL, IMSL_NO_SPACE_FOR_INVERSE);

				++error;

			    }

			}

			if (!error) {

			    Mfloat      div;

			    for (j = 0; j < nca; j++) {

				sset (nra, F_ZERO, (ginva + j * ginva_col_dim), 1);

			    }

			    for (j = 0; j < irank; j++) {

				for (i = 0; i < nca; i++) {

				    div = v[j * nca + i] / s[j];

				    saxpy (nra, div, (u + j * nra), 1, (ginva + i * ginva_col_dim), 1);

				}

			    }

			}

		    }

		}

		imsl_f_m1ran (a_col_dim, nra, a, a);

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

	v = NULL;

    }

    else {

	if (user_v) {

	    if (v_user_memory) {

		imsl_f_m1ran (v_col_dim, nca, v, v);

	    }

	    else {

		for (i = 0; i < nca; i++) {

		    scopy (min_dim, (v + i), nca, (v_vector + i * v_col_dim), 1);

		}

		imsl_free(v);

	    }

	}

	else {

	    if (free_v) {

		Mfloat     *v_new;

		v_new = (Mfloat *) imsl_malloc (v_col_dim * nca * sizeof (Mfloat));

		for (i = 0; i < nca; i++) {

		    scopy (min_dim, (v + i), nca, (v_new + i * v_col_dim), 1);

		}

		imsl_free (v);

		v = v_new;

	    }

	    else {

		imsl_f_m1ran (v_col_dim, nca, v, v);

	    }

	    for (i = min_dim; i < v_col_dim; i++) {

		sset (nca, F_ZERO, (v + i), v_col_dim);

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

	    imsl_f_m1ran (u_col_dim, nra, u, u);

	    for (i = min_dim; i < u_col_dim; i++) {

		sset (nra, F_ZERO, (u + i), u_col_dim);

	    }

	    *p_u = u;

	}

    }

    else {

	imsl_f_m1ran (u_col_dim, nra, u, u);

    }

    u = NULL;



    if (!user_ginva && return_ginva) {

	if (error) {

	    if (ginva)

		imsl_free (ginva);

	}

	else {

	    for (i = nra; i < ginva_col_dim; i++) {

		sset (nca, F_ZERO, (ginva + i), ginva_col_dim);

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

	scopy (min_dim, s, 1, lv_s, 1);

	imsl_free (s);

    }

    s = NULL;



    return argptr;

}

/*

  -----------------------------------------------------------------------

    IMSL Name:  L2VRR/DL2VRR (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 1, 1985



    Purpose:    Compute the singular value decomposition of a real

                matrix.



    Usage:      CALL L2VRR (NRA, NCA, A, LDA, IPATH, TOL, IRANK,

                            S, U, LDU, V, LDV, WKA, WK)



    Arguments:  See LSVRR/DLSVRR.



    Remarks:    See LSVRR/DLSVRR.



    Chapter:    MATH/LIBRARY Linear Systems



    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_l2vrr (Mint *nra, Mint *nca, Mfloat *a, Mint *lda, Mint *ipath,

                Mfloat *tol, Mint *irank, Mfloat *s, Mfloat *u, Mint *ldu,

                Mfloat *v, Mint *ldv, Mfloat *wka, Mfloat *wk)

#else

static void l_l2vrr (nra, nca, a, lda, ipath, tol, irank, s, u,

                ldu, v, ldv, wka, wk)

    Mint       *nra, *nca;

    Mfloat     *a;

    Mint       *lda, *ipath;

    Mfloat     *tol;

    Mint       *irank;

    Mfloat      s[], u[];

    Mint       *ldu;

    Mfloat      v[];

    Mint       *ldv;

    Mfloat      wka[], wk[];

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

    Mint        i, inde, indw, jobu, jobv;





    imsl_e1psh ("l_l2vrr");



    if (*nra <= 0 || *nca <= 0) {

	imsl_e1sti (1, *nra);

	imsl_e1sti (2, *nca);



/*		(5, 1, "Both the number of rows and the number of columns of the input matrix have to be positive while NRA = %(i1) and NCA = %(i2) are given.");

*/

	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NRA_AND_NCA_GT_ZERO);

    }

    else {

	if (*nra > *lda) {

	    imsl_e1sti (1, *nra);

	    imsl_e1sti (2, *lda);



/*			(5, 2, "The number of rows of A must be less than or equal to its leading dimension while NRA = %(i1) and LDA = %(i2) are given.");

*/

	    imsl_ermes (IMSL_TERMINAL,

		IMSL_NEED_SMALLER_NRA_VALUE);

	}

	if (*ldu <= 0) {

	    imsl_e1sti (1, *ldu);



/*			(5, 3, "The leading dimension of U must be greater than zero.  LDU = %(i1) is given.");

*/

	    imsl_ermes (IMSL_TERMINAL, IMSL_NEED_LDU_GT_ZERO);

	}

	if (*ldv <= 0) {

	    imsl_e1sti (1, *ldv);



/*	    (5, 4, "The leading dimension of V must be greater than zero.  LDV = %(i1) is given.");*/

	    imsl_ermes (IMSL_TERMINAL, IMSL_NEED_LDV_GT_ZERO);

	}

	if (imsl_n1rcd (0) == 0) {

	    jobu = mod (*ipath, 100) / 10;

	    jobv = mod (*ipath, 10);

	    if (((jobu != 0 && jobu != 1) && jobu != 2) || (jobv !=

		    0 && jobv != 1)) {

		imsl_e1sti (1, jobu);

		imsl_e1sti (2, jobv);



/*				(5, 7, "Error in computation control flag. The IJ decimal expansion of IPATH is I = %(i1) and J = %(i2).  I must be either 0, 1 or 2 and J must be either 0 or 1.");

*/

		imsl_ermes (IMSL_TERMINAL,

		    IMSL_CONTROL_FLAG_ERROR);

	    }

	    else {

		/*

		 * MAKE A COPY OF A IN WKA AND WORK WITH WKA ONLY.

		 */

		for (i = 1; i <= *nca; i++) {

		    scopy (*nra, A (i - 1, 0), 1, &wka[(i - 1) ** nra],

			1);

		}



		inde = 1;

		indw = *nca + 1;

		l_l3vrr (nra, nca, wka, nra, ipath, tol, irank, s, u,

		    ldu, v, ldv, &wk[inde - 1], &wk[indw - 1]);

	    }

	}

    }



    imsl_e1pop ("l_l2vrr");

    return;

}				/* end of function */

#undef  A



#if 0 /* old l3vrr */

/*

  -----------------------------------------------------------------------

    IMSL Name:  L3VRR/DL3VRR (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    January 27, 1988



    Purpose:    Compute the singular value decomposition of a real

                matrix.



    Usage:      CALL L3VRR (NRA, NCA, A, LDA, IPATH, TOL, IRANK,

                            S, U, LDU, V, LDV, E, WORK)



    Arguments:  See LSVRR/DLSVRR



    Remarks:    See LSVRR/DLSVRR.



    Chapter:    MATH/LIBRARY Linear Systems



    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_l3vrr (Mint *nra, Mint *nca, Mfloat *a, Mint *lda, Mint *ipath,

                Mfloat *tol, Mint *irank, Mfloat *s, Mfloat *u, Mint *ldu,

                Mfloat *v, Mint *ldv, Mfloat *e, Mfloat *work)

#else

static void l_l3vrr (nra, nca, a, lda, ipath, tol, irank, s, u,

                ldu, v, ldv, e, work)

    Mint       *nra, *nca;

    Mfloat     *a;

    Mint       *lda, *ipath;

    Mfloat     *tol;

    Mint       *irank;

    Mfloat      s[], *u;

    Mint       *ldu;

    Mfloat     *v;

    Mint       *ldv;

    Mfloat      e[], work[];

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

#define U(I_,J_)	(u+(I_)*(*ldu)+(J_))

#define V(I_,J_)	(v+(I_)*(*ldv)+(J_))

#ifdef COMPUTER_DECOSF

    int         wantu, wantv;

    int        IMSLFALSE = 0, IMSLTRUE = 1;

#else

    long        wantu, wantv;

    long        IMSLFALSE = 0, IMSLTRUE = 1;

#endif

    Mint        _l0, _l1, _l2, _l3, i, iend, info, istart, iter, j, jobu,

                k, kase, kk, l, ll, lls, lm1, lp1, ls, lu, m, maxit, minmn,

                mm, mm1, mp1, nct, nctp1, ncu, nrt, nrtp1;

    Mfloat      _f0, _f1, anorm, b, c, cs, el, emm1, eps, f, g, scale, shift,

                sl, sm, smm1, sn, stol, t, t1, test, ztest;





    imsl_e1psh ("l_l3vrr");

    eps = imsl_amach (4);

    /*

     * IF TOL.LT.0, COMPUTE THE INFINITY NORM OF A AND SINGULAR VALUE

     * TOLERANCE FOR RANK ESTIMATION

     */

    if (*tol < F_ZERO) {

	anorm = F_ZERO;

	for (i = 1; i <= *nra; i++) {

	    anorm = imsl_f_max (anorm, imsl_sasum (*nca, A (0, i - 1), *lda));

	}

	stol = fabs (*tol) * anorm;

    }

    else {

	stol = *tol;

    }

    /*

     * SET THE MAXIMUM NUMBER OF ITERATIONS.

     */

    maxit = 30;

    /* DETERMINE WHAT IS TO BE COMPUTED. */

    wantu = IMSLFALSE;

    wantv = IMSLFALSE;

    jobu = mod (*ipath, 100) / 10;

    ncu = *nra;

    if (jobu == 2)

	ncu = imsl_i_min (*nra, *nca);

    if (jobu != 0)

	wantu = IMSLTRUE;

    if (mod (*ipath, 10) != 0)

	wantv = IMSLTRUE;

    /*

     * REDUCE A TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS IN S AND

     * THE SUPER-DIAGONAL ELEMENTS IN E.

     */

    info = 0;

    nct = imsl_i_min (*nra - 1, *nca);

    nrt = imsl_i_max (0, imsl_i_min (*nca - 2, *nra));

    lu = imsl_i_max (nct, nrt);

    if (lu >= 1) {

	for (l = 1; l <= lu; l++) {

	    s[l - 1] = F_ZERO;

	    lp1 = l + 1;

	    if (l <= nct) {

		/*

		 * COMPUTE THE TRANSFORMATION FOR THE L-TH NCA AND PLACE THE

		 * L-TH DIAGONAL IN S(L).

		 */

		s[l - 1] = imsl_snrm2 (*nra - l + 1, A (l - 1, l - 1), 1);

		if (s[l - 1] != F_ZERO) {

		    if (*A (l - 1, l - 1) != F_ZERO)

			s[l - 1] = sign (s[l - 1], *A (l - 1, l - 1));

		    sscal (*nra - l + 1, F_ONE / s[l - 1], A (l - 1, l - 1),

			1);

		    *A (l - 1, l - 1) += F_ONE;

		}

		s[l - 1] = -s[l - 1];

	    }

	    if (*nca >= lp1) {

		if (l <= nct && s[l - 1] != F_ZERO) {

		    /* APPLY THE TRANSFORMATION. */

		    _l0 = *nra - l + 1;

		    _l1 = *nca - lp1 + 1;

		    _f0 = -F_ONE / *A (l - 1, l - 1);

		    _l2 = 1;

		    _f1 = F_ZERO;

		    _l3 = 1;

		    imsl_sgemv ("T", sizeof ("T"), &_l0,

			&_l1, &_f0, A (lp1 - 1, l - 1), lda, A (l - 1, l - 1),

			&_l2, &_f1, &work[*nra], &_l3);

		    imsl_sger (*nra - l + 1, *nca - lp1 + 1, F_ONE, A (l - 1, l - 1),

			1, &work[*nra], 1, A (lp1 - 1, l - 1), *lda);

		}

		/*

		 * PLACE THE L-TH ROW OF A INTO E FOR THE SUBSEQUENT

		 * CALCULATION OF THE ROW TRANSFORMATION.

		 */

		scopy (*nca - lp1 + 1, A (lp1 - 1, l - 1), *lda, &e[lp1 - 1],

		    1);

	    }

	    /*

	     * PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK

	     * MULTIPLICATION.

	     */

	    if (wantu && l <= nct)

		scopy (*nra - l + 1, A (l - 1, l - 1), 1, U (l - 1, l - 1),

		    1);

	    if (l <= nrt) {

		/*

		 * COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE L-TH

		 * SUPER-DIAGONAL IN E(L).

		 */

		e[l - 1] = imsl_snrm2 (*nca - l, &e[lp1 - 1], 1);

		if (e[l - 1] != F_ZERO) {

		    if (e[lp1 - 1] != F_ZERO)

			e[l - 1] = sign (e[l - 1], e[lp1 - 1]);

		    sscal (*nca - l, F_ONE / e[l - 1], &e[lp1 - 1],

			1);

		    e[lp1 - 1] += F_ONE;

		}

		e[l - 1] = -e[l - 1];

		if (lp1 <= *nra && e[l - 1] != F_ZERO) {

		    /* APPLY THE TRANSFORMATION. */

		    _l0 = *nra - l;

		    _l1 = *nca - lp1 + 1;

		    _f0 = F_ONE;

		    _l2 = 1;

		    _f1 = F_ZERO;

		    _l3 = 1;

		    imsl_sgemv ("N", sizeof ("N"), &_l0, &_l1, &_f0, A (lp1 - 1, lp1 - 1),

			lda, &e[lp1 - 1], &_l2, &_f1, &work[lp1 - 1], &_l3);

		    imsl_sger (*nra - l, *nca - lp1 + 1, -F_ONE / e[lp1 - 1],

			&work[lp1 - 1], 1, &e[lp1 - 1], 1, A (lp1 - 1, lp1 - 1),

			*lda);

		}

		/*

		 * PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT BACK

		 * MULTIPLICATION.

		 */

		if (wantv)

		    scopy (*nca - lp1 + 1, &e[lp1 - 1], 1, V (l - 1, lp1 - 1),

			1);

	    }

	}

    }

    /*

     * SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.

     */

    m = imsl_i_min (*nca, *nra + 1);

    nctp1 = nct + 1;

    nrtp1 = nrt + 1;

    if (nct < *nca)

	s[nctp1 - 1] = *A (nctp1 - 1, nctp1 - 1);

    if (*nra < m)

	s[m - 1] = F_ZERO;

    if (nrtp1 < m)

	e[nrtp1 - 1] = *A (m - 1, nrtp1 - 1);

    e[m - 1] = F_ZERO;

    /* IF REQUIRED, GENERATE U. */

    if (wantu) {

	if (ncu >= nctp1) {

	    for (j = nctp1; j <= ncu; j++) {

		sset (*nra, F_ZERO, U (j - 1, 0), 1);

		*U (j - 1, j - 1) = F_ONE;

	    }

	}

	if (nct >= 1) {

	    for (ll = 1; ll <= nct; ll++) {

		l = nct - ll + 1;

		if (s[l - 1] != F_ZERO) {

		    lp1 = l + 1;

		    if (ncu >= lp1) {

			_l0 = *nra - l + 1;

			_l1 = ncu - lp1 + 1;

			_f0 = -F_ONE / *U (l - 1, l - 1);

			_l2 = 1;

			_f1 = F_ZERO;

			_l3 = 1;

			imsl_sgemv ("T", sizeof ("T"), &_l0, &_l1, &_f0, U (lp1 - 1, l - 1),

			    ldu, U (l - 1, l - 1), &_l2, &_f1, &work[*nra], &_l3);

			imsl_sger (*nra - l + 1, ncu - lp1 + 1, F_ONE,

			    U (l - 1, l - 1), 1, &work[*nra], 1, U (lp1 - 1, l - 1),

			    *ldu);

		    }

		    sscal (*nra - l + 1, -F_ONE, U (l - 1, l - 1), 1);

		    *U (l - 1, l - 1) += F_ONE;

		    lm1 = l - 1;

		    if (lm1 >= 1)

			sset (lm1, F_ZERO, U (l - 1, 0), 1);

		}

		else {

		    sset (*nra, F_ZERO, U (l - 1, 0), 1);

		    *U (l - 1, l - 1) = F_ONE;

		}

	    }

	}

    }

    /* IF IT IS REQUIRED, GENERATE V. */

    if (wantv) {

	for (ll = 1; ll <= *nca; ll++) {

	    l = *nca - ll + 1;

	    lp1 = l + 1;

	    if (l <= nrt && e[l - 1] != F_ZERO) {

		_l0 = *nca - 1;

		_l1 = *nca - lp1 + 1;

		_f0 = -F_ONE / *V (l - 1, lp1 - 1);

		_l2 = 1;

		_f1 = F_ZERO;

		_l3 = 1;

		imsl_sgemv ("T", sizeof ("T"), &_l0, &_l1, &_f0, V (lp1 - 1, lp1 - 1),

		    ldv, V (l - 1, lp1 - 1), &_l2, &_f1, &work[*nra], &_l3);

		imsl_sger (*nca - l, *nca - lp1 + 1, F_ONE, V (l - 1, lp1 - 1),

		    1, &work[*nra], 1, V (lp1 - 1, lp1 - 1), *ldv);

	    }

	    sset (*nca, F_ZERO, V (l - 1, 0), 1);

	    *V (l - 1, l - 1) = F_ONE;

	}

    }

    /*

     * MAIN ITERATION LOOP FOR THE SINGULAR VALUES.

     */

    mm = m;

    iter = 0;

L_60:

    ;

    /*

     * QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.

     */

    if (m != 0 && iter <= maxit) {

	/*

	 * THIS SECTION OF THE PROGRAM INSPECTS FOR NEGLIGIBLE ELEMENTS IN

	 * THE S AND E ARRAYS. ON COMPLETION THE VARIABLES KASE AND L ARE SET

	 * AS FOLLOWS.

	 * 

	 * KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M

	 * 

	 * KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M

	 * 

	 * KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND S(L), ..., S(M) ARE

	 * NOT NEGLIGIBLE (QR STEP).

	 * 

	 * KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).

	 */

	for (ll = 1; ll <= m; ll++) {

	    l = m - ll;

	    if (l == 0)

		goto L_80;

	    test = fabs (s[l - 1]) + fabs (s[l]);

	    ztest = test + fabs (e[l - 1]);

	    if (fabs (ztest - test) <= eps * ztest) {

		e[l - 1] = F_ZERO;

		goto L_80;

	    }

	}

L_80:

	;

	if (l == m - 1) {

	    kase = 4;

	}

	else {

	    lp1 = l + 1;

	    mp1 = m + 1;

	    for (lls = lp1; lls <= mp1; lls++) {

		ls = m - lls + lp1;

		if (ls == l)

		    goto L_100;

		test = F_ZERO;

		if (ls != m)

		    test += fabs (e[ls - 1]);

		if (ls != l + 1)

		    test += fabs (e[ls - 2]);

		ztest = test + fabs (s[ls - 1]);

		if (ztest == test) {

		    s[ls - 1] = F_ZERO;

		    goto L_100;

		}

	    }

    L_100:

	    ;

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

	    mm1 = m - 1;

	    f = e[m - 2];

	    e[m - 2] = F_ZERO;

	    for (kk = l; kk <= mm1; kk++) {

		k = mm1 - kk + l;

		t1 = s[k - 1];

		imsl_srotg (&t1, &f, &cs, &sn);

		s[k - 1] = t1;

		if (k != l) {

		    f = -sn * e[k - 2];

		    e[k - 2] *= cs;

		}

		_l0 = 1;

		_l1 = 1;

		if (wantv)

		    imsl_srot (*nca, V (k - 1, 0), _l0, V (m - 1, 0),

			_l1, cs, sn);

	    }

	    /* SPLIT AT NEGLIGIBLE S(L). */

	}

	else if (kase == 2) {

	    f = e[l - 2];

	    e[l - 2] = F_ZERO;

	    for (k = l; k <= m; k++) {

		t1 = s[k - 1];

		imsl_srotg (&t1, &f, &cs, &sn);

		s[k - 1] = t1;

		f = -sn * e[k - 1];

		e[k - 1] *= cs;

		_l0 = 1;

		_l1 = 1;

		if (wantu)

		    imsl_srot (*nra, U (k - 1, 0), _l0, U (l - 2, 0),

			_l1, cs, sn);

	    }

	    /* PERFORM ONE QR STEP. */

	}

	else if (kase == 3) {

	    /* CALCULATE THE SHIFT. */

/*			scale = imsl_f_vmax(fabs(s[m - 1]), fabs(s[m - 2]), fabs(e[m - 2]),

				      fabs(s[l - 1]), fabs(e[l - 1]), FEND);

*/

	    scale = imsl_f_vmax (5, fabs (s[m - 1]), fabs (s[m - 2]), fabs (e[m - 2]),

		fabs (s[l - 1]), fabs (e[l - 1]));

	    sm = s[m - 1] / scale;

	    smm1 = s[m - 2] / scale;

	    emm1 = e[m - 2] / scale;

	    sl = s[l - 1] / scale;

	    el = e[l - 1] / scale;

	    b = ((smm1 + sm) * (smm1 - sm) + imsl_fi_power (emm1, 2)) / F_TWO;

	    c = imsl_fi_power (sm * emm1, 2);

	    shift = F_ZERO;

	    if (b != F_ZERO || c != F_ZERO) {

		shift = sqrt (imsl_fi_power (b, 2) + c);

		if (b < F_ZERO)

		    shift = -shift;

		shift = c / (b + shift);

	    }

	    f = (sl + sm) * (sl - sm) - shift;

	    g = sl * el;

	    /* CHASE ZEROS. */

	    mm1 = m - 1;

	    for (k = l; k <= mm1; k++) {

		imsl_srotg (&f, &g, &cs, &sn);

		if (k != l)

		    e[k - 2] = f;

		f = cs * s[k - 1] + sn * e[k - 1];

		e[k - 1] = cs * e[k - 1] - sn * s[k - 1];

		g = sn * s[k];

		s[k] *= cs;

		_l0 = 1;

		_l1 = 1;

		if (wantv)

		    imsl_srot (*nca, V (k - 1, 0), _l0, V (k, 0), _l1, cs, sn);

		imsl_srotg (&f, &g, &cs, &sn);

		s[k - 1] = f;

		f = cs * e[k - 1] + sn * s[k];

		s[k] = -sn * e[k - 1] + cs * s[k];

		g = sn * e[k];

		e[k] *= cs;

		_l0 = 1;

		_l1 = 1;

		if (wantu && k < *nra)

		    imsl_srot (*nra, U (k - 1, 0), _l0, U (k, 0), _l1,

			cs, sn);

	    }

	    e[m - 2] = f;

	    iter += 1;

	    /* CONVERGENCE. */

	}

	else if (kase == 4) {

	    /* MAKE THE SINGULAR VALUE POSITIVE. */

	    if (s[l - 1] < F_ZERO) {

		s[l - 1] = -s[l - 1];

		if (wantv)

		    sscal (*nca, -F_ONE, V (l - 1, 0), 1);

	    }

	    /* ORDER THE SINGULAR VALUE. */

    L_140:

	    if (l == mm)

		goto L_150;

	    if (s[l - 1] >= s[l])

		goto L_150;

	    t = s[l - 1];

	    s[l - 1] = s[l];

	    s[l] = t;

	    if (wantv && l < *nca)

		sswap (*nca, V (l - 1, 0), 1, V (l, 0), 1);

	    if (wantu && l < *nra)

		sswap (*nra, U (l - 1, 0), 1, U (l, 0), 1);

	    l += 1;

	    goto L_140;

    L_150:

	    ;

	    iter = 0;

	    m -= 1;

	}

	goto L_60;

    }

    if (iter > maxit)

	info = m;

    if (info == 0) {

	/* ESTIMATE THE RANK OF A */

	ls = imsl_i_min (*nca + 1, *nra);

	minmn = imsl_i_min (*nca, *nra);

	*irank = minmn;

	for (i = 1; i <= minmn; i++) {

	    if (fabs (s[i - 1]) <= stol) {

		*irank = i - 1;

		goto L_170;

	    }

	}

L_170:

	;

    }

    else {

	istart = info + 1;

	iend = imsl_i_min (*nra, *nca);

	imsl_e1sti (1, istart);

	imsl_e1sti (2, iend);



/*		(4, 1, "Convergence can only be obtained for the %(i1), ... , %(i2) singular values and their corresponding singular vectors.");

*/

	imsl_ermes (IMSL_FATAL, IMSL_SLOWCONVERGENT_MATRIX);

    }

    imsl_e1pop ("l_l3vrr");



    return;

}				/* end of function */

#undef A

#undef U

#undef V



#endif /*old l3vrr */



/*Translated by FOR_C++, v0.1, on 11/29/91 at 12:54:47 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 11/29/91 at 12:54:43

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  L3VRR/DL3VRR (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    November 24, 1991



    Purpose:    Compute the singular value decomposition of a real

                matrix.



    Usage:      CALL L3VRR (NRA, NCA, A, LDA, IPATH, TOL, IRANK,

                            S, U, LDU, V, LDV, E, WORK)



    Arguments:  See LSVRR/DLSVRR



    Remarks:    See LSVRR/DLSVRR.



    Chapter:    MATH/LIBRARY Linear Systems



    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_l3vrr (Mint *nra, Mint *nca, Mfloat *a, Mint *lda, Mint *ipath,

                Mfloat *tol, Mint *irank, Mfloat *s, Mfloat *u, Mint *ldu,

                Mfloat *v, Mint *ldv, Mfloat *e, Mfloat *work)

#else

static void l_l3vrr (nra, nca, a, lda, ipath, tol, irank, s, u,

                ldu, v, ldv, e, work)

    Mint        *nra, *nca;

    Mfloat      *a;

    Mint        *lda, *ipath;

    Mfloat      *tol;

    Mint        *irank;

    Mfloat       s[], *u;

    Mint        *ldu;

    Mfloat      *v;

    Mint        *ldv;

    Mfloat       e[], work[];

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

#define U(I_,J_)	(u+(I_)*(*ldu)+(J_))

#define V(I_,J_)	(v+(I_)*(*ldv)+(J_))

    Mint   usean, wantu, wantv;

    Mint         _l0, _l1, _l2, _l3, i, iend, info, istart, iter, j, jobu,

                k, kase, kk, l, ll, lls, lm1, lp1, ls, lu, m, maxit, minmn,

                mm, mm1, mp1, nct, nctp1, ncu, nrt, nrtp1;

    Mfloat       _f0, anorm, b, big, c, cs, el, emm1, eps, f, g, one, scale, shift,

                sl, sm, small, smm1, sn, stol, t, t1, test, zero, ztest;





    imsl_e1psh ("L3VRR");

    zero = 0.0e0;

    one = 1.0e0;

    small = imsl_amach (1);

    big = imsl_amach (2);

    if (small * big < one)

	small = one / big;

    eps = imsl_amach (4);

    /*

     * IF TOL.LT.0, COMPUTE THE INFINITY NORM OF A AND SINGULAR VALUE

     * TOLERANCE FOR RANK ESTIMATION

     */

    anorm = zero;

    for (i = 1; i <= *nra; i++) {

	anorm = imsl_f_max (anorm, imsl_sasum (*nca, A (0, i - 1), *lda));

    }

    usean = FALSE;

    if (*tol < zero) {

	stol = fabs (*tol) * anorm;

    }

    else {

	stol = *tol;

    }

    /*

     * SET THE MAXIMUM NUMBER OF ITERATIONS.

     */

    maxit = 30;

    /* DETERMINE WHAT IS TO BE COMPUTED. */

    wantu = FALSE;

    wantv = FALSE;

    jobu = mod (*ipath, 100) / 10;

    ncu = *nra;

    if (jobu == 2)

	ncu = imsl_i_min (*nra, *nca);

    if (jobu != 0)

	wantu = TRUE;

    if (mod (*ipath, 10) != 0)

	wantv = TRUE;

    /*

     * REDUCE A TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS IN S AND

     * THE SUPER-DIAGONAL ELEMENTS IN E.

     */

    info = 0;

    nct = imsl_i_min (*nra - 1, *nca);

    nrt = imsl_i_max (0, imsl_i_min (*nca - 2, *nra));

    lu = imsl_i_max (nct, nrt);

    if (lu >= 1) {

	for (l = 1; l <= lu; l++) {

	    s[l - 1] = zero;

	    lp1 = l + 1;

	    if (l <= nct) {

		/*

		 * COMPUTE THE TRANSFORMATION FOR THE L-TH NCA AND PLACE THE

		 * L-TH DIAGONAL IN S(L).

		 */

		s[l - 1] = imsl_snrm2 (*nra - l + 1, A (l - 1, l - 1), 1);

		if (s[l - 1] >= small) {

		    if (*A (l - 1, l - 1) != zero)

			s[l - 1] = sign (s[l - 1], *A (l - 1, l - 1));

		    sscal (*nra - l + 1, one / s[l - 1], A (l - 1, l - 1),

			1);

		    *A (l - 1, l - 1) += one;

		}

		s[l - 1] = -s[l - 1];

	    }



	    if (*nca >= lp1) {

		if (l <= nct && fabs (s[l - 1]) >= small) {

		    /* APPLY THE TRANSFORMATION. */

		    _l0 = *nra - l + 1;

		    _l1 = *nca - lp1 + 1;

		    _f0 = -one / *A (l - 1, l - 1);

		    _l2 = 1;

		    _l3 = 1;

		    imsl_sgemv ("T", sizeof ("T"), &_l0,

			&_l1, &_f0,

			A (lp1 - 1, l - 1), lda, A (l - 1, l - 1), &_l2,

			&zero, &work[*nra], &_l3);

		    imsl_sger (*nra - l + 1, *nca - lp1 + 1, one, A (l - 1, l - 1),

			1, &work[*nra], 1, A (lp1 - 1, l - 1), *lda);

		}

		/*

		 * PLACE THE L-TH ROW OF A INTO E FOR THE SUBSEQUENT

		 * CALCULATION OF THE ROW TRANSFORMATION.

		 */

		scopy (*nca - lp1 + 1, A (lp1 - 1, l - 1), *lda, &e[lp1 - 1],

		    1);

	    }

	    /*

	     * PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK

	     * MULTIPLICATION.

	     */

	    if (wantu && l <= nct)

		scopy (*nra - l + 1, A (l - 1, l - 1), 1, U (l - 1, l - 1),

		    1);

	    if (l <= nrt) {

		/*

		 * COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE L-TH

		 * SUPER-DIAGONAL IN E(L).

		 */

		e[l - 1] = imsl_snrm2 (*nca - l, &e[lp1 - 1], 1);

		if (e[l - 1] >= small) {

		    if (e[lp1 - 1] != zero)

			e[l - 1] = sign (e[l - 1], e[lp1 - 1]);

		    sscal (*nca - l, one / e[l - 1], &e[lp1 - 1], 1);

		    e[lp1 - 1] += one;

		}

		e[l - 1] = -e[l - 1];

		if (lp1 <= *nra && fabs (e[l - 1]) >= small) {

		    /* APPLY THE TRANSFORMATION. */

		    _l0 = *nra - l;

		    _l1 = *nca - lp1 + 1;

		    _l2 = 1;

		    _l3 = 1;

		    imsl_sgemv ("N", sizeof ("N"), &_l0, &_l1, 

			    &one, A (lp1 - 1, lp1 - 1), lda, &e[lp1 - 1],

			&_l2, &zero, &work[lp1 - 1], &_l3);

		    imsl_sger (*nra - l, *nca - lp1 + 1, -one / e[lp1 - 1],

			&work[lp1 - 1], 1, &e[lp1 - 1], 1, A (lp1 - 1, lp1 - 1),

			*lda);

		}

		/*

		 * PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT BACK

		 * MULTIPLICATION.

		 */

		if (wantv)

		    scopy (*nca - lp1 + 1, &e[lp1 - 1], 1, V (l - 1, lp1 - 1),

			1);

	    }

	}

    }

    /*

     * SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.

     */

    m = imsl_i_min (*nca, *nra + 1);

    nctp1 = nct + 1;

    nrtp1 = nrt + 1;

    if (nct < *nca)

	s[nctp1 - 1] = *A (nctp1 - 1, nctp1 - 1);

    if (*nra < m)

	s[m - 1] = zero;

    if (nrtp1 < m)

	e[nrtp1 - 1] = *A (m - 1, nrtp1 - 1);

    e[m - 1] = zero;

    /* IF REQUIRED, GENERATE U. */

    if (wantu) {

	if (ncu >= nctp1) {

	    for (j = nctp1; j <= ncu; j++) {

		sset (*nra, zero, U (j - 1, 0), 1);

		*U (j - 1, j - 1) = one;

	    }

	}

	if (nct >= 1) {

	    for (ll = 1; ll <= nct; ll++) {

		l = nct - ll + 1;

		if (fabs (s[l - 1]) >= small) {

		    lp1 = l + 1;

		    if (ncu >= lp1) {

			_l0 = *nra - l + 1;

			_l1 = ncu - lp1 + 1;

			_f0 = -one / *U (l - 1, l - 1);

			_l2 = 1;

			_l3 = 1;

			imsl_sgemv ("T", sizeof ("T"), &_l0, 

				&_l1, &_f0, 

				U (lp1 - 1, l - 1), ldu,

			    U (l - 1, l - 1), &_l2, &zero, &work[*nra],

			    &_l3);

			imsl_sger (*nra - l + 1, ncu - lp1 + 1, one, U (l - 1, l - 1),

			    1, &work[*nra], 1, U (lp1 - 1, l - 1), *ldu);

		    }

		    sscal (*nra - l + 1, -one, U (l - 1, l - 1), 1);

		    *U (l - 1, l - 1) += one;

		    lm1 = l - 1;

		    if (lm1 >= 1)

			sset (lm1, zero, U (l - 1, 0), 1);

		}

		else {

		    sset (*nra, zero, U (l - 1, 0), 1);

		    *U (l - 1, l - 1) = one;

		}

	    }

	}

    }

    /* IF IT IS REQUIRED, GENERATE V. */

    if (wantv) {

	for (ll = 1; ll <= *nca; ll++) {

	    l = *nca - ll + 1;

	    lp1 = l + 1;

	    if (l <= nrt && fabs (e[l - 1]) >= small) {

		_l0 = *nca - l;

		_l1 = *nca - lp1 + 1;

		_f0 = -one / *V (l - 1, lp1 - 1);

		_l2 = 1;

		_l3 = 1;

		imsl_sgemv ("T", sizeof ("T"), &_l0, &_l1, 

			&_f0, V (lp1 - 1, lp1 - 1),

		    ldv, V (l - 1, lp1 - 1), &_l2, &zero, &work[*nra],

		    &_l3);

		imsl_sger (*nca - l, *nca - lp1 + 1, one, V (l - 1, lp1 - 1),

		    1, &work[*nra], 1, V (lp1 - 1, lp1 - 1), *ldv);

	    }

	    sset (*nca, zero, V (l - 1, 0), 1);

	    *V (l - 1, l - 1) = one;

	}

    }

    /*

     * MAIN ITERATION LOOP FOR THE SINGULAR VALUES.

     */

    mm = m;

    iter = 0;

L_60:

    ;

    /*

     * QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.

     */

    if (m != 0 && iter <= maxit) {

	/*

	 * THIS SECTION OF THE PROGRAM INSPECTS FOR NEGLIGIBLE ELEMENTS IN

	 * THE S AND E ARRAYS. ON COMPLETION THE VARIABLES KASE AND L ARE SET

	 * AS FOLLOWS.

	 * 

	 * KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M

	 * 

	 * KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M

	 * 

	 * KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND S(L), ..., S(M) ARE

	 * NOT NEGLIGIBLE (QR STEP).

	 * 

	 * KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).

	 */

	for (ll = 1; ll <= m; ll++) {

	    l = m - ll;

	    if (l == 0)

		goto L_80;

	    test = fabs (s[l - 1]) + fabs (s[l]);

	    ztest = test + fabs (e[l - 1]);

	    /*

	     * This is a shift from the default test for small to a less

	     * stringent test whenever there are many iterations (maxit-1)

	     * without convergence.  From that point the less stringent test

	     * (eps*matrix norm) is used.

	     */

	    usean = usean || (iter == maxit - 1);

	    if (fabs (ztest - test) <= eps * ztest || (usean && fabs (e[l - 1]) <=

		    eps * anorm)) {

		e[l - 1] = zero;

		goto L_80;

	    }

	}

L_80:

	;

	if (l == m - 1) {

	    kase = 4;

	}

	else {

	    lp1 = l + 1;

	    mp1 = m + 1;

	    for (lls = lp1; lls <= mp1; lls++) {

		ls = m - lls + lp1;

		if (ls == l)

		    goto L_100;

		test = zero;

		if (ls != m)

		    test += fabs (e[ls - 1]);

		if (ls != l + 1)

		    test += fabs (e[ls - 2]);

		ztest = test + fabs (s[ls - 1]);

		if (ztest == test) {

		    s[ls - 1] = zero;

		    goto L_100;

		}

	    }

    L_100:

	    ;

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

	    mm1 = m - 1;

	    f = e[m - 2];

	    e[m - 2] = zero;

	    for (kk = l; kk <= mm1; kk++) {

		k = mm1 - kk + l;

		t1 = s[k - 1];

		imsl_srotg (&t1, &f, &cs, &sn);

		s[k - 1] = t1;

		if (k != l) {

		    f = -sn * e[k - 2];

		    e[k - 2] *= cs;

		}

		if (wantv) {

		    _l0 = 1;

		    _l1 = 1;

		    imsl_srot (*nca, V (k - 1, 0), _l0, V (m - 1, 0),

			_l1, cs, sn);

		}

	    }

	    /* SPLIT AT NEGLIGIBLE S(L). */

	}

	else if (kase == 2) {

	    f = e[l - 2];

	    e[l - 2] = zero;

	    for (k = l; k <= m; k++) {

		t1 = s[k - 1];

		imsl_srotg (&t1, &f, &cs, &sn);

		s[k - 1] = t1;

		f = -sn * e[k - 1];

		e[k - 1] *= cs;

		if (wantu) {

		    _l0 = 1;

		    _l1 = 1;

		    imsl_srot (*nra, U (k - 1, 0), _l0, U (l - 2, 0),

			_l1, cs, sn);

		}

	    }

	    /* PERFORM ONE QR STEP. */

	}

	else if (kase == 3) {

	    /* CALCULATE THE SHIFT. */

	    scale = imsl_f_vmax (5, fabs (s[m - 1]), fabs (s[m - 2]), fabs (e[m - 2]),

		fabs (s[l - 1]), fabs (e[l - 1]));

	    sm = s[m - 1] / scale;

	    smm1 = s[m - 2] / scale;

	    emm1 = e[m - 2] / scale;

	    sl = s[l - 1] / scale;

	    el = e[l - 1] / scale;

	    b = ((smm1 + sm) * (smm1 - sm) + imsl_fi_power (emm1, 2)) / 2.0e0;

	    c = imsl_fi_power (sm * emm1, 2);

	    shift = zero;

	    if (b != zero || c != zero) {

		shift = sqrt (imsl_fi_power (b, 2) + c);

		if (b < zero)

		    shift = -shift;

		shift = c / (b + shift);

	    }

	    f = (sl + sm) * (sl - sm) - shift;

	    g = sl * el;

	    /* CHASE ZEROS. */

	    mm1 = m - 1;

	    for (k = l; k <= mm1; k++) {

		imsl_srotg (&f, &g, &cs, &sn);

		if (k != l)

		    e[k - 2] = f;

		f = cs * s[k - 1] + sn * e[k - 1];

		e[k - 1] = cs * e[k - 1] - sn * s[k - 1];

		g = sn * s[k];

		s[k] *= cs;

		if (wantv) {

		    _l0 = 1;

		    _l1 = 1;

		    imsl_srot (*nca, V (k - 1, 0), _l0, V (k, 0), _l1,

			cs, sn);

		}

		imsl_srotg (&f, &g, &cs, &sn);

		s[k - 1] = f;

		f = cs * e[k - 1] + sn * s[k];

		s[k] = -sn * e[k - 1] + cs * s[k];

		g = sn * e[k];

		e[k] *= cs;

		if (wantu && k < *nra) {

		    _l0 = 1;

		    _l1 = 1;

		    imsl_srot (*nra, U (k - 1, 0), _l0, U (k, 0), _l1,

			cs, sn);

		}

	    }

	    e[m - 2] = f;

	    iter += 1;

	    /* CONVERGENCE. */

	}

	else if (kase == 4) {

	    /* MAKE THE SINGULAR VALUE POSITIVE. */

	    if (s[l - 1] < zero) {

		s[l - 1] = -s[l - 1];

		if (wantv)

		    sscal (*nca, -one, V (l - 1, 0), 1);

	    }

	    /* ORDER THE SINGULAR VALUE. */

    L_140:

	    if (l == mm)

		goto L_150;

	    if (s[l - 1] >= s[l])

		goto L_150;

	    t = s[l - 1];

	    s[l - 1] = s[l];

	    s[l] = t;

	    if (wantv && l < *nca)

		sswap (*nca, V (l - 1, 0), 1, V (l, 0), 1);

	    if (wantu && l < *nra)

		sswap (*nra, U (l - 1, 0), 1, U (l, 0), 1);

	    l += 1;

	    goto L_140;

    L_150:

	    ;

	    iter = 0;

	    m -= 1;

	}

	goto L_60;

    }

    if (iter > maxit)

	info = m;

    if (info == 0) {

	/* ESTIMATE THE RANK OF A */

	ls = imsl_i_min (*nca + 1, *nra);

	minmn = imsl_i_min (*nca, *nra);

	*irank = minmn;

	for (i = 1; i <= minmn; i++) {

	    if (fabs (s[i - 1]) <= stol) {

		*irank = i - 1;

		goto L_170;

	    }

	}

L_170:

	;

    }

    else {

	istart = info + 1;

	iend = imsl_i_min (*nra, *nca);

	imsl_e1sti (1, istart);

	imsl_e1sti (2, iend);

/*

	(4, 1, "Convergence can only be obtained for the %(i1), ... , 

	%(i2) singular values and their corresponding singular vectors.");

*/

	imsl_ermes (IMSL_FATAL, IMSL_SLOWCONVERGENT_MATRIX);

    }

    imsl_e1pop ("L3VRR");



    return;

}				/* end of function */



#undef A

#undef U

#undef V

