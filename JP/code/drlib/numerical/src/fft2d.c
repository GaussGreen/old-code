#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
void        imsl_f3tcf (Mint *n, Mfloat *c, Mfloat *ch, Mfloat *wa, Mfloat
               *imsl_fac);
void        imsl_f3tcb (Mint *n, Mfloat *c, Mfloat *ch, Mfloat *wa, Mfloat
               *imsl_fac);
#else
void        imsl_f3tcf ();
void        imsl_f3tcb ();
#endif
#ifdef ANSI
static VA_LIST_HACK l_fft_2d_complex (Mint n, Mint m, Mf_complex *seq, va_list argptr);
static void l_f2t2d (Mint *nra, Mint *nca, Mf_complex *a, Mint *lda,
                Mf_complex *coef, Mint *ldcoef, Mfloat *wff1,
                Mfloat *wff2, Mf_complex *cwk, Mfloat *cpy);
static void l_f2t2b (Mint *nrcoef, Mint *nccoef, Mf_complex *a, Mint *lda,
                Mf_complex *coef, Mint *ldcoef, Mfloat *wff1,
                Mfloat *wff2, Mf_complex *cwk, Mfloat *cpy);
#else
static VA_LIST_HACK l_fft_2d_complex ();
static void l_f2t2d ();
static void l_f2t2b ();
#endif

static Mf_complex *lv_coef;
#ifdef ANSI
Mf_complex *imsl_c_fft_2d_complex (Mint n, Mint m, Mf_complex *seq,...)
#else
Mf_complex *imsl_c_fft_2d_complex (n, m, seq, va_alist)
    Mint        n;
    Mint        m;
    Mf_complex *seq;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, seq);
#ifdef DOUBLE
    imsl_e1psh ("imsl_z_fft_2d_complex");
#else
    imsl_e1psh ("imsl_c_fft_2d_complex");
#endif
    lv_coef = NULL;
    IMSL_CALL (l_fft_2d_complex (n, m, seq, argptr));
    va_end (argptr);
#ifdef DOUBLE
    imsl_e1pop ("imsl_z_fft_2d_complex");
#else
    imsl_e1pop ("imsl_c_fft_2d_complex");
#endif
    return lv_coef;
}


#ifdef ANSI
static VA_LIST_HACK l_fft_2d_complex (Mint n, Mint m, Mf_complex *seq, va_list argptr)
#else
static VA_LIST_HACK l_fft_2d_complex (n, m, seq, argptr)
    Mint        n;
    Mint        m;
    Mf_complex *seq;
    va_list     argptr;
#endif
{
    Mint        code;
    Mint        arg_number = 3;
    Mint        forward = 1;
    Mint        result_user = 0;
    Mint        backward = 0;
    Mint        p_col_dim = m;
    Mint        q_col_dim = m;
    Mfloat     *params1 = NULL;
    Mfloat     *params2 = NULL;
    Mfloat     *copy = NULL;
    Mf_complex *work = NULL;
    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_P_COL_DIM:
	    p_col_dim = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_Q_COL_DIM:
	    q_col_dim = va_arg (argptr, Mint);
	    arg_number++;
	    break;
	case IMSL_BACKWARD:
	    forward = 0;
	    backward = 1;
	    break;
	case IMSL_RETURN_USER:
	    result_user = 1;
	    lv_coef = va_arg (argptr, Mf_complex *);
	    arg_number++;
	    break;
	case 0:
	    break;
	default:
	    /* Argument number %(I2) is an unknown */
	    /* optional argument %(I1). */
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, arg_number);
	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    break;
	}
    }

    if (imsl_n1rty (0))
	goto RETURN;

    if (m <= 0) {
	imsl_e1sti (1, m);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NCA_GT_ZERO);
    }

    if (n <= 0) {
	imsl_e1sti (1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NRA_GT_ZERO);
    }


    if (result_user && (lv_coef == NULL)) {
	imsl_e1stl (1, "result");
	imsl_e1stl (2, "IMSL_RETURN_USER");
	imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
	goto RETURN;
    }

    if (seq == NULL) {
	imsl_e1stl (1, "seq");
	imsl_ermes (IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
	goto RETURN;
    }

    if (imsl_n1rty (0))
	goto RETURN;

    params1 = imsl_c_fft_complex_init (m);
    params2 = imsl_c_fft_complex_init (n);
    if (imsl_n1rty (0))
	goto FREE_SPACE;

    work = (Mf_complex *) imsl_malloc (n * sizeof (*work));
    copy = (Mfloat *) imsl_malloc (2 * ((n > m) ? n : m) * sizeof (*copy));


    if (!result_user) {
	lv_coef = (Mf_complex *) imsl_malloc (n * m * sizeof (*lv_coef));
    }

    if (lv_coef == NULL || work == NULL || copy == NULL) {
	/* Not enough memory, with %(L1) = %(I1). */
	imsl_e1sti (1, n);
	imsl_e1stl (1, "n");
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto FREE_SPACE;
    }

    if (forward) {
	l_f2t2d (&m, &n, seq, &p_col_dim, lv_coef, &q_col_dim,
	    params1, params2, work, copy);
    }

    if (backward) {
	l_f2t2b (&m, &n, lv_coef, &p_col_dim, seq, &q_col_dim,
	    params1, params2, work, copy);
    }
FREE_SPACE:
    if (params1 != NULL)
	imsl_free (params1);
    if (params2 != NULL)
	imsl_free (params2);
    if (work != NULL)
	imsl_free (work);
    if (copy != NULL)
	imsl_free (copy);
RETURN:
    return (argptr);
}

/* Structured by FOR_STRUCT, v0.2, on 08/30/90 at 11:53:35
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F2T2D/DF2T2D (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    May 1, 1990

    Purpose:    Compute Fourier coefficients of a d_complex periodic
                two-dimensional array.

    Usage:      CALL F2T2D (NRA, NCA, A, LDA, COEF, LDCOEF, WFF1, WFF2,
                            CWK, CPY)

    Arguments:  (See FFT2D)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* CWK is not used here, but leave the calling sequence intact */
#ifdef ANSI
static void l_f2t2d (Mint *nra, Mint *nca, Mf_complex *a, Mint *lda,
                Mf_complex *coef, Mint *ldcoef, Mfloat *wff1,
                Mfloat *wff2, Mf_complex *cwk, Mfloat *cpy)
#else
static void l_f2t2d (nra, nca, a, lda, coef, ldcoef, wff1, wff2,
                cwk, cpy)
    Mint       *nra, *nca;
    Mf_complex *a;
    Mint       *lda;
    Mf_complex *coef;
    Mint       *ldcoef;
    Mfloat      wff1[], wff2[];
    Mf_complex  cwk[];
    Mfloat      cpy[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define COEF(I_,J_)	(coef+(I_)*(*ldcoef)+(J_))
    Mint        i, j;


    imsl_e1psh ("l_f2t2d");
    /* CHECK NRA */
    if (*nra <= 0) {
	imsl_e1sti (1, *nra);

	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NRA_GT_ZERO);
    }
    /* CHECK NCA */
    if (*nca <= 0) {
	imsl_e1sti (1, *nca);
	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NCA_GT_ZERO);
    }
    if (imsl_n1rty (0) != 0)
	goto L_9000;

    /* NOTE: The folowing is not a pure translation, but is */
    /* modified to follow the order of matrix access */
    /* exactly as per th e fortran version.  This does */
    /* not affect the final results, but does make   */
    /* debugging possible.                           */
    /*
     * TAKE THE TRANSFORM OF THE ROWS OF A, PUTTING THE RESULTS IN THE ROWS
     * OF COEF
     */
    for (i = 1; i <= *nca; i++) {
	for (j = 1; j <= *nra; j++) {
	    cpy[j * 2 - 2] = imsl_fc_convert (*A (i - 1, j - 1));
	    cpy[j * 2 - 1] = imsl_c_aimag (*A (i - 1, j - 1));
	}
	if (*nra > 1) {
	    imsl_f3tcf (nra, cpy, wff1, &wff1[*nra * 2], &wff1[*nra * 4]);
	}
	for (j = 1; j <= *nra; j++) {
	    *COEF (i - 1, j - 1) = imsl_cf_convert (cpy[j * 2 - 2], cpy[j * 2 - 1]);
	}
    }
    /* TAKE TRANSFORM OF THE COLUMNS OF COEF */
    for (i = 1; i <= *nra; i++) {
	for (j = 1; j <= *nca; j++) {
	    cpy[j * 2 - 2] = imsl_fc_convert (*COEF (j - 1, i - 1));
	    cpy[j * 2 - 1] = imsl_c_aimag (*COEF (j - 1, i - 1));
	}
	if (*nca > 1) {
	    imsl_f3tcf (nca, cpy, wff2, &wff2[*nca * 2], &wff2[*nca * 4]);
	}
	for (j = 1; j <= *nca; j++) {
	    *COEF (j - 1, i - 1) = imsl_cf_convert (cpy[j * 2 - 2], cpy[j * 2 - 1]);
	}
    }

L_9000:
    ;
    imsl_e1pop ("l_f2t2d");
    return;
}				/* end of function */
#undef  A
#undef  COEF
/*----------------------------------------------------------------------- */

/*  IMSL Name:  F2T2B/DF2T2B (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 24, 1990

    Purpose:    Compute the inverse Fourier transform of a d_complex
                periodic two-dimensional array.

    Usage:      CALL F2T2B (NRCOEF, NCCOEF, A, LDA, COEF, LDCOEF, WFF1,
                            WFF2, CWK, CPY)

    Arguments:  (See FFT2B)

    Chapter:    MATH/LIBRARY Transforms

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* CWK is not used here, but leave the calling sequence intact. */
#ifdef ANSI
static void l_f2t2b (Mint *nrcoef, Mint *nccoef, Mf_complex *a, Mint *lda,
                Mf_complex *coef, Mint *ldcoef, Mfloat *wff1,
                Mfloat *wff2, Mf_complex *cwk, Mfloat *cpy)
#else
static void l_f2t2b (nrcoef, nccoef, a, lda, coef, ldcoef, wff1,
                wff2, cwk, cpy)
    Mint       *nrcoef, *nccoef;
    Mf_complex *a;
    Mint       *lda;
    Mf_complex *coef;
    Mint       *ldcoef;
    Mfloat      wff1[], wff2[];
    Mf_complex  cwk[];
    Mfloat      cpy[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define COEF(I_,J_)	(coef+(I_)*(*ldcoef)+(J_))
    Mint        i, j;


    imsl_e1psh ("l_f2t2b");
    /* CHECK NRCOEF */
    if (*nrcoef < 1) {
	imsl_e1sti (1, *nrcoef);

	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NRCOEF_GT_ZERO);
	goto L_9000;
    }
    /* CHECK NCCOEF */
    if (*nccoef < 1) {
	imsl_e1sti (1, *nccoef);

	imsl_ermes (IMSL_TERMINAL, IMSL_NEED_NCCOEF_GT_ZERO);
	goto L_9000;
    }
    /* CHECK LDA */
    if (*lda < *nrcoef) {
	imsl_e1sti (1, *lda);
	imsl_e1sti (2, *nrcoef);

	imsl_ermes (IMSL_TERMINAL, IMSL_LDA_IS_LT_NRCOEF);
	goto L_9000;
    }
    /* CHECK LDCOEF */
    if (*ldcoef < *nrcoef) {
	imsl_e1sti (1, *ldcoef);
	imsl_e1sti (2, *nrcoef);

	imsl_ermes (IMSL_TERMINAL, IMSL_LDCOEF_IS_LT_NRCOEF);
	goto L_9000;
    }


    /*
     * TAKE THE TRANSFORM OF THE ROWS OF COEF, PUTTING THE RESULTS IN THE
     * ROWS OF A
     */
    if (*nccoef > 1) {
	for (i = 1; i <= *nrcoef; i++) {
	    /* COPY COEF TO CPY */
	    for (j = 1; j <= *nccoef; j++) {
		cpy[j * 2 - 2] = imsl_fc_convert (*COEF (j - 1, i - 1));
		cpy[j * 2 - 1] = imsl_c_aimag (*COEF (j - 1, i - 1));
	    }
	    imsl_f3tcb (nccoef, cpy, wff2, &wff2[*nccoef * 2], &wff2[*nccoef * 4]);
	    /* COPY CPY TO A */
	    for (j = 1; j <= *nccoef; j++) {
		*A (j - 1, i - 1) = imsl_cf_convert (cpy[j * 2 - 2], cpy[j * 2 - 1]);
	    }
	}
    }
    else {
	for (i = 1; i <= *nrcoef; i++) {
	    *A (0, i - 1) = *COEF (0, i - 1);
	}
    }
    /*
     * TAKE THE TRANSFORM OF THE COLUMNS OF COEF
     */
    if (*nrcoef > 1) {
	for (i = 1; i <= *nccoef; i++) {
	    /* COPY COEF TO CPY */
	    for (j = 1; j <= *nrcoef; j++) {
		cpy[j * 2 - 2] = imsl_fc_convert (*A (i - 1, j - 1));
		cpy[j * 2 - 1] = imsl_c_aimag (*A (i - 1, j - 1));
	    }

	    imsl_f3tcb (nrcoef, cpy, wff1, &wff1[*nrcoef * 2], &wff1[*nrcoef * 4]);
	    /* COPY CPY TO A */
	    for (j = 1; j <= *nrcoef; j++) {
		*A (i - 1, j - 1) = imsl_cf_convert (cpy[j * 2 - 2], cpy[j * 2 - 1]);
	    }
	}
    }
    else {
	for (i = 1; i <= *nccoef; i++) {
	    *A (i - 1, 0) = *COEF (i - 1, 0);
	}
    }
L_9000:
    ;
    imsl_e1pop ("l_f2t2b");
    return;
}				/* end of function */
#undef A
#undef COEF
