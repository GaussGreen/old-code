#define ADR(t,x)    ( t = x, &t )



#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef TRUE
#undef TRUE
#endif
#ifdef FALSE
#undef FALSE
#endif

#ifdef ANSI
static VA_LIST_HACK l_eig_gen(Mint n, Mf_complex *a, va_list argptr);
static void l_e3ccg(Mint *n, Mf_complex *a, Mint *lda, Mint *low,
                        Mint *igh, Mfloat scale[]);
static void l_e3cch (Mint *n, Mint *low, Mint *igh, Mf_complex ort[],
                        Mf_complex *a, Mint *lda, Mf_complex eval[],
                        Mint *vector, Mf_complex *evec, Mint *ldevec,
                        Mf_complex work[]);
static void l_e3lcg (Mint *n, Mf_complex *a, Mint *lda, Mf_complex eval[],
                        Mf_complex *acopy, Mfloat rwk[], Mf_complex *cwk,
                        Mint iwk[]);
static void l_e4ccg (Mint *n, Mint *low, Mint *igh, Mf_complex *a,
                        Mint *lda, Mf_complex ort[], Mf_complex work[]);
static void l_e5ccg (Mint *n, Mint *low, Mint *igh, Mfloat scale[],
                        Mint *nvec, Mf_complex *evec, Mint *ldevec);
static void l_ccgcg (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *b,
                        Mint *ldb);
static void l_e6ccg (Mint *n, Mf_complex *a, Mint *lda, Mf_complex eval[],
                        Mf_complex *evec, Mint *ldevec, Mf_complex *acopy,
                        Mfloat rwk[], Mf_complex *cwk, Mint iwk[]);
#else
static VA_LIST_HACK	l_eig_gen();
static void 	l_e3ccg();
static void 	l_e3cch();
static void	l_e3lcg();
static void 	l_e4ccg();
static void 	l_e5ccg();
static void 	l_ccgcg();
static void 	l_e6ccg();
#endif

static Mf_complex	*lv_eval;

#ifdef ANSI
Mf_complex *imsl_c_eig_gen(Mint n, Mf_complex *a, ...)
#else
Mf_complex *imsl_c_eig_gen(n, a, va_alist)
    Mint	n;
    Mf_complex	*a;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr,a);

#ifdef DOUBLE
    imsl_e1psh("imsl_z_eig_gen");
#else
    imsl_e1psh("imsl_c_eig_gen");
#endif
    lv_eval = NULL;
    IMSL_CALL(l_eig_gen(n, a, argptr));
    va_end(argptr);
#ifdef DOUBLE
    imsl_e1pop("imsl_z_eig_gen");
#else
    imsl_e1pop("imsl_c_eig_gen");
#endif
    return lv_eval;
}


#ifdef ANSI
static VA_LIST_HACK l_eig_gen(Mint n, Mf_complex *a, va_list argptr)
#else
static VA_LIST_HACK l_eig_gen(n, a, argptr)
    Mint	n;
    Mf_complex	*a;
    va_list	argptr;
#endif
{
    Mint	    code;
    Mint	    arg_number  = 2;
    Mint	    a_col_dim   = n;
    Mf_complex	    **evec	= NULL;
    Mf_complex	    *evecu	= NULL;
    Mint	    evecu_col_dim = n;
    Mint	    vectors	= 0;
    Mint	    vectors_user = 0;
    Mint	    evals_user	= 0;
    Mf_complex	    *cwork	= NULL;
    Mf_complex	    *acopy	= NULL;
    Mfloat	    *work	= NULL;
    Mint	    *iwork	= NULL;

    code = 1;
    while (code > 0) {
	code = va_arg(argptr, Mint);
	arg_number++;
	switch (code) {
	    case IMSL_VECTORS:
		evec = va_arg(argptr, Mf_complex**);
		arg_number++;
		vectors = 1;
		break;
	    case IMSL_VECTORS_USER:
		evecu = va_arg(argptr, Mf_complex*);
		arg_number++;
		vectors_user = 1;
		break;
	    case IMSL_RETURN_USER:
		lv_eval = va_arg(argptr, Mf_complex*);
		arg_number++;
		evals_user = 1;
		break;
	    case IMSL_A_COL_DIM:
		a_col_dim = va_arg(argptr, Mint);
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
    }
    if (imsl_n1rty(0)) goto RETURN;

    iwork   = (Mint *) imsl_malloc(n*sizeof(*iwork));
    work    = (Mfloat *) imsl_malloc(n*sizeof(*work));
    cwork	= (Mf_complex *) imsl_malloc(2*n*sizeof(*cwork));
    if (vectors)
	*evec	= (Mf_complex *) imsl_malloc(n*n*sizeof(**evec));
    acopy = (Mf_complex *) imsl_malloc(n*n*sizeof(*acopy));

    if (iwork == NULL || work==NULL || cwork==NULL || acopy==NULL || 
		(vectors && (*evec==NULL)) ) {
	imsl_e1stl(1, "n");
	imsl_e1sti(1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto FREE_SPACE;
    }
    if (lv_eval == NULL) {
	lv_eval = (Mf_complex *)imsl_malloc (n*sizeof(*lv_eval));
	if (lv_eval == NULL) {
	    imsl_e1stl(1, "n");
	    imsl_e1sti(1, n);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE;
	}
    }

    if (!vectors && !vectors_user) {
	l_e3lcg(&n, a, &a_col_dim, lv_eval, acopy, work, cwork, iwork);
	}

    if (vectors) { 
	imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
        l_e6ccg(&n, a, &a_col_dim, lv_eval, *evec, &n, acopy,
			 work, cwork, iwork);
	imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
	imsl_trncr(n, n, *evec, n, n, n, *evec, n); 
	}
    if (vectors_user) {
	imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
        l_e6ccg(&n, a, &a_col_dim, lv_eval, evecu, &evecu_col_dim,
			 acopy, work, cwork, iwork);
	imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
	imsl_trncr(n, n, evecu, evecu_col_dim, n, n, evecu,
			 evecu_col_dim); 
	}


FREE_SPACE:
    if (iwork != NULL) imsl_free(iwork);
    if (work  != NULL) imsl_free(work);
    if (cwork != NULL) imsl_free(cwork);
    if (acopy != NULL) imsl_free(acopy);
RETURN:
    if (imsl_n1rty(0)>3) {
       if (!evals_user) {
           if (lv_eval != NULL) imsl_free(lv_eval);
       }
       lv_eval = NULL;  
    }
    return (argptr);
}





/* Structured by FOR_STRUCT, v0.2, on 11/08/90 at 13:35:32
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E3CCG/DE3CCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 19, 1985

    Purpose:    Balance a d_complex general matrix.

    Usage:      CALL E3CCG (N, A, LDA, LOW, IGH, SCALE)

    Arguments:
       N      - Order of the matrix.  (Input)
       A      - Complex general N by N matrix.  (Input/Output)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       LOW    - The lower boundary of the balanced matrix.  (Output)
       IGH    - The upper boundary of the balanced matrix.  (Output)
       SCALE  - Real vector of length N containing information about
                the similarity transformation.  (Output)

    Remark:
       E3CCG is based on the EISPACK routine CBAL.

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e3ccg(Mint *n, Mf_complex *a, Mint *lda, Mint *low,
			Mint *igh, Mfloat scale[])
#else
static void l_e3ccg (n, a, lda, low, igh, scale)
    Mint        *n;
    Mf_complex  *a;
    Mint        *lda, *low, *igh;
    Mfloat       scale[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
    Mint   noconv;
    Mint         _l0, _l1, i, iexc, j, k, l, m;
    Mfloat       b2, c, f, g, r, radix, s, temp;
    Mf_complex   zdum;


#define CABS1(zdum)	(Mfloat)(fabs( imsl_fc_convert( (zdum) ) ) + fabs( imsl_c_aimag( (zdum) ) ))

    imsl_e1psh ("E3CCG ");

    radix = nint (pow (10.0e0, imsl_amach (5)));

    b2 = radix * radix;
    k = 1;
    l = *n;
    goto L_40;

L_10:
    scale[m - 1] = j;
    if (j != m) {
	imsl_cswap (&l, A (j - 1, 0), ADR (_l0, 1), A (m - 1, 0), ADR (_l1, 1));
	imsl_cswap (ADR (_l0, *n - k + 1), A (k - 1, j - 1), lda, A (k - 1, m - 1),
	    lda);
    }

L_20:
    if (iexc == 2)
	goto L_70;
    /*
     * Search for rows isolating an eigenvalue and push them down
     */
L_30:
    if (l == 1)
	goto L_170;
    l -= 1;

L_40:
    for (j = l; j >= 1; j--) {
	for (i = 1; i <= l; i++) {
	    if (i == j)
		goto L_50;
	    if (imsl_fc_convert (*A (i - 1, j - 1)) != 0.0e0 || imsl_c_aimag (*A (i - 1, j - 1)) !=
		0.0e0)
		goto L_60;
    L_50:
	    ;
	}

	m = l;
	iexc = 1;
	goto L_10;
L_60:
	;
    }

    goto L_80;
    /*
     * Search for columns isolating an eigenvalue and push them left
     */
L_70:
    k += 1;

L_80:
    for (j = k; j <= l; j++) {

	for (i = k; i <= l; i++) {
	    if (i == j)
		goto L_90;
	    if (imsl_fc_convert (*A (j - 1, i - 1)) != 0.0e0 || imsl_c_aimag (*A (j - 1, i - 1)) !=
		0.0e0)
		goto L_100;
    L_90:
	    ;
	}

	m = k;
	iexc = 2;
	goto L_10;
L_100:
	;
    }
    /*
     * Now balance the submatrix in rows K to L
     */
    sset (l - k + 1, 1.0e0, &scale[k - 1], 1);
    /* Iterative loop for norm reduction */
L_110:
    noconv = 0;

    for (i = k; i <= l; i++) {
	c = imsl_scasum (ADR (_l0, l - k + 1), A (i - 1, k - 1), ADR (_l1, 1));
	temp = CABS1 (*A (i - 1, i - 1));
	c -= temp;
	r = imsl_scasum (ADR (_l0, l - k + 1), A (k - 1, i - 1), lda);
	r -= temp;
	/*
	 * Guard against zero C or R due to underflow
	 */
	if (c == 0.0e0 || r == 0.0e0)
	    goto L_160;
	g = r / radix;
	f = 1.0e0;
	s = c + r;
L_120:
	if (c >= g)
	    goto L_130;
	f *= radix;
	c *= b2;
	goto L_120;
L_130:
	g = r * radix;
L_140:
	if (c < g)
	    goto L_150;
	f /= radix;
	c /= b2;
	goto L_140;
	/* Now balance */
L_150:
	if ((c + r) / f < 0.95e0 * s) {
	    g = 1.0e0 / f;
	    scale[i - 1] *= f;
	    noconv = 1;

	    imsl_csscal (ADR (_l0, *n - k + 1), &g, A (k - 1, i - 1), lda);
	    imsl_csscal (&l, &f, A (i - 1, 0), ADR (_l0, 1));
	}
L_160:
	;
    }

    if (noconv)
	goto L_110;

L_170:
    *low = k;
    *igh = l;

    imsl_e1pop ("E3CCG ");
    return;
#undef	CABS1
}				/* end of function */



























/* Structured by FOR_STRUCT, v0.2, on 11/08/90 at 13:59:46
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E3CCH/DE3CCH (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 9, 1986

    Purpose:    Determine the eigenvalues and (optionally) the
                eigenvectors of a d_complex upper Hessenberg matrix.

    Usage:      CALL E3CCH (N, LOW, IGH, ORT, A, LDA, EVAL, VECTOR, EVEC,
               &   LDEVEC, WORK)

    Arguments:
       N      - Order of the matrix.  (Input)
       LOW    - First eigenvalue not already isolated in A.  (Input)
       IGH    - Last eigenvalue not already isolated in A.  (Input)
       ORT    - Complex vector of length N containing information about
                the reduction of the matrix to upper Hessenberg form.
                ORT is destroyed by this routine.  (Input/Output)
                If there is no such transformation then ORT should
                initially be set to zero.
       A      - Complex upper Hessenberg matrix of order N.
                (Input/Output)
                It is destroyed by this routine.
       LDA    - The leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       EVAL   - Complex vector of length N containing the eigenvalues.
                of A.  (Output)
       VECTOR - Logical variable, .TRUE. if the eigenvectors are
                computed.  (Input)
       EVEC   - Complex N by N matrix containing on input the
                transformation matrix used to reduce the matrix to upper
                Hessenberg form.  On output, it contains the
                eigenvectors.  (Input/Output)
                If there is no such transformation then EVEC should
                initially be set to the identity matrix.
       LDEVEC - The leading dimension of EVEC exactly as specified in the
                dimensions statement of the calling program.  (Input)
       WORK   - Complex work vector of length N.  (Output)

    Remarks:
       E3CCH is based on the EISPACK routine COMQR2.

    Copyright:  1986 by IMSL, Inc.  All rights reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e3cch (Mint *n, Mint *low, Mint *igh, Mf_complex ort[],
			Mf_complex *a, Mint *lda, Mf_complex eval[],
			Mint *vector, Mf_complex *evec, Mint *ldevec,
			Mf_complex work[])
#else
static void l_e3cch (n, low, igh, ort, a, lda, eval, vector, evec,
                ldevec, work)
    Mint        *n, *low, *igh;
    Mf_complex   ort[], *a;
    Mint        *lda;
    Mf_complex   eval[];
    Mint  *vector;
    Mf_complex  *evec;
    Mint        *ldevec;
    Mf_complex   work[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))
    Mint         _l0, _l1, _l2, _l3, i, ien, itn, its, j,
                k, l, m, tmp_itn;
    Mfloat       _f0, anorm, sr, temp3, tiny, tmp1, tmp2,
                tol, tr, tst1, tst2, yr;
    Mf_complex   _cx0, _cx1, s, t, temp1, temp2, tmp3, tmp4,
                x, y, zdum, zz;
    Mf_complex ctemp1;
    /* *SQRT */

#define CABS1(zdum)	(Mfloat)(fabs( imsl_fc_convert( (zdum) ) ) + fabs( imsl_c_aimag( (zdum) ) ))

    imsl_e1psh ("E3CCH ");
    /* Set constants */
    tol = imsl_amach (4);
    tiny = 100.0 * imsl_f_max (imsl_amach (1), 1.0 / imsl_amach (2));
    /* Initialize EVEC */
    if (*vector) {
	for (j = 1; j <= *n; j++) {
	    imsl_cset (n, ADR (_cx0, imsl_cf_convert (0.0, 0.0)), EVEC (j - 1, 0), ADR (_l0, 1));
	}
	imsl_cset (n, ADR (_cx0, imsl_cf_convert (1.0, 0.0)), evec, ADR (_l0, *ldevec +
		1));
    }
    /*
     * Form matrix of accumulated transformations from info in ORT
     */
    if (*vector) {
	for (i = *igh - 1; i >= (*low + 1); i--) {
	    if ((imsl_fc_convert (ort[i - 1]) != 0.0 || imsl_c_aimag (ort[i - 1]) !=
		    0.0) && (imsl_fc_convert (*A (i - 2, i - 1)) != 0.0 || imsl_c_aimag (*A (i - 2, i - 1)) !=
		    0.0)) {
		anorm = imsl_fc_convert (*A (i - 2, i - 1)) * imsl_fc_convert (ort[i - 1]) +
		    imsl_c_aimag (*A (i - 2, i - 1)) * imsl_c_aimag (ort[i - 1]);
		/*
		 * ANORM is 0 if A(I,I-1) and ORT(I) are both very small,
		 * even if not zero
		 */
		if (anorm != 0.0) {
		    imsl_ccopy (ADR (_l0, *igh - i), A (i - 2, i), ADR (_l1, 1),
			&ort[i], ADR (_l2, 1));
		    imsl_cgemv ("CONJG-TRANS", sizeof ("CONJG-TRANS"), ADR (_l0, *igh -
			    i + 1), ADR (_l1, *igh - i + 1), ADR (_cx0, imsl_cf_convert (1.0, 0.0)),
			EVEC (i - 1, i - 1), ldevec, &ort[i - 1], ADR (_l2, 1),
			ADR (_cx1, imsl_cf_convert (0.0, 0.0)), work, ADR (_l3, 1));
		    imsl_cgerc (ADR (_l0, *igh - i + 1), ADR (_l1, *igh - i +
			    1), ADR (_cx0, imsl_cf_convert (1.0 / anorm, 0.0)), &ort[i - 1],
			ADR (_l2, 1), work, ADR (_l3, 1), EVEC (i - 1, i - 1),
			ldevec);
		}
	    }
	}
    }
    /* Create real subdiagonal elements */
    for (i = *low + 1; i <= *igh; i++) {
	if (imsl_c_aimag (*A (i - 2, i - 1)) != 0.0e0) {
	    anorm = imsl_c_abs (*A (i - 2, i - 1));
	    y = imsl_c_div (*A (i - 2, i - 1), imsl_cf_convert (anorm, 0.));
	    *A (i - 2, i - 1) = imsl_cf_convert (anorm, 0.0);
	    imsl_cscal (ADR (_l0, *n - i + 1), ADR (_cx0, imsl_c_conjg (y)), A (i - 1, i - 1),
		lda);
	    imsl_cscal (ADR (_l0, imsl_i_min (i + 1, *igh)), &y, A (i - 1, 0), ADR (_l1, 1));
	    if (*vector)
		imsl_cscal (ADR (_l0, *igh - *low + 1), &y, EVEC (i - 1, *low - 1),
		    ADR (_l1, 1));
	}
    }
    /* Store roots isolated by balancing */
    for (i = 1; i <= *n; i++) {
	if (i < *low || i > *igh)
	    eval[i - 1] = *A (i - 1, i - 1);
    }

    ien = *igh;
    t = imsl_cf_convert (0.0e0, 0.);
    itn = 30 ** n;
    tmp_itn = itn;
    /* Search for next eigenvalue */
L_50:
    if (ien < *low)
	goto L_140;
    its = 0;
    /*
     * Look for single small sub-diagonal element
     */
L_60:
    for (l = ien; l >= *low; l--) {
	if (l == *low)
	    goto L_80;
	if (fabs (imsl_fc_convert (*A (l - 2, l - 1))) <= imsl_f_max (tol * (CABS1 (*A (l - 2, l - 2)) +
		    CABS1 (*A (l - 1, l - 1))), tiny))
	    goto L_80;
    }
    /* Form shift */
L_80:
    if (l == ien) {
	/* A root found */
	*A (ien - 1, ien - 1) = imsl_c_add (*A (ien - 1, ien - 1), t);
	eval[ien - 1] = *A (ien - 1, ien - 1);
	ien -= 1;
	goto L_50;
    }
    if (itn == 0) {
	imsl_e1sti(1, tmp_itn);
        imsl_ermes (IMSL_FATAL, IMSL_SLOW_CONVERGENCE_GEN);
	goto L_9000;
    }
    if (its == 10 || its == 20) {
	/* Form exceptional shift */
	s = imsl_cf_convert (fabs (imsl_fc_convert (*A (ien - 2, ien - 1))) + fabs (imsl_fc_convert (*A (ien - 3, ien - 2))), 0.);
    }
    else {
	s = *A (ien - 1, ien - 1);
	x = imsl_c_mul (*A (ien - 1, ien - 2), imsl_cf_convert (imsl_fc_convert (*A (ien - 2, ien - 1)), 0.));
	if (imsl_fc_convert (x) != 0.0e0 || imsl_c_aimag (x) != 0.0e0) {
	    y = imsl_c_div ((imsl_c_sub (*A (ien - 2, ien - 2), s)), imsl_cf_convert (2.0e0, 0.));
	    zz = imsl_c_sqrt (imsl_c_add (imsl_cf_power (y, (Mdouble) 2), x));
	    temp1 = imsl_c_mul (y, imsl_c_conjg (zz));
	    zz = imsl_c_mul (imsl_cf_convert (sign (1.0, imsl_fc_convert (temp1)), 0.), zz);
	    s = imsl_c_sub (s, imsl_c_div (x, (imsl_c_add (y, zz))));
	}
    }

    imsl_cadd (ADR (_l0, ien - *low + 1), ADR (_cx0, imsl_c_neg (s)), A (*low - 1, *low - 1),
	ADR (_l1, *lda + 1));

    t = imsl_c_add (t, s);
    its += 1;
    itn -= 1;
    /* Reduce to triangle (rows) */
    for (i = l + 1; i <= ien; i++) {
#if 0
	sr = imsl_fc_convert (*A (i - 2, i - 1));
	tmp1 = imsl_c_abs (*A (i - 2, i - 2));
	anorm = hypot (tmp1, sr);
	x = imsl_c_div (*A (i - 2, i - 2), imsl_cf_convert (anorm, 0.));
	eval[i - 2] = x;
	*A (i - 2, i - 2) = imsl_cf_convert (anorm, 0.0);
	tmp2 = sr / anorm;
	*A (i - 2, i - 1) = imsl_cf_convert (0.0e0, tmp2);
#endif
	sr = (A (i - 2, i - 1))->re;
	tmp1 = hypot((A(i-2, i-2))->re, (A (i-2, i-2))->im);
        anorm = hypot (tmp1, sr);
        x = imsl_c_div (*A (i - 2, i - 2), imsl_cf_convert (anorm, 0.));
        eval[i - 2] = x;  
	(A(i-2, i-2))->re = anorm;
	(A(i-2, i-2))->im = 0.0;
	tmp2 = sr / anorm;
	(A(i-2, i-1))->re = 0.0;
	(A(i-2, i-1))->im = tmp2;

	for (j = i; j <= *n; j++) {
	    y = *A (j - 1, i - 2);
	    zz = *A (j - 1, i - 1);
#if 0
	    tmp3 = imsl_c_conjg (x);
	    tmp3 = imsl_c_mul (tmp3, y);
	    tmp1 = imsl_c_aimag (*A (i - 2, i - 1));
	    tmp4 = imsl_c_mul (imsl_cf_convert (tmp1, 0.), zz);
	    *A (j - 1, i - 2) = imsl_c_add (tmp3, tmp4);
	    tmp3 = imsl_c_mul (x, zz);
	    tmp4 = imsl_c_mul (imsl_cf_convert (tmp1, 0.), y);
	    tmp4 = imsl_c_neg (tmp4);
	    *A (j - 1, i - 1) = imsl_c_add (tmp4, tmp3);
#endif
	    tmp3.re = x.re*y.re + x.im*y.im;
	    tmp3.im = x.re*y.im - x.im*y.re;
	    tmp1 = (A(i-2, i-1))->im;
	    tmp4.re = tmp1*zz.re;
	    tmp4.im = tmp1*zz.im;
	    (A(j-1, i-2))->re = tmp3.re + tmp4.re;
	    (A(j-1, i-2))->im = tmp3.im + tmp4.im;
	    tmp3.re = x.re*zz.re - x.im*zz.im;
	    tmp3.im = x.re*zz.im + x.im*zz.re;
	    tmp4.re = -tmp1*y.re;
	    tmp4.im = -tmp1*y.im;
	    (A(j-1, i-1))->re = tmp3.re + tmp4.re;
	    (A(j-1, i-1))->im = tmp3.im + tmp4.im;
	}

    }

    s = imsl_cf_convert (sr, imsl_c_aimag (*A (ien - 1, ien - 1)));
    if (imsl_c_aimag (s) != 0.0e0) {
	anorm = imsl_c_abs (*A (ien - 1, ien - 1));
	s = imsl_c_div (*A (ien - 1, ien - 1), imsl_cf_convert (anorm, 0.));
	*A (ien - 1, ien - 1) = imsl_cf_convert (anorm, 0.0);
	if (ien != *n)
	    imsl_cscal (ADR (_l0, *n - ien), ADR (_cx0, imsl_c_conjg (s)), A (ien, ien - 1),
		lda);
    }
    /* Inverse operation (columns) */
    for (j = l + 1; j <= ien; j++) {
	x = eval[j - 2];
	for (i = 1; i <= j; i++) {
	    y = *A (j - 2, i - 1);
	    zz = *A (j - 1, i - 1);
	    if (i == j) {
#if 0
		y = imsl_cf_convert (imsl_fc_convert (y), 0.);
		temp2 = imsl_c_add (imsl_c_mul (x, y), imsl_c_mul (imsl_cf_convert (imsl_c_aimag (*A (j - 2, j - 1)), 0.),
			zz));
		temp3 = imsl_c_aimag (*A (j - 2, i - 1));
		*A (j - 2, i - 1) = imsl_cf_convert (imsl_fc_convert (temp2), temp3);

#endif

		y.im = 0.0;
		ctemp1.re = x.re*y.re - x.im*y.im;
		ctemp1.im = x.re*y.im + x.im*y.re;
		temp2.re = ctemp1.re + (A(j-2, j-1))->im*zz.re;
		temp2.im = ctemp1.im + (A(j-2, j-1))->im*zz.im;
		temp3 = (A(j-2, i-1))->im;
		(A(j-2, i-1))->re = temp2.re;
		(A(j-2, i-1))->im = temp3;

	    }
	    else {

#if 0
		*A (j - 2, i - 1) = imsl_c_add (imsl_c_mul (x, y), imsl_c_mul (imsl_cf_convert (imsl_c_aimag (*A (j - 2, j - 1)), 0.),
			zz));
#endif

		ctemp1.re = x.re*y.re - x.im*y.im;
		ctemp1.im = x.re*y.im + x.im*y.re;
		(A(j-2, i-1))->re = ctemp1.re + (A(j-2, j-1))->im*zz.re;
		(A(j-2, i-1))->im = ctemp1.im + (A(j-2, j-1))->im*zz.im;
	    }
#if 0
	    *A (j - 1, i - 1) = imsl_c_sub (imsl_c_mul (imsl_c_conjg (x), zz), imsl_c_mul (imsl_cf_convert (imsl_c_aimag (*A (j - 2, j - 1)), 0.),
		    y));
#endif

		ctemp1.re = x.re*zz.re + x.im*zz.im;
		ctemp1.im = x.re*zz.im - x.im*zz.re;
		(A(j-1, i-1))->re = ctemp1.re - (A(j-2, j-1))->im*y.re;
		(A(j-1, i-1))->im = ctemp1.im - (A(j-2, j-1))->im*y.im;
		
	}

	if (*vector) {
	    for (i = *low; i <= *igh; i++) {
		y = *EVEC (j - 2, i - 1);
		zz = *EVEC (j - 1, i - 1);

#if 0
		*EVEC (j - 2, i - 1) = imsl_c_add (imsl_c_mul (x, y), imsl_c_mul (imsl_cf_convert (imsl_c_aimag (*A (j - 2, j - 1)), 0.),
			zz));
		*EVEC (j - 1, i - 1) = imsl_c_sub (imsl_c_mul (imsl_c_conjg (x), zz), imsl_c_mul (imsl_cf_convert (imsl_c_aimag (*A (j - 2, j - 1)), 0.),
			y));
#endif

		ctemp1.re = x.re*y.re - x.im*y.im;
		ctemp1.im = x.re*y.im + x.im*y.re;
		(EVEC(j-2, i-1))->re = ctemp1.re + (A(j-2, j-1))->im * zz.re;
		(EVEC(j-2, i-1))->im = ctemp1.im + (A(j-2, j-1))->im * zz.im;

		ctemp1.re = x.re*zz.re + x.im*zz.im;
		ctemp1.im = x.re*zz.im - x.im*zz.re;
		(EVEC(j-1, i-1))->re = ctemp1.re - (A(j-2, j-1))->im * y.re;
		(EVEC(j-1, i-1))->im = ctemp1.im - (A(j-2, j-1))->im * y.im;

	    }
	}

    }

    if (imsl_c_aimag (s) != 0.0e0) {
	imsl_cscal (&ien, &s, A (ien - 1, 0), ADR (_l0, 1));
	if (*vector)
	    imsl_cscal (ADR (_l0, *igh - *low + 1), &s, EVEC (ien - 1, *low - 1),
		ADR (_l1, 1));
    }
    goto L_60;
    /* All roots found */
L_140:
    ;
    if (!*vector)
	goto L_9000;
    /*
     * Backsubstitute to find vectors of upper triangular form
     */
    anorm = 0.0e0;
    for (i = 1; i <= *n; i++) {
	k = imsl_icamax (ADR (_l0, *n - i + 1), A (i - 1, i - 1), lda);
	j = (k - 1) / *lda + i;
	anorm = imsl_f_max (anorm, CABS1 (*A (j - 1, i - 1)));
    }

    if (*n == 1 || anorm == 0.0e0)
	goto L_9000;

    for (ien = *n; ien >= 2; ien--) {
	x = eval[ien - 1];
	*A (ien - 1, ien - 1) = imsl_cf_convert (1.0e0, 0.);

	for (i = ien - 1; i >= 1; i--) {
	    zz = imsl_cdotu (ADR (_l0, ien - i), A (i, i - 1), lda, A (ien - 1, i),
		ADR (_l1, 1));
	    y = imsl_c_sub (x, eval[i - 1]);
	    if (CABS1 (zz) <= tiny) {
		*A (ien - 1, i - 1) = imsl_cf_convert (0.0, 0.);
	    }
	    else if (CABS1 (y) <= tiny * imsl_f_max (CABS1 (zz), 1.0e0)) {
		tst1 = anorm;
		yr = tst1;
	L_160:
		yr *= 0.01e0;
		tst2 = anorm + yr;
		if (tst2 > tst1)
		    goto L_160;
		*A (ien - 1, i - 1) = imsl_c_div (zz, imsl_cf_convert (yr, 0.));
	    }
	    else {
		*A (ien - 1, i - 1) = imsl_c_div (zz, y);
	    }
	    /* Overflow control */
	    tr = CABS1 (*A (ien - 1, i - 1));
	    if (tr > 1.0 / sqrt (tol)) {
		imsl_csscal (ADR (_l0, ien - i + 1), ADR (_f0, 1.0e0 / tr), A (ien - 1, i - 1),
		    ADR (_l1, 1));
	    }
	}

    }
    /*
     * End backsubstitution Vectors of isolated roots
     */
    for (i = 1; i <= (*n - 1); i++) {
	if (i < *low || i > *igh) {
	    imsl_ccopy (ADR (_l0, *n - i), A (i, i - 1), lda, EVEC (i, i - 1),
		ldevec);
	}
    }
    /*
     * Multiply by transformation matrix to give vectors of original full
     * matrix
     */
    for (j = *n; j >= (*low + 1); j--) {
	m = imsl_i_min (j, *igh);
	imsl_cgemv ("NOT-TRANS", sizeof ("NOT-TRANS"), ADR (_l0, *igh - *low +
		1), ADR (_l1, m - *low + 1), ADR (_cx0, imsl_cf_convert (1.0, 0.0)), EVEC (*low - 1, *low - 1),
	    ldevec, A (j - 1, *low - 1), ADR (_l2, 1), ADR (_cx1, imsl_cf_convert (0.0, 0.0)),
	    work, ADR (_l3, 1));
	imsl_ccopy (ADR (_l0, *igh - *low + 1), work, ADR (_l1, 1), EVEC (j - 1, *low - 1),
	    ADR (_l2, 1));
    }

L_9000:
    ;
    imsl_e1pop ("E3CCH ");
    return;
#undef	CABS1
}				/* end of function */



























/* Structured by FOR_STRUCT, v0.2, on 11/08/90 at 13:33:40
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E3LCG/DE3LCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 11, 1990

    Purpose:    Compute all of the eigenvalues of a d_complex matrix.

    Usage:      CALL E3LCG (N, A, LDA, EVAL, ACOPY, RWK, CWK, IWK)

    Arguments:  (See EVLCG)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e3lcg (Mint *n, Mf_complex *a, Mint *lda, Mf_complex eval[],
			Mf_complex *acopy, Mfloat rwk[], Mf_complex *cwk,
			Mint iwk[])
#else
static void l_e3lcg (n, a, lda, eval, acopy, rwk, cwk, iwk)
    Mint        *n;
    Mf_complex  *a;
    Mint        *lda;
    Mf_complex   eval[], *acopy;
    Mfloat       rwk[];
    Mf_complex  *cwk;
    Mint         iwk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define ACOPY(I_,J_)	(acopy+(I_)*(*n)+(J_))
#define CWK(I_,J_)	(cwk+(I_)*(*n)+(J_))
    Mint         _l0, _l1, _l2, i, igh, j, k, low;
    Mf_complex   *evec = NULL;


    imsl_e1psh ("E3LCG ");
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
	goto L_9000;
    }
    /* Copy A to ACOPY */
    l_ccgcg (n, a, lda, acopy, n);
    /* Balance the matrix */
    l_e3ccg (n, acopy, n, &low, &igh, rwk);
    /*
     * Reduce to upper Hessenberg form by elementary similarity
     * transformations
     */
    l_e4ccg (n, &low, &igh, acopy, n, cwk, CWK (1, 0));
    /* Find eigenvalues and eigenvectors */
	_l0 = 0;
    l_e3cch (n, &low, &igh, cwk, acopy, n, eval, &_l0, evec,
	ADR (_l1, 1), CWK (1, 0));
    if (imsl_n1rty (0) != 0)
	goto L_9000;
    /*
     * Sort the eigenvalue magnitudes. Ultimately want eigenvalues in
     * descending magnitude.
     */
    for (i = 1; i <= *n; i++) {
	iwk[i - 1] = i;
	rwk[i - 1] = -imsl_c_abs (eval[i - 1]);
    }


    imsl_svrgp (*n, &rwk[0], &rwk[0], iwk);
    /* Move the eigenvalues */
    for (i = 1; i <= *n; i++) {
	for (j = i; j <= *n; j++) {
	    if (iwk[j - 1] == i) {
		k = iwk[i - 1];
		iwk[i - 1] = j;
		iwk[j - 1] = k;
		goto L_30;
	    }
	}
L_30:
	;
    }
    /* Take product of cycles and permute */
    for (i = *n - 1; i >= 1; i--) {
	imsl_cswap (ADR (_l0, 1), &eval[i - 1], ADR (_l1, 1), &eval[iwk[i - 1] -
		1], ADR (_l2, 1));
    }

L_9000:
    imsl_e1pop ("E3LCG ");
    return;
}				/* end of function */


#undef A
#undef ACOPY
#undef CWK


















/* Structured by FOR_STRUCT, v0.2, on 11/08/90 at 13:37:31
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E4CCG/DE4CCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 9, 1986

    Purpose:    Reduce a d_complex general matrix to an upper Hessenberg
                matrix by orthogonal similarity transformations.

    Usage:      CALL E4CCG (N, LOW, IGH, A, LDA, ORT, WORK)

    Arguments:
       N      - Order of the matrix.  (Input)
       LOW    - Lower boundary of the balanced matrix.  (Input)
       IGH    - Upper boundary of the balanced matrix.  (Input)
                If the matrix is not balanced set LOW=1 and IGH=N.
       A      - Complex N by N matrix to be reduced.  (Input/Output)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       ORT    - Complex vector of length IGH containing information about
                the transformation.  (Output)
       WORK   - Complex work vector of length N.  (Output)

    Remark:
       E4CCG is based on the EISPACK routine CORTH.

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e4ccg (Mint *n, Mint *low, Mint *igh, Mf_complex *a,
			Mint *lda, Mf_complex ort[], Mf_complex work[])
#else
static void l_e4ccg (n, low, igh, a, lda, ort, work)
    Mint        *n, *low, *igh;
    Mf_complex  *a;
    Mint        *lda;
    Mf_complex   ort[], work[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
    Mint         _l0, _l1, _l2, _l3, i, m;
    Mfloat       f, g, h, scale;
    Mf_complex   _cx0, _cx1;


    imsl_e1psh ("E4CCG ");

    for (m = *low + 1; m <= (*igh - 1); m++) {
	ort[m - 1] = imsl_cf_convert (0., 0.);
	/* Scale column (TOL then not needed) */
	scale = imsl_scasum (ADR (_l0, *igh - m + 1), A (m - 2, m - 1), ADR (_l1, 1));
	if (scale == 0.0e0)
	    goto L_20;

	for (i = m; i <= *igh; i++) {
	    ort[i - 1] = imsl_c_div (*A (m - 2, i - 1), imsl_cf_convert (scale, 0.));
	}
	g = imsl_scnrm2 (ADR (_l0, *igh - m + 1), &ort[m - 1], ADR (_l1, 1));

	h = imsl_fi_power (g, 2);
	f = imsl_c_abs (ort[m - 1]);
	if (f == 0.0e0) {
	    ort[m - 1] = imsl_cf_convert (g, 0.);
	    *A (m - 2, m - 1) = imsl_cf_convert (scale, 0.);
	}
	else {
	    h += f * g;
	    g /= f;
	    ort[m - 1] = imsl_c_mul (imsl_cf_convert (1.0e0 + g, 0.), ort[m - 1]);
	}
	/* Form (I-(U*UT)/H) * A */
	imsl_cgemv ("CONJ-TRANS", sizeof ("CONJ-TRANS"), ADR (_l0, *igh - m +
		1), ADR (_l1, *n - m + 1), ADR (_cx0, imsl_cf_convert (1.0, 0.0)), A (m - 1, m - 1),
	    lda, &ort[m - 1], ADR (_l2, 1), ADR (_cx1, imsl_cf_convert (0.0, 0.0)), work,
	    ADR (_l3, 1));
	imsl_cgerc (ADR (_l0, *igh - m + 1), ADR (_l1, *n - m + 1), ADR (_cx0, imsl_cf_convert (-1.0 /
		    h, 0.0)), &ort[m - 1], ADR (_l2, 1), work, ADR (_l3, 1), A (m - 1, m - 1),
	    lda);
	/* Form (I-(U*UT)/H)*A*(I-(U*UT)/H) */
	imsl_cgemv ("NOT-TRANS", sizeof ("NOT-TRANS"), igh, ADR (_l0, *igh -
		m + 1), ADR (_cx0, imsl_cf_convert (1.0, 0.0)), A (m - 1, 0), lda, &ort[m - 1],
	    ADR (_l1, 1), ADR (_cx1, imsl_cf_convert (0.0, 0.0)), work, ADR (_l2, 1));
	imsl_cgerc (igh, ADR (_l0, *igh - m + 1), ADR (_cx0, imsl_cf_convert (-1.0 / h,
		    0.0)), work, ADR (_l1, 1), &ort[m - 1], ADR (_l2, 1), A (m - 1, 0),
	    lda);

	ort[m - 1] = imsl_c_mul (imsl_cf_convert (scale, 0.), ort[m - 1]);
	*A (m - 2, m - 1) = imsl_c_neg (imsl_c_mul (imsl_cf_convert (g, 0.), *A (m - 2, m - 1)));
L_20:
	;
    }

L_9000:
    imsl_e1pop ("E4CCG ");
    return;
}				/* end of function */


#undef A











/* Structured by FOR_STRUCT, v0.2, on 11/08/90 at 13:34:39
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CCGCG/DCCGCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 5, 1985

    Purpose:    Copy a d_complex general matrix.

    Usage:      CALL CCGCG (N, A, LDA, B, LDB)

    Arguments:
       N      - Order of the matrices A and B.  (Input)
       A      - Complex matrix of order N.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       B      - Complex matrix of order N containing a copy of A.
                (Output)
       LDB    - Leading dimension of B exactly as specified in the
                dimension statement of the calling program.  (Input)

    Keyword:    Basic matrix operation

    GAMS:       D1b8

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_ccgcg (Mint *n, Mf_complex *a, Mint *lda, Mf_complex *b,
			Mint *ldb)
#else
static void l_ccgcg (n, a, lda, b, ldb)
    Mint        *n;
    Mf_complex  *a;
    Mint        *lda;
    Mf_complex  *b;
    Mint        *ldb;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define B(I_,J_)	(b+(I_)*(*ldb)+(J_))
    Mint         _l0, _l1, _l2, j;


    imsl_e1psh ("CCGCG ");
    /* Check N */
    if (*n < 1) {
	imsl_e1sti (1, *n);

	imsl_ermes (IMSL_TERMINAL, IMSL_LARGER_N_VALUE_NEEDED);
	goto L_9000;
    }
    /* Check LDA */
    if (*lda < *n) {
	imsl_e1sti (1, *lda);
	imsl_e1sti (2, *n);

	imsl_ermes (IMSL_TERMINAL, IMSL_LARGER_LDA_VALUE_NEEDED);
	goto L_9000;
    }
    /* Check LDB */
    if (*ldb < *n) {
	imsl_e1sti (1, *ldb);
	imsl_e1sti (2, *n);

	imsl_ermes (IMSL_TERMINAL, IMSL_LARGER_LDB_VALUE_NEEDED);
	goto L_9000;
    }
    /* Copy */
    if (*lda == *n && *ldb == *n) {
	imsl_ccopy (ADR (_l0, *n ** n), a, ADR (_l1, 1), b, ADR (_l2, 1));
    }
    else if (*lda >= *ldb) {
	for (j = 1; j <= *n; j++) {
	    imsl_ccopy (n, A (j - 1, 0), ADR (_l0, 1), B (j - 1, 0), ADR (_l1, 1));
	}
    }
    else {
	for (j = *n; j >= 1; j--) {
	    imsl_ccopy (n, A (j - 1, 0), ADR (_l0, -1), B (j - 1, 0), ADR (_l1, -1));
	}
    }

L_9000:
    ;
    imsl_e1pop ("CCGCG ");
    return;
}				/* end of function */

#undef A
#undef B




/* Structured by FOR_STRUCT, v0.2, on 11/14/90 at 13:40:51
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E6CCG/DE6CCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 10, 1990

    Purpose:    Compute all of the eigenvalues and eigenvectors of a
                d_complex matrix. ( Deprecate October 9, 1990)

    Usage:      CALL E6CCG (N, A, LDA, EVAL, EVEC, LDEVEC, ACOPY, RWK,
                            CWK, IWK )

    Arguments:  (See EVCCG)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e6ccg (Mint *n, Mf_complex *a, Mint *lda, Mf_complex eval[],
			Mf_complex *evec, Mint *ldevec, Mf_complex *acopy,
			Mfloat rwk[], Mf_complex *cwk, Mint iwk[])
#else
static void l_e6ccg (n, a, lda, eval, evec, ldevec, acopy, rwk,
                cwk, iwk)
    Mint        *n;
    Mf_complex  *a;
    Mint        *lda;
    Mf_complex   eval[], *evec;
    Mint        *ldevec;
    Mf_complex  *acopy;
    Mfloat       rwk[];
    Mf_complex  *cwk;
    Mint         iwk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))
#define ACOPY(I_,J_)	(acopy+(I_)*(*n)+(J_))
#define CWK(I_,J_)	(cwk+(I_)*(*n)+(J_))
    Mint         _l0, _l1, _l2, i, igh, j, k, low;
    Mfloat       _f0, ri;
    Mf_complex   w;


    imsl_e1psh ("E6CCG ");
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
    /* Check LDEVEC */
    if (*ldevec < *n) {
	imsl_e1sti (1, *ldevec);
	imsl_e1sti (2, *n);

	imsl_ermes (IMSL_TERMINAL, IMSL_LDEVEC_VALUE_TOO_SMALL);
    }
    if (imsl_n1rty (0) != 0)
	goto L_9000;
    /* Copy A to ACOPY */
    l_ccgcg (n, a, lda, acopy, n);
    /* Balance the matrix */
    l_e3ccg (n, acopy, n, &low, &igh, rwk);
    /*
     * Reduce to upper Hessenberg form by elementary similarity
     * transformations
     */
    l_e4ccg (n, &low, &igh, acopy, n, cwk, CWK (1, 0));
    /* Find eigenvalues and eigenvectors */
    _l0 = 1;
    l_e3cch (n, &low, &igh, cwk, acopy, n, eval, &_l0, evec,
	ldevec, CWK (1, 0));
    if (imsl_n1rty (0) != 0)
	goto L_9000;
    /* Backtransform the eigenvectors */
    l_e5ccg (n, &low, &igh, rwk, n, evec, ldevec);
    /*
     * Sort the eigenvalue magnitudes. Ultimately want eigenvalues in
     * descending magnitude.
     */
    for (i = 1; i <= *n; i++) {
	iwk[i - 1] = i;
	rwk[i - 1] = -imsl_c_abs (eval[i - 1]);
    }

    imsl_svrgp (*n, &rwk[0], &rwk[0], iwk);
    /* Move the eigenvalues */
    for (i = 1; i <= *n; i++) {
	for (j = i; j <= *n; j++) {
	    if (iwk[j - 1] == i) {
		k = iwk[i - 1];
		iwk[i - 1] = j;
		iwk[j - 1] = k;
		goto L_30;
	    }
	}
L_30:
	;
    }
    /* Take product of cycles and permute */
    for (i = *n - 1; i >= 1; i--) {
	imsl_cswap (n, EVEC (i - 1, 0), ADR (_l0, 1), EVEC (iwk[i - 1] - 1, 0),
	    ADR (_l1, 1));
	imsl_cswap (ADR (_l0, 1), &eval[i - 1], ADR (_l1, 1), &eval[iwk[i - 1] -
		1], ADR (_l2, 1));
    }

    /*
     * Normalize eigenvectors using the Euclidean Norm
     */
    for (j = 1; j <= *n; j++) {
	ri = imsl_scnrm2 (n, EVEC (j - 1, 0), ADR (_l0, 1));
	if (ri > 0.0) {
	    imsl_csscal (n, ADR (_f0, 1.0 / ri), EVEC (j - 1, 0), ADR (_l0, 1));
	}
    }
    /*
     * Normalize each eigenvector so that its biggest component is real and
     * positive. The reason for this is that the eigenvectors then form a
     * right- hand system.
     */
    for (j = 1; j <= *n; j++) {
	for (i = 1; i <= *n; i++) {


#if 0
	    rwk[i - 1] = imsl_fc_convert (imsl_c_mul (*EVEC (j - 1, i - 1), imsl_c_conjg (*EVEC (j - 1, i - 1))));
#endif


	    rwk[i-1] = (EVEC (j - 1, i - 1))->re * (EVEC (j - 1, i - 1))->re
		+ (EVEC (j - 1, i - 1))->im * (EVEC (j - 1, i - 1))->im;

	}
	i = imsl_isamax (*n, &rwk[0], 1);
	if (imsl_c_abs (*EVEC (j - 1, i - 1)) == 0.0)
	    goto L_70;
	/* w = conjugate(evec)/imsl_c_abs(evec) */
	w = imsl_c_div (imsl_c_conjg (*EVEC (j - 1, i - 1)), imsl_cf_convert (imsl_c_abs (*EVEC (j - 1, i - 1)), 0.));
	/* x = x*w */
	imsl_cscal (n, &w, EVEC (j - 1, 0), ADR (_l0, 1));
	*EVEC (j - 1, i - 1) = imsl_cf_convert (imsl_fc_convert (*EVEC (j - 1, i - 1)), 0.);
L_70:
	;
    }

L_9000:
    imsl_e1pop ("E6CCG ");
    return;
}				/* end of function */

#undef A
#undef EVEC
#undef ACOPY
#undef CWK







/* Structured by FOR_STRUCT, v0.2, on 11/14/90 at 13:54:01
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E5CCG/DE5CCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 19, 1985

    Purpose:    Back transform eigenvectors of a matrix transformed
                by E3CCG.

    Usage:      CALL E5CCG (N, LOW, IGH, SCALE, NVEC, EVEC, LDEVEC)

    Arguments:
       N      - Order of the matrix.  (Input)
       LOW    - The lower boundary of the balanced matrix.  (Input)
       IGH    - The upper boundary of the balanced matrix.  (Input)
       SCALE  - Real vector of length N containing information about the
                transformation (usually computed by E3CCG).  (Input)
       NVEC   - Number of columns in EVEC to be back transformed.
                (Input)
       EVEC   - Complex N by N matrix containing the eigenvectors to
                be back transformed.  (Input/Output)
       LDEVEC - Leading dimension of EVEC exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remark:
       E5CCG is based on the EISPACK routine CBABK2.

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e5ccg (Mint *n, Mint *low, Mint *igh, Mfloat scale[], 
			Mint *nvec, Mf_complex *evec, Mint *ldevec)
#else
static void l_e5ccg (n, low, igh, scale, nvec, evec, ldevec)
    Mint        *n, *low, *igh;
    Mfloat       scale[];
    Mint        *nvec;
    Mf_complex  *evec;
    Mint        *ldevec;
#endif
{
#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))
    Mint         i, k;


    imsl_e1psh ("E5CCG ");

    if (*nvec == 0)
	goto L_9000;

    for (i = *low; i <= *igh; i++) {
	imsl_csscal (nvec, &scale[i - 1], EVEC (0, i - 1), ldevec);
    }

    for (i = *low - 1; i >= 1; i--) {
	k = scale[i - 1];
	if (k != i) {
	    imsl_cswap (nvec, EVEC (0, i - 1), ldevec, EVEC (0, k - 1), ldevec);
	}
    }

    for (i = *igh + 1; i <= *n; i++) {
	k = scale[i - 1];
	if (k != i) {
	    imsl_cswap (nvec, EVEC (0, i - 1), ldevec, EVEC (0, k - 1), ldevec);
	}
    }

L_9000:
    imsl_e1pop ("E5CCG ");
    return;
}				/* end of function */

#undef EVEC
