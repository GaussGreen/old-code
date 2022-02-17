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
static VA_LIST_HACK l_eig_herm(Mint n, Mf_complex *a, va_list argptr);
static void l_e3chf(Mint *n, Mf_complex *a, Mint *lda, Mfloat *d, Mfloat *e,
                    Mf_complex *tau, Mf_complex *w);
static void l_e4chf(Mint *n, Mint *nevec, Mf_complex *a, Mint *lda,
                    Mf_complex *tau, Mf_complex *evec, Mint *ldevec,
                    Mf_complex *work);
static void l_e4csf(Mint *n, Mfloat d[], Mfloat e[], Mint *vector,
                    Mfloat *evec, Mint *ldevec);
static void l_chfcg(Mint *n, Mf_complex *a, Mint *lda);
static void l_crgcg(Mint *n, Mfloat *a, Mint *lda, Mf_complex *b, Mint
*ldb);
static void l_chemv(Mchar *uplo, unsigned uplo_s, Mint *n,
                    Mf_complex *alpha, Mf_complex *a, Mint *lda,
                    Mf_complex *x, Mint *incx, Mf_complex *imsl_beta,
                    Mf_complex *y, Mint *incy);
static void l_cher(Mchar *uplo, unsigned uplo_s, Mint *n, Mfloat *alpha,
                   Mf_complex *x, Mint *incx, Mf_complex *a, Mint *lda);
static void l_cher2(Mchar *uplo, unsigned uplo_s, Mint *n,
                    Mf_complex *alpha, Mf_complex *x, Mint *incx,
                    Mf_complex *y, Mint *incy, Mf_complex *a, Mint *lda);
static void l_e3bsf(Mint *n, Mint *mxeval, Mfloat *elow, Mfloat *ehigh,
			Mint *neval, Mfloat eval[], Mfloat d[],
			Mfloat e[], Mfloat e2[], Mfloat wk1[],
			Mfloat wk2[], Mint ind[]);
static void l_crrcr(Mint *nra, Mint *nca, Mfloat *a, Mint *lda,
			Mint *nrb, Mint *ncb, Mf_complex *b,
			Mint *ldb);
static void l_e3esf(Mint *n, Mint *nevec, Mfloat eval[], Mfloat *evec,
			Mint *ldevec, Mfloat d[], Mfloat e[], Mfloat e2[],
			Mfloat wk1[], Mfloat wk2[], Mfloat wk3[], Mfloat wk4[],
			Mfloat wk5[], Mint ind[]);
static void l_e3bhf(Mint *n, Mint *mxeval, Mf_complex *a, Mint *lda,
                Mfloat *elow, Mfloat *ehigh, Mint *neval, Mfloat eval[],
                Mf_complex *acopy, Mfloat *rwk, Mf_complex *cwk,
                Mint iwk[]);
static void l_e3lhf(Mint *n, Mf_complex *a, Mint *lda, Mfloat eval[],
                Mf_complex *acopy, Mfloat rwk[], Mf_complex cwk[],
                Mint iwk[]);
static void l_e5chf(Mint *n, Mf_complex *a, Mint *lda, Mfloat eval[],
                Mf_complex *evec, Mint *ldevec, Mf_complex *acopy,
                Mfloat rwk[], Mf_complex cwk[], Mint iwk[]);
static void l_e3fhf(Mint *n, Mint *mxeval, Mf_complex *a, Mint *lda,
                Mfloat *elow, Mfloat *ehigh, Mint *neval, Mfloat eval[],
                Mf_complex *evec, Mint *ldevec, Mf_complex *acopy,
                Mfloat *ecopy, Mfloat *rwk, Mf_complex *cwk, Mint iwk[]);
#else
static VA_LIST_HACK	l_eig_herm();
static void l_e3chf();
static void l_e3bhf();
static void l_e5chf();
static void l_e3fhf();
static void l_e4chf();
static void l_e4csf();
static void l_chfcg();
static void l_crgcg();
static void l_chemv();
static void l_cher();
static void l_cher2();
static void l_e3lhf();
static void l_e3bsf();
static void l_crrcr();
static void l_e3esf();
#endif

static Mfloat	*lv_eval;

#ifdef ANSI
Mfloat *imsl_c_eig_herm(Mint n, Mf_complex *a, ...)
#else
Mfloat *imsl_c_eig_herm(n, a, va_alist)
    Mint	n;
    Mf_complex	*a;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr,a);

    E1PSH("imsl_c_eig_herm", "imsl_z_eig_herm");
    lv_eval = NULL;
    IMSL_CALL(l_eig_herm(n, a, argptr));
    va_end(argptr);
    E1POP("imsl_c_eig_herm", "imsl_z_eig_herm"); 
    return lv_eval;
}


#ifdef ANSI
static VA_LIST_HACK l_eig_herm(Mint n, Mf_complex *a, va_list argptr)
#else
static VA_LIST_HACK l_eig_herm(n, a, argptr)
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
    Mf_complex	    *non_square	= NULL;
    Mdouble         elow        = -1000000.0;
    Mdouble         ehigh       =  1000000.0;
    Mfloat          elow_float;
    Mfloat          ehigh_float;
    Mint	    mxeval	= n;
    Mint	    evecu_col_dim = n;
    Mint            user_range  = 0;
    Mint	    vectors	= 0;
    Mint	    vectors_user = 0;
    Mint	    evals_user	= 0;
    Mint            user_number = 0;
    Mint 	    *neval	= NULL;
    Mint            neval2      = 0;
    Mf_complex	    *cwork	= NULL;
    Mf_complex	    *acopy	= NULL;
    Mfloat	    *ecopy	= NULL;
    Mfloat	    *work	= NULL;
    Mint	    *iwork	= NULL;
    Mint	    *iwk	= NULL;
    Mint	    count	= 0;
    Mint 	    i		= 0;
    Mint 	    j		= 0;
    Mint            ln, l1;

    code = 1;
    while (code > 0) {
	code = va_arg(argptr, int);
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
		lv_eval = va_arg(argptr, Mfloat*);
		arg_number++;
		evals_user = 1;
		break;
            case IMSL_RANGE:
                elow = va_arg(argptr, Mdouble);
                ehigh = va_arg(argptr, Mdouble);
                elow_float = (Mfloat) elow;
                ehigh_float = (Mfloat) ehigh;
                user_range = 1;
                arg_number += 2;
                break;
	    case IMSL_RANGE_ADR:
		elow_float  = *(va_arg(argptr, Mfloat *));
		ehigh_float = *(va_arg(argptr, Mfloat *));
		elow        =  elow_float;
		ehigh       =  ehigh_float;
		user_range = 1;
		arg_number += 2;
		break;
	    case IMSL_A_COL_DIM:
		a_col_dim = va_arg(argptr, Mint);
		arg_number++;
		break;
	    case IMSL_EVECU_COL_DIM:
		evecu_col_dim = va_arg(argptr, Mint);
		arg_number++;
		break;
            case IMSL_RETURN_NUMBER:
                neval = va_arg(argptr, Mint *);
                user_number = 1;
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

    if (user_range && vectors_user) {
        imsl_ermes(IMSL_TERMINAL, IMSL_USE_IMSL_VECTORS_OPTION);
    }

    if (a==NULL) {
        imsl_e1stl (1, "a");
        imsl_ermes(IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
    }
    if (vectors_user && (evecu==NULL)) {
        imsl_e1stl (1, "evecu");
        imsl_e1stl (2, "IMSL_VECTORS_USER");
        imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
    }

    if (n <= 0) {
	imsl_e1sti(1, n);
	imsl_ermes(IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);
    } else {
        if (n > a_col_dim) {
	  imsl_e1sti(1, n);
	  imsl_e1sti(2, a_col_dim);
	  imsl_e1stl(1, "a");
	  imsl_ermes(IMSL_TERMINAL, IMSL_COL_DIM_LESS_ORDER);
	}
        if (n > evecu_col_dim) {
          imsl_e1sti(1, n);
          imsl_e1sti(2, evecu_col_dim);
          imsl_e1stl(1, "evecu");
          imsl_ermes(IMSL_TERMINAL, IMSL_COL_DIM_LESS_ORDER);
        }
    }
    if (imsl_n1rty(0)) goto RETURN;

    if (user_range) {
	if (!vectors && !vectors_user) {
		work	= (Mfloat *) imsl_malloc(5*n*sizeof(*work));
		iwk 	= (Mint *) imsl_malloc(mxeval*sizeof(*iwk));
		}
	else {
		work 	= (Mfloat *) imsl_malloc(8*n*sizeof(*work));
		iwk 	= (Mint *) imsl_malloc(mxeval*sizeof(*iwk));
		ecopy 	= (Mfloat *) imsl_malloc(n*mxeval*sizeof(*ecopy));
		}
	}
    else {
    	if (vectors || vectors_user)
    		work    = (Mfloat *) imsl_malloc((n*n+n)*sizeof(*work));
    	else
		work 	= (Mfloat *) imsl_malloc(n*sizeof(*work));
	}

    iwork	= (Mint *) imsl_malloc(n*sizeof(*iwork));
    cwork	= (Mf_complex *) imsl_malloc(2*n*sizeof(*cwork));
    acopy = (Mf_complex *) imsl_malloc(n*n*sizeof(*acopy));

    if (vectors)
	*evec	= (Mf_complex *) imsl_malloc(n*n*sizeof(**evec));

    if (work==NULL || cwork==NULL || acopy==NULL || 
		(vectors && (*evec==NULL)) || iwork==NULL ||
		(user_range && !vectors && !vectors_user && (iwk==NULL)) ) {
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

    if (!user_range && !vectors && !vectors_user) {
	l_e3lhf(&n, a, &a_col_dim, lv_eval, acopy, work, cwork, iwork);
	}

    if (!user_range && vectors) { 
	imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
        l_e5chf(&n, a, &a_col_dim, lv_eval, *evec, &n, acopy,
			 work, cwork, iwork);
	imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
        imsl_trncr(n, n, *evec, n, n, n, *evec, n);
    }
    if (!user_range && vectors_user) {
	imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
        l_e5chf(&n, a, &a_col_dim, lv_eval, evecu, &evecu_col_dim,
			 acopy, work, cwork, iwork);
	imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
	imsl_trncr(n, n, evecu, evecu_col_dim, n, n, evecu,
			 evecu_col_dim); 
	}

    if (user_range && !vectors && !vectors_user) {
        l_e3bhf(&n, &mxeval, a, &a_col_dim, &elow_float, &ehigh_float,
                         &neval2, lv_eval, acopy, work, cwork, iwk);
        /*&neval2 is needed because sending a (possibly) nil pointer, */
        /*neval, can cause an error.                                  */ 
        if (user_number) *neval = neval2;
	}

    if (user_range && vectors) {
        imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
        l_e3fhf(&n, &mxeval, a, &a_col_dim, &elow_float, &ehigh_float,
                         &neval2, lv_eval, *evec, &n, acopy, ecopy, work,
			 cwork, iwk);
        /*&neval2 is needed because sending a (possibly) nil pointer, */
        /*neval, can cause an error.                                  */
        if (user_number)  *neval = neval2;
        imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
	if (neval2 > 0) {
            non_square = (Mf_complex *) imsl_malloc(n*neval2
				*sizeof(*non_square));
            if (non_square==NULL) {
                imsl_e1stl(1, "n");
                imsl_e1sti(1, n);
                imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
                goto FREE_SPACE;
            }
            count = 0;
            for (i=0; i < n; ++i){
                for (j = 0; j < neval2; ++j){
                        if(j == 0)
                                *(non_square+count) = *(*evec+i);
                        else
                                *(non_square+count) = *(*evec+i+n);
                        count++;
                }
            }
            if (*evec != NULL) imsl_free(*evec);
            *evec = (Mf_complex *) imsl_malloc(n*neval2*sizeof(**evec));
            ln = n*neval2; l1 = 1;
            imsl_ccopy(&ln, non_square, &l1, *evec, &l1);
            imsl_free(non_square);
	    }
	}

/*This option is not allowed because of difficulty with leading */
/*dimension checks by low level subroutines. 
    if (user_range && vectors_user) {
        imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
        l_e3fhf(&n, &mxeval, a, &a_col_dim, &elow_float, &ehigh_float,
                         neval, lv_eval, evecu, &n, acopy, ecopy, work,
			 cwork, iwk);
        imsl_trncr(n, n, a, a_col_dim, n, n, a, a_col_dim);
        imsl_trncr(n, n, evecu, evecu_col_dim, n, n, evecu,
                         evecu_col_dim);
	}
*/

FREE_SPACE:
    if (work != NULL) imsl_free(work);
    if (cwork != NULL) imsl_free(cwork);
    if (acopy != NULL) imsl_free(acopy);
    if (iwk != NULL) imsl_free(iwk);
    if (iwork != NULL) imsl_free(iwork);
    if (ecopy != NULL) imsl_free(ecopy);
RETURN:
    if (imsl_n1rty(0)>3) {
        if (!evals_user && (lv_eval != NULL)) imsl_free(lv_eval);
	lv_eval = NULL;
    }
    return (argptr);
}

#if 0
#ifdef ANSI
void imsl_cgerc(Mint *m, Mint *n, Mf_complex *alpha, Mf_complex x[],
		Mint *incx, Mf_complex y[], Mint *incy, Mf_complex a[],
		Mint *lda);
void imsl_cgemv(Mchar *trans, unsigned trans_s, Mint *m, Mint *n,
		Mf_complex *alpha, Mf_complex a[], Mint *lda,
		Mf_complex x[], Mint *incx, Mf_complex *imsl_beta,
		Mf_complex y[], Mint *incy);
void imsl_caxpy(Mint *n, Mf_complex *ca, Mf_complex cx[], Mint *incx, 
		Mf_complex cy[], Mint *incy);
void imsl_cset(Mint *n, Mf_complex *ca, Mf_complex cx[], Mint *incx);
Mint imsl_icamax(Mint *n, Mf_complex *cx, Mint *incx);
#else
void imsl_cgerc();
void imsl_cgemv();
void imsl_caxpy();
void imsl_cset();
Mint imsl_icamax();
#endif
#endif

/*----------------------------------------------------------------------- */

/*  IMSL Name:  E3CHF/DE3CHF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 9, 1986

    Purpose:    Reduce a Hermitian matrix in full storage mode to
                a symmetric tridiagonal matrix using unitary
                transformations.

    Usage:      CALL E3CHF (N, A, LDA, D, E, TAU, W)

    Arguments:
       N      - Order of the matrix.  (Input)
       A      - Hermitian matrix.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       D      - Vector of length N containing the diagonal of the
                tridiagonal matrix.  (Output)
       E      - Vector of length N containing the subdiagonal of the
                tridiagonal matrix in E(2:N).  E(1) = 0.  (Output)
       TAU    - Complex vector of length N containing diagonal of the
                unitary matrix T.  See reference Mueller (1966).
                (Output)
       W      - Complex vector of length N, used as workspace.  (Output)

    Remark:
       E3CHF is based on the routine 'householder hermitian' in
       Mueller, Dennis (1966), Householder's method for d_complex matrices
       and eigensystems of Hermitian matrices, Numerische Mathematik,
       (8), p72-92 and the EISPACK routine HTRIDI.

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e3chf(Mint *n, Mf_complex *a, Mint *lda, Mfloat *d, Mfloat *e,
                    Mf_complex *tau, Mf_complex *w)
#else
static void l_e3chf(n, a, lda, d, e, tau, w)
	Mint            *n;
	Mf_complex      *a;
	Mint            *lda;
	Mfloat           d[], e[];
	Mf_complex       tau[], w[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint             _l0, _l1, _l2, i, k, l;
	Mfloat           _f0, bb, delta, ratio, rho, root,
	                scale, vr, zero;
	Mf_complex       _cx0, _cx1, temp1;


	zero = F_ZERO;
	/* Perform N-2 similarity transforms */
	for (k = 2; k <= (*n - 1); k++) {
		tau[k - 1] = imsl_cf_convert(zero, F_ZERO);
		l = k - 1;
		scale = F_ZERO;
                _l0 = *n - l; _l1 = 1; _l2 = 1;
		imsl_ccopy(&_l0, A(l - 1, k - 1), &_l1, &w[k - 1],
			   &_l2);
		for (i = k; i <= *n; i++) {
			scale += fabs(imsl_fc_convert(w[i - 1])) + fabs(imsl_c_aimag(w[i - 1]));
		}
		if (scale <= zero)
			goto L_30;
                _l0 = *n - l; _l1 = 1;
                _cx0 = imsl_cf_convert(1 / scale, zero);
		imsl_cscal(&_l0, &_cx0,
			   &w[k - 1], &_l1);
                _l0 = *n - l; _l1 = 1; _l2 = 1;
		imsl_ccopy(&_l0, &w[k - 1], &_l1, A(l - 1, k - 1),
			   &_l2);
                _l0 = *n - l; _l1 = 1; _l2 = 1;
	/*	vr = imsl_fc_convert(imsl_cdotc(&_l0, &w[k - 1], &_l1,
					     &w[k - 1], &_l2)); */

	/* problem with imsl_c_mul() in imsl_cdotc() call above
	   replaced with in-line code below ... (K.C. 10-19-90) */

		vr = F_ZERO;
		for (i=0; i <= _l0-1; ++i)
			vr += w[k-1+i].re * w[k-1+i].re +
					w[k-1+i].im * w[k-1+i].im;


		root = sqrt(vr);
		if ((fabs(imsl_fc_convert(w[k - 1])) + fabs(imsl_c_aimag(w[k - 1]))) ==
		    zero) {
			delta = vr;
			tau[0] = imsl_cf_convert(-root, F_ZERO);
			w[k - 1] = imsl_cf_convert(root, F_ZERO);
			*A(l - 1, k - 1) = w[k - 1];
		} else {
			temp1 = w[k - 1];
			root *= imsl_c_abs(temp1);
			delta = vr + root;
			ratio = vr / root;
			tau[0] = imsl_c_neg(imsl_c_mul(imsl_cf_convert(ratio, F_ZERO), imsl_c_conjg(temp1)));
			w[k - 1] = imsl_c_mul(imsl_cf_convert(ratio + F_ONE, F_ZERO), temp1);
			*A(l - 1, k - 1) = w[k - 1];
		}
		/*
		 * The matrix to be used in the similarity transformation has
		 * been determined.  Transformation follows
		 */
		for (i = k; i <= *n; i++) {
			tau[i - 1] = w[i - 1];
		}

                _l0 = *n - l; _l1 = 1; _l2 = 1;
                _cx0 = imsl_cf_convert(F_ONE, F_ZERO);
                _cx1 = imsl_cf_convert(F_ZERO, F_ZERO);
		l_chemv("LOWER", sizeof("LOWER"), &_l0, &_cx0,
			   A(k - 1, k - 1), lda, &tau[k - 1], &_l1, &_cx1,
			   &w[k - 1], &_l2);
                _l0 = *n - l; _l1 = 1;
                _cx0 = imsl_cf_convert(1 / delta, zero);
		imsl_cscal(&_l0, &_cx0,
			   &tau[k - 1], &_l1);
		/* rho = u*Nv */
                _l0 = *n - l; _l1 = 1; _l2 = 1;
		temp1 = imsl_cdotc(&_l0, &w[k - 1], &_l1, &tau[k - 1],
				   &_l2);
		rho = imsl_fc_convert(temp1);
                _l0 = *n - l; _l1 = 1; _l2 = 1;
                _cx0 = imsl_cf_convert(-F_ONE, F_ZERO);
		l_cher2("LOWER", sizeof("LOWER"), &_l0, &_cx0,
			   &tau[k - 1], &_l1, &w[k - 1], &_l2, A(k - 1, k - 1),
			   lda);
                _l0 = *n - l; _l1 = 1;
                _f0 = rho * delta;
		l_cher("LOWER", sizeof("LOWER"), &_l0, &_f0, &tau[k - 1], &_l1, A(k - 1, k - 1), lda);
		tau[k - 1] = tau[0];
		e[k - 1] = scale;
L_30:
		;
	}
	/*
	 * The matrix has been reduced to tri- diagonal Hermitian form. The
	 * sub- diagonal has been temporarily stored in vector TAU. Store the
	 * diagonal of the reduced matrix in D.
	 */
	for (i = 1; i <= *n; i++) {
		d[i - 1] = imsl_fc_convert(*A(i - 1, i - 1));
	}
	/*
	 * Obtain diagonal of matrix T by performing a diagonal unitary
	 * similarity transformation
	 */
	tau[0] = imsl_cf_convert(F_ONE, F_ZERO);
	if (*n > 1)
		tau[*n - 1] = imsl_c_conjg(*A(*n - 2, *n - 1));
	e[0] = F_ZERO;
	e[*n - 1] = F_ONE;
	/*
	 * Calculate subdiagonal E of the real symmetric tridiagonal matrix
	 * and TAU, the diagonal of the diagonal unitary matrix
	 */
	for (i = 2; i <= *n; i++) {
		bb = imsl_c_abs(tau[i - 1]);
		e[i - 1] *= bb;
		*A(i - 1, i - 1) = imsl_cf_convert(imsl_fc_convert(*A(i - 1, i - 1)), bb);
		if (bb == F_ZERO) {
			tau[i - 1] = imsl_cf_convert(F_ONE, F_ZERO);
			bb = F_ONE;
		}
		tau[i - 1] = imsl_c_div(imsl_c_mul(tau[i - 2], tau[i - 1]), imsl_cf_convert(bb, F_ZERO));
	}

	return;
}				/* end of function */
#undef  A
/*----------------------------------------------------------------------- */

/*  IMSL Name:  E4CHF/DE4CHF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 9, 1986

    Purpose:    Reduce a Hermitian matrix in full storage mode to
                a symmetric tridiagonal matrix using unitary
                transformations.

    Usage:      CALL E4CHF (N, NEVEC, A, LDA, TAU, EVEC, LDEVEC, WORK)

    Arguments:
       N      - Order of the matrix.  (Input)
       NEVEC  - Number of eigenvectors to be determined.  (Input)
       A      - Complex matrix of order N.  Its lower triangle
                contains information on the Householder reduction
                from E3CHF.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       TAU    - Complex vector of length N containing diagonal of the
                unitary matrix T.  See reference Mueller (1966).
                (Input)
       EVEC   - Complex matrix of order N.  On input it contains the
                eigenvectors of the tridiagonal matrix.  On output it
                contains the eigenvectors of the original matrix.
                (Input/Output)
       LDEVEC - Leading dimension of EVEC exactly as specified in the
                dimension statement of the calling program. (Input)
       WORK   - Complex work array of length N.  (Output)

    Remark:
       E4CHF is based on the routine 'reverse' in
       Mueller, Dennis (1966), Householder's method for d_complex matrices
       and eigensystems of Hermitian matrices, Numerische Mathematik,
       (8), p72-92.

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e4chf(Mint *n, Mint *nevec, Mf_complex *a, Mint *lda, 
                    Mf_complex *tau, Mf_complex *evec, Mint *ldevec,
                    Mf_complex *work)
#else
static void l_e4chf(n, nevec, a, lda, tau, evec, ldevec, work)
	Mint            *n, *nevec;
	Mf_complex      *a;
	Mint            *lda;
	Mf_complex       tau[], *evec;
	Mint            *ldevec;
	Mf_complex       work[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))
	Mint             _l0, _l1, _l2, j, nr;
	Mfloat           delta;
	Mf_complex       _cx0, _cx1;


	for (j = 2; j <= *n; j++) {
                _cx0 = imsl_c_conjg(tau[j - 1]);
		imsl_cscal(nevec, &_cx0, EVEC(0, j - 1),
			   ldevec);
	}
	/*
	 * Recover the Householder matrices in reverse order
	 */
	for (nr = *n - 1; nr >= 2; nr--) {
		delta = imsl_c_aimag(*A(nr - 1, nr - 1)) * imsl_c_abs(*A(nr - 2, nr - 1));
		if (delta != F_ZERO) {
                        _l0 = *n - nr + 1; _l1 = 1; _l2 = 1;
                        _cx0 = imsl_cf_convert(F_ONE, F_ZERO);
                        _cx1 = imsl_cf_convert(F_ZERO, F_ZERO);
			imsl_cgemv("CONJUGATE", sizeof("CONJUGATE"), &_l0
			, nevec, &_cx0, EVEC(0, nr - 1),
				   ldevec, A(nr - 2, nr - 1), &_l1, &_cx1,
				   work, &_l2);
                        _l0 = *n - nr + 1; _l1 = 1; _l2 = 1;
                        _cx0 = imsl_cf_convert(-F_ONE / delta, F_ZERO);
			imsl_cgerc(&_l0, nevec, &_cx0, A(nr - 2, nr - 1), &_l1, work,
			&_l2,
				   EVEC(0, nr - 1), ldevec);
		}
	}
	return;
}				/* end of function */
#undef  A
#undef  EVEC
/*----------------------------------------------------------------------- */

/*  IMSL Name:  E4CSF/DE4CSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 1,1990

    Purpose:    Compute eigenvalues and (optionally) eigenvectors of a
                symmetric triadiagonal matrix.

    Usage:      CALL E4CSF (N, D, E, VECTOR, EVEC, LDEVEC)

    Arguments:
       N      - Order of the matrix.  (Input)
       D      - Real vector of length N.  On input, the diagonal elements
                of the matrix.  On output, the eigenvalues in increasing
                order.  (Input/Output)
       E      - Real vector of length N.  On input, the elemets of the
                off diagonal.  E(1) is arbitrary.  On output, E is
                destroyed.  (Input/Output)
       VECTOR - Logical variable, .TRUE. if eigenvectors are to be
                computed.  (Input)
       EVEC   - Real matrix of order N.  On input it contains the
                transformation matrix used to transform the original
                matrix to a tridiagonal matrix.  If the eigenvectors
                of a tridiagonal matrix are desired, then it contains
                the identity matrix.  On output the eigenvector
                corresponding to EVAL(J) is stored in the J-th
                column.  (Input/Output)
       LDEVEC - Leading dimension of EVEC exactly as specified in the
                dimension statement of the calling program.  (Input)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e4csf(Mint *n, Mfloat d[], Mfloat e[], Mint *vector,
                    Mfloat *evec, Mint *ldevec)
#else
static void l_e4csf(n, d, e, vector, evec, ldevec)
	Mint            *n;
	Mfloat           d[], e[];
	Mint            *vector;
	Mfloat          *evec;
	Mint            *ldevec;
#endif
{
#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))
	Mint             _l0, _l1, i, iter, j, k, l, m;
	Mfloat           _f0, b, c, f, g, p, r, s, tiny, tol;


	imsl_e1psh("l_e4csf");

	if (*n == 1)
		goto L_9000;

	scopy(*n - 1, &e[1], 1, &e[0], 1);
	e[*n - 1] = F_ZERO;

	tiny = 100.0 * imsl_amach(1);
	tol = imsl_amach(4);
	iter = 0;
	for (l = 1; l <= *n; l++) {
		/* Look for small sub-diagonal element */
L_10:
		for (m = l; m <= *n; m++) {
			if (m == *n)
				goto L_30;
			if (fabs(e[m - 1]) <= imsl_f_max(tol * (fabs(d[m - 1]) +
							 fabs(d[m])), tiny))
				goto L_30;
		}

L_30:
		p = d[l - 1];
		if (m == l)
			goto L_60;
		if (iter == 30 ** n) {

                imsl_e1sti(1, 30 ** n);
/*			imsl_ermes(4, 1, "The iteration for the eigenvalues did not converge.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_SLOW_CONVERGENCE_GEN);
			goto L_9000;
		}
		iter += 1;
		/* Form shift */
		g = (d[l] - p) / (F_TWO * e[l - 1]);
                _f0 = F_ONE;
		r = hypot(g, _f0);
		g = d[m - 1] - p + e[l - 1] / (g + sign(r, g));
		s = F_ONE;
		c = F_ONE;
		p = F_ZERO;

		for (i = m - 1; i >= l; i--) {
			f = s * e[i - 1];
			b = c * e[i - 1];
			imsl_srotg(&g, &f, &c, &s);
			e[i] = g;
			if (g == F_ZERO)
				goto L_50;
			g = d[i] - p;
			r = (d[i - 1] - g) * s + F_TWO * c * b;
			p = s * r;
			d[i] = g + p;
			g = c * r - b;
			/* Form vector */
                        _l0 = 1; _l1 = 1;
			if (*vector)
				imsl_srot(*n, EVEC(i, 0), _l0, EVEC(i - 1, 0), _l1,
					  c, s);

		}

		d[l - 1] -= p;
		e[l - 1] = g;
		e[m - 1] = F_ZERO;
		goto L_10;
		/* Recover from underflow */
L_50:
		d[i] -= p;
		e[m - 1] = F_ZERO;
		goto L_10;
L_60:
		;
	}
	/* Order eigenvalues and eigenvectors */
	for (i = 1; i <= (*n - 1); i++) {
		k = i;
		p = d[i - 1];

		for (j = i + 1; j <= *n; j++) {
			if (d[j - 1] < p) {
				k = j;
				p = d[j - 1];
			}
		}

		if (k != i) {
			d[k - 1] = d[i - 1];
			d[i - 1] = p;
			if (*vector)
				sswap(*n, EVEC(i - 1, 0), 1, EVEC(k - 1, 0), 1);
		}
L_80:
		;

	}

L_9000:
	imsl_e1pop("l_e4csf");
	return;
}				/* end of function */
#undef  EVEC
/*----------------------------------------------------------------------- */

/*  IMSL Name:  CHFCG/DCHFCG (Single/Double precision version)

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
static void l_chfcg(Mint *n, Mf_complex *a, Mint *lda)
#else
static void l_chfcg(n, a, lda)
	Mint            *n;
	Mf_complex      *a;
	Mint            *lda;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint             i, j;
	Mfloat           eps;


	imsl_e1psh("l_chfcg");
	/* Check N */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "The argument N = %(i1).  It must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);
		goto L_9000;
	}
	/* Check LDA */
	if (*lda < *n) {
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 2, "The argument LDA = %(i1).  It must be at least as large as N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);
		goto L_9000;
	}
	/* Check that A is Hermitian */
	eps = F_TEN * imsl_amach(4);
	for (i = 1; i <= *n; i++) {
		if (fabs(imsl_c_aimag(*A(i - 1, i - 1))) != F_ZERO) {
			if (fabs(imsl_c_aimag(*A(i - 1, i - 1))) > eps * fabs(imsl_fc_convert(*A(i - 1, i - 1)))) {
				imsl_e1sti(1, i);
				imsl_e1stc(1, *A(i - 1, i - 1));
/*				imsl_ermes(4, 2, "The matrix element A(%(i1),%(i1)) = %(c1).  The diagonal of a Hermitian matrix must be real.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_HERMITIAN_DIAG_REAL);
				goto L_9000;
			} else {
				imsl_e1sti(1, i);
				imsl_e1stc(1, *A(i - 1, i - 1));
                                imsl_ermes(IMSL_FATAL, IMSL_HERMITIAN_DIAG_REAL);
				*A(i - 1, i - 1) = imsl_cf_convert(imsl_fc_convert(*A(i - 1, i - 1)), F_ZERO);
			}
		}
	}

	for (j = 1; j <= (*n - 1); j++) {
		for (i = j + 1; i <= *n; i++) {
			*A(j - 1, i - 1) = imsl_c_conjg(*A(i - 1, j - 1));
		}
	}

L_9000:
	;
	imsl_e1pop("l_chfcg");
	return;
}				/* end of function */
#undef  A
/*----------------------------------------------------------------------- */

/*  IMSL Name:  CRGCG/DCRGCG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 8, 1986

    Purpose:    Copy a real general matrix to a d_complex general matrix.

    Usage:      CALL CRGCG (N, A, LDA, B, LDB)

    Arguments:
       N      - Order of the matrices A and B.  (Input)
       A      - Real matrix of order N.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       B      - Complex matrix of order N containing a copy of A.
                (Output)
       LDB    - Leading dimension of B exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remark:
       The matrices A and B may be the same.

    Keyword:    Basic matrix operation

    GAMS:       D1b9

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_crgcg(Mint *n, Mfloat *a, Mint *lda, Mf_complex *b, Mint *ldb)
#else
static void l_crgcg(n, a, lda, b, ldb)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda;
	Mf_complex      *b;
	Mint            *ldb;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define B(I_,J_)	(b+(I_)*(*ldb)+(J_))
	Mint             i, j;


	imsl_e1psh("l_crgcg");
	/* Check N */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "The order of the matrix must be at least 1 while N = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);
		goto L_9000;
	}
	/* Check LDA */
	if (*lda < *n) {
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 2, "The leading dimension of A must be at least as large as N while LDA = %(i1) and N = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDA_GE_N);
		goto L_9000;
	}
	/* Check LDB */
	if (*ldb < *n) {
		imsl_e1sti(1, *ldb);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 3, "The leading dimension of B must be at least as large as N while LDB = %(i1) and N = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDB_GE_N);
		goto L_9000;
	}
	/* Copy A to B */
	for (j = *n; j >= 1; j--) {
		for (i = *n; i >= 1; i--) {
			*B(j - 1, i - 1) = imsl_cf_convert(*A(j - 1, i - 1), F_ZERO);
		}
	}

L_9000:
	;
	imsl_e1pop("l_crgcg");
	return;
}				/* end of function */
#undef  A
#undef  B
/*----------------------------------------------------------------------- */

/*  IMSL Name:  CHEMV (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 28, 1989

    Purpose:    Perform the matrix-vector multiplication
                y = alpha*A*x + imsl_beta*y, where A is a Hermitian matrix.

    Usage:      CALL CHEMV (UPLO, N, ALPHA, A, LDA, X, INCX, BETA,
                            Y, INCY)

    Arguments:
       UPLO   - Character specifing the storage structure.
                (Input)
                   UPLO              Structure
                'U' or 'u'      Hermitian matrix A is referenced using
                                its upper triangular part.
                'L' or 'l'      Hermitian matrix A is referenced using
                                its lower triangular part.
       N      - Order of the matrix A.  (Input)
       ALPHA  - Complex scalar multiplier for the matrix-vector product.
                (Input)
       A      - Complex array of size LDA by N containing the matrix of
                order N in Hermitian form.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       X      - Complex vector of length (N-1)*IABS(INCX)+1.  (Input)
       INCX   - Displacement between elements of X.  (Input)
       BETA   - Complex scalar multiplier for Y.  (Input)
                When BETA is zero, Y is not referenced.  In that case,
                BETA*Y is defined as the zero vector.
       Y      - Complex vector of length (N-1)*IABS(INCY)+1.  (Input)
       INCY   - Displacement between elements of Y.  (Input)

    GAMS:       D1b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_chemv(Mchar *uplo, unsigned uplo_s, Mint *n,
                    Mf_complex *alpha, Mf_complex *a, Mint *lda, 
                    Mf_complex *x, Mint *incx, Mf_complex *imsl_beta,
	            Mf_complex *y, Mint *incy)
#else
static void l_chemv(uplo, uplo_s, n, alpha, a, lda, x, incx, imsl_beta,
         	   y, incy)
	Mchar           *uplo;
	unsigned        uplo_s;
	Mint            *n;
	Mf_complex      *alpha, *a;
	Mint            *lda;
	Mf_complex       x[];
	Mint            *incx;
	Mf_complex      *imsl_beta, y[];
	Mint            *incy;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
        Mint             lower, upper;
	Mint             _l0, _l1, i, ix, iy, j, ky;
	Mf_complex       _cx0, temp;


	upper = imsl_l1ame(uplo, uplo_s, "U", sizeof("U"));
	lower = imsl_l1ame(uplo, uplo_s, "L", sizeof("L"));
	/*
	 * Test the input parameters.
	 */
	if (*n < 0) {
		imsl_e1psh("l_chemv");
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "N must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
		imsl_e1pop("l_chemv");
		goto L_9000;
	} else if ((*lda < *n) || (*lda == 0)) {
		imsl_e1psh("l_chemv");
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 2, "LDA must be greater than or equal to N and greater than zero while LDA = %(i1) and N = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_LDA_VALUE_GIVEN);
		imsl_e1pop("l_chemv");
		goto L_9000;
	} else if (*incx == 0) {
		imsl_e1psh("l_chemv");
		imsl_e1sti(1, *incx);

/*		imsl_ermes(5, 3, "INCX must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
		imsl_e1pop("l_chemv");
		goto L_9000;
	} else if (*incy == 0) {
		imsl_e1psh("l_chemv");
		imsl_e1sti(1, *incy);

/*		imsl_ermes(5, 4, "INCY must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCY_EQUALS_ZERO);
		imsl_e1pop("l_chemv");
		goto L_9000;
	} else if ((!upper) && (!lower)) {
		imsl_e1psh("l_chemv");
		imsl_e1stl(1, uplo);

/*		imsl_ermes(5, 5, "UPLO must be set equal to U or L while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_UPLO_VALUE);
		imsl_e1pop("l_chemv");
		goto L_9000;
	}
	/* Quick return if possible */
	if (*n == 0 || (imsl_c_eq(*alpha, imsl_cf_convert(F_ZERO, F_ZERO)) && imsl_c_eq(*imsl_beta, imsl_cf_convert(F_ONE, F_ZERO))))
		goto L_9000;

	ix = 1;
	iy = 1;
	if (*incx < 0)
		ix = (-*n + 1) ** incx + 1;
	if (*incy < 0)
		iy = (-*n + 1) ** incy + 1;

	if (imsl_c_eq(*imsl_beta, imsl_cf_convert(F_ZERO, F_ZERO))) {
                _cx0 = imsl_cf_convert(F_ZERO, F_ZERO);
                _l0 = abs(*incy);
		imsl_cset(n, &_cx0, y, &_l0);
	} else if (!imsl_c_eq(*imsl_beta, imsl_cf_convert(F_ONE, F_ZERO))) {
                _l0 = abs(*incy);
		imsl_cscal(n, imsl_beta, y, &_l0);
	}
	if (imsl_c_eq(*alpha, imsl_cf_convert(F_ZERO, F_ZERO)))
		goto L_9000;

	if (upper) {
		for (j = 1; j <= *n; j++) {
			temp = imsl_c_mul(*alpha, x[ix - 1]);
			ky = iy + (j - 2) * imsl_i_min(*incy, 0);
                        _l0 = j - 1; _l1 = 1;
			imsl_caxpy(&_l0, &temp, A(j - 1, 0), &_l1,
				   &y[ky - 1], incy);
			ky = iy + (j - 1) ** incy;
			y[ky - 1] = imsl_c_add(y[ky - 1], imsl_c_mul(temp, imsl_cf_convert(imsl_fc_convert(*A(j - 1, j - 1)), F_ZERO)));
			for (i = j + 1; i <= *n; i++) {
				ky += *incy;
				y[ky - 1] = imsl_c_add(y[ky - 1], imsl_c_mul(temp, imsl_c_conjg(*A(i - 1, j - 1))));
			}
			ix += *incx;
		}
	} else {
		for (j = 1; j <= *n; j++) {
			temp = imsl_c_mul(*alpha, x[ix - 1]);
			ky = iy;
			for (i = 1; i <= (j - 1); i++) {
				y[ky - 1] = imsl_c_add(y[ky - 1], imsl_c_mul(temp, imsl_c_conjg(*A(i - 1, j - 1))));
				ky += *incy;
			}
			y[ky - 1] = imsl_c_add(y[ky - 1], imsl_c_mul(temp, imsl_cf_convert(imsl_fc_convert(*A(j - 1, j - 1)), F_ZERO)));
			ky += *incy + (*n - j - 1) * imsl_i_min(*incy, 0);
                        _l0 = *n - j;  _l1 = 1;
 			imsl_caxpy(&_l0, &temp, A(j - 1, j), &_l1,
				   &y[ky - 1], incy);
			ix += *incx;
		}
	}

L_9000:
	return;
}				/* end of function */
#undef  A
/*----------------------------------------------------------------------- */

/*  IMSL Name:  CHER  (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 12, 1989

    Purpose:    Perform the rank-one matrix update A = A +
                alpha*x*ctrans(x) to the Hermitian matrix A, with A and x
                d_complex and alpha real.  The vector ctrans(x) is the
                conjugate transpose of x.

    Usage:      CALL CHER (UPLO, N, ALPHA, X, INCX, A, LDA)

    Arguments:
       UPLO   - Character string indicating the matrix storage.  (Input)
                If UPLO is 'U' of 'u' then only the upper half of A is
                used.  If UPLO is 'L' of 'l' then only the lower half of
                A is used.
       N      - Order of the matrix A.  (Input)
       ALPHA  - Scalar multiplier for the matrix-vector product.
                (Input)
       X      - Complex vector of length (N-1)*IABS(INCX)+1.  (Input)
       INCX   - Displacement between elements of X.  (Input)
       A      - Complex array of order N.  (Input/Output)
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
static void l_cher(Mchar *uplo, unsigned uplo_s, Mint *n, Mfloat *alpha,
                   Mf_complex *x, Mint *incx, Mf_complex *a, Mint *lda)
#else
static void l_cher(uplo, uplo_s, n, alpha, x, incx, a, lda)
	Mchar           *uplo;
	unsigned        uplo_s;
	Mint            *n;
	Mfloat          *alpha;
	Mf_complex       x[];
	Mint            *incx;
	Mf_complex       a[];
	Mint            *lda;
#endif
{
	Mint             lower, upper;
	Mint             _l0, _l1, ix, j;
	Mf_complex       temp;


	upper = imsl_l1ame(uplo, uplo_s, "U", sizeof("U"));
	lower = imsl_l1ame(uplo, uplo_s, "L", sizeof("L"));
	/*
	 * Test the input parameters.
	 */
	if (*n < 0) {
		imsl_e1psh("l_cher");
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "N must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
		imsl_e1pop("l_cher");
		goto L_9000;
	} else if ((*lda < *n) || (*lda == 0)) {
		imsl_e1psh("l_cher");
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 2, "LDA must be greater than or equal to N and greater than zero while LDA = %(i1) and N = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_LDA_VALUE_GIVEN);
		imsl_e1pop("l_cher");
		goto L_9000;
	} else if (*incx == 0) {
		imsl_e1psh("l_cher");
		imsl_e1sti(1, *incx);

/*		imsl_ermes(5, 3, "INCX must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
		imsl_e1pop("l_cher");
		goto L_9000;
	} else if ((!upper) && (!lower)) {
		imsl_e1psh("l_cher");
		imsl_e1stl(1, uplo);

/*		imsl_ermes(5, 4, "UPLO must be set equal to U or L while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_UPLO_VALUE);
		imsl_e1pop("l_cher");
		goto L_9000;
	}
	/* Quick return if possible */
	if (*n == 0 || *alpha == F_ZERO)
		goto L_9000;

	ix = 1;
	if (*incx < 0)
		ix = (-*n + 1) ** incx + 1;

	for (j = 1; j <= *n; j++) {
		temp = imsl_c_mul(imsl_cf_convert(*alpha, F_ZERO), imsl_c_conjg(x[ix - 1]));
		if (upper) {
			if (*incx >= 0) {
                                _l0 = j - 1; _l1 = 1;
				imsl_caxpy(&_l0, &temp, x, incx, &a[*lda * (j - 1)],
					   &_l1);
			} else {
                                _l0 = j - 1; _l1 = 1;
				imsl_caxpy(&_l0, &temp, &x[ix - *incx - 1],
				     incx, &a[*lda * (j - 1)], &_l1);
			}
		} else {
			if (*incx >= 0) {
                                _l0 = *n - j; _l1 = 1;
				imsl_caxpy(&_l0, &temp, &x[ix + *incx - 1],
				 incx, &a[*lda * (j - 1) + j], &_l1);
			} else {
                                _l0 = *n - j; _l1 = 1;
				imsl_caxpy(&_l0, &temp, x, incx, &a[*lda * (j - 1) + j],
					   &_l1);
			}
		}
		a[*lda * (j - 1) + j - 1] = imsl_cf_convert(imsl_fc_convert(a[*lda * (j - 1) + j - 1]) +
			     imsl_fc_convert(imsl_c_mul(x[ix - 1], temp)), F_ZERO);
		ix += *incx;
	}

L_9000:
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  CHER2  (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 14, 1989

    Purpose:    Perform a rank-two matrix update to the Hermitian
                matrix A, A = A + alpha*x*ctrans(y) +
                imsl_c_conjg(alpha)*y*ctrans(x).  Here ctrans( ) represents the
                conjugate transpose of the vectors.

    Usage:      CALL CHER2 (UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA)

    Arguments:
       UPLO   - Character specifying the matrix storage.  (Input)
                If UPLO is 'U' of 'u', then the upper half of A is used.
                If UPLO is 'L' of 'l', then the lower half of A is used.
       N      - Order of the matrix A.  (Input)
       ALPHA  - Complex scalar multiplier for the vector-vector product.
                (Input)
       X      - Complex vector of length (N-1)*IABS(INCX)+1.  (Input)
       INCX   - Displacement between elements of X.  (Input)
       Y      - Complex vector of length (N-1)*IABS(INCY)+1.  (Input)
       INCY   - Displacement between elements of Y.  (Input)
       A      - Complex array of order N.  (Input/Output)
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
#ifdef  ANSI
static void l_cher2(Mchar *uplo, unsigned uplo_s, Mint *n,
                    Mf_complex *alpha, Mf_complex *x, Mint *incx, 
                    Mf_complex *y, Mint *incy, Mf_complex *a, Mint *lda)
#else
static void l_cher2(uplo, uplo_s, n, alpha, x, incx, y, incy,
	   a, lda)
	Mchar           *uplo;
	unsigned        uplo_s;
	Mint            *n;
	Mf_complex      *alpha, x[];
	Mint            *incx;
	Mf_complex       y[];
	Mint            *incy;
	Mf_complex      *a;
	Mint            *lda;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
        Mint	        lower, upper;
	Mint             _l0, _l1, ix, iy, j;
	Mf_complex       tempx, tempy;


	upper = imsl_l1ame(uplo, uplo_s, "U", sizeof("U"));
	lower = imsl_l1ame(uplo, uplo_s, "L", sizeof("L"));
	/*
	 * Test the input parameters.
	 */
	if (*n < 0) {
		imsl_e1psh("l_cher2");
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "N must be greater than or equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
		imsl_e1pop("l_cher2");
		goto L_9000;
	} else if ((*lda < *n) || (*lda == 0)) {
		imsl_e1psh("l_cher2");
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 2, "LDA must be greater than or equal to N and greater than zero while LDA = %(i1) and N = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_LDA_VALUE_GIVEN);
		imsl_e1pop("l_cher2");
		goto L_9000;
	} else if (*incx == 0) {
		imsl_e1psh("l_cher2");
		imsl_e1sti(1, *incx);

/*		imsl_ermes(5, 3, "INCX must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
		imsl_e1pop("l_cher2");
		goto L_9000;
	} else if (*incy == 0) {
		imsl_e1psh("l_cher2");
		imsl_e1sti(1, *incy);

/*		imsl_ermes(5, 4, "INCY must not be equal to zero while %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INCY_EQUALS_ZERO);
		imsl_e1pop("l_cher2");
		goto L_9000;
	} else if ((!upper) && (!lower)) {
		imsl_e1psh("l_cher2");
		imsl_e1stl(1, uplo);

/*		imsl_ermes(5, 5, "UPLO must be set equal to U or L while %(l1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_UPLO_VALUE);
		imsl_e1pop("l_cher2");
		goto L_9000;
	}
	/* Quick return if possible */
	if (*n == 0 || imsl_c_eq(*alpha, imsl_cf_convert(F_ZERO, F_ZERO)))
		goto L_9000;

	ix = 1;
	iy = 1;
	if (*incx < 0)
		ix = 1 - (*n - 1) ** incx;
	if (*incy < 0)
		iy = 1 - (*n - 1) ** incy;

	for (j = 1; j <= *n; j++) {
		tempx = imsl_c_conjg(imsl_c_mul(*alpha, x[ix - 1]));
		tempy = imsl_c_mul(*alpha, imsl_c_conjg(y[iy - 1]));
		if (upper) {
			if (*incx >= 0) {
                                _l0 = j - 1; _l1 = 1;
				imsl_caxpy(&_l0, &tempy, x, incx, A(j - 1, 0),
					   &_l1);
			} else {
                                _l0 = j - 1; _l1 = 1;
				imsl_caxpy(&_l0, &tempy, &x[ix - *incx - 1],
					   incx, A(j - 1, 0), &_l1);
			}
			if (*incy >= 0) {
                                _l0 = j - 1; _l1 = 1;
				imsl_caxpy(&_l0, &tempx, y, incy, A(j - 1, 0),
					   &_l1);
			} else {
                                _l0 = j - 1; _l1 = 1;
				imsl_caxpy(&_l0, &tempx, &y[iy - *incy - 1],
					   incy, A(j - 1, 0), &_l1);
			}
		} else {
			if (*incx >= 0) {
                                _l0 = *n - j; _l1 = 1;
				imsl_caxpy(&_l0, &tempy, &x[ix + *incx - 1],
					   incx, A(j - 1, j), &_l1);
			} else {
                                _l0 = *n - j; _l1 = 1;
				imsl_caxpy(&_l0, &tempy, x, incx, A(j - 1, j),
					   &_l1);
			}
			if (*incy >= 0) {
                                _l0 = *n - j; _l1 = 1;
				imsl_caxpy(&_l0, &tempx, &y[iy + *incy - 1],
					   incy, A(j - 1, j), &_l1);
			} else {
                                _l0 = *n - j; _l1 = 1;
				imsl_caxpy(&_l0, &tempx, y, incy, A(j - 1, j),
					   &_l1);
			}
		}
		*A(j - 1, j - 1) = imsl_cf_convert(imsl_fc_convert(imsl_c_add(imsl_c_add(*A(j - 1, j - 1),
										   imsl_c_mul(y[iy - 1], tempx)), imsl_c_mul(x[ix - 1], tempy))), F_ZERO);
		ix += *incx;
		iy += *incy;
	}

L_9000:
	return;
}				/* end of function */
#undef  A

/* Structured by FOR_STRUCT, v0.2, on 09/19/90 at 16:29:40
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E3BSF/DE3BSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 31, 1985

    Purpose:    Compute the eigenvalues in a given range of a symmetric
                tridiagonal matrix.

    Usage:      CALL E3BSF (N, MXEVAL, ELOW, EHIGH, NEVAL, EVAL, D, E,
                            E2, WK1, WK2, IND)

    Arguments:
       N      - Order of the matrix.  (Input)
       MXEVAL - Maximum number of eigenvalues to be computed.  (Input)
       ELOW   - Lower limit of the interval in which the eigenvalues are
                sought.  (Input)
       EHIGH  - Upper limit of the interval in which the eigenvalues are
                sought.  (Input)
       NEVAL  - Number of eigenvalues found.  (Output)
       EVAL   - Real vector of length at least MXEVAL containing the
                eigenvalues of A in the interval (ELOW,EHIGH) in
                increasing order.  (Output)
                Only the first NEVAL elements of EVAL are significant.
       D      - Real vector of length N containing the main diagonal
                of the symmetric tridiagonal matrix.  (Input)
       E      - Real vector of length N containing, in its last N-1
                positions, the subdiagonal elements of the symmetric
                tridiagonal matrix.  (Input)
       E2     - Real vector of length N containing, in the last N-1
                positions, the squares of the subdiagonal elements
                of the symmetric tridiagonal matrix.  (Output)
       WK1    - Real work vector of length N.  (Output)
       WK2    - Real work vector of length N.  (Output)
       IND    - Integer vector of length at least MXEVAL containing the
                submatrix indices associated with the corresponding
                NEVAL eigenvalues in EVAL.  (Output)

    Remark:
       E3BSF is associated with EISPACK routine BISECT

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e3bsf(Mint *n, Mint *mxeval, Mfloat *elow, Mfloat *ehigh,
			Mint *neval, Mfloat eval[], Mfloat d[],
			Mfloat e[], Mfloat e2[], Mfloat wk1[],
			Mfloat wk2[], Mint ind[])
#else
static void l_e3bsf(n, mxeval, elow, ehigh, neval, eval, d, e,
	   e2, wk1, wk2, ind)
	Mint            *n, *mxeval;
	Mfloat          *elow, *ehigh;
	Mint            *neval;
	Mfloat           eval[], d[], e[], e2[], wk1[], wk2[];
	Mint             ind[];
#endif
{
	Mint             i, ip, iq, ir, is, isturm, itag, j, k, l, m1, m2;
	Mfloat           elb, eps, eps1, t1, t2, u, ub, v, x0, x1, xu;


	imsl_e1psh("l_e3bsf");
	/* INITIALIZATIONS */
	eps = imsl_amach(4);
	elb = *elow;
	ub = *ehigh;
	itag = 0;
	t1 = *elow;
	t2 = *ehigh;
	eps1 = F_ZERO;
	/*
	 * LOOK FOR SMALL SUB-DIAGONAL ENTRIES
	 */
	e2[0] = F_ZERO;
	for (i = 2; i <= *n; i++) {
		e2[i - 1] = e[i - 1] * e[i - 1];
		if (fabs(e[i - 1]) <= eps * (fabs(d[i - 1]) + fabs(d[i - 2])))
			e2[i - 1] = F_ZERO;
	}
	/*
	 * DETERMINE THE NUMBER OF EIGENVALUES IN THE INTERVAL
	 */
	ip = 1;
	iq = *n;
	x1 = ub;
	isturm = 1;
	goto L_130;
L_20:
	*neval = is;
	x1 = elb;
	isturm = 2;
	goto L_130;
L_30:
	*neval -= is;
	if (*neval > *mxeval) {
		/*
		 * UNDERESTIMATE OF NUMBER OF EIGENVALUES IN INTERVAL
		 */
		imsl_e1sti(1, *neval);
		imsl_e1sti(2, *mxeval);
		imsl_e1str(3, *elow);
		imsl_e1str(4, *ehigh);

/*		imsl_ermes(3, 1, "The determined number of eigenvalues in the interval (%(r3),%(r4)) is NEVAL = %(i1).  However, the input value for the maximum number of eigenvalues in this interval is MXEVAL = %(i2).");
*/
                imsl_ermes(IMSL_WARNING, IMSL_NEVAL_MXEVAL_MISMATCH);
		goto L_9000;
	}
	iq = 0;
	ir = 0;
	/*
	 * ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING INTERVAL BY THE
	 * GERSCHGORIN BOUNDS
	 */
L_40:
	if (ir == *neval)
		goto L_9000;
	itag += 1;
	ip = iq + 1;
	xu = d[ip - 1];
	x0 = d[ip - 1];
	u = F_ZERO;

	for (iq = ip; iq <= *n; iq++) {
		x1 = u;
		u = F_ZERO;
		v = F_ZERO;
		if (iq != *n) {
			u = fabs(e[iq]);
			v = e2[iq];
		}
		xu = imsl_f_min(d[iq - 1] - (x1 + u), xu);
		x0 = imsl_f_max(d[iq - 1] + (x1 + u), x0);
		if (v == F_ZERO)
			goto L_60;
	}

L_60:
	x1 = imsl_f_max(fabs(xu), fabs(x0)) * eps;
	if (eps1 <= F_ZERO)
		eps1 = -x1;
	if (ip == iq) {
		/*
		 * CHECK FOR ISOLATED ROOT WITHIN INTERVAL
		 */
		if (d[ip - 1] <= t1 || d[ip - 1] >= t2)
			goto L_200;
		m1 = ip;
		m2 = ip;
		wk2[ip - 1] = d[ip - 1];
		goto L_180;
	}
	x1 *= iq - ip + 1;
	elb = imsl_f_max(t1, xu - x1);
	ub = imsl_f_min(t2, x0 + x1);
	x1 = elb;
	isturm = 3;
	goto L_130;
L_70:
	m1 = is + 1;
	x1 = ub;
	isturm = 4;
	goto L_130;
L_80:
	m2 = is;
	if (m1 > m2)
		goto L_200;
	/* FIND ROOTS BY BISECTION */
	x0 = ub;
	isturm = 5;

	sset(m2 - m1 + 1, ub, &wk2[m1 - 1], 1);
	sset(m2 - m1 + 1, elb, &wk1[m1 - 1], 1);
	/*
	 * LOOP FOR K-TH EIGENVALUE FOR K=M2 STEP -1 UNTIL M1 DO -- (-DO- NOT
	 * USED TO LEGALIZE -COMPUTED GO TO-)
	 */
	k = m2;
L_90:
	xu = elb;
	for (i = k; i >= m1; i--) {
		if (xu < wk1[i - 1]) {
			xu = wk1[i - 1];
			goto L_110;
		}
	}

L_110:
	if (x0 > wk2[k - 1])
		x0 = wk2[k - 1];
	/* NEXT BISECTION STEP */
L_120:
	x1 = (xu + x0) * F_HALF;
	if ((x0 - xu) <= (F_TWO * eps * (fabs(xu) + fabs(x0)) + fabs(eps1)))
		goto L_170;
	/* IN-LINE PROCEDURE FOR STURM SEQUENCE */
L_130:
	is = ip - 1;
	u = F_ONE;

	for (i = ip; i <= iq; i++) {
		if (u == F_ZERO) {
			v = fabs(e[i - 1]) / eps;
			if (e2[i - 1] == F_ZERO)
				v = F_ZERO;
		} else {
			v = e2[i - 1] / u;
		}
		u = d[i - 1] - x1 - v;
		if (u < F_ZERO)
			is += 1;
	}

	if (isturm == 1)
		goto L_20;
	if (isturm == 2)
		goto L_30;
	if (isturm == 3)
		goto L_70;
	if (isturm == 4)
		goto L_80;
	/* REFINE INTERVALS */
L_150:
	if (is < k) {
		xu = x1;
		if (is < m1) {
			wk1[m1 - 1] = x1;
			goto L_120;
		}
		wk1[is] = x1;
		if (wk2[is - 1] > x1)
			wk2[is - 1] = x1;
		goto L_120;
	}
L_160:
	x0 = x1;
	goto L_120;
	/* K-TH EIGENVALUE FOUND */
L_170:
	wk2[k - 1] = x1;
	k -= 1;
	if (k >= m1)
		goto L_90;
	/*
	 * ORDER EIGENVALUES TAGGED WITH THEIR SUBMATRIX ASSOCIATIONS
	 */
L_180:
	is = ir;
	ir += m2 - m1 + 1;
	j = 1;
	k = m1;

	for (l = 1; l <= ir; l++) {
		if (j > is) {
			eval[l - 1] = wk2[k - 1];
			ind[l - 1] = itag;
			k += 1;
		} else {
			if (k > m2)
				goto L_200;
			if (wk2[k - 1] < eval[l - 1]) {
				scopy(is - j + 1, &eval[l - 1], -1, &eval[l], -1);
				icopy(is - j + 1, &ind[l - 1], -1, &ind[l], -1);
				eval[l - 1] = wk2[k - 1];
				ind[l - 1] = itag;
				k += 1;
			} else {
				j += 1;
			}
		}
	}
L_200:
	if (iq < *n)
		goto L_40;

L_9000:
	imsl_e1pop("l_e3bsf");

	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 09/19/90 at 17:15:23
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CRRCR/DCRRCR (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 8, 1986

    Purpose:    Copy a real rectangular matrix to a d_complex rectangular
                matrix.

    Usage:      CALL CRRCR (NRA, NCA, A, LDA, NRB, NCB, B, LDB)

    Arguments:
       NRA    - Number of rows in A.  (Input)
       NCA    - Number of columns in A.  (Input)
       A      - Real NRA by NCA rectangular matrix.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       NRB    - Number of rows in B.  (Input)
                It must be the same as NRA.
       NCB    - Number of columns in B.  (Input)
                It must be the same as NCA.
       B      - Complex NRB by NCB rectangular matrix containing a copy
                of A.  (Output)
       LDB    - Leading dimension of B exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remark:
       The matrices A and B may be the same.

    Keyword:    Basic matrix operation

    GAMS:       D1b9

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_crrcr(Mint *nra, Mint *nca, Mfloat *a, Mint *lda,
			Mint *nrb, Mint *ncb, Mf_complex *b,
			Mint *ldb)
#else
static void l_crrcr(nra, nca, a, lda, nrb, ncb, b, ldb)
	Mint            *nra, *nca;
	Mfloat          *a;
	Mint            *lda, *nrb, *ncb;
	Mf_complex      *b;
	Mint            *ldb;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define B(I_,J_)	(b+(I_)*(*ldb)+(J_))
	Mint             i, j;


	imsl_e1psh("CRRCR ");
	/* Check NRA */
	if (*nra < 1) {
		imsl_e1sti(1, *nra);

/*		imsl_ermes(5, 1, "The number of rows in A must be at least 1 while NRA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NRA_MUST_BE_AT_LEAST_1);
	}
	/* Check NCA */
	if (*nca < 1) {
		imsl_e1sti(1, *nca);

/*		imsl_ermes(5, 2, "The number of columns in A must be at least 1 while NCA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NCA_MUST_BE_AT_LEAST_1);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Check NRB */
	if (*nra != *nrb) {
		imsl_e1sti(1, *nra);
		imsl_e1sti(2, *nrb);

/*		imsl_ermes(5, 3, "The number of rows in B must be equal to the number of rows in A while NRB = %(i1) and NRA = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_NRB_EQUAL_TO_NRA);
	}
	/* Check NCB */
	if (*nca != *ncb) {
		imsl_e1sti(1, *nca);
		imsl_e1sti(2, *ncb);

/*		imsl_ermes(5, 4, "The number of columns in B must be equal to the number of columns in A while NCB = %(i1) and NCA = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_NCB_EQUAL_TO_NCA);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Check LDA */
	if (*lda < *nra) {
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *nra);

/*		imsl_ermes(5, 5, "The leading dimension of A must be at least as large as the number of rows in A while LDA = %(i1) and NRA = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDA_GE_NRA);
	}
	/* Check LDB */
	if (*ldb < *nrb) {
		imsl_e1sti(1, *ldb);
		imsl_e1sti(2, *nrb);

/*		imsl_ermes(5, 6, "The leading dimension of B must be at least as large as the number of rows in B while LDB = %(i1) and NRB = %(i2) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDB_GE_NRB);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Copy A to B */
	for (j = *nca; j >= 1; j--) {
		for (i = *nra; i >= 1; i--) {
			*B(j - 1, i - 1) = imsl_cf_convert(*A(j - 1, i - 1), F_ZERO);
		}
	}

L_9000:
	;
	imsl_e1pop("CRRCR ");
	return;
}				/* end of function */

#undef A
#undef B
/* Structured by FOR_STRUCT, v0.2, on 09/19/90 at 17:16:20
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E3ESF/DE3ESF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 15, 1985

    Purpose:    Computes the eigenvectors of a symmetric tridiagonal
                matrix corresponding to a set of ordered approximate
                eigenvalues, using inverse iteration.

    Usage:      CALL E3ESF (N, NEVEC, EVAL, EVEC, LDEVEC, D, E, E2,
                            WK1, WK2, WK3, WK4, WK5, IND)

    Arguments:
       N      - Order of the matrix.  (Input)
       NEVEC  - Number of eigenvectors to be computed.  (Input)
       EVAL   - Real vector of length NEVEC containing the eigenvalues
                the eigenvector computation is based on.  (Input)
       EVEC   - Real matrix of order N containing, in its first NEVEC
                columns, the eigenvectors associated with the eigenvalues
                in EVAL.  Eigenvector J corresponds to eigenvalue J for
                J = 1 through NEVEC.  (Output)
       LDEVEC - Leading dimension of EVEC exactly as specified in
                the dimension statement of the calling program.  (Input)
       D      - Real vector of length containig the diagonal elements
                of the tri-diagonal matrix.  (Input)
       E      - Real vector of length N.  On input, the elemets of the
                off diagonal.  E(1) is arbitrary.  (Input)
       E2     - Real vector of length N containing the elements of the
                off diagonal squared.  (Output)
       WK1    - Real work vector of length N.  (Output)
       WK2    - Real work vector of length N.  (Output)
       WK3    - Real work vector of length N.  (Output)
       WK4    - Real work vector of length N.  (Output)
       WK5    - Real work vector of length N.  (Output)
       IND    - Integer vector of length NEVEC containing the submatrix
                indices associated with the corresponding NEVEC
                eigenvalues in EVAL.  (Input)

    Remark:
      This routine is based on the EISPACK routine TINVIT.

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e3esf(Mint *n, Mint *nevec, Mfloat eval[], Mfloat *evec,
			Mint *ldevec, Mfloat d[], Mfloat e[], Mfloat e2[],
			Mfloat wk1[], Mfloat wk2[], Mfloat wk3[], Mfloat wk4[],
			Mfloat wk5[], Mint ind[])
#else
static void l_e3esf(n, nevec, eval, evec, ldevec, d, e, e2, wk1,
	   wk2, wk3, wk4, wk5, ind)
	Mint            *n, *nevec;
	Mfloat           eval[], *evec;
	Mint            *ldevec;
	Mfloat           d[], e[], e2[], wk1[], wk2[], wk3[], wk4[], wk5[];
	Mint             ind[];
#endif
{
#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))
	Mint       err;
	Mint             _l0, _l1, i, igroup, iq, ir, is, itag, its, j,
	                kp;
	Mfloat           anorm, eps, eps2, eps3, eps4, u,
	                uk, v, vnorm, x0, x1, xu;


	imsl_e1psh("E3ESF ");

	eps = imsl_amach(4);
	e2[0] = F_ZERO;
	for (i = 2; i <= *n; i++) {
		e2[i - 1] = imsl_fi_power(e[i - 1], 2);
		if (fabs(e[i - 1]) <= eps * (fabs(d[i - 1]) + fabs(d[i - 2])))
			e2[i - 1] = F_ZERO;
	}
	itag = 0;
	iq = 0;
/*	err = FALSE;  */
	err = 0;
	/* Compute norm of A */
	anorm = fabs(d[0]);
	for (i = 2; i <= *n; i++) {
		anorm = imsl_f_max(anorm, fabs(d[i - 1]) + fabs(e[i - 1]));
	}
	/* Establish and process next submatrix */
L_30:
	kp = iq + 1;

	for (iq = kp; iq <= *n; iq++) {
		if (iq == *n)
			goto L_50;
		if (e2[iq] == F_ZERO)
			goto L_50;
	}
	/* Find vectors by inverse iteration */
L_50:
	;
	itag += 1;
	is = 0;

	for (ir = 1; ir <= *nevec; ir++) {
		if (ind[ir - 1] != itag)
			goto L_140;
		its = 1;
		x1 = eval[ir - 1];
		if (is == 0) {
			/* Check for isolated root */
			xu = F_ONE;
			if (iq - kp + 1 == 1) {
				wk5[kp - 1] = F_ONE;
				goto L_130;
			}
			vnorm = fabs(d[kp - 1]);
			for (i = kp + 1; i <= iq; i++) {
				vnorm = imsl_f_max(vnorm, fabs(d[i - 1]) + fabs(e[i - 1]));
			}
			/*
			 * EPS2 is the criterion for grouping, EPS3 replaces
			 * zero pivots and equal roots are modified by EPS3,
			 * EPS4 is taken very small to avoid overflow The
			 * tolerance criterion for grouping has been changed
			 * from eps2=sqrt(eps) vnorm to eps*vnorm.
			 */
			eps2 = eps * vnorm;
			eps3 = eps * vnorm;
			uk = iq - kp + 1;
			eps4 = uk * eps3;
			eps4 *= 1000.0e0;
			uk = eps4 / sqrt(uk);
			is = kp;
			igroup = 0;
		} else {
			/* Look for close or coincident roots */
			if (fabs(x1 - x0) < eps2) {
				igroup += 1;
				if (x1 <= x0)
					x1 = x0 + eps3;
			} else {
				igroup = 0;
			}
		}
		/*
		 * Elimination with interchanges and initialization of vector
		 */
		v = F_ZERO;

		for (i = kp; i <= iq; i++) {
			wk5[i - 1] = uk;
			if (i == kp) {
				u = d[i - 1] - x1 - xu * v;
				if (i != iq)
					v = e[i];
			} else {
				if (fabs(e[i - 1]) < fabs(u)) {
					xu = e[i - 1] / u;
					wk4[i - 1] = xu;
					wk1[i - 2] = u;
					wk2[i - 2] = v;
					wk3[i - 2] = F_ZERO;
					u = d[i - 1] - x1 - xu * v;
					if (i != iq)
						v = e[i];
				} else {
					xu = u / e[i - 1];
					wk4[i - 1] = xu;
					wk1[i - 2] = e[i - 1];
					wk2[i - 2] = d[i - 1] - x1;
					wk3[i - 2] = F_ZERO;
					if (i != iq)
						wk3[i - 2] = e[i];
					u = v - xu * wk2[i - 2];
					v = -xu * wk3[i - 2];
				}
			}
		}

		if (fabs(u) < eps * eps3)
			u = eps3;
		wk1[iq - 1] = u;
		wk2[iq - 1] = F_ZERO;
		wk3[iq - 1] = F_ZERO;
		/* Back substitution for I=IQ */
L_80:
		for (i = iq; i >= kp; i--) {
			wk5[i - 1] = (wk5[i - 1] - u * wk2[i - 1] - v * wk3[i - 1]) /
				wk1[i - 1];
			v = u;
			u = wk5[i - 1];
		}
		/*
		 * Orthogonalize with respect to previous members of IGROUP
		 */
		if (igroup != 0) {
			for (j = ir - 1; j <= (ir - igroup); j++) {
				if (ind[j - 1] == itag) {
					xu = imsl_sdot(iq - kp + 1, &wk5[kp - 1], 1, EVEC(j - 1, kp - 1),
						       1);
					saxpy(iq - kp + 1, -xu, EVEC(j - 1, kp - 1), 1,
						   &wk5[kp - 1], 1);
				}
			}
		}
		_l0 = iq-kp+1;
		_l1 = 1;
		vnorm = imsl_sasum(_l0, &wk5[kp - 1], _l1);
		if (vnorm < F_ONE) {
			/* Forward substitution */
			if (its < 5) {
				if (vnorm <= eps3) {
					wk5[is - 1] = eps4;
					is += 1;
					if (is > iq)
						is = kp;
				} else {
					xu = eps4 / vnorm;
					sscal(iq - kp + 1, xu, &wk5[kp - 1], 1);
				}
				/*
				 * Elimination operations on next vector
				 * iterate
				 */
		L_110:
				for (i = kp + 1; i <= iq; i++) {
					u = wk5[i - 1];
					/*
					 * If WK1(I-1) .EQ. E(I), a row
					 * interchange was performed earlier
					 * in the triangularization process
					 */
					if (wk1[i - 2] == e[i - 1]) {
						u = wk5[i - 2];
						wk5[i - 2] = wk5[i - 1];
					}
					wk5[i - 1] = u - wk4[i - 1] * wk5[i - 2];
				}
				its += 1;
				goto L_80;
			}
			/* Set error --non-converged eigenvector */
			if (!err) {

/*				imsl_ermes(3, 2, "The iteration for at least one eigenvector failed to converge.  Some of the eigenvectors may be inaccurate.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_LOST_ORTHOGONALITY);
			/*	err = TRUE; */
				err = 1;
				xu = F_ZERO;
			}
		}
		/* Store eigenvector in EVEC */
L_130:
		sset(*n, F_ZERO, EVEC(ir - 1, 0), 1);
		saxpy(iq - kp + 1, xu, &wk5[kp - 1], 1, EVEC(ir - 1, kp - 1),
			   1);
		x0 = x1;
L_140:
		;
	}

	if (iq < *n)
		goto L_30;

L_9000:
	imsl_e1pop("E3ESF ");

	return;
}				/* end of function */

#undef EVEC

/* Structured by FOR_STRUCT, v0.2, on 10/18/90 at 15:45:58
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E3BHF/DE3BHF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 12, 1990

    Purpose:    Compute the eigenvalues in a given range of a d_complex
                Hermitian matrix.

    Usage:      CALL E3BHF (N, MXEVAL, A, LDA, ELOW, EHIGH, NEVAL, EVAL,
                            ACOPY, RWK, CWK, IWK)

    Arguments:  (See EVBHF)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e3bhf(Mint *n, Mint *mxeval, Mf_complex *a, Mint *lda,
		Mfloat *elow, Mfloat *ehigh, Mint *neval, Mfloat eval[],
		Mf_complex *acopy, Mfloat *rwk, Mf_complex *cwk,
		Mint iwk[])
#else
static void l_e3bhf(n, mxeval, a, lda, elow, ehigh, neval, eval,
      acopy, rwk, cwk, iwk)
	Mint            *n, *mxeval;
	Mf_complex      *a;
	Mint            *lda;
	Mfloat          *elow, *ehigh;
	Mint            *neval;
	Mfloat           eval[];
	Mf_complex      *acopy;
	Mfloat          *rwk;
	Mf_complex      *cwk;
	Mint             iwk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define ACOPY(I_,J_)	(acopy+(I_)*(*n)+(J_))
#define RWK(I_,J_)	(rwk+(I_)*(*n)+(J_))
#define CWK(I_,J_)	(cwk+(I_)*(*n)+(J_))
	Mint             i, j, k;


	imsl_e1psh("E3BHF ");
	/* Check N */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "The argument N = %(i1).  The order of the matrix must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);
		goto L_9000;
	}
	/* Check MXEVAL */
	if (*mxeval <= 0 || *mxeval > *n) {
		imsl_e1sti(1, *mxeval);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 2, "The argument MXEVAL = %(i1).  The maximum number of eigenvalues to be calculated must be greater than 0 and less than or equal to N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NUMBER_MAX_EIGENVALUES);
	}
	/* Check LDA */
	if (*lda < *n) {
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 3, "The argument LDA = %(i1).  The leading dimension of the matrix must be at least equal to the order of the matrix, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);
	}
	/* Check ELOW against EHIGH */
	if (*elow >= *ehigh) {
		imsl_e1str(1, *elow);
		imsl_e1str(2, *ehigh);

/*		imsl_ermes(5, 4, "The lower limit of the interval in which the eigenvalues are to be sought, ELOW = %(r1), must be strictly less than the upper limit, EHIGH = %(r2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_ELOW_LESS_THAN_EHIGH);
	}
	if (imsl_n1rty(0) > 0)
		goto L_9000;
	/* Copy A to ACOPY */
	imsl_ccgcg(n, a, lda, acopy, n);
	/* Extend ACOPY to its lower half */
	l_chfcg(n, acopy, n);
	if (imsl_n1rty(1) == 3) {
		imsl_e1mes(-1, 2, " ");
		goto L_9000;
	}
	if (imsl_n1rty(0) > 0)
		goto L_9000;
	/* Reduce to symmetric triagonal */
	l_e3chf(n, acopy, n, RWK(0, 0), RWK(1, 0), CWK(0, 0), CWK(1, 0));
	/* Find eigenvalues */
	l_e3bsf(n, mxeval, elow, ehigh, neval, eval, RWK(0, 0), RWK(1, 0),
		   RWK(2, 0), RWK(3, 0), RWK(4, 0), iwk);
        if (*neval == 0) goto L_9000;
	/*
	 * Sort the eigenvalue magnitudes. Ultimately want eigenvalues in
	 * descending magnitude.
	 */
	for (i = 1; i <= *neval; i++) {
		iwk[i - 1] = i;
		*RWK(0, i - 1) = -fabs(eval[i - 1]);
	}


	imsl_svrgp(*neval, RWK(0, 0), RWK(0, 0), iwk);
	/* Move the eigenvalues */
	for (i = 1; i <= *neval; i++) {
		for (j = i; j <= *neval; j++) {
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
	for (i = *neval - 1; i >= 1; i--) {
		sswap(1, &eval[i - 1], 1, &eval[iwk[i - 1] - 1], 1);
	}

L_9000:
	imsl_e1pop("E3BHF ");
	return;
}				/* end of function */

#undef A
#undef ACOPY
#undef RWK
#undef CWK


/* Structured by FOR_STRUCT, v0.2, on 10/18/90 at 15:44:30
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E3LHF/DE3LHF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 12, 1990

    Purpose:    Compute all of the eigenvalues of a d_complex Hermitian
                matrix.

    Usage:      CALL E3LHF (N, A, LDA, EVAL, ACOPY, RWK, CWK, IWK)

    Arguments:  (See EVLHF)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e3lhf(Mint *n, Mf_complex *a, Mint *lda, Mfloat eval[],
		Mf_complex *acopy, Mfloat rwk[], Mf_complex cwk[],
		Mint iwk[])
#else
static void l_e3lhf(n, a, lda, eval, acopy, rwk, cwk, iwk)
	Mint            *n;
	Mf_complex      *a;
	Mint            *lda;
	Mfloat           eval[];
	Mf_complex      *acopy;
	Mfloat           rwk[];
	Mf_complex       cwk[];
	Mint             iwk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define ACOPY(I_,J_)	(acopy+(I_)*(*n)+(J_))
	Mint             _l0, _l1, i, j, k;
	Mfloat           evec[1][1];


	imsl_e1psh("E3LHF ");
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

/*		imsl_ermes(5, 2, "The argument LDA = %(i1).  The leading dimension of the matrix must be at least equal to the order of the matrix, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);
		goto L_9000;
	}
	/* Copy A to ACOPY */
	imsl_ccgcg(n, a, lda, acopy, n);
	/* Extend ACOPY to its lower half */
	l_chfcg(n, acopy, n);
	if (imsl_n1rty(1) == 4) {
		imsl_e1mes(-1, 2, " ");
		goto L_9000;
	}
	if (imsl_n1rty(0) > 0)
		goto L_9000;
	/* Reduce to symmetric triagonal */
	l_e3chf(n, acopy, n, eval, rwk, &cwk[0], &cwk[*n]);
	/* Find eigenvalues */
	_l0 = 0;
	_l1 = 1;
	l_e4csf(n, eval, rwk, &_l0, &evec[0][0], &_l1);
	/*
	 * Sort the eigenvalue magnitudes. Ultimately want eigenvalues in
	 * descending magnitude.
	 */
	for (i = 1; i <= *n; i++) {
		iwk[i - 1] = i;
		rwk[i - 1] = -fabs(eval[i - 1]);
	}


	imsl_svrgp(*n, &rwk[0], &rwk[0], iwk);
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
		sswap(1, &eval[i - 1], 1, &eval[iwk[i - 1] - 1], 1);
	}

L_9000:
	;
	imsl_e1pop("E3LHF ");
	return;
}				/* end of function */

#undef A
#undef ACOPY
/* Structured by FOR_STRUCT, v0.2, on 10/18/90 at 15:45:15
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E5CHF/DE5CHF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 11, 1990

    Purpose:    Compute all of the eigenvalues and eigenvectors of a
                d_complex Hermitian matrix.

    Usage:      CALL E5CHF (N, A, LDA, EVAL, EVEC, LDEVEC, ACOPY, RWK,
                            CWK, IWK)

    Arguments:  (See EVCHF)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e5chf(Mint *n, Mf_complex *a, Mint *lda, Mfloat eval[],
		Mf_complex *evec, Mint *ldevec, Mf_complex *acopy,
		Mfloat rwk[], Mf_complex cwk[], Mint iwk[])
#else
static void l_e5chf(n, a, lda, eval, evec, ldevec, acopy, rwk,
      cwk, iwk)
	Mint            *n;
	Mf_complex      *a;
	Mint            *lda;
	Mfloat           eval[];
	Mf_complex      *evec;
	Mint            *ldevec;
	Mf_complex      *acopy;
	Mfloat           rwk[];
	Mf_complex       cwk[];
	Mint             iwk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))
#define ACOPY(I_,J_)	(acopy+(I_)*(*n)+(J_))
	Mint             _l0, _l1, i, j, k;
	Mfloat           _f0, ri;
	Mf_complex       w;


	imsl_e1psh("E5CHF ");
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

/*		imsl_ermes(5, 2, "The argument LDA = %(i1).  The leading dimension of the matrix must be at least equal to the order, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);
	}
	/* Check LDEVEC */
	if (*ldevec < *n) {
		imsl_e1sti(1, *ldevec);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 3, "The argument LDEVEC = %(i1).  The leading dimension of the eigenvector matrix must be at least equal to the order, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDEVEC_VALUE_TOO_SMALL);
	}
	if (imsl_n1rty(0) > 0)
		goto L_9000;
	/* Copy A to ACOPY */
	for (i = 1; i <= *n; i++) {
		imsl_ccopy(&i, A(i - 1, 0), ADR(_l0, 1), ACOPY(i - 1, 0), ADR(_l1, 1));
	}
	l_chfcg(n, acopy, n);
	/* Extend ACOPY to its lower half */
	if (imsl_n1rty(1) == 4) {
		imsl_e1mes(-1, 2, " ");
		goto L_9000;
	}
	if (imsl_n1rty(0) > 0)
		goto L_9000;
	/* Set real work array to the identity */
	sset(*n ** n, F_ZERO, &rwk[*n], 1);
	sset(*n, F_ONE, &rwk[*n], *n + 1);
	/* Reduce to symmetric tridiagonal */
	l_e3chf(n, acopy, n, eval, rwk, cwk, &cwk[*n]);
	/* Find eigenvalues and eignevectors */
	_l0 = 1;
	l_e4csf(n, eval, rwk, &_l0, &rwk[*n], n);
	if (imsl_n1rty(1) == 4)
		goto L_9000;
	/* Copy vectors to d_complex matrix */
	l_crgcg(n, &rwk[*n], n, evec, ldevec);
	/* Backtransform eigenvectors */
	l_e4chf(n, n, acopy, n, cwk, evec, ldevec, &cwk[*n]);
	/*
	 * Sort the eigenvalue magnitudes. Ultimately want eigenvalues in
	 * descending magnitude.
	 */
	for (i = 1; i <= *n; i++) {
		iwk[i - 1] = i;
		rwk[i - 1] = -fabs(eval[i - 1]);
	}


	imsl_svrgp(*n, &rwk[0], &rwk[0], iwk);
	/* Move the eigenvalues */
	for (i = 1; i <= *n; i++) {
		for (j = i; j <= *n; j++) {
			if (iwk[j - 1] == i) {
				k = iwk[i - 1];
				iwk[i - 1] = j;
				iwk[j - 1] = k;
				goto L_40;
			}
		}
L_40:
		;
	}
	/* Take product of cycles and permute */
	for (i = *n - 1; i >= 1; i--) {
		imsl_cswap(n, EVEC(i - 1, 0), ADR(_l0, 1), EVEC(iwk[i - 1] - 1, 0),
			   ADR(_l1, 1));
		sswap(1, &eval[i - 1], 1, &eval[iwk[i - 1] - 1], 1);
	}

	/*
	 * Normalize eigenvectors using the Euclidean Norm
	 */
	for (j = 1; j <= *n; j++) {
		ri = imsl_scnrm2(n, EVEC(j - 1, 0), ADR(_l0, 1));
		if (ri > F_ZERO) {
			imsl_csscal(n, ADR(_f0, F_ONE / ri), EVEC(j - 1, 0), ADR(_l0, 1));
		}
	}
	/*
	 * Normalize each eigenvector so that its biggest component is real
	 * and positive. The reason for this is that the eigenvectors then
	 * form a right- hand system.
	 */
	for (j = 1; j <= *n; j++) {
		for (i = 1; i <= *n; i++) {
			rwk[i - 1] = imsl_fc_convert(imsl_c_mul(*EVEC(j - 1, i - 1), imsl_c_conjg(*EVEC(j - 1, i - 1))));
		}
		i = imsl_isamax(*n, &rwk[0], 1);
		if (imsl_c_abs(*EVEC(j - 1, i - 1)) == F_ZERO)
			goto L_80;
		/* w = conjugate(evec)/imsl_c_abs(evec) */
		w = imsl_c_div(imsl_c_conjg(*EVEC(j - 1, i - 1)), imsl_cf_convert(imsl_c_abs(*EVEC(j - 1, i - 1)), F_ZERO));
		/* x = x*w */
		imsl_cscal(n, &w, EVEC(j - 1, 0), ADR(_l0, 1));
		*EVEC(j - 1, i - 1) = imsl_cf_convert(imsl_fc_convert(*EVEC(j - 1, i - 1)), F_ZERO);
L_80:
		;
	}

L_9000:
	imsl_e1pop("E5CHF ");
	return;
}				/* end of function */

#undef A
#undef EVEC
#undef ACOPY
/* Structured by FOR_STRUCT, v0.2, on 10/18/90 at 16:14:14
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E3FHF/DE3FHF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 12, 1990

    Purpose:    Compute the eigenvalues in a given range and the
                corresponding eigenvectors of a d_complex Hermitian matrix.

    Usage:      CALL E3FHF (N, MXEVAL, A, LDA, ELOW, EHIGH, NEVAL, EVAL,
                            EVEC, LDEVEC, ACOPY, ECOPY, RWK, CWK, IWK)

    Arguments:  (See EVFHF)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_e3fhf(Mint *n, Mint *mxeval, Mf_complex *a, Mint *lda,
		Mfloat *elow, Mfloat *ehigh, Mint *neval, Mfloat eval[],
		Mf_complex *evec, Mint *ldevec, Mf_complex *acopy,
		Mfloat *ecopy, Mfloat *rwk, Mf_complex *cwk, Mint iwk[])
#else
static void l_e3fhf(n, mxeval, a, lda, elow, ehigh, neval, eval,
      evec, ldevec, acopy, ecopy, rwk, cwk, iwk)
	Mint            *n, *mxeval;
	Mf_complex      *a;
	Mint            *lda;
	Mfloat          *elow, *ehigh;
	Mint            *neval;
	Mfloat           eval[];
	Mf_complex      *evec;
	Mint            *ldevec;
	Mf_complex      *acopy;
	Mfloat          *ecopy, *rwk;
	Mf_complex      *cwk;
	Mint             iwk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))
#define ACOPY(I_,J_)	(acopy+(I_)*(*n)+(J_))
#define ECOPY(I_,J_)	(ecopy+(I_)*(*n)+(J_))
#define RWK(I_,J_)	(rwk+(I_)*(*n)+(J_))
#define CWK(I_,J_)	(cwk+(I_)*(*n)+(J_))
	Mint             _l0, _l1, i, j, k;
	Mfloat           _f0, ri;
	Mf_complex       w;


	imsl_e1psh("E3FHF ");
	/* Check N */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "The argument N = %(i1).  The order of the matrix must be at least 1.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);
		goto L_9000;
	}
	/* Check MXEVAL */
	if (*mxeval <= 0 || *mxeval > *n) {
		imsl_e1sti(1, *mxeval);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 2, "The argument MXEVAL = %(i1).  The maximum number of eigenvalues to be calculated must be greater than 0 and less than or equal to N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NUMBER_MAX_EIGENVALUES);
	}
	/* Check LDA */
	if (*lda < *n) {
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 3, "The argument LDA = %(i1).  The leading dimension of the matrix must be at least equal to the order of the matrix, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);
	}
	/* Check ELOW against EHIGH */
	if (*elow >= *ehigh) {
		imsl_e1str(1, *elow);
		imsl_e1str(2, *ehigh);

/*		imsl_ermes(5, 4, "The lower limit of the interval in which the eigenvalues are to be sought, ELOW = %(r1), must be strictly less than the upper limit, EHIGH = %(r2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_ELOW_LESS_THAN_EHIGH);
	}
	/* Check LDEVEC */
	if (*ldevec < *n) {
		imsl_e1sti(1, *ldevec);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 5, "The argument LDEVEC = %(i1).  The leading dimension of the matrix must be at least equal to the order of the matrix, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDEVEC_VALUE_TOO_SMALL);
	}
	if (imsl_n1rty(0) > 0)
		goto L_9000;
	/* Copy A to ACOPY */
	imsl_ccgcg(n, a, lda, acopy, n);
	/* Extend ACOPY to its lower half */
	l_chfcg(n, acopy, n);
	if (imsl_n1rty(1) == 3) {
		imsl_e1mes(-1, 3, " ");
		goto L_9000;
	}
	if (imsl_n1rty(0) > 0)
		goto L_9000;
	/* Reduce to symmetric triagonal */
	l_e3chf(n, acopy, n, RWK(0, 0), RWK(1, 0), CWK(0, 0), CWK(1, 0));
	/* Find eigenvalues */
	l_e3bsf(n, mxeval, elow, ehigh, neval, eval, RWK(0, 0), RWK(1, 0),
		   RWK(2, 0), RWK(3, 0), RWK(4, 0), iwk);
        if (*neval == 0) goto L_9000;
	if (*neval <= *mxeval) {
		/* Set real EVEC matrix to the identity */
		sset(*n ** neval, F_ZERO, ecopy, 1);
		sset(*neval, F_ONE, ecopy, *n + 1);
		/* Find eigenvectors */
		l_e3esf(n, neval, eval, ecopy, n, RWK(0, 0), RWK(1, 0), RWK(2, 0),
		RWK(3, 0), RWK(4, 0), RWK(5, 0), RWK(6, 0), RWK(7, 0), iwk);
		/* Copy vectors to d_complex matrix */
		l_crrcr(n, neval, ecopy, n, n, neval, evec, ldevec);
		/* Backtransform eigenvectors */
		l_e4chf(n, neval, acopy, n, CWK(0, 0), evec, ldevec, CWK(1, 0));
		/*
		 * Sort the eigenvalue magnitudes. Ultimately want
		 * eigenvalues in descending magnitude.
		 */
		for (i = 1; i <= *neval; i++) {
			iwk[i - 1] = i;
			*RWK(0, i - 1) = -fabs(eval[i - 1]);
		}


		imsl_svrgp(*neval, RWK(0, 0), RWK(0, 0), iwk);
		/* Move the eigenvalues */
		for (i = 1; i <= *neval; i++) {
			for (j = i; j <= *neval; j++) {
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
		for (i = *neval - 1; i >= 1; i--) {
			imsl_cswap(n, EVEC(i - 1, 0), ADR(_l0, 1), EVEC(iwk[i - 1] -
							1, 0), ADR(_l1, 1));
			sswap(1, &eval[i - 1], 1, &eval[iwk[i - 1] - 1], 1);
		}

		/*
		 * Normalize eigenvectors using the Euclidean Norm
		 */
		for (j = 1; j <= *neval; j++) {
			ri = imsl_scnrm2(n, EVEC(j - 1, 0), ADR(_l0, 1));
			if (ri > F_ZERO) {
				imsl_csscal(n, ADR(_f0, F_ONE / ri), EVEC(j - 1, 0), ADR(_l0, 1));
			}
		}
		/*
		 * Normalize each eigenvector so that its biggest component
		 * is real and positive. The reason for this is that the
		 * eigenvectors then form a right- hand system.
		 */
		for (j = 1; j <= *neval; j++) {
			for (i = 1; i <= *neval; i++) {
				*RWK(0, i - 1) = imsl_fc_convert(imsl_c_mul(*EVEC(j - 1, i - 1), imsl_c_conjg(*EVEC(j - 1, i - 1))));
			}
			i = imsl_isamax(*n, RWK(0, 0), 1);
			if (imsl_c_abs(*EVEC(j - 1, i - 1)) == F_ZERO)
				goto L_70;
			/* w = conjugate(evec)/imsl_c_abs(evec) */
			w = imsl_c_div(imsl_c_conjg(*EVEC(j - 1, i - 1)), imsl_cf_convert(imsl_c_abs(*EVEC(j - 1, i - 1)), F_ZERO));
			/* x = x*w */
			imsl_cscal(n, &w, EVEC(j - 1, 0), ADR(_l0, 1));
			*EVEC(j - 1, i - 1) = imsl_cf_convert(imsl_fc_convert(*EVEC(j - 1, i - 1)), F_ZERO);
	L_70:
			;
		}
	}
L_9000:
	imsl_e1pop("E3FHF ");
	return;
}				/* end of function */

#undef A
#undef EVEC
#undef ACOPY
#undef ECOPY
#undef RWK
#undef CWK
