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

#ifdef ANSI
static VA_LIST_HACK l_eig_gen(Mint n, Mfloat *a, va_list argptr);
static void 	l_e3lrg(Mint *n, Mfloat *a, Mint *lda, Mf_complex eval[],
			Mfloat *acopy, Mfloat *wk, Mint *iwk);
static void	l_e4lrg(Mint *nm, Mint *n, Mint *low, Mint *igh,
			Mfloat *a, Mint int_[]);
static void	l_e5crh(Mint *n, Mint *low, Mint *igh, Mfloat *a,
			Mint *lda, Mfloat wr[], Mfloat wi[], Mfloat *ecopy);
static void	l_e5lrg(Mint *lda, Mint *n, Mfloat *a, Mfloat wr[],
			Mfloat wi[]);
static void	l_e6crg(Mint *n, Mint *low, Mint *igh, Mfloat scale[],
			Mint *nvec, Mf_complex *evec, Mint *ldevec);
static void	l_e8crg(Mint *n, Mfloat *a, Mint *lda, Mf_complex eval[],
			Mf_complex *evec, Mint *ldevec, Mfloat *acopy,
			Mfloat *ecopy, Mfloat *wk, Mint iwk[]);
static void	l_e9crg(Mint *n, Mint *low, Mint *igh, Mfloat *a,
			Mint *lda, Mfloat scale[]);
static void	l_eacrg(Mint *n, Mint *low, Mint *igh, Mfloat *a,
			Mint *lda, Mfloat ort[]);
static void	l_ebcrg(Mint *n, Mint *low, Mint *igh, Mfloat *a,
			Mint *lda, Mfloat ort[], Mfloat *evec);
static void	l_eccrg(Mfloat *ar, Mfloat *imsl_ai, Mfloat *br,
			Mfloat *imsl_bi, Mfloat *cr, Mfloat *imsl_ci);
#else
static VA_LIST_HACK	l_eig_gen();
static void 	l_e3lrg();
static void	l_e4lrg();
static void	l_e5crh();
static void	l_e5lrg();
static void	l_e6crg();
static void	l_e8crg();
static void 	l_e9crg();
static void	l_eacrg();
static void	l_ebcrg();
static void	l_eccrg();
#endif

static Mf_complex	*lv_eval;

#ifdef ANSI
Mf_complex *imsl_f_eig_gen(Mint n, Mfloat *a, ...)
#else
Mf_complex *imsl_f_eig_gen(n, a, va_alist)
    Mint	n;
    Mfloat	*a;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr,a);

    E1PSH("imsl_f_eig_gen", "imsl_d_eig_gen");
    lv_eval = NULL;
    IMSL_CALL(l_eig_gen(n, a, argptr));
    va_end(argptr);
    E1POP("imsl_f_eig_gen", "imsl_d_eig_gen");
    return lv_eval;
}


#ifdef ANSI
static VA_LIST_HACK l_eig_gen(Mint n, Mfloat *a, va_list argptr)
#else
static VA_LIST_HACK l_eig_gen(n, a, argptr)
    Mint	n;
    Mfloat	*a;
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
    Mint	    *iwork	= NULL;
    Mfloat	    *acopy	= NULL;
    Mfloat	    *ecopy	= NULL;
    Mfloat	    *work	= NULL;

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
    if (a==NULL) {
        imsl_e1stl (1, "a");
        imsl_ermes(IMSL_TERMINAL, IMSL_REQ_ARGUMENT_IS_NULL);
    }
    if (vectors_user && (evecu==NULL)) {
        imsl_e1stl (1, "evecu");
        imsl_e1stl (2, "IMSL_VECTORS_USER");
        imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);
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

    if (!vectors && !vectors_user) {
	work    = (Mfloat *) imsl_malloc(4*n*sizeof(*work));
	iwork	= (Mint *) imsl_malloc(2*n*sizeof(*iwork));
    }
    if (vectors || vectors_user) {
	ecopy	= (Mfloat *) imsl_malloc((n*n+n)*sizeof(*ecopy));
	work    = (Mfloat *) imsl_malloc(6*n*sizeof(*work));
	iwork	= (Mint *) imsl_malloc(n*sizeof(*iwork));
    }
    if (vectors)
	*evec	= (Mf_complex *) imsl_malloc(n*n*sizeof(**evec));

    acopy = (Mfloat *) imsl_malloc(n*n*sizeof(*acopy));
    if (work==NULL || iwork==NULL || acopy==NULL || 
		(vectors && (*evec==NULL)) || ((vectors || vectors_user)
		&& ecopy==NULL)) {
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

    if (!vectors && !vectors_user)
	l_e3lrg(&n, a, &a_col_dim, lv_eval, acopy, work, iwork);

    if (vectors) { 
        imsl_f_m1ran(n,a_col_dim, a, a);
        l_e8crg(&n, a, &n, lv_eval, *evec, &n, acopy, ecopy,
			 work, iwork);
	imsl_f_m1ran(a_col_dim, n, a, a);
	imsl_c_m1ran(n, evecu_col_dim, *evec, *evec); 
	}
    if (vectors_user) {
        imsl_f_m1ran(n,a_col_dim, a, a);
        l_e8crg(&n, a, &n, lv_eval, evecu, &evecu_col_dim,
			 acopy, ecopy, work, iwork);
	imsl_f_m1ran(a_col_dim, n, a, a);
	imsl_c_m1ran(n, evecu_col_dim, evecu, evecu); 
	}


FREE_SPACE:
    if (work != NULL) imsl_free(work);
    if (iwork != NULL) imsl_free(iwork);
    if (acopy != NULL) imsl_free(acopy);
    if (ecopy != NULL) imsl_free(ecopy);
RETURN:
    if (imsl_n1rty(0)>3) {
       if (!evals_user) {
	   if (lv_eval != NULL) imsl_free(lv_eval);
       }
       lv_eval = NULL;  
    }
    return (argptr);
}



/*Translated by FOR_C++, v0.1, on 08/17/90 at 13:48:09 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/17/90 at 13:48:07
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E3LRG/DE3LRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 1,1990

    Purpose:    Compute all of the eigenvalues of a real matrix.

    Usage:      CALL E3LRG (N, A, LDA, EVAL, ACOPY, WK, IWK)

    Arguments:  (See EVLRG)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void     l_e3lrg(Mint *n, Mfloat *a, Mint *lda, Mf_complex eval[],
                        Mfloat *acopy, Mfloat *wk, Mint *iwk)
#else
static void l_e3lrg(n, a, lda, eval, acopy, wk, iwk)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda;
	Mf_complex       eval[];
	Mfloat          *acopy, *wk;
	Mint            *iwk;
#endif
{
#define A(I_,J_)	(a+(I_)*(alda)+(J_))
#define ACOPY(I_,J_)	(acopy+(I_)*(an)+(J_))
#define WK(I_,J_)	(wk+(I_)*(an)+(J_))
#define IWK(I_,J_)	(iwk+(I_)*(an)+(J_))
	Mint alda = *lda;
	Mint an = *n;

	Mint             _l0, i, igh, j, k, low, nprime;


	imsl_e1psh("E3LRG ");
	/* Check N */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "The argument N = %(i1).  The order of the matrix must be at least 1. ");
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
	/* Exit if any errors have occurred. */
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Copy A to ACOPY */
	imsl_crgrg(*n, a, *lda, acopy, *n);
	/* Balance the matrix. */
	l_e9crg(n, &low, &igh, acopy, n, WK(0, 0));
	/*
	 * Reduce to upper Hessenberg using elementary transformations.
	 */
	l_e4lrg(n, n, &low, &igh, acopy, IWK(0, 0));
	/*
	 * Zero out the multipliers stored in the remaining triangle under
	 * the Hessenberg matrix. E5LRG requires a Hessenberg matrix.
	 */
	for (j = 1; j <= (*n - 2); j++) {
		for (i = j + 2; i <= *n; i++) {
			*ACOPY(j - 1, i - 1) = F_ZERO;
		}
	}
	/*
	 * After balancing, the order of the Hessenberg matrix is nprime.
	 * Input acopy(low,low)
	 */
	nprime = igh - low + 1;
	/* Find the eigenvalues. */
	l_e5lrg(n, &nprime, ACOPY(low - 1, low - 1), WK(1, low - 1), WK(2, low - 1));
	/*
	 * Gather the eigenvalues i = 1,low -1 and i = igh+1,n are the real
	 * eigenvalues i = low , igh are the d_complex
	 */
	for (i = 1; i <= (low - 1); i++) {
		*WK(1, i - 1) = *ACOPY(i - 1, i - 1);
		*WK(2, i - 1) = F_ZERO;
		eval[i - 1] = imsl_cf_convert(*WK(1, i - 1), *WK(2, i - 1));
	}
	for (i = igh + 1; i <= *n; i++) {
		*WK(1, i - 1) = *ACOPY(i - 1, i - 1);
		*WK(2, i - 1) = F_ZERO;
		eval[i - 1] = imsl_cf_convert(*WK(1, i - 1), *WK(2, i - 1));
	}
	for (i = low; i <= igh; i++) {
		eval[i - 1] = imsl_cf_convert(*WK(1, i - 1), *WK(2, i - 1));
	}

	/*
	 * Final sort of eigenvalues and eigenvectors.
	 */
	for (i = 1; i <= *n; i++) {
		*IWK(1, i - 1) = i;
		*WK(3, i - 1) = -imsl_c_abs(eval[i - 1]);
	}

	/*
	 * Eigenvalues are sorted into ascending order by magnitude
	 */
	imsl_svrgp(*n, WK(3, 0), WK(3, 0), IWK(1, 0));

	/* Resort the record of the pivots. */
	for (i = 1; i <= *n; i++) {
		for (j = i; j <= *n; j++) {
			if (*IWK(1, j - 1) == i) {
				k = *IWK(1, i - 1);
				*IWK(1, i - 1) = j;
				*IWK(1, j - 1) = k;
				goto L_80;
			}
		}
L_80:
		;
	}

	/*
	 * Take the product of cycles and permute. This will leave the
	 * eigenvalues and eigenvectors sorted in descending order of
	 * magnitude.
	 */
	_l0 = 1;
	for (i = *n - 1; i >= 1; i--) {
		imsl_cswap(&_l0, &eval[i - 1], &_l0, &eval[*IWK(1, i - 1) -
							   1], &_l0);
	}
	/*
	 * Normalize system so that d_complex values have imaginary part
	 * positive, then negative. Sorting may lose this order.
	 */
	j = 1;
L_100:
	;
	if (j >= *n)
		goto L_110;
	if (imsl_c_aimag(eval[j - 1]) != F_ZERO) {
		if (imsl_c_aimag(eval[j - 1]) < F_ZERO) {
			imsl_cswap(&_l0, &eval[j - 1], &_l0, &eval[j],
				   &_l0);
		}
		j += 2;
	} else {
		j += 1;
	}
	goto L_100;
L_110:
	;
L_9000:
	imsl_e1pop("E3LRG ");
	return;
}				/* end of function */

#undef A
#undef ACOPY
#undef WK
#undef IWK




/*Translated by FOR_C++, v0.1, on 08/17/90 at 13:50:06 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/17/90 at 13:50:04
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E4LRG/DE4LRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 1,1990

    Purpose:    Reduce a real general matrix to upper Hessenberg
                form by stabilized elementary similarity transformations.

    Usage:      CALL E4LRG (N, LOW, IGH, A, LDA, ORT, WORK)

    Arguments:
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       N      - Order of the matrix.  (Input)
       LOW    - Lower boundary of the balanced matrix.  (Input)
       IGH    - Upper boundary of the balanced matrix.  (Input)
                 If imsl_e9crg is not used set low=1 and igh=n
       A      - Real N by N matrix to be reduced.  (Input/Output)
                On output, the multipliers which were used in the
                reduction are stored in the remaining triangle under
                the Hessenberg matrix.
       INT    - Integer vector of length IGH containing information on
                the rows and columns interchanged in the
                reduction. Only elements low through igh are used.
                (Output)
    Remark:
       this subroutine is a translation of the algol procedure elmhes,
       num. math. 12, 349-368(1968) by martin and wilkinson.
       handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  ----------------------------------------------------------------------- */
#ifdef ANSI
static void     l_e4lrg(Mint *nm, Mint *n, Mint *low, Mint *igh,
                        Mfloat *a, Mint int_[])
#else
static void l_e4lrg(nm, n, low, igh, a, int_)
	Mint            *nm, *n, *low, *igh;
	Mfloat          *a;
	Mint             int_[];
#endif
{
#define A(I_,J_)	(a+(I_)*(anm)+(J_))
	Mint anm = *nm;

	Mint             i, ip, j, m, mm1, mp1;
	Mfloat           x, y;


	imsl_e1psh("imsl_e4lrg ");
	for (m = *low + 1; m <= (*igh - 1); m++) {
		mm1 = m - 1;
		x = F_ZERO;
		ip = m;

		for (j = m; j <= *igh; j++) {
			if (fabs(*A(mm1 - 1, j - 1)) <= fabs(x))
				goto L_10;
			x = *A(mm1 - 1, j - 1);
			ip = j;
	L_10:
			;
		}

		int_[m - 1] = ip;
		if (x == F_ZERO)
			goto L_100;

		/* interchange columns of a */
		for (j = 1; j <= *igh; j++) {
			y = *A(ip - 1, j - 1);
			*A(ip - 1, j - 1) = *A(m - 1, j - 1);
			*A(m - 1, j - 1) = y;
		}

		y = *A(mm1 - 1, ip - 1);
		*A(mm1 - 1, ip - 1) = *A(mm1 - 1, m - 1);
		*A(mm1 - 1, m - 1) = y;

		mp1 = m + 1;

		for (i = mp1; i <= *igh; i++) {
			*A(mm1 - 1, i - 1) /= x;
		}
		/*
		 * To avoid passing through the address space of A two
		 * times(as would be done if the LHS and RHS updates were
		 * performed one after the other), these updates are
		 * interleaved. Note that the inner loop, with label 50, must
		 * be executed for j = m,n, while the other inner loop, 60
		 * ,is executed only for j=m+1,igh).  To accomplish this, the
		 * outer loop(on j=m,n below) must be broken into three
		 * separate cases: j=m , j=m+1,igh, j=igh+1,n.
		 */
		y = *A(m - 1, ip - 1);
		*A(m - 1, ip - 1) = *A(m - 1, m - 1);
		*A(m - 1, m - 1) = y;
		y = -y;
		for (i = mp1; i <= *igh; i++) {
			*A(m - 1, i - 1) += y ** A(mm1 - 1, i - 1);
		}

		for (j = mp1; j <= *igh; j++) {
			y = *A(j - 1, ip - 1);
			*A(j - 1, ip - 1) = *A(j - 1, m - 1);
			*A(j - 1, m - 1) = y;
			y = -y;
			for (i = mp1; i <= *igh; i++) {
				*A(j - 1, i - 1) += y ** A(mm1 - 1, i - 1);
			}

			/*
			 * compiler directive (ignore recrdeps)
			 */
			for (i = 1; i <= *igh; i++) {
				*A(m - 1, i - 1) += *A(mm1 - 1, j - 1) ** A(j - 1, i - 1);
			}
		}

		for (j = *igh + 1; j <= *n; j++) {
			y = *A(j - 1, ip - 1);
			*A(j - 1, ip - 1) = *A(j - 1, m - 1);
			*A(j - 1, m - 1) = y;
			y = -y;
			for (i = mp1; i <= *igh; i++) {
				*A(j - 1, i - 1) += y ** A(mm1 - 1, i - 1);
			}
		}

L_100:
		;
	}

	imsl_e1pop("imsl_e4lrg ");
	return;
}				/* end of function */


#undef A



/*Translated by FOR_C++, v0.1, on 08/17/90 at 13:51:39 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/17/90 at 13:51:36
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E5LRG/DE5LRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 1,1990

    Purpose:    Determine the eigenvalues of a real upper Hessenberg
                matrix using an heuristically stable combination of
                the (double step,shifted,implicit)LR and QR
                algorithms.

    Usage:      CALL E5LRG (LDA, N , A, WR, WI)

    Arguments:
       LDA    - The leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       N      - Order of the matrix.  (Input)
       A      - Real upper Hessenberg matrix.  (Input/Output)
                 On output a is destroyed.
       WR and WI - are the real and imaginary parts
                   of the eigenvalues. If the algorithm fails to
                   converge after 100*N iterations, only the
                   eigenvalues IQUIT+1, IQUIT+2, ..., N are correct.
                   (Output)

    Remark :    This routine was develped from work done by David
                Watkins and Jeff Haag at Washington State University,
                1989 - 1990.
    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  ----------------------------------------------------------------------- */
#define	ONE	1.0
#define	ZERO	0.0

#ifdef ANSI
static void     l_e5lrg(Mint *lda, Mint *n, Mfloat *a, Mfloat wr[],
                        Mfloat wi[])
#else
static void l_e5lrg(lda, n, a, wr, wi)
	Mint            *lda, *n;
	Mfloat          *a, wr[], wi[];
#endif
{
#define A(I_,J_)	(a+(I_)*(alda)+(J_))
	Mint alda = *lda;
	Mint             i, ibot, ich, imx, inc, itop, its, j, jtop,
	                k, ktop, max_, nits;
	Mfloat           adotw, amx, anorm, big, cpar, det, dif, disc, radical,
	                rdiag, save, sig, small, t1, t2, temp, trace, v[3],
	                w[3], wdota, z;


	/*
	 * ICH chooses the method.
	 * 
	 * ICH .NE. 1  gives the fastest method, the "almost" LR algorithm with
	 * partial pivoting.  Orthogonal transformations (QR) are used only
	 * at the beginning and end of each step.
	 * 
	 * ICH .EQ. 1  gives a hybrid LR-QR algorithm that combines orthogonal
	 * transformations (QR) with Gaussian elimination (LR) without
	 * pivoting. The mix is determined by the "caution parameter" CPAR.
	 * CPAR is a caution parameter that need not be set unless ICH=1. If
	 * ICH=1, CPAR must be set to a value between 0 and 1, inclusive.
	 * Values outside this range are illegal and will be reset.  The
	 * choice CPAR=0 gives the straight QR algorithm, which is the
	 * slowest (and some would say safest) of the options.  The choice
	 * CPAR=1 gives an algorithm that uses an orthogonal transformation
	 * in place of Gaussian elimination only when pivoting would have
	 * been needed for stability.  This is the fastest option under the
	 * choice ICH= Intermediate choices of CPAR give algorithms of
	 * intermediate speed.  The exact use of CPAR is as follows.  At each
	 * step a element p is used to eliminate two other elements q and r.
	 * |p|*CPAR > imsl_i_max{|q|,|r|}, Gaussian elimination is used.
	 * Otherwise, an orthogonal transformation is used.
	 * 
	 * Although other possibilities exist, only ich = 1 and cpar = 1.0 are
	 * used.
	 */
	imsl_e1psh("E5LRG ");
	ich = 1;
	cpar = F_ONE;
	anorm = ZERO;
	for (j = 1; j <= *n; j++) {
		for (i = 1; i <= imsl_i_min(j + 1, *n); i++) {
			anorm += fabs(*A(j - 1, i - 1));
		}
	}
	ibot = *n;
	itop = *n;
	nits = 0;
	/* There are two branches to label 30. */
L_30:
	its = 0;

	/*
	 * Here control is returned to the driving program upon normal
	 * completion
	 */
L_40:
	if (ibot == 0)
		goto L_9000;
	if (ibot == 1) {
		wr[ibot - 1] = *A(ibot - 1, ibot - 1);
		wi[ibot - 1] = ZERO;
		goto L_9000;
	} else if (nits > 100 ** n) {

		/*
		 * The method is not converging rapidly enough. IQUIT is the
		 * error flag. A positive value of IQUIT indicates failure.
		 * IQUIT = 0 on a successful return. Only the eigenvalues
		 * IQUIT+1, IQUIT+2, ..., N are correct.
		 */
		imsl_e1sti(1, 100 ** n);

/*		imsl_ermes(4, 1, "The iteration for the eigenvalues did not converge after %(i1) iterations.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_SLOW_CONVERGENCE_GEN);
		goto L_9000;
	}
	/*
	 * Search for zero on subdiagonal
	 */
L_50:
	big = fabs(*A(itop - 2, itop - 2)) + fabs(*A(itop - 1, itop - 1));
	if (big == ZERO)
		big = anorm;
	small = *A(itop - 2, itop - 1);
	if (big + small != big) {
		itop -= 1;
		if (itop > 1)
			goto L_50;
	} else {
		*A(itop - 2, itop - 1) = ZERO;
	}
	ktop = itop;

	/*
	 * Deflate or isolate irreducible block
	 */
	if (ibot == itop) {

		/*
		 * one (real) eigenvalue
		 */
		wr[ibot - 1] = *A(ibot - 1, ibot - 1);
		wi[ibot - 1] = ZERO;
		ibot -= 1;
		itop = ibot;
		goto L_30;
	} else {
		trace = *A(ibot - 2, ibot - 2) + *A(ibot - 1, ibot - 1);
		if (ibot == itop + 1) {

			/*
			 * two eigenvalues
			 */
			rdiag = *A(ibot - 1, ibot - 2) ** A(ibot - 2, ibot - 1);
			dif = (*A(ibot - 2, ibot - 2) - *A(ibot - 1, ibot - 1)) /
				F_TWO;
			disc = imsl_fi_power(dif, 2) + rdiag;
			radical = sqrt(fabs(disc));
			if (disc >= ZERO) {

				/*
				 * both eigenvalues real
				 */
				radical = dif + sign(radical, dif);
				wr[ibot - 2] = *A(ibot - 1, ibot - 1) + radical;
				wr[ibot - 1] = wr[ibot - 2];
				if (radical != ZERO) {
					wr[ibot - 1] = *A(ibot - 1, ibot - 1) - rdiag / radical;
				}
				wi[ibot - 2] = ZERO;
				wi[ibot - 1] = ZERO;
			} else {

				/*
				 * a d_complex pair of eigenvalues
				 */
				wr[ibot - 2] = *A(ibot - 1, ibot - 1) + dif;
				wr[ibot - 1] = *A(ibot - 1, ibot - 1) + dif;
				wi[ibot - 2] = radical;
				wi[ibot - 1] = -radical;
			}
			ibot -= 2;
			itop = ibot;
			goto L_30;
		} else {

			/*
			 * Search for two consecutive small subdiagonal
			 * entries
			 */
			if (its > 6) {
				jtop = ibot - 2;
				if (jtop > itop) {
			L_60:
					v[0] = ((*A(ibot - 1, ibot - 1) - *A(jtop - 1, jtop - 1)) *
						(*A(ibot - 2, ibot - 2) - *A(jtop - 1, jtop - 1)) -
						*A(ibot - 2, ibot - 1) ** A(ibot - 1, ibot - 2)) /
						*A(jtop - 1, jtop) + *A(jtop, jtop - 1);
					v[1] = *A(jtop - 1, jtop - 1) + *A(jtop, jtop) -
						trace;
					v[2] = *A(jtop, jtop + 1);
					big = v[0] * (fabs(*A(jtop - 2, jtop - 2)) + fabs(*A(jtop - 1, jtop - 1)) +
						      fabs(*A(jtop, jtop)));
					if (big == ZERO)
						big = anorm;
					small = fabs(*A(jtop - 2, jtop - 1)) * (fabs(v[1]) +
								fabs(v[2]));
					if (big + small != big) {
						jtop -= 1;
						/*
						 * This is the only branch to
						 * 60
						 */
						if (jtop > itop)
							goto L_60;
					}
				}
				if (jtop > itop) {
					itop = jtop;
				}
			}
			/*
			 * Perform double GR step on an irreducible block.
			 * The initial bulge is created here by a reflector
			 * (Householder transformation).
			 */
			v[0] = ((*A(ibot - 1, ibot - 1) - *A(itop - 1, itop - 1)) *
				(*A(ibot - 2, ibot - 2) - *A(itop - 1, itop - 1)) - *A(ibot - 2, ibot - 1) *
				*A(ibot - 1, ibot - 2)) / *A(itop - 1, itop) + *A(itop, itop - 1);
			v[1] = *A(itop - 1, itop - 1) + *A(itop, itop) - trace;
			v[2] = *A(itop, itop + 1);

			/*
			 * If necessary, perform an exceptional shift
			 */
			if (mod(its, 10) == 0 && its != 0) {
				t1 = fabs(*A(ibot - 2, ibot - 1)) + fabs(*A(ibot - 3, ibot - 2));
				det = imsl_fi_power(t1, 2);
				trace = 1.5 * t1;
				v[0] = (((*A(itop - 1, itop - 1) - trace) ** A(itop - 1, itop - 1) +
					 det) / *A(itop - 1, itop)) + *A(itop, itop - 1);
				v[1] = *A(itop - 1, itop - 1) + *A(itop, itop) - trace;
			}
			sig = ZERO;
			z = imsl_f_vmax(3, fabs(v[0]), fabs(v[1]), fabs(v[2]));
			if (z != ZERO) {
				for (k = 1; k <= 3; k++) {
					v[k - 1] /= z;
				}
				sig = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
				if (v[0] < ZERO)
					sig = -sig;
				v[0] += sig;
				save = v[0];
				for (k = 1; k <= 3; k++) {
					w[k - 1] = v[k - 1] / save;
					v[k - 1] /= sig;
				}
				sig *= z;

				/*
				 * Multiply on left, (row operations),
				 * forming Q...A...
				 */
				if (ktop != itop) {
					*A(itop - 2, itop - 1) *= ONE - v[0];
				}
				for (j = itop; j <= ibot; j++) {
					wdota = *A(j - 1, itop - 1) + w[1] ** A(j - 1, itop) +
						w[2] ** A(j - 1, itop + 1);
					*A(j - 1, itop - 1) += -wdota * v[0];
					*A(j - 1, itop) += -wdota * v[1];
					*A(j - 1, itop + 1) += -wdota * v[2];
				}

				/*
				 * Multiply on right, (column operations),
				 * forming Q...A...Q to complete similarity
				 * transformation
				 */
				max_ = imsl_i_min(itop + 3, ibot);
				for (j = ktop; j <= max_; j++) {
					adotw = *A(itop - 1, j - 1) + w[1] ** A(itop, j - 1) +
						w[2] ** A(itop + 1, j - 1);
					*A(itop - 1, j - 1) += -adotw * v[0];
					*A(itop, j - 1) += -adotw * v[1];
					*A(itop + 1, j - 1) += -adotw * v[2];
				}
			}
			/*
			 * Now chase bulge, using either reflectors or
			 * Gaussian elimination.
			 */
			if (ibot - itop > 2) {
				for (i = itop; i <= (ibot - 3); i++) {
					if (ich == 1 && cpar * fabs(*A(i - 1, i)) <= imsl_f_max(fabs(*A(i - 1, i + 1)),
						  fabs(*A(i - 1, i + 2)))) {

						/*
						 * Perform an orthogonal
						 * similarity transformation
						 * using a reflector.
						 */
						v[0] = *A(i - 1, i);
						v[1] = *A(i - 1, i + 1);
						v[2] = *A(i - 1, i + 2);
						sig = ZERO;
						z = imsl_f_vmax(3, fabs(v[0]), fabs(v[1]), fabs(v[2]));
						if (z != ZERO) {
							for (k = 1; k <= 3; k++) {
								v[k - 1] /= z;
							}
							sig = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] *
								      v[2]);
							if (v[0] < ZERO)
								sig = -sig;
							v[0] += sig;
							save = v[0];
							for (k = 1; k <= 3; k++) {
								w[k - 1] = v[k - 1] / save;
								v[k - 1] /= sig;
							}
							sig *= z;

							/*
							 * Multiply on left,
							 * forming Q...A
							 */
							for (j = i + 1; j <= ibot; j++) {
								wdota = *A(j - 1, i) + w[1] ** A(j - 1, i + 1) +
									w[2] ** A(j - 1, i + 2);
								*A(j - 1, i) += -wdota * v[0];
								*A(j - 1, i + 1) += -wdota * v[1];
								*A(j - 1, i + 2) += -wdota * v[2];
							}
							*A(i - 1, i) = -sig;
							*A(i - 1, i + 1) = ZERO;
							*A(i - 1, i + 2) = ZERO;

							/*
							 * Multiply on right,
							 * forming Q...A...Q
							 */
							max_ = imsl_i_min(i + 4, ibot);
							for (j = ktop; j <= max_; j++) {
								adotw = *A(i, j - 1) + w[1] ** A(i + 1, j - 1) +
									w[2] ** A(i + 2, j - 1);
								*A(i, j - 1) += -adotw * v[0];
								*A(i + 1, j - 1) += -adotw * v[1];
								*A(i + 2, j - 1) += -adotw * v[2];
							}
						}
					} else {

						/*
						 * Gaussian elimination with
						 * or without partial
						 * pivoting
						 * 
						 * Check whether interchange is
						 * needed.
						 */
						amx = fabs(*A(i - 1, i));
						imx = i + 1;
						if (ich != 1) {
							for (inc = 2; inc <= 3; inc++) {
								if (fabs(*A(i - 1, i + inc - 1)) >
								    amx) {
									amx = fabs(*A(i - 1, i + inc - 1));
									imx = i + inc;
								}
							}
						}
						if (amx != F_ZERO) {

							/*
							 * Interchange if
							 * necessary.
							 */
							if (imx != i + 1) {
								for (j = i; j <= ibot; j++) {
									temp = *A(j - 1, i);
									*A(j - 1, i) = *A(j - 1, imx - 1);
									*A(j - 1, imx - 1) = temp;
								}
								max_ = imsl_i_min(i + 4, ibot);
								for (j = ktop; j <= max_; j++) {
									temp = *A(i, j - 1);
									*A(i, j - 1) = *A(imx - 1, j - 1);
									*A(imx - 1, j - 1) = temp;
								}
							}
							/*
							 * Perform the
							 * elimination
							 */
							t1 = *A(i - 1, i + 1) / *A(i - 1, i);
							t2 = *A(i - 1, i + 2) / *A(i - 1, i);

							/*
							 * Multiply on left,
							 * forming
							 * (L^-1)...A...
							 */
							for (j = i + 1; j <= ibot; j++) {
								*A(j - 1, i + 1) += -t1 ** A(j - 1, i);
								*A(j - 1, i + 2) += -t2 ** A(j - 1, i);
							}
							*A(i - 1, i + 1) = ZERO;
							*A(i - 1, i + 2) = ZERO;

							/*
							 * Multiply on right,
							 * forming
							 * (L^-1)...A...L
							 */
							max_ = imsl_i_min(i + 4, ibot);
							for (j = ktop; j <= max_; j++) {
								*A(i, j - 1) += t1 ** A(i + 1, j - 1) +
									t2 ** A(i + 2, j - 1);
							}
						}
					}
				}
			}
			/*
			 * Reflect away the last little bulge
			 */
			v[0] = *A(ibot - 3, ibot - 2);
			v[1] = *A(ibot - 3, ibot - 1);
			if (v[1] != ZERO) {
				z = fabs(v[0]) + fabs(v[1]);
				for (i = 1; i <= 2; i++) {
					v[i - 1] /= z;
				}
				sig = sqrt(v[0] * v[0] + v[1] * v[1]);
				if (v[0] < ZERO)
					sig = -sig;
				v[0] += sig;
				save = v[0];
				for (i = 1; i <= 2; i++) {
					w[i - 1] = v[i - 1] / save;
					v[i - 1] /= sig;
				}
				sig *= z;

				/*
				 * Multiply on left, forming Q...A...
				 */
				for (j = ibot - 1; j <= ibot; j++) {
					wdota = *A(j - 1, ibot - 2) + w[1] ** A(j - 1, ibot - 1);
					*A(j - 1, ibot - 2) += -wdota * v[0];
					*A(j - 1, ibot - 1) += -wdota * v[1];
				}
				*A(ibot - 3, ibot - 2) = -sig;
				*A(ibot - 3, ibot - 1) = ZERO;

				/*
				 * Multiply on right, forming Q...A...Q
				 */
				for (j = ktop; j <= ibot; j++) {
					adotw = *A(ibot - 2, j - 1) + w[1] ** A(ibot - 1, j - 1);
					*A(ibot - 2, j - 1) += -adotw * v[0];
					*A(ibot - 1, j - 1) += -adotw * v[1];
				}
			}
			/*
			 * Count iterations and repeat the process
			 */
			its += 1;
			nits += 1;
			itop = ibot;
			/* This is the only branch to label 40 */
			goto L_40;
		}
	}
L_9000:
	imsl_e1pop("E5LRG ");
	return;
}				/* end of function */


#undef ONE
#undef ZERO
#undef A






/*Translated by FOR_C++, v0.1, on 08/17/90 at 13:59:29 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/17/90 at 13:59:28
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E6CRG/DE6CRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1984

    Purpose:    Back transform eigenvectors of a matrix transformed
                by E3CRG.

    Usage:      CALL E6CRG (N, LOW, IGH, SCALE, NVEC, EVEC, LDEVEC)

    Arguments:
       N      - Order of the matrix.  (Input)
       LOW    - Lower boundary index of the balanced matrix.  (Input)
       IGH    - Upper boundary index of the balanced matrix.  (Input)
       SCALE  - Real vector of length N containing information about the
                transformation (usually computed by E3CRG).  (Input)
       NVEC   - Number of columns in EVEC to be back transformed.
                (Input)
       EVEC   - Complex N by N matrix containing the eigenvectors to
                be back transformed.  (Input/Output)
       LDEVEC - Leading dimension of EVEC exactly as specified in the
                dimension statement of the calling program.  (Input)

    Remark:
       E6CRG is based on the EISPACK routine BALBAK.

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void     l_e6crg(Mint *n, Mint *low, Mint *igh, Mfloat scale[],
                        Mint *nvec, Mf_complex *evec, Mint *ldevec)
#else
static void l_e6crg(n, low, igh, scale, nvec, evec, ldevec)
	Mint            *n, *low, *igh;
	Mfloat           scale[];
	Mint            *nvec;
	Mf_complex      *evec;
	Mint            *ldevec;
#endif
{
#define EVEC(I_,J_)	(evec+(I_)*(aldevec)+(J_))
	Mint aldevec = *ldevec;
	Mint             i, k;


	imsl_e1psh("E6CRG ");

	if (*nvec == 0)
		goto L_9000;

	for (i = *low; i <= *igh; i++) {
		imsl_csscal(nvec, &scale[i - 1], EVEC(0, i - 1), ldevec);
	}

	for (i = *low - 1; i >= 1; i--) {
		k = scale[i - 1];
		if (k != i) {
			imsl_cswap(nvec, EVEC(0, i - 1), ldevec, EVEC(0, k - 1), ldevec);
		}
	}

	for (i = *igh + 1; i <= *n; i++) {
		k = scale[i - 1];
		if (k != i) {
			imsl_cswap(nvec, EVEC(0, i - 1), ldevec, EVEC(0, k - 1), ldevec);
		}
	}

L_9000:
	imsl_e1pop("E6CRG ");
	return;
}				/* end of function */


#undef EVEC




/*Translated by FOR_C++, v0.1, on 08/17/90 at 13:57:31 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/17/90 at 13:57:29
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E8CRG/DE8CRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 1, 1990

    Purpose:    Compute all of the eigenvalues and eigenvectors of a real
                matrix.

    Usage:      CALL E8CRG (N, A, LDA, EVAL, EVEC, LDEVEC, ACOPY, ECOPY,
                            RWK,IWK)

    Arguments:  (See EVCRG)

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */



#ifdef ANSI
static void     l_e8crg(Mint *n, Mfloat *a, Mint *lda, Mf_complex eval[],
                        Mf_complex *evec, Mint *ldevec, Mfloat *acopy,
                        Mfloat *ecopy, Mfloat *wk, Mint iwk[])
#else
static void l_e8crg(n, a, lda, eval, evec, ldevec, acopy, ecopy,
	   wk, iwk)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda;
	Mf_complex       eval[], *evec;
	Mint            *ldevec;
	Mfloat          *acopy, *ecopy, *wk;
	Mint             iwk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(alda)+(J_))
#define EVEC(I_,J_)	(evec+(I_)*(aldevec)+(J_))
#define ACOPY(I_,J_)	(acopy+(I_)*(an)+(J_))
#define ECOPY(I_,J_)	(ecopy+(I_)*(an)+(J_))
#define WK(I_,J_)	(wk+(I_)*(an)+(J_))
	Mint alda = *lda;
	Mint aldevec = *ldevec;
	Mint an = *n;

	Mint             _l0, i, igh, j, k, low;
	Mfloat           _f0, ri;
	Mf_complex       w;


	imsl_e1psh("E8CRG ");
	/* Check N */
	if (*n < 1) {
		imsl_e1sti(1, *n);

/*		imsl_ermes(5, 1, "The argument N = %(i1).  The order of the matrix must be at least 1. ");
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
	/* Check LDEVEC */
	if (*ldevec < *n) {
		imsl_e1sti(1, *ldevec);
		imsl_e1sti(2, *n);

/*		imsl_ermes(5, 3, "The argument LDEVEC = %(i1).  The leading dimension of the eigenvector matrix must be at least equal to the order of the matrix, N = %(i2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);
	}
	/* Check for errors */
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Copy A to ACOPY */
	imsl_crgrg(*n, a, *lda, acopy, *n);
	/* Balance */
	l_e9crg(n, &low, &igh, acopy, n, WK(0, 0));
	/* Transform to Hessenberg */
	l_eacrg(n, &low, &igh, acopy, n, WK(1, 0));
	l_ebcrg(n, &low, &igh, acopy, n, WK(1, 0), ecopy);
	/*
	 * Find vectors, values of Hessenberg matrix
	 */
	l_e5crh(n, &low, &igh, acopy, n, WK(2, 0), WK(3, 0), ecopy);

	/*
	 * wk(1,3) contains the real components of the eigenvalues wk(1,4)
	 * contains the imaginary components of the eigenvaues
	 */
	for (i = 1; i <= *n; i++) {
		eval[i - 1] = imsl_cf_convert(*WK(2, i - 1), *WK(3, i - 1));
	}
	/* Unpack eigenvectors(originally l_e4crh) */
	j = *n;
L_20:
	;
	if (j < 1)
		goto L_50;
	if (*WK(3, j - 1) == F_ZERO) {
		for (i = 1; i <= *n; i++) {
			*EVEC(j - 1, i - 1) = imsl_cf_convert(*ECOPY(j - 1, i - 1), F_ZERO);
		}
		j -= 1;
	} else {
		for (i = 1; i <= *n; i++) {
			*EVEC(j - 2, i - 1) = imsl_cf_convert(*ECOPY(j - 2, i - 1), *ECOPY(j - 1, i - 1));
			*EVEC(j - 1, i - 1) = imsl_c_conjg(*EVEC(j - 2, i - 1));
		}
		j -= 2;
	}
	goto L_20;
L_50:
	;
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Back transform */
	l_e6crg(n, &low, &igh, WK(0, 0), n, evec, ldevec);
	/*
	 * Sort the eigenvalue magnitudes. Ultimately want eigenvalues in
	 * descending magnitude.
	 */
	for (i = 1; i <= *n; i++) {
		iwk[i - 1] = i;
		*WK(4, i - 1) = -imsl_c_abs(eval[i - 1]);
	}


	imsl_svrgp(*n, WK(4, 0), WK(4, 0), iwk);
	/* Move the eigenvalues */
	for (i = 1; i <= *n; i++) {
		for (j = i; j <= *n; j++) {
			if (iwk[j - 1] == i) {
				k = iwk[i - 1];
				iwk[i - 1] = j;
				iwk[j - 1] = k;
				goto L_80;
			}
		}
L_80:
		;
	}
	/* Take product of cycles and permute */
	_l0 = 1;
	for (i = *n - 1; i >= 1; i--) {
		imsl_cswap(n, EVEC(i - 1, 0), &_l0, EVEC(iwk[i - 1] - 1, 0),
			   &_l0);
		imsl_cswap(&_l0, &eval[i - 1], &_l0, &eval[iwk[i - 1] -
							   1], &_l0);
	}

	/*
	 * Normalize eigenvectors using the Euclidean Norm
	 */
	for (j = 1; j <= *n; j++) {
		ri = imsl_scnrm2(n, EVEC(j - 1, 0), &_l0);
		if (ri > F_ZERO) {
			_f0 = F_ONE/ri;
			imsl_csscal(n, &_f0, EVEC(j - 1, 0), &_l0);
		}
	}
	/*
	 * Normalize system so that d_complex values have imaginary part
	 * positive, then negative. Sorting may lose this order.
	 */
	j = 1;
L_110:
	;
	if (j >= *n)
		goto L_120;
	if (imsl_c_aimag(eval[j - 1]) != F_ZERO) {
		if (imsl_c_aimag(eval[j - 1]) < F_ZERO) {
			imsl_cswap(n, EVEC(j - 1, 0), &_l0, EVEC(j, 0), &_l0);
			imsl_cswap(&_l0, &eval[j - 1], &_l0, &eval[j],
				   &_l0);
		}
		j += 2;
	} else {
		j += 1;
	}
	goto L_110;
L_120:
	;
	/*
	 * Normalize each eigenvector so that its biggest component is real
	 * and positive. The reason for this is that the eigenvectors then
	 * form a right- hand system.
	 */
	for (j = 1; j <= *n; j++) {
		for (i = 1; i <= *n; i++) {
			*WK(5, i - 1) = imsl_fc_convert(imsl_c_mul(*EVEC(j - 1, i - 1), imsl_c_conjg(*EVEC(j - 1, i - 1))));
		}
		i = imsl_isamax(*n, WK(5, 0), 1);
		if (imsl_c_abs(*EVEC(j - 1, i - 1)) == F_ZERO)
			goto L_140;
		/* w = conjugate(evec)/imsl_c_abs(evec) */
		w = imsl_c_div(imsl_c_conjg(*EVEC(j - 1, i - 1)), imsl_cf_convert(imsl_c_abs(*EVEC(j - 1, i - 1)), F_ZERO));
		/* x = x*w */
		imsl_cscal(n, &w, EVEC(j - 1, 0), &_l0);
		*EVEC(j - 1, i - 1) = imsl_cf_convert(imsl_fc_convert(*EVEC(j - 1, i - 1)), F_ZERO);
L_140:
		;
	}

L_9000:
	imsl_e1pop("E8CRG ");
	return;
}				/* end of function */



#undef A
#undef EVEC
#undef ACOPY
#undef ECOPY
#undef WK





/*Translated by FOR_C++, v0.1, on 08/17/90 at 13:54:14 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/17/90 at 13:54:12
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E9CRG/DE9CRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 1, 1990

    Purpose:    Balance a real general matrix.

    Usage:      CALL E9CRG (N, LOW, IGH, A, LDA, SCALE)

    Arguments:
       N      - Order of the matrix.  (Input)
       LOW    - Lower boundary of the balanced matrix.  (Output)
       IGH    - Upper boundary of the balanced matrix.  (Output)
       A      - Real general N by N matrix.  (Input/Output)
       LDA    - Leading dimension of A exactly as specified in the
                the dimension statement of the calling program.  (Input)
       SCALE  - Real vector of length N containing information about
                the similarity transformation.  (Output)

    Remark:
       this subroutine is a translation of the algol procedure balance,
       num. math. 13, 293-304(1969) by parlett and reinsch.
       handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	CONST1	0.95
#define	ONE	1.0
#define	RADIX	16.0
#define	ZERO	0.0

#ifdef ANSI
static void     l_e9crg(Mint *n, Mint *low, Mint *igh, Mfloat *a,
                        Mint *lda, Mfloat scale[])
#else
static void l_e9crg(n, low, igh, a, lda, scale)
	Mint            *n, *low, *igh;
	Mfloat          *a;
	Mint            *lda;
	Mfloat           scale[];
#endif
{
#define A(I_,J_)	(a+(I_)*(alda)+(J_))
	Mint alda = *lda;

/*	LOGICAL32       noconv; */
	Mlong            noconv;
	Mint             i, iexc, j, k, l, m;
	Mfloat           b2, c, f, g, r, s;



	imsl_e1psh("imsl_e9crg ");
	/*
	 * To avoid introducing rounding errors during the balancing process,
	 * the elements of D(diagonal matrices) are restricted to be exact
	 * powers of the radix base employed for floating point arithmetic.
	 */
	b2 = RADIX * RADIX;
	/*
	 * Similarity transformations are used to make corresponding rows and
	 * columns of the matrix being balanced have comparable norms. This
	 * reduces the overall norm of the matrix while leaving the
	 * eigenvalues unchanged.
	 */
	k = 1;
	l = *n;
	/* This is the only branch to 60 */
	goto L_50;
	/*
	 * in-line procedure for row and column exchange The label 10 is
	 * referenced 2 times
	 */
L_10:
	scale[m - 1] = j;
	if (j != m) {

		for (i = 1; i <= l; i++) {
			f = *A(j - 1, i - 1);
			*A(j - 1, i - 1) = *A(m - 1, i - 1);
			*A(m - 1, i - 1) = f;
		}

		for (i = k; i <= *n; i++) {
			f = *A(i - 1, j - 1);
			*A(i - 1, j - 1) = *A(i - 1, m - 1);
			*A(i - 1, m - 1) = f;
		}
	}
	/* This is the only branch to 90 */
	if (iexc == 1) {
		goto L_40;
	} else if (iexc == 2) {
		goto L_80;
	}
	/*
	 * Search for rows isolating an eigenvalue and push them down This is
	 * the only branch to 230
	 */
L_40:
	if (l == 1)
		goto L_200;
	l -= 1;
	/* for j=l step -1 until 1 do -- */
L_50:
	for (j = l; j >= 1; j--) {

		for (i = 1; i <= l; i++) {
			if (i != j) {
				if (*A(i - 1, j - 1) != ZERO)
					goto L_70;
			}
		}

		m = l;
		iexc = 1;
		goto L_10;
L_70:
		;
	}

	goto L_90;
	/*
	 * Search for columns isolating an eigenvalue and push them left
	 */
L_80:
	k += 1;

L_90:
	for (j = k; j <= l; j++) {

		for (i = k; i <= l; i++) {
			if (i != j) {
				if (*A(j - 1, i - 1) != ZERO)
					goto L_110;
			}
		}

		m = k;
		iexc = 2;
		goto L_10;
L_110:
		;
	}
	/*
	 * Now balance the submatrix in rows K to L
	 */
	for (i = k; i <= l; i++) {
		scale[i - 1] = ONE;
	}
	/* Iterative loop for norm reduction */
L_130:
	noconv = FALSE;

	for (i = k; i <= l; i++) {
		c = ZERO;
		r = ZERO;

		for (j = k; j <= l; j++) {
			if (j != i) {
				c += fabs(*A(i - 1, j - 1));
				r += fabs(*A(j - 1, i - 1));
			}
		}
		/*
		 * Guard against zero C or R due to underflow
		 */
		if (c == ZERO || r == ZERO)
			goto L_190;
		g = r / RADIX;
		f = ONE;
		s = c + r;
L_150:
		if (c < g) {
			f *= RADIX;
			c *= b2;
			goto L_150;
		}
		g = r * RADIX;
L_160:
		if (c >= g) {
			f /= RADIX;
			c /= b2;
			goto L_160;
		}
		/* Now balance */
		if ((c + r) / f < CONST1 * s) {
			g = ONE / f;
			scale[i - 1] *= f;
			noconv = TRUE;

			for (j = k; j <= *n; j++) {
				*A(j - 1, i - 1) *= g;
			}

			for (j = 1; j <= l; j++) {
				*A(i - 1, j - 1) *= f;
			}
		}
L_190:
		;
	}

	if (noconv)
		goto L_130;

L_200:
	*low = k;
	*igh = l;
	/*
	 * Output is a matrix that is balanced in the norm given by summing
	 * the absolute magnitudes of the matrix elements.
	 */
	imsl_e1pop("imsl_e9crg ");
	return;
}				/* end of function */


#undef CONST1
#undef ONE
#undef RADIX
#undef ZERO
#undef A



/*Translated by FOR_C++, v0.1, on 08/17/90 at 14:00:27 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/17/90 at 14:00:24
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  EACRG/DEACRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 1, 1990
    Purpose:    Reduce a real general matrix to an upper Hessenberg
                matrix by orthogonal similarity transformations.

    Usage:      CALL EACRG (N, LOW, IGH, A, LDA, ORT)

    Arguments:
       N      - Order of the matrix.  (Input)
       LOW    - Lower boundary of the balanced matrix.  (Input)
       IGH    - Upper boundary of the balanced matrix.  (Input)
       A      - Real N by N matrix to be reduced.  (Input/Output)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       ORT    - Real vector of length IGH containing information about
                the transformation. Only elements low through igh
                are used. (Output)

    Remark:
       this subroutine is a translation of the algol procedure orthes,
       num. math. 12, 349-368(1968) by martin and wilkinson.
       handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	ONE	1.0
#define	ZERO	0.0

#ifdef ANSI
static void     l_eacrg(Mint *n, Mint *low, Mint *igh, Mfloat *a,
                        Mint *lda, Mfloat ort[])
#else
static void l_eacrg(n, low, igh, a, lda, ort)
	Mint            *n, *low, *igh;
	Mfloat          *a;
	Mint            *lda;
	Mfloat           ort[];
#endif
{
#define A(I_,J_)	(a+(I_)*(alda)+(J_))
	Mint alda = *lda;
	Mint             i, j, m, mp;
	Mfloat           f, g, h, scale, small, big ;


	imsl_e1psh("imsl_eacrg");
        /*
         * Small and Big are used to compute a bound on the smallest number
         * that can be reciprocated.
         */
        small = imsl_amach(1);
        big   = imsl_amach(2);
        if (small*big <= 1.0) 
            small = 1.0/big;

	for (m = *low + 1; m <= (*igh - 1); m++) {
		h = ZERO;
		ort[m - 1] = ZERO;
		scale = ZERO;
		/*
		 * Scale column (ALGOL TOL then not needed)
		 */
		for (i = m; i <= *igh; i++) {
			scale += fabs(*A(m - 2, i - 1));
		}

		if (scale > small){
		    
	

		scale = ONE / scale;
		/*
		 * for i=igh step -1 until m do -- This has been reversed
		 */
		for (i = m; i <= *igh; i++) {
			ort[i - 1] = *A(m - 2, i - 1) * scale;
			h += ort[i - 1] * ort[i - 1];
		}

		g = -sign(sqrt(h), ort[m - 1]);
		h = ONE / (h - ort[m - 1] * g);
		ort[m - 1] -= g;
		/*
		 * Form (I-(U*UT)/H) * A
		 * 
		 * Compiler directive(prefer vector)
		 */
		for (j = m; j <= *n; j++) {
			f = ZERO;

			for (i = m; i <= *igh; i++) {
				f += ort[i - 1] ** A(j - 1, i - 1);
			}

			f *= -h;

			for (i = m; i <= *igh; i++) {
				*A(j - 1, i - 1) += f * ort[i - 1];
			}

		}
		/*
		 * Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
		 * 
		 * Compiler directive(prefer vector)
		 */
		for (i = 1; i <= *igh; i++) {
			f = ZERO;
			/*
			 * for j=igh step -1 until m do -- This has been
			 * reversed
			 */
			for (j = m; j <= *igh; j++) {
				f += ort[j - 1] ** A(j - 1, i - 1);
			}

			f *= -h;

			for (j = m; j <= *igh; j++) {
				*A(j - 1, i - 1) += f * ort[j - 1];
			}

		}

		ort[m - 1] /= scale;
		*A(m - 2, m - 1) = g / scale;

                } else {
                        goto L_90;
                }

L_90:       
		;
	}

	imsl_e1pop("imsl_eacrg");
	return;
}				/* end of function */

#undef ONE
#undef ZERO
#undef A








/*Translated by FOR_C++, v0.1, on 08/17/90 at 14:01:08 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/17/90 at 14:01:07
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  EBCRG/DEBCRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    June 1, 1990

    Purpose:    Accumulate the transformations in reducing a real general
                matrix to an upper Hessenberg matrix by orthogonal
                similarity transformations.

    Usage:      CALL EBCRG (N, LOW, IGH, A, LDA, ORT, EVEC)

    Arguments:
       N      - Order of the matrix.  (Input)
       LOW    - Lower boundary index of the balanced matrix.  (Input)
       IGH    - Upper boundary index of the balanced matrix.  (Input)
       A      - Real N by N matrix.  (Input)
                The strict lower triangular part contains information
                about the transformation.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       ORT    - Real vector of length IGH containing information about
                the transformation.  (Input)
       EVEC   - Real N by N matrix containing the transformation matrix.
                (Output)

    Remark:
       this subroutine is a translation of the algol procedure ortrans,
       num. math. 16, 181-204(1970) by peters and wilkinson.
       handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	ONE	1.0
#define	ZERO	0.0

#ifdef ANSI
static void     l_ebcrg(Mint *n, Mint *low, Mint *igh, Mfloat *a,
                        Mint *lda, Mfloat ort[], Mfloat *evec)
#else
static void l_ebcrg(n, low, igh, a, lda, ort, evec)
	Mint            *n, *low, *igh;
	Mfloat          *a;
	Mint            *lda;
	Mfloat           ort[], *evec;
#endif
{
#define A(I_,J_)	(a+(I_)*(alda)+(J_))
#define EVEC(I_,J_)	(evec+(I_)*(alda)+(J_))
	Mint alda = *lda;

	Mint             i, j, mp;
	Mfloat           g, small, big;



	imsl_e1psh("imsl_ebcrg ");
        /*
         * Small and Big are used to compute a bound on the smallest number
         * that can be reciprocated.
         */
        small = imsl_amach(1);
        big = imsl_amach(2);
        if (small * big < 1.0)
                small = 1.0 / big;

	for (j = 1; j <= *n; j++) {

		/* Initialize EVEC to identity matrix */
		for (i = 1; i <= *n; i++) {
			*EVEC(j - 1, i - 1) = ZERO;
		}

		*EVEC(j - 1, j - 1) = ONE;
	}

	/* for mp=igh-1 step -1 until low+1 do - */
	for (mp = *igh - 1; mp >= (*low + 1); mp--) {
		/* This is the only branch to 70 */
		if ( fabs(*A(mp - 2, mp - 1)) > small){
		    
		for (i = mp + 1; i <= *igh; i++) {
			ort[i - 1] = *A(mp - 2, i - 1);
		}
		/*
		 * Compiler directive(prefer vector)
		 */
		for (j = mp; j <= *igh; j++) {
			g = ZERO;

			for (i = mp; i <= *igh; i++) {
				g += ort[i - 1] ** EVEC(j - 1, i - 1);
			}
			/*
			 * Divisor below is negative of H formed in
			 * imsl_eacrg.  Double division avoids possible
			 * underflow.
			 */
			g = (g / ort[mp - 1]) / *A(mp - 2, mp - 1);

			for (i = mp; i <= *igh; i++) {
				*EVEC(j - 1, i - 1) += g * ort[i - 1];
			}

		}


                } else {
                        goto L_70;
                }

L_70:
		;
	}

	imsl_e1pop("imsl_ebcrg ");
	return;
}				/* end of function */



#undef ONE
#undef ZERO
#undef A
#undef EVEC




/*Translated by FOR_C++, v0.1, on 08/17/90 at 15:58:37 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/17/90 at 15:58:36
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  ECCRG/DECCRG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Purpose:    COMPLEX DIVIDE

    Usage:      CALL ECCRG (AR,AI,BR,BI,CR,CI)

    Arguments:
       AR    - Real part of the numerator.  (Input)
       AI    - Imaginary part of the numerator.  (Input)
       BR    - Real part of the denominator.  (Input)
       BI    - Imaginary part of the denominator. ( Input)
       CR    - Real part of the solution.(Output)
       CI    - Imaginary part of the solution.(Output)

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  ----------------------------------------------------------------------- */
#ifdef ANSI
static void     l_eccrg(Mfloat *ar, Mfloat *imsl_ai, Mfloat *br,
                        Mfloat *imsl_bi, Mfloat *cr, Mfloat *imsl_ci)
#else
static void l_eccrg(ar, imsl_ai, br, imsl_bi, cr, imsl_ci)
	Mfloat          *ar, *imsl_ai, *br, *imsl_bi, *cr, *imsl_ci;
#endif
{
	Mfloat           ais, ars, bis, brs, s;


	s = fabs(*br) + fabs(*imsl_bi);
	ars = *ar / s;
	ais = *imsl_ai / s;
	brs = *br / s;
	bis = *imsl_bi / s;
	s = imsl_fi_power(brs, 2) + imsl_fi_power(bis, 2);
	*cr = (ars * brs + ais * bis) / s;
	*imsl_ci = (ais * brs - ars * bis) / s;
	return;
}				/* end of function */
/*Translated by FOR_C++, v0.1, on 08/05/91 at 17:20:18 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/05/91 at 17:20:13
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  E5CRH/DE5CRH (Single/Double precision version)

    Computer:   SUN4/SINGLE

    Revised:    June 1, 1990

    Purpose:    Determines the eigenvalues and eigenvectors
                of a real upper Hessenberg matrix by the QR method.
                The eigenvectors of a real general matrix can also
                be found if imsl_eacrg  and  imsl_ebcrg  have been used to reduce
                this general matrix to Hessenberg form and to
                accumulate the similarity transformations.

    Usage:      CALL E5CRH (N, LOW, IGH, A, LDA, WR,WI, ECOPY)

       N      - Order of the matrix.  (Input)
       LOW    - Index of the first location in the balanced matrix.
                (Input)
       IGH    - Index of the last location in the balanced matrix.
                (Input)
       A      - Real upper Hessenberg matrix.  (Input)

       LDA    - The leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)

       WR and WI - contain the real and imaginary parts,
                   respectively, of the eigenvalues.  The eigenvalues
                   are unordered except that d_complex conjugate imsl_pairs
                   of values appear consecutively with the eigenvalue
                   having the positive imaginary part first. Each of
                   WR and WI are real vectors of length N. (Output)
       ECOPY   -  On input ecopy contains the transformation matrix.
                  On output ecopy contains the real and imaginary
                  parts of the eigenvectors. If the i-th eigenvalue is
                  real, the i-th column of ecopy contains its eigenvector
                  If the i-th eigenvalue is d_complex with positive
                  imaginary part, the i-th and (i+1)-th columns of ecopy
                  contain the real and imaginary parts of its eigenvector
                  If an error exit is made, none of the eigenvectors
                  have been found. Ecopy is a real N*N matrix.
                  ( Input/Output)

     Remarks:
       calls imsl_eccrg for d_complex division.
       this subroutine is a translation of the algol procedure hqr2,
       num. math. 16, 181-204(1970) by peters and wilkinson.
       handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).

    Chapter:    MATH/LIBRARY Eigensystem Analysis

    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

       this subroutine is a translation of the algol procedure hqr2,
       num. math. 16, 181-204(1970) by peters and wilkinson.
       handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).

  -----------------------------------------------------------------------
 */
#define	CONST	0.75
#define	CONST1	(-0.4375)
#define	CONST2	0.01
#define	ONE	1.0
#define	TWO	2.0
#define	ZERO	0.0

#ifdef ANSI
static void l_e5crh (Mint *n, Mint *low, Mint *igh, Mfloat *a,
		Mint *lda, Mfloat wr[], Mfloat wi[], Mfloat *ecopy)
#else
static void l_e5crh (n, low, igh, a, lda, wr, wi, ecopy)
    Mint        *n, *low, *igh;
    Mfloat      *a;
    Mint        *lda;
    Mfloat       wr[], wi[], *ecopy;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define ECOPY(I_,J_)	(ecopy+(I_)*(*lda)+(J_))
    Mint   notlas, usebla;
    Mint         _l0, _l1, enm2, i, ien, itn, its, j, k, l, m, na;
    Mfloat       _f0, _f1, _f2, big, cc, norm, p, p1, q, r, ra, s, s1, s2,
                sa, small, ss, t, t1, tol, tst1, tst2, u, vi, vr, w, x,
                x1, y, zz;


    imsl_e1psh ("imsl_e5crh ");
    /*
     * Small and Big are used to compute a bound on the smallest number that
     * can be reciprocated.
     */
    small = imsl_amach (1);
    big = imsl_amach (2);
    if (small * big < 1.0)
	small = 1.0 / big;
    tol = imsl_amach (4);
    norm = ZERO;
    /*
     * All constants have been put in a parameter statement for efficiency.
     * 
     * The QR transformation(with shifts) preserves the  upper Hessenberg form
     * of the original matrix. Shifting is neccessary for rapid convergence.
     * Ultimately, the eigenvalues are eithe isolated on the diagonal or are
     * eigenvalues of a 2 X 2 submatrix on the diagonal. Note that a
     * nonsymmetric real matrix can have d_complex eigenvalues.
     * 
     * Compute matrix norm
     */

        for (j = 1; j <= *n; j++) {
                for (i = 1; i <= imsl_i_min(*n, j + 1); i++) {
                        norm += fabs(*A(j - 1, i - 1));
                }
        }

        /*
         * If imsl_e9crg is not used, set low = 1 and igh = n.
         * 
         * 
         * Store eigenvalues isolated by imsl_e9crg
         */
        for (i = 1; i <= (*low - 1); i++) {
                wr[i - 1] = *A(i - 1, i - 1);
                wi[i - 1] = ZERO;
        }

        for (i = *igh + 1; i <= *n; i++) {
                wr[i - 1] = *A(i - 1, i - 1);
                wi[i - 1] = ZERO;
        }

        ien = *igh;
        /*
         * T will only get changed by an exceptional shift
         */
        t = ZERO;
        /* Search for next eigenvalue */
L_50:
        if (ien < *low)
                goto L_330;
        itn = 50;
        its = 1;
        na = ien - 1;
        enm2 = na - 1;
        /*
         * Begin iteration look for single small sub-diagonal element
         */
L_60:
        for (l = ien; l >= (*low + 1); l--) {
                /* Local scaling. */
                s = fabs(*A(l - 2, l - 2)) + fabs(*A(l - 1, l - 1));
                /*
                 * Relax scale if no convergence. Reset iteration counter for
                 * eigenvalue
                 */
                if (itn == 0) {
                        s = norm;
                }
                if (fabs(*A(l - 2, l - 1)) <= tol * s)
                        goto L_80;
        }
        /* Form shift */
L_80:
        x = *A(ien - 1, ien - 1);

        if (l == ien)
                goto L_270;
        y = *A(na - 1, na - 1);
        w = *A(na - 1, ien - 1) ** A(ien - 1, na - 1);

        if (l == na)
                goto L_280;
        if (itn == 0) {
                imsl_e1sti(1, 50);


                imsl_e1mes(4, 1, "The iteration for the eigenvalues did not converge after %(i1) iterations.");

                goto L_9000;
        }
        if (mod(its, 10) != 0)
                goto L_100;
        /* Form exceptional shift */
        t += x;

        for (i = *low; i <= ien; i++) {
                *A(i - 1, i - 1) -= x;
        }

        s = fabs(*A(na - 1, ien - 1)) + fabs(*A(enm2 - 1, na - 1));
        x = CONST * s;
        y = x;
        w = CONST1 * s * s;
L_100:
        its += 1;
        itn -= 1;
        /*
         * The shifts of origin at each stage are taken to be the eigenvalues
         * of the 2 X 2 matrix in the bottom right-hand corner of the current
         * matrix.
         * 
         * Form shift and then look for two consecutive small sub-diagonal
         * elements in the submatrix
         */
        for (m = ien - 2; m >= l; m--) {
                zz = *A(m - 1, m - 1);
                r = x - zz;
                s = y - zz;
                /*
                 * orig        p = (r*s-w)/a(m+1,m) + a(m,m+1) orig        q
                 * = a(m+1,m+1) - zz - r - s orig        r = a(m+2,m+1)
                 * Commented code was changed to the following because of
                 * overflow problems on the VAX.
                 */
                p = r * s - w;
                q = *A(m, m) - zz - r - s;
                r = *A(m, m + 1);
                if (fabs(*A(m - 1, m)) > 1.0e0) {
                        p = p / *A(m - 1, m) + *A(m, m - 1);
                } else {
                        p += *A(m, m - 1) ** A(m - 1, m);
                        q *= *A(m - 1, m);
                        r *= *A(m - 1, m);
                }
                if (m == l)
                        goto L_120;
                tst1 = fabs(p) * (fabs(*A(m - 2, m - 2)) + fabs(zz) + fabs(*A(m, m)));
                tst2 = fabs(*A(m - 2, m - 1)) * (fabs(q) + fabs(r));
                if (tst2 <= tol * tst1) {
                        s = fabs(p) + fabs(q) + fabs(r);
                        if (s > small) {
                                p /= s;
                                q /= s;
                                r /= s;
                                goto L_120;
                        }
                }
        }


L_120:
        *A(m - 1, m + 1) = ZERO;

        for (i = m + 3; i <= ien; i++) {
                *A(i - 3, i - 1) = ZERO;
                *A(i - 4, i - 1) = ZERO;
        }

        /*
         * Double QR step involving rows L to EN and columns M to EN
         */
        for (k = m; k <= na; k++) {
                notlas = k != na;
                if (k != m) {
                        /* Begin setup of Householder vector */
                        p = *A(k - 2, k - 1);
                        q = *A(k - 2, k);
                        r = ZERO;
                        if (notlas)
                                r = *A(k - 2, k + 1);
                        x = fabs(p) + fabs(q) + fabs(r);

                        if (x > small) {
                                /*
                                 * x1 is used to avoid multiple divides
                                 */
                                x1 = ONE / x;
                                p *= x1;
                                q *= x1;
                                r *= x1;
                                s = sign(sqrt(p * p + q * q + r * r), p);
                                s2 = ONE / s;
                                if (k == m)
                                        goto L_140;
                                *A(k - 2, k - 1) = -s * x;
                                goto L_150;
                L_140:
                                if (l != m)
                                        *A(k - 2, k - 1) = -*A(k - 2, k - 1);
                L_150:
                                p += s;
                                p1 = ONE / p;
                                /*
                                 * s2 and p1 are used to avoid multiple
                                 * divides
                                 */
                                x = p * s2;
                                y = q * s2;
                                zz = r * s2;
                                q *= p1;
                                r *= p1;

                                if (notlas)
                                        goto L_210;
                        } else {
                                goto L_260;
                        }
                } else {
                        s = sign(sqrt(p * p + q * q + r * r), p);
                        /*
                         * test s so not to divide by zero later.
                         */
                        if (fabs(s) > small) {
                                if (k == m)
                                        goto L_160;
                                *A(k - 2, k - 1) = -s * x;
                                goto L_170;
                L_160:
                                if (l != m)
                                        *A(k - 2, k - 1) = -*A(k - 2, k - 1);
                L_170:
                                p += s;
                                x = p / s;
                                y = q / s;
                                zz = r / s;
                                q /= p;
                                r /= p;

                                if (notlas)
                                        goto L_210;
                        } else {
                                goto L_260;
                        }
                }
                /*
                 * Row modification The logical flag below is arbitrarily set
                 * as shown. The BLAS are used to apply transformations only
                 * when all dimensions are .gt. 16.
                 */
                usebla = *n - k > 16 && k > 16;
                if (usebla) {
                        cc = ONE - x;
                        ss = -y;
                        _l0 = *n - k + 1;
                        imsl_srot(_l0, A(k - 1, k - 1), *lda, A(k - 1, k),
                                  *lda, cc, ss);
                        _l0 = imsl_i_min (ien, k + 3);
                        _l1 = 1;
                        imsl_srot(_l0,A(k - 1, 0), _l1,
                                  A(k, 0), _l1, cc, ss);
                        _l0 = *igh - *low + 1;
                        imsl_srot(_l0,  ECOPY(k - 1, *low - 1),
                                  _l1, ECOPY(k, *low - 1), _l1, cc, ss);
                } else {

                        for (j = k; j <= *n; j++) {
                                u = *A(j - 1, k - 1) + q ** A(j - 1, k);
                                *A(j - 1, k - 1) += -u * x;
                                *A(j - 1, k) += -u * y;
                        }

                        j = imsl_i_min(ien, k + 3);
                        /* Column modification */
                        for (i = 1; i <= j; i++) {
                                u = (x ** A(k - 1, i - 1)) + (y ** A(k, i - 1));                                *A(k - 1, i - 1) -= u;
                                *A(k, i - 1) += -u * q;
                        }
                        /* Accumulate transformations */
                        for (i = *low; i <= *igh; i++) {
                                u = (x ** ECOPY(k - 1, i - 1)) + (y ** ECOPY(k, i - 1));
                                *ECOPY(k - 1, i - 1) -= u;
                                *ECOPY(k, i - 1) += -u * q;
                        }
                }
                goto L_250;
L_210:
                ;
                /*
                 * Row modification
                 * 
                 */
                for (j = k; j <= *n; j++) {
                        u = *A(j - 1, k - 1) + q ** A(j - 1, k) + r ** A(j - 1, k + 1);
                        *A(j - 1, k - 1) += -u * x;
                        *A(j - 1, k) += -u * y;
                        *A(j - 1, k + 1) += -u * zz;
                }

                j = imsl_i_min(ien, k + 3);
                /* Column modification */
                for (i = 1; i <= j; i++) {
                        u = (x ** A(k - 1, i - 1)) + (y ** A(k, i - 1)) + (zz ** A(k + 1, i - 1));
                        *A(k - 1, i - 1) -= u;
                        *A(k, i - 1) += -u * q;
                        *A(k + 1, i - 1) += -u * r;
                }
                /* Accumulate transformations */
                for (i = *low; i <= *igh; i++) {
                        u = (x ** ECOPY(k - 1, i - 1)) + (y ** ECOPY(k, i - 1)) + (zz *
                                                      *ECOPY(k + 1, i - 1));
                        *ECOPY(k - 1, i - 1) -= u;
                        *ECOPY(k, i - 1) += -u * q;
                        *ECOPY(k + 1, i - 1) += -u * r;
                }
L_250:
                ;

L_260:
                ;

        }

        /*
         * for next iteration on current eigenvalue
         */
        goto L_60;
        /* One root found */
L_270:
        *A(ien - 1, ien - 1) = x + t;
        wr[ien - 1] = *A(ien - 1, ien - 1);
        wi[ien - 1] = ZERO;
        ien = na;

        goto L_50;
        /* Two roots found */
L_280:
        p = (y - x) / TWO;
        q = p * p + w;
        zz = sqrt(fabs(q));
        *A(ien - 1, ien - 1) = x + t;
        x = *A(ien - 1, ien - 1);
        *A(na - 1, na - 1) = y + t;
        if (q >= ZERO) {
                /* Real pair */
                zz = p + sign(zz, p);
                wr[na - 1] = x + zz;
                wr[ien - 1] = wr[na - 1];
                if (zz != ZERO)
                        wr[ien - 1] = x - w / zz;
                wi[na - 1] = ZERO;
                wi[ien - 1] = ZERO;
                x = *A(na - 1, ien - 1);
                s = fabs(x) + fabs(zz);
                p = x / s;
                q = zz / s;
                r = sqrt(p * p + q * q);
                p /= r;
                q /= r;
                /*
                 * Row modification
                 */
                for (j = na; j <= *n; j++) {
                        u = *A(j - 1, na - 1);
                        *A(j - 1, na - 1) = q * u + p ** A(j - 1, ien - 1);
                        *A(j - 1, ien - 1) = q ** A(j - 1, ien - 1) - p * u;
                }
                /* Column modification */
                for (i = 1; i <= ien; i++) {
                        u = *A(na - 1, i - 1);
                        *A(na - 1, i - 1) = q * u + p ** A(ien - 1, i - 1);
                        *A(ien - 1, i - 1) = q ** A(ien - 1, i - 1) - p * u;
                }
                /* Accumulate transformations */
                for (i = *low; i <= *igh; i++) {
                        u = *ECOPY(na - 1, i - 1);
                        *ECOPY(na - 1, i - 1) = q * u + p ** ECOPY(ien - 1, i - 1);
                        *ECOPY(ien - 1, i - 1) = q ** ECOPY(ien - 1, i - 1) - p * u;
                }
                goto L_320;
        } else {
                /* Complex pair */
                wr[na - 1] = x + p;
                wr[ien - 1] = x + p;
                wi[na - 1] = zz;
                wi[ien - 1] = -zz;
                goto L_320;
        }
L_320:
        ien = enm2;
        goto L_50;
        /*
         * All roots found. Backsubstitute to find vectors of upper
         * triangular form Only
         */
L_330:
        if (norm == ZERO)
                goto L_9000;
        /* for ien=n step -1 until 1 do -- */
        for (ien = *n; ien >= 1; ien--) {
                p = wr[ien - 1];
                q = wi[ien - 1];
                na = ien - 1;
                /* Aritmetic if was taken out */
                if (q == ZERO) {
                        /* Real vector */
        L_340:
                        m = ien;
                        *A(ien - 1, ien - 1) = ONE;
                        if (na == 0)
                                goto L_560;

                        for (i = ien - 1; i >= 1; i--) {
                                w = *A(i - 1, i - 1) - p;
                                r = ZERO;

                                for (j = m; j <= ien; j++) {
                                        r += *A(j - 1, i - 1) ** A(ien - 1, j - 1);
                                }
                                if (wi[i - 1] >= ZERO)
                                        goto L_360;
                                zz = w;
                                s = r;
                                goto L_430;
                L_360:
                                m = i;
                                if (wi[i - 1] != ZERO)
                                        goto L_390;
                                t = w;
                                if (t != ZERO)
                                        goto L_380;
                                tst1 = norm;
                                t = tst1;
                L_370:
                                t *= CONST2;
                                tst2 = norm + t;
                                if (tst2 > tst1)
                                        goto L_370;
                L_380:
                                *A(ien - 1, i - 1) = -r / t;
                                goto L_410;
                                /* Solve real equations */
                L_390:
                                x = *A(i, i - 1);
                                y = *A(i - 1, i);
                                q = (wr[i - 1] - p) * (wr[i - 1] - p) + wi[i - 1] * wi[i - 1];
                                t = (x * s - zz * r) / q;
                                *A(ien - 1, i - 1) = t;
                                if (fabs(x) <= fabs(zz))
                                        goto L_400;
                                *A(ien - 1, i) = (-r - w * t) / x;
                                goto L_410;
                L_400:
                                *A(ien - 1, i) = (-s - y * t) / zz;

                                /* Overflow control */
                L_410:
                                t = fabs(*A(ien - 1, i - 1));
                                if (t > small) {
                                        tst1 = t;
                                        tst2 = tst1 + ONE / tst1;
                                        if (tst2 > tst1)
                                                goto L_430;
                                        t1 = ONE / t;
                                        /*
                                         * t1 is used to avoid multiple
                                         * divides
                                         */
                                        for (j = i; j <= ien; j++) {
                                                *A(ien - 1, j - 1) *= t1;
                                        }
                                } else {
                                        goto L_430;
                                }
                L_430:
                                ;
                        }
                        /* End real vector */
                        goto L_560;
                        /* Complex vector */
                } else if (q < 0) {
        L_440:
                        m = na;
                        /*
                         * Last vector component chosen imaginary so that
                         * eigenvector matrix is triangular
                         */
                        if (fabs(*A(na - 1, ien - 1)) <= fabs(*A(ien - 1, na - 1)))
                                goto L_450;
                        *A(na - 1, na - 1) = q / *A(na - 1, ien - 1);
                        *A(ien - 1, na - 1) = -(*A(ien - 1, ien - 1) - p) / *A(na - 1, ien - 1);

                        goto L_460;
                        /* Call to a d_complex divide routine */
        L_450:
                        _f0 = 0.0;
                        _f1 = -*A (ien - 1, na - 1);
                        _f2 = *A (na - 1, na - 1) - p;
                        l_eccrg(&_f0, &_f1, &_f2, &q, A(na - 1, na - 1), A(ien - 1, na - 1));
        L_460:
                        *A(na - 1, ien - 1) = ZERO;
                        *A(ien - 1, ien - 1) = ONE;
                        enm2 = na - 1;
                        if (enm2 == 0)
                                goto L_560;

                        for (i = ien - 2; i >= 1; i--) {
                                w = *A(i - 1, i - 1) - p;
                                ra = ZERO;
                                sa = ZERO;

                                for (j = m; j <= ien; j++) {
                                        ra += *A(j - 1, i - 1) ** A(na - 1, j - 1);
                                        sa += *A(j - 1, i - 1) ** A(ien - 1, j - 1);
                                }

                                if (wi[i - 1] >= ZERO)
                                        goto L_480;
                                zz = w;
                                r = ra;
                                s = sa;

                                goto L_550;
                L_480:
                                m = i;
                                if (wi[i - 1] != ZERO)
                                        goto L_490;
                                /* Call to d_complex divide routine */
                                _f0 = -ra;
                                _f1 = -sa;
                                l_eccrg(&_f0, &_f1, &w, &q, A(na - 1, i - 1),
                                           A(ien - 1, i - 1));

                                goto L_530;
                                /* Solve d_complex equations */
                L_490:
                                x = *A(i, i - 1);
                                y = *A(i - 1, i);
                                vr = (wr[i - 1] - p) * (wr[i - 1] - p) + wi[i - 1] * wi[i - 1] -
                                        q * q;
                                vi = (wr[i - 1] - p) * TWO * q;
                                if (vr != ZERO || vi != ZERO)
                                        goto L_510;
                                tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) +
                                               fabs(zz));
                                vr = tst1;

                L_500:
                                vr *= CONST2;
                                tst2 = tst1 + vr;

                                if (tst2 > tst1)
                                        goto L_500;
                                /* Call to d_complex divide routine */
                L_510:
                                _f0 = x * r - zz * ra + q * sa;
                                 _f1 = x * s - zz * sa - q * ra;

                                l_eccrg(&_f0,  &_f1, &vr, &vi, A(na - 1, i - 1), A(ien - 1, i - 1));
                                if (fabs(x) <= fabs(zz) + fabs(q))
                                        goto L_520;
                                *A(na - 1, i) = (-ra - w ** A(na - 1, i - 1) + q ** A(ien - 1, i - 1)) /
                                        x;
                                *A(ien - 1, i) = (-sa - w ** A(ien - 1, i - 1) - q ** A(na - 1, i - 1)) /
                                        x;

                                goto L_530;
                                /* Call to d_complex divide routine */
                L_520:
                                _f0 = -r - y ** A (na - 1, i - 1);
                                _f1 = -s - y ** A (ien - 1, i - 1);
                                l_eccrg(&_f0,  &_f1, &zz, &q, A(na - 1, i), A(ien - 1, i));
                                /* Overflow control */
                L_530:
                                t = imsl_f_max(fabs(*A(na - 1, i - 1)), fabs(*A(ien - 1, i - 1)));

                                /*
                                 * if (t .ne. zero) then Commented code was
                                 * replaced because of potential underflow/
                                 * overflow errors.
                                 */
                                if (t > small) {
                                        t1 = ONE / t;
                                        tst1 = t;
                                        tst2 = tst1 + t1;
                                        if (tst2 <= tst1) {
                                                /*
                                                 * t1 is used to avoid
                                                 * multiple divides
                                                 */
                                                for (j = i; j <= ien; j++) {
                                                        *A(na - 1, j - 1) *= t1;                                                        *A(ien - 1, j - 1) *= t1;
                                                }
                                        }

                                }
                L_550:
                                ;
                        }

                }
                /* End d_complex vector */
L_560:
                ;
        }
        /*
         * End back substitution. Vectors of isolated roots
         */
        for (i = 1; i <= (*low - 1); i++) {

                for (j = i; j <= *n; j++) {
                        *ECOPY(j - 1, i - 1) = *A(j - 1, i - 1);
                }
        }

        for (i = *igh + 1; i <= *n; i++) {

                for (j = i; j <= *n; j++) {
                        *ECOPY(j - 1, i - 1) = *A(j - 1, i - 1);
                }
        }

        /*
         * Multiply by transformation matrix to give vectors of original full
         * matrix for j=n step -1 until low do --
         */
        for (j = *n; j >= *low; j--) {
                for (i = *low; i <= *igh; i++) {

                        zz = ZERO;

                        for (k = *low; k <= imsl_i_min(j, *igh); k++) {
                                zz += *ECOPY(k - 1, i - 1) ** A(j - 1, k - 1);
                        }
                        *ECOPY(j - 1, i - 1) = zz;
                }
        }

L_9000:
        imsl_e1pop("imsl_e5crh ");
        return;
}                               /* end of function */


#undef CONST
#undef CONST1
#undef CONST2
#undef ONE
#undef TWO
#undef ZERO
#undef A
#undef ECOPY
