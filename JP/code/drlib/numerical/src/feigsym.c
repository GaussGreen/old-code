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

#define ADR(t,x)    ( t = x, &t )



#ifdef ANSI

static VA_LIST_HACK l_eig_sym(Mint n, Mfloat *a, va_list argptr);

static void l_e5csf(Mint *n, Mfloat *a, Mint *lda, Mfloat *eval, Mfloat *evec,

            Mint *ldevec, Mfloat *wk, Mint *iwk);

static void l_e7csf(Mint *lb, Mint *nb, Mint *n, Mfloat *eval, Mfloat *dwk,

            Mfloat *e, Mfloat *fnorm, Mfloat *evec, Mint *ldevec, Mfloat *wk);

static void l_e3fsf(Mint *n, Mint *mxeval, Mfloat *a, Mint *lda, Mfloat *elow,

            Mfloat *ehigh, Mint *neval, Mfloat eval[], Mfloat *evec,

            Mint *ldevec, Mfloat wk[], Mint iwk[]);

static void l_e4bsf(Mint *n, Mint *mxeval, Mfloat *elow, Mfloat *ehigh,

            Mint *neval, Mfloat eval[], Mfloat d[], Mfloat e[],

            Mfloat e2[], Mfloat wk1[], Mfloat wk2[]);

static void l_e5bsf (Mint *n, Mint *mxeval, Mfloat *a, Mint *lda,

            Mfloat *elow, Mfloat *ehigh, Mint *neval, Mfloat eval[],

            Mfloat *wk, Mint iwk[]);

static void l_e6esf(Mint *n, Mint *nevec, Mfloat eval[], Mfloat *evec,

            Mint *ldevec, Mfloat d[], Mfloat e[],

            Mfloat *wk);

static void l_e4lsf(Mint *n, Mfloat *a, Mint *lda, Mfloat *eval, Mfloat *wk,

            Mint *iwk);

static void l_e5lsf(Mint *n, Mfloat *d, Mfloat *e2, Mint *iwk);

static void l_e6csf(Mint *n, Mfloat *a, Mint *lda, Mfloat *eval, Mfloat *e,

            Mfloat *e2, Mfloat *evec, Mint *ldevec, Mint *eigvec, 

	    Mfloat *scale);

static void l_scopy(Mint n, Mfloat sx[], Mint incx, Mfloat sy[], Mint incy);

#else

static VA_LIST_HACK	l_eig_sym();

static void l_e4lsf();

static void l_e5csf();

static void l_e5lsf();

static void l_e6csf();

static void l_e7csf();

static void l_e3fsf();

static void l_e4bsf();

static void l_e5bsf();

static void l_e6esf();

static void l_scopy();

#endif



static Mfloat	*lv_eval;



#ifdef ANSI

Mfloat *imsl_f_eig_sym(Mint n, Mfloat *a, ...)

#else

Mfloat *imsl_f_eig_sym(n, a, va_alist)

    Mint	n;

    Mfloat	*a;

    va_dcl

#endif

{

    va_list	argptr;



    VA_START(argptr,a);



#ifdef DOUBLE

    imsl_e1psh("imsl_d_eig_sym");

#else

    imsl_e1psh("imsl_f_eig_sym");

#endif

    lv_eval = NULL;

    IMSL_CALL(l_eig_sym(n, a, argptr));

    va_end(argptr);

#ifdef DOUBLE

    imsl_e1pop("imsl_d_eig_sym");

#else

    imsl_e1pop("imsl_f_eig_sym");

#endif

    return lv_eval;

}





#ifdef ANSI

static VA_LIST_HACK l_eig_sym(Mint n, Mfloat *a, va_list argptr)

#else

static VA_LIST_HACK l_eig_sym(n, a, argptr)

    Mint	n;

    Mfloat	*a;

    va_list	argptr;

#endif

{

    Mint	    code;

    Mint	    arg_number  = 2;

    Mint	    a_col_dim   = n;

    Mfloat	    **evec	= NULL;

    Mfloat	    *evecu	= NULL;

    Mfloat          *a_temp     = NULL;

    Mdouble	    elow	= -1000000.0;

    Mdouble	    ehigh	=  1000000.0;

    Mfloat	    elow_float;

    Mfloat	    ehigh_float;

    Mint	    i;

    Mint	    j;

    Mint	    count;

    Mint	    evecu_col_dim = n;

    Mint	    user_range  = 0;

    Mint	    vectors	= 0;

    Mint	    vectors_user = 0;

    Mint            user_number = 0;

    Mint	    evals_user	= 0;

    Mint	    mxeval	= n;

    Mint	    *neval	= 0;

    Mint            neval2      = 0;

    Mint	    *iwork	= NULL;

    Mfloat	    *work	= NULL;

    Mfloat	    *non_square	= NULL;



    code = 1;

    while (code > 0) {

	code = va_arg(argptr, int);

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

/*	(5,1, "Optional arguments IMSL_RANGE and IMSL_VECTORS_USER are incompatible.  Use IMSL_VECTORS instead.");

*/

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



    if (!user_range && (!vectors && !vectors_user)) {

	work    = (Mfloat *) imsl_malloc(2*n*sizeof(*work));

	iwork	= (Mint *) imsl_malloc(n*sizeof(*iwork));

    }

    if (!user_range && (vectors || vectors_user)) {

	work    = (Mfloat *) imsl_malloc(3*n*sizeof(*work));

	iwork	= (Mint *) imsl_malloc(n*sizeof(*iwork));

    }

    if (!user_range && vectors)

	*evec	= (Mfloat *) imsl_malloc(n*n*sizeof(**evec));



    if (user_range && (!vectors && !vectors_user)) {

	work 	= (Mfloat *) imsl_malloc(5*n*sizeof(*work));

	iwork 	= (Mint *) imsl_malloc(n*sizeof(*iwork));

    }

    if (user_range && vectors) {

	work    = (Mfloat *) imsl_malloc(9*n*sizeof(*work));

	*evec	= (Mfloat *) imsl_malloc(n*mxeval*sizeof(**evec));

	iwork	= (Mint *) imsl_malloc(n*sizeof(*iwork));

    }



/*This option is not allowed because of difficulty with leading */

/*dimension checks by low level subroutines.   

    if (user_range && vectors_user) {

	work    = (Mfloat *) imsl_malloc(9*n*sizeof(*work));

	iwork	= (Mint *) imsl_malloc(n*sizeof(*iwork));

    }

*/



    a_temp      = (Mfloat *) imsl_malloc(n*n*sizeof(*a));

    

    if ( (work==NULL) || (a_temp==NULL) || (iwork == NULL) || 

		(user_range && vectors && *evec==NULL) ) {

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

    scopy(n*n, a, 1, a_temp, 1);



    if (!user_range && (!vectors && !vectors_user))

	l_e4lsf(&n, a, &a_col_dim, lv_eval, work, iwork);

    if (!user_range && vectors) { 

        l_e5csf(&n, a, &a_col_dim, lv_eval, *evec, &n,

			 work, iwork);

	imsl_trnrr(n, n, *evec, n, n, n, *evec, n);

	}

    if (!user_range && vectors_user) {

        l_e5csf(&n, a, &a_col_dim, lv_eval, evecu, &evecu_col_dim,

			 work, iwork);

	imsl_trnrr(n, n, evecu, evecu_col_dim, n, n, evecu,

			 evecu_col_dim);

	}

    if (user_range && (!vectors && !vectors_user)) {

	l_e5bsf(&n, &mxeval, a, &a_col_dim, &elow_float, &ehigh_float,

			&neval2, lv_eval, work, iwork);

        /*&neval2 is needed because sending a (possibly) nil pointer, */

        /*neval, can cause an error.                                  */ 

        if (user_number) *neval = neval2;

    }

    if (user_range && vectors) {

	l_e3fsf(&n, &mxeval, a, &a_col_dim, &elow_float, &ehigh_float,

			&neval2, lv_eval, *evec, &n, work, iwork);

        /*&neval2 is needed because sending a (possibly) nil pointer, */

        /*neval, can cause an error.                                  */

        if (user_number)  *neval = neval2;

	if (neval2 > 0) {

        non_square = (Mfloat *) imsl_malloc(n*neval2*sizeof(*non_square));

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

        *evec = (Mfloat *) imsl_malloc(n*neval2*sizeof(**evec));

        if (*evec==NULL) {

          imsl_e1stl(1, "n"); 

          imsl_e1sti(1, n);

          imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

          goto FREE_SPACE;

        }

        scopy(n*neval2, non_square, 1, *evec, 1);

        if (non_square != NULL) imsl_free(non_square);

	}

    }



/*This option is not allowed because of difficulty with leading */

/*dimension checks by low level subroutines.                

    if (user_range && vectors_user) {

	l_e3fsf(&n, &mxeval, a, &a_col_dim, &elow_float, &ehigh_float,

			neval, lv_eval, evecu, &evecu_col_dim,

			work, iwork);

	imsl_trnrr(n, n, evecu, evecu_col_dim, n, n, evecu,

			 evecu_col_dim);

    }

*/



FREE_SPACE:

    if (work != NULL) imsl_free(work);

    if (iwork != NULL) imsl_free(iwork);

    if (a_temp != NULL) {

      scopy(n*n, a_temp, 1, a, 1);

      imsl_free(a_temp);

    }

RETURN:

    if (imsl_n1rty(0)>3) {

        if (!evals_user && (lv_eval != NULL)) imsl_free(lv_eval);

	lv_eval = NULL;

    }

    return (argptr);

}









/*Translated by FOR_C++, v0.1, on 07/30/90 at 17:05:54 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 07/30/90 at 17:05:52

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

       DWK    - The eigenvalues as returned by l_e5lsf. (Input)

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



/*WK is not used, but leave the calling sequence intact. */

#ifdef ANSI

static void l_e7csf(Mint *lb, Mint *nb, Mint *n, Mfloat *eval, Mfloat *dwk,

            Mfloat *e, Mfloat *fnorm, Mfloat *evec, Mint *ldevec, Mfloat *wk)

#else

static void l_e7csf(lb, nb, n, eval, dwk, e, fnorm, evec, ldevec, wk)

	Mint            *lb, *nb, *n;

	Mfloat           eval[], dwk[], e[], *fnorm, *evec;

	Mint            *ldevec;

	Mfloat           wk[];

#endif

{

#define EVEC(I_,J_)	(evec+(I_)*(sldevec)+(J_))

	Mint		sldevec = *ldevec;

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

	if (l == *nb) {

                *lb = l;

		goto L_9000;

	      }

	/*

	 * j represents how many iterations the QR was applied in order to

	 * obtain the l-th eigenvalue.

	 */

	j = 0;

	/*

	 * First test for convergence. dwk(1:n) are the perfect shifts as

	 * returned by l_e5lsf. l-th eigenvalue found ?

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

		sigma = eval[l - 1] - e[l] / sign(sqrt(powl(p, 2) + powl(e[l], 2)) +

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









/* Structured by FOR_STRUCT, v0.2, on 10/10/90 at 14:59:13

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  E4BSF/DE4BSF (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    August 27, 1990



    Purpose:    Compute the eigenvalues in a given interval of a

                symmetri

                tridiagonal matrix.



    Usage:      CALL E4BSF (N, MXEVAL, ELOW, EHIGH, NEVAL, EVAL, D, E,

                            E2, WK1, WK2)



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



    Remark:

       E4BSF is associated with EISPACK routine BISECT



    Chapter:    MATH/LIBRARY Eigensystem Analysis



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_e4bsf(Mint *n, Mint *mxeval, Mfloat *elow, Mfloat *ehigh,

		Mint *neval, Mfloat eval[], Mfloat d[], Mfloat e[],

		Mfloat e2[], Mfloat wk1[], Mfloat wk2[])

#else

static void l_e4bsf(n, mxeval, elow, ehigh, neval, eval, d, e,

      e2, wk1, wk2)

	Mint            *n, *mxeval;

	Mfloat          *elow, *ehigh;

	Mint            *neval;

	Mfloat           eval[], d[], e[], e2[], wk1[], wk2[];

#endif

{

	Mint             i, ip, iq, ir, is, isturm, j, k, l, m1, m2;

	Mfloat           elb, eps, eps1, t1, t2, u, ub, v, x0, x1, xu;



	imsl_e1psh("E4BSF ");

	/* INITIALIZATIONS */

	eps = imsl_amach(4);

	elb = *elow;

	ub = *ehigh;

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



/*		(3, 1, " The number of eigenvalues in the specified interval exceeds MXEVAL. NEVAL contains the number of eigenvalues in the .interval. No eigenvalues will be returned.");

*/

                imsl_ermes(IMSL_WARNING, IMSL_NO_EIGENVALUES_RETURNED);

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

			k += 1;

		} else {

			if (k > m2)

				goto L_200;

			if (wk2[k - 1] < eval[l - 1]) {

				l_scopy(is - j + 1, &eval[l - 1], -1, &eval[l], -1);

				eval[l - 1] = wk2[k - 1];

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

	imsl_e1pop("E4BSF ");



	return;

}				/* end of function */







/* Structured by FOR_STRUCT, v0.2, on 10/17/90 at 09:25:52

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  E6ESF/DE6ESF (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    August 1, 1990



    Purpose:    Computes the eigenvectors of a symmetric tridiagonal

                matrix corresponding to a set of ordered approximate

                eigenvalues, using inverse iteration.



    Usage:      CALL E6ESF (N, NEVEC, EVAL, EVEC, LDEVEC, D, E,

                            WK )



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

       D      - Real vector of length n containing the diagonal elements

                of the tri-diagonal matrix.  (Input)

       E      - Real vector of length N.  On input, the elements of the

                off diagonal.  E(1) is arbitrary.  (Input)

       WK     - Real work vector of length 7*N.  (Output)



    Chapter:    MATH/LIBRARY Eigensystem Analysis



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_e6esf(Mint *n, Mint *nevec, Mfloat eval[], Mfloat *evec,

			Mint *ldevec, Mfloat d[], Mfloat e[], 

			Mfloat *wk)

#else



static void l_e6esf(n, nevec, eval, evec, ldevec, d, e, wk)

	Mint            *n, *nevec;

	Mfloat           eval[], *evec;

	Mint            *ldevec;

	Mfloat           d[], e[], *wk;

#endif

{

#define EVEC(I_,J_)	(evec+(I_)*(sldevec)+(J_))

#define WK(I_,J_)	(wk+(I_)*(sn)+(J_))

	Mint		sldevec = *ldevec;

	Mint		sn = *n;

	Mint             _l0, _l1, i, ir, isave, iseed, its, j, k,

	                kb, kp1, nm1, nm2;

	Mfloat           _f0, _f1, anorm, eps, eval1, sc,

	                t, temp, temp1, temp2, temp3, temp4, vnorm, wnorm;





	imsl_e1psh("E6ESF ");



	eps = imsl_amach(4);

	temp4 = F_ONE / eps;

	/*

	 * Retrieve the current value of the random number generator.

	 */

/*	imsl_rnget(&isave); */

	isave = imsl_random_seed_get();

	iseed = 123457;

/*	imsl_rnset(&iseed); */

	imsl_random_seed_set(iseed);

	/*

	 * Compute the norm of the symmetric tridiagonal matrix A. This norm

	 * will be used later to check the eigenvector. The norm of the

	 * tridiagonal matrix A is a good estimate of the size of the largest

	 * eigenvalue of the original symmetric matrix A. Temp1 contains the

	 * sum of the diagona terms of the tridiagonal matrix A . Temp2

	 * contains the sum of the subdiagonal terms of the tridiagonal

	 * matrix A.

	 */

	e[0] = F_ZERO;

	temp1 = imsl_sasum(*n, d, 1);

	temp2 = imsl_sasum(*n, &e[0], 1);

	anorm = (temp1 + 2 * temp2) / *n;

	/*

	 * Find vectors by inverse iteration. This method finds the solution

	 * to the inhomogeneous set of equations ( A - eval*I)*x = b.

	 */

	for (ir = 1; ir <= *nevec; ir++) {

		/*

		 * Initialize the starting vector. Randomly get the

		 * right-hand side vector, b. 'B' is a function of eval, so

		 * we are effectively using a different right-hand side for

		 * each eval.

		 */

		for (i = 1; i <= *n; i++) {

/*			*WK(3, i - 1) = imsl_rnunf(); */

			imsl_f_random_uniform(1, IMSL_RETURN_USER,

						WK(3,i-1), 0);

		}

		/*

		 * Two iterations of inverse iteration. Using the random

		 * vector b makes 2 iterations of inverse iteration

		 * sufficient.

		 */

		for (its = 1; its <= 2; its++) {

			/*

			 * Copy e to wk(1,3). Wk(1,3) is a vector containing

			 * the superdiagonal of the tridiagonal matrix in

			 * wk(1,3) through wk(n-1,3).

			 */

			for (i = 1; i <= (*n - 1); i++) {

				*WK(2, i - 1) = e[i];

			}

			/* Save the subdiagonal. */

			scopy(*n, e, 1, WK(0, 0), 1);

			/*

			 * Define (A-eval(ir)*Identity) by subtracting

			 * eval(ir) from the diagonal. ( A is the tridiagonal

			 * matrix )

			 */

			for (i = 1; i <= *n; i++) {

				*WK(1, i - 1) = d[i - 1] - eval[ir - 1];

			}

			/*

			 * Solve symmetric tridiagonal system. WK(1,1)

			 * represents E(the subdiagonal) WK(1,2) represents

			 * (D-eval(ir)*I) WK(1,3) represents the

			 * superdiagonal WK(1,4) represents the right-hand

			 * side of the system. WK(1,4) will ultimately be the

			 * eigenvector for the specified eigenvalue.

			 */

			*WK(0, 0) = *WK(1, 0);

			nm1 = *n - 1;

			if (nm1 >= 1) {

				*WK(1, 0) = *WK(2, 0);

				*WK(2, 0) = F_ZERO;

				*WK(2, *n - 1) = F_ZERO;



				for (k = 1; k <= nm1; k++) {

					kp1 = k + 1;

					/* Find the largest of the two rows */

					if (fabs(*WK(0, kp1 - 1)) > fabs(*WK(0, k - 1))) {

						/* Interchange row */

						t = *WK(0, kp1 - 1);

						*WK(0, kp1 - 1) = *WK(0, k - 1);

						*WK(0, k - 1) = t;

						t = *WK(1, kp1 - 1);

						*WK(1, kp1 - 1) = *WK(1, k - 1);

						*WK(1, k - 1) = t;

						t = *WK(2, kp1 - 1);

						*WK(2, kp1 - 1) = *WK(2, k - 1);

						*WK(2, k - 1) = t;

						t = *WK(3, kp1 - 1);

						*WK(3, kp1 - 1) = *WK(3, k - 1);

						*WK(3, k - 1) = t;

					}

					/* Zero elements */

					if (*WK(0, k - 1) != F_ZERO) {

						t = -*WK(0, kp1 - 1) / *WK(0, k - 1);

					} else {

						t = F_ZERO;

					}

					*WK(0, kp1 - 1) = *WK(1, kp1 - 1) + t ** WK(1, k - 1);

					*WK(1, kp1 - 1) = *WK(2, kp1 - 1) + t ** WK(2, k - 1);

					*WK(2, kp1 - 1) = F_ZERO;

					*WK(3, kp1 - 1) += t ** WK(3, k - 1);

				}

			}

			/*

			 * Back solve If wk(n,1) .eq. 0.0, then set right

			 * hand side to a random number/eps*wk(n,4).

			 */

			if (*WK(0, *n - 1) != F_ZERO) {

				*WK(3, *n - 1) /= *WK(0, *n - 1);

			} else {

/*				*WK(3, *n - 1) *= imsl_rnunf() * (temp4); */

				imsl_f_random_uniform(1, IMSL_RETURN_USER,

							&temp, 0);

				*WK(3,*n-1) *= temp*temp4;

			}

			nm2 = *n - 2;

			if (*n != 1) {

				if (*WK(0, nm1 - 1) != F_ZERO) {

					*WK(3, nm1 - 1) = (*WK(3, nm1 - 1) - *WK(1, nm1 - 1) *

					  *WK(3, *n - 1)) / *WK(0, nm1 - 1);

				} else {

/*					*WK(3, nm1 - 1) *= imsl_rnunf() * (temp4); */

					imsl_f_random_uniform(1, IMSL_RETURN_USER, &temp, 0);

					*WK(3,nm1-1) *= temp*temp4;

				}

				if (nm2 >= 1) {

					for (kb = 1; kb <= nm2; kb++) {

						k = nm2 - kb + 1;

						if (*WK(0, k - 1) != F_ZERO) {

							*WK(3, k - 1) = (*WK(3, k - 1) - *WK(1, k - 1) *

									 *WK(3, k) - *WK(2, k - 1) ** WK(3, k + 1)) /

								*WK(0, k - 1);

						} else {

/*							*WK(3, k - 1) *= imsl_rnunf() * (temp4); */

							imsl_f_random_uniform(1,

IMSL_RETURN_USER, &temp, 0);

							*WK(3, k-1) *=

							temp*temp4;

						}

					}

				}

			}

			/*

			 * Make sure wk(*,4) is not getting too big.

			 */

			wnorm = imsl_snrm2(*n, WK(3, 0), 1);

			if (wnorm > temp4)

				goto L_70;

		}

		/*

		 * Set normalized value of EVEC WK4 is the stored eigenvector

		 */

L_70:

		sc = F_ONE / wnorm;

		sscal(*n, sc, WK(3, 0), 1);

		scopy(*n, WK(3, 0), 1, EVEC(ir - 1, 0), 1);

		/*

		 * Check the eigenvector by using a Rayleigh quotient. This

		 * involves computing (trans(X)*A*X)/ trans(X)*X. But, the

		 * previous normalization of the eigenvector takes into

		 * account trans(X)*X, so the numerator need only to be

		 * calculated. Finally, computing the distance between

		 * trans(X)*A*X and the true eigenvalue will reveal how good

		 * the eigenvector is. Note: tridiagonal = U + trans(U) Thus,

		 * trans(X)*A*X = trans(X)*U*X + trans(X)*trans(U)*X =

		 * trans(X)*U*X+trans(trans(X)*U*X) = 2*trans(X)*U*X.

		 * 

		 * First step: Multiply U*X

		 */

		for (i = 1; i <= (*n - 1); i++) {

			/*

			 * For j=1,n-1, multiply diag(j)*eigenvector(j)

			 */

			for (j = 1; j <= (*n - 1); j++) {

				*WK(4, j - 1) = d[j - 1] ** WK(3, j - 1);

			}

			/*

			 * For k=2,n, multiply the

			 * superdiagonal(k)*eigenvector(k)

			 */

			for (k = 2; k <= *n; k++) {

				*WK(5, k - 1) = e[k - 1] ** WK(3, k - 1);

			}

			/*

			 * For i = 1,n-1, compute the sum of

			 * .5*diag(i)*eigenvector(i) +

			 * superdiagonal(i+1)*eigenvector(i+1).

			 */

			*WK(6, i - 1) = (F_HALF) ** WK(4, i - 1) + *WK(5, i);

		}

		/* Compute .5*diag(n)*eigenvector(n) */

		*WK(6, *n - 1) = (F_HALF) * d[*n - 1] ** WK(3, *n - 1);

		/*

		 * Second Step: Multiply trans(X)*(U*X)

		 */

		temp3 = imsl_sdot(*n, WK(3, 0), 1, WK(6, 0), 1);

		/* Multiply 2*eval1 */

		eval1 = 2 * temp3;

		/*

		 * Third Step: Check abs(eval - eval1) .gt.

		 * norm(tridiagonal)*imsl_amach(4).

		 */

		if (fabs(eval[ir - 1] - eval1) > anorm * sqrt(eps)) {



/*			(3, 2, "Inverse iteration did not converge. Eigenvector is not correct for the specified eigenvalue.");

*/

                imsl_ermes(IMSL_WARNING, IMSL_SLOW_CONVERGENCE_2);

		}

	}

	/*

	 * Orthogonalize all of the vectors using HMGS algorithm.

	 */

	for (j = 1; j <= *nevec; j++) {

		vnorm = imsl_snrm2(*n, EVEC(j - 1, 0), 1);

		*WK(3, j - 1) = F_ZERO;

		if (vnorm == F_ZERO)

			goto L_120;

		*WK(3, j - 1) = F_ONE / vnorm;

		sscal(*n, *WK(3, j - 1), EVEC(j - 1, 0), 1);

		if (j == *nevec)

			goto L_120;

		_l0 = *nevec-j;

		_f0 = F_ONE;

		_l1 = 1;

		_f1 = F_ZERO;

		imsl_sgemv("T", sizeof("T"), n, &_l0, &_f0,

			   EVEC(j, 0), ldevec, EVEC(j - 1, 0), &_l1, &_f1,

			   WK(3, j), &_l1);

		imsl_sger(*n, *nevec - j, -F_ONE, EVEC(j - 1, 0), 1, WK(3, j), 1,

			  EVEC(j, 0), *ldevec);

L_120:

		;

	}

	/*

	 * Complete the orthogonalization by applying the transformations and

	 * check for loss of orthogonality.

	 */

	vnorm = F_ZERO;

	for (j = *nevec; j >= 1; j--) {

		if (*WK(3, j - 1) == F_ZERO)

			vnorm += F_ONE;

		if (j == *nevec)

			goto L_130;

		_l0 = *nevec-j;

		_f0 = F_ONE;

		_l1 = 1;

		_f1 = F_ZERO;

		imsl_sgemv("T", sizeof("T"), n, &_l0, &_f0,

			   EVEC(j, 0), ldevec, EVEC(j - 1, 0), &_l1, &_f1,

			   WK(3, j), &_l1);

		/*

		 * c                                   Compute norm of error

		 * matrix

		 */

		vnorm += imsl_sdot(*nevec - j, WK(3, j), 1, WK(3, j), 1);

		imsl_sger(*n, *nevec - j, -F_ONE, EVEC(j - 1, 0), 1, WK(3, j), 1,

			  EVEC(j, 0), *ldevec);

L_130:

		;

	}

	/*

	 * The eigenvectors may not be orthogonal. A perfect test is vnorm

	 * .le. imsl_amach(4)**2.

	 */

	if (vnorm > eps) {



/*		(3, 3, "The eigenvectors have lost orthogonality.");

*/

                imsl_ermes(IMSL_WARNING, IMSL_LOST_ORTHOGONALITY_2);

	}

	/*

	 * Initialize the seed to the number imsl_rnget retrieved initially.

	 */

/*	imsl_rnset(&isave); */

	imsl_random_seed_set(isave);

	imsl_e1pop("E6ESF ");



	return;

}				/* end of function */



#undef EVEC

#undef WK

/* Structured by FOR_STRUCT, v0.2, on 10/17/90 at 15:48:49

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  SCOPY (Single precision version)



    Computer:   FORC/SINGLE



    Revised:    August 9, 1986



    Purpose:    Copy a vector X to a vector Y, both single precision.



    Usage:      CALL SCOPY (N, SX, INCX, SY, INCY)



    Arguments:

       N      - Length of vectors X and Y.  (Input)

       SX     - Real vector of length MAX(N*IABS(INCX),1).  (Input)

       INCX   - Displacement between elements of SX.  (Input)

                X(I) is defined to be.. SX(1+(I-1)*INCX) if INCX .GE. 0

                or SX(1+(I-N)*INCX) if INCX .LT. 0.

       SY     - Real vector of length MAX(N*IABS(INCY),1).  (Output)

                SCOPY copies X(I) to Y(I) for I=1,...,N. X(I) and Y(I)

                refer to specific elements of SX and SY, respectively.

                See INCX and INCY argument descriptions.

       INCY   - Displacement between elements of SY.  (Input)

                Y(I) is defined to be.. SY(1+(I-1)*INCY) if INCY .GE. 0

                or SY(1+(I-N)*INCY) if INCY .LT. 0.



    Keywords:   Level 1 BLAS; SCOPY



    GAMS:       D1a



    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations

                STAT/LIBRARY Mathematical Support



    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */



#ifdef ANSI

static void l_scopy (Mint n, Mfloat sx[], Mint incx, Mfloat sy[], Mint incy)

#else

static void l_scopy(n, sx, incx, sy, incy)

	Mint             n;

	Mfloat           sx[];

	Mint             incx;

	Mfloat           sy[];

	int             incy;

#endif

{

	Mint             i, ix, iy;





	if (n > 0) {

		if (incx != 1 || incy != 1) {

			/* CODE FOR UNEQUAL INCREMENTS */

			ix = 1;

			iy = 1;

			if (incx < 0)

				ix = (-n + 1) * incx + 1;

			if (incy < 0)

				iy = (-n + 1) * incy + 1;

			for (i = 1; i <= n; i++) {

				sy[iy - 1] = sx[ix - 1];

				ix += incx;

				iy += incy;

			}

		} else {

			for (i = 1; i <= n; i++) {

				sy[i - 1] = sx[i - 1];

			}

		}

	}

	return;

}				/* end of function */



/*Translated by FOR_C++, v0.1, on 11/21/91 at 15:33:19 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 11/21/91 at 15:33:17

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  E4LSF/DE4LSF (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    May 31, 1991



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

static void l_e4lsf(Mint *n, Mfloat *a, Mint *lda, Mfloat *eval, Mfloat *wk, 

	Mint *iwk)

#else

static void l_e4lsf (n, a, lda, eval, wk, iwk)

    Mint        *n;

    Mfloat      *a;

    Mint        *lda;

    Mfloat       eval[], *wk;

    Mint         iwk[];

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

#define WK(I_,J_)	(wk+(I_)*(*n)+(J_))

    Mint         _l0, i;

    Mfloat       scale;





    /* Check N */

    if (*n < 1) {

	imsl_e1sti (1, *n);

/*

	(5, 1, "The argument N = %(i1).  The order of the matrix must be at least 1.");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);

    }

    if (imsl_n1rty (0) != 0)

	goto L_9000;

    /* Check LDA */

    if (*lda < *n) {

	imsl_e1sti (1, *lda);

	imsl_e1sti (2, *n);

/*

	(5, 2, "The argument LDA = %(i1).  The leading dimension of the matrix must be at least equal to the order of the matrix, N = %(i2).");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);

    }

    if (imsl_n1rty (0) != 0)

	goto L_9000;

    /*

     * The lower triangle is used. The upper part of the matrix is preserved.

     * A copy of 'A' is not needed. Reduce to tridiagonal Do not accumulate

     * transformations

     */

    _l0 = 0;

    l_e6csf (n, a, lda, eval, WK (0, 0), WK (1, 0), a, lda, &_l0,

	&scale);

    /* Matrix is zero; exit */

    if (scale == 0.0e0)

	goto L_9000;

    /* Find eigenvalues */

    l_e5lsf (n, eval, WK (1, 0), iwk);

    /* Scale back the eigenvalues. */

    if (scale != 1.0e0)

	sscal (*n, scale, eval, 1);

    /*

     * Sort eigenvalues into ascending order of magnitude.

     */

    imsl_svrbn (n, eval, WK (1, 0));

    /*

     * Resort into descending order of magnitude

     */

    for (i = 1; i <= *n; i++) {

	eval[i - 1] = *WK (1, *n - i);

    }

    /*

     * Retrieve the original symmetric matrix. Copy the upper triangle to

     * lower so that a is restored.

     */

    for (i = 1; i <= (*n - 1); i++) {

	scopy (*n - i, A (i, i - 1), *lda, A (i - 1, i), 1);

    }



L_9000:

    return;

}				/* end of function */





#undef A

#undef WK



/*----------------------------------------------------------------------- */



/*  IMSL Name:  E5LSF/DE5LSF (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    May 31, 1991



    Purpose:    Compute eigenvalues of a symmetric tridiagonal matrix.



    Usage:      CALL E5LSF (N, D, E2,  IWK)



    Arguments:

       N      - Order of the matrix.  (Input)

       D      - Real vector of length N.  On input, the diagonal elements

                of the matrix. On output, the eigenvalues of the

                matrix.  (Input/Output)

       E2     - Real vector of length N.  On input, the squares of

                the elements of the off diagonal.  E2(0) is arbitrary.

                On output, E2 is destroyed.  (Input/Output)

       IWK    - Integer work array.  (Input)



    Chapter:    MATH/LIBRARY Eigensystem Analysis



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  ----------------------------------------------------------------------- */

#ifdef ANSI

static void l_e5lsf(Mint *n, Mfloat *d, Mfloat *e2, Mint *iwk)

#else

static void l_e5lsf (n, d, e2, iwk)

    Mint        *n;

    Mfloat       d[], e2[];

    Mint         iwk[];

#endif

{

    Mint         i, iter, j, k, l, m;

    Mfloat       cossq, e, epslon, gama, h, p, prevc, prevga, recipe, sigma, sinsq;





    imsl_e1psh ("l_e5lsf ");



    if (*n == 1)

	goto L_9000;

    /*

     * The following code is based on the PWK algorithm, pages 164-169,

     * "Symmetric Eigenvalue Problem" B. Parlett.

     */

    epslon = imsl_fi_power (imsl_amach (4), 2);

    e2[0] = 0.0e0;

    iter = 0;

    /*

     * Bound on smallest number that can be reciprocated.

     */

    h = imsl_amach (1);

    if (h * imsl_amach (2) < 1.0e0)

	h = 1.0e0 / imsl_amach (2);

    /*

     * iwk() is used for the tracking when the current block decouples j is

     * the pointer.

     */

    j = 1;

    iwk[j - 1] = 1;

    /*

     * k is the start of the current block and m is the end of the current

     * block.

     */

    k = 1;

    m = *n;

    /*

     * Deflate the current block i times if possible.

     */

L_10:

    e = epslon * (0.5e0 * e2[m - 1] + d[m - 1] * d[m - 1]);

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

     * Check if current block is empty, i.e. no nonzero subdiagonal elements

     */

    if ((k == m) && (j > 1)) {

	j -= 1;

	k = iwk[j - 1];

	goto L_10;

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

     * Initialize the implicit shift sigma The Wilkinson shift is used. This

     * is the eigenvalue of the SE 2 by 2 nearest the SE entry.

     */

    e = (d[m - 2] - d[m - 1]) * 0.5e0;

    sigma = d[m - 1] - e2[m - 1] / sign (sqrt (imsl_fi_power (e, 2) + e2[m - 1]) +

	fabs (e), e);

    /* Initialize variables for the QR */

    cossq = 1.0e0;

    sinsq = 0.0e0;

    gama = d[k - 1] - sigma;

    p = gama * gama;

    /*

     * The next 13 lines are executed outside the main loop to avoid a

     * conditional statement. Start the QR factorization. Pre and Post

     * multiply the current block by the first plane rotation that attempts

     * to eliminate the off diagonal element b(k) of the block.

     */

    e = p + e2[k];

    recipe = 1.0e0 / e;

    prevc = cossq;

    cossq = p * recipe;

    sinsq = e2[k] * recipe;

    prevga = gama;

    gama = cossq * (d[k] - sigma) - sinsq * prevga;

    d[k - 1] = prevga + (d[k] - gama);

    p = e2[k];

    if (cossq != 0.0e0)

	p = imsl_fi_power (gama, 2) / cossq;

    /*

     * Finish the QR factorization of the current block. The following loop

     * "chases the bulge".

     */

    for (i = k + 1; i <= (m - 1); i++) {

	e = p + e2[i];

	recipe = 1.0e0 / e;

	e2[i - 1] = sinsq * e;

	prevc = cossq;

	/*

	 * cossq and sinsq repesent the i-th plane rotation, where cossq +

	 * sinsq = 1.

	 */

	cossq = p * recipe;

	sinsq = e2[i] * recipe;

	prevga = gama;

	gama = cossq * (d[i] - sigma) - sinsq * prevga;

	/*

	 * update i-th diagonal element that represents the shifted i-th

	 * eigenvalue.

	 */

	d[i - 1] = prevga + (d[i] - gama);

	if (cossq == 0.0e0) {

	    p = prevc * e2[i];

	}

	else {

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

     * Iterate back to the top The iteration counter(.le. 100) is arbitrary.

     * Deflation is done if 100 iterations were taken. Usually only 2 - 6

     * iterations are required.

     */

    if (iter <= 100) {

	goto L_10;

    }

    else {

/*

	(3, 1, "The iteration for the eigenvalue failed  to converge in 100 iterations before               deflating ");

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

    imsl_e1pop ("l_e5lsf ");



    return;

}				/* end of function */

#define trunc( f )       (double)((Mint)(f))





/*Translated by FOR_C++, v0.1, on 11/21/91 at 16:01:44 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 11/21/91 at 16:01:41

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  E6CSF/DE6CSF (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    May 31, 1991



    Purpose:    Reduce a real symmetric matrix to a tridiagonal

                matrix by orthogonal similarity transformations.



    Usage:      CALL E6CSF (N, A, LDA, EVAL, E, E2, EVEC, LDEVEC,

                            EIGVEC, SCALE)



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

       EVEC   - Accumulated product of transformations.  (Output)



       LDEVEC - Leading dimension of EVEC exactly as specified in the

                dimension statement of the calling program.  (Input)

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

static void l_e6csf(Mint *n, Mfloat *a, Mint *lda, Mfloat *eval, Mfloat *e,

            Mfloat *e2, Mfloat *evec, Mint *ldevec, Mint *eigvec, Mfloat *scale)

#else

static void l_e6csf (n, a, lda, eval, e, e2, evec, ldevec, eigvec,

                scale)

    Mint        *n;

    Mfloat      *a;

    Mint        *lda;

    Mfloat       eval[], e[], e2[], *evec;

    Mint        *ldevec;

    Mint  *eigvec;

    Mfloat      *scale;

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))

    Mint         _l0, i, j, k;

    Mfloat       c, p, s, t, tol;





    imsl_e1psh ("l_e6csf ");



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

     * Tol = imsl_amach(1)/imsl_amach(4) is used to ensure that a transformation

     * is skipped whenever there is a danger that it might not be accurately

     * orthogonal.

     */

    tol = imsl_amach (1) / imsl_amach (4);

    /*

     * Preserve the diagonal elements and find the element of maximum

     * magnitude

     */

    *scale = 0.0e0;

    for (j = 1; j <= *n; j++) {

	eval[j - 1] = *A (j - 1, j - 1);

	for (i = j; i <= *n; i++) {

	    if (fabs (*A (j - 1, i - 1)) > *scale)

		*scale = fabs (*A (j - 1, i - 1));

	}

    }

    /* Matrix is all zero. */

    if (*scale == 0.0e0) {

	/* Set the vector eval to zero */

	sset (*n, 0.0e0, eval, 1);

	/* Set evec to the identity matrix */

	if (*eigvec) {

	    for (i = 1; i <= *n; i++) {

		sset (*n, 0.0e0, EVEC (i - 1, 0), 1);

		*EVEC (i - 1, i - 1) = 1.0e0;

	    }

	}

	/* Exit */

	goto L_220;

    }

    /* Reduction is trival */

    if (*n == 1) {

	*scale = 1.0e0;

	if (*eigvec)

	    *EVEC (0, 0) = 1.0e0;

	goto L_220;

    }

    /*

     * Form the ratio imsl_alog10(scale) to imsl_alog(base of machine)

     */

  

    if ((Mint) (*scale) == imsl_i_machine(6))

        k = 1;

    else

        k = trunc (log10 (*scale) / imsl_amach (5));

    /* imsl_imach(3) is the base of the machine */

    s = 1.0e0 / (Mfloat) (imsl_imach (3));

    /* Adjust for scale < 1.0 */

    if (*scale < 1.0e0) {

	s = 1.0e0 / s;

	k = -k;

    }

    /*

     * Compute scale. This product is typically exact.

     */

    *scale = 1.0e0;

    for (i = 1; i <= k; i++) {

	*scale *= s;

    }

    /*

     * Scale the matrix (lower triangle) This product is typically exact.

     */

    if (*scale != 1.0e0) {

	for (j = 1; j <= *n; j++) {

	    for (i = j; i <= *n; i++) {

		*A (j - 1, i - 1) *= *scale;

	    }

	}

    }



    for (k = 1; k <= (*n - 2); k++) {

	/*

	 * swap of the original and transformed diagonal elements

	 */

	s = eval[k - 1];

	eval[k - 1] = *A (k - 1, k - 1);

	*A (k - 1, k - 1) = s;

	/*

	 * computation of the square of the elements of a(k+2:n,k)

	 */

	e2[k] = 0.0e0;

	e[k] = *A (k - 1, k);

	for (i = k + 2; i <= *n; i++) {

	    e[i - 1] = *A (k - 1, i - 1);

	    e2[k] += e[i - 1] * e[i - 1];

	}

	/*

	 * Skip the k-th transformation if the sum of squares is too small

	 * for orthogonality to be guaranteed.

	 */

	if ((e2[k] + imsl_fi_power (e[k], 2)) <= tol) {

	    e[k] = 0.0e0;

	    *A (k - 1, k) = 0.0e0;

	    e2[k] = 0.0e0;

	    /*

	     * If e2(k+1) is zero then column k is in tridiagonal form, no

	     * reduction necesssary, skip to column k+1.

	     */

	}

	else if (e2[k] <= tol) {

	    e[k] = *A (k - 1, k);

	    *A (k - 1, k) = 0.0e0;

	    e2[k] = e[k] * e[k];

	    /*

	     * Reduce column k to upper hessenberg form using Householder

	     * reflectors.

	     */

	}

	else {

	    /*

	     * compute the l2 norm of a(k+1:n, k) and the codiagonal element.

	     */

	    e2[k] += imsl_fi_power (e[k], 2);

	    p = -sign (sqrt (e2[k]), e[k]);

	    t = 1.0e0 / p;

	    /*

	     * form u   ; store in a(k+1:n, k) k+1 Note that t*e(k+1) <= 0.0.

	     */

	    *A (k - 1, k) = t * e[k] - 1.0e0;

	    for (i = k + 2; i <= *n; i++) {

		*A (k - 1, i - 1) = t * e[i - 1];

	    }

	    /* Normalizing factor of the reflection */

	    t = 1.0e0 / *A (k - 1, k);

	    /*

	     * ----------------------------------- The active submatrix is

	     * now updated. k An equivalent form for A can be employed: k

	     * k-1           t          t A  = A    +   u  * f   +  f * u k+1

	     * k+1

	     * 

	     * k-1 f = imsl_gamma*( A * u   +  sigma*u ) k+1          k+1

	     * 

	     * t    k-1 sigma = 0.5*u(1)*( u  * A  * u ) k+1       k+1

	     * ----------------------------------- k-1 Compute A * u ; store

	     * in e(k+1,...,) k+1

	     */

	    e[*n - 1] = *A (k - 1, *n - 1) ** A (*n - 1, *n - 1);

	    for (j = *n - 1; j >= (k + 1); j--) {

		s = *A (k - 1, j - 1) ** A (j - 1, j - 1);

		c = *A (k - 1, j - 1);

		for (i = *n; i >= (j + 1); i--) {

		    s += *A (k - 1, i - 1) ** A (j - 1, i - 1);

		    e[i - 1] += c ** A (j - 1, i - 1);

		}

		e[j - 1] = s;

	    }

	    /*

	     * t    k-1 Compute  u  * A  * u k+1       k+1

	     */

	    s = 0.0e0;

	    for (i = k + 1; i <= *n; i++) {

		s += e[i - 1] ** A (k - 1, i - 1);

	    }

	    /* Compute sigma */

	    s *= 0.5e0 * t;

	    /* Compute f */

	    for (i = k + 1; i <= *n; i++) {

		e[i - 1] = t * (e[i - 1] + s ** A (k - 1, i - 1));

	    }

	    /*

	     * k Compute A ;transformation of the unreduced submatrix .

	     */

	    for (j = k + 1; j <= *n; j++) {

		/*

		 * ******* IBM uses directive ( IGNORE RECDEPS ******

		 */

		s = e[j - 1];

		c = *A (k - 1, j - 1);

		for (i = j; i <= *n; i++) {

		    *A (j - 1, i - 1) += s ** A (k - 1, i - 1) + c * e[i - 1];

		}

	    }

	    e[k] = p;

	}

    }

    /*

     * Take care of the last elements of the matrix.

     */

    s = eval[*n - 2];

    eval[*n - 2] = *A (*n - 2, *n - 2);

    *A (*n - 2, *n - 2) = s;

    e2[*n - 1] = imsl_fi_power (*A (*n - 2, *n - 1), 2);

    e[*n - 1] = *A (*n - 2, *n - 1);

    s = eval[*n - 1];

    eval[*n - 1] = *A (*n - 1, *n - 1);

    *A (*n - 1, *n - 1) = s;

    /* Reciprocate scale for output */

    *scale = 1.0e0 / *scale;

    /*

     * Construct right operator of the similarity transformation by

     * accumulating the Householder matricies.

     */

L_160:

    if (*eigvec) {

	/* Set the matrix Z to zero */

	for (i = 1; i <= *n; i++) {

	    sset (*n, 0.0e0, EVEC (i - 1, 0), 1);

	}

	if (*n > 1)

	    *EVEC (*n - 2, *n - 2) = 1.0e0;

	*EVEC (*n - 1, *n - 1) = 1.0e0;

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

	    *EVEC (k - 1, k - 1) = 1.0e0;

	    /*

	     * u(1) associated with u   is the k+1 reciprocal of a(k+1,k).

	     */

	    if (*A (k - 1, k) != 0.0e0) {

		t = 1.0e0 / *A (k - 1, k);

		for (j = k + 1; j <= *n; j++) {

		    s = 0.0e0;

		    /*

		     * Form the dot product of a(k+1:n,k) with evec(k+1:n,j)

		     */

		    for (i = k + 1; i <= *n; i++) {

			s += *A (k - 1, i - 1) ** EVEC (j - 1, i - 1);

		    }

		    /*

		     * Normalize the dot products by a(k+1,k)

		     */

		    s *= t;

		    /*

		     * Form the saxpy of evec(k+1:n,j) <- evec(k+1:n,j) + s *

		     * a(k+1:n,k)

		     */

		    for (i = k + 1; i <= *n; i++) {

			*EVEC (j - 1, i - 1) += s ** A (k - 1, i - 1);

		    }

		}

	    }

	}

    }



L_220:

    e2[0] = 0.0e0;

    e[0] = 0.0e0;

L_9000:

    imsl_e1pop ("l_e6csf ");

    return;

}				/* end of function */









/*Translated by FOR_C++, v0.1, on 01/06/92 at 10:42:43 */

/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */

/* Structured by FOR_STRUCT, v0.2, on 01/06/92 at 10:42:38

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  E5CSF/DE5CSF (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    May 31, 1991



    Purpose:    Compute all of the eigenvalues and eigenvectors of a real

                symmetric matrix.



    Usage:      CALL E5CSF (N, A, LDA, EVAL, EVEC, LDEVEC, WK, IWK)



    Arguments:  (See EVCSF)



    Chapter:    MATH/LIBRARY Eigensystem Analysis



    Copyright:  1990 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_e5csf (Mint *n, Mfloat *a, Mint *lda, Mfloat *eval,

		Mfloat *evec, Mint *ldevec, Mfloat *wk, Mint *iwk)

#else

static void l_e5csf (n, a, lda, eval, evec, ldevec, wk, iwk)

    Mint        *n;

    Mfloat      *a;

    Mint        *lda;

    Mfloat       eval[], *evec;

    Mint        *ldevec;

    Mfloat      *wk;

    Mint         iwk[];

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))

#define WK(I_,J_)	(wk+(I_)*(*n)+(J_))

    Mint   first;

    Mint         _l0, i, j, k, lb, nb;

    Mfloat       fnorm, h, scale, test;





    imsl_e1psh ("imsl_e5csf ");

    /*

     * wk(1,1) is used in l_e5lsf and imsl_e7csf wk(1,2) and wk(1,3) are

     * the subdiagonal and the squares of the subdiagonal elements

     * respectively Check N

     */

    if (*n < 1) {

	imsl_e1sti (1, *n);

/*

	(5, 1, "The argument N = %(i1).  The order of the matrix must be at least 1.");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);

	goto L_9000;

    }

    /* Check LDA */

    if (*lda < *n) {

	imsl_e1sti (1, *lda);

	imsl_e1sti (2, *n);

/*

	(5, 2, "The argument LDA = %(i1).  The leading dimension of the matrix must be at least equal to the order of the matrix, N = %(i2).");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);

    }

    /* Check LDEVEC */

    if (*ldevec < *n) {

	imsl_e1sti (1, *ldevec);

	imsl_e1sti (2, *n);

/*

	(5, 3, "The argument LDEVEC = %(i1).  The leading dimension of the eigenvector matrix must be at least the order of the matrix, N = %(i2).");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_LDEVEC_VALUE_TOO_SMALL);

    }

    if (imsl_n1rty (0) > 0)

	goto L_9000;



    /* Reduction is trivial. */

    if (*n == 1) {

	*EVEC (0, 0) = 1.0e0;

	eval[0] = *A (0, 0);

	goto L_9000;

    }

    /*

     * The logical variable first is is used to prevent us from forming the

     * square of the off-diagonal elements since they are returned in wk(1,3)

     * from l_e6csf.

     * 

     * The lower half of 'a' is used. The upper half will remain intact. Thus, a

     * copy of 'A' is not needed.

     */

    first = TRUE;

    /*

     * Reduce matrix to tridiagonal form and accumulate the transformation

     * returned in evec(*,*).

     */

    l_e6csf (n, a, lda, eval, WK (1, 0), WK (2, 0), evec, ldevec, ADR (_l0, TRUE),

	&scale);

    /*

     * If scale is equal to zero, then evec is the identity matrix and eval

     * is set to the zero vector. Exit.

     */

    if (scale == 0.0e0)

	goto L_9000;

    /*

     * Form the Frobenius norm of the matrix used for a scale factor

     */

    fnorm = 0.0e0;

    for (i = 1; i <= *n; i++) {

	fnorm += imsl_fi_power (eval[i - 1], 2) + *WK (2, i - 1);

    }

    fnorm = sqrt (fnorm * scale);

    /*

     * Test = 1.0e0/imsl_amach(2) is the bound on smallest number that can be

     * reciprocated.

     */

    test = imsl_amach (1);

    if (test * imsl_amach (2) < 1.0e0)

	test = 1.0e0 / imsl_amach (2);

    /*

     * lb and nb are the start and end of current block respectively.

     */

    lb = 1;

    nb = *n;

    /* Find the start of the block. */

L_20:

    for (i = lb; i <= (nb - 1); i++) {

	h = imsl_amach (4) * (fabs (eval[i - 1]) + fabs (*WK (1, i)));

	if (test < h)

	    test = h;

	if (fabs (*WK (1, i)) > test)

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

	if (fabs (*WK (1, i)) <= test) {

	    nb = i;

	    goto L_60;

	}

    }

    /*

     * Make a copy of the diagonal elements of the tridiagonal matrix eval

     */

L_60:

    scopy (nb - lb + 1, &eval[lb - 1], 1, WK (0, lb - 1), 1);

    /*

     * Square off-diagonal elements if not the first time.

     */

    if (!first) {

	for (k = 1; k <= (nb - lb); k++) {

	    /*

	     * lb + k is used instead of lb because wk(lb,3) is a dummy

	     * element

	     */

	    *WK (2, lb + k - 1) = imsl_fi_power (*WK (1, lb + k - 1), 2);

	}

    }

    /*

     * Find the eigenvalues of the tridiagonal block.

     */

    l_e5lsf (ADR (_l0, nb - lb + 1), WK (0, lb - 1), WK (2, lb - 1), iwk);

    /*

     * Sort eigenvalues into ascending order of magnitude.

     */

    imsl_svrbn (ADR (_l0, nb - lb + 1), WK (0, lb - 1), WK (0, lb - 1));

    /*

     * Find the eigenvectors of the tri- diagonal block and update evec.

     * Employing the explicit perfect shift, accumulate the remaining

     * similarity transforms to find the eigenvectors.

     */





    l_e7csf (&lb, &nb, n, eval, WK (0, 0), WK (1, 0), &fnorm, evec, ldevec,

	WK (2, 0));

    /*

     * If lb less than n implies the matrix has split or that a Wilkinson

     * shift was used.

     */

    if (lb < nb) {

	/*

	 * This is to ensure that the off-diagon elements are squared for use

	 * in l_e5lsf

	 */

	first = FALSE;

	goto L_20;

	/*

	 * lb equal to nb implies that we are done with the current block.

	 */

    }

    else if (lb == nb) {

	lb = nb + 1;

	nb = *n;

	goto L_20;

    }

    /* Scale back the eigenvalues. */

L_80:

    if (scale != 1.0e0)

	sscal (*n, scale, eval, 1);

    /*

     * Final sort of eigenvalues and eigenvectors.

     */

    for (i = 1; i <= *n; i++) {

	iwk[i - 1] = i;

	*WK (2, i - 1) = -fabs (eval[i - 1]);

    }

    /*

     * Eigenvalues are sorted into ascending order by magnitude

     */

    imsl_svrgp (*n, WK (2, 0), WK (2, 0), iwk);

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

     * Finally leave the eigenvalues and eigenvectors sorted in descending

     * order of magnitude.

     */

    for (i = *n - 1; i >= 1; i--) {

	if (i == iwk[i - 1])

	    goto L_120;

	sswap (*n, EVEC (i - 1, 0), 1, EVEC (iwk[i - 1] - 1, 0), 1);

	sswap (1, &eval[i - 1], 1, &eval[iwk[i - 1] - 1], 1);

L_120:

	;

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

     * Retrieve the original symmetric matrix. Copy the upper triangle to

     * lower.

     */

    for (i = 1; i <= (*n - 1); i++) {

	scopy (*n - i, A (i, i - 1), *lda, A (i - 1, i), 1);

    }



L_9000:

    imsl_e1pop ("imsl_e5csf ");

    return;

}				/* end of function */

/*----------------------------------------------------------------------- */



/*  IMSL Name:  E3FSF/DE3FSF (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    February 19, 1991



    Purpose:    Compute eigenvalues and eigenvectors in a given range

                of a real symmetric matrix.



    Usage:      CALL E3FSF (N, MXEVAL, A, LDA, ELOW, EHIGH, NEVAL, EVAL,

                            EVEC, LDEVEC, WK, IWK)



    Arguments:  (See EVESF)



    Chapter:    MATH/LIBRARY Eigensystem Analysis



    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_e3fsf (Mint *n, Mint *mxeval, Mfloat *a, Mint *lda,

		Mfloat *elow, Mfloat *ehigh, Mint *neval, Mfloat eval[],

		Mfloat *evec, Mint *ldevec, Mfloat wk[], Mint iwk[])

#else

static void l_e3fsf (n, mxeval, a, lda, elow, ehigh, neval, eval,

                evec, ldevec, wk, iwk)

    Mint        *n, *mxeval;

    Mfloat      *a;

    Mint        *lda;

    Mfloat      *elow, *ehigh;

    Mint        *neval;

    Mfloat       eval[], *evec;

    Mint        *ldevec;

    Mfloat      *wk;

    Mint         iwk[];

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))

#define WK(I_,J_)	(wk+(I_)*(*n)+(J_))

    Mint         _l0, _l1, _l2, i, j, k;

    Mfloat       _f0, _f1, scale;





    imsl_e1psh ("E3FSF ");

    /* Check N */

    if (*n < 1) {

	imsl_e1sti (1, *n);

/*

	(5, 1, "The argument N = %(i1).  The order of the matrix must be at least 1.");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);

	goto L_9000;

    }

    /* Check MXEVAL */

    if (*mxeval <= 0 || *mxeval > *n) {

	imsl_e1sti (1, *mxeval);

	imsl_e1sti (2, *n);

/*

	(5, 2, "The argument MXEVAL = %(i1).  The maximum number of eigenvalues to be calculated must be greater than 0 and less than or equal to N = %(i2).");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_NUMBER_MAX_EIGENVALUES);

    }

    /* Check LDA */

    if (*lda < *n) {

	imsl_e1sti (1, *lda);

	imsl_e1sti (2, *n);

/*

	(5, 3, "The argument LDA = %(i1).  The leading dimension of the matrix must be at least equal to the order of the matrix, N = %(i2).");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);

    }

    /* Check ELOW against EHIGH */

    if (*elow >= *ehigh) {

	imsl_e1str (1, *elow);

	imsl_e1str (2, *ehigh);

/*

	(5, 4, "The lower limit of the interval in which the eigenvalues are to be sought, ELOW = %(r1), must be strictly less than the upper limit, EHIGH = %(r2).");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_NEED_ELOW_LESS_THAN_EHIGH);

    }

    /* Check LDEVEC */

    if (*ldevec < *n) {

	imsl_e1sti (1, *ldevec);

	imsl_e1sti (2, *n);

/*

	(5, 5, "The argument LDEVEC = %(i1).  The leading dimension of the eigenvector matrix must be at least equal to the order of the matrix, N = %(i2).");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_LDEVEC_VALUE_TOO_SMALL);

    }

    if (imsl_n1rty (0) > 0)

	goto L_9000;

    /* Reduce to tridiagonal */

    l_e6csf (n, a, lda, WK (0, 0), WK (1, 0), WK (2, 0), a, lda, ADR (_l0, FALSE),

	&scale);

    /*

     * Save the diagonal for use later in E6ESF.

     */

    scopy (*n, WK (0, 0), 1, WK (3, 0), 1);

    /* Find all eigenvalues */

    l_e5lsf (n, WK (0, 0), WK (2, 0), iwk);



    /*

     * Count the number of eigenvalues in the scaled interval. Wk(1,5) needs

     * to be initialized to zero. Otherwise extraneous data gets left after

     * placing the eigenvalues in the vector.

     */

    sset (*n, 0.0, WK (4, 0), 1);

    *neval = 0;

    for (i = 1; i <= *n; i++) {

	if ((*WK (0, i - 1) >= *elow / scale) && (*WK (0, i - 1) <= *ehigh /

		scale)) {

	    *neval += 1;

	    *WK (4, i - 1) = *WK (0, i - 1);

	}

    }

    /*

     * Sort wk(1,5) into ascending order of magnitude.

     */

    imsl_svrbn (n, WK (4, 0), WK (4, 0));

    /*

     * Resort into descending order of magnitude

     */

    for (i = 1; i <= *n; i++) {

	*WK (5, i - 1) = *WK (4, *n - i);

    }

    /* Copy the neval eigenvalues to eval. */

    scopy (*neval, WK (5, 0), 1, eval, 1);





    if (*neval <= *mxeval) {

	/* Restore the diagonal */

	scopy (*n, WK (3, 0), 1, WK (0, 0), 1);

	/*

	 * Find corresponding eigenvectors by using Inverse Iteration.

	 */

	l_e6esf (n, neval, eval, evec, ldevec, WK (0, 0), WK (1, 0), WK (2, 0));

	/*

	 * Back transform eigenvectors of tri-diagonal to original matrix

	 * This involves applying the householder transformations(stored in

	 * A) to evec.

	 */

	for (j = *n - 2; j >= 1; j--) {

	    if (*A (j - 1, j) != 0.0e0) {

		imsl_sgemv ("transpose", sizeof ("transpose"), ADR (_l0, *n -

			j), neval, ADR (_f0, 1.0e0 / *A (j - 1, j)), EVEC (0, j),

		    ldevec, A (j - 1, j), ADR (_l1, 1), ADR (_f1, 0.0e0), WK (0, 0),

		    ADR (_l2, 1));

		imsl_sger (*n - j, *neval, 1.0e0, A (j - 1, j), 1, WK (0, 0),

		    1, EVEC (0, j), *ldevec);

	    }

	}

    }

    /*

     * Scale back the eigenvalues. Use the scale factor from E6CSF.

     */

    if (scale != 1.0e0)

	sscal (*neval, scale, eval, 1);

    /*

     * Final sort of eigenvalues and eigenvectors.

     */

    for (i = 1; i <= *neval; i++) {

	iwk[i - 1] = i;

	*WK (0, i - 1) = -fabs (eval[i - 1]);

    }

    /*

     * Eigenvalues are sorted into ascending order by magnitude

     */

    

    if (*neval > 0)

   	 imsl_svrgp (*neval, WK (0, 0), WK (0, 0), iwk);

    /* Resort the record of the pivots. */

    for (i = 1; i <= *neval; i++) {

	for (j = i; j <= *neval; j++) {

	    if (iwk[j - 1] == i) {

		k = iwk[i - 1];

		iwk[i - 1] = j;

		iwk[j - 1] = k;

		goto L_60;

	    }

	}

L_60:

	;

    }

    /*

     * Finally leave the eigenvalues and eigenvectors sorted in descending

     * order of magnitude.

     */

    for (i = *neval - 1; i >= 1; i--) {

	if (i == iwk[i - 1])

	    goto L_70;

	sswap (*n, EVEC (i - 1, 0), 1, EVEC (iwk[i - 1] - 1, 0), 1);

	sswap (1, &eval[i - 1], 1, &eval[iwk[i - 1] - 1], 1);

L_70:

	;

    }

    /*

     * Normalize each eigenvector so that its biggest component is positive.

     * The eigenvectors then form a right-hand system.

     */

    for (j = 1; j <= *neval; j++) {

	i = imsl_isamax (*n, EVEC (j - 1, 0), 1);

	if (*EVEC (j - 1, i - 1) < 0.0) {

	    for (k = 1; k <= *n; k++) {

		*EVEC (j - 1, k - 1) = -*EVEC (j - 1, k - 1);

	    }

	}

    }

    /*

     * Copy the upper triangle to lower so that a is restored.

     */

    for (i = 1; i <= (*n - 1); i++) {

	scopy (*n - i, A (i, i - 1), *lda, A (i - 1, i), 1);

    }



L_9000:

    imsl_e1pop ("E3FSF ");

    return;

}				/* end of function */

/*----------------------------------------------------------------------- */



/*  IMSL Name:  E5BSF/DE5BSF (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    February 20,1991



    Purpose:    Compute the eigenvalues in a given range of a real

                symmetric matrix.



    Usage:      CALL E5BSF (N, MXEVAL, A, LDA, ELOW, EHIGH, NEVAL, EVAL,

                            WK, IWK)



    Arguments:  (See EVBSF)



    Chapter:    MATH/LIBRARY Eigensystem Analysis



    Copyright:  1991 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_e5bsf (Mint *n, Mint *mxeval, Mfloat *a, Mint *lda,

		Mfloat *elow, Mfloat *ehigh, Mint *neval, Mfloat eval[],

		Mfloat *wk, Mint iwk[])

#else

static void l_e5bsf (n, mxeval, a, lda, elow, ehigh, neval, eval,

                wk, iwk)

    Mint        *n, *mxeval;

    Mfloat      *a;

    Mint        *lda;

    Mfloat      *elow, *ehigh;

    Mint        *neval;

    Mfloat       eval[], *wk;

    Mint         iwk[];

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

#define WK(I_,J_)	(wk+(I_)*(*n)+(J_))

    Mint         _l0, i;

    Mfloat       scale;





    imsl_e1psh ("E5BSF ");

    /* Check N */

    if (*n < 1) {

	imsl_e1sti (1, *n);

/*

	(5, 1, "The argument N = %(i1).  The order of the matrix must be at least 1.");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_MATRIX_ORDER_TOO_SMALL);

	goto L_9000;

    }

    /* Check MXEVAL */

    if (*mxeval <= 0 || *mxeval > *n) {

	imsl_e1sti (1, *mxeval);

	imsl_e1sti (2, *n);

/*

	(5, 2, "The argument MXEVAL = %(i1).  The maximum number of eigenvalues to be calculated must be greater than 0 and less than or equal to N = %(i2).");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_NUMBER_MAX_EIGENVALUES);

	goto L_9000;

    }

    /* Check LDA */

    if (*lda < *n) {

	imsl_e1sti (1, *lda);

	imsl_e1sti (2, *n);

/*

	(5, 3, "The argument LDA = %(i1).  The leading dimension of the matrix must be at least equal to the order of the matrix, N = %(i2).");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_LDA_VALUE_TOO_SMALL);

    }

    /* Check ELOW against EHIGH */

    if (*elow >= *ehigh) {

	imsl_e1str (1, *elow);

	imsl_e1str (2, *ehigh);

/*

	(5, 4, "The lower limit of the interval in which the eigenvalues are to be sought, ELOW = %(r1), must be strictly less than the upper limit, EHIGH = %(r2).");

*/

	imsl_ermes(IMSL_TERMINAL, IMSL_NEED_ELOW_LESS_THAN_EHIGH);

    }

    if (imsl_n1rty (0) > 0)

	goto L_9000;



    /* Reduce to tridiagonal */

    l_e6csf (n, a, lda, WK (0, 0), WK (1, 0), WK (2, 0), a, lda, ADR (_l0, FALSE),

	&scale);

    /* Find all eigenvalues */

    l_e5lsf (n, WK (0, 0), WK (2, 0), iwk);



    /*

     * Count the number of eigenvalues in the scaled interval. Wk(1,4) needs

     * to be initialized to zero. Otherwise extraneous data gets left after

     * placing the neval eigenvalues in the vector(wk(1,4)).

     */

    sset (*n, 0.0, WK (3, 0), 1);

    *neval = 0;

    for (i = 1; i <= *n; i++) {

	if ((*WK (0, i - 1) >= *elow / scale) && (*WK (0, i - 1) <= *ehigh /

		scale)) {

	    *neval += 1;

	    *WK (3, i - 1) = *WK (0, i - 1);

	}

    }



    /* Copy wk(1,4) to wk(1,5). */

    scopy (*n, WK (3, 0), 1, WK (4, 0), 1);



    /*

     * Sort wk(1,5) into ascending order of magnitude.

     */

    imsl_svrbn (n, WK (4, 0), WK (4, 0));

    /*

     * Resort into descending order of magnitude and put into eval.

     */

    for (i = 1; i <= *n; i++) {

	*WK (3, i - 1) = *WK (4, *n - i);

    }



    scopy (*neval, WK (3, 0), 1, eval, 1);



    /*

     * Pick off the neval eigenvalues and scale them.

     */

    if (scale != 1.0e0)

	sscal (*neval, scale, eval, 1);

    /*

     * Retrieve the original symmetric matrix. Copy the upper triangle to

     * lower so that a is restored.

     */

    for (i = 1; i <= (*n - 1); i++) {

	scopy (*n - i, A (i, i - 1), *lda, A (i - 1, i), 1);

    }



L_9000:

    imsl_e1pop ("E5BSF ");

    return;

}				/* end of function */