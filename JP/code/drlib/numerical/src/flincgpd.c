#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif



#define ADR(t,x)    ( t = x, &t )

#define SIGN(A) ((A > 0.0) ? 1.0 : -1.0)

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



static VA_LIST_HACK PROTO (l_lin_sol_def_cg, (Mint n, void (*amultp) (Mfloat *, Mfloat *), Mfloat *b,

	            va_list));

static void PROTO (l_p2grc, (Mint *ido, Mint *n, Mfloat *x, Mfloat *p, Mfloat *r,

	            Mfloat *z, Mfloat *relerr, Mint *itmax, Mfloat *tri,

	            Mfloat *work, Mint *iwork));

    static void PROTO (l_p3grc, (Mint *comput, Mint *iter, Mfloat *alpha, Mfloat *imsl_beta,

	            Mfloat *tri, Mfloat *eigest, Mfloat *work, Mint

	           *iwork));

    static void PROTO (l_e3csb, (Mint *n, Mfloat *a, Mint *lda, Mint *ncoda, Mfloat *d,

	            Mfloat *e, Mint *vector, Mfloat *evec, Mint *ldevec));



    static void PROTO (l_e5asf, (Mint *, Mint *, Mfloat *,

	            Mfloat[], Mfloat[], Mfloat[], Mint *,

	            Mint *));

    static void PROTO (l_e3asb, (Mint *, Mint *, Mfloat *, Mint *,

	            Mint *, Mint *, Mfloat[], Mfloat *,

	            Mfloat[]));



    static Mfloat *lv_x;

    static Mint lv_iter;

    static Mint lv_jump;

    static Mfloat lv_relerr;

    static short int lv_return_relerr;

#ifdef ANSI

    Mfloat     *imsl_f_lin_sol_def_cg (Mint n, void (*amultp) (Mfloat *, Mfloat *), Mfloat *b,...)

#else

    Mfloat     *imsl_f_lin_sol_def_cg (n, amultp, b, va_alist)

    Mint        n;

    void        (*amultp) ();

    Mfloat     *b;

va_dcl

#endif

{

    va_list     argptr;

    VA_START (argptr, b);



    E1PSH ("imsl_f_lin_sol_def_cg", "imsl_d_lin_sol_def_cg");



    lv_x = NULL;

    IMSL_CALL (l_lin_sol_def_cg (n, amultp, b, argptr));

    va_end (argptr);



    E1POP ("imsl_f_lin_sol_def_cg", "imsl_d_lin_sol_def_cg");



    return lv_x;

}





#ifdef ANSI

static VA_LIST_HACK l_lin_sol_def_cg (Mint n, void (*amultp) (Mfloat *, Mfloat *), Mfloat *b, va_list

                argptr)

#else

static VA_LIST_HACK l_lin_sol_def_cg (n, amultp, b, argptr)

    Mint        n;

    void        (*amultp) ();

    Mfloat     *b;

    va_list     argptr;

#endif

{

#ifdef ANSI

    typedef void (*FUNC) (Mfloat *, Mfloat *);

#else

    typedef void (*FUNC) ();

#endif

    FUNC        precond;

    Mfloat     *diagonal = NULL;

    Mint        maxiter;

    Mfloat     *relative_error = NULL;



    int         code = 1;

    int         arg_number = 3;

    short int   user_solution = 0;

    short int   precondition = 0;

    short int   user_jacobi = 0;

    short int   user_maxiter = 0;

    short int   user_relerr = 0;

    short int   strict_pos_neg = 0;



    Mfloat     *p = NULL;

    Mfloat     *r = NULL;

    Mfloat     *z = NULL;

    Mfloat     *tri = NULL;

    Mfloat     *work = NULL;

    Mint       *iwork = NULL;



    Mint        ido, i;

    Mint        error = 0;

    Mint       *maxit;

    Mfloat      relerr = F_ZERO;

    relerr = sqrt (imsl_amach (4));

    lv_return_relerr = 0;



    if (n < 1) {

	imsl_e1sti (1, n);

	imsl_ermes (IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);

	++error;

    }



    if (!error) {

	maxiter = (Mint) sqrt ((Mfloat) n);

	if (maxiter < 1000)

	    maxiter = 1000;

    }



    while (code > 0) {

	code = va_arg (argptr, int);

	++arg_number;

	switch (code) {

	case IMSL_RETURN_USER:

	    lv_x = va_arg (argptr, Mfloat *);

	    ++arg_number;

	    user_solution = 1;

	    break;



	case IMSL_MAX_ITER:

	    maxit = va_arg (argptr, Mint *);

	    user_maxiter = 1;

	    ++arg_number;

	    break;



	case IMSL_REL_ERR:

	    relative_error = va_arg (argptr, Mfloat *);

	    user_relerr = 1;

	    lv_return_relerr = 1;

	    ++arg_number;

	    break;



	case IMSL_PRECOND:

	    precond = va_arg (argptr, FUNC);

	    precondition = 1;

	    ++arg_number;

	    break;



	case IMSL_JACOBI:

	    diagonal = va_arg (argptr, Mfloat *);

	    user_jacobi = 1;

	    ++arg_number;

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



    if (!b) {

	imsl_e1stl (1, "b");

	imsl_ermes (IMSL_TERMINAL, IMSL_UNEXPECTED_NULL_POINTER);

	++error;

    }



    if (user_solution && !lv_x) {

	imsl_e1stl (1, "x");

	imsl_e1stl (2, "IMSL_RETURN_USER");

	imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);

	++error;

    }

    if (user_jacobi) {

	if (!diagonal) {

	    imsl_e1stl (1, "diagonal");

	    imsl_e1stl (2, "IMSL_JACOBI");

	    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);

	    ++error;

	}

	else {

	    strict_pos_neg = SIGN (diagonal[0]);

	    i = -1;

	    while (strict_pos_neg && ++i < n) {

		if (diagonal[i] == F_ZERO || (SIGN (diagonal[i]) != strict_pos_neg)) {

		    imsl_ermes (IMSL_TERMINAL, IMSL_BAD_JACOBI_PRECONDITION);

		    error = 1;

		    strict_pos_neg = 0;

		}

	    }

	}

    }

    if (user_maxiter) {

	if (!maxit) {

	    imsl_e1stl (1, "maxiter");

	    imsl_e1stl (2, "IMSL_MAX_ITER");

	    imsl_ermes (IMSL_TERMINAL, IMSL_OPTIONAL_ARG_NULL_1);

	    ++error;

	}

	else {

	    maxiter = *maxit;

	    if (maxiter < 1) {

		imsl_e1stl (1, "maxiter");

		imsl_e1sti (1, maxiter);

		imsl_ermes (IMSL_TERMINAL, IMSL_BAD_MAX_ITERATIONS);

		++error;

	    }

	}

    }



    if (error)

	return argptr;



    p = (Mfloat *) imsl_malloc (n * sizeof (Mfloat));

    r = (Mfloat *) imsl_malloc (n * sizeof (Mfloat));

    z = (Mfloat *) imsl_malloc (n * sizeof (Mfloat));

    tri = (Mfloat *) calloc (2 * maxiter, sizeof (Mfloat));

    lv_jump = maxiter;

    work = (Mfloat *) imsl_malloc (5 * maxiter * sizeof (Mfloat));

    iwork = (Mint *) imsl_malloc (maxiter * sizeof (Mint));

    if (!iwork || !tri || !z || !r || !p) {

	imsl_e1sti (1, n);

	imsl_e1sti (2, maxiter);

	imsl_e1stl (1, "n");

	imsl_e1stl (2, "maxiter");

	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);

	error = 1;

    }

    else {

	if (!user_solution) {

	    lv_x = (Mfloat *) imsl_malloc (n * sizeof (Mfloat));

	    if (!lv_x) {

		imsl_e1sti (1, n);

		imsl_e1stl (1, "n");

		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);

		error = 1;

	    }

	    else {

		imsl_scopy (n, b, 1, lv_x, 1);

	    }

	}

	if (!error) {

	    ido = 0;

	    imsl_scopy (n, b, 1, r, 1);

	    while (ido < 3 && !error) {

		l_p2grc (&ido, &n, lv_x, p, r, z, &relerr, &maxiter, tri, work, iwork);





		if (!(error = imsl_n1rty (1) > 3 ? 1 : 0)) {

		    if (ido == 1) {

			amultp (p, z);

		    }

		    else if (ido == 2) {

			if (precondition) {

			    precond (r, z);

			}

			else if (user_jacobi) {

			    for (i = 0; i < n; i++) {

				z[i] = r[i] / diagonal[i];

			    }

			}

			else {

			    for (i = 0; i < n; i++)

				z[i] = r[i];

			}

		    }

		}

	    }

	}

    }

    imsl_free (iwork);

    imsl_free (work);

    imsl_free (tri);

    imsl_free (z);

    imsl_free (r);

    imsl_free (p);

    if (error) {

	if (!user_solution) {

	    imsl_free (lv_x);

	    lv_x = NULL;

	}

	if (user_maxiter) {

	    *maxit = F_ZERO;

	}

	if (user_relerr) {

	    *relative_error = F_ZERO;

	}

    }

    else {

	if (user_maxiter) {

	    *maxit = lv_iter;

	}

	if (user_relerr) {

	    *relative_error = lv_relerr;

	}

    }

    return argptr;

}



/* Structured by FOR_STRUCT, v0.2, on 08/31/90 at 09:46:28

    Options SET: fmt=t s=n

  -----------------------------------------------------------------------

    IMSL Name:  P2GRC/DP2GRC (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    February 10, 1986



    Purpose:    Solve a real symmetric definite linear system using a

                preconditioned conjugate gradient method with reverse

                communication.



    Usage:      CALL P2GRC (IDO, N, X, P, R, Z, RELERR, ITMAX, TRI,

                            WORK, IWORK)



    Arguments:  See PCGRC/DPCGRC



    Chapter:    MATH/LIBRARY Linear Systems



    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_p2grc (Mint *ido, Mint *n, Mfloat *x, Mfloat *p, Mfloat *r,

                Mfloat *z, Mfloat *relerr, Mint *itmax, Mfloat *tri,

                Mfloat *work, Mint *iwork)

#else

static void l_p2grc (ido, n, x, p, r, z, relerr, itmax, tri, work, iwork)

    Mint       *ido, *n;

    Mfloat      x[], p[], r[], z[], *relerr;

    Mint       *itmax;

    Mfloat      tri[], work[];

    Mint        iwork[];

#endif

{

    Mint        convrg;

    Mint        i;

    static Mint jump;

    Mfloat      t, xnorm, znorm;

    static Mfloat alpha, imsl_beta, defa, defm, eigest, pz, s;





    imsl_e1psh ("l_p2grc");

    /* Check */

    if (*ido == 0) {

	/* Check N */

	if (*n < 1) {

	    imsl_e1sti (1, *n);



	    imsl_ermes (IMSL_TERMINAL,

		IMSL_MATRIX_ORDER_NOT_POSITIVE);

	    goto L_9000;

	}

	/* Check ITMAX */

	if (*itmax < 1) {

	    imsl_e1sti (1, *itmax);

	    imsl_e1stl (1, "ITMAX");

	    imsl_ermes (IMSL_TERMINAL, IMSL_BAD_MAX_ITERATIONS);

	    goto L_9000;

	}

    }

    else if (*ido < 0 || *ido > 2) {

	imsl_e1sti (1, *ido);



	imsl_ermes (IMSL_TERMINAL, IMSL_WRONG_IDO_VALUE);

	goto L_9000;

    }

    else {

	/* Jump to correct place */

	if (jump == 1)

	    goto L_10;

	if (jump == 2)

	    goto L_30;

	if (jump == 3)

	    goto L_50;

    }

    /* Initialize */

    eigest = -F_ONE;

    /* Set p = x */

    imsl_scopy (*n, x, 1, p, 1);

    /* Compute z = A*p */

    *ido = 1;

    jump = 1;

    goto L_9000;

L_10:

    ;

    /* Compute residual */

    saxpy (*n, -F_ONE, z, 1, r, 1);

    /* For ITER = 1, 2, ... */

    lv_iter = 1;

L_20:

    ;

    /* Compute z = M**(-1)*r */

    *ido = 2;

    jump = 2;

    goto L_9000;

L_30:

    ;

    /* T = z . r */

    t = imsl_sdot (*n, z, 1, r, 1);

    /* Check the precondintioning matrix */

    if (t == F_ZERO) {

	/* Check for nonzero residual */

	if (imsl_sasum (*n, r, 1) != F_ZERO) {

/*	          (4, 1, "The preconditioning matrix is singular");*/

	    imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_M_MATRIX);

	}

	*ido = 3;

	goto L_9000;

    }

    else if (lv_iter == 1) {

	defm = sign (F_ONE, t);

    }

    else if (defm * t < F_ZERO) {

/*		(4, 2, "The preconditioning matrix is not definite.");*/

	imsl_ermes (IMSL_FATAL, IMSL_NOT_POSITIVE_DEFINITE_M);

	*ido = 3;

	goto L_9000;

    }

    if (lv_iter == 1) {

	imsl_beta = F_ONE;

	imsl_scopy (*n, z, 1, p, 1);

    }

    else {

	imsl_beta = t / s;

	for (i = 1; i <= *n; i++) {

	    p[i - 1] = z[i - 1] + imsl_beta * p[i - 1];

	}

    }

    s = t;

    /* Compute z = A*P */

    *ido = 1;

    jump = 3;

    goto L_9000;

L_50:

    ;

    /* PZ = p . z */

    pz = imsl_sdot (*n, p, 1, z, 1);



    /* Check linear system */

    if (pz == F_ZERO) {

/*		(4, 3, "The linear system is singular."); <already done>*/

/*	        imsl_e1stl(1,"linear system");  */

	imsl_ermes (IMSL_FATAL, IMSL_SINGULAR_A_MATRIX);

	*ido = 3;

	goto L_9000;

    }

    else if (lv_iter == 1) {

	defa = sign (F_ONE, pz);

    }

    else if (defa * pz < F_ZERO) {

/*            (4, 4, "The linear system is not definite."); <already done>*/

	imsl_ermes (IMSL_FATAL, IMSL_NOT_POSITIVE_DEFINITE_A);

	*ido = 3;

	goto L_9000;

    }

    alpha = s / pz;

    saxpy (*n, alpha, p, 1, x, 1);

    saxpy (*n, -alpha, z, 1, r, 1);

    /* Preliminary convergence check */

    znorm = imsl_snrm2 (*n, z, 1);

    xnorm = imsl_snrm2 (*n, x, 1);



/*

printf("%e  %e  %e  %e\n",

	znorm, xnorm, eigest, *relerr * (F_ONE - eigest) * xnorm);

*/



    convrg = lv_iter > 3 && znorm <= *relerr * (F_ONE - eigest) * xnorm;

    /*

     * Compute estimate of the largest eigenvalue of the iteration matrix

     */

    l_p3grc (&convrg, &lv_iter, &alpha, &imsl_beta, tri, &eigest, work, iwork);

    /* Check for convergence */

    convrg = convrg && znorm <= *relerr * (F_ONE - eigest) * xnorm;

    if (convrg) {

	*ido = 3;

	if (lv_return_relerr) {

	    lv_relerr = xnorm * (F_ONE - eigest);

	    if (lv_relerr != F_ZERO) {

		lv_relerr = znorm / lv_relerr;

	    }

	    else {

		lv_relerr = imsl_amach (6);

	    }

	}

    }

    else {

	lv_iter += 1;

	if (lv_iter > *itmax) {

	    imsl_e1sti (1, *itmax);



/*			(4, 5, "The iteration did not converge in MAXITER = %(i1) iterations.");*/

	    imsl_ermes (IMSL_FATAL, IMSL_NO_CONVERGENCE);

	    *ido = 3;

	    goto L_9000;

	}

	goto L_20;

    }



L_9000:

    ;

    imsl_e1pop ("l_p2grc");

    return;

}			     /* end of function */

/*----------------------------------------------------------------------- */



/*  IMSL Name:  P3GRC/DP3GRC (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    February 10, 1986



    Purpose:    Compute an estimate of the largest eigenvalue of the

                precondition conjugate gradient iteration matrix.



    Usage:      CALL P3GRC (COMPUT, ITER, ALPHA, BETA, TRI, EIGEST,

                            WORK, IWORK)



    Arguments:

       COMPUT - Logical variable.  (Input)

                If COMPUT is true, then a new EIGEST is computed.

       ITER   - Iteration number.  (Input)

       ALPHA  - Conjugate gradient parameter for the ITER iteration.

                (Input)

       BETA   - Conjugate gradient parameter for the ITER iteration.

                (Input)

       TRI    - Band symmetric tridiagonal matrix similar to the

                iteration matrix.  (Input/Output)

                TRI is set only by this routine and must not be changed

                between calls.  See Hageman and Young (1981), imsl_page 150,

                equation 7-5.9 for the matrix.

       EIGEST - Estimate of the largest eigenvalue.  (Output)

       WORK   - Workspace of length 4*ITER.  (Output)

       IWORK  - Integer workspace.  (Output)



    Remark:

       It two successive values of EIGEST produce the same value (within

       one percent) for the ratio used in the error tolerance,

       1/(1-EIGEST), then the iteration for the eigenvalue is assumed to

       have converged and no more computations for it are made.



    Chapter:    MATH/LIBRARY Linear Systems



    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_p3grc (Mint *comput, Mint *iter, Mfloat *alpha, Mfloat *imsl_beta,

                Mfloat *tri, Mfloat *eigest, Mfloat *work, Mint *iwork)

#else

static void l_p3grc (comput, iter, alpha, imsl_beta, tri, eigest, work,

                iwork)

    Mint       *comput;

    Mint       *iter;

    Mfloat     *alpha, *imsl_beta, tri[], *eigest, work[];

    Mint        iwork[];

#endif

{

    Mint        _l0, _l1, _l2, _l3, iacopy, irwk;

    Mfloat      eig[1], imsl_gamma, omega;

    static Mfloat alpold, eigold, gamold;

    Mfloat     *acopy;

    Mfloat     *wk;

    /* Initialize */

    if (*iter == 1) {

	*eigest = -F_ONE;

	imsl_gamma = *alpha;

	tri[0] = F_ZERO;

	tri[lv_jump + *iter - 1] = F_ONE - F_ONE / imsl_gamma;

    }

    else {

	/* Add new terms to tridiagonal matrix */

	imsl_gamma = F_ONE / (*imsl_beta / alpold + F_ONE / *alpha);

	omega = sqrt (*imsl_beta) / alpold;

	tri[lv_jump + *iter - 1] = F_ONE - F_ONE / imsl_gamma;

	tri[*iter - 1] = omega;

	if (*comput) {

	    /* Partition workspace */

	    iacopy = 1;

	    irwk = iacopy + 2 ** iter;

	    /*

	     * Compute largest Evalue of tridiagonal matrix

	     */

	    _l0 = 1;

	    _l1 = 2;

	    _l2 = 1, _l3 = 0;

	    imsl_f_m1ran (2, lv_jump, tri, tri);



	    acopy = (Mfloat *) malloc (2 ** iter * sizeof (*acopy));

	    wk = (Mfloat *) malloc (3 ** iter * sizeof (*wk));

	    l_e3asb (iter, &_l0, tri, &_l1, &_l2,

		&_l3, eig, acopy, wk);

	    if (*eig > 1.0) {

		printf ("recomputing\n");

		(*iter)--;

		l_e3asb (iter, &_l0, tri, &_l1, &_l2,

		    &_l3, eig, acopy, wk);

		(*iter)++;

	    }

	    free (acopy);

	    free (wk);



	    imsl_f_m1ran (lv_jump, 2, tri, tri);

	    /*

	     * If errors in Evalue computation then ignore and use of Evalue

	     */

	    if (imsl_n1rty (1) >= 4) {

		*eigest = eigold;

		imsl_e1mes (-1, 0, " ");

	    }

	    else {

		*eigest = eig[0];

	    }

	}

	else {

	    /* If COMPUT is false return old value */

	    *eigest = eigold;

	}

    }

    gamold = imsl_gamma;

    alpold = *alpha;

    eigold = *eigest;



    return;

}



/*----------------------------------------------------------------------- */



/*  IMSL Name:  E3CSB/DE3CSB (Single/Double precision version)



    Computer:   FORC/SINGLE



    Revised:    June 3, 1985



    Purpose:    Reduce a real band symmetric matrix to a symmetric

                tridiagonal matrix using and optionally accumulating

                orthogonal transformations.



    Usage:      CALL E3CSB (N, A, LDA, NCODA, D, E, VECTOR, EVEC, LDEVEC)



    Arguments:

       N      - Order of the matrix.  (Input)

       A      - Band symmetric matrix.  (Input/Output)

                On input, A is the band symmetric matrix.  On output,

                A is destroyed.

       LDA    - Leading dimension of A exactly as specified in the

                dimension statement of the calling program.  (Input)

       NCODA  - Number of codiagonals.  (Input)

       D      - Real vector containing the diagonal elements of the

                tridiagonal matrix.  (Output)

       E      - Real vector containing the subdiagonal elements of the

                tridiagonal matrix in E(2:N).  (Output)

       VECTOR - Logical variable, .TRUE. if both the eigenvectors and the

                eigenvalues are to be stored.  (Input)

       EVEC   - Real N by N matrix containing the orthogonal

                transformation matrix produced in the reduction to

                tridiagonal form.  If VECTOR = .FALSE. then it can be a

                dummy.  (Output)

       LDEVEC - Leading dimension of EVEC exactly as specified in the

                dimension statement of the calling program.  (Output)



    Remark:

       E3CSB is based on the EISPACK routine BANDR.



    Chapter:    MATH/LIBRARY Eigensystem Analysis



    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.



    Warranty:   IMSL warrants only that IMSL testing has been applied

                to this code.  No other warranty, expressed or implied,

                is applicable.



  -----------------------------------------------------------------------

 */

#ifdef ANSI

static void l_e3csb (Mint *n, Mfloat *a, Mint *lda, Mint *ncoda, Mfloat *d,

                Mfloat *e, Mint *vector, Mfloat *evec, Mint *ldevec)

#else

static void l_e3csb (n, a, lda, ncoda, d, e, vector, evec, ldevec)

    Mint       *n;

    Mfloat     *a;

    Mint       *lda, *ncoda;

    Mfloat      d[], e[];

    Mint       *vector;

    Mfloat     *evec;

    Mint       *ldevec;

#endif

{

#define A(I_,J_)	(a+(I_)*(*lda)+(J_))

#define EVEC(I_,J_)	(evec+(I_)*(*ldevec)+(J_))

    Mint        _d_l, _d_m, _do0, _do1, _l0, _l1, _l2, iugl, j, k, l, r;

    Mfloat      b1, b2, c2, dmin, dminrt, f1, f2, g, s2, sparam[5], u;





    dmin = imsl_fi_power (F_TWO, -64);

    dminrt = imsl_fi_power (F_TWO, -32);

    /* Initialize diagonal scaling matrix */

    sset (*n, F_ONE, d, 1);



    if (*vector) {

	for (j = 1; j <= *n; j++) {

	    sset (*n, F_ZERO, EVEC (j - 1, 0), 1);

	}

	sset (*n, F_ONE, evec, *ldevec + 1);

    }

    if (*ncoda == 0) {

	imsl_scopy (*n, A (0, 0), *lda, d, 1);

	sset (*n, F_ZERO, e, 1);

	goto L_9000;

    }

    else if (*ncoda == 1) {

	goto L_60;

    }

    for (k = 1; k <= (*n - 2); k++) {

	for (r = imsl_i_min (*ncoda, *n - k); r >= 2; r--) {

	    g = *A (k + r - 1, *ncoda - r);

	    *A (k + r - 2, 0) = *A (k + r - 2, *ncoda - r + 1);

	    iugl = k;



	    for (j = k + r, _do0 = DOCNT (k + r, *n, _do1 = *ncoda); _do0 > 0; j += _do1, _do0--) {

		if (g == F_ZERO)

		    goto L_30;

		b1 = *A (j - 2, 0) / g;

		b2 = b1 * d[j - 2] / d[j - 1];

		s2 = F_ONE / (F_ONE + b1 * b2);

		if (s2 < F_HALF) {

		    b1 = g / *A (j - 2, 0);

		    b2 = b1 * d[j - 1] / d[j - 2];

		    c2 = F_ONE - s2;

		    d[j - 2] *= c2;

		    d[j - 1] *= c2;

		    f1 = F_TWO ** A (j - 1, *ncoda - 1);

		    f2 = b1 ** A (j - 2, *ncoda);

		    *A (j - 1, *ncoda - 1) += -b2 * (b1 ** A (j - 1, *ncoda - 1) -

			*A (j - 1, *ncoda)) - f2;

		    *A (j - 2, *ncoda) += b2 * (b2 ** A (j - 1, *ncoda) +

			f1);

		    *A (j - 1, *ncoda) += b1 * (f2 - f1);



		    sparam[0] = F_ZERO;

		    sparam[2] = b2;

		    sparam[3] = -b1;

		    _l0 = j - iugl - 1;

		    _l1 = 1;

		    _l2 = 1;

		    imsl_srotm (_l0, A (j - 1, *ncoda - j + iugl),

			_l1, A (j - 2, *ncoda - j + iugl + 1), _l2,

			sparam);



		    iugl = j;

		    *A (j - 2, 0) += b2 * g;

		    _l0 = imsl_i_min (*ncoda, *n - j + 1)

			- 1;

		    _l1 = *lda - 1;

		    _l2 = *lda - 1;

		    if (j != *n) {

			imsl_srotm (_l0, A (j, *ncoda - 1), _l1, A (j, *ncoda - 2),

			    _l2, sparam);

			if (j + *ncoda <= *n)

			    g = b2 ** A (j + *ncoda - 1, 0);

		    }

		    _l0 = 1;

		    _l1 = 1;

		    if (*vector) {

			imsl_srotm (*n, EVEC (j - 1, 0), _l0, EVEC (j - 2, 0),

			    _l1, sparam);

		    }

		}

		else {

		    u = d[j - 2];

		    d[j - 2] = s2 * d[j - 1];

		    d[j - 1] = s2 * u;

		    f1 = F_TWO ** A (j - 1, *ncoda - 1);

		    f2 = b1 ** A (j - 1, *ncoda);

		    u = b1 * (f2 - f1) + *A (j - 2, *ncoda);

		    *A (j - 1, *ncoda - 1) = b2 * (b1 ** A (j - 1, *ncoda - 1) -

			*A (j - 2, *ncoda)) + f2 - *A (j - 1, *ncoda - 1);

		    *A (j - 2, *ncoda) = b2 * (b2 ** A (j - 2, *ncoda) + f1) +

			*A (j - 1, *ncoda);

		    *A (j - 1, *ncoda) = u;



		    sparam[0] = F_ONE;

		    sparam[1] = b2;

		    sparam[4] = b1;

		    _l0 = j - iugl - 1;

		    _l1 = 1;

		    _l2 = 1;

		    imsl_srotm (_l0, A (j - 2, *ncoda - j + iugl + 1),

			_l1, A (j - 1, *ncoda - j + iugl), _l2,

			sparam);



		    iugl = j;

		    *A (j - 2, 0) = b2 ** A (j - 2, 0) + g;

		    _l0 = imsl_i_min (*ncoda, *n - j + 1)

			- 1;

		    _l1 = *lda - 1;

		    _l2 = *lda - 1;

		    if (j != *n) {

			imsl_srotm (_l0, A (j, *ncoda - 2), _l1, A (j, *ncoda - 1),

			    _l2, sparam);



			if (j + *ncoda <= *n) {

			    g = *A (j + *ncoda - 1, 0);

			    *A (j + *ncoda - 1, 0) *= b1;

			}

		    }

		    _l0 = 1;

		    _l1 = 1;

		    if (*vector) {

			imsl_srotm (*n, EVEC (j - 2, 0), _l0, EVEC (j - 1, 0),

			    _l1, sparam);

		    }

		}



	    }



    L_30:

	    ;

	}



	if (mod (k, 64) == 0) {

	    /*

	     * Rescale to avoid underflow or overflow

	     */

	    for (j = k; j <= *n; j++) {

		if (d[j - 1] < dmin) {

		    l = imsl_i_max (1, *ncoda + 2 - j);

		    sscal (*ncoda - l + 1, dminrt, A (j - 1, l - 1),

			1);

		    sscal (imsl_i_min (*ncoda, *n - j), dminrt, A (j, *ncoda - 1),

			*lda - 1);

		    if (*vector)

			sscal (*n, dminrt, EVEC (j - 1, 0), 1);

		    *A (j - 1, *ncoda) *= dmin;

		    d[j - 1] /= dmin;

		}

	    }

	}

    }

L_60:

    ;

    /* Form square root of scaling matrix */

    for (j = 2; j <= *n; j++) {

	e[j - 1] = sqrt (d[j - 1]);

    }



    if (*vector) {

	for (k = 2; k <= *n; k++) {

	    sscal (*n, e[k - 1], EVEC (k - 1, 0), 1);

	}

    }

    u = F_ONE;



    for (j = 2; j <= *n; j++) {

	*A (j - 1, *ncoda - 1) *= u * e[j - 1];

	u = e[j - 1];

	*A (j - 1, *ncoda) *= d[j - 1];

	d[j - 1] = *A (j - 1, *ncoda);

	e[j - 1] = *A (j - 1, *ncoda - 1);

    }



    d[0] = *A (0, *ncoda);

    e[0] = F_ZERO;



L_9000:

    ;

    return;

}			     /* end of function */

#undef  A

#undef  EVEC



/*Translated by FOR_C, v3.4 (P), on 09/16/93 at 08:16:12 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn pf=,/home/elpaso/clib/include/imsl_int.h s=dvowa str=ln - prototypes */

/* Structured by FOR_STRUCT, v2.0, on 09/16/93 at 08:16:11

 *  Options SET: fmt=t s=n

 *-----------------------------------------------------------------------

 *  IMSL Name:  E3ASB/DE3ASB (Single/Double precision version)

 *

 *  Computer:   FORC/SINGLE

 *

 *  Revised:    August 28, 1990

 *

 *  Purpose:    Compute the largest or smallest eigenvalues of a real

 *              symmetric matrix in band symmetric storage mode.

 *

 *  Usage:      CALL E3ASB (N, NEVAL, A, LDA, NCODA, SMALL, EVAL, ACOPY,

 *                          WK)

 *

 *  Arguments:  (See EVASB)

 *

 *  Chapter:    MATH/LIBRARY Eigensystem Analysis

 *

 *  Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

 *

 *  Warranty:   IMSL warrants only that IMSL testing has been applied

 *              to this code.  No other warranty, expressed or implied,

 *              is applicable.

 *

 *-----------------------------------------------------------------------

 * */

#ifdef ANSI

static void l_e3asb (Mint *n, Mint *neval, Mfloat *a, Mint *lda,

                Mint *ncoda, Mint *small, Mfloat eval[], Mfloat *acopy,

                Mfloat wk[])

#else

static void l_e3asb (n, neval, a, lda, ncoda, small, eval, acopy, wk)

    int        *n;

    int        *neval;

    Mfloat     *a;

    Mint       *lda;

    Mint       *ncoda;

    Mint       *small;

    Mfloat      eval[];

    Mfloat     *acopy;

    Mfloat      wk[];

#endif

{

#define A(I_,J_)	(*(a+(I_)*(*lda)+(J_)))

#define ACOPY(I_,J_)	(*(acopy+(I_)*(*ncoda + 1)+(J_)))

    Mint        _l0, _l1, _l2, i, i_, id, idef, ie, ie2;

    Mfloat      eps1, evec[1][1];





    imsl_e1psh ("l_e3asb");



    /* Copy A to ACOPY */

    imsl_crbrb (n, a, lda, ADR (_l0, 0), ncoda, acopy, ADR (_l1, *ncoda + 1),

	ADR (_l2, 0), ncoda);

    /*

     * Reduce to symmetric tridiagonal matrix

     */

    id = 1;

    ie = id + *n;

    ie2 = ie + *n;

    l_e3csb (n, acopy, ADR (_l0, *ncoda + 1), ncoda, &wk[id - 1], &wk[ie - 1],

	ADR (_l1, FALSE), (Mfloat *) evec, ADR (_l2, 1));

    /* Find eigenvalues */

    eps1 = 0.0e0;

    idef = 0;

    l_e5asf (n, neval, &eps1, &wk[id - 1], &wk[ie - 1], &wk[ie2 - 1],

	small, &idef);

    imsl_scopy (*neval, &wk[id - 1], 1, eval, 1);

    /*

     * Sort eigenvalues into ascending order of magnitude.

     */

    imsl_svrbn (neval, eval, wk);

    /*

     * Resort into descending order of magnitude

     */

    for (i = 1; i <= *neval; i++) {

	i_ = i - 1;

	eval[i_] = wk[*neval - i];

    }



    imsl_e1pop ("l_e3asb");

    return;

}

#undef	ACOPY

#undef	A



/*Translated by FOR_C, v3.4 (P), on 09/16/93 at 08:19:58 */

/*FOR_C Options SET: c do=r ftn=ln io=p ndnt=4 op=aimn pf=,/home/elpaso/clib/include/imsl_int.h s=dvowa str=ln - prototypes */

/* Structured by FOR_STRUCT, v2.0, on 09/16/93 at 08:19:56

 *  Options SET: fmt=t s=n

 *-----------------------------------------------------------------------

 *  IMSL Name:  E5ASF/DE5ASF (Single/Double precision version)

 *

 *  Computer:   FORC/SINGLE

 *

 *  Revised:    July 1, 1990

 *

 *  Purpose:    Compute the smallest or largest eigenvalues of a

 *              symmetric tridiagonal matrix.

 *

 *  Usage:      CALL E5ASF (N, NEVAL, EPS1, D, E, E2, SMALL, IDEF)

 *

 *  Arguments:

 *     N      - Order of the matrix.  (Input)

 *     NEVAL  - Number of eigenvalues to be computed.  (Input)

 *     EPS1   - Tolerance for the theoretical absolute error of the

 *              computed eigenvalues.  (Input)

 *     D      - Real vector of length N.  On input, the diagonal elements

 *              of the matrix.  On output, if SMALL=.TRUE. the

 *              eigenvalues are in increasing order.  Otherwise they are

 *              in decreasing order.  (Input/Output)

 *     E      - Real vector of length N.  On input, the elements of the

 *              off diagonal.  E(1) is arbitrary.  On output, E is

 *              destroyed.  (Input/Output)

 *     E2     - Work vector of length N.  (Output)

 *     SMALL  - Logical variable.  (Input)

 *              If .TRUE. the smallest NEVAL eigenvalues are computed.

 *              If .FALSE. the largest NEVAL eigenvalues are computed.

 *     IDEF   - Flag indicating if the the matrix is definite.  (Input)

 *              If IDEF=+1 then the matrix is assumed to be positive

 *              definite. If IDEF=-1 then it is assumed to be negative

 *              definite.  If IDEF=0 then no such assumption is made.

 *

 *  Remark:

 *     This routine is based on the EISPACK routine RATQR.

 *

 *  Chapter:    MATH/LIBRARY Eigensystem Analysis

 *

 *  Copyright:  1990 by IMSL, Inc.  All Rights Reserved.

 *

 *  Warranty:   IMSL warrants only that IMSL testing has been applied

 *              to this code.  No other warranty, expressed or implied,

 *              is applicable.

 *

 *-----------------------------------------------------------------------

 * */

#ifdef ANSI

static void l_e5asf (Mint *n, Mint *neval, Mfloat *eps1,

                Mfloat d[], Mfloat e[], Mfloat e2[], Mint *small,

                Mint *idef)

#else

static void l_e5asf (n, neval, eps1, d, e, e2, small, idef)

    int        *n;

    int        *neval;

    Mfloat     *eps1;

    Mfloat      d[];

    Mfloat      e[];

    Mfloat      e2[];

    Mint       *small;

    Mint       *idef;

#endif

{

    Mint        conv;

    Mint        i, i_, j, j_, jdef, k, k_, l, l_;

    Mfloat      delta, ep, eps, err, f, p, q, qp, r, s, tot;

    imsl_e1psh ("E5ASF ");



    eps = imsl_amach (4);

    conv = FALSE;

    /* Set E2 to E squared */

    for (i = 2; i <= *n; i++) {

	i_ = i - 1;

	e2[i_] = e[i_] * e[i_];

    }

    jdef = *idef;



    if (!*small) {

	/*

	 * Negate elements of D for largest values

	 */

	sscal (*n, -1.0, d, 1);

	jdef = -jdef;

    }

    err = 0.0e0;

    s = 0.0e0;

    /*

     * Look for small sub-diagonal entries and define intial shift from lower

     * GERSCHGORIN bound. Copy E2 array into E2

     */

    tot = d[0];

    q = 0.0e0;

    j = 0;



    for (i = 1; i <= *n; i++) {

	i_ = i - 1;

	p = q;

	if (i == 1) {

	    e2[i_] = 0.0e0;

	}

	else {

	    if (p <= eps * (fabs (d[i_]) + fabs (d[i_ - 1])))

		e2[i_] = 0.0e0;

	}

	/*

	 * Count also if element of E2 has underflowed

	 */

	if (e2[i_] <= imsl_amach (1))

	    j += 1;

	q = 0.0e0;

	if (i != *n)

	    q = fabs (e[i_ + 1]);

	tot = imsl_f_min (d[i_] - p - q, tot);

    }



    if (jdef == 1 && tot < 0.0e0) {

	tot = 0.0e0;

    }

    else {

	sadd (*n, -tot, d, 1);

    }



    for (k = 1; k <= *neval; k++) {

	k_ = k - 1;

	/* Next QR transformation */

L_30:

	tot += s;

	delta = d[*n - 1] - s;

	i = *n;

	f = eps * fabs (tot);

	if (*eps1 < f)

	    *eps1 = f;

	if (delta <= *eps1) {

	    if (delta < (-*eps1)) {

		imsl_e1sti (1, *idef);

		if (*small) {

		    imsl_e1stl (1, "positive");

		}

		else {

		    imsl_e1stl (1, "negative");

		}



		imsl_e1mes (4, 1, "The argument IDEF = %(i1) and the matrix is not %(l1) definite.");

		goto L_9000;

	    }

	    goto L_70;

	}

	/*

	 * Replace small sub-diagonal squares by zero to reduce the incidence

	 * of underflows

	 */

	if (k != *n) {

	    for (j = k + 1; j <= *n; j++) {

		j_ = j - 1;

#if 0

		if (e2[j_] <= powi (eps * fabs (d[j_] + d[j_ - 1]), 2))

		    e2[j_] = 0.0e0;

#endif

		if (e2[j_] <= (eps * fabs (d[j_] + d[j_ - 1]) * (eps * fabs (d[j_] + d[j_ - 1]))))

		    e2[j_] = 0.0e0;

	    }

	}



	f = e2[*n - 1] / delta;

	qp = delta + f;

	p = 1.0e0;

	for (i = *n - 1; i >= k; i--) {

	    i_ = i - 1;

	    q = d[i_] - s - f;

	    r = q / qp;

	    p = p * r + 1.0e0;

	    ep = f * r;

	    d[i_ + 1] = qp + ep;

	    delta = q - ep;

	    if (delta <= *eps1) {

		if (delta < (-*eps1)) {

		    imsl_e1sti (1, *idef);

		    if (*small) {

			imsl_e1stl (1, "positive");

		    }

		    else {

			imsl_e1stl (1, "negative");

		    }



		    imsl_e1mes (4, 1, "The argument IDEF = %(i1) and the matrix is not %(l1) definite.");

		    goto L_9000;

		}

		goto L_70;

	    }

	    f = e2[i_] / q;

	    qp = delta + f;

	    e2[i_ + 1] = qp * ep;

	}



	d[k_] = qp;

	s = qp / p;

	if (tot + s > tot)

	    goto L_30;

	/*

	 * Set error -- irregular end of iteration

	 */

	if (!conv) {

	    imsl_e1sti (1, k);

/*

	    imsl_e1mes (3, 1, "The iteration for an eigenvalue failed to converge.  The best estimate will be returned."); */

	    imsl_ermes (IMSL_WARNING, IMSL_BEST_ESTIMATE_RETURNED);

	    conv = TRUE;

	}

	/* Deflate minimum diagonal element */

	s = 0.0e0;

	delta = qp;



	for (j = k; j <= *n; j++) {

	    j_ = j - 1;

	    if (d[j_] <= delta) {

		i = j;

		delta = d[j_];

	    }

	}

	/* Convergence */

L_70:

	if (i < *n)

	    e2[i] = e2[i - 1] * f / qp;



	for (j = i - 1; j >= k; j--) {

	    j_ = j - 1;

	    d[j_ + 1] = d[j_] - s;

	}

	/*

	 * Scopy was replaced with a do loop because with optimizing

	 * compilers the overlapping of storage is a problem. call imsl_scopy

	 * (i-k, e2(k), -1, e2(k+1), -1)

	 */

	for (l = i - k; l >= 1; l--) {

	    l_ = l - 1;

	    e2[k + l_] = e2[k + l_ - 1];

	}



	d[k_] = tot;

	err += fabs (delta);

	e2[k_] = err;

    }

    /*

     * Negate elements of D for largest values

     */

    if (!*small)

	sscal (*n, -1.0, d, 1);



L_9000:

    imsl_e1pop ("E5ASF ");

    return;

}

