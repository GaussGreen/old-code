#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK    l_zeros_sys_eqn(void (*fcn)(Mint, Mfloat[], Mfloat[]),
                                  Mint n, va_list argptr);
static void l_n2qnf(void (*fcn) (Mint, Mfloat[], Mfloat[]), Mfloat*,
                    Mint*, Mint*, Mfloat[], Mfloat[], Mfloat*, Mfloat[],
                    Mfloat*, Mfloat[], Mfloat[], Mfloat[]);
static void l_n2qnj(void (*fcn)(Mint, Mfloat[], Mfloat[]),
                    void (*lsjac)(Mint, Mfloat[], Mfloat[]), Mfloat*,
                    Mint*, Mint*, Mfloat[], Mfloat[], Mfloat*, Mfloat[],
                    Mfloat*, Mfloat[], Mfloat[], Mfloat[]);
static void l_n3qnf(void (*fcn) (Mint, Mfloat[], Mfloat[]), Mfloat*,
                    Mint*, Mfloat[], Mfloat[], Mfloat*, Mfloat[],
                    Mfloat[], Mint*, Mint*, Mint*, Mfloat*, Mint*,
                    Mfloat*, Mint*, Mint*, Mint*, Mint*, Mfloat[],
                    Mfloat[], Mfloat[], Mfloat[], Mfloat[]);
static void l_n3qnj(void (*fcn)(Mint, Mfloat[], Mfloat[]),
                    void (*lsjac)( Mint, Mfloat[], Mfloat[]), Mfloat *,
                    Mint*, Mfloat[], Mfloat[], Mfloat*, Mfloat[],
                    Mfloat[], Mint*, Mint*, Mint*, Mint*, Mfloat*,
                    Mint*, Mint*, Mint*, Mint*, Mfloat[], Mfloat[],
                    Mfloat[], Mfloat[], Mfloat[]);
static void l_n4qnf(void (*fcn) (Mint, Mfloat[], Mfloat[]), Mint*,
                    Mfloat[], Mfloat[], Mfloat*, Mint*, Mint*, Mint*,
                    Mfloat*, Mfloat[], Mfloat[]);
static void l_n5qnf(Mint*, Mint*, Mfloat*, Mlong*, Mint[], Mint*,
                    Mfloat[], Mfloat[], Mfloat[]);
static void l_n6qnf(Mint*, Mint*, Mfloat*, Mfloat[]);
static void l_n7qnf(Mint*, Mfloat[], Mint*, Mfloat[], Mfloat[],
                    Mfloat[], Mfloat[], Mfloat[], Mfloat[]);
static void l_n8qnf(Mint*, Mint*, Mfloat[], Mfloat[], Mfloat[],
                    Mfloat[], Mlong*);
static void l_n9qnf(Mint*, Mint*, Mfloat*, Mint*, Mfloat[], Mfloat[]);
#else
static VA_LIST_HACK      l_zeros_sys_eqn();
static void l_n2qnf();
static void l_n2qnj();
static void l_n3qnf();
static void l_n3qnj();
static void l_n4qnf();
static void l_n5qnf();
static void l_n6qnf();
static void l_n7qnf();
static void l_n8qnf();
static void l_n9qnf();
#endif

static Mfloat       *lv_value;

#ifdef ANSI
#if defined(COMPUTER_HP97C)
Mfloat *imsl_f_zeros_sys_eqn(void (*fcn) (Mint, Mfloat*, Mfloat*),
                             Mint n, ...)
#else
Mfloat *imsl_f_zeros_sys_eqn(void (*fcn) (Mint, Mfloat[], Mfloat[]),
                             Mint n, ...)
#endif
#else
Mfloat *imsl_f_zeros_sys_eqn(fcn, n, va_alist)
    void (*fcn) ();
    Mint n;
    va_dcl
#endif
{
    va_list argptr;

    VA_START(argptr, n);
    E1PSH("imsl_f_zeros_sys_eqn", "imsl_d_zeros_sys_eqn");
    lv_value = NULL;
    IMSL_CALL(l_zeros_sys_eqn(fcn, n, argptr));
    va_end(argptr);
    E1POP("imsl_f_zeros_sys_eqn", "imsl_d_zeros_sys_eqn"); 
    return lv_value;
}


#ifdef ANSI
#if defined(COMPUTER_HP97C)
static VA_LIST_HACK l_zeros_sys_eqn(void (*fcn) (Mint, Mfloat*, Mfloat*),
                               Mint n, va_list argptr)
#else
static VA_LIST_HACK l_zeros_sys_eqn(void (*fcn) (Mint, Mfloat[], Mfloat[]),
                               Mint n, va_list argptr)
#endif
#else
static VA_LIST_HACK l_zeros_sys_eqn(fcn, n, argptr)
    void (*fcn) ();
    Mint n;
    va_list       argptr;
#endif
{
    Mint          i, k, code;
    Mint          arg_number    = 2;
    Mint          maxitn        = 200;
    Mfloat        *xguess       = NULL;
    Mfloat        *xguess_float = NULL;
    Mfloat        *x_float      = NULL;
    Mfloat        err_rel;
    Mfloat        sqrteps;
    Mint          user_xguess   = 0;
    Mfloat        f;
    Mfloat        *fnorm        = NULL;
    Mfloat        *fvec         = NULL;
    Mfloat        *r            = NULL;
    Mfloat        *qtf          = NULL;
    Mfloat        *wk           = NULL;
    Mfloat        *fjac         = NULL;
    void          (*lsjac)();
    Mint          user_jacobian = 0;
    Mint          return_user   = 0;
    Mint          fnorm_user    = 0;

    sqrteps  = sqrt(imsl_amach(4));
    err_rel  = sqrteps;

    code = 1;
    while (code > 0) {
        code = va_arg(argptr, Mint);
        arg_number++;
        switch(code) {
            case IMSL_XGUESS:
                arg_number++;
                xguess = va_arg(argptr, Mfloat*);
                user_xguess = 1;
                break;
            case IMSL_JACOBIAN:
                arg_number++;
                lsjac = va_arg(argptr, void*);
                user_jacobian = 1;
                break;
            case IMSL_ERR_REL:
                arg_number++; 
                err_rel = (Mfloat) va_arg(argptr, Mdouble);
                break;
            case IMSL_ERR_REL_ADR:
                arg_number++; 
                err_rel = *(va_arg(argptr, Mfloat *));
                break;
            case IMSL_MAX_ITN:
                arg_number++; 
                maxitn = va_arg(argptr, Mint);
                break;
            case IMSL_RETURN_USER:
                arg_number++;
                lv_value = va_arg(argptr, Mfloat*);
                return_user = 1;
                break;
            case IMSL_FNORM:
                arg_number++;
                fnorm = va_arg(argptr, Mfloat*);
                fnorm_user = 1;
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

    if (n < 1) {
        imsl_e1sti(1, n);
        imsl_ermes(IMSL_TERMINAL, IMSL_NUMBER_EQUN_UNKN_LT_1); 
    }

    if (err_rel < 0) {
        imsl_e1str(1, err_rel);
        imsl_ermes(IMSL_TERMINAL, IMSL_ERRREL_LESS_THAN_ZERO);
    }

    if (maxitn < 1) {
        imsl_e1sti(1, maxitn);
        imsl_ermes(IMSL_TERMINAL, IMSL_ITMAX_LESS_THAN_ZERO); 
    }

    if (imsl_n1rty(0)) goto RETURN;

    fvec    = (Mfloat *) imsl_malloc (n*sizeof(*fvec));
    fjac    = (Mfloat *) imsl_malloc (n*n*sizeof(*fjac));
    r       = (Mfloat *) imsl_malloc ((n*(n+1)/2)*sizeof(*r));
    qtf     = (Mfloat *) imsl_malloc (n*sizeof(*qtf));
    wk      = (Mfloat *) imsl_malloc (5*n*sizeof(*wk));
    x_float = (Mfloat *) imsl_malloc (n*sizeof(*x_float));
    xguess_float = (Mfloat *) imsl_malloc (n*sizeof(*xguess_float)); 

    if (!user_xguess) 
       for (i=0; i<n; i++)  xguess_float[i] = F_ZERO;
    else
       for (i=0; i<n; i++)  xguess_float[i] = (Mfloat) xguess[i];

    if (lv_value == NULL) {
       lv_value = (Mfloat *) imsl_malloc (n*sizeof(*lv_value));
       if (lv_value == NULL) {
          imsl_e1sti(1, n);
          imsl_e1stl(1, "n");
          imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
          goto FREE_SPACE;
       }
    }
    if (fvec == NULL || r == NULL || qtf == NULL || wk == NULL ||
        fjac == NULL || x_float == NULL || xguess_float == NULL) {
          imsl_e1sti(1, n);
          imsl_e1stl(1, "n");
          imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
          goto FREE_SPACE;
    }

    if (user_jacobian)
       l_n2qnj(fcn, lsjac, &err_rel, &n, &maxitn, xguess_float, 
               x_float, &f, fvec, fjac, r, qtf, wk); 
    else
       l_n2qnf(fcn, &err_rel, &n, &maxitn, xguess_float, x_float,
               &f, fvec, fjac, r, qtf, wk);

    for (i=0; i<n; i++) {
        lv_value[i] = x_float[i];
    }
    
    if (fnorm_user) *fnorm = f;

FREE_SPACE:
    if (fvec != NULL)          imsl_free (fvec);
    if (qtf != NULL)           imsl_free (qtf);
    if (wk != NULL)            imsl_free (wk);
    if (x_float != NULL)       imsl_free (x_float);
    if (xguess_float != NULL)  imsl_free (xguess_float);
    if (fjac != NULL)          imsl_free (fjac);
    if (r != NULL)             imsl_free (r);

RETURN:
    if (imsl_n1rty(0) == 5) {
        if (!return_user && lv_value != NULL)  imsl_free(lv_value);
        lv_value = NULL;
    }
    return (argptr);
}


/* -----------------------------------------------------------------------
    IMSL Name:  N2QNF/DN2QNF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 16, 1985

    Purpose:    Solve a system of nonlinear equations using the
                Levenberg-Marquardt algorithm and a finite difference
                Jacobian.

    Usage:      CALL N2QNF (FCN, ERRREL, N, ITMAX, XINIT, X,
                            FNORM, FVEC, FJAC, R, QTF, WK)

    Arguments:  (See NEQNF)

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_n2qnf(void (*fcn) (Mint, Mfloat*, Mfloat*), Mfloat *errrel,
                    Mint *n, Mint *itmax, Mfloat xinit[], Mfloat x[], 
                    Mfloat *fnorm, Mfloat fvec[], Mfloat *fjac, Mfloat r[],
                    Mfloat qtf[], Mfloat wk[])
#else
static void l_n2qnf(void (*fcn) (Mint, Mfloat[], Mfloat[]), Mfloat *errrel,
                    Mint *n, Mint *itmax, Mfloat xinit[], Mfloat x[], 
                    Mfloat *fnorm, Mfloat fvec[], Mfloat *fjac, Mfloat r[],
                    Mfloat qtf[], Mfloat wk[])
#endif
#else
static void l_n2qnf(fcn, errrel, n, itmax, xinit, x, fnorm, fvec,
	   fjac, r, qtf, wk)
	void            (*fcn) ();
	Mfloat          *errrel;
	Mint            *n, *itmax;
	Mfloat          xinit[], x[], *fnorm, fvec[], *fjac, r[], qtf[],
	                wk[];
#endif
{
#define FJAC(I_,J_)	(fjac+(I_)*(*n)+(J_))
	Mint             info, iter, lr, maxfev, ml, mode, mu, nfev, nprint;
	Mfloat           epsfcn;
	static Mfloat    factor = 1.0e2;



	imsl_e1psh("l_n2qnf");

	if (imsl_n1rty(0) != 0)
		goto L_9000;

	info = 0;
	maxfev = *itmax * (*n + 1);
	ml = *n - 1;
	mu = *n - 1;
	epsfcn = F_ZERO;
	mode = 2;
	sset(*n, F_ONE, wk, 1);
	nprint = 0;
	lr = (*n * (*n + 1)) / 2;
	/* Copy initial guesses into X */
	scopy(*n, xinit, 1, x, 1);

	l_n3qnf(fcn, errrel, n, x, fvec, fjac, r, qtf, &maxfev, &ml, &mu,
		&epsfcn, &mode, &factor, &nprint, &info, &nfev, &lr, &wk[0],
		   &wk[*n], &wk[*n * 2], &wk[*n * 3], &wk[*n * 4]);

	if (info == 5)
		info = 4;
	*fnorm = imsl_sdot(*n, fvec, 1, fvec, 1);

	if (info == 2) {
		iter = *itmax * (*n + 1);
		imsl_e1sti(1, iter);
/*		(4, 1, "The number of calls to the function has    */
/*                      exceeded ITMAX*(N+1) = %(i1).  The user may */
/*                      try a new initial guess.");                 */
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_FCN_EVALS);
	} else if (info == 3) {
		imsl_e1str(1, *errrel);
/*		(4, 2, "The bound for the relative error, ERR_REL =  */
/*                      %(r1), is too small.  No further improvement */
/*                      in the approximate solution is possible.     */
/*                      The user should increase ERR_REL.");          */
                imsl_ermes(IMSL_WARNING, IMSL_NO_BETTER_POINT);
	} else if (info == 4) {
/*		(4, 3, "The iteration has not made good progress.   */
/*                       The user may try a new initial guess.");    */
                imsl_ermes(IMSL_WARNING, IMSL_NO_PROGRESS);
	}
L_9000:
	imsl_e1pop("l_n2qnf");

	return;
}				/* end of function */
#undef  FJAC
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N3QNF/DN3QNF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 30, 1985

    Purpose:

    Usage:      CALL N3QNF (FCN, ERRREL, N, X, FVEC, FJAC, R, QTF,
                            MAXFEV, ML, MU, EPSFCN, MODE, FACTOR, NPRINT,
                            INFO, NFEV, LR, DIAG, WK1, WK2, WK3, WK4)

    Arguments:
       FCN    - A real function subroutine supplied by the
                user.  FCN must be declared EXTERNAL in the
                CALLING PROGRAM.  FCN specifies the system of
                equations to be solved and should be of the
                following form
                              SUBROUTINE FCN(X,F,N)
                              REAL X(*),F(*)
                              F(1)=
                               .
                              F(N)=
                              RETURN
                              END
                Where X is given.  FCN must not alter X.
       ERRREL - Stopping criterion.  The root is accepted if the
                relative error between two successive approximations
                to this root is within ERRREL.  (Input)
       N      - The number of equations to be solved and the number
                of unknowns.  (Input)
       X      - A vector of length N.  X contains the best estimate
                of the root found by NEQNF.  (Output)
       FVEC   - A vector of length N.  FVEC contains the functions
                evaluated at the point X.
       FJAC   - An N by N matrix.  FJAC contains the orthogonal
                matrix Q produced by the QR factorization of the
                final approximate Jacobian.
       R      - A vector of length N*(N+1)/2.  R contains the upper
                triangular matrix produced by the QR factorization
                of the final approximation Jacobian.  R is stored
                rowwise.
       QTF    - A vector of length N.  QTF contains the vector
                (Q transpose)*FVEC.
       MAXFEV - Maximum number of calls to FCN. (Input)
       ML     -
       MU     -
       EPSFCN -
       MODE   -
       FACTOR -
       NPRINT - Number of iterates to be printed. (
       INFO   -
       NFEV   - Number of calls to FCN.  (Output)
       LR     - Length of the vector R.  (Input)
       DIAG   -
       WK1    - Real work vector of length N.  (Output)
       WK2    - Real work vector of length N.  (Output)
       WK3    - Real work vector of length N.  (Output)
       WK4    - Real work vector of length N.  (Output)

    Copyright:  1985 by IMSL, Inc.  All rights reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_n3qnf(void (*fcn) (Mint, Mfloat*, Mfloat*), Mfloat *errrel,
                    Mint *n, Mfloat x[], Mfloat fvec[], Mfloat *fjac,
                    Mfloat r[], Mfloat qtf[], Mint *maxfev, Mint *ml, Mint *mu,
                    Mfloat *epsfcn, Mint *mode, Mfloat *factor, Mint *nprint,
                    Mint *info, Mint *nfev, Mint *lr, Mfloat diag[],
                    Mfloat wk1[], Mfloat wk2[], Mfloat wk3[], Mfloat wk4[])
#else
static void l_n3qnf(void (*fcn) (Mint, Mfloat[], Mfloat[]), Mfloat *errrel,
                    Mint *n, Mfloat x[], Mfloat fvec[], Mfloat *fjac,
                    Mfloat r[], Mfloat qtf[], Mint *maxfev, Mint *ml, Mint *mu,
                    Mfloat *epsfcn, Mint *mode, Mfloat *factor, Mint *nprint,
                    Mint *info, Mint *nfev, Mint *lr, Mfloat diag[],
                    Mfloat wk1[], Mfloat wk2[], Mfloat wk3[], Mfloat wk4[])
#endif
#else
static void l_n3qnf(fcn, errrel, n, x, fvec, fjac, r, qtf, maxfev, ml, mu,
                    epsfcn, mode, factor, nprint, info, nfev, lr, diag, wk1,
	            wk2, wk3, wk4)
	void            (*fcn) ();
	Mfloat          *errrel;
	Mint            *n;
	Mfloat           x[], fvec[], *fjac, r[], qtf[];
	Mint            *maxfev, *ml, *mu;
	Mfloat          *epsfcn;
	Mint            *mode;
	Mfloat          *factor;
	Mint            *nprint, *info, *nfev, *lr;
	Mfloat           diag[], wk1[], wk2[], wk3[], wk4[];
#endif
{
#define FJAC(I_,J_)	(fjac+(I_)*(*n)+(J_))
	Mlong	        jeval, sing;
	Mint             _l0, _l1, i, iflag, iter, iwa[1], j, jm1, l, msum,
	                ncfail, ncsuc, nslow1, nslow2;
	Mfloat           actred, delta, epsmch, fnorm, fnorm1, pnorm, prered,
	                ratio, sum, temp, xnorm;
        Mlong            IMSLFALSE = 0, IMSLTRUE = 1;

	epsmch = imsl_amach(4);
	*info = 0;
	iflag = 0;
	*nfev = 0;
	/*
	 * Check the input parameters for errors
	 */
	if (*mode == 2) {
		for (j = 1; j <= *n; j++) {
			if (diag[j - 1] <= F_ZERO)
				goto L_150;
		}
	}
	/*
	 * Evaluate the function at the starting point and calculate its norm
	 */
	iflag = 1;
	imsl_e1usr("ON");
	(*fcn) (*n, x, fvec);
	imsl_e1usr("OFF");
	*nfev = 1;
	if (iflag < 0)
		goto L_150;
	fnorm = imsl_snrm2(*n, fvec, 1);
	/*
	 * Determine the number of calls to FCN Needed to compute the
	 * jacobian matrix
	 */
	msum = imsl_i_min(*ml + *mu + 1, *n);

	/*
	 * Initialize iteration counter and monitors
	 */
	iter = 1;
	ncsuc = 0;
	ncfail = 0;
	nslow1 = 0;
	nslow2 = 0;
	/* Beginning of the outer loop */
L_20:
	jeval = IMSLTRUE;
	/* Calculate the jacobian matrix */
	iflag = 2;
	l_n4qnf(fcn, n, x, fvec, fjac, &iflag, ml, mu, epsfcn, wk1, wk2);
	*nfev += msum;
	if (iflag < 0)
		goto L_150;
	/*
	 * Compute the QR factorization of the jacobian
	 */
	_l0 = IMSLFALSE;
	_l1 = 1;
	{	/* fix Mlong/Mint mismatch for _l0; */ 
	 	Mlong _l0_tmp = _l0;
		l_n5qnf(n, n, fjac, (Mlong *)&_l0_tmp, iwa, &_l1, wk1, wk2, wk3);
		_l0 = _l0_tmp;
	}

	/*
	 * On the first iteration and if MODE is 1, scale according to the
	 * norms of the columns of the intial Jacobian
	 */
	if (iter == 1) {
		if (*mode != 2) {
			scopy(*n, wk2, 1, diag, 1);
			for (j = 1; j <= *n; j++) {
				if (wk2[j - 1] == F_ZERO)
					diag[j - 1] = F_ONE;
			}
		}
		/*
		 * On the first iteration, calculate the norm of the scaled X
		 * and initialize the step bound delta
		 */
		for (j = 1; j <= *n; j++) {
			wk3[j - 1] = diag[j - 1] * x[j - 1];
		}
		xnorm = imsl_snrm2(*n, wk3, 1);
		delta = *factor * xnorm;
		if (delta == F_ZERO)
			delta = *factor;
	}
	/*
	 * Form (Q transpose)*FVEC and store in QTF.
	 */
	scopy(*n, fvec, 1, qtf, 1);
	for (j = 1; j <= *n; j++) {
		if (*FJAC(j - 1, j - 1) != F_ZERO) {
			sum = imsl_sdot(*n - j + 1, FJAC(j - 1, j - 1), 1, &qtf[j - 1],
					1);
			temp = -sum / *FJAC(j - 1, j - 1);
			saxpy(*n - j + 1, temp, FJAC(j - 1, j - 1), 1, &qtf[j - 1],
				   1);
		}
	}
	/*
	 * Copy the triangular factor of the QR factorization into R
	 */
	sing = IMSLFALSE;
	for (j = 1; j <= *n; j++) {
		l = j;
		jm1 = j - 1;
		if (jm1 >= 1) {
			for (i = 1; i <= jm1; i++) {
				r[l - 1] = *FJAC(j - 1, i - 1);
				l += *n - i;
			}
		}
		r[l - 1] = wk1[j - 1];
		if (wk1[j - 1] == F_ZERO)
			sing = IMSLTRUE;
	}
	/*
	 * Accumulate the orthogonal factor in FJAC
	 */
	l_n6qnf(n, n, fjac, wk1);
	/* Rescale if necessary */
	if (*mode != 2) {
		for (j = 1; j <= *n; j++) {
			diag[j - 1] = imsl_f_max(diag[j - 1], wk2[j - 1]);
		}
	}
	/*
	 * Beginning of the inner loop If requested, call FCN to enable
	 * printing of iterates
	 */
L_90:
	if (*nprint <= 0)
		goto L_100;
	if (iflag < 0)
		goto L_150;
L_100:
	;
	/* Determine the direction P */
	l_n7qnf(n, r, lr, diag, qtf, &delta, wk1, wk2, wk3);
	/*
	 * Store the direction P and X + P Calculate the norm of P
	 */
	sscal(*n, -F_ONE, wk1, 1);
	for (j = 1; j <= *n; j++) {
		wk2[j - 1] = x[j - 1] + wk1[j - 1];
		wk3[j - 1] = diag[j - 1] * wk1[j - 1];
	}
	pnorm = imsl_snrm2(*n, wk3, 1);
	/*
	 * On the first iteration, adjust the initial step bound
	 */
	if (iter == 1)
		delta = imsl_f_min(delta, pnorm);
	/*
	 * Evaluate the function at X + P and calculate its norm
	 */
	imsl_e1usr("ON");
	(*fcn) (*n, wk2, wk4);
	imsl_e1usr("OFF");
	*nfev += 1;
	fnorm1 = imsl_snrm2(*n, wk4, 1);
	/* Compute the scaled actual reduction */
	actred = -F_ONE;
	if (fnorm1 < fnorm)
		actred = F_ONE - imsl_fi_power(fnorm1 / fnorm, 2);
	/*
	 * Compute the scaled predicted reduction
	 */
	l = 1;
	for (i = 1; i <= *n; i++) {
		sum = imsl_sdot(*n - i + 1, &r[l - 1], 1, &wk1[i - 1], 1);
		l += *n - i + 1;
		wk3[i - 1] = qtf[i - 1] + sum;
	}
	temp = imsl_snrm2(*n, wk3, 1);
	prered = F_ONE;
	if (temp < fnorm)
		prered = F_ONE - imsl_fi_power(temp / fnorm, 2);
	/*
	 * Compute the ratio of the actual to the prdeicted reduction
	 */
	ratio = F_ZERO;
	if (prered > F_ZERO)
		ratio = actred / prered;
	/* Update the step bound */
	if (ratio >= 0.1) {
		ncfail = 0;
		ncsuc += 1;
		if (ratio >= F_HALF || ncsuc > 1)
			delta = imsl_f_max(delta, pnorm / F_HALF);
		if (fabs(ratio - F_ONE) <= 0.1)
			delta = pnorm / F_HALF;
	} else {
		ncsuc = 0;
		ncfail += 1;
		delta *= F_HALF;
	}
	/* Test for successful iteration */
	if (ratio >= 0.0001) {
		/*
		 * Successful iteration. Update X, FVEC, and their norms
		 */
		scopy(*n, wk2, 1, x, 1);
		scopy(*n, wk4, 1, fvec, 1);
		for (j = 1; j <= *n; j++) {
			wk2[j - 1] = diag[j - 1] * x[j - 1];
		}
		xnorm = imsl_snrm2(*n, wk2, 1);
		fnorm = fnorm1;
		iter += 1;
	}
	/*
	 * Determine the progress of the iteration
	 */
	nslow1 += 1;
	if (actred >= 0.001)
		nslow1 = 0;
	if (jeval)
		nslow2 += 1;
	if (actred >= 0.1)
		nslow2 = 0;
	/* Test for convergence */
	if (delta <= *errrel * xnorm || fnorm == F_ZERO)
		*info = 1;
	if (*info != 0)
		goto L_150;
	/*
	 * Tests for termination and stringent tolerances
	 */
	if (*nfev >= *maxfev)
		*info = 2;
	if (0.1 * imsl_f_max(0.1 * delta, pnorm) <= epsmch * xnorm)
		*info = 3;
	if (nslow2 == 5)
		*info = 4;
	if (nslow1 == 10)
		*info = 5;
	if (*info != 0)
		goto L_150;
	/*
	 * Criterion for recalculating Jacobian approximation by forward
	 * differences
	 */
	if (ncfail != 2) {
		/*
		 * Calculate the rank one modification to the Jacobian and
		 * update QTF if necessary
		 */
		for (j = 1; j <= *n; j++) {
			sum = imsl_sdot(*n, FJAC(j - 1, 0), 1, wk4, 1);
			wk2[j - 1] = (sum - wk3[j - 1]) / pnorm;
			wk1[j - 1] = diag[j - 1] * ((diag[j - 1] * wk1[j - 1]) / pnorm);
			if (ratio >= 0.0001)
				qtf[j - 1] = sum;
		}
		/*
		 * Compute the QR factoracation of the updated Jacobian
		 */
		l_n8qnf(n, n, r, wk1, wk2, wk3, &sing);
		l_n9qnf(n, n, fjac, n, wk2, wk3);
                _l0 = 1;
                _l1 = 1;
		l_n9qnf(&_l0, n, qtf, &_l1, wk2, wk3);
		/* End of the inner loop */
		jeval = IMSLFALSE;
		goto L_90;
	}
	/* End of the outer loop */
	goto L_20;
	/*
	 * Termination, either normal or user imposed
	 */
L_150:
	if (iflag < 0)
		*info = iflag;
	iflag = 0;
	return;
}				/* end of function */
#undef  FJAC
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N4QNF/DN4QNF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1985

    Purpose:    Calculate the Jacobian in order to solve a system of
                nonlinear equations

    Usage:      CALL N4QNF (FCN, N, X, FVEC, FJAC, IFLAG, ML, MU, EPSFCN,
                            WK1, WK2)

    Arguments:
       FCN    - A real function subroutine supplied by the user.  FCN
                must be declared EXTERNAL in the CALLING PROGRAM.
                FCN specifies the system of equations to be solved and
                should be of the following form
                               SUBROUTINE FCN(X,F,N)
                               REAL X(*),F(*)
                               F(1)=
                               .
                               F(N)=
                               RETURN
                               END
                Where X is given.  FCN must not alter X.
       N      - The number of equations to be solbed and the number
                of unknowns.  (Input)
       X      - A vector of length N.  X contains the best estimate
                of the root found by NEQNF.  (Output)
       FVEC   - A vector of length N.  FVEC contains the functions
                evaluated at the point X.  (Input)
       FJAC   - An N by N matrix.  FJAC contains the orthogonal
                matrix Q produced by the QR factorization of the
                final approximate Jacobian.  (Output)
       IFLAG  -
       ML     -
       MU     -
       EPSFCN -
       WK1    - Real work array of length N.  (Output)
       WK2    - Real work array of length N.  (Output)

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_n4qnf(void (*fcn) (Mint, Mfloat*, Mfloat*), Mint *n,
                    Mfloat x[], Mfloat fvec[], Mfloat *fjac, Mint *iflag,
                    Mint *ml, Mint *mu, Mfloat *epsfcn, Mfloat wk1[],
                    Mfloat wk2[])
#else
static void l_n4qnf(void (*fcn) (Mint, Mfloat[], Mfloat[]), Mint *n,
                    Mfloat x[], Mfloat fvec[], Mfloat *fjac, Mint *iflag,
                    Mint *ml, Mint *mu, Mfloat *epsfcn, Mfloat wk1[],
                    Mfloat wk2[])
#endif
#else
static void l_n4qnf(fcn, n, x, fvec, fjac, iflag, ml, mu, epsfcn, wk1, wk2)
	void            (*fcn) ();
	Mint            *n;
	Mfloat           x[], fvec[], *fjac;
	Mint            *iflag, *ml, *mu;
	Mfloat          *epsfcn, wk1[], wk2[];
#endif
{
#define FJAC(I_,J_)	(fjac+(I_)*(*n)+(J_))
	Mint             _d_l, _d_m, _do0, _do1, _do2, _do3, i, j, k, msum;
	Mfloat           eps, epsmch, h, temp;


	epsmch = imsl_amach(4);
	eps = sqrt(imsl_f_max(*epsfcn, epsmch));
	msum = *ml + *mu + 1;
	if (msum >= *n) {
		/*
		 * Computation of dense approxiamte JACOBIAN
		 */
		for (j = 1; j <= *n; j++) {
			temp = x[j - 1];
			h = eps * fabs(temp);
			if (h == F_ZERO)
				h = eps;
			x[j - 1] = temp + h;
			imsl_e1usr("ON");
			(*fcn) (*n, x, wk1);
			imsl_e1usr("OFF");
			if (*iflag < 0)
				goto L_9000;
			x[j - 1] = temp;
			for (i = 1; i <= *n; i++) {
				*FJAC(j - 1, i - 1) = (wk1[i - 1] - fvec[i - 1]) / h;
			}
		}
		goto L_9000;
	}
	/*
	 * Computation of banded approximate JACOBIAN
	 */
	for (k = 1; k <= msum; k++) {
		for (j = k, _do0 = DOCNT(k, *n, _do1 = msum); _do0 > 0; j += _do1, _do0--) {
			wk2[j - 1] = x[j - 1];
			h = eps * fabs(wk2[j - 1]);
			if (h == F_ZERO)
				h = eps;
			x[j - 1] = wk2[j - 1] + h;
		}
		imsl_e1usr("ON");
		(*fcn) (*n, x, wk1);
		imsl_e1usr("OFF");
		if (*iflag < 0)
			goto L_9000;
		for (j = k, _do2 = DOCNT(k, *n, _do3 = msum); _do2 > 0; j += _do3, _do2--) {
			x[j - 1] = wk2[j - 1];
			h = eps * fabs(wk2[j - 1]);
			if (h == F_ZERO)
				h = eps;
			for (i = 1; i <= *n; i++) {
				*FJAC(j - 1, i - 1) = F_ZERO;
				if (i >= j - *mu && i <= j + *ml)
					*FJAC(j - 1, i - 1) = (wk1[i - 1] - fvec[i - 1]) /
						h;
			}
		}
	}
L_9000:
	return;
}				/* end of function */
#undef   FJAC
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N5QNF/DN5QNF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1985

    Purpose:

    Usage:      CALL N5QNF (M, N, A, PIVOT, IPVT, LIPVT, RDIAG, ACNORM,
                            WK)

    Arguments:
       M      -
       N      - The number of equations to be solbed and the number
                of unknowns.  (Input)
       A      -
       PIVOT  -
       IPVT   -
       LIPVT  -
       RDIAG  -
       ACNORM -
       WK     - Real work array of length N.  (Output)

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* LIPVT is not used here, but leave the calling sequence intact. */
#ifdef ANSI
static void l_n5qnf(Mint *m, Mint *n, Mfloat *a, Mlong *pivot, Mint ipvt[],
                    Mint *lipvt, Mfloat rdiag[], Mfloat acnorm[], Mfloat wk[])
#else
static void l_n5qnf(m, n, a, pivot, ipvt, lipvt, rdiag, acnorm,
	   wk)
	Mint            *m, *n;
	Mfloat          *a;
	Mlong            *pivot;
	Mint             ipvt[], *lipvt;
	Mfloat           rdiag[], acnorm[], wk[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*n)+(J_))
	Mint             j, jp1, k, kmax, minmn;
	Mfloat           ajnorm, big, epsmch, small, sum, temp;


	epsmch = imsl_amach(4);
	/*
	 * Compute the initial column norms and initialize several arrays
	 */
	for (j = 1; j <= *n; j++) {
		acnorm[j - 1] = imsl_snrm2(*m, A(j - 1, 0), 1);
		if (*pivot)
			ipvt[j - 1] = j;
	}
	scopy(*n, acnorm, 1, rdiag, 1);
	scopy(*n, rdiag, 1, wk, 1);
	/*
	 * Reduce A to R with Householder transformations
	 */
	minmn = imsl_i_min(*m, *n);
	for (j = 1; j <= minmn; j++) {
		if (*pivot) {
			/*
			 * Bring the column of largest norm into the pivot
			 * position
			 */
			kmax = j;
			for (k = j; k <= *n; k++) {
				if (rdiag[k - 1] > rdiag[kmax - 1])
					kmax = k;
			}
			if (kmax != j) {
				sswap(*m, A(j - 1, 0), 1, A(kmax - 1, 0), 1);
				rdiag[kmax - 1] = rdiag[j - 1];
				wk[kmax - 1] = wk[j - 1];
				k = ipvt[j - 1];
				ipvt[j - 1] = ipvt[kmax - 1];
				ipvt[kmax - 1] = k;
			}
		}
		/*
		 * Compute the Householder transformation to reduce the J-TH
		 * column of A to a multiple of the J-TH unit vector
		 */
		ajnorm = imsl_snrm2(*m - j + 1, A(j - 1, j - 1), 1);
		small = imsl_amach(1);
		big = imsl_amach(2);
		if (small * big < F_ONE)
			small = F_ONE / big;
		if (ajnorm != F_ZERO) {
			if (*A(j - 1, j - 1) < F_ZERO)
				ajnorm = -ajnorm;
			sscal(*m - j + 1, F_ONE / ajnorm, A(j - 1, j - 1), 1);
			*A(j - 1, j - 1) += F_ONE;
			/*
			 * Apply the transformation to the remaining columns
			 * and update the norms
			 */
			jp1 = j + 1;
			if (*n >= jp1) {
				for (k = jp1; k <= *n; k++) {
					sum = imsl_sdot(*m - j + 1, A(j - 1, j - 1), 1, A(k - 1, j - 1),
							1);
					temp = sum / *A(j - 1, j - 1);
					saxpy(*m - j + 1, -temp, A(j - 1, j - 1), 1, A(k - 1, j - 1),
						   1);
					if (*pivot && rdiag[k - 1] != F_ZERO) {
						temp = *A(k - 1, j - 1) / rdiag[k - 1];
						rdiag[k - 1] *= sqrt(imsl_f_max(F_ZERO, F_ONE - imsl_fi_power(temp, 2)));
						if (0.05 * imsl_fi_power(rdiag[k - 1] / wk[k - 1], 2) <=
						    epsmch) {
							rdiag[k - 1] = imsl_snrm2(*m - j, A(k - 1, jp1 - 1),
									 1);
							wk[k - 1] = rdiag[k - 1];
						}
					}
				}
			}
		}
		rdiag[j - 1] = -ajnorm;
	}
	return;
}				/* end of function */
#undef  A
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N6QNF/DN6QNF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1985

    Purpose:    Calculate the orthogonal factor in to solve a system of
                nonlinear equations

    Usage:      CALL N6QNF (M, N, Q, WK)

    Arguments:
       M      -
       N      - The number of equations to be solbed and the number
                of unknowns.  (Input)
       Q      -
       WK     - Real work array of length N. (

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_n6qnf(Mint *m, Mint *n, Mfloat *q, Mfloat wk[])
#else
static void l_n6qnf(m, n, q, wk)
	Mint            *m, *n;
	Mfloat          *q, wk[];
#endif
{
#define Q(I_,J_)	(q+(I_)*(*n)+(J_))
	Mint             j, jm1, k, l, minmn, np1;
	Mfloat           big, small, sum, temp;


	minmn = imsl_i_min(*m, *n);
	if (minmn >= 2) {
		for (j = 2; j <= minmn; j++) {
			jm1 = j - 1;
			sset(jm1, F_ZERO, Q(j - 1, 0), 1);
		}
	}
	/*
	 * Initialize remaining columns to those of the identity matrix
	 */
	np1 = *n + 1;
	if (*m >= np1) {
		for (j = np1; j <= *m; j++) {
			sset(*m, F_ZERO, Q(j - 1, 0), 1);
			*Q(j - 1, j - 1) = F_ONE;
		}
	}
	/* Accumulate Q from its factored form */
	for (l = 1; l <= minmn; l++) {
		k = minmn - l + 1;
		scopy(*m - k + 1, Q(k - 1, k - 1), 1, &wk[k - 1], 1);
		sset(*m - k + 1, F_ZERO, Q(k - 1, k - 1), 1);
		*Q(k - 1, k - 1) = F_ONE;
		small = imsl_amach(1);
		big = imsl_amach(2);
		if (big * small < F_ONE)
			small = F_ONE / big;
		if (wk[k - 1] != F_ZERO) {
			for (j = k; j <= *m; j++) {
				sum = imsl_sdot(*m - k + 1, Q(j - 1, k - 1), 1, &wk[k - 1],
						1);
				temp = sum / wk[k - 1];
				saxpy(*m - k + 1, -temp, &wk[k - 1], 1, Q(j - 1, k - 1),
					   1);
			}
		}
	}
	return;
}				/* end of function */
#undef  Q
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N7QNF/DN7QNF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1985

    Purpose:

    Usage:      CALL N7QNF (N, R, LR, DIAG, QTB, DELTA, X, WK1, WK2)

    Arguments:
       N      - The number of equations to be solbed and the number
                of unknowns.  (Input)
       R      -
       LR     -
       DIAG   -
       QTB    -
       DELTA  -
       X      -
       WK1    - Real work array of length N.
       WK2    - Real work array of length N.

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* LR is not used, but leave the calling sequence intact. */
#ifdef ANSI
static void l_n7qnf(Mint *n, Mfloat r[], Mint *lr, Mfloat diag[], Mfloat qtb[],
                    Mfloat delta[], Mfloat x[], Mfloat wk1[], Mfloat wk2[])
#else
static void l_n7qnf(n, r, lr, diag, qtb, delta, x, wk1, wk2)
	Mint            *n;
	Mfloat           r[];
	Mint            *lr;
	Mfloat           diag[], qtb[], *delta, x[], wk1[], wk2[];
#endif
{
	Mint             i, j, jj, jp1, k, l;
	Mfloat           alpha, big, bnorm, epsmch, gnorm, qnorm, sgnorm,
	                small, sum, temp;


	epsmch = imsl_amach(4);
	/*
	 * First, calculate the GAUSS-NEWTON direction
	 */
	jj = (*n * (*n + 1)) / 2 + 1;
	for (k = 1; k <= *n; k++) {
		j = *n - k + 1;
		jp1 = j + 1;
		jj -= k;
		l = jj + 1;
		sum = imsl_sdot(*n - jp1 + 1, &r[l - 1], 1, &x[jp1 - 1], 1);
		temp = r[jj - 1];
		if (temp == F_ZERO) {
			l = j;
			for (i = 1; i <= j; i++) {
				temp = imsl_f_max(temp, fabs(r[l - 1]));
				l += *n - i;
			}
			temp *= epsmch;
			if (temp == F_ZERO)
				temp = epsmch;
		}
		x[j - 1] = (qtb[j - 1] - sum) / temp;
	}
	/*
	 * Test whether the GAUSS-NEWTON direction is acceptable
	 */
	sset(*n, F_ZERO, wk1, 1);
	for (j = 1; j <= *n; j++) {
		wk2[j - 1] = diag[j - 1] * x[j - 1];
	}
	qnorm = imsl_snrm2(*n, wk2, 1);
	if (qnorm > *delta) {
		/*
		 * The GAUSS-NEWTON direction is not acceptabel. Next,
		 * calculate the scaled gradient direction
		 */
		l = 1;
		for (j = 1; j <= *n; j++) {
			temp = qtb[j - 1];
			saxpy(*n - j + 1, temp, &r[l - 1], 1, &wk1[j - 1], 1);
			l += *n - j + 1;
			wk1[j - 1] /= diag[j - 1];
		}
		/*
		 * Calculate the norm of the scaled gradient and test for the
		 * special case in wjich the scaled gradient is zero
		 */
		gnorm = imsl_snrm2(*n, wk1, 1);
		sgnorm = F_ZERO;
		alpha = *delta / qnorm;
		small = imsl_amach(1);
		big = imsl_amach(2);
		if (big * small < F_ONE)
			small = F_ONE / big;
		if (gnorm != F_ZERO) {
			/*
			 * Calculate the point along the scaled gradient at
			 * which the quadratic is minimized
			 */
			for (j = 1; j <= *n; j++) {
				wk1[j - 1] = (wk1[j - 1] / gnorm) / diag[j - 1];
			}
			l = 1;
			for (j = 1; j <= *n; j++) {
				sum = imsl_sdot(*n - j + 1, &r[l - 1], 1, &wk1[j - 1],
						1);
				l += *n - j + 1;
				wk2[j - 1] = sum;
			}
			temp = imsl_snrm2(*n, wk2, 1);
			sgnorm = (gnorm / temp) / temp;
			/*
			 * Test whether the scaled gradient direction is
			 * acceptable
			 */
			alpha = F_ZERO;
			if (sgnorm < *delta) {
				/*
				 * The scaled gradient direction is not
				 * acceptable. Finally, calculate the point
				 * along the dogleg at which the quadratic is
				 * minimized
				 */
				bnorm = imsl_snrm2(*n, qtb, 1);
				temp = (bnorm / gnorm) * (bnorm / qnorm) * (sgnorm / *delta);
				temp += -(*delta / qnorm) * imsl_fi_power(sgnorm / *delta, 2) +
					sqrt(imsl_fi_power(temp - (*delta / qnorm), 2) + (F_ONE - imsl_fi_power(*delta /
														   qnorm, 2)) * (F_ONE - imsl_fi_power(sgnorm / *delta, 2)));
				alpha = ((*delta / qnorm) * (F_ONE - imsl_fi_power(sgnorm / *delta, 2))) /
					temp;
			}
		}
		/*
		 * Form appropriate convex combination of the GAUSS-NEWTON
		 * direction and the scaled gradient direction
		 */
		temp = (F_ONE - alpha) * imsl_f_min(sgnorm, *delta);
		for (j = 1; j <= *n; j++) {
			x[j - 1] = temp * wk1[j - 1] + alpha * x[j - 1];
		}
	}
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N8QNF/DN8QNF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1985

    Purpose:

    Usage:      CALL N8QNF (M, N, S, U, V, W, SING)

    Arguments:
       M      -
       N      -
       S      -
       U      -
       V      -
       W      -
       SING   -

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_n8qnf(Mint *m, Mint *n, Mfloat s[], Mfloat u[], Mfloat v[],
                    Mfloat w[], Mlong *sing)
#else
static void l_n8qnf(m, n, s, u, v, w, sing)
	Mint            *m, *n;
	Mfloat           s[], u[], v[], w[];
	Mlong            *sing;
#endif
{
	Mint             i, j, jj, l, nm1, nmj;
	Mfloat           giant, tau, temp, temp1, temp2, temp3, temp4;
        Mlong             IMSLFALSE = 0, IMSLTRUE = 1;

	giant = imsl_amach(2);
	jj = (*n * (2 ** m - *n + 1)) / 2 - (*m - *n);
	/*
	 * Move the nontrivial part of the last column of S into W
	 */
	scopy(*m - *n + 1, &s[jj - 1], 1, &w[*n - 1], 1);
	/*
	 * Rotate the vector V into a multiple of the N-th unit vector in
	 * such a way that a spike is introduced into W
	 */
	nm1 = *n - 1;
	if (nm1 >= 1) {
		for (nmj = 1; nmj <= nm1; nmj++) {
			j = *n - nmj;
			jj -= *m - j + 1;
			w[j - 1] = F_ZERO;
			if (v[j - 1] != F_ZERO) {
				/*
				 * Determine a GIVENS rotation which
				 * eliminates the J-th element of V
				 */
				if (fabs(v[*n - 1]) < fabs(v[j - 1])) {
					temp2 = v[*n - 1] / v[j - 1];
					temp3 = F_HALF / sqrt(0.25 + 0.25 * imsl_fi_power(temp2, 2));
					temp1 = temp3 * temp2;
					tau = F_ONE;
					if (fabs(temp1) * giant > F_ONE)
						tau = F_ONE / temp1;
				} else {
					temp4 = v[j - 1] / v[*n - 1];
					temp1 = F_HALF / sqrt(0.25 + 0.25 * imsl_fi_power(temp4, 2));
					temp3 = temp1 * temp4;
					tau = temp3;
				}
				/*
				 * Apply the transformation to V and store
				 * the information necessary to recover the
				 * GIVENS rotation
				 */
				v[*n - 1] = temp3 * v[j - 1] + temp1 * v[*n - 1];
				v[j - 1] = tau;
				/*
				 * Apply the transformation to S and extend
				 * the spike in W
				 */
				l = jj;
				for (i = j; i <= *m; i++) {
					temp = temp1 * s[l - 1] - temp3 * w[i - 1];
					w[i - 1] = temp3 * s[l - 1] + temp1 * w[i - 1];
					s[l - 1] = temp;
					l += 1;
				}
			}
		}
	}
	/*
	 * Add the spike from the rank 1 update to W
	 */
	for (i = 1; i <= *m; i++) {
		w[i - 1] += v[*n - 1] * u[i - 1];
	}
	/* Eliminate the spike */
	*sing = IMSLFALSE;
	if (nm1 >= 1) {
		for (j = 1; j <= nm1; j++) {
			if (w[j - 1] != F_ZERO) {
				/*
				 * Determine a GIVENS rotation which
				 * eliminates the J-th element of the spike
				 */
				if (fabs(s[jj - 1]) < fabs(w[j - 1])) {
					temp2 = s[jj - 1] / w[j - 1];
					temp3 = F_HALF / sqrt(0.25 + 0.25 * imsl_fi_power(temp2, 2));
					temp1 = temp3 * temp2;
					tau = F_ONE;
					if (fabs(temp1) * giant > F_ONE)
						tau = F_ONE / temp1;
				} else {
					temp4 = w[j - 1] / s[jj - 1];
					temp1 = F_HALF / sqrt(0.25 + 0.25 * imsl_fi_power(temp4, 2));
					temp3 = temp1 * temp4;
					tau = temp3;
				}
				/*
				 * Apply the transformation to S and reduce
				 * the spike in W
				 */
				l = jj;
				for (i = j; i <= *m; i++) {
					temp = temp1 * s[l - 1] + temp3 * w[i - 1];
					w[i - 1] = -temp3 * s[l - 1] + temp1 * w[i - 1];
					s[l - 1] = temp;
					l += 1;
				}
				/*
				 * Sotre the information necessary to recover
				 * the GIVENS rotation
				 */
				w[j - 1] = tau;
			}
			/*
			 * Test for zero diagonal elements in the output S
			 */
			if (s[jj - 1] == F_ZERO)
				*sing = IMSLTRUE;
			jj += *m - j + 1;
		}
	}
	/*
	 * Move W back into the last column of the output S
	 */
	scopy(*m - *n + 1, &w[*n - 1], 1, &s[jj - 1], 1);
	if (s[jj - 1] == F_ZERO)
		*sing = IMSLTRUE;
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N9QNF/DN9QNF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1985

    Purpose:

    Usage:      CALL N9QNF (M, N, A, V, W)

    Arguments:
       M      -
       N      - The number of equations to be solved and the number
                of unknowns.  (Input)
       A      -
       V      -
       W      -

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_n9qnf(Mint *m, Mint *n, Mfloat *a, Mint *lda, Mfloat v[], 
                    Mfloat w[])
#else
static void l_n9qnf(m, n, a, lda, v, w)
	Mint            *m, *n;
	Mfloat          *a;
	Mint            *lda;
	Mfloat           v[], w[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint             i, j, nm1, nmj;
	Mfloat           temp, temp1, temp2;


	nm1 = *n - 1;
	if (nm1 >= 1) {
		for (nmj = 1; nmj <= nm1; nmj++) {
			j = *n - nmj;
			if (fabs(v[j - 1]) > F_ONE) {
				temp1 = F_ONE / v[j - 1];
				temp2 = sqrt(F_ONE - imsl_fi_power(temp1, 2));
			} else {
				temp2 = v[j - 1];
				temp1 = sqrt(F_ONE - imsl_fi_power(temp2, 2));
			}
			for (i = 1; i <= *m; i++) {
				temp = temp1 ** A(j - 1, i - 1) - temp2 ** A(*n - 1, i - 1);
				*A(*n - 1, i - 1) = temp2 ** A(j - 1, i - 1) + temp1 ** A(*n - 1, i - 1);
				*A(j - 1, i - 1) = temp;
			}
		}
		/*
		 * Apply the second set of GIVENS rotations to A
		 */
		for (j = 1; j <= nm1; j++) {
			if (fabs(w[j - 1]) > F_ONE) {
				temp1 = F_ONE / w[j - 1];
				temp2 = sqrt(F_ONE - imsl_fi_power(temp1, 2));
			} else {
				temp2 = w[j - 1];
				temp1 = sqrt(F_ONE - imsl_fi_power(temp2, 2));
			}
			for (i = 1; i <= *m; i++) {
				temp = temp1 ** A(j - 1, i - 1) + temp2 ** A(*n - 1, i - 1);
				*A(*n - 1, i - 1) = -temp2 ** A(j - 1, i - 1) + temp1 *
					*A(*n - 1, i - 1);
				*A(j - 1, i - 1) = temp;
			}
		}
	}
	return;
}				/* end of function */
#undef  A


/* -----------------------------------------------------------------------
    IMSL Name:  N2QNJ/DN2QNJ (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 16, 1985

    Purpose:    Solve a system of nonlinear equations using the
                Levenberg-Marquardt algorithm and a user-supplied
                Jacobian.

    Usage:      CALL N2QNJ (FCN, LSJAC, ERRREL, N, ITMAX, XINIT, X,
                            FNORM, FVEC, FJAC, R, QTF, WK)

    Arguments:  (See NEQNF)

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_n2qnj(void (*fcn)(Mint, Mfloat*, Mfloat*), void (*lsjac)(Mint,
                    Mfloat*, Mfloat*), Mfloat *errrel, Mint *n, Mint *itmax,
                    Mfloat xinit[], Mfloat x[], Mfloat *fnorm, Mfloat fvec[],
                    Mfloat *fjac, Mfloat r[], Mfloat qtf[], Mfloat wk[])
#else
static void l_n2qnj(void (*fcn)(Mint, Mfloat[], Mfloat[]), void (*lsjac)(Mint,
                    Mfloat[], Mfloat[]), Mfloat *errrel, Mint *n, Mint *itmax,
                    Mfloat xinit[], Mfloat x[], Mfloat *fnorm, Mfloat fvec[],
                    Mfloat *fjac, Mfloat r[], Mfloat qtf[], Mfloat wk[])
#endif
#else
static void l_n2qnj(fcn, lsjac, errrel, n, itmax, xinit, x, fnorm, fvec, fjac,
                    r, qtf, wk)
	void            (*fcn) (), (*lsjac) ();
	Mfloat          *errrel;
	Mint            *n, *itmax;
	Mfloat           xinit[], x[], *fnorm, fvec[], *fjac, r[], qtf[],
	                wk[];
#endif
{
#define FJAC(I_,J_)	(fjac+(I_)*(*n)+(J_))
	Mint             info, iter, lr, maxfev, ml, mode, mu, nfev, nprint;
	static Mfloat    factor = 1.0e2;



	imsl_e1psh("N2QNJ ");

	if (imsl_n1rty(0) != 0)
		goto L_9000;

	info = 0;
	maxfev = *itmax;
	ml = *n - 1;
	mu = *n - 1;
	mode = 2;
	sset(*n, F_ONE, wk, 1);
	nprint = 0;
	lr = (*n * (*n + 1)) / 2;
	/* Copy initial guesses into X */
	scopy(*n, xinit, 1, x, 1);

	l_n3qnj(fcn, lsjac, errrel, n, x, fvec, fjac, r, qtf, &maxfev,
	       &ml, &mu, &mode, &factor, &nprint, &info, &nfev, &lr, &wk[0],
		   &wk[*n], &wk[*n * 2], &wk[*n * 3], &wk[*n * 4]);

	if (info == 5)
		info = 4;
        *fnorm = imsl_sdot(*n, fvec, 1, fvec, 1);

	if (info == 2) {
		iter = *itmax;
		imsl_e1sti(1, iter);
/*		(4, 1, "The number of calls to the function has    */
/*                      exceeded MAX_ITN = %(i1).  The user may try  */
/*                      a new initial guess.");                    */
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_FCN_EVALS);
	} else if (info == 3) {
		imsl_e1str(1, *errrel);
/*		(4, 2, "The bound for the relative error, ERRREL =  */
/*                      %(r1), is too small.  No further improvement*/
/*                      in the approximate solution is possible.    */
/*                      The user should increase ERRREL.");         */
                imsl_ermes(IMSL_WARNING, IMSL_NO_BETTER_POINT);
	} else if (info == 4) {
/*		(4, 3, "The iteration has not made good progress.  */
/*                      The user may try a new initial guess.");   */
                imsl_ermes(IMSL_WARNING, IMSL_NO_PROGRESS);
	}
L_9000:
	imsl_e1pop("N2QNJ ");

	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N3QNJ/DN3QNJ (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 30, 1985

    Purpose:    Solve a system of nonlinear equations using the
                Levenberg-Marquardt algorithm with a user-supplied
                Jacobian.

    Usage:      CALL N3QNJ (FCN, LSJAC, ERRREL, N, X, FVEC, FJAC, R,
                            QTF, MAXFEV, ML, MU, MODE, FACTOR, NPRINT,
                            INFO, NFEV, LR, DIAG, WK1, WK2, WK3, WK4)

    Arguments:
       FCN    - A real function subroutine supplied by the
                user.  FCN must be declared EXTERNAL in the
                CALLING PROGRAM.  FCN specifies the system of
                equations to be solved and should be of the
                following form
                              SUBROUTINE FCN(X,F,N)
                              REAL X(*),F(*)
                              F(1)=
                               .
                              F(N)=
                              RETURN
                              END
                Where X is given.  FCN must not alter X.
       LSJAC  - A real function subroutine supplied by the user.
                LSJAC must be declared txternal in the CALLING
                PROGRAM.  LSJAC specifies the Jacobian matrix
                of the system of equations to be solved and should
                be of the following form
                              SUBROUTINE LSJAC (N, X, FJAC)
                              INTEGER N
                              REAL X, FJAC(N,N)
                              FJAC(1,1) =
                                .
                                .
                                .
                              FJAC(N,N) =
                              RETURN
                              END
       ERRREL - Stopping criterion.  The root is accepted if the
                relative error between two successive approximations
                to this root is within ERRREL.  (Input)
       N      - The number of equations to be solved and the number
                of unknowns.  (Input)
       X      - A vector of length N.  X contains the best estimate
                of the root found by NEQNF.  (Output)
       FVEC   - A vector of length N.  FVEC contains the functions
                evaluated at the point X.
       FJAC   - An N by N matrix.  FJAC contains the orthogonal
                matrix Q produced by the QR factorization of the
                final approximate Jacobian.
       R      - A vector of length N*(N+1)/2.  R contains the upper
                triangular matrix produced by the QR factorization
                of the final approximation Jacobian.  R is stored
                rowwise.
       QTF    - A vector of length N.  QTF contains the vector
                (Q transpose)*FVEC.
       MAXFEV - Maximum number of calls to FCN. (Input)
       ML     -
       MU     -
       MODE   -
       FACTOR -
       NPRINT - Number of iterates to be printed. (
       INFO   -
       NFEV   - Number of calls to FCN.  (Output)
       LR     - Length of the vector R.  (Input)
       DIAG   -
       WK1    - Real work vector of length N.  (Output)
       WK2    - Real work vector of length N.  (Output)
       WK3    - Real work vector of length N.  (Output)
       WK4    - Real work vector of length N.  (Output)

    Copyright:  1985 by IMSL, Inc.  All rights reserved

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* ML and MU are not used here, but leave in calling sequence. */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_n3qnj(void (*fcn)(Mint, Mfloat*, Mfloat*), void (*lsjac)(
                    Mint, Mfloat*, Mfloat*), Mfloat *errrel, Mint *n,
                    Mfloat x[], Mfloat fvec[], Mfloat *fjac, Mfloat r[],
                    Mfloat qtf[], Mint *maxfev, Mint *ml, Mint *mu, Mint *mode, 
                    Mfloat *factor, Mint *nprint, Mint *info, Mint *nfev,
                    Mint *lr, Mfloat diag[], Mfloat wk1[], Mfloat wk2[],
                    Mfloat wk3[], Mfloat wk4[])
#else
static void l_n3qnj(void (*fcn)(Mint, Mfloat[], Mfloat[]), void (*lsjac)(
                    Mint, Mfloat[], Mfloat[]), Mfloat *errrel, Mint *n,
                    Mfloat x[], Mfloat fvec[], Mfloat *fjac, Mfloat r[],
                    Mfloat qtf[], Mint *maxfev, Mint *ml, Mint *mu, Mint *mode, 
                    Mfloat *factor, Mint *nprint, Mint *info, Mint *nfev,
                    Mint *lr, Mfloat diag[], Mfloat wk1[], Mfloat wk2[],
                    Mfloat wk3[], Mfloat wk4[])
#endif
#else
static void l_n3qnj(fcn, lsjac, errrel, n, x, fvec, fjac, r, qtf, maxfev,
                    ml, mu, mode, factor, nprint, info, nfev, lr, diag, wk1,
	            wk2, wk3, wk4)
	void            (*fcn) (), (*lsjac) ();
	Mfloat          *errrel;
	Mint            *n;
	Mfloat           x[], fvec[], *fjac, r[], qtf[];
	Mint            *maxfev, *ml, *mu, *mode;
	Mfloat          *factor;
	Mint            *nprint, *info, *nfev, *lr;
	Mfloat           diag[], wk1[], wk2[], wk3[], wk4[];
#endif
{
#define FJAC(I_,J_)	(fjac+(I_)*(*n)+(J_))
	Mlong            jeval, sing;
	Mlong            _l0;
	Mint             i, iflag, iter, iwa[1], j, jm1, l, msum,
	                ncfail, ncsuc, nslow1, nslow2, ntmp1;
	Mfloat           actred, delta, epsmch, fnorm, fnorm1, pnorm, prered,
	                ratio, sum, temp, xnorm;
        Mlong            IMSLFALSE = 0, IMSLTRUE = 1;


	epsmch = imsl_amach(4);
	*info = 0;
	iflag = 0;
	*nfev = 0;
	/*
	 * Check the input parameters for errors
	 */
	if (*mode == 2) {
		for (j = 1; j <= *n; j++) {
			if (diag[j - 1] <= F_ZERO)
				goto L_150;
		}
	}
	/*
	 * Evaluate the function at the starting point and calculate its norm
	 */
	iflag = 1;
	imsl_e1usr("ON");
	(*fcn) (*n, x, fvec);
	imsl_e1usr("OFF");
	*nfev = 1;
	if (iflag < 0)
		goto L_150;
	fnorm = imsl_snrm2(*n, fvec, 1);
	/*
	 * Determine the number of calls to FCN Needed to compute the
	 * jacobian matrix
	 * 
	 * NOTE!!!  In NEQNJ, unlike NEQNF, the Jacobian is evaluated via the
	 * subroutine LSJAC, rather than by repeated calls to the function
	 * subroutine FCN!  Accordingly, MSUM has been set to 1 rather than
	 * MIN0(ML+MU+1,N). MSUM = MIN0(ML+MU+1,N)
	 */
	msum = 1;

	/*
	 * Initialize iteration counter and monitors
	 */
	iter = 1;
	ncsuc = 0;
	ncfail = 0;
	nslow1 = 0;
	nslow2 = 0;
	/* Beginning of the outer loop */
L_20:
	jeval = IMSLTRUE;
	/* Calculate the jacobian matrix */
	iflag = 2;
	(*lsjac) (*n, x, fjac);
	*nfev += msum;
	if (iflag < 0)
		goto L_150;
	/*
	 * Compute the QR factorization of the jacobian
	 */
        _l0 = IMSLFALSE;
        ntmp1 = 1;
	l_n5qnf(n, n, fjac, (Mlong *)&_l0, iwa, &ntmp1, wk1, wk2, wk3);
	/*
	 * On the first iteration and if MODE is 1, scale according to the
	 * norms of the columns of the Mintial Jacobian
	 */
	if (iter == 1) {
		if (*mode != 2) {
			scopy(*n, wk2, 1, diag, 1);
			for (j = 1; j <= *n; j++) {
				if (wk2[j - 1] == F_ZERO)
					diag[j - 1] = F_ONE;
			}
		}
		/*
		 * On the first iteration, calculate the norm of the scaled X
		 * and initialize the step bound delta
		 */
		for (j = 1; j <= *n; j++) {
			wk3[j - 1] = diag[j - 1] * x[j - 1];
		}
		xnorm = imsl_snrm2(*n, wk3, 1);
		delta = *factor * xnorm;
		if (delta == F_ZERO)
			delta = *factor;
	}
	/*
	 * Form (Q transpose)*FVEC and store in QTF.
	 */
	scopy(*n, fvec, 1, qtf, 1);
	for (j = 1; j <= *n; j++) {
		if (*FJAC(j - 1, j - 1) != F_ZERO) {
			sum = imsl_sdot(*n - j + 1, FJAC(j - 1, j - 1), 1, &qtf[j - 1],
					1);
			temp = -sum / *FJAC(j - 1, j - 1);
			saxpy(*n - j + 1, temp, FJAC(j - 1, j - 1), 1, &qtf[j - 1],
				   1);
		}
	}
	/*
	 * Copy the triangular factor of the QR factorization into R
	 */
	sing = IMSLFALSE;
	for (j = 1; j <= *n; j++) {
		l = j;
		jm1 = j - 1;
		if (jm1 >= 1) {
			for (i = 1; i <= jm1; i++) {
				r[l - 1] = *FJAC(j - 1, i - 1);
				l += *n - i;
			}
		}
		r[l - 1] = wk1[j - 1];
		if (wk1[j - 1] == F_ZERO)
			sing = IMSLTRUE;
	}
	/*
	 * Accumulate the orthogonal factor in FJAC
	 */
	l_n6qnf(n, n, fjac, wk1);
	/* Rescale if necessary */
	if (*mode != 2) {
		for (j = 1; j <= *n; j++) {
			diag[j - 1] = imsl_f_max(diag[j - 1], wk2[j - 1]);
		}
	}
	/*
	 * Beginning of the inner loop If requested, call FCN to enable
	 * prMinting of iterates
	 */
L_90:
	if (*nprint <= 0)
		goto L_100;
	if (iflag < 0)
		goto L_150;
L_100:
	;
	/* Determine the direction P */
	l_n7qnf(n, r, lr, diag, qtf, &delta, wk1, wk2, wk3);
	/*
	 * Store the direction P and X + P Calculate the norm of P
	 */
	sscal(*n, -F_ONE, wk1, 1);
	for (j = 1; j <= *n; j++) {
		wk2[j - 1] = x[j - 1] + wk1[j - 1];
		wk3[j - 1] = diag[j - 1] * wk1[j - 1];
	}
	pnorm = imsl_snrm2(*n, wk3, 1);
	/*
	 * On the first iteration, adjust the initial step bound
	 */
	if (iter == 1)
		delta = imsl_f_min(delta, pnorm);
	/*
	 * Evaluate the function at X + P and calculate its norm
	 */
	imsl_e1usr("ON");
	(*fcn) (*n, wk2, wk4);
	imsl_e1usr("OFF");
	*nfev += 1;
	fnorm1 = imsl_snrm2(*n, wk4, 1);
	/* Compute the scaled actual reduction */
	actred = -F_ONE;
	if (fnorm1 < fnorm)
		actred = F_ONE - imsl_fi_power(fnorm1 / fnorm, 2);
	/*
	 * Compute the scaled predicted reduction
	 */
	l = 1;
	for (i = 1; i <= *n; i++) {
		sum = imsl_sdot(*n - i + 1, &r[l - 1], 1, &wk1[i - 1], 1);
		l += *n - i + 1;
		wk3[i - 1] = qtf[i - 1] + sum;
	}
	temp = imsl_snrm2(*n, wk3, 1);
	prered = F_ONE;
	if (temp < fnorm)
		prered = F_ONE - imsl_fi_power(temp / fnorm, 2);
	/*
	 * Compute the ratio of the actual to the prdeicted reduction
	 */
	ratio = F_ZERO;
	if (prered > F_ZERO)
		ratio = actred / prered;
	/* Update the step bound */
	if (ratio >= 0.1) {
		ncfail = 0;
		ncsuc += 1;
		if (ratio >= F_HALF || ncsuc > 1)
			delta = imsl_f_max(delta, pnorm / F_HALF);
		if (fabs(ratio - F_ONE) <= 0.1)
			delta = pnorm / F_HALF;
	} else {
		ncsuc = 0;
		ncfail += 1;
		delta *= F_HALF;
	}
	/* Test for successful iteration */
	if (ratio >= 0.0001) {
		/*
		 * Successful iteration. Update X, FVEC, and their norms
		 */
		scopy(*n, wk2, 1, x, 1);
		scopy(*n, wk4, 1, fvec, 1);
		for (j = 1; j <= *n; j++) {
			wk2[j - 1] = diag[j - 1] * x[j - 1];
		}
		xnorm = imsl_snrm2(*n, wk2, 1);
		fnorm = fnorm1;
		iter += 1;
	}
	/*
	 * Determine the progress of the iteration
	 */
	nslow1 += 1;
	if (actred >= 0.001)
		nslow1 = 0;
	if (jeval)
		nslow2 += 1;
	if (actred >= 0.1)
		nslow2 = 0;
	/* Test for convergence */
	if (delta <= *errrel * xnorm || fnorm == F_ZERO)
		*info = 1;
	if (*info != 0)
		goto L_150;
	/*
	 * Tests for termination and stringent tolerances
	 */
	if (*nfev >= *maxfev)
		*info = 2;
	if (0.1 * imsl_f_max(0.1 * delta, pnorm) <= epsmch * xnorm)
		*info = 3;
	if (nslow2 == 5)
		*info = 4;
	if (nslow1 == 10)
		*info = 5;
	if (*info != 0)
		goto L_150;
	/*
	 * Criterion for recalculating Jacobian approximation by forward
	 * differences
	 */
	if (ncfail != 2) {
		/*
		 * Calculate the rank one modification to the Jacobian and
		 * update QTF if necessary
		 */
		for (j = 1; j <= *n; j++) {
			sum = imsl_sdot(*n, FJAC(j - 1, 0), 1, wk4, 1);
			wk2[j - 1] = (sum - wk3[j - 1]) / pnorm;
			wk1[j - 1] = diag[j - 1] * ((diag[j - 1] * wk1[j - 1]) / pnorm);
			if (ratio >= 0.0001)
				qtf[j - 1] = sum;
		}
		/*
		 * Compute the QR factoracation of the updated Jacobian
		 */
                ntmp1 = 1;
		l_n8qnf(n, n, r, wk1, wk2, wk3, &sing);
		l_n9qnf(n, n, fjac, n, wk2, wk3);
		l_n9qnf(&ntmp1, n, qtf, &ntmp1, wk2, wk3);
		/* End of the inner loop */
		jeval = IMSLFALSE;
		goto L_90;
	}
	/* End of the outer loop */
	goto L_20;
	/*
	 * Termination, either normal or user imposed
	 */
L_150:
	if (iflag < 0)
		*info = iflag;
	iflag = 0;
	return;
}				/* end of function */
