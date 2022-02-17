#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif 

#ifdef ANSI
static VA_LIST_HACK	l_zeros_fcn(Mfloat (*fcn)(Mfloat), va_list argptr);
static void	l_zreal(Mfloat (*f)(Mfloat), Mfloat errabs, Mfloat errrel,
			Mfloat eps, Mfloat eta, Mint nroot, Mint itmax,
			Mfloat xguess[], Mfloat x[], Mint infer[]);
#else
static VA_LIST_HACK	l_zeros_fcn();
static void	l_zreal();
#endif

static Mfloat	*lv_x;


#ifdef ANSI
Mfloat  *imsl_f_zeros_fcn(Mfloat (*fcn)(Mfloat), ...)
#else
Mfloat  *imsl_f_zeros_fcn(fcn, va_alist)
    Mfloat (*fcn) ();
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr, fcn);
    E1PSH("imsl_f_zeros_fcn", "imsl_d_zeros_fcn");
    lv_x = NULL;
    IMSL_CALL(l_zeros_fcn(fcn, argptr));
    va_end(argptr);
    E1POP("imsl_f_zeros_fcn", "imsl_d_zeros_fcn");
    return lv_x;
}

#ifdef ANSI
static VA_LIST_HACK l_zeros_fcn(Mfloat (*fcn)(Mfloat), va_list argptr)
#else
static VA_LIST_HACK l_zeros_fcn(fcn, argptr)
    Mfloat (*fcn) ();
    va_list	argptr;
#endif
{
    Mint	    code;
    Mint	    arg_number  = 1;
    Mfloat	    err_abs;
    Mfloat	    err_rel;
    Mfloat	    eps;
    Mfloat	    eta         = 0.01;
    Mfloat	    *xguess     = NULL;
    Mint 	    num_roots   = 1;
    Mint	    itmax       = 100;
    Mint	    *info       = NULL;
    Mint	    **pinfo     = 0;
    Mint	    user_guess  = 0;
    Mint	    user_info   = 0;
    Mint            user_pinfo  = 0;
    Mint            return_user = 0;
    Mint	    i;

    eps     = sqrt(imsl_amach(4));
    err_abs = eps; 
    err_rel = eps;

    code = 1;
    while (code > 0) {
	code = va_arg(argptr, Mint);
        arg_number++;
	switch (code) { 
	    case IMSL_RETURN_USER:
		arg_number++;
		lv_x = va_arg(argptr, Mfloat*);
                return_user = 1;
		break;
	    case IMSL_ERR_ABS:
		arg_number++;
		err_abs = (Mfloat) va_arg(argptr, Mdouble);
		break;
	    case IMSL_ERR_REL:
		arg_number++;
		err_rel = (Mfloat) va_arg(argptr, Mdouble);
		break;
	    case IMSL_ETA:
		arg_number++;
		eta = (Mfloat) va_arg(argptr, Mdouble);
		break;
	    case IMSL_EPS:
		arg_number++;
		eps = (Mfloat) va_arg(argptr, Mdouble);
		break;
	    case IMSL_ERR_ABS_ADR:
		arg_number++;
		err_abs = *(va_arg(argptr, Mfloat *));
		break;
	    case IMSL_ERR_REL_ADR:
		arg_number++;
		err_rel = *(va_arg(argptr, Mfloat *));
		break;
	    case IMSL_ETA_ADR:
		arg_number++;
		eta = *(va_arg(argptr, Mfloat *));
		break;
	    case IMSL_EPS_ADR:
		arg_number++;
		eps = *(va_arg(argptr, Mfloat *));
		break;
	    case IMSL_MAX_ITN:
		arg_number++;
		itmax = va_arg(argptr, Mint);
		break;
	    case IMSL_INFO_USER:
		arg_number++;
		user_info = 1;
		info = va_arg(argptr, Mint*);
		break;
	    case IMSL_INFO:
		arg_number++;
		user_pinfo = 1;
		pinfo = va_arg(argptr, Mint**);
		break;
	    case IMSL_NUM_ROOTS:
		arg_number++;
		num_roots = va_arg(argptr, Mint);
		break;
	    case IMSL_XGUESS:
		arg_number++;
		user_guess = 1;
		xguess = va_arg(argptr, Mfloat*);
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

    if (info == NULL) 
        info = (Mint *) imsl_malloc(num_roots*sizeof(*info));

    if (!user_guess)
        xguess  = (Mfloat *) imsl_malloc(num_roots*sizeof(*xguess));

    if (xguess==NULL) {
        imsl_e1sti(1, num_roots);
        imsl_e1stl(1, "num_roots");
        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	goto FREE_SPACE;
    }
    else if (!user_guess) {
        for (i=0; i<num_roots; ++i)
             *(xguess+i) = F_ZERO;
    }

    if (lv_x == NULL) {
	lv_x = (Mfloat *) imsl_malloc (num_roots*sizeof(*lv_x));
	if (lv_x == NULL) {
	    imsl_e1sti(1, num_roots);
	    imsl_e1stl(1, "num_roots");
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE;
	}
    }

    l_zreal (fcn, err_abs, err_rel, eps, eta, num_roots, itmax, xguess,
             lv_x, info);

    if (user_pinfo) *pinfo = info;

FREE_SPACE:
    if (!user_info && !user_pinfo && info != NULL)  imsl_free(info);
    if (!user_guess && xguess != NULL)              imsl_free(xguess);

RETURN:
    if (imsl_n1rty(0) > 4) {
        if (!return_user && lv_x != NULL)  imsl_free(lv_x);
        lv_x = NULL;
    }
    return (argptr);
}



/* -----------------------------------------------------------------------
    IMSL Name:  ZREAL/DZREAL (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 16, 1985

    Purpose:    Find the real zeros of a real function using Muller's
                method.

    Usage:      CALL ZREAL (F, ERRABS, ERRREL, EPS, ETA, NROOT,
                            ITMAX, XGUESS, X, INFO)

    Arguments:
       F      - User-supplied FUNCTION to compute the value of the
                function of which the zeros will be found.
                The form is
                F(X), where
                X      - The point at which the function is evaluated.
                         (Input)
                         X should not be changed by F.
                F      - The computed function value at the point X.
                         (Output)
                F must be declared EXTERNAL in the calling program.
       ERRABS - First stopping criterion.  (Input)
                A zero X(I) is accepted if
                    ABS(F(X(I)) .LT. ERRABS.
       ERRREL - Second stopping criterion is the relative error.  (Input)
                A zero X(I) is accepted if the relative change
                of two successive approximations to X(I) is less
                than ERRREL.
       EPS    - See ETA.  (Input)
       ETA    - Spread criteria for multiple zeros.  (Input)
                If the zero X(I) has been computed and
                ABS(X(I)-X(J)) .LT. EPS, where X(J) is a previously
                computed zero, then the computation is restarted with
                a guess equal to X(I) + ETA.
       NROOT  - The number of zeros to be found by ZREAL.  (Input)
       ITMAX  - The maximum allowable number of iterations per
                zero.  (Input)
       XGUESS - A vector of length NROOT.  (Input)
                XGUESS contains the initial guesses for the zeros.
       X      - A vector of length NROOT.  (Output)
                X contains the computed zeros.
       INFO   - An integer vector of length NROOT.  (Output)
                INFO(J) contains the number of iterations used in
                finding the J-th zero when convergence was achieved.
                If convergence was not obtained in ITMAX iterations,
                INFO(J) will be greater than ITMAX.

    Remarks:
    1. Informational error
       Type Code
         3   1  Failure to converge within ITMAX iterations for at
                least one of the NROOT roots.

    2. ZREAL always returns the last approximation for zero J in
       X(J).  If the convergence criterion is satisfied, then INFO(J)
       is less than or equal to ITMAX.  If the convergence criterion
       is not satisfied, then INFO(J) is set to ITMAX+1.

    3. ZREAL assumes that there exist NROOT distinct real zeros for
       the function F and that they can be reached from the initial
       guesses supplied.  The routine is designed so that convergence
       to any single zero cannot be obtained from two different
       initial guesses.

    4. Scaling the X vector in the function F may be required, if
       any of the zeros are known to be less than one.

    Keywords:   Muller; Roots; Zeros; Nonlinear equation; Twice
                continuously differentiable

    GAMS:       F1a

    Chapter:    MATH/LIBRARY Nonlinear Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_zreal(Mfloat (*f)(Mfloat), Mfloat errabs, Mfloat errrel, 
                    Mfloat eps, Mfloat eta, Mint nroot, Mint itmax,
		    Mfloat xguess[], Mfloat x[], Mint infer[])
#else
static void l_zreal(f, errabs, errrel, eps, eta, nroot, itmax, xguess, x, infer)
	Mfloat           (*f) (), errabs, errrel, eps, eta;
	Mint             nroot, itmax;
	Mfloat           xguess[], x[];
	int             infer[];
#endif
{
	Mint             i, jk, l;
	Mfloat           bi, d, dd, den, di, digt, dm, dn, fprt, frt, h, p;
	Mfloat		 p1, p2, rt, tem, x0, x1, x2;


/*	imsl_e1psh("ZREAL");    */
	/* COPY GUESSES INTO VECTOR X */
	scopy(nroot, xguess, 1, x, 1);

	digt = errrel;
	p = -F_ONE;
	p1 = F_ONE;
	p2 = F_ZERO;
	h = F_ZERO;
	for (l = 1; l <= nroot; l++) {
		jk = 0;
		if (x[l - 1] != F_ZERO) {
			p = 0.9 * x[l - 1];
			p1 = 1.1 * x[l - 1];
			p2 = x[l - 1];
		}
		rt = p;
		goto L_60;
L_10:
		if (jk != 1)
			goto L_20;
		rt = p1;
		x0 = fprt;
		goto L_60;
L_20:
		if (jk != 2)
			goto L_30;
		rt = p2;
		x1 = fprt;
		goto L_60;
L_30:
		if (jk != 3)
			goto L_50;
		x2 = fprt;
		d = -F_HALF;
		if (x[l - 1] == F_ZERO) {
			h = -F_ONE;
		} else {
			h = -0.1 * x[l - 1];
		}
L_40:
		dd = F_ONE + d;
		bi = x0 * imsl_fi_power(d, 2) - x1 * imsl_fi_power(dd, 2) + x2 * (dd + d);
		den = imsl_fi_power(bi, 2) - F_FOUR * x2 * d * dd * (x0 * d - (x1 * dd) + x2);
		if (den <= F_ZERO) {
			den = F_ZERO;
		} else {
			den = sqrt(den);
		}
		dn = bi + den;
		dm = bi - den;
		if (fabs(dn) <= fabs(dm)) {
			den = dm;
		} else {
			den = dn;
		}
		if (den == F_ZERO)
			den = F_ONE;
		di = -dd * (x2 + x2) / den;
		h *= di;
		rt += h;
		/* TEST FOR CONVERGENCE */
		if (fabs(h) < fabs(rt) * digt)
			goto L_80;
		goto L_60;
L_50:
		if (fabs(fprt) < fabs(x2 * F_TEN)) {
			x0 = x1;
			x1 = x2;
			x2 = fprt;
			d = di;
			goto L_40;
		}
		di *= F_HALF;
		h *= F_HALF;
		rt -= h;
L_60:
		jk += 1;
		if (jk >= itmax) {
			/* WARNING ERROR: ITERATIONS = MAXIMUM.  Failure to */
                        /* converge within  ITMAX = %(I1) iterations for at */
                        /* least one of the NROOT = %(i2) roots. */
			x[l - 1] = rt;
			infer[l - 1] = itmax + 1;
			imsl_e1sti(1, itmax);
			imsl_e1sti(2, nroot);
			imsl_ermes(IMSL_WARNING, IMSL_NO_CONVERGE_MAX_ITER);
			goto L_90;
		}
		imsl_e1usr("ON");
		frt = (*f) (rt);
		imsl_e1usr("OFF");
		fprt = frt;
		if (l >= 2) {
			for (i = 2; i <= l; i++) {
				tem = rt - x[i - 2];
				if (fabs(tem) < eps) {
					rt += eta;
					jk -= 1;
					goto L_60;
				}
				fprt /= tem;
			}
		}
		/* TEST FOR CONVERGENCE */
		if ((fabs(frt) >= errabs) || (fabs(fprt) >= errabs))
			goto L_10;

L_80:
		x[l - 1] = rt;
		infer[l - 1] = jk;

L_90:
		;
	}

/*	imsl_e1pop("ZREAL");    */

	return;
}				/* end of function */
