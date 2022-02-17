#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK   l_min_uncon_deriv(Mfloat (*fcn)(Mfloat), Mfloat (*grad)(Mfloat),
                                  Mfloat a, Mfloat b, va_list argptr);
static void      l_uvmid(Mfloat (*f)(Mfloat), Mfloat (*g)(Mfloat),
                         Mfloat *xguess, Mfloat *errrel, Mfloat *gtol, 
                         Mint *maxfn, Mfloat *a, Mfloat *b, Mfloat *x,
                         Mfloat *fx, Mfloat *gx);
static void      l_u3mid(Mfloat *x1, Mfloat *x2, Mfloat *f1, Mfloat *f2,
                         Mfloat *g1, Mfloat *g2, Mfloat *xmin);

#else
static VA_LIST_HACK   l_min_uncon_deriv();
static void      l_uvmid();
static void      l_u3mid();
#endif

static Mfloat    lv_minpt; 

#ifdef ANSI
Mfloat imsl_f_min_uncon_deriv(Mfloat (*fcn) (Mfloat), Mfloat (*grad) (Mfloat),
                             Mfloat a, Mfloat b, ...)
#else
Mfloat imsl_f_min_uncon_deriv(fcn, grad, a, b, va_alist)
    Mfloat (*fcn) ();
    Mfloat (*grad) ();
    Mfloat a;
    Mfloat b;
    va_dcl
#endif
{
    va_list argptr;

    VA_START(argptr, b);
    E1PSH("imsl_f_min_uncon_deriv", "imsl_d_min_uncon_deriv"); 
    lv_minpt = F_ZERO;
    IMSL_CALL(l_min_uncon_deriv(fcn, grad, a, b, argptr));
    va_end(argptr);
    E1POP("imsl_f_min_uncon_deriv", "imsl_d_min_uncon_deriv");  
    return lv_minpt;
}




#ifdef ANSI
static VA_LIST_HACK l_min_uncon_deriv(Mfloat (*fcn) (Mfloat), Mfloat (*grad)(Mfloat),
                                Mfloat a, Mfloat b, va_list argptr)
#else
static VA_LIST_HACK l_min_uncon_deriv(fcn, grad, a, b, argptr)
    Mfloat (*fcn) ();
    Mfloat (*grad) ();
    Mfloat a;
    Mfloat b;
    va_list       argptr;
#endif
{
    Mfloat	  a_float;
    Mfloat	  b_float;
    Mint          code;
    Mint          arg_number  = 4;
    Mint          maxfn       = 1000;
    Mfloat        xguess;
    Mfloat        err_rel     = 0.0001;
    Mfloat        grad_tol;
    Mfloat        sqrteps;
    Mfloat        f           = F_ZERO;
    Mfloat        g           = F_ZERO;
    Mfloat        *fvalue     = NULL;
    Mfloat        *gvalue     = NULL;
    Mfloat        *x          = NULL;
    Mint          fvalue_user = 0;
    Mint          gvalue_user = 0;
    Mint          x_user      = 0;

    xguess   = (a+b)/F_TWO;
    sqrteps  = sqrt(imsl_amach(4));
    grad_tol = sqrteps;
    err_rel  = sqrteps;

    code = 1;
    while (code > 0) {
        code = va_arg(argptr, Mint);
        arg_number++;
        switch(code) {

            case IMSL_XGUESS:
                arg_number++;
                xguess = (Mfloat) va_arg(argptr, Mdouble);
                break;
            case IMSL_ERR_REL:
                arg_number++; 
                err_rel = (Mfloat) va_arg(argptr, Mdouble);
                break;
            case IMSL_GRAD_TOL:
                arg_number++; 
                grad_tol = (Mfloat) va_arg(argptr, Mdouble);
                break;
            case IMSL_XGUESS_ADR:
                arg_number++;
                xguess = *(va_arg(argptr, Mfloat *));
                break;
            case IMSL_ERR_REL_ADR:
                arg_number++; 
                err_rel = *(va_arg(argptr, Mfloat *));
                break;
            case IMSL_GRAD_TOL_ADR:
                arg_number++; 
                grad_tol = *(va_arg(argptr, Mfloat *));
                break;
            case IMSL_MAX_FCN:
                arg_number++; 
                maxfn = va_arg(argptr, Mint);
                break;
            case IMSL_FVALUE:
                arg_number++;
                fvalue = va_arg(argptr, Mfloat*);
                fvalue_user = 1;
                break;
            case IMSL_GVALUE:
                arg_number++;
                gvalue = va_arg(argptr, Mfloat*);
                gvalue_user = 1;
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

    a_float = a; 
    b_float = b;

    l_uvmid(fcn, grad, &xguess, &err_rel, &grad_tol, &maxfn, &a_float, 
            &b_float, &lv_minpt, &f, &g);

    if (x_user) *x = lv_minpt;
    if (fvalue_user) *fvalue = f;
    if (gvalue_user) *gvalue = g;

RETURN:
    if (imsl_n1rty(0) > 4) lv_minpt = imsl_amach(6);  /* NaN */
    return (argptr);
}
/* -----------------------------------------------------------------------
    IMSL Name:  UVMID/DUVMID (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 11, 1985

    Purpose:    Find the minimum point of a smooth function of a single
                variable using both function evaluations and first
                derivative evaluations.

    Usage:      CALL UVMID (F, G, XGUESS, ERRREL, GTOL, MAXFN, A, B, X,
                            FX, GX)

    Arguments:
       F      - User-supplied FUNCTION to define the function to be
                minimized.  The form is
                F(X), where
                X      - The point at which the function is to be
                         evaluated.  (Input)
                F      - The computed value of the function at X.
                         (Output)
                F must be declared EXTERNAL in the calling program.
       G      - User-supplied FUNCTION to compute the derivative of
                the function.  The form is
                G(X), where
                X      - The point at which the derivative is to be
                         computed.  (Input)
                G      - The computed value of the derivative at X.
                         (Output)
                G must be declared EXTERNAL in the calling program.
       XGUESS - An initial guess of the minimum point of F.  (Input)
       ERRREL - The required relative accuracy in the final value of X.
                (Input)
                This is first stopping criterion.  On a normal return
                the solution X is in an interval which contains a local
                minimum and is less than or equal to MAX(1.0,ABS(X)) *
                ERRREL.  When the given ERRREL is less than machine
                epsilon, SQRT(machine epsilon) is used as ERRREL.
       GTOL   - The derivative tolerance used to decide if the current
                point is a local minimum.  (Input)
                This is second stopping criterion.  X is returned as a
                solution when GX is less than or equal to GTOL.
                GTOL should be nonnegative, otherwise zero would be used.
       MAXFN  - Maximum number of function evaluations allowed.
                (Input)
       A      - A is the lower endpoint of the interval in which the
                minimum point of F is to be located.  (Input)
       B      - B is the upper endpoint of the interval in which the
                minimum point of F is to be located.  (Input)
       X      - The point at which a minimum value of F is found.
                (Output)
       FX     - The function value at point X.  (Output)
       GX     - The derivative value at point X.  (Output)

    Remark:
       Informational errors
       Type Code
         3   1  The final value of X is at the lower bound.  The minimum
                is probably beyond the bound.
         3   2  The final value of X is at the upper bound.  The minimum
                is probably beyond the bound.
         4   3  The maximum number of function evaluations has been
                exceeded.

    Keywords:   Unconstrained univariate minimization; Safeguarded
                cubic interpolation; Secant method; Optimization;
                Univariate function

    GAMS:       G1a1b

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_uvmid(Mfloat (*f)(Mfloat), Mfloat (*g)(Mfloat), Mfloat *xguess,
                    Mfloat *errrel, Mfloat *gtol, Mint *maxfn, Mfloat *a,
                    Mfloat *b, Mfloat *x, Mfloat *fx, Mfloat *gx)
#else
static void l_uvmid(f, g, xguess, errrel, gtol, maxfn, a, b, x, fx, gx)
	Mfloat           (*f) (), (*g) (), *xguess, *errrel, *gtol;
	Mint            *maxfn;
	Mfloat          *a, *b, *x, *fx, *gx;
#endif
{
	Mint             nf, ng;
	Mfloat           d, eps, fa, fb, fxl, fxmin, gxl, gxmin, st, tmp,
	                 tol, tolg, tot, xl, xmin, xtol;


	imsl_e1psh("UVMID ");

	*x = *xguess;
	if (*a > *x || *b < *x) {
		imsl_e1str(1, *a);
		imsl_e1str(2, *x);
		imsl_e1str(3, *b);

/*		imsl_ermes(5, 1, "The initial guess must lie between A and B, while A = %(r1), XGUESS = %(r2) and B = %(r3) are given.");
*/              imsl_ermes(IMSL_TERMINAL, IMSL_NEED_A_LE_XGUESS_LE_B);
	} else {
		eps = imsl_amach(4);
		if (*errrel < F_ZERO) {
			tol = sqrt(eps);
		} else {
			tol = *errrel;
		}
		tolg = imsl_f_max(*gtol, F_ZERO);
		nf = 0;
		ng = 0;
		imsl_e1usr("ON");
		fa = (*f) (*a);
		imsl_e1usr("OFF");

		imsl_e1usr("ON");
		fb = (*f) (*b);
		imsl_e1usr("OFF");

		imsl_e1usr("ON");
		*fx = (*f) (*x);
		imsl_e1usr("OFF");
		nf += 3;

		if (fa < *fx) {
			if (fa < fb) {
				imsl_e1usr("ON");
				*gx = (*g) (*a);
				imsl_e1usr("OFF");
				ng += 1;
				*x = *a;
				*fx = fa;
				if (*gx > F_ZERO) {

/*					imsl_ermes(3, 1, "The final value of x is at the lower bound.");
*/
                                        imsl_ermes(IMSL_WARNING, IMSL_MIN_AT_LOWERBOUND);
					goto L_9000;
				}
			} else {
				imsl_e1usr("ON");
				*gx = (*g) (*b);
				imsl_e1usr("OFF");
				ng += 1;
				*x = *b;
				*fx = fb;
				if (*gx < F_ZERO) {

/*					imsl_ermes(3, 2, "The final value of x is at the upper bound. ");
*/
                                        imsl_ermes(IMSL_WARNING, IMSL_MIN_AT_UPPERBOUND);
					goto L_9000;
				}
			}
		} else if (fb < *fx) {
			imsl_e1usr("ON");
			*gx = (*g) (*b);
			imsl_e1usr("OFF");
			ng += 1;
			*x = *b;
			*fx = fb;
			if (*gx < F_ZERO) {

/*				imsl_ermes(3, 2, "The final value of x is at the upper bound.");
*/
                                imsl_ermes(IMSL_WARNING, IMSL_MIN_AT_UPPERBOUND);
				goto L_9000;
			}
		} else {

			imsl_e1usr("ON");
			*gx = (*g) (*x);
			imsl_e1usr("OFF");
			ng += 1;
		}
		/* Search for initial interval. */
		st = fabs(*gx);
L_10:
		if (fabs(*gx) <= tolg)
			goto L_9000;
		xtol = imsl_f_max(F_ONE, fabs(*x)) * tol;
		if (st < xtol)
			st = xtol;
		if (*gx > F_ZERO) {
			xl = *x - st;
			if (xl <= *a) {
				xl = *a;
				fxl = fa;
				goto L_20;
			}
		} else {
			xl = *x + st;
			if (xl >= *b) {
				xl = *b;
				fxl = fb;
				goto L_20;
			}
		}

		imsl_e1usr("ON");
		fxl = (*f) (xl);
		imsl_e1usr("OFF");
		nf += 1;
		if (nf > *maxfn) {
/*			imsl_ermes(4, 3, "The maximum number of function evaluations has been exceeded.");
*/
                        imsl_e1sti(1, *maxfn);
                        imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_FCN_EVAL);
			goto L_9000;
		}
L_20:
		imsl_e1usr("ON");
		gxl = (*g) (xl);
		imsl_e1usr("OFF");
		ng += 1;

		if (fxl < *fx) {
			if (gxl ** gx > F_ZERO) {
				tmp = (*gx - gxl) / st;
				*x = xl;
				*fx = fxl;
				*gx = gxl;
				st = fabs(*gx / tmp);
				goto L_10;
			} else {
				tmp = *x;
				*x = xl;
				xl = tmp;
				tmp = *fx;
				*fx = fxl;
				fxl = tmp;
				tmp = *gx;
				*gx = gxl;
				gxl = tmp;
			}
		}
		/*
		 * Initial interval found, check the convergence.
		 */
		tot = fabs(xl - *x);
		if (fabs(tot - xtol) <= F_TEN * eps)
			goto L_9000;
		if (fabs(*gx) <= tolg)
			goto L_9000;
		/* Call cubic fit routine. */
L_30:
		l_u3mid(x, &xl, fx, &fxl, gx, &gxl, &xmin);
		d = xmin - *x;
		tot = xl - *x;
		if (fabs(d) < 0.25e0 * fabs(tot)) {
			d = 0.25e0 * tot;
		} else if (fabs(d) > 0.75e0 * fabs(tot)) {
			d = 0.75e0 * tot;
		}
		xtol = imsl_f_max(F_ONE, fabs(*x)) * tol;
		if (fabs(d) < xtol)
			d = sign(xtol, d);
		xmin = *x + d;
		imsl_e1usr("ON");
		fxmin = (*f) (xmin);
		imsl_e1usr("OFF");
		nf += 1;
		if (nf > *maxfn) {
/*			imsl_ermes(4, 3, "The maximum number of function evaluations has been exceeded.");
*/
                        imsl_e1sti(1, *maxfn);
                        imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_FCN_EVAL);
			goto L_9000;
		}
		imsl_e1usr("ON");
		gxmin = (*g) (xmin);
		imsl_e1usr("OFF");
		ng += 1;

		if (fxmin < *fx) {
			/* Decide next interval */
			if (*gx * gxmin < F_ZERO) {
				xl = *x;
				fxl = *fx;
				gxl = *gx;
				*x = xmin;
				*fx = fxmin;
				*gx = gxmin;
			} else {
				*x = xmin;
				*fx = fxmin;
				*gx = gxmin;
			}
		} else {
			xl = xmin;
			fxl = fxmin;
			gxl = gxmin;
		}
		/* Check for convergence. */
		tot = fabs(xl - *x);
		if (fabs(tot - xtol) <= F_TEN * eps)
			goto L_9000;
		if (fabs(*gx) <= tolg)
			goto L_9000;
		goto L_30;
	}

L_9000:
	imsl_e1pop("UVMID ");
	return;
}				/* end of function */



/* -----------------------------------------------------------------------
    IMSL Name:  U3MID/DU3MID (Single/Double precision version)

    Computer:   $COMPUTER/$PRECISION

    Revised:    February 25, 1986

    Purpose:    Perform a cubic fit given two points, their function
                values and their gradient values.

    Usage:      CALL U3MID (X1, X2, F1, F2, G1, G2, XMIN)

    Arguments:
       X1     - The point yielding the smaller function value of the
                two points.  (Input)
       X2     - The point yielding the larger function value of the
                two points.  (Input)
       F1     - The value of the function evaluated at X1.  (Input)
       F2     - The value of the function evaluated at X2.  (Input)
       G1     - The value of the gradient evaluated at X1.  (Input)
       G2     - The value of the gradient evaluated at X2.  (Input)
       XMIN   - The point at which the cubic fit has a minimum.
                (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u3mid(Mfloat *x1, Mfloat *x2, Mfloat *f1, Mfloat *f2, Mfloat *g1,
                    Mfloat *g2, Mfloat *xmin)
#else
static void l_u3mid(x1, x2, f1, f2, g1, g2, xmin)
	Mfloat          *x1, *x2, *f1, *f2, *g1, *g2, *xmin;
#endif
{
	Mfloat           dx, tmp1, tmp2;


	dx = *x1 - *x2;
	tmp1 = *g2 + *g1 + F_THREE * (*f2 - *f1) / dx;
	tmp2 = sqrt(tmp1 * tmp1 - *g2 ** g1);
	if (dx > F_ZERO)
		tmp2 = -tmp2;

	*xmin = *x1 - dx * ((*g1 - tmp2 - tmp1) / (*g1 - *g2 - 2 * tmp2));
	if (*xmin < *x1) {
		if (*xmin < *x2)
			*xmin = *x1 - dx * ((*g1 + tmp2 - tmp1) / (*g1 - *g2 + 2 * tmp2));
	} else if (*xmin > *x1) {
		if (*xmin > *x2)
			*xmin = *x1 - dx * ((*g1 + tmp2 - tmp1) / (*g1 - *g2 + 2 * tmp2));
	}
	return;
}				/* end of function */
