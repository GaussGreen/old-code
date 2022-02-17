#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK   l_min_uncon(Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                             va_list argptr);
static void      l_uvmif(Mfloat (*f) (Mfloat), Mfloat *a, Mfloat *b, 
                         Mfloat *xguess, Mfloat *step, Mfloat *xacc, 
                         Mint *maxfn, Mfloat *x);

#else
static VA_LIST_HACK   l_min_uncon();
static void      l_uvmif();
#endif

static Mfloat    lv_minpt;    

#ifdef ANSI
Mfloat imsl_f_min_uncon(Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b, ...)
#else
Mfloat imsl_f_min_uncon(fcn, a, b, va_alist)
    Mfloat (*fcn) ();
    Mfloat a;
    Mfloat b;
    va_dcl
#endif
{
    va_list argptr;

    VA_START(argptr, b);
    E1PSH("imsl_f_min_uncon", "imsl_f_min_uncon");
    lv_minpt = F_ZERO;
    IMSL_CALL(l_min_uncon(fcn, a, b, argptr));
    va_end(argptr);
    E1POP("imsl_f_min_uncon", "imsl_f_min_uncon"); 
    return lv_minpt;
}




#ifdef ANSI
static VA_LIST_HACK l_min_uncon(Mfloat (*fcn) (Mfloat), Mfloat a, Mfloat b,
                           va_list argptr)
#else
static VA_LIST_HACK l_min_uncon(fcn, a, b, argptr)
    Mfloat (*fcn) ();
    Mfloat a;
    Mfloat b;
    va_list       argptr;
#endif
{
    Mint          code;
    Mint          arg_number = 3;
    Mint          max_fcn = 1000;
    Mfloat        xguess;
    Mfloat        step = F_ONE;
    Mfloat        err_abs = 0.0001;
    Mfloat        a_float, b_float;

    xguess  = (a+b)/F_TWO;
#ifdef DOUBLE
    err_abs = sqrt(imsl_dmach(4));
#endif

    code = 1;
    while (code > 0) {
        code = va_arg(argptr, int);
        arg_number++;
        switch(code) {
            case IMSL_XGUESS:
                arg_number++;
                xguess = (Mfloat) va_arg(argptr, Mdouble);
                break;
            case IMSL_STEP:
                arg_number++; 
                step = (Mfloat) va_arg(argptr, Mdouble);
                break;
            case IMSL_ERR_ABS:
                arg_number++; 
                err_abs = (Mfloat) va_arg(argptr, Mdouble);
                break;
            case IMSL_XGUESS_ADR:
                arg_number++;
                xguess = *(va_arg(argptr, Mfloat *));
                break;
            case IMSL_STEP_ADR:
                arg_number++; 
                step = *(va_arg(argptr, Mfloat *));
                break;
            case IMSL_ERR_ABS_ADR:
                arg_number++; 
                err_abs = *(va_arg(argptr, Mfloat *));
                break;
            case IMSL_MAX_FCN:
                arg_number++; 
                max_fcn = va_arg(argptr, Mint);
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

    a_float = (Mfloat) a;
    b_float = (Mfloat) b;
    l_uvmif(fcn, &a_float, &b_float, &xguess, &step, &err_abs, &max_fcn, 
            &lv_minpt);

RETURN:
    if (imsl_n1rty(0) > 4)  lv_minpt = imsl_amach(6);    /* NaN */
    return (argptr);
}

/* -----------------------------------------------------------------------
    IMSL Name:  UVMIF/DUVMIF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    July 11, 1985

    Purpose:    Find the minimum of a smooth univariate function within
                an interval.

    Usage:      CALL UVMRF (F, A, B, XGUESS, STEP, XACC, MAXFN, X)

    Arguments:
       F      - Function value evaluated at X.  (Input)
       A, B   - Bounds to define the interval in which the minimum has
                to be found.  (Input)
                B has to be greater than A.
       XGUESS - An initial guess of the minimum.  (Input)
       STEP   - An order of magnitude estimate of the required
                change in X.  (Input)
       XACC   - The required absolute accuracy in the final
                value of X.  On a normal return there are
                points on either side of X within a distance
                XACC at which F is no less than F at X.  (Input)
       MAXFN  - A limit, which must be set to a positive
                integer, on the number of function evaluations.  (Input)
       X      - The scalar variable of the calculation.  It gives
                the least calculated value of F.  (Output)

    Remarks:
       Informational errors
       Type Code
         3   1  The final value of X is at a bound.  The minimum
                is probably beyond the bound.
         3   2  Computer rounding errors prevent further
                refinement of X.
         4   3  The number of function evaluations has exceeded MAXFN.

    Keywords:   Unconstrained univariate minimization; Safeguarded
                quadratic interpolation

    GAMS:       G1a1a

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_uvmif (Mfloat (*f)(Mfloat), Mfloat *a, Mfloat *b, Mfloat *xguess,
                     Mfloat *step, Mfloat *xacc, Mint *maxfn, Mfloat *x)
#else
static void l_uvmif(f, a, b, xguess, step, xacc, maxfn, x)
	Mfloat          (*f) (), *a, *b, *xguess, *step, *xacc;
	Mint            *maxfn;
	Mfloat          *x;
#endif
{
	Mint             info, is, nf;
	Mfloat           acc, bl, bound, bu, da, db, dc, fa, fb, fc, fd, fx, h,
                         rho, sigma, st, temp, tol, xa, xb, xc, xd, xx;


	imsl_e1psh("UVMIF ");

	*x = *xguess;
	nf = 0;
	goto L_70;
L_10:
	xb = *x;
	fb = fx;
	info = 0;
	/*
	 * ALTER ANY UNSUITABLE INITIAL PARAMETERS, INCLUDING INCREASING THE
	 * STEP IF NECESSARY IN ORDER THAT THE COMPUTER PRECISION GIVES A
	 * CHANGE IN X.
	 */
	if (*x < *a || *x > *b) {
		info = 6;
		goto L_270;
	}
	acc = imsl_f_max(F_ZERO, *xacc);
	bl = *x - *a;
	bu = *b - *x;
	bound = imsl_f_max(bl, bu);
	if (bound <= acc) {
		info = 2;
		goto L_270;
	}
	st = *step;
L_15:
	if (st < F_ZERO) {
		if (fabs(st) > bl)
			st = -bl;
	} else {
		if (fabs(st) > bu)
			st = bu;
	}
	if (st == F_ZERO) {
		st = 0.01 * bound;
		if (bl > bu)
			st = -st;
	}
	if (fabs(st) < acc)
		st = sign(acc, st);
	goto L_30;
L_20:
	st *= F_FIVE;
L_30:
	xa = xb + F_HALF * fabs(st);
	xc = xb + fabs(st);
	if (xa <= xb)
		goto L_20;
	if (xc <= xa)
		goto L_20;

	if (*x + st < *a || *x + st > *b) {
		if (info == 2) {
			goto L_270;
		} else {
			info = 2;
			st = -st;
			goto L_15;
		}
	}
	/*
	 * CALCULATE THE NEXT TRIAL VALUE OF X FROM XB AND ST.
	 */
L_40:
	is = 1;
L_50:
	*x = xb + st;
	h = xb + 1.5 * st;
	if (h <= *a) {
		*x = *a;
		is = 2;
	} else if (h >= *b) {
		*x = *b;
		is = 2;
	}
	/* CALCULATE THE NEXT VALUE OF F. */
L_60:
	if (nf >= *maxfn) {
		info = 4;
		goto L_270;
	}
L_70:
	nf += 1;
	imsl_e1usr("ON ");
	fx = (*f)(*x);
	imsl_e1usr("OFF");
	if (nf <= 1)
		goto L_10;
	/*
	 * REVERSE ST IF THE INITIAL STEP SEEMS TO BE UPHILL.
	 */
	if (nf >= 3)
		goto L_80;
	if (fx < fb)
		goto L_90;
	xc = *x;
	fc = fx;
	if (xb == *a || xb == *b)
		goto L_95;
	st = -st;
	goto L_40;
	/*
	 * ENSURE THAT FB IS THE LEAST CALCULATED VALUE OF F.
	 */
L_80:
	xd = xc;
	fd = fc;
	if (fx >= fb)
		goto L_120;
L_90:
	xc = xb;
	fc = fb;
	xb = *x;
	fb = fx;
	/*
	 * CALCULATE AN EXTRA FUNCTION VALUE IF X IS AT A BOUND.
	 */
	if (is >= 4)
		goto L_140;
	if (is <= 1)
		goto L_110;
	if (is == 3)
		goto L_100;
L_95:
	is = 3;
	h = imsl_f_max(0.9e0 * acc, 0.01e0 * fabs(xb - xc));
	*x = xb + sign(h, xc - xb);
	if (fabs(*x - xc) < fabs(*x - xb))
		*x = F_HALF * (xb + xc);
	temp = (xb - *x) / (xb - xc);
	if (temp <= F_ZERO) {
		info = 3;
		goto L_270;
	}
	goto L_60;
	/*
	 * THIS STAGE IS REACHED WHEN A BRACKET IS FOUND NEAR A BOUND.
	 */
L_100:
	xa = xd;
	fa = fd;
	is = 4;
	goto L_170;
	/*
	 * CALCULATE THE NEXT STEP IN THE SEARCH FOR A BRACKET, TRYING TO
	 * OVERSHOOT THE PREDICTED MINIMUM BY THE FACTOR THREE.
	 */
L_110:
	if (nf <= 2)
		goto L_50;
	rho = (xb - xc) / (xb - xd);
	sigma = (fb - fc) / (fb - fd);
	h = F_NINE;
	if (sigma < rho)
		h = 1.5e0 * (rho - sigma / rho) / (sigma - rho);
	h = imsl_f_min(h, F_NINE);
	h = imsl_f_max(h, F_TWO);
	st *= h;
	goto L_50;
	/*
	 * RETURN IF THE MINIMUM SEEMS TO BE BEYOND A BOUND.
	 */
L_120:
	if (is >= 4)
		goto L_130;
	if (is == 3) {
		info = 3;
		goto L_270;
	}
	/*
	 * A BRACKET HAS BEEN FOUND SO BRANCH TO PREDICT THE MINIMUM.
	 */
	xa = *x;
	fa = fx;
	is = 4;
	goto L_170;
	/*
	 * ENSURE THAT XA, XB, XC AND XD ARE ORDERED MONOTONICALLY.
	 */
L_130:
	xc = *x;
	fc = fx;
L_140:
	temp = (xb - xc) / (xa - xd);
	if (temp > F_ZERO)
		goto L_160;
L_150:
	h = xa;
	xa = xd;
	xd = h;
	h = fa;
	fa = fd;
	fd = h;
	/*
	 * IF THERE ARE THREE CONSECUTIVE EQUAL VALUES OF F, ENSURE THAT FB
	 * IS THE MIDDLE ONE.
	 */
L_160:
	if (fa == fb)
		goto L_170;
	if (fd != fb)
		goto L_170;
	if (fc != fb)
		goto L_170;
	xc = xb;
	xb = *x;
	fc = fb;
	fb = fx;
	goto L_150;
	/*
	 * USE THE MINIMA OF TWO QUADRATICS TO CALCULATE A TOLERANCE ON THE
	 * NEXT
	 */
L_170:
	da = (fb - fa) / (xb - xa);
	db = (fc - fb) / (xc - xb);
	if (is >= 5)
		goto L_180;
	is = 5;
	tol = 0.01e0 * fabs(xa - xc);
	goto L_190;
L_180:
	dc = (fd - fc) / (xd - xc);
	tol = F_ZERO;
	if (db == F_ZERO)
		goto L_190;
	tol = fabs(xa - xc);
	temp = (dc - db) / (xa - xc);
	if (temp >= F_ZERO)
		goto L_190;
	h = F_HALF * fabs(db * ((xd - xb) / (dc - db) + (xa - xc) / (db - da)));
	if (h < tol)
		tol = h;
L_190:
	tol = imsl_f_max(tol, 0.9e0 * acc);
	/*
	 * SET X TO THE VALUE THAT MINIMIZES THE INTERPOLATING QUADRATIC.
	 */
	*x = xb;
	if (da != db)
		*x = F_HALF * (da * (xb + xc) - db * (xa + xb)) / (da - db);
	/* ENSURE THAT ABS(XA-XB).GE.ABS(XB-XC). */
	temp = (xa - xb) / (xb - xc);
	if (temp >= F_ONE)
		goto L_200;
	h = xa;
	xa = xc;
	xc = h;
	h = fa;
	fa = fc;
	fc = h;
	/* TEST FOR CONVERGENCE. */
L_200:
	if (fabs(xa - xb) <= acc)
		goto L_280;
	/*
	 * IF ABS(XA-XB).LE.2.9*ACC, CHOOSE THE NEXT X TO AVOID AN
	 * UNNECESSARY FUNCTION EVALUATION.
	 */
	temp = (xa - xb) / (xb - xc);
	if (temp > F_TEN)
		goto L_220;
	if (fabs(xa - xb) > 2.9 * acc)
		goto L_210;
	*x = F_HALF * (xa + xb);
	if (fabs(*x - xb) <= acc)
		goto L_230;
	*x = 0.67 * xb + 0.33 * xa;
	goto L_230;
	/*
	 * IF (XA-XB)/(XB-XC).LE.10, ENSURE THAT THE DISTANCE FROM X TO XB IS
	 * NORMALLY AT LEAST TOL.
	 */
L_210:
	if (fabs(*x - xb) >= tol)
		goto L_230;
	*x = xb + sign(tol, *x - xb);
	if (fabs(*x - xc) < fabs(*x - xb))
		*x = xb + sign(tol, xa - xb);
	if (fabs(*x - xa) < fabs(*x - xb))
		*x = F_HALF * (xa + xb);
	goto L_230;
	/*
	 * WHEN (XA-XB)/(XB-XC).GT.10, ENSURE THAT X IS IN THE LONGER
	 * INTERVAL, AND TRY TO OVERSHOOT THE MINIMUM.
	 */
L_220:
	h = (*x - xb) * sign(F_THREE, xa - xb);
	h = imsl_f_vmax(3, h, tol, fabs(xb - xc));
	h = imsl_f_min(h, 0.1 * fabs(xa - xb));
	*x = xb + sign(h, xa - xb);
	/*
	 * TEST WHETHER ROUNDING ERRORS MAKE THE NEW X THE SAME AS A PREVIOUS
	 * X.
	 */
L_230:
	if (*x == xb)
		goto L_240;
	if ((xa - *x) / (xa - xb) <= F_ZERO)
		goto L_240;
	if ((*x - xc) / (xa - xb) > F_ZERO)
		goto L_60;
	/*
	 * SEEK A NEW VALUE OF X BETWEEN XA AND XB, BUT ROUNDING ERRORS MAY
	 * NOT ALLOW ONE.
	 */
L_240:
	*x = xa;
L_250:
	xx = F_HALF * (*x + xb);
	temp = (xx - xb) / (xa - xb);
	if (temp <= F_ZERO)
		goto L_260;
	temp = (*x - xx) / (xa - xb);
	if (temp <= F_ZERO)
		goto L_260;
	*x = xx;
	goto L_250;
L_260:
	if (*x != xa)
		goto L_60;
	info = 5;
	/* RETURN FROM THE SUBROUTINE. */
L_270:
	if (info == 2) {
             /* imsl_ermes(5, 1, "The interval between bounds is too small.  No change in X can be made."); */
                imsl_ermes(IMSL_TERMINAL, IMSL_INTERVAL_TOO_SMALL);

	} else if (info == 6) {
             /* imsl_ermes(5, 2, "The initial guess is out of bounds.  Rerun with a new guess."); */
                imsl_ermes(IMSL_TERMINAL, IMSL_INIT_GUESS_OUT_BOUNDS);

	} else if (info == 5) {
             /* imsl_ermes(3, 1, "Computer rounding errors prevent further refinement of X.");*/
                imsl_ermes(IMSL_WARNING, IMSL_NO_MORE_PROGRESS);

	} else if (info == 3) {
             /* imsl_ermes(3, 2, "The final value for X is at a bound.\n"); */
                imsl_ermes(IMSL_WARNING, IMSL_MIN_AT_BOUND);

	} else if (info == 4) {
		imsl_e1sti(1, *maxfn);
             /* imsl_ermes(4, 3, "The maximum number of function evaluations has been exceeded."); */
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_FCN_EVAL);
	}
L_280:
	*x = xb;

	imsl_e1pop("UVMIF ");
	return;
}				/* end of function */
