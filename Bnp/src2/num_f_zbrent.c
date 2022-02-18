/*--------------------------------------------------------------
        FILE: num_f_zbrent.c
        PURPOSE: Brent method one-dimensional root finding from NR (Chapter 9.3)
        AUTHOR: Dimitri Mayevski
        DATE: 11/07/2002
  --------------------------------------------------------------*/

#include "math.h"
#include "num_h_zbrent.h"
#include "stddef.h"
#include "uterror.h"

#define FACTOR 1.6
#define NTRY 50

Err num_f_zbrac(
    Err (*func)(double, double*, void*),
    double* x1,
    double* x2,
    double* f1,
    double* f2,
    void*   static_data)
/*	Given a function func and an initial guessed range x1 to x2, the routine expands the range
        geometrically until a root is bracketed by the returned values x1 and x2
        or until the range becomes unacceptably large. The static_data parameter is passed
        as the last func's argument. */
{
    Err err = NULL;
    int j;
    if (*x1 == *x2)
        return serror("Bad initial range in num_f_zbrac");
    err = (*func)(*x1, f1, static_data);
    if (err)
        return err;
    err = (*func)(*x2, f2, static_data);
    if (err)
        return err;

    for (j = 1; j <= NTRY; j++)
    {
        if ((*f1) * (*f2) < 0.0)
            return NULL;
        if (fabs(*f1) < fabs(*f2))
            err = (*func)(*x1 += FACTOR * (*x1 - *x2), f1, static_data);
        else
            err = (*func)(*x2 += FACTOR * (*x2 - *x1), f2, static_data);
        if (err)
            return err;
    }
    return serror("Max number of iterations exceeded in num_f_zbrac");
}

#define ITMAX 50   /* Maximum allowed number of iterations */
#define EPS 3.0e-8 /* Machine floating-point precision */

Err num_f_zbrent(
    Err (*func)(double, double*, void*),
    double  x1,
    double  x2,
    double  f1,
    double  f2,
    double  tol,
    void*   static_data,
    double* answer)
/*	Using Brent's method, find the root of a function func known to lie between x1 and x2. The
        root, returned as zbrent, will be refined until its accuracy is tol. */
{
    Err    err = NULL;
    int    iter;
    double a = x1, b = x2, c = x2, d, e, min1, min2;
    double fa = f1, fb = f2, fc, p, q, r, s, tol1, xm;
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        return serror("Root must be bracketed in num_f_zbrent");
    fc = fb;
    for (iter = 1; iter <= ITMAX; iter++)
    {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
        {
            c  = a; /* Rename a, b, c and adjust bounding interval d */
            fc = fa;
            e = d = b - a;
        }
        if (fabs(fc) < fabs(fb))
        {
            a  = b;
            b  = c;
            c  = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol; /* Convergence check */
        xm   = 0.5 * (c - b);
        if (fabs(xm) <= tol1 || fb == 0.0)
        {
            *answer = b;
            return NULL;
        }
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
        {
            s = fb / fa; /* Attempt inverse quadratic interpolation */
            if (a == c)
            {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            }
            else
            {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0)
                q = -q; /* Check whether in bounds */
            p    = fabs(p);
            min1 = 3.0 * xm * q - fabs(tol1 * q);
            min2 = fabs(e * q);
            if (2.0 * p < (min1 < min2 ? min1 : min2))
            {
                e = d; /* Accept interpolation */
                d = p / q;
            }
            else
            {
                d = xm; /* Interpolation failed, use bisection */
                e = d;
            }
        }
        else
        { /* Bounds decreasing too slowly, use bisection */
            d = xm;
            e = d;
        }
        a  = b; /* Move last best guess to a */
        fa = fb;
        if (fabs(d) > tol1) /* Evaluate new trial root */
            b += d;
        else
            b += (xm >= 0.0 ? fabs(tol1) : -fabs(tol1));
        err = (*func)(b, &fb, static_data);
        if (err)
            return err;
    }
    return serror("Maximum number of iterations exceeded in num_f_zbrent");
}
