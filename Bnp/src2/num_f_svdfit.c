/* ====================================================================
   FILENAME		: num_f_svdfit.c

   PURPOSE		: Singular Value Decomposition
                  Find a linear combination of predefined functions
                                  that best matches a set ( x, y=F(x) ) of points

   DESCRIPTION:	Given a set of data points x[1..ndata], y[1..ndata] with
                individual standard deviations sig[1..ndata], uses
                chi square minimisation to determine the coefficients a[1..ma]
                of the fitting function y = Sum(i) a_i * afunc_i(x)
                Arrays u[1..ndata][1..ma] , v[1..ma][1..ma], w[1..ma] provide
                workspace on input; on output they define the singular value
                decomposition and can be used to obtain the covariance matrix.

   RETURN:	The program returns values for the ma fit parameters a and the
                associated chi square.

   INPUT:	The user supplies a routine funcs(x,afunc,ma) that returns the
                ma basis functions evaluated at x in the array afunc[1..ma]

  ========================================================================= */

#include "utallhdr.h"

#define NRANSI
#define TOL 1.0e-5

Err svdfit(
    double   x[],
    double   y[],
    double   sig[],
    int      ndata,
    double   a[],
    int      ma,
    double** u,
    double** v,
    double   w[],
    double*  chisq,
    void (*funcs)(double, double[], int))
{
    void   svbksb(double** u, double w[], double** v, int m, int n, double b[], double x[]);
    Err    err, svdcmp(double**a, int m, int n, double w[], double**v);
    int    j, i;
    double wmax, tmp, thresh, sum, *b, *afunc;

    b     = vector(1, ndata);
    afunc = vector(1, ma);
    for (i = 1; i <= ndata; i++)
    {
        (*funcs)(x[i], afunc, ma);
        tmp = 1.0 / sig[i];
        for (j = 1; j <= ma; j++)
            u[i][j] = afunc[j] * tmp;
        b[i] = y[i] * tmp;
    }
    err = svdcmp(u, ndata, ma, w, v);
    if (err)
        return err;
    wmax = 0.0;
    for (j = 1; j <= ma; j++)
        if (w[j] > wmax)
            wmax = w[j];
    thresh = TOL * wmax;
    for (j = 1; j <= ma; j++)
        if (w[j] < thresh)
            w[j] = 0.0;
    svbksb(u, w, v, ndata, ma, b, a);
    *chisq = 0.0;
    for (i = 1; i <= ndata; i++)
    {
        (*funcs)(x[i], afunc, ma);
        for (sum = 0.0, j = 1; j <= ma; j++)
            sum += a[j] * afunc[j];
        *chisq += (tmp = (y[i] - sum) / sig[i], tmp * tmp);
    }
    free_vector(afunc, 1, ma);
    free_vector(b, 1, ndata);
    return NULL;
}
#undef TOL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software koV219. */
