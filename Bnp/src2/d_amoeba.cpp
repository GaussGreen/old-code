// Sample page from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
// Copyright (C) 1988-1992 by Cambridge University Press. Programs Copyright (C) 1988-1992 by
// Numerical Recipes Software.

#include "d_amoeba.h"

#include "math.h"
#include "srt_h_all.h"

#define TINY 1.0e-8  // A small number.
#define NMAX 5000    // Maximum allowed number of function evaluations.
#define NDIMMAX 50   // Maximum number of dimensions

#define GET_PSUM                              \
    for (j = 0; j < ndim; j++)                \
    {                                         \
        for (sum = 0.0, i = 0; i < mpts; i++) \
            sum += p[i][j];                   \
        psum[j] = sum;                        \
    }

#define SWAP(a, b)   \
    {                \
        swap = (a);  \
        (a)  = (b);  \
        (b)  = swap; \
    }

static Err amotry(
    double** p,
    double*  y,
    double*  psum,
    int      ndim,
    Err (*funk)(double*, double*, void*),
    int     ihi,
    double  fac,
    double* res,
    void*   comm)
// Extrapolates by a factor fac through the face of the simplex across from the high point, tries
// it, and replaces the high point if the new point is better.
{
    Err    err = NULL;
    int    j;
    double fac1, fac2, ytry, ptry[NDIMMAX];
    fac1 = (1.0 - fac) / ndim;
    fac2 = fac1 - fac;
    for (j = 0; j < ndim; j++)
        ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
    err = (*funk)(ptry, &ytry, comm);  // Evaluate the function at the trial point.
    if (err)
        return err;
    if (ytry < y[ihi])
    {  // If it’s better than the highest, then replace the highest.
        y[ihi] = ytry;
        for (j = 0; j < ndim; j++)
        {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
    }
    *res = ytry;
    return NULL;
}

Err d_amoeba(
    double** p,
    double*  y,
    int      ndim,
    double   ftol,
    Err (*funk)(double*, double*, void*),
    int*  nfunk,
    void* comm)
/* Multidimensional minimization of the function funk(x) where x[0..ndim-1] is a vector in ndim
dimensions, by the downhill simplex method of Nelder and Mead. The matrix p[0..ndim][0..ndim-1] is
input. Its ndim+1 rows are ndim-dimensional vectors which are the vertices of the starting simplex.
Also input is the vector y[0..ndim], whose components must be preinitialized to the values of funk
evaluated at the ndim+1 vertices (rows) of p; and ftol the fractional convergence tolerance to be
achieved in the function value (n.b.!). On output, p and y will have been reset to ndim+1 new points
all within ftol of a minimum function value, and nfunk gives the number of function evaluations
taken. comm is a pointer to a communication structure (external parameters that will be passed to
the objective function (3rd parameter)) */
{
    Err    err = NULL;
    int    i, ihi, ilo, inhi, j, mpts = ndim + 1;
    double rtol, sum, swap, ysave, ytry, psum[NDIMMAX];
    *nfunk = 0;

    GET_PSUM
    for (;;)
    {
        ilo = 0;
        // First we must determine which point is the highest (worst), next-highest, and lowest
        // (best), by looping over the points in the simplex.
        ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
        for (i = 0; i < mpts; i++)
        {
            if (y[i] <= y[ilo])
                ilo = i;
            if (y[i] > y[ihi])
            {
                inhi = ihi;
                ihi  = i;
            }
            else if (y[i] > y[inhi] && i != ihi)
                inhi = i;
        }
        rtol = 2.0 * fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]) + TINY);
        // Compute the fractional range from highest to lowest and return if satisfactory.
        if (rtol < ftol)
        {  // If returning, put best point and value in slot 0.
            SWAP(y[0], y[ilo])
            for (i = 0; i < ndim; i++)
                SWAP(p[0][i], p[ilo][i])
            break;
        }
        if (*nfunk >= NMAX)
            return serror("NMAX exceeded in d_amoeba");
        *nfunk += 2;
        // Begin a new iteration. First extrapolate by a factor -1 through the face of the simplex
        // across from the high point, i.e., reflect the simplex from the high point.

        err = amotry(p, y, psum, ndim, funk, ihi, -1.0, &ytry, comm);
        if (err)
            return err;

        if (ytry <= y[ilo])
        // Gives a result better than the best point, so try an additional extrapolation by a
        // factor 2.
        {
            err = amotry(p, y, psum, ndim, funk, ihi, 2.0, &ytry, comm);
            if (err)
                return err;
        }
        else if (ytry >= y[inhi])
        {
            // The reflected point is worse than the second-highest, so look for an intermediate
            // lower point, i.e., do a one-dimensional contraction.
            ysave = y[ihi];
            err   = amotry(p, y, psum, ndim, funk, ihi, 0.5, &ytry, comm);
            if (err)
                return err;

            if (ytry >= ysave)
            {   // Can’t seem to get rid of that high point. Better
                // contract around the lowest (best) point.
                for (i = 0; i < mpts; i++)
                {
                    if (i != ilo)
                    {
                        for (j = 0; j < ndim; j++)
                            p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
                        err = (*funk)(psum, &y[i], comm);
                        if (err)
                            return err;
                    }
                }
                *nfunk += ndim;  // Keep track of function evaluations.
                GET_PSUM         // Recompute psum.
            }
        }
        else
            --(*nfunk);  // Correct the evaluation count.
    }                    // Go back for the test of doneness and the next iteration.
    return NULL;
}
