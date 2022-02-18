/* ============================================================

  FILENAME:	    num_f_fmin.c


  PURPOSE:	    diagonalise_symmetric_matrix
  ============================================================= */
#define NRANSI
#include "utallhdr.h"

EXTERN int     nn;
EXTERN double* fvec;
EXTERN void (*nrfuncv)(int n, double v[], double f[]);

double fmin_c(double x[])
{
    int    i;
    double sum;

    (*nrfuncv)(nn, x, fvec);
    for (sum = 0.0, i = 1; i <= nn; i++)
        sum += SQR(fvec[i]);
    return 0.5 * sum;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software +135[)6=. */
