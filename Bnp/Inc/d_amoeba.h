// Sample page from NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
// Copyright (C) 1988-1992 by Cambridge University Press. Programs Copyright (C)
// 1988-1992 by Numerical Recipes Software.

#ifndef __D_AMOEBA_H__
#define __D_AMOEBA_H__

Err d_amoeba(double **p, double *y, int ndim, double ftol,
             Err (*funk)(double *, double *, void *), int *nfunk, void *comm);

#endif // #ifndef __D_AMOEBA_H__
