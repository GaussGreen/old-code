/* ==========================================================
   FILENAME :        num_f_mrqmin.c

   PURPOSE:          for Levenberg-Marquardt
   ========================================================== */

#include "utallhdr.h"
#include <num_h_gaussj.h"

#define NRANSI

/* return err from other NR functions - K Chau 5/5/95 */

Err mrqmin(double x[],     /* Vector of data ( with y_i=f(x_i  ,a) ) */
           double y[],     /* Vector of target values to reach for f(x) */
           double sig[],   /* Weighting of each value */
           int ndata,      /* Number of data available (size of x and y )*/
           double a[],     /* Vector of parameters used to optimise */
           int ia[],       /* Flag to know whether to use the param or not */
           int ma,         /* Number of parameters (size of a) */
           double **covar, /* Covariance matrix for a ([1..ma][1..ma]) */
           double **alpha, /* Hessian matrix for F(...  ,a) ([1..ma][1..ma]) */
           double *chisq,  /* Value of the criteria after minimisation */
           Err (*funcs)(double, double[], double *, double[], int),
           double *alamda) /* Length of the search on the gradient line */
{
  void covsrt(double **covar, int ma, int ia[], int mfit);
  Err mrqcof(double x[], double y[], double sig[], int ndata, double a[],
             int ia[], int ma, double **alpha, double beta[], double *chisq,
             Err (*funcs)(double, double[], double *, double[], int));
  int j, k, l;
  Err err = NULL;
  static int mfit;
  static double ochisq, *atry, *beta, *da, **oneda;

  /* Initialization */
  if (*alamda < 0.0) {
    atry = dvector(1, ma);
    beta = dvector(1, ma);
    da = dvector(1, ma);
    for (mfit = 0, j = 1; j <= ma; j++)
      if (ia[j])
        mfit++;
    oneda = dmatrix(1, mfit, 1, 1);
    *alamda = 0.001;
    err = mrqcof(x, y, sig, ndata, a, ia, ma, alpha, beta, chisq, funcs);
    if (err) {
      free_dmatrix(oneda, 1, mfit, 1, 1);
      free_dvector(da, 1, ma);
      free_dvector(beta, 1, ma);
      free_dvector(atry, 1, ma);
      return err;
    }
    ochisq = (*chisq);
    for (j = 1; j <= ma; j++)
      atry[j] = a[j];
  }
  /* Alter linearized fitting matrix by augmenting diagonal elements */
  for (j = 1; j <= mfit; j++) {
    /* Alteration in the new version of the software */
    for (k = 1; k <= mfit; k++) {
      covar[j][k] = alpha[j][k];
    }
    covar[j][j] = alpha[j][j] * (1.0 + (*alamda));
    oneda[j][1] = beta[j];
  }

  /* Add a check to prevent useless Gaussj-2 when a coeff has no impact */
  for (j = 1; j <= mfit; j++) {
    /* All the chisqr partial derivatives for this coeff are null */
    if (covar[j][j] == 0.0) {
      /* Replace covar by 1.0 to solve 1.0 * Dx = 0 rather than : 0 * Dx = 0 */
      covar[j][j] = 1.0;
      oneda[j][1] = 0;
    }
  }

  /* Matrix solution */
  if (err = gaussj(covar, mfit, oneda, 1)) {
    free_dmatrix(oneda, 1, mfit, 1, 1);
    free_dvector(da, 1, ma);
    free_dvector(beta, 1, ma);
    free_dvector(atry, 1, ma);
    return err;
  }

  /* Once converged  , evaluate covariance matrix */
  for (j = 1; j <= mfit; j++)
    da[j] = oneda[j][1];
  if (*alamda == 0.0) {
    covsrt(covar, ma, ia, mfit);
    free_dmatrix(oneda, 1, mfit, 1, 1);
    free_dvector(da, 1, ma);
    free_dvector(beta, 1, ma);
    free_dvector(atry, 1, ma);
    return NULL;
  }

  /* Did the trial succeed ?*/
  for (j = 0, l = 1; l <= ma; l++)
    if (ia[l])
      atry[l] = a[l] + da[++j];
  err = mrqcof(x, y, sig, ndata, atry, ia, ma, covar, da, chisq, funcs);
  if (err) {
    free_dmatrix(oneda, 1, mfit, 1, 1);
    free_dvector(da, 1, ma);
    free_dvector(beta, 1, ma);
    free_dvector(atry, 1, ma);
    return err;
  }
  /* Success: accept the new solution */
  if (*chisq < ochisq) {
    *alamda *= 0.1;
    ochisq = (*chisq);
    /* Alteration in the new version of the software */
    for (j = 1; j <= mfit; j++) {
      for (k = 1; k <= mfit; k++)
        alpha[j][k] = covar[j][k];
      beta[j] = da[j];
    }
    for (l = 1; l <= ma; l++)
      a[l] = atry[l];

  } else
  /* Failure: reject the new solution and increase alambda */
  {
    *alamda *= 10.0;
    *chisq = ochisq;
  }
  return NULL;
}

#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 0>)"?. */
