/* ======================================================
   FILENAME:  num_h_sobenberg.h

   PURPOSE:   Minimalisation using Sobenberg method
   ====================================================== */

#ifndef NUM_H_SOBENBERG_H
#define NUM_H_SOBENBERG_H

Err sobenberg(double *data, double *target, double *weight, long ndata,
              double *param, double *min_param, double *max_param, long nparam,
              long tot_pts, long bst_pts, long nclust, char *rsc_mth,
              long niter, Err (*funcs)(double, double[], double *, int),
              Err (*dfuncs)(double, double[], double *, double[], int),
              double *chisq);

#endif
