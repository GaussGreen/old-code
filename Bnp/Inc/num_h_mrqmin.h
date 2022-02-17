/* ======================================================
   FILENAME:  num_h_mrqmin.h

   PURPOSE:   for LevenebergMarquardt optimisation
   ====================================================== */
#ifndef NUM_H_MRQMIN_H
#define NUM_H_MRQMIN_H

Err mrqmin(double x[], double y[], double sig[], int ndata, double a[],
           int ia[], int ma, double **covar, double **alpha, double *chisq,
           Err (*funcs)(double, double[], double *, double[], int),
           double *alamda);

#endif
