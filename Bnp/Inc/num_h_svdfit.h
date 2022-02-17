/* ======================================================
   FILENAME:  num_h_svdfit.h

   PURPOSE:   ???
   ====================================================== */

#ifndef NUM_H_SVDFIT_H
#define NUM_H_SVDFIT_H

Err svdfit(double x[], double y[], double sig[], int ndata, double a[], int ma,
           double **u, double **v, double w[], double *chisq,
           void (*funcs)(double, double[], int));

#endif
