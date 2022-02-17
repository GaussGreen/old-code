/* ==============================================================================
   FILENAME:  num_h_bestfit.h

   PURPOSE:   multi-dimensional quadratic fit
   ==============================================================================
 */

#ifndef NUM_H_BESTFIT_H
#define NUM_H_BESTFIT_H

Err quadr_best_fit(double **x, double *y, long n, long dim, double *a,
                   double *grad, double **hess, double *min);

#endif
