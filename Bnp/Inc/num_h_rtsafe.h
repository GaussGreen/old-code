/* ======================================================
   FILENAME:  num_h_rtsafe.h

   PURPOSE:   Safe Newton
   ====================================================== */

#ifndef NUM_H_RTSAFE_H
#define NUM_H_RTSAFE_H

SrtErr rtsafe(
    Err (*func)(double, double*, double*),
    double  x1,
    double  x2,
    double  xacc,
    int     num_iter,
    double* answer);

SrtErr rtsafe_with_par(
    Err (*func)(double, double*, double*, double*),
    double  x1,
    double  x2,
    double  xacc,
    int     num_iter,
    double* answer,
    double* par);

#endif
