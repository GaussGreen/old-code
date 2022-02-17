/* ======================================================
   FILENAME:  num_h_chebyshev.h
   
   PURPOSE:   Broydn 
   ====================================================== */

#ifndef NUM_H_CHEBYSHEV_H
#define NUM_H_CHEBYSHEV_H

/* --------------------------------------------------------------------------
   Chebyshev evaluation : all arguments are input c[0..m-1] is an array of 
   Chebyshev coefficients.
   -------------------------------------------------------------------------- */
double chebev(double a, double b, double c[], int m, double x);

#endif



