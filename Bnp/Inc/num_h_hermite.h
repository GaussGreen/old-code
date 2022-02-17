/* -------------------------------------------------------------
   FILENAME:   num_h_hermite.h

   PURPOSE:    provide the roots of the Hermite polynomials for
               a quick integration of a function f on the
                           gaussian density as a discrete sum on the w[]
               (cf Numerical Recipes in C p. 154)
   ------------------------------------------------------------- */

#ifndef NUM_H_HERMITE_H
#define NUM_H_HERMITE_H

Err gauss_hermite(double x[], double w[], int n);

Err store_hermite_in_file(int n);

Err hermite_gauss_quick(int n, double x[], double w[]);

Err HermiteStandard(double *x, double *w, int n);

#endif