/* ===============================================================
   FILE_NAME:	num_f_spline.cxx

   PURPOSE:     spline a function
   =============================================================== */

#include "utallhdr.h"

#define NRANSI

void spline(double x[], double y[], int n, double yp1, double ypn,
            double y2[]) {
  int i;
  int k;
  double p;
  double qn;
  double sig;
  double un;
  double *u;

  u = vector(1, n - 1);

  if (yp1 > 0.99e30) {
    y2[1] = u[1] = 0.0;
  } else {
    y2[1] = -0.5;
    u[1] = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1);
  }

  for (i = 2; i <= (n - 1); i++) {
    sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
           (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
  }

  if (ypn > 0.99e30) {
    qn = un = 0.0;
  } else {
    qn = 0.5;
    un = (3.0 / (x[n] - x[n - 1])) *
         (ypn - (y[n] - y[n - 1]) / (x[n] - x[n - 1]));
  }

  y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0);

  for (k = (n - 1); k >= 1; k--) {
    y2[k] = y2[k] * y2[k + 1] + u[k];
  }

  free_vector(u, 1, n - 1);
}
/******************************************************
 * Routine to compute the second derivatives
 * in in the x (or y) direction for a given matrix of 2D points.
 * C. Godart: 6/01/99
 * See also: splin2 to actually perform the interpolation
 *
 *splie2 return y2a a matrix of doubles
 *representing the second partial derivatives in one direction
 *for a matrix of point representing a surface.
 *This is a m*n call in term of complexity.
 ******************************************************/

void splie2(double x1a[], double x2a[], double **ya, int m, int n,
            double **y2a) {
  int j;
  /* remember that setting 1.0e30 actually sets the second derivative
   at the boundary to be null */
  for (j = 1; j <= m; j++)
    spline(x2a, ya[j], n, 1.0e30, 1.0e30, y2a[j]);
}

#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 0>)"?. */
