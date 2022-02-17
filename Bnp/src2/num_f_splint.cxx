/* ===============================================================
   FILE_NAME:	num_f_splint.cxx

   PURPOSE:     spline a function
   =============================================================== */

#include "num_h_spline.h"
#include "utallhdr.h"

/**
  modified E.Auld to return Err
  (nrerror replaced with serror).
  The reason is that srt functions cannot be allowed
  to call exit()
**/
Err splint(double xa[], double ya[], double y2a[], int n, double x, double *y) {

  int klo;
  int khi;
  int k;

  double h;
  double b;
  double a;

  klo = 1;
  khi = n;

  while (khi - klo > 1) {
    k = (khi + klo) >> 1;

    if (xa[k] > x) {
      khi = k;
    } else {
      klo = k;
    }
  }

  h = xa[khi] - xa[klo];

  if (h == 0.0) {
    return serror("Bad xa input to routine splint");
  }

  a = (xa[khi] - x) / h;
  b = (x - xa[klo]) / h;
  *y =
      a * ya[klo] + b * ya[khi] +
      ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
  return NULL;
}
/******************************************************
 * splin2. Routine to perform bi-dimensional spline interpolation
 * C. Godart: 6/01/99
 *  See also: splie2.
 ******************************************************/

void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n,
            double x1, double x2, double *y) {
  int j;
  double *ytmp, *yytmp;

  ytmp = dvector(1, m);
  yytmp = dvector(
      1,
      m); /* Perform m evaluations of the row splines constructed by splie2 */
  for (j = 1; j <= m; j++) /* using the one-dimensional spline evaluator: */
    splint(x2a, ya[j], y2a[j], n, x2, &yytmp[j]); /* splint */
  /* So I get a column of yytmp */
  /* on which I perform a spline interpolation */
  /* first compute 1d ytmp coefficients */
  spline(x1a, yytmp, m, 1.0e30, 1.0e30, ytmp);
  /* Then interpolate */
  /* Note that this time for each spline interpolation
     I perform one spline computation: the algorithm in two dimension
     jumps from log(n) to m*log(n) */
  splint(x1a, yytmp, ytmp, m, x1, y);
  free_dvector(yytmp, 1, m);
  free_dvector(ytmp, 1, m);
}
/* (C) Copr. 1986-92 Numerical Recipes Software 0>)"?. */
