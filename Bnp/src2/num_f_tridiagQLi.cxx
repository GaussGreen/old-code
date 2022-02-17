/* ========================================================================

  MODULE:	num_f_tridiagQLi.cxx

  LIBRARY:	NUM_LIB

  FUNCTION:	tridiagonal_QL_implicit

  DESCRIPTION:	From Numerical Recipes in C page 480
                QL algorithm with implicit shifts        , to determine the
eigen values and eigen vectors of a real        , symmetric , tridiagonal matrix
, previously reduced by the Householder method (cf uthouseholder) On input ,
d[0..n-1][0..n-1] contains the diagonal elements of a tridiagonal matrix. On
output        , it returns the iegen values. The vector e[0..n-1] inputs the sub
diagonal elements of the matrix        , with e[0] arbitrary. On output        ,
e[...] is destroyed If the eigenvectors of a tridiagonal matrix are desired ,
the matrix z[0..n-1][0..n-1] is input as the identity matrix. If the eigen
vectors of a matrix that has been reduced by Householder method
(householder_tridiagonalisation) are required        , than z is input as the
matrix output by this function. The kth column of z (i.e. z[...][k] returns the
                                normalised eigenvector corresponding to d[k]).

========================================================================== */

#include "math.h"
#include "num_h_pythag.h"
#include "utallhdr.h"

Err tridiagonal_QL_implicit(double d[], double e[], int n, double **z) {
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;
  Err err = NULL;

  for (i = 1; i < n; i++)
    e[i - 1] = e[i];
  e[n - 1] = 0.0;
  for (l = 0; l < n; l++) {
    iter = 0;
    /* Look for a single small subdiagonal element to split the matrix */
    do {
      for (m = l; m <= n - 2; m++) {
        dd = fabs(d[m]) + fabs(d[m + 1]);
        if ((double)(fabs(e[m]) + dd) == dd)
          break;
      }
      if (m != l) {
        if (iter++ == 30)
          return serror("Too many iterations in tridiagonal_QL_implicit");
        g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        r = pythag(g, 1.0);
        /* this is dm - ks */
        g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
        s = c = 1.0;
        p = 0.0;
        /* A plane rotation        , as in the original QL        , followed by
           Givens rotation to restore tridiagonal form */
        for (i = m - 1; i >= l; i--) {
          f = s * e[i];
          b = c * e[i];
          e[i + 1] = (r = pythag(f, g));
          if (r == 0.0) {
            d[i + 1] -= p;
            e[m] = 0.0;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          d[i + 1] = g + (p = s * r);
          g = c * r - b;
          for (k = 0; k < n; k++) {
            f = z[k][i + 1];
            z[k][i + 1] = s * z[k][i] + c * f;
            z[k][i] = c * z[k][i] - s * f;
          }
        }
        if (r == 0.0 && i >= l)
          continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }

    } while (m != l);

  } /* END for (l = 0 ; l < n ; l++ ) loop */

  /* Return a success message */
  return NULL;

} /* END Err tridiagonal_QL_implicit(...) */

/* (C) Copr. 1986-92 Numerical Recipes Software +135[)6=. */

Err tridag(double a[], double b[], double c[], double r[], double u[],
           unsigned long n) {
  unsigned long j;
  double bet, *gam;

  gam = dvector(1, n);
  if (b[1] == 0.0)
    serror("Error 1 in tridag");
  u[1] = r[1] / (bet = b[1]);
  for (j = 2; j <= n; j++) {
    gam[j] = c[j - 1] / bet;
    bet = b[j] - a[j] * gam[j];
    if (bet == 0.0)
      serror("Error 2 in tridag");
    u[j] = (r[j] - a[j] * u[j - 1]) / bet;
  }
  for (j = (n - 1); j >= 1; j--)
    u[j] -= gam[j + 1] * u[j + 1];
  free_dvector(gam, 1, n);

  return NULL;
}
