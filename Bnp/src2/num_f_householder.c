/* ========================================================================

  MODULE:	num_f_householder.C

  LIBRARY:	NUM_LIB

  FUNCTION:	householder_tridiagonalisation

  INPUTS:


  DESCRIPTION:	From Numerical Recipes in C page 474
                This function performs the Householder reduction of a real  ,
                                symmetric matrix a[0..n-1][0..n-1]. On output  ,
the matrix a is replaced by the orthogonal matrix effecting the transformation
                                d[0..n-1] returns the diagonal elements of the
tridiagonal matrix  , and e[0..n-1] the off diagonal elements  , with e[0] = 0

========================================================================== */

#include "math.h"
#include "utallhdr.h"

Err householder_tridiagonalisation(double **a, int n, double d[], double e[]) {
  int l, k, j, i;
  double scale, hh, h, g, f;

  for (i = n - 1; i >= 1; i--) {
    l = i - 1;
    h = scale = 0.0;
    if (l > 0) {
      for (k = 0; k <= l; k++)
        scale += fabs(a[i][k]);

      /* Skip trnasformation */
      if (scale == 0.0)
        e[i] = a[i][l];
      else {
        for (k = 0; k <= l; k++) {
          /* Use scaled a's for transformation */
          a[i][k] /= scale;
          /* Form sigma in h */
          h += a[i][k] * a[i][k];
        }

        f = a[i][l];
        g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i] = scale * g;
        /* Now h is equation 11.2.4 */
        h -= f * g;
        /* Store u in the ith row of a */
        a[i][l] = f - g;
        f = 0.0;
        for (j = 0; j <= l; j++) {
          /* Store u/H in the ith column of a */
          a[j][i] = a[i][j] / h;
          /* Form an element of A.u in g */
          g = 0.0;
          for (k = 0; k <= j; k++)
            g += a[j][k] * a[i][k];

          for (k = j + 1; k <= l; k++)
            g += a[k][j] * a[i][k];
          /* Form element of p in temporarily unused element of e */
          e[j] = g / h;
          f += e[j] * a[i][j];
        }
        /* Form K  , equation 11.2.11 */
        hh = f / (h + h);
        /* Form q and store in e overwriting p */
        for (j = 0; j <= l; j++) {
          f = a[i][j];
          e[j] = g = e[j] - hh * f;
          /* Reduce a  , equation 11.2.13 */
          for (k = 0; k <= j; k++)
            a[j][k] -= (f * e[k] + g * a[i][k]);
        }
      } /* END if scale == 0.00 */

    } /* END if l > 0 */

    else
      e[i] = a[i][l];
    d[i] = h;
  } /* END for i loop  */

  d[0] = 0.0;
  e[0] = 0.0;

  /* Get eigenvectors */
  for (i = 0; i < n; i++) {
    /* Begin accumulation of transformation matrices */
    l = i - 1;
    if (d[i]) {
      for (j = 0; j <= l; j++) {
        g = 0.0;
        /* Use u and u/H stored in a to form P.Q */
        for (k = 0; k <= l; k++)
          g += a[i][k] * a[k][j];
        for (k = 0; k <= l; k++)
          a[k][j] -= g * a[k][i];
      }
    }
    d[i] = a[i][i];
    /* Reset row and column of a to identity matrix for next iteration */
    a[i][i] = 1.0;
    for (j = 0; j <= l; j++)
      a[j][i] = a[i][j] = 0.0;
  }

  /* Return a success message */
  return NULL;

} /* END Err householder_tridiagonalisation(...) */

/* (C) Copr. 1986-92 Numerical Recipes Software +135[)6=. */
/* ------------------------------------------------------------------------- */