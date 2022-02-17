/* ==========================================================================
   FILE_NAME:	num_f_chebyshev.c

   PURPOSE:     a few utilities using Chebyshev polynomial
                (cf NRC pages 190ff)
   ========================================================================== */

#include "utallhdr.h"

#define CHEBYSHEV_BIG 1.0e+50

/* --------------------------------------------------------------------------
   Chebyshev evaluation : all arguments are input c[0..m-1] is an array of
   Chebyshev coefficients.
   -------------------------------------------------------------------------- */
double chebev(double a, double b, double c[], int m, double x) {
  double d = 0.0, dd = 0.0, sv, y, y2;
  int j;

  if ((x - a) * (x - b) > 0.0)
    return CHEBYSHEV_BIG;

  /* Change of variable */
  y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a));

  /* Clenshaw's recurrence */
  for (j = m - 1; j >= 1; j--) {
    sv = d;
    d = y2 * d - dd + c[j];
    dd = sv;
  }

  return y * d - dd + 0.5 * c[0];
}

#undef CHEBYSHEV_BIG

/* (C) Copr. 1986-92 Numerical Recipes Software. */
