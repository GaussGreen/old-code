/*******************************************************************************
**
**	lightable_swap_fct.c
**
**      Authors:        K.Chau
**			J.Malhi
*******************************************************************************/

/* ==========================================================================
   include files
   ========================================================================== */

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/* ==========================================================================
   lightable_swap_fct(...)
   ========================================================================== */

double lightable_swap_fct(int num_path, double y20, double sy2, double y30,
                          double sy3, double rho, int compd, int nfp,
                          double coupon, double spread, double t, double df,
                          int condition) {

  long seed;
  double mc = 0;
  double sy2t, sy3t, sy22t, sy32t, disc, e1, e2, y2, y3;
  int i, j;
  double lp;

  seed = RANDINIT;

  sy2t = sy2 * sqrt(t);
  sy3t = sy3 * sqrt(t);
  sy22t = sy2t * sy2t / 2;
  sy32t = sy3t * sy3t / 2;
  for (i = 0; i < num_path; i++) {
    e1 = gauss_sample(&seed);
    e2 = rho * e1 + gauss_sample(&seed) * sqrt(1 - rho * rho);
    y2 = y20 * exp(-sy22t + sy2t * e1);
    y3 = y30 * exp(-sy32t + sy3t * e2);
    lp = 1.0;
    if (((y3 - y2 > spread) && (condition == 1)) ||
        ((y3 - y2 < spread) && (condition == 0))) {
      disc = 1 / (1 + y2 / compd);
      for (j = 1; j < nfp; j++) {
        lp = 1 + disc * lp;
      }
      lp *= disc;
      mc += (coupon - y2) * lp;
    }
  }

  return df * mc / (compd * num_path);
}
