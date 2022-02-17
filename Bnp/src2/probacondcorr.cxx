
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"

/*	Numerical calculation of conditional correlation        ,
        i.e. correl ( X         , Y / abs ( X * Y ) > k
        where X and Y are Gaussians with correl rho */
Err proba_cond_corr(double mean1, double mean2, double std1, double std2,
                    double rho, double k, int npth, double *empave1,
                    double *empave2, double *empexpx1sqr, double *empexpx2sqr,
                    double *empexpx1x2) {
  Err err = NULL;
  int i, n;
  double rho2, g1, g2;
  double ***g = NULL;
  double mu1, mu2, covar, var1, var2;

  g = f3tensor(0, npth - 1, 0, 0, 0, 1);
  if (!g) {
    goto FREE_RETURN;
  }

  /*	Initialise Sobol sequences */
  sobol_init(0, npth - 1, 0, 0, 0, 1);
  sobol_cube(g, 0, npth - 1, 0, 0, 0, 1);

  mu1 = mu2 = covar = var1 = var2 = 0.0;
  n = 0;
  rho2 = sqrt(1.0 - rho * rho);
  for (i = 0; i < npth - 1; i++) {
    g1 = g[i][0][0];
    g2 = rho * g1 + rho2 * g[i][0][1];

    g1 = mean1 + std1 * g1;
    g2 = mean2 + std2 * g2;

    if (fabs(g1 * g2) > k) {
      n++;
      mu1 += g1;
      mu2 += g2;
      covar += g1 * g2;
      var1 += g1 * g1;
      var2 += g2 * g2;
    }
  }

  if (!n) {
    err = "No path in conditional event";
    goto FREE_RETURN;
  }

  mu1 /= n;
  mu2 /= n;
  var1 /= n;
  var2 /= n;
  covar /= n;

  *empave1 = mu1;
  *empave2 = mu2;
  *empexpx1sqr = var1;
  *empexpx2sqr = var2;
  *empexpx1x2 = covar;

FREE_RETURN:

  if (g)
    free_f3tensor(g, 0, npth - 1, 0, 0, 0, 1);
  sobol_free();
  return err;
}
