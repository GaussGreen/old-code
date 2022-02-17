/*	Delivery options */

#include "math.h"
#include "utallhdr.h"
#include <NUM_H_ALLHDR.H>
#include <OPFNCTNS.H>

double opdeliv(int num_bonds,      /*	Number of bonds in the basket        ,
                                      index 0 = cheapest */
               double *fwds,       /*	Forwards */
               double *conv_facts, /*	Conversion factors */
               double **cov,       /*	Covariance matrix */
               double mat,         /*	Maturity */
               int log_or_norm,    /*	0: norm        , 1: log */
               long npth)          /*	Number of Sobol points */
{
  int *dir = (int *)calloc(num_bonds, sizeof(int));
  int i, j;
  double ***g;
  double *val = (double *)calloc(num_bonds, sizeof(double)),
         *var = (double *)calloc(num_bonds, sizeof(double));
  double min, res;

  for (i = 0; i < num_bonds; i++) {
    dir[i] = 1;
    for (j = 0; j < num_bonds; j++) {
      cov[i][j] *= mat;
    }
    var[i] = 0.5 * cov[i][i];
  }

  g = gen_sobol_mat(num_bonds, cov, dir, npth);

  res = 0.0;
  for (j = 0; j < npth; j++) {
    if (log_or_norm) {
      for (i = 0; i < num_bonds; i++) {
        val[i] = fwds[i] * exp(g[j][0][i] - var[i]);
      }
    } else {
      for (i = 0; i < num_bonds; i++) {
        val[i] = fwds[i] + g[j][0][i];
      }
    }

    min = val[0] / conv_facts[0];
    for (i = 1; i < num_bonds; i++) {
      if (val[i] / conv_facts[i] < min) {
        min = val[i] / conv_facts[i];
      }
    }

    res += min / npth;
  }

  free(dir);
  free(val);
  free(var);
  free_f3tensor(g, 0, npth - 1, 0, 0, 0, num_bonds - 1);

  return res;
}
