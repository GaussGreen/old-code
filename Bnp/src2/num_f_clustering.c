/* ===============================================================
   FILENAME : num_f_clustering.c

   PURPOSE:   clustering algorithm
   =============================================================== */

#include "utallhdr.h"

#define BIG 1.0e+16

double dist_sqr(double *x, double *y, long DIM) {

  long i;
  double res;

  res = 0.0;
  for (i = 1; i <= DIM; i++) {
    res += (x[i] - y[i]) * (x[i] - y[i]);
  }

  return res;
}

void find_initial_bars(double **pts, double **bar, long npts, long nbar,
                       long dim) {
  long i, j, k, pts_per_cluster = (long)(npts / nbar);

  for (j = 1; j <= nbar; j++) {
    memset(&(bar[j][1]), 0, dim * sizeof(double));
  }

  for (i = 1; i <= nbar; i++) {
    for (j = 1; j <= pts_per_cluster; j++) {
      for (k = 1; k <= dim; k++) {
        bar[i][k] += pts[(i - 1) * pts_per_cluster + j][k] / pts_per_cluster;
      }
    }
  }
}

/* do not end iterations before return = 1*/
int clustering(double **pts, double **bar, long npts, long nbar, long dim) {

  /* declarations */
  long i, j, k, test, *assoc_bar = lngvector(1, npts),
                      *npts_in_bar = lngvector(1, nbar);

  static long pt = 0;

  double temp, dist;

  /* compute distances and associate bars */

  memset(&(npts_in_bar[1]), 0, nbar * sizeof(long));

  for (i = 1; i <= npts; i++) {
    dist = BIG;
    for (j = 1; j <= nbar; j++) {
      temp = dist_sqr(pts[i], bar[j], dim);
      if (temp < dist) {
        dist = temp;
        assoc_bar[i] = j;
      }
    }
    npts_in_bar[assoc_bar[i]]++;
  }

  test = 0;
  for (j = 1; j <= nbar; j++) {
    if (!npts_in_bar[j]) {
      test = 1;
      pt++;
      if (pt > npts) {
        pt = 1;
      }

      for (k = 1; k <= dim; k++) {
        bar[j][k] = pts[pt][k];
      }
    }
  }

  if (test) {
    free_lngvector(assoc_bar, 1, npts);
    free_lngvector(npts_in_bar, 1, nbar);
    return 0;
  }

  for (j = 1; j <= nbar; j++) {
    memset(&(bar[j][1]), 0, dim * sizeof(double));
  }

  for (i = 1; i <= npts; i++) {
    for (k = 1; k <= dim; k++) {
      bar[assoc_bar[i]][k] += pts[i][k] / npts_in_bar[assoc_bar[i]];
    }
  }

  /* free memory */

  free_lngvector(assoc_bar, 1, npts);
  free_lngvector(npts_in_bar, 1, nbar);
  return 1;
}
