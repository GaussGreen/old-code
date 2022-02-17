/* ===============================================================
   FILENAME : num_h_clustering.h

   PURPOSE:   clustering algorithm
   =============================================================== */

#ifndef NUM_H_CLUSTERING_H
#define NUM_H_CLUSTERING_H

double dist_sqr(double *x, double *y, long DIM);

int clustering(double **pts, double **bar, long npts, long nbar, long dim);

void find_initial_bars(double **pts, double **bar, long npts, long nbar,
                       long dim);

#endif
