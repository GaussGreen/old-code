/* ======================================================
   FILENAME:  num_h_spline.h

   PURPOSE:   Spline
   ====================================================== */

#ifndef NUM_H_SPLINE_H
#define NUM_H_SPLINE_H

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splie2(double x1a[], double x2a[], double **ya, int m, int n,
            double **y2a);

#endif
