/* ======================================================
   FILENAME:  num_h_gausslaguerre.h

   PURPOSE:   integration by Gauss Laguerre (NUMC p152)
   ====================================================== */

#ifndef NUM_H_GAUSSLAGUERRE_H
#define NUM_H_GAUSSLAGUERRE_H

void gaulag(double x[], double w[], int n, double alpha);

void GaussLag(double *x, double *w, int n, double alpha);

#endif