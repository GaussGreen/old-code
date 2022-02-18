/* ======================================================
   FILENAME:  num_h_splint.h

   PURPOSE:   SPlint
   ====================================================== */

#ifndef SPLINT_H
#define SPLINT_H

Err splint(double xa[], double ya[], double y2a[], int n, double x, double* y);

void splin2(
    double   x1a[],
    double   x2a[],
    double** ya,
    double** y2a,
    int      m,
    int      n,
    double   x1,
    double   x2,
    double*  y);

#endif
