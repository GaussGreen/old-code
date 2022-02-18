/* ======================================================
   FILENAME:  num_h_interp.h

   PURPOSE:   linear interpolation functions
   ====================================================== */

#ifndef NUM_H_INTERP_H
#define NUM_H_INTERP_H

double interp(double* x, double* y, int l, double xt, int method, double* q);

double interp_columns(double** x, double** y, int l, double xt, int method, double* q, const int n);

double lin_interp_2d(
    double x, double y, double* xa, double* ya, double** za, long num_xa, long num_ya);

void hunt(double xx[], unsigned long n, double x, unsigned long* jlo);

double InterpLagrange1D(double* Xdata, double* Ydata, unsigned int n, double x);

double InterpLagrangeColumns1D(double** Xdata, double** Ydata, unsigned int n, double x, int m);

int BinarySearch(double* xx, int n, double x);

#endif
