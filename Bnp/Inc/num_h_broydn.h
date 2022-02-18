/* ======================================================
   FILENAME:  num_h_broydn.h

   PURPOSE:   Broydn
   ====================================================== */

#ifndef NUM_H_BROYDN_H
#define NUM_H_BROYDN_H

char* broydn(double x[], int n, int* check, void (*vecfunc)(int, double[], double[]));

int zbrac(double (*func)(double), double* x1, double* x2);

#endif
