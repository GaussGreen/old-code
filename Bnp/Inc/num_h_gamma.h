/* ======================================================
   FILENAME:  num_h_gamma.h

   PURPOSE:   Broydn
   ====================================================== */
#ifndef NUM_H_GAMMA_H
#define NUM_H_GAMMA_H

/* Header file for gamma functions */

double gammln(double xx);

double gammp(double a, double x);

void gcf(double *gammcf, double a, double x, double *gln);

void gser(double *gamser, double a, double x, double *gln);

double gammq(double a, double x);

double betaincmp(double a, double b, double x);

double betacmp(double alpha, double beta);

double fact(long n);
#endif