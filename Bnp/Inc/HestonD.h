// HestonD.h : New Heston implementation

#ifndef __HESTOND_H__
#define __HESTOND_H__

Err HestonDOptions(double f0, double sigma, double alpha, double rho,
                   double lam, double shift, double mat, int nK, double *K,
                   char **rec_pay_str, int *want, double **res);

Err HestonDCalibrate(double f0, double volf0, double K1, double vol1, double K2,
                     double vol2, double mat, double *sigma, double *alpha,
                     double *rho, double lam, double shift, int calib_smile);

Err HestonDDensityFT(double sigma, double alpha, double rho, double lam,
                     double mat, double u_re, double u_im, double *h_re,
                     double *h_im);

#endif // #ifndef __HESTOND_H__