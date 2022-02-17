// FxHestonD.h : FX Heston model

#ifndef __FXHESTOND_H__
#define __FXHESTOND_H__

typedef struct _SFxHestonUnd {
  char yc_dom[80], yc_for[80];
  double S0, shift;
  int nt;
  long D_star, today;
  double *t, *sigma, *alpha, *lam, *rho;
} SFxHestonUnd;

char *FreeFxHestonUnd(void *ptr);

Err FxHestonDOptions(SFxHestonUnd *und, long fix_date, long pay_date, int nK,
                     double *K, char **rec_pay_str, double *res);

Err NewtonD(double x, double x_min, double x_max,
            Err (*func)(double, double *, void *), double tol, int maxiter,
            void *comm);

Err FxHestonDCalibrate(char *yc_dom, char *yc_for, double spot_fx, long D_star,
                       double beta, int ndates, long *dates, double *lam,
                       double *volf0, double *K1, double *volK1, double *K2,
                       double *volK2, int *calib_smile,
                       // Output:
                       double *sigma, double *alpha, double *rho);

#endif // #ifndef __FXHESTOND_H__