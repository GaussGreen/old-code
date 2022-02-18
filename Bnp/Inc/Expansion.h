// Expansion.h : Automatic expansion pricing for stoch vol models

#ifndef __EXPANSION_H__
#define __EXPANSION_H__

#include "uterror.h"

Err sabr_test_get_d(double sigma, double alpha, double rho, double t, double** d, int nd);
Err sabr_test4_get_d(double sigma, double alpha, double rho, double t, double** d, int nd);
Err heston_test_get_d(
    double sigma, double gamma, double alpha, double rho, double t, double* d, int nd);
Err sabr_test_get_all_d(double sigma, double alpha, double rho, double t, double** d, int nd);
Err heston_test_get_all_d(
    double sigma, double gamma, double alpha, double rho, double t, double** d, int nd);

Err calc_impvols(
    double** d,
    int*     uplus,
    int*     bplus,
    int*     splus,
    int      m,
    int      n,
    double   eps,
    double*  yk,
    int      nyk,
    double** vols);

Err sabr_test_mc(
    double  sigma,
    double  alpha,
    double  rho,
    double  t,
    double  fwd,
    double* K,
    int     nK,
    int     nstp,
    long    npaths,
    int     cf_ord,
    double* pv,
    double* stddev);

Err sabr_beta_test_mc(
    double  sigma,
    double  alpha,
    double  beta,
    double  rho,
    double  t,
    double  fwd,
    double* K,
    int     nK,
    int     nstp,
    long    npaths,
    int     cf_ord,
    long    seed,
    double* pv,
    double* stddev);

Err sabr_beta_test_mc4(
    double  sigma,
    double  alpha,
    double  beta,
    double  rho,
    double  t,
    double  fwd,
    double* K,
    int     nK,
    int     nstp,
    long    npaths,
    int     cf_ord,
    long    seed,
    double* pv,
    double* stddev);

Err sabr_beta_newexp(
    double  sigma,
    double  alpha,
    double  beta,
    double  rho,
    double  mat,
    double  fwd,
    double* K,
    int     nK,
    double* res);

#endif  // #ifndef __EXPANSION_H__