#ifndef OPSABRNUMER_H
#define OPSABRNUMER_H

#include "utError.h"
#include "utTypes.h"

Err op_sabr_adi(double forward, double strike, double maturity, double disc,
                int call_put_var, /*	0:	(F - K)+
                                                          1:	(K - F)+
                                                          2:	F^2	*/
                double sigma_beta, double alpha, double beta, double rho,
                int nt, int nx, int nz, double *result);

double
op_sabr_adi_simple(double forward, double strike, double maturity, double disc,
                   int call_put_var, /*	0:	(F - K)+
                                                             1:	(K - F)+
                                                             2:	F^2	*/
                   double sigma_beta, double alpha, double beta, double rho,
                   int nt, int nx, int nz);

Err op_sabr_mc_ln(double forward, int nb_strike, double *strike,
                  double maturity, double disc, double sigma, double alpha,
                  double rho, int nt, long npaths, double *res, double *volres,
                  double *fwdres);

Err op_sabr_mc_nor(double forward, int nb_strike, double *strike,
                   double maturity, double disc, double sigma, double alpha,
                   double rho, int nt, long npaths, double *res,
                   double *volres);

Err op_sabr_mc_beta_Euler(double forward, int nb_strike,
                          double *strike,  // [in]
                          double maturity, // Option tenor in years
                          int optionType,  // 1 for call  , -1 for put
                          double df,       // discount factor
                          double sigma,    // "sigma beta"
                          double alpha,    // annualized lognormal vol of sigma
                          double beta,     // in [0  ,1]
                          double rho,      // in (-1  , 1)
                          int nt,          // number of time points
                          long npaths,     // number of paths
                          long seed, double *value, double *sterr);

Err op_sabr_calib_adi(double forward, double strike, double maturity,
                      double tgt_vol, SrtDiffusionType input_vol_type,
                      double alpha, double beta, double rho, int nt, int nx,
                      int nz, int nbiter, double precision, double *res);

Err srt_ADINewRhoSABR(double forward, double maturity, double ATMVol,
                      double alpha, double beta, double rho, double newbeta,
                      int nt, int nx, int nz, int iNumStrikes,
                      double *pdNewRho);

Err srt_AlphaFromADI(double forward, double maturity, double ATMVol,
                     double alpha, double beta, double rho, double Disc,

                     int nt, int nx, int nz, int iNumStrikes,

                     double strikeMin, double strikeMax,

                     double *pdNewAlpha);

#endif