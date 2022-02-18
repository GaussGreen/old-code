#ifndef OPSABRCALIB_H
#define OPSABRCALIB_H

/* Calibrates the SABR parameters to a market smile */

Err opsabrcalib(
    double  Fwd,
    double  maturity,
    int     nbr_strikes,
    double* strikes,
    double* marketvols,
    double* ATMVol,
    double* alpha,
    int     freeze_alpha, /* if 0, alpha is not calibrated */
    double* beta,
    int     freeze_beta, /* If 0, beta is not calibrated */
    double* rho,
    int     freeze_rho, /* if 0, rho is not calibrated */
    double* fitting_error);

/**********************************************************************
The following functions compute SABR vol and some derivatives as an
explicit function of the ATM volatility. They come from a cpp file
***********************************************************************/

/* SABR lognormal Vol */
double SabrLognormalVol(
    double        vol,
    double        alpha,
    double        beta,
    double        rho,
    double        forward,
    double        strike,
    double        maturity,
    SABR_VOL_TYPE cst_vol);

/* Derivative with respect to the vovol alpha */
double SabrDsigmaDalpha(
    double        vol,
    double        alpha,
    double        beta,
    double        rho,
    double        forward,
    double        strike,
    double        maturity,
    SABR_VOL_TYPE cst_vol);

/* Derivative with respect to the power parameter beta */
double SabrDsigmaDbeta(
    double        vol,
    double        alpha,
    double        beta,
    double        rho,
    double        forward,
    double        strike,
    double        maturity,
    SABR_VOL_TYPE cst_vol);

/* Derivative with respect to the correlation parameter rho */
double SabrDsigmaDrho(
    double        vol,
    double        alpha,
    double        beta,
    double        rho,
    double        forward,
    double        strike,
    double        maturity,
    SABR_VOL_TYPE cst_vol);

/* Derivative with respect to the ATM lognormal volatility */
double SabrDsigmaDvolatm(
    double vol,
    double alpha,
    double beta,
    double rho,
    double forward,
    double strike,
    double maturity);

Err compute_ATMVol(
    double  forward,
    double* strikes,
    double* market_vols,
    int     nbr_strikes,
    int*    freeze_ATMVol,
    double* ATMVol);

Err compute_weights(
    double* strikes,
    int     nbr_strikes,
    double  forward,
    double  ATMVol,
    double  maturity,
    double* weights_of_strikes);

/**************************************************************************
*************************************************************************/

#endif
