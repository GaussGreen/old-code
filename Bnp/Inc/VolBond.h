#ifndef __VOLBOND_H
#define __VOLBOND_H

#include "VolBondProdStruct.h"

Err Price_VolBondCoupon(long today,
                        //	The underlying
                        char *yc,      //	yc
                        TermStruct *l, // underlying
                        //	SABR Parameters
                        double sabrBetavol, double sabrAlpha, double sabrBeta,
                        double sabrRho,
                        // Complementary Swap Sabr parameters
                        double swapRho, double ComplementarySabrbeta,
                        double ComplementarySabrbetaVol,
                        // Numerical Parameter
                        int NGaussQuadrature, double *x, double *w,
                        volBondCoupon *volBondCpn, double *VolBondCouponPrice);

Err Price_VolBond(long today,
                  //	The underlying
                  char *yc,      //	yc
                  TermStruct *l, // underlying
                  //	SABR Parameters
                  double *sabrBetaVol, double *sabrAlpha, double *sabrBeta,
                  double *sabrRho,
                  // Complementary Swap Sabr parameters
                  double *swapRho, double *CompSabrBeta,
                  double *CompSabrBetaVol,
                  // Numerical Parameter
                  int NGaussQuadrature, volBondStruct *volBond,
                  double *VolBondCpnPrices, double *VolBondPrice);

Err VolBond_Price(
    long today, char *yc, TermStruct *l,
    // VolBond Parameters//
    long NbCoupons, double *notional, char **first_SwapFreq,
    char **first_SwapBasis, char **first_SwapRefRate, char **second_SwapFreq,
    char **second_SwapBasis, char **second_SwapRefRate, long *first_fixingDate,
    long *first_startDate, long *first_endDate, long *second_fixingDate,
    long *second_startDate, long *second_endDate, long *payment_Date,
    double *alphaC, double *alphaP, double *betaC, double *betaP,
    double *strikeC, double *strikeP,
    // Main Swap Sabr Parameters//
    double *sabrBetaVol, double *sabrAlpha, double *sabrBeta, double *sabrRho,
    // Complementary Swap Sabr Parameters//
    double *swapRho, double *CompSabrBeta, double *CompSabrBetaVol,
    // Numerical Parameter//
    int NGaussQuadrature, double *cpnPrices, double *priceValue);

Err VolBond_PriceAutoCal(
    long today, char *yc, char *vc, char *refrate, double defaultLambda,
    double lgm2f_alpha, double lgm2f_gamma, double lgm2f_rho,
    // VolBond Parameters//
    long NbCoupons, double *notional, char **first_SwapFreq,
    char **first_SwapBasis, char **first_SwapRefRate, char **second_SwapFreq,
    char **second_SwapBasis, char **second_SwapRefRate, long *first_fixingDate,
    long *first_startDate, long *first_endDate, long *second_fixingDate,
    long *second_startDate, long *second_endDate, long *payment_Date,
    double *alphaC, double *alphaP, double *betaC, double *betaP,
    double *strikeC, double *strikeP,
    // Main Swap Sabr Parameters//
    double *sabrBetaVol, double *sabrAlpha, double *sabrBeta, double *sabrRho,
    // Complementary Swap Sabr Parameters//
    double *swapRho, double *CompSabrBeta, double *CompSabrBetaVol,
    // Numerical Parameter//
    int NGaussQuadrature, double *lambdaVec, double *sigmaVec,
    double *cpnPrices, double *priceValue);

#endif //__VOLBOND_H