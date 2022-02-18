#ifndef __COPULAGAUSSIAN_H
#define __COPULAGAUSSIAN_H

#include "srt_h_types.h"

typedef struct
{
    long         lNumPaths;
    long         lNbPoints;
    SrtMCSamType eMCType;

    int    iInterpMethod;
    int    iAdaptGrid;
    long   lAdaptNbPoints;
    int    iUseOldGridMethod;
    double dNbStd;
    double dMaxError;

    /* All the numerical params */
    int      iHasSavedGaussian;
    int      iFreeSavedGaussian;
    int      iNbFwd;
    double** dSavedGaussian;

    int    iForceLogBeta;
    int    iUseSABR;
    double dPi;
    double dCalibBetaStd;

    int iIntegSkipFirst;

    char cFileName[140];

} COPULAGAUSSIAN_Params, *COPULAGAUSSIAN_PARAMS;

//void copula_gaussian_set_default_num_params(COPULAGAUSSIAN_PARAMS sParams){};
//
//Err copula_gaussian_init_gaussian(int iNbFwd, COPULAGAUSSIAN_PARAMS sParams)
//{
//    return "";
//};
//
//void copula_gaussian_free_gaussian(COPULAGAUSSIAN_PARAMS sParams){};
//
///* Compute the Cumulative of a SABR distribution !!! */
//Err copula_gaussian_get_SABR_cumulative(
//    double dMaturity,
//    double dForward,
//    double dSigmaBeta,
//    double dAlpha,
//    double dBeta,
//    double dRho,
//
//    int  iAdaptGrid,
//    long lAdaptPoints,
//    int  iUsePoints,
//
//    long    lNumPoints,
//    double* dPoints,
//    double* dCumulative)
//{
//    return "";
//};
//
///* Compute the Cumulative of a SABR distribution !!!	*/
///* Recurcive version with control of interpolation		*/
//Err copula_gaussian_get_BMM_linterp_cumulative(
//    double dMaturity,
//    double dForward,
//    double dForward1,
//    double dSigmaBeta1,
//    double dForward2,
//    double dSigmaBeta2,
//    double dBeta,
//    double dPi,
//
//    double dMaxError,
//    double dNbStdLeft,
//    double dNbStdRight,
//
//    long*    lNumPoints,
//    double** dPoints,
//    double** dCumulative)
//{
//    return "";
//};
//
///* Pricing of max(dWeights[0->n-1] * dFwds[0->n-1] - K, 0) */
///* ******************************************************* */
//
//Err copula_gaussian_basket_SABR(/* Marginales Distributions */
//                                double           dMaturity,
//                                int              iNbFwd,
//                                double*          dFwds,
//                                double*          dVols,
//                                double*          dAlpha,
//                                double*          dBeta,
//                                double*          dRho,
//                                SrtDiffusionType eTypeInput,
//
//                                /* Payoff Definition */
//                                double*        dWeights,
//                                int            iNbStrikes,
//                                double*        dStrike,
//                                SrtCallPutType eCallPut,
//
//                                /* Model Parameters */
//                                double**              dCorrelation,
//                                COPULAGAUSSIAN_PARAMS sParams,
//
//                                /* Results */
//                                double* dPremium)
//{
//    return "";
//};
//
//Err copula_gaussian_numer(/* Marginales Distributions */
//                          int      iNbFwd,
//                          double** dX,
//                          double** dCumulative,
//                          long*    lNbPoints,
//
//                          /* Copula Parameters */
//                          int      iIsChoMatrix,
//                          double** dCorrMatrix,
//
//                          /* Numerical Parameters */
//                          COPULAGAUSSIAN_PARAMS sParams,
//
//                          /* Result */
//                          double** dResult)
//{
//    return "";
//};
//
//Err copula_gaussian_GearedOption_SABR(/* Marginales Distributions */
//                                      double           dFloatMaturity,
//                                      double           dCMSMaturity,
//                                      int              iNbFwd,
//                                      double*          dFwds,
//                                      double*          dVols,
//                                      double*          dAlpha,
//                                      double*          dBeta,
//                                      double*          dRho,
//                                      SrtDiffusionType eTypeInput,
//
//                                      /* Payoff Definition */
//                                      double*        dWeights,
//                                      double         dMargin,
//                                      int            iNbStrike,
//                                      double*        dStrike,
//                                      SrtCallPutType eCallPut,
//
//                                      /* Model Parameters */
//                                      double**              dCorrMatrix,
//                                      COPULAGAUSSIAN_PARAMS sParams,
//
//                                      /* Results */
//                                      double* dPremium)
//{
//    return "";
//};

#endif
