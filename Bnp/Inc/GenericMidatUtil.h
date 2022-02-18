#ifndef GENERICMIDATUTIL
#define GENERICMIDATUTIL

#include "srt_h_all.h"

typedef struct
{
    int    iNbFactor;
    double dNumeraire;

    double dAlpha;
    double dGamma;
    double dRho;

    double dAlpha2;
    double dRhoAlpha;

    int     iNbVols;
    double* dTimes;
    double* dForward;
    double* dBeta;
    double* dLambda;
    double* dSigma;

    /* 2 Factor Specific */
    double dStartBeta;
    double dStartBeta2;

    double* dLambda2;
    double* dBeta2;
    double* dSigma2;
    double* dCorrel;

} genmidat_model, *GENMIDAT_MODEL;

Err genmidat_alloc_model(int iNbExe, int iNbFactor, GENMIDAT_MODEL sModel);

void genmidat_init_NULL_model(GENMIDAT_MODEL sModel);

Err genmidat_init_model(
    double         dAlpha,
    double         dGamma,
    double         dRho,
    double*        dTimes,
    double*        dForward,
    double*        dSigma,
    double*        dBeta,
    double*        dSigma2,
    double*        dCorrel,
    double         dStartBeta2,
    double         dNumeraire,
    GENMIDAT_MODEL sModel);

void genmidat_free_model(GENMIDAT_MODEL sModel);

Err genmidat_copy_model(GENMIDAT_MODEL sSourceModel, GENMIDAT_MODEL sDestModel);

void genmidat_convert_lambda_into_beta(GENMIDAT_MODEL sModel);
void genmidat_convert_beta_into_lambda(GENMIDAT_MODEL sModel);

Err genmidat_fill_2factor(GENMIDAT_MODEL sModel);

double genmidat_Swap_correlation(
    double dTime1, double dTime2, double dAlpha, double dGamma, double dRho);

Err genmidat_free_und_struct(SrtUndPtr pUndDesc);

Err SrtInitGenericMidatUnd(
    char*   undName, /* und name */
    int     iNbFactor,
    double  dNumeraire,
    double  dAlpha,
    double  dGamma,
    double  dRho,
    int     iNbVols,
    double* dTimes,
    double* dForward,
    double* dSigma,
    double* dBeta,
    double* dSigma2,
    double* dCorrel,
    double  dStartBeta2);

Err genmidat_get_model(char* cUndName, GENMIDAT_MODEL sModel);

Err genmidat_implied_volatility_and_price_model(
    int            iIntegStartIndex,
    int            iIntegEndIndex,
    int            iUndIndex,
    int            iNbPeriod,
    GENMIDAT_MODEL sModel,
    double*        dVolatility,
    double*        dPrice);

Err genmidat_implied_volatility_and_price_und(
    char*   cUndName,
    int     iIntegStartIndex,
    int     iIntegEndIndex,
    int     iUndIndex,
    int     iNbPeriod,
    double* dVolatility,
    double* dPrice);

Err genmidat_implied_correlation_model(
    int            iIntegStartIndex,
    int            iIntegEndIndex,
    int            iUndIndex1,
    int            iNbPeriod1,
    int            iUndIndex2,
    int            iNbPeriod2,
    GENMIDAT_MODEL sModel,
    double*        dCorrelation);

Err genmidat_implied_correlation_und(
    char*   cUndName,
    int     iIntegStartIndex,
    int     iIntegEndIndex,
    int     iUndIndex1,
    int     iNbPeriod1,
    int     iUndIndex2,
    int     iNbPeriod2,
    double* dCorrelation);

#endif