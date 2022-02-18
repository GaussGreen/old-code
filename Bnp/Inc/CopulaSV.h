#ifndef __COPULASV_H
#define __COPULASV_H

#include "srt_h_types.h"

typedef struct
{
    long   lNumPaths;
    double dMinTime;
    int    iMCType;

    long         lNbPoints;
    int          iStudDegree;
    int          iNConv;
    SrtMCSamType eMCType;

    int    iAdaptGrid;
    long   lAdaptNbPoints;
    int    iUseOldGridMethod;
    double dNbStd;
    double dMaxError;

    int iCalcImpliedRho;

    int    iForceLogBeta;
    int    iUseSABR;
    double dPi;
    double dCalibBetaStd;

    int      iHasSavedGaussian;
    int      iFreeSavedGaussian;
    int      iNbFwd;
    double** dSavedGaussian;

    int iIntegSkipFirst;

    char cFileName[140];

} COPULASV_Params, *COPULASV_PARAMS;

void copula_sv_set_default_num_params(COPULASV_PARAMS sParams);

Err copula_sv_single_pricer(/* Product Parameters */
                            double  dTime,
                            int     iNbStrike,
                            double* dStrikes,

                            /* Model Parameters */
                            double dForward,
                            double dShift,
                            double dSigma,
                            double dAlpha,
                            double dRho,
                            double dVar0,
                            double dVarInf,
                            double dLambda,

                            /* Numerical Params */
                            COPULASV_PARAMS sParams,

                            /* Outputs */
                            double* dValue);

Err copula_sv_numer(/* Marginales Distributions */
                    int      iNbFwd,
                    double** xa,
                    double** ya,
                    long*    lNbPoints,

                    /* Copula Parameters */
                    double   dTime,
                    double*  dForward,
                    double*  dShift,
                    double*  dSigma,
                    double*  dAlpha,
                    double*  dRho,
                    double*  dVar0,
                    double*  dVarInf,
                    double*  dLambda,
                    double** dCorrMatrix,

                    /* Numerical Parameters */
                    COPULASV_PARAMS sParams,

                    /* Result */
                    double** dResult);

Err copula_sv_from_sabr_to_heston(
    double           dMaturity,
    double           dFwd,
    double           dSigma,
    double           dAlpha,
    double           dBeta,
    double           dRho,
    SrtDiffusionType eTypeInput,

    double* dSigmaH,
    double* dAlphaH,
    double* dShiftH,
    double* dRhoH,
    double* dVar0H,
    double* dVarInfH,
    double* dLambdaH);

/* Pricing of max(dWeights[0->n-1] * dFwds[0->n-1] - K, 0) */
Err copula_sv_basket_SABR(
    int              nFwds,
    double*          dFwds,
    double*          dWeights,
    int              iNbStrikes,
    double*          dStrike,
    double*          dVols,
    double*          dAlpha,
    double*          dBeta,
    double*          dRho,
    double**         dCorrelation,
    double           dMaturity,
    SrtCallPutType   eCallPut,
    SrtDiffusionType eTypeInput,
    COPULASV_PARAMS  sParams,
    double*          dPremium);

/* Pricing of (dWeigths[0]*dFwds[0] + dMargin) * max(dWeights[1->n-1] * dFwds[1->n-1] - K, 0) */
Err copula_sv_GearedOption_SABR(
    int              nFwds,
    double*          dFwds,
    double*          dWeights,
    double           dMargin,
    int              iNbStrikes,
    double*          dStrike,
    double*          dVols,
    double*          dAlpha,
    double*          dBeta,
    double*          dRho,
    double**         dCorrelation,
    double           dMaturity,
    SrtCallPutType   eCallPut,
    SrtDiffusionType eTypeInput,
    COPULASV_PARAMS  sParams,
    double*          dPremium);

/* Monte Carlo to check results */
Err copula_sv_basket_SABR_MC(
    int              nFwds,
    double*          dFwds,
    double*          dWeights,
    int              iNbStrikes,
    double*          dStrike,
    double*          dVols,
    double*          dAlpha,
    double*          dBeta,
    double*          dRho,
    double**         dCorrelation,
    double           dMaturity,
    SrtCallPutType   eCallPut,
    SrtDiffusionType eTypeInput,
    COPULASV_PARAMS  sParams,
    double*          dPremium);

Err copula_sv_quadrant_dependence_SABR(
    int              nFwds,
    double*          dFwds,
    double*          dWeights,
    int              iNbFract,
    double*          dFract,
    double*          dVols,
    double*          dAlpha,
    double*          dBeta,
    double*          dRho,
    double**         dCorrelation,
    double           dMaturity,
    SrtDiffusionType eTypeInput,
    COPULASV_PARAMS  sParams,
    double*          dResult);

Err copula_student_quadrant_dependence(
    int             iNbFwd,
    double*         dWeights,
    int             iNbFract,
    double*         dFract,
    double**        dCorrMatrix,
    COPULASV_PARAMS sParams,
    double*         dResult);

Err copula_sv_get_implied_linear_correl(
    int      iNbFwd,
    int      iIndex1,
    int      iIndex2,
    double*  dCoef1,
    double*  dCoef2,
    double*  dCross,
    double** dCorrelation,
    double*  dRho);

#endif
