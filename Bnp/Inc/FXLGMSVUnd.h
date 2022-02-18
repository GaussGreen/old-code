#ifndef FXLGMSVUND
#define FXLGMSVUND

#include "LGMSVUtil.h"
#include "math.h"
#include "srt_h_all.h"

void LGMSVSolveODE(
    double  FreqRe,
    double  FreqIm,
    double  Fwd,
    int     endi,
    double* dt,
    double* CoefInt,
    double* Coef1ReT,
    double* Coef1ImT,
    double* Coef2ReT,
    double* Coef2ImT,
    double* Coef3ReT,
    double* Coef3ImT,

    /* Outputs */
    double* LogTFRe,
    double* LogTFIm);

void HestonClosedFormApprox(
    /* Parameters of diffusion */
    int     iNbPWTime, /* Piece Wise Term Structures  */
    double* dPWTime,
    double* dSigma,
    double* dAlpha,
    double* dLevelEps,
    double* dLambdaEps,
    double* dRho,

    /* Product description */
    double dFwd,
    double dStrike,
    double dExTime,

    /* Numerical parameters */
    int    iNbX,
    double dIntegParam,

    /* Output */
    double* Price);

Err LGMSVBondVolApprox(
    // Parameters of diffusion
    int     iNbPWTime,  // Piece Wise Constant Term Structures
    double* dPWTime,
    double* dSigma,
    double* dAlpha,
    double* dLevelEps,
    double* dLambdaEps,
    double* dRho,
    double* dRho2,

    double dLambda,
    double dTStar,

    double* dAlphaLGM,
    double* dRhoLGM,
    double  dGammaLGM,

    // Product description
    int    sign,
    double dFixTime,
    double dBondMat,

    /* Numerical parameters */
    int    iNbX,
    double dIntegParam,

    // Outputs
    double* dBondATMVol);

typedef struct
{
    /* Today */
    long today;

    /* TS Dates and Times */
    int     num_dates;
    long*   dates;
    double* times;

    /* Domestic Underlying */
    irm_sv* dom_irmsv;

    /* Foreign Underlying */
    irm_sv* for_irmsv;

    /* FX Underlying */
    double fx_spot;

    double* fx_sigma;

    /* Correlations */
    double*** correlation;  // correlation[TimeIndex][VariableIndex1][VariableIndex2]
                            // 0 => Dom 1
                            // 1 => Dom 2
                            // 2 => Dom SV
                            // 3 => For 1
                            // 4 => For 2
                            // 5 => For SV
                            // 6 => FX

} fxlgmsv_str, *FXLGMSV_STR;

Err fxlgmsv_free_und_struct(SrtUndPtr pUndDesc);

Err fxlgmsv_get_struct_from_und(char* und, fxlgmsv_str** fxlgmsv);

Err fxlgmsv_get_fx_and_correl_ts(
    char*      und,
    double*    fx_spot,
    int*       num_time,
    double**   fx_time,
    double**   fx_vol,
    int*       num_rho,
    double**   rho_time,
    double**** correlation);

Err SrtInitFXLGMSVUnd(
    char* undName, /* und name */

    char* dom_undName, /* domestic underlying name */
    char* for_undName, /* foreign underlying name */

    /*	FX Underlying	*/
    double  fx_spot,
    int     fx_n_sigma,
    long*   fx_sigma_date,
    double* fx_sigma,

    /*	Correlations	*/
    long*     rho_date,
    int       rho_n,
    double*** correlation);

Err fxlgmsv_get_ts_from_und(
    char* undName, /* fxlgmsv und name */

    /* TS Dates and Times */
    int*     num_dates,
    long**   dates,
    double** times,

    /* Domestic Underlying */
    int* dom_one2F,

    double** dom_sigma,

    double** dom_tau,

    double* dom_alpha,
    double* dom_gamma,
    double* dom_rho,

    double** dom_alphaSV,
    double** dom_lamSV,
    double** dom_rhoSV,
    double** dom_rho2SV,

    double* dom_tstar,

    /* Foreign Underlying */
    int* for_one2F,

    double** for_sigma,

    double** for_tau,

    double* for_alpha,
    double* for_gamma,
    double* for_rho,

    double** for_alphaSV,
    double** for_lamSV,
    double** for_rhoSV,
    double** for_rho2SV,

    double* for_tstar,

    /* FX Underlying */
    double* fx_spot,

    double** fx_sigma,

    /* Correlations */
    double**** correlation);

typedef struct
{
    /* TS Dates and Times */
    int     num_times;
    double* times;

    /* Domestic Underlying */
    LGMSV_model* dom_lgmsv_model;

    /* Foreign Underlying */
    LGMSV_model* for_lgmsv_model;

    /* FX Underlying */
    double  fx_spot;
    double* fx_sigma;

    /* Correlations */
    double*** correlation;  // correlation[TimeIndex][VariableIndex1][VariableIndex2]
                            // 0 => Dom 1
                            // 1 => Dom 2
                            // 2 => Dom SV
                            // 3 => For 1
                            // 4 => For 2
                            // 5 => For SV
                            // 6 => FX

} fxlgmsv_model, *FXLGMSV_MODEL;

void fxlgmsv_init_NULL_model(fxlgmsv_model* model);

Err fxlgmsv_fill_model_from_struct(fxlgmsv_str* str, fxlgmsv_model* model);

Err fxlgmsv_get_model_from_und(char* undName, fxlgmsv_model** fxlgmsv);

Err qtolgmsv_check_und(char* und, int* for_one2F);

Err FXLGMSV_FXImpliedVolApprox(
    // Underlying Name
    char* UndName,

    // Product description
    double dFixTime,
    double dBondMat,

    // Numerical parameters
    int    iNbX,
    double dIntegParam,

    // Outputs
    double* FXImpliedVol);

Err FXLGMSV_FXImpliedVolPert(
    // Underlying Name
    char* UndName,

    // Product description
    double dFixTime,
    double dBondMat,

    // Numerical parameters
    int NHermiteQuad,
    int NLegendreQuad,

    // Outputs
    double* FXImpliedVol);

#endif
