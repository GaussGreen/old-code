///* ========================================================================== */
//
#ifndef BGMUTILS_H
#define BGMUTILS_H

#include "BGMTypes.h"
#include "uterror.h"
#include "uttypes.h"
//
//Err interp_STRIKETYPE(char* str, STRIKE_TYPE* val);
//
//Err interp_SMILEMOVE(char* str, SMILE_MOVE_TYPE* val);
//
//Err interp_DISPLAYTYPE(char* str, DISPLAY_TYPE* val);
//
//Err srt_BGMSetUpSchedule(
//    long            MaxNumPeriod,
//    char*           szRefRate,
//    char*           szYieldCurveName,
//    long*           tFixing,
//    long*           tStart,
//    long*           tPay,
//    double*         tCvg,
//    double*         tExp,
//    SrtCompounding* freq,
//    SrtBasisCode*   bas);
//
//Err srt_BGMSetUpScheduleFixLeg(
//    long            MaxFixPeriod, /*The number of periods for the fixed leg schedule */
//    char*           szYieldCurveName,
//    long*           tFixing,
//    long*           tStart,
//    long*           tPay,
//    double*         tCvg,
//    double*         tExp,
//    SrtCompounding* freq,
//    SrtBasisCode*   bas);
//
//Err srt_BGM_FwdRate(
//    long    start,
//    long    theoend,
//    char*   szRefRate,
//    char*   szYieldCurveName,
//    int     CashOrLibor, /* 0 Cash or 1 Libor */
//    double* FwdRate);
//
//Err srt_BGM_FwdSwapRate(
//    long    start,
//    long    theoend,
//    char*   szRefRate,
//    char*   szYieldCurveName,
//    char*   szFrequency,
//    char*   szBasis,
//    int     CashOrLibor, /* 0 cash, 1 Libor */
//    double* FwdRate);
//
//Err srt_BGM_Level(long start, long theoend, char* szRefRate, char* szYieldCurveName, double* Level);
//
//Err srt_f_BGMGetVol(
//    long        today, /* only use in case of striketype = STDEV or Cash = 0 */
//    long        tStart,
//    long        tEnd,
//    double      Fwd, /* Libor Forward */
//    double      Spread,
//    int         CashOrLibor,
//    char*       szVolCurveName,
//    double      strike, /* Assuming Libor Strike */
//    STRIKE_TYPE striketype,
//    double*     Vol);
//
//Err srt_f_BGMGetSABRVol(
//    long             today, /* only use in case of striketype = STDEV or Cash = 0 */
//    long             tStart,
//    long             tEnd,
//    double           Fwd, /* Libor Forward */
//    double           Spread,
//    int              CashOrLibor,
//    char*            szVolCurveName,
//    double           strike, /* Assuming Libor Strike */
//    STRIKE_TYPE      StrikeOrStdev,
//    SABRVolComponent ResultType,
//    double*          Result);
//
//Err srt_BGMFillLiborTS(
//    long MaxNumPeriod, double** ppdMarketVols, SMILE_MOVE_TYPE ModelSmileMove, double** ppdLiborTS);
//
//Err srt_f_BGMFromBetaVoltoShiftedVol(
//    double ATMVol, double Forward, double Beta, double maturity, double* ShiftedVol, double* Shift);
//
//Err srt_f_BGMFromShiftedVoltoBetaVol(
//    double  Forward,
//    double  Maturity,
//    double  ShiftedVol,
//    double  Shift,
//    double* Beta,
//    double* BetaVol);
//
//Err srt_f_matching_moments(double* shift, double* varshift, double M1, double M2, double M3);
//
//Err srt_f_sumshiftedlognormal(
//    double*  Forwards,
//    double*  Shifts,
//    double*  VolShifts,
//    double*  Weights,
//    double** Correl,
//    int      Num_assets,
//    double   Maturity,
//    double*  Forward_bskt,
//    double*  Shift_bskt,
//    double*  VolShift_bskt);
//
//Err srt_f_BGMAdjustVolandGetShifts(
//    long      MaxNumPeriod,
//    long      MaturityInPeriod,
//    double**  ppdMarketVols,
//    double*   pdBeta,
//    double*** pppdWeights,
//    double*   tExp,
//    double**  tFwdSwap,
//    double**  tFwdSwapSpread,
//    double**  tFwdSwapLibor,
//    double**  tFwdSwapShift,
//    double*   tLiborShift,
//    double**  Strikes,
//    int       CashOrLibor,
//    double*   pdshift);
//Err srt_f_BGM_GetStrikes(
//    long        MaxNumPeriod,
//    long        MaturityInPeriod,
//    long*       tStart,
//    long*       tPay,
//    double**    tFwdSwap,
//    char*       szVolCurveName,
//    double      CapStrike,
//    STRIKE_TYPE CapStrikeType,
//    double      SwapStrike,
//    STRIKE_TYPE SwapStrikeType,
//    double**    Strikes);
//
//Err srt_f_BGMSABRtoFudge(
//    double  forward,
//    double  maturity,
//    double* alphaSABR,
//    double* betaSABR,
//    double* rhoSABR,
//    double* betavolSABR,
//    double* shift,
//    double* volup,
//    double* voldown,
//    int     way /* 0 for from SABR to Fudge, 1 for from Fudge to SABR */
//);
//
//Err srt_f_BGMFillsCapletforSABRRisk(
//    int     iStartCaplet,
//    long    iStartSwaption,
//    long    IEndSwaption,
//    double* tLibor,
//    double* tLiborVols,
//    double* tExp,
//    double* shifts,
//    double* shiftVols,
//    double  alphas,
//    double  betas,
//    double  rhos,
//    double* betavols,
//    double* smileAlphas,
//    double* smileBetas,
//    double* smilebetaVols,
//    double* SV_shifts,
//    double* up_shiftVols,
//    double* down_shiftVols);
//
//Err srt_f_BGMInterpolSABR(
//    double   undCaplet,
//    char*    szRefRate,
//    char*    szYieldCurveName,
//    char*    szVolCurveName,
//    double*  matLiquidOptions,
//    double*  undLiquidOptions,
//    int      n_LiquidOptions,
//    double*  liquidalphas,
//    double*  liquidbetas,
//    double*  liquidrhos,
//    double** Correl,
//    int      dimCorrel,
//    double*  matOutputs,
//    double*  undOutputs,
//    int      n_matOutputs,
//    int      n_undOutputs,
//    double** atmvolOutputs,
//    double** alphaOutputs,
//    double** betaOutputs,
//    double** rhoOutputs,
//    double   LastMaturity);
//
//Err srt_f_BGMSensitivitySABR(
//    double   undCaplet,
//    char*    szRefRate,
//    char*    szYieldCurveName,
//    char*    szVolCurveName,
//    double*  inputalphas,
//    double*  inputbetas,
//    double*  inputrhos,
//    double** Correl,
//    int      dimCorrel,
//    double   matSwaption,
//    double   undSwaption,
//    double   shiftatmvol,
//    double   shiftalpha,
//    double   shiftbeta,
//    double   shiftrho,
//    double** greeksOutputs);
//
//Err srt_bgm_FudgeForDecreasingCorrelation(
//    double** Correl, int NbrMaturities, double** DecreasingCorrel);
//
//Err srt_bgm_SymmetricToPositive(double** OldCorrel, double** NewCorrel, int NewNbrMaturities);
//
//Err srt_bgm_InterpolExtrapolCorrel(
//    double** OldCorrel,
//    double** NewCorrel,
//    double*  OldMaturities,
//    double*  NewMaturities,
//    int      NbrOldMaturities,
//    int      NewNbrMaturities);
//
//Err srt_bgm_InterpolCorrelation(
//    int      OldNbrMaturities,
//    int      NewNbrMaturities,
//    double*  OldMaturities,
//    double*  NewMaturities,
//    double** OldCorrel,
//    double** NewCorrel);
//
//Err srt_bgm_GetShiftedLogModelFrom2Strikes(
//    double  forward,
//    double  maturity,
//    double* strikes,
//    double* marketvols,
//    double* shift,
//    double* volshift,
//    double* error);
//
//Err srt_BGMNewRhoSABR(
//    double  forward,
//    double  maturity,
//    double  ATMVol,
//    double  alpha,
//    double  beta,
//    double  rho,
//    double  newbeta,
//    double* Newrho);
#endif