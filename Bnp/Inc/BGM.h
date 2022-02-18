

#include <BGMUTILS.H>
#include <NUM_H_ALLHDR.H>

#include "utallhdr.h"

/* ========================================================================== */

#ifndef BGM_H
#define BGM_H

Err shiftmatrix1to0(int dim, double** matrix, double** matrixshifted);

Err shiftmatrix0to1(int dim, double** matrixshifted, double** matrix);

Err get_bgm_fra(
    long Fradate, int spotlag, char* cRefRateCode, char* szYieldCurveName, double* dFra);

Err get_bgm_SwapFraWeights(
    double  maturity,
    double  underlying,
    char*   szYieldCurveName,
    char*   cRefRateCode,
    double* weights);

Err get_bgm_SwapFraHatWeights(
    double  maturity,
    double  underlying,
    char*   szYieldCurveName,
    char*   cRefRateCode,
    double* weights,
    double* hatweights);

Err srt_BGMF_Weight(
    long      MaxNumPeriod,
    long      MaturityInPeriod,
    char*     szRefRate,
    char*     szYieldCurveName,
    long      today,
    long*     tStart,
    long*     tPay,
    double*   tCvg,
    double*** pppdWeights);

Err srt_BGMF_WeightNew(
    long      MaxNumPeriod,
    long      MaturityInPeriod,
    char*     szRefRate,
    char*     szYieldCurveName,
    char*     szFrequency,
    char*     szBasis,
    long      today,
    long*     tStart,
    long*     tPay,
    double*   tCvg,
    double*** pppdWeights);

Err build_instrument(double* vecs, double** matres, int dim);

Err get_bgm_AllCoverages(
    double** instruments,
    int      num_of_fras,
    char*    szYieldCurveName,
    char*    cRefRateCode,
    double*  cvgs);

Err get_bgm_SwaptionOMEGAMatrix(
    double* hatweights, int num_of_tenors, double** OmegaMatrix, double* coverages);

Err get_bgm_CapletOMEGAMatrix(
    int CapletIndex, int num_of_tenors, double** OmegaMatrix, double* coverages);

Err srt_bgm_OMEGA_Constraints(
    double**  instruments,
    int       num_of_instruments,
    int       num_of_caplets,
    double*** OMEGA,
    char*     szYieldCurveName,
    char*     cRefRateCode);

Err get_bgm_SwapNumOfPeriods(
    double maturity,
    double underlying,
    char*  szYieldCurveName,
    char*  cRefRateCode,
    int*   VectorSize);

/* Get the Market Vols for a given set of instruments : the instruments are described by : a
 * maturity in number of year, an underlying in number of year */
Err get_bgm_MarketVols(
    double** instruments,
    int      num_of_instruments,
    double   strike,
    double*  VolsMark,
    char*    cRefRateCode,
    char*    szYieldCurveName,
    char*    szVolCurveName);

Err srt_sdp(
    double*** const_mat,
    double*   const_val,
    double**  cmat,
    int       dim,
    int       num_const,
    double**  xmat,
    double*   obj_val,
    double*   error,
    int*      return_error,
    double    toler,
    int       printlevel,
    int       max_iterations);

Err get_bgm_HistoricalVarianceCovarianceMatrix(
    double** capletsinstruments,
    int      num_of_fras,
    double   strike,
    double** historicalcorrelation,
    double** histvarcovar,
    char*    szYieldCurveName,
    char*    cRefRateCode,
    char*    szVolCurveName);

Err srt_bgm_Calibration(
    double** instruments,
    int      num_of_instruments,
    double** capletsinstruments,
    int      num_of_fras,
    double   strike,
    double** histcorrelation,
    char*    cRefRateCode,
    char*    szYieldCurveName,
    char*    szVolCurveName,
    double** calibratedmatrix,
    double*  result,
    double*  error,
    int*     return_err,
    double   toler,
    int      printlevel,
    int      niter);

Err srt_bgm_Calibration_with_PCA(
    double** instruments,
    int      num_of_instruments,
    double** capletsinstruments,
    int      num_of_fras,
    double   strike,
    double** histcorrelation,
    char*    cRefRateCode,
    char*    szYieldCurveName,
    char*    szVolCurveName,
    int*     num_of_factors,
    double** volfunc,
    double*  result,
    double*  error,
    int*     return_err,
    double   toler,
    int      printlevel,
    int      niter);

Err srt_bgm_MatrixReduction(double* eigen_val, int dim, int* numfactors, double** VolMatrix);

Err srt_bgm_VarCovarRebuild(double** VolMatrix, double*** NewVarCovar, int dimrow, int dimcol);

Err srt_bgm_MCEvolve(
    double** LiborMatrix,
    int      ndiscretization,
    int*     Indices,
    double*  maturity,
    double*  coverages,
    int LastTenor, /* Last tenor necessary to compute all the discount factors, numeraire of the
                      measure in which we price*/
    int      FirstTenor,
    double** Brownianincs,
    int      num_of_factors,
    double** VolMatrix,
    double** VarCovar);

Err srt_bgm_MCEvolveCasc(
    double** LiborMatrix,
    int      ndiscretization,
    int*     Indices,
    double*  maturity,
    double*  coverages,
    double*  shifts,
    int LastTenor,  /* LastTenor-1 is the Last Tenor necessary to compute all the discount factors,
                       numeraire of the measure in which we price*/
    int FirstTenor, /*FirstTenor for which we actually need to compute the Libor rate, typically,
                       the Tenor just before the first event date*/
    double**  Brownianincs,
    int*      factors,
    double*** VolMatrices,
    double*** VarCovarCasc);

Err srt_bgm_MCEvolveJumping(
    double** LiborMatrix,
    int      ndiscretization,
    int*     Indices,
    int*     JumpingIndices,
    double*  maturity,
    double*  coverages,
    int LastTenor,  /* LastTenor-1 is the Last Tenor necessary to compute all the discount factors,
                       numeraire of the measure in which we price*/
    int FirstTenor, /*FirstTenor for which we actually need to compute the Libor rate, typically,
                       the Tenor just before the first event date*/
    double** Brownianincs,
    int      num_of_factors,
    double** VolMatrix,
    double** VarCovar);

Err srt_bgm_MCEvolveJumpingCasc(
    double** LiborMatrix,
    int      ndiscretization,
    int*     Indices,
    int*     JumpingIndices,
    double*  maturity,
    double*  coverages,
    int LastTenor,  /* LastTenor-1 is the Last Tenor necessary to compute all the discount factors,
                       numeraire of the measure in which we price*/
    int FirstTenor, /*FirstTenor for which we actually need to compute the Libor rate, typically,
                       the Tenor just before the first event date*/
    double**  Brownianincs,
    int*      factors,
    double*** VolMatrices,
    double*** VarCovarCasc);

Err srt_bgm_ComputeRollingNumeraire(
    double** LiborMatrix,
    int*     TenorIndicesinDisc,
    int      LastTenorIndex,
    double*  RollingNumeraire,
    double*  TenorCoverages);

Err srt_bgm_ComputeDiscountFactors(
    double* Libors,
    long*   dfdates, /*discount factor dates for this event date augmented by the LastTenor date*/
    int     numdfdates,
    double* TenorCoverages,
    double* EndCoverages,
    double  StartCoverage,
    int     StartIndex,
    int*    EndIndices,
    double* dfs);

Err srt_f_BGMVolMatrix(
    long     MaxNumPeriod,
    long     MaturityInPeriod,
    long     nCorrelRow,
    double** ppdHistCorrel,
    double*  pdBeta,
    char*    szRefRate,
    char*    szYieldCurveName,
    char*    szVolCurveName,
    char*    szFrequency,
    char*    szBasis,
    double   CapStrike,
    char*    szCapStrikeType,
    double   SwapStrike,
    char*    szSwaptionStrikeType,
    int      CashOrLibor, /* 0 cash, 1 Libor */
    char*    szCorrelMode,
    double** answer,
    double*  answershift);

Err srt_f_BGMImpliedCorrel(
    long     MaxNumPeriod,
    char*    szRefRate,
    char*    szYieldCurveName,
    char*    szVolCurveName,
    char*    szCorrelMode, /* SLIDING or CONVERGING */
    char*    szSmileMode,  /* SLIDING or CONVERGING */
    char*    szFrequency,
    char*    szBasis,
    double   CapStrike,
    char*    szCapStrikeType,
    double   SwapStrike,
    char*    szSwaptionStrikeType,
    int      CashOrLibor,      /* 0 cash, 1 Libor */
    int      InputLiborTSFlag, /* 0 noinput 1 LiborTS already defined */
    double** answerCorrel,
    double** answerLiborTS,
    double*  Frequency);

Err srt_f_FwdVolMatrix(
    long     nCorrelRow,
    double** histcorrel,
    long     nInstr,
    long**   Instr,
    long**   InstrDate,
    char*    szRefRate,
    long     NumPeriodEnd,
    char*    szYieldCurveName,
    char*    szVolCurveName,
    double** VolMatrix,
    long     nVolMatrixTime,
    long     nVolMatrixLibor,
    char*    szFrequency,
    char*    szBasis,
    long*    MaxNumPeriod,
    long*    FreqFloatRate,
    long*    FreqFixedRate,
    char* szDisplaytype, /* NULL or TS == term structure, CUM == cumulative vol upto expiry of the
                            option */
    double*** answer);

Err srt_f_BGMGetATMSwaptionVol(
    char*    BGMUnd,
    double** ppdCorrelInput,
    long     nCorrelSize,
    char*    szCorrelType,
    double** ppdVolMatrix,
    long     nVolInput,
    char*    szRefRateInput,
    char*    szYieldCurveNameInput,
    int      UndOrMatrix, /* 0 if BGMUnd is input, 1 if ppdCorrel and ppd VolMatrix are input */
    long*    Expiry,
    long     nExpiry,
    long*    Underlying,
    long     nUnderlying,
    char*    SwapFrequency,
    char*    SwapBasis,
    double*  Frequency, /* Output, the frequency of the Swaption */
    double** SwaptionATMVol);

Err srt_f_BGMGetATMSwaptionVolNew(
    double**         ppdCorrelInput,
    long             nCorrelSize,
    char*            szCorrelType,
    double**         ppdVolMatrix, /*This matrix is in Normal Vols!*/
    long             nVolInput,
    char*            szRefRateInput,
    char*            szYieldCurveNameInput,
    double*          Expiry,
    long             nExpiry,
    double*          Underlying,
    long             nUnderlying,
    double           VolMat, /* The Cumulative Vol up to that maturity is what is returned*/
    char*            SwapFrequency,
    char*            SwapBasis,
    SrtDiffusionType VolType,
    double*          Frequency, /* Output, the frequency of the Swaption */
    double**         SwaptionATMVol);

Err srt_bgm_GetVolMatrixSlidingSliding(
    double** VolMatrix,
    double** correl,
    int      num_tenors,
    int*     factorsofmatrix,
    double** VolMatrixPCA);

/* Needed for BGMStripATM.c */
Err srt_BGMFwdVolComputeTenors(
    long     MaxNumPeriod,
    long     MaturityInPeriod,
    long*    tStart,
    long*    tPay,
    char*    szRefRate,
    char*    szBasis,
    int      CashOrLibor, /* 0 Cash, 1 Libor */
    double*  tLibor,
    double** tFwdSwap,
    double** tFwdSwapSpread);

Err srt_BGMFwdVolComputeTenorsNew(
    long     MaxNumPeriod,
    long     MaturityInPeriod,
    long*    tStart,
    long*    tPay,
    char*    szRefRate,
    char*    szYieldCurveName,
    char*    szFrequency,
    char*    szBasis,
    int      CashOrLibor, /* 0 Cash, 1 Libor */
    double*  tLibor,
    double** tFwdSwap,
    double** tFwdSwapSpread);

/* Needed for BGMStripATM.c */
Err srt_BGMFwdVol_computeCoeff(
    long            indexTime,
    long            indexTenor,
    double*         tLibor,
    double*         tExp,
    double**        tFwdSwap,
    double***       pppdWeights,
    double**        ppdHistCorrel,
    SMILE_MOVE_TYPE ModelCorrelMove,
    double**        ppdTempVol,
    double*         A,
    double*         B,
    double*         C);

/* Needed for BGMStripATM.c */
Err BGM_FillMarketVol(
    long        MaxNumPeriod,
    long        MaturityInPeriod,
    long*       tStart,
    long*       tPay,
    double**    tFwdSwap,
    double**    tFwdSwapSpread,
    int         CashOrLibor,
    double      CapStrike,
    STRIKE_TYPE CapStrikeType,
    double      SwapStrike,
    STRIKE_TYPE SwapStrikeType,
    char*       szRefRate,
    char*       szVolCurveName,
    double**    ppdMarketVols);

/* Needed for BGMStripATM.c */
double srt_BGMFwdVol_computecovar(
    long            indexTime,
    long            indexStart,
    long            indexEnd,   /* For the Swap */
    long            indexTenor, /* on the libor indexStart up to indexTenor - 1*/
    double*         tLibor,
    double***       pppdWeights,
    double**        ppdHistCorrel,
    SMILE_MOVE_TYPE CorrelMode,
    double**        ppdTempVol);
#endif