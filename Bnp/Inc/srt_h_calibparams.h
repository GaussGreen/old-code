#ifndef SRT_H_CALIBPARAMS_H
#define SRT_H_CALIBPARAMS_H

/* ------------------------------------------------------------------------
   STRUCTURE: SrtCalibParam

   FIELDS:
                - eAlgoType: LEVENBERG_MARQUARDT  , SIMPLEX  ,
  SIMULATED_ANNEALING  , SOBENBERG Which algorithm to use in the calibration: LM
  , SimuAnn  , Simplex...
                - eCalibType: GLOBAL  , FIXED
                        Which correlation calibration to go for: fixed corr  ,
  global  , blind on all...
                - lNumIter: any strictly positive long
                        Number of succeeded iterations in the algorithm
                - lNumSobenPoints: any strictly positive long
                        Number of points in SOBENBERG
                - lNumBestSobenPoints: any strictly positive long
                        Number of best points in SOBENBERG
                - lNumClusters: any strictly positive long
                        Number of clusters in SOBENBERG
                - szSobenMethod: "BESTFIT" or "BARC"
                        Method for minimum research within the clusters
                - dSimplexInitScale: any strictly positive double
                        Initial scaling for simplex and annealing
                - dSimplexTol: any strictly positive double
                        Fractional tolerance for simplex and annealing
                - dMaxTemp: any strictly positive double
                        Initial temperature for annealing
                - dMinTemp: any strictly positive double
                        Final temperature for annealing
                - dDecreaseFactor: any strictly positive double
                        decreasing factor for annealing temperature
                - bOneTau: SRT_NO or SRT_YES
                        Use a TermStruct of tau (NO) or just one single value
  (YES)
                - bFreezeTau: SRT_NO or SRT_YES
                        Freeze the tau structure
                - lNumSobolPaths: 0 or any positive long
                        Number of paths used in the presampling with Sobol (0==
  do not use)
                - bLgmStartPoint:  SRT_YES or SRT_NO
                        Use LGM as a starting point for the Calibration
  (Cheyette)
                - lLgmNumIter:
                        Number of iterations used for the LGM starting point
  calibration
                - bSmoothSigma: add an extra criteria for smoothness of sigma
                - dSigmaMin: >=0.0
                        Minimum value authorised for sigma in the Sobol
  pre-sampling
                - dSigmaMax;
                        Maximum value authorised for sigma in the Sobol
  pre-sampling
                - dLambdaMin;
                        Minimum value authorised for lambda in the Sobol
  pre-sampling
                - dLambdaMax;
                        Maximum value authorised for lambda in the Sobol
  pre-sampling
                - dOptionsWeight;
                        Weight of the Option Prices matrix in the minimisation
  criteria
                - dBetaMin;
                        Minimum value authorised for beta
                - dBetaMax;
                        Maximum value authorised for beta
                - dOmegaMin;
                        Minimum value authorised for omega (= beta 2 - beta 1)
                - dOmegaMax;
                        Maximum value authorised for omega
                - dAlphaMin;
                        Minimum value authorised for alpha
                - dAlphaMax;
                        Maximum value authorised for alpha
                - dGammaMin;
                        Minimum value authorised for gamma
                - dGammaMax;
                        Maximum value authorised for gamma
                - dRhoMin;
                        Minimum value authorised for rho
                - dRhoMax;
                        Maximum value authorised for rho
                - bAggregate
                        Flag between old and new calibration
  ------------------------------------------------------------------------ */
typedef struct {
  SrtCalibAlgoType eAlgoType;
  SrtCalibType eCalibType;
  long lNumIter;
  long lNumSobenPoints;
  long lNumBestSobenPoints;
  long lNumClusters;
  char *szSobenMethod;
  double dSimplexInitScale;
  double dSimplexTol;
  double dMaxTemp;
  double dMinTemp;
  double dDecreaseFactor;
  SRT_Boolean bOneTau;
  SRT_Boolean bFreezeTau;
  SRT_Boolean bOneBeta;
  SRT_Boolean bFreezeBeta;
  SRT_Boolean bOneOmega;
  SRT_Boolean bFreezeOmega;
  long lNumSobolPaths;
  SRT_Boolean bLgmStartPoint;
  long lLgmNumIter;
  SrtCalibAlgoType eLgmAlgoType;
  SRT_Boolean bSmoothSigma;
  double dSmoothSigmaWeight;
  double dSigmaMin;
  double dSigmaMax;
  double dLambdaMin;
  double dLambdaMax;
  double dOptionsWeight;
  double dBetaMin;
  double dBetaMax;
  double dAlphaMin;
  double dAlphaMax;
  double dGammaMin;
  double dGammaMax;
  double dRhoMin;
  double dRhoMax;
  double dOmegaMin;
  double dOmegaMax;
#ifndef PVMPI
  double bAggregate;
#else
  SRT_Boolean bAggregate;
#endif
} SrtCalibParam;

Err srt_f_set_CalibParams(String *pszParamStrings, String *pszValueStrings,
                          int iNumParams, SrtCalibParam *psCalibParam,
                          SrtMdlType eModelType);

typedef struct {

  SrtCalibAlgoType eAlgoType;
  long lNumIter;
  SRT_Boolean bSmoothFXVol;
  double dSmoothFXVolWeight;
  SRT_Boolean bCalibCorrel;

} SrtFXCalibParam;

Err srt_f_set_FXCalibParams(String *pszCalibParamNames,
                            String *pszCalibParamValues, int iNumCalibParamss,
                            SrtFXCalibParam *psFXCalibParams);

#endif
