#ifndef __CMS_SPREAD_OPTIONS_H
#define __CMS_SPREAD_OPTIONS_H

#include "CopulaGaussian.h"

#define XL_MAX_ID_LEN 40
#define CMSS_MAXSL 3

typedef struct {
  long lToday;
  char cYcName[XL_MAX_ID_LEN];
  char cVolName[XL_MAX_ID_LEN];
  char cFixedFreq[XL_MAX_ID_LEN];
  char cFixedBasis[XL_MAX_ID_LEN];
  char cFloatBasis[XL_MAX_ID_LEN];

  Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                      char *vol_curve_name, double start_date, double end_date,
                      double cash_strike, int zero, char *ref_rate_name,
                      double *vol, double *power);

  long lNumStrikesInVol; /*	Array of strikes in vol matrix */
  double *dStrikesInVol;
  SrtDiffusionType eVolType; /*	Type of vol in matrix      , SRT_NORMAL or
                                SRT_LOGNORMAL */
  int iCashVol; /*	1: matrix is a cash vol 0: matrix is a swap vol */
  int iIsSABR;
  int iIsSABRAF;

} cms_spread_mkt, *CMS_SPREAD_MKT;

void free_cms_spread_mkt(cms_spread_mkt *market);

typedef enum CMSSpreadModelType_ {
  CMSS_NORMAL = 0,
  CMSS_LOGNORMAL,
  CMSS_SHIFTEDLOG,
  CMSS_SLUSER,
  CMSS_MIXEDSL,
  CMSS_COPULA,
  CMSS_COPULAMIXEDSL,
  CMSS_COPULASV

} CMSSpreadModelType;

Err interp_cms_spread_model(const char *constStr, CMSSpreadModelType *val);

typedef struct {
  /* For Option / Digital part */
  CMSSpreadModelType eModelChoice;
  double dCorrel;
  double dBeta1;
  double dBeta2;
  int iUseCMSSmile;
  int iUseCMSImpliedVol;
  int iUseCMSImpliedSmile;
  int iUseCMSTEC;

  double dDefaultBeta;
  double dNbStdBetaCalib;
  int iUseShiftOrBeta;

  int iNbCorrel;
  double *dCorrelMat;
  double *dCorrelMin;
  double *dCorrelMax;
  double dATMShift;

  /* For Float part */
  CMSSpreadModelType eFloatModelChoice;
  int iFloatNbCorrel;
  double *dFloatCorrelMat;
  double *dFloatCorrel1;
  double *dFloatCorrel2;
  double dFloatBeta;
  int iUseFullCopula;

  /* For Quanto Part */
  int iNbQuantoCorrel;
  double *dQuantoMat;
  double *dQuantoCorrel1;
  double *dQuantoCorrel2;
  double *dFloatQuantoCorrel;

  int iNbFxVol;
  double *dFxMat;
  double *dFxVols;

  /* For Digital */
  double dDigitCallSpread;

  /* For Mixed SL */
  long lNbIter;

  /* For Copula */
  long lNbPaths;
  long lNbPoints;
  int iStudDegree;
  int iNConv;
  SrtMCSamType eMCType;

  /* For Copula SV */
  int iUseMC;
  double dMinTime;
  double dVolCorrel;
  double dCross1Shift;
  double dCross2Shift;

  /* For Copula MSL */
  long lNbPointsCopulaSL;
  double dNbSTDCopulaSL;

  /* For ATM adjstment */
  int iAdjustCorrelATM;
  int iSolveOnVol;
  double dATMPrecision;

  /* For Gaussian Copula */
  int iUseGaussianCopula;
  int iFreeSavedGaussian;
  COPULAGAUSSIAN_PARAMS sGaussianParams;

} cms_spread_model_params, *CMS_SPREAD_MODEL_PARAMS;

void cms_spread_set_default_model_params(CMS_SPREAD_MODEL_PARAMS sModelParams);

void free_cms_spread_model_params(cms_spread_model_params *params);

typedef struct {
  /* Forward Infos */
  double dCMS1;
  double dCMS2;
  double dFloatCMS;

  /* Correlation */
  double dCorrel;

  /* Normal Case */
  double dNormalVol1;
  double dNormalVol2;
  double dNormalSpreadVol;

  /* Lognormal Case */
  double dLogVol1;
  double dLogVol2;

  /* Shifter Log Case */
  double dSLVol1;
  double dBetaSL1;
  double dSLVol2;
  double dBetaSL2;

  /* Copula Case */
  double dAlpha1;
  double dBeta1;
  double dRho1;
  double dAlpha2;
  double dBeta2;
  double dRho2;

} cms_spread_infos, *CMS_SPREAD_INFOS;

typedef struct {
  long lFixingDate;
  long lStartDate;
  long lTheoEndDate;
  long lActEndDate;
  long lPayDate;
  int iIsCaplet;
  int iIsQuanto;

  double dMaturity;
  double dDelay;
  double dRateConv;
  double dNumPeriods;
  double dFwdZC;

  char *cFreq;
  char *cBasis;
  char *cRefRate;

  SrtCompounding sFreq;
  SrtBasisCode sBasis;

  double dLevel;
  double dForward;
  double dSpread;
  double dCMS;

  double dATMNormalVol;
  double dATMLogVol;

  int iNbSL;
  double dProba[CMSS_MAXSL];
  double dBeta[CMSS_MAXSL];
  double dATMSLVol[CMSS_MAXSL];
  double dShift[CMSS_MAXSL];
  double dConvFwd[CMSS_MAXSL];
  double dShiftedCMS[CMSS_MAXSL];

  double dSigmaBeta;
  double dSABRAlpha;
  double dSABRBeta;
  double dSABRRho;

} cms_struct, *CMS_STRUCT;

Err cms_init_cms_struct(long lFixingDate, long lStartDate, long lTheoEndDate,
                        long lPayDate, char *cFreq, char *cBasis,
                        char *cRefRate, int iIsQuanto, CMS_SPREAD_MKT cCmsMkt,
                        CMS_STRUCT cCmsStruct);

Err calc_cms_struct(CMS_STRUCT cCmsStruct, CMS_SPREAD_MODEL_PARAMS sModelParams,
                    CMS_SPREAD_MKT cCmsMkt);

/* Pricing of (gf * CMSf + mf) times an option or a digital on g1 * CMS1 + g2 *
 * CMS2 */
Err cms_spread_option_digital(
    /* Type */
    int iIsOption,

    long lFixingDate, long lStartDate, long lPayDate,

    /* CMS1 Description */
    long lTheoEndDate1, char *cFreq1, char *cBasis1, char *cRefRate1,
    int iIsQuanto1,

    /* CMS2 Description */
    long lTheoEndDate2, char *cFreq2, char *cBasis2, char *cRefRate2,
    int iIsQuanto2,

    /* Floating CMS Description */
    long lFloatFixingDate, long lFloatStartDate, long lFloatTheoEndDate,
    char *cFloatFreq, char *cFloatBasis, char *cFloatRefRate,
    int iFloatIsQuanto,

    /* Payoff Description */

    /* Floating part */
    double dFloatGear, double dFloatMargin,

    /* Option part */
    double dGear1, double dGear2, int iNbStrike, double *dStrike,
    int iCallPut,   /* 0: Call      , 1: Put */
    double dPayRec, /* Pay: -1.0      , Rec: +1.0 */

    /* Model Description */
    CMS_SPREAD_MODEL_PARAMS sModelParams,

    /* Market Informations */
    CMS_SPREAD_MKT sDomMarket, CMS_SPREAD_MKT sQuantoMarket,

    /* Outputs */
    double *dValue, CMS_SPREAD_INFOS sInfos);

Err cms_spread_option_struct(
    /* CMS */
    CMS_STRUCT sCms1, CMS_STRUCT sCms2,

    /* Payoff Description */
    double dGear1, double dGear2, int iNbStrike, double *dStrike,
    int iCallPut,   /* 0: Call      , 1: Put */
    double dPayRec, /* Pay: -1.0      , Rec: +1.0 */

    /* Model Description */
    CMS_SPREAD_MODEL_PARAMS sModelParams,

    /* Outputs */
    double *dValue, CMS_SPREAD_INFOS sInfos);

Err cms_spread_digital_struct(
    /* CMS */
    CMS_STRUCT sCms1, CMS_STRUCT sCms2,

    /* Payoff Description */
    double dGear1, double dGear2, int iNbStrike, double *dStrike,
    int iCallPut,   /* 0: Call      , 1: Put */
    double dPayRec, /* Pay: -1.0      , Rec: +1.0 */

    /* Model Description */
    CMS_SPREAD_MODEL_PARAMS sModelParams,

    /* Outputs */
    double *dValue, CMS_SPREAD_INFOS sInfos);

Err cms_spread_copula_mixed_sl(double dFwd1, double dShift1, double dProba1,
                               double dStd11, double dStd12,

                               double dFwd2, double dShift2, double dProba2,
                               double dStd21, double dStd22,

                               double dRho,

                               double dGear1, double dGear2, int iNbStrikes,
                               double *dStrikes,

                               CMS_SPREAD_MODEL_PARAMS sParams,

                               double *dValues);

Err cms_geared_spread_option_struct(
    /* CMS */
    CMS_STRUCT sFloatCms, CMS_STRUCT sCms1, CMS_STRUCT sCms2,

    /* Payoff Description */
    double dFloatGear, double dMargin, double dGear1, double dGear2,
    int iNbStrike, double *dStrike, int iCallPut, /* 0: Call      , 1: Put */
    double dPayRec, /* Pay: -1.0      , Rec: +1.0 */

    /* Model Description */
    CMS_SPREAD_MODEL_PARAMS sModelParams,

    /* Outputs */
    double *dValue, CMS_SPREAD_INFOS sInfos);

Err cms_geared_spread_digital_struct(
    /* CMS */
    CMS_STRUCT sFloatCms, CMS_STRUCT sCms1, CMS_STRUCT sCms2,

    /* Payoff Description */
    double dFloatGear, double dMargin, double dGear1, double dGear2,
    int iNbStrike, double *dStrike, int iCallPut, /* 0: Call      , 1: Put */
    double dPayRec, /* Pay: -1.0      , Rec: +1.0 */

    /* Model Description */
    CMS_SPREAD_MODEL_PARAMS sModelParams,

    /* Outputs */
    double *dValue, CMS_SPREAD_INFOS sInfos);

typedef struct {
  CMS_STRUCT sCms1;
  CMS_STRUCT sCms2;

  CMS_SPREAD_MODEL_PARAMS sModelParams;

  double dGear1;
  double dGear2;
  double dForward;
  double dPrice;
  double dImpliedVol;

  int iSolveOnVol;

  double dInitialCorrel;

} cmsatmcalib_params, *CMSATMCALIBPARAMS;

Err cms_adjust_correl_atm(CMS_STRUCT sCms1, CMS_STRUCT sCms2, double dGear1,
                          double dGear2, int iMinMax,
                          /* Model Description */
                          CMS_SPREAD_MODEL_PARAMS sModelParams);

#endif