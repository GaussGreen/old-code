#ifndef __CMSCTS_PROD_STRUCT_H
#define __CMSCTS_PROD_STRUCT_H

#include "CMSSpreadOptions.h"
#include "LGMSVCalibApprox.h"
#include "MCEBOptimisation.h"
#include "srt_h_all.h"

/* Market Structure */
typedef struct {
  long lToday;

  /* Market Curves and defaults */
  char cYcName[XL_MAX_ID_LEN];
  char cVcName[XL_MAX_ID_LEN];
  char cRefRateName[XL_MAX_ID_LEN];
  char cFixedFreq[XL_MAX_ID_LEN];
  char cFixedBasis[XL_MAX_ID_LEN];
  char cFloatBasis[XL_MAX_ID_LEN];

  /* Function to get volatility */
  char *(*GetCashVol)(char *cVolCurve, double dStartDate, double dEndDate,
                      double dCashStrike, int iZero, char *cRefRate,
                      double *dVol, double *dPower);

  /* CMS Informations */
  long lNumStrikesInVol; /*	Array of strikes in vol matrix */
  double *dStrikesInVol;
  SrtDiffusionType eVolType; /*	Type of vol in matrix      , SRT_NORMAL or
                                SRT_LOGNORMAL */
  int iCashVol; /*	1: matrix is a cash vol 0: matrix is a swap vol */
  int iIsSABR;

  /* Foreign Funding Informations */
  char cFundYcName[XL_MAX_ID_LEN];
  double dFundSpotFx;
  long lFundFxSpotDate;

  /* EOD flags 0: Intraday      , 1: End of Day */
  int iEODFixFlag;
  int iEODPayFlag;
  int iEODExeFlag;

} cmscts_market, *CMSCTS_MARKET;

void cmscts_free_market(CMSCTS_MARKET sMarket);

/* Deal Structures */
/* *************** */

/* Funding */
typedef struct {
  int iStartIndex;
  int iEndIndex;
  int iNbUsedCpn;

  int iNumStartIndex;
  int iNbUsedNumCpn;

  double dCashFx;
  char cFundYcName[XL_MAX_ID_LEN];
  double *dNot;
  double *dMargin;
  double *dFixCoupon;
  double *dInitSpread;
  double *dNewSpread;
  long lStartDate;
  double dInitEquiEx;
  double dFinalEquiEx;

  double *dCouponPlusEx;
  double *dCouponPlusExPartial;

  double *dMarketValue;

  int *iFundingControll;

} cmscts_fund_aux, *CMSCTS_FUND_AUX;

typedef struct {
  int iCCY; /*	0: domestic      , 1: other */
  char cRefRate[20];

  int iNbCpn;
  double *dNot;
  long *lFixDate;
  long *lStartDate;
  long *lEndDate;
  long *lPayDate;
  double *dCoverage;
  double *dMargin;
  double
      *dFixCoupon; /*	Past coupon fixing if relevant      ,
                      includes spr      , but not mrg      , cvg and notional */

  /* calculation aux */
  CMSCTS_FUND_AUX sAux;

} cmscts_funding, *CMSCTS_FUNDING;

void cmscst_funding_free_aux(CMSCTS_FUNDING sFunding);
void cmscts_free_funding(CMSCTS_FUNDING sFunding);

/* Exotic Coupons */
typedef enum CMSCTS_COUPON_TYPE_ {
  CMSCTS_FLOATER = 0,
  CMSCTS_TIMESWAP

} CMSCTS_COUPON_TYPE;

Err cmscts_interp_coupon_type(char *cCouponType,
                              CMSCTS_COUPON_TYPE *eCouponType);

typedef struct {
  int iStartIndex;
  int iEndIndex;
  int iNbUsedCpn;

  int iNumStartIndex;
  int iNbUsedNumCpn;

  /* Accrual Def */
  int *iFloatingType; /* 0: no floating      , 1: floating in advance      , 2:
                         floating in arrears */

  int *iIVNbTrimFix;
  int **iIVTrimIdx;
  double **dIVTrimWeights;

  int *iIndexFixed;
  int *iFullyFixed;

  int *iCallNbTrimFix;
  int **iCallTrimIdx;
  double **dCallTrimWeights;

  /* Coupon Payoff = dAlpha + dBeta * CMSS + sum(nbOpt * (CMSS - K)+) */
  double *dAlpha;
  double *dBeta;
  int *iNbStrikes;
  double **dStrikes;
  double **dNbOpt;

  /* Floor / Cap for Time Swap Coupon */
  int iFirstFloorCap;
  int *iHasFloorCap;

  /* Controll */
  int *iExoticControll;

  /* Market Value */
  double *dMarketValue;

  /* Partial PV for Target Note */
  double dFixedPV;
  double dPartialPV;
  double dExtraTarget;

  /* Adjusted Floating Fixing Date */
  long *lNewPaidFixingDate;

} cmscts_exotic_aux, *CMSCTS_EXOTIC_AUX;

typedef struct {
  /* Schedule */
  int iNbCpn;
  long *lStartDate;
  long *lEndDate;
  long *lPayDate;
  double *dCoverage;
  double *dNot;

  /* Coupons */
  CMSCTS_COUPON_TYPE *eCouponType;

  /* Paid CMS */
  int iFixingLag;
  int iNbPaidCMS;
  long *lPaidFixingDate;
  double **dPaidGearing;
  char ***cPaidTenor;
  char **cPaidFreq;
  char **cPaidBasis;
  char **cPaidRef;
  double *dPaidMargin;
  double *dFloor;
  double *dCap;

  double *dPastPaidSpreadCMS;

  /* Accrual Def */
  int iAccuralFixLag;
  char cAccrualFixFreq[5];
  int iCountWeekEnd;

  int *iNbFix;
  double **dWeights;
  long **lFixDates;

  int iNbCMS;
  double **dGear;
  char ***cTenor;
  char **cCMSFreq;
  char **cCMSBasis;
  char **cCMSRef;
  double *dLowBarrier;
  double *dHighBarrier;

  double **dPastAccrualSpread;

  /* calculation aux */
  CMSCTS_EXOTIC_AUX sAux;

} cmscts_exotic, *CMSCTS_EXOTIC;

void cmscts_free_exotic(CMSCTS_EXOTIC sExotic);

/* Right to Call */
typedef struct {
  int iStartIndex;
  int iNbUsedCall;

  /* Calling Details */
  int *iStartCpnIdx;
  int *iStartFundIdx;

  /* Exercised Details */
  int iIndexExercised;

} cmscts_call_aux, *CMSCTS_CALL_AUX;

typedef struct {
  double dPayRec; /* -1: Pay      , +1: Rec */

  int iNbCall;
  long *lExeDate;
  long *lSettlDate;
  double *dFee;

  /* Exercised */
  int iIsExercised;
  long lExercisedDate;

  /* aux for calculation */
  CMSCTS_CALL_AUX sAux;

} cmscts_call, *CMSCTS_CALL;

void cmscts_free_call(CMSCTS_CALL sCall);

/* Exotic Coupons */
typedef enum CMSCTS_KO_TYPE_ {
  CMSCTS_KO_CMS = 0,
  CMSCTS_KO_COUPON

} CMSCTS_KO_TYPE;

Err cmscts_interp_KO_type(const char *constStr, CMSCTS_KO_TYPE *eKOType);

typedef enum CMSCTS_TARGET_TYPE_ {
  CMSCTS_TARGET_FULL_COUPON = 0,
  CMSCTS_TARGET_FULL_GUARANTED,
  CMSCTS_TARGET_LAST_GUARANTED,

} CMSCTS_TARGET_TYPE;

Err cmscts_interp_TARGET_type(const char *constStr,
                              CMSCTS_TARGET_TYPE *eTargetType);

typedef struct {
  int iStartIndex;
  int iNbUsedKO;

  /* KO Details */
  long *lNewKODates;
  int *iIsReallyKO;
  int *iStartCpnIdx;
  int *iStartFundIdx;

  /* Coupon Payoff = dAlpha + dBeta * CMSS + sum(nbOpt * (CMSS - K)+) */
  double *dAlpha;
  double *dBeta;
  int *iNbStrikes;
  double **dStrikes;
  double **dNbOpt;

  /* KOed ? */
  int iHasKO;
  int iKOIndex;
  double dCumulKOTarget;

} cmscts_ko_aux, *CMSCTS_KO_AUX;

typedef struct {
  int iNbKO;

  CMSCTS_KO_TYPE eKOType;

  int iNbKOCMS;
  int iFixingLag;
  long *lFixingDate;
  long *lSettlDate;

  double **dKOGearing;
  char ***cKOTenor;
  char **cKOFreq;
  char **cKOBasis;
  char **cKORef;

  double *dLowBarrier;
  double *dHighBarrier;
  double *dTarget;

  int iTargetGuaranted;
  int iTargetCapped;

  double *dPastKOSpreadCMS;

  CMSCTS_KO_AUX sAux;

} cmscts_ko, *CMSCTS_KO;

void cmscts_free_ko(CMSCTS_KO sKO);

/* Deal structure */
typedef struct {
  long lStartDate;
  long lTheoEndDate;
  double dPayRec;

  /* Funding */
  CMSCTS_FUNDING sFunding;

  /* Exotic */
  CMSCTS_EXOTIC sExotic;

  /* Call */
  CMSCTS_CALL sCall;

  /* KO */
  CMSCTS_KO sKO;

} cmscts_deal, *CMSCTS_DEAL;

Err cmscts_allocate_deal(CMSCTS_DEAL sDeal);
void cmscts_free_deal(CMSCTS_DEAL sDeal);

// added by J.B
Err cmscts_funding_allocate_aux(CMSCTS_FUNDING sFunding);
Err cmscts_check_funding(CMSCTS_FUNDING sFunding);

typedef enum CMSCTS_CALIB_STRAT_ {
  CMSCTS_DIAGCAP = 0,
  CMSCTS_DIAGLONG,
  CMSCTS_DIAGSHORT,
  CMSCTS_LONGSHORT,
  CMSCTS_SHORTLONG,
  CMSCTS_CMS1CMS2,
  CMSCTS_CMS2CMS1,
  CMSCTS_DIAGCMS1,
  CMSCTS_DIAGCMS2,
  CMSCTS_CMS1DIAG,
  CMSCTS_CMS2DIAG,
  CMSCTS_CUSTOM

} CMSCTS_CALIB_STRAT;

Err cmscts_interp_calib_strat(const char *constStr, CMSCTS_CALIB_STRAT *val);

typedef struct {
  char cCalibTenor[10];
  long lEndDate;
  char cFreq[10];
  char cBasis[10];
  char cRefRate[10];

} cmscts_custom_calib, *CMSCTS_CUSTOM_CALIB;

/* Calibration / Model Parameters */
typedef struct {
  /* Precalibrated Und */
  int iUseCalib;
  char *cUndName;
  int iUseSV;

  /* LGM Parameters */
  int iNbFactor;
  double dLambda;
  double dLGMAlpha;
  double dLGMGamma;
  double dLGMRho;

  /* LGMSV Parameters */
  int iNbTimesSV;
  double *dTimesSV;
  double *dAlphaSV;
  double *dLamSV;
  double *dRhoSV;
  double *dRho2SV;
  double dTStar;

  /* Calibration Parameters */
  int iFixLambda;
  CPD_DIAG_CALIB_PARAM sPrimCalibParams;
  CPD_DIAG_CALIB_PARAM sSecCalibParams;
  DIAG_CALIB_LM_PARAMS sLMParams;
  LGMSV_CALIBPARAMS sLGMSVCalibParams;

  /* Numerical Parameters */
  LGMSV_NUMERPARAMS sLGMSVNumerParams;

  /* Calibration strategy */
  CMSCTS_CALIB_STRAT eCalibStrat;
  char cMaxCalibTenor[5];
  double dCalibSmileSTD;
  int iUseInArrearsFixing;

  /* Custom Calibration */
  CMSCTS_CUSTOM_CALIB sPrimInst;
  CMSCTS_CUSTOM_CALIB sSecInst;

} cmscts_calib, *CMSCTS_CALIB;

Err cmscts_allocate_calib(CMSCTS_CALIB sCalib);
void cmscts_free_calib(CMSCTS_CALIB sCalib);

/* All Pricing Params */
/* ****************** */

/* For Coupon Pricing */
typedef struct {
  /* Triming for Accrual */
  /* 0: no trim      , 1: x fixings max      , 2: x time min between two fixings
   */
  int iIVTrimType;
  int iIVMaxFix;
  double dIVMinFixTime;

  int iCallTrimType;
  int iCallMaxFix;
  double dCallMinFixTime;

  /* Digital Params */
  int iValueZero;
  double dIVCallSpread;
  int iIVBuySell;
  double dCallCallSpread;
  int iCallBuySell;
  int iAccrueOnBarrier;

  double dMinBarrier;
  double dMaxBarrier;

  /* Call Params */
  double dMCEBCallSpread;

  /* KO Params */
  int iKOBuySell;
  double dKOCallSpread;

  /* Pricing */
  CMS_SPREAD_MODEL_PARAMS sCMSSpreadModel;

  /* Correlations */
  int iUseLGMCorrel;
  double dLGMCorrelAlpha;
  double dLGMCorrelGamma;
  double dLGMCorrelRho;

  int iNbCorrPair;
  char **cPairNames;

  int iNbCorrTimes;
  double *dCorrTimes;
  double **dPairCorrels;

} cmscts_cpn_numparams, *CMSCTS_CPN_NUMPARAMS;

/* For Monte Carlo */
typedef struct {
  double dMCMaxTime;
  long lNbSteps;
  long lNbPaths;

  int iAdjustIV;
  int iNumeraireType;
  int iEstimateCallOnly;

  LGMSVPARAM sLGMSVParams;

  /* Params for CMS */
  int iTrimFreq;

} cmscts_simul_params, *CMSCTS_SIMUL_PARAMS;

typedef struct {
  int iCalcSwitchAdjust;
  int iGMANbFactor;
  double dGMAAlpha;
  double dGMAGamma;
  double dGMARho;

  CMSCTS_CALIB_STRAT eCalibStrat;
  char cLGMModelName[256];
  int iLGMNbFactor;
  double dLGMLambda;
  int iLGMFixLambda;
  double dLGMAlpha;
  double dLGMGamma;
  double dLGMRho;

  double dLGMLambdaShift;

} cmscts_reserve_params, *CMSCTS_RESERVE_PARAMS;

typedef struct {
  CMSCTS_CPN_NUMPARAMS sCouponPricingParams;
  CMSCTS_SIMUL_PARAMS sSimulParams;
  CMSCTS_RESERVE_PARAMS sReserveParams;

} cmscts_pricing_params, *CMSCTS_PRICING_PARAMS;

Err cmscts_allocate_pricing_params(CMSCTS_PRICING_PARAMS sPricingParams);
void cmscts_free_pricing_params(CMSCTS_PRICING_PARAMS sPricingParams);

/* Function to check and fill the structures */
Err cmscts_fill_check_all_struct(CMSCTS_DEAL sDeal, CMSCTS_MARKET sMarket,
                                 CMSCTS_CPN_NUMPARAMS sNumerParams);

typedef struct {
  long lNbCall;
  long *lExeDates;
  double *dMarketIV;
  double *dModelIV;
  double *dNewModelIV;
  double *dFees;
  double *dOTCall;
  double *dOTCallPartial;

} cmscts_fwdiv, *CMSCTS_FWDIV;

Err cnscts_allocate_fwdiv(CMSCTS_FWDIV sFwdIVInfos, long lNbCall);

void cmscts_free_fwdiv(CMSCTS_FWDIV sFwdIVInfos);

typedef struct {
  double dFunding;
  double dExotic;
  double dCall;
  double dStd;

  double dNewIV;
  double dCallInit;
  double dSwitchAdjust;
  double dGMAInitCall;

  CPD_CALIB_INST_DATA sInstDatas;
  LGMSV_MODEL sLGMSVModel;

  int iCalcExeProbas;

  CMSCTS_FWDIV sFwdIVInfos;
  CPD_CALIB_INST_DATA sSwitchInstDatas;

} cmscts_outputs, *CMSCTS_OUTPUTS;

Err cmscts_allocate_outputs(CMSCTS_OUTPUTS sOutputs);
void cmscts_free_outputs(CMSCTS_OUTPUTS sOutputs);

#endif