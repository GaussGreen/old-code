#ifndef _CR_IF_UTIL_
#define _CR_IF_UTIL_

#include "CRIFProdStruct.h"

#define CRIF_MAX_FUND_CPN 780
#define CRIF_MAX_FLT_CPN 780
#define CRIF_MAX_ALL_CPN 1560

typedef struct {
  int iIndexCoupon;

  int iSavedCoupon;

  /* For payment */
  double tpay_tstar_alpha;
  double tpay_tstar_beta;
  double tpay_tstar_gamma;

  double tpay_tstar_beta_2;
  double tpay_tstar_gamma_2;
  double tpay_tstar_gamma_12;

  /* All DF needed by CMS */
  int iNbDF;
  double tdf_tstar_alpha[CRIF_MAX_ALL_CPN];
  double tdf_tstar_beta[CRIF_MAX_ALL_CPN];
  double tdf_tstar_gamma[CRIF_MAX_ALL_CPN];

  double tdf_tstar_beta_2[CRIF_MAX_ALL_CPN];
  double tdf_tstar_gamma_2[CRIF_MAX_ALL_CPN];
  double tdf_tstar_gamma_12[CRIF_MAX_ALL_CPN];

  double dAllDF[CRIF_MAX_ALL_CPN];

  /* Each CMS reconstruction */
  int iNbCMSDF;
  int iIndexCMSDF[CRIF_MAX_ALL_CPN];
  double dCMSCoverage[CRIF_MAX_ALL_CPN];
  double dSpread;

  /* For Payoff */
  double dCoef;
  double dAlpha;
  double dBeta;
  double dGamma;

  double dFloor;
  double dCap;

  /* For IV control Variate  */
  double dCVStrike;

} crif_couponarg, *CRIF_COUPONARG;

typedef struct {
  int iIndexCall;

  /* for funding */
  int iNbFundCoupon;
  int iNbFundCouponPartial;
  double tpay_tstar_alpha[CRIF_MAX_FUND_CPN];
  double tpay_tstar_beta[CRIF_MAX_FUND_CPN];
  double tpay_tstar_gamma[CRIF_MAX_FUND_CPN];

  double tpay_tstar_beta_2[CRIF_MAX_FUND_CPN];
  double tpay_tstar_gamma_2[CRIF_MAX_FUND_CPN];
  double tpay_tstar_gamma_12[CRIF_MAX_FUND_CPN];

  /* for settlement */
  double tset_tstar_alpha;
  double tset_tstar_beta;
  double tset_tstar_gamma;

  double tset_tstar_beta_2;
  double tset_tstar_gamma_2;
  double tset_tstar_gamma_12;

  /* Cpn infos for regression*/
  CRIF_COUPONARG mceb_sCouponArg;

} crif_callarg, *CRIF_CALLARG;

typedef struct {
  /* Disctretisation */
  int iNbTimes;
  double *dTimes;
  double *dDates;
  double *dSigma;
  double *dAlphaSV;
  double *dLambdaSV;
  double *dLevelSV;
  double *dRhoSV;
  double *dRho2SV;
  double *dLGMAlpha;
  double *dLGMRho;

  /* Payoff infos */
  int iNbEvent;
  int *iEvalEvent;
  int *iIndexEvent;
  void **vPayoffParams;
  double *dModelValue;

  /* Path Dependent Infos */
  double *dPathInfos;
  int iHasFloorCap;
  double *dSavedCoupon;
  double *dSavedDF;

  int iHasCallFee;

  /* MCEB Parameters */
  int iAdjustIV;
  MCEBPARAMS sMCEBParams;
  int *iOptimise;
  double *dMarketValue;

} crif_mcarg, *CRIF_MCARG;

Err crif_allocate_mc_arg(CRIF_MCARG sMCArg);
void crif_free_mc_arg(CRIF_MCARG sMCArg);

typedef struct {
  /* Model used for pricing */
  LGMSV_MODEL sModel;

  /* Model used for Tau Ts Adjustment */
  char cLGM_fixedTau_ModelName[256];
  char cLGM_TauTs_ModelName[256];

  /* All the arguments for MC */
  CRIF_MCARG sMCArg;

} crif_simularg, *CRIF_SIMULARG;

Err crif_allocate_simularg(CRIF_SIMULARG sSimulArg);
void crif_free_simularg(CRIF_SIMULARG sSimulArg);

typedef struct {
  int iIsCoupon;
  int iIsFixing;
  int iIsCall;

  CRIF_DEAL sDeal;
  CRIF_SIMULARG sSimulArg;
  int iNumeraireType;

  CRIF_COUPONARG sCouponArg;
  CRIF_CALLARG sCallArg;

  /* All Path dependent infos */
  double *dPathInfos;
  double *dSavedCoupon;
  double *dSavedDF;

  int iAdjustPayoff;

} crif_payarg, *CRIF_PAYARG;

Err crif_calibrate_model(CMSCTS_MARKET sMarket, CRIF_DEAL sDeal,
                         CRIF_SIMULARG sSimulArg, CMSCTS_CALIB sCalibration,
                         CRIF_OUTPUTS sOutputs);

Err crif_calculate_option(CMSCTS_MARKET sMarket, CRIF_DEAL sDeal,
                          CRIF_SIMULARG sSimulArg,
                          CRIF_PRICING_PARAMS sPricingParams,
                          CRIF_OUTPUTS sOutputs);

Err crif_calculate_switch_adjust(CMSCTS_MARKET sMarket, CRIF_DEAL sDeal,
                                 CRIF_SIMULARG sSimulArg,
                                 CMSCTS_CALIB sCalibration,
                                 CRIF_PRICING_PARAMS sPricingParams,
                                 CRIF_OUTPUTS sOutputs);

Err crif_fill_outputs(CRIF_DEAL sDeal, CRIF_SIMULARG sSimulArg,
                      CRIF_OUTPUTS sOutputs);

#endif