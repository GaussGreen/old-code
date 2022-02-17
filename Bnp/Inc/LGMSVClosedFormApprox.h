
#ifndef LGMSVCLOSEDFORMAPPROX_H
#define LGMSVCLOSEDFORMAPPROX_H

#include "DiagCalibDLM.h"
#include "DiagCalibGen.h"
#include "LGMSVUtil.h"

#define LGMSV_FIRSTBREAK 10
#define LGMSV_LASTBREAK 1000000

void lgmsv_app_set_default_params(int *iNbX, double *iNbSigmaXLeft,
                                  double *iNbSigmaXRight, double *dIntegParam,
                                  int *iIntegMethod, double *dVolLimit,
                                  int *iCalibLGM, double *dMinStd,
                                  double *dMaxStd);

typedef struct {
  int iNbX;
  double dParam1;
  double dParam2;
  double dIntegParam;
  double dVolLimit;
  int iCalibLGM;
  int iIntegMethod;

  double dMinStd;
  double dMaxStd;

} LGMSV_NumerParams, *LGMSV_NUMERPARAMS;

void lgmsv_app_set_default_params_struct(LGMSV_NUMERPARAMS NumerParams);

void lgmsv_init_numer_params(LGMSV_NUMERPARAMS NumerParams, int iNbX,
                             double dParam1, double dParam2, double dIntegParam,
                             double dVolLimit, int iCalibLGM, int iIntegMethod,
                             double dMinStd, double dMaxStd);

Err LGMSVCmsRateApprox(
    char *und_name,      /*	Name of the underlying */
    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lFixDate, long lStartDate, long lEndDate, long lPayDate,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *dFwd);

Err LGMSVOptionApprox(
    char *und_name,      /*	Name of the underlying */
    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lExDate, long lStartDate, long lEndDate, long lPayDate, int IsCMS,
    long lRStartDate, long lREndDate, int iIsRatioNum, int iNbStrike,
    double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice);

Err LGMSVLiborOptionApprox(
    char *und_name,               /*	Name of the underlying */
    char *yc_name,                /*	Name of the yield curve */
    char *swaption_ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,
    char *libor_ref_rate_name, /*	Name of the reference rate */
    char *libor_freq,          /*	Frequency and basis of underlying swaptions */
    char *libor_basis,

    long lExDate, long lStartDate, long lEndDate, long lPayDate, int IsCMS,
    long lLiborFixDate, long lLiborStartDate, long lLiborEndDate, int iNbStrike,
    double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice);

Err LGMSVOptionApprox_TS(

    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lExDate, long lStartDate, long lEndDate, long lPayDate, int IsCMS,
    int iNbStrike, double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    double dLambdaX, int iNbPWTime, /* Piece Wise Term Structures  */
    double *dPWTime, double *dSigmaTS, double *dAlphaTS, double *dLambdaEpsTS,
    double *dRhoTS, double dTStar,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice);

Err LGMSVCmsRateApprox_TS(

    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lFixingDate, long lStartDate, long lEndDate, long lPayDate,

    double dLambdaX, int iNbPWTime, /* Piece Wise Term Structures  */
    double *dPWTime, double *dSigmaTS, double *dAlphaTS, double *dLambdaEpsTS,
    double *dRhoTS, double dTStar,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice);

Err LGMSVLiborOptionApprox_TS(

    char *yc_name,                /*	Name of the yield curve */
    char *swaption_ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis, char *libor_ref_rate_name, char *libor_freq,
    char *libor_basis,

    long lExDate, long lStartDate, long lEndDate, long lPayDate, int iIsCMS,
    long lLiborFixDate, long lLiborStartDate, long lLiborEndDate,

    int iNbStrike, double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    double dLambdaX, int iNbPWTime, /* Piece Wise Term Structures  */
    double *dPWTime, double *dSigmaTS, double *dAlphaTS, double *dLambdaEpsTS,
    double *dRhoTS, double dTStar,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice);

Err Calculate_HestonEquivalent_LGMSV(
    int iNbCoupon, long *lCouponDate, double *dCouponTime, double *dCouponCvg,
    double *dCouponDf, double dTPay, int iIsCMS, long lToday, char *cYcName,
    double dTau, double dTStar, double *dFwdSwapCash, double *dLevel,
    double *dShift, double *dCoefVol, double *dCoefMeanRev, double *dCoefCMS);

void LGMSVClosedFormApprox(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    double dLambdaX,

    int iNbPWTime, /* Piece Wise Term Structures  */
    double *dPWTime, double *dSigma, double *dAlpha, double *dLambdaEps,
    double *dRho,

    double dTStar, /* Tstar in years from today */

    /* Product description */
    double dFwd, int iNbStrike, double *dStrike, double dExTime,
    double dCoefVol,

    double dSwitchTime, double dCoefMeanRev1, double dCoefCMS1,
    double dCoefMeanRev2, double dCoefCMS2,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Outputs */
    double *Price);

Err Construct_Swap_Schedule(long lTodayDate, char *cYcName, long lStartDate,
                            long lEndDate, SrtCompounding iFreq,
                            SrtBasisCode iBasis, int *iNbCoupon,
                            long **lCouponDate, double **dCouponTime,
                            double **dCouponCvg, double **dCouponDf);

typedef struct {
  int iNbCoef;

  double *dt;
  double *CoefIntT;
  double *Coef1ReT;
  double *Coef2ReT;
  double *Coef2ImT;
  double *Coef3ReT;
  double *Coef3ImT;

  int iNbX;
  double *X;
  double *W;
  double *IntRe;
  double *IntIm;

  /* special case for low volatility */
  double *XX;
  double *WW;

  LGMSV_NUMERPARAMS NumerParams;

} LGMSV_PricingConst, *LGMSV_PRICINGCONST;

Err Initialise_PricingConst(LGMSV_MODEL model, LGMSV_NUMERPARAMS numer_params,
                            LGMSV_PRICINGCONST PricingConst);

void Free_PricingConst(LGMSV_PRICINGCONST PricingConst);

typedef struct {
  double dShift;
  double dCoefVol;
  int iNegVol;

  double dCoefMeanRev;
  double dCoefCMS;
  double dNewCoefCMS;

  double dSwitchTime;
  double dCoefMeanRev2;
  double dCoefCMS2;
  double dNewCoefCMS2;

  double dMinStd;
  double dMaxStd;
  double dATMNormalStd;
  double dMinStrike;
  double dMaxStrike;
  double dMinVol;
  double dMaxVol;

  /* for LGM2F */
  double dCoefVol_2F;

  double dCoefMeanRev_2F;
  double dCoefCMS_2F;
  double dNewCoefCMS_2F;

  double dCoefMeanRev2_2F;
  double dCoefCMS2_2F;
  double dNewCoefCMS2_2F;

  double dNewCoefCMS_cross_2F;
  double dNewCoefCMS2_cross_2F;

  /* Precalculations for LGM */
  double dCpnG1[MAX_CPN];
  double dCpnG2[MAX_CPN];

  double dExG1;
  double dExG2;

} LGMSV_HestonInst, *LGMSV_HESTONINST;

typedef struct {
  int iEndIndex;
  int iDoSwitch;
  int iSwitchIndex;
  double dNewSwitchTime;
  int iLiborIndex;
  int iNbStrike;
  double dLogShiftFwd;
  double *dLogShiftStrike;

} LGMSV_NumerInst, *LGMSV_NUMERINST;

struct InstParams;
typedef struct InstParams InstParams;
typedef InstParams *INSTPARAMS;

struct InstParams {
  CALIBCPNSCHEDULEDLM sCalibCpnSchedule;
  LGMSV_HESTONINST HestonInst;
  LGMSV_NUMERINST NumerInst;

  InstParams *NextInstParams;
};

Err Initialise_NumerInst(LGMSV_MODEL model, CALIBINSTRUMENTDLM CalibInst,
                         LGMSV_HESTONINST HestonParam,
                         LGMSV_NUMERINST NumerInst);

void Update_NumerInst(CALIBINSTRUMENTDLM CalibInst,
                      LGMSV_HESTONINST HestonParam, LGMSV_NUMERINST NumerInst);

Free_NumerInst(LGMSV_NUMERINST NumerInst);

Err Calculate_ShiftAndVol_LGMSV_FromLGM_struct(/* product information */
                                               CALIBINSTRUMENTDLM CalibInst,
                                               CALIBCPNSCHEDULEDLM
                                                   CalibCpnSchedule,

                                               /* model information */
                                               LGMSV_MODEL model,

                                               /* output */
                                               LGMSV_HESTONINST HestonParam);

Err Calculate_HestonEquivalent_LGMSV_FromLGM_struct(
    /*	Market information */
    /* product information */
    CALIBINSTRUMENTDLM CalibInst, CALIBCPNSCHEDULEDLM CalibCpnSchedule,

    /* model information */
    LGMSV_MODEL model,

    /* output */
    LGMSV_HESTONINST HestonParam);

void LGMSVClosedFormApprox_struct(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    LGMSV_MODEL model,

    /* Product description */
    CALIBINSTRUMENTDLM CalibInst, LGMSV_HESTONINST HestonParam,

    /* Parameter of grids */
    LGMSV_PRICINGCONST NumerParams, LGMSV_NUMERINST InstNumer,

    /* Outputs */
    double *Price);

void LGMSVFwdClosedFormApprox_struct(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    LGMSV_MODEL model,

    /* Product description */
    CALIBINSTRUMENTDLM CalibInst, LGMSV_HESTONINST HestonParam,

    /* Parameter of grids */
    LGMSV_PRICINGCONST NumerParams, LGMSV_NUMERINST InstNumer,

    /* Outputs */
    double *Price);

Err ConstructSingleInstrumentSchedule(char *sYcName, long lToday,
                                      char *cInstFreq, char *cInstBasis,
                                      long lExeDate, long lStartDate,
                                      long lTheoEndDate,
                                      CALIBCPNSCHEDULEDLM CalibCpnSchedule,
                                      CALIBEXESCHEDULEDLM CalibExeSchedule);

Err LGMSVOptionApprox_new(
    char *und_name,      /*	Name of the underlying */
    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lExDate, long lStartDate, long lTheoEndDate, long lPayDate, int IsCMS,
    long lRStartDate, long lREndDate, int iIsRatioNum, int iNbStrike,
    double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice);

Err LGMSVOptionApproxFromModel(
    LGMSV_model *sModel,

    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lExDate, long lStartDate, long lTheoEndDate, long lPayDate, int IsCMS,
    long lRStartDate, long lREndDate, int iIsRatioNum, int iNbStrike,
    double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice);

Err LGMSVCmsRateApprox_new(
    char *und_name,      /*	Name of the underlying */
    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lFixDate, long lStartDate, long lTheoEndDate, long lPayDate,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *dCmsRate);

Err Calculate_CoefShiftAndVol_LGMSV(/* product information */
                                    double dExeTime, int iNbCoupon,
                                    double *dCouponTime, double *dCouponCvg,
                                    double *dCouponDf,

                                    /* model information */
                                    double dTau,
                                    int iNbPWTime, /* Piece Wise Term Structures
                                                    */
                                    double *dPWTime, double *dSigmaTS,
                                    double dTStar,

                                    /* output */
                                    double *dLevel, double *dSwpCash,
                                    double *dShift, double *dCoefVol);

void Calculate_CoefCMSandMeanRev_LGMSV(/* product information */
                                       int iNbCoupon, double *dCouponTime,
                                       double *dCouponCvg, double *dCouponDf,
                                       double dTPay, int iIsCMS,
                                       double dTRStart, double dTREnd,
                                       int iIsRatioNum,

                                       /* model information */
                                       double dTau, double dTStar,
                                       double dCoefVol,

                                       /* output */
                                       double *dCoefMeanRev, double *dCoefCMS);

void LGMSVConstructLegendreGrid(int iNbInt, double *dBreakPoints,
                                int *iNbPoints, double *dX, double *dW);

#endif