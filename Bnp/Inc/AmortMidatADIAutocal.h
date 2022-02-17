// prevent multiple inclusions
#pragma once

//////////////////////
//	warnings
#pragma warning(disable : 4786) //"identifier was truncated to '255' characters
                                // in the debug information"
// NB: force warnings for unused arguments and local parameters
#pragma warning(1 : 4100 4101)

#include "AmortMidatADICalib.h"
#include "AmortMidatADIUtils.h"

/////////////////////////////////////////////////
/// structure that contains all necessary arguments
///	to be passed to MidAt_AutoCal
typedef struct _Comm_MidAtAutocal {
  ///////Market
  const long m_lToday;
  const char *m_szYC;
  const char *m_szVC;
  Err (*m_get_cash_vol)(char *, double, double, double, int, char *, double *,
                        double *); //	functor to get cash vol from the market
  /// MidAt General Specs
  const long m_lTheoEnd;
  const char *m_szRefRate; // Refrate
  const char
      *m_szFreq; // frequency of the fixed leg -"M"      ,"Q"      ,"S" ,"A"
  const char *m_szBasis;          // basis of the fixed leg
  const double *m_pdCoupon_Begin; // coupon
  const char *m_szPayRec;         //"REC" or "PAY"
  // MidAt Exercise
  const long *m_plEx_Begin; // exercise dates
  const long *m_plEx_End;
  const long *m_plExStart_Begin; // exercise premium dats
  // MidAt Fixed Leg
  const long *m_plFixPay_Begin;
  const long *m_plFixPay_End;
  const long *m_plFixStart_Begin;
  const long *m_plFixEnd_Begin;
  const double *m_pdFixCvg_Begin;
  // MidAt Funding/Floating leg
  const long *m_plFltPay_Begin;
  const long *m_plFltPay_End;
  const long *m_plFltStart_Begin;
  const long *m_plFltEnd_Begin;
  const double *m_pdFltCvg_Begin;
  const double *m_pdMargin_Begin;
  const double *m_pdSpread_Begin;
  // term structure
  const int m_nCalibTauToCap; // 1 - calibrate Tau to cap      , 0 otherwise
  const double *m_pdLamT_Begin;
  const double *m_pdLamT_End;
  const double *m_pdLam_Begin;
  const double *m_pdAlpha;
  const double *m_pdGamma;
  const double *m_pdRho;
  // calibration params
  const int m_nOnefequi;
  // grid
  const int m_nNumPointT;
  const int m_nNumPointX;
  /// parameters needed by diagcalib
  const cpd_diag_calib_param *m_pParam_Swaption;
  const cpd_diag_calib_param *m_pParam_Caplet;
  const diag_calib_lm_params *m_pCalibParam;
  /// m_results;
  double *m_pdOptionPV;
  double *m_pdFixLegPV;
  double *m_pdFltLegPV;
} Comm_MidAtAutocal;

const char *MIDAT_AUTOCAL(const Comm_MidAtAutocal *pComm);

const char *GenMidAt_AutoCal(
    /// Market
    long lToday, const char *szYC, const char *szVC,
    Err (*get_cash_vol)(char *, double, double, double, int, char *, double *,
                        double *),
    const char *szCorrel,
    Err (*get_correl)(char *, double, double, double, double *),
    /// MidAt Specs
    const char *szRefRate,
    int eod_flag, //	EOD Flag 0: I      , 1: E
    long lTheoEnd,
    const char *szPayRec, //"REC" or "PAY"
    // MidAt Exercise
    const long *plEx_Begin, // exercise dates
    const long *plEx_End,
    const long *plExStart_Begin, // exercise premium dats
    const double *pdExFee_Begin, // exercise fee
    // MidAt Fixed Leg
    const char *szFreq, const char *szBasis, const long *plFixStart_Begin,
    const long *plFixStart_End, const long *plFixEnd_Begin,
    const double *pdFixCvg_Begin, const long *plFixPay_Begin,
    const double *pdFixNotional_Begin,
    const double *pdFixCoupon_Begin, // coupon
    const double *pdFixFee_Begin,
    // MidAt Funding/Floating leg
    const long *plFltStart_Begin, const long *plFltStart_End,
    const long *plFltEnd_Begin, const double *pdFltCvg_Begin,
    const long *plFltPay_Begin, const double *pdMargin_Begin,
    const double *pdSpread_Begin, const double *pdFltNotional_Begin,
    // Given Term Structure
    const long *plLamDate_Begin, const long *plLamDate_End,
    const double *pdLam_Begin, const double *pdAlpha, const double *pdGamma,
    const double *pdRho,
    /// Calibration Params
    const char *default_ref,        //	ref rate
    const char *default_swap_freq,  // swap freq
    const char *default_swap_basis, // swap basis
    double dMinTime, double dMinInterval, long lNotice, double dMaxStdShort,
    int nFixLam,   // 0: calib lambda to cap      , 1: fix lambda calib	to
                   // diagonal
    int n1FEqui,   // 1: 2F lambda will calibrate to the cap priced within
                   // calibrated 1F	with the given lambda
    int nSkipLast, // If 1      , the last option is disregardedand the forward
                   // volatility is flat from option	n-1
    int nUseJump, double dMaxVarJump, int nStrikeType, int nEuroModel,
    // grid
    int nNumPointT, int nNumPointX,
    /// Results      ,
    // NB: Set the sig related ptrs to 0 to calibrate sigma; else use the sigma
    // passed in !!!
    long **pplSigDate_Begin, long **pplSigDate_End, double **ppdSig_Begin,
    double *pdOptionPV, double *pdFixLegPV, double *pdFltLegPV,
    double *pdExProb, double *pdExBoundary, int nUse_Backward,
    int *pnBool_CalibInst, double *pdPrice_CalibInst);

const char *MidAt_AutoCal(
    ///////Market
    long lToday, const char *szYC, const char *szVC,
    Err (*get_cash_vol)(char *, double, double, double, int, char *, double *,
                        double *), //	functor to get cash vol from the market
    /// MidAt General Specs
    long lTheoEnd,
    const char *szRefRate, // Refrate
    const char
        *szFreq, // frequency of the fixed leg -"M"      ,"Q"      ,"S" ,"A"
    const char *szBasis,          // basis of the fixed leg
    const double *pdCoupon_Begin, // coupon
    const char *szPayRec,         //"REC" or "PAY"
    // MidAt Exercise
    const long *plEx_Begin, // exercise dates
    const long *plEx_End,
    const long *plExStart_Begin, // exercise premium dats
    // MidAt Fixed Leg
    const long *plFixPay_Begin, const long *plFixPay_End,
    const long *plFixStart_Begin, const long *plFixEnd_Begin,
    const double *pdFixCvg_Begin,
    // MidAt Funding/Floating leg
    const long *plFltPay_Begin, const long *plFltPay_End,
    const long *plFltStart_Begin, const long *plFltEnd_Begin,
    const double *pdFltCvg_Begin, const double *pdMargin_Begin,
    const double *pdSpread_Begin,
    // term structure
    int nCalibTauToCap, // 1 - calibrate Tau to cap      , 0 otherwise
    double *pdLamT_Begin, double *pdLamT_End, double *pdLam_Begin,
    const double *pdAlpha, const double *pdGamma, const double *pdRho,
    // calibration params
    int nOnefequi,
    // grid
    int nNumPointT, int nNumPointX,
    /// parameters needed by diagcalib
    const cpd_diag_calib_param *pParam_Swaption,
    const cpd_diag_calib_param *pParam_Caplet,
    const diag_calib_lm_params *pCalibParam,
    /// results      ,
    double *pdOptionPV, double *pdFixLegPV, double *pdFltLegPV,
    double *pdExProb, double *pdExBoundary, int nUse_Backward);

const char *MidAt_AutoCal_Diag(
    ///////Market
    long lToday, const char *szYC, const char *szVC,
    Err (*get_cash_vol)(char *, double, double, double, int, char *, double *,
                        double *), //	functor to get cash vol from the market
    /// MidAt General Specs
    long lTheoEnd,
    const char *szRefRate, // Refrate
    const char
        *szFreq, // frequency of the fixed leg -"M"      ,"Q"      ,"S" ,"A"
    const char *szBasis,          // basis of the fixed leg
    const double *pdCoupon_Begin, // coupon
    const char *szPayRec,         //"REC" or "PAY"
    // MidAt Exercise
    const long *plEx_Begin, // exercise dates
    const long *plEx_End,
    const long *plExStart_Begin, // exercise premium dats
    // MidAt Fixed Leg
    const long *plFixPay_Begin, const long *plFixPay_End,
    const long *plFixStart_Begin, const long *plFixEnd_Begin,
    const double *pdFixCvg_Begin,
    // MidAt Funding/Floating leg
    const long *plFltPay_Begin, const long *plFltPay_End,
    const long *plFltStart_Begin, const long *plFltEnd_Begin,
    const double *pdFltCvg_Begin, const double *pdMargin_Begin,
    const double *pdSpread_Begin,
    // term structure
    int nCalibTauToCap, // 1 - calibrate Tau to cap      , 0 otherwise
    const double *pdLam_Begin, const double *pdAlpha, const double *pdGamma,
    const double *pdRho,
    // calibration params
    int nOnefequi,
    // grid
    int nNumPointT, int nNumPointX,
    /// parameters needed by diagcalib
    const cpd_diag_calib_param *pParam_Swaption,
    const cpd_diag_calib_param *pParam_Caplet,
    const diag_calib_lm_params *pCalibParam,
    /// results      ,
    double *pdOptionPV_Begin);

const char *Euro_AutoCal_Diag(
    ///////Market
    long lToday, const char *szYC, const char *szVC,
    Err (*get_cash_vol)(char *, double, double, double, int, char *, double *,
                        double *), //	functor to get cash vol from the market
    /// MidAt General Specs
    long lTheoEnd,
    const char *szRefRate, // Refrate
    const char
        *szFreq, // frequency of the fixed leg -"M"      ,"Q"      ,"S" ,"A"
    const char *szBasis,          // basis of the fixed leg
    const double *pdCoupon_Begin, // coupon
    const char *szPayRec,         //"REC" or "PAY"
    // MidAt Exercise
    const long *plEx_Begin, // exercise dates
    const long *plEx_End,
    const long *plExStart_Begin, // exercise premium dats
    // MidAt Fixed Leg
    const long *plFixPay_Begin, const long *plFixPay_End,
    const long *plFixStart_Begin, const long *plFixEnd_Begin,
    const double *pdFixCvg_Begin,
    // MidAt Funding/Floating leg
    const long *plFltPay_Begin, const long *plFltPay_End,
    const long *plFltStart_Begin, const long *plFltEnd_Begin,
    const double *pdFltCvg_Begin, const double *pdMargin_Begin,
    const double *pdSpread_Begin,
    // term structure
    int nCalibTauToCap, // 1 - calibrate Tau to cap      , 0 otherwise
    const double *pdLam_Begin, const double *pdAlpha, const double *pdGamma,
    const double *pdRho,
    // calibration params
    int nOnefequi,
    // grid
    int nNumPointT, int nNumPointX,
    /// parameters needed by diagcalib
    const cpd_diag_calib_param *pParam_Swaption,
    const cpd_diag_calib_param *pParam_Caplet,
    const diag_calib_lm_params *pCalibParam,
    /// results      ,
    double *pdOptionPV_Begin);

const char *_MidAt_AutoCal(
    ///////Market
    long lToday, const char *szYC, const char *szVC,
    Err (*get_cash_vol)(char *, double, double, double, int, char *, double *,
                        double *), //	functor to get cash vol from the market
    /// MidAt General Specs
    long lTheoEnd,
    const char *szRefRate, // Refrate
    const char
        *szFreq, // frequency of the fixed leg -"M"      ,"Q"      ,"S" ,"A"
    const char *szBasis,          // basis of the fixed leg
    const double *pdCoupon_Begin, // coupon
    const char *szPayRec,         //"REC" or "PAY"
    // MidAt Exercise
    const long *plEx_Begin, // exercise dates
    const long *plEx_End,
    const long *plExStart_Begin, // exercise premium dats
    // MidAt Fixed Leg
    const long *plFixPay_Begin, const long *plFixPay_End,
    const long *plFixStart_Begin, const long *plFixEnd_Begin,
    const double *pdFixCvg_Begin,
    // MidAt Funding/Floating leg
    const long *plFltPay_Begin, const long *plFltPay_End,
    const long *plFltStart_Begin, const long *plFltEnd_Begin,
    const double *pdFltCvg_Begin, const double *pdMargin_Begin,
    const double *pdSpread_Begin,
    // term structure
    int nCalibTauToCap, // 1 - calibrate Tau to cap      , 0 otherwise
    double *pdLamT_Begin, double *pdLamT_End, double *pdLam_Begin,
    const double *pdAlpha, const double *pdGamma, const double *pdRho,
    // calibration params
    const char **pszSwapTenor_Begin, const double *pdSwapStrike_Begin,
    const char **pszCapTenor_Begin, const double *pdCapStrike_Begin,
    const cpd_diag_calib_param *pParam_Swaption,
    const cpd_diag_calib_param *pParam_Caplet,
    const diag_calib_lm_params *pCalibParam, int nOnefequi,
    // grid
    int nNumPointT, int nNumPointX,
    /// results      ,
    double *pdOptionPV, double *pdFixLegPV, double *pdFltLegPV,
    double *pdExProb, double *pdExBoundary, int nUse_Backward);
