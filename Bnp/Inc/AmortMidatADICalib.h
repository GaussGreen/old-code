// prevent multiple inclusions
#pragma once

////////////////////
//	warnings
#pragma warning(disable : 4786) //"identifier was truncated to '255' characters
                                // in the debug information"

// NB: force warnings for unused arguments and local parameters
#pragma warning(1 : 4100 4101)

#include "AmortMidatADIPrice.h"
#include "CPDCalib.h"
#include "srt_h_all.h"

typedef struct _Calibrate_LamTS_DiagEuro_Common_ {
  const char *m_szYC;
  const char *m_szVC;
  Err (*m_get_cash_vol)(char *, double, double, double, int, char *, double *,
                        double *);
  const long m_lTheoEnd;
  const char *m_szRefRate;
  const char *m_szFreq;
  const char *m_szBasis;
  const long *const m_plEx_Begin;
  const long *const m_plEx_End;
  const char **m_pszSwapTenor;
  const double *m_pdSwapStrike;
  const char **m_pszCapTenor;
  const double *m_pdCapStrike;
  const cpd_diag_calib_param *m_pParam_Swaption;
  const cpd_diag_calib_param *m_pParam_Caplet;
  const diag_calib_lm_params *m_pCalibParam;
  const int m_nOnefequi;

} _Calib_LamTS_DiagEuro_Comm;

void _free_Calib_LamTS_DiagEuro_Comm(
    _Calib_LamTS_DiagEuro_Comm *pCalib_LamTS_DiagEuro);

void _init_Calib_LamTS_DiagEuro_Comm(
    long lToday, const char *szYC, const char *szVC,
    Err (*get_cash_vol)(char *, double, double, double, int, char *, double *,
                        double *),
    long lTheoEnd, const char *szRefRate, const char *szFreq,
    const char *szBasis, const long *plEx_Begin, const long *plEx_End,
    const int nOnefequi, const double *pdCoupon_Begin,
    const long *plFixStart_Begin, const long *plFixStart_End,
    const long *plFixEnd_Begin, const long *plFixPay_Begin,
    const long *plFltStart_Begin, const long *plFltStart_End,
    const long *plFltEnd_Begin, const long *plFltPay_Begin,
    const double *pdMargin_Begin, const double *pdSpread_Begin,
    const cpd_diag_calib_param *pParam_Swaption,
    const cpd_diag_calib_param *pParam_Caplet,
    const diag_calib_lm_params *pCalibParam,
    /// output
    _Calib_LamTS_DiagEuro_Comm *pCalib_LamTS_DiagEuro);

const char *
_Calibrate_TS_(int nNumPointT, int nNumPointX, const _GenMidAt *pOption_Begin,
               const _GenMidAt *pOption_End, const _PCQ_Seq *pNumeraire_Begin,
               const double *pOptionPrice_Begin,
               const _Calib_LamTS_DiagEuro_Comm *pCalib_LamTS_DiagEuro_Begin,
               // model
               const double *pdLamT_Begin, const double *pdLamT_End,
               const double *pdAlpha, const double *pdGamma,
               const double *pdRho,
               /// output
               double *pdLam_InitGuess_Begin // lambda ts
);

void Calibrate_MidAt_Preprocess(
    ///////Market
    long lToday, const char *szYC,
    /// MidAt General Specs
    const char *szRefRate, // Refrate
    const char
        *szFreq, // frequency of the fixed leg -"M"      ,"Q"      ,"S" ,"A"
    const double *pdCoupon_Begin, // coupon
    // MidAt Exercise
    const long *plEx_Begin,
    const long *plEx_End, // exercise dates
    // MidAt Fixed Leg
    const long *plFixStart_Begin, const long *plFixStart_End,
    const long *plFixEnd_Begin, const long *plFixPay_Begin,
    // MidAt Funding/Floating leg
    const long *plFltStart_Begin, const long *plFltStart_End,
    const long *plFltEnd_Begin, const long *plFltPay_Begin,
    const double *pdMargin_Begin, const double *pdSpread_Begin,
    // results
    char **pszSwapTenor, double *pdSwapStrike, char **pszCapTenor,
    double *pdCapStrike);

const char *_Calibrate_MidAt(
    ///////Market
    const char *szYC, const char *szVC,
    Err (*get_cash_vol)(char *, double, double, double, int, char *, double *,
                        double *), //	functor to get cash vol from the market
    /// MidAt General Specs
    long lTheoEnd,
    const char *szRefRate, // Refrate
    const char
        *szFreq, // frequency of the fixed leg -"M"      ,"Q"      ,"S" ,"A"
    const char *szBasis, // basis of the fixed leg
    // MidAt Exercise
    const long *plEx_Begin,
    const long *plEx_End, // exercise dates
    // calib params
    const char **pszSwapTenor_Begin, const double *pdSwapStrike_Begin,
    const char **pszCapTenor_Begin, const double *pdCapStrike_Begin,
    int nCalibTauToCap, // 1 - calibrate Tau to cap      , 0 otherwise
    const double *pdAlpha, const double *pdGamma, const double *pdRho,
    const cpd_diag_calib_param *pParam_Swaption,
    const cpd_diag_calib_param *pParam_Caplet,
    const diag_calib_lm_params *pCalibParam, int nOnefequi,
    const double *pdLamT_Begin, // NB: might be changed!
    const double *pdLamT_End,   // NB: might be changed!
    double *pdLam_Begin,        // NB: might be changed!
    // results
    double **ppdSigT_Begin, double **ppdSigT_End, double **ppdSig);

const char *Calibrate_MidAt(
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
    // MidAt Exercise
    const long *plEx_Begin,
    const long *plEx_End, // exercise dates
    // MidAt Fixed Leg
    const long *plFixStart_Begin, const long *plFixStart_End,
    const long *plFixEnd_Begin, const long *plFixPay_Begin,
    // MidAt Funding/Floating leg
    const long *plFltStart_Begin, const long *plFltStart_End,
    const long *plFltEnd_Begin, const long *plFltPay_Begin,
    const double *pdMargin_Begin, const double *pdSpread_Begin,
    // term structure
    int nCalibTauToCap, // 1 - calibrate Tau to cap      , 0 otherwise
    const double *pdAlpha, const double *pdGamma, const double *pdRho,
    // calibration Param_Caplet
    const cpd_diag_calib_param *pParam_Swaption,
    const cpd_diag_calib_param *pParam_Caplet,
    const diag_calib_lm_params *pCalibParam, int nOnefequi,
    double *pdLamT_Begin, // NB: might be changed!
    double *pdLamT_End,   // NB: might be changed!
    double *pdLam_Begin,  // NB: might be changed!
    // results
    double **ppdSigT_Begin, double **ppdSigT_End, double **ppdSig_Begin);

const char *Calibrate_GenMidAt(
    /// Market
    long lToday, const char *szYC, const char *szVC,
    Err (*get_cash_vol)(char *, double, double, double, int, char *, double *,
                        double *),
    const char *szCorrel,
    Err (*get_correl)(char *, double, double, double, double *),
    // GenMidAt general info
    int eod_flag, //	EOD Flag 0: I      , 1: E
    long lTheoEnd,
    const char *szRefRate, // ref rate
    const char *szFreq,    // swap_freq      , //swap freq
    const char *szBasis,   //	swap basis
    // Exercise info
    const long *plEx_Begin, const long *plEx_End,
    // Fixed leg
    const long *plFixStart_Begin, const long *plFixStart_End,
    const long *plFixEnd_Begin, const long *plFixPay_Begin,
    const double *pdFixRate_Begin, const double *pdFixNotional_Begin,
    const double *pdFixFee_Begin,
    // Floating leg
    const long *plFltStart_Begin, const long *plFltStart_End,
    const long *plFltEnd_Begin, const long *plFltPay_Begin,
    const double *pdMargin_Begin, const double *pdSpread_Begin,
    const double *pdFltNotional_Begin,
    // Given term structure
    const double *pdLamT_Begin, const double *pdLamT_End,
    const double *pdLam_Begin, const double *pdAlpha, const double *pdGamma,
    const double *pdRho,
    //	Calib params
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
    // Calibration results
    double **ppdSigT_Begin, double **ppdSigT_End, double **ppdSig_Begin,
    // market prices
    int *pnExBool, double *ppdDiagPrice);

const char *Calibrate_TS_Diag_Euro_MidAt(
    // Grid
    int nNumPointT, int nNumPointX,
    // market
    const char *szYC, const char *szVC,
    Err (*get_cash_vol)(char *, double, double, double, int, char *, double *,
                        double *), //	functor to get cash vol from the market
    long lToday, const char *szFreq, const char *szBasis, const char *szRefRate,
    long lTheoEnd, int nOnefequi,
    // model
    const double *pdLambda_Begin, // lambda from the matrices
    const double *pdAlpha, const double *pdGamma, const double *pdRho,
    /// Ex
    const long *plEx_Begin, const long *plEx_End, const long *plExFee_Begin,
    const double *pdExFee_Begin,
    /// Ex Start
    const long *plExStart_Begin,
    //// fixed leg
    const long *plFixPay_Begin, const long *plFixPay_End,
    const long *plFixStart_Begin, const long *plFixEnd_Begin,
    const long *plFixFee_Begin, const double *pdFixCoverage_Begin,
    const double *pdFixNotional_Begin, const double *pdFixFee_Begin,
    //// floating leg
    const long *plFltPay_Begin, const long *plFltPay_End,
    const long *plFltStart_Begin, const long *plFltEnd_Begin,
    const double *pdFltCoverage_Begin, const double *pdFltNotional_Begin,
    const double *pdMargin_Begin, const double *pdSpread_Begin,
    // PayRec
    double dPayRec, const cpd_diag_calib_param *pParam_Swaption,
    const cpd_diag_calib_param *pParam_Caplet,
    const diag_calib_lm_params *pCalibParam,
    // results
    double *pdOptionPV_Begin, double *pdCoupon_Begin,
    double *pdLam_InitGuess_Begin // term struture returned
);
