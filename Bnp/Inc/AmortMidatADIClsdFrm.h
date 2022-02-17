// prevent multiple inclusions
#pragma once

////////////////////
//	warnings
#pragma warning(disable : 4786) //"identifier was truncated to '255' characters
                                //in the debug information"

// NB: force warnings for unused arguments and local parameters
#pragma warning(1 : 4100 4101)

#include "AmortMidatADIUtils.h"

const char *GenMidAt_ClsdFrm_1F(
    // market
    const char *szYC, long lToday,
    // model
    const double *pdSigT_Begin, const double *pdSigT_End,
    const double *pdSig_Begin, const double *pdLamT_Begin,
    const double *pdLamT_End, const double *pdLam_Begin,
    /// Call
    long lEx, double dCoupon,
    //// fixed leg
    const long *plFixPay_Begin, const long *plFixPay_End,
    const double *pdFixCoverage_Begin, const double *pdFixNotional_Begin,
    //// floating leg
    long lFltStart, const long *plFltPay_Begin, const long *plFltPay_End,
    const double *pdFltCoverage_Begin, const double *pdFltNotional_Begin,
    const double *pdMargin_Begin, const double *pdSpread_Begin,
    /// pay or receive
    double dPayRec,
    /// optional params
    const long *plExFeePay, const double *pdExFee, const long *plFixFeePay,
    const double *pdFixFee,
    /// result
    double *pdOptionPV);

const char *GenMidAt_ClsdFrm_2F(
    // market
    const char *szYC, long lToday,
    // model
    const double *pdSigT_Begin, const double *pdSigT_End,
    const double *pdSig_Begin, const double *pdLamT_Begin,
    const double *pdLamT_End, const double *pdLam_Begin, double dAlpha,
    double dGamma, double dRho,
    /// Call
    long lEx, double dCoupon,
    //// fixed leg
    const long *plFixPay_Begin, const long *plFixPay_End,
    const double *pdFixCoverage_Begin, const double *pdFixNotional_Begin,
    //// floating leg
    long lFltStart, const long *plFltPay_Begin, const long *plFltPay_End,
    const double *pdFltCoverage_Begin, const double *pdFltNotional_Begin,
    const double *pdMargin_Begin, const double *pdSpread_Begin,
    /// pay or receive
    double dPayRec, const long *plExFeePay, const double *pdExFee,
    const long *plFixFeePay, const double *pdFixFee,
    /// result
    double *pdOptionPV);

void Convert_LamToG_1F(const double *pdLamT_Begin, const double *pdLamT_End,
                       const double *pdLam_Begin, const double *pdGT_Begin,
                       const double *pdGT_End,
                       // results
                       double *pdG_Begin);

void Convert_LamToG_2F(const double *pdLamT_Begin, const double *pdLamT_End,
                       const double *pdLam_Begin, const double *pdGT_Begin,
                       const double *pdGT_End, double dGamma,
                       // results
                       double *pdG1_Begin, double *pdG2_Begin);

const char *
Convert_SigToZeta_1F(const double *pdSigT_Begin, const double *pdSigT_End,
                     const double *pdSig_Begin, const double *pdLamT_Begin,
                     const double *pdLamT_End, const double *pdLam_Begin,
                     const double *pdZetaT_Begin, const double *pdZetaT_End,
                     // results
                     double *pdZeta_Begin);

const char *
Convert_SigToZeta_2F(const double *pdSigT_Begin, const double *pdSigT_End,
                     const double *pdSig_Begin, const double *pdLamT_Begin,
                     const double *pdLamT_End, const double *pdLam_Begin,
                     const double *pdZetaT_Begin, const double *pdZetaT_End,
                     double dAlpha, double dGamma, double dRho,
                     // results
                     double *pdZeta1_Begin, double *pdZeta2_Begin,
                     double *pdZeta12_Begin);
