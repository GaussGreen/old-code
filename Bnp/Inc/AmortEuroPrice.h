// prevent multiple inclusions
#pragma once

//////////////////////
//	warnings
#pragma warning(disable : 4786) //"identifier was truncated to '255' characters
                                // in the debug information"

// NB: force warnings for unused arguments and local parameters
#pragma warning(1 : 4100 4101)

#include "AmortMidatADIUtils.h"

double FltLeg(long lToday, const char *szYC, const long *plStart_Begin,
              const long *plStart_End, const long *plEnd_Begin,
              const long *plPay_Begin, const double *pdMargin_Begin,
              const double *pdSpread_Begin, const double *pdNotional_Begin,
              SrtBasisCode eBasis, double *pdLevel,
              double *pdLevelMultCoupon_Margin,
              double *pdLevelMultCoupon_Spread);

double FixLeg(long lToday, const char *szYC, const long *plStart_Begin,
              const long *plStart_End, const long *plEnd_Begin,
              const long *plPay_Begin, const double *pdCoupon_Begin,
              const double *pdNotional_Begin, SrtBasisCode eBasis,
              double *pdLevel, double *pdLevelMultCoupon);

double Swap_Rate(long lToday, const char *szYC, const long *plFixStart_Begin,
                 const long *plFixStart_End, const long *plFixEnd_Begin,
                 const long *plFixPay_Begin, SrtBasisCode eFixBasis,
                 const long *plFltStart_Begin, const long *plFltStart_End,
                 const long *plFltEnd_Begin, const long *plFltPay_Begin,
                 const double *pdMargin_Begin, const double *pdSpread_Begin,
                 SrtBasisCode eFltBasis);

const char *Price_DiagGenEuro(
    // Market
    long lToday, const char *szYC, const char *szVC,
    Err (*get_correl)(char *, double, double, double, double *),
    const char *szCorrel,
    /// GenMidat general Info
    const char *szRefRate, const char *szFreq, const char *szBasis,
    long lTheoEnd,
    // Exercise
    long lFirstEx,
    // Fix Leg
    const long *plFixStart_Begin, const long *plFixStart_End,
    const long *plFixEnd_Begin, const double *pdFixNotional_Begin,
    const double *pdFixRate_Begin, const double *pdFixFee_Begin,
    // Floating Leg
    const long *plFltStart_Begin, const long *plFltStart_End,
    const long *plFltEnd_Begin, const double *pdFltNotional_Begin,
    const double *pdMargin_Begin, const double *pdSpread_Begin,
    // Calibratino specs
    double dMinTime, double dMinInterval, int nUseVol,
    // Output
    int *pnExBool, double *pdDiagPrice);