#ifndef EURO_AMORT_SWAPTION_h
#define EURO_AMORT_SWAPTION_h

// ------------------------------------------------------------------------------------------------------------------
// // EuroAmortSwaption.h
//
// Declarations for functions to price european amortizing swaptions
//
// ------------------------------------------------------------------------------------------------------------------
// // C files
#include "math.h"

// Files included to allow compilation
#include "swp_h_curve_struct.h"  // needed by srt_h_lgmprotos.h

// Utility files
#include "utCurrency.h"  // needed by utTypes.h
#include "utDates.h"
#include "utError.h"
#include "utNrutil.h"
#include "utTypes.h"

// Sort files
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"
#include "swp_h_df.h"
#include "swp_h_spread.h"
#include "swp_h_swaption.h"
#include "swp_h_vol.h"

// Functions defined elsewhere

// Function to price a european amortizing swaption
Err EuropeanAmortizingSwaption(
    long    today,
    long    StartDate,
    long    TheoEndDate,
    long    NumPeriod,
    long*   lvFixedStartDates,
    long*   lvFixedEndDates,
    long    NoticePeriod,
    double  Coupon,
    double* dvNotionals,  // 0,..,NumPeriod  ( dvNotionals[NumPeriod]=0 )
    char*   szYieldCurveName,
    char*   szFreq,
    char*   szBasis,
    char*   szPayRec,
    char*   szRefRate,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*),
    SrtDiffusionType srt_vol_type,
    double*          dvPayDates,
    double*          dvEquivalentStrikes,
    double*          dvEquivalentNotionals,
    double*          dPrice,
    double*          dFixedLeg,
    double*          dFloatLeg,
    double*          dSwapRate,
    SigKapTS**       lgmSigKapTSptr,
    double*          dLGMVol,
    double*          dvReplicatingSwaptions,
    double           dMargin);

Err EuropeanAmortizingSwaption2(
    long AsOfDate,     // today or calculation date
    long FixDate,      // common fix date of underlying co-initial swaption
    long TheoEndDate,  // last end date of underlying
    long NumPeriods,  // total number of european swaptions that amortizing swaption decomposes into
    long*   plAcrlStartDates,  // vector(1,NumPeriods),
    long*   plAcrlEndDates,    // vector(1,NumPeriods),
    long*   plPayDates,        // vector(1,NumPeriods),
    double* pFixCvgsAdjusted,  // vector(1,NumPeriods), the same as pFixCvgs
    double* pFixCvgs,
    double  dCoupon,
    double* pdFixNotionals,  // vector(1,NumPeriods),
    char*   pszYieldCurveName,
    char*   szVolCurveName, /*	vc */
    char*   pszFreq,
    char*   pszBasis,
    char*   pszPayRec,
    char*   pszRefRate,
    Err (*get_cash_vol)(
        char*   vol_curve_name,
        double  start_date,
        double  end_date,
        double  cash_strike,
        int     zero,
        char*   ref_rate_name,
        double* vol,
        double* power),
    /*								LGMErr  (*GetVol)( Date, Date, double, SRT_Boolean, double
       *), SrtDiffusionType srt_vol_type,*/
    // Output
    double* dvPayDates,  // caller must pass in a vector(0,NumPeriods-1); returned zero bond values
                         // to pay dates
    double* dvReplicatingStrikes,  // caller must pass in a vector(0,NumPeriods-1);returned strikes
                                   // of replicating swaptions
    double* dvReplicatingNotionals,     // caller must pass in a vector(0,NumPeriods-1); returned
                                        // notionals of replicating swaptions
    double*    dPrice,                  //  returned amortizing swaption price
    double*    dFixedPV,                // returned fixed leg pv
    double*    dFloatPV,                // returned floating leg pv
    double*    dSwapRate,               // returned amortizing swap rate
    SigKapTS** lgmSigKapTSPtrPtr,       // returned term structure of vol and mean reversion
    double*    dLGMVol,                 // returned lgm vol
    double*    dvReplicatingSwaptions,  // caller must pass in a vector(0,NumPeriods-1); returned
                                        // values of underlying replicating european swaptions

    // input
    double* pdMargins,       // vector(1,NumPeriods),
    double* pdCoupon,        // optional; vector(1,NumPeriods),
    double* pdFltNotionals,  // optional; vector(1,NumPeriods),
    double  dExerciseFee);

// function solves for Coupon given BondMult
// also see function XLZeroCouponSwaption(...)
char* ZeroCouponSwaption_ComputeCoupon(
    long    Today,
    char*   pcYieldCurveName,
    double  dBondMult,
    long    lNumPeriods,
    long*   plFixedDates,
    long    lNumFltPeriods,
    long*   plFltDates,
    double* pdMargins,
    double* pdSpread,
    double* pdFixCov,
    double* pdFltCov,
    double* pdCoupon);

void ComputeCoverages(
    long  Today,
    char* pcYCName,

    // fixed leg info
    double* pdFixedCvgs,   // 0, NumFixPeriods -1
    long*   plFixedDates,  // 0,NumFixPeriods
    long    NumFixPeriods,

    // flt leg info
    double* pdFltCvgs,     // 0, NumFltPeriods -1
    double* pdFltSpreads,  // 0, NumFltPeriods -1
    double* pdMargins,     // 0, NumFltPeriods -1
    // double *pdFltNotionals,// 0, NumFltPeriods -1
    long* plFltDates,  // 0, NumFltPeriods
    long  NumFltPeriods,

    // output
    double* pdFixedCvgs_adjusted,

    double* pdCoupon,
    double* pdFltAdj);

void ComputeCoverages2(
    long  Today,
    char* pcYCName,

    // fixed leg info
    double* pdFixedCvgs,   // 0, NumFixPeriods -1
    long*   plFixedDates,  // 0,NumFixPeriods
    long    NumFixPeriods,
    double* pdFixedNotionals,

    // flt leg info
    double* pdFltCvgs,     // 0, NumFltPeriods -1
    double* pdFltSpreads,  // 0, NumFltPeriods -1
    double* pdMargins,     // 0, NumFltPeriods -1
    // double *pdFltNotionals,// 0, NumFltPeriods -1
    long* plFltDates,  // 0, NumFltPeriods
    long  NumFltPeriods,

    // output
    double* pdFixedCvgs_adjusted,

    double* pdCoupon,
    double* pdFltAdj);

typedef struct _EASData
{
    long    m_Today;
    char*   m_pYieldCurve;
    long    m_lNumFltPeriods;
    long    m_NumPeriods;
    long    m_AsOfDate;
    long    m_FutureDate;
    long    m_AcrlStartDate;
    long*   m_plPayDates;  // vector(1,NumPeriods)
    char*   m_pszYieldCurveName;
    double* m_pdA;  // vector(0,NumPeriods)
    double  m_G_FutureDate;
    double  m_G_AcrlStartDate;
    double* m_pdG;  // vector(1,NumPeriods)
    double  m_zeta;
    // double *m_pdASFV;
    double* m_pdASFLevel;
    double* m_pdASFZeroBond;
    double* m_pdCoverages;
    double  m_dExerciseFee;
} EASData, *pEASData;

typedef struct _EASData_ZC_BondMult
{
    long    m_Today;
    char*   m_pYieldCurve;
    double* m_pdFixCov_adjusted;

    long m_lNumPeriods;
    long m_lNumFltPeriods;

    long* m_plFixedDates;

    long* m_plFltDates;

    double  m_BondMult;
    double* m_pdMargins;
    double* m_pdSpread;
    double* m_pdFixCov;
    double* m_pdFltCov;

} EASData_ZC_BondMult, *pEASData_ZC_BondMult;

typedef struct _EASData_ZC_Price
{
    char*               m_und_name;
    char*               m_mkt_name;
    int                 m_notperiod;
    char*               m_ref_name;
    char*               m_freq_name;
    char*               m_basis_name;
    double              m_start_date;
    double              m_end_date;
    double              m_notional;
    unsigned short int* m_num_margins;
    unsigned short int* m_should_be_1_margins;
    double*             m_margins;
    double              m_lpxCoupon;
    double              m_lpxBondMult;

} EASData_ZC_Price, *pEASData_ZC_Price;

#endif  // EURO_AMORT_SWAPTION_h