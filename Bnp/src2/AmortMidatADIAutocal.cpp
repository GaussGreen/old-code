#include "AmortMidatADIAutocal.h"

#include "AmortMidatADICalib.h"
#include "AmortMidatADIPrice.h"
#include "AmortMidatAutocal.h"
#include "swp_h_all.h"

const char* MidAt_AutoCal_Diag(
    ///////Market
    long        lToday,
    const char* szYC,
    const char* szVC,
    Err (*get_cash_vol)(
        char*, double, double, double, int, char*, double*, double*),  //	functor to get cash
                                                                       // vol from the market
    /// MidAt General Specs
    long          lTheoEnd,
    const char*   szRefRate,       // Refrate
    const char*   szFreq,          // frequency of the fixed leg -"M","Q","S","A"
    const char*   szBasis,         // basis of the fixed leg
    const double* pdCoupon_Begin,  // coupon
    const char*   szPayRec,        //"REC" or "PAY"
    // MidAt Exercise
    const long* plEx_Begin,  // exercise dates
    const long* plEx_End,
    const long* plExStart_Begin,  // exercise premium dats
    // MidAt Fixed Leg
    const long*   plFixPay_Begin,
    const long*   plFixPay_End,
    const long*   plFixStart_Begin,
    const long*   plFixEnd_Begin,
    const double* pdFixCvg_Begin,
    // MidAt Funding/Floating leg
    const long*   plFltPay_Begin,
    const long*   plFltPay_End,
    const long*   plFltStart_Begin,
    const long*   plFltEnd_Begin,
    const double* pdFltCvg_Begin,
    const double* pdMargin_Begin,
    const double* pdSpread_Begin,
    // term structure
    int           nCalibTauToCap,  // 1 - calibrate Tau to cap, 0 otherwise
    const double* pdLam_Begin,
    const double* pdAlpha,
    const double* pdGamma,
    const double* pdRho,
    // calibration params
    int nOnefequi,
    // grid
    int nNumPointT,
    int nNumPointX,
    /// parameters needed by diagcalib
    const cpd_diag_calib_param* pParam_Swaption,
    const cpd_diag_calib_param* pParam_Caplet,
    const diag_calib_lm_params* pCalibParam,
    /// results,
    double* pdOptionPV_Begin)
{
    const int    nSize_Ex  = plEx_End - plEx_Begin;
    const int    nSize_Fix = plFixPay_End - plFixPay_Begin;
    const int    nSize_Flt = plFltPay_End - plFltPay_Begin;
    const double dLamT     = 0.1;  // 1m

    for (; plEx_Begin < plEx_End;
         ++plEx_Begin, ++plExStart_Begin, ++pdOptionPV_Begin, ++pdCoupon_Begin, ++pdLam_Begin)
    {
        const char* szErr = "GenMidAt_AutoCal_Diag(...): internal dates error!";
        // locate indeces
        const long* plFixPay_nI_Begin = (const long*)find_min_greater_than(
            plExStart_Begin, plFixPay_Begin, nSize_Fix, sizeof(long), lless_v);
        const long* plFltPay_nI_Begin = (const long*)find_min_greater_than(
            plExStart_Begin, plFltPay_Begin, nSize_Flt, sizeof(long), lless_v);

        const int nOffSet_Fix = plFixPay_nI_Begin - plFixPay_Begin;
        const int nOffSet_Flt = plFltPay_nI_Begin - plFltPay_Begin;

        // check that indeces are located correctly
        if (nOffSet_Fix == nSize_Fix)
            return szErr;
        if (nOffSet_Flt == nSize_Flt)
            return szErr;

        szErr = MidAt_AutoCal(
            lToday,
            szYC,
            szVC,
            get_cash_vol,
            lTheoEnd,
            szRefRate,       // Refrate
            szFreq,          // frequency of the fixed leg -"M","Q","S","A"
            szBasis,         // basis of the fixed leg
            pdCoupon_Begin,  // coupon
            szPayRec,        //"REC" or "PAY"
            plEx_Begin,      // exercise dates
            plEx_End,
            plExStart_Begin,  // exercise premium dats
            // MidAt Fixed Leg
            plFixPay_Begin + nOffSet_Fix,
            plFixPay_End,
            plFixStart_Begin + nOffSet_Fix,
            plFixEnd_Begin + nOffSet_Fix,
            pdFixCvg_Begin + nOffSet_Fix,
            // MidAt Funding/Floating leg
            plFltPay_Begin + nOffSet_Flt,
            plFltPay_End,
            plFltStart_Begin + nOffSet_Flt,
            plFltEnd_Begin + nOffSet_Flt,
            pdFltCvg_Begin + nOffSet_Flt,
            pdMargin_Begin + nOffSet_Flt,
            pdSpread_Begin + nOffSet_Flt,
            // term structure
            nCalibTauToCap,  // 1 - calibrate Tau to cap, 0 otherwise
            (double*)&dLamT,
            1 + (double*)&dLamT,
            (double*)pdLam_Begin,
            pdAlpha,
            pdGamma,
            pdRho,
            // calibration params
            nOnefequi,
            // grid
            nNumPointT,
            nNumPointX,
            /// parameters needed by diagcalib
            pParam_Swaption,
            pParam_Caplet,
            pCalibParam,
            /// results,
            pdOptionPV_Begin,
            0,  // double *pdFixLegPV,
            0,  // double *pdFltLegPV,
            0,  // double *pdExProb,
            0,  // double *pdExBoundary,
            0   // i0nt nUse_Backward
        );

        if (szErr)
            return szErr;
    }

    return 0;
}

const char* Euro_AutoCal_Diag(
    ///////Market
    long        lToday,
    const char* szYC,
    const char* szVC,
    Err (*get_cash_vol)(
        char*, double, double, double, int, char*, double*, double*),  //	functor to get cash
                                                                       // vol from the market
    /// MidAt General Specs
    long          lTheoEnd,
    const char*   szRefRate,       // Refrate
    const char*   szFreq,          // frequency of the fixed leg -"M","Q","S","A"
    const char*   szBasis,         // basis of the fixed leg
    const double* pdCoupon_Begin,  // coupon
    const char*   szPayRec,        //"REC" or "PAY"
    // MidAt Exercise
    const long* plEx_Begin,  // exercise dates
    const long* plEx_End,
    const long* plExStart_Begin,  // exercise premium dats
    // MidAt Fixed Leg
    const long*   plFixPay_Begin,
    const long*   plFixPay_End,
    const long*   plFixStart_Begin,
    const long*   plFixEnd_Begin,
    const double* pdFixCvg_Begin,
    // MidAt Funding/Floating leg
    const long*   plFltPay_Begin,
    const long*   plFltPay_End,
    const long*   plFltStart_Begin,
    const long*   plFltEnd_Begin,
    const double* pdFltCvg_Begin,
    const double* pdMargin_Begin,
    const double* pdSpread_Begin,
    // term structure
    int           nCalibTauToCap,  // 1 - calibrate Tau to cap, 0 otherwise
    const double* pdLam_Begin,
    const double* pdAlpha,
    const double* pdGamma,
    const double* pdRho,
    // calibration params
    int nOnefequi,
    // grid
    int nNumPointT,
    int nNumPointX,
    /// parameters needed by diagcalib
    const cpd_diag_calib_param* pParam_Swaption,
    const cpd_diag_calib_param* pParam_Caplet,
    const diag_calib_lm_params* pCalibParam,
    /// results,
    double* pdOptionPV_Begin)
{
    const int    nSize_Ex  = plEx_End - plEx_Begin;
    const int    nSize_Fix = plFixPay_End - plFixPay_Begin;
    const int    nSize_Flt = plFltPay_End - plFltPay_Begin;
    const double dLamT     = 0.1;  // 1m

    for (; plEx_Begin < plEx_End;
         ++plEx_Begin, ++plExStart_Begin, ++pdOptionPV_Begin, ++pdCoupon_Begin, ++pdLam_Begin)
    {
        const char* szErr = "GenMidAt_AutoCal_Diag(...): internal dates error!";
        // locate indeces
        const long* plFixPay_nI_Begin = (const long*)find_min_greater_than(
            plExStart_Begin, plFixPay_Begin, nSize_Fix, sizeof(long), lless_v);
        const long* plFltPay_nI_Begin = (const long*)find_min_greater_than(
            plExStart_Begin, plFltPay_Begin, nSize_Flt, sizeof(long), lless_v);

        const int nOffSet_Fix = plFixPay_nI_Begin - plFixPay_Begin;
        const int nOffSet_Flt = plFltPay_nI_Begin - plFltPay_Begin;

        // check that indeces are located correctly
        if (nOffSet_Fix == nSize_Fix)
            return szErr;
        if (nOffSet_Flt == nSize_Flt)
            return szErr;

        szErr = MidAt_AutoCal(
            lToday,
            szYC,
            szVC,
            get_cash_vol,
            lTheoEnd,
            szRefRate,       // Refrate
            szFreq,          // frequency of the fixed leg -"M","Q","S","A"
            szBasis,         // basis of the fixed leg
            pdCoupon_Begin,  // coupon
            szPayRec,        //"REC" or "PAY"
            plEx_Begin,      // exercise dates
            1 + plEx_Begin,
            plExStart_Begin,  // exercise premium dats
            // MidAt Fixed Leg
            plFixPay_Begin + nOffSet_Fix,
            plFixPay_End,
            plFixStart_Begin + nOffSet_Fix,
            plFixEnd_Begin + nOffSet_Fix,
            pdFixCvg_Begin + nOffSet_Fix,
            // MidAt Funding/Floating leg
            plFltPay_Begin + nOffSet_Flt,
            plFltPay_End,
            plFltStart_Begin + nOffSet_Flt,
            plFltEnd_Begin + nOffSet_Flt,
            pdFltCvg_Begin + nOffSet_Flt,
            pdMargin_Begin + nOffSet_Flt,
            pdSpread_Begin + nOffSet_Flt,
            // term structure
            nCalibTauToCap,  // 1 - calibrate Tau to cap, 0 otherwise
            (double*)&dLamT,
            1 + (double*)&dLamT,
            (double*)pdLam_Begin,
            pdAlpha,
            pdGamma,
            pdRho,
            // calibration params
            nOnefequi,
            // grid
            nNumPointT,
            nNumPointX,
            /// parameters needed by diagcalib
            pParam_Swaption,
            pParam_Caplet,
            pCalibParam,
            /// results,
            pdOptionPV_Begin,
            0,  // double *pdFixLegPV,
            0,  // double *pdFltLegPV,
            0,  // double *pdExProb,
            0,  // double *pdExBoundary,
            0   // i0nt nUse_Backward
        );

        if (szErr)
            return szErr;
    }

    return 0;
}

const char* GenMidAt_AutoCal(
    /// Market
    long        lToday,
    const char* szYC,
    const char* szVC,
    Err (*get_cash_vol)(char*, double, double, double, int, char*, double*, double*),
    const char* szCorrel,
    Err (*get_correl)(char*, double, double, double, double*),
    /// MidAt Specs
    const char* szRefRate,
    int         eod_flag,  //	EOD Flag 0: I, 1: E
    long        lTheoEnd,
    const char* szPayRec,  //"REC" or "PAY"
    // MidAt Exercise
    const long*   plEx_Begin,  // exercise dates
    const long*   plEx_End,
    const long*   plExStart_Begin,  // exercise premium dats
    const double* pdExFee_Begin,    // exercise fee
    // MidAt Fixed Leg
    const char*   szFreq,
    const char*   szBasis,
    const long*   plFixStart_Begin,
    const long*   plFixStart_End,
    const long*   plFixEnd_Begin,
    const double* pdFixCvg_Begin,
    const long*   plFixPay_Begin,
    const double* pdFixNotional_Begin,
    const double* pdFixCoupon_Begin,  // coupon
    const double* pdFixFee_Begin,
    // MidAt Funding/Floating leg
    const long*   plFltStart_Begin,
    const long*   plFltStart_End,
    const long*   plFltEnd_Begin,
    const double* pdFltCvg_Begin,
    const long*   plFltPay_Begin,
    const double* pdMargin_Begin,
    const double* pdSpread_Begin,
    const double* pdFltNotional_Begin,
    // Given Term Structure
    const long*   plLamDate_Begin,
    const long*   plLamDate_End,
    const double* pdLam_Begin,
    const double* pdAlpha,
    const double* pdGamma,
    const double* pdRho,
    /// Calibration Params
    const char* default_ref,         //	ref rate
    const char* default_swap_freq,   // swap freq
    const char* default_swap_basis,  // swap basis
    double      dMinTime,
    double      dMinInterval,
    long        lNotice,
    double      dMaxStdShort,
    int         nFixLam,  // 0: calib lambda to cap, 1: fix lambda calib	to diagonal
    int n1FEqui,    // 1: 2F lambda will calibrate to the cap priced within calibrated 1F	with the
                    // given lambda
    int nSkipLast,  // If 1, the last option is disregardedand the forward volatility is flat from
                    // option	n-1
    int    nUseJump,
    double dMaxVarJump,
    int    nStrikeType,
    int    nEuroModel,
    // grid
    int nNumPointT,
    int nNumPointX,
    /// Results,
    // NB: Set the sig related ptrs to 0 to calibrate sigma; else use the sigma passed in !!!
    long**   pplSigDate_Begin,
    long**   pplSigDate_End,
    double** ppdSig_Begin,
    double*  pdOptionPV,
    double*  pdFixLegPV,
    double*  pdFltLegPV,
    double*  pdExProb,
    double*  pdExBoundary,
    int      nUse_Backward,
    int*     pnBool_CalibInst,
    double*  pdPrice_CalibInst)
{
    const char* szErr = 0;

    // sigma can be calibrated
    int     nSize_Sig    = *pplSigDate_End - *pplSigDate_Begin;
    double *pdSigT_Begin = (double*)malloc(nSize_Sig * sizeof(double)),
           *pdSigT_End   = pdSigT_Begin + nSize_Sig;

    // lambda fixed
    const int nSize_Lam    = plLamDate_End - plLamDate_Begin;
    double *  pdLamT_Begin = (double*)malloc(nSize_Lam * sizeof(double)),
           *pdLamT_End     = pdLamT_Begin + nSize_Lam;

    _fill_tenor(lToday, plLamDate_Begin, plLamDate_End, pdLamT_Begin);
    _fill_tenor(lToday, *pplSigDate_Begin, *pplSigDate_End, pdSigT_Begin);

    // Calibrate if any of the 3 sigma related ptrs is 0
    // else use the sigma ts passed in
    if (!(*pplSigDate_Begin) || !(*pplSigDate_End) || !(*ppdSig_Begin))
    {
        szErr = Calibrate_GenMidAt(
            lToday,
            szYC,
            szVC,
            get_cash_vol,
            szCorrel,
            get_correl,
            eod_flag,
            lTheoEnd,
            szRefRate,
            szFreq,
            szBasis,
            plEx_Begin,
            plEx_End,
            plFixStart_Begin,
            plFixStart_End,
            plFixEnd_Begin,
            plFixPay_Begin,
            pdFixCoupon_Begin,
            pdFixNotional_Begin,
            pdFixFee_Begin,
            plFltStart_Begin,
            plFltStart_End,
            plFltEnd_Begin,
            plFltPay_Begin,
            pdMargin_Begin,
            pdSpread_Begin,
            pdFltNotional_Begin,
            pdLamT_Begin,
            pdLamT_Begin + nSize_Lam,
            pdLam_Begin,
            pdAlpha,
            pdGamma,
            pdRho,
            default_ref,
            default_swap_freq,
            default_swap_basis,
            dMinTime,
            dMinInterval,
            lNotice,
            dMaxStdShort,
            nFixLam,
            n1FEqui,
            nSkipLast,
            nUseJump,
            dMaxVarJump,
            nStrikeType,
            nEuroModel,
            &pdSigT_Begin,
            &pdSigT_End,
            ppdSig_Begin,
            pnBool_CalibInst,
            pdPrice_CalibInst);

        nSize_Sig = pdSigT_End - pdSigT_Begin;
    }

#ifdef _DEBUG
    _print_vector(pdSigT_Begin, nSize_Sig);
    _print_vector(*ppdSig_Begin, nSize_Sig);
    _print_vector(pdLamT_Begin, nSize_Lam);
    _print_vector(pdLam_Begin, nSize_Lam);
#endif

    /// Price ....
    szErr = _Price_GenMidAt(
        nNumPointT,
        nNumPointX,
        szYC,
        lToday,
        pdSigT_Begin,
        pdSigT_Begin + nSize_Sig,
        *ppdSig_Begin,
        pdLamT_Begin,
        pdLamT_Begin + nSize_Lam,
        pdLam_Begin,
        pdAlpha,
        pdGamma,
        pdRho,
        plEx_Begin,
        plEx_End,
        plExStart_Begin,
        pdExFee_Begin,
        plExStart_Begin,
        plFixPay_Begin,
        plFixPay_Begin + (plFixStart_End - plFixStart_Begin),
        plFixStart_Begin,
        plFixEnd_Begin,
        plFixPay_Begin,  //&vec_lFixFeePay[0],
        pdFixCvg_Begin,
        pdFixNotional_Begin,
        pdFixCoupon_Begin,
        pdFixFee_Begin,
        plFltPay_Begin,
        plFltPay_Begin + (plFltStart_End - plFltStart_Begin),
        plFltStart_Begin,
        plFltEnd_Begin,
        pdFltCvg_Begin,
        pdFltNotional_Begin,
        pdMargin_Begin,
        pdSpread_Begin,
        _PayRec(szPayRec),
        pdFixLegPV,
        pdFltLegPV,
        pdOptionPV,
        pdExProb,
        pdExBoundary,
        nUse_Backward);

    if (!(*pplSigDate_Begin) || !(*pplSigDate_End) || !(*ppdSig_Begin))
    {
        *pplSigDate_Begin = (long*)calloc((pdSigT_End - pdSigT_Begin), sizeof(long));
        *pplSigDate_End   = *pplSigDate_Begin + (pdSigT_End - pdSigT_Begin);
        _fill_date(lToday, pdSigT_Begin, pdSigT_End, *pplSigDate_Begin);
        free(pdSigT_Begin);
    }

    return szErr;
}

/// the part of MidAt_AutoCal that does not
/// have the calibration pre-processing
const char* _MidAt_AutoCal(
    ///////Market
    long        lToday,
    const char* szYC,
    const char* szVC,
    Err (*get_cash_vol)(
        char*, double, double, double, int, char*, double*, double*),  //	functor to get cash
                                                                       // vol from the market
    /// MidAt General Specs
    long          lTheoEnd,
    const char*   szRefRate,       // Refrate
    const char*   szFreq,          // frequency of the fixed leg -"M","Q","S","A"
    const char*   szBasis,         // basis of the fixed leg
    const double* pdCoupon_Begin,  // coupon
    const char*   szPayRec,        //"REC" or "PAY"
    // MidAt Exercise
    const long* plEx_Begin,  // exercise dates
    const long* plEx_End,
    const long* plExStart_Begin,  // exercise premium dats
    // MidAt Fixed Leg
    const long*   plFixPay_Begin,
    const long*   plFixPay_End,
    const long*   plFixStart_Begin,
    const long*   plFixEnd_Begin,
    const double* pdFixCvg_Begin,
    // MidAt Funding/Floating leg
    const long*   plFltPay_Begin,
    const long*   plFltPay_End,
    const long*   plFltStart_Begin,
    const long*   plFltEnd_Begin,
    const double* pdFltCvg_Begin,
    const double* pdMargin_Begin,
    const double* pdSpread_Begin,
    // term structure
    int           nCalibTauToCap,  // 1 - calibrate Tau to cap, 0 otherwise
    double*       pdLamT_Begin,
    double*       pdLamT_End,
    double*       pdLam_Begin,
    const double* pdAlpha,
    const double* pdGamma,
    const double* pdRho,
    // calibration params
    const char**                pszSwapTenor_Begin,
    const double*               pdSwapStrike_Begin,
    const char**                pszCapTenor_Begin,
    const double*               pdCapStrike_Begin,
    const cpd_diag_calib_param* pParam_Swaption,
    const cpd_diag_calib_param* pParam_Caplet,
    const diag_calib_lm_params* pCalibParam,
    int                         nOnefequi,
    // grid
    int nNumPointT,
    int nNumPointX,
    /// results,
    double* pdOptionPV,
    double* pdFixLegPV,
    double* pdFltLegPV,
    double* pdExProb,
    double* pdExBoundary,
    int     nUse_Backward)
{
    const int nSize_Ex  = plEx_End - plEx_Begin;
    const int nSize_Fix = plFixPay_End - plFixPay_Begin;
    const int nSize_Flt = plFltPay_End - plFltPay_Begin;

    const double dFee      = 0.;
    const double dNotional = 1.;

    const size_t  size = sizeof(double);
    const double* pdExFee =
        (const double*)memset((double*)malloc(nSize_Ex * size), dFee, nSize_Ex * size);

    const double* pdFixFee =
        (const double*)memset((double*)malloc(nSize_Ex * size), dFee, nSize_Ex * size);
    // const double *pdCoupon_Begin = 	_memset(&dCoupon,_alloca(nSize_Fix*size),nSize_Fix,size);
    const double* pdFixNotional =
        (const double*)memset((double*)malloc(nSize_Ex * size), 1, nSize_Ex * size);

    const double* pdFltNotional =
        (const double*)memset((double*)malloc(nSize_Ex * size), 1, nSize_Ex * size);

    double *pdSigT_Begin = 0, *pdSigT_End = 0, *pdSig_Begin = 0;

    // Step 1. Calibrate Midat
    const char* szErr = _Calibrate_MidAt(
        szYC,
        szVC,
        get_cash_vol,
        lTheoEnd,
        szRefRate,
        szFreq,
        szBasis,
        plEx_Begin,
        plEx_End,
        pszSwapTenor_Begin,
        pdSwapStrike_Begin,
        pszCapTenor_Begin,
        pdCapStrike_Begin,
        nCalibTauToCap,
        pdAlpha,
        pdGamma,
        pdRho,
        pParam_Swaption,
        pParam_Caplet,
        pCalibParam,
        nOnefequi,
        pdLamT_Begin,
        pdLamT_End,
        pdLam_Begin,
        &pdSigT_Begin,
        &pdSigT_End,
        &pdSig_Begin);

    if (szErr)
        return szErr;

    // Step 2: pricing
    return _Price_GenMidAt(
        nNumPointT,
        nNumPointX,
        szYC,
        lToday,
        pdSigT_Begin,
        pdSigT_End,
        pdSig_Begin,
        pdLamT_Begin,
        pdLamT_End,
        pdLam_Begin,
        pdAlpha,
        pdGamma,
        pdRho,
        plEx_Begin,
        plEx_End,
        plExStart_Begin,
        pdExFee,  //&vec_dExFee[0],
        plExStart_Begin,
        plFixPay_Begin,
        plFixPay_End,
        plFixStart_Begin,
        plFixEnd_Begin,
        plFixPay_Begin,  //&vec_lFixFeePay[0],
        pdFixCvg_Begin,
        pdFixNotional,   //&vec_dFixNotional[0],
        pdCoupon_Begin,  //&vec_dExCoupon[0],
        pdFixFee,        //&vec_dFixFee[0],
        plFltPay_Begin,
        plFltPay_End,
        plFltStart_Begin,
        plFltEnd_Begin,
        pdFltCvg_Begin,
        pdFltNotional,  //&vec_dFltNotional[0],
        pdMargin_Begin,
        pdSpread_Begin,
        _PayRec(szPayRec),
        pdFixLegPV,
        pdFltLegPV,
        pdOptionPV,
        pdExProb,
        pdExBoundary,
        nUse_Backward);
}

const char* MIDAT_AUTOCAL(const Comm_MidAtAutocal* pComm)
{
    return MidAt_AutoCal(
        pComm->m_lToday,
        pComm->m_szYC,
        pComm->m_szVC,
        pComm->m_get_cash_vol,
        pComm->m_lTheoEnd,
        pComm->m_szRefRate,       // Refrate
        pComm->m_szFreq,          // frequency of the fixed leg -"M","Q","S","A"
        pComm->m_szBasis,         // basis of the fixed leg
        pComm->m_pdCoupon_Begin,  // coupon
        pComm->m_szPayRec,        //"REC" or "PAY"
        pComm->m_plEx_Begin,      // exercise dates
        pComm->m_plEx_End,
        pComm->m_plExStart_Begin,  // exercise premium dats
        pComm->m_plFixPay_Begin,
        pComm->m_plFixPay_End,
        pComm->m_plFixStart_Begin,
        pComm->m_plFixEnd_Begin,
        pComm->m_pdFixCvg_Begin,
        pComm->m_plFltPay_Begin,
        pComm->m_plFltPay_End,
        pComm->m_plFltStart_Begin,
        pComm->m_plFltEnd_Begin,
        pComm->m_pdFltCvg_Begin,
        pComm->m_pdMargin_Begin,
        pComm->m_pdSpread_Begin,
        pComm->m_nCalibTauToCap,  // 1 - calibrate Tau to cap, 0 otherwise
        (double*)pComm->m_pdLamT_Begin,
        (double*)pComm->m_pdLamT_End,
        (double*)pComm->m_pdLam_Begin,
        pComm->m_pdAlpha,
        pComm->m_pdGamma,
        pComm->m_pdRho,
        pComm->m_nOnefequi,
        pComm->m_nNumPointT,
        pComm->m_nNumPointX,
        pComm->m_pParam_Swaption,
        pComm->m_pParam_Caplet,
        pComm->m_pCalibParam,
        pComm->m_pdOptionPV,
        pComm->m_pdFixLegPV,
        pComm->m_pdFltLegPV,
        0,
        0,
        0);
}

const char* MidAt_AutoCal(
    ///////Market
    long        lToday,
    const char* szYC,
    const char* szVC,
    Err (*get_cash_vol)(
        char*, double, double, double, int, char*, double*, double*),  //	functor to get cash
                                                                       // vol from the market
    /// MidAt General Specs
    long          lTheoEnd,
    const char*   szRefRate,       // Refrate
    const char*   szFreq,          // frequency of the fixed leg -"M","Q","S","A"
    const char*   szBasis,         // basis of the fixed leg
    const double* pdCoupon_Begin,  // coupon
    const char*   szPayRec,        //"REC" or "PAY"
    // MidAt Exercise
    const long* plEx_Begin,  // exercise dates
    const long* plEx_End,
    const long* plExStart_Begin,  // exercise premium dats
    // MidAt Fixed Leg
    const long*   plFixPay_Begin,
    const long*   plFixPay_End,
    const long*   plFixStart_Begin,
    const long*   plFixEnd_Begin,
    const double* pdFixCvg_Begin,
    // MidAt Funding/Floating leg
    const long*   plFltPay_Begin,
    const long*   plFltPay_End,
    const long*   plFltStart_Begin,
    const long*   plFltEnd_Begin,
    const double* pdFltCvg_Begin,
    const double* pdMargin_Begin,
    const double* pdSpread_Begin,
    // term structure
    int           nCalibTauToCap,  // 1 - calibrate Tau to cap, 0 otherwise
    double*       pdLamT_Begin,
    double*       pdLamT_End,
    double*       pdLam_Begin,
    const double* pdAlpha,
    const double* pdGamma,
    const double* pdRho,
    // calibration params
    int nOnefequi,
    // grid
    int nNumPointT,
    int nNumPointX,
    /// parameters needed by diagcalib
    const cpd_diag_calib_param* pParam_Swaption,
    const cpd_diag_calib_param* pParam_Caplet,
    const diag_calib_lm_params* pCalibParam,
    /// results,
    double* pdOptionPV,
    double* pdFixLegPV,
    double* pdFltLegPV,
    double* pdExProb,
    double* pdExBoundary,
    int     nUse_Backward)
{
    const int nSize_Ex  = plEx_End - plEx_Begin;
    const int nSize_Fix = plFixPay_End - plFixPay_Begin;
    const int nSize_Flt = plFltPay_End - plFltPay_Begin;

    // const double dFee = 0.;
    // const double dNotional = 1.;

    // double *pdSigT_Begin=0,*pdSigT_End=0,*pdSig_Begin=0;

    // Swaption tenor
    char**  pszSwapTenor = (char**)malloc(nSize_Ex * sizeof(char*));
    double* pdSwapStrike = (double*)malloc(nSize_Ex * sizeof(double));

    // cap tenor
    char**  pszCapTenor = (char**)malloc(nSize_Ex * sizeof(char*));
    double* pdCapStrike = (double*)malloc(nSize_Ex * sizeof(double));

    // delegate to prepare for the call to diagcalibdlm
    Calibrate_MidAt_Preprocess(
        lToday,
        szYC,
        szRefRate,       // Refrate
        szFreq,          // frequency of the fixed leg -"M","Q","S","A"
        pdCoupon_Begin,  // coupon
        plEx_Begin,
        plEx_End,
        plFixStart_Begin,
        plFixStart_Begin + nSize_Fix,
        plFixEnd_Begin,
        plFixPay_Begin,
        plFltStart_Begin,
        plFltStart_Begin + nSize_Flt,
        plFltEnd_Begin,
        plFltPay_Begin,
        pdMargin_Begin,
        pdSpread_Begin,
        pszSwapTenor,
        pdSwapStrike,
        pszCapTenor,
        pdCapStrike);

    return _MidAt_AutoCal(
        lToday,
        szYC,
        szVC,
        get_cash_vol,
        lTheoEnd,
        szRefRate,
        szFreq,
        szBasis,
        pdCoupon_Begin,
        szPayRec,
        plEx_Begin,
        plEx_End,
        plExStart_Begin,
        plFixPay_Begin,
        plFixPay_End,
        plFixStart_Begin,
        plFixEnd_Begin,
        pdFixCvg_Begin,
        plFltPay_Begin,
        plFltPay_End,
        plFltStart_Begin,
        plFltEnd_Begin,
        pdFltCvg_Begin,
        pdMargin_Begin,
        pdSpread_Begin,
        nCalibTauToCap,
        pdLamT_Begin,
        pdLamT_End,
        pdLam_Begin,
        pdAlpha,
        pdGamma,
        pdRho,
        (const char**)pszSwapTenor,
        pdSwapStrike,
        (const char**)pszCapTenor,
        pdCapStrike,
        pParam_Swaption,
        pParam_Caplet,
        pCalibParam,
        nOnefequi,
        nNumPointT,
        nNumPointX,
        pdOptionPV,
        pdFixLegPV,
        pdFltLegPV,
        pdExProb,
        pdExBoundary,
        nUse_Backward);
}
