#include "AmortEuroPrice.h"

#include "swp_h_amortswaption.h"
#include "swp_h_convslidingcorrel.h"

static void _FixPV(
    long          lToday,
    const char*   szYC,
    const long*   plStart_Begin,
    const long*   plStart_End,
    const long*   plEnd_Begin,
    const long*   plPay_Begin,
    const double* pdCoupon_Begin,
    const double* pdNotional_Begin,
    SrtBasisCode  eBasis,
    // results
    double* pdPV,
    double* pdLevel,
    double* pdLevel_Mult_Coupon)
{
    _ASSERTE(pdPV && pdLevel && pdLevel_Mult_Coupon);

    *pdPV                = 0.;
    *pdLevel             = 0.;
    *pdLevel_Mult_Coupon = 0.;

    for (; plStart_Begin < plStart_End;
         ++plStart_Begin, ++plEnd_Begin, ++plPay_Begin, ++pdCoupon_Begin)
    {
        double dCvg = 0., dDF = 0., dInc = 0.;
        _ASSERTE(plStart_Begin && plEnd_Begin && pdCoupon_Begin && pdNotional_Begin);
        dCvg = coverage(*plStart_Begin, *plEnd_Begin, eBasis);
        _ASSERTE(*plPay_Begin >= lToday);
        dDF  = swp_f_df(lToday, *plPay_Begin, szYC);
        dInc = dCvg * dDF;
        // level
        *pdLevel += dInc;
        // level multiplied by coupon
        dInc *= (*pdCoupon_Begin);
        *pdLevel_Mult_Coupon += dInc;
        // level multiplied by coupon and notional
        dInc *= (*pdNotional_Begin);
        *pdPV += dInc;
    }
}

double Swap_Rate(
    long          lToday,
    const char*   szYC,
    const long*   plFixStart_Begin,
    const long*   plFixStart_End,
    const long*   plFixEnd_Begin,
    const long*   plFixPay_Begin,
    SrtBasisCode  eFixBasis,
    const long*   plFltStart_Begin,
    const long*   plFltStart_End,
    const long*   plFltEnd_Begin,
    const long*   plFltPay_Begin,
    const double* pdMargin_Begin,
    const double* pdSpread_Begin,
    SrtBasisCode  eFltBasis)
{
    double dFltPV, dLevel_Flt, dLevel_Mult_Margin, dLevel_Mult_Spread;
    double dFixPV, dLevel_Fix, dLevel_Mult_Coupon;

    const int nSize_Fix = plFixStart_End - plFixStart_Begin;
    const int nSize_Flt = plFltStart_End - plFltStart_Begin;

    double*      pdFixNotional_Begin = (double*)malloc(nSize_Fix * sizeof(double));
    double*      pdFltNotional_Begin = (double*)malloc(nSize_Flt * sizeof(double));
    const double dEqui_Notional      = 1.;

    _ASSERTE(nSize_Fix >= 0 && nSize_Flt >= 0);
    _memset(&dEqui_Notional, pdFixNotional_Begin, nSize_Fix, sizeof(dEqui_Notional));
    _memset(&dEqui_Notional, pdFltNotional_Begin, nSize_Flt, sizeof(dEqui_Notional));

    dFltPV = FltLeg(
        lToday,
        szYC,
        plFltStart_Begin,
        plFltStart_End,
        plFltEnd_Begin,
        plFltPay_Begin,
        pdMargin_Begin,
        pdSpread_Begin,
        pdFltNotional_Begin,
        eFltBasis,
        &dLevel_Flt,
        &dLevel_Mult_Margin,
        &dLevel_Mult_Spread);

    dFixPV = FixLeg(
        lToday,
        szYC,
        plFixStart_Begin,
        plFixStart_End,
        plFixEnd_Begin,
        plFixPay_Begin,
        (double*)memset(malloc(nSize_Fix * sizeof(double)), 0, nSize_Fix * sizeof(double)),
        pdFixNotional_Begin,
        eFixBasis,
        &dLevel_Fix,
        &dLevel_Mult_Coupon);

    _ASSERTE(dLevel_Fix != 0.);
    return dFltPV / dLevel_Fix;
}

double FltLeg(
    long          lToday,
    const char*   szYC,
    const long*   plStart_Begin,
    const long*   plStart_End,
    const long*   plEnd_Begin,
    const long*   plPay_Begin,
    const double* pdMargin_Begin,
    const double* pdSpread_Begin,
    const double* pdNotional_Begin,
    SrtBasisCode  eBasis,
    double*       pdLevel,
    double*       pdLevelMultCoupon_Margin,
    double*       pdLevelMultCoupon_Spread)
{
    const int nSize = plStart_End - plStart_Begin;
    double    dPV_Margin, dPV_Spread;

    double dResult = 0., df_nI, df_nI_1;

    // margins and spreads
    _FixPV(
        lToday,
        szYC,
        plStart_Begin,
        plStart_End,
        plEnd_Begin,
        plPay_Begin,
        pdMargin_Begin,
        pdNotional_Begin,
        eBasis,
        &dPV_Margin,
        pdLevel,
        pdLevelMultCoupon_Margin);

    _FixPV(
        lToday,
        szYC,
        plStart_Begin,
        plStart_End,
        plEnd_Begin,
        plPay_Begin,
        pdSpread_Begin,
        pdNotional_Begin,
        eBasis,
        &dPV_Spread,
        pdLevel,
        pdLevelMultCoupon_Spread);

    dResult = dPV_Margin + dPV_Spread;

    // first period
    _ASSERTE(*plStart_Begin >= lToday);
    df_nI_1 = swp_f_df(lToday, *plStart_Begin, szYC);

    for (; plStart_Begin < plStart_End; ++plStart_Begin, ++plPay_Begin, ++pdNotional_Begin)
    {
        _ASSERTE(*plPay_Begin >= lToday);
        df_nI = swp_f_df(lToday, *plPay_Begin, szYC);
        dResult += *pdNotional_Begin * (df_nI_1 - df_nI);
        df_nI_1 = df_nI;
    }
    return dResult;
}

double FixLeg(
    long          lToday,
    const char*   szYC,
    const long*   plStart_Begin,
    const long*   plStart_End,
    const long*   plEnd_Begin,
    const long*   plPay_Begin,
    const double* pdCoupon_Begin,
    const double* pdNotional_Begin,
    SrtBasisCode  eBasis,
    double*       pdLevel,
    double*       pdLevelMultCoupon)
{
    double dPV = 0.;

    _FixPV(
        lToday,
        szYC,
        plStart_Begin,
        plStart_End,
        plEnd_Begin,
        plPay_Begin,
        pdCoupon_Begin,
        pdNotional_Begin,
        eBasis,
        &dPV,
        pdLevel,
        pdLevelMultCoupon);

    return dPV;
}

const char* Price_DiagGenEuro(
    // Market
    long        lToday,
    const char* szYC,
    const char* szVC,
    Err (*get_correl)(char*, double, double, double, double*),
    const char* szCorrel,
    /// GenMidat general Info
    const char* szRefRate,
    const char* szFreq,
    const char* szBasis,
    long        lTheoEnd,
    // Exercise
    long lFirstEx,
    // Fix Leg
    const long*   plFixStart_Begin,
    const long*   plFixStart_End,
    const long*   plFixEnd_Begin,
    const double* pdFixNotional_Begin,
    const double* pdFixRate_Begin,
    const double* pdFixFee_Begin,
    // Floating Leg
    const long*   plFltStart_Begin,
    const long*   plFltStart_End,
    const long*   plFltEnd_Begin,
    const double* pdFltNotional_Begin,
    const double* pdMargin_Begin,
    const double* pdSpread_Begin,
    // Calibratino specs
    double dMinTime,
    double dMinInterval,
    int    nUseVol,
    // Output
    int*    pnExBool,
    double* pdDiagPrice)
{
    int        i, k, l;
    Err        err = NULL;
    long       lMinDate;
    const long lMinInterval = (long)(dMinInterval / (double)YEARS_IN_DAY + 0.5);
    const int  lNFixNot     = plFixStart_End - plFixStart_Begin;
    const int  lNFltNot     = plFltStart_End - plFltStart_Begin;
    const int  nex          = lNFixNot;
    long       StartDate, SpotDate;

    double** Correl = NULL;
    int      NDimCorrel;

    double         exer_fee;
    double         dPrice;
    SrtCompounding srtFreq;
    SrtBasisCode   srtBasis;

    char*  cRec = NULL;
    double lastcaltime;

    _fill_date(lToday, &dMinTime, 1 + &dMinTime, &lMinDate);

    cRec = "REC";

    err = interp_compounding(szFreq, &srtFreq);
    if (err)
        return err;

    err = interp_basis(szBasis, &srtBasis);
    if (err)
        return err;

    SpotDate = add_unit(lToday, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    // lNFixNot = am->fix_leg->num_cpn;
    // lNFloatNot = am->fund_leg->num_cpn;
    //*nex = am->fix_leg->num_cpn;

    if (nex < 1)
        return 0;

    //	Skip calls to be exercised before today
    i = 0;
    l = 0;

    //*firstex = 0;
    memset(pnExBool, 0, nex * sizeof(int));
    memset(pdDiagPrice, 0, nex * sizeof(double));

    // EndDate = am->theoEndDate;

    lastcaltime = *plFixStart_Begin;
    for (k = 0; k < nex; ++k)
    {
        // Check if kth exercise is worth calibrating ..
        if ((plFixStart_Begin[k] - lastcaltime >= lMinInterval) &&
            (plFixStart_Begin[k] >= DMAX(lMinDate, lFirstEx)))
        {
            StartDate   = plFixStart_Begin[k];
            lastcaltime = StartDate;  // am->fix_leg->cpn[k].start_time;

            // NB: because of the assumption that fixed leg freq is the same as exercise freq
            /// but this exercise fee is in fact fixed fee.
            exer_fee = pdFixFee_Begin[k];

            while ((l < lNFltNot) && (plFltStart_Begin[l] < plFixStart_Begin[k] - 10.0))
            {
                ++l;
            }

            NDimCorrel = lNFixNot - k;
            err        = Compute_CoInitalSwaps_Correl2(
                SpotDate,
                (char*)szCorrel,
                get_correl,
                StartDate,
                lTheoEnd,
                srtFreq,
                srtBasis,
                NDimCorrel,
                &Correl);
            if (err)
                return err;

            err = AmortizedSwaptionShiftedLogForMAD(
                (char*)szYC,
                (char*)szVC,
                (char*)szRefRate,

                srtFreq,
                srtBasis,

                exer_fee,

                (long*)plFixStart_Begin + k,
                (long*)plFixEnd_Begin + k,
                lNFixNot - k,
                (double*)pdFixNotional_Begin + k,
                (double*)pdFixRate_Begin + k,

                (long*)plFltStart_Begin + l,
                (long*)plFltEnd_Begin + l,
                lNFltNot - l,
                (double*)pdFltNotional_Begin + l,
                (double*)pdMargin_Begin + l,
                (double*)pdSpread_Begin + l,

                SRT_PUT,

                10000,

                0,
                0,

                nUseVol,

                Correl,
                &dPrice);

            if (Correl)
            {
                free_dmatrix(Correl, 0, NDimCorrel - 1, 0, NDimCorrel - 1);
                Correl = NULL;
            }

            if (err)
                return err;

            pdDiagPrice[k] = dPrice;
            pnExBool[k]    = 1;
        }
    }

    // free_dmatrix(Correl, 0, NDimCorrel-1,0, NDimCorrel-1);
    // Correl = NULL;
    return err;
}
