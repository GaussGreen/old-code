#include "AmortMidatAutocal.h"

#include "AmortMidatADIUtils.h"
#include "AmortMidatCalib.h"
#include "EuroAmortSwaption.h"
#include "LGM2Fgrfn.h"
#include "LGMSVUtil.h"
#include "amortmidatcaller.h"
#include "amortmidatprodstruct.h"
#include "crtdbg.h"
#include "malloc.h"
#include "math.h"
#include "memory.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srt_h_lgmtypes.h"
#include "srtaccess.h"
#include "swp_h_amortswaption.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "swp_h_convslidingcorrel.h"

#define NSTD_LGM 7.0

Err amortMidat_compute_diagonal_prices_new(
    char*    yc,
    char*    vc,
    char*    ref,
    char*    cFreq,
    char*    cBasis,
    long     today,
    int      eod_flag,
    double   mintime,
    double   mininterval,
    int      UseVol,
    int*     nex,
    int*     firstex,
    int**    ex_bool,
    double** diag_prices,
    int      model,
    Err (*get_correl)(
        char* vol_curve_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,
    Err (*get_cash_vol)(
        char*   vol_curve_name,
        double  start_date,
        double  end_date,
        double  cash_strike,
        int     zero,
        char*   ref_rate_name,
        double* vol,
        double* power),

    AM_STR am);

static double _multiply_level_coupon(
    long          lToday,
    const char*   szYC,
    const long*   plStart_Begin,
    const long*   plStart_End,
    const long*   plEnd_Begin,
    const long*   plPay_Begin,
    SrtBasisCode  eBasis,
    const double* pdCoupon_Begin  // 0 to return level only
)
{
    const int nSize = plStart_End - plStart_Begin;
    int       nI;
    double    dResult = 0., df, cvg, cpn;
    for (nI = 0; nI < nSize; ++nI)
    {
        cvg = coverage(plStart_Begin[nI], plEnd_Begin[nI], eBasis);
        df  = swp_f_df(lToday, plPay_Begin[nI], szYC);
        cpn = pdCoupon_Begin ? pdCoupon_Begin[nI] : 1.;
        dResult += cpn * cvg * df;
    }
    return dResult;
}

static double _leg_coupon_pv_(
    long          lToday,
    const char*   szYC,
    const long*   plStart_Begin,
    const long*   plStart_End,
    const long*   plEnd_Begin,
    const long*   plPay_Begin,
    const double* pdCoupon_Begin,
    const double* pdNotional_Begin,
    SrtBasisCode  eBasis)
{
    const int nSize = plStart_End - plStart_Begin;
    int       nI;
    double    dResult = 0., df, cvg;
    for (nI = 0; nI < nSize; ++nI)
    {
        cvg = coverage(plStart_Begin[nI], plEnd_Begin[nI], eBasis);
        df  = swp_f_df(lToday, plPay_Begin[nI], szYC);
        dResult += pdNotional_Begin[nI] * pdCoupon_Begin[nI] * cvg * df;
    }
    return dResult;
}

static double _flt_leg(
    long          lToday,
    const char*   szYC,
    const long*   plStart_Begin,
    const long*   plStart_End,
    const long*   plEnd_Begin,
    const long*   plPay_Begin,
    const double* pdMargin_Begin,
    const double* pdSpread_Begin,
    const double* pdNotional_Begin,
    SrtBasisCode  eBasis)
{
    const int nSize = plStart_End - plStart_Begin;
    int       nI;
    double    dResult = 0., df_nI, df_nI_1;

    // margins and spreads
    dResult = _leg_coupon_pv_(
        lToday,
        szYC,
        plStart_Begin,
        plStart_End,
        plEnd_Begin,
        plPay_Begin,
        pdMargin_Begin,
        pdNotional_Begin,
        eBasis);
    dResult += _leg_coupon_pv_(
        lToday,
        szYC,
        plStart_Begin,
        plStart_End,
        plEnd_Begin,
        plPay_Begin,
        pdSpread_Begin,
        pdNotional_Begin,
        eBasis);

    // first period
    df_nI_1 = swp_f_df(lToday, plStart_Begin[0], szYC);
    for (nI = 0; nI < nSize; ++nI)
    {
        df_nI = swp_f_df(lToday, plPay_Begin[nI], szYC);
        dResult += pdNotional_Begin[nI] * (df_nI_1 - df_nI);
        df_nI_1 = df_nI;
    }
    return dResult;
}

static double _fix_leg(
    long          lToday,
    const char*   szYC,
    const long*   plStart_Begin,
    const long*   plStart_End,
    const long*   plEnd_Begin,
    const long*   plPay_Begin,
    const double* pdCoupon_Begin,
    const double* pdNotional_Begin,
    SrtBasisCode  eBasis)
{
    return _leg_coupon_pv_(
        lToday,
        szYC,
        plStart_Begin,
        plStart_End,
        plEnd_Begin,
        plPay_Begin,
        pdCoupon_Begin,
        pdNotional_Begin,
        eBasis);
}

static double _equi_strike(
    long          lToday,
    const char*   szYC,
    const long*   plFixStart_Begin,
    const long*   plFixStart_End,
    const long*   plFixEnd_Begin,
    const long*   plFixPay_Begin,
    const double* pdFixCoupon_Begin,
    const double* pdFixNotional_Begin,
    SrtBasisCode  eFixBasis,
    const long*   plFltStart_Begin,
    const long*   plFltStart_End,
    const long*   plFltEnd_Begin,
    const long*   plFltPay_Begin,
    const double* pdMargin_Begin,
    const double* pdSpread_Begin,
    const double* pdFltNotional_Begin,
    SrtBasisCode  eFltBasis,
    double        dEqui_Notional)
{
    const double dFltPV = _flt_leg(
        lToday,
        szYC,
        plFltStart_Begin,
        plFltStart_End,
        plFltEnd_Begin,
        plFltPay_Begin,
        pdMargin_Begin,
        pdSpread_Begin,
        pdFltNotional_Begin,
        eFltBasis);

    const double dFixPV = _fix_leg(
        lToday,
        szYC,
        plFixStart_Begin,
        plFixStart_End,
        plFixEnd_Begin,
        plFixPay_Begin,
        pdFixCoupon_Begin,
        pdFixNotional_Begin,
        eFixBasis);

    const double dMarginLevel = _multiply_level_coupon(
        lToday,
        szYC,
        plFltStart_Begin,
        plFltStart_End,
        plFltEnd_Begin,
        plFltPay_Begin,
        eFltBasis,
        pdMargin_Begin);

    const double dSpreadLevel = _multiply_level_coupon(
        lToday,
        szYC,
        plFltStart_Begin,
        plFltStart_End,
        plFltEnd_Begin,
        plFltPay_Begin,
        eFltBasis,
        pdSpread_Begin);

    const double dFixLevel = _multiply_level_coupon(
        lToday,
        szYC,
        plFixStart_Begin,
        plFixStart_End,
        plFixEnd_Begin,
        plFixPay_Begin,
        eFixBasis,
        0);

    // floating size
    const int nFltSize = plFltStart_End - plFltStart_Begin;
    // floating leg first df
    const double dDF_Front = swp_f_df(lToday, plFltStart_Begin[0], szYC);
    // floating leg last df
    const double dDF_Back = swp_f_df(lToday, plFltPay_Begin[nFltSize - 1], szYC);

    // match cash flow to return equivalent strike
    double dResult = (dFixPV - dFltPV) / dEqui_Notional;
    dResult += (dDF_Front - dDF_Back + dSpreadLevel + dMarginLevel);
    return dResult / dFixLevel;
}

static double _equi_notional(
    long          lToday,
    const char*   szYC,
    const long*   plStart_Begin,
    const long*   plStart_End,
    const long*   plEnd_Begin,
    const long*   plPay_Begin,
    const double* pdMargin_Begin,
    const double* pdSpread_Begin,
    const double* pdNotional_Begin,
    SrtBasisCode  eBasis)
{
    const int nSize = plStart_End - plStart_Begin;
    // floating leg first df
    const double dDF_Front = swp_f_df(lToday, plStart_Begin[0], szYC);
    // floating leg last df
    const double dDF_Back = swp_f_df(lToday, plPay_Begin[nSize - 1], szYC);
    const double dFltPV   = _flt_leg(
        lToday,
        szYC,
        plStart_Begin,
        plStart_End,
        plEnd_Begin,
        plPay_Begin,
        pdMargin_Begin,
        pdSpread_Begin,
        pdNotional_Begin,
        eBasis);

    const double dLevelMargin = _multiply_level_coupon(
        lToday, szYC, plStart_Begin, plStart_End, plEnd_Begin, plPay_Begin, eBasis, pdMargin_Begin);

    const double dLevelSpread = _multiply_level_coupon(
        lToday, szYC, plStart_Begin, plStart_End, plEnd_Begin, plPay_Begin, eBasis, pdSpread_Begin);

    _ASSERTE((dDF_Front - dDF_Back + dLevelMargin + dLevelSpread) > 0.);

    return dFltPV / (dDF_Front - dDF_Back + dLevelMargin + dLevelSpread);
}

static Err _LognormalSwapVol(
    long        lToday,
    const char* szYC,
    const char* szVC,
    const char* szFreq,
    const char* szBasis,
    const char* szRefRate,
    long        lEx,
    long        lStart,
    long        lEnd,
    double      dStrike,
    double*     pdVol)
{
    double       dPower, dFRA, dPrice;
    const double dEx   = (lEx - lToday) / 365.;
    char*        szErr = swp_f_vol((char*)szVC, lStart, lEnd, dStrike, pdVol, &dPower);
    if (szErr)
        return szErr;

    if (dPower < 1.)
    {
        szErr = swp_f_ForwardRate(
            lStart, lEnd, (char*)szFreq, (char*)szBasis, (char*)szYC, (char*)szRefRate, &dFRA);
        if (szErr)
            return szErr;
        dPrice = srt_f_optblknrm(dFRA, dStrike, *pdVol, dEx, 1., SRT_CALL, SRT_PREMIUM);
        szErr  = srt_f_optimpvol(dPrice, dFRA, dStrike, dEx, 1, SRT_CALL, SRT_LOGNORMAL, pdVol);
        if (szErr)
            return szErr;
    }
    return 0;
}

Err amortMidat_compute_diagonal_prices_new_equistrike(
    char*    yc,
    char*    vc,
    char*    ref,
    char*    cFreq,
    char*    cBasis,
    long     today,
    int      eod_flag,
    double   mintime,
    double   mininterval,
    int      UseVol,
    int*     nex,
    int*     firstex,
    int**    ex_bool,
    double** diag_prices,
    double*  pdequinotional,
    double*  pdlong_strikes,
    int      model,
    Err (*get_correl)(
        char* vol_curve_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,
    Err (*get_cash_vol)(
        char*   vol_curve_name,
        double  start_date,
        double  end_date,
        double  cash_strike,
        int     zero,
        char*   ref_rate_name,
        double* vol,
        double* power),

    AM_STR am)
{
    int i, j, k, l;
    Err err = NULL;

    long StartDate, EndDate, SpotDate, ExDate;

    int     lNFixNot;
    double* dFixNotionals = NULL;
    double* dFixRates     = NULL;
    double* dFee          = NULL;
    long *  lPayDates = NULL, *lStartDates = NULL, *lEndDates = NULL;

    long *lFixStartDates = NULL, *lFixEndDates = NULL;
    long *lFloatStartDates = NULL, *lFloatEndDates = NULL;

    double* dCoverages = NULL;

    int     lNFloatNot;
    double* dFloatNotionals = NULL;
    double* dMargins        = NULL;
    double* dSpreads        = NULL;

    double** Correl = NULL;
    int      NDimCorrel;

    double         exer_fee;
    double         dPrice;
    SrtCompounding srtFreq;
    SrtBasisCode   srtBasis;

    // Output of Badr function
    double*   dvPayDates             = NULL;
    double*   dvReplicatingStrikes   = NULL;
    double*   dvReplicatingNotionals = NULL;
    double*   dvReplicatingSwaptions = NULL;
    double    dFixedPV;
    double    dFloatPV;
    double    dSwapRate;
    SigKapTS* lgmSigKapTSPtrPtr = NULL;
    double    dLGMVol;

    char* cRec = NULL;

    double* pdMargins = NULL;
    double* dFixNots  = NULL;

    double lastcaltime, dEqui_Strike, dVol;
    double dEquiv_Notional;

    /*	SrtDiffusionType srt_vol_type; */

    cRec = "REC";

    err = interp_compounding(cFreq, &srtFreq);
    if (err)
    {
        smessage("Wrong Frequency");
        err = "Wrong Frequency";
        goto FREE_RETURN;
    }
    err = interp_basis(cBasis, &srtBasis);
    if (err)
    {
        smessage("Wrong Basis");
        err = "Wrong Basis";
        goto FREE_RETURN;
    }

    SpotDate = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    lNFixNot       = am->fix_leg->num_cpn;
    dFixNotionals  = (double*)calloc(lNFixNot, sizeof(double));
    dFixRates      = (double*)calloc(lNFixNot, sizeof(double));
    lFixStartDates = (long*)calloc(lNFixNot, sizeof(long));
    lFixEndDates   = (long*)calloc(lNFixNot, sizeof(long));
    dFee           = (double*)calloc(lNFixNot, sizeof(double));
    if ((!dFixNotionals) || (!dFixRates) || (!dFee))
    {
        err = "Allocation failed in amortMidat_compute_diagonal_prices";
        smessage("Allocation failed in amortMidat_compute_diagonal_prices");
        goto FREE_RETURN;
    }

    // Inputs for Badr function
    lStartDates = (long*)calloc(lNFixNot + 1, sizeof(long));
    lEndDates   = (long*)calloc(lNFixNot + 1, sizeof(long));
    lPayDates   = (long*)calloc(lNFixNot + 1, sizeof(long));
    dCoverages  = (double*)calloc(lNFixNot + 1, sizeof(double));
    dFixNots    = (double*)calloc(lNFixNot + 1, sizeof(double));
    pdMargins   = (double*)calloc(lNFixNot + 1, sizeof(double));

    // Outputs of Badr function
    dvPayDates             = (double*)calloc(lNFixNot + 1, sizeof(double));
    dvReplicatingStrikes   = (double*)calloc(lNFixNot + 1, sizeof(double));
    dvReplicatingNotionals = (double*)calloc(lNFixNot + 1, sizeof(double));
    dvReplicatingSwaptions = (double*)calloc(lNFixNot + 1, sizeof(double));
    if ((!dvPayDates) || (!dvReplicatingStrikes) || (!dvReplicatingNotionals) ||
        (!dvReplicatingSwaptions) || (!pdMargins) || (!lStartDates) || (!lEndDates) ||
        (!lPayDates) || (!dCoverages))
    {
        err = "Allocation failed in amortMidat_compute_diagonal_prices";
        smessage("Allocation failed in amortMidat_compute_diagonal_prices");
        goto FREE_RETURN;
    }

    pdMargins[0]   = 0;
    lEndDates[0]   = 0;
    lStartDates[0] = 0;
    lStartDates[0] = 0;
    dFixNots[0]    = 0;
    j              = 0;
    for (i = 0; i < lNFixNot; ++i)
    {
        dFixNotionals[i]  = am->fix_leg->notional[i];
        lFixStartDates[i] = am->fix_leg->cpn[i].start_date;
        lFixEndDates[i]   = am->fix_leg->cpn[i].end_date;

        dFixRates[i] = am->fix_leg->rate[i];
        dFee[i]      = am->fix_leg->fee[i];

        // Badr function
        lStartDates[i + 1] = am->fix_leg->cpn[i].start_date;
        lEndDates[i + 1]   = am->fix_leg->cpn[i].end_date;
        lPayDates[i + 1]   = am->fix_leg->cpn[i].end_date;
        dCoverages[i + 1]  = am->fix_leg->cpn[i].cvg;
        dFixNots[i + 1]    = am->fix_leg->notional[i];
        while ((j < am->fund_leg->num_cpn - 1) &&
               (am->fund_leg->cpn[j].start_date < am->fix_leg->cpn[i].start_date - 10))
        {
            ++j;
        }
        pdMargins[i + 1] = am->fund_leg->margin[j];
    }

    lNFloatNot = am->fund_leg->num_cpn;

    lFloatStartDates = (long*)calloc(lNFloatNot, sizeof(long));
    lFloatEndDates   = (long*)calloc(lNFloatNot, sizeof(long));
    dFloatNotionals  = (double*)calloc(lNFloatNot, sizeof(double));

    dMargins = (double*)calloc(lNFloatNot, sizeof(double));
    dSpreads = (double*)calloc(lNFloatNot, sizeof(double));
    if ((!dFloatNotionals) || (!dMargins))
    {
        err = "Allocation failed in amortMidat_compute_diagonal_prices";
        smessage("Allocation failed in amortMidat_compute_diagonal_prices");
        goto FREE_RETURN;
    }

    for (i = 0; i < lNFloatNot; ++i)
    {
        dFloatNotionals[i] = am->fund_leg->notional[i];

        lFloatStartDates[i] = am->fund_leg->cpn[i].start_date;
        lFloatEndDates[i]   = am->fund_leg->cpn[i].end_date;

        dMargins[i] = am->fund_leg->margin[i];
        dSpreads[i] = am->fund_leg->spread[i];
    }

    /*	Skip calls to be exercised before today */
    i = 0;
    j = 0;

    l = 0;

    *nex     = am->fix_leg->num_cpn - j;
    *firstex = j;

    *ex_bool     = (int*)calloc(*nex, sizeof(int));
    *diag_prices = (double*)calloc(*nex, sizeof(double));
    if ((!(*ex_bool)) || (!(*diag_prices)))
    {
        err = "allocation failed in amortMidat_compute_diagonal_prices";
        smessage("allocation failed in amortMidat_compute_diagonal_prices");
        goto FREE_RETURN;
    }

    EndDate = am->theoEndDate;

    lastcaltime     = 0.0;
    dEquiv_Notional = _equi_notional(
        today,
        yc,
        lFloatStartDates,
        lFloatStartDates + lNFloatNot,
        lFloatEndDates,
        lFloatEndDates,
        dMargins,
        dSpreads,
        dFloatNotionals,
        srtBasis  // eFltBasis; wrong !!!
    );

    *pdequinotional = dEquiv_Notional;

    if (am->fix_leg->cpn[j].start_time < DMAX(mintime, am->call[0].ex_time))
    {
        (*ex_bool)[0]     = 0;
        (*diag_prices)[0] = 0;
    }
    else
    {
        (*ex_bool)[0] = 1;
        exer_fee      = dFee[j];
        StartDate     = am->fix_leg->cpn[j].start_date;
        lastcaltime   = am->fix_leg->cpn[j].start_time;

        while ((l < am->fund_leg->num_cpn) &&
               (am->fund_leg->cpn[l].start_time < am->fix_leg->cpn[j].start_time - 10.0 / 365.0))
        {
            ++l;
        }

        /// preivously  ... european amortizing swaption prices
        /// currently, equivalent strike
        dEqui_Strike = _equi_strike(
            today,
            yc,
            //// fixed leg
            lFixStartDates + j,
            lFixStartDates + lNFixNot,
            lFixEndDates + j,
            lFixEndDates + j,
            dFixRates + j,
            dFixNotionals + j,
            srtBasis,  // eFixBasis,
            /// floating leg
            lFloatStartDates + l,
            lFloatStartDates + lNFloatNot,
            lFloatEndDates + l,
            lFloatEndDates + l,
            dMargins + l,
            dSpreads + l,
            dFloatNotionals + l,
            srtBasis,  // eFltBasis; wrong !!!
            dEquiv_Notional);

        pdlong_strikes[0] = dEqui_Strike;

        /// get vol
        err = _LognormalSwapVol(
            today,
            yc,
            vc,
            cFreq,
            cBasis,
            ref,
            *(lFloatStartDates + l),
            *(lFloatStartDates + l),
            *(lFloatEndDates + lNFloatNot - 1),
            dEqui_Strike,
            &dVol);

        /// get pricer
        err = swp_f_Swaption(
            *(lFloatStartDates + l),
            *(lFloatEndDates + lNFloatNot - 1),
            cFreq,
            cBasis,
            dVol,
            dEqui_Strike,  //// wrong ! must adjust for exercise fee
            cRec,          // szPayRec,
            ref,
            yc,
            "PREMIUM",
            "LOGNORMAL",
            &dPrice);

        if (err)
            return err;

        /*		if(model==1)
                        {
                                NDimCorrel = lNFixNot-j;
                                err = Compute_CoInitalSwaps_Correl2(SpotDate,
                                                                                                CorrelName,
                                                                                                get_correl,
                                                                                                StartDate,
                                                                                                EndDate,
                                                                                                srtFreq,
                                                                                                srtBasis,
                                                                                                NDimCorrel,
                                                                                                &Correl);
                                if(err)
                                {
                                        goto FREE_RETURN;
                                }

                                err = AmortizedSwaptionShiftedLogForMAD(
                                                         yc,
                                                         vc,
                                                         ref,

                                                         srtFreq,
                                                         srtBasis,

                                                         exer_fee,

                                                         lFixStartDates+j,
                                                         lFixEndDates+j,
                                                         lNFixNot-j,
                                                         dFixNotionals+j,
                                                         dFixRates+j,

                                                         lFloatStartDates+l,
                                                         lFloatEndDates+l,
                                                         lNFloatNot-l,
                                                         dFloatNotionals+l,
                                                         dMargins+l,
                                                         dSpreads+l,

                                                         SRT_PUT,

                                                         10000,

                                                         0,
                                                         0,

                                                         UseVol,

                                                         Correl,
                                                         &dPrice
                                            );
                                if(err)
                                {
                                        smessage("Pb in pricing European with Shifted Log SMM");
                                        goto FREE_RETURN;
                                }

                                if(Correl)
                                {
                                        free_dmatrix(Correl, 0, NDimCorrel-1,0, NDimCorrel-1);
                                        Correl = NULL;
                                }
                        }
                        else
                        {
                                ExDate = add_unit (StartDate, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
                                err = EuropeanAmortizingSwaption2(today,
                                                                        ExDate,
                                                                        EndDate,
                                                                        lNFixNot-j,
                                                                        lStartDates+j,
                                                                        lEndDates+j,
                                                                        lPayDates+j,
                                                                        dCoverages+j,
                                                                        dCoverages+j,
                                                                        dFixRates[0],
                                                                        dFixNots+j,
                                                                        yc,
                                                                        vc,
                                                                        cFreq,
                                                                        cBasis,
                                                                        cRec,
                                                                        ref,
                                                                        get_cash_vol,

                                                                        // Output
                                                                        dvPayDates,
                                                                        dvReplicatingStrikes,
                                                                        dvReplicatingNotionals,
                                                                        &dPrice,
                                                                        &dFixedPV,
                                                                        &dFloatPV,
                                                                        &dSwapRate,
                                                                        &lgmSigKapTSPtrPtr,
                                                                        &dLGMVol,
                                                                        dvReplicatingSwaptions,

                                                                        //Inputs
                                                                        pdMargins+j,
                                                                        0,//dFixRates+j,
                                                                        0,
                                                                        exer_fee);
                                if(err)
                                {
                                        smessage("Pb in pricing European with 1F");
                                        goto FREE_RETURN;
                                }
                        }
                        */

        (*diag_prices)[0] = dPrice;
    }

    for (k = 1; k < *nex; ++k)
    {
        if ((am->fix_leg->cpn[j + k].start_time - lastcaltime < mininterval) ||
            (am->fix_leg->cpn[j + k].start_time < DMAX(mintime, am->call[0].ex_time)))
        {
            (*ex_bool)[k]     = 0;
            (*diag_prices)[k] = 0;
        }
        else
        {
            StartDate   = am->fix_leg->cpn[j + k].start_date;
            lastcaltime = am->fix_leg->cpn[j + k].start_time;

            exer_fee = dFee[j + k];

            (*ex_bool)[k] = 1;

            while ((l < am->fund_leg->num_cpn) &&
                   (am->fund_leg->cpn[l].start_time <
                    am->fix_leg->cpn[j + k].start_time - 10.0 / 365.0))
            {
                ++l;
            }

            /// preivously  ... european amortizing swaption prices
            /// currently, equivalent strike
            dEqui_Strike = _equi_strike(
                today,
                yc,
                //// fixed leg
                lFixStartDates + k + j,
                lFixStartDates + lNFixNot,
                lFixEndDates + k + j,
                lFixEndDates + k + j,
                dFixRates + k + j,
                dFixNotionals + k + j,
                srtBasis,  // eFixBasis,
                /// floating leg
                lFloatStartDates + l,
                lFloatStartDates + lNFloatNot,
                lFloatEndDates + l,
                lFloatEndDates + l,
                dMargins + l,
                dSpreads + l,
                dFloatNotionals + l,
                srtBasis,  // eFltBasis; wrong !!!
                dEquiv_Notional);

            pdlong_strikes[k] = dEqui_Strike;

            /// get vol
            err = _LognormalSwapVol(
                today,
                yc,
                vc,
                cFreq,
                cBasis,
                ref,
                *(lFloatStartDates + l),
                *(lFloatStartDates + l),
                *(lFloatEndDates + lNFloatNot - 1),
                dEqui_Strike,
                &dVol);

            /// get pricer
            err = swp_f_Swaption(
                *(lFloatStartDates + l),
                *(lFloatEndDates + lNFloatNot - 1),
                cFreq,
                cBasis,
                dVol,
                dEqui_Strike,  //// wrong ! must adjust for exercise fee
                cRec,          // szPayRec,
                ref,
                yc,
                "PREMIUM",
                "LOGNORMAL",
                &dPrice);

            if (err)
                return err;

            /*if(model==1)
            {
                    NDimCorrel = lNFixNot-(k+j);
                    err = Compute_CoInitalSwaps_Correl2(SpotDate,
                                                                                    CorrelName,
                                                                                    get_correl,
                                                                                    StartDate,
                                                                                    EndDate,
                                                                                    srtFreq,
                                                                                    srtBasis,
                                                                                    NDimCorrel,
                                                                                    &Correl);
                    if(err)
                    {
                            goto FREE_RETURN;
                    }

                    err = AmortizedSwaptionShiftedLogForMAD(
                                     yc,
                                     vc,
                                     ref,

                                     srtFreq,
                                     srtBasis,

                                     exer_fee,


                                     lFixStartDates+k+j,
                                     lFixEndDates+k+j,
                                     lNFixNot-(k+j),
                                     dFixNotionals+k+j,
                                     dFixRates+k+j,

                                     lFloatStartDates+l,
                                     lFloatEndDates+l,
                                     lNFloatNot-l,
                                     dFloatNotionals+l,
                                     dMargins+l,
                                     dSpreads+l,

                                     SRT_PUT,

                                     10000,

                                     0,
                                     0,

                                     UseVol,

                                     Correl,
                                     &dPrice
                        );
                    if(err)
                    {
                            smessage("Pb in pricing European with Shifted Log SMM");
                            goto FREE_RETURN;
                    }

                    if(Correl)
                    {
                            free_dmatrix(Correl, 0, NDimCorrel-1,0, NDimCorrel-1);
                            Correl = NULL;
                    }
            }
            else
            {
                    ExDate = add_unit (StartDate, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
                    err = EuropeanAmortizingSwaption2(today,
                                                    ExDate,
                                                    EndDate,
                                                    lNFixNot-(k+j),
                                                    lStartDates+(k+j),
                                                    lEndDates+(k+j),
                                                    lPayDates+(k+j),
                                                    dCoverages+(k+j),
                                                    dCoverages+(k+j),
                                                    dFixRates[0],
                                                    dFixNots+(k+j),
                                                    yc,
                                                    vc,
                                                    cFreq,
                                                    cBasis,
                                                    cRec,
                                                    ref,
                                                    get_cash_vol,

                                                    // Output
                                                    dvPayDates,
                                                    dvReplicatingStrikes,
                                                    dvReplicatingNotionals,
                                                    &dPrice,
                                                    &dFixedPV,
                                                    &dFloatPV,
                                                    &dSwapRate,
                                                    &lgmSigKapTSPtrPtr,
                                                    &dLGMVol,
                                                    dvReplicatingSwaptions,

                                                    //Inputs
                                                    pdMargins+k+j,
                                                    0,//dFixRates+k+j,
                                                    0,
                                                    exer_fee);
                    if(err)
                    {
                            smessage("Pb in pricing European with 1F");
                            goto FREE_RETURN;
                    }
            }*/

            (*diag_prices)[k] = dPrice;
        }
    }

FREE_RETURN:

    if (lStartDates)
        free(lStartDates);
    if (lEndDates)
        free(lEndDates);
    if (lPayDates)
        free(lPayDates);
    if (dCoverages)
        free(dCoverages);
    if (dFixNots)
        free(dFixNots);
    if (pdMargins)
        free(pdMargins);

    if (dvPayDates)
        free(dvPayDates);
    if (dvReplicatingStrikes)
        free(dvReplicatingStrikes);
    if (dvReplicatingNotionals)
        free(dvReplicatingNotionals);
    if (dvReplicatingSwaptions)
        free(dvReplicatingSwaptions);

    if (dFixNotionals)
        free(dFixNotionals);
    if (lFixStartDates)
        free(lFixStartDates);
    if (lFixEndDates)
        free(lFixEndDates);
    if (dFixRates)
        free(dFixRates);
    if (dFee)
        free(dFee);
    if (dFloatNotionals)
        free(dFloatNotionals);
    if (lFloatStartDates)
        free(lFloatStartDates);
    if (lFloatEndDates)
        free(lFloatEndDates);
    if (dMargins)
        free(dMargins);
    if (dSpreads)
        free(dSpreads);

    if (Correl)
    {
        // free_dmatrix(Correl, 0, NDimCorrel-1,0, NDimCorrel-1);
        Correl = NULL;
    }

    if (err)
    {
        if (*ex_bool)
        {
            free(*ex_bool);
            *ex_bool = NULL;
        }
        if (*diag_prices)
        {
            free(*diag_prices);
            *diag_prices = NULL;
        }
    }

    /// unreferenced variables
    get_cash_vol;
    CorrelName;
    model;
    UseVol;
    eod_flag;
    get_correl;
    dFixedPV;
    dSwapRate;
    ExDate;
    dLGMVol;
    NDimCorrel;
    dFloatPV;

    return err;
}

Err am_calib_und_new(
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I, 1: E */

    char* yc,                 /*	yc */
    char* vc,                 /*	vc */
    char* default_ref,        /*	ref rate */
    char* default_swap_freq,  /*	swap freq */
    char* default_swap_basis, /*	swap basis */

    char* ref,        /*	ref rate */
    char* swap_freq,  /*	swap freq */
    char* swap_basis, /*	swap basis */

    double lambda, /*	lambda if unique */
    double alpha,  /*	alpha */
    double gamma,  /*	gamma */
    double rho,    /*	rho */
    /*	Calib params */
    double mintime,
    double mininterval,

    int    notperiod,
    double max_std_short,
    int    one2F,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                       to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */
    int    use_jump,
    double max_var_jump,

    int strike_type,
    int european_model,

    Err (*get_correl)(
        char* correl_cube_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,
    /*
                    Err  (*GetVolForBadr)( Date, Date, double, SRT_Boolean, double *),
                    char *cVolType,
    */
    /*	End of calib params */
    AM_STR am,          /*	structure */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    AM_UND und);

/*	Fill underlying structure from calibration instruments */
Err am_calib_und_new_reduc(
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I, 1: E */

    char* yc,                 /*	yc */
    char* vc,                 /*	vc */
    char* default_ref,        /*	ref rate */
    char* default_swap_freq,  /*	swap freq */
    char* default_swap_basis, /*	swap basis */

    char* ref,        /*	ref rate */
    char* swap_freq,  /*	swap freq */
    char* swap_basis, /*	swap basis */

    double lambda, /*	lambda if unique */
    double alpha,  /*	alpha */
    double gamma,  /*	gamma */
    double rho,    /*	rho */
    /*	Calib params */
    double mintime,
    double mininterval,

    int    notperiod,
    double max_std_short,
    int    one2F,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                       to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */
    int    use_jump,
    double max_var_jump,

    int strike_type,
    int european_model,

    Err (*get_correl)(
        char* correl_cube_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,
    /*
                    Err  (*GetVolForBadr)( Date, Date, double, SRT_Boolean, double *),
                    char *cVolType,
    */
    /*	End of calib params */
    AM_STR am,          /*	structure */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    AM_UND und)
{
    int   i, nex, firstex;
    long* lex = NULL;

    long end_struct;  //,
                      // end_swap;
    double *long_strikes = NULL, *short_strikes = NULL;

    int     Nfix;
    int*    ex_bool         = NULL;
    double* diag_prices     = NULL;
    double* ex_fee          = NULL;
    double* fixNotional     = NULL;
    long*   fix_start_dates = NULL;
    long*   fix_end_dates   = NULL;
    long*   fix_pay_dates   = NULL;

    int     Nfloat;
    double* floatNotional     = NULL;
    long*   float_start_dates = NULL;
    long*   float_end_dates   = NULL;
    long*   float_pay_dates   = NULL;
    double* float_margin      = NULL;
    double* float_spread      = NULL;

    long start_date;
    long end_date;

    int UseVol;

    SrtCompounding floatfreq;
    SrtBasisCode   floatbasis;

    char* float_freq  = NULL;
    char* float_basis = NULL;

    AM_FUND_LEG fund_leg;
    AM_FIX_LEG  fix_leg;
    double*     pdlong_strikes  = 0;
    double*     pdfixNotional   = 0;
    double*     pdfloatNotional = 0;
    double      dequinotional;

    Err err = NULL;

    /*	Initialise */
    fix_leg  = am->fix_leg;
    fund_leg = am->fund_leg;

    und->sigma_date = NULL;
    und->sigma_time = NULL;
    und->sigma      = NULL;

    und->has_inst_data = 0;
    cpd_init_calib_inst_data(&(und->inst_data));

    und->num_prices     = 0;
    und->exercise_dates = NULL;
    und->mkt_prices     = NULL;
    und->mdl_prices     = NULL;

    und->lambda = lambda;

    /* write default values */
    und->lambda_n      = 1;
    und->pdlambda      = 0;
    und->pdlambda_date = 0;
    und->pdlambda_time = 0;

    und->alpha = alpha;
    und->gamma = gamma;
    und->rho   = rho;
    und->today = today;

    strcpy(und->yc, yc);
    strcpy(und->vc, vc);
    strcpy(und->ref, ref);
    strcpy(und->swap_freq, swap_freq);
    strcpy(und->swap_basis, swap_basis);
    strcpy(und->name, "CALIB");

    /*	Exercise dates for calibration */
    end_struct = fund_leg->cpn[fund_leg->num_cpn - 1].pay_date;
    //	end_swap = end_struct;
    //	if (fix_leg->cpn[fix_leg->num_cpn-1].pay_date > end_swap)
    //	{
    //		end_swap = fix_leg->cpn[fix_leg->num_cpn-1].pay_date;
    //	}

    if (am->num_calls > 0 && !(am->num_calls == 1 && am->call[0].ex_date <= today + eod_flag))
    {
        /*	If call dates, choose call dates as option expiries for calibration */
        nex = am->num_calls;
        lex = (long*)calloc(nex, sizeof(long));
        if (!lex)
        {
            err = "Allocation error (2) in am_calib_und";
            goto FREE_RETURN;
        }

        for (i = 0; i < nex; i++)
        {
            lex[i] = am->call[i].ex_date;
        }
    }
    else
    {
        /*	If no call dates, exit */
        am->num_calls = 0;
        goto FREE_RETURN;
    }

    UseVol = 0;
    if (strike_type == 5)
    {
        UseVol = 1;
    }

    Nfix = am->fix_leg->num_cpn;

    //// Changed by Albert Wang on 12/19/03
    //// previously
    // short_strikes = NULL;
    short_strikes = (double*)calloc(am->fix_leg->num_cpn, sizeof(double));

    fix_start_dates = (long*)calloc(am->fix_leg->num_cpn, sizeof(long));
    fix_end_dates   = (long*)calloc(am->fix_leg->num_cpn, sizeof(long));
    fix_pay_dates   = (long*)calloc(am->fix_leg->num_cpn, sizeof(long));
    fixNotional     = (double*)calloc(am->fix_leg->num_cpn, sizeof(double));
    long_strikes    = (double*)calloc(am->fix_leg->num_cpn, sizeof(double));
    ex_fee          = (double*)calloc(am->fix_leg->num_cpn, sizeof(double));
    if ((!fixNotional) || (!long_strikes) || (!ex_fee) || (!fix_start_dates) || (!fix_end_dates) ||
        (!fix_pay_dates))
    {
        err = "memory allocation failed in am_calib_und";
        smessage("memory allocation failed in am_calib_und");
        goto FREE_RETURN;
    }
    for (i = 0; i < am->fix_leg->num_cpn; ++i)
    {
        fix_start_dates[i] = am->fix_leg->cpn[i].start_date;
        fix_end_dates[i]   = am->fix_leg->cpn[i].end_date;
        fix_pay_dates[i]   = am->fix_leg->cpn[i].end_date;
        fixNotional[i]     = am->fix_leg->notional[i];
        long_strikes[i]    = am->fix_leg->rate[i];
        short_strikes[i]   = am->fix_leg->rate[i];
        ex_fee[i]          = am->fix_leg->fee[i];
    }

    Nfloat            = am->fund_leg->num_cpn;
    float_start_dates = (long*)calloc(am->fund_leg->num_cpn, sizeof(long));
    float_end_dates   = (long*)calloc(am->fund_leg->num_cpn, sizeof(long));
    float_pay_dates   = (long*)calloc(am->fund_leg->num_cpn, sizeof(long));
    floatNotional     = (double*)calloc(am->fund_leg->num_cpn, sizeof(double));
    float_margin      = (double*)calloc(am->fund_leg->num_cpn, sizeof(double));
    float_spread      = (double*)calloc(am->fund_leg->num_cpn, sizeof(double));
    if ((!floatNotional) || (!float_margin) || (!float_spread) || (!float_start_dates) ||
        (!float_end_dates) || (!float_pay_dates))
    {
        err = "memory allocation failed in am_calib_und";
        smessage("memory allocation failed in am_calib_und");
        goto FREE_RETURN;
    }
    for (i = 0; i < am->fund_leg->num_cpn; ++i)
    {
        float_start_dates[i] = am->fund_leg->cpn[i].start_date;
        float_end_dates[i]   = am->fund_leg->cpn[i].end_date;
        float_pay_dates[i]   = am->fund_leg->cpn[i].end_date;
        floatNotional[i]     = am->fund_leg->notional[i];
        float_margin[i]      = am->fund_leg->margin[i];
        float_spread[i]      = am->fund_leg->spread[i];
    }

    start_date = am->fix_leg->cpn[0].start_date;
    end_date   = am->theoEndDate;

    err = swp_f_get_ref_rate_details(ref, &floatbasis, &floatfreq);
    err = translate_compounding(&float_freq, floatfreq);
    err = translate_basis(&float_basis, floatbasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    pdlong_strikes  = _alloca(am->fix_leg->num_cpn * sizeof(double));
    pdfixNotional   = _alloca(Nfix * sizeof(double));
    pdfloatNotional = _alloca(Nfloat * sizeof(double));

    /// Previously err = amortMidat_compute_diagonal_prices_new(yc, vc, ref, swap_freq, swap_basis,
    /// 1. Compute equivalent notional and strikes for the european amortizing underlyings and
    /// 2. Get Market price for these equivalent strike eurpean swaptions
    err = amortMidat_compute_diagonal_prices_new_equistrike(
        yc,
        vc,
        ref,
        swap_freq,
        swap_basis,
        today,
        eod_flag,
        mintime,
        mininterval,
        UseVol,
        &nex,
        &firstex,
        &ex_bool,
        &diag_prices,
        &dequinotional,
        pdlong_strikes,
        european_model,
        get_correl,
        CorrelName,
        get_cash_vol,
        am);

    /// overwrite the notionals
    _memset(&dequinotional, pdlong_strikes, am->fix_leg->num_cpn, sizeof(*pdlong_strikes));

    und->num_prices     = nex;
    und->exercise_dates = (double*)calloc(nex, sizeof(double));
    und->mkt_prices     = (double*)calloc(nex, sizeof(double));
    und->mdl_prices     = (double*)calloc(nex, sizeof(double));

    for (i = 0; i < nex; ++i)
    {
        und->exercise_dates[i] = am->fix_leg->cpn[firstex + i].start_date;
        und->mkt_prices[i]     = diag_prices[i];
    }

    /// Calibrate term structure again based on the combination of amortizing and equivalent strike
    /// swaptions
    /// 3. Calibrate term structure based on these equivalent strike swaptions
    err = amortMidat_cpd_calib_diagonal_new(
        notperiod,
        yc,
        vc,
        default_ref,
        default_swap_basis,
        default_swap_freq,

        get_cash_vol,  // GetCpdAutocalCashVol

        0,
        1,  /// shift type

        ex_bool,
        diag_prices,
        ex_fee,

        swap_freq,
        swap_basis,
        Nfix,
        fix_start_dates,
        fix_end_dates,
        fix_pay_dates,
        pdlong_strikes,  // long_strikes,
        pdfixNotional,   // fixNotional,

        float_freq,
        float_basis,
        Nfloat,
        float_start_dates,
        float_end_dates,
        float_pay_dates,
        float_margin,
        float_spread,
        pdfloatNotional,  // floatNotional

        pdlong_strikes,  // short_strikes

        strike_type,

        max_std_short,

        fix_lambda,
        one_f_equi,

        skip_last,

        use_jump,
        max_var_jump,

        &lambda,

        one2F,

        alpha,
        gamma,
        rho,
        &(und->sigma_n),
        &(und->sigma_time),
        &(und->sigma));
    if (err)
    {
        goto FREE_RETURN;
    }

    und->lambda = lambda;

    /* write default values */
    und->lambda_n      = 1;
    und->pdlambda      = 0;
    und->pdlambda_date = 0;
    und->pdlambda_time = 0;

    und->has_inst_data = 1;

    und->sigma_date = (double*)calloc(und->sigma_n, sizeof(double));
    if (!und->sigma_date)
    {
        err = "Allocation error (6) in am_calib_und";
        goto FREE_RETURN;
    }

    for (i = 0; i < und->sigma_n; i++)
    {
        und->sigma_date[i] = today + und->sigma_time[i] * DAYS_IN_YEAR + 1.0e-08;
    }

FREE_RETURN:

    if (err)
        am_free_und(und);

    //	if (float_freq) free(float_freq);
    //	if (float_basis) free(float_basis);

    if (fix_start_dates)
        free(fix_start_dates);
    if (fix_end_dates)
        free(fix_end_dates);
    if (fix_pay_dates)
        free(fix_pay_dates);
    if (fixNotional)
        free(fixNotional);
    if (long_strikes)
        free(long_strikes);
    if (ex_fee)
        free(ex_fee);

    if (float_start_dates)
        free(float_start_dates);
    if (float_end_dates)
        free(float_end_dates);
    if (float_pay_dates)
        free(float_pay_dates);
    if (floatNotional)
        free(floatNotional);
    if (float_margin)
        free(float_margin);
    if (float_spread)
        free(float_spread);

    if (lex)
        free(lex);
    if (short_strikes)
        free(short_strikes);

    if (ex_bool)
        free(ex_bool);
    if (diag_prices)
        free(diag_prices);

    return err;
}

Err am_calib_und_new_ts(
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I, 1: E */

    char* yc,                 /*	yc */
    char* vc,                 /*	vc */
    char* default_ref,        /*	ref rate */
    char* default_swap_freq,  /*	swap freq */
    char* default_swap_basis, /*	swap basis */

    char* ref,        /*	ref rate */
    char* swap_freq,  /*	swap freq */
    char* swap_basis, /*	swap basis */

    int     nlambda,
    double* pdlambda_time,
    double* pdlambda,

    double alpha, /*	alpha */
    double gamma, /*	gamma */
    double rho,   /*	rho */
    /*	Calib params */
    double mintime,
    double mininterval,

    int    notperiod,
    double max_std_short,
    int    one2F,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                       to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */
    int    use_jump,
    double max_var_jump,

    int strike_type,
    int european_model,

    Err (*get_correl)(
        char* correl_cube_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,
    /*
                    Err  (*GetVolForBadr)( Date, Date, double, SRT_Boolean, double *),
                    char *cVolType,
    */
    /*	End of calib params */
    AM_STR am,          /*	structure */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    AM_UND und);

typedef struct _am_adi_arg_amortmidatautocal_
{
    AM_PAY_ARG    m_pBase;
    int*          pnum_fund_cpn_done;
    int*          pnum_fix_cpn_done;
    double*       pdInit_Prev;
    double*       pdIV;
    const double* pdStrike;
    long          lToday;

    /// exercise probability related variables
    int nNumEx;
    int nCompExProb;

    double*** pppdExIndicator;  // [nth of exercises][state i] [ state - j]
    double*** pppdCondiDF;      // [nth of exercises][state i] [ state - j]
    double*   pdAvgSwapRate;

    int* pdUB_D1;  //[nth of exercises]
    int* pdLB_D1;  //[nth of exercises]
    int* pdUB_D2;  //[nth of exercises]
    int* pdLB_D2;  //[nth of exercises]
    int* pnEx;

} genmidat_arg, *ptr_genmidat_arg;

Err am_fill_check_all_struct_AmortMidatAutocal(
    /*	Today's date */
    long today,
    long theoEndDate,

    /*	The underlying */
    int use_calib, /*	0: use lgm2dund, 1: calibrate */

    /*		if calib */
    char*  yc,         /*	yc */
    char*  vc,         /*	vc (only if calib) */
    char*  ref,        /*	ref rate (only if calib) */
    char*  swap_freq,  /*	swap freq (only if calib) */
    char*  swap_basis, /*	swap basis (only if calib) */
    double lambda,     /*	lambda if unique */
    double alpha,      /*	alpha */
    double gamma,      /*	gamma */
    double rho,        /*	rho */
    /*	End of calib params */

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */

    /*		funding */
    char*   fund_ref,
    double* fund_not,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,

    /*		cf */
    double* fix_not,
    int     fix_ncpn,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee,

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    double mintime,
    double mininterval,

    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    Err (*get_correl)(
        char* vol_curve_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,

    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                               to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Results */
    AM_STR am,
    AM_UND und,
    int*   call_feat, /*	0: No callable feature to be valued
                              1: Callable feature to be valued through adi */
    AM_ADI_ARG adi_arg)
{
    Err  szErr = NULL;
    char swap_freq_loc[256];
    char swap_basis_loc[256];
    int  intfreq;

    /*	Initialisation */
    am->fund_leg    = NULL;
    am->fix_leg     = NULL;
    am->call        = NULL;
    am->theoEndDate = theoEndDate;

    und->sigma_n    = 0;
    und->sigma_date = NULL;
    und->sigma_time = NULL;
    und->sigma      = NULL;

    strcpy(swap_basis_loc, fix_basis[0]);

    intfreq = (int)(12 * (fix_pay[fix_ncpn - 1] - fix_start[fix_ncpn - 1]) / 365.0 + 0.5);
    if (intfreq == 1)
    {
        strcpy(swap_freq_loc, "M");
    }
    else if (intfreq == 3)
    {
        strcpy(swap_freq_loc, "Q");
    }
    else if (intfreq == 6)
    {
        strcpy(swap_freq_loc, "S");
    }
    else if (intfreq == 12)
    {
        strcpy(swap_freq_loc, "A");
    }

    cpd_init_calib_inst_data(&(und->inst_data));
    und->has_inst_data = 0;
    und->today         = today;

    adi_arg->time     = NULL;
    adi_arg->date     = NULL;
    adi_arg->sig1     = NULL;
    adi_arg->sig_time = NULL;
    adi_arg->void_prm = NULL;
    adi_arg->is_event = NULL;
    adi_arg->ifr      = NULL;

    /*	Funding leg */

    am->fund_leg = (AM_FUND_LEG)malloc(sizeof(am_fund_leg));
    if (!am->fund_leg)
    {
        szErr = "Memory allocation error (1) in am_fill_check_all_AmortMidatAutocal";
        goto FREE_RETURN;
    }

    szErr = am_fill_fund_leg(
        und->today,
        eod_fix_flag,
        fund_not,
        fund_ncpn,
        fund_fix,
        fund_start,
        fund_end,
        fund_pay,
        fund_basis,
        fund_spr,
        fund_mrg,
        am->fund_leg);
    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	fix leg */

    am->fix_leg = (AM_FIX_LEG)malloc(sizeof(am_fix_leg));
    if (!am->fix_leg)
    {
        szErr = "Memory allocation error (2) in am_fill_check_all_struct_AmortMidatAutocal";
        goto FREE_RETURN;
    }

    szErr = am_fill_fix_leg(
        und->today,
        eod_fix_flag,
        fix_not,
        fix_ncpn,
        //			fix_fix,
        fix_start,
        fix_end,
        fix_pay,
        fix_basis,
        fix_rate,
        fix_fee,
        am->fix_leg);

    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	Calls */

    if (ncall > 0 && ex_date[ncall - 1] >= und->today + eod_ex_flag)
    {
        szErr = am_fill_calls(und->today, eod_ex_flag, ncall, pay_rec, ex_date, set_date, fee, am);
        if (szErr)
        {
            goto FREE_RETURN;
        }
    }
    else
    {
        am->num_calls = 0;
        am->call      = NULL;
    }

    /*	Underlying */
    if (am->num_calls > 0 && am->call)
    {
        if (use_calib)
        {
            // previously ... if reversion is needed.
            szErr = am_calib_und_new(
                // szErr = am_calib_und_new_AmortMidatAutocal(
                today,
                eod_ex_flag,

                yc,
                vc,
                ref,
                swap_freq,
                swap_basis,

                fund_ref,
                swap_freq_loc,
                swap_basis_loc,
                lambda,
                alpha,
                gamma,
                rho,

                mintime,
                mininterval,

                notperiod,
                max_std_short,
                one2F,
                fix_lambda,

                one_f_equi,

                skip_last,

                use_jump,
                max_var_jump,
                strike_type,
                european_model,

                get_correl,
                CorrelName,
                /*
                                                                GetVolForBadr,
                                                                cVolType,
                */
                am,

                get_cash_vol,

                und);
        }
        else
        {
            szErr = am_fill_und(lgm2dund, vc, ref, swap_freq, swap_basis, und);
        }
    }

    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	Tree */

    if (am->num_calls > 0 && am->call)
    {
        szErr = am_fill_adi_arg(und, am, get_cash_vol, req_stp, req_stpx, adi_arg);
        if (szErr)
        {
            goto FREE_RETURN;
        }

        *call_feat = 1;
    }
    else
    {
        *call_feat = 0;
    }

FREE_RETURN:

    if (szErr)
    {
        am_free_all_struct(am, und, *call_feat, adi_arg);
    }

    return szErr;
}

Err am_fill_check_all_struct_AmortMidatAutocal_reduc(
    /*	Today's date */
    long today,
    long theoEndDate,

    /*	The underlying */
    int use_calib, /*	0: use lgm2dund, 1: calibrate */

    /*		if calib */
    char*  yc,         /*	yc */
    char*  vc,         /*	vc (only if calib) */
    char*  ref,        /*	ref rate (only if calib) */
    char*  swap_freq,  /*	swap freq (only if calib) */
    char*  swap_basis, /*	swap basis (only if calib) */
    double lambda,     /*	lambda if unique */
    double alpha,      /*	alpha */
    double gamma,      /*	gamma */
    double rho,        /*	rho */
    /*	End of calib params */

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */

    /*		funding */
    char*   fund_ref,
    double* fund_not,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,

    /*		cf */
    double* fix_not,
    int     fix_ncpn,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee,

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    double mintime,
    double mininterval,

    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    Err (*get_correl)(
        char* vol_curve_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,

    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                               to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Results */
    AM_STR am,
    AM_UND und,
    int*   call_feat, /*	0: No callable feature to be valued
                              1: Callable feature to be valued through adi */
    AM_ADI_ARG adi_arg)
{
    Err  szErr = NULL;
    char swap_freq_loc[256];
    char swap_basis_loc[256];
    int  intfreq;

    /*	Initialisation */
    am->fund_leg    = NULL;
    am->fix_leg     = NULL;
    am->call        = NULL;
    am->theoEndDate = theoEndDate;

    und->sigma_n    = 0;
    und->sigma_date = NULL;
    und->sigma_time = NULL;
    und->sigma      = NULL;

    strcpy(swap_basis_loc, fix_basis[0]);

    intfreq = (int)(12 * (fix_pay[fix_ncpn - 1] - fix_start[fix_ncpn - 1]) / 365.0 + 0.5);
    if (intfreq == 1)
    {
        strcpy(swap_freq_loc, "M");
    }
    else if (intfreq == 3)
    {
        strcpy(swap_freq_loc, "Q");
    }
    else if (intfreq == 6)
    {
        strcpy(swap_freq_loc, "S");
    }
    else if (intfreq == 12)
    {
        strcpy(swap_freq_loc, "A");
    }

    cpd_init_calib_inst_data(&(und->inst_data));
    und->has_inst_data = 0;
    und->today         = today;

    adi_arg->time     = NULL;
    adi_arg->date     = NULL;
    adi_arg->sig1     = NULL;
    adi_arg->sig_time = NULL;
    adi_arg->void_prm = NULL;
    adi_arg->is_event = NULL;
    adi_arg->ifr      = NULL;

    /*	Funding leg */

    am->fund_leg = (AM_FUND_LEG)malloc(sizeof(am_fund_leg));
    if (!am->fund_leg)
    {
        szErr = "Memory allocation error (1) in am_fill_check_all_AmortMidatAutocal";
        goto FREE_RETURN;
    }

    szErr = am_fill_fund_leg(
        und->today,
        eod_fix_flag,
        fund_not,
        fund_ncpn,
        fund_fix,
        fund_start,
        fund_end,
        fund_pay,
        fund_basis,
        fund_spr,
        fund_mrg,
        am->fund_leg);
    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	fix leg */

    am->fix_leg = (AM_FIX_LEG)malloc(sizeof(am_fix_leg));
    if (!am->fix_leg)
    {
        szErr = "Memory allocation error (2) in am_fill_check_all_struct_AmortMidatAutocal";
        goto FREE_RETURN;
    }

    szErr = am_fill_fix_leg(
        und->today,
        eod_fix_flag,
        fix_not,
        fix_ncpn,
        //			fix_fix,
        fix_start,
        fix_end,
        fix_pay,
        fix_basis,
        fix_rate,
        fix_fee,
        am->fix_leg);

    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	Calls */

    if (ncall > 0 && ex_date[ncall - 1] >= und->today + eod_ex_flag)
    {
        szErr = am_fill_calls(und->today, eod_ex_flag, ncall, pay_rec, ex_date, set_date, fee, am);
        if (szErr)
        {
            goto FREE_RETURN;
        }
    }
    else
    {
        am->num_calls = 0;
        am->call      = NULL;
    }

    /*	Underlying */
    if (am->num_calls > 0 && am->call)
    {
        if (use_calib)
        {
            szErr = am_calib_und_new_reduc(
                today,
                eod_ex_flag,

                yc,
                vc,
                ref,
                swap_freq,
                swap_basis,

                fund_ref,
                swap_freq_loc,
                swap_basis_loc,
                lambda,
                alpha,
                gamma,
                rho,

                mintime,
                mininterval,

                notperiod,
                max_std_short,
                one2F,
                fix_lambda,

                one_f_equi,

                skip_last,

                use_jump,
                max_var_jump,
                strike_type,
                european_model,

                get_correl,
                CorrelName,
                /*
                                                                GetVolForBadr,
                                                                cVolType,
                */
                am,

                get_cash_vol,

                und);
        }
        else
        {
            szErr = am_fill_und(lgm2dund, vc, ref, swap_freq, swap_basis, und);
        }
    }

    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	Tree */

    if (am->num_calls > 0 && am->call)
    {
        szErr = am_fill_adi_arg(und, am, get_cash_vol, req_stp, req_stpx, adi_arg);
        if (szErr)
        {
            goto FREE_RETURN;
        }

        *call_feat = 1;
    }
    else
    {
        *call_feat = 0;
    }

FREE_RETURN:

    if (szErr)
    {
        am_free_all_struct(am, und, *call_feat, adi_arg);
    }

    return szErr;
}

Err am_fill_check_all_struct_ts_AmortMidatAutocal(
    /*	Today's date */
    long today,
    long theoEndDate,

    /*	The underlying */
    int use_calib, /*	0: use lgm2dund, 1: calibrate */

    /*		if calib */
    char* yc,         /*	yc */
    char* vc,         /*	vc (only if calib) */
    char* ref,        /*	ref rate (only if calib) */
    char* swap_freq,  /*	swap freq (only if calib) */
    char* swap_basis, /*	swap basis (only if calib) */

    int     nlambda,
    double* pdlambda_time,
    double* pdlambda,

    double alpha, /*	alpha */
    double gamma, /*	gamma */
    double rho,   /*	rho */
    /*	End of calib params */

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */

    /*		funding */
    char*   fund_ref,
    double* fund_not,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,

    /*		cf */
    double* fix_not,
    int     fix_ncpn,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee,

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    double mintime,
    double mininterval,

    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    Err (*get_correl)(
        char* vol_curve_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,

    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                               to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Results */
    AM_STR am,
    AM_UND und,
    int*   call_feat, /*	0: No callable feature to be valued
                              1: Callable feature to be valued through adi */
    AM_ADI_ARG adi_arg)
{
    Err  szErr = NULL;
    char swap_freq_loc[256];
    char swap_basis_loc[256];
    int  intfreq;

    /*	Initialisation */
    am->fund_leg    = NULL;
    am->fix_leg     = NULL;
    am->call        = NULL;
    am->theoEndDate = theoEndDate;

    und->sigma_n    = 0;
    und->sigma_date = NULL;
    und->sigma_time = NULL;
    und->sigma      = NULL;

    strcpy(swap_basis_loc, fix_basis[0]);

    intfreq = (int)(12 * (fix_pay[fix_ncpn - 1] - fix_start[fix_ncpn - 1]) / 365.0 + 0.5);
    if (intfreq == 1)
    {
        strcpy(swap_freq_loc, "M");
    }
    else if (intfreq == 3)
    {
        strcpy(swap_freq_loc, "Q");
    }
    else if (intfreq == 6)
    {
        strcpy(swap_freq_loc, "S");
    }
    else if (intfreq == 12)
    {
        strcpy(swap_freq_loc, "A");
    }

    cpd_init_calib_inst_data(&(und->inst_data));
    und->has_inst_data = 0;
    und->today         = today;

    adi_arg->time     = NULL;
    adi_arg->date     = NULL;
    adi_arg->sig1     = NULL;
    adi_arg->sig_time = NULL;
    adi_arg->void_prm = NULL;
    adi_arg->is_event = NULL;
    adi_arg->ifr      = NULL;

    /*	Funding leg */

    am->fund_leg = (AM_FUND_LEG)malloc(sizeof(am_fund_leg));
    if (!am->fund_leg)
    {
        szErr = "Memory allocation error (1) in am_fill_check_all_AmortMidatAutocal";
        goto FREE_RETURN;
    }

    szErr = am_fill_fund_leg(
        und->today,
        eod_fix_flag,
        fund_not,
        fund_ncpn,
        fund_fix,
        fund_start,
        fund_end,
        fund_pay,
        fund_basis,
        fund_spr,
        fund_mrg,
        am->fund_leg);
    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	fix leg */

    am->fix_leg = (AM_FIX_LEG)malloc(sizeof(am_fix_leg));
    if (!am->fix_leg)
    {
        szErr = "Memory allocation error (2) in am_fill_check_all_AmortMidatAutocal";
        goto FREE_RETURN;
    }

    szErr = am_fill_fix_leg(
        und->today,
        eod_fix_flag,
        fix_not,
        fix_ncpn,
        //			fix_fix,
        fix_start,
        fix_end,
        fix_pay,
        fix_basis,
        fix_rate,
        fix_fee,
        am->fix_leg);

    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	Calls */

    if (ncall > 0 && ex_date[ncall - 1] >= und->today + eod_ex_flag)
    {
        szErr = am_fill_calls(und->today, eod_ex_flag, ncall, pay_rec, ex_date, set_date, fee, am);
        if (szErr)
        {
            goto FREE_RETURN;
        }
    }
    else
    {
        am->num_calls = 0;
        am->call      = NULL;
    }

    /*	Underlying */
    if (am->num_calls > 0 && am->call)
    {
        if (use_calib)
        {
            szErr = am_calib_und_new_ts(
                today,
                eod_ex_flag,

                yc,
                vc,
                ref,
                swap_freq,
                swap_basis,

                fund_ref,
                swap_freq_loc,
                swap_basis_loc,
                nlambda,
                pdlambda_time,
                pdlambda,
                alpha,
                gamma,
                rho,

                mintime,
                mininterval,

                notperiod,
                max_std_short,
                one2F,
                fix_lambda,

                one_f_equi,

                skip_last,

                use_jump,
                max_var_jump,
                strike_type,
                european_model,

                get_correl,
                CorrelName,
                /*
                                                                GetVolForBadr,
                                                                cVolType,
                */
                am,

                get_cash_vol,

                und);
        }
        else
        {
            szErr = am_fill_und(lgm2dund, vc, ref, swap_freq, swap_basis, und);
        }
    }

    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	Tree */

    if (am->num_calls > 0 && am->call)
    {
        szErr = am_fill_adi_arg_ts(
            nlambda, pdlambda, pdlambda_time, und, am, get_cash_vol, req_stp, req_stpx, adi_arg);
        if (szErr)
        {
            goto FREE_RETURN;
        }

        *call_feat = 1;
    }
    else
    {
        *call_feat = 0;
    }

FREE_RETURN:

    if (szErr)
    {
        am_free_all_struct(am, und, *call_feat, adi_arg);
    }

    return szErr;
}

Err am_fill_check_all_struct_ts_AmortMidatAutocal_reduc(
    /*	Today's date */
    long today,
    long theoEndDate,

    /*	The underlying */
    int use_calib, /*	0: use lgm2dund, 1: calibrate */

    /*		if calib */
    char* yc,         /*	yc */
    char* vc,         /*	vc (only if calib) */
    char* ref,        /*	ref rate (only if calib) */
    char* swap_freq,  /*	swap freq (only if calib) */
    char* swap_basis, /*	swap basis (only if calib) */

    int     nlambda,
    double* pdlambda_time,
    double* pdlambda,

    double alpha, /*	alpha */
    double gamma, /*	gamma */
    double rho,   /*	rho */
    /*	End of calib params */

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */

    /*		funding */
    char*   fund_ref,
    double* fund_not,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,

    /*		cf */
    double* fix_not,
    int     fix_ncpn,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee,

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    double mintime,
    double mininterval,

    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    Err (*get_correl)(
        char* vol_curve_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,

    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                               to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Results */
    AM_STR am,
    AM_UND und,
    int*   call_feat, /*	0: No callable feature to be valued
                              1: Callable feature to be valued through adi */
    AM_ADI_ARG adi_arg)
{
    Err  szErr = NULL;
    char swap_freq_loc[256];
    char swap_basis_loc[256];
    int  intfreq;

    /*	Initialisation */
    am->fund_leg    = NULL;
    am->fix_leg     = NULL;
    am->call        = NULL;
    am->theoEndDate = theoEndDate;

    und->sigma_n    = 0;
    und->sigma_date = NULL;
    und->sigma_time = NULL;
    und->sigma      = NULL;

    strcpy(swap_basis_loc, fix_basis[0]);

    intfreq = (int)(12 * (fix_pay[fix_ncpn - 1] - fix_start[fix_ncpn - 1]) / 365.0 + 0.5);
    if (intfreq == 1)
    {
        strcpy(swap_freq_loc, "M");
    }
    else if (intfreq == 3)
    {
        strcpy(swap_freq_loc, "Q");
    }
    else if (intfreq == 6)
    {
        strcpy(swap_freq_loc, "S");
    }
    else if (intfreq == 12)
    {
        strcpy(swap_freq_loc, "A");
    }

    cpd_init_calib_inst_data(&(und->inst_data));
    und->has_inst_data = 0;
    und->today         = today;

    adi_arg->time     = NULL;
    adi_arg->date     = NULL;
    adi_arg->sig1     = NULL;
    adi_arg->sig_time = NULL;
    adi_arg->void_prm = NULL;
    adi_arg->is_event = NULL;
    adi_arg->ifr      = NULL;

    /*	Funding leg */

    am->fund_leg = (AM_FUND_LEG)malloc(sizeof(am_fund_leg));
    if (!am->fund_leg)
    {
        szErr = "Memory allocation error (1) in am_fill_check_all_AmortMidatAutocal";
        goto FREE_RETURN;
    }

    szErr = am_fill_fund_leg(
        und->today,
        eod_fix_flag,
        fund_not,
        fund_ncpn,
        fund_fix,
        fund_start,
        fund_end,
        fund_pay,
        fund_basis,
        fund_spr,
        fund_mrg,
        am->fund_leg);
    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	fix leg */

    am->fix_leg = (AM_FIX_LEG)malloc(sizeof(am_fix_leg));
    if (!am->fix_leg)
    {
        szErr = "Memory allocation error (2) in am_fill_check_all_AmortMidatAutocal";
        goto FREE_RETURN;
    }

    szErr = am_fill_fix_leg(
        und->today,
        eod_fix_flag,
        fix_not,
        fix_ncpn,
        //			fix_fix,
        fix_start,
        fix_end,
        fix_pay,
        fix_basis,
        fix_rate,
        fix_fee,
        am->fix_leg);

    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	Calls */

    if (ncall > 0 && ex_date[ncall - 1] >= und->today + eod_ex_flag)
    {
        szErr = am_fill_calls(und->today, eod_ex_flag, ncall, pay_rec, ex_date, set_date, fee, am);
        if (szErr)
        {
            goto FREE_RETURN;
        }
    }
    else
    {
        am->num_calls = 0;
        am->call      = NULL;
    }

    /*	Underlying */
    if (am->num_calls > 0 && am->call)
    {
        if (use_calib)
        {
            szErr = am_calib_und_new_ts(
                today,
                eod_ex_flag,

                yc,
                vc,
                ref,
                swap_freq,
                swap_basis,

                fund_ref,
                swap_freq_loc,
                swap_basis_loc,
                nlambda,
                pdlambda_time,
                pdlambda,
                alpha,
                gamma,
                rho,

                mintime,
                mininterval,

                notperiod,
                max_std_short,
                one2F,
                fix_lambda,

                one_f_equi,

                skip_last,

                use_jump,
                max_var_jump,
                strike_type,
                european_model,

                get_correl,
                CorrelName,
                /*
                                                                GetVolForBadr,
                                                                cVolType,
                */
                am,

                get_cash_vol,

                und);
        }
        else
        {
            szErr = am_fill_und(lgm2dund, vc, ref, swap_freq, swap_basis, und);
        }
    }

    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	Tree */

    if (am->num_calls > 0 && am->call)
    {
        szErr = am_fill_adi_arg_ts(
            nlambda, pdlambda, pdlambda_time, und, am, get_cash_vol, req_stp, req_stpx, adi_arg);
        if (szErr)
        {
            goto FREE_RETURN;
        }

        *call_feat = 1;
    }
    else
    {
        *call_feat = 0;
    }

FREE_RETURN:

    if (szErr)
    {
        am_free_all_struct(am, und, *call_feat, adi_arg);
    }

    return szErr;
}

//// optimized payoff function
//// that utilizes MidAt "put-call parity"
Err _genmidat_nocall_payoff_(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    void*   yc,
    double* lam,
    double* ts_time,
    int     nb_ts,
    double  gamma,
    double  rho,
    double  phi1,
    double  phi2,
    double  phi12,
    /* Nodes data */
    int      l1,
    int      u1,
    int      l2,
    int      u2,
    double*  r1,
    double** r2,
    int      nprod,
    /* Vector of results to be updated */
    double*** prod_val)
{
    /// alias
    const ptr_genmidat_arg func_params    = (ptr_genmidat_arg)func_parm;
    const AM_PAY_ARG       am_arg         = (AM_PAY_ARG)(func_params->m_pBase);
    const int              call_idx       = am_arg->call_idx;
    const AM_EVAL_CONST    eval_const     = (AM_EVAL_CONST)(&(am_arg->eval_const));
    const AM_STR           am             = am_arg->am;
    const AM_CALL          call           = am->call + call_idx;
    const AM_UND           und            = (AM_UND)(am_arg->und);
    const AM_FUND_LEG      fund_leg       = am->fund_leg;
    const AM_FIX_LEG       fix_leg        = am->fix_leg;
    const int              nCompExProb    = func_params->nCompExProb;
    double **              ppdExIndicator = 0, *pdAvgSwapRate = 0;
    int                    nNum_Boundary = 0;

    const double dDelta_R1 = r1[l1 + 1] - r1[l1];
    const double _R1       = r1[l1] - dDelta_R1;  // offset by 1
    const double dDelta_R2 = r2[l1][l2 + 1] - r2[l1][l2];
    double       _R2;  // = r2[l1][l2] - dDelta_R2; // offset by 1

    AM_FUND_CPN fund_cpn;
    AM_FIX_CPN  fix_cpn;

    int    i, j, l, ndf;
    double fund_leg_pv, fix_leg_pv, iv;
    double fee, fee_total;

    double dFltInitCoeff_R1, dFltInitCoeff_R2, dFltInitBase, dFltInitPayment;
    double dFltTermCoeff_R1, dFltTermCoeff_R2, dFltTermBase, dFltTermPayment;
    double dCallFeeCoeff_R1, dCallFeeCoeff_R2, dCallFeeBase;
    double dFixFeeBase, dFixFee;

    //// NB: Optimzation in the time direction ; compute DF only between 2 exercise dates
    /// total number of funding leg coupons
    const int num_fund_cpn = call->num_fund_cpn - (*func_params->pnum_fund_cpn_done);
    const int num_fix_cpn  = call->num_fix_cpn - (*func_params->pnum_fix_cpn_done);

    // allocate mem
    double* pdFixCoeff_R1 = _alloca(num_fix_cpn * sizeof(double));
    double* pdFixCoeff_R2 = _alloca(num_fix_cpn * sizeof(double));
    double* pdFltCoeff_R1 = _alloca(num_fund_cpn * sizeof(double));
    double* pdFltCoeff_R2 = _alloca(num_fund_cpn * sizeof(double));
    double* pdFixPayment  = _alloca(num_fix_cpn * sizeof(double));
    double* pdFltPayment  = _alloca(num_fund_cpn * sizeof(double));
    double* pdFixBase     = _alloca(num_fix_cpn * sizeof(double));
    double* pdFltBase     = _alloca(num_fund_cpn * sizeof(double));

    // precompute constants for fund_coupons
    for (l = 0; l < num_fund_cpn; l++)
    {
        ndf              = eval_const->fund_idx[l];
        pdFltCoeff_R1[l] = exp(-eval_const->df_beta[ndf] * dDelta_R1);
        pdFltCoeff_R2[l] = exp(-eval_const->df_gamma[ndf] * dDelta_R2);
        pdFltBase[l]     = exp(-eval_const->df_alpha[ndf] - eval_const->df_beta[ndf] * _R1);
        fund_cpn         = am->fund_leg->cpn + call->fund_idx + l;
        pdFltBase[l] *= fund_cpn->cpn;
    }

    // precompute constants for "fund initial"
    ndf              = eval_const->start_idx;
    dFltInitCoeff_R1 = exp(-eval_const->df_beta[ndf] * dDelta_R1);
    dFltInitCoeff_R2 = exp(-eval_const->df_gamma[ndf] * dDelta_R2);
    dFltInitBase     = exp(-eval_const->df_alpha[ndf] - eval_const->df_beta[ndf] * _R1);
    dFltInitBase *= fund_leg->notional[call->fund_idx];

    // precompute constants for "fund terminal"
    ndf              = eval_const->fund_idx[num_fund_cpn - 1];
    dFltTermCoeff_R1 = exp(-eval_const->df_beta[ndf] * dDelta_R1);
    dFltTermCoeff_R2 = exp(-eval_const->df_gamma[ndf] * dDelta_R2);
    dFltTermBase     = exp(-eval_const->df_alpha[ndf] - eval_const->df_beta[ndf] * _R1);
    dFltTermBase *= (*func_params->pdInit_Prev);

    // precompute constants for "fix fee"
    ndf         = eval_const->start_idx;
    dFixFeeBase = exp(-eval_const->df_alpha[ndf] - eval_const->df_beta[ndf] * _R1);
    dFixFeeBase *= fix_leg->fee[call->fix_idx];
    if (call->pay_rec)       /// factor Fix fee into call fee
        dFixFeeBase *= -1.;  /// receiving floating

    // precompute constants for call fee
    ndf              = eval_const->fee_idx;
    dCallFeeCoeff_R1 = exp(-eval_const->df_beta[ndf] * dDelta_R1);
    dCallFeeCoeff_R2 = exp(-eval_const->df_gamma[ndf] * dDelta_R2);
    dCallFeeBase     = exp(-eval_const->df_alpha[ndf] - eval_const->df_beta[ndf] * _R1);

    dCallFeeBase *= (fabs(call->fee) > 1.0e-08 ? call->fee : 0.);

    // precompute constants for fund_coupons
    for (l = 0; l < num_fix_cpn; l++)
    {
        ndf              = eval_const->fix_disc_idx[l];
        pdFixCoeff_R1[l] = exp(-eval_const->df_beta[ndf] * dDelta_R1);
        pdFixCoeff_R2[l] = exp(-eval_const->df_gamma[ndf] * dDelta_R2);
        pdFixBase[l]     = exp(-eval_const->df_alpha[ndf] - eval_const->df_beta[ndf] * _R1);

        fix_cpn = fix_leg->cpn + call->fix_idx + l;
        pdFixBase[l] *= fix_cpn->cpn;
    }

    //	Eval payoff
#ifdef _DEBUG
    _RPT0(_CRT_WARN, "\n _genmidat_nocall_payoff_(...): iv \n");
#endif  ///_DEBUG

    // allocate mem for exercie indicators
    if (nCompExProb)
    {
        *func_params->pdUB_D1 = u1;
        *func_params->pdLB_D1 = l1;
        *func_params->pdUB_D2 = u2;
        *func_params->pdLB_D2 = l2;

        *func_params->pppdExIndicator = dmatrix(l1, u1, l2, u2);
        ppdExIndicator                = *func_params->pppdExIndicator;
        pdAvgSwapRate                 = func_params->pdAvgSwapRate;
        _ASSERTE(pdAvgSwapRate != 0);
        *pdAvgSwapRate = 0.;
    }

    for (i = l1; i <= u1; i++)
    {
        _R2 = r2[i][l2] - dDelta_R2;

        // update funding legs
        for (l = 0; l < num_fund_cpn; l++)
        {
            ndf             = eval_const->fund_idx[l];
            pdFltPayment[l] = (pdFltBase[l] *= pdFltCoeff_R1[l]);
            pdFltPayment[l] *= exp(-eval_const->df_gamma[ndf] * _R2);
        }

        /// fuding leg init
        ndf             = eval_const->start_idx;
        dFltInitPayment = (dFltInitBase *= dFltInitCoeff_R1);
        dFltInitPayment *= exp(-eval_const->df_gamma[ndf] * _R2);

        /// funding leg terminal
        ndf             = eval_const->fund_idx[num_fund_cpn - 1];
        dFltTermPayment = (dFltTermBase *= dFltTermCoeff_R1);
        dFltTermPayment *= exp(-eval_const->df_gamma[ndf] * _R2);

        for (l = 0; l < num_fix_cpn; l++)
        {
            ndf             = eval_const->fix_disc_idx[l];
            pdFixPayment[l] = (pdFixBase[l] *= pdFixCoeff_R1[l]);
            pdFixPayment[l] *= exp(-eval_const->df_gamma[ndf] * _R2);
        }

        /// fix fee
        ndf     = eval_const->start_idx;
        dFixFee = (dFixFeeBase *= dFltInitCoeff_R1);
        dFixFee *= exp(-eval_const->df_gamma[ndf] * _R2);

        // call fee
        ndf = eval_const->fee_idx;
        fee = (dCallFeeBase *= dCallFeeCoeff_R1);
        fee *= exp(-eval_const->df_gamma[ndf] * _R2);

        for (j = l2; j <= u2; j++)
        {
            //	Funding leg
            fund_leg_pv = (dFltInitPayment *= dFltInitCoeff_R2);
            fund_leg_pv -= (dFltTermPayment *= dFltTermCoeff_R2);
            for (l = 0; l < num_fund_cpn; l++)
                fund_leg_pv += (pdFltPayment[l] *= pdFltCoeff_R2[l]);

            // Fixed leg
            fix_leg_pv = 0.0;
            for (l = 0; l < num_fix_cpn; l++)
                fix_leg_pv += (pdFixPayment[l] *= pdFixCoeff_R2[l]);

            //	Intrinsic value
            iv = call->pay_rec ? (fund_leg_pv - fix_leg_pv) : (fix_leg_pv - fund_leg_pv);

            fee *= dCallFeeCoeff_R2;                          /// call fee
            fee_total = fee + (dFixFee *= dFltInitCoeff_R2);  // factor fixed fee into call fee

            // update indicator function - exercise 1. and not exercise 0.
            if (nCompExProb)
            {
                ppdExIndicator[i][j] = (prod_val[i][j][0] - iv) > (-fee_total) ? 0. : 1.;

                // detect flipping/boundary
                if (i > l1 && j > l2 && j < u2)
                {
                    const double dExInd = ppdExIndicator[i][j];
                    const double dA     = ppdExIndicator[i - 1][j - 1];
                    const double dB     = ppdExIndicator[i - 1][j];
                    const double dC     = ppdExIndicator[i - 1][j + 1];
                    const double dD     = ppdExIndicator[i][j - 1];

                    // if dExInd is different from the ajacent 4 cells (i-1,j-1),
                    // (i-1,j),(i-1,j+1),(i,j-1) then (i,j) sit right on the boundary
                    int nNumID = 0;
                    if (dExInd == dA)
                        ++nNumID;
                    if (dExInd == dB)
                        ++nNumID;
                    if (dExInd == dC)
                        ++nNumID;
                    if (dExInd == dD)
                        ++nNumID;
                    if (nNumID <= 2)
                    {
                        double dFixLeg = 0.;
                        double dFltLeg = 0.;
                        {
                            const int num_fund_cpn = call->num_fund_cpn;
                            const int num_fix_cpn  = call->num_fix_cpn;

                            double _dFltInitCoeff_R1, _dFltInitCoeff_R2, _dFltInitBase,
                                _dFltInitPayment;
                            double _dFltTermCoeff_R1, _dFltTermCoeff_R2, _dFltTermBase,
                                _dFltTermPayment;
                            const double __R1 = r1[i] - dDelta_R1;     // offset by 1
                            const double __R2 = r2[i][j] - dDelta_R2;  // offset by 1

                            double* _pdFixCoeff_R1 = _alloca(num_fix_cpn * sizeof(double));
                            double* _pdFixCoeff_R2 = _alloca(num_fix_cpn * sizeof(double));
                            double* _pdFltCoeff_R1 = _alloca(num_fund_cpn * sizeof(double));
                            double* _pdFltCoeff_R2 = _alloca(num_fund_cpn * sizeof(double));
                            double* _pdFixPayment  = _alloca(num_fix_cpn * sizeof(double));
                            double* _pdFltPayment  = _alloca(num_fund_cpn * sizeof(double));
                            double* _pdFixBase     = _alloca(num_fix_cpn * sizeof(double));
                            double* _pdFltBase     = _alloca(num_fund_cpn * sizeof(double));

                            // precompute constants for fund_coupons
                            for (l = 0; l < num_fund_cpn; l++)
                            {
                                ndf               = eval_const->fund_idx[l];
                                _pdFltCoeff_R1[l] = exp(-eval_const->df_beta[ndf] * dDelta_R1);
                                _pdFltCoeff_R2[l] = exp(-eval_const->df_gamma[ndf] * dDelta_R2);
                                _pdFltBase[l]     = exp(
                                    -eval_const->df_alpha[ndf] - eval_const->df_beta[ndf] * __R1);
                                _pdFltBase[l] *= (am->fund_leg->cpn + call->fund_idx + l)->cpn;
                            }

                            // precompute constants for "fund initial"
                            ndf               = eval_const->start_idx;
                            _dFltInitCoeff_R1 = exp(-eval_const->df_beta[ndf] * dDelta_R1);
                            _dFltInitCoeff_R2 = exp(-eval_const->df_gamma[ndf] * dDelta_R2);
                            _dFltInitBase =
                                exp(-eval_const->df_alpha[ndf] - eval_const->df_beta[ndf] * __R1);
                            _dFltInitBase *= fund_leg->notional[call->fund_idx];

                            // precompute constants for "fund terminal"
                            ndf               = eval_const->fund_idx[num_fund_cpn - 1];
                            _dFltTermCoeff_R1 = exp(-eval_const->df_beta[ndf] * dDelta_R1);
                            _dFltTermCoeff_R2 = exp(-eval_const->df_gamma[ndf] * dDelta_R2);
                            _dFltTermBase =
                                exp(-eval_const->df_alpha[ndf] - eval_const->df_beta[ndf] * __R1);
                            _dFltTermBase *= (*func_params->pdInit_Prev);

                            // precompute constants for fix coupons
                            for (l = 0; l < num_fix_cpn; l++)
                            {
                                ndf               = eval_const->fix_disc_idx[l];
                                _pdFixCoeff_R1[l] = exp(-eval_const->df_beta[ndf] * dDelta_R1);
                                _pdFixCoeff_R2[l] = exp(-eval_const->df_gamma[ndf] * dDelta_R2);
                                _pdFixBase[l]     = exp(
                                    -eval_const->df_alpha[ndf] - eval_const->df_beta[ndf] * __R1);

                                // fix_cpn = fix_leg->cpn + call->fix_idx + l;
                                _pdFixBase[l] *= (fix_leg->cpn + call->fix_idx + l)->cpn;
                            }

                            // update funding legs
                            for (l = 0; l < num_fund_cpn; l++)
                            {
                                ndf              = eval_const->fund_idx[l];
                                _pdFltPayment[l] = (_pdFltBase[l] *= _pdFltCoeff_R1[l]);
                                _pdFltPayment[l] *= exp(-eval_const->df_gamma[ndf] * __R2);
                            }

                            /// fuding leg init
                            ndf              = eval_const->start_idx;
                            _dFltInitPayment = (_dFltInitBase *= _dFltInitCoeff_R1);
                            _dFltInitPayment *= exp(-eval_const->df_gamma[ndf] * __R2);

                            /// funding leg terminal
                            ndf              = eval_const->fund_idx[num_fund_cpn - 1];
                            _dFltTermPayment = (_dFltTermBase *= _dFltTermCoeff_R1);
                            _dFltTermPayment *= exp(-eval_const->df_gamma[ndf] * __R2);

                            for (l = 0; l < num_fix_cpn; l++)
                            {
                                ndf              = eval_const->fix_disc_idx[l];
                                _pdFixPayment[l] = (_pdFixBase[l] *= _pdFixCoeff_R1[l]);
                                _pdFixPayment[l] *= exp(-eval_const->df_gamma[ndf] * __R2);
                            }
                            //	Funding leg
                            dFltLeg = (_dFltInitPayment *= _dFltInitCoeff_R2);
                            // dFltLeg -= (_dFltTermPayment*= _dFltTermCoeff_R2);
                            for (l = 0; l < num_fund_cpn; l++)
                                dFltLeg += (_pdFltPayment[l] *= _pdFltCoeff_R2[l]);

                            // Fixed leg
                            dFixLeg = 0.0;
                            for (l = 0; l < num_fix_cpn; l++)
                                dFixLeg += (_pdFixPayment[l] *= _pdFixCoeff_R2[l]);
                        }

                        // const double dSwapRate = fix_cpn->cpn* fund_leg_pv/fix_leg_pv;
                        *pdAvgSwapRate += dFltLeg / dFixLeg;
                        ++nNum_Boundary;
                    }
                }
            }

            prod_val[i][j][0] =
                (prod_val[i][j][0] - iv) > (-fee_total) ? (prod_val[i][j][0] - iv) : (-fee_total);
        }

        // check that state R2 flip no more than once given each R1
#ifdef _DEBUG
        if (nCompExProb)
        {
            const double* pdExInd = 0;
            int           nFlip   = 0;
            for (pdExInd = ppdExIndicator[i] + l2 + 1; pdExInd < ppdExIndicator[i] + u2 + 1;
                 ++pdExInd)
            {
                if (*pdExInd != *(pdExInd - 1))
                    ++nFlip;
            }
            _ASSERTE(nFlip <= 1);
        }
#endif
    }

    if (nCompExProb)
    {
        _ASSERTE(nNum_Boundary != 0);
        _ASSERTE(nNum_Boundary <= (u1 - l1 - 1) > (u2 - l2 - 2) ? (u1 - l1 - 1) : (u2 - l2 - 2));
        _ASSERTE(func_params->pdStrike != 0);

        *pdAvgSwapRate /= nNum_Boundary;
        *pdAvgSwapRate *= *func_params->pdStrike;
        ++func_params->pdAvgSwapRate;
    }

    //	update time optimzation variables
    (*func_params->pnum_fund_cpn_done) += num_fund_cpn;
    (*func_params->pnum_fix_cpn_done) += num_fix_cpn;
    (*func_params->pdInit_Prev) = fund_leg->notional[call->fund_idx];

    // if this is the first exercise date, update dIV
    if (call_idx == 0)
    {
        fund_leg_pv =
            swp_f_df(func_params->lToday, call->set_date, yc) * fund_leg->notional[call->fund_idx];

        //	Fund Spread Coupons */
        for (l = 0; l < call->num_fund_cpn; l++)
        {
            fund_cpn = am->fund_leg->cpn + call->fund_idx + l;
            fund_leg_pv += swp_f_df(func_params->lToday, fund_cpn->pay_date, yc) * fund_cpn->cpn;
        }

        //	PV of FIX leg
        fix_leg_pv = 0.0;
        for (l = 0; l < call->num_fix_cpn; l++)
        {
            fix_cpn = fix_leg->cpn + call->fix_idx + l;
            fix_leg_pv += swp_f_df(func_params->lToday, fix_cpn->pay_date, yc) * fix_cpn->cpn;
        }

        //	Intrinsic value
        if (call->pay_rec == 0)
        {
            (*func_params->pdIV) = fix_leg_pv - fund_leg_pv;
        }
        else
        {
            (*func_params->pdIV) = fund_leg_pv - fix_leg_pv;
        }
    }

    /// unreferenced variable
    nprod;
    phi12;
    phi1;
    phi2;
    rho;
    gamma;
    nb_ts;
    ts_time;
    lam;
    yc;
    evt_time;
    evt_date;

    return 0;
}

// copied
static double H_func(double lam, double t1, double t2)
{
    return exp(lam * t2) - exp(lam * t1);
}

// copied
static void LGM2FExpectations(
    int     nstept,
    double* time,
    double  lam1,
    double* sig_time,
    double* sig1,
    int     nb_sig,
    double  alpha,
    double  gamma,
    double  rho,
    double* fwd1,
    double* fwd1_mid,
    double* fwd3,
    double* fwd3_mid,
    double* var1,
    double* var2,
    double* phi1,
    double* phi2,
    double* phi12)
{
    double lam2, vol2, lam21, lam22, lam12;
    double t1, t2, ta, tb;
    double coef11, coef12, coef13, coef21, coef22, coef23, fact1, fact2, rhoalpha;

    double H1, H2, H21, H22, H12;
    double st1;

    int i, j, nb_sig_minus1;

    lam2          = lam1 + gamma;
    lam21         = 2.0 * lam1;
    lam22         = 2.0 * lam2;
    lam12         = (lam1 + lam2);
    rhoalpha      = rho * alpha;
    nb_sig_minus1 = nb_sig - 1;

    fact1 = 1.0 / alpha / sqrt(1.0 - rho * rho);
    fact2 = rho / sqrt(1.0 - rho * rho);

    /* constant for expectation r1_dim1 */
    coef11 = (lam2 + rho * alpha * lam1) / (lam1 * lam1 * lam2);
    coef12 = -1.0 / (2.0 * lam1 * lam1);
    coef13 = -rhoalpha / (lam2 * lam12);

    /* constant for expectation r3 */
    coef21 = alpha * (rho * lam2 + alpha * lam1) / (lam1 * lam2 * lam2);
    coef22 = -alpha * alpha / (2.0 * lam2 * lam2);
    coef23 = -rhoalpha / (lam1 * lam12);

    /* initialisation */
    H1 = H2 = H21 = H22 = H12 = 0.0;
    t1                        = 0.0;
    j                         = 0;
    vol2                      = sig1[0] * sig1[0];

    phi1[0]  = 0.0;
    phi2[0]  = 0.0;
    phi12[0] = 0.0;

    fwd1[0] = 0;
    fwd3[0] = 0;

    for (i = 1; i < nstept; i++)
    {
        t2 = 0.5 * (t1 + time[i]);

        /* first from t1 to tmid */

        ta = t1;
        tb = sig_time[j];

        st1 = 0.0;

        while (tb < t2 && j < nb_sig_minus1)
        {
            H1 += vol2 * H_func(lam1, ta, tb);
            H2 += vol2 * H_func(lam2, ta, tb);
            H21 += vol2 * H_func(lam21, ta, tb);
            H22 += vol2 * H_func(lam22, ta, tb);
            H12 += vol2 * H_func(lam12, ta, tb);
            st1 += vol2 * (tb - ta);

            j++;
            vol2 = sig1[j] * sig1[j];
            ta   = tb;
            tb   = sig_time[j];
        }

        H1 += vol2 * H_func(lam1, ta, t2);
        H2 += vol2 * H_func(lam2, ta, t2);
        H21 += vol2 * H_func(lam21, ta, t2);
        H22 += vol2 * H_func(lam22, ta, t2);
        H12 += vol2 * H_func(lam12, ta, t2);
        st1 += vol2 * (t2 - ta);

        var1[i - 1] = st1;
        var2[i - 1] = st1;

        phi1[i]  = exp(-lam21 * t2) * H21;
        phi2[i]  = exp(-lam22 * t2) * H22;
        phi12[i] = exp(-lam12 * t2) * H12;

        fwd1_mid[i - 1] = coef11 * exp(-lam1 * t2) * H1 + coef12 * phi1[i] + coef13 * phi12[i];
        fwd3_mid[i - 1] =
            fact1 * (coef21 * exp(-lam2 * t2) * H2 + coef22 * phi2[i] + coef23 * phi12[i]) -
            fact2 * fwd1_mid[i - 1];

        phi1[i] /= lam21;
        phi2[i] *= alpha * alpha / lam22;
        phi12[i] *= alpha / lam12;

        t1 = t2;

        t2 = time[i];

        /* first from tmid to t2 */

        ta = t1;
        tb = sig_time[j];

        while (tb < t2 && j < nb_sig_minus1)
        {
            H1 += vol2 * H_func(lam1, ta, tb);
            H2 += vol2 * H_func(lam2, ta, tb);
            H21 += vol2 * H_func(lam21, ta, tb);
            H22 += vol2 * H_func(lam22, ta, tb);
            H12 += vol2 * H_func(lam12, ta, tb);
            st1 += vol2 * (tb - ta);

            j++;
            vol2 = sig1[j] * sig1[j];
            ta   = tb;
            tb   = sig_time[j];
        }

        H1 += vol2 * H_func(lam1, ta, t2);
        H2 += vol2 * H_func(lam2, ta, t2);
        H21 += vol2 * H_func(lam21, ta, t2);
        H22 += vol2 * H_func(lam22, ta, t2);
        H12 += vol2 * H_func(lam12, ta, t2);
        st1 += vol2 * (t2 - ta);

        var1[i - 1] = st1;
        var2[i - 1] = st1;

        phi1[i]  = exp(-lam21 * t2) * H21;
        phi2[i]  = exp(-lam22 * t2) * H22;
        phi12[i] = exp(-lam12 * t2) * H12;

        fwd1[i] = coef11 * exp(-lam1 * t2) * H1 + coef12 * phi1[i] + coef13 * phi12[i];
        fwd3[i] = fact1 * (coef21 * exp(-lam2 * t2) * H2 + coef22 * phi2[i] + coef23 * phi12[i]) -
                  fact2 * fwd1[i];

        phi1[i] /= lam21;
        phi2[i] *= alpha * alpha / lam22;
        phi12[i] *= alpha / lam12;

        t1 = t2;
    }
}

//// payoff function used to compute exercise probabilities
static Err am_exerprob_4_lgm2f_adi_fixgrid(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    void*   yc,
    double* lam,
    double* ts_time,
    int     nb_ts,
    double  gamma,
    double  rho,
    double  phi1,
    double  phi2,
    double  phi12,
    /* Nodes data */
    int      l1,
    int      u1,
    int      l2,
    int      u2,
    double*  r1,
    double** r2,
    int      nprod,
    /* Vector of results to be updated */
    double*** prod_val)
{
    /// alias
    const ptr_genmidat_arg func_params    = (ptr_genmidat_arg)func_parm;
    const AM_PAY_ARG       am_arg         = (AM_PAY_ARG)(func_params->m_pBase);
    const int              call_idx       = am_arg->call_idx;
    const double**         ppdExIndicator = *(func_params->pppdExIndicator);
    const double**         ppdCondiDF     = *(func_params->pppdCondiDF);
    const int              nEx            = *(func_params->pnEx);
    int                    i, j;

    _ASSERTE(l1 == *func_params->pdLB_D1);
    _ASSERTE(u1 == *func_params->pdUB_D1);
    _ASSERTE(l2 == *func_params->pdLB_D2);
    _ASSERTE(u2 == *func_params->pdUB_D2);

#ifdef _DEBUG
    for (i = l1; i <= u1; i++)
        for (j = l2; j <= u2; j++)
            _ASSERTE(ppdExIndicator[i][j] == 1. || ppdExIndicator[i][j] == 0.);
#endif

    // if computing the probability of exercising at call_idxTH exercise date
    if (nEx == call_idx)
    {
        for (i = l1; i <= u1; i++)
        {
            for (j = l2; j <= u2; j++)
            {
                prod_val[i][j][0] = ppdExIndicator[i][j];
                /// TEST TEST TEST TEST TEST
                // prod_val[i][j][0] = 1.;
                /// TEST TEST TEST TEST TEST
            }
        }
    }
    else
    {
        for (i = l1; i <= u1; i++)
        {
            for (j = l2; j <= u2; j++)
            {
                _ASSERTE(ppdCondiDF[i][j] > 0.);
                prod_val[i][j][0] =
                    ppdExIndicator[i][j] == 1. ? 0. : prod_val[i][j][0] / ppdCondiDF[i][j];
                /// TEST TEST TEST TEST TEST
                // prod_val[i][j][0] = prod_val[i][j][0] / ppdCondiDF[i][j];
                /// TEST TEST TEST TEST TEST
            }
        }
    }

    /// unreferenced variable
    nprod;
    phi12;
    phi1;
    phi2;
    rho;
    gamma;
    nb_ts;
    ts_time;
    lam;
    yc;
    evt_time;
    evt_date;
    r1;
    r2;

    return 0;
}

//// payoff function used to compute exercise probabilities
static Err am_notexerprob_4_lgm2f_adi_fixgrid(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    void*   yc,
    double* lam,
    double* ts_time,
    int     nb_ts,
    double  gamma,
    double  rho,
    double  phi1,
    double  phi2,
    double  phi12,
    /* Nodes data */
    int      l1,
    int      u1,
    int      l2,
    int      u2,
    double*  r1,
    double** r2,
    int      nprod,
    /* Vector of results to be updated */
    double*** prod_val)
{
    /// alias
    const ptr_genmidat_arg func_params    = (ptr_genmidat_arg)func_parm;
    const AM_PAY_ARG       am_arg         = (AM_PAY_ARG)(func_params->m_pBase);
    const int              call_idx       = am_arg->call_idx;
    const double**         ppdExIndicator = *(func_params->pppdExIndicator);
    const double**         ppdCondiDF     = *(func_params->pppdCondiDF);
    const int              nEx            = *(func_params->pnEx);
    int                    i, j;

    _ASSERTE(l1 == *func_params->pdLB_D1);
    _ASSERTE(u1 == *func_params->pdUB_D1);
    _ASSERTE(l2 == *func_params->pdLB_D2);
    _ASSERTE(u2 == *func_params->pdUB_D2);

#ifdef _DEBUG
    for (i = l1; i <= u1; i++)
        for (j = l2; j <= u2; j++)
            _ASSERTE(ppdExIndicator[i][j] == 1. || ppdExIndicator[i][j] == 0.);
#endif

    // if computing the probability of exercising at call_idxTH exercise date
    if (nEx == call_idx)
    {
        for (i = l1; i <= u1; i++)
        {
            for (j = l2; j <= u2; j++)
            {
                prod_val[i][j][0] = 1. - ppdExIndicator[i][j];
                ///// TEST TEST TEST
                // prod_val[i][j][0] = 1.;
                ///// TEST TEST TEST
            }
        }
    }
    else
    {
        for (i = l1; i <= u1; i++)
        {
            for (j = l2; j <= u2; j++)
            {
                prod_val[i][j][0] =
                    (ppdExIndicator[i][j] == 1.) ? 0. : (prod_val[i][j][0] / ppdCondiDF[i][j]);
                ///// TEST TEST TEST
                /// prod_val[i][j][0] = prod_val[i][j][0]/ppdCondiDF[i][j];
                ///// TEST TEST TEST
            }
        }
    }

    /// unreferenced variable
    nprod;
    phi12;
    phi1;
    phi2;
    rho;
    gamma;
    nb_ts;
    ts_time;
    lam;
    yc;
    evt_time;
    evt_date;
    r1;
    r2;

    return 0;
}

//// pay off function used to compute conditional discount factors
Err am_CondiDF_4_lgm2f_adi_fixgrid(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    void*   yc,
    double* lam,
    double* ts_time,
    int     nb_ts,
    double  gamma,
    double  rho,
    double  phi1,
    double  phi2,
    double  phi12,
    /* Nodes data */
    int      l1,
    int      u1,
    int      l2,
    int      u2,
    double*  r1,
    double** r2,
    int      nprod,
    /* Vector of results to be updated */
    double*** prod_val)
{
    /// alias
    const ptr_genmidat_arg func_params = (ptr_genmidat_arg)func_parm;
    double**               ppdCondiDF  = 0;
    int                    i, j;

    _ASSERTE(l1 == *func_params->pdLB_D1);
    _ASSERTE(u1 == *func_params->pdUB_D1);
    _ASSERTE(l2 == *func_params->pdLB_D2);
    _ASSERTE(u2 == *func_params->pdUB_D2);

    (*func_params->pppdCondiDF) = dmatrix(l1, u1, l2, u2);
    ppdCondiDF                  = (*func_params->pppdCondiDF);

    for (i = l1; i <= u1; i++)
    {
        for (j = l2; j <= u2; j++)
        {
            /// store result of df
            ppdCondiDF[i][j] = prod_val[i][j][0];
            /// overwrite to compute conditional df at the next time slice
            prod_val[i][j][0] = 1.;
        }
    }

    /// unreferenced variable
    nprod;
    phi12;
    phi1;
    phi2;
    rho;
    gamma;
    nb_ts;
    ts_time;
    lam;
    yc;
    evt_time;
    evt_date;
    r1;
    r2;

    return 0;
}

static Err _lgm2f_adi_shrinkgrid(
    const int* plx,
    const int* pux,
    const int* plz,
    const int* puz,
    /*	Time data		*/
    int     nstept,
    double* time,
    double* date,

    /*	Discretisation	*/
    int nsteps,

    /*	Model data		*/
    double  lam1,
    double* sig_time,
    double* sig1,
    int     nb_sig,
    double  alpha,
    double  gamma,
    double  rho,

    /*	Product data */
    void** func_parm_tab,
    int*   eval_evt,

    /*	Market data */
    double* ifr,
    char*   yc,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date,
                       double evt_time,
                       void*  func_parm,

                       /* Market data	*/
                       void* yc,

                       /* Model data	*/
                       double* lam,
                       double* ts_time,
                       int     nb_ts,
                       double  gamma,
                       double  rho,
                       double  phi1,
                       double  phi2,
                       double  phi12,

                       /* Gride data	*/
                       int      l1,
                       int      u1,
                       int      l2,
                       int      u2,
                       double*  r1_dim2,
                       double** r2,

                       /* Vector of results to be updated */
                       int       nprod,
                       double*** prod_val),
    /*	Result */
    int     nprod,
    double* res)
{
    Err szErr = NULL;
    int i, j, step;
    // int				nstepx, nstepz;
    double   meshx, meshz, dt;
    double   mu_r1i, mu_r3i;
    double   r_temp;
    double   std1, std3, r2_temp;
    double **r1_dim2 = NULL, **r2 = NULL, ***values = NULL, ***values_p1 = NULL,
           ***values_temp = NULL, **mu_r1 = NULL, **mu_r3 = NULL, **var_r1 = NULL, **var_r3 = NULL,
           **r = NULL, **r_init = NULL;

    CNPDE_TEMP_2D_ADI pdestr, *pde = &pdestr;

    const int nstepx = pux[nstept - 1] + 1;
    const int nstepz = puz[nstept - 1] + 1;

    //	corresponding index to the 0 value of x and z
    const int index_x = (nstepx - 1) / 2;
    const int index_z = (nstepz - 1) / 2;

    /*
    lx = 0;
    ux = nstepx - 1;
    lz = 0;
    uz = nstepz - 1;
    */

    //	Allocations	of time vectors
    double* phi1     = _alloca(nstept * sizeof(double));
    double* phi2     = _alloca(nstept * sizeof(double));
    double* phi12    = _alloca(nstept * sizeof(double));
    double* var1     = _alloca(nstept * sizeof(double));
    double* var3     = _alloca(nstept * sizeof(double));
    double* fwd1     = _alloca(nstept * sizeof(double));
    double* fwd3     = _alloca(nstept * sizeof(double));
    double* fwd1_mid = _alloca(nstept * sizeof(double));
    double* fwd3_mid = _alloca(nstept * sizeof(double));
    double *r1_bar, *r3_bar, *r1_dim1;

    // Constant
    const double lam2           = lam1 + gamma;
    const double rho2           = sqrt(1.0 - rho * rho);
    const double rhoalpha       = rho * alpha;
    const double rhoalpha_plus1 = 1.0 + rhoalpha;
    const double rho2alpha      = alpha * rho2;
    const double fact1          = -rho * gamma / rho2;
    const double rhoq           = rho / rho2;

    if (!pde)
    {
        szErr = "Memory allocation error (1) in lgm2fpde";
        goto FREE_RETURN;
    }

    /*	Calculate the expecations of r1 and r3 */
    LGM2FExpectations(
        nstept,
        time,
        lam1,
        sig_time,
        sig1,
        nb_sig,
        alpha,
        gamma,
        rho,
        fwd1,
        fwd1_mid,
        fwd3,
        fwd3_mid,
        var1,
        var3,
        phi1,
        phi2,
        phi12);

    //	Calculation of the number of steps in each direction: since local volatility
    //	is the same, the mesh has to be the same, but the number of steps has to be adjusted
    std1 = sqrt(phi1[nstept - 1]);
    std3 = 1.0 / rho2 *
           sqrt(
               phi2[nstept - 1] / alpha / alpha + rho * rho * phi1[nstept - 1] -
               2.0 * rho * rho / alpha * phi12[nstept - 1]);

#if 0
		nstepx = (int) (nsteps * sqrt(std1 / std3) + 0.5);
		nstepz = (int) (nsteps * sqrt(std3 / std1) + 0.5);
	
		//	nstep has to be a odd nuber		
		nstepx = ((int) (nstepx / 2)) * 2 + 1;
		nstepz = ((int) (nstepz / 2)) * 2 + 1;
	
	
	//	we want at least three points in each directions										*/
	if (nstepx < 3)
	{
		nstepx = 3;
		nstepz = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepz = ((int) (nstepz / 2)) * 2 + 1;

		if (nstepz < 3)
		{
			nstepz = 3;		
		}
	}

	if (nstepz < 3)
	{
		nstepz = 3;
		nstepx = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepx = ((int) (nstepx / 2)) * 2 + 1;

		if (nstepx < 3)
		{
			nstepx = 3;		
		}
	}
#endif

    meshx = 2.0 * NSTD_LGM * sqrt(phi1[nstept - 1]) / (nstepx - 1);
    meshz = 2.0 * NSTD_LGM * std3 / (nstepz - 1);

    //	Allocations of space vectors
    r1_bar  = _alloca(nstepx * sizeof(double));
    r3_bar  = _alloca(nstepz * sizeof(double));
    r1_dim1 = _alloca(nstepx * sizeof(double));

    mu_r1   = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    mu_r3   = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    var_r1  = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    var_r3  = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    r       = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    r_init  = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    r1_dim2 = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    r2      = dmatrix(0, nstepx - 1, 0, nstepz - 1);

    values    = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
    values_p1 = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);

    if (!r1_dim2 || !r2 || !values || !values_p1 || !mu_r1 || !mu_r3 || !var_r1 || !var_r3 || !r ||
        !r_init || !r1_bar || !r1_dim1 || !r3_bar)
    {
        szErr = "Memory allocation error (2) in lgm2fpde";
        goto FREE_RETURN;
    }

    //	Then discretise space in the orthogonal system r1 / r3
    r1_bar[0] = -(nstepx - 1) / 2.0 * meshx;
    r3_bar[0] = -(nstepz - 1) / 2.0 * meshz;

    for (i = 1; i < nstepx; i++)
    {
        r1_bar[i] = r1_bar[i - 1] + meshx;
    }

    for (j = 1; j < nstepz; j++)
    {
        r3_bar[j] = r3_bar[j - 1] + meshz;
    }

    //	Corresponding r1 and r2
    for (i = 0; i < nstepx; i++)
    {
        r1_dim1[i] = r1_bar[i] + fwd1[nstept - 1];
        r2_temp    = rhoalpha * r1_dim1[i] + rho2alpha * fwd3[nstept - 1];
        r_temp     = rhoalpha_plus1 * r1_bar[i];

        for (j = 0; j < nstepz; j++)
        {
            r1_dim2[i][j] = r1_dim1[i];
            r2[i][j]      = r2_temp + rho2alpha * r3_bar[j];
            r_init[i][j]  = r_temp + rho2alpha * r3_bar[j];
        }
    }

    //	Final payoff valuation
    if (!eval_evt[nstept - 1])
    {
        szErr = "No event at last step in lgm2f_pde";
        goto FREE_RETURN;
    }

    //	Eval payoff
    szErr = payoff_func(
        date[nstept - 1],
        time[nstept - 1],
        func_parm_tab[nstept - 1],
        yc,
        &lam1,
        date,
        1,
        gamma,
        rho,
        phi1[nstept - 1],
        phi2[nstept - 1],
        phi12[nstept - 1],
        plx[nstept - 1],
        pux[nstept - 1],
        plz[nstept - 1],
        puz[nstept - 1],
        r1_dim1,
        r2,
        nprod,
        values_p1);

    if (szErr)
    {
        goto FREE_RETURN;
    }

    //	Initialize the CNPDE_TEMP_2D
    num_f_pde_init_2d_adi(pde, nstepx, nstepz, nprod);

    if (!pde)
    {
        szErr = "Memory allocation error (3) in lgm2fpde";
        goto FREE_RETURN;
    }

    //	now do the backward pde
    for (step = nstept - 2; step >= 0; step--)
    {
        dt = time[step + 1] - time[step];

        r_temp = ifr[step] + rhoalpha_plus1 * fwd1_mid[step] + rho2alpha * fwd3_mid[step];

        for (i = plx[step]; i <= pux[step]; i++)
        {
            mu_r1i = -lam1 * r1_bar[i] * dt;
            mu_r3i = fact1 * r1_bar[i];

            for (j = plz[step]; j <= puz[step]; j++)
            {
                mu_r1[i][j] = mu_r1i;
                mu_r3[i][j] = (mu_r3i - lam2 * r3_bar[j]) * dt;

                var_r1[i][j] = var1[step];
                var_r3[i][j] = var3[step];

                r[i][j] = (r_init[i][j] + r_temp) * dt;
            }
        }

        //	convolve
        num_f_pde_one_step_backward_2f_adi(
            pde,
            nstepx,
            r1_bar,
            nstepz,
            r3_bar,
            0,
            nprod - 1,
            values_p1,
            mu_r1,
            mu_r3,
            var_r1,
            var_r3,
            r,
            values,
            plx[step],
            pux[step],
            plz[step],
            puz[step]);

        /*	Eval payoff */
        if (eval_evt[step])
        {
            /*	Corresponding r1 and r2					*/
            for (i = plx[step]; i <= pux[step]; i++)
            {
                r1_dim1[i] = r1_bar[i] + fwd1[step];
                r2_temp    = rhoalpha * r1_dim1[i] + rho2alpha * fwd3[step];

                for (j = plz[step]; j <= puz[step]; j++)
                {
                    r2[i][j] = r2_temp + rho2alpha * r3_bar[j];
                }
            }

            szErr = payoff_func(
                date[step],
                time[step],
                func_parm_tab[step],
                yc,
                &lam1,
                date,
                1,
                gamma,
                rho,
                phi1[step],
                phi2[step],
                phi12[step],
                plx[step],
                pux[step],
                plz[step],
                puz[step],
                r1_dim1,
                r2,
                nprod,
                values);
            if (szErr)
            {
                goto FREE_RETURN;
            }
        }

        values_temp = values_p1;
        values_p1   = values;
        values      = values_temp;
    }

    /* copy the result					*/
    for (i = 0; i < nprod; i++)
    {
        res[i] = values_p1[index_x][index_z][i];
    }

FREE_RETURN:

    /*	allocation (2)		*/
    if (pde)
        num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);
    if (r1_dim2)
        free_dmatrix(r1_dim2, 0, nstepx - 1, 0, nstepz - 1);
    if (r2)
        free_dmatrix(r2, 0, nstepx - 1, 0, nstepz - 1);
    if (values)
        free_f3tensor(values, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
    if (values_p1)
        free_f3tensor(values_p1, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
    if (mu_r1)
        free_dmatrix(mu_r1, 0, nstepx - 1, 0, nstepz - 1);
    if (mu_r3)
        free_dmatrix(mu_r3, 0, nstepx - 1, 0, nstepz - 1);
    if (var_r1)
        free_dmatrix(var_r1, 0, nstepx - 1, 0, nstepz - 1);
    if (var_r3)
        free_dmatrix(var_r3, 0, nstepx - 1, 0, nstepz - 1);
    if (r)
        free_dmatrix(r, 0, nstepx - 1, 0, nstepz - 1);
    if (r_init)
        free_dmatrix(r_init, 0, nstepx - 1, 0, nstepz - 1);

    return szErr;

    // unreferenced
    nsteps;
}

static Err _lgm2f_adi_fixgrid(
    /*	Time data		*/
    int     nstept,
    double* time,
    double* date,
    /*	Discretisation	*/
    int nsteps,
    /*	Model data		*/
    double  lam1,
    double* sig_time,
    double* sig1,
    int     nb_sig,
    double  alpha,
    double  gamma,
    double  rho,
    /*	Product data */
    void** func_parm_tab,
    int*   eval_evt,

    /*	Market data */
    double* ifr,
    char*   yc,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date,
                       double evt_time,
                       void*  func_parm,

                       /* Market data	*/
                       void* yc,

                       /* Model data	*/
                       double* lam,
                       double* ts_time,
                       int     nb_ts,
                       double  gamma,
                       double  rho,
                       double  phi1,
                       double  phi2,
                       double  phi12,

                       /* Gride data	*/
                       int      l1,
                       int      u1,
                       int      l2,
                       int      u2,
                       double*  r1_dim2,
                       double** r2,

                       /* Vector of results to be updated */
                       int       nprod,
                       double*** prod_val),
    /*	Result */
    int     nprod,
    double* res,
    int*    plx,
    int*    pux,
    int*    plz,
    int*    puz)
{
    Err szErr = NULL;

    int    i, j, step, index_x, index_z;
    int    nstepx, nstepz;
    double meshx, meshz, dt;
    double mu_r1i, mu_r3i;
    double lam2;
    double fact1, rho2, rhoq, rhoalpha, rhoalpha_plus1, rho2alpha;
    double r_temp;

    double std1, std3, r2_temp;

    double **r1_dim2 = NULL, **r2 = NULL, ***values = NULL, ***values_p1 = NULL,
           ***values_temp = NULL, **mu_r1 = NULL, **mu_r3 = NULL, **var_r1 = NULL, **var_r3 = NULL,
           **r = NULL, **r_init = NULL;

    int lx, ux, lz, uz;

    clock_t t1, t2;

    CNPDE_TEMP_2D_ADI pdestr, *pde = &pdestr;

    //	Allocations	of time vectors
    double* phi1     = _alloca(nstept * sizeof(double));
    double* phi2     = _alloca(nstept * sizeof(double));
    double* phi12    = _alloca(nstept * sizeof(double));
    double* var1     = _alloca(nstept * sizeof(double));
    double* var3     = _alloca(nstept * sizeof(double));
    double* fwd1     = _alloca(nstept * sizeof(double));
    double* fwd3     = _alloca(nstept * sizeof(double));
    double* fwd1_mid = _alloca(nstept * sizeof(double));
    double* fwd3_mid = _alloca(nstept * sizeof(double));
    double *r1_bar, *r3_bar, *r1_dim1;

    t1 = clock();

    /* Constant	*/
    lam2 = lam1 + gamma;

    rho2           = sqrt(1.0 - rho * rho);
    rhoalpha       = rho * alpha;
    rhoalpha_plus1 = 1.0 + rhoalpha;
    rho2alpha      = alpha * rho2;
    fact1          = -rho * gamma / rho2;
    rhoq           = rho / rho2;

    /*	Allocations	of time vectors						*/
    phi1     = dvector(0, nstept - 1);
    phi2     = dvector(0, nstept - 1);
    phi12    = dvector(0, nstept - 1);
    var1     = dvector(0, nstept - 1);
    var3     = dvector(0, nstept - 1);
    fwd1     = dvector(0, nstept - 1);
    fwd3     = dvector(0, nstept - 1);
    fwd1_mid = dvector(0, nstept - 1);
    fwd3_mid = dvector(0, nstept - 1);

    if (!pde || !var1 || !var3 || !phi1 || !phi2 || !phi12 || !fwd1 || !fwd3 || !fwd1_mid ||
        !fwd3_mid)
    {
        szErr = "Memory allocation error (1) in lgm2fpde";
        goto FREE_RETURN;
    }

    /*	Calculate the expecations of r1 and r3 */
    LGM2FExpectations(
        nstept,
        time,
        lam1,
        sig_time,
        sig1,
        nb_sig,
        alpha,
        gamma,
        rho,
        fwd1,
        fwd1_mid,
        fwd3,
        fwd3_mid,
        var1,
        var3,
        phi1,
        phi2,
        phi12);

    /*	Calculation of the number of steps in each direction: since local volatility
     */
    /*	is the same, the mesh has to be the same, but the number of steps has to be adjusted	*/

    std1 = sqrt(phi1[nstept - 1]);
    std3 = 1.0 / rho2 *
           sqrt(
               phi2[nstept - 1] / alpha / alpha + rho * rho * phi1[nstept - 1] -
               2.0 * rho * rho / alpha * phi12[nstept - 1]);

    nstepx = (int)(nsteps * sqrt(std1 / std3) + 0.5);
    nstepz = (int)(nsteps * sqrt(std3 / std1) + 0.5);

    /*	nstep has to be a odd nuber			*/
    nstepx = ((int)(nstepx / 2)) * 2 + 1;
    nstepz = ((int)(nstepz / 2)) * 2 + 1;

    /*	we want at least three points in each directions
     */
    if (nstepx < 3)
    {
        nstepx = 3;
        nstepz = (int)(nsteps * nsteps / 3.0 + 0.5);
        nstepz = ((int)(nstepz / 2)) * 2 + 1;

        if (nstepz < 3)
        {
            nstepz = 3;
        }
    }

    if (nstepz < 3)
    {
        nstepz = 3;
        nstepx = (int)(nsteps * nsteps / 3.0 + 0.5);
        nstepx = ((int)(nstepx / 2)) * 2 + 1;

        if (nstepx < 3)
        {
            nstepx = 3;
        }
    }

    /*	corresponding index to the 0 value of x and z	*/
    index_x = (nstepx - 1) / 2;
    index_z = (nstepz - 1) / 2;

    meshx = 2.0 * NSTD_LGM * sqrt(phi1[nstept - 1]) / (nstepx - 1);
    meshz = 2.0 * NSTD_LGM * std3 / (nstepz - 1);

    //	Allocations of space vectors
    r1_bar  = _alloca(nstepx * sizeof(double));
    r3_bar  = _alloca(nstepz * sizeof(double));
    r1_dim1 = _alloca(nstepx * sizeof(double));

    values    = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
    values_p1 = f3tensor(0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
    mu_r1     = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    mu_r3     = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    var_r1    = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    var_r3    = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    r         = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    r_init    = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    r1_dim2   = dmatrix(0, nstepx - 1, 0, nstepz - 1);
    r2        = dmatrix(0, nstepx - 1, 0, nstepz - 1);

    if (!r1_dim2 || !r2 || !values || !values_p1 || !mu_r1 || !mu_r3 || !var_r1 || !var_r3 || !r ||
        !r_init)
    {
        szErr = "Memory allocation error (2) in lgm2fpde";
        goto FREE_RETURN;
    }

    t2 = clock();

    smessage("Phase 1 -preprocessing, time in sec: %.2f", (double)(t2 - t1) / CLOCKS_PER_SEC);
    smessage("Phase 2 -convolution, stept: %d stepx: %d stepz: %d", nstept, nstepx, nstepz);

    t1 = clock();

    /*	Then discretise space in the orthogonal system r1 / r3	*/

    r1_bar[0] = -(nstepx - 1) / 2.0 * meshx;
    r3_bar[0] = -(nstepz - 1) / 2.0 * meshz;

    for (i = 1; i < nstepx; i++)
    {
        r1_bar[i] = r1_bar[i - 1] + meshx;
    }

    for (j = 1; j < nstepz; j++)
    {
        r3_bar[j] = r3_bar[j - 1] + meshz;
    }

    /*	Corresponding r1 and r2					*/
    for (i = 0; i < nstepx; i++)
    {
        r1_dim1[i] = r1_bar[i] + fwd1[nstept - 1];
        r2_temp    = rhoalpha * r1_dim1[i] + rho2alpha * fwd3[nstept - 1];
        r_temp     = rhoalpha_plus1 * r1_bar[i];

        for (j = 0; j < nstepz; j++)
        {
            r1_dim2[i][j] = r1_dim1[i];
            r2[i][j]      = r2_temp + rho2alpha * r3_bar[j];
            r_init[i][j]  = r_temp + rho2alpha * r3_bar[j];
        }
    }

    /*	Final payoff valuation					*/
    if (!eval_evt[nstept - 1])
    {
        szErr = "No event at last step in lgm2f_pde";
        goto FREE_RETURN;
    }

    /*	Eval payoff */
    szErr = payoff_func(
        date[nstept - 1],
        time[nstept - 1],
        func_parm_tab[nstept - 1],
        yc,
        &lam1,
        date,
        1,
        gamma,
        rho,
        phi1[nstept - 1],
        phi2[nstept - 1],
        phi12[nstept - 1],
        0,
        nstepx - 1,
        0,
        nstepz - 1,
        r1_dim1,
        r2,
        nprod,
        values_p1);

    if (szErr)
    {
        goto FREE_RETURN;
    }

    /*	Initialize the CNPDE_TEMP_2D		*/

    num_f_pde_init_2d_adi(pde, nstepx, nstepz, nprod);

    if (!pde)
    {
        szErr = "Memory allocation error (3) in lgm2fpde";
        goto FREE_RETURN;
    }

    lx = 0;
    ux = nstepx - 1;
    lz = 0;
    uz = nstepz - 1;

    *plx = lx;
    *pux = ux;
    *plz = lz;
    *puz = uz;

    /*	now do the backward pde					*/
    for (step = nstept - 2; step >= 0; step--)
    {
        dt = time[step + 1] - time[step];

        r_temp = ifr[step] + rhoalpha_plus1 * fwd1_mid[step] + rho2alpha * fwd3_mid[step];

        for (i = lx; i <= ux; i++)
        {
            mu_r1i = -lam1 * r1_bar[i] * dt;
            mu_r3i = fact1 * r1_bar[i];

            for (j = lz; j <= uz; j++)
            {
                mu_r1[i][j] = mu_r1i;
                mu_r3[i][j] = (mu_r3i - lam2 * r3_bar[j]) * dt;

                var_r1[i][j] = var1[step];
                var_r3[i][j] = var3[step];

                r[i][j] = (r_init[i][j] + r_temp) * dt;
            }
        }

        /*	convolve							*/

        num_f_pde_one_step_backward_2f_adi(
            pde,
            nstepx,
            r1_bar,
            nstepz,
            r3_bar,
            0,
            nprod - 1,
            values_p1,
            mu_r1,
            mu_r3,
            var_r1,
            var_r3,
            r,
            values,
            lx,
            ux,
            lz,
            uz);

        /*	Eval payoff */
        if (eval_evt[step])
        {
            /*	Corresponding r1 and r2					*/
            for (i = lx; i <= ux; i++)
            {
                r1_dim1[i] = r1_bar[i] + fwd1[step];
                r2_temp    = rhoalpha * r1_dim1[i] + rho2alpha * fwd3[step];

                for (j = lz; j <= uz; j++)
                {
                    r2[i][j] = r2_temp + rho2alpha * r3_bar[j];
                }
            }

            szErr = payoff_func(
                date[step],
                time[step],
                func_parm_tab[step],
                yc,
                &lam1,
                date,
                1,
                gamma,
                rho,
                phi1[step],
                phi2[step],
                phi12[step],
                lx,
                ux,
                lz,
                uz,
                r1_dim1,
                r2,
                nprod,
                values);
            if (szErr)
            {
                goto FREE_RETURN;
            }
        }

        values_temp = values_p1;
        values_p1   = values;
        values      = values_temp;

        /* new indexes: we cut the PDE after NSTD_LGM number of standard deviations */
        /* but we need at least three points to do the pde
         */

        /// do not reduce the number of points when computing the exercise probablities
#if 0
		tem = (int) (NSTD_LGM * sqrt(phi1[step]) / meshx + 0.5);
		ux = min(nstepx-1, index_x + tem);
		lx = max(0, index_x - tem);
		
		if (ux - lx < 2)
		{
			tem += 1;
			ux = min(nstepx-1, index_x + tem);
			lx = max(0, index_x - tem);
		}

		tem = (int) (NSTD_LGM / rho2 * sqrt(phi2[step] / alpha / alpha
									  + rho * rho * phi1[step]
									  - 2.0 * rho * rho / alpha * phi12[step]) / meshx + 0.5);			
		uz = min(nstepz-1, index_z + tem);
		lz = max(0, index_z - tem);

		if (uz - lz < 2)
		{
			tem += 1;
			uz = min(nstepz-1, index_z + tem);
			lz = max(0, index_z - tem);
		}
#endif
    }

    /* copy the result					*/
    for (i = 0; i < nprod; i++)
    {
        res[i] = values_p1[index_x][index_z][i];
    }

    t2 = clock();

    smessage("Phase 2 -convolution, time in sec: %.2f", (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

    /*	allocation (2)		*/
    if (pde)
        num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

    if (r1_dim2)
        free_dmatrix(r1_dim2, 0, nstepx - 1, 0, nstepz - 1);
    if (r2)
        free_dmatrix(r2, 0, nstepx - 1, 0, nstepz - 1);
    if (values)
        free_f3tensor(values, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
    if (values_p1)
        free_f3tensor(values_p1, 0, nstepx - 1, 0, nstepz - 1, 0, nprod - 1);
    if (mu_r1)
        free_dmatrix(mu_r1, 0, nstepx - 1, 0, nstepz - 1);
    if (mu_r3)
        free_dmatrix(mu_r3, 0, nstepx - 1, 0, nstepz - 1);
    if (var_r1)
        free_dmatrix(var_r1, 0, nstepx - 1, 0, nstepz - 1);
    if (var_r3)
        free_dmatrix(var_r3, 0, nstepx - 1, 0, nstepz - 1);
    if (r)
        free_dmatrix(r, 0, nstepx - 1, 0, nstepz - 1);
    if (r_init)
        free_dmatrix(r_init, 0, nstepx - 1, 0, nstepz - 1);

    return szErr;
}

static void _init_genmidat_arg(
    genmidat_arg*     pAmortArg,
    ptr_genmidat_arg* ppAmortArg,
    int*              pdUB_D1,
    int*              pdLB_D1,
    int*              pdUB_D2,
    int*              pdLB_D2,
    int*              pnum_fix_cpn_done,
    int*              pnum_fund_cpn_done,
    int*              pnEx,
    double*           pdInit_Prev,
    double*           pdIV,
    const double*     pdStrike,
    double***         pppdExIndicator,
    double***         pppdCondiDF,
    double*           pdAvgSwapRate,
    int*              pnExToTime,
    AM_ADI_ARG        adi_arg,
    long              lToday,
    int               nCompExProb)
{
    static int _accumulate_(const int* piBegin, const int* piEnd);
    const int  nNumEx = _accumulate_(adi_arg->is_event, (adi_arg->is_event) + adi_arg->nstp);
    int        nI, nJ;

    //// initialize the parameters
    for (nI = 0, nJ = -1; nI < adi_arg->nstp; ++nI)
    {
        pAmortArg[nI].m_pBase            = (AM_PAY_ARG)(adi_arg->void_prm[nI]);
        pAmortArg[nI].pnum_fix_cpn_done  = pnum_fix_cpn_done;
        pAmortArg[nI].pnum_fund_cpn_done = pnum_fund_cpn_done;
        pAmortArg[nI].pdInit_Prev        = pdInit_Prev;
        pAmortArg[nI].pdIV               = pdIV;
        pAmortArg[nI].pdStrike           = pdStrike;
        pAmortArg[nI].lToday             = lToday;
        pAmortArg[nI].nNumEx             = nNumEx;
        pAmortArg[nI].pnEx               = pnEx;
        pAmortArg[nI].nCompExProb        = nCompExProb;

        if (adi_arg->is_event[nI])
        {
            ++nJ;
            pAmortArg[nI].pppdExIndicator = &pppdExIndicator[nJ];
            pAmortArg[nI].pppdCondiDF     = &pppdCondiDF[nJ];
            pAmortArg[nI].pdAvgSwapRate   = pdAvgSwapRate ? &pdAvgSwapRate[nJ] : 0;

            pAmortArg[nI].pdUB_D1 = &(pdUB_D1[nJ]);
            pAmortArg[nI].pdLB_D1 = &(pdLB_D1[nJ]);
            pAmortArg[nI].pdUB_D2 = &(pdUB_D2[nJ]);
            pAmortArg[nI].pdLB_D2 = &(pdLB_D2[nJ]);
            pnExToTime[nJ]        = nI;
        }
        else
        {
            pAmortArg[nI].pppdExIndicator = 0;
            pAmortArg[nI].pppdCondiDF     = 0;
            pAmortArg[nI].pdAvgSwapRate   = 0;

            pAmortArg[nI].pdUB_D1 = 0;
            pAmortArg[nI].pdLB_D1 = 0;
            pAmortArg[nI].pdUB_D2 = 0;
            pAmortArg[nI].pdLB_D2 = 0;
        }

        ppAmortArg[nI] = &pAmortArg[nI];
    }
}

static void _copy_(int* pnBegin, int* pnEnd, int nValue)
{
    for (; pnBegin < pnEnd; ++pnBegin)
        *pnBegin = nValue;
}

static Err _compute_exercise_prob_(
    ptr_genmidat_arg* ppAmortArg,
    AM_ADI_ARG        adi_arg,
    int               lx,
    int               ux,
    int               lz,
    int               uz,
    const int*        pnExToTime,
    double*           pdExProb  // returned exercise probablilities
)
{
    static int  _accumulate_(const int* piBegin, const int* piEnd);
    static void _copy_(int* pnBegin, int* pnEnd, int nValue);

    const int nNumStep = adi_arg->nstp;
    const int nNumEx   = _accumulate_(adi_arg->is_event, (adi_arg->is_event) + nNumStep);

    // allocate mem
    int* plx  = _alloca(nNumStep * sizeof(int));
    int* pux  = _alloca(nNumStep * sizeof(int));
    int* plz  = _alloca(nNumStep * sizeof(int));
    int* puz  = _alloca(nNumStep * sizeof(int));
    int* pnEx = ppAmortArg[0]->pnEx;

    /// first compute conditional discount factors
    double dDummy;
    Err    szErr = 0;
    int    nJ;  //, nEx;

    /// set mem
    _copy_(plx, plx + nNumStep, lx);
    _copy_(pux, pux + nNumStep, ux);
    _copy_(plz, plz + nNumStep, lz);
    _copy_(puz, puz + nNumStep, uz);

    //// _lgm2f_adi_shrinkgrid is idential to lgm2f_adi ,except
    //// that the bounds are determined exogenously;  this variation is needed for
    //// computing exercise probablities; bounds of the 2D grid are stored in plx,pux,plz,puz
    if (szErr = _lgm2f_adi_shrinkgrid(
            plx,
            pux,
            plz,
            puz,
            adi_arg->nstp,
            adi_arg->time,
            adi_arg->date,
            adi_arg->nstpx,
            adi_arg->lam,
            adi_arg->sig_time,
            adi_arg->sig1,
            adi_arg->nb_sig,
            adi_arg->alpha,
            adi_arg->gamma,
            adi_arg->rho,
            (void**)ppAmortArg,  // adi_arg->void_prm,//
            adi_arg->is_event,
            adi_arg->ifr,
            adi_arg->yc,
            am_CondiDF_4_lgm2f_adi_fixgrid,
            1,
            &dDummy))
    {
        return szErr;
    }

    /// probablility of not exercising the deal at all
    nJ    = nNumEx - 1;
    *pnEx = nJ;
    if (szErr = _lgm2f_adi_shrinkgrid(
            plx,
            pux,
            plz,
            puz,
            pnExToTime[nJ] + 1,  // adi_arg->nstp,
            adi_arg->time,
            adi_arg->date,
            adi_arg->nstpx,
            adi_arg->lam,
            adi_arg->sig_time,
            adi_arg->sig1,
            adi_arg->nb_sig,
            adi_arg->alpha,
            adi_arg->gamma,
            adi_arg->rho,
            (void**)ppAmortArg,  // adi_arg->void_prm,//
            adi_arg->is_event,
            adi_arg->ifr,
            adi_arg->yc,
            am_notexerprob_4_lgm2f_adi_fixgrid,
            1,
            &pdExProb[nNumEx]))
    {
        return szErr;
    }

    pdExProb[nNumEx] /= dDummy;

    for (nJ = nNumEx - 1; nJ >= 0; --nJ)
    {
        *pnEx = nJ;

        if (szErr = _lgm2f_adi_shrinkgrid(
                plx,
                pux,
                plz,
                puz,
                pnExToTime[nJ] + 1,  // adi_arg->nstp,
                adi_arg->time,
                adi_arg->date,
                adi_arg->nstpx,
                adi_arg->lam,
                adi_arg->sig_time,
                adi_arg->sig1,
                adi_arg->nb_sig,
                adi_arg->alpha,
                adi_arg->gamma,
                adi_arg->rho,
                (void**)ppAmortArg,  // adi_arg->void_prm,//
                adi_arg->is_event,
                adi_arg->ifr,
                adi_arg->yc,
                am_exerprob_4_lgm2f_adi_fixgrid,
                1,
                &pdExProb[nJ]))
        {
            return szErr;
        }

        pdExProb[nJ] /= dDummy;
    }

    // resample to 1.
    {
        double dSumProb = 0.;

        for (nJ = 0; nJ <= nNumEx; ++nJ)
            dSumProb += pdExProb[nJ];

        for (nJ = 0; nJ <= nNumEx; ++nJ)
            pdExProb[nJ] /= dSumProb;
    }

    return szErr;
}

//// prices generic midat only
static Err _price_genmidat_(
    AM_ADI_ARG        adi_arg,  // structure
    ptr_genmidat_arg* ppAmortArg,
    double*           pdOptionPV)
{
    double _pdOptionPV[2] = {0};
    Err    szErr          = 0;

    /// delegate ....
    szErr = lgm2f_adi(
        adi_arg->nstp,
        adi_arg->time,
        adi_arg->date,
        adi_arg->nstpx,
        adi_arg->lam,
        adi_arg->sig_time,
        adi_arg->sig1,
        adi_arg->nb_sig,
        adi_arg->alpha,
        adi_arg->gamma,
        adi_arg->rho,
        (void**)ppAmortArg,
        adi_arg->is_event,
        adi_arg->ifr,
        adi_arg->yc,
        _genmidat_nocall_payoff_,
        1,
        _pdOptionPV);

    /// NB: because *pdOptionPV = option value - swap(first exercise)
    /// add Swap value to the option
    *pdOptionPV = _pdOptionPV[0];
    (*pdOptionPV) += *(ppAmortArg[0]->pdIV);

    return szErr;
}

//// prices generic midat only
static Err _price_genmidat_ts_(
    AM_ADI_ARG        adi_arg,  // structure
    ptr_genmidat_arg* ppAmortArg,
    const double*     pdLambdaValue,
    const double*     pdLambdaTime,  //// NB: calendar time, not year fraction !!!!!!!
    int               nLambdaSize,
    double*           pdOptionPV)
{
    double    _pdOptionPV[2] = {0};
    const int ndisc_method   = 1;
    Err       szErr          = 0;

    /// delegate ....
    szErr = lgm2fTau_adi(
        adi_arg->nstp,
        adi_arg->time,
        adi_arg->date,

        adi_arg->nstpx,
        ndisc_method,

        adi_arg->sig1,
        adi_arg->sig_time,
        adi_arg->nb_sig,

        (double*)pdLambdaValue,
        (double*)pdLambdaTime,  //// NB: calendar time, not year fraction !!!!!!!
        nLambdaSize,

        adi_arg->alpha,
        adi_arg->gamma,
        adi_arg->rho,
        (void**)ppAmortArg,  // adi_arg->void_prm,
        adi_arg->is_event,
        adi_arg->ifr,
        adi_arg->yc,
        _genmidat_nocall_payoff_,  // payoff_lgm2fTau_pde,//am_payoff_4_lgm2f_adi,
        1,
        _pdOptionPV);

    /// NB: because *pdOptionPV = option value - swap(first exercise)
    /// add Swap value to the option
    *pdOptionPV = _pdOptionPV[0];
    (*pdOptionPV) += *(ppAmortArg[0]->pdIV);

    return szErr;
}

static void _free_genmidat_arg(
    ptr_genmidat_arg* ppAmortArg, const int* pnExToTime, int nNumEx, int lx, int ux, int lz, int uz)
{
    int      nJ, nExToTime;
    double** pdmatrix = 0;
    if (ppAmortArg[0]->nCompExProb)
    {
        for (nJ = 0; nJ < nNumEx; ++nJ)
        {
            nExToTime = pnExToTime[nJ];
            pdmatrix  = *(ppAmortArg[nExToTime]->pppdExIndicator);
            if (pdmatrix)
                free_dmatrix(pdmatrix, lx, ux, lz, uz);
            pdmatrix = *(ppAmortArg[nExToTime]->pppdCondiDF);
            if (pdmatrix)
                free_dmatrix(pdmatrix, lx, ux, lz, uz);
        }
    }
}

//// prices both generic midat and exercise probabilites
static Err _price_genmidat_exprob_(
    AM_ADI_ARG        adi_arg,  // structure
    ptr_genmidat_arg* ppAmortArg,
    const int*        pnExToTime,
    double*           pdOptionPV,
    double*           pdExProb)
{
    const int nNumEx         = _accumulate_(adi_arg->is_event, (adi_arg->is_event) + adi_arg->nstp);
    double    _pdOptionPV[2] = {0};
    int       lx, ux, lz, uz;

    Err szErr;

    //// delegate to price the option
    //// _lgm2f_adi_fixgrid is idential to lgm2f_adi ,except
    //// that once the grid is set up - by the solver - at the termination, it remains fixed  as the
    ///solver / advances backwards in time; this variation is needed for / computing exercise
    ///probablities; bounds of the 2D grid are stored in plx,pux,plz,puz
    if (szErr = _lgm2f_adi_fixgrid(
            adi_arg->nstp,
            adi_arg->time,
            adi_arg->date,
            adi_arg->nstpx,
            adi_arg->lam,
            adi_arg->sig_time,
            adi_arg->sig1,
            adi_arg->nb_sig,
            adi_arg->alpha,
            adi_arg->gamma,
            adi_arg->rho,
            (void**)ppAmortArg,
            adi_arg->is_event,
            adi_arg->ifr,
            adi_arg->yc,
            _genmidat_nocall_payoff_,
            1,
            _pdOptionPV,
            &lx,
            &ux,
            &lz,
            &uz))
    {
        _free_genmidat_arg(ppAmortArg, pnExToTime, nNumEx, lx, ux, lz, uz);
        return szErr;
    }

    /// NB: because *pdOptionPV = option value - swap(first exercise)
    /// add Swap value to the option
    *pdOptionPV = _pdOptionPV[0];
    (*pdOptionPV) += *(ppAmortArg[0]->pdIV);

    /// delegate to compute exercise probabilities
    szErr = _compute_exercise_prob_(ppAmortArg, adi_arg, lx, ux, lz, uz, pnExToTime, pdExProb);

    _free_genmidat_arg(ppAmortArg, pnExToTime, nNumEx, lx, ux, lz, uz);

    return szErr;
}

static Err _price_genmidat(
    AM_STR     am,       // structure
    AM_UND     und,      // structure
    AM_ADI_ARG adi_arg,  // structure
    long       lToday,
    double     dStrike,
    double*    pdOptionPV,
    double**   ppdExProb,
    double**   ppdAvgStrike,
    int*       pnEx)
{
    // total number of exercise ;
    // NB: is_event contains {0,1 sequence}, indicating occurrence of an event/exercise
    const int nNumEx = _accumulate_(adi_arg->is_event, (adi_arg->is_event) + adi_arg->nstp);

    const int nCompExProb = ppdExProb && ppdAvgStrike ? 1 : 0;

    /// optimized midat payoff function utilizes the midat "put-call" parity
    /// num_fund_cpn_done and num_fix_cpn_done keep track of the number
    /// of coupons already accounted for by later exercises
    int    num_fund_cpn_done = 0, num_fix_cpn_done = 0;
    double dInit_Prev = 0.;
    double dIV_FirstEx;  // dIV_FirstEx stores the PV of the swap with no cash flow until the first
                         // exercise date
    Err szErr;

    /// user-provided payoff function arguments
    genmidat_arg* pAmortArg = memset(
        _alloca(adi_arg->nstp * sizeof(genmidat_arg)), 0, adi_arg->nstp * sizeof(genmidat_arg));
    ptr_genmidat_arg* ppAmortArg = memset(
        _alloca(adi_arg->nstp * sizeof(ptr_genmidat_arg)),
        0,
        adi_arg->nstp * sizeof(ptr_genmidat_arg));

    // global exercise indicator
    double*** pppdExIndicator =
        memset(_alloca(nNumEx * sizeof(double**)), 0, nNumEx * sizeof(double**));

    // global conditional discount factor
    // i.e. discount factors observed at time tnEx and paid at tnEx+1 (>tnEx)
    double*** pppdCondiDF =
        memset(_alloca(nNumEx * sizeof(double**)), 0, nNumEx * sizeof(double**));

    // global average strike
    // double* pdAvgStrike = memset(_alloca(nNumEx * sizeof(double)),0,nNumEx*sizeof(double));

    /// bounds of the 2D grid at exercises
    int* pdUB_D1 = _alloca(nNumEx * sizeof(int));
    int* pdLB_D1 = _alloca(nNumEx * sizeof(int));
    int* pdUB_D2 = _alloca(nNumEx * sizeof(int));
    int* pdLB_D2 = _alloca(nNumEx * sizeof(int));

    int nEx;
    int one_lam = 1;

    /// time step correponding to each exercise
    int* pnExToTime =
        memset(_alloca(nNumEx * sizeof(int)), 0, nNumEx * sizeof(int));  /// initialize to 0

    *pnEx = nNumEx;

    // NB: exercise probablities; +1 to include
    // probability of not exercising at all
    if (nCompExProb)
    {
        *ppdExProb    = calloc(nNumEx + 1, sizeof(double));
        *ppdAvgStrike = calloc(nNumEx, sizeof(double));
    }

    //// initialize payoff function parameters
    _init_genmidat_arg(
        pAmortArg,
        ppAmortArg,
        pdUB_D1,
        pdLB_D1,
        pdUB_D2,
        pdLB_D2,
        &num_fix_cpn_done,
        &num_fund_cpn_done,
        &nEx,
        &dInit_Prev,
        &dIV_FirstEx,
        &dStrike,
        pppdExIndicator,
        pppdCondiDF,
        ppdAvgStrike ? *ppdAvgStrike : 0,
        pnExToTime,
        adi_arg,
        lToday,
        nCompExProb);

    /// compute option value only
    if (!nCompExProb)
        return _price_genmidat_(adi_arg, ppAmortArg, pdOptionPV);
    // else compute both option value and exercise probabilities
    else
        szErr = _price_genmidat_exprob_(adi_arg, ppAmortArg, pnExToTime, pdOptionPV, *ppdExProb);

    und;
    am;

    return szErr;
}

static Err _price_genmidat_ts(
    AM_STR        am,       // structure
    AM_UND        und,      // structure
    AM_ADI_ARG    adi_arg,  // structure
    long          lToday,
    double        dStrike,
    const double* pdLambdaValue,
    const double* pdLambdaTime,  //// NB: calendar time, not year fraction !!!!!!!
    int           nLambdaSize,
    double*       pdOptionPV,
    double**      ppdExProb,
    double**      ppdAvgSwapRate,
    int*          pnEx)
{
    // total number of exercise ;
    // NB: is_event contains {0,1 sequence}, indicating occurrence of an event/exercise
    const int nNumEx = _accumulate_(adi_arg->is_event, (adi_arg->is_event) + adi_arg->nstp);

    const int nCompExProb = ppdExProb && ppdAvgSwapRate ? 1 : 0;

    /// optimized midat payoff function utilizes the midat "put-call" parity
    /// num_fund_cpn_done and num_fix_cpn_done keep track of the number
    /// of coupons already accounted for by later exercises
    int    num_fund_cpn_done = 0, num_fix_cpn_done = 0;
    double dInit_Prev = 0.;
    double dIV_FirstEx;  // dIV_FirstEx stores the PV of the swap with no cash flow until the first
                         // exercise date
    Err szErr;

    /// user-provided payoff function arguments
    genmidat_arg* pAmortArg = memset(
        _alloca(adi_arg->nstp * sizeof(genmidat_arg)), 0, adi_arg->nstp * sizeof(genmidat_arg));
    ptr_genmidat_arg* ppAmortArg = memset(
        _alloca(adi_arg->nstp * sizeof(ptr_genmidat_arg)),
        0,
        adi_arg->nstp * sizeof(ptr_genmidat_arg));

    // global exercise indicator
    double*** pppdExIndicator =
        memset(_alloca(nNumEx * sizeof(double**)), 0, nNumEx * sizeof(double**));

    // global conditional discount factor
    // i.e. discount factors observed at time tnEx and paid at tnEx+1 (>tnEx)
    double*** pppdCondiDF =
        memset(_alloca(nNumEx * sizeof(double**)), 0, nNumEx * sizeof(double**));

    /// bounds of the 2D grid at exercises
    int* pdUB_D1 = _alloca(nNumEx * sizeof(int));
    int* pdLB_D1 = _alloca(nNumEx * sizeof(int));
    int* pdUB_D2 = _alloca(nNumEx * sizeof(int));
    int* pdLB_D2 = _alloca(nNumEx * sizeof(int));

    int nEx;
    int one_lam = 1;

    /// time step correponding to each exercise
    int* pnExToTime =
        memset(_alloca(nNumEx * sizeof(int)), 0, nNumEx * sizeof(int));  /// initialize to 0

    // NB: exercise probablities; +1 to include
    // probability of not exercising at all
    if (nCompExProb)
    {
        *ppdExProb      = calloc(nNumEx + 1, sizeof(double));
        *ppdAvgSwapRate = calloc(nNumEx + 1, sizeof(double));
    }

    *pnEx = nNumEx;

    //// initialize payoff function parameters
    _init_genmidat_arg(
        pAmortArg,
        ppAmortArg,
        pdUB_D1,
        pdLB_D1,
        pdUB_D2,
        pdLB_D2,
        &num_fix_cpn_done,
        &num_fund_cpn_done,
        &nEx,
        &dInit_Prev,
        &dIV_FirstEx,
        &dStrike,
        pppdExIndicator,
        pppdCondiDF,
        ppdAvgSwapRate ? *ppdAvgSwapRate : 0,
        pnExToTime,
        adi_arg,
        lToday,
        nCompExProb);

    /// compute option value only
    if (!nCompExProb)
        return _price_genmidat_ts_(
            adi_arg, ppAmortArg, pdLambdaValue, pdLambdaTime, nLambdaSize, pdOptionPV);
    // else compute both option value and exercise probabilities
    else
        szErr =
            "exercise probability not available with Tau TS yet !";  // _price_genmidat_exprob_(adi_arg,ppAmortArg,pnExToTime,pdOptionPV,*ppdExProb);

    und;
    am;

    return szErr;
}

static int _accumulate_(const int* piBegin, const int* piEnd)
{
    int iSum = 0;
    for (; piBegin != piEnd; ++piBegin)
        iSum += (*piBegin);
    return iSum;
}

Err convert_funding_to_domestic_amort(
    //	Inputs
    long today,        //	Today
    long not_ex_date,  //	Date at which the
                       //		initial notional
                       //		exchange takes place
                       //		(or has taken place)
    int eod_fix_flag,  //	0: I, 1: E
    int eod_pay_flag,  //	0: I, 1: E

    double fx_fund_dom,            //	Fx fund/dom, 2bd fwd
    long   fx_fund_dom_spot_date,  //	Spot date for Fx

    char*   dom_yc,      //	Domestic discount curve
    int     nb_dom_not,  //	Number of Domestic notional
    double* dom_not,     //	Domestic notional

    int    fund_ncpn,   //	Number of coupon
    long*  fund_fix,    //	Fixing dates
    long*  fund_start,  //	Start dates
    long*  fund_pay,    //	Pay dates
    char** fund_basis,  //	Basis
    //	The following are modified
    char* fund_yc,         //	Funding discount curve,
                           //	changed to domestic
    double* fund_not,      //	Funding notional
                           //		in funding ccy,
                           //	converted to domestic
    double* fund_spr,      //	Spread in funding ccy,
                           //	put to 0
    double* fund_mrg,      //	Margin in funding ccy,
                           //	converted to margin over
                           //	cash libor in domestic currency
    double* fund_fix_cpn,  //	Fixing: contains spread
                           //		but not margin,
                           //	converted to equivalent
                           //	domestic cash-flows
    //	The following are returned
    long*   fund_start_date,  //	Start date of the funding
    double* eq_final_ex,      //	Domestic cash-flow equivalent
                              //		to final exchange
                              //	(to be delivered
                              //	at funding start date)
    double* eq_init_ex);      //	Domestic cash-flow equivalent
                              //		to initial exchange
                              //	(to be delivered
                              //	at initial exchange date)

Err am_caller_AmortMidatAutocal(
    /*	Today's date */
    long today,

    /*	The underlying */
    int use_calib, /*	0: use lgm2fund, 1: calibrate */

    /*		if calib */
    char*  yc,         /*	yc */
    char*  vc,         /*	vc */
    char*  ref,        /*	ref rate (only if calib) */
    char*  swap_freq,  /*	swap freq (only if calib) */
    char*  swap_basis, /*	swap basis (only if calib) */
    double lambda,     /*	lambda if unique */
    double alpha,      /*	alpha */
    double gamma,      /*	gamma */
    double rho,        /*	rho */

    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */
    long start_date,  /*	Date at which initial notional exchange occurs */
    long theoEndDate, /* End Date of the structure */

    /*		funding */
    char*   fund_ref,
    int     fund_ccy,    /*	0: domestic, 1: other */
    double* fund_not,    /*	If different from domestic or foreign (fund_ccy = 1) */
    char*   fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 1) */
    double  fx_fund_dom, /*	If different from domestic or foreign (fund_ccy = 1) 2 bd fwd */
    long    fx_fund_dom_spot_date,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    double* fund_fix_cpn, /*	Past coupon fixing if relevant,
                                          includes spr, but not mrg, cvg and notional */
    /*		fix */
    double* fix_not,
    int     fix_ncpn,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee, /*	Exercise Fee */

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    double mintime,
    double mininterval,
    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    /*	Function to get correlation stored from a vol cube */
    Err (*get_correl)(
        char* correl_cube_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,

    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                       to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_pay_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Exercised flag */
    int    exercised,   /*	Flag */
    long   ex_date_ex,  /*	Date when exercised */
    long   ex_date_set, /*	Corresponding settlement date */
    double ex_fee,      /*	Corresponding fee */

    int nReduCalibPoints,
    int nOutPutExProb,

    /*	Results */
    double*  fund_val, /*	Value of the funding leg */
    double*  fix_val,  /*	Value of the Power Dual leg */
    double*  call_val, /*	Value of the callable feature */
    double** ppdExProb,
    double** ppdAvgSwapRate,
    int*     pnEx,
    int      export_ts, /*	1: Export TS, 0: don't */
    AM_UND   und_exp)
{
    am_str*     am      = NULL;
    am_und*     und     = NULL;
    am_adi_arg* adi_arg = NULL;

    AM_FUND_LEG fund_leg;
    AM_FUND_CPN fund_cpn;
    AM_FIX_LEG  fix_leg;
    AM_FIX_CPN  fix_cpn;

    int          call_feat;
    double       fund_leg_pv, fix_leg_pv;
    SrtBasisCode bas;

    double df, coupon;

    double call;  //,iv;
    int    i, j;
    int    free_struct = 0;

    double temp;

    int    for_fund, nEx;
    long   fund_start_date, fin_not_date;
    double eq_final_ex, eq_init_ex;

    Err szErr = NULL;

    /*	If exercised */
    if (exercised)
    {
        i = 0;
        while (i < fix_ncpn && fix_start[i] < ex_date_ex)
        {
            i++;
        }
        fix_ncpn = i;

        /*	Structure is called before start: return 0 */
        if (fix_ncpn == 0)
        {
            *fund_val = *fix_val = *call_val = 0.0;
            return NULL;
        }

        i = 0;
        while (i < fund_ncpn && fund_start[i] < ex_date_ex)
        {
            i++;
        }
        fund_ncpn = i;

        ncall     = 0;
        exercised = 0;

        szErr = am_caller_AmortMidatAutocal(
            today,
            use_calib,
            yc,
            vc,
            ref,
            swap_freq,
            swap_basis,
            lambda,
            alpha,
            gamma,
            rho,
            get_cash_vol,
            lgm2dund,
            start_date,
            theoEndDate,
            fund_ref,
            fund_ccy,
            fund_not,
            fund_ccy_yc,
            fx_fund_dom,
            fx_fund_dom_spot_date,

            fund_ncpn,
            fund_fix,
            fund_start,
            fund_end,
            fund_pay,
            fund_basis,
            fund_spr,
            fund_mrg,
            fund_fix_cpn,

            fix_not,
            fix_ncpn,
            fix_start,
            fix_end,
            fix_pay,
            fix_basis,
            fix_rate,
            fix_fee,

            ncall,
            pay_rec,
            ex_date,
            set_date,
            fee,

            req_stp,
            req_stpx,

            mintime,
            mininterval,
            notperiod,
            one2F,
            use_jump,
            max_var_jump,
            strike_type,
            european_model,

            get_correl,
            CorrelName,

            max_std_short,
            fix_lambda,
            one_f_equi,
            skip_last,

            eod_fix_flag,
            eod_pay_flag,
            eod_ex_flag,

            exercised,
            ex_date_ex,
            ex_date_set,
            ex_fee,

            nReduCalibPoints,
            nOutPutExProb,

            fund_val,
            fix_val,
            call_val,
            ppdExProb,
            ppdAvgSwapRate,
            &nEx,
            export_ts,
            und_exp);

        if (szErr)
        {
            goto FREE_RETURN;
        }

        if (ex_date_set >= today + eod_pay_flag)
        {
            *call_val = -ex_fee * swp_f_df(today, ex_date_set, yc);
        }

        goto FREE_RETURN;
    }

    if (fund_ccy == 1)
    {
        fund_ccy = 0;
        for_fund = 1;
        szErr    = convert_funding_to_domestic_amort(
            today,
            start_date,
            eod_fix_flag,
            eod_pay_flag,
            fx_fund_dom,
            fx_fund_dom_spot_date,
            yc,
            fix_ncpn,
            fix_not,
            fund_ncpn,
            fund_fix,
            fund_start,
            fund_pay,
            fund_basis,
            fund_ccy_yc,
            fund_not,
            fund_spr,
            fund_mrg,
            fund_fix_cpn,
            &fund_start_date,
            &eq_final_ex,
            &eq_init_ex);
        if (szErr)
        {
            return szErr;
        }
    }
    else
    {
        for_fund = 0;
    }

    am      = calloc(1, sizeof(am_str));
    und     = calloc(1, sizeof(am_und));
    adi_arg = calloc(1, sizeof(am_adi_arg));

    if (!am || !und || !adi_arg)
    {
        szErr = "memory allocation failure in am_caller";
        goto FREE_RETURN;
    }

    /*	Initialise structures */
    free_struct = 0;
    szErr       = am_fill_check_all_struct_AmortMidatAutocal(
        today,
        theoEndDate,

        use_calib,
        yc,
        vc,
        ref,
        swap_freq,
        swap_basis,
        lambda,
        alpha,
        gamma,
        rho,

        lgm2dund,

        fund_ref,
        fund_not,
        fund_ncpn,
        fund_fix,
        fund_start,
        fund_end,
        fund_pay,
        fund_basis,
        fund_spr,
        fund_mrg,

        fix_not,
        fix_ncpn,
        fix_start,
        fix_end,
        fix_pay,
        fix_basis,
        fix_rate,
        fix_fee,

        ncall,
        pay_rec,
        ex_date,
        set_date,
        fee,

        req_stp,
        req_stpx,

        get_cash_vol,

        mintime,
        mininterval,
        notperiod,
        one2F,
        use_jump,
        max_var_jump,
        strike_type,
        european_model,

        get_correl,
        CorrelName,
        /*
                        GetVolForBadr,
                        cVolType,
        */
        max_std_short,
        fix_lambda,
        one_f_equi,
        skip_last,

        eod_fix_flag,
        eod_ex_flag,

        am,
        und,

        &call_feat,
        adi_arg);

    if (szErr)
    {
        goto FREE_RETURN;
    }
    free_struct = 1;

    /*	1)	Value funding leg */
    fund_leg    = am->fund_leg;
    fund_leg_pv = 0.0;

    /*	Cash libor */
    if (fund_leg->num_cpn > 0)
    {
        if (for_fund)
        {
            for (j = 0; j < fund_leg->num_cpn - 1; ++j)
            {
                if (fund_leg->cpn[j].pay_date >= today + eod_pay_flag)
                {
                    temp = swp_f_df(today, fund_leg->cpn[j].pay_date, yc) *
                           (fund_not[j] - fund_not[j + 1]);
                    fund_leg_pv += temp;
                }
            }

            if (fund_leg->cpn[fund_leg->num_cpn - 1].pay_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fund_leg->cpn[fund_leg->num_cpn - 1].pay_date, yc) *
                       fund_not[fund_leg->num_cpn - 1];
                fund_leg_pv += temp;
            }
        }
        else
        {
            fund_leg_pv += swp_f_df(today, fund_leg->cpn[0].start_date, yc) * fund_leg->notional[0];
        }
    }

    /*	Coupons: spread + margin */
    for (i = 0; i < fund_leg->num_cpn; i++)
    {
        fund_cpn = fund_leg->cpn + i;
        temp     = swp_f_df(today, fund_cpn->pay_date, yc) * fund_cpn->cpn;
        fund_leg_pv += temp;
    }

    /*	PV of coupons fixed in the past and not yet paid */
    i = 0;
    while (i < fund_ncpn && fund_fix[i] < today + eod_fix_flag)
    {
        if (fund_pay[i] >= today + eod_pay_flag)
        {
            szErr = interp_basis(fund_basis[i], &bas);
            if (szErr)
            {
                goto FREE_RETURN;
            }

            fund_leg_pv += (fund_fix_cpn[i] + fund_mrg[i]) *
                           coverage(fund_start[i], fund_pay[i], bas) * fund_not[0] *
                           swp_f_df(today, fund_pay[i], yc);
        }

        i++;
    }

    /*	2)	Value fix leg */

    fix_leg    = am->fix_leg;
    fix_leg_pv = 0.0;
    // fix_fee = 0.;

    /*	Coupons */
    for (i = 0; i < fix_leg->num_cpn; i++)
    {
        fix_cpn = fix_leg->cpn + i;

        /*	Discount */
        df = swp_f_df(today, fix_cpn->pay_date, yc);

        coupon = fix_cpn->cpn;

        /*	Coupon pv */
        temp = df * coupon;
        fix_leg_pv += temp;
    }

    /*	PV of coupons of past periods and not yet paid */
    i = 0;
    //	while (i < fix_ncpn && fix_fix[i] < today + eod_fix_flag)
    while (i < fix_ncpn && fix_start[i] < today + eod_fix_flag)
    {
        if (fix_pay[i] >= today + eod_pay_flag)
        {
            szErr = interp_basis(fix_basis[i], &bas);
            if (szErr)
            {
                goto FREE_RETURN;
            }

            //			fix_leg_pv += fix_fix_cpn[i]
            fix_leg_pv += fix_rate[i] * coverage(fix_start[i], fix_pay[i], bas) * fix_not[i] *
                          swp_f_df(today, fix_pay[i], yc);
        }

        i++;
    }

    /*	Initial and final exchange */
    if (for_fund)
    {
        /*	Final */
        if (fix_leg->num_cpn > 0)
        {
            for (j = i; j < fix_leg->num_cpn - 1; ++j)
            {
                temp =
                    swp_f_df(today, fix_leg->cpn[j].pay_date, yc) * (fix_not[j] - fix_not[j + 1]);
                fix_leg_pv += temp;
            }

            temp = swp_f_df(today, fix_leg->cpn[fix_leg->num_cpn - 1].pay_date, yc) *
                   fix_not[fix_leg->num_cpn - 1];
            fix_leg_pv += temp;
        }
        else
        {
            fin_not_date = fund_pay[fund_ncpn - 1];
            if (fin_not_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fin_not_date, yc) * fix_not[fund_leg->num_cpn - 1];
                fix_leg_pv += temp;
            }
        }

        /*	Initial */
        if (start_date >= today + eod_pay_flag)
        {
            fix_leg_pv -= swp_f_df(today, start_date, yc) * fix_not[0];
        }
    }

    /*	4)	If there is at least one call after today, value call feature */
    if (call_feat == 1)
    {
        smessage("Launching adi, time steps requested: %d, actual: %d", req_stp, adi_arg->nstp);
        // NB: most cases fix rates are the same! but might not be correct for some cases!!! shall
        // come back later ...
        szErr = _price_genmidat(
            am,
            und,
            adi_arg,
            today,
            fix_rate[0],
            &call,
            (nOutPutExProb ? ppdExProb : 0),
            (nOutPutExProb ? ppdAvgSwapRate : 0),
            pnEx);

        if (szErr)
        {
            goto FREE_RETURN;
        }
    }
    else
    {
        call = 0.0;
    }

    *fund_val = fund_leg_pv;
    *fix_val  = fix_leg_pv;
    *call_val = call;

    if (export_ts)
    {
        am_copy_und(und, und_exp);
    }

FREE_RETURN:

    if (free_struct)
    {
        am_free_all_struct(am, und, call_feat, adi_arg);
    }
    if (am)
        free(am);
    if (und)
        free(und);
    if (adi_arg)
        free(adi_arg);

    return szErr;
}

Err am_caller_AmortMidatAutocal_reduc(
    /*	Today's date */
    long today,

    /*	The underlying */
    int use_calib, /*	0: use lgm2fund, 1: calibrate */

    /*		if calib */
    char*  yc,         /*	yc */
    char*  vc,         /*	vc */
    char*  ref,        /*	ref rate (only if calib) */
    char*  swap_freq,  /*	swap freq (only if calib) */
    char*  swap_basis, /*	swap basis (only if calib) */
    double lambda,     /*	lambda if unique */
    double alpha,      /*	alpha */
    double gamma,      /*	gamma */
    double rho,        /*	rho */

    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */
    long start_date,  /*	Date at which initial notional exchange occurs */
    long theoEndDate, /* End Date of the structure */

    /*		funding */
    char*   fund_ref,
    int     fund_ccy,    /*	0: domestic, 1: other */
    double* fund_not,    /*	If different from domestic or foreign (fund_ccy = 1) */
    char*   fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 1) */
    double  fx_fund_dom, /*	If different from domestic or foreign (fund_ccy = 1) 2 bd fwd */
    long    fx_fund_dom_spot_date,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    double* fund_fix_cpn, /*	Past coupon fixing if relevant,
                                          includes spr, but not mrg, cvg and notional */
    /*		fix */
    double* fix_not,
    int     fix_ncpn,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee, /*	Exercise Fee */

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    double mintime,
    double mininterval,
    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    /*	Function to get correlation stored from a vol cube */
    Err (*get_correl)(
        char* correl_cube_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,

    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                       to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_pay_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Exercised flag */
    int    exercised,   /*	Flag */
    long   ex_date_ex,  /*	Date when exercised */
    long   ex_date_set, /*	Corresponding settlement date */
    double ex_fee,      /*	Corresponding fee */

    int nReduCalibPoints,
    int nOutPutExProb,

    /*	Results */
    double*  fund_val, /*	Value of the funding leg */
    double*  fix_val,  /*	Value of the Power Dual leg */
    double*  call_val, /*	Value of the callable feature */
    double** ppdExProb,
    double** ppdAvgSwapRate,
    int*     pnEx,
    int      export_ts, /*	1: Export TS, 0: don't */
    AM_UND   und_exp)
{
    am_str*     am      = NULL;
    am_und*     und     = NULL;
    am_adi_arg* adi_arg = NULL;

    AM_FUND_LEG fund_leg;
    AM_FUND_CPN fund_cpn;
    AM_FIX_LEG  fix_leg;
    AM_FIX_CPN  fix_cpn;

    int          call_feat;
    double       fund_leg_pv, fix_leg_pv;
    SrtBasisCode bas;

    double df, coupon;

    double call;  //,iv;
    int    i, j;
    int    free_struct = 0;

    double temp;

    int    for_fund, nEx;
    long   fund_start_date, fin_not_date;
    double eq_final_ex, eq_init_ex;

    Err szErr = NULL;

    /*	If exercised */
    if (exercised)
    {
        i = 0;
        while (i < fix_ncpn && fix_start[i] < ex_date_ex)
        {
            i++;
        }
        fix_ncpn = i;

        /*	Structure is called before start: return 0 */
        if (fix_ncpn == 0)
        {
            *fund_val = *fix_val = *call_val = 0.0;
            return NULL;
        }

        i = 0;
        while (i < fund_ncpn && fund_start[i] < ex_date_ex)
        {
            i++;
        }
        fund_ncpn = i;

        ncall     = 0;
        exercised = 0;

        szErr = am_caller_AmortMidatAutocal(
            today,
            use_calib,
            yc,
            vc,
            ref,
            swap_freq,
            swap_basis,
            lambda,
            alpha,
            gamma,
            rho,
            get_cash_vol,
            lgm2dund,
            start_date,
            theoEndDate,
            fund_ref,
            fund_ccy,
            fund_not,
            fund_ccy_yc,
            fx_fund_dom,
            fx_fund_dom_spot_date,

            fund_ncpn,
            fund_fix,
            fund_start,
            fund_end,
            fund_pay,
            fund_basis,
            fund_spr,
            fund_mrg,
            fund_fix_cpn,

            fix_not,
            fix_ncpn,
            fix_start,
            fix_end,
            fix_pay,
            fix_basis,
            fix_rate,
            fix_fee,

            ncall,
            pay_rec,
            ex_date,
            set_date,
            fee,

            req_stp,
            req_stpx,

            mintime,
            mininterval,
            notperiod,
            one2F,
            use_jump,
            max_var_jump,
            strike_type,
            european_model,

            get_correl,
            CorrelName,

            max_std_short,
            fix_lambda,
            one_f_equi,
            skip_last,

            eod_fix_flag,
            eod_pay_flag,
            eod_ex_flag,

            exercised,
            ex_date_ex,
            ex_date_set,
            ex_fee,

            nReduCalibPoints,
            nOutPutExProb,

            fund_val,
            fix_val,
            call_val,
            ppdExProb,
            ppdAvgSwapRate,
            &nEx,
            export_ts,
            und_exp);

        if (szErr)
        {
            goto FREE_RETURN;
        }

        if (ex_date_set >= today + eod_pay_flag)
        {
            *call_val = -ex_fee * swp_f_df(today, ex_date_set, yc);
        }

        goto FREE_RETURN;
    }

    if (fund_ccy == 1)
    {
        fund_ccy = 0;
        for_fund = 1;
        szErr    = convert_funding_to_domestic_amort(
            today,
            start_date,
            eod_fix_flag,
            eod_pay_flag,
            fx_fund_dom,
            fx_fund_dom_spot_date,
            yc,
            fix_ncpn,
            fix_not,
            fund_ncpn,
            fund_fix,
            fund_start,
            fund_pay,
            fund_basis,
            fund_ccy_yc,
            fund_not,
            fund_spr,
            fund_mrg,
            fund_fix_cpn,
            &fund_start_date,
            &eq_final_ex,
            &eq_init_ex);
        if (szErr)
        {
            return szErr;
        }
    }
    else
    {
        for_fund = 0;
    }

    am      = calloc(1, sizeof(am_str));
    und     = calloc(1, sizeof(am_und));
    adi_arg = calloc(1, sizeof(am_adi_arg));

    if (!am || !und || !adi_arg)
    {
        szErr = "memory allocation failure in am_caller";
        goto FREE_RETURN;
    }

    /*	Initialise structures */
    free_struct = 0;
    szErr       = am_fill_check_all_struct_AmortMidatAutocal_reduc(
        today,
        theoEndDate,

        use_calib,
        yc,
        vc,
        ref,
        swap_freq,
        swap_basis,
        lambda,
        alpha,
        gamma,
        rho,

        lgm2dund,

        fund_ref,
        fund_not,
        fund_ncpn,
        fund_fix,
        fund_start,
        fund_end,
        fund_pay,
        fund_basis,
        fund_spr,
        fund_mrg,

        fix_not,
        fix_ncpn,
        fix_start,
        fix_end,
        fix_pay,
        fix_basis,
        fix_rate,
        fix_fee,

        ncall,
        pay_rec,
        ex_date,
        set_date,
        fee,

        req_stp,
        req_stpx,

        get_cash_vol,

        mintime,
        mininterval,
        notperiod,
        one2F,
        use_jump,
        max_var_jump,
        strike_type,
        european_model,

        get_correl,
        CorrelName,
        /*
                        GetVolForBadr,
                        cVolType,
        */
        max_std_short,
        fix_lambda,
        one_f_equi,
        skip_last,

        eod_fix_flag,
        eod_ex_flag,

        am,
        und,

        &call_feat,
        adi_arg);

    if (szErr)
    {
        goto FREE_RETURN;
    }
    free_struct = 1;

    /*	1)	Value funding leg */
    fund_leg    = am->fund_leg;
    fund_leg_pv = 0.0;

    /*	Cash libor */
    if (fund_leg->num_cpn > 0)
    {
        if (for_fund)
        {
            for (j = 0; j < fund_leg->num_cpn - 1; ++j)
            {
                if (fund_leg->cpn[j].pay_date >= today + eod_pay_flag)
                {
                    temp = swp_f_df(today, fund_leg->cpn[j].pay_date, yc) *
                           (fund_not[j] - fund_not[j + 1]);
                    fund_leg_pv += temp;
                }
            }

            if (fund_leg->cpn[fund_leg->num_cpn - 1].pay_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fund_leg->cpn[fund_leg->num_cpn - 1].pay_date, yc) *
                       fund_not[fund_leg->num_cpn - 1];
                fund_leg_pv += temp;
            }
        }
        else
        {
            fund_leg_pv += swp_f_df(today, fund_leg->cpn[0].start_date, yc) * fund_leg->notional[0];
        }
    }

    /*	Coupons: spread + margin */
    for (i = 0; i < fund_leg->num_cpn; i++)
    {
        fund_cpn = fund_leg->cpn + i;
        temp     = swp_f_df(today, fund_cpn->pay_date, yc) * fund_cpn->cpn;
        fund_leg_pv += temp;
    }

    /*	PV of coupons fixed in the past and not yet paid */
    i = 0;
    while (i < fund_ncpn && fund_fix[i] < today + eod_fix_flag)
    {
        if (fund_pay[i] >= today + eod_pay_flag)
        {
            szErr = interp_basis(fund_basis[i], &bas);
            if (szErr)
            {
                goto FREE_RETURN;
            }

            fund_leg_pv += (fund_fix_cpn[i] + fund_mrg[i]) *
                           coverage(fund_start[i], fund_pay[i], bas) * fund_not[0] *
                           swp_f_df(today, fund_pay[i], yc);
        }

        i++;
    }

    /*	2)	Value fix leg */

    fix_leg    = am->fix_leg;
    fix_leg_pv = 0.0;
    // fix_fee = 0.;

    /*	Coupons */
    for (i = 0; i < fix_leg->num_cpn; i++)
    {
        fix_cpn = fix_leg->cpn + i;

        /*	Discount */
        df = swp_f_df(today, fix_cpn->pay_date, yc);

        coupon = fix_cpn->cpn;

        /*	Coupon pv */
        temp = df * coupon;
        fix_leg_pv += temp;
    }

    /*	PV of coupons of past periods and not yet paid */
    i = 0;
    //	while (i < fix_ncpn && fix_fix[i] < today + eod_fix_flag)
    while (i < fix_ncpn && fix_start[i] < today + eod_fix_flag)
    {
        if (fix_pay[i] >= today + eod_pay_flag)
        {
            szErr = interp_basis(fix_basis[i], &bas);
            if (szErr)
            {
                goto FREE_RETURN;
            }

            //			fix_leg_pv += fix_fix_cpn[i]
            fix_leg_pv += fix_rate[i] * coverage(fix_start[i], fix_pay[i], bas) * fix_not[i] *
                          swp_f_df(today, fix_pay[i], yc);
        }

        i++;
    }

    /*	Initial and final exchange */
    if (for_fund)
    {
        /*	Final */
        if (fix_leg->num_cpn > 0)
        {
            for (j = i; j < fix_leg->num_cpn - 1; ++j)
            {
                temp =
                    swp_f_df(today, fix_leg->cpn[j].pay_date, yc) * (fix_not[j] - fix_not[j + 1]);
                fix_leg_pv += temp;
            }

            temp = swp_f_df(today, fix_leg->cpn[fix_leg->num_cpn - 1].pay_date, yc) *
                   fix_not[fix_leg->num_cpn - 1];
            fix_leg_pv += temp;
        }
        else
        {
            fin_not_date = fund_pay[fund_ncpn - 1];
            if (fin_not_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fin_not_date, yc) * fix_not[fund_leg->num_cpn - 1];
                fix_leg_pv += temp;
            }
        }

        /*	Initial */
        if (start_date >= today + eod_pay_flag)
        {
            fix_leg_pv -= swp_f_df(today, start_date, yc) * fix_not[0];
        }
    }

    /*	4)	If there is at least one call after today, value call feature */
    if (call_feat == 1)
    {
        smessage("Launching adi, time steps requested: %d, actual: %d", req_stp, adi_arg->nstp);
        {
            /// 4. Compute exercise probabilities based on the "equivalent" term structure
            double* pdEquiExProb = _alloca((*pnEx + 1) * sizeof(double));
            double* pdAvgStrike  = _alloca((*pnEx) * sizeof(double));
            szErr                = _price_genmidat(
                am, und, adi_arg, today, fix_rate[0], &call, &pdEquiExProb, &pdAvgStrike, pnEx);
            /// 5. Determine the "most probable" equivalent strike swaptions

            /// 6. Get Market price for the corresponding "most probable" amortizing swaptions,
            /// while keeping the rest equivalent strikes swaptions untouched
            /// amortMidat_compute_diagonal_prices_new_equistrike

            /// 7. Calibrate again with this new set of prices
            // amortMidat_cpd_calib_diagonal_new

            /// 8. price again with then new term structure
            szErr = _price_genmidat(
                am, und, adi_arg, today, fix_rate[0], &call, &pdEquiExProb, &pdAvgStrike, pnEx);
        }

        if (szErr)
        {
            goto FREE_RETURN;
        }
    }
    else
    {
        call = 0.0;
    }

    *fund_val = fund_leg_pv;
    *fix_val  = fix_leg_pv;
    *call_val = call;

    if (export_ts)
    {
        am_copy_und(und, und_exp);
    }

FREE_RETURN:

    if (free_struct)
    {
        am_free_all_struct(am, und, call_feat, adi_arg);
    }
    if (am)
        free(am);
    if (und)
        free(und);
    if (adi_arg)
        free(adi_arg);

    return szErr;
}

Err am_caller_ts_AmortMidatAutocal(
    /*	Today's date */
    long today,

    /*	The underlying */
    int use_calib, /*	0: use lgm2fund, 1: calibrate */

    /*		if calib */
    char*   yc,            /*	yc */
    char*   vc,            /*	vc */
    char*   ref,           /*	ref rate (only if calib) */
    char*   swap_freq,     /*	swap freq (only if calib) */
    char*   swap_basis,    /*	swap basis (only if calib) */
    double* pdLambdaValue, /*	lambda if unique */
    double* pdLambdaTime,
    int     nLambdaSize,
    double  alpha, /*	alpha */
    double  gamma, /*	gamma */
    double  rho,   /*	rho */

    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */
    long start_date,  /*	Date at which initial notional exchange occurs */
    long theoEndDate, /* End Date of the structure */

    /*		funding */
    char*   fund_ref,
    int     fund_ccy,    /*	0: domestic, 1: other */
    double* fund_not,    /*	If different from domestic or foreign (fund_ccy = 1) */
    char*   fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 1) */
    double  fx_fund_dom, /*	If different from domestic or foreign (fund_ccy = 1) 2 bd fwd */
    long    fx_fund_dom_spot_date,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    double* fund_fix_cpn, /*	Past coupon fixing if relevant,
                                          includes spr, but not mrg, cvg and notional */
    /*		fix */
    double* fix_not,
    int     fix_ncpn,
    //			long		*fix_fix,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee, /*	Exercise Fee */
                     //			double		*fix_fix_cpn,			/*	Past coupon fixing if relevant
                     //*/

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    double mintime,
    double mininterval,
    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    /*	Function to get correlation stored from a vol cube */
    Err (*get_correl)(
        char* correl_cube_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,
    /*
                            Err (*GetVolForBadr)( Date, Date, double, SRT_Boolean, double *),
                            char *cVolType,
    */
    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                       to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_pay_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Exercised flag */
    int    exercised,   /*	Flag */
    long   ex_date_ex,  /*	Date when exercised */
    long   ex_date_set, /*	Corresponding settlement date */
    double ex_fee,      /*	Corresponding fee */

    int nReduCalibPoints,
    int nOutPutExProb,

    /*	Results */
    double*  fund_val, /*	Value of the funding leg */
    double*  fix_val,  /*	Value of the Power Dual leg */
    double*  call_val, /*	Value of the callable feature */
    double** ppdExProb,
    double** ppdAvgSwapRate,  // averag swap rate at the boundary
    int*     pnEx,
    int      export_ts, /*	1: Export TS, 0: don't */
    AM_UND   und_exp)
{
    am_str*     am      = NULL;
    am_und*     und     = NULL;
    am_adi_arg* adi_arg = NULL;

    AM_FUND_LEG fund_leg;
    AM_FUND_CPN fund_cpn;
    AM_FIX_LEG  fix_leg;
    AM_FIX_CPN  fix_cpn;

    int          call_feat;
    double       fund_leg_pv, fix_leg_pv;
    SrtBasisCode bas;

    double df, coupon;

    double call;
    int    i, j;
    int    free_struct = 0;

    double temp;

    int    for_fund, nEx;
    long   fund_start_date, fin_not_date;
    double eq_final_ex, eq_init_ex;

    Err szErr = NULL;

    /*	If exercised */
    if (exercised)
    {
        i = 0;
        while (i < fix_ncpn && fix_start[i] < ex_date_ex)
        {
            i++;
        }
        fix_ncpn = i;

        /*	Structure is called before start: return 0 */
        if (fix_ncpn == 0)
        {
            *fund_val = *fix_val = *call_val = 0.0;
            return NULL;
        }

        i = 0;
        while (i < fund_ncpn && fund_start[i] < ex_date_ex)
        {
            i++;
        }
        fund_ncpn = i;

        ncall     = 0;
        exercised = 0;

        szErr = am_caller_ts_AmortMidatAutocal(
            today,
            use_calib,
            yc,
            vc,
            ref,
            swap_freq,
            swap_basis,
            pdLambdaValue,
            pdLambdaTime,
            nLambdaSize,

            alpha,
            gamma,
            rho,
            get_cash_vol,
            lgm2dund,
            start_date,
            theoEndDate,
            fund_ref,
            fund_ccy,
            fund_not,
            fund_ccy_yc,
            fx_fund_dom,
            fx_fund_dom_spot_date,

            fund_ncpn,
            fund_fix,
            fund_start,
            fund_end,
            fund_pay,
            fund_basis,
            fund_spr,
            fund_mrg,
            fund_fix_cpn,

            fix_not,
            fix_ncpn,
            fix_start,
            fix_end,
            fix_pay,
            fix_basis,
            fix_rate,
            fix_fee,

            ncall,
            pay_rec,
            ex_date,
            set_date,
            fee,

            req_stp,
            req_stpx,

            mintime,
            mininterval,
            notperiod,
            one2F,
            use_jump,
            max_var_jump,
            strike_type,
            european_model,

            get_correl,
            CorrelName,

            max_std_short,
            fix_lambda,
            one_f_equi,
            skip_last,

            eod_fix_flag,
            eod_pay_flag,
            eod_ex_flag,

            exercised,
            ex_date_ex,
            ex_date_set,
            ex_fee,
            nReduCalibPoints,
            nOutPutExProb,

            fund_val,
            fix_val,
            call_val,
            ppdExProb,
            ppdAvgSwapRate,
            &nEx,
            export_ts,
            und_exp);

        if (szErr)
        {
            goto FREE_RETURN;
        }

        if (ex_date_set >= today + eod_pay_flag)
        {
            *call_val = -ex_fee * swp_f_df(today, ex_date_set, yc);
        }

        goto FREE_RETURN;
    }

    if (fund_ccy == 1)
    {
        fund_ccy = 0;
        for_fund = 1;
        szErr    = convert_funding_to_domestic_amort(
            today,
            start_date,
            eod_fix_flag,
            eod_pay_flag,
            fx_fund_dom,
            fx_fund_dom_spot_date,
            yc,
            fix_ncpn,
            fix_not,
            fund_ncpn,
            fund_fix,
            fund_start,
            fund_pay,
            fund_basis,
            fund_ccy_yc,
            fund_not,
            fund_spr,
            fund_mrg,
            fund_fix_cpn,
            &fund_start_date,
            &eq_final_ex,
            &eq_init_ex);
        if (szErr)
        {
            return szErr;
        }
    }
    else
    {
        for_fund = 0;
    }

    am      = calloc(1, sizeof(am_str));
    und     = calloc(1, sizeof(am_und));
    adi_arg = calloc(1, sizeof(am_adi_arg));

    if (!am || !und || !adi_arg)
    {
        szErr = "memory allocation failure in am_caller";
        goto FREE_RETURN;
    }

    /*	Initialise structures */
    free_struct = 0;
    szErr       = am_fill_check_all_struct_ts_AmortMidatAutocal(
        today,
        theoEndDate,

        use_calib,
        yc,
        vc,
        ref,
        swap_freq,
        swap_basis,

        nLambdaSize,
        pdLambdaTime,
        pdLambdaValue,

        alpha,
        gamma,
        rho,

        lgm2dund,

        fund_ref,
        fund_not,
        fund_ncpn,
        fund_fix,
        fund_start,
        fund_end,
        fund_pay,
        fund_basis,
        fund_spr,
        fund_mrg,

        fix_not,
        fix_ncpn,
        fix_start,
        fix_end,
        fix_pay,
        fix_basis,
        fix_rate,
        fix_fee,

        ncall,
        pay_rec,
        ex_date,
        set_date,
        fee,

        req_stp,
        req_stpx,

        get_cash_vol,

        mintime,
        mininterval,
        notperiod,
        one2F,
        use_jump,
        max_var_jump,
        strike_type,
        european_model,

        get_correl,
        CorrelName,
        /*
                        GetVolForBadr,
                        cVolType,
        */
        max_std_short,
        fix_lambda,
        one_f_equi,
        skip_last,

        eod_fix_flag,
        eod_ex_flag,

        am,
        und,

        &call_feat,
        adi_arg);

    if (szErr)
    {
        goto FREE_RETURN;
    }
    free_struct = 1;

    /*	1)	Value funding leg */

    fund_leg    = am->fund_leg;
    fund_leg_pv = 0.0;

    /*	Cash libor */
    if (fund_leg->num_cpn > 0)
    {
        if (for_fund)
        {
            for (j = 0; j < fund_leg->num_cpn - 1; ++j)
            {
                if (fund_leg->cpn[j].pay_date >= today + eod_pay_flag)
                {
                    temp = swp_f_df(today, fund_leg->cpn[j].pay_date, yc) *
                           (fund_not[j] - fund_not[j + 1]);
                    fund_leg_pv += temp;
                }
            }

            if (fund_leg->cpn[fund_leg->num_cpn - 1].pay_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fund_leg->cpn[fund_leg->num_cpn - 1].pay_date, yc) *
                       fund_not[fund_leg->num_cpn - 1];
                fund_leg_pv += temp;
            }
        }
        else
        {
            fund_leg_pv += swp_f_df(today, fund_leg->cpn[0].start_date, yc) * fund_leg->notional[0];
        }
    }

    /*	Coupons: spread + margin */
    for (i = 0; i < fund_leg->num_cpn; i++)
    {
        fund_cpn = fund_leg->cpn + i;
        temp     = swp_f_df(today, fund_cpn->pay_date, yc) * fund_cpn->cpn;
        fund_leg_pv += temp;
    }

    /*	PV of coupons fixed in the past and not yet paid */
    i = 0;
    while (i < fund_ncpn && fund_fix[i] < today + eod_fix_flag)
    {
        if (fund_pay[i] >= today + eod_pay_flag)
        {
            szErr = interp_basis(fund_basis[i], &bas);
            if (szErr)
            {
                goto FREE_RETURN;
            }

            fund_leg_pv += (fund_fix_cpn[i] + fund_mrg[i]) *
                           coverage(fund_start[i], fund_pay[i], bas) * fund_not[0] *
                           swp_f_df(today, fund_pay[i], yc);
        }

        i++;
    }

    /*	2)	Value fix leg */

    fix_leg    = am->fix_leg;
    fix_leg_pv = 0.0;

    /*	Coupons */
    for (i = 0; i < fix_leg->num_cpn; i++)
    {
        fix_cpn = fix_leg->cpn + i;

        /*	Discount */
        df = swp_f_df(today, fix_cpn->pay_date, yc);

        coupon = fix_cpn->cpn;

        /*	Coupon pv */
        temp = df * coupon;
        fix_leg_pv += temp;
    }

    /*	PV of coupons of past periods and not yet paid */
    i = 0;
    //	while (i < fix_ncpn && fix_fix[i] < today + eod_fix_flag)
    while (i < fix_ncpn && fix_start[i] < today + eod_fix_flag)
    {
        if (fix_pay[i] >= today + eod_pay_flag)
        {
            szErr = interp_basis(fix_basis[i], &bas);
            if (szErr)
            {
                goto FREE_RETURN;
            }

            //			fix_leg_pv += fix_fix_cpn[i]
            fix_leg_pv += fix_rate[i] * coverage(fix_start[i], fix_pay[i], bas) * fix_not[i] *
                          swp_f_df(today, fix_pay[i], yc);
        }

        i++;
    }

    /*	Initial and final exchange */
    if (for_fund)
    {
        /*	Final */
        if (fix_leg->num_cpn > 0)
        {
            for (j = i; j < fix_leg->num_cpn - 1; ++j)
            {
                temp =
                    swp_f_df(today, fix_leg->cpn[j].pay_date, yc) * (fix_not[j] - fix_not[j + 1]);
                fix_leg_pv += temp;
            }

            temp = swp_f_df(today, fix_leg->cpn[fix_leg->num_cpn - 1].pay_date, yc) *
                   fix_not[fix_leg->num_cpn - 1];
            fix_leg_pv += temp;
        }
        else
        {
            fin_not_date = fund_pay[fund_ncpn - 1];
            if (fin_not_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fin_not_date, yc) * fix_not[fund_leg->num_cpn - 1];
                fix_leg_pv += temp;
            }
        }

        /*	Initial */
        if (start_date >= today + eod_pay_flag)
        {
            fix_leg_pv -= swp_f_df(today, start_date, yc) * fix_not[0];
        }
    }

    /*	4)	If there is at least one call after today, value call feature */

    if (call_feat == 1)
    {
        smessage("Launching adi, time steps requested: %d, actual: %d", req_stp, adi_arg->nstp);
        szErr = _price_genmidat_ts(
            am,
            und,
            adi_arg,
            today,
            fix_rate[0],
            pdLambdaValue,
            pdLambdaTime,
            nLambdaSize,
            &call,
            (nOutPutExProb ? ppdExProb : 0),
            (nOutPutExProb ? ppdAvgSwapRate : 0),
            pnEx);

        if (szErr)
        {
            goto FREE_RETURN;
        }
    }
    else
    {
        call = 0.0;
    }

    *fund_val = fund_leg_pv;
    *fix_val  = fix_leg_pv;
    *call_val = call;

    if (export_ts)
    {
        am_copy_und(und, und_exp);
    }

FREE_RETURN:

    if (free_struct)
    {
        am_free_all_struct(am, und, call_feat, adi_arg);
    }
    if (am)
        free(am);
    if (und)
        free(und);
    if (adi_arg)
        free(adi_arg);

    return szErr;
}

Err am_caller_ts_AmortMidatAutocal_reduc(
    /*	Today's date */
    long today,

    /*	The underlying */
    int use_calib, /*	0: use lgm2fund, 1: calibrate */

    /*		if calib */
    char*   yc,            /*	yc */
    char*   vc,            /*	vc */
    char*   ref,           /*	ref rate (only if calib) */
    char*   swap_freq,     /*	swap freq (only if calib) */
    char*   swap_basis,    /*	swap basis (only if calib) */
    double* pdLambdaValue, /*	lambda if unique */
    double* pdLambdaTime,
    int     nLambdaSize,
    double  alpha, /*	alpha */
    double  gamma, /*	gamma */
    double  rho,   /*	rho */

    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    /*		if no calilb */
    char* lgm2dund,

    /*	The structure */
    long start_date,  /*	Date at which initial notional exchange occurs */
    long theoEndDate, /* End Date of the structure */

    /*		funding */
    char*   fund_ref,
    int     fund_ccy,    /*	0: domestic, 1: other */
    double* fund_not,    /*	If different from domestic or foreign (fund_ccy = 1) */
    char*   fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 1) */
    double  fx_fund_dom, /*	If different from domestic or foreign (fund_ccy = 1) 2 bd fwd */
    long    fx_fund_dom_spot_date,
    int     fund_ncpn,
    long*   fund_fix,
    long*   fund_start,
    long*   fund_end,
    long*   fund_pay,
    char**  fund_basis,
    double* fund_spr,
    double* fund_mrg,
    double* fund_fix_cpn, /*	Past coupon fixing if relevant,
                                          includes spr, but not mrg, cvg and notional */
    /*		fix */
    double* fix_not,
    int     fix_ncpn,
    //			long		*fix_fix,
    long*   fix_start,
    long*   fix_end,
    long*   fix_pay,
    char**  fix_basis,
    double* fix_rate,
    double* fix_fee, /*	Exercise Fee */
                     //			double		*fix_fix_cpn,			/*	Past coupon fixing if relevant
                     //*/

    /*		calls */
    int     ncall,
    int     pay_rec, /*	0: rec pd, 1: pay pd */
    long*   ex_date,
    long*   set_date,
    double* fee,

    /*	Numerical params */
    int req_stp,
    int req_stpx,

    /*	Calib params */
    double mintime,
    double mininterval,
    int    notperiod,
    int    one2F,
    int    use_jump,
    double max_var_jump,
    int    strike_type,
    int    european_model,

    /*	Function to get correlation stored from a vol cube */
    Err (*get_correl)(
        char* correl_cube_name, double start_date, double end_date, double strike, double* vol),
    char* CorrelName,
    /*
                            Err (*GetVolForBadr)( Date, Date, double, SRT_Boolean, double *),
                            char *cVolType,
    */
    double max_std_short,
    int    fix_lambda, /*	0: calib lambda to cap, 1: fix lambda calib
                                                       to diagonal */
    int one_f_equi,    /*	1F equivalent flag:
                                                       if set to 1, then 2F lambda will calibrate
                                                       to the cap priced within calibrated 1F
                                                       with the given lambda */
    int skip_last,     /*	If 1, the last option is disregarded
                                                       and the forward volatility is flat from option
                                                       n-1 */

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I, 1: E */
    int eod_pay_flag, /*	0: I, 1: E */
    int eod_ex_flag,  /*	0: I, 1: E */

    /*	Exercised flag */
    int    exercised,   /*	Flag */
    long   ex_date_ex,  /*	Date when exercised */
    long   ex_date_set, /*	Corresponding settlement date */
    double ex_fee,      /*	Corresponding fee */

    int nReduCalibPoints,
    int nOutPutExProb,

    /*	Results */
    double*  fund_val, /*	Value of the funding leg */
    double*  fix_val,  /*	Value of the Power Dual leg */
    double*  call_val, /*	Value of the callable feature */
    double** ppdExProb,
    double** ppdAvgSwapRate,
    int*     pnEx,
    int      export_ts, /*	1: Export TS, 0: don't */
    AM_UND   und_exp)
{
    am_str*     am      = NULL;
    am_und*     und     = NULL;
    am_adi_arg* adi_arg = NULL;

    AM_FUND_LEG fund_leg;
    AM_FUND_CPN fund_cpn;
    AM_FIX_LEG  fix_leg;
    AM_FIX_CPN  fix_cpn;

    int          call_feat;
    double       fund_leg_pv, fix_leg_pv;
    SrtBasisCode bas;

    double df, coupon;

    double call;
    int    i, j;
    int    free_struct = 0;

    double temp;

    int    for_fund, nEx;
    long   fund_start_date, fin_not_date;
    double eq_final_ex, eq_init_ex;

    Err szErr = NULL;

    /*	If exercised */
    if (exercised)
    {
        i = 0;
        while (i < fix_ncpn && fix_start[i] < ex_date_ex)
        {
            i++;
        }
        fix_ncpn = i;

        /*	Structure is called before start: return 0 */
        if (fix_ncpn == 0)
        {
            *fund_val = *fix_val = *call_val = 0.0;
            return NULL;
        }

        i = 0;
        while (i < fund_ncpn && fund_start[i] < ex_date_ex)
        {
            i++;
        }
        fund_ncpn = i;

        ncall     = 0;
        exercised = 0;

        szErr = am_caller_ts_AmortMidatAutocal(
            today,
            use_calib,
            yc,
            vc,
            ref,
            swap_freq,
            swap_basis,
            pdLambdaValue,
            pdLambdaTime,
            nLambdaSize,

            alpha,
            gamma,
            rho,
            get_cash_vol,
            lgm2dund,
            start_date,
            theoEndDate,
            fund_ref,
            fund_ccy,
            fund_not,
            fund_ccy_yc,
            fx_fund_dom,
            fx_fund_dom_spot_date,

            fund_ncpn,
            fund_fix,
            fund_start,
            fund_end,
            fund_pay,
            fund_basis,
            fund_spr,
            fund_mrg,
            fund_fix_cpn,

            fix_not,
            fix_ncpn,
            fix_start,
            fix_end,
            fix_pay,
            fix_basis,
            fix_rate,
            fix_fee,

            ncall,
            pay_rec,
            ex_date,
            set_date,
            fee,

            req_stp,
            req_stpx,

            mintime,
            mininterval,
            notperiod,
            one2F,
            use_jump,
            max_var_jump,
            strike_type,
            european_model,

            get_correl,
            CorrelName,

            max_std_short,
            fix_lambda,
            one_f_equi,
            skip_last,

            eod_fix_flag,
            eod_pay_flag,
            eod_ex_flag,

            exercised,
            ex_date_ex,
            ex_date_set,
            ex_fee,
            nReduCalibPoints,
            nOutPutExProb,

            fund_val,
            fix_val,
            call_val,
            ppdExProb,
            ppdAvgSwapRate,
            &nEx,
            export_ts,
            und_exp);

        if (szErr)
        {
            goto FREE_RETURN;
        }

        if (ex_date_set >= today + eod_pay_flag)
        {
            *call_val = -ex_fee * swp_f_df(today, ex_date_set, yc);
        }

        goto FREE_RETURN;
    }

    if (fund_ccy == 1)
    {
        fund_ccy = 0;
        for_fund = 1;
        szErr    = convert_funding_to_domestic_amort(
            today,
            start_date,
            eod_fix_flag,
            eod_pay_flag,
            fx_fund_dom,
            fx_fund_dom_spot_date,
            yc,
            fix_ncpn,
            fix_not,
            fund_ncpn,
            fund_fix,
            fund_start,
            fund_pay,
            fund_basis,
            fund_ccy_yc,
            fund_not,
            fund_spr,
            fund_mrg,
            fund_fix_cpn,
            &fund_start_date,
            &eq_final_ex,
            &eq_init_ex);
        if (szErr)
        {
            return szErr;
        }
    }
    else
    {
        for_fund = 0;
    }

    am      = calloc(1, sizeof(am_str));
    und     = calloc(1, sizeof(am_und));
    adi_arg = calloc(1, sizeof(am_adi_arg));

    if (!am || !und || !adi_arg)
    {
        szErr = "memory allocation failure in am_caller";
        goto FREE_RETURN;
    }

    /*	Initialise structures */
    free_struct = 0;
    szErr       = am_fill_check_all_struct_ts_AmortMidatAutocal(
        today,
        theoEndDate,

        use_calib,
        yc,
        vc,
        ref,
        swap_freq,
        swap_basis,

        nLambdaSize,
        pdLambdaTime,
        pdLambdaValue,

        alpha,
        gamma,
        rho,

        lgm2dund,

        fund_ref,
        fund_not,
        fund_ncpn,
        fund_fix,
        fund_start,
        fund_end,
        fund_pay,
        fund_basis,
        fund_spr,
        fund_mrg,

        fix_not,
        fix_ncpn,
        fix_start,
        fix_end,
        fix_pay,
        fix_basis,
        fix_rate,
        fix_fee,

        ncall,
        pay_rec,
        ex_date,
        set_date,
        fee,

        req_stp,
        req_stpx,

        get_cash_vol,

        mintime,
        mininterval,
        notperiod,
        one2F,
        use_jump,
        max_var_jump,
        strike_type,
        european_model,

        get_correl,
        CorrelName,
        /*
                        GetVolForBadr,
                        cVolType,
        */
        max_std_short,
        fix_lambda,
        one_f_equi,
        skip_last,

        eod_fix_flag,
        eod_ex_flag,

        am,
        und,

        &call_feat,
        adi_arg);

    if (szErr)
    {
        goto FREE_RETURN;
    }
    free_struct = 1;

    /*	1)	Value funding leg */

    fund_leg    = am->fund_leg;
    fund_leg_pv = 0.0;

    /*	Cash libor */
    if (fund_leg->num_cpn > 0)
    {
        if (for_fund)
        {
            for (j = 0; j < fund_leg->num_cpn - 1; ++j)
            {
                if (fund_leg->cpn[j].pay_date >= today + eod_pay_flag)
                {
                    temp = swp_f_df(today, fund_leg->cpn[j].pay_date, yc) *
                           (fund_not[j] - fund_not[j + 1]);
                    fund_leg_pv += temp;
                }
            }

            if (fund_leg->cpn[fund_leg->num_cpn - 1].pay_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fund_leg->cpn[fund_leg->num_cpn - 1].pay_date, yc) *
                       fund_not[fund_leg->num_cpn - 1];
                fund_leg_pv += temp;
            }
        }
        else
        {
            fund_leg_pv += swp_f_df(today, fund_leg->cpn[0].start_date, yc) * fund_leg->notional[0];
        }
    }

    /*	Coupons: spread + margin */
    for (i = 0; i < fund_leg->num_cpn; i++)
    {
        fund_cpn = fund_leg->cpn + i;
        temp     = swp_f_df(today, fund_cpn->pay_date, yc) * fund_cpn->cpn;
        fund_leg_pv += temp;
    }

    /*	PV of coupons fixed in the past and not yet paid */
    i = 0;
    while (i < fund_ncpn && fund_fix[i] < today + eod_fix_flag)
    {
        if (fund_pay[i] >= today + eod_pay_flag)
        {
            szErr = interp_basis(fund_basis[i], &bas);
            if (szErr)
            {
                goto FREE_RETURN;
            }

            fund_leg_pv += (fund_fix_cpn[i] + fund_mrg[i]) *
                           coverage(fund_start[i], fund_pay[i], bas) * fund_not[0] *
                           swp_f_df(today, fund_pay[i], yc);
        }

        i++;
    }

    /*	2)	Value fix leg */

    fix_leg    = am->fix_leg;
    fix_leg_pv = 0.0;

    /*	Coupons */
    for (i = 0; i < fix_leg->num_cpn; i++)
    {
        fix_cpn = fix_leg->cpn + i;

        /*	Discount */
        df = swp_f_df(today, fix_cpn->pay_date, yc);

        coupon = fix_cpn->cpn;

        /*	Coupon pv */
        temp = df * coupon;
        fix_leg_pv += temp;
    }

    /*	PV of coupons of past periods and not yet paid */
    i = 0;
    //	while (i < fix_ncpn && fix_fix[i] < today + eod_fix_flag)
    while (i < fix_ncpn && fix_start[i] < today + eod_fix_flag)
    {
        if (fix_pay[i] >= today + eod_pay_flag)
        {
            szErr = interp_basis(fix_basis[i], &bas);
            if (szErr)
            {
                goto FREE_RETURN;
            }

            //			fix_leg_pv += fix_fix_cpn[i]
            fix_leg_pv += fix_rate[i] * coverage(fix_start[i], fix_pay[i], bas) * fix_not[i] *
                          swp_f_df(today, fix_pay[i], yc);
        }

        i++;
    }

    /*	Initial and final exchange */
    if (for_fund)
    {
        /*	Final */
        if (fix_leg->num_cpn > 0)
        {
            for (j = i; j < fix_leg->num_cpn - 1; ++j)
            {
                temp =
                    swp_f_df(today, fix_leg->cpn[j].pay_date, yc) * (fix_not[j] - fix_not[j + 1]);
                fix_leg_pv += temp;
            }

            temp = swp_f_df(today, fix_leg->cpn[fix_leg->num_cpn - 1].pay_date, yc) *
                   fix_not[fix_leg->num_cpn - 1];
            fix_leg_pv += temp;
        }
        else
        {
            fin_not_date = fund_pay[fund_ncpn - 1];
            if (fin_not_date >= today + eod_pay_flag)
            {
                temp = swp_f_df(today, fin_not_date, yc) * fix_not[fund_leg->num_cpn - 1];
                fix_leg_pv += temp;
            }
        }

        /*	Initial */
        if (start_date >= today + eod_pay_flag)
        {
            fix_leg_pv -= swp_f_df(today, start_date, yc) * fix_not[0];
        }
    }

    /*	4)	If there is at least one call after today, value call feature */

    if (call_feat == 1)
    {
        smessage("Launching adi, time steps requested: %d, actual: %d", req_stp, adi_arg->nstp);
        szErr = _price_genmidat_ts(
            am,
            und,
            adi_arg,
            today,
            fix_rate[0],
            pdLambdaValue,
            pdLambdaTime,
            nLambdaSize,
            &call,
            (nOutPutExProb ? ppdExProb : 0),
            (nOutPutExProb ? ppdAvgSwapRate : 0),
            pnEx);

        if (szErr)
        {
            goto FREE_RETURN;
        }
    }
    else
    {
        call = 0.0;
    }

    *fund_val = fund_leg_pv;
    *fix_val  = fix_leg_pv;
    *call_val = call;

    if (export_ts)
    {
        am_copy_und(und, und_exp);
    }

FREE_RETURN:

    if (free_struct)
    {
        am_free_all_struct(am, und, call_feat, adi_arg);
    }
    if (am)
        free(am);
    if (und)
        free(und);
    if (adi_arg)
        free(adi_arg);

    return szErr;
}
