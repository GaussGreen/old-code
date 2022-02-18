/****************************************************************************/
/* Purpose */
/****************************************************************************/

/****************************************************************************/
/* EXTERNAL ROUTINES USED */
/****************************************************************************/
/*

*/
/****************************************************************************/
/* Headers */
/****************************************************************************/
#include "grf_h_all.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_lgmUSprotos.h"
#include "srt_h_lgmUStypes.h"
#include "srt_h_lgmprotos.h"
#include "swp_h_vol.h"

/****************************************************************************/
/* Prototypes of private functions */
/****************************************************************************/
static void   TestOverride(LGMCalParm* CalReqPtr);
static LGMErr copyExerBdry(LGMSwptns* ExerBdryPtr, SrtLgmExerBdryData* lgmExerBdryData);
static LGMErr copyTSData(SigKapTS* SKtsPtr, SrtLgmTSData* TSData);

/* Free arrays that will be used for output, if already allocated */
static void free_outputs(
    int* convertTS,    /* 1=compute new sigs and taus; 0=don't bother */
    int* findExerBdry, /* 1=find swap rates at exercise boundary; 0=don't bother */
    SrtLgmExerBdryData*
        lgmExerBdryData, /* ptr to exercise boundary data structure (NULL => not req'd) */
    SrtLgmRefSwptnData*
        lgmRefSwptnData,     /* ptr to reference swaption data structure (NULL => not req'd) */
    SrtLgmTSData* lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr*    atcTSData);   /* ptr to zeta/G data (NULL => not req'd) */

/* Gets a sabr param for the autocal time swap */
LGMErr GetSABRParam(long dstart, long dend, char* vol_id, SABR_VOL_TYPE component, double* result)
{
    LGMErr err = NULL;
    int    icomponent;
    double strike = 0.05;
    double power;

    switch (component)
    {
    case SABR_ATM_LOG:
        icomponent = 0;
        break;
    case SABR_ATM_NORM:
        icomponent = 1;
        break;
    case SABR_ATM_BETA:
        icomponent = 2;
        break;
    case ALPHA:
        icomponent = 5;
        break;
    case BETA:
        icomponent = 6;
        break;
    case RHO:
        icomponent = 7;
        break;
    }

    err = swp_f_SABRvol(vol_id, dstart, dend, strike, result, &power, icomponent);
    return err;
}

static struct
{
    LGMErr (*GetVol)(long, long, double, SRT_Boolean, double*); /* A cash volatility function */
    char*       YieldCurveName;
    long        lSpotLag;
    SrtCurvePtr YieldCurve;
    LGMMarkConv Conventions;
    long        CalcDate;
    double      MaxStrikeStd;
    long        IsLognormal;
} StaticMarketInfo;

LGMErr SetUpStaticMarketIfo(
    LGMErr (*GetVol)(long, long, double, SRT_Boolean, double*),
    char*  YCName,
    double dMaxStrikeStd,
    char*  char_vol_type)
{
    SrtDiffusionType srt_vol_type;
    LGMErr           err;

    StaticMarketInfo.GetVol         = GetVol;
    StaticMarketInfo.YieldCurveName = YCName;
    StaticMarketInfo.YieldCurve     = lookup_curve(YCName);
    StaticMarketInfo.CalcDate       = get_clcndate_from_yldcrv(StaticMarketInfo.YieldCurve);
    StaticMarketInfo.MaxStrikeStd   = dMaxStrikeStd;

    if (err = interp_diffusion_type(char_vol_type, &srt_vol_type))
        return err;

    if (srt_vol_type == SRT_LOGNORMAL)
        StaticMarketInfo.IsLognormal = 1;
    else
        StaticMarketInfo.IsLognormal = 0;

    return LGMCcyDefaults(StaticMarketInfo.YieldCurve, &StaticMarketInfo.Conventions);
}

/* Wrapper to the input volatility function that cuts the volatility at +/- StdDev above and below
 * the ATM */
/* Only affects the LGM calibration */
LGMErr GetTruncVol(long lStart, long lEnd, double dStrike, SRT_Boolean bUnknown, double* pdVol)
{
    /* local variables */
    SwapDP swap_dp;
    double dForward, dATMvol, dExpiry, dBound, dStdDev;
    LGMErr err;
    long   lFix;

    /* Set up the swap DP */
    if (err = swp_f_setSwapDP(
            lStart,
            lEnd,
            StaticMarketInfo.Conventions.sfreq,
            StaticMarketInfo.Conventions.sbasis,
            &swap_dp))
        return err;
    swap_dp.spot_lag = StaticMarketInfo.Conventions.lag;

    /* Calculate the CASH forward swap rate */
    if (err =
            swp_f_ForwardRate_SwapDP(&swap_dp, StaticMarketInfo.YieldCurveName, "CASH", &dForward))
        return err;

    /* Get the CASH ATM vol */
    if (err = StaticMarketInfo.GetVol(lStart, lEnd, dForward, bUnknown, &dATMvol))
        return err;

    /* Convert the ATM vol to lognormal if necessary */
    lFix    = add_unit(swap_dp.start, -swap_dp.spot_lag, SRT_BDAY, SUCCEEDING);
    dExpiry = (lFix - StaticMarketInfo.CalcDate) * YEARS_IN_DAY;
    if (!StaticMarketInfo.IsLognormal)
        if (err = srt_f_optbetavoltoblkvol(dForward, dForward, dATMvol, dExpiry, 0.0, &dATMvol))
            return err;

    /* truncate the strike if necessary */
    dStdDev = exp(StaticMarketInfo.MaxStrikeStd * dATMvol * dExpiry);
    if (dStrike > (dBound = dForward * dStdDev))
        dStrike = dBound;
    else if (dStrike < (dBound = dForward / dStdDev))
        dStrike = dBound;

    /* return the volatility from the adjusted strike */
    return StaticMarketInfo.GetVol(lStart, lEnd, dStrike, bUnknown, pdVol);
}

// -------------------------------------------------------------------------------------------------------------------------------------------------
// //
//
// LGMCallableTimeSwapCaller
//
LGMErr LGMCallableTimeSwapCaller(
    //	Exercise
    long  nEx,             /* number of exercise dates		*/
    Date* tEx,             /* [0,1,...,nEx-1] notification dates:
                                           first coupon to be called is the first coupon with a start date
                                           on or after the notification date */
    Date* tStart,          /* [0,1,...,nEx-1] settlement dates for each exercise
                                                   i.e. dates on which the fee is paid */
    double* exerFee,       /* [0,1,...,nEx-1] fee paid by option holder to exercise
                                                   most frequently, the bond is called at par, i.e.
                              exerFee = 1.0 */
    const char* PayRecStr, /* RECEIVER (call on the bond) or PAYER (put on the bond) */

    /*	Coupons: the structure must be decomposed into a swap that delivers:
- on the funding side, Libor CASH FLAT
- on the exotic side, (cpn + gear*accrued_index)*accrued_ratio*coverage */
    double* cpn,                /* [0,1,...,nCpn-1]  coupons*/
    double* tgear_floatdigital, /*[0,1,...,nCpn-1] gearings for floating digitals */
    double* cvgCpn,             /* [0,1,...,nCpn-1]  coverages  */
    long    nCpn,               /* number of coupon periods */
    Date*   tCpnStart,          /* [0,1,...,nCpn-1] coupon start date */
    Date*   tCpnPay,            /* [0,1,...,nCpn-1] coupon pay date */
    double* tFundingPayment,    /* [0,1,...,nCpn-1] payment of the floating leg */

    // Spread:  not currently supported
    double* dvSpread,

    // floating digital
    SrtCompounding freqAccrued,
    double         CorrelStart,
    double         CorrelEnd,

    // Barriers details
    double   call_spread_up,
    double   call_spread_low,
    double** barriers,
    int      iAccrueOnBarrier,

    // Index details
    char*          szIndexRefRate,
    SrtCompounding freqIndex,

    //	Market details
    String ycName, /* yield curve name */
    String vcName, /* Vol curve name */

    // Stuff needed for the generic autocal code
    LGMErr (*GetVol)(long, long, double, SRT_Boolean, double*),
    /* Autocal volatility function to get market cash vol */
    char* char_vol_type, /* normal or log normal swaption vols */

    //	Calibration
    double tau,               /* if fixed tau, use this value for tau (in years) */
    int    calibrationmethod, /* 1 = Swaptions ATM, 2 = Swaptions Strike eq, 3 = Caplets ATM  */
    int    timeswapvolmethod, /* 1 = Model, 2 = Sliding, 3 = Converging  */
    int    atmvolmethod,      /* 1 = Lognormal, 2 = Normal, 3 = SigmaBeta */
    long   observation_freq,
    long   num_subdiscret_obsfreq,
    double shiftvol,
    long   nx,     /* number of state variable steps in the Pat's integral */
    double maxstd, /* Maximum number of std between forward and strike */
    int    nofunding,
    double dStrikeVolCap, /* Maximum lognormal std dev for smile around ATM (i.e., vol is flat
                             outside) */

    // Exercise Stuff
    int  isExercised,  // 0 not exercised; 1 exercised
    long lExDate,      // date when exercised
    int  eod_fix_flag, /*	EOD Fixing Flag 0: I, 1: E */
    int  eod_pay_flag, /*	EOD Payment Flag 0: I, 1: E */
    int  eod_ex_flag,  /*	EOD Exercise Flag 0: I, 1: E */

    // For period started in the past
    int     nDaysInFirstCoupon,
    long*   lvFirstCouponFixingDates,
    double* dvFirstCouponAccrualWeights,
    double* dvFirstCouponFixingOrSpread,
    double  dFloatCouponFixing,
    long    lPayDate,

    /* Accrual Period Day Rule -- not currently used:  1: only use business days, 3: count holidays
       using past setting */
    int iAccrualPeriodDayRule,

    /* Coupon Index Lag -- not currently used:  number of business days before start that each index
       sets */
    int iCouponIndexLag,

    /* End Lag -- not currently used:  number of business days before the end of the period that we
       stop checking the index */
    int iEndLag,

    /* calibrated term structure */
    SrtLgmTSData* lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr*    atcTSData, /* ptr to zeta/G data (NULL => not req'd) */

    /*	Results	*/
    double* Option,
    double* FixedLeg,
    //	double  *dFloatLeg,
    double** forwardTimeSwap, /*[0..1]*[0..nEx-1]  the first is the forward computed */
    double*  dvStrikes        /* the strikes used in the calibration. */
)
{
    /* Arguments needed in the generical Autocal code */
    int convertTS    = 1; /* 1 = cmpute new sigmas, taus & store in ts; 0 = don't bother */
    int findExerBdry = 1; /* 1 = find swap rates at exercise boundary; 0 = don't bother */
    SrtLgmExerBdryData* lgmExerBdryData = NULL;
    SrtLgmRefSwptnData* lgmRefSwptnData = NULL;
    FWD_VOL_STR*        fwdVolStr       = NULL;
    long                skipCal;
    Date                tNow, tfirst;
    SrtReceiverType     srt_p_r_flag;
    SrtCurvePtr         yldcrv = NULL;
    LGMErr (*GetBeta)(Date, Date, double*);
    SrtDiffusionType srt_vol_type;
    int              usefixtau, usecaps;
    LGMCalParm*      CalReqPtr = NULL;
    LGMErr (*RecVol)(Date, Date, double);
    ConvParams*      EvalParamsPtr = NULL;
    SrtCallTimeSwap* dealPtr       = NULL;
    SrtCallTimeSwap* dealTSPtr     = NULL;
    LGMDealType      dealtype;
    LGM_TS*          LGMtsPtr    = NULL;
    SigKapTS*        SigKaptsPtr = NULL;
    LGMSwptns*       ExerIntoPtr = NULL;
    LGMSwptns*       ExerBdryPtr = NULL;
    double           LGMVal;
    LGMErr           error;
    long             nExLeft;
    int              i, j, iEx, iLastStart;
    double           intrinsic_value;
    double           pv_fixedleg, dStubPV;

    // Initialize the output to 0
    *Option   = 0.0;
    *FixedLeg = 0.0;
    dStubPV   = 0.0;

    /* Set up the static market info */
    if (SetUpStaticMarketIfo(GetVol, ycName, dStrikeVolCap, char_vol_type))
        return "Error in setting up static market info in CallableTimeSwapCaller";
    tNow = StaticMarketInfo.CalcDate;

    // Check to see if the option is exercised, at which point we can return the current PV
    if (isExercised)
    {
        // Find the last coupon start date before the exercise date
        iLastStart = 0;
        while (iLastStart < nCpn && tCpnStart[iLastStart] < tNow)
            iLastStart++;
        // If the deal is exercised before the first start, then there is nothing to do
        if (iLastStart == 0)
        {
            return 0;
        }
        // Otherwise, calculate the PV of the fixed period that we are on
        return calcRangeAccrual(
            lPayDate,
            cpn[iLastStart - 1] + tgear_floatdigital[iLastStart - 1] * dFloatCouponFixing,
            eod_fix_flag,
            iAccrueOnBarrier,
            nDaysInFirstCoupon,
            lvFirstCouponFixingDates,
            dvFirstCouponAccrualWeights,
            dvFirstCouponFixingOrSpread,
            barriers[iLastStart - 1][1],
            barriers[iLastStart - 1][0],
            call_spread_up,
            call_spread_low,
            szIndexRefRate,
            cvgCpn[iLastStart - 1],
            FixedLeg);
    }

    // If Number of Days = 0, then skip the stub calculation
    if (nDaysInFirstCoupon > 0)
    {
        // Find the first start date in the future
        iLastStart = 0;
        while (iLastStart < nCpn && tCpnStart[iLastStart] < tNow)
            iLastStart++;
        // If there are no start dates in the future, then the PV is 0
        if (iLastStart == nCpn)
        {
            *Option   = 0.0;
            *FixedLeg = 0.0;
            return 0;
        }
        // Calculate the PV of the current period (0 if we haven't started yet)
        if (iLastStart == 0)
            *FixedLeg = 0.0;
        else if (
            error = calcRangeAccrual(
                lPayDate,
                cpn[iLastStart - 1] + tgear_floatdigital[iLastStart - 1] * dFloatCouponFixing,
                eod_fix_flag,
                iAccrueOnBarrier,
                nDaysInFirstCoupon,
                lvFirstCouponFixingDates,
                dvFirstCouponAccrualWeights,
                dvFirstCouponFixingOrSpread,
                barriers[iLastStart - 1][1],
                barriers[iLastStart - 1][0],
                call_spread_up,
                call_spread_low,
                szIndexRefRate,
                cvgCpn[iLastStart - 1],
                &dStubPV))
            return error;
        // Reset the dates and number of coupons to only include the future values
        cpn                = &cpn[iLastStart];
        tgear_floatdigital = &tgear_floatdigital[iLastStart];
        cvgCpn             = &cvgCpn[iLastStart];
        tCpnStart          = &tCpnStart[iLastStart];
        tCpnPay            = &tCpnPay[iLastStart];
        tFundingPayment    = &tFundingPayment[iLastStart];
        barriers           = &barriers[iLastStart];
        nCpn -= iLastStart;
        // Find the first future exercise date
        iEx = 0;
        while (iEx < nEx && tEx[iEx] < tNow)
            iEx++;
        // Move forward in time if we have already passed the exercise today
        if (tEx[iEx] == tNow && eod_ex_flag)
            iEx++;
        // Reset the dates and number of coupons to only include the future values
        tEx     = &tEx[iEx];
        tStart  = &tStart[iEx];
        exerFee = &exerFee[iEx];
        nEx -= iEx;
    }

    /*	1.- Initialise generic data
            ---------------------------	*/

    /* Free arrays that will be used for output, if already allocated */
    free_outputs(&convertTS, &findExerBdry, lgmExerBdryData, lgmRefSwptnData, lgmTSData, atcTSData);

    /* Set operations to be done by LGMnewautocal */
    skipCal = 0; /* always calibrate */

    /* Set information about today, and today's market */
    yldcrv = lookup_curve(ycName);
    //		tNow = get_clcndate_from_yldcrv (yldcrv);    /* get calculation date */

    /* Get the exponent beta to use to interpret market vols. Useful only to cope with the generica
     * autocal code */
    /* Transform char_vol_type to SrtDiffusionType */
    error = interp_diffusion_type(char_vol_type, &srt_vol_type);
    if (error)
        return error;

    if (srt_vol_type == SRT_LOGNORMAL)
    {
        GetBeta = LGMReturnExp1; /* This will return beta=1 for all swaptions */
    }
    else
    {
        GetBeta = LGMReturnExp0; /* This will return beta=0 for all swaptions */
    }

    /* Determine the first possible day that an exercise can take place */
    tfirst = add_unit(tNow, eod_ex_flag, SRT_BDAY, SUCCEEDING);

    /* Interpret pay-receive flag and basis */
    error = interp_rec_pay(PayRecStr, &srt_p_r_flag);
    if (error)
        return error;

    /* Always calibrate to swaptions  */
    usecaps = 0;
    if (tau == 0)
    {
        error = "No calibration to the full cap for the moment";
        return error;
    }
    else
    {
        usefixtau = 1; /* use fixed tau */
    }

    CalReqPtr =
        LGMSetCalibMeth(1, usefixtau, usecaps, tau, 0, 0, 0, maxstd, NULL, NULL, NULL, NULL, NULL);

    if (CalReqPtr == NULL)
    {
        return ("alloc failed in LGMcaller");
    }

    if (usefixtau != 1) /* calibrate tau to the full cap  */
    {
        CalReqPtr->calmeth = TenorAndDiag;
    }

    /* Set numerical parameters for convolution to their default values */
    EvalParamsPtr = LGMSetDefaultEvalParms();
    if (!EvalParamsPtr)
    {
        LGMFreeCalParm(&CalReqPtr);
        return "alloc failed in LGMcaller";
    }

    /* Record Vol */
    RecVol = LGMRecVolDummy;

    /*	2.- Initialise specific callable inverse floater data
            ----------------------------------------------------- */

    /* Set deal type */
    dealtype = CallTimeSwap;

    error = LGMFillCallableTimeSwap(
        &dealPtr,   /* output: the deal */
        &dealTSPtr, /* output: the underlying TS deal */
        &nExLeft,
        tfirst,
        ycName,
        nCpn,
        tCpnStart,
        tCpnPay,
        cpn,
        tgear_floatdigital,
        cvgCpn,
        tFundingPayment,
        nEx,
        tEx,
        tStart,
        exerFee,
        call_spread_up,
        call_spread_low,
        observation_freq,
        barriers,
        freqIndex,
        freqAccrued,
        CorrelStart,
        CorrelEnd,
        srt_p_r_flag,
        timeswapvolmethod,
        atmvolmethod,
        calibrationmethod,
        num_subdiscret_obsfreq,
        shiftvol,
        nx,
        nofunding,
        vcName,
        GetSABRParam,
        1,
        eod_pay_flag,
        eod_fix_flag);

    if (error)
    {
        LGMFreeCalParm(&CalReqPtr);
        srt_free(EvalParamsPtr);
        LGMFreeCallableTimeSwap(&dealPtr);
        return error;
    }

    if (nExLeft == 0)
    {
        *Option = 0.0;
        LGMFreeCalParm(&CalReqPtr);
        srt_free(EvalParamsPtr);
        LGMFreeCallableTimeSwap(&dealPtr);
        return NULL;
    }

    /*	3.- Call to LGMAutocal
            ---------------------- */

    /* allocate memory for the lgmRefSwpnData so that we can return the strikes */
    lgmRefSwptnData                   = srt_malloc(sizeof(SrtLgmRefSwptnData));
    lgmRefSwptnData->refCapSwptnArr   = 0;
    lgmRefSwptnData->refFixSwptnArr   = 0;
    lgmRefSwptnData->refLongSwptnArr  = 0;
    lgmRefSwptnData->refShortSwptnArr = 0;

    if (!(error = LGMnewautocal(
              0,
              skipCal,
              convertTS,
              findExerBdry, /* task list */
              tNow,
              eod_ex_flag,
              ycName, /* info about today's market */
              GetTruncVol,
              GetBeta,
              RecVol, /* swaption & caplet price info */
              dealtype,
              (void*)dealPtr, /* the deal */
              CalReqPtr,
              EvalParamsPtr, /* calib & eval methods */
              &LGMVal,
              &intrinsic_value, /* value & intrinsic value of deal */
              &ExerIntoPtr,
              &ExerBdryPtr,    /* swaptions representing the underlying & exer boundary */
              lgmRefSwptnData, /* reference swaptions */
              &LGMtsPtr,
              &SigKaptsPtr, /* calibrated term structures */
              NULL)))
    {
        /* Price of the underlying time swap */
        pv_fixedleg = 0;
        if (dealTSPtr)
        {
            error = CompVal_TimeSwap(dealTSPtr, tNow, ycName, &pv_fixedleg, 0, 0, 0);
        }
        pv_fixedleg += (dealPtr->pv_fixedleg) + dStubPV;
        if (dealPtr->PayRec == SRT_RECEIVER)
        {
            pv_fixedleg = -pv_fixedleg;
        }

        *FixedLeg = pv_fixedleg;

        /* to check pricing of forward time swap and the strikes used for the calibration*/
        for (j = 0; j < nEx; j++)
        {
            dvStrikes[j] = lgmRefSwptnData->refLongSwptnArr[j].strike;
            for (i = 0; i <= 1; i++)
                forwardTimeSwap[i][j] = (dealPtr->tForwardTS)[i][j];
        }
    }

    /*	4.- Free and return
            ------------------- */

    /* Free unneeded structures */ /* these are not needed to get output */
    LGMFreeCalParm(&CalReqPtr);    /* contains the calibration request */
    srt_free(EvalParamsPtr);       /* contains the convolution parameters */
    LGMFreeSwptns(&ExerIntoPtr);   /* structure containing underlying swaptions of MidAt */
    if (atcTSData)
    {
        *atcTSData = LGMtsPtr;
    }
    else
    {
        LGMFreeLGM_TS(&LGMtsPtr);
    }

    if (lgmRefSwptnData)
        free(lgmRefSwptnData);

    if (error)
    {
        LGMFreeCallableTimeSwap(&dealPtr);
        LGMFreeSwptns(&ExerBdryPtr);
        LGMFreeSigKapTS(&SigKaptsPtr);
        if (atcTSData)
        {
            *atcTSData = NULL;
            LGMFreeLGM_TS(&LGMtsPtr);
        }
        return error;
    }

    /* Unpack output */
    *Option = LGMVal;

    /* If exercise boundary has been found, copy it into tExBdry & rExBdry */
    if (findExerBdry != 0 && ExerBdryPtr != NULL && ExerBdryPtr->n > 0 && lgmExerBdryData != NULL)
    {
        error = copyExerBdry(ExerBdryPtr, lgmExerBdryData);
    }

    if (error)
    {
        LGMFreeCallableTimeSwap(&dealPtr);
        LGMFreeSwptns(&ExerBdryPtr);
        LGMFreeSigKapTS(&SigKaptsPtr);
        if (atcTSData)
        {
            *atcTSData = NULL;
            LGMFreeLGM_TS(&LGMtsPtr);
        }
        if (lgmExerBdryData)
        {
            lgmExerBdryData->NexerBdry = 0;
            if (lgmExerBdryData->exerBdryArr)
            {
                srt_free(lgmExerBdryData->exerBdryArr);
            }
        }
        return error;
    }
    LGMFreeSwptns(&ExerBdryPtr); /* no longer needed */

    /* If sigma-kappa term structure has been found, copy it into output arrays */
    if (convertTS != 0 && SigKaptsPtr != NULL && SigKaptsPtr->numS > 0 && SigKaptsPtr->numK > 0 &&
        SigKaptsPtr->sdate != NULL && SigKaptsPtr->sig != NULL && SigKaptsPtr->kdate != NULL &&
        SigKaptsPtr->kap != NULL && lgmTSData != NULL)
    {
        error = copyTSData(SigKaptsPtr, lgmTSData);
        if (error)
        {
            LGMFreeCallableTimeSwap(&dealPtr);
            LGMFreeSigKapTS(&SigKaptsPtr);
            if (atcTSData)
            {
                *atcTSData = NULL;
                LGMFreeLGM_TS(&LGMtsPtr);
            }
            return error;
        }
    }
    LGMFreeSigKapTS(&SigKaptsPtr); /* no longer needed */

    /* Copy forward vol */
    if (fwdVolStr)
    {
        fwdVolStr->fwdVolCpn = dealPtr->fwdVolCpn;
        fwdVolStr->fwdVolEx  = dealPtr->fwdVolEx;
        fwdVolStr->fwdVolMat = dmatrix(0, dealPtr->fwdVolCpn, 0, dealPtr->fwdVolEx);
        if (!fwdVolStr->fwdVolMat)
        {
            LGMFreeCallableTimeSwap(&dealPtr);
            return "Allocation of forward volatility matrix falied in LGMCallInvFloaterCaller";
        }

        for (i = 0; i <= dealPtr->fwdVolCpn; i++)
        {
            for (j = 0; j <= dealPtr->fwdVolEx; j++)
            {
                fwdVolStr->fwdVolMat[i][j] = dealPtr->fwdVolMat[i][j];
            }
        }
    }

    LGMFreeCallableTimeSwap(&dealPtr);
    LGMFreeCallableTimeSwap(&dealTSPtr);

    return error;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------
// //
//
// LGMNewCallableTimeSwapCaller
//
/* Create two structures, one for the option to call, and one for the time swap itself  */
LGMErr LGMFillCallableTimeSwap(
    SrtCallTimeSwap** dealPtrPtr,   /* output: the deal */
    SrtCallTimeSwap** dealTSPtrPtr, /* output: the TS deal */
    long*             nExleftPtr,   /* output: effective number of exer */
    Date              tfirst,
    char*             ycName, /* info about today */
    long              n,
    Date*             tStartCpn,
    Date*             tPayCpn, /* coupon leg */
    double*           tCpn,
    double*           tgear_floatdigital,
    double*           tCvgCpn,
    double*           tFundingPayment,
    long              nEx,
    Date*             tEx,
    Date*             tExStart,
    double*           exFee, /* exercise info */
    double            call_spread_up,
    double            call_spread_low,
    long              observation_freq,
    double**          barriers,
    SrtCompounding    freqIndex,
    SrtCompounding    freqAccrued,
    double            CorrelStart,
    double            CorrelEnd,
    SrtReceiverType   payrec, /* PAYER or RECEIVER */
    int               timeswapvolmethod,
    int               atmvolmethod,
    int               calibrationmethod,
    long              num_subdiscret_obsfreq,
    double            shiftvol,
    long              nx,
    int               nofunding,
    char*             vol_id,
    LGMErr (*GetSABRParam)(long, long, char*, SABR_VOL_TYPE, double*),
    double exFeeTimeSwap,
    int    eod_pay_flag, /*	0: I, 1: E */
    int    eod_fix_flag)    /*	0: I, 1: E */
{
    LGMErr           error = NULL;
    LGMMarkConv      conventions; /* standard swaption conventions */
    Date             tNow, tSpot;
    long             i, j, iFirst, iFirstTS;
    long             jEff, nExEff, iPrev;
    SrtCallTimeSwap* ptr    = NULL;
    SrtCallTimeSwap* ptrTS  = NULL;
    SrtCurvePtr      yldcrv = NULL;
    long             numdays_in_period;
    long             remaining_numdays_in_period;
    long             first_future_fixing;
    *dealPtrPtr = NULL;

    /* get yield curve */
    yldcrv = lookup_curve(ycName);
    tNow   = get_clcndate_from_yldcrv(yldcrv);     /* get calculation date */
    error  = LGMCcyDefaults(yldcrv, &conventions); /* get currency conventions */
    if (error != NULL)
        return (error);
    tSpot = add_unit(tNow, conventions.lag, SRT_BDAY, SUCCEEDING);

    /* Step 2: Determine the relevant exercise and coupon dates */
    /* Figure out the effective number of exercise dates */

    nExEff = 0;
    iPrev  = n;
    i      = n - 1;
    for (j = nEx - 1; j >= 0 && tEx[j] >= tfirst; j--)
    {
        while (i >= 0 && tStartCpn[i] >= tEx[j])
            i--;
        i++; /* i is first cpn date on or after tEx[j] */

        if (i < iPrev)
        {
            iPrev = i; /* exercise has value */
            nExEff++;
        }
    }
    iFirst = i; /* Index of the first coupon after first exercise */

    *nExleftPtr = nExEff; /* Number of exercises after today */
    if (nExEff == 0)
    {
        return NULL; /* no exercise dates left; deal is worth zero */
    }
    if (nExEff != nEx)
        return "Inconsistent number of exercises";

    /* Computes the number of periods IN THE TIME SWAP DEAL */
    i = n - 1;
    while (i >= 0 && tPayCpn[i] >= tNow + eod_pay_flag)
        i--;
    i++;
    if (i < n)
    {
        iFirstTS = i;
    }
    else
    {
        return NULL;
    }

    if (iFirstTS < iFirst) /* then we create a time swap deal */
    {
        ptrTS = LGMCreateCallableTimeSwap(1, iFirst - iFirstTS, observation_freq);
        switch (freqIndex)
        {
        case SRT_ANNUAL:
            ptrTS->undIndex = 12;
            break;
        case SRT_SEMIANNUAL:
            ptrTS->undIndex = 6;
            break;
        case SRT_QUARTERLY:
            ptrTS->undIndex = 3;
            break;
        case SRT_MONTHLY:
            ptrTS->undIndex = 1;
            break;
        }

        switch (freqAccrued)
        {
        case SRT_ANNUAL:
            ptrTS->undAccrued = 12;
            break;
        case SRT_SEMIANNUAL:
            ptrTS->undAccrued = 6;
            break;
        case SRT_QUARTERLY:
            ptrTS->undAccrued = 3;
            break;
        case SRT_MONTHLY:
            ptrTS->undAccrued = 1;
            break;
        }

        /* Fills the coupon, barriers information for the TS deal */
        for (i = iFirstTS; i < iFirst; i++)
        {
            ptrTS->tCpnStart[i - iFirstTS]       = tStartCpn[i];
            ptrTS->tCpnPay[i - iFirstTS]         = tPayCpn[i];
            ptrTS->tCpn[i - iFirstTS]            = tCpn[i];
            ptrTS->tFundingPayment[i - iFirstTS] = tFundingPayment[i];
            ptrTS->gear[i - iFirstTS]            = tgear_floatdigital[i];
            ptrTS->tCvgCpn[i - iFirstTS]         = tCvgCpn[i];
            ptrTS->barriers[i - iFirstTS][0]     = barriers[i][0];
            ptrTS->barriers[i - iFirstTS][1]     = barriers[i][1];

            /* Computes the observation days and related info */
            numdays_in_period = tPayCpn[i] - tStartCpn[i]; /* Assumes paydate = enddate */
            for (j = 1; j <= observation_freq; j++)
            {
                /* for each subperiod ... */
                ptrTS->observationdays[i - iFirstTS][j - 1][0] =
                    (int)(tStartCpn[i] + (j - 1) * (numdays_in_period + 0.0) / (observation_freq - 1.0)); /* Start Date */
                ptrTS->observationdays[i - iFirstTS][j - 1][1] = add_unit(
                    ptrTS->observationdays[i - iFirstTS][j - 1][0],
                    ptrTS->undIndex,
                    SRT_MONTH,
                    NO_BUSDAY_CONVENTION);
                ptrTS->observationdays[i - iFirstTS][j - 1][1] = add_unit(
                    ptrTS->observationdays[i - iFirstTS][j - 1][1],
                    0,
                    SRT_BDAY,
                    MODIFIED_SUCCEEDING);
            }
            if (observation_freq == 2) /* HARD CODE FOR THE PERIODIC CAP */
            {
                if (nx == 1)
                {
                    ptrTS->ratiodays_for_subperiod[i - iFirstTS][0] = 0;
                    ptrTS->ratiodays_for_subperiod[i - iFirstTS][1] = 1;
                }
                else
                {
                    ptrTS->ratiodays_for_subperiod[i - iFirstTS][0] = 1;
                    ptrTS->ratiodays_for_subperiod[i - iFirstTS][1] = 0;
                }
            }
            else
            {
                for (j = 1; j <= observation_freq; j++)
                {
                    if (j < observation_freq)
                    {
                        ptrTS->ratiodays_for_subperiod[i - iFirstTS][j - 1] =
                            (double)(ptrTS->observationdays[i - iFirstTS][j][0] - ptrTS->observationdays[i - iFirstTS][j - 1][0] + 0.0) /
                            (numdays_in_period + 0.0);
                    }
                    else
                    {
                        ptrTS->ratiodays_for_subperiod[i - iFirstTS][observation_freq - 1] =
                            (double)(tPayCpn[i] - (ptrTS->observationdays)[i - iFirstTS][observation_freq - 1][0] + 0.0) /
                            (numdays_in_period + 0.0);
                    }
                }
            }
        }

        ptrTS->PayRec                 = payrec;
        ptrTS->call_spread_up         = call_spread_up;
        ptrTS->call_spread_low        = call_spread_low;
        ptrTS->CorrelStart            = CorrelStart;
        ptrTS->CorrelEnd              = CorrelEnd;
        ptrTS->num_subdiscret_obsfreq = num_subdiscret_obsfreq;
        ptrTS->calibrationmethod      = calibrationmethod;
        ptrTS->timeswapvolmethod      = timeswapvolmethod;
        ptrTS->atmvolmethod           = atmvolmethod;
        ptrTS->num_subdiscret_obsfreq = num_subdiscret_obsfreq;
        ptrTS->shiftvol               = shiftvol;
        ptrTS->nofunding              = nofunding;
        ptrTS->vol_id                 = vol_id;
        ptrTS->GetSABRParam           = GetSABRParam;
    }

    /* Step 3: Allocate Callable Time Swap with nExEff exer dates and n-i cpn periods */
    ptr = LGMCreateCallableTimeSwap(nExEff, n - iFirst, observation_freq);
    /*	ptrTS = LGMCreateCallableTimeSwap(1,n-iFirstTS,observation_freq); */

    /* Step 4: Fill the Callable Time Swap structure with the input data */
    /* Fills und of index and accrued */
    switch (freqIndex)
    {
    case SRT_ANNUAL:
        ptr->undIndex = 12;
        break;
    case SRT_SEMIANNUAL:
        ptr->undIndex = 6;
        break;
    case SRT_QUARTERLY:
        ptr->undIndex = 3;
        break;
    case SRT_MONTHLY:
        ptr->undIndex = 1;
        break;
    }

    /*		switch (freqIndex) {
            case SRT_ANNUAL:	ptrTS->undIndex = 12;
                                                    break;
            case SRT_SEMIANNUAL:	ptrTS->undIndex = 6;
                                                            break;
            case SRT_QUARTERLY:		ptrTS->undIndex = 3;
                                                            break;
            case SRT_MONTHLY:		ptrTS->undIndex = 1;
                                                            break;
            }
            */

    switch (freqAccrued)
    {
    case SRT_ANNUAL:
        ptr->undAccrued = 12;
        break;
    case SRT_SEMIANNUAL:
        ptr->undAccrued = 6;
        break;
    case SRT_QUARTERLY:
        ptr->undAccrued = 3;
        break;
    case SRT_MONTHLY:
        ptr->undAccrued = 1;
        break;
    }
    /*
                    switch (freqAccrued) {
            case SRT_ANNUAL:	ptrTS->undAccrued = 12;
                                                    break;
            case SRT_SEMIANNUAL:	ptrTS->undAccrued = 6;
                                                            break;
            case SRT_QUARTERLY:		ptrTS->undAccrued = 3;
                                                            break;
            case SRT_MONTHLY:		ptrTS->undAccrued = 1;
                                                            break;
                    }
                    */

    /* Fills the coupon, barriers information for the option deal */
    for (i = iFirst; i < n; i++)
    {
        ptr->tCpnStart[i - iFirst]       = tStartCpn[i];
        ptr->tCpnPay[i - iFirst]         = tPayCpn[i];
        ptr->tCpn[i - iFirst]            = tCpn[i];
        ptr->tFundingPayment[i - iFirst] = tFundingPayment[i];
        ptr->gear[i - iFirst]            = tgear_floatdigital[i];
        ptr->tCvgCpn[i - iFirst]         = tCvgCpn[i];
        ptr->barriers[i - iFirst][0]     = barriers[i][0];
        ptr->barriers[i - iFirst][1]     = barriers[i][1];

        /* Computes the observation days and related info */
        numdays_in_period = tPayCpn[i] - tStartCpn[i]; /* Assumes paydate = enddate */
        if (tStartCpn[i] < tNow + eod_fix_flag)        /* Fixing in the past ? */
            first_future_fixing = tNow + eod_fix_flag;
        else
            first_future_fixing = tStartCpn[i];
        remaining_numdays_in_period = tPayCpn[i] - first_future_fixing;
        for (j = 1; j <= observation_freq; j++)
        {
            /* for each subperiod ... */
            ptr->observationdays[i - iFirst][j - 1][0] =
                (int)(first_future_fixing + (j - 1) * (remaining_numdays_in_period + 0.0) / (observation_freq + 0.0)); /* Start Date */
            ptr->observationdays[i - iFirst][j - 1][1] = add_unit(
                ptr->observationdays[i - iFirst][j - 1][0],
                ptr->undIndex,
                SRT_MONTH,
                NO_BUSDAY_CONVENTION);
            ptr->observationdays[i - iFirst][j - 1][1] = add_unit(
                ptr->observationdays[i - iFirst][j - 1][1], 0, SRT_BDAY, MODIFIED_SUCCEEDING);
        }
        if (observation_freq == 2) /* HARD CODE FOR THE PERIODIC CAP */
        {
            if (nx == 1) /* HARD CODE IN ARREARS */
            {
                ptr->ratiodays_for_subperiod[i - iFirst][0] = 0;
                ptr->ratiodays_for_subperiod[i - iFirst][1] = 1;
            }
            else
            {
                ptr->ratiodays_for_subperiod[i - iFirst][0] = 1;
                ptr->ratiodays_for_subperiod[i - iFirst][1] = 0;
            }
        }
        else
        {
            for (j = 1; j <= observation_freq; j++)
            {
                if (j < observation_freq)
                {
                    ptr->ratiodays_for_subperiod[i - iFirst][j - 1] =
                        (double)(ptr->observationdays[i - iFirst][j][0] - ptr->observationdays[i - iFirst][j - 1][0] + 0.0) /
                        (numdays_in_period + 0.0);
                }
                else
                {
                    ptr->ratiodays_for_subperiod[i - iFirst][observation_freq - 1] =
                        (double)(tPayCpn[i] - (ptr->observationdays)[i - iFirst][observation_freq - 1][0] + 0.0) /
                        (numdays_in_period + 0.0);
                }
            }
        }
    }

    /* Fills the exercise information */
    iPrev = n;
    i     = n - 1;
    jEff  = nExEff;

    for (j = nEx - 1; j >= 0 && tEx[j] >= tfirst; j--)
    {
        while (i >= 0 && tStartCpn[i] >= tEx[j])
            i--;
        i++; /* i is first cpn date on or after tEx[j] */
        if (i < iPrev)
        {
            iPrev = i;
            jEff--;
            ptr->tEx[jEff]    = tEx[j];
            ptr->tSet[jEff]   = tExStart[j];
            ptr->iSet[jEff]   = i - iFirst;
            ptr->strike[jEff] = exFee[j];
        }
    }

    /*
    ptrTS->tEx[0] = tNow;
    ptrTS->tSet[0] = tStartCpn[iFirstTS];
    ptrTS->iSet[0] = 0;
    ptrTS->strike[0] = exFeeTimeSwap;
    */

    if (jEff != 0)
    {
        LGMFreeCallableTimeSwap(&ptr);
        LGMFreeCallableTimeSwap(&ptrTS);
        return "Serious error";
    }

    ptr->PayRec = payrec;
    /*	ptrTS->PayRec = payrec;	*/

    /* Fills  and index details */
    ptr->call_spread_up  = call_spread_up;
    ptr->call_spread_low = call_spread_low;
    ptr->CorrelStart     = CorrelStart;
    ptr->CorrelEnd       = CorrelEnd;

    /*
            ptrTS->call_spread_up = call_spread_up;
            ptrTS->call_spread_low = call_spread_low;
            ptrTS->CorrelStart = CorrelStart;
            ptrTS->CorrelEnd = CorrelEnd;
            */

    /* Fills the time swap forward pricing info data */
    ptr->calibrationmethod      = calibrationmethod;
    ptr->timeswapvolmethod      = timeswapvolmethod;
    ptr->atmvolmethod           = atmvolmethod;
    ptr->num_subdiscret_obsfreq = num_subdiscret_obsfreq;
    ptr->shiftvol               = shiftvol;

    ptr->nx = nx;

    /*	ptr->nx = 192;	*/
    /* HARD CODE OR THE NUMBER OF INTEGRATION POINTS */

    ptr->nofunding    = nofunding;
    ptr->vol_id       = vol_id;
    ptr->GetSABRParam = GetSABRParam;

    /*
    ptrTS->num_subdiscret_obsfreq = num_subdiscret_obsfreq;
    ptrTS->vol_id = vol_id;
    ptrTS->GetSABRParam = GetSABRParam;
    ptrTS->timeswapvolmethod = timeswapvolmethod;
    ptrTS->atmvolmethod = atmvolmethod;
    */

    *dealPtrPtr   = ptr;
    *dealTSPtrPtr = ptrTS;
    return NULL;
}

/********************/
/* copying routines */
/********************/
static LGMErr copyExerBdry(LGMSwptns* ExerBdryPtr, SrtLgmExerBdryData* lgmExerBdryData)
{
    long i;

    if (lgmExerBdryData == NULL)
        return ("no lgmExerBdryData structure");

    lgmExerBdryData->NexerBdry = ExerBdryPtr->n;

    if (lgmExerBdryData->exerBdryArr != NULL)
        srt_free(lgmExerBdryData->exerBdryArr);
    lgmExerBdryData->exerBdryArr =
        (SrtLgmExerBdry*)srt_calloc(lgmExerBdryData->NexerBdry, sizeof(SrtLgmExerBdry));

    if (lgmExerBdryData->exerBdryArr == NULL)
        return ("alloc failed in LGM, exerbdry");

    for (i = 0; i < lgmExerBdryData->NexerBdry; i++)
    {
        lgmExerBdryData->exerBdryArr[i].parRate = ExerBdryPtr->Rfix[i];
        lgmExerBdryData->exerBdryArr[i].bgnDate = ExerBdryPtr->tEx[i];
    }
    return (NULL);
}

static LGMErr copyTSData(SigKapTS* SKtsPtr, SrtLgmTSData* TSData)
{
    long j;

    if (TSData == NULL)
        return (" no TSData structure");
    TSData->NTS   = max(SKtsPtr->numK, SKtsPtr->numS);
    TSData->TSArr = (SrtLgmTS*)srt_calloc(TSData->NTS, sizeof(SrtLgmTS));
    if (TSData->TSArr == NULL)
        return ("alloc failed in LGM, termstruct");

    for (j = 0; j < TSData->NTS; j++)
    {
        if (j < SKtsPtr->numS)
        {
            TSData->TSArr[j].sigma   = SKtsPtr->sig[j];
            TSData->TSArr[j].sigmaDt = SKtsPtr->sdate[j];
        }
        else
        {
            TSData->TSArr[j].sigma   = 0;
            TSData->TSArr[j].sigmaDt = 0;
        }

        if (j < SKtsPtr->numK)
        {
            TSData->TSArr[j].tauDt = SKtsPtr->kdate[j];
            TSData->TSArr[j].tau   = SKtsPtr->kap[j];
            if (fabs(TSData->TSArr[j].tau) < 0.000001)
                TSData->TSArr[j].tau = 1000000.;
            else
                TSData->TSArr[j].tau = 1 / TSData->TSArr[j].tau;
        }
        else
        {
            TSData->TSArr[j].tau   = 0;
            TSData->TSArr[j].tauDt = 0;
        }
    }
    return (NULL);
}

/* Free arrays that will be used for output, if already allocated */
static void free_outputs(
    int* convertTS,    /* 1=compute new sigs and taus; 0=don't bother */
    int* findExerBdry, /* 1=find swap rates at exercise boundary; 0=don't bother */
    SrtLgmExerBdryData*
        lgmExerBdryData, /* ptr to exercise boundary data structure (NULL => not req'd) */
    SrtLgmRefSwptnData*
        lgmRefSwptnData,     /* ptr to reference swaption data structure (NULL => not req'd) */
    SrtLgmTSData* lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr*    atcTSData)    /* ptr to zeta/G data (NULL => not req'd) */
{
    if (*findExerBdry != 0 && lgmExerBdryData != NULL)
    {
        lgmExerBdryData->NexerBdry = 0;
        if (lgmExerBdryData->exerBdryArr != NULL)
        {
            srt_free(lgmExerBdryData->exerBdryArr);
        }
    }
    else
    {
        *findExerBdry = 0;
    }

    if (*convertTS != 0 && lgmTSData != NULL)
    {
        lgmTSData->NTS = 0;
        if (lgmTSData->TSArr != NULL)
        {
            srt_free(lgmTSData->TSArr);
        }
    }
    else
        *convertTS = 0;

    if (*convertTS != 0 && atcTSData != NULL && *atcTSData != NULL)
    {
        if ((*atcTSData)->zdate)
            srt_free((*atcTSData)->zdate);
        if ((*atcTSData)->zeta)
            srt_free((*atcTSData)->zeta);
        if ((*atcTSData)->Gdate)
            srt_free((*atcTSData)->Gdate);
        if ((*atcTSData)->G)
            srt_free((*atcTSData)->G);
        (*atcTSData)->zdate = NULL;
        (*atcTSData)->zeta  = NULL;
        (*atcTSData)->Gdate = NULL;
        (*atcTSData)->G     = NULL;
    }
}

static int IsAccrued(int iAccrueOnBarrier, double dFixing, double dUpBarrier, double dLowerBarrier)
{
    if (iAccrueOnBarrier)
    {
        return (dFixing <= dUpBarrier) && (dFixing >= dLowerBarrier);
    }
    else
    {
        return (dFixing < dUpBarrier) && (dFixing > dLowerBarrier);
    }
}

static char* BS(
    long           lFixing,
    double         dStrike,
    char*          szRefRate,
    double         dSpread,
    SrtCallPutType srtCallPut,
    double*        out_pdPrice)
{
    // Variable declarations
    double         dVol, dFwd, dExpiry;
    char*          err;
    long           lStart, lEnd;
    SrtCompounding Comp;
    SrtBasisCode   Basis;

    // Get the expiry
    dExpiry = (lFixing - StaticMarketInfo.CalcDate) / 365.0;

    // If the expiry is less than 0, the value is 0
    if (dExpiry < 0)
    {
        *out_pdPrice = 0.0;
        return 0;
    }

    // Get the refrate details
    swp_f_get_ref_rate_details(szRefRate, &Basis, &Comp);

    // Get the start and end dates
    lStart = add_unit(
        lFixing, StaticMarketInfo.Conventions.lag, SRT_BDAY, StaticMarketInfo.Conventions.cbdconv);
    lEnd = add_unit(lStart, (int)(12 / Comp), SRT_MONTH, StaticMarketInfo.Conventions.cbdconv);

    // Get the forward
    dFwd =
        swp_f_swapcashrate(lStart, lEnd, Basis, Comp, StaticMarketInfo.YieldCurveName, szRefRate);

    // If the expiry is 0, return the intrinsic value
    if (dExpiry == 0.0)
    {
        *out_pdPrice =
            srt_f_optblksch(dFwd, dStrike - dSpread, 0.0, dExpiry, 1.0, srtCallPut, PREMIUM);
        return 0;
    }

    // Get the volatility
    if (err = StaticMarketInfo.GetVol(lStart, lEnd, dStrike, 0, &dVol))
        return err;

    // Price the option
    if (StaticMarketInfo.IsLognormal)
        *out_pdPrice =
            srt_f_optblksch(dFwd, dStrike - dSpread, dVol, dExpiry, 1.0, srtCallPut, PREMIUM);
    else
        *out_pdPrice =
            srt_f_optblknrm(dFwd, dStrike - dSpread, dVol, dExpiry, 1.0, srtCallPut, PREMIUM);

    // return
    return 0;
}

static char* putSpread(
    long    lFixing,
    double  dBarrier,
    double  dSpread,
    char*   szRefRate,
    double  dPutSpread,
    double* out_pdProb)
{
    // Variable declaration
    double dUpPV, dDownPV;
    char*  err;

    // Check that the barrier is not zero
    if (dBarrier == 0.0)
    {
        *out_pdProb = 0.0;
        return 0;
    }

    // Check that the put spread is non-zero
    if (dPutSpread == 0.0)
        return "Cannot have a negative put spread";

    // Calculate the put spread
    if (err = BS(lFixing, dBarrier + dPutSpread, szRefRate, dSpread, SRT_PUT, &dUpPV))
        return err;
    if (err = BS(lFixing, dBarrier, szRefRate, dSpread, SRT_PUT, &dDownPV))
        return err;

    // Return
    *out_pdProb = (dUpPV - dDownPV) / dPutSpread;
    return 0;
}

static char* getSinglePV(
    long    lFixing,
    double  dUpperBarrier,
    double  dLowerBarrier,
    double  dUpperSpread,
    double  dLowerSpread,
    char*   szRefRate,
    double  dSpread,
    double* out_pdProb)
{
    // Variable declarations
    double dUpperProb, dLowerProb;
    char*  err;

    // Calculate the put spread prices
    if (err = putSpread(lFixing, dUpperBarrier, dSpread, szRefRate, dUpperSpread, &dUpperProb))
        return err;

    if (err = putSpread(lFixing, dLowerBarrier, dSpread, szRefRate, dLowerSpread, &dLowerProb))
        return err;

    // return
    *out_pdProb = dUpperProb - dLowerProb;
    return 0;
}

LGMErr calcRangeAccrual(
    long    lPayDate,
    double  dCoupon,
    int     iEODFixFlag,
    int     iAccrueOnBarrier,
    int     nDaysInCoupon,
    long*   lvFixingDates,
    double* dvWeights,
    double* dvFixingOrSpread,
    double  dUpperBarrier,
    double  dLowerBarrier,
    double  dUpperSpread,
    double  dLowerSpread,
    char*   szRefRate,
    double  dCvg,
    double* out_pdFixed)
{
    // Variable declaration
    int    i;
    double dProb, dAccrual, dTotalDays;
    long   lToday = StaticMarketInfo.CalcDate;
    char*  err;

    // Loop until today, counting the number of accrual dates
    i          = 0;
    dAccrual   = 0.0;
    dTotalDays = 0.0;
    while (i < nDaysInCoupon && lvFixingDates[i] < lToday)
    {
        if (IsAccrued(iAccrueOnBarrier, dvFixingOrSpread[i], dUpperBarrier, dLowerBarrier))
            dAccrual += dvWeights[i];
        dTotalDays += dvWeights[i];
        i++;
    }

    // Test today, depending on the fixing flag
    if (i < nDaysInCoupon && iEODFixFlag)
    {
        if (IsAccrued(iAccrueOnBarrier, dvFixingOrSpread[i], dUpperBarrier, dLowerBarrier))
            dAccrual += dvWeights[i];
        dTotalDays += dvWeights[i];
        i++;
    }

    // Loop over the remaining dates, calculating the PV of each
    while (i < nDaysInCoupon)
    {
        if (err = getSinglePV(
                lvFixingDates[i],
                dUpperBarrier,
                dLowerBarrier,
                dUpperSpread,
                dLowerSpread,
                szRefRate,
                dvFixingOrSpread[i],
                &dProb))
            return err;
        dAccrual += dvWeights[i] * dProb;
        dTotalDays += dvWeights[i];
        i++;
    }

    // return the value
    if (dTotalDays == 0)
        return ("Total number of days is zero!");

    *out_pdFixed = dCvg * dAccrual * dCoupon *
                   swp_f_df(lToday, lPayDate, StaticMarketInfo.YieldCurveName) / dTotalDays;

    return 0;
}

/* ------------------------------------------------------------------------------------------------------------------
 */
/*
        New midat interface with fixing date in vol function
*/
/* ------------------------------------------------------------------------------------------------------------------
 */

/* Static pointer to 3 date vol function */
static LGMErr (*STATIC_ptrFixingGetVol)(Date, Date, Date, double, SRT_Boolean, double*);
static int           STATIC_lag;
static SrtBusDayConv STATIC_BusDayConv;

/* Conversion routine between 2 date and 3 date vol function */
static LGMErr STATIC_GetVol(
    Date lStartDate, Date lEndDate, double dStrike, SRT_Boolean b, double* Vol)
{
    /* Calculate the fixing date */
    Date lFixingDate;
    lFixingDate = add_unit(lStartDate, -STATIC_lag, SRT_BDAY, STATIC_BusDayConv);
    return (*STATIC_ptrFixingGetVol)(lFixingDate, lStartDate, lEndDate, dStrike, b, Vol);
}

LGMErr TestLGMautocalCallerFix(
    long    nEx,         /* nEx is number of exercises */
    Date*   tEx,         /* notification (exercise) dates, [0,1,...,nEx-1] */
    Date*   tStart,      /* start dates for each exercise, [0,1,...,nEx-1] */
    double* Strike,      /* total value paid at tStart[j] for fixed leg, [0,1,...,nEx-1] */
    long    nPay,        /* nPay is number of fixed leg coupons */
    Date*   tPay,        /* pay dates for period i, [0,1, ...,nPay-1] */
    double* Payment,     /* total fixed leg payment(last includes notional), [0,...,nPay-1] */
    double* RedFirstPay, /* reduction in 1rst payment after exercise, [0,...,nEx-1] */
    char*   PayRecStr,   /* RECEIVER or PAYER */
                         /* information about today */
    String ycName,       /* pointer to market structures */
    int    endofday,     /* 1=too late to exercise deals today, 0=not too late */
                         /* today's volatilities */
    LGMErr (*GetVol)(
        Date, Date, Date, double, SRT_Boolean, double*), /* function to get swaption vols */
    char* char_vol_type, /* determines whether vol is normal or log normal */
                         /* calibration method to use */
    int    LGMOneTwoFactor,
    int    usefixtau, /* 1=calibrate with fixed tau */
    int    usecaps,   /* 1=use caplets for calibration, 0=use only swaptions */
    double tau,       /* if fixed tau, use this value for tau (in years) */
    double alpha,
    double gamma,
    double rho,
    int    calibrationmeth, /* 1 = fixexp, 2=backboot, 3=fixed sigma */
    int    strikechoice,    /* 1 = IRR, 2 = dIRR */
    double maxstd,          /* Maximum number of std between forward and strike */
                            /* requested operation */
    int      skipEval,      /* 1=calibrate only, 0=calibrate & evaluate deal */
    int      convertTS,     /* 1=compute new sigs and taus; 0=don't bother */
    int      findExerBdry,  /* 1=find swap rates at exercise boundary; 0=don't bother */
    long*    Zeta1Dates,
    double*  StartZeta1s,
    long*    TauDates,
    double*  StartTaus,
    double** HybridShortInstrsIndex,

    /* outputs */
    String  outfile,   /* output file name for log file (unused) */
    double* LGMvalPtr, /* LGM value of mid-atlantic */
    SrtLgmExerBdryData*
        lgmExerBdryData, /* ptr to exercise boundary data structure (NULL => not req'd) */
    SrtLgmRefSwptnData*
        lgmRefSwptnData,     /* ptr to reference swaption data structure (NULL => not req'd) */
                             /* calibrated term structure */
    SrtLgmTSData* lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr*    atcTSData, /* ptr to zeta/G data (NULL => not req'd) */
                             /* miscellaneous */
    double* intrinsicValPtr, /* ptr to intrinsic value (NULL => not req'd) */
    char*   szRefRate)
{
    /* Variable declaration */
    char*       err;
    SrtCurvePtr yldcrv = NULL;
    LGMMarkConv conventions; /* standard swaption conventions */

    /* set up the static vol function ptr */
    STATIC_ptrFixingGetVol = GetVol;
    yldcrv                 = lookup_curve(ycName);

    /* Find the currency defaults and set up the static data */
    if (err = LGMCcyDefaults(yldcrv, &conventions))
        return err;

    STATIC_lag        = conventions.lag;
    STATIC_BusDayConv = conventions.sbdconv;

    // Call the regular autocal
    return TestLGMautocalCallerRR(
        nEx,           /* nEx is number of exercises */
        tEx,           /* notification (exercise) dates, [0,1,...,nEx-1] */
        tStart,        /* start dates for each exercise, [0,1,...,nEx-1] */
        Strike,        /* total value paid at tStart[j] for fixed leg, [0,1,...,nEx-1] */
        nPay,          /* nPay is number of fixed leg coupons */
        tPay,          /* pay dates for period i, [0,1, ...,nPay-1] */
        Payment,       /* total fixed leg payment(last includes notional), [0,...,nPay-1] */
        RedFirstPay,   /* reduction in 1rst payment after exercise, [0,...,nEx-1] */
        PayRecStr,     /* RECEIVER or PAYER */
        ycName,        /* pointer to market structures */
        endofday,      /* 1=too late to exercise deals today, 0=not too late */
        STATIC_GetVol, /* function to get swaption vols */
        char_vol_type, /* determines whether vol is normal or log normal */
        LGMOneTwoFactor,
        usefixtau, /* 1=calibrate with fixed tau */
        usecaps,   /* 1=use caplets for calibration, 0=use only swaptions */
        tau,       /* if fixed tau, use this value for tau (in years) */
        alpha,
        gamma,
        rho,
        calibrationmeth, /* 1 = fixexp, 2=backboot, 3=fixed sigma */
        strikechoice,    /* 1 = IRR, 2 = dIRR */
        maxstd,          /* Maximum number of std between forward and strike */
        skipEval,        /* 1=calibrate only, 0=calibrate & evaluate deal */
        convertTS,       /* 1=compute new sigs and taus; 0=don't bother */
        findExerBdry,    /* 1=find swap rates at exercise boundary; 0=don't bother */
        Zeta1Dates,
        StartZeta1s,
        TauDates,
        StartTaus,
        HybridShortInstrsIndex,
        outfile,         /* output file name for log file (unused) */
        LGMvalPtr,       /* LGM value of mid-atlantic */
        lgmExerBdryData, /* ptr to exercise boundary data structure (NULL => not req'd) */
        lgmRefSwptnData, /* ptr to reference swaption data structure (NULL => not req'd) */
        lgmTSData,       /* ptr to tau/sigma data (NULL => not req'd) */
        atcTSData,       /* ptr to zeta/G data (NULL => not req'd) */
        intrinsicValPtr, /* ptr to intrinsic value (NULL => not req'd) */
        szRefRate);
}
