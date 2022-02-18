#include "opfnctns.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lgmUStypes.h"
#include "srt_h_lgmtypes.h"
#include "swp_h_spread.h"
#include "swp_h_vol.h"

#define INTRINSICDEFAULT -333333
#include "srt_h_all.h"
#include "srt_h_lgmUSprotos.h"
#include "srt_h_lgmprotos.h"
#include "swp_h_all.h"
#include "swp_h_swap_pricing.h"

#define MAX_CPN 600

static LGMErr FitAllZetaAmortisingSwaptions(
    LGM_TS*    tsPtr,
    LGMCalSet* CSPtr,
    double     NumPeriod,
    double*    index,
    double*    d,
    double*    dstart,
    double*    mktprice);

LGMErr LGMCalFixKapamortising(
    LGM_TS**   LGMtsPtrPtr, /* Return: Calibrated term structure */
    double     kap,         /* Fixed kappa to use to construct G */
    int        usestarts,   /* construct G on start dates or pay dates */
    long       nDate,       /* number of dates for G if usestarts=0 */
    Date*      GDate,       /* dates for G if usestarts=0 */
    int        usecaps,     /* calibrate on caplets or swaptions */
    LGMCalSet* CSPtr,
    double     NumPeriod,
    double*    index,
    double*    d,
    double*    dstart,
    double*    mktprice); /* Reference caplets */

static double RecSqZetaDer(
    double  ystar,
    long    n,
    double* a,
    double  DStart, /* receiver swaption */
    double  sqzeta,
    double* Gpay,
    double  Gst); /* term structure */

double RecGDer(
    double  ystar,
    long    n,
    double* a,
    double  DStart, /* receiver swaption */
    double  sqzeta,
    double* Gpay,
    double  Gst, /* term structure */
    double  dGst,
    double* dGpay);

Err lgm_2f_autocal_calibrate(
    void*               dealPtr,
    LGM_TS**            lgm_ts_ptr_ptr,
    SrtLgmRefSwptnData* lgm_ref_swp_data,
    LGMCalParm*         cal_req,
    Err (*srt_f_get_vol)(Date, Date, double, SRT_Boolean, double*),
    Err (*srt_f_get_beta)(Date, Date, double*),
    String yc_name);

LGMErr LGMamortisingCalibration(
    /* info about today's market */
    Date   tNow,   /* first possible exercise date (today+eod) */
    String ycname, /* yield curve name for discount factors */
    LGMErr (*GetVol)(
        Date, Date, double, SRT_Boolean, double*), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date, double*), /* function to provide exponent (beta) to interpret vols */
    LGMErr (*RecVol)(Date, Date, double), /* function called whenever LGM uses a vol for pricing */
    double  NumPeriod,
    double* index,
    double* d,
    double* dstart,
    double* mktprice,
    /* info about the deal */
    LGMDealType dealType,
    void*       dealptr,
    long        nEx,    /* number of exercise dates to calibrate to */
    Date        tlast,  /* last date that needs to be calibrated */
    Date*       TauArr, /* calibrate on these exercise dates */
    double*     FVArr,  /* use swaptions with these forward vals of fixed legs */
    /* calibration methodology */
    LGMCalParm* CalReq,        /* structure defining calibration method to be used */
                               /* calibrated term structure */
    LGM_TS**    LGMtsPtrPtr,   /* calibrated zeta-G term structure */
    LGMCalSet** RefDealPtrPtr, /* Reference instruments used for calibration */
    SrtLgmRefSwptnData*
        lgmRefSwptnData); /* ptr to reference swaption data structure (NULL => not req'd) */

static double FitOneCap(
    double price, double cvg, double Rfix, double DStart, double Dend, SrtReceiverType recpay);

static LGMErr FitAllZetaCaps(
    LGM_TS*    tsPtr,  /* Return: zeta's in LGM term structure ts */
    LGMCalSet* CSPtr); /* Reference caplets */

static LGM_TSPtr CreateTSwithKappa(
    double     kap,       /* Fixed kappa to use to construct G */
    int        usestarts, /* construct G on start dates or pay dates */
    long       ndate,     /* number of dates for G if usestarts=0 */
    Date*      Gdate,     /* dates for G if usestarts=0 */
    int        usecaps,   /* final date from caplet or swaption */
    LGMCalSet* CSPtr);    /* Reference instruments */

LGMErr LGMExtractInfoFromDeal(
    /* information about today's market */
    Date   tNow,   /* evaluation date */
    Date   tfirst, /* first possible exercise date (eval date + eod) */
    String ycName, /* yield curve for discount factors */
    LGMErr (*GetVol)(
        Date, Date, double, SRT_Boolean, double*), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date, double*), /* function to provide exponent (beta) to interpret vols */
                              /* information about the deal */
    LGMDealType dealtype,     /* type of deal */
    void*       dealPtr,      /* pointer to the deal */
                              /* output */
    long*       EffnEx,       /* Effective number of exercise dates */
    Date*       LasttPay,     /* Last pay date of deal */
    Date**      TauArrPtr,    /* Array of effective exercise dates */
    double**    FVArrPtr,     /* PV of fixed leg/PV of strike for these dates */
    double*     intValPtr,    /* intrinsic value of the deal */
    LGMSwptns** ExerIntoPtr); /* closest European swaptions underlying deal */

LGMErr TestLGMamortisingautocalCaller(
    long    nEx,         /* nEx is number of exercises */
    Date*   tEx,         /* notification (exercise) dates, [0,1,...,nEx-1] */
    Date*   tStart,      /* start dates for each exercise, [0,1,...,nEx-1] */
    double* Strike,      /* total value paid at tStart[j] for fixed leg, [0,1,...,nEx-1] */
    long    nPay,        /* nPay is number of fixed leg coupons */
    Date*   tPay,        /* pay dates for period i, [0,1, ...,nPay-1] */
    double* Payment,     /* total fixed leg payment(last includes notional), [0,...,nPay-1] */
    double* RedFirstPay, /* reduction in 1rst payment after exercise, [0,...,nEx-1] */
    char*   PayRecStr,   /* RECEIVER or PAYER */
    double  NumPeriod,
    double* index,
    double* d,
    double* dstart,
    double* mktprice,
    /* information about today */
    String ycName,   /* pointer to market structures */
    int    endofday, /* 1=too late to exercise deals today, 0=not too late */
                     /* today's volatilities */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* function to get swaption vols */
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
        lgmRefSwptnData,      /* ptr to reference swaption data structure (NULL => not req'd) */
                              /* calibrated term structure */
    SrtLgmTSData* lgmTSData,  /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr*    atcTSData,  /* ptr to zeta/G data (NULL => not req'd) */
                              /* miscellaneous */
    double* intrinsicValPtr); /* ptr to intrinsic value (NULL => not req'd) */

/* Fits value of zeta at exercise date for a swaption */
static double FitOneZetaSwap(
    long            n,
    double*         a,
    double          DStart,
    SrtReceiverType recpay,
    double*         Gpay,
    double          Gst,
    double          price);

LGMErr LGMamortisingnewautocal(
    /* task list */
    int skipEval,     /* 0=calibrate & value deal, 1=calibrate only */
    int skipCalib,    /* 0=calibrate & value deal, 1=value deal only */
    int convertTS,    /* 1=compute & output sig-kappa ts equivalent to LGM ts */
    int findExerBdry, /* 1=find swap rates at exer boundary  */
                      /* information about today & today's market prices */
    Date   tNow,      /* eval as if today is tNow */
    int    eod,       /* 0=can exercise on tNow, 1=cannot exercise on tNow */
    String ycName,    /* market pointer for discount factors */
    LGMErr (*GetVol)(
        Date, Date, double, SRT_Boolean, double*), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date, double*), /* function to provide exponent (beta) to interpret vols */
    LGMErr (*RecVol)(Date, Date, double), /* function called whenever LGM uses a vol for pricing */
    double  NumPeriod,
    double* index,
    double* d,
    double* dstart,
    double* mktprice,
    /* information about the deal */
    LGMDealType dealtype,       /* type of deal */
    void*       dealPtr,        /* pointer to the deal */
                                /* calibration and eval parameters */
    LGMCalParm* CalReq,         /* structure defining calibration method to be used */
    ConvParams* EvalParms,      /* structure containing numerical convolution parameters */
                                /* outputs */
    double*     LGMValPtr,      /* value of the deal */
    double*     intrinValPtr,   /* intrinsic value of the deal */
    LGMSwptns** ExerIntoPtrPtr, /* European swaptions most similar to underlying */
    LGMSwptns** ExerBdryPtrPtr, /* European swaptions struck at the exercise boundary */
    SrtLgmRefSwptnData*
               lgmRefDeals,    /* ptr to reference swaption data structure (NULL => not req'd) */
    LGM_TS**   LGMtsPtrPtr,    /* calibrated zeta-G term structure */
    SigKapTS** SigKaptsPtrPtr, /* an equivalent sigma-kappa term structure */
    String     outfile);           /* output file name for log file (unused) */

/* Convert xExBdry[i] = x(t) at tExBdry[i] into equivalent swaptions which
would be ATM at the exercise boundary, for i = 0,...,nExBdry-1 */
static LGMSwptnsPtr FindBdry(
    int     iLGMOneTwoFactor,
    Date    tNow,    /* Evaluation date */
    String  ycName,  /* information about today */
    LGM_TS* tsPtr,   /* term structure */
    long    nExBdry, /* number of boundary points*/
    Date*   tExBdry, /* exercise boundary dates */
    double* xExBdry, /* exercise boundary (in x) */
    Date    tEnd);      /* last pay date */

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

/* Check calibration methods and data in the  request CR */
static void checkcalibrationmethod(LGMCalParm* CRptr);

/* Routine to construct the "1 into k" swaptions */
/* Note: RefDealsCommon must be called first */
static LGMErr RefDeals1intok(
    long         nArr,
    Date*        TauArr,
    double*      FVArr, /* exercise date and relative PV arrays */
    LGMMarkConv* conv,
    String       ycname, /* market conventions, discount curve */
    LGMCalParm*  CalReq, /* calibration methodology */
    LGMDealType  DealType,
    void*        dealPtr, /* deal type and params */
    LGMErr (*GetVol)(
        Date, Date, double, SRT_Boolean, double*), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date, double*), /* function to provide exponent (beta) to interpret vols */
    LGMCalSet* CSPtr);        /* Output */

/* Routine to construct the caplets */
/* Note: RefDealsCommon must be called first */
static LGMErr RefDealsCap(
    LGMDealType DealType,
    void*       DealPtr,
    long        nArr,
    Date*       TauArr,
    double*     FVArr,
    LGMCalParm* CalReq,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*),
    LGMErr (*GetBeta)(Date, Date, double*),
    String       YcName,
    LGMMarkConv* conv,
    LGMCalSet*   CSPtr); /*OUTPUT*/

/* Routine to construct the long swaptions */
/* Note: RefDealsCommon must be called first */
static LGMErr RefDealsLong(
    long         nArr,
    Date*        TauArr,
    double*      FVArr, /* exercise date and relative PV arrays */
    LGMMarkConv* conv,
    String       ycname, /* market conventions, discount curve */
    LGMCalParm*  CalReq, /* calibration methodology */
    LGMDealType  DealType,
    void*        dealPtr, /* deal type and params */
    LGMErr (*GetVol)(
        Date, Date, double, SRT_Boolean, double*), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date, double*), /* function to provide exponent (beta) to interpret vols */
    LGMCalSet* CSPtr);
/* Routine to construct the common parts of all reference instrumetns */
static LGMErr RefDealsCommon(
    Date         tNow,
    long         nEx, /* evaluation date, number of exercises */
    Date         tLast,
    Date*        TauArr, /* last pay date, exer dates */
    LGMCalParm*  CalReq, /* Methodology */
    LGMMarkConv* conv,
    String       ycname, /* market conventions, discount curve */
    LGMCalSet*   CSPtr);
/* Check n, the array of exercise dates and the array of forward values
Return the first i with t[i]>tNow, and the last t with t[i]<tLast */
static LGMErr checkinputarrays(
    long        n,
    Date        tNow,
    Date        tLast,
    Date*       t,
    double*     RatArr,
    LGMCalParm* CR,
    long*       firstiPtr,
    long*       lastiPtr);

/* Routine to find the fixed rates of the reference instruments */
static LGMErr FindFixRateRefInst(
    double       Rsw, /* swap rate */
    double       Lvl, /* level, adjusted for gearing */
    double       DStart,
    double       DEnd,
    double       DLast, /* df's to start date and end dates */
    Date         tEx,
    Date         tStart,
    Date         tEnd,   /* dates of refence swaption */
    long         islong, /* 1=use Rdata1, 0=use Rdata2 */
    long         nArr,
    double*      FVarr,
    Date*        TauArr, /* PV ratio & exer dates of orig deal */
    Date         tNow,
    Date         tLast,                                         /* end date of original deal */
    LGMMarkConv* conv,                                          /* market conventions */
    LGMCalParm*  CalReq,                                        /* calibration request */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* gets swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*),                     /* gets exp (beta) for vols */
    double*    rfix, /* returns fix rate of ref.swaption */
    double*    std,
    LGMCalSet* CSPtr);

/* Use linear interpolation (linear at ends) to find the PV ratio at thedate
with FVArr=1 at tLast */
double InterpPVratio(Date thedate, long n, double* FV, Date* Tau, Date tLast);

static double find_gear(Date thedate, Date* d, double* gear, int n);

static double LinearFlatFlat(Date thedate, long n, Date* t, double* Rdata);

static LGMSwptnsPtr FindBdry(
    int     iLGMOneTwoFactor,
    Date    tNow,    /* Evaluation date */
    String  ycName,  /* information about today */
    LGM_TS* tsPtr,   /* term structure */
    long    nEx,     /* number of boundary points*/
    Date*   tExBdry, /* exercise boundary dates */
    double* xExBdry, /* exercise boundary (in x) */
    Date    tEnd)       /* last pay date */

{
    LGMErr      error = NULL;
    LGMSwptns*  ExerBdryPtr;
    SrtCurvePtr yldcrv = NULL;
    LGMMarkConv conv;
    Date        tEx, tStart, tPay;
    double      Rfix;
    double      FloatRate;
    long        j, k;

    yldcrv = lookup_curve(ycName);
    error  = LGMCcyDefaults(yldcrv, &conv);
    if (error)
        return (NULL);

    /* allocate space */
    ExerBdryPtr = LGMCreateSwptns(nEx);
    if (ExerBdryPtr == NULL)
        return (NULL);

    k = -1;
    for (j = 0; j < nEx; j++)
    {
        tEx = tExBdry[j];
        if (tEx >= tNow && tEx < tEnd - 30)
        {
            if (iLGMOneTwoFactor == 1)
            {
                // error = LGMGetRfixFromX(tEx, tEnd, xExBdry[j], tNow, ycName, tsPtr, &Rfix);
                // if(error) return NULL;

                tStart = add_unit(tEx, conv.lag, SRT_BDAY, SUCCEEDING);
                tPay   = add_unit(tStart, 12 / (int)conv.cfreq, SRT_MONTH, conv.cbdconv);

                error = LGMGetFloatRateFromX(tStart, tPay, xExBdry[j], ycName, tsPtr, &FloatRate);
                if (error)
                    return NULL;
            }
            else
                Rfix = xExBdry[j];
            if (error == NULL)
            {
                k++;
                ExerBdryPtr->tEx[k]  = tEx;
                ExerBdryPtr->tEnd[k] = tEnd;
                // ExerBdryPtr->Rfix[k] = Rfix;
                ExerBdryPtr->Rfix[k] = FloatRate;
            }
        }
    }
    ExerBdryPtr->n = k + 1;
    if (k < 1)
        LGMFreeSwptns(&ExerBdryPtr);

    return (ExerBdryPtr);
}
/*************************************************************************************************************/
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

/************************************************************************************************************/
static LGMErr CompValCapFlt(
    Date   tNow,
    String ycName,                                              /* yield curve */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* cap vols */
    LGMErr (*GetBeta)(Date, Date, double*),                     /* cap exponent (beta) */
    double*         ratioPtr,                                   /* ratio */
    double*         intValPtr,                                  /* intrinsic value */
    SrtReceiverType PayRec,                                     /* pay/rec */
    Date            tStart,
    double          Strike, /* settlement, strike */
    long            n,
    Date*           t, /* num or periods, cpn dates */
    double*         amax,
    double*         amin,
    double*         marg,     /* handle max and min, margin*/
    double          lvg,      /* leverage */
    SrtBasisCode    cpnbasis, /* coupon basis */
    SrtBasisCode    basis)       /* basis for rates */
{
    SrtCurvePtr yldcrv = NULL; /* yield curve pointer */
    LGMMarkConv conventions;   /* standard market conventions */
    LGMErr      error = NULL;
    Date        tEx;
    long        i, lag;
    double      CEVvolmax, CEVvolmin, beta, Vcapmax, Vcapmin, pvCpn;
    double      df, dfSt, dfEnd, cvgCap, cvgCpn;
    double      Level, Forward, Yrtoexp;
    double      sum, ratio, intVal;

    yldcrv = lookup_curve(ycName);
    error  = LGMCcyDefaults(yldcrv, &conventions); /* get currency conventions */
    if (error != NULL)
        return (error);
    lag = conventions.lag; /* get lag */
    if (n < 1)
        return ("no coupons");

    sum = 0.;
    for (i = 1; i <= n; i++) /* find PV of coupon leg */
    {
        dfEnd = swp_f_df(tNow, t[i], ycName); /* calculate value of floor */
        dfSt  = swp_f_df(tNow, t[i - 1], ycName);
        if (dfSt == SRT_DF_ERROR || dfEnd == SRT_DF_ERROR)
            return ("no discount factor");
        cvgCap = coverage(t[i - 1], t[i], basis);
        tEx    = add_unit(t[i - 1], -lag, SRT_BDAY, SUCCEEDING);

        error = GetVol(t[i - 1], t[i], amax[i], SRT_TRUE, &CEVvolmax); /* Get CEVvol for cap max*/
        if (error)
            return (error);
        error = GetVol(t[i - 1], t[i], amin[i], SRT_TRUE, &CEVvolmin); /* Get CEVvol for cap min */
        if (error)
            return (error);
        error = GetBeta(t[i - 1], t[i], &beta); /* Get CEV exponent */
        if (error)
            return (error);

        Yrtoexp = (tEx - tNow) * YEARS_IN_DAY;
        Level   = cvgCap * dfEnd;
        Forward = (dfSt - dfEnd) / (Level);

        Vcapmax =
            LGMCEVCapletPrice(tNow, tEx, SRT_PAYER, amax[i], cvgCap, dfSt, dfEnd, CEVvolmax, beta);

        if (Vcapmax < 0.0)
            return ("can't get cap max price");
        Vcapmin =
            LGMCEVCapletPrice(tNow, tEx, SRT_PAYER, amin[i], cvgCap, dfSt, dfEnd, CEVvolmin, beta);

        cvgCpn = coverage(t[i - 1], t[i], cpnbasis);
        pvCpn  = cvgCpn * (marg[i] * dfEnd + lvg * (Vcapmin - Vcapmax) / cvgCap);
        sum    = sum + pvCpn;
    }
    sum = sum + dfEnd; /* add notional to end of coupon leg */

    df = swp_f_df(tNow, tStart, ycName);
    if (df == SRT_DF_ERROR)
        return ("no discount factor");

    intVal = sum - df * Strike; /* intrinsic value */
    if (PayRec == SRT_PAYER)
        intVal = -intVal;
    intVal = max(0., intVal);

    if (Strike == 0.)
        ratio = 0.0;
    else
        ratio = sum / (df * Strike); /* ratio of PV's */

    if (ratio > 5.0 || ratio < 0.2) /* clearly, one sign deal */
        ratio = 0.;

    *ratioPtr  = ratio;
    *intValPtr = intVal;
    return (NULL);
}

/*****************************************************************************************************************/
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

/*****************************************************************************************************************/

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

/********************************************************************************/
/*	Returns intrinsic value and the ratio of the PV of fixed leg to the PV of
the strike. If deal is strictly positive or strictly negative, ratio is set
to 0. Note that tNow is needed only for computing the intrinsic value */
static LGMErr CompValOpt(
    Date            tNow,
    String          ycName,
    double*         ratioPtr,
    double*         intValPtr,
    SrtReceiverType PayRec,
    Date            tStart,
    double          Strike,
    double          redfirstPay,
    long            nPay,
    Date*           tPay,
    double*         Payment)
{
    long   j;
    double sum, df, ratio, intVal;
    double signleg;

    signleg = 1.; /* intrin val = signleg*(fix leg - strike) */
    if (PayRec == SRT_PAYER)
        signleg = -1.;

    intVal = INTRINSICDEFAULT;
    ratio  = 0.;

    df = swp_f_df(tNow, tStart, ycName);
    if (df == SRT_DF_ERROR)
        return ("no discount factor");
    Strike = df * Strike;

    if (nPay < 1) /* no payments - escape */
        return ("no payments");

    sum = 0.;
    for (j = 0; j < nPay; j++) /* find PV of fixed leg */
    {
        df = swp_f_df(tNow, tPay[j], ycName);
        if (df == SRT_DF_ERROR)
            return ("no discount factor");
        sum = sum + df * Payment[j];
        if (j == 0)
            sum = sum - df * redfirstPay; /* account for reduced first payment */
    }
    intVal = max(0, signleg * (sum - Strike)); /* intrinsic value */

    if (Strike == 0.)
        ratio = 0.0;
    else
        ratio = sum / Strike; /* ratio of PV's */

    if (ratio > 5.0 || ratio < 0.2) /* clearly, one sign deal */
        ratio = 0.;

    *ratioPtr  = ratio;
    *intValPtr = intVal;
    return (NULL);
}

/********************************************************************************/
/*	Returns intrinsic value and the ratio of the PV of fixed leg to the PV of
the strike (strike is assumed to be the first payment). If deal is strictly
positive or strictly negative, ratio is set to 0. Note that tNow is needed only
for computing the intrinsic value */
static LGMErr CompValGen(
    Date    tNow,
    String  ycName,
    double* ratioPtr,
    double* intValPtr,
    long    nPay,
    Date*   tPay,
    double* Payment)
{
    long   j;
    double Strike, df, intVal, ratio;

    intVal = INTRINSICDEFAULT;
    ratio  = 0.;

    if (nPay < 1) /* no payments - escape */
        return ("no payments");

    for (j = 1; j < nPay; j++) /* find PV of all payments after strike */
    {
        df = swp_f_df(tNow, tPay[j], ycName);
        if (df == SRT_DF_ERROR)
            return ("no discount factor");
        intVal = intVal + df * Payment[j];
    }

    df = swp_f_df(tNow, tPay[0], ycName);
    if (df == SRT_DF_ERROR)
        return ("no discount factor");
    Strike = df * Payment[0];

    if (Strike == 0.) /* deal has one sign - escape */
        ratio = 0.0;
    else
        ratio = -(intVal / Strike); /* ratio of PV's */
    intVal = intVal + Strike;       /* intrinsic val */

    if (ratio > 5.0 || ratio < 0.2) /* clearly, one sign deal */
        ratio = 0.;

    *ratioPtr  = ratio;
    *intValPtr = intVal;
    return (NULL);
}

static Date StFromEx(Date tEx, int lag, int calbus, BusDayConv bdr)
{
    Date tStart;
    if (calbus == 1)
        tStart = add_unit(tEx, lag, SRT_BDAY, bdr);
    else
        tStart = add_unit(tEx, lag, SRT_DAY, bdr);
    return (tStart);
}

static Date ExFromSt(Date tStart, int lag, int calbus, BusDayConv bdr)
{
    Date tEx;
    if (calbus == 1)
        tEx = add_unit(tStart, -lag, SRT_BDAY, NO_BUSDAY_CONVENTION);
    else
    {
        tEx = add_unit(tStart, -lag, SRT_DAY, NO_BUSDAY_CONVENTION);
        while (add_unit(tEx, lag, SRT_DAY, bdr) > tStart)
            tEx--;
    }
    return (tEx);
}
#undef INTRINSICDEFAULT

/* Check calibration methods and data in the calibration
request CR. Change to default methods if necessary */
static void checkcalibrationmethod(LGMCalParm* CR)
{
    double lbound, ubound, slope, dt, isign;
    long   i, n;
    long   needdata1, needdata2;

    if (CR->calmeth >= LGMLastCalMeth)
        CR->calmeth = FixExp; /* default method */

    switch (CR->calmeth)
    {
    case FixKappa:
        if (CR->kap > 4.0) /* mean reversion timescale greater than 3 months */
            CR->kap = 4.0;
        else if (CR->kap < -4.0)
            CR->kap = -4.0;
        break;

    case GivenG:
        n = CR->numG;
        if (n < 2 || CR->Gdate == NULL || CR->G == NULL)
        {
            CR->calmeth = FixSigma; /* default method */
            break;
        }
        for (i = 1; i < n; i++) /* check date order */
        {
            if (CR->Gdate[i] <= CR->Gdate[i - 1])
            {
                CR->calmeth = FixSigma;
                break;
            }
        }

        isign = 1.0;
        if (CR->G[0] < CR->G[n - 1])
            isign = -1.0;
        CR->G[0] = isign * CR->G[0];
        for (i = 1; i < n; i++)
        {
            CR->G[i] = isign * CR->G[i];
            if (CR->G[i] > CR->G[i - 1]) /* ensure G(t) is decreasing */
                CR->G[i] = CR->G[i - 1];
        }

        dt    = (double)(CR->Gdate[n - 1] - CR->Gdate[0]);
        slope = (CR->G[0] - CR->G[n - 1]) * 365.0 / dt; /* re-scale */
        if (fabs(slope) < 1.e-8)                        /* G values too small? */
        {
            CR->calmeth = FixSigma;
            break;
        }
        for (i = 0; i < n; i++)
            CR->G[i] = (CR->G[i] - CR->G[n - 1]) / slope;

        for (i = 1; i < n; i++)
        {
            dt    = (double)(CR->Gdate[i] - CR->Gdate[i - 1]);
            slope = (CR->G[i - 1] - CR->G[i]) * 365.0 / dt;
            if (slope > 10.0)
            {
                CR->calmeth = FixSigma; /* G too ragged */
                break;
            }
        }
        break;

    case FixSigma: /* nothing to check */
        break;

    case GivenZeta:
        if (CR->numZ < 2 || CR->zdate == NULL || CR->zeta == NULL)
        {
            CR->calmeth = FixSigma; /* default method */
            break;
        }
        for (i = 1; i < CR->numZ; i++) /* check dates */
        {
            if (CR->zdate[i] <= CR->zdate[i - 1])
            {
                CR->calmeth = FixSigma;
                break;
            }
        }
        for (i = 1; i < CR->numZ; i++) /* ensure that zeta values are increasing */
        {
            if (CR->zeta[i] < CR->zeta[i - 1])
                CR->zeta[i] = CR->zeta[i - 1];
        }
        break;

    case FixExp:       /* checked later */
    case TenorAndDiag: /* nothing to check */
        break;
    } /* end switch */

    if (CR->MinMonToEx < 3)  /* minimum time to exercise at least three months */
        CR->MinMonToEx = 3;  /* for the 1 into k swaptions */
    if (CR->MinMonToEx > 12) /* minimum time to exercise less than twelve months */
        CR->MinMonToEx = 12; /* for the 1 into k swaptions */

    if (CR->keep1intok > 1)
        CR->keep1intok = 1;
    if (CR->keep1intok < -1)
        CR->keep1intok = -1;

    /* check strike data */
    if (CR->Rmeth >= LastRMeth)
    {
        CR->Rmeth = dIRR;
        return;
    }
    if (CR->Rmeth == IRR || CR->Rmeth == dIRR || CR->Rmeth == EMK1 ||
        CR->Rmeth == EMK2) /* Nothing to check */
        return;

    /* determine which data is needed */
    /* Rdata1 is used to construct strikes (Rfix) for long swaptions */
    /* Rdata2 is used to construct strikes (Rfix) for all other instruments */
    needdata1 = 0;
    needdata2 = 0;

    if (CR->usecaps == 1) /* All methods except Tenor & Diagonal and Fixed Expiry */
        needdata2 = 1;    /* use either long "k into n-k" swaptions */
    else                  /* or the "k into 3mo" caplets			*/
        needdata1 = 1;

    if (CR->calmeth == FixExp || CR->calmeth == TenorAndDiag)
    {
        needdata1 = 1; /* These methods use two sets of swaptions */
        needdata2 = 1;
    }

    if (CR->Rmeth == addshift)
    {
        lbound = -0.03;
        ubound = 0.03;
    } /* Rswap - 3% < Rfix < Rswap + 3% */
    if (CR->Rmeth == propshift)
    {
        lbound = 0.5;
        ubound = 1.0;
    } /* 0.5*Rswap < Rfix < 2.0*Rfix */
    if (CR->Rmeth == fracstd)
    {
        lbound = -2.0;
        ubound = 2.0;
    } /* Rswap - 2stddev < Rfix < Rswap + 2stddev */
    if (CR->Rmeth == givenR)
    {
        lbound = 0.0010;
        ubound = 0.99;
    } /* 10bps < Rfix < 99% */

    /* check data for long swaptions */
    if (needdata1 == 1)
    {
        if (CR->numR1 == 0 || CR->Rdate1 == NULL || CR->Rdata1 == NULL)
        {
            CR->Rmeth = dIRR;
            return;
        }
        for (i = 1; i < CR->numR1; i++)
        {
            if (CR->Rdate1[i] < CR->Rdate1[i - 1])
            {
                CR->Rmeth = dIRR;
                return;
            }
        }
        for (i = 0; i < CR->numR1; i++)
        {
            if (CR->Rdata1[i] < lbound)
                CR->Rdata1[i] = lbound;
            if (CR->Rdata1[i] > ubound)
                CR->Rdata1[i] = ubound;
        }
    }

    /* check data for all other swaptions/caps */
    if (needdata2 == 1)
    {
        if (CR->numR2 == 0 || CR->Rdate2 == NULL || CR->Rdata2 == NULL)
        {
            CR->Rmeth = dIRR;
            return;
        }
        for (i = 1; i < CR->numR2; i++)
        {
            if (CR->Rdate2[i] < CR->Rdate2[i - 1])
            {
                CR->Rmeth = dIRR;
                return;
            }
        }
        for (i = 0; i < CR->numR2; i++)
        {
            if (CR->Rdata2[i] < lbound)
                CR->Rdata2[i] = lbound;
            if (CR->Rdata2[i] > ubound)
                CR->Rdata2[i] = ubound;
        }
    }
    return;
}

/* Routine to construct the "1 into k" swaptions. For all these swaptions,
the exercise date is ftEx, which is the latter of the first exercise date tEx[1] or
MinMoEx months from today. The start dates of all these swaptions are ftStart, and
the fixed leg pay dates of swaption j are tPay[i], for i=ffirst,...,j. The swaptions run
from j=nfirst to j=nlast.
The coverage between paydates ftStart and ftPay[ffirst] is cvgfix, and
the discount factor to ftStart is DfixStart. */
/* Note: RefDealsCommon must be called first */

static LGMErr RefDeals1intok(
    long         nArr,
    Date*        TauArr,
    double*      FVArr, /* exercise date and relative PV arrays */
    LGMMarkConv* conv,
    String       ycname, /* market conventions, discount curve */
    LGMCalParm*  CalReq, /* calibration methodology */
    LGMDealType  DealType,
    void*        dealPtr, /* deal type and params */
    LGMErr (*GetVol)(
        Date, Date, double, SRT_Boolean, double*), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date, double*), /* function to provide exponent (beta) to interpret vols */
    LGMCalSet* CSPtr)         /* Output */
{
    long             MinMonEx, swmonth, fullper;
    Date             spot, ftStart, ftEx;
    Date             tLast;
    long             i, j, islong, lag, nfirst;
    double           reflevel, level, DLast, refsw, rfix;
    SrtCallInvFltPtr cinf;
    LGMErr           error;
    Date             tNow, *tPay;                       /* alias to simplify coding */
    long             ffirst, n, nPay;                   /* alias to simplify coding */
    double           cvgfix, DfixStart, *cvgpay, *Dpay; /* alias to simplify coding */

    /* initialize */
    error          = NULL;
    CSPtr->fixflag = 1;             /* flag 1 into k swaptions as existing */
    tNow           = CSPtr->tNow;   /* alias to simplify coding */
    n              = CSPtr->n;      /* alias to simplify coding */
    nPay           = CSPtr->nPay;   /* alias to simplify coding */
    tPay           = CSPtr->tPay;   /* alias to simplify coding */
    cvgpay         = CSPtr->cvgpay; /* alias to simplify coding */
    Dpay           = CSPtr->Dpay;   /* alias to simplify coding */
    lag            = conv->lag;     /* alias to simplify coding */

    swmonth = 6;                   /* swmonth is number of months per fixed leg period */
    fullper = 155;                 /* if period has more than fullper days, count it as full*/
    if (conv->sfreq == SRT_ANNUAL) /* period for fulfilling minimum swap length requirement */
    {
        swmonth = 12;
        fullper = 310;
    }
    MinMonEx = CalReq->MinMonToEx; /* minimum months to exercise */
    if (MinMonEx < 3)
        MinMonEx = 3;
    if (MinMonEx > 12)
        MinMonEx = 12;

    if (DealType == CallInvFlt)
    {
        cinf = (SrtCallInvFltPtr)dealPtr;
    }

    /* Find the exercise date for the 1 into k swaptions */
    spot    = add_unit(tNow, lag, SRT_BDAY, SUCCEEDING);       /* spot of today */
    ftStart = add_unit(spot, MinMonEx, SRT_MONTH, SUCCEEDING); /* plus MinMoEx months */
    ftEx    = add_unit(ftStart, -lag, SRT_BDAY, SUCCEEDING);   /* less lag bus days */

    if (ftEx < CSPtr->tEx[1])
        ftEx = CSPtr->tEx[1]; /* exer date for the 1 into k */

    ftStart = add_unit(ftEx, lag, SRT_BDAY, SUCCEEDING); /* start for 1 into k */

    CSPtr->ftEx    = ftEx;
    CSPtr->ftStart = ftStart;

    /* find first pay date for the 1 into k swaptions */
    for (i = 1; i <= nPay && tPay[i] <= ftStart; i++)
        ;
    ffirst    = i;
    cvgfix    = coverage(ftStart, tPay[i], conv->sbasis);
    DfixStart = swp_f_df(tNow, ftStart, ycname);
    if (DfixStart == SRT_DF_ERROR)
        return ("no discount factor");

    CSPtr->ffirst    = ffirst;
    CSPtr->cvgfix    = cvgfix;
    CSPtr->DfixStart = DfixStart;

    /* find last pay date of shortest and longest "1 into k" swaptions */
    if (tPay[i] - ftStart > fullper)
        nfirst = i + conv->minswap - 1;
    else
        nfirst = i + conv->minswap;

    if (nfirst > nPay)
        nfirst = nPay;

    CSPtr->nfirst = nfirst;
    if (nfirst > n)
        CSPtr->nlast = nfirst;
    else
        CSPtr->nlast = n;

    /* Allocate space for strikes and market price */
    CSPtr->Rfix = (double*)srt_calloc(nPay + 1, sizeof(double));
    CSPtr->Vfix = (double*)srt_calloc(nPay + 1, sizeof(double));
    if (CSPtr->Rfix == NULL || CSPtr->Vfix == NULL)
        return ("allocation failed in RefDeals1intok");

    /* Find fixed rates Rfix of 1 into k swaptions */
    tLast = tPay[n];
    DLast = Dpay[n];

    reflevel = cvgfix * Dpay[ffirst];
    for (i = ffirst + 1; i < nfirst; i++)
    {
        reflevel += cvgpay[i] * Dpay[i];
    }

    if (DealType == CallInvFlt)
    {
        level = cvgfix * Dpay[ffirst] *
                (1.0 + find_gear(tPay[ffirst], cinf->tCpnPay, cinf->gear, cinf->nCpn));
        for (i = ffirst + 1; i < nfirst; i++)
        {
            level += cvgpay[i] * Dpay[i] *
                     (1.0 + find_gear(tPay[i], cinf->tCpnPay, cinf->gear, cinf->nCpn));
        }
    }

    for (j = nfirst; j <= CSPtr->nlast; j++) /* for each 1 into k swaption */
    {
        reflevel += cvgpay[j] * Dpay[j];          /* compute forward level and swap rate */
        refsw = (DfixStart - Dpay[j]) / reflevel; /* cash swap rate */

        if (DealType == CallInvFlt)
        {
            level += cvgpay[j] * Dpay[j] *
                     (1.0 + find_gear(tPay[j], cinf->tCpnPay, cinf->gear, cinf->nCpn));
            level /= DfixStart;
        }
        else
        {
            level = reflevel / DfixStart;
        }

        islong = 0; /* use second set of Rdata if other set are long swaptions */
        if (CalReq->usecaps == 1)
            islong = 1; /* use first set of Rdata if other set are caps */

        error = FindFixRateRefInst(
            refsw, /* swap rate of ref swaption*/
            level,
            DfixStart,
            Dpay[j],
            DLast, /* df's for ref start and end, actual end */
            ftEx,
            ftStart,
            tPay[j], /* exer, start, and end date of ref swaption */
            islong,  /* which set of Rdata for strikes */
            nArr,
            FVArr,
            TauArr,
            tNow,
            tLast, /* PV ratio, exer, and end date of orig deal */
            conv,
            CalReq,
            GetVol,
            GetBeta, /* Calibration method, market vol and smile */
            &rfix,   /* strike of the reference swaption */
            NULL,
            CSPtr);

        if (error != NULL)
            return (error);

        CSPtr->Rfix[j] = rfix;
    }
    return (error);
}

/*******************************************/
/* Routine to construct the caplets */
/* Construct the caplets j, j=1,2, ..., nEx. Caplet j has
the fixing date tEx[j], the start date as tStart[j], and the end date tEnd[j].
Discount factors to the start date and end date are DStart[j] and Dcap[j], and
coverage from the start to end date is cvgcap[j].
The  new quantities we need are tEnd[], Dcap[], and cvgcap[]. The rest are
calculated in RefDealsCommon */
/* Note: RefDealsCommon must be called first to set nEx, tEx[], tStart[], and DStart[] */

static LGMErr RefDealsCap(
    LGMDealType DealType,
    void*       DealPtr,
    long        nArr,
    Date*       TauArr,
    double*     FVArr,
    LGMCalParm* CalReq,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*),
    LGMErr (*GetBeta)(Date, Date, double*),
    String       YcName,
    LGMMarkConv* conv,
    LGMCalSet*   CSPtr) /*OUTPUT*/
{
    LGMErr         error = NULL;
    long           i0, j0, j;
    long           islong;
    SrtCallInvFlt* BerInvFlt;
    SrtCallCapFlt* BerCapFlt;
    SrtSimMidAt*   MidAt;

    /* Allocate space for tEnd[], Dcap[], cvgcap[] */
    CSPtr->tEnd    = (Date*)srt_calloc(CSPtr->nEx + 1, sizeof(Date));
    CSPtr->cvgcap  = (double*)srt_calloc(CSPtr->nEx + 1, sizeof(double));
    CSPtr->Dcap    = (double*)srt_calloc(CSPtr->nEx + 1, sizeof(double));
    CSPtr->Rfcap   = (double*)srt_calloc(CSPtr->nEx + 1, sizeof(double));
    CSPtr->Vcap    = (double*)srt_calloc(CSPtr->nEx + 1, sizeof(double));
    CSPtr->VegaCap = (double*)srt_calloc(CSPtr->nEx + 1, sizeof(double));

    if (CSPtr->tEnd == NULL || CSPtr->cvgcap == NULL || CSPtr->Dcap == NULL ||
        CSPtr->Rfcap == NULL || CSPtr->Vcap == NULL)
        return ("Allocation Failed in RefDealsCap");

    CSPtr->capflag = 1;

    switch (DealType)
    {
    case CallCapFlt:
        BerCapFlt = (SrtCallCapFlt*)DealPtr;
        j0        = BerCapFlt->FirstEx;
        i0        = BerCapFlt->iSet[j0];
        break;

    case CallInvFlt:
        BerInvFlt = (SrtCallInvFlt*)DealPtr;
        break;

    case SimpleMidAt:
        MidAt = (SrtSimMidAt*)DealPtr;
        break;
    }

    /* Update Midat Structure for special methods of calibration */
    if (CalReq->LGMOneTwoFactor == 2 && DealType == SimpleMidAt &&
        ((CalReq->Rmeth == EMK1) || (CalReq->Rmeth == EMK2)))
    {
        error = update_autocal_midat_struct(MidAt, YcName, GetVol, GetBeta);
        if (error)
            return error;
    }

    /* Fill in CSPtr->tEnd[], CSPtr->cvgcap[], CSPtr->Dcap[] */
    for (j = 1; j <= CSPtr->nEx; j++) /*nEx IS ALSO THE COUPON NBR */
    {
        CSPtr->tEnd[j] = add_unit(
            CSPtr->tStart[j],
            12 / (int)conv->cfreq,
            SRT_MONTH,
            conv->cbdconv); /* NO_BUSDAY_CONVENTION or BUSDAY_CONVENTION...... */

        CSPtr->cvgcap[j] = coverage(
            CSPtr->tStart[j],
            CSPtr->tEnd[j],
            conv->cbasis); /* NO_BUSDAY_CONVENTION or BUSDAY_CONVENTION...... */

        CSPtr->Dcap[j] = swp_f_df(CSPtr->tNow, CSPtr->tEnd[j], YcName);
        if (CSPtr->Dcap[j] == SRT_DF_ERROR)
            return ("No Discount Factor");

        switch (DealType)
        {
        case CallCapFlt:
            CSPtr->Rfcap[j] = BerCapFlt->amax[j + i0];
            break;

        case CallInvFlt:
            islong = 0;
            error  = FindFixRateRefInst(
                (CSPtr->DStart[j] - CSPtr->Dcap[j]) / (CSPtr->cvgcap[j] * CSPtr->Dcap[j]),
                CSPtr->cvgcap[j] * CSPtr->Dcap[j] *
                    (1.0 +
                     find_gear(
                         CSPtr->tEnd[j], BerInvFlt->tCpnPay, BerInvFlt->gear, BerInvFlt->nCpn)),
                CSPtr->DStart[j],
                CSPtr->Dcap[j],
                CSPtr->Dpay[CSPtr->n],
                CSPtr->tEx[j],
                CSPtr->tStart[j],
                CSPtr->tEnd[j],
                islong,
                nArr,
                FVArr,
                TauArr,
                CSPtr->tNow,
                CSPtr->tPay[CSPtr->n],
                conv,
                CalReq,
                GetVol,
                GetBeta,
                &CSPtr->Rfcap[j],
                NULL,
                CSPtr);

            if (error != NULL)
                return (error);
            break;

        case SimpleMidAt:
        default:
            islong = 0;
            error  = FindFixRateRefInst(
                (CSPtr->DStart[j] - CSPtr->Dcap[j]) / (CSPtr->cvgcap[j] * CSPtr->Dcap[j]),
                CSPtr->cvgcap[j] * CSPtr->Dcap[j],
                CSPtr->DStart[j],
                CSPtr->Dcap[j],
                CSPtr->Dpay[CSPtr->n],
                CSPtr->tEx[j],
                CSPtr->tStart[j],
                CSPtr->tEnd[j],
                islong,
                nArr,
                FVArr,
                TauArr,
                CSPtr->tNow,
                CSPtr->tPay[CSPtr->n],
                conv,
                CalReq,
                GetVol,
                GetBeta,
                &CSPtr->Rfcap[j],
                NULL,
                CSPtr);

            if (error != NULL)
                return (error);
            break;
        }
    }

    return error;
}

/*******************************************/
/* Routine to construct the long swaptions */
/* Construct the long swaption j. j=1,2, ..., nEx. Long swaption j has
exercise date tEx[j], start date tStart[j], and fixed leg pay dates
tPay[i], i=ifirst[j] to i=nlong[j].
The discount factor to the start date is DStart[j], and
the discount factors to the pay dates tPay[i] are Dpay[i].
The coverage between start and the first pay date is cvgfirst[j], and
the coverages for the period between tPay[i-1] and tPay[i] is cvgpay[i] for i>ifirst[j]
*/
/* Note: RefDealsCommon must be called first to set nEx, n, nPay, tEx[], tStart[], tPay, cvgpay,
   Dpay, n, and nPay */

static LGMErr RefDealsLong(
    long         nArr,
    Date*        TauArr,
    double*      FVArr, /* exercise date and relative PV arrays */
    LGMMarkConv* conv,
    String       ycname, /* market conventions, discount curve */
    LGMCalParm*  CalReq, /* calibration methodology */
    LGMDealType  DealType,
    void*        dealPtr, /* For EMKi Method */
    LGMErr (*GetVol)(
        Date, Date, double, SRT_Boolean, double*), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date, double*), /* function to provide exponent (beta) to interpret vols */
    LGMCalSet* CSPtr)         /* Output */
{
    long               i, j, iend, islong;
    long               fullper;
    double             reflevel, level, DLast, refsw, rfix, std;
    Date               tNow, tLast;
    SrtCallInvFltPtr   cinf;
    SrtCallTimeSwapPtr ctimeswap;
    SrtSimMidAt*       MidAt;
    LGMErr             error;
    long               nEx, n, *ifirst, *nlong;           /* aliases to simplify code */
    Date *             tEx, *tStart, *tPay;               /* aliases to simplify code */
    double *           cvgpay, *Dpay, *cvgfirst, *DStart; /* aliases to simplify code */

    /* initialize */
    error  = NULL;          /* aliases to simplify code */
    nEx    = CSPtr->nEx;    /* aliases to simplify code */
    n      = CSPtr->n;      /* aliases to simplify code */
    tNow   = CSPtr->tNow;   /* aliases to simplify code */
    tEx    = CSPtr->tEx;    /* aliases to simplify code */
    tStart = CSPtr->tStart; /* aliases to simplify code */
    DStart = CSPtr->DStart; /* aliases to simplify code */
    cvgpay = CSPtr->cvgpay; /* aliases to simplify code */
    tPay   = CSPtr->tPay;   /* aliases to simplify code */
    Dpay   = CSPtr->Dpay;   /* aliases to simplify code */

    CSPtr->longflag = 1;           /* flag long swaptions as existing */
    fullper         = 155;         /* if period has more than fullper days, count it as full*/
    if (conv->sfreq == SRT_ANNUAL) /* period for fulfilling minimum swap length requirement */
        fullper = 310;

    if (DealType == CallInvFlt)
    {
        cinf = (SrtCallInvFltPtr)dealPtr;
    }
    else if (DealType == SimpleMidAt)
    {
        MidAt = (SrtSimMidAt*)dealPtr;
    }
    else if (DealType == CallTimeSwap)
    {
        ctimeswap = (SrtCallTimeSwapPtr)dealPtr;
    }

    /* Allocate space */
    CSPtr->ifirst   = (long*)srt_calloc(nEx + 1, sizeof(long));
    CSPtr->nlong    = (long*)srt_calloc(nEx + 1, sizeof(long));
    CSPtr->cvgfirst = (double*)srt_calloc(nEx + 1, sizeof(double));
    CSPtr->Rflong   = (double*)srt_calloc(nEx + 1, sizeof(double));
    CSPtr->Vlong    = (double*)srt_calloc(nEx + 1, sizeof(double));

    CSPtr->VegaLong = (double*)srt_calloc(nEx + 1, sizeof(double));
    CSPtr->StDLong  = (double*)srt_calloc(nEx + 1, sizeof(double));
    CSPtr->CEVLong  = (double*)srt_calloc(nEx + 1, sizeof(double));
    CSPtr->FrLong   = (double*)srt_calloc(nEx + 1, sizeof(double));
    CSPtr->BetaLong = (double*)srt_calloc(nEx + 1, sizeof(double));

    if (CSPtr->ifirst == NULL || CSPtr->nlong == NULL || CSPtr->cvgfirst == NULL ||
        CSPtr->Rflong == NULL || CSPtr->Vlong == NULL || CSPtr->CEVLong == NULL ||
        CSPtr->BetaLong == NULL || CSPtr->FrLong == NULL)

        return ("allocation failed in RefDealsLong");

    ifirst   = CSPtr->ifirst;   /* aliases to simplify code */
    nlong    = CSPtr->nlong;    /* aliases to simplify code */
    cvgfirst = CSPtr->cvgfirst; /* aliases to simplify code */

    /* for each tEx[j], find the first pay date tPay[i] which has tPay[i]>tEx[j] */
    /* if a pay date tPay[i] lies between tEx[j] and tStart[j], change tPay[i] to tStart[j] */
    i = 1;
    for (j = 1; j <= nEx; j++)
    {
        for (; i <= n && tPay[i] < tEx[j]; i++)
            ;
        if (tPay[i] <= tStart[j])
        {
            tPay[i] = tStart[j];
            i++;
        }
        ifirst[j] = i; /* first pay date strictly after tStart[j] */
        if (i > n && j > 1)
        {
            nEx        = j - 1; /* eliminate exercise dates whose start date occurs after tLast */
            CSPtr->nEx = nEx;
        }
    }

    /* find coverage from start to first pay date for each swaption */
    for (j = 1; j <= nEx; j++)
    {
        i           = ifirst[j];
        cvgfirst[j] = coverage(tStart[j], tPay[i], conv->sbasis);
    }

    /* last pay date for long swaptions */
    for (j = 1; j <= nEx; j++)
    {
        i = ifirst[j];
        if (tPay[i] - tStart[j] > fullper)
            nlong[j] = i + conv->minswap - 1; /* start to ifirst counts as full period */
        else
            nlong[j] = i + conv->minswap; /* start to ifirst too short to be full period */

        if (nlong[j] < n) /* long swaption ends at paydate n, unless can't */
            nlong[j] = n; /* fit minimum length swaption in */
    }

    /* Find fixed rate for long swaptions */
    tLast = tPay[n];
    DLast = Dpay[n];

    /* Update Midat Structure for special methods of calibration */
    if (CalReq->LGMOneTwoFactor == 2 && DealType == SimpleMidAt &&
        ((CalReq->Rmeth == EMK1) || (CalReq->Rmeth == EMK2)))
    {
        error = update_autocal_midat_struct(MidAt, ycname, GetVol, GetBeta);
        if (error)
            return error;
    }

    for (j = 1; j <= nEx; j++)
    {
        i    = ifirst[j];
        iend = nlong[j];

        reflevel = cvgfirst[j] * Dpay[i]; /* first period may be broken */
        i++;
        for (; i <= iend; i++)
        {
            reflevel = reflevel + cvgpay[i] * Dpay[i];
        }
        refsw = (DStart[j] - Dpay[iend]) / reflevel; /* cash swap rate */

        if (DealType == CallInvFlt)
        {
            i     = ifirst[j];
            level = cvgfirst[j] * Dpay[i] *
                    (1.0 + find_gear(tPay[i], cinf->tCpnPay, cinf->gear, cinf->nCpn));
            i++;
            for (; i <= iend; i++)
            {
                level += cvgpay[i] * Dpay[i] *
                         (1.0 + find_gear(tPay[i], cinf->tCpnPay, cinf->gear, cinf->nCpn));
            }
            level /= DStart[j];
        }
        else
        {
            level = reflevel / DStart[j];
        }

        islong = 1; /* use first set of Rdata if relevant */

        /* CALIBRATION AT THE MONEY SWAPTIONS FOR CALLABLE TIME SWAPS. */

        if (DealType == CallTimeSwap && ctimeswap->calibrationmethod == 1)
        {
            rfix = refsw;
            std  = 0;
        }
        else /* CALIBRATION AT THE EQUIVALENT STRIKE */
        {
            error = FindFixRateRefInst(
                refsw, /* swap rate of ref swaption */
                level,
                DStart[j],
                Dpay[iend],
                DLast, /* df's for ref start and end, actual end */
                tEx[j],
                tStart[j],
                tPay[iend], /* exer, start, and end date of ref swaption */
                islong,     /* which set of Rdata to use for strikes */
                nArr,
                FVArr,
                TauArr,
                tNow,
                tLast, /* PV ratio, exer, and end date of orig deal */
                conv,
                CalReq,
                GetVol,
                GetBeta, /* Calibration method, market vol and smile */
                &rfix,   /* strike of the reference swaption */
                &std,
                CSPtr);
        }

        if (error != NULL)
            return (error);

        CSPtr->Rflong[j]  = rfix;
        CSPtr->StDLong[j] = std;
    }

    return error;
}

static LGMErr RefDealsCommon(
    Date         tNow,
    long         nEx, /* evaluation date, number of exercises */
    Date         tLast,
    Date*        TauArr, /* last pay date, exer dates */
    LGMCalParm*  CalReq, /* Methodology */
    LGMMarkConv* conv,
    String       ycname, /* market conventions, discount curve */
    LGMCalSet*   CSPtr)    /* Output */
{
    Date *tEx, *tStart, *tPay;
    long  extraperiods, j, n, nPay;
    int   ResExer;

    /* initialize */
    CSPtr->longflag  = 0; /* flag all reference instruments as not existing */
    CSPtr->shortflag = 0;
    CSPtr->capflag   = 0;
    CSPtr->fixflag   = 0;
    CSPtr->tNow      = tNow;
    ResExer          = CalReq->respectExer;

    /* Set up fix leg for a standard swap with exercise date tfirstEx and last
    pay date tLast. For safety, add extraperiods periods to end of leg */
    /* Find fixed leg pay dates tPay[1,...,nPay]; start date is tPay[0] */
    extraperiods = conv->minswap;
    tPay         = LGMFixLegSched(TauArr[1], tLast, extraperiods, conv, &nPay);
    n            = nPay - extraperiods;

    CSPtr->n    = n;
    CSPtr->nPay = nPay;
    CSPtr->tPay = tPay;

    if (tPay == NULL || n < 1)
        return ("failure in RefDealsCommon");

    CSPtr->cvgpay = LGMFixLegCvg(nPay, tPay, conv);        /* coverages */
    CSPtr->Dpay   = LGMFixLegDF(nPay, tPay, tNow, ycname); /* discount factors to pay dates */

    if (CSPtr->cvgpay == NULL || CSPtr->Dpay == NULL)
        return ("allocation failed in RefDealsCommon");

    /* Allocate space for common exer dates, start dates, and dis factors to start dates */
    if (ResExer == 0)
        CSPtr->nEx = n;
    else
        CSPtr->nEx = nEx;

    CSPtr->tEx = tEx = (Date*)srt_calloc(CSPtr->nEx + 1, sizeof(Date));
    CSPtr->tStart = tStart = (Date*)srt_calloc(CSPtr->nEx + 1, sizeof(Date));
    CSPtr->DStart          = (double*)srt_calloc(CSPtr->nEx + 1, sizeof(double));

    if (tEx == NULL || tStart == NULL || CSPtr->DStart == NULL)
        return ("allocation failed in RefDealsCommon");

    /* Fill in common exercise dates, start dates, and dis factors to start dates */
    if (ResExer != 0)
    {
        for (j = 1; j <= nEx; j++)
        {
            tEx[j]    = TauArr[j];
            tStart[j] = add_unit(tEx[j], conv->lag, SRT_BDAY, SUCCEEDING);
        }
    }
    else
    {
        for (j = 1; j <= n; j++)
        {
            tStart[j] = tPay[j - 1];
            tEx[j]    = add_unit(tStart[j], -conv->lag, SRT_BDAY, SUCCEEDING);
        }
        tEx[1] = TauArr[1]; /* shouldn't make a difference */
    }

    for (j = 1; j <= CSPtr->nEx; j++)
    {
        CSPtr->DStart[j] = swp_f_df(tNow, tStart[j], ycname);
        if (CSPtr->DStart[j] == SRT_DF_ERROR)
            return ("no discount factor");
    }

    return (NULL); /* done */
}

static LGMErr checkinputarrays(
    long        n,
    Date        tNow,
    Date        tLast,
    Date*       t,
    double*     RatArr,
    LGMCalParm* CR,
    long*       firstiPtr,
    long*       lastiPtr)
{
    long i, needratios;
    long firsti, lasti;

    needratios = 0;
    if (CR->Rmeth == IRR || CR->Rmeth == dIRR) /* these methods choose the strikes Rfix */
        needratios = 1;                        /* based on the fwd value ratios of the
                                                               /* underlying deal */
    if (n < 1 || t == NULL)
        return ("no exer dates for calibration");

    for (i = 1; i < n; i++)
    {
        if (t[i] > t[i + 1])
            return ("exer dates out of order in calib");
    }

    /* find first date in t which is after tNow */
    for (firsti = 1; firsti <= n && t[firsti] <= tNow; firsti++)
        ;
    if (firsti > n)
        return ("no exercise dates left");
    *firstiPtr = firsti;

    /* find the last date in t which is before tLast */
    for (lasti = n; lasti > 0 && t[lasti] >= tLast; lasti--)
        ;
    if (lasti < firsti)
        return ("no exercise dates before tLast");
    *lastiPtr = lasti;

    if (needratios > 0)
    {
        if (RatArr == NULL)
            return ("no fwd values in calibrator");
        for (i = firsti; i <= lasti; i++)
        {
            if (RatArr[i] < 0.25 || RatArr[i] > 4.0)
                return ("bad fwd values");
        }
    }
    return (NULL);
}

static LGMErr FindFixRateRefInst(
    double       Rsw, /* swap rate for ref deal */
    double       Lvl,
    double       DStart,
    double       DEnd,
    double       DLast, /* df's to start date and end dates */
    Date         tEx,
    Date         tStart,
    Date         tEnd,   /* dates of refence swaption */
    long         islong, /* 1 = use Rdata1, 0 = use Rdata2 */
    long         nArr,
    double*      FVArr,
    Date*        TauArr, /* PV ratio & exer dates of orig deal */
    Date         tNow,
    Date         tLast,                                         /* end date of original deal */
    LGMMarkConv* conv,                                          /* market conventions */
    LGMCalParm*  CalReq,                                        /* calibration request */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* gets swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*),                     /* gets exp (beta) for vols */
    double*    rfix, /* returns fix rate of ref.swaption */
    double*    std,
    LGMCalSet* CSPtr)
{
    LGMErr           error = NULL;
    double           stddev, CEVvol, ATMVol, StDLong, beta;
    double           gamma, theta1, theta2;
    Date             tEx2;
    SrtDiffusionType input_vol;
    long             j;
    double           numer, MidatStrike;

    /*	Calculate standard deviation */
    error = GetVol(tStart, tEnd, Rsw, SRT_FALSE, &CEVvol);
    if (error)
        return error;
    error = GetBeta(tStart, tEnd, &beta);
    if (error)
        return error;

    if (beta == 0)
    {
        input_vol = SRT_NORMAL;
    }
    else
    {
        input_vol = SRT_LOGNORMAL;
    }

    stddev = fabs(CEVvol) * pow(Rsw, beta);
    stddev *= sqrt(((double)(tEx - tNow)) / 365.);

    if (CalReq->Rmeth != IRR && CalReq->Rmeth != dIRR && CalReq->Rmeth != EMK1 &&
        CalReq->Rmeth != EMK2)
    {
        if (islong == 1)
            gamma = LinearFlatFlat(tEx, CalReq->numR1, CalReq->Rdate1, CalReq->Rdata1);
        else
            gamma = LinearFlatFlat(tEx, CalReq->numR2, CalReq->Rdate2, CalReq->Rdata2);
    }

    /*	Calculate midat strike */
    theta1      = InterpPVratio(tEx, nArr, FVArr, TauArr, tLast);
    tEx2        = add_unit(tEnd, -conv->lag, SRT_BDAY, SUCCEEDING);
    theta2      = InterpPVratio(tEx2, nArr, FVArr, TauArr, tLast);
    numer       = theta1 - theta2 * DEnd / DStart;
    MidatStrike = numer / Lvl;

    switch (CalReq->Rmeth)
    {
    case IRR:

        /*
                        theta = InterpPVratio(tEx, nArr, FVArr, TauArr, tLast);
                        if (DLast > DStart/(1.0 + 0.25*Rsw))
                                DLast = DStart/(1.0 + 0.25*Rsw);
                        *rfix = Rsw*(1.0 + (theta-1.0)/(1.0-DLast/DStart));
                        break;
        */
        return "IRR Not available any more - use strike choices 1 (IRR), 2 (EMK1) or 3 (EMK2)";
        break;

    case dIRR:

        *rfix = MidatStrike;
        break;

    case addshift:
        *rfix = Rsw + gamma;
        break;

    case propshift:
        *rfix = Rsw * (1.0 + gamma);
        break;

    case fracstd:
        *rfix = Rsw + gamma * stddev;
        break;

    case givenR:
        *rfix = gamma;
        break;

    case EMK1:

        if (islong == 1)
        {
            *rfix = MidatStrike;
            *std  = (MidatStrike - Rsw) / stddev;
        }
        else if (islong == 0)
        {
            j = BasicDichotomie(CSPtr->tEx, 1, CSPtr->nEx, tEx);

            StDLong = (CSPtr->StDLong[j]);

            error = GetVol(tStart, tEnd, Rsw, SRT_TRUE, &ATMVol);

            *rfix = Rsw + StDLong * stddev;
        }

        break;

    case EMK2:

        if (islong == 1)
        {
            *rfix = MidatStrike;
        }
        else if (islong == 0)
        {
            *rfix = Rsw;
        }

        break;
    }

    /* ensure that the fixed rate is within 2 stddev's of par */
    if (*rfix < Rsw - CalReq->maxstd * stddev)
    {
        *rfix = Rsw - CalReq->maxstd * stddev;
    }
    else if (*rfix > Rsw + CalReq->maxstd * stddev)
    {
        *rfix = Rsw + CalReq->maxstd * stddev;
    }

    return error;
}

static double find_gear(Date thedate, Date* d, double* gear, int n)
{
    int i;

    i = 0;
    while (i < n - 1 && d[i] < thedate)
    {
        i++;
    }

    return gear[i];
}

/* Helper routine  for picking strikes */
/* Use interpolation to find the PV ratio at thedate */
/* Linear at front, linear at back with the proviso that
        FVArr=1 at tLast
*/
double InterpPVratio(Date thedate, long n, double* FV, Date* Tau, Date tLast)
{
    double slope, ans;
    long   i;

    if (n == 1 || thedate >= Tau[n])
    {
        slope = (1.0 - FV[n]) / ((double)(tLast - Tau[n]));
        ans   = 1.0 + slope * ((double)(thedate - tLast));
        return (ans);
    }
    else if (thedate <= Tau[1])
    {
        slope = (FV[2] - FV[1]) / ((double)(Tau[2] - Tau[1]));
        ans   = FV[1] + slope * ((double)(thedate - Tau[1]));
        return (ans);
    }
    else
        ;
    for (i = 2; Tau[i] <= thedate; i++)
        ;
    slope = (FV[i] - FV[i - 1]) / ((double)(Tau[i] - Tau[i - 1]));
    ans   = FV[i] + slope * ((double)(thedate - Tau[i]));
    return (ans);
}

static double LinearFlatFlat(Date thedate, long n, Date* t, double* Rdata)
{
    double slope, ans;
    long   i;

    if (n == 1 || thedate >= t[n - 1])
    {
        ans = Rdata[n - 1];
        return (ans);
    }
    else if (thedate <= t[0])
    {
        ans = Rdata[0];
        return (ans);
    }

    for (i = 1; t[i] <= thedate; i++)
        ;
    slope = (Rdata[i] - Rdata[i - 1]) / ((double)(t[i] - t[i - 1]));
    ans   = Rdata[i] + slope * ((double)(thedate - t[i]));
    return (ans);
}

static double FitOneSwap(
    long            n,
    double*         a,
    double          DStart,
    SrtReceiverType recpay, /* the swaption */
    double          zeta,   /* zeta */
    double*         Gstptr,
    double*         Gpay, /* the existing term structure */
    double          dGst,
    double*         dGpay, /* direction of change */
    double          price)          /* market price of swaption */
{
    SrtCallPutType callput;
    double         Gst;
    double         fwd, sens, dsens, theta, favg, stddev, ratio;
    double         arg, value, newvalue, der, ystar, sqzeta, coef, dcoef;
    double         timevalue, isign;
    long           i, iter, itermax, success;

    /* initialize */
    itermax = 25;
    callput = SRT_CALL;
    zeta    = fabs(zeta);
    sqzeta  = sqrt(zeta);
    Gst     = *Gstptr;

    /* check for intrinsic */
    /*	if (fabs (LGMRecVal(&ystar, n, a, DStart, 1.0e-15, Gpay, Gst) - price) < 1.0e-08)
            {
                    return 1.0e-15;
            }
    */
    /* compute forward value and sensitivity of swaption */
    fwd   = 0.;
    sens  = 0.0;
    dsens = 0.0;
    for (i = 1; i <= n; i++)
    {
        fwd   = fwd + a[i];
        sens  = sens + a[i] * (Gst - Gpay[i]);
        dsens = dsens + a[i] * (dGst - dGpay[i]);
    }
    theta = fwd - DStart;
    if (recpay == SRT_PAYER) /* need only consider receivers now */
        price = price + theta;

    timevalue = price - theta;
    if (timevalue <= 0.0) /* best we can do is C=0 */
        return (zeta);

    if (sens < -1.e-6) /* slopes of G are wrong */
        return (-1.0);

    if (dsens < 0.0)
    {
        for (i = 1; i <= n; i++)
            dGpay[i] = -dGpay[i];
        dGst  = -dGst;
        dsens = -dsens;
    }

    if (dsens < 1.e-05)
        return (-1.); /* can't fit, no moment */

    /* see if the swaption can be fit */
    if (sens >= 1.e-04)
    {
        value = LGMRecVal(&ystar, n, a, DStart, sqzeta, Gpay, Gst);
        if (price <= value) /* best fit is coef=0 */
            return (zeta);
    }

    /* get initial guess by linear approximation (neglects convexity) */
    favg   = 0.5 * (fwd + DStart);
    stddev = LGMNormImpVol(price / favg, fwd / favg, DStart / favg, 1.0, 1.0, callput);
    stddev = stddev * favg;
    if (stddev < 0.0)
        return (-1.0);

    ratio = sqzeta / stddev;

    /* global newton to get initial guess */
    coef    = 0.;
    dcoef   = (1.0 - ratio * sens) / (ratio * dsens);
    success = -1;
    value   = -1.0;
    for (iter = 1; iter <= itermax && success <= 0; iter++)
    {
        newvalue = -1.0; /* calculate func & der at C+dC */
        der      = 0.;
        for (i = 1; i <= n; i++)
        {
            arg = ratio * (Gst - Gpay[i] + (coef + dcoef) * (dGst - dGpay[i]));
            if (fabs(theta) <= 1.e-5 || fabs(arg) <= 1.e-4)
                newvalue =
                    newvalue + a[i] * arg * (1. - 0.5 * arg * theta * (1.0 - arg * theta / 3.0));
            else
                newvalue = newvalue + a[i] * (1.0 - exp(-theta * arg)) / theta;
            der = der + a[i] * ratio * (dGst - dGpay[i]) * exp(-theta * arg);
        }

        if (fabs(newvalue) < fabs(value)) /* if improvement, accept step */
        {
            coef  = coef + dcoef;
            value = newvalue;
            dcoef = -newvalue / der; /* new Newton step */
        }
        else
            dcoef = 0.5 * dcoef; /* else cut step in half */

        if (fabs(newvalue) <= 1.e-06) /* check convergence */
            success++;                /* go one step past OK convergence */
    }
    coef = coef + dcoef;

    /* Use global Newton on coef to match LGM price to market price */
    dcoef   = coef;
    coef    = 0.;
    value   = -price;
    success = -2;
    isign   = 1.0;

    for (iter = 1; iter <= itermax && success <= 0; iter++)
    {
        coef = coef + isign * dcoef; /* update term structure */
        Gst  = Gst + isign * dcoef * dGst;
        for (i = 1; i <= n; i++)
            Gpay[i] = Gpay[i] + isign * dcoef * dGpay[i];

        newvalue = LGMRecVal(&ystar, n, a, DStart, sqzeta, Gpay, Gst) - price;

        if (fabs(newvalue) / (timevalue + 0.001) < 1.e-5 ||
            fabs(dcoef) < 1.e-6) /* check convergence */
            success++;

        if (fabs(newvalue) < fabs(value))
        {
            value = newvalue; /* step successful */
            isign = 1.0;      /* take new step */
            der   = RecGDer(ystar, n, a, DStart, sqzeta, Gpay, Gst, dGst, dGpay);
            dcoef = -value / der; /* Newton step */
        }
        else
        {
            isign = -1.0;        /* step unsuccessful */
            dcoef = 0.5 * dcoef; /* back up half step */
        }
    }

    if (success >= 0)
    {
        coef = coef + isign * dcoef;       /* update term structure */
        Gst  = Gst + isign * dcoef * dGst; /* for last step */
        for (i = 1; i <= n; i++)
            Gpay[i] = Gpay[i] + isign * dcoef * dGpay[i];

        *Gstptr = Gst;
        return (zeta);
    }
    else
        return (-1.);
}

LGMErr LGMamortisingnewautocal(
    /* task list */
    int skipEval,     /* 0=calibrate & value deal, 1=calibrate only */
    int skipCalib,    /* 0=calibrate & value deal, 1=value deal only */
    int convertTS,    /* 1=compute & output sig-kappa ts equivalent to LGM ts */
    int findExerBdry, /* 1=find swap rates at exer boundary  */
                      /* information about today & today's market prices */
    Date   tNow,      /* eval as if today is tNow */
    int    eod,       /* 0=can exercise on tNow, 1=cannot exercise on tNow */
    String ycName,    /* market pointer for discount factors */
    LGMErr (*GetVol)(
        Date, Date, double, SRT_Boolean, double*), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date, double*), /* function to provide exponent (beta) to interpret vols */
    LGMErr (*RecVol)(Date, Date, double), /* function called whenever LGM uses a vol for pricing */
    double  NumPeriod,
    double* index,
    double* d,
    double* dstart,
    double* mktprice,
    /* information about the deal */
    LGMDealType dealtype,       /* type of deal */
    void*       dealPtr,        /* pointer to the deal */
                                /* calibration and eval parameters */
    LGMCalParm* CalReq,         /* structure defining calibration method to be used */
    ConvParams* EvalParms,      /* structure containing numerical convolution parameters */
                                /* outputs */
    double*     LGMValPtr,      /* value of the deal */
    double*     intrinValPtr,   /* intrinsic value of the deal */
    LGMSwptns** ExerIntoPtrPtr, /* European swaptions most similar to underlying */
    LGMSwptns** ExerBdryPtrPtr, /* European swaptions struck at the exercise boundary */
    SrtLgmRefSwptnData*
               lgmRefDeals,    /* ptr to reference swaption data structure (NULL => not req'd) */
    LGM_TS**   LGMtsPtrPtr,    /* calibrated zeta-G term structure */
    SigKapTS** SigKaptsPtrPtr, /* an equivalent sigma-kappa term structure */
    String     outfile)            /* output file name for log file (unused) */
{
    LGMErr       error = NULL;
    LGMMarkConv  conventions;        /* standard swaption conventions */
    LGMCalSetPtr RefDealsPtr = NULL; /* set of reference instruments */
    Date         tfirst;             /* first possible exercise date = tNow + eod */
    long         nEx;                /* number of exercise dates */

    long       EffnEx;                      /* effective number of exercise dates */
    Date       LasttPay, *TauArr = NULL;    /* info extracted from deal for calibrator */
    double*    FVArr  = NULL;               /* info extracted from deal for calibrator */
    double     LGMVal = 0., intrinVal = 0.; /* the value and intrinsic value of the deal */
    LGMSwptns* ExerIntoPtr = NULL;          /* European swaptions most similar to underlying */
    LGMSwptns* ExerBdryPtr = NULL;          /* European swaptions struck at the exercise boundary */
    LGM_TS*    LGMtsPtr    = NULL;          /* calibrated zeta-G term structure */
    SigKapTS*  SigKaptsPtr = NULL;          /* an equivalent sigma-kappa term structure */

    long        nExBdry;          /* number of points in the exercise boundary */
    Date*       tExBdry   = NULL; /* exer bdry is x(t)=xExBdry[i] at t=tExBdry[i] */
    double*     xExBdry   = NULL; /* for i=0,1,...,nExBdry-1 */
    SrtCurvePtr yldcrv    = NULL; /* yield curve pointer */
    long        startover = 0;    /* to speed up debugging */

    double kappa_term, yr_to_exp, prev_yr_to_exp, one_kappa;
    long   i;

    yldcrv = lookup_curve(ycName);

    if (findExerBdry != 0 && *ExerBdryPtrPtr != NULL) /* clear for output */
        LGMFreeSwptns(ExerBdryPtrPtr);

    if (skipCalib == 1)
    {
        error = LGMVerifyLGM_TS(*LGMtsPtrPtr);
        if (error)
            skipCalib = 0; /* can't skip if no valid term structue */
    }

    if (skipCalib != 1 && *LGMtsPtrPtr != NULL)
        if (CalReq->LGMOneTwoFactor == 1)
            LGMFreeLGM_TS(LGMtsPtrPtr);
        else if (CalReq->LGMOneTwoFactor == 2)
            LGMFreeLGM2F_TS(LGMtsPtrPtr);
    if (skipCalib != 1 && *ExerIntoPtrPtr != NULL)
        LGMFreeSwptns(ExerIntoPtrPtr); /* clear for output */
    if (skipCalib != 1 && lgmRefDeals != NULL)
        LGMFreeRefSwptnData(lgmRefDeals); /* clear for output */
    if (convertTS != 0 && *SigKaptsPtrPtr != NULL)
        LGMFreeSigKapTS(SigKaptsPtrPtr); /* clear for output */

    error = LGMCcyDefaults(yldcrv, &conventions); /* get currency code */
    if (error)
        return (error);

    tfirst = tNow;
    if (eod == 1)
        tfirst = add_unit(tNow, 1, SRT_BDAY, SUCCEEDING /*, conv.ccy*/);

    /* Step 1: Check validity of deal */
    switch (dealtype)
    {
    case GenEur:
        error = LGMValidGenEur((SrtGenEur*)dealPtr, tfirst, &nEx);
        break;
    case SimpleMidAt:
        error = LGMValidSimMidAt((SrtSimMidAt*)dealPtr, tfirst, &nEx);
        break;
    case GenMidAt:
        error = LGMValidGenMidAt((SrtGenMidAt*)dealPtr, tfirst, &nEx);
        break;
    case CallInvFlt:
        error = LGMValidCallInvFlt((SrtCallInvFlt*)dealPtr, tfirst, &nEx);
        break;
    case CallTimeSwap:
        error = LGMValidCallTimeSwap((SrtCallTimeSwap*)dealPtr, tfirst, &nEx);
        break;
    case CallCapFlt:
        error = LGMValidCallCapFlt((SrtCallCapFlt*)dealPtr, tfirst, &nEx);
        break;
    }
    if (error)
        return (error); /* problem with the deal */

    if (nEx <= 0)
    {
        *LGMValPtr    = 0.0; /* tfirst is strictly after last exer Date */
        *intrinValPtr = 0.0;
        return (NULL);
    }

    /* Step 2: Extract the relevant characteristics of the deal for calibration
    The effective number of exercise dates and the exercise dates - EffnEx, TauArr[]
    The PV of fixed leg/PV of strike at each exercise date - FVArr[]
    The last pay date of the deal - LasttPay

    Also, compute the intrinsic value of the deal, and store the underlying
    exercise dates, enddates, and the fixed rate Rfix of the vanilla swaption
    that comes closest to describing the underlying in ExerInto structure,
    and return these to the calling routine as "ExerIntoPtrPtr" */

    error = LGMExtractInfoFromDeal(
        tNow,
        tfirst,
        ycName, /* info about today's market */
        GetVol,
        GetBeta, /* function to provide swaption/cap vols */
        dealtype,
        dealPtr, /* the deal */
        &EffnEx,
        &LasttPay,
        /* NBR OF EXERCISE DATES, THE LAST PAYT DATE  */ /* THE MID-AT DEAL THEO END DATE. */
        &TauArr,
        &FVArr, /* info for calibrator */
        &intrinVal,
        &ExerIntoPtr); /* other info to return as output */

    if (error)
    {
        srt_free(TauArr);
        srt_free(FVArr);
        LGMFreeSwptns(&ExerIntoPtr);
        return (error);
    }

    if (EffnEx == 1 && tfirst == TauArr[1]) /* only exercise date is on tfirst */
    {
        *LGMValPtr = intrinVal; /* LGM value is the intrinsic value */
        srt_free(TauArr);
        srt_free(FVArr);
        LGMFreeSwptns(&ExerIntoPtr);
        return (NULL);
    }

    /* Step 3: Create and calibrate a LGM ts */
    if (skipCalib != 1)
    {
        error = LGMamortisingCalibration(
            tNow,
            ycName, /* info about today's market */
            GetVol,
            GetBeta,
            RecVol, /* info about today's swaption/caplet prices */
            NumPeriod,
            index,
            d,
            dstart,
            mktprice,
            dealtype,
            dealPtr,
            EffnEx,
            LasttPay, /* info about the deal */
            TauArr,
            FVArr,
            /* info about the deal */
            CalReq,       /* the requested calibration method */
            &LGMtsPtr,    /* the calibrated term structure in LGM */
            &RefDealsPtr, /* reference deals used for calibration */
            lgmRefDeals); /* Reference swaptions stock structure */

        LGMFreeCalSet(&RefDealsPtr); /* free unneeded structures */
    }
    else
        LGMtsPtr = *LGMtsPtrPtr;

    srt_free(FVArr); /* free unneeded structures */

    if (error)
    {
        LGMFreeSwptns(&ExerIntoPtr);
        if (skipCalib != 1)
        {
            if (CalReq->LGMOneTwoFactor == 1)
                LGMFreeLGM_TS(LGMtsPtrPtr);
            else if (CalReq->LGMOneTwoFactor == 2)
                LGMFreeLGM2F_TS(LGMtsPtrPtr);
            LGMFreeRefSwptnData(lgmRefDeals);
        }
        return (error);
    }

    /* test calibration against component prices (debug - remove for no debug) */
    /*	testcompprice(LGM1FtsPtr, tNow, ycName, (SrtSimMidAt*)dealPtr); */

    /* Step 4: Create the equivalent sigma_kappa term structure if requested */
    /* On error, returns a NULL pointer. We're treating this as non-fatal */
    if (skipEval != 1 && convertTS != 0)
    {
        if (CalReq->LGMOneTwoFactor == 1)
            SigKaptsPtr = LGMConvertZGtoSigKap(tNow, LGMtsPtr);
        else
        {
            SigKaptsPtr = LGMCreateSigKapTS(nEx, 1);

            SigKaptsPtr->numS     = nEx;
            SigKaptsPtr->numK     = 1;
            SigKaptsPtr->kdate[0] = LasttPay;
            SigKaptsPtr->kap[0]   = LGMtsPtr->one_kappa;

            one_kappa = LGMtsPtr->one_kappa;

            for (i = 1; i <= nEx; i++)
            {
                SigKaptsPtr->sdate[i - 1] = TauArr[i];
                if (i == 1)
                {
                    yr_to_exp  = (double)(TauArr[i] - tNow) * YEARS_IN_DAY;
                    kappa_term = (exp(2 * one_kappa * yr_to_exp) - 1) / (2 * one_kappa);
                    SigKaptsPtr->sig[i - 1] = sqrt(LGMtsPtr->Zeta1[0] / kappa_term);
                }
                else
                {
                    yr_to_exp      = (double)(TauArr[i] - tNow) * YEARS_IN_DAY;
                    prev_yr_to_exp = (double)(TauArr[i - 1] - tNow) * YEARS_IN_DAY;
                    kappa_term =
                        (exp(2 * one_kappa * yr_to_exp) - exp(2 * one_kappa * prev_yr_to_exp)) /
                        (2 * one_kappa);
                    SigKaptsPtr->sig[i - 1] =
                        sqrt(LGMtsPtr->Zeta1[i - 1] - LGMtsPtr->Zeta1[i - 2]) / kappa_term;
                }
            }
        }
    }

    srt_free(TauArr);
    /* Step 5: Evaluate deal, create xExBdry array and fill it with exer bdry */
    if (skipEval != 1 || findExerBdry != 0)
        error = NewMidAtEval(
            EvalParms, /* Convolution parameters */
            tNow,
            eod,
            ycName,
            CalReq,
            LGMtsPtr, /* Information about today */
            GetVol,
            GetBeta, /* swaption/cap vols */
            dealtype,
            dealPtr, /* the deal */
            &LGMVal, /* value of the deal */
            &nExBdry,
            &tExBdry,
            (findExerBdry ? &xExBdry : NULL)); /* array of exercise points (in x) */

    if (error)
    {
        LGMFreeSwptns(&ExerIntoPtr);
        if (skipCalib != 1)
        {
            if (CalReq->LGMOneTwoFactor == 1)
                LGMFreeLGM_TS(LGMtsPtrPtr);
            else if (CalReq->LGMOneTwoFactor == 2)
                LGMFreeLGM2F_TS(LGMtsPtrPtr);
            LGMFreeRefSwptnData(lgmRefDeals);
        }
        if (convertTS != 0)
        {
            if (CalReq->LGMOneTwoFactor == 1)
            {
                LGMFreeSigKapTS(&SigKaptsPtr);
                srt_free(xExBdry);
                srt_free(tExBdry);
            }
        }
        return (error);
    }

    /* Step 6: Find exer boundary & store as set of swaptions in ExerBdryPtr */
    /* On error, returns a NULL pointer. We're treating this as non-fatal */
    if (skipEval != 1 && findExerBdry != 0 && nExBdry > 0 && tExBdry != NULL && xExBdry != NULL)
    {
        ExerBdryPtr = FindBdry(
            CalReq->LGMOneTwoFactor, tNow, ycName, LGMtsPtr, nExBdry, tExBdry, xExBdry, LasttPay);
    }

    if (skipEval != 1 || findExerBdry != 0)
    {
        srt_free(xExBdry);
        srt_free(tExBdry); /* free unneeded structure */
    }
    /* Set outputs and return */
    if (skipEval != 1)
    {
        *LGMValPtr    = LGMVal;
        *intrinValPtr = intrinVal;
    }
    if (skipCalib != 1)
    {
        *ExerIntoPtrPtr = ExerIntoPtr;
        *LGMtsPtrPtr    = LGMtsPtr;
    }
    else
        LGMFreeSwptns(&ExerIntoPtr);
    if (skipEval != 1 && findExerBdry != 0)
        *ExerBdryPtrPtr = ExerBdryPtr;
    if (skipEval != 1 && convertTS != 0)
        *SigKaptsPtrPtr = SigKaptsPtr;

    return (error);
}

LGMErr TestLGMamortisingautocalCaller(
    long    nEx,         /* nEx is number of exercises */
    Date*   tEx,         /* notification (exercise) dates, [0,1,...,nEx-1] */
    Date*   tStart,      /* start dates for each exercise, [0,1,...,nEx-1] */
    double* Strike,      /* total value paid at tStart[j] for fixed leg, [0,1,...,nEx-1] */
    long    nPay,        /* nPay is number of fixed leg coupons */
    Date*   tPay,        /* pay dates for period i, [0,1, ...,nPay-1] */
    double* Payment,     /* total fixed leg payment(last includes notional), [0,...,nPay-1] */
    double* RedFirstPay, /* reduction in 1rst payment after exercise, [0,...,nEx-1] */
    char*   PayRecStr,   /* RECEIVER or PAYER */
    double  NumPeriod,
    double* index,
    double* d,
    double* dstart,
    double* mktprice,
    /* information about today */
    String ycName,   /* pointer to market structures */
    int    endofday, /* 1=too late to exercise deals today, 0=not too late */
                     /* today's volatilities */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* function to get swaption vols */
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
    double* intrinsicValPtr) /* ptr to intrinsic value (NULL => not req'd) */
{
    LGMDealType      ItsAMidAt;
    SrtSimMidAt*     MidAtPtr = NULL;
    SrtReceiverType  srt_p_r_flag;
    SrtDiffusionType srt_vol_type;
    LGMErr (*GetBeta)(Date, Date, double*);
    LGMErr (*RecVol)(Date, Date, double);
    LGMCalParm* CalReqPtr     = NULL;
    ConvParams* EvalParamsPtr = NULL;
    double      intrinsicVal  = 0;
    double      LGMVal        = 0;
    LGMSwptns*  ExerIntoPtr   = NULL;
    LGMSwptns*  ExerBdryPtr   = NULL;
    LGM_TS*     LGMtsPtr      = NULL;
    SigKapTS*   SigKaptsPtr   = NULL;
    Date        tNow, tfirst;
    SrtCurvePtr yldcrv = NULL;

    LGMErr error = NULL;
    long   skipCal;

    /*	1.- Initialise generic data
            ---------------------------	*/

    /* Free arrays that will be used for output, if already allocated */
    free_outputs(&convertTS, &findExerBdry, lgmExerBdryData, lgmRefSwptnData, lgmTSData, atcTSData);

    /* Set operations to be done by LGMnewautocal */
    skipCal = 0; /* always calibrate */

    /* Set information about today, and today's market */
    yldcrv = lookup_curve(ycName);
    tNow   = get_clcndate_from_yldcrv(yldcrv); /* get calculation date */

    /* Endofday and GetVol available from input */

    /* Get the exponent beta to use to interpret market vols */
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
    tfirst = tNow;
    if (endofday == 1)
    {
        tfirst = add_unit(tNow, 1, SRT_BDAY, SUCCEEDING);
    }

    /* Interpret pay-receive flag and basis */
    error = interp_rec_pay(PayRecStr, &srt_p_r_flag);
    if (error)
        return error;

    /* Set calibration method */
    CalReqPtr = LGMSetCalibMeth(
        LGMOneTwoFactor,
        usefixtau,
        usecaps,
        tau,
        alpha,
        gamma,
        rho,
        maxstd,
        Zeta1Dates,
        StartZeta1s,
        TauDates,
        StartTaus,
        HybridShortInstrsIndex);

    if (CalReqPtr == NULL)
    {
        return ("alloc failed in LGMcaller");
    }

    if (LGMOneTwoFactor == 2)
    {
        if (usefixtau == 0)
        {
            if ((calibrationmeth == 1) && (HybridShortInstrsIndex != NULL))
                CalReqPtr->calmeth = HybridShortAndDiag;
            else if ((calibrationmeth == 2) || (HybridShortInstrsIndex == NULL))
                CalReqPtr->calmeth = FullCapAndDiag;
            else
            {
                return serror(
                    "Calibration Method should be 1: HybridShortAndDiag or 2: FullCapAndDiag ");
                /* CalReqPtr->calmeth = TenorAndDiag; */
            }
        }
    }

    /* DEFINE INSTRUMENTS TO CALIBRATE ON */
    if (LGMOneTwoFactor == 1)
    {
        if (calibrationmeth == 1 && usefixtau != 1)
            /* CALIBRATION ON FIXED MATURITY SWAPTIONS + DIAG SWAPTIONS OR CAPLETS IF usecaps==1 */
            CalReqPtr->calmeth = FixExp;

        else if (calibrationmeth == 2 && usefixtau != 1)
            /* CALIBRATION ON LONG & SHORT SWAPTIONS  */
            CalReqPtr->calmeth = TenorAndDiag;

        else if (calibrationmeth == 3 && usefixtau != 1)
            CalReqPtr->calmeth = FixSigma;
    }

    /* DEFINE STRIKE METHOD */
    if (strikechoice == 1)
        CalReqPtr->Rmeth = IRR;
    else if (strikechoice == 2)
        CalReqPtr->Rmeth = dIRR;
    else if ((strikechoice == 3))
        CalReqPtr->Rmeth = EMK1;
    else if ((strikechoice == 4))
        CalReqPtr->Rmeth = EMK2;

    /* Set numerical parameters for convolution to their default values */
    EvalParamsPtr = LGMSetDefaultEvalParms();
    if (!EvalParamsPtr)
    {
        LGMFreeCalParm(&CalReqPtr);
        return "alloc failed in LGMcaller";
    }

    /* Record Vol */
    RecVol = LGMRecVolDummy;

    /*	2.- Initialise specific Midat data
            ---------------------------------- */

    /* Set deal type */
    ItsAMidAt = SimpleMidAt;

    /* allocate MidAtlantic stucture, and copy input data into it */
    /* This routine ASSUMES that the option holder will receive all payments
    whose pay dates tPay[i] are strictly after the settlement date. This shoud be
    replaced by user input */
    error = LGMFillSimMidAt(
        tfirst,
        &MidAtPtr,
        nPay,
        Payment,
        tPay,
        nEx,
        tEx,
        tStart,
        Strike,
        RedFirstPay,
        srt_p_r_flag);

    if (error != NULL)
    {
        LGMFreeCalParm(&CalReqPtr);
        srt_free(EvalParamsPtr);
        LGMFreeSimMidAt(&MidAtPtr);
        return (error);
    }

    /*	3.- Call to LGMAutocal
            ---------------------- */

    error = LGMamortisingnewautocal(
        skipEval,
        skipCal,
        convertTS,
        findExerBdry, /* task list */
        tNow,
        endofday,
        ycName, /* info about today's market */
        GetVol,
        GetBeta,
        RecVol, /* swaption & caplet price info */
        NumPeriod,
        index,
        d,
        dstart,
        mktprice,
        ItsAMidAt,
        (void*)MidAtPtr, /* the deal */
        CalReqPtr,
        EvalParamsPtr, /* calib & eval methods */
        &LGMVal,
        &intrinsicVal, /* value & intrinsic value of deal */
        &ExerIntoPtr,
        &ExerBdryPtr,    /* swaptions representing the underlying & exer boundary */
        lgmRefSwptnData, /* reference swaptions */
        &LGMtsPtr,
        &SigKaptsPtr, /* calibrated term structures */
        outfile);     /* unused */

    /*	4.- Free and return
            ------------------- */

    /* free unneeded structures */ /* these are not needed to get output */
    LGMFreeSimMidAt(&MidAtPtr);    /* structure containing the deal */
    LGMFreeCalParm(&CalReqPtr);    /* structure containing the requested cal method */
    srt_free(EvalParamsPtr);       /* structure containing the convolution parameters */
    LGMFreeSwptns(&ExerIntoPtr);   /* structure containing underlying swaptions of MidAt */

    /* keep track of Autocal Term Struct if required, otherwise free it */
    if (LGMOneTwoFactor == 1)
    {
        if (atcTSData)
        {
            *atcTSData = LGMtsPtr;
        }
        else
        {
            LGMFreeLGM_TS(&LGMtsPtr);
        }
    }
    else if (LGMOneTwoFactor == 2)
    {
        if (atcTSData)
        {
            *atcTSData = LGMtsPtr;
        }
        else if (LGMtsPtr)
            LGMFreeLGM2F_TS(&LGMtsPtr);
    }

    if (error != NULL)
    {
        LGMFreeSwptns(&ExerBdryPtr);
        LGMFreeSigKapTS(&SigKaptsPtr);
        if (LGMOneTwoFactor == 1 && atcTSData)
        {
            *atcTSData = NULL;
            LGMFreeLGM_TS(&LGMtsPtr);
        }
        return (error);
    }

    /* Unpack output */
    *LGMvalPtr = LGMVal;
    if (intrinsicValPtr != NULL)
        *intrinsicValPtr = intrinsicVal;

    /* If exercise boundary has been found, copy it into tExBdry & rExBdry */
    if (skipEval != 1 && findExerBdry != 0 && ExerBdryPtr != NULL && ExerBdryPtr->n > 0 &&
        lgmExerBdryData != NULL)
        error = copyExerBdry(ExerBdryPtr, lgmExerBdryData);
    if (error != NULL)
    {
        LGMFreeSwptns(&ExerBdryPtr); /* no longer needed */
        LGMFreeSigKapTS(&SigKaptsPtr);
        if (LGMOneTwoFactor == 1 && atcTSData)
        {
            *atcTSData = NULL;
            LGMFreeLGM_TS(&LGMtsPtr);
        }
        if (lgmExerBdryData != NULL)
        {
            lgmExerBdryData->NexerBdry = 0;
            if (lgmExerBdryData->exerBdryArr != NULL)
                srt_free(lgmExerBdryData->exerBdryArr);
        }
        return (error);
    }
    LGMFreeSwptns(&ExerBdryPtr); /* no longer needed */

    /* If sigma-kappa term structure has been found, copy it into output arrays */
    if (skipEval != 1 && convertTS != 0 && SigKaptsPtr != NULL && SigKaptsPtr->numS > 0 &&
        SigKaptsPtr->numK > 0 && SigKaptsPtr->sdate != NULL && SigKaptsPtr->sig != NULL &&
        SigKaptsPtr->kdate != NULL && SigKaptsPtr->kap != NULL && lgmTSData != NULL)
    {
        error = copyTSData(SigKaptsPtr, lgmTSData);
        if (error != NULL)
        {
            LGMFreeSigKapTS(&SigKaptsPtr);
            if (LGMOneTwoFactor == 1 && atcTSData)
            {
                *atcTSData = NULL;
                LGMFreeLGM_TS(&LGMtsPtr);
            }
            return (error);
        }
    }
    LGMFreeSigKapTS(&SigKaptsPtr); /* no longer needed */

    return error;
}

static LGMErr FitAllZetaAmortisingSwaptions(
    LGM_TS*    tsPtr,
    LGMCalSet* CSPtr,
    double     NumPeriod,
    double*    index,
    double*    d,
    double*    dstart,
    double*    mktprice)
{
    LGMErr          error;
    SrtReceiverType recpay;
    long            i, ifirst, j, n, nlong, nEx, k, ratio;
    double *        dispayment, *Gpay, Gst, noswaption;
    double*         coupon_G;
    double*         zeta;

    /* initialize */
    error    = NULL;
    nEx      = CSPtr->nEx;
    coupon_G = dvector(0, nEx);
    zeta     = dvector(0, nEx);
    if (tsPtr->numZ < nEx + 1)
    {
        srt_free(tsPtr->zeta);
        srt_free(tsPtr->zdate);
        tsPtr->zeta  = (double*)srt_calloc(nEx + 1, sizeof(double));
        tsPtr->zdate = (Date*)srt_calloc(nEx + 1, sizeof(Date));
        if (tsPtr->zeta == NULL || tsPtr->zdate == NULL)
            return ("allocation failed in FitAllZetaSwaps");
    }
    tsPtr->numZ = nEx + 1;

    for (k = tsPtr->numG - 1; k > 0; k--)
    {
        if (tsPtr->G[k - 1] < tsPtr->G[k])
            tsPtr->G[k - 1] = tsPtr->G[k];
    }

    dispayment = (double*)srt_calloc((unsigned int)NumPeriod, sizeof(double));
    Gpay       = (double*)srt_calloc(nEx, sizeof(double));

    if (dispayment == NULL || Gpay == NULL)
    {
        srt_free(dispayment);
        srt_free(Gpay);
        return ("allocation failed in FitAllZetaSwaps");
    }

    for (i = 1; i <= CSPtr->nPay; i++)
        Gpay[i] = LGMGFromTS(CSPtr->tPay[i], tsPtr); /* get G(t) at pay dates */

    if (dstart[0] == 0)
    {
        recpay = SRT_PAYER;
    }
    else
    {
        recpay = SRT_RECEIVER;
    }

    tsPtr->zdate[0] = CSPtr->tNow;
    tsPtr->zeta[0]  = 0.0;
    noswaption      = index[nEx];

    /*Find the ratio (call frequency/fixed frequency) */
    ratio = 0;

    while (index[ratio] == index[0])
    {
        ratio++;
    }

    /* For each j, find the zeta(tEx[j]) that fits the LGM value of
long swaption j to the market price Vlong[j] */
    for (j = 1; j <= nEx; j++)
    {
        ifirst = (long)index[j - 1];
        nlong  = (long)NumPeriod;

        for (i = ifirst + 2; i <= nlong; i++)
        {
            dispayment[i - ifirst - 1] = d[i];
        }
        n = nlong - ifirst - 1;
        for (i = 1; i < n + 1; i++)
        {
            coupon_G[i] = tsPtr->G[j - 1 + i * ratio];
        }

        Gst = tsPtr->G[j - 1];

        if ((j == 1) || ((index[j - 1] > index[j - 2]) && (j <= nEx - noswaption + 1)))
        {
            zeta[j] = FitOneZetaSwap(
                n, dispayment, dstart[ifirst + 2], recpay, coupon_G, Gst, mktprice[ifirst + 1]);
        }
        else
        {
            zeta[j] = zeta[j - 1];
        }

        if (zeta[j] < 0)
        {
            srt_free(Gpay);
            srt_free(dispayment);
            return ("failed to fit zeta for swaption");
        }
        tsPtr->zeta[j]  = zeta[j]; /* update term structure */
        tsPtr->zdate[j] = CSPtr->tEx[j];
    }
    srt_free(dispayment);
    srt_free(Gpay);
    return (error);
}

LGMErr LGMCalFixKapamortising(
    LGM_TS**   LGMtsPtrPtr, /* Return: Calibrated term structure */
    double     kap,         /* Fixed kappa to use to construct G */
    int        usestarts,   /* construct G on start dates or pay dates */
    long       ndate,       /* number of dates for G if usestarts=0 */
    Date*      Gdate,       /* dates for G if usestarts=0 */
    int        usecaps,     /* calibrate on caplets or swaptions */
    LGMCalSet* CSPtr,
    double     NumPeriod,
    double*    index,
    double*    d,
    double*    dstart,
    double*    mktprice) /* Reference caplets */
{
    LGMErr  error;
    LGM_TS* tsPtr;

    /* initialize */
    error        = NULL;
    *LGMtsPtrPtr = NULL;

    /* Create term structure, fill in  G dates, calculate and fill in G values from kappa */
    tsPtr = CreateTSwithKappa(kap, usestarts, ndate, Gdate, usecaps, CSPtr);
    if (tsPtr == NULL)
        return ("error in CreateTSwithKappa");

    /* Pick the zeta's in the term structure *ts by calibrating to the caplets or swaptions
            in the calibration set CS */
    if (usecaps == 1)
    {
        error = FitAllZetaCaps(tsPtr, CSPtr);
    }
    else
    {
        error = FitAllZetaAmortisingSwaptions(tsPtr, CSPtr, NumPeriod, index, d, dstart, mktprice);
    }

    if (error != NULL)
    {
        LGMFreeLGM_TS(&tsPtr);
        return (error);
    }

    /* Normalize the calibrated term structure */
    if (usecaps == 0)
    {
        error = LGMVerifyLGM_TS_For_DiagSwpts(CSPtr, tsPtr);

        /*error = LGMVerifyLGM_TS(tsPtr);*/
        if (error != NULL)
            LGMFreeLGM_TS(&tsPtr);
    }

    *LGMtsPtrPtr = tsPtr;
    return (error);
}

LGMErr LGMamortisingCalibration(
    /* info about today's market */
    Date   tNow,   /*  evaluation date */
    String ycname, /* market pointer for discount factors */
    LGMErr (*GetVol)(
        Date, Date, double, SRT_Boolean, double*), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date, double*), /* function to provide exponent (beta) to interpret vols */
    LGMErr (*RecVol)(Date, Date, double), /* function called whenever LGM uses a vol for pricing */
    double  NumPeriod,
    double* index,
    double* d,
    double* dstart,
    double* mktprice,
    /* info about the deal */
    LGMDealType dealType,
    void*       dealptr,
    long        nEx, /* number of exercise dates to calibrate to */
    Date        tlast,
    Date*       TauArr,        /* calibrate on these exercise dates */
    double*     FVArr,         /* use swaptions with these forward vals of fixed legs */
                               /* calibration methodology */
    LGMCalParm* CalReq,        /* structure defining calibration method to be used */
                               /* calibrated term structure */
    LGM_TS**    LGMtsPtrPtr,   /* CALIBRATED ZETA-G  TS*/
    LGMCalSet** RefDealPtrPtr, /* Reference instruments used for calibration */
    SrtLgmRefSwptnData*
        lgmRefSwptnData) /* ptr to reference swaption data structure (NULL => not req'd) */
{
    LGMErr       error;
    LGM_TS*      LGMtsPtr = NULL;
    LGMCalSet*   RefDealPtr;
    SrtSimMidAt* MidAt;
    LGMMarkConv  conv;
    SrtCrvPtr    yldcrv;
    long         noPairsFlag = 0;

    /* initialize */
    error          = NULL;
    *RefDealPtrPtr = NULL;
    *LGMtsPtrPtr   = NULL;
    yldcrv         = lookup_curve(ycname);
    error          = LGMCcyDefaults(yldcrv, &conv); /* get currency defaults */
    if (error != NULL)
        return (error);

    /* Check calibration request for sensibility data */

    checkcalibrationmethod(CalReq);

    if (CalReq->LGMOneTwoFactor == 1)
    {
        /* Step 1: Construct reference instruments */
        error = LGMMakeRefDeals(
            dealType,
            dealptr,
            tNow,
            tlast,
            nEx,
            TauArr,
            FVArr, /* exer dates and relative PV of fixed leg of deal */
            ycname,
            &conv,  /* discount curve & market conventions */
            CalReq, /* method and data for choosing strikes */
            GetVol,
            GetBeta,      /* information for choosing strikes */
            &RefDealPtr); /* reference instruments */
        if (error != NULL)
        {
            LGMFreeCalSet(&RefDealPtr);
            return (error);
        }

        /* Step 2: Price reference instruments, and record swaptions used */
        error = LGMPriceRefDeals(GetVol, GetBeta, RecVol, RefDealPtr, lgmRefSwptnData);
        if (error != NULL)
        {
            LGMFreeCalSet(&RefDealPtr);
            return (error);
        }
    }

    /*	Step 3: Calibrate */
    if (CalReq->LGMOneTwoFactor == 1)
    {
        switch (CalReq->calmeth)
        {
        /* Use input kappa to construct G(t); then get zeta(t) by calibrating on
        the long "k into n-k" swaptions or the caplets, depending on usecaps */
        case FixKappa:
            error = LGMCalFixKapamortising(
                &LGMtsPtr,         /* Return: Calibrated term structure */
                CalReq->kap,       /* Fixed kappa to use to construct G */
                CalReq->usestarts, /* construct G on start dates or pay dates */
                RefDealPtr->n,     /* number of dates for G if usestarts=0 */
                RefDealPtr->tPay,  /* dates for G if usestarts=0 */
                CalReq->usecaps,   /* calibrate on caplets or swaptions */
                RefDealPtr,
                NumPeriod,
                index,
                d,
                dstart,
                mktprice);

            /* Reference instruments */
            break;

        /* Use input G(t) for G(t); then get zeta(t) by calibrating on the long
        "k into n-k" swaptions or the caplets, depending on usecaps */
        case GivenG:
            error = LGMCalGivenG(
                &LGMtsPtr,       /* Return: Calibrated term structure */
                CalReq->numG,    /* number of G(t) values */
                CalReq->Gdate,   /* dates t for G(t) */
                CalReq->G,       /* G(t) at these dates */
                CalReq->usecaps, /* calibrate on caplets or swaptions */
                RefDealPtr);     /* Reference instruments */
            break;

        /* Simultaneously calibrate on the long "k into n-k" swaptions and either the
        short "k into 1" swaptions or caplets depending on usecaps */
        /* This is the calibration used by the old LGM autocal */
        case TenorAndDiag:
            noPairsFlag = 0;
            error       = LGMCalTenorDiag(
                &LGMtsPtr,       /* Return: Calibrated term structure */
                CalReq->usecaps, /* calibrate on caplets or swaptions */
                RefDealPtr,      /* Reference instruments */
                &noPairsFlag);   /* signals if there are no pairs to calibrate on */
            if (noPairsFlag != 1)      /* if successful, continue */
                break;
            else                     /* else try the the Fixed sigma calibration */
                CalReq->usecaps = 0; /* calibrating on the diagonal */
        /* Use constant sigma to construct zeta(t), then get G(t) by calibrating on
        the long "k into n-k" swaptions or the caplets, depending on usecaps */
        case FixSigma:
            error = LGMCalFixSig(
                &LGMtsPtr,       /* Return: Calibrated term structure */
                CalReq->usecaps, /* calibrate on caplets or swaptions */
                RefDealPtr);     /* Reference instruments */
            break;

        /* Use input zeta(t) for zeta(t), then get G(t) by calibrating on
        the long "k into n-k" swaptions or the caplets, depending on usecaps */
        case GivenZeta:
            error = LGMCalGivenZeta(
                &LGMtsPtr,       /* Return: Calibrated term structure */
                CalReq->numZ,    /* number of zeta(t) values */
                CalReq->zdate,   /* dates t for zeta(t) */
                CalReq->zeta,    /* zeta(t) at these dates */
                CalReq->usecaps, /* calibrate on caplets or swaptions */
                RefDealPtr);     /* Reference instruments */
            break;

        /* Calibrate on the "1 into k" swaptions to get G(t), then get zeta(t) by calibrating
        on the long "k into n-k" swaptions or the caplets, depending on usecaps */
        case FixExp:
            error = LGMCalFixExp(
                &LGMtsPtr,          /* Return: Calibrated term structure */
                CalReq->usecaps,    /* calibrate on caplets or swaptions */
                CalReq->keep1intok, /* keep 1intok zeta */
                RefDealPtr);        /* Reference instruments */
            break;
        } /* end switch */

        if (error != NULL)
        {
            if (CalReq->LGMOneTwoFactor == 1)
                LGMFreeLGM_TS(&LGMtsPtr);
            if (CalReq->LGMOneTwoFactor == 2)
                LGMFreeLGM2F_TS(&LGMtsPtr);
            LGMFreeCalSet(&RefDealPtr);
        }

        *RefDealPtrPtr = RefDealPtr;
        *LGMtsPtrPtr   = LGMtsPtr;
    }
    else if (CalReq->LGMOneTwoFactor == 2)
    {
        MidAt        = (SrtSimMidAt*)dealptr;
        MidAt->tLast = tlast;

        error = lgm_2f_autocal_calibrate(
            dealptr, LGMtsPtrPtr, lgmRefSwptnData, CalReq, GetVol, GetBeta, ycname);
        if (error != NULL)
        {
            if (LGMtsPtr)
                LGMFreeLGM2F_TS(&LGMtsPtr);
            LGMFreeCalSet(&RefDealPtr);
        }
    }
    return (error);
}

static LGMErr FitAllZetaCaps(
    LGM_TS*    tsPtr, /* Return: zeta's in LGM term structure ts */
    LGMCalSet* CSPtr) /* Reference caplets */
{
    LGMErr          error;
    SrtReceiverType recpay;
    long            j, k, nEx;
    double          Gstart, Gend, deltaG, stddev, zeta;

    /* initialize */
    error = NULL;
    nEx   = CSPtr->nEx;

    if (tsPtr->numZ < nEx + 1)
    {
        srt_free(tsPtr->zeta);
        srt_free(tsPtr->zdate);
        tsPtr->zeta  = (double*)srt_calloc(nEx + 1, sizeof(double));
        tsPtr->zdate = (Date*)srt_calloc(nEx + 1, sizeof(Date));
        if (tsPtr->zeta == NULL || tsPtr->zdate == NULL)
            return ("allocation failed in FitAllZetaCaps");
    }
    tsPtr->numZ = nEx + 1;

    for (k = tsPtr->numG - 1; k > 0; k--)
    {
        if (tsPtr->G[k - 1] < tsPtr->G[k])
            tsPtr->G[k - 1] = tsPtr->G[k];
    }

    tsPtr->zdate[0] = CSPtr->tNow;
    tsPtr->zeta[0]  = 0.0;

    recpay = SRT_RECEIVER;
    for (j = 1; j <= nEx; j++)
    {
        stddev = FitOneCap(
            CSPtr->Vcap[j],
            CSPtr->cvgcap[j],
            CSPtr->Rfcap[j],
            CSPtr->DStart[j],
            CSPtr->Dcap[j],
            recpay);
        if (stddev < 0.0) /* error */
            return ("failed to fit cap");

        Gstart = LGMGFromTS(CSPtr->tStart[j], tsPtr);
        Gend   = LGMGFromTS(CSPtr->tEnd[j], tsPtr);
        deltaG = Gstart - Gend;
        zeta   = (stddev / deltaG) * (stddev / deltaG);

        tsPtr->zeta[j]  = zeta;
        tsPtr->zdate[j] = CSPtr->tEx[j];
    }
    return (error);
}

static LGM_TSPtr CreateTSwithKappa(
    double     kap,       /* Fixed kappa to use to construct G */
    int        usestarts, /* construct G on start dates or pay dates */
    long       ndate,     /* number of dates for G if usestarts=0 */
    Date*      Gdate,     /* dates for G if usestarts=0 */
    int        usecaps,   /* final date from caplet or swaption */
    LGMCalSet* CSPtr)     /* Reference instruments */
{
    LGM_TS* tsPtr;
    long    i, k, n, nPay, nEx, numG;
    Date    tLast;
    double  dt, dtn, arg, argn, factor;

    nEx = CSPtr->nEx;

    /* Create term structure, fill in numG and Gdates */
    if (usestarts != 0)
    {
        numG  = nEx + 1;
        tsPtr = LGMCreateLGM_TS(nEx + 1, numG);
        if (tsPtr == NULL)
            return (NULL);

        tsPtr->numG = numG;
        for (i = 0; i < nEx; i++)
            tsPtr->Gdate[i] = CSPtr->tStart[i + 1];

        if (usecaps == 1) /* get last date */
            tsPtr->Gdate[nEx] = CSPtr->tEnd[nEx];
        else
        {
            n    = CSPtr->n;
            nPay = CSPtr->nPay;
            for (i = n; i <= nPay && CSPtr->tPay[i] <= tsPtr->Gdate[nEx - 1]; i++)
                ;
            tsPtr->Gdate[nEx] = CSPtr->tPay[i];
            if (tsPtr->Gdate[nEx] <= tsPtr->Gdate[nEx - 1])
                tsPtr->numG = numG = nEx;
        }
    }
    else /* use input G dates */
    {
        if (usecaps == 1) /* find last date to use */
            tLast = CSPtr->tEnd[nEx];
        else
            tLast = CSPtr->tPay[CSPtr->n];

        for (i = 1; i <= ndate && Gdate[i] <= CSPtr->tStart[1];
             i++) /* find number of dates to use */
            ;
        for (n = ndate; n >= 1 && Gdate[n] >= tLast; n--)
            ;
        numG = 2;
        if (n >= i)
            numG = 3 + n - i;

        tsPtr = LGMCreateLGM_TS(nEx + 1, numG); /* create term structure */
        if (tsPtr == NULL)
            return (NULL);

        tsPtr->numG            = numG; /* fill in G dates */
        tsPtr->Gdate[0]        = CSPtr->tStart[1];
        tsPtr->Gdate[numG - 1] = tLast;
        for (k = i; k <= n; k++)
            tsPtr->Gdate[k + 1 - i] = Gdate[k];
    }

    /* fill in G values */
    dtn  = ((double)(tsPtr->Gdate[numG - 1] - tsPtr->Gdate[0])) / 365.0;
    argn = kap * dtn;
    if (fabs(argn) < 1.0e-6)
    {
        for (i = 0; i < numG; i++)
        {
            dtn         = ((double)(tsPtr->Gdate[numG - 1] - tsPtr->Gdate[i])) / 365.0;
            dt          = ((double)(tsPtr->Gdate[i] - tsPtr->Gdate[0])) / 365.0;
            tsPtr->G[i] = dtn * (1 - 0.5 * kap * dt);
        }
    }
    else
    {
        factor = dtn / (1 - exp(-argn));
        for (i = 0; i < numG; i++)
        {
            arg         = kap * ((double)(tsPtr->Gdate[i] - tsPtr->Gdate[0])) / 365.0;
            tsPtr->G[i] = (exp(-arg) - exp(-argn)) * factor;
        }
    }
    return (tsPtr);
}

static double FitOneZetaSwap(
    long            n,
    double*         a,
    double          DStart,
    SrtReceiverType recpay,
    double*         Gpay,
    double          Gst,
    double          price)
{
    SrtCallPutType callput;
    double         fwd, sens, ystar, favg, stddev;
    double         value, newvalue, der, zeta, sqzeta, dsqzeta;
    long           i, iter, itermax, success;

    /* check for intrinsic */
    if (fabs(LGMRecVal(&ystar, n, a, DStart, 1.0e-15, Gpay, Gst) - price) < 1.0e-08)
    {
        return 1.0e-15;
    }

    /* initialize */
    itermax = 25;
    callput = SRT_CALL;

    /* compute forward value and sensitivity of swaption */
    fwd  = 0.;
    sens = 0.0;
    for (i = 1; i <= n; i++)
    {
        fwd  = fwd + a[i];
        sens = sens + a[i] * (Gst - Gpay[i]);
    }

    /* compute better sensitivity */
    if (fabs((fwd - DStart) / sens) > 1.e-5)
    {
        ystar = LGMFindystar(n, a, DStart, 0., Gpay, Gst);
        sens  = -(fwd - DStart) / ystar;
    }

    /* need only consider receivers now */
    if (recpay == SRT_PAYER)
        price = price + fwd - DStart;

    /* get initial guess by linear approximation (neglects convexity) */
    favg   = 0.5 * (fwd + DStart);
    stddev = LGMNormImpVol(price / favg, fwd / favg, DStart / favg, 1.0, 1.0, callput);

    /* Use global Newton on sqzeta = sqrt(zeta(tEx)) to match LGM price to market price */
    success = -2; /* set up for global newton */
    sqzeta  = 0.;
    dsqzeta = stddev * favg / fabs(sens); /* initial guess */
    value   = -price;

    for (iter = 1; iter <= itermax && success <= 0; iter++)
    {
        newvalue = LGMRecVal(&ystar, n, a, DStart, sqzeta + dsqzeta, Gpay, Gst) - price;
        if (fabs(newvalue) / price < 1.e-7 ||
            fabs(dsqzeta) < 1.e-7) /* take three steps after OK convergence */
            success++;
        if (iter == 1 || fabs(newvalue) < fabs(value)) /* step successful */
        {
            sqzeta  = sqzeta + dsqzeta; /* take new step */
            value   = newvalue;
            der     = RecSqZetaDer(ystar, n, a, DStart, sqzeta, Gpay, Gst);
            dsqzeta = -value / der; /* Newton step */
        }
        else
            dsqzeta = 0.5 * dsqzeta;
    }
    sqzeta = sqzeta + dsqzeta;
    zeta   = sqzeta * sqzeta;
    if (success >= 0) /* convergence */
        return (zeta);
    else
        return (-1.); /* failure */
}

static double FitOneCap(
    double price, double cvg, double Rfix, double DStart, double Dend, SrtReceiverType recpay)
{
    double         fwd, stddev;
    SrtCallPutType callput;

    if (recpay == SRT_RECEIVER)
        callput = SRT_CALL;
    else
        callput = SRT_PUT;

    fwd    = (1.0 + cvg * Rfix) * Dend;
    stddev = LGMBlackImpVol(price, fwd, DStart, 1.0, 1.0, callput);
    return (stddev);
}

static double RecSqZetaDer(
    double  ystar,
    long    n,
    double* a,
    double  DStart, /* receiver swaption */
    double  sqzeta,
    double* Gpay,
    double  Gst) /* term structure */
{
    double der, arg, deltaG;
    long   i;

    if (sqzeta <= 1.e-12)
        sqzeta = 1.e-12;

    arg = ystar / sqzeta;
    der = 0.0;
    for (i = 1; i <= n; i++)
    {
        deltaG = Gst - Gpay[i];
        der    = der + deltaG * a[i] * LGMsafeGauss(arg + deltaG * sqzeta);
    }

    if (Gst < Gpay[n])
        der = -der; /* this should never occur */

    return (der);
}