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
#include "LGMCalib.h"

#include "LGM2FRefInstrs.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"

/****************************************************************************/
/* Prototypes of private functions */
/****************************************************************************/
/* Check calibration methods and data in the  request CR */
static void checkcalibrationmethod(LGMCalParm* CRptr);

/* After the G(t) values have been filled in the term structure, this
routine finds the zeta values by calibrating to the caplets in CSPtr */
static LGMErr FitAllZetaCaps(
    LGM_TS*    tsPtr,  /* Return: zeta's in LGM term structure ts */
    LGMCalSet* CSPtr); /* Reference caplets */

/* After the zeta values have been filled in the term structure, this
routine finds the G values by calibrating to the caplets */
static LGMErr FitAllGCaps(
    LGM_TS*    tsPtr,  /* Return: calibrated LGM term structure tsPtr */
    LGMCalSet* CSPtr); /* Reference caplets */

/* Fits the LGM price of a caplet to the market price */
static double FitOneCap(
    double price, double cvg, double Rfix, double DStart, double Dend, SrtReceiverType recpay);

/* After the G values have been filled in the term structure, this
routine finds the zeta values by calibrating to the long swaptions */
static LGMErr FitAllZetaSwaps(
    LGM_TS*    tsPtr,  /* Return: zeta's in LGM term structure ts */
    LGMCalSet* CSPtr); /* Reference caplets */

/* Fits value of zeta at exercise date for a swaption */
static double FitOneZetaSwap(
    long            n,
    double*         a,
    double          DStart,
    SrtReceiverType recpay,
    double*         Gpay,
    double          Gst,
    double          price);

/* After the zeta values have been filled in the term structure, this
routine finds the G(t) values by calibrating to the long swaptions */
static LGMErr FitAllGBkwdSwap(LGM_TS* tsPtr, LGMCalSet* CSPtr);

/* Determine swaption/caplet pairs to use to find G(t) */
static long findswapcappairs(long* jpair, LGMCalSet* CSPtr);

/* Determine long/short swaption pairs to use to find G(t) */
static long findswappairs(long* jpair, LGMCalSet* CSPtr);

/* Determine G(t) so that the LGM value of long swaption j and either
the LGM value of the caplet or short swaption match their market prices */
static LGMErr FitAllGPairs(long npair, long* jpair, long usecaps, LGM_TS* tsPtr, LGMCalSet* CSPtr);

/* Picks coef so the LGM price of a long swaption matches the market price when
        Gst = Gst + coef*dGst							(G at start date)
        Gpay[i] = Gpay[i] + coef*dGpay[i], i=1, ...,n	(G at pay dates) */
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
    double          price);         /* market price of swaption */

/* Picks coef so the LGM price of a long swaption matches the market price when
        Gst = Gst + coef*dGst							(G at start date)
        Gpay[i] = Gpay[i] + coef*dGpay[i], i=1, ...,n	(G at pay dates)
and sqrt(zeta)*Gst is held constant at product */
static double FitOneSwapCap(
    long            n,
    double*         a,
    double          DStart,
    SrtReceiverType recpay, /* the swaption */
    double          zeta,
    double          product, /* zeta, product */
    double*         Gstptr,
    double*         Gpay, /* the existing term structure */
    double          dGst,
    double*         dGpay, /* direction of change */
    double          price);         /* market price of swaption */

/* Picks coef and sqzeta so the LGM prices of the long and short
swaptions match their market price when
        Gst = Gst + coef*dGst							(G at start date)
        Gpay[i] = Gpay[i] + coef*dGpay[i], i=1, ...,n	(G at pay dates)
        zeta(tEx) = sqzeta*sqzeta */
static double FitOneLongShort(
    long            n1,
    double*         a1, /* swaption 1 */
    long            n2,
    double*         a2, /* swaption 2 */
    double          DStart,
    SrtReceiverType recpay, /* PV of strike for both swaptions */
    double*         Gstptr,
    double*         Gpay, /* the existing term structure */
    double          dGst,
    double*         dGpay, /* direction of change */
    double          price1,
    double          price2); /* market prices of swaptions */

/* Derivative of LGM value with resp. to sqrt(zeta) for a receiver swaption */
/* NOTE: LGMRecVal must be run first so that we have ystar */
static double RecSqZetaDer(
    double  ystar,
    long    n,
    double* a,
    double  DStart, /* receiver swaption */
    double  sqzeta,
    double* Gpay,
    double  Gst); /* term structure */

/* Derivative of  LGM value with respect to C where
Gst -> Gst + C*dGst,   Gpay[i] -> Gpay[i] + CdGpay[i]
for a receiver swaption */
/* NOTE: LGMRecVal must be ran first so that we have ystar */
static double RecGDer(
    double  ystar,
    long    n,
    double* a,
    double  DStart, /* receiver swaption */
    double  sqzeta,
    double* Gpay,
    double  Gst, /* term structure */
    double  dGst,
    double* dGpay); /* direction of differentiation */

/* Create a term structure and calculate G(t) from fixed kappa */
static LGM_TSPtr CreateTSwithKappa(
    double     kap,       /* Fixed kappa to use to construct G */
    int        usestarts, /* construct G on start dates or pay dates */
    long       ndate,     /* number of dates for G if usestarts=0 */
    Date*      Gdate,     /* dates for G if usestarts=0 */
    int        usecaps,   /* final date from caplet or swaption */
    LGMCalSet* CSPtr);    /* Reference instruments */

/* Create a term structure and fills G(t) from input arrays */
static LGM_TSPtr CreateTSwithGs(
    long       numG,   /* number of values of G */
    Date*      Gdate,  /* construct G on these dates */
    double*    G,      /* values of G(t) */
    LGMCalSet* CSPtr); /* Reference instruments */

/* Create a term structure and compute zeta(t) from constant sigma */
static LGM_TSPtr CreateTSwithSigma(LGMCalSet* CSPtr); /* Reference instruments */

/* Create a term structure and fills in the input zeta(t) */
static LGM_TSPtr CreateTSfromZetas(
    long       nZ,     /* number of zeta values */
    Date*      Zdate,  /* zeta given on these dates */
    double*    Zval,   /* zeta values */
    LGMCalSet* CSPtr); /* Reference instruments */

/****************************************************************************/
/**** Comments ****/
/**** HOLIDAYS ***/
/* Wherever we need to use the add_unit function,
the currency code ccy is available */

/****************************************************************************/
/* Main LGM calibration routine */
LGMErr LGMCalibration(
    /* info about today's market */
    Date   tNow,   /*  evaluation date */
    String ycname, /* market pointer for discount factors */
    LGMErr (*GetVol)(
        Date, Date, double, SRT_Boolean, double*), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date, double*), /* function to provide exponent (beta) to interpret vols */
    LGMErr (*RecVol)(Date, Date, double), /* function called whenever LGM uses a vol for pricing */
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
            error = LGMCalFixKap(
                &LGMtsPtr,         /* Return: Calibrated term structure */
                CalReq->kap,       /* Fixed kappa to use to construct G */
                CalReq->usestarts, /* construct G on start dates or pay dates */
                RefDealPtr->n,     /* number of dates for G if usestarts=0 */
                RefDealPtr->tPay,  /* dates for G if usestarts=0 */
                CalReq->usecaps,   /* calibrate on caplets or swaptions */
                RefDealPtr);       /* Reference instruments */
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

/***************************************************/
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

/**************************************************************/
/******************** calibration routines ********************/
/**************************************************************/
/* The following routine first creates a term structure  and constructs G(t) from
the incoming fixed kappa:
        If usestarts != 0, it will construct G(t) on each start date of the caplets or
swaptions, and on the end date of the last caplet (if usecaps=1) or tPay[n] (if usecaps!=1).
        If usestarts = 0, it will construct G(t) on the start date of the first caplet
or swaption and on the end date of the last caplet or tPay[n] , and on the dates
Gdates[j], j=1,..., which fall between the first start date and the last end date.
Note that Gdate[0] is ignored.  In the call from LGMCalibration, Gdate[] are the pay dates
 of the fixed leg of the longest swaption used to calibrate
        It then finds zeta(t) at tNow, and at each exercise date of the set of caplets or
swaptions by calibrating to the set of caplets  (if usecaps=1) or long swaptions
(if usecaps!=1) in the calibration set CSPtr. */

/* If usecaps == 1, the following information must be set in CSPtr:
tNow							evaluation date
nEx>=1							number of reference caplets
tEx[j] for j=1,...,nEx			exercise date for caplet j
tStart[j] for j=1,...,nEx		start date for caplet j
DStart[j] for j=1,...,nEx		dis factor from tNow to the start date
tEnd[j] for j=1,...,nEx			end date for caplet j
Dcap[j] for j=1,...,nEx			dis factor to the end date of caplet j
cvgcap[j] for j=1,...,nEx		coverage from start date to end date for caplet j
Rfcap[j] for j=1,...,nEx		fixed rate for caplet j
Vcap[j]	for j=1,...,nEx			market value at tNow for caplet j
        Note the index j=0 is ignored in the above arrays. */

/* If usecaps != 1, the following information must be set in CSPtr:
tNow							evaluation date
nEx>=1							number of reference swaptions
tEx[j] for j=1,...,nEx			exercise date for swaption j
tStart[j] for j=1,...,nEx		start date for swaption j
DStart[j] for j=1,...,nEx		dis factor from tNow to the start date
nPay							total number of fixed leg pay dates
tPay[i] for i=1,...,nPay		fixed leg pay dates
cvgpay[i] for i=1,...,nPay		coverage from tPay[i-1] to tPay[i]
Dpay[i] for i=1,...,nPay		dis factor from tNow to tPay[i]
ifirst [j] for j=1,...,nEx		fixed leg pay dates are i = ifirst[j] to nlong[j]
nlong [j] for j=1,...,nEx
cvgfirst t[j] for j=1,...,nEx	first period coverage (from tStart[j] to tPay[ifirst[j]])
Rlong [j] for j=1,...,nEx		fixed rate for swaption j
Vlong [j] for j=1,...,nEx		market price of swaption j
        Note the indices j=0 and i=0 are ignored in the above arrays */
LGMErr LGMCalFixKap(
    LGM_TS**   LGMtsPtrPtr, /* Return: Calibrated term structure */
    double     kap,         /* Fixed kappa to use to construct G */
    int        usestarts,   /* construct G on start dates or pay dates */
    long       ndate,       /* number of dates for G if usestarts=0 */
    Date*      Gdate,       /* dates for G if usestarts=0 */
    int        usecaps,     /* calibrate on caplets or swaptions */
    LGMCalSet* CSPtr)       /* Reference caplets */
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
        error = FitAllZetaCaps(tsPtr, CSPtr);
    else
        error = FitAllZetaSwaps(tsPtr, CSPtr);

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

/*************************************************************************/
/* This routine is the same as above, except that it uses the input arrays
        G(t) = G[i] at t=Gdate[i], i=1,...,ndate
to define G(t) in the term structure */
LGMErr LGMCalGivenG(
    LGM_TS**   LGMtsPtrPtr, /* Return: Calibrated term structure */
    long       ndate,       /* number of dates for G */
    Date*      Gdate,       /* dates for G(t) */
    double*    G,           /* values for G(t) */
    int        usecaps,     /* calibrate on caplets or swaptions */
    LGMCalSet* CSPtr)       /* Reference instruments */
{
    LGMErr  error;
    LGM_TS* tsPtr;

    /* initialize */
    error        = NULL;
    *LGMtsPtrPtr = NULL;

    /* Create term structure, fill in  G dates, calculate and fill in G values from kappa */
    tsPtr = CreateTSwithGs(ndate, Gdate, G, CSPtr);
    if (tsPtr == NULL)
        return ("error in CreateTSwithG");

    /* Pick the zeta's in the term structure *ts by calibrating to the caplets in CS */
    if (usecaps == 1)
        error = FitAllZetaCaps(tsPtr, CSPtr);
    else
        error = FitAllZetaSwaps(tsPtr, CSPtr);

    if (error != NULL)
    {
        LGMFreeLGM_TS(&tsPtr);
        return (error);
    }

    /* Normalize the calibrated term structure */
    error = LGMVerifyLGM_TS_For_DiagSwpts(CSPtr, tsPtr);
    if (error != NULL)
        LGMFreeLGM_TS(&tsPtr);

    *LGMtsPtrPtr = tsPtr;
    return (error);
}

/**************************************************/
/* Fixed sigma and given zeta calibration		  */
/**************************************************/
/* This routine creates a term structure and calibrates it to the set of caplets or the
set of long (k into n-k) swaptions calibration set CSPtr. It first constructs zeta(t)
using a fixed value of sigma: zeta[j] = sigma*sigma*(tEx[j]-tNow) on the exercise
dates tEx[j], j=1,...,nEx. The particular value of sigma used is irrelevent.
/* If usecaps ==1, the following information must be set in CSPtr:
tNow						evaluation date
nEx>=1						number of reference caplets
tEx[j] for j=1,...,nEx		exercise date for caplet j
tStart[j] for j=1,...,nEx		start date for caplet j
DStart[j] for j=1,...,nEx		dis factor from tNow to the start date
tEnd[j] for j=1,...,nEx		end date for caplet j
Dcap[j] for j=1,...,nEx		dis factor to the end date of caplet j
cvgcap[j] for j=1,...,nEx	coverage from start date to end date for caplet j
Rfcap[j] for j=1,...,nEx	fixed rate for caplet j
Vcap[j]	for j=1,...,nEx		market value at tNow for caplet j
        Note the index j=0 is ignored in the above arrays */

/* If usecaps !=1, the following information must be set in CSPtr:
tNow						evaluation date
nEx>=1						number of reference swaptions
tEx[j] for j=1,...,nEx		exercise date for swaption j
tStart[j] for j=1,...,nEx		start date for swaption j
DStart[j] for j=1,...,nEx		dis factor from tNow to the start date
n							tPay[n] is the last date we need to
calibrate to tPay[i] for i=1,...,nPay	fixed leg pay dates cvgpay[i] for i=1,...,nPay	coverage
from tPay[i-1] to tPay[i] Dpay[i] for i=1,...,nPay	dis factor from tNow to tPay[i] ifirst [j]
for j=1,...,nEx	fixed leg pay dates are i = ifirst[j] to nlong[j] nlong [j] for j=1,...,nEx cvgfirst
t[j] for j=1,...,nEx  first period coverage (from tStart[j] to tPay[ifirst[j]]) Rlong [j] for
j=1,...,nEx	fixed rate for swaption j Vlong [j] for j=1,...,nEx	market price of swaption j
        Note the indices j=0 and i=0 are ignored in the above arrays */
LGMErr LGMCalFixSig(
    LGM_TS**   LGMtsPtrPtr, /* Return: Calibrated term structure */
    int        usecaps,     /* calibrate on caplets or long k into n-k swaptions */
    LGMCalSet* CSPtr)       /* Reference instruments */
{
    LGMErr  error;
    LGM_TS* tsPtr;

    /* initialize */
    error        = NULL;
    *LGMtsPtrPtr = NULL;

    /* Create term structure, fill in  zeta dates, calculate and fill in zeta values */
    tsPtr = CreateTSwithSigma(CSPtr);
    if (tsPtr == NULL)
        return ("error in CreateTSwithSigma");

    /* Pick the G(t) values in the term structure *tsPtr by calibrating to the caplets in *CSPtr */
    if (usecaps == 1)
        error = FitAllGCaps(tsPtr, CSPtr);
    else
        error = FitAllGBkwdSwap(tsPtr, CSPtr);

    if (error != NULL)
    {
        LGMFreeLGM_TS(&tsPtr);
        return (error);
    }

    /* Normalize the calibrated term structure */
    error = LGMVerifyLGM_TS_For_DiagSwpts(CSPtr, tsPtr);
    if (error != NULL)
        LGMFreeLGM_TS(&tsPtr);

    *LGMtsPtrPtr = tsPtr;
    return (error);
}

/*************************************************************************/
/* This routine is the same as above, except that it uses the input arrays
        zeta(t) = Zval[i] at t=Zdate[i], i=1,...,numZ
to define zeta(t) in the term structure */
LGMErr LGMCalGivenZeta(
    LGM_TS**   LGMtsPtrPtr, /* Return: Calibrated term structure */
    long       numZ,        /* number of zeta(t) values */
    Date*      Zdate,       /* dates for zeta(t) */
    double*    Zval,        /* zeta(t) values */
    int        usecaps,     /* calibrate on caplets or long k into n-k swaptions */
    LGMCalSet* CSPtr)       /* Reference caplets */
{
    LGMErr  error;
    LGM_TS* tsPtr;

    /* initialize */
    error        = NULL;
    *LGMtsPtrPtr = NULL;

    /* Create term structure, fill in  numZ, zeta dates, and zeta values */
    tsPtr = CreateTSfromZetas(numZ, Zdate, Zval, CSPtr);
    if (tsPtr == NULL)
        return ("error in CreateTSfromZetas");

    /* Pick the zeta's in the term structure *tsPtr by calibrating to the caplets in CSPtr */
    if (usecaps == 1)
        error = FitAllGCaps(tsPtr, CSPtr);
    else
        error = FitAllGBkwdSwap(tsPtr, CSPtr);

    if (error != NULL)
    {
        LGMFreeLGM_TS(&tsPtr);
        return (error);
    }

    /* Normalize the calibrated term structure */
    error = LGMVerifyLGM_TS_For_DiagSwpts(CSPtr, tsPtr);
    if (error != NULL)
        LGMFreeLGM_TS(&tsPtr);

    *LGMtsPtrPtr = tsPtr;
    return (error);
}

/*****************************************************************************/
/* After the term structure has been created, and the zeta values filled in, this
routine finds the G values by calibrating to the caplets, and then
fills them into the term structure */
static LGMErr FitAllGCaps(
    LGM_TS*    tsPtr, /* Return: calibrated LGM term structure ts */
    LGMCalSet* CSPtr) /* Reference caplets */
{
    LGMErr          error;
    SrtReceiverType recpay;
    long            j, k, nEx;
    double          deltaG, Gend, stddev, zetaex;
    double          dt1, dt2, slope;

    /* initialize */
    error = NULL;
    nEx   = CSPtr->nEx;

    if (tsPtr->numG < nEx + 1)
    {
        srt_free(tsPtr->G);
        srt_free(tsPtr->Gdate);
        tsPtr->G     = (double*)srt_calloc(nEx + 1, sizeof(double));
        tsPtr->Gdate = (Date*)srt_calloc(nEx + 1, sizeof(Date));
        if (tsPtr->G == NULL || tsPtr->Gdate == NULL)
            return ("allocation failed in FitAllGCaps");
    }
    tsPtr->numG = nEx + 1;

    /* fill in Gdates */
    for (j = 0; j < nEx; j++)
        tsPtr->Gdate[j] = CSPtr->tStart[j + 1];
    tsPtr->Gdate[nEx] = CSPtr->tEnd[nEx];

    recpay        = SRT_RECEIVER;
    tsPtr->G[nEx] = 0;

    for (j = nEx; j > 0; j--)
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

        zetaex = LGMZetaFromTS(CSPtr->tEx[j], CSPtr->tNow, tsPtr);
        deltaG = stddev / sqrt(zetaex); /* drop in G from tStart[j] to tEnd[j] */

        if (CSPtr->tEnd[j] <= tsPtr->Gdate[j]) /* tEnd[j] betweeen Gdate[j-1] and Gdate[j] */
        {
            dt1             = (double)(tsPtr->Gdate[j] - tsPtr->Gdate[j - 1]);
            dt2             = (double)(CSPtr->tEnd[j] - CSPtr->tStart[j]);
            tsPtr->G[j - 1] = tsPtr->G[j] + deltaG * dt1 / dt2;
        }
        else
        {
            for (k = j; tsPtr->Gdate[k] < CSPtr->tEnd[j] && k < tsPtr->numG; k++)
                ;
            dt1 = (double)(tsPtr->Gdate[k] - tsPtr->Gdate[k - 1]); /* tEnd[j] between Gdate[k-1] and
                                                                      Gdate[k] */
            dt2             = (double)(CSPtr->tEnd[j] - tsPtr->Gdate[k]);
            slope           = (tsPtr->G[k] - tsPtr->G[k - 1]) / dt1;
            Gend            = tsPtr->G[k] + slope * dt2;
            tsPtr->G[j - 1] = deltaG + Gend;
            if (tsPtr->G[j - 1] < tsPtr->G[j]) /* change to feasible set if needed */
                tsPtr->G[j - 1] = tsPtr->G[j];
        }
    }

    return (error);
}

/******************************************************/
/* Fits the LGM price of a caplet to the market price */
/* price is market price, cvg = coverage of period, Rfix = fixed rate,
   DStart, Dend are dis factors to start and end dates
The output is stddev = sqrt(zeta(tEx))*|G(tStart) - G(tEnd)| */
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

/*****************************************************************************/
/* After the term structure has been created, and the G values filled in, this
routine finds the zeta values by calibrating to the caplets, and then
fills them into the term structure */
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

/*************************************************************************************/
/* After the term structure has been created, and the G values filled in, this
routine finds the zeta values by calibrating to the long swaptions, and then
fills them into the term structure */
static LGMErr FitAllZetaSwaps(LGM_TS* tsPtr, LGMCalSet* CSPtr)
{
    LGMErr          error;
    SrtReceiverType recpay;
    long            i, ifirst, ind, j, n, nlong, nEx, k;
    double *        dispayment, *Gpay, Gst, zeta;

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
            return ("allocation failed in FitAllZetaSwaps");
    }
    tsPtr->numZ = nEx + 1;

    for (k = tsPtr->numG - 1; k > 0; k--)
    {
        if (tsPtr->G[k - 1] < tsPtr->G[k])
            tsPtr->G[k - 1] = tsPtr->G[k];
    }

    dispayment = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));
    Gpay       = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));

    if (dispayment == NULL || Gpay == NULL)
    {
        srt_free(dispayment);
        srt_free(Gpay);
        return ("allocation failed in FitAllZetaSwaps");
    }

    for (i = 1; i <= CSPtr->nPay; i++)
        Gpay[i] = LGMGFromTS(CSPtr->tPay[i], tsPtr); /* get G(t) at pay dates */

    recpay          = SRT_RECEIVER;
    tsPtr->zdate[0] = CSPtr->tNow;
    tsPtr->zeta[0]  = 0.0;

    /* For each j, find the zeta(tEx[j]) that fits the LGM value of
    long swaption j to the market price Vlong[j] */
    for (j = 1; j <= nEx; j++)
    {
        ifirst        = CSPtr->ifirst[j]; /* discounted payments of fixed leg of swaption j */
        nlong         = CSPtr->nlong[j];
        dispayment[1] = CSPtr->Rflong[j] * CSPtr->cvgfirst[j] * CSPtr->Dpay[ifirst];

        for (i = ifirst + 1; i <= nlong; i++)
            dispayment[i + 1 - ifirst] = CSPtr->Rflong[j] * CSPtr->cvgpay[i] * CSPtr->Dpay[i];

        n             = nlong + 1 - ifirst;
        dispayment[n] = dispayment[n] + CSPtr->Dpay[nlong]; /* add notional to final payment */

        Gst = LGMGFromTS(CSPtr->tStart[j], tsPtr);

        ind  = ifirst - 1;
        zeta = FitOneZetaSwap(
            n, dispayment, CSPtr->DStart[j], recpay, &(Gpay[ind]), Gst, CSPtr->Vlong[j]);
        if (zeta < 0)
        {
            srt_free(Gpay);
            srt_free(dispayment);
            return ("failed to fit zeta for swaption");
        }
        tsPtr->zeta[j]  = zeta; /* update term structure */
        tsPtr->zdate[j] = CSPtr->tEx[j];
    }
    srt_free(dispayment);
    srt_free(Gpay);
    return (error);
}

/********************************************************************************/
/* Fits value of zeta at exercise date for a swaption */
/* The fixed leg has the disc payments a[i] = payment[i]*df(tPay[i]), i=1,...,n.
The strike is worth DStart = df(tStart), and  Gpay[i] = G(t) at tPay[i], i=1,...,n
and Gst = G(tStart) */
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

/******************************************************************************/
/* After the term structure has been created and the zeta values filled in, this
routine finds the G(t) values by calibrating to the long swaptions and
fills them into the term structure */
static LGMErr FitAllGBkwdSwap(
    LGM_TS*    tsPtr, /* Return: zeta's in LGM term structure tsPtr */
    LGMCalSet* CSPtr) /* Reference instruments */
{
    LGMErr          error;
    SrtReceiverType recpay;
    long            i, ifirst, previfirst, nlong;
    long            ind, j, n, nEx;
    double          zeta, Gst, dGst;
    double *        dispayment, *Gpay, *dGpay;
    double          dt1, dt2;

    /* initialize */
    error = NULL;
    nEx   = CSPtr->nEx; /* number of swaptions to fit */
    n     = CSPtr->n;   /* we need to fit term structure only for t<tPay[n] */

    if (tsPtr->numG < nEx + 1)
    {
        srt_free(tsPtr->G);
        srt_free(tsPtr->Gdate);
        tsPtr->G     = (double*)srt_calloc(nEx + 1, sizeof(double));
        tsPtr->Gdate = (Date*)srt_calloc(nEx + 1, sizeof(Date));
        if (tsPtr->G == NULL || tsPtr->Gdate == NULL)
            return ("allocation failed in FitAllGBckwdSwap");
    }
    tsPtr->numG = nEx + 1;

    dispayment = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));
    Gpay       = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));
    dGpay      = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));

    if (dispayment == NULL || Gpay == NULL || dGpay == NULL)
    {
        LGMFreeLGM_TS(&tsPtr);
        srt_free(dispayment);
        srt_free(Gpay);
        srt_free(dGpay);
        return ("allocation failed in FitAllGBkwd");
    }

    for (j = 0; j < nEx; j++)
        tsPtr->Gdate[j] = CSPtr->tStart[j + 1];
    tsPtr->Gdate[nEx] = CSPtr->tPay[n];
    tsPtr->G[nEx]     = 0.0;

    recpay = SRT_RECEIVER;

    for (i = 0; i <= CSPtr->nPay; i++)
    {
        Gpay[i]  = 0.;
        dGpay[i] = 0.;
    }

    /* For each j, find G(tStart[j]) so that the LGM value of long swaption j
    agrees with the market price Vlong[j] */
    previfirst = CSPtr->nPay + 1;
    for (j = nEx; j > 0; j--)
    {
        ifirst = CSPtr->ifirst[j]; /* compute discounted payments of fixed leg of swaption j */
        nlong  = CSPtr->nlong[j];
        dispayment[1] = CSPtr->Rflong[j] * CSPtr->cvgfirst[j] * CSPtr->Dpay[ifirst];

        for (i = ifirst + 1; i <= nlong; i++)
            dispayment[i + 1 - ifirst] = CSPtr->Rflong[j] * CSPtr->cvgpay[i] * CSPtr->Dpay[i];

        n             = nlong + 1 - ifirst;
        dispayment[n] = dispayment[n] + CSPtr->Dpay[nlong]; /* add notional to final payment */

        /* update up Gpay and dGpay arrays */
        Gst  = tsPtr->G[j];
        dGst = 1.0;
        for (i = ifirst; i <= nlong; i++)
            dGpay[i] = 0.0;

        for (i = ifirst; i < previfirst; i++)
        {
            Gpay[i]  = tsPtr->G[j];
            dt1      = (double)(tsPtr->Gdate[j] - CSPtr->tPay[i]);
            dt2      = (double)(tsPtr->Gdate[j] - tsPtr->Gdate[j - 1]);
            dGpay[i] = dGst * dt1 / dt2;
        }
        zeta = LGMZetaFromTS(CSPtr->tEx[j], CSPtr->tNow, tsPtr);
        ind  = ifirst - 1;
        zeta = FitOneSwap(
            n,
            dispayment,
            CSPtr->DStart[j],
            recpay, /* the swaption */
            zeta,   /* zeta(tEx) */
            &Gst,
            &(Gpay[ind]), /* the term structure */
            dGst,
            &(dGpay[ind]),    /* direction G */
            CSPtr->Vlong[j]); /* market price */

        if (zeta < 0.0)
        {
            LGMFreeLGM_TS(&tsPtr);
            srt_free(Gpay);
            srt_free(dGpay);
            srt_free(dispayment);
            return ("failed to fit G for swaption");
        }

        if (Gst > tsPtr->G[j]) /* replace by feasible value if needed */
            tsPtr->G[j - 1] = Gst;
        else
            tsPtr->G[j - 1] = tsPtr->G[j];

        for (i = ifirst; i < previfirst; i++)
        {
            dt1     = (double)(tsPtr->Gdate[j] - CSPtr->tPay[i]);
            dt2     = (double)(tsPtr->Gdate[j] - tsPtr->Gdate[j - 1]);
            Gpay[i] = tsPtr->G[j] + (tsPtr->G[j - 1] - tsPtr->G[j]) * dt1 / dt2;
        }
        previfirst = CSPtr->ifirst[j];
    }
    srt_free(dispayment);
    srt_free(Gpay);
    srt_free(dGpay);
    return (error);
}

/******************************************************************/
/* Fixed expiry date, calibrate on 1 into k swaptions and caplets */
/******************************************************************/
/* This routine creates a term structure and calibrates it to the set of fixed expiry
(1 into k) swaptions and the set of caplets in the calibration set CSPtr.
It first constructs G(t) from by calibrating the 1 into k swaptions to their
market value. It then finds zeta(t) at tNow, and at each exercise date tEx by
calibrating the caplets to their market values */
/* The following information must be set in CSPtr:
        tNow						evaluation date
Definition of the 1 into k swaptions:
        ftEx							exercise date for all 1 into k
swaptions ftStart							start date for all 1 into k
swaptions
        DfixStart							dis factor from tNow to
common start date fstt tPay[i] for i=ffirst,...,nlast	common set of pay dates for 1 into k
swaptions
        ffirst							pay dates of swaption j are tPay[i]
for nfirst							i = ffirst, ..., j nlast
for j=nfirst, ..., nlast cvgfix							coverage from start
date to first pay date cvgpay[i], i=ffirst+1,...,nlast	coverage from tPay[i-1] to tPay[i] Dpay[i],
i=ffirst, ..., nlast	dis factor from tNow to tPay[i] Rfix[j] for j=nfirst,...,nlast	fixed rate
for caplet j Vfix[j]	for j=nfirst,...,nlast	market value at tNow for caplet j
*/
/* If usecaps == 1, the caplets must also be defined:
        nEx>1						number of reference caplets
        tEx[j] for j=1,...,nEx		exercise date for caplet j
        tStart[j] for j=1,...,nEx		start date for caplet j
        DStart[j] for j=1,...,nEx		dis factor from tNow to the start date
        tEnd[j] for j=1,...,nEx		end date for caplet j
        Dcap[j] for j=1,...,nEx		dis factor to the end date of caplet j
        cvgcap[j] for j=1,...,nEx	coverage from start date to end date for caplet j
        Rfcap[j] for j=1,...,nEx	fixed rate for caplet j
        Vcap[j]	for j=1,...,nEx		market value at tNow for caplet j
Note the indices i=0 and j=0 is ignored in the above arrays
*/
/* If usecaps != 1, the long "k into n-k" swaptions must also be defined:
        nEx>1						number of reference swaptions
        tEx[j] for j=1,...,nEx		exercise date for swaption j
        tStart[j] for j=1,...,nEx		start date for swaption j
        DStart[j] for j=1,...,nEx		dis factor from tNow to the start date
        nPay						total number of fixed leg pay dates
        tPay[i] for i=1,...,nPay	fixed leg pay dates
        cvgpay[i] for i=1,...,nPay	coverage from tPay[i-1] to tPay[i]
        Dpay[i] for i=1,...,nPay	dis factor from tNow to tPay[i]

        ifirst [j] for j=1,...,nEx	fixed leg pay dates are i = ifirst[j] to nlong[j]
        nlong [j] for j=1,...,nEx
        cvgfirst t[j] for j=1,...,nEx  first period coverage (from tStart[j] to tPay[ifirst[j]])
        Rlong [j] for j=1,...,nEx	fixed rate for swaption j
        Vlong [j] for j=1,...,nEx	market price of swaption j
        Note the indices j=0 and i=0 are ignored in the above arrays
*/
LGMErr LGMCalFixExp(
    LGM_TS**   LGMtsPtrPtr, /* Return: Calibrated term structure */
    int        usecaps,     /* calibrate on caplets or long "k into n-k" swaptions */
    int        keep,        /* keep 1intok zeta */
    LGMCalSet* CSPtr)       /* Reference caplets */
{
    LGMErr  error;
    LGM_TS* tsPtr = NULL;
    long    numG, numZ, j0, j, nEx;
    Date    zdatefix; /* tEx for "1 into k" */
    double  zetafix;  /* zeta(tEx) for "1 into k" */
    double  dint, dt, slope, zmax, zmin;

    /* initialize */
    error = NULL;
    nEx   = CSPtr->nEx;

    /* Create term structure */
    numG  = CSPtr->nlast - CSPtr->nfirst + 2;
    numZ  = nEx + 2;
    tsPtr = LGMCreateLGM_TS(numZ, numG);
    if (tsPtr == NULL)
    {
        LGMFreeLGM_TS(&tsPtr);
        return ("Memory Allocation Error");
    }

    /* Find the G(t) values by calibrating on the 1 into k swaptions
            and fill numG, the G dates, and the G values */
    error = FitAllGFwdSwap(tsPtr, CSPtr, &zdatefix, &zetafix);
    if (error != NULL)
    {
        LGMFreeLGM_TS(&tsPtr);
        return (error);
    }

    /* Pick the zeta's in the term structure *tsPtr by calibrating to the caplets */
    /* or to the long swaptions in CSPtr */
    if (usecaps == 1)
        error = FitAllZetaCaps(tsPtr, CSPtr);
    else
        error = FitAllZetaSwaps(tsPtr, CSPtr);

    if (error != NULL)
    {
        LGMFreeLGM_TS(&tsPtr);
        return (error);
    }

    /* Slip the zeta at the exercise date of the "1 into k" into the
    term structure */
    if (keep < 0 || zdatefix <= tsPtr->zdate[0]) /* do nothing */
        tsPtr->numZ = nEx + 1;

    else if (zdatefix >= tsPtr->zdate[nEx] + 365) /* append to end */
    {
        tsPtr->zdate[nEx + 1] = zdatefix;
        tsPtr->zeta[nEx + 1]  = zetafix;
        tsPtr->numZ           = nEx + 2;
    }

    else if (zdatefix >= tsPtr->zdate[nEx])
    {
        if (keep > 0) /* replace last */
        {
            tsPtr->zdate[nEx] = zdatefix;
            tsPtr->zeta[nEx]  = zetafix;
            tsPtr->numZ       = nEx + 1;
        }
        else
            tsPtr->numZ = nEx + 1; /* do nothing */
    }

    else
    {
        for (j0 = 0; tsPtr->zdate[j0] <= zdatefix; j0++)
            ;
        if (keep > 0)
        {
            if (j0 > 1 && (abs(tsPtr->zdate[j0 - 1] - zdatefix) < abs(tsPtr->zdate[j0] - zdatefix)))
                j0 = j0 - 1;
            tsPtr->zdate[j0] = zdatefix; /* replace j0 */
            tsPtr->zeta[j0]  = zetafix;
            tsPtr->numZ      = nEx + 1;
        }

        else if (tsPtr->zdate[j0 - 1] == zdatefix)
            tsPtr->numZ = nEx + 1; /* do nothing */

        else /* insert at j0 */
        {
            dint  = ((double)(tsPtr->zdate[j0] - tsPtr->zdate[j0 - 1])) / 365.0;
            dt    = ((double)(zdatefix - tsPtr->zdate[j0 - 1])) / 365.0;
            slope = (tsPtr->zeta[j0] - tsPtr->zeta[j0 - 1]) / dint;
            zmax =
                min(tsPtr->zeta[j0 - 1] + 2.5 * slope * dt,
                    tsPtr->zeta[j0] - 0.4 * slope * (dint - dt));
            zmin =
                max(tsPtr->zeta[j0 - 1] + 0.4 * slope * dt,
                    tsPtr->zeta[j0] - 2.5 * slope * (dint - dt));
            if (zetafix > zmax)
                zetafix = zmax;
            if (zetafix < zmin)
                zetafix = zmin;

            for (j = nEx; j >= j0; j--)
            {
                tsPtr->zdate[j + 1] = tsPtr->zdate[j];
                tsPtr->zeta[j + 1]  = tsPtr->zeta[j];
            }
            tsPtr->zdate[j0] = zdatefix;
            tsPtr->zeta[j0]  = zetafix;
            tsPtr->numZ      = nEx + 2;
        }
    }

    /* Normalize the calibrated term structure */
    error = LGMVerifyLGM_TS_For_DiagSwpts(CSPtr, tsPtr);
    if (error != NULL)
        LGMFreeLGM_TS(&tsPtr);

    *LGMtsPtrPtr = tsPtr;
    return (error);
}

/**********************************************************************************/
/* After the term structure has been created this routine finds the G(t) values
by calibrating to the "1 into k" swaptions and fills them into the term structure */
LGMErr FitAllGFwdSwap(LGM_TS* tsPtr, LGMCalSet* CSPtr, Date* zdatePtr, double* zetaPtr)
{
    LGMErr          error;
    SrtReceiverType recpay;
    long            i, j, ifirst, ind, nPay;
    long            nfirst, nlast, n, numG;
    double          zeta, Gst, dGst;
    double *        dispayment, *Gpay, *dGpay;
    double          slope, dt, dummy;
    Date            zdate;

    error = NULL;

    /* allocate space */
    nPay       = CSPtr->nPay;
    dispayment = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));
    Gpay       = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));
    dGpay      = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));

    if (dispayment == NULL || Gpay == NULL || dGpay == NULL)
    {
        srt_free(dispayment);
        srt_free(Gpay);
        srt_free(dGpay);
        return ("allocation failed in FitAllGFwd");
    }

    /* initialize */
    ifirst = CSPtr->ffirst;
    nfirst = CSPtr->nfirst;
    nlast  = CSPtr->nlast;

    numG = nlast - nfirst + 2;
    if (tsPtr->numG < numG)
    {
        srt_free(tsPtr->G);
        srt_free(tsPtr->Gdate);
        tsPtr->G     = (double*)srt_calloc(numG, sizeof(double));
        tsPtr->Gdate = (Date*)srt_calloc(numG, sizeof(Date));
        if (tsPtr->G == NULL || tsPtr->Gdate == NULL)
            return ("allocation failed in FitAllGFwsSwaps");
    }
    tsPtr->numG = numG;

    ind   = ifirst - 1;
    zdate = CSPtr->ftEx;
    zeta  = 0.0001 * ((double)(zdate - CSPtr->tNow)) / 365.0;

    recpay = SRT_RECEIVER;

    tsPtr->Gdate[0] = CSPtr->ftStart;
    tsPtr->G[0]     = 0.0;
    for (j = CSPtr->nfirst; j <= CSPtr->nlast; j++)
    {
        i               = j + 1 - CSPtr->nfirst;
        tsPtr->Gdate[i] = CSPtr->tPay[j];
    }

    for (i = 0; i <= nPay; i++)
    {
        Gpay[i]  = 0.;
        dGpay[i] = 0.;
    }

    /* For each j, find G(tStart[j]) so that the LGM value of long swaption j
    agrees with the market price Vlong[j] */
    for (j = nfirst; j <= nlast; j++)
    {
        dispayment[1] =
            CSPtr->Rfix[j] * CSPtr->cvgfix * CSPtr->Dpay[ifirst]; /* compute discounted payments */
        for (i = ifirst + 1; i <= j; i++)                         /* of fixed leg of swaption j */
            dispayment[i + 1 - ifirst] = CSPtr->Rfix[j] * CSPtr->cvgpay[i] * CSPtr->Dpay[i];
        n             = j + 1 - ifirst;
        dispayment[n] = dispayment[n] + CSPtr->Dpay[j]; /* add notional to final payment */

        Gst  = 0;
        dGst = 0.;

        if (j == nfirst)
        {
            for (i = ifirst; i <= CSPtr->nfirst; i++)
                dGpay[i] = ((double)(CSPtr->ftStart - CSPtr->tPay[i])) / 365.0;
        }
        else
        {
            for (i = ifirst; i < j; i++)
                dGpay[i] = 0.0;
            dGpay[j] = -1.0;
            Gpay[j]  = Gpay[j - 1];
        }
        dummy = FitOneSwap(
            n,
            dispayment,
            CSPtr->DfixStart,
            recpay, /* the swaption */
            zeta,   /* zeta(tEx) */
            &Gst,
            &(Gpay[ind]), /* the term structure */
            dGst,
            &(dGpay[ind]),   /* direction G */
            CSPtr->Vfix[j]); /* market price */

        if (dummy <= -1.0)
        {
            srt_free(Gpay);
            srt_free(dGpay);
            srt_free(dispayment);
            return ("failed to fit G for swaption");
        }
        tsPtr->G[j + 1 - nfirst] = Gpay[j];
    }

    /* Normalize G(t) in term structure */
    tsPtr->numG = nlast + 2 - nfirst;
    dt          = ((double)(tsPtr->Gdate[tsPtr->numG - 1] - tsPtr->Gdate[0])) / 365.0;
    slope       = (tsPtr->G[0] - tsPtr->G[tsPtr->numG - 1]) / dt;
    for (j = 0; j < tsPtr->numG; j++)
        tsPtr->G[j] = (tsPtr->G[j] - tsPtr->G[tsPtr->numG - 1]) / slope;

    zeta      = zeta * slope * slope;
    *zdatePtr = zdate;
    *zetaPtr  = zeta;

    /* free and leave */
    srt_free(dispayment);
    srt_free(Gpay);
    srt_free(dGpay);
    return (error);
}

/************************************************************************/
/* Fixed tenor and diag, calibrate on k into n-k and caplets			*/
/************************************************************************/
/* This routine creates a term structure and calibrates it to the set of long
(k into n-k) swaptions and the set of caplets in the calibration set CSPtr to
their market values.It constructs both G(t) and zeta(t) simultaneously by
a backwards bootstrap calibration method. */
/* The following information must be set in CSPtr:
        tNow						evaluation date
Definition of the long "k into n-k" swaptions:
        nEx>1						number of reference swaptions
        tEx[j] for j=1,...,nEx		exercise date for swaption j
        tStart[j] for j=1,...,nEx		start date for swaption j
        DStart[j] for j=1,...,nEx		dis factor from tNow to the start date
        nPay						total number of fixed leg pay dates
        tPay[i] for i=1,...,nPay	fixed leg pay dates
        cvgpay[i] for i=1,...,nPay	coverage from tPay[i-1] to tPay[i]
        Dpay[i] for i=1,...,nPay	dis factor from tNow to tPay[i]

        ifirst [j] for j=1,...,nEx	fixed leg pay dates are i = ifirst[j] to nlong[j]
        nlong [j] for j=1,...,nEx
        cvgfirst t[j] for j=1,...,nEx  first period coverage (from tStart[j] to tPay[ifirst[j]])
        Rlong [j] for j=1,...,nEx	fixed rate for swaption j
        Vlong [j] for j=1,...,nEx	market price of swaption j
Note the indices j=0 and i=0 are ignored in the above arrays
Definition of the caplets also requires:
        tEnd[j] for j=1,...,nEx		end date for caplet j
        Dcap[j] for j=1,...,nEx		dis factor to the end date of caplet j
        cvgcap[j] for j=1,...,nEx	coverage from start date to end date for caplet j
        Rfcap[j] for j=1,...,nEx	fixed rate for caplet j
        Vcap[j]	for j=1,...,nEx		market value at tNow for caplet j */
/* Note the indices i=0 and j=0 is ignored in the above arrays */
LGMErr LGMCalTenorDiag(
    LGM_TS**   LGMtsPtrPtr, /* Return: Calibrated term structure */
    int        usecaps,     /* calibrate on caplets or short "k into 1" swaptions */
    LGMCalSet* CSPtr,       /* Reference instruments */
    long*      noPairsFlag)      /* signals whether there are any pairs to work with */
{
    LGMErr  error;
    LGM_TS* tsPtr;
    long    npair, *jpair;

    /* initialize */
    error = NULL;

    /* Create term structure */
    /* This is larger than the eventual term structure */
    tsPtr = LGMCreateLGM_TS(CSPtr->nEx + 1, CSPtr->nEx + 1);

    jpair = (long*)srt_calloc(CSPtr->nEx + 1, sizeof(long));

    if (tsPtr == NULL || jpair == NULL)
    {
        LGMFreeLGM_TS(&tsPtr);
        srt_free(jpair);
        return ("allocation error");
    }

    /* Determine which pairs of swaptions to use */
    if (usecaps == 1)
        npair = findswapcappairs(jpair, CSPtr);
    else
        npair = findswappairs(jpair, CSPtr);
    if (npair <= 0)
    {
        LGMFreeLGM_TS(&tsPtr);
        srt_free(jpair);
        *noPairsFlag = 1;
        return (NULL);
    }
    tsPtr->numG = npair + 1;

    /* Find the G(t) values by calibrating on the long/short
            swaption pairs or the long/caplet swaption pairs,
            and fill numG, the G dates, and the G values */
    error = FitAllGPairs(npair, jpair, usecaps, tsPtr, CSPtr);
    srt_free(jpair);

    if (error != NULL)
    {
        LGMFreeLGM_TS(&tsPtr);
        return (error);
    }

    /* Pick the zeta's in term structure *tsPtr by calibrating
            to long swaptions in CSPtr */
    error = FitAllZetaSwaps(tsPtr, CSPtr);

    if (error != NULL)
    {
        LGMFreeLGM_TS(&tsPtr);
        return (error);
    }

    /* Normalize the calibrated term structure */
    /*error = LGMVerifyLGM_TS(tsPtr);*/
    error = LGMVerifyLGM_TS_For_DiagSwpts(CSPtr, tsPtr);
    if (error != NULL)
        LGMFreeLGM_TS(&tsPtr);

    *LGMtsPtrPtr = tsPtr;
    return (error);
}

/**************************************/
/* Determine swaption/caplet pairs to use to find G(t) */
/* Use only caplets whose start date is on or after the previously used
 caplets end date.
Note that the code doesn't use the last pair found */
static long findswapcappairs(long* jpair, LGMCalSet* CSPtr)
{
    long j, jmax, k;
    Date tprev;

    /* use only pairs whose long swaptions end date is tPay[n] or earlier */
    for (jmax = 1; jmax <= CSPtr->nEx && CSPtr->nlong[jmax] <= CSPtr->n; jmax++)
        ;
    jmax--;
    /* use only caplets whose start date is on or after the previously used
     caplet's end date */
    tprev = CSPtr->tStart[1];
    k     = 0;
    for (j = 1; j <= jmax; j++)
    {
        if (CSPtr->tStart[j] >= tprev)
        {
            k++;
            jpair[k] = j;
            tprev    = CSPtr->tEnd[j];
        }
    }
    return (k); /* number of pairs to use for backward bootstrap */
}

/**********************************************************/
/* Determine long/short swaption pairs to use to find G(t) */
/* Use only pairs whose last paydate for short swaption is
greater than previous pair's, and whose short swaption is shorter
than the long swaption.
Note: the calibration code does not use the last pair in jpair[k] */
static long findswappairs(long* jpair, LGMCalSet* CSPtr)
{
    long j, k;

    k = 0;
    for (j = 1; j <= CSPtr->nEx && CSPtr->nshort[j] <= CSPtr->n; j++)
    {
        if (j == CSPtr->nEx || CSPtr->nshort[j] < CSPtr->nshort[j + 1])
        {
            k++; /* new last short paydate. Use this pair */
            jpair[k] = j;
        }
    }

    return (k); /* number of pairs to use for backward bootstrap */
}

/***************************************************************/
static LGMErr FitAllGPairs(long npair, long* jpair, long usecaps, LGM_TS* tsPtr, LGMCalSet* CSPtr)
{
    LGMErr          error;
    SrtReceiverType recpay;
    long            i, j, k, ifirst, previfirst, jprev, nlong;
    long            ind, n, nsw, nsw2, nshort;
    long            numG;
    double          Gst, dGst;
    double *        dispaylong, *dispayshort, *Gpay, *dGpay;
    double          dt1, dt2, stddev, product, zeta;
    double          scale;

    /* initialize */
    error = NULL;
    n     = CSPtr->n; /* we need to fit term structure only for t<tPay[n] */

    dispaylong  = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));
    dispayshort = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));
    Gpay        = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));
    dGpay       = (double*)srt_calloc(CSPtr->nPay + 1, sizeof(double));

    if (dispaylong == NULL || dispayshort == NULL || Gpay == NULL || dGpay == NULL)
    {
        srt_free(dispaylong);
        srt_free(dispayshort);
        srt_free(Gpay);
        srt_free(dGpay);
        return ("allocation failed in FitAllGPairs");
    }

    numG = npair + 1;
    if (tsPtr->numG < numG)
    {
        srt_free(tsPtr->G);
        srt_free(tsPtr->Gdate);
        tsPtr->G     = (double*)srt_calloc(numG, sizeof(double));
        tsPtr->Gdate = (Date*)srt_calloc(numG, sizeof(Date));
        if (tsPtr->G == NULL || tsPtr->Gdate == NULL)
            return ("allocation failed in FitAllGFwsSwaps");
    }
    tsPtr->numG = numG;

    for (k = 1; k <= npair; k++)
    {
        j                   = jpair[k];
        tsPtr->Gdate[k - 1] = CSPtr->tStart[j];
    }
    tsPtr->Gdate[npair] = CSPtr->tPay[n];
    tsPtr->G[npair]     = 0.0;
    tsPtr->G[npair - 1] = ((double)(tsPtr->Gdate[npair] - tsPtr->Gdate[npair - 1])) / 365.0;

    recpay = SRT_RECEIVER;

    for (i = 0; i <= CSPtr->nPay; i++)
    {
        Gpay[i]  = 0.;
        dGpay[i] = 0.;
    }

    /* Fill in G at last start date and later */
    j      = jpair[npair];
    ifirst = CSPtr->ifirst[j];
    for (i = ifirst; i <= CSPtr->nPay; i++)
        Gpay[i] = ((double)(CSPtr->tPay[n] - CSPtr->tPay[i])) / 365.0;

    /* For k = npair-1, npair-2, ..., 1, for each j=jpair[k] find G(tStart[j])
    so that the LGM value of long swaption j and either the LGM value of the caplet
    or short swaption agrees with their market prices */

    previfirst = CSPtr->ifirst[j];
    for (k = npair - 1; k >= 1; k--)
    {
        j     = jpair[k];
        jprev = jpair[k + 1];

        /* compute discounted payments of fixed leg of long swaption j */
        ifirst        = CSPtr->ifirst[j];
        nlong         = CSPtr->nlong[j];
        dispaylong[1] = CSPtr->Rflong[j] * CSPtr->cvgfirst[j] * CSPtr->Dpay[ifirst];

        for (i = ifirst + 1; i <= nlong; i++)
            dispaylong[i + 1 - ifirst] = CSPtr->Rflong[j] * CSPtr->cvgpay[i] * CSPtr->Dpay[i];
        nsw             = nlong + 1 - ifirst;
        dispaylong[nsw] = dispaylong[nsw] + CSPtr->Dpay[nlong];

        /* update up Gpay and dGpay arrays */
        Gst  = tsPtr->G[k];
        dGst = 1.0;
        for (i = ifirst; i <= nlong; i++)
            dGpay[i] = 0.0;

        for (i = ifirst; i < previfirst; i++)
        {
            Gpay[i]  = tsPtr->G[k];
            dt1      = (double)(tsPtr->Gdate[k] - CSPtr->tPay[i]);
            dt2      = (double)(tsPtr->Gdate[k] - tsPtr->Gdate[k - 1]);
            dGpay[i] = dGst * dt1 / dt2;
        }

        /* if usecaps==1, fit caplet j */
        if (usecaps == 1)
        {
            stddev = FitOneCap(
                CSPtr->Vcap[j],
                CSPtr->cvgcap[j],
                CSPtr->Rfcap[j],
                CSPtr->DStart[j],
                CSPtr->Dcap[j],
                recpay);
            dt1     = (double)(CSPtr->tEnd[j] - CSPtr->tStart[j]);
            dt2     = (double)(CSPtr->tStart[jprev] - CSPtr->tStart[j]);
            product = stddev * dt2 / dt1;
            ind     = ifirst - 1;
            zeta    = FitOneSwapCap(
                nsw,
                dispaylong,
                CSPtr->DStart[j],
                recpay, /* the swaption */
                zeta,
                product, /* unused, product of sqzeta*dGst */
                &Gst,
                &(Gpay[ind]), /* the existing term structure */
                dGst,
                &(dGpay[ind]),    /* direction of changes in G */
                CSPtr->Vlong[j]); /* market price */
        }
        /* if usecaps!=1, compute discounted payments of short swaption j */
        else
        {
            nshort         = CSPtr->nshort[j];
            dispayshort[1] = CSPtr->Rfshort[j] * CSPtr->cvgfirst[j] * CSPtr->Dpay[ifirst];
            for (i = ifirst + 1; i <= nshort; i++)
                dispayshort[i + 1 - ifirst] = CSPtr->Rfshort[j] * CSPtr->cvgpay[i] * CSPtr->Dpay[i];
            nsw2              = nshort + 1 - ifirst;
            dispayshort[nsw2] = dispayshort[nsw2] + CSPtr->Dpay[nshort];

            /* fit pair of long & short swaptions*/
            ind  = ifirst - 1;
            zeta = FitOneLongShort(
                nsw,
                dispaylong, /* the swaption 1 */
                nsw2,
                dispayshort, /* the swaption 2 */
                CSPtr->DStart[j],
                recpay, /* the swaptions */
                &Gst,
                &(Gpay[ind]), /* the existing term structure */
                dGst,
                &(dGpay[ind]), /* direction of changes in G */
                CSPtr->Vlong[j],
                CSPtr->Vshort[j]); /* market price */
        }

        if (zeta < 0.)
        {
            srt_free(dispaylong);
            srt_free(dispayshort);
            srt_free(Gpay);
            srt_free(dGpay);
            return ("failed to fit pair in FitAllGPairs");
        }

        tsPtr->G[k - 1] = Gst;
        for (i = ifirst; i < previfirst; i++)
        {
            dt1     = (double)(tsPtr->Gdate[k] - CSPtr->tPay[i]);
            dt2     = (double)(tsPtr->Gdate[k] - tsPtr->Gdate[k - 1]);
            Gpay[i] = tsPtr->G[k] + (tsPtr->G[k - 1] - tsPtr->G[k]) * dt1 / dt2;
        }
        previfirst = CSPtr->ifirst[j];
    }

    /* finished calibrating G[k] ... re-scale G */
    scale = (tsPtr->G[0] - tsPtr->G[numG - 1]) * 365.0 /
            ((double)(tsPtr->Gdate[numG - 1] - tsPtr->Gdate[0]));
    for (j = 0; j < numG; j++)
        tsPtr->G[j] = (tsPtr->G[j] - tsPtr->G[numG - 1]) / scale;

    srt_free(dispaylong);
    srt_free(dispayshort);
    srt_free(Gpay);
    srt_free(dGpay);
    return (error);
}

/******************************************************************/
/* Picks the value of C so that the LGM price of a long swaption
matches the market price when
        Gst = Gst + C*dGst							(G at start date)
        Gpay[i] = Gpay[i] + C*dGpay[i], i=1, ...,n	(G at pay dates)

Zeta is held fixed at its input value.

The routine sets Gpay[i] to its new value Gpay[i] + C*dGpay[i]
sets *Gstptr to Gst's new value, Gst+C*dGst, and returns the
 input zeta(tEx).

On error, returns -1. Note that zeta must be positive */
/* The fixed leg has discounted payments a[i] = payment[i]*df(tPay[i]), i=1,...,n.
The strike is worth DStart = dis factor(tStart) */
double FitOneSwap(
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
    itermax = 256;
    callput = SRT_CALL;
    zeta    = fabs(zeta);
    sqzeta  = sqrt(zeta);
    Gst     = *Gstptr;

    /* check for intrinsic */
    /*	if (fabs (LGMRecVal(&ystar, n, a, DStart, 1.0e-15, Gpay, Gst) - price) < 1.0e-08)
            {
                    return 1.0e-15;
            }*/

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

/******************************************************************/
/* Picks the value of C so that the LGM price of a long swaption
matches the market price when
        Gst = Gst + C*dGst							(G at start date)
        Gpay[i] = Gpay[i] + C*dGpay[i], i=1, ...,n	(G at pay dates)

zeta(tEx) is varied so that the product
        C*dGst*sqrt(zeta)
is held fixed. Provided product has been properly
chosen, this keeps the caplet with fixing date tEx, start date tStart and end date
calibrated to its market value

The routine sets Gpay[i] to its new value Gpay[i] + C*dGpay[i]
sets *Gstptr to Gst's new value, Gst+C*dGst, and returns the
calibrated zeta(tEx).

If both swaption and caplet cannot be fit, it will fit only the swaption,
using the input value of zeta.

On error, returns -1. Note that zeta must be positive */
/* The fixed leg has discounted payments a[i] = payment[i]*df(tPay[i]), i=1,...,n.
The strike is worth DStart = dis factor(tStart) */
static double FitOneSwapCap(
    long            n,
    double*         a,
    double          DStart,
    SrtReceiverType recpay, /* the swaption */
    double          zeta,
    double          product, /* zeta, product */
    double*         Gstptr,
    double*         Gpay, /* the existing term structure */
    double          dGst,
    double*         dGpay, /* direction of change */
    double          price)          /* market price of swaption */
{
    SrtCallPutType callput;
    double         Gst = *Gstptr;
    double         fwd, sens, dsens, favg, stddev, ratio;
    double         value, newvalue, timevalue, der, ystar, sqzeta, coef, dcoef;
    double         theta, arg, z, dz, isign, expfac;
    long           i, iter, itermax, success;

    /* check for intrinsic */
    if (fabs(LGMRecVal(&ystar, n, a, DStart, 1.0e-15, Gpay, Gst) - price) < 1.0e-08)
    {
        return 1.0e-15;
    }

    /* initialize */
    itermax = 25;
    success = 0;
    callput = SRT_CALL;

    if (product <= 0.0)
        return (-1.0);

    for (i = 1; i <= n; i++)
        dGpay[i] = dGpay[i] / dGst;
    product = product / dGst;
    dGst    = 1.0;

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

    if (sens < -1.e-6) /* slopes of G are wrong */
        return (-1.0);

    if (sens < 1.e-04 || dsens < 1.e-05)
        return (zeta); /* can't fit both, best fit is C=0 */

    /* see if the swaption can be fit */
    value = LGMRecVal(&ystar, n, a, DStart, product, dGpay, dGst);
    if (value >= price)
        return (zeta); /* can't fit both, best fit is C=0 */

    /* get initial guess by linear approximation (neglects convexity) */
    favg   = 0.5 * (fwd + DStart);
    stddev = LGMNormImpVol(price / favg, fwd / favg, DStart / favg, 1.0, 1.0, callput);
    stddev = stddev * favg;
    if (stddev < 0.0)
        return (-1.0);

    ratio = product / (stddev * dGst);

    /* global newton to get initial guess */
    z       = 0.;
    dz      = (1.0 - ratio * dsens) / sens;
    success = -1;
    value   = -1.0;
    for (iter = 1; iter <= itermax && success <= 0; iter++)
    {
        newvalue = -1.0; /* calculate func & der at z+dz */
        der      = 0.;
        for (i = 1; i <= n; i++)
        {
            arg    = (z + dz) * (Gst - Gpay[i]) + ratio * (dGst - dGpay[i]);
            expfac = exp(-theta * arg);
            if (fabs(theta) <= 1.e-5 || fabs(arg) <= 1.e-4)
                newvalue =
                    newvalue + a[i] * arg * (1. - 0.5 * arg * theta * (1.0 - arg * theta / 3.0));
            else
                newvalue = newvalue + a[i] * (1.0 - expfac) / theta;
            der = der + a[i] * (Gst - Gpay[i]) * expfac;
        }

        if (fabs(newvalue) < fabs(value)) /* if improvement, accept step */
        {
            z     = z + dz;
            value = newvalue;
            dz    = -newvalue / der; /* new Newton step */
        }
        else
            dz = 0.5 * dz; /* else cut step in half */

        if (fabs(newvalue) <= 1.e-06) /* check convergence */
            success++;                /* go one step past OK convergence */
    }
    z = z + dz;

    /* Use global Newton on coef to match LGM price to market price */
    coef    = 0.;
    dcoef   = ratio / z;
    value   = -price;
    success = -3;
    isign   = 1.0;

    for (iter = 1; iter <= itermax && success <= 0; iter++)
    {
        coef = coef + isign * dcoef; /* update term structure */
        Gst  = Gst + isign * dcoef * dGst;
        for (i = 1; i <= n; i++)
            Gpay[i] = Gpay[i] + isign * dcoef * dGpay[i];
        sqzeta = product / (coef * dGst);

        newvalue = LGMRecVal(&ystar, n, a, DStart, sqzeta, Gpay, Gst) - price;

        if (fabs(newvalue) / (timevalue + .001) < 1.e-5 ||
            fabs(dcoef) < 1.e-6) /* check convergence */
            success++;

        if (fabs(newvalue) < fabs(value))
        {
            value = newvalue; /* step successful */
            isign = 1.0;      /* take new step */
            der   = RecGDer(ystar, n, a, DStart, sqzeta, Gpay, Gst, dGst, dGpay);
            der   = der - (product / (coef * coef * dGst)) *
                            RecSqZetaDer(ystar, n, a, DStart, sqzeta, Gpay, Gst);
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
        sqzeta = product / (coef * dGst);

        *Gstptr = Gst;
        zeta    = sqzeta * sqzeta;
        return (zeta);
    }
    else
        return (-1.0);
}

/***************************************************************************/
/* Picks the value of coef and sqzeta so that the LGM prices of the long and
short swaptions match their market price when
        Gst = Gst + coef*dGst							(G at start date)
        Gpay[i] = Gpay[i] + coef*dGpay[i], i=1, ...,n	(G at pay dates)
        zeta(tEx) = sqzeta*sqzeta

The routine sets Gpay[i] to its new value Gpay[i] + coef*dGpay[i]
sets *Gstptr to Gst's new value, Gst+coef*dGst, and returns the
calibrated zeta(tEx).

On error, returns -1. Note that zeta must be positive */
/* The fixed legs of the swaptions have the discounted payments
         a1[i] = payment1[i]*df(tPay[i]), i=1,...,n1
         a2[i] = payment2[i]*df(tPay[i]), i=1,...,n2
respectively, and the strike is worth DStart = df(tStart) */
static double FitOneLongShort(
    long            n1,
    double*         a1, /* swaption 1 */
    long            n2,
    double*         a2, /* swaption 2 */
    double          DStart,
    SrtReceiverType recpay, /* PV of strike for both swaptions */
    double*         Gstptr,
    double*         Gpay, /* the existing term structure */
    double          dGst,
    double*         dGpay, /* direction of change */
    double          price1,
    double          price2) /* market prices of swaptions */
{
    SrtCallPutType callput;
    double         Gst, favg, arg, z, dz, cz, dcz;
    double         expfac, olderror, newerror;
    double         fwd1, theta1, stddev1, f1, df1z, df1cz;
    double         fwd2, theta2, stddev2, f2, df2z, df2cz;
    double         value1, value2, derG1, derG2, derSqz1, derSqz2;
    double         ystar1, ystar2;
    double         isign, zeta, sqzeta, coef, deltaSqZ, det;
    long           i, n;
    long           iter, itermax, success;

    /* initialize */
    itermax = 25;
    n       = (n1 > n2 ? n1 : n2);
    callput = SRT_CALL;
    Gst     = *Gstptr;

    /* compute forward value of the swaptions */
    fwd1 = 0.;
    for (i = 1; i <= n1; i++)
        fwd1 = fwd1 + a1[i];
    theta1 = fwd1 - DStart;

    fwd2 = 0.;
    for (i = 1; i <= n2; i++)
        fwd2 = fwd2 + a2[i];
    theta2 = fwd2 - DStart;

    if (recpay == SRT_PAYER) /* need only consider receivers now */
    {
        price1 = price1 + theta1;
        price2 = price2 + theta2;
    }

    /* check for intrinsic (1) */

    if (fabs(price1 - max(theta1, 0)) < 1.0e-08)
    {
        return 1.0e-15;
    }

    if (fabs(price2 - max(theta2, 0)) < 1.0e-08)
    {
        return 1.0e-15;
    }

    if (price1 < max(theta1, 0.0) + 1.e-15)
        return (-1.0); /* cannot fit */
    if (price2 < max(theta2, 0.0) + 1.e-15)
        return (-1.0); /* cannot fit */

    /* get initial guess by linear approximation (neglects convexity) */
    favg    = 0.5 * (fwd1 + DStart);
    stddev1 = favg * LGMNormImpVol(price1 / favg, fwd1 / favg, DStart / favg, 1.0, 1.0, callput);

    favg    = 0.5 * (fwd2 + DStart);
    stddev2 = favg * LGMNormImpVol(price2 / favg, fwd2 / favg, DStart / favg, 1.0, 1.0, callput);

    /* check for intrinsic (2) */

    if (fabs(stddev1) < 1.0e-08 || fabs(stddev2) < 1.0e-08)
    {
        return 1.0e-15;
    }

    if (stddev1 < 0 || stddev2 < 0)
        return (-1.0); /* cannot fit */

    success  = -2;
    olderror = 2.0;
    z        = 0.;
    dz       = 0.;
    cz       = 0.;
    dcz      = 0.;
    for (iter = 1; iter <= itermax && success <= 0; iter++)
    {
        /* calculate f1, f2, and their derivatives */
        f1    = 0.;
        df1z  = 0.;
        df1cz = 0.;
        for (i = 1; i <= n1; i++)
        {
            arg = (Gst - Gpay[i]) * stddev2 * (z + dz) + (dGst - dGpay[i]) * stddev2 * (cz + dcz);
            expfac = exp(-arg * theta1);
            if (fabs(arg * theta1) < 1.e-5)
                f1 = f1 +
                     a1[i] * arg * (1.0 - 0.5 * arg * theta1 + arg * arg * theta1 * theta1 / 6.0);
            else
                f1 = f1 + a1[i] * (1.0 - expfac) / theta1;
            df1z  = df1z + a1[i] * (Gst - Gpay[i]) * stddev2 * expfac;
            df1cz = df1cz + a1[i] * (dGst - dGpay[i]) * stddev2 * expfac;
        }

        f2    = 0.;
        df2z  = 0.;
        df2cz = 0;
        for (i = 1; i <= n2; i++)
        {
            arg = (Gst - Gpay[i]) * stddev1 * (z + dz) + (dGst - dGpay[i]) * stddev1 * (cz + dcz);
            expfac = exp(-arg * theta2);
            if (fabs(arg * theta2) < 1.e-5)
                f2 = f2 +
                     a2[i] * arg * (1.0 - 0.5 * arg * theta2 + arg * arg * theta2 * theta2 / 6.0);
            else
                f2 = f2 + a2[i] * (1.0 - expfac) / theta2;
            df2z  = df2z + a2[i] * (Gst - Gpay[i]) * stddev1 * expfac;
            df2cz = df2cz + a2[i] * (dGst - dGpay[i]) * stddev1 * expfac;
        }
        /* check convergence */
        newerror = (f1 - 1) * (f1 - 1) + (f2 - 1) * (f2 - 1);
        if (newerror < 1.e-10)
            success++;
        /* if improvement, accept step and calculate new Newton step */
        if (iter == 1 || newerror < olderror)
        {
            z        = z + dz;
            cz       = cz + dcz;
            olderror = newerror;
            det      = df1z * df2cz - df1cz * df2z;
            dz       = (df2cz * (1.0 - f1) - df1cz * (1 - f2)) / det;
            dcz      = (df1z * (1.0 - f2) - df2z * (1 - f1)) / det;
        }
        else
        {
            dz  = 0.5 * dz;
            dcz = 0.5 * dcz;
        }
    }
    z  = z + dz;
    cz = cz + dcz;

    coef   = cz / z; /* initial guess */
    sqzeta = stddev1 * stddev2 * z;
    if (coef <= 0.0) /* cannot improve on coef=0 */
        return (sqzeta * sqzeta);

    /* update term structure */
    Gst = Gst + coef * dGst;
    for (i = 1; i <= n; i++)
        Gpay[i] = Gpay[i] + coef * dGpay[i];

    /* Use global Newton to match LGM prices to market prices */
    success = -3;
    isign   = 1.0;
    for (iter = 1; iter <= itermax && success <= 0; iter++)
    {
        value1   = LGMRecVal(&ystar1, n1, a1, DStart, sqzeta, Gpay, Gst) - price1;
        value2   = LGMRecVal(&ystar2, n2, a2, DStart, sqzeta, Gpay, Gst) - price2;
        newerror = value1 * value1 / (price1 * price1) + value2 * value2 / (price2 * price2);
        if (newerror < 1.e-10)
            success++;

        if (iter == 1 || newerror < olderror)
        {
            isign    = 1.0;
            olderror = newerror;
            derG1    = RecGDer(ystar1, n1, a1, DStart, sqzeta, Gpay, Gst, dGst, dGpay);
            derSqz1  = RecSqZetaDer(ystar1, n1, a1, DStart, sqzeta, Gpay, Gst);

            derG2   = RecGDer(ystar2, n2, a2, DStart, sqzeta, Gpay, Gst, dGst, dGpay);
            derSqz2 = RecSqZetaDer(ystar2, n2, a2, DStart, sqzeta, Gpay, Gst);

            det      = derSqz1 * derG2 - derSqz2 * derG1; /* Newton step */
            deltaSqZ = (-value1 * derG2 + value2 * derG1) / det;
            coef     = (value1 * derSqz2 - value2 * derSqz1) / det;
        }
        else
        {
            isign    = -1.0;
            coef     = 0.5 * coef;
            deltaSqZ = 0.5 * deltaSqZ;
        }
        /* update term structure */
        sqzeta = sqzeta + isign * deltaSqZ;
        Gst    = Gst + isign * coef * dGst;
        for (i = 1; i <= n; i++)
            Gpay[i] = Gpay[i] + isign * coef * dGpay[i];
    }

    if (success >= -1)
    {
        *Gstptr = Gst;
        zeta    = sqzeta * sqzeta;
        return (zeta);
    }
    else
        return (-1.0);
}

/*****************************************************/
/****** LGM value and derivatives of swaptions *******/
/*****************************************************/
/* Calculate the LGM value of a receiver with
        a[i] = discounted fixed	leg payments, i=1,...,n
        DStart = discount factor to start date
        sqzeta = sqrt(zeta(t)) at the exercise date
        Gpay[i] = G at i-th pay date, i=1,...,n
        Gst = G(t) at the start date */
double LGMRecVal(
    double* ystar,
    long    n,
    double* a,
    double  DStart, /* receiver swaption */
    double  sqzeta,
    double* Gpay,
    double  Gst) /* term structure */
{
    double price, arg, deltaG;
    long   i;

    /* find value of y where swaption is ATM */
    *ystar = LGMFindystar(n, a, DStart, sqzeta * sqzeta, Gpay, Gst);

    if (sqzeta <= 1.e-12)
        sqzeta = 1.e-12;
    arg   = (*ystar) / sqzeta;
    price = -DStart * LGMsafeNorm(arg);
    for (i = 1; i <= n; i++)
    {
        deltaG = Gst - Gpay[i];
        price  = price + a[i] * LGMsafeNorm(arg + deltaG * sqzeta);
    }

    if (Gst < Gpay[n])
    {
        price = -(DStart + price); /* this should never occur */
        for (i = 1; i <= n; i++)
            price = price + a[i];
    }
    return (price);
}

/*******************************************************************/
/* returns the derivative of the LGM price with respect to sqrt(zeta)
        for the same swaption as above */
/* NOTE: LGMRecVal must be ran first so that we have ystar */
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

/*****************************************************************/
/* returns the derivative of the LGM price with respect to C (at C=0) where
        Gst -> Gst + C*dGst
        Gpay[i] -> Gpay[i] + C*dGpay[i]
for the same swaption as above */
/* NOTE: LGMRecVal must be ran first so that we have ystar */
static double RecGDer(
    double  ystar,
    long    n,
    double* a,
    double  DStart, /* receiver swaption */
    double  sqzeta,
    double* Gpay,
    double  Gst, /* term structure */
    double  dGst,
    double* dGpay) /* direction of differentiation */
{
    double der, arg;
    long   i;

    if (sqzeta <= 1.e-12)
        sqzeta = 1.e-12;

    arg = ystar / sqzeta;
    der = 0.0;
    for (i = 1; i <= n; i++)
        der = der + (dGst - dGpay[i]) * a[i] * LGMsafeGauss(arg + (Gst - Gpay[i]) * sqzeta);

    der = der * sqzeta;

    if (Gst < Gpay[n])
        der = -der; /* this should never occur */

    return (der);
}

/**************************************************************/
/* finds ystar, the state variable at which a swaption is ATM */
/* Same swaption as above */
double LGMFindystar(long n, double* a, double DStart, double zeta, double* Gpay, double Gst)
{
    double y, term, func, der, dfunc, dy, target;
    long   iter, itermax, i;
    long   success;

    itermax = 25;
    success = -2;
    target  = log(DStart);
    y       = 0.; /* initial guess */

    for (iter = 1; iter < itermax && success <= 0; iter++)
    {
        func = 0.;
        der  = 0.;
        for (i = 1; i <= n; i++)
        {
            term = a[i] * exp(-(Gst - Gpay[i]) * (y + 0.5 * (Gst - Gpay[i]) * zeta));
            func = func + term;
            der  = der - term * (Gst - Gpay[i]);
        }
        dfunc = target - log(func);
        dy    = dfunc * func / der;
        y     = y + dy;
        if (fabs(dy) < 1.e-9 || fabs(dfunc) < 1.e-9)
            success++;
    }
    return (y);
}

/*****************************************************************************/
/* Create a term structure, calculate G(t) from fixed kappa, fill in numG,
Gdates, G values */
/*****************************************************************************/
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

/*****************************************************************************/
/* Create a term structure, check out G(t) from user input, fill in numG,
Gdates, G values */
/*****************************************************************************/
static LGM_TSPtr CreateTSwithGs(
    long       numG,  /* number of values of G */
    Date*      Gdate, /* construct G on these dates */
    double*    G,     /* values of G(t) */
    LGMCalSet* CSPtr) /* Reference instruments */
{
    LGM_TS* tsPtr;
    long    i, nEx;
    double  slope, dt, isign;

    nEx = CSPtr->nEx;
    if (numG < 2 || Gdate == NULL || G == NULL) /* too little data */
        return (NULL);

    /* Create term structure, fill in numG and Gdates */
    tsPtr = LGMCreateLGM_TS(nEx + 1, numG);
    if (tsPtr == NULL)
        return (NULL);

    /* fill in G dates */
    isign = 1.0;
    if (G[0] < G[numG - 1])
        isign = -1.;
    tsPtr->numG = numG;
    for (i = 0; i < numG; i++)
    {
        tsPtr->Gdate[i] = Gdate[i];
        tsPtr->G[i]     = isign * (G[i] - G[numG - 1]);
    }

    /* ensure that G(t) is decreasing */
    for (i = 1; i < numG; i++)
    {
        if (tsPtr->G[i] > tsPtr->G[i - 1])
            tsPtr->G[i] = tsPtr->G[i - 1];
    }

    /* check out Gdates and Gvalues */
    for (i = 1; i < numG; i++) /* check date order */
    {
        if (tsPtr->Gdate[i] < tsPtr->Gdate[i - 1])
        {
            LGMFreeLGM_TS(&tsPtr);
            return (NULL);
        }
    }

    dt    = (double)(tsPtr->Gdate[numG - 1] - tsPtr->Gdate[0]);
    slope = tsPtr->G[0] * 365.0 / dt; /* re-scale */
    if (fabs(slope) < 1.e-8)          /* G values too small? */
    {
        LGMFreeLGM_TS(&tsPtr);
        return (NULL);
    }

    for (i = 0; i < numG; i++)
        tsPtr->G[i] = tsPtr->G[i] / slope;

    for (i = 1; i < numG; i++)
    {
        dt    = (double)(tsPtr->Gdate[i] - tsPtr->Gdate[i - 1]); /* G values too ragged? */
        slope = (tsPtr->G[i - 1] - tsPtr->G[i]) * 365.0 / dt;
        if (fabs(slope) > 10.0)
        {
            LGMFreeLGM_TS(&tsPtr);
            return (NULL);
        }
    }
    return (tsPtr);
}

/*****************************************************************************/
/* Create a term structure, compute zeta from constant sigma at the exercise
 dates, and fill in numZ, zdate[], and zeta[]   */
/*****************************************************************************/
static LGM_TSPtr CreateTSwithSigma(LGMCalSet* CSPtr) /* Reference instruments */
{
    LGM_TS* tsPtr;
    long    j, nEx, numG;
    Date    tNow;
    double  sigma;

    nEx  = CSPtr->nEx;
    numG = nEx + 1;
    tNow = CSPtr->tNow;

    /* Create term structure */
    tsPtr = LGMCreateLGM_TS(nEx + 1, numG);
    if (tsPtr == NULL)
        return (NULL);

    sigma           = 0.01;
    tsPtr->numZ     = nEx + 1;
    tsPtr->zdate[0] = tNow;
    tsPtr->zeta[0]  = 0.0;

    for (j = 1; j <= nEx; j++)
    {
        tsPtr->zdate[j] = CSPtr->tEx[j];
        tsPtr->zeta[j]  = sigma * sigma * ((double)(tsPtr->zdate[j] - tNow)) / 365.0;
    }
    return (tsPtr);
}

/*****************************************************************************/
/* Create a term structure, check out input zeta(t) and fill in numZ,
   zdate[], and zeta[]   */
/*****************************************************************************/
static LGM_TSPtr CreateTSfromZetas(
    long       nZ,    /* number of zeta values */
    Date*      Zdate, /* zeta given on these dates */
    double*    Zval,  /* zeta values */
    LGMCalSet* CSPtr) /* Reference instruments */
{
    LGM_TS* tsPtr;
    long    j0, j, numG, numZ;
    Date    tNow;
    double  z0, slope, dt;

    numG = CSPtr->nEx + 1;
    tNow = CSPtr->tNow;

    if (nZ < 2 || Zdate == NULL || Zval == NULL) /* too few values */
        return (NULL);

    /* check dates */
    for (j = 1; j < nZ; j++)
    {
        if (Zdate[j] <= Zdate[j - 1])
            return (NULL);
    }

    /* find all dates after tNow */
    for (j0 = 0; j0 < nZ && Zdate[j0] <= tNow; j0++)
        ;
    numZ = nZ + 1 - j0;
    if (numZ < 2) /* no date after tNow */
        return (NULL);

    /* create term structure */
    tsPtr = LGMCreateLGM_TS(numZ, numG);
    if (tsPtr == NULL)
        return (NULL);

    /* interpolate to get Zval(tNow) */
    if (j0 > 0)
        slope = (Zval[j0] - Zval[j0 - 1]) / ((double)(Zdate[j0] - Zdate[j0 - 1]));
    else
        slope = (Zval[1] - Zval[0]) / ((double)(Zdate[1] - Zdate[0]));

    z0 = Zval[j0] + slope * ((double)(tNow - Zdate[j0]));

    /* fill in the zeta dates & values, change to feasible set if needed */
    tsPtr->numZ     = numZ;
    tsPtr->zdate[0] = tNow;
    tsPtr->zeta[0]  = 0.0;

    for (j = j0; j < nZ; j++)
    {
        tsPtr->zdate[j + 1 - j0] = Zdate[j];
        tsPtr->zeta[j + 1 - j0]  = Zval[j] - z0;
        if (tsPtr->zeta[j + 1 - j0] <= tsPtr->zeta[j - j0])
            tsPtr->zeta[j + 1 - j0] = tsPtr->zeta[j - j0];
    }

    /* Re-scale term structure */
    dt    = (double)(tsPtr->zdate[numZ - 1] - tNow);
    slope = (tsPtr->zeta[numZ - 1]) * 365.0 / dt;
    if (fabs(slope) < 1.e-16)
    {
        LGMFreeLGM_TS(&tsPtr); /* zeta values too small */
        return (NULL);
    }

    for (j = 0; j < numZ; j++)
        tsPtr->zeta[j] = 0.0001 * (tsPtr->zeta[j]) / slope;

    return (tsPtr);
}
