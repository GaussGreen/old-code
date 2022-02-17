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
#include "OPFNCTNS.H>
#include "SRT_H_ALL.H>
#include "math.h"
#include "srt_h_lgmUStypes.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"
#define CEVvolShift 1e-2

/****************************************************************************/
/* Prototypes of private functions */
/****************************************************************************/
/* Check n  , the array of exercise dates and the array of forward values
Return the first i with t[i]>tNow  , and the last t with t[i]<tLast */
static LGMErr checkinputarrays(long n, Date tNow, Date tLast, Date *t,
                               double *RatArr, LGMCalParm *CR, long *firstiPtr,
                               long *lastiPtr);

/* Routine to construct the common parts of all reference instrumetns */
static LGMErr
RefDealsCommon(Date tNow, long nEx, /* evaluation date  , number of exercises */
               Date tLast, Date *TauArr, /* last pay date  , exer dates */
               LGMCalParm *CalReq,       /* Methodology */
               LGMMarkConv *conv,
               String ycname,     /* market conventions  , discount curve */
               LGMCalSet *CSPtr); /* Output */

/* Routine to construct the long swaptions */
/* Note: RefDealsCommon must be called first */
static LGMErr RefDealsLong(
    long nArr, Date *TauArr,
    double *FVArr,                    /* exercise date and relative PV arrays */
    LGMMarkConv *conv, String ycname, /* market conventions  , discount curve */
    LGMCalParm *CalReq,               /* calibration methodology */
    LGMDealType DealType, void *dealPtr, /* deal type and params */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* function to provide exponent (beta) to interpret vols */
    LGMCalSet *CSPtr); /* Output */

/* Routine to construct the short swaptions */
/* Note: RefDealslong must be called first */
static LGMErr RefDealsShort(
    long nArr, Date *TauArr,
    double *FVArr,                    /* exercise date and relative PV arrays */
    LGMMarkConv *conv, String ycname, /* market conventions  , discount curve */
    LGMCalParm *CalReq,               /* calibration methodology */
    LGMDealType DealType, void *dealPtr, /* deal type and params */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* function to provide exponent (beta) to interpret vols */
    LGMCalSet *CSPtr); /* Output */

/* Routine to construct the caplets */
/* Note: RefDealsCommon must be called first */
static LGMErr
RefDealsCap(LGMDealType DealType, void *DealPtr, long nArr, Date *TauArr,
            double *FVArr, LGMCalParm *CalReq,
            LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double *),
            LGMErr (*GetBeta)(Date, Date, double *), String YcName,
            LGMMarkConv *conv, LGMCalSet *CSPtr); /*OUTPUT*/

/* Routine to construct the "1 into k" swaptions */
/* Note: RefDealsCommon must be called first */
static LGMErr RefDeals1intok(
    long nArr, Date *TauArr,
    double *FVArr,                    /* exercise date and relative PV arrays */
    LGMMarkConv *conv, String ycname, /* market conventions  , discount curve */
    LGMCalParm *CalReq,               /* calibration methodology */
    LGMDealType DealType, void *dealPtr, /* deal type and params */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* function to provide exponent (beta) to interpret vols */
    LGMCalSet *CSPtr); /* Output */

/* Routine to find the fixed rates of the reference instruments */
static LGMErr FindFixRateRefInst(
    double Rsw, /* swap rate */
    double Lvl, /* level  , adjusted for gearing */
    double DStart, double DEnd,
    double DLast,                     /* df's to start date and end dates */
    Date tEx, Date tStart, Date tEnd, /* dates of refence swaption */
    long islong,                      /* 1=use Rdata1  , 0=use Rdata2 */
    long nArr, double *FVarr,
    Date *TauArr,          /* PV ratio & exer dates of orig deal */
    Date tNow, Date tLast, /* end date of original deal */
    LGMMarkConv *conv,     /* market conventions */
    LGMCalParm *CalReq,    /* calibration request */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *),              /* gets swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double *), /* gets exp (beta) for vols */
    double *rfix, /* returns fix rate of ref.swaption */
    double *std, LGMCalSet *CSPtr);

/* Use linear interpolation (flat at ends) to obtain Rdata at thedate */
static double LinearFlatFlat(Date thedate, long n, Date *t, double *Rdata);

/* Use linear interpolation (linear at ends) to find the PV ratio at thedate
with FVArr=1 at tLast */
static double InterpPVratio(Date thedate, long n, double *FV, Date *Tau,
                            Date tLast);

static double find_gear(Date thedate, Date *d, double *gear, int n);

/****************************************************************************/
/**** Comments ****/
/**** HOLIDAYS ***/
/* Wherever we need to use the add_unit function  ,
the currency code ccy is available */

/****************************************************************************/
/* The following routines construct the reference instruments */
/* Main routine to construct the reference instruments */
LGMErr LGMMakeRefDeals(
    LGMDealType DealType, void *dealptr, Date tNow, /*EVALUATION DATE */
    Date tlast, long nArr, /*NUMBER OF EXERCISE DATES IN DEAL */
    Date *TauArr, /* [*  ,1  ,2  ,...  ,nEx] ARRAY OF EXERCISE DATES AFTER TODAY
                   */
    double *RatArr, /* [*  ,1  ,2  ,...  ,nEx] RATIO OF PV OF FIXED LEG TO PV OF
                       STRIKE */
    String YcName,  /* DISCOUNT CURVE */
    LGMMarkConv *conv,  /* MARKET PLACE CONVENTIONS */
    LGMCalParm *CalReq, /* METHODS AND DATA FOR CHOOSING STRIKES*/
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* SWAPTION CAP/VOL */
    LGMErr (*GetBeta)(Date, Date,
                      double *), /* EXPONENT BETA TO INTERPRET VOLS */
    LGMCalSet **RefDealsPtrPtr)  /* REFERENCE INSTRUMENTS (OUTPUT)*/
{
  LGMErr error;
  LGMCalSet *RefDealsPtr;
  int makelong;   /* 1=need long k into n-k swaptions */
  int makeshort;  /* 1=need short k into 1 swaptions */
  int makecap;    /* 1=need long k into 3 mo caplets */
  int make1intok; /* 1=need fixed exercise date "1 into k" swaptions */
  long firsti, lasti, nEx;
  Date *tExArr;
  double *FVArr;

  /* initialize */
  error = NULL;
  *RefDealsPtrPtr = NULL;

  /* check the arrays TauArr  , RatArr and find first date in TauArr
  after tNow and the last date before tLast */

  error = checkinputarrays(nArr, tNow, tlast, TauArr, RatArr, CalReq, &firsti,
                           &lasti);
  if (error != NULL)
    return (error);

  RefDealsPtr = (LGMCalSet *)srt_calloc(1, sizeof(LGMCalSet));
  if (RefDealsPtr == NULL)
    return ("alloc failed in MakeRefDeals");

  nEx = lasti - firsti + 1;
  tExArr =
      &(TauArr[firsti - 1]); /* tExArr[*  ,1  ,...  ,nEx] are relevent dates */
  FVArr = &(RatArr[firsti -
                   1]); /* FVArr[*  ,1  ,...  ,nEx] are relevent fwd values */

  /* Create reference instruments */
  /* Set tNow  , n  , nPay  , the long swaption pay dates tPay  , the coverages
     cvgpay  , and the dis factors Dpay from tNow to the pay dates. Set nEx  ,
     and common exercise dates (tEx)  , start dates (tStart)  ,
          and df to start date (DStart) */

  error = RefDealsCommon(tNow, nEx, tlast, tExArr, CalReq, conv, YcName,
                         RefDealsPtr);
  if (error != NULL) {
    LGMFreeCalSet(&RefDealsPtr);
    return (error);
  }

  /* determine which deals are needed */
  makelong = 0;
  makeshort = 0;
  makecap = 0;
  make1intok = 0;

  if (CalReq->LGMOneTwoFactor == 1) {
    if (CalReq->usecaps == 1) /* All methods except Tenor & Diagonal */
      makecap = 1;            /* use either long "k into n-k" swaptions */
    else                      /* or the "k into 3mo" caplets			*/
      makelong = 1;

    if (CalReq->calmeth ==
        FixExp)       /* The fixed expiry method also uses the "1 into k" */
      make1intok = 1; /* swaptions to find G(t) */

    if (CalReq->calmeth == TenorAndDiag) /* TenorAndDiag uses long "k into n-k"
                                            swaptions for zeta(t) */
    {
      makelong = 1; /* It uses long "k into n-k" swaptions and either */
      if (CalReq->usecaps == 1)
        makecap = 1; /* "k into 3mo caplets  , or */
      else
        makeshort = 1; /* short "k into 1" swaptions */
    }
  } else if (CalReq->LGMOneTwoFactor == 2) {
    makelong = 1; /* It uses long "k into n-k" swaptions and either */
    makecap = 1;  /* "k into 3mo caplets  , or */
  }

  /* create long swaptions */
  if (makelong != 0) {
    error = RefDealsLong(nEx, tExArr, FVArr, conv, YcName, CalReq, DealType,
                         dealptr, GetVol, GetBeta, RefDealsPtr);
    if (error != NULL) {
      LGMFreeCalSet(&RefDealsPtr);
      return (error);
    }
  }

  /* create the short swaptions */
  if (makeshort != 0) {
    error = RefDealsShort(nEx, tExArr, FVArr, conv, YcName, CalReq, DealType,
                          dealptr, GetVol, GetBeta, RefDealsPtr);
    if (error != NULL) {
      LGMFreeCalSet(&RefDealsPtr);
      return (error);
    }
  }

  /* create caplets */
  if (makecap != 0) {
    error = RefDealsCap(DealType, dealptr, nEx, tExArr, FVArr, CalReq, GetVol,
                        GetBeta, YcName, conv, RefDealsPtr);

    if (error != NULL) {
      LGMFreeCalSet(&RefDealsPtr);
      return (error);
    }
  }

  /* create the "1 into k swaptions" */
  if (make1intok != 0) {
    error = RefDeals1intok(nEx, tExArr, FVArr, conv, YcName, CalReq, DealType,
                           dealptr, GetVol, GetBeta, RefDealsPtr);
    if (error != NULL) {
      LGMFreeCalSet(&RefDealsPtr);
      return (error);
    }
  }

  *RefDealsPtrPtr = RefDealsPtr;
  return (error);
}

/*******************************************/
/* Check n  , the array of exercise dates  , and (if it's needed)
the array of forward values of the fixed leg
Return the first i with t[i]>tNow  , and the last t with t[i]<tLast */
static LGMErr checkinputarrays(long n, Date tNow, Date tLast, Date *t,
                               double *RatArr, LGMCalParm *CR, long *firstiPtr,
                               long *lastiPtr) {
  long i, needratios;
  long firsti, lasti;

  needratios = 0;
  if (CR->Rmeth == IRR ||
      CR->Rmeth == dIRR) /* these methods choose the strikes Rfix */
    needratios = 1;      /* based on the fwd value ratios of the
                                         /* underlying deal */
  if (n < 1 || t == NULL)
    return ("no exer dates for calibration");

  for (i = 1; i < n; i++) {
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

  if (needratios > 0) {
    if (RatArr == NULL)
      return ("no fwd values in calibrator");
    for (i = firsti; i <= lasti; i++) {
      if (RatArr[i] < 0.25 || RatArr[i] > 4.0)
        return ("bad fwd values");
    }
  }
  return (NULL);
}

/*******************************************/
/* Routine to constuct the common elements of all the reference instruments */
/* First construct the long fixed leg tPay[]  , which runs from the last paydate
tLast = tPay[n] to the spot of the first exercise date  , tPay[0]. For
simplicity  , also extend the pay dates to either tPay[n+1] (for annual
compounding) or tPay[n+2] (for semiannual compounding) to ensure we can fit in a
minimum length swap whose start date is before tPay[n]. Construct the fixed leg
pay dates tPay[i]  , i = 1  ,...  ,n. the discount factor Dpay[i] from tNow to
the pay date the coverage cvgpay[i] between tPay[i-1] and tPay[i]. If ResExer =
1  , then set the common exercise dates to be tEx[i]  , and set the start dates
tStart[i] to be spot of tEx[i]. If use starts=0  , then let the start dates be
tPay[i]  , i=1 to n-1  , and set tEx[i] to be lag business days earlier. In
either case set DStart[i] to be the discount factor from tNow to tStart.
*/
static LGMErr
RefDealsCommon(Date tNow, long nEx, /* evaluation date  , number of exercises */
               Date tLast, Date *TauArr, /* last pay date  , exer dates */
               LGMCalParm *CalReq,       /* Methodology */
               LGMMarkConv *conv,
               String ycname,    /* market conventions  , discount curve */
               LGMCalSet *CSPtr) /* Output */
{
  Date *tEx, *tStart, *tPay;
  long extraperiods, j, n, nPay;
  int ResExer;

  /* initialize */
  CSPtr->longflag = 0; /* flag all reference instruments as not existing */
  CSPtr->shortflag = 0;
  CSPtr->capflag = 0;
  CSPtr->fixflag = 0;
  CSPtr->tNow = tNow;
  ResExer = CalReq->respectExer;

  /* Set up fix leg for a standard swap with exercise date tfirstEx and last
  pay date tLast. For safety  , add extraperiods periods to end of leg */
  /* Find fixed leg pay dates tPay[1  ,...  ,nPay]; start date is tPay[0] */
  extraperiods = conv->minswap;
  tPay = LGMFixLegSched(TauArr[1], tLast, extraperiods, conv, &nPay);
  n = nPay - extraperiods;

  CSPtr->n = n;
  CSPtr->nPay = nPay;
  CSPtr->tPay = tPay;

  if (tPay == NULL || n < 1)
    return ("failure in RefDealsCommon");

  CSPtr->cvgpay = LGMFixLegCvg(nPay, tPay, conv); /* coverages */
  CSPtr->Dpay =
      LGMFixLegDF(nPay, tPay, tNow, ycname); /* discount factors to pay dates */

  if (CSPtr->cvgpay == NULL || CSPtr->Dpay == NULL)
    return ("allocation failed in RefDealsCommon");

  /* Allocate space for common exer dates  , start dates  , and dis factors to
   * start dates */
  if (ResExer == 0)
    CSPtr->nEx = n;
  else
    CSPtr->nEx = nEx;

  CSPtr->tEx = tEx = (Date *)srt_calloc(CSPtr->nEx + 1, sizeof(Date));
  CSPtr->tStart = tStart = (Date *)srt_calloc(CSPtr->nEx + 1, sizeof(Date));
  CSPtr->DStart = (double *)srt_calloc(CSPtr->nEx + 1, sizeof(double));

  if (tEx == NULL || tStart == NULL || CSPtr->DStart == NULL)
    return ("allocation failed in RefDealsCommon");

  /* Fill in common exercise dates  , start dates  , and dis factors to start
   * dates */
  if (ResExer != 0) {
    for (j = 1; j <= nEx; j++) {
      tEx[j] = TauArr[j];
      tStart[j] = add_unit(tEx[j], conv->lag, SRT_BDAY, SUCCEEDING);
    }
  } else {
    for (j = 1; j <= n; j++) {
      tStart[j] = tPay[j - 1];
      tEx[j] = add_unit(tStart[j], -conv->lag, SRT_BDAY, SUCCEEDING);
    }
    tEx[1] = TauArr[1]; /* shouldn't make a difference */
  }

  for (j = 1; j <= CSPtr->nEx; j++) {
    CSPtr->DStart[j] = swp_f_df(tNow, tStart[j], ycname);
    if (CSPtr->DStart[j] == SRT_DF_ERROR)
      return ("no discount factor");
  }

  return (NULL); /* done */
}

/*******************************************/
/* Routine to construct the long swaptions */
/* Construct the long swaption j. j=1  ,2  , ...  , nEx. Long swaption j has
exercise date tEx[j]  , start date tStart[j]  , and fixed leg pay dates
tPay[i]  , i=ifirst[j] to i=nlong[j].
The discount factor to the start date is DStart[j]  , and
the discount factors to the pay dates tPay[i] are Dpay[i].
The coverage between start and the first pay date is cvgfirst[j]  , and
the coverages for the period between tPay[i-1] and tPay[i] is cvgpay[i] for
i>ifirst[j]
*/
/* Note: RefDealsCommon must be called first to set nEx  , n  , nPay  , tEx[]  ,
   tStart[]  , tPay  , cvgpay  , Dpay  , n  , and nPay */

static LGMErr RefDealsLong(
    long nArr, Date *TauArr,
    double *FVArr,                    /* exercise date and relative PV arrays */
    LGMMarkConv *conv, String ycname, /* market conventions  , discount curve */
    LGMCalParm *CalReq,               /* calibration methodology */
    LGMDealType DealType, void *dealPtr, /* For EMKi Method */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* function to provide exponent (beta) to interpret vols */
    LGMCalSet *CSPtr) /* Output */
{
  long i, j, iend, islong;
  long fullper;
  double reflevel, level, DLast, refsw, rfix, std;
  Date tNow, tLast;
  SrtCallInvFltPtr cinf;
  SrtCallTimeSwapPtr ctimeswap;
  SrtSimMidAt *MidAt;
  LGMErr error;
  long nEx, n, *ifirst, *nlong;              /* aliases to simplify code */
  Date *tEx, *tStart, *tPay;                 /* aliases to simplify code */
  double *cvgpay, *Dpay, *cvgfirst, *DStart; /* aliases to simplify code */

  /* initialize */
  error = NULL;           /* aliases to simplify code */
  nEx = CSPtr->nEx;       /* aliases to simplify code */
  n = CSPtr->n;           /* aliases to simplify code */
  tNow = CSPtr->tNow;     /* aliases to simplify code */
  tEx = CSPtr->tEx;       /* aliases to simplify code */
  tStart = CSPtr->tStart; /* aliases to simplify code */
  DStart = CSPtr->DStart; /* aliases to simplify code */
  cvgpay = CSPtr->cvgpay; /* aliases to simplify code */
  tPay = CSPtr->tPay;     /* aliases to simplify code */
  Dpay = CSPtr->Dpay;     /* aliases to simplify code */

  CSPtr->longflag = 1; /* flag long swaptions as existing */
  fullper = 155; /* if period has more than fullper days  , count it as full*/
  if (conv->sfreq ==
      SRT_ANNUAL) /* period for fulfilling minimum swap length requirement */
    fullper = 310;

  if (DealType == CallInvFlt) {
    cinf = (SrtCallInvFltPtr)dealPtr;
  } else if (DealType == SimpleMidAt) {
    MidAt = (SrtSimMidAt *)dealPtr;
  } else if (DealType == CallTimeSwap) {
    ctimeswap = (SrtCallTimeSwapPtr)dealPtr;
  }

  /* Allocate space */
  CSPtr->ifirst = (long *)srt_calloc(nEx + 1, sizeof(long));
  CSPtr->nlong = (long *)srt_calloc(nEx + 1, sizeof(long));
  CSPtr->cvgfirst = (double *)srt_calloc(nEx + 1, sizeof(double));
  CSPtr->Rflong = (double *)srt_calloc(nEx + 1, sizeof(double));
  CSPtr->Vlong = (double *)srt_calloc(nEx + 1, sizeof(double));

  CSPtr->VegaLong = (double *)srt_calloc(nEx + 1, sizeof(double));
  CSPtr->StDLong = (double *)srt_calloc(nEx + 1, sizeof(double));
  CSPtr->CEVLong = (double *)srt_calloc(nEx + 1, sizeof(double));
  CSPtr->FrLong = (double *)srt_calloc(nEx + 1, sizeof(double));
  CSPtr->BetaLong = (double *)srt_calloc(nEx + 1, sizeof(double));

  if (CSPtr->ifirst == NULL || CSPtr->nlong == NULL ||
      CSPtr->cvgfirst == NULL || CSPtr->Rflong == NULL ||
      CSPtr->Vlong == NULL || CSPtr->CEVLong == NULL ||
      CSPtr->BetaLong == NULL || CSPtr->FrLong == NULL)

    return ("allocation failed in RefDealsLong");

  ifirst = CSPtr->ifirst;     /* aliases to simplify code */
  nlong = CSPtr->nlong;       /* aliases to simplify code */
  cvgfirst = CSPtr->cvgfirst; /* aliases to simplify code */

  /* for each tEx[j]  , find the first pay date tPay[i] which has tPay[i]>tEx[j]
   */
  /* if a pay date tPay[i] lies between tEx[j] and tStart[j]  , change tPay[i]
   * to tStart[j] */
  i = 1;
  for (j = 1; j <= nEx; j++) {
    for (; i <= n && tPay[i] < tEx[j]; i++)
      ;
    if (tPay[i] <= tStart[j]) {
      tPay[i] = tStart[j];
      i++;
    }
    ifirst[j] = i; /* first pay date strictly after tStart[j] */
    if (i > n && j > 1) {
      nEx =
          j -
          1; /* eliminate exercise dates whose start date occurs after tLast */
      CSPtr->nEx = nEx;
    }
  }

  /* find coverage from start to first pay date for each swaption */
  for (j = 1; j <= nEx; j++) {
    i = ifirst[j];
    cvgfirst[j] = coverage(tStart[j], tPay[i], conv->sbasis);
  }

  /* last pay date for long swaptions */
  for (j = 1; j <= nEx; j++) {
    i = ifirst[j];
    if (tPay[i] - tStart[j] > fullper)
      nlong[j] =
          i + conv->minswap - 1; /* start to ifirst counts as full period */
    else
      nlong[j] =
          i + conv->minswap; /* start to ifirst too short to be full period */

    if (nlong[j] < n) /* long swaption ends at paydate n  , unless can't */
      nlong[j] = n;   /* fit minimum length swaption in */
  }

  /* Find fixed rate for long swaptions */
  tLast = tPay[n];
  DLast = Dpay[n];

  /* Update Midat Structure for special methods of calibration */
  if (CalReq->LGMOneTwoFactor == 2 && DealType == SimpleMidAt &&
      ((CalReq->Rmeth == EMK1) || (CalReq->Rmeth == EMK2))) {
    error = update_autocal_midat_struct(MidAt, ycname, GetVol, GetBeta);
    if (error)
      return error;
  }

  for (j = 1; j <= nEx; j++) {
    i = ifirst[j];
    iend = nlong[j];

    reflevel = cvgfirst[j] * Dpay[i]; /* first period may be broken */
    i++;
    for (; i <= iend; i++) {
      reflevel = reflevel + cvgpay[i] * Dpay[i];
    }
    refsw = (DStart[j] - Dpay[iend]) / reflevel; /* cash swap rate */

    if (DealType == CallInvFlt) {
      i = ifirst[j];
      level = cvgfirst[j] * Dpay[i] *
              (1.0 + find_gear(tPay[i], cinf->tCpnPay, cinf->gear, cinf->nCpn));
      i++;
      for (; i <= iend; i++) {
        level +=
            cvgpay[i] * Dpay[i] *
            (1.0 + find_gear(tPay[i], cinf->tCpnPay, cinf->gear, cinf->nCpn));
      }
      level /= DStart[j];
    } else {
      level = reflevel / DStart[j];
    }

    islong = 1; /* use first set of Rdata if relevant */

    /* CALIBRATION AT THE MONEY SWAPTIONS FOR CALLABLE TIME SWAPS. */

    if (DealType == CallTimeSwap && ctimeswap->calibrationmethod == 1) {
      rfix = refsw;
      std = 0;
    } else /* CALIBRATION AT THE EQUIVALENT STRIKE */
    {
      error = FindFixRateRefInst(
          refsw, /* swap rate of ref swaption */
          level, DStart[j], Dpay[iend],
          DLast, /* df's for ref start and end  , actual end */
          tEx[j], tStart[j],
          tPay[iend], /* exer  , start  , and end date of ref swaption */
          islong,     /* which set of Rdata to use for strikes */
          nArr, FVArr, TauArr, tNow,
          tLast, /* PV ratio  , exer  , and end date of orig deal */
          conv, CalReq, GetVol,
          GetBeta, /* Calibration method  , market vol and smile */
          &rfix,   /* strike of the reference swaption */
          &std, CSPtr);
    }

    if (error != NULL)
      return (error);

    CSPtr->Rflong[j] = rfix;
    CSPtr->StDLong[j] = std;
  }

  return error;
}

/*******************************************/
/*******************************************/
/* Routine to construct the short swaptions */
/* Short swaption j has the same exercise dates and start dates as the long
swaptions. The fixed leg date pay dates are tPay[i]  , i=ifirst[j] to
i=nshort[j]. The discount factors to the pay dates tPay[i] are Dpay[i] like the
long swaption  , the coverage between start and the first pay date is
cvgfirst[j] like the long swaption  , and the coverages for the period between
tPay[i-1] and tPay[i] is cvgpay[i] for i>ifirst[j]. The only new variables we
need are nshort[j]. */
/* Note: RefDealslong must be called first to set n  , nPay  , tPay[]  ,
 * ifirst[] */

static LGMErr RefDealsShort(
    long nArr, Date *TauArr,
    double *FVArr,                    /* exercise date and relative PV arrays */
    LGMMarkConv *conv, String ycname, /* market conventions  , discount curve */
    LGMCalParm *CalReq,               /* calibration methodology */
    LGMDealType DealType, void *dealPtr, /* deal type and params */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* function to provide exponent (beta) to interpret vols */
    LGMCalSet *CSPtr) /* Output */
{
  long i, j, iend, islong;
  long fullper;
  double reflevel, level, DLast, refsw, rfix;
  Date tNow, tLast;
  SrtCallInvFltPtr cinf;
  SrtSimMidAt *MidAt;
  LGMErr error;
  long nEx, *ifirst, *nshort;                /* aliases to simplify code */
  Date *tEx, *tStart, *tPay;                 /* aliases to simplify code */
  double *cvgpay, *Dpay, *cvgfirst, *DStart; /* aliases to simplify code */

  /* initialize */
  error = NULL;
  CSPtr->shortflag = 1; /* flag short swaptions as existing */

  nEx = CSPtr->nEx;           /* aliases to simplify code */
  tNow = CSPtr->tNow;         /* aliases to simplify code */
  tEx = CSPtr->tEx;           /* aliases to simplify code */
  ifirst = CSPtr->ifirst;     /* aliases to simplify code */
  cvgfirst = CSPtr->cvgfirst; /* aliases to simplify code */
  tStart = CSPtr->tStart;     /* aliases to simplify code */
  DStart = CSPtr->DStart;     /* aliases to simplify code */
  cvgpay = CSPtr->cvgpay;     /* aliases to simplify code */
  tPay = CSPtr->tPay;         /* aliases to simplify code */
  Dpay = CSPtr->Dpay;         /* aliases to simplify code */

  fullper = 155; /* if period has more than fullper days  , count it as full*/
  if (conv->sfreq ==
      SRT_ANNUAL) /* period for fulfilling minimum swap length requirement */
    fullper = 310;

  if (DealType == CallInvFlt) {
    cinf = (SrtCallInvFltPtr)dealPtr;
  } else if (DealType == SimpleMidAt) {
    MidAt = (SrtSimMidAt *)dealPtr;
  }

  /* Allocate storage for nshort[j] */
  CSPtr->nshort = (long *)srt_calloc(nEx + 1, sizeof(long));
  CSPtr->Rfshort = (double *)srt_calloc(nEx + 1, sizeof(double));
  CSPtr->Vshort = (double *)srt_calloc(nEx + 1, sizeof(double));

  if (CSPtr->nshort == NULL || CSPtr->Rfshort == NULL || CSPtr->Vshort == NULL)
    return ("allocation failed in RefDealsShort");

  nshort = CSPtr->nshort; /* aliases to simplify code */

  /* Find end period nshort[j] */
  for (j = 1; j <= nEx; j++) {
    i = ifirst[j];
    if (tPay[i] - tStart[j] > fullper)
      nshort[j] =
          i + conv->minswap - 1; /* start to ifirst counts as full period */
    else
      nshort[j] =
          i + conv->minswap; /* start to ifirst too short to be full period */
  }

  /* Find fixed rates Rsh of short swaptions */
  tLast = tPay[CSPtr->n];
  DLast = Dpay[CSPtr->n];

  /* Update Midat Structure for special methods of calibration */
  if (CalReq->LGMOneTwoFactor == 2 && DealType == SimpleMidAt &&
      ((CalReq->Rmeth == EMK1) || (CalReq->Rmeth == EMK2))) {
    error = update_autocal_midat_struct(MidAt, ycname, GetVol, GetBeta);
    if (error)
      return error;
  }

  for (j = 1; j <= nEx; j++) { /* for each short swaption used */
    i = ifirst[j];             /* compute forward level and swap rate */
    iend = nshort[j]; /* cap dates are i = ifirst[j] to iend=nshort[j] */

    reflevel = cvgfirst[j] * Dpay[i]; /* first period may be broken */
    i++;
    for (; i <= iend; i++) {
      reflevel = reflevel + cvgpay[i] * Dpay[i];
    }
    refsw = (DStart[j] - Dpay[iend]) / reflevel; /* cash swap rate */

    if (DealType == CallInvFlt) {
      i = ifirst[j];
      level = cvgfirst[j] * Dpay[i] *
              (1.0 + find_gear(tPay[i], cinf->tCpnPay, cinf->gear, cinf->nCpn));
      i++;
      for (; i <= iend; i++) {
        level +=
            cvgpay[i] * Dpay[i] *
            (1.0 + find_gear(tPay[i], cinf->tCpnPay, cinf->gear, cinf->nCpn));
      }
      level /= DStart[j];
    } else {
      level = reflevel / DStart[j];
    }

    islong = 0; /* use second set of Rdata if relevant */

    error = FindFixRateRefInst(
        refsw, /* swap rate of ref swaption*/
        level, DStart[j], Dpay[iend],
        DLast, /* df's for ref start and end  , actual end */
        tEx[j], tStart[j],
        tPay[iend], /* exer  , start  , and end date of ref swaption */
        islong,     /* which set of Rdata for strikes */
        nArr, FVArr, TauArr, tNow,
        tLast, /* PV ratio  , exer  , and end date of orig deal */
        conv, CalReq, GetVol,
        GetBeta, /* Calibration method  , market vol and smile */
        &rfix,   /* strike of the reference swaption */
        NULL, CSPtr);

    if (error != NULL)
      return (error);

    CSPtr->Rfshort[j] = rfix;
  }

  return error;
}

/*******************************************/
/* Routine to construct the caplets */
/* Construct the caplets j  , j=1  ,2  , ...  , nEx. Caplet j has
the fixing date tEx[j]  , the start date as tStart[j]  , and the end date
tEnd[j]. Discount factors to the start date and end date are DStart[j] and
Dcap[j]  , and coverage from the start to end date is cvgcap[j]. The  new
quantities we need are tEnd[]  , Dcap[]  , and cvgcap[]. The rest are
calculated in RefDealsCommon */
/* Note: RefDealsCommon must be called first to set nEx  , tEx[]  , tStart[]  ,
 * and DStart[] */

static LGMErr
RefDealsCap(LGMDealType DealType, void *DealPtr, long nArr, Date *TauArr,
            double *FVArr, LGMCalParm *CalReq,
            LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double *),
            LGMErr (*GetBeta)(Date, Date, double *), String YcName,
            LGMMarkConv *conv, LGMCalSet *CSPtr) /*OUTPUT*/
{
  LGMErr error = NULL;
  long i0, j0, j;
  long islong;
  SrtCallInvFlt *BerInvFlt;
  SrtCallCapFlt *BerCapFlt;
  SrtSimMidAt *MidAt;

  /* Allocate space for tEnd[]  , Dcap[]  , cvgcap[] */
  CSPtr->tEnd = (Date *)srt_calloc(CSPtr->nEx + 1, sizeof(Date));
  CSPtr->cvgcap = (double *)srt_calloc(CSPtr->nEx + 1, sizeof(double));
  CSPtr->Dcap = (double *)srt_calloc(CSPtr->nEx + 1, sizeof(double));
  CSPtr->Rfcap = (double *)srt_calloc(CSPtr->nEx + 1, sizeof(double));
  CSPtr->Vcap = (double *)srt_calloc(CSPtr->nEx + 1, sizeof(double));
  CSPtr->VegaCap = (double *)srt_calloc(CSPtr->nEx + 1, sizeof(double));

  if (CSPtr->tEnd == NULL || CSPtr->cvgcap == NULL || CSPtr->Dcap == NULL ||
      CSPtr->Rfcap == NULL || CSPtr->Vcap == NULL)
    return ("Allocation Failed in RefDealsCap");

  CSPtr->capflag = 1;

  switch (DealType) {
  case CallCapFlt:
    BerCapFlt = (SrtCallCapFlt *)DealPtr;
    j0 = BerCapFlt->FirstEx;
    i0 = BerCapFlt->iSet[j0];
    break;

  case CallInvFlt:
    BerInvFlt = (SrtCallInvFlt *)DealPtr;
    break;

  case SimpleMidAt:
    MidAt = (SrtSimMidAt *)DealPtr;
    break;
  }

  /* Update Midat Structure for special methods of calibration */
  if (CalReq->LGMOneTwoFactor == 2 && DealType == SimpleMidAt &&
      ((CalReq->Rmeth == EMK1) || (CalReq->Rmeth == EMK2))) {
    error = update_autocal_midat_struct(MidAt, YcName, GetVol, GetBeta);
    if (error)
      return error;
  }

  /* Fill in CSPtr->tEnd[]  , CSPtr->cvgcap[]  , CSPtr->Dcap[] */
  for (j = 1; j <= CSPtr->nEx; j++) /*nEx IS ALSO THE COUPON NBR */
  {
    CSPtr->tEnd[j] = add_unit(
        CSPtr->tStart[j], 12 / (int)conv->cfreq, SRT_MONTH,
        conv->cbdconv); /* NO_BUSDAY_CONVENTION or BUSDAY_CONVENTION...... */

    CSPtr->cvgcap[j] = coverage(
        CSPtr->tStart[j], CSPtr->tEnd[j],
        conv->cbasis); /* NO_BUSDAY_CONVENTION or BUSDAY_CONVENTION...... */

    CSPtr->Dcap[j] = swp_f_df(CSPtr->tNow, CSPtr->tEnd[j], YcName);
    if (CSPtr->Dcap[j] == SRT_DF_ERROR)
      return ("No Discount Factor");

    switch (DealType) {
    case CallCapFlt:
      CSPtr->Rfcap[j] = BerCapFlt->amax[j + i0];
      break;

    case CallInvFlt:
      islong = 0;
      error = FindFixRateRefInst(
          (CSPtr->DStart[j] - CSPtr->Dcap[j]) /
              (CSPtr->cvgcap[j] * CSPtr->Dcap[j]),
          CSPtr->cvgcap[j] * CSPtr->Dcap[j] *
              (1.0 + find_gear(CSPtr->tEnd[j], BerInvFlt->tCpnPay,
                               BerInvFlt->gear, BerInvFlt->nCpn)),
          CSPtr->DStart[j], CSPtr->Dcap[j], CSPtr->Dpay[CSPtr->n],
          CSPtr->tEx[j], CSPtr->tStart[j], CSPtr->tEnd[j], islong, nArr, FVArr,
          TauArr, CSPtr->tNow, CSPtr->tPay[CSPtr->n], conv, CalReq, GetVol,
          GetBeta, &CSPtr->Rfcap[j], NULL, CSPtr);

      if (error != NULL)
        return (error);
      break;

    case SimpleMidAt:
    default:
      islong = 0;
      error = FindFixRateRefInst(
          (CSPtr->DStart[j] - CSPtr->Dcap[j]) /
              (CSPtr->cvgcap[j] * CSPtr->Dcap[j]),
          CSPtr->cvgcap[j] * CSPtr->Dcap[j], CSPtr->DStart[j], CSPtr->Dcap[j],
          CSPtr->Dpay[CSPtr->n], CSPtr->tEx[j], CSPtr->tStart[j],
          CSPtr->tEnd[j], islong, nArr, FVArr, TauArr, CSPtr->tNow,
          CSPtr->tPay[CSPtr->n], conv, CalReq, GetVol, GetBeta,
          &CSPtr->Rfcap[j], NULL, CSPtr);

      if (error != NULL)
        return (error);
      break;
    }
  }

  return error;
}

/*******************************************/
/* Routine to construct the "1 into k" swaptions. For all these swaptions  ,
the exercise date is ftEx  , which is the latter of the first exercise date
tEx[1] or MinMoEx months from today. The start dates of all these swaptions are
ftStart  , and the fixed leg pay dates of swaption j are tPay[i]  , for i=ffirst
,...  ,j. The swaptions run from j=nfirst to j=nlast. The coverage between
paydates ftStart and ftPay[ffirst] is cvgfix  , and the discount factor to
ftStart is DfixStart. */
/* Note: RefDealsCommon must be called first */

static LGMErr RefDeals1intok(
    long nArr, Date *TauArr,
    double *FVArr,                    /* exercise date and relative PV arrays */
    LGMMarkConv *conv, String ycname, /* market conventions  , discount curve */
    LGMCalParm *CalReq,               /* calibration methodology */
    LGMDealType DealType, void *dealPtr, /* deal type and params */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* function to provide exponent (beta) to interpret vols */
    LGMCalSet *CSPtr) /* Output */
{
  long MinMonEx, swmonth, fullper;
  Date spot, ftStart, ftEx;
  Date tLast;
  long i, j, islong, lag, nfirst;
  double reflevel, level, DLast, refsw, rfix;
  SrtCallInvFltPtr cinf;
  LGMErr error;
  Date tNow, *tPay;                         /* alias to simplify coding */
  long ffirst, n, nPay;                     /* alias to simplify coding */
  double cvgfix, DfixStart, *cvgpay, *Dpay; /* alias to simplify coding */

  /* initialize */
  error = NULL;
  CSPtr->fixflag = 1;     /* flag 1 into k swaptions as existing */
  tNow = CSPtr->tNow;     /* alias to simplify coding */
  n = CSPtr->n;           /* alias to simplify coding */
  nPay = CSPtr->nPay;     /* alias to simplify coding */
  tPay = CSPtr->tPay;     /* alias to simplify coding */
  cvgpay = CSPtr->cvgpay; /* alias to simplify coding */
  Dpay = CSPtr->Dpay;     /* alias to simplify coding */
  lag = conv->lag;        /* alias to simplify coding */

  swmonth = 6;   /* swmonth is number of months per fixed leg period */
  fullper = 155; /* if period has more than fullper days  , count it as full*/
  if (conv->sfreq ==
      SRT_ANNUAL) /* period for fulfilling minimum swap length requirement */
  {
    swmonth = 12;
    fullper = 310;
  }
  MinMonEx = CalReq->MinMonToEx; /* minimum months to exercise */
  if (MinMonEx < 3)
    MinMonEx = 3;
  if (MinMonEx > 12)
    MinMonEx = 12;

  if (DealType == CallInvFlt) {
    cinf = (SrtCallInvFltPtr)dealPtr;
  }

  /* Find the exercise date for the 1 into k swaptions */
  spot = add_unit(tNow, lag, SRT_BDAY, SUCCEEDING); /* spot of today */
  ftStart =
      add_unit(spot, MinMonEx, SRT_MONTH, SUCCEEDING); /* plus MinMoEx months */
  ftEx = add_unit(ftStart, -lag, SRT_BDAY, SUCCEEDING); /* less lag bus days */

  if (ftEx < CSPtr->tEx[1])
    ftEx = CSPtr->tEx[1]; /* exer date for the 1 into k */

  ftStart = add_unit(ftEx, lag, SRT_BDAY, SUCCEEDING); /* start for 1 into k */

  CSPtr->ftEx = ftEx;
  CSPtr->ftStart = ftStart;

  /* find first pay date for the 1 into k swaptions */
  for (i = 1; i <= nPay && tPay[i] <= ftStart; i++)
    ;
  ffirst = i;
  cvgfix = coverage(ftStart, tPay[i], conv->sbasis);
  DfixStart = swp_f_df(tNow, ftStart, ycname);
  if (DfixStart == SRT_DF_ERROR)
    return ("no discount factor");

  CSPtr->ffirst = ffirst;
  CSPtr->cvgfix = cvgfix;
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
  CSPtr->Rfix = (double *)srt_calloc(nPay + 1, sizeof(double));
  CSPtr->Vfix = (double *)srt_calloc(nPay + 1, sizeof(double));
  if (CSPtr->Rfix == NULL || CSPtr->Vfix == NULL)
    return ("allocation failed in RefDeals1intok");

  /* Find fixed rates Rfix of 1 into k swaptions */
  tLast = tPay[n];
  DLast = Dpay[n];

  reflevel = cvgfix * Dpay[ffirst];
  for (i = ffirst + 1; i < nfirst; i++) {
    reflevel += cvgpay[i] * Dpay[i];
  }

  if (DealType == CallInvFlt) {
    level =
        cvgfix * Dpay[ffirst] *
        (1.0 + find_gear(tPay[ffirst], cinf->tCpnPay, cinf->gear, cinf->nCpn));
    for (i = ffirst + 1; i < nfirst; i++) {
      level +=
          cvgpay[i] * Dpay[i] *
          (1.0 + find_gear(tPay[i], cinf->tCpnPay, cinf->gear, cinf->nCpn));
    }
  }

  for (j = nfirst; j <= CSPtr->nlast; j++) /* for each 1 into k swaption */
  {
    reflevel += cvgpay[j] * Dpay[j]; /* compute forward level and swap rate */
    refsw = (DfixStart - Dpay[j]) / reflevel; /* cash swap rate */

    if (DealType == CallInvFlt) {
      level +=
          cvgpay[j] * Dpay[j] *
          (1.0 + find_gear(tPay[j], cinf->tCpnPay, cinf->gear, cinf->nCpn));
      level /= DfixStart;
    } else {
      level = reflevel / DfixStart;
    }

    islong = 0; /* use second set of Rdata if other set are long swaptions */
    if (CalReq->usecaps == 1)
      islong = 1; /* use first set of Rdata if other set are caps */

    error = FindFixRateRefInst(
        refsw, /* swap rate of ref swaption*/
        level, DfixStart, Dpay[j],
        DLast, /* df's for ref start and end  , actual end */
        ftEx, ftStart,
        tPay[j], /* exer  , start  , and end date of ref swaption */
        islong,  /* which set of Rdata for strikes */
        nArr, FVArr, TauArr, tNow,
        tLast, /* PV ratio  , exer  , and end date of orig deal */
        conv, CalReq, GetVol,
        GetBeta, /* Calibration method  , market vol and smile */
        &rfix,   /* strike of the reference swaption */
        NULL, CSPtr);

    if (error != NULL)
      return (error);

    CSPtr->Rfix[j] = rfix;
  }
  return (error);
}

/*******************************************/
/* Routine to find the fixed rates of the reference instruments */
/* Use method in CalReq to pick fixed rate of a swaption with
exercise date tEx  , start date tStart  , last pay date tEnd.
Match to the ratio of the PV of the fixed leg to the PV of the strike at
the same exercise date  , where the ratios are stored in FVarr for dates TauArr
for n=1  ,...  ,nArr. By definition  , the ratio is 1 at date tLast (last pay
date of the original deal) Use GetVol and GetBeta to ensure that the strike is
within two std deviations of par */
static LGMErr FindFixRateRefInst(
    double Rsw, /* swap rate for ref deal */
    double Lvl, double DStart, double DEnd,
    double DLast,                     /* df's to start date and end dates */
    Date tEx, Date tStart, Date tEnd, /* dates of refence swaption */
    long islong,                      /* 1 = use Rdata1  , 0 = use Rdata2 */
    long nArr, double *FVArr,
    Date *TauArr,          /* PV ratio & exer dates of orig deal */
    Date tNow, Date tLast, /* end date of original deal */
    LGMMarkConv *conv,     /* market conventions */
    LGMCalParm *CalReq,    /* calibration request */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *),              /* gets swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double *), /* gets exp (beta) for vols */
    double *rfix, /* returns fix rate of ref.swaption */
    double *std, LGMCalSet *CSPtr) {
  LGMErr error = NULL;
  double stddev, CEVvol, ATMVol, StDLong, beta;
  double gamma, theta1, theta2;
  Date tEx2;
  SrtDiffusionType input_vol;
  long j;
  double numer, MidatStrike;

  /*	Calculate standard deviation */
  error = GetVol(tStart, tEnd, Rsw, SRT_FALSE, &CEVvol);
  if (error)
    return error;
  error = GetBeta(tStart, tEnd, &beta);
  if (error)
    return error;

  if (beta == 0) {
    input_vol = SRT_NORMAL;
  } else {
    input_vol = SRT_LOGNORMAL;
  }

  stddev = fabs(CEVvol) * pow(Rsw, beta);
  stddev *= sqrt(((double)(tEx - tNow)) / 365.);

  if (CalReq->Rmeth != IRR && CalReq->Rmeth != dIRR && CalReq->Rmeth != EMK1 &&
      CalReq->Rmeth != EMK2) {
    if (islong == 1)
      gamma =
          LinearFlatFlat(tEx, CalReq->numR1, CalReq->Rdate1, CalReq->Rdata1);
    else
      gamma =
          LinearFlatFlat(tEx, CalReq->numR2, CalReq->Rdate2, CalReq->Rdata2);
  }

  /*	Calculate midat strike */
  theta1 = InterpPVratio(tEx, nArr, FVArr, TauArr, tLast);
  tEx2 = add_unit(tEnd, -conv->lag, SRT_BDAY, SUCCEEDING);
  theta2 = InterpPVratio(tEx2, nArr, FVArr, TauArr, tLast);
  numer = theta1 - theta2 * DEnd / DStart;
  MidatStrike = numer / Lvl;

  switch (CalReq->Rmeth) {
  case IRR:

    /*
                    theta = InterpPVratio(tEx  , nArr  , FVArr  , TauArr  ,
       tLast); if (DLast > DStart/(1.0 + 0.25*Rsw)) DLast = DStart/(1.0 +
       0.25*Rsw); *rfix = Rsw*(1.0 + (theta-1.0)/(1.0-DLast/DStart)); break;
    */
    return "IRR Not available any more - use strike choices 1 (IRR)  , 2 "
           "(EMK1) or 3 (EMK2)";
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

    if (islong == 1) {
      *rfix = MidatStrike;
      *std = (MidatStrike - Rsw) / stddev;
    } else if (islong == 0) {
      j = BasicDichotomie(CSPtr->tEx, 1, CSPtr->nEx, tEx);

      StDLong = (CSPtr->StDLong[j]);

      error = GetVol(tStart, tEnd, Rsw, SRT_TRUE, &ATMVol);

      *rfix = Rsw + StDLong * stddev;
    }

    break;

  case EMK2:

    if (islong == 1) {
      *rfix = MidatStrike;
    } else if (islong == 0) {
      *rfix = Rsw;
    }

    break;
  }

  /* ensure that the fixed rate is within 2 stddev's of par */
  if (*rfix < Rsw - CalReq->maxstd * stddev) {
    *rfix = Rsw - CalReq->maxstd * stddev;
  } else if (*rfix > Rsw + CalReq->maxstd * stddev) {
    *rfix = Rsw + CalReq->maxstd * stddev;
  }

  return error;
}

/* Helper routine  for picking strikes */
/* Use linear interpolation (flat at ends) to obtain Rdata at thedate */
static double LinearFlatFlat(Date thedate, long n, Date *t, double *Rdata) {
  double slope, ans;
  long i;

  if (n == 1 || thedate >= t[n - 1]) {
    ans = Rdata[n - 1];
    return (ans);
  } else if (thedate <= t[0]) {
    ans = Rdata[0];
    return (ans);
  }

  for (i = 1; t[i] <= thedate; i++)
    ;
  slope = (Rdata[i] - Rdata[i - 1]) / ((double)(t[i] - t[i - 1]));
  ans = Rdata[i] + slope * ((double)(thedate - t[i]));
  return (ans);
}

/* Helper routine  for picking strikes */
/* Use interpolation to find the PV ratio at thedate */
/* Linear at front  , linear at back with the proviso that
        FVArr=1 at tLast
*/
static double InterpPVratio(Date thedate, long n, double *FV, Date *Tau,
                            Date tLast) {
  double slope, ans;
  long i;

  if (n == 1 || thedate >= Tau[n]) {
    slope = (1.0 - FV[n]) / ((double)(tLast - Tau[n]));
    ans = 1.0 + slope * ((double)(thedate - tLast));
    return (ans);
  } else if (thedate <= Tau[1]) {
    slope = (FV[2] - FV[1]) / ((double)(Tau[2] - Tau[1]));
    ans = FV[1] + slope * ((double)(thedate - Tau[1]));
    return (ans);
  } else
    ;
  for (i = 2; Tau[i] <= thedate; i++)
    ;
  slope = (FV[i] - FV[i - 1]) / ((double)(Tau[i] - Tau[i - 1]));
  ans = FV[i] + slope * ((double)(thedate - Tau[i]));
  return (ans);
}

/******************************************************************************/
/*** Routines for finding the market prices of the reference instruments ******/
/******************************************************************************/
/* Main routine to find the market prices of the reference instruments */
LGMErr LGMPriceRefDeals(
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* function to provide exponent (beta) to interpret vols */
    LGMErr (*RecVol)(Date, Date,
                     double), /* function to provide swaption/cap vols */
    LGMCalSet *CSPtr,         /* Output */
    SrtLgmRefSwptnData *lgmRefSwptnData) /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
{
  SrtReceiverType recpay;
  LGMErr error;
  double CEVvol, beta;
  long i, j, iend, nEx;
  Date tNow, *tEx, *tStart, *tPay, *tEnd;
  double *DStart, *Dpay, *cvgpay;
  double SwapLevel, SwapRate, Yrtoexp;
  long CpnIndex;

  double PVTest = 0.;

  error = NULL;
  recpay = SRT_RECEIVER;

  nEx = CSPtr->nEx;
  tNow = CSPtr->tNow;
  tEx = CSPtr->tEx;
  tStart = CSPtr->tStart;
  DStart = CSPtr->DStart;
  tPay = CSPtr->tPay;
  Dpay = CSPtr->Dpay;
  cvgpay = CSPtr->cvgpay;
  tEnd = CSPtr->tEnd;

  /* find prices of long swaptions if needed */
  if (CSPtr->longflag != 0) {
    if (lgmRefSwptnData != NULL) {
      lgmRefSwptnData->NrefLongSwptn = nEx;
      lgmRefSwptnData->refLongSwptnArr = (SrtLgmRefSwptn *)srt_calloc(
          lgmRefSwptnData->NrefLongSwptn, sizeof(SrtLgmRefSwptn));
      if (lgmRefSwptnData->refLongSwptnArr == NULL) {
        lgmRefSwptnData->NrefLongSwptn = 0;
      }
    }

    for (j = 1; j <= nEx; j++) {
      i = CSPtr->ifirst[j];   /* FIRST PAYMENT */
      iend = CSPtr->nlong[j]; /* LAST PAYMENT */

      /* GET VOL  , RECORD VOL & GET CEV EXPONENT */
      error =
          GetVol(tStart[j], tPay[iend], CSPtr->Rflong[j], SRT_FALSE, &CEVvol);
      error = RecVol(tStart[j], tPay[iend], CSPtr->Rflong[j]);
      error = GetBeta(tStart[j], tPay[iend], &beta);
      if (error != NULL)
        return (error);

      /* COMPUTATION OF THE REFS SWAPTION VEGA */

      SwapLevel = 0.;
      SwapLevel = CSPtr->cvgfirst[j] * Dpay[i]; /* i=CSPtr->ifirst[j] */
      for (CpnIndex = (i + 1); CpnIndex <= iend; CpnIndex++)
        SwapLevel = SwapLevel + cvgpay[CpnIndex] * Dpay[CpnIndex];

      SwapRate = (DStart[j] - Dpay[iend]) / SwapLevel;
      Yrtoexp = (double)(tEx[j] - tNow) * YEARS_IN_DAY;

      CSPtr->FrLong[j] = SwapRate;
      CSPtr->CEVLong[j] = CEVvol;
      CSPtr->BetaLong[j] = beta;

      CSPtr->Vlong[j] = LGMCEVSwaptionPrice(
          tNow, tEx[j], recpay, CSPtr->Rflong[j], DStart[j], CSPtr->cvgfirst[j],
          Dpay[i], (iend - i), &(cvgpay[i + 1]), &(Dpay[i + 1]), CEVvol, beta);

      CSPtr->VegaLong[j] = LGMCEVSwaptionPrice(
          tNow, tEx[j], recpay, CSPtr->Rflong[j], DStart[j], CSPtr->cvgfirst[j],
          Dpay[i], (iend - i), &(cvgpay[i + 1]), &(Dpay[i + 1]),
          CEVvol + CEVvolShift, beta);
      CSPtr->VegaLong[j] -= CSPtr->Vlong[j];
      CSPtr->VegaLong[j] /= CEVvolShift;

      CSPtr->StDLong[j] =
          (CSPtr->Rflong[j] - SwapRate) / (CEVvol * SwapRate * sqrt(Yrtoexp));
      CSPtr->StDLong[j] *= (recpay == SRT_RECEIVER ? 1 : -1);

      if (CSPtr->Vlong[j] < 0.)
        return ("can't get long swaption price");

      /* Store swaptions descriptions if needed */
      if (lgmRefSwptnData != NULL) {
        lgmRefSwptnData->refLongSwptnArr[j - 1].bgnDt = tStart[j];
        lgmRefSwptnData->refLongSwptnArr[j - 1].endDt = tPay[iend];
        lgmRefSwptnData->refLongSwptnArr[j - 1].strike = CSPtr->Rflong[j];
        lgmRefSwptnData->refLongSwptnArr[j - 1].bsPv = CSPtr->Vlong[j];
        lgmRefSwptnData->refLongSwptnArr[j - 1].bsVol = CEVvol;
      }
    }
  }

  /* find prices of short swaptions if needed */
  if (CSPtr->shortflag != 0) {
    if (lgmRefSwptnData != NULL) {
      lgmRefSwptnData->NrefShortSwptn = nEx;
      lgmRefSwptnData->refShortSwptnArr = (SrtLgmRefSwptn *)srt_calloc(
          lgmRefSwptnData->NrefShortSwptn, sizeof(SrtLgmRefSwptn));
      if (lgmRefSwptnData->refShortSwptnArr == NULL) {
        lgmRefSwptnData->NrefShortSwptn = 0;
      }
    }

    for (j = 1; j <= nEx; j++) {
      i = CSPtr->ifirst[j];    /* first payment */
      iend = CSPtr->nshort[j]; /* last payment */

      error = GetVol(tStart[j], tPay[iend], CSPtr->Rfshort[j], SRT_FALSE,
                     &CEVvol);                                  /* Get CEVvol */
      error = RecVol(tStart[j], tPay[iend], CSPtr->Rfshort[j]); /* Record vol */
      error = GetBeta(tStart[j], tPay[iend], &beta); /* Get CEV exponent */
      if (error != NULL)
        return (error);

      /* COMPUTE THE LEVEL AND THE SWAP RATE */
      SwapLevel = 0.;
      SwapLevel = CSPtr->cvgfirst[j] * Dpay[i]; /* i=CSPtr->ifirst[j] */
      for (CpnIndex = (i + 1); CpnIndex <= iend; CpnIndex++)
        SwapLevel = SwapLevel + cvgpay[CpnIndex] * Dpay[CpnIndex];

      SwapRate = (DStart[j] - Dpay[iend]) / SwapLevel;
      Yrtoexp = (double)(tEx[j] - tNow) * YEARS_IN_DAY;

      CSPtr->Vshort[j] = LGMCEVSwaptionPrice(
          tNow, tEx[j],                /* eval date  , exer date */
          recpay,                      /* receiver or payer */
          CSPtr->Rfshort[j],           /* strike */
          DStart[j],                   /* dis fac to start */
          CSPtr->cvgfirst[j], Dpay[i], /* stub cvg  , dis fac to stub paydate */
          (iend - i), &(cvgpay[i + 1]), &(Dpay[i + 1]), /* rest of fixed leg */
          CEVvol, beta);

      if (CSPtr->Vshort[j] < 0.)
        return ("can't get short swaption price");

      /* Store swaptions descriptions if needed */
      if (lgmRefSwptnData != NULL) {
        lgmRefSwptnData->refShortSwptnArr[j - 1].bgnDt = tStart[j];
        lgmRefSwptnData->refShortSwptnArr[j - 1].endDt = tPay[iend];
        lgmRefSwptnData->refShortSwptnArr[j - 1].strike = CSPtr->Rfshort[j];
        lgmRefSwptnData->refShortSwptnArr[j - 1].bsPv = CSPtr->Vshort[j];
        lgmRefSwptnData->refShortSwptnArr[j - 1].bsVol = CEVvol;
      }
    }
  }

  /* find prices of caplets if needed */
  if (CSPtr->capflag != 0) {
    if (lgmRefSwptnData != NULL) {
      lgmRefSwptnData->NrefCapSwptn = nEx;
      lgmRefSwptnData->refCapSwptnArr = (SrtLgmRefSwptn *)srt_calloc(
          lgmRefSwptnData->NrefCapSwptn, sizeof(SrtLgmRefSwptn));
      if (lgmRefSwptnData->refCapSwptnArr == NULL) {
        lgmRefSwptnData->NrefCapSwptn = 0;
      }
    }

    for (j = 1; j <= nEx; j++) {
      /* GET CEVvol  , RECORD VOL AND GET CEV EXPONENT */
      error = GetVol(tStart[j], tEnd[j], CSPtr->Rfcap[j], SRT_TRUE, &CEVvol);
      error = RecVol(tStart[j], tEnd[j], CSPtr->Rfcap[j]);
      error = GetBeta(tStart[j], tEnd[j], &beta);
      if (error != NULL)
        return (error);

      /*
      SwapLevel = CSPtr->cvgcap[j]*CSPtr->Dcap[j];
      SwapRate = (DStart[j]-CSPtr->Dcap[j])/SwapLevel;
      Yrtoexp = (double)(tEx[j]-tNow)*YEARS_IN_DAY;

      CSPtr->Vcap[j] = LGMCEVCapletGreek(recpay  ,SwapRate  ,CSPtr->Rfcap[j]
      ,SwapLevel  ,Yrtoexp  ,CEVvol  ,beta  ,PREMIUM); if (CSPtr->Vcap[j] < 0.)
      return("can't get cap price"); CSPtr->VegaCap[j] =
      LGMCEVCapletGreek(recpay  ,SwapRate  ,CSPtr->Rfcap[j]  ,SwapLevel ,Yrtoexp
      ,CEVvol  ,beta  ,VEGA);

      */

      CSPtr->Vcap[j] = LGMCEVCapletPrice(tNow, tEx[j], recpay, CSPtr->Rfcap[j],
                                         CSPtr->cvgcap[j], DStart[j],
                                         CSPtr->Dcap[j], CEVvol, beta);

      CSPtr->VegaCap[j] = LGMCEVCapletPrice(
          tNow, tEx[j], recpay, CSPtr->Rfcap[j], CSPtr->cvgcap[j], DStart[j],
          CSPtr->Dcap[j], (CEVvol + CEVvolShift), beta);

      CSPtr->VegaCap[j] -= CSPtr->Vcap[j];
      CSPtr->VegaCap[j] /= CEVvolShift;

      /* Store swaptions descriptions if needed */
      if (lgmRefSwptnData != NULL) {
        lgmRefSwptnData->refCapSwptnArr[j - 1].bgnDt = tStart[j];
        lgmRefSwptnData->refCapSwptnArr[j - 1].endDt = tEnd[j];
        lgmRefSwptnData->refCapSwptnArr[j - 1].strike = CSPtr->Rfcap[j];
        lgmRefSwptnData->refCapSwptnArr[j - 1].bsPv = CSPtr->Vcap[j];
        lgmRefSwptnData->refCapSwptnArr[j - 1].bsVol = CEVvol;
      }
    }
  }

  /* find prices of 1 into k swaptions if needed */
  if (CSPtr->fixflag != 0) {
    if (lgmRefSwptnData != NULL) {
      lgmRefSwptnData->NrefFixSwptn = CSPtr->nlast - CSPtr->nfirst + 1;
      lgmRefSwptnData->refFixSwptnArr = (SrtLgmRefSwptn *)srt_calloc(
          lgmRefSwptnData->NrefFixSwptn, sizeof(SrtLgmRefSwptn));
      if (lgmRefSwptnData->refFixSwptnArr == NULL) {
        lgmRefSwptnData->NrefFixSwptn = 0;
      }
    }

    i = CSPtr->ffirst;
    for (j = CSPtr->nfirst; j <= CSPtr->nlast; j++) {
      error = GetVol(CSPtr->ftStart, tPay[j], CSPtr->Rfix[j], SRT_FALSE,
                     &CEVvol);                                 /* Get CEVvol */
      error = RecVol(CSPtr->ftStart, tPay[j], CSPtr->Rfix[j]); /* Record vol */
      error = GetBeta(CSPtr->ftStart, tPay[j], &beta); /* Get CEV exponent */
      if (error != NULL)
        return (error);

      /* COMPUTE THE LEVEL AND THE SWAP RATE  */
      SwapLevel = 0.;
      SwapLevel = CSPtr->cvgfix * Dpay[i];
      for (CpnIndex = (i + 1); CpnIndex <= j; CpnIndex++)
        SwapLevel = SwapLevel + cvgpay[CpnIndex] * Dpay[CpnIndex];

      SwapRate = (CSPtr->DfixStart - Dpay[iend]) / SwapLevel;
      Yrtoexp = (double)(CSPtr->ftEx - tNow) * YEARS_IN_DAY;

      CSPtr->Vfix[j] =
          LGMCEVSwaptionPrice(tNow, CSPtr->ftEx, recpay, CSPtr->Rfix[j],
                              CSPtr->DfixStart, CSPtr->cvgfix, Dpay[i], (j - i),
                              &(cvgpay[i + 1]), &(Dpay[i + 1]), CEVvol, beta);

      if (CSPtr->Vfix[j] < 0.)
        return ("can't get 1 into k swaption price");

      /* Store swaptions descriptions if needed */
      if (lgmRefSwptnData != NULL) {
        lgmRefSwptnData->refFixSwptnArr[j - CSPtr->nfirst].bgnDt =
            CSPtr->ftStart;
        lgmRefSwptnData->refFixSwptnArr[j - CSPtr->nfirst].endDt = tPay[j];
        lgmRefSwptnData->refFixSwptnArr[j - CSPtr->nfirst].strike =
            CSPtr->Rfix[j];
        lgmRefSwptnData->refFixSwptnArr[j - CSPtr->nfirst].bsPv =
            CSPtr->Vfix[j];
        lgmRefSwptnData->refFixSwptnArr[j - CSPtr->nfirst].bsVol = CEVvol;
      }
    }
  }

  return (NULL);
}

double LGMCEVSwaptionPrice(
    Date tNow,              /* evaluation date */
    Date tExer,             /* exercise date */
    SrtReceiverType recpay, /* receiver or payer */
    double Rfix,            /*fixed rate of swaption */
    double DStart,  /* dis fac from tNow to start date (PV of floating leg) */
    double stubcvg, /* cvg for initial stub period (if any) */
    double Dstub,   /* dis factor to the stub's pay date */
    long n,         /* number of other fixed leg pay dates */
    double *cvg,    /* [0  ,...  ,n-1] coverages of rest of periods */
    double *Dpay,   /* [0  ,...  ,n-1] dis factor to rest of pay dates */
    double CEVvol,  /* CEV vol */
    double beta)    /* CEV exponent */
{
  double level, Rsw, yrtoexp, price;
  double normvol, disfac;
  SrtCallPutType callput;
  long i;

  if (recpay == SRT_RECEIVER)
    callput = SRT_PUT;
  else
    callput = SRT_CALL;

  /* compute level and swap rate */
  level = stubcvg * Dstub;
  for (i = 0; i < n; i++)
    level = level + cvg[i] * Dpay[i];

  Rsw = (DStart - Dpay[n - 1]) / level;
  yrtoexp = (tExer - tNow) / 365.0;

  /* get equivalent vol for normal (absolute diffusion) model */
  normvol = LGMNormalVolFromCEV(Rsw, Rfix, yrtoexp, CEVvol, beta);

  disfac = 1.0;
  price = level * LGMNormOptPrice(Rsw, Rfix, normvol, yrtoexp, disfac, callput);
  return (price);
}

/******************************************************/
/* Returns the CEV model's price of a caplet/floorlet */
double LGMCEVCapletPrice(
    Date tNow,              /* evaluation date */
    Date tExer,             /* exercise date */
    SrtReceiverType recpay, /* receiver or payer (floorlet/caplet) */
    double Rfix,            /* fixed rate of swaption */
    double cvg,             /* coverage over the period */
    double DStart, /* dis fac from tNow to start date (PV of floating leg) */
    double Dend,   /* dis factor to the end date */
    double CEVvol, /* CEV vol */
    double beta)   /* CEV exponent */
{
  double level, Rsw, yrtoexp, price;
  double normvol, disfac;
  SrtCallPutType callput;

  if (recpay == SRT_RECEIVER)
    callput = SRT_PUT;
  else
    callput = SRT_CALL;

  /* compute level and swap rate */
  level = cvg * Dend;
  Rsw = (DStart - Dend) / level;
  yrtoexp = (tExer - tNow) / 365.0;

  /* get equivalent vol for normal (absolute diffusion) model */
  normvol = LGMNormalVolFromCEV(Rsw, Rfix, yrtoexp, CEVvol, beta);

  disfac = 1.0;
  price = level * LGMNormOptPrice(Rsw, Rfix, normvol, yrtoexp, disfac, callput);
  return (price);
}

static double find_gear(Date thedate, Date *d, double *gear, int n) {
  int i;

  i = 0;
  while (i < n - 1 && d[i] < thedate) {
    i++;
  }

  return gear[i];
}

#undef CEVvolShitf
