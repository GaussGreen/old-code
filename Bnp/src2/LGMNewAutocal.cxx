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
#include "srt_h_lgmUSprotos.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"

/****************************************************************************/
/* Headers */
/****************************************************************************/
#define INTRINSICDEFAULT -333333

/****************************************************************************/
/* Prototypes of private functions */
/****************************************************************************/
/****************************************************************************/

static Date StFromEx(Date tEx, int lag, int calbus, BusDayConv bdr);

static Date ExFromSt(Date tStart, int lag, int calbus, BusDayConv bdr);

/*	Returns intrinsic value and the ratio of the PV of fixed leg to the PV
of the strike (strike is assumed to be the first payment). If deal is strictly
positive or strictly negative        , ratio is set to 0. Note that tNow is
needed only for computing the intrinsic value */
static LGMErr CompValGen(Date tNow, String ycName, double *ratio,
                         double *intValPtr, long nPay, Date *tPay,
                         double *Payment);

/*	Returns intrinsic value and the ratio of the PV of fixed leg to the PV
of the strike. If deal is strictly positive or strictly negative        , ratio
is set to 0. Note that tNow is needed only for computing the intrinsic value */
static LGMErr CompValOpt(Date tNow, String ycName, double *ratio,
                         double *intValPtr, SrtReceiverType PayRec, Date tStart,
                         double Strike, double redfirstPay, long nPay,
                         Date *tPay, double *Payment);

/*	Returns intrinsic value and the ratio PV of cpn leg to  PV of the
strike. If deal is strictly positive or strictly negative        , ratio is set
to 0
*/
static LGMErr CompValInvFlt(
    Date tNow, String ycName, /* yield curve */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double *), /* cap vols */
    LGMErr (*GetBeta)(Date, Date, double *), /* cap exponent (beta) */
    double *ratioPtr,                        /* ratio */
    double *intValPtr,                       /* intrinsic value */
    SrtReceiverType PayRec,                  /* pay/rec */
    Date tStart, double Strike,              /* settlement        , strike */
    long n, Date *t, /* num or periods        , cpn dates */
    double *a, double *marg,
    double lvg,            /* handle        , margin        , leverage */
    SrtBasisCode cpnbasis, /* coupon basis */
    SrtBasisCode basis);   /* basis for rates */

static LGMErr CompValCapFlt(
    Date tNow, String ycName, /* yield curve */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double *), /* cap vols */
    LGMErr (*GetBeta)(Date, Date, double *), /* cap exponent (beta) */
    double *ratioPtr,                        /* ratio */
    double *intValPtr,                       /* intrinsic value */
    SrtReceiverType PayRec,                  /* pay/rec */
    Date tStart, double Strike,              /* settlement        , strike */
    long n, Date *t, /* num or periods        , cpn dates */
    double *amax, double *amin,
    double *marg,          /* handle max and min        , margin*/
    double lvg,            /* leverage */
    SrtBasisCode cpnbasis, /* coupon basis */
    SrtBasisCode basis);   /* basis for rates */

/* Convert xExBdry[i] = x(t) at tExBdry[i] into equivalent swaptions which
would be ATM at the exercise boundary        , for i = 0        ,...
,nExBdry-1 */
static LGMSwptnsPtr FindBdry(int iLGMOneTwoFactor,
                             Date tNow,       /* Evaluation date */
                             String ycName,   /* information about today */
                             LGM_TS *tsPtr,   /* term structure */
                             long nExBdry,    /* number of boundary points*/
                             Date *tExBdry,   /* exercise boundary dates */
                             double *xExBdry, /* exercise boundary (in x) */
                             Date tEnd);      /* last pay date */

/* computes the LGM prices of the deal components (debug) */
static void testcompprice(LGM_TS *tsPtr, Date tNow, String ycname,
                          SrtSimMidAt *dealPtr);

/***************************************************************************************/
/**** Comments ****/
/* Outputs */
/* double LGMValue        , double intrinsicValue be allocated in the calling
routine and their addresses passed to LGMnewautocal as LGMValPtr and
intrinValPtr to receive the value and intrinsic value of the deal

The pointers *ExerIntoPtrPtr        , *ExerBdryPtrPtr        , *LGMtsPtrPtr ,
*SigKaptsPtrPtr must be allocated in the calling routine and their addresses
passed to LGMnewautocal On exit        , *ExerIntoPtrPtr and *ExerBdryPtrPtr
will be the first element of arrays of LGMSwptns (exer date        , final
paydate        , fixed rate) which are most similar to the underlying
midatlantic and define the exercise boundary        , respectively. *LGMtsPtrPtr
will point to the calibrated zeta-G term structure        , and if convertTS is
on        , *SigKaptsPtrPtr will point to the equivalent sigma-kappa term
structure
*/

/**** HOLIDAYS ***/
/* Wherever we need to use the add_unit function        ,
the currency code ccy is available */

/***************************************************************************************/
/* MAIN PROGRAM to calib zeta(t) & G(t) and evaluate an American swaption */
/***************************************************************************************/
/* NOTE: Notional must be 1   */
LGMErr LGMNewAmerican(
    /* task list */
    int skipEval,     /* 0=calibrate & value deal        , 1=calibrate only */
    int skipCalib,    /* 0=calibrate & value deal        , 1=value deal only */
    int convertTS,    /* 1=compute & output sig-kappa ts equivalent to LGM ts */
    int findExerBdry, /* 1=find swap rates at exer boundary  */
                      /* information about today & today's yield curve */
    Date tNow,        /* eval as if today is tNow */
    int eod, /* 0=can exercise on tNow        , 1=cannot exercise on tNow */
    String ycName,       /* market pointer for discount factors */
    double rflt_current, /* current floating rate (if fixed); ignored if
                      non-positive */
                         /* information about today's option prices */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* function to provide exponent (beta) to interpret vols */
    LGMErr (*RecVol)(
        Date, Date,
        double), /* function called whenever LGM uses a vol for pricing */
                 /* the deal */
    LGMAmerType dealtype, /* type of deal */
    void *dealPtr,        /* pointer to the deal */
                          /* calibration and eval parameters */
    LGMCalParm *CalReq,   /* structure defining calibration method to be used */
    ConvParams
        *EvalParms, /* structure containing numerical convolution parameters */
                    /* outputs */
    double *LGMValPtr,    /* value of the deal */
    double *intrinValPtr, /* intrinsic value of the deal */
    LGMSwptns *
        *ExerIntoPtrPtr, /* European swaptions most similar to underlying */
    LGMSwptns **
        ExerBdryPtrPtr, /* European swaptions struck at the exercise boundary */
    SrtLgmRefSwptnData *lgmRefDeals, /* ptr to reference swaption data structure
                                  (NULL => not req'd) */
    LGM_TS **LGMtsPtrPtr,            /* calibrated zeta-G term structure */
    SigKapTS **SigKaptsPtrPtr, /* an equivalent sigma-kappa term structure */
    String szRefRate) /* previously output file name for log file (unused):  now
                         refrate for MUNI case */
{
  /* declarations */
  LGMErr error = NULL;
  SrtCurvePtr yldcrv = NULL;    /* yield curve pointer */
  LGMMarkConv conventions;      /* standard swaption conventions */
  Date tfirst;                  /* first possible exercise date = tNow + eod */
  int k, exfreqmax;             /* maximum months between exercise dates */
  double MidAtValue[7];         /* value of the mid-atlantics */
  int skipCal, convert, findEx; /* task list for LGMnewAutocal */
  LGMDealType ItsAMidAt;
  SrtSimMidAt *MidAtPtr;
  long startover = 0; /* to speed up debugging */

starthere: /* for debugging */
           /* initialize */
  exfreqmax = 2;
  yldcrv = lookup_curve(ycName);

  /* if a non-null refrate is passed in        , set the defaults based on it ,
   * otherwise use the currency defaults */
  if ((szRefRate != 0) && (strcmp(szRefRate, "UTFBM3") == 0))
    error = LGMMuniDefaults(szRefRate, &conventions);
  else
    error = LGMCcyDefaults(yldcrv, &conventions); /* get swap conventions */
  if (error != NULL)
    return (error);

  tfirst = tNow;
  if (eod == 1)
    tfirst = add_unit(tNow, 1, SRT_BDAY, SUCCEEDING);

  /* Step 1: Check validity of deal */
  if (dealtype != SimAmer)
    return ("wrong deal type");
  dealPtr = (SrtSimAmer *)dealPtr;
  LGMValidSimAmer(dealPtr, tfirst);
  if (error != NULL)
    return (error); /* problem with the deal */

  /* Step 2: Create & value mid-atlantic with 2 months between exercise dates ,
  then with one month between exercise dates        , and then extrapolate to
  find the value of the American option */
  ItsAMidAt = SimpleMidAt;
  for (k = exfreqmax; k >= 1; k--) {
    if (k == exfreqmax) {
      skipCal = skipCalib; /* only calibrate first time through */
      convert = convertTS; /* only convert ts first time through */
    } else {
      skipCal = 1; /* only calibrate first time through */
      convert = 0; /* only convert ts first time through */
    }
    if (k != 1)
      findEx = 0; /* only find exer bdry last time through */
    else
      findEx = findExerBdry; /* only find exer bdry last time through */

    MidAtPtr =
        LGMCreateMidAtFromAmer(/* create MidAt with k months between exercises
                                */
                               tfirst, rflt_current, ycName, &conventions, k,
                               dealPtr);
    if (MidAtPtr != NULL) {
      error = LGMnewautocal(
          skipEval, skipCal, convert, findEx, /* task list */
          tNow, eod, ycName,                  /* info about today's market */
          GetVol, GetBeta, RecVol,            /* swaption & caplet price info */
          ItsAMidAt, (void *)MidAtPtr,        /* the deal */
          CalReq, EvalParms,                  /* calib & eval methods */
          &(MidAtValue[k]), intrinValPtr, /* value & intrinsic value of deal */
          ExerIntoPtrPtr, ExerBdryPtrPtr, /* swaptions representing underlying &
                                 exer boundary */
          lgmRefDeals,                    /* reference swaptions */
          LGMtsPtrPtr, SigKaptsPtrPtr,    /* calibrated term structures */
          szRefRate);                     /* unused */

      LGMFreeSimMidAt(&MidAtPtr); /* free midAt structure containing deal */
    } else
      error = "failed to get MidAt from American";

    if (error != NULL) {
      if (skipCal != 1) {
        LGMFreeSwptns(ExerIntoPtrPtr);
        if (CalReq->LGMOneTwoFactor == 1)
          LGMFreeLGM_TS(LGMtsPtrPtr);
        else if (CalReq->LGMOneTwoFactor == 2)
          LGMFreeLGM2F_TS(LGMtsPtrPtr);
      }
      if (convert != 0)
        LGMFreeSigKapTS(SigKaptsPtrPtr);
      if (findEx != 0)
        LGMFreeSwptns(ExerBdryPtrPtr);
      return (error);
    }
  }
  /* extrapolate to get value of the american */
  if (skipEval != 1) {
    if (exfreqmax == 1)
      *LGMValPtr = MidAtValue[1]; /* use monthly exercise for american value */
    else                          /* extrapolate (error goes like h*h */
      *LGMValPtr = (4.0 * MidAtValue[1] - MidAtValue[2]) / 3.0;
  }
  /*	for debugging purposes */
  if (startover > 0)
    goto starthere;

  return (NULL);
}

/*****************************************************************************************/
/* MAIN PROGRAM to calib zeta(t) & G(t) and evaluate a mid_atlantic or Bermudan
 * Inv Flter*/
/*****************************************************************************************/
LGMErr LGMnewautocal(
    /* task list */
    int skipEval,     /* 0=calibrate & value deal        , 1=calibrate only */
    int skipCalib,    /* 0=calibrate & value deal        , 1=value deal only */
    int convertTS,    /* 1=compute & output sig-kappa ts equivalent to LGM ts */
    int findExerBdry, /* 1=find swap rates at exer boundary  */
                      /* information about today & today's market prices */
    Date tNow,        /* eval as if today is tNow */
    int eod, /* 0=can exercise on tNow        , 1=cannot exercise on tNow */
    String ycName, /* market pointer for discount factors */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* function to provide exponent (beta) to interpret vols */
    LGMErr (*RecVol)(
        Date, Date,
        double), /* function called whenever LGM uses a vol for pricing */
                 /* information about the deal */
    LGMDealType dealtype, /* type of deal */
    void *dealPtr,        /* pointer to the deal */
                          /* calibration and eval parameters */
    LGMCalParm *CalReq,   /* structure defining calibration method to be used */
    ConvParams
        *EvalParms, /* structure containing numerical convolution parameters */
                    /* outputs */
    double *LGMValPtr,    /* value of the deal */
    double *intrinValPtr, /* intrinsic value of the deal */
    LGMSwptns *
        *ExerIntoPtrPtr, /* European swaptions most similar to underlying */
    LGMSwptns **
        ExerBdryPtrPtr, /* European swaptions struck at the exercise boundary */
    SrtLgmRefSwptnData *lgmRefDeals, /* ptr to reference swaption data structure
                                  (NULL => not req'd) */
    LGM_TS **LGMtsPtrPtr,            /* calibrated zeta-G term structure */
    SigKapTS **SigKaptsPtrPtr, /* an equivalent sigma-kappa term structure */
    String szRefRate)          /* output file name for log file (unused) */
{
  LGMErr error = NULL;
  LGMMarkConv conventions;         /* standard swaption conventions */
  LGMCalSetPtr RefDealsPtr = NULL; /* set of reference instruments */
  Date tfirst; /* first possible exercise date = tNow + eod */
  long nEx;    /* number of exercise dates */

  long EffnEx;                   /* effective number of exercise dates */
  Date LasttPay, *TauArr = NULL; /* info extracted from deal for calibrator */
  double *FVArr = NULL;          /* info extracted from deal for calibrator */
  double LGMVal = 0.,
         intrinVal = 0.; /* the value and intrinsic value of the deal */
  LGMSwptns *ExerIntoPtr =
      NULL; /* European swaptions most similar to underlying */
  LGMSwptns *ExerBdryPtr =
      NULL; /* European swaptions struck at the exercise boundary */
  LGM_TS *LGMtsPtr = NULL;      /* calibrated zeta-G term structure */
  SigKapTS *SigKaptsPtr = NULL; /* an equivalent sigma-kappa term structure */

  long nExBdry;           /* number of points in the exercise boundary */
  Date *tExBdry = NULL;   /* exer bdry is x(t)=xExBdry[i] at t=tExBdry[i] */
  double *xExBdry = NULL; /* for i=0        ,1        ,...        ,nExBdry-1 */
  SrtCurvePtr yldcrv = NULL; /* yield curve pointer */
  long startover = 0;        /* to speed up debugging */

  double kappa_term, yr_to_exp, prev_yr_to_exp, one_kappa;
  long i;

  yldcrv = lookup_curve(ycName);

  if (findExerBdry != 0 && *ExerBdryPtrPtr != NULL) /* clear for output */
    LGMFreeSwptns(ExerBdryPtrPtr);

  if (skipCalib == 1) {
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

  /* if a non-null refrate is passed in        , set the defaults based on it ,
   * otherwise use the currency defaults */
  if ((szRefRate != 0) && (strcmp(szRefRate, "UTFBM3") == 0))
    error = LGMMuniDefaults(szRefRate, &conventions);
  else
    error = LGMCcyDefaults(yldcrv, &conventions); /* get currency code */
  if (error)
    return (error);

  tfirst = tNow;
  if (eod == 1)
    tfirst = add_unit(tNow, 1, SRT_BDAY, SUCCEEDING /*        , conv.cxxcy*/);

  /* Step 1: Check validity of deal */
  switch (dealtype) {
  case GenEur:
    error = LGMValidGenEur((SrtGenEur *)dealPtr, tfirst, &nEx);
    break;
  case SimpleMidAt:
    error = LGMValidSimMidAt((SrtSimMidAt *)dealPtr, tfirst, &nEx);
    break;
  case GenMidAt:
    error = LGMValidGenMidAt((SrtGenMidAt *)dealPtr, tfirst, &nEx);
    break;
  case CallInvFlt:
    error = LGMValidCallInvFlt((SrtCallInvFlt *)dealPtr, tfirst, &nEx);
    break;
  case CallTimeSwap:
    error = LGMValidCallTimeSwap((SrtCallTimeSwap *)dealPtr, tfirst, &nEx);
    break;
  case CallCapFlt:
    error = LGMValidCallCapFlt((SrtCallCapFlt *)dealPtr, tfirst, &nEx);
    break;
  }
  if (error)
    return (error); /* problem with the deal */

  if (nEx <= 0) {
    *LGMValPtr = 0.0; /* tfirst is strictly after last exer Date */
    *intrinValPtr = 0.0;
    return (NULL);
  }

  /* Step 2: Extract the relevant characteristics of the deal for calibration
  The effective number of exercise dates and the exercise dates - EffnEx ,
  TauArr[] The PV of fixed leg/PV of strike at each exercise date - FVArr[] The
  last pay date of the deal - LasttPay

  Also        , compute the intrinsic value of the deal        , and store the
  underlying exercise dates        , enddates        , and the fixed rate Rfix
  of the vanilla swaption that comes closest to describing the underlying in
  ExerInto structure        , and return these to the calling routine as
  "ExerIntoPtrPtr" */

  error = LGMExtractInfoFromDeal(
      tNow, tfirst, ycName, /* info about today's market */
      GetVol, GetBeta,      /* function to provide swaption/cap vols */
      dealtype, dealPtr,    /* the deal */
      &EffnEx, &LasttPay,
      /* NBR OF EXERCISE DATES        , THE LAST PAYT DATE  */ /* THE MID-AT
                                                            DEAL THEO
                                                            END DATE. */
      &TauArr, &FVArr,           /* info for calibrator */
      &intrinVal, &ExerIntoPtr); /* other info to return as output */

  if (error) {
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
  if (skipCalib != 1) {
    error = LGMCalibration(
        tNow, ycName,            /* info about today's market */
        GetVol, GetBeta, RecVol, /* info about today's swaption/caplet prices */
        dealtype, dealPtr, EffnEx, LasttPay, /* info about the deal */
        TauArr, FVArr,                       /* info about the deal */
        CalReq,       /* the requested calibration method */
        &LGMtsPtr,    /* the calibrated term structure in LGM */
        &RefDealsPtr, /* reference deals used for calibration */
        lgmRefDeals); /* Reference swaptions stock structure */

    LGMFreeCalSet(&RefDealsPtr); /* free unneeded structures */
  } else
    LGMtsPtr = *LGMtsPtrPtr;

  srt_free(FVArr); /* free unneeded structures */

  if (error) {
    LGMFreeSwptns(&ExerIntoPtr);
    if (skipCalib != 1) {
      if (CalReq->LGMOneTwoFactor == 1)
        LGMFreeLGM_TS(LGMtsPtrPtr);
      else if (CalReq->LGMOneTwoFactor == 2)
        LGMFreeLGM2F_TS(LGMtsPtrPtr);
      LGMFreeRefSwptnData(lgmRefDeals);
    }
    return (error);
  }

  /* test calibration against component prices (debug - remove for no debug) */
  /*	testcompprice(LGM1FtsPtr        , tNow        , ycName        ,
   * (SrtSimMidAt*)dealPtr); */

  /* Step 4: Create the equivalent sigma_kappa term structure if requested */
  /* On error        , returns a NULL pointer. We're treating this as non-fatal
   */
  if (skipEval != 1 && convertTS != 0) {
    if (CalReq->LGMOneTwoFactor == 1)
      SigKaptsPtr = LGMConvertZGtoSigKap(tNow, LGMtsPtr);
    else {
      SigKaptsPtr = LGMCreateSigKapTS(nEx, 1);

      SigKaptsPtr->numS = nEx;
      SigKaptsPtr->numK = 1;
      SigKaptsPtr->kdate[0] = LasttPay;
      SigKaptsPtr->kap[0] = LGMtsPtr->one_kappa;

      one_kappa = LGMtsPtr->one_kappa;

      for (i = 1; i <= nEx; i++) {
        SigKaptsPtr->sdate[i - 1] = TauArr[i];
        if (i == 1) {
          yr_to_exp = (double)(TauArr[i] - tNow) * YEARS_IN_DAY;
          kappa_term = (exp(2 * one_kappa * yr_to_exp) - 1) / (2 * one_kappa);
          SigKaptsPtr->sig[i - 1] = sqrt(LGMtsPtr->Zeta1[0] / kappa_term);
        } else {
          yr_to_exp = (double)(TauArr[i] - tNow) * YEARS_IN_DAY;
          prev_yr_to_exp = (double)(TauArr[i - 1] - tNow) * YEARS_IN_DAY;
          kappa_term = (exp(2 * one_kappa * yr_to_exp) -
                        exp(2 * one_kappa * prev_yr_to_exp)) /
                       (2 * one_kappa);
          SigKaptsPtr->sig[i - 1] =
              sqrt(LGMtsPtr->Zeta1[i - 1] - LGMtsPtr->Zeta1[i - 2]) /
              kappa_term;
        }
      }
    }
  }

  srt_free(TauArr);
  /* Step 5: Evaluate deal        , create xExBdry array and fill it with exer
   * bdry */
  if (skipEval != 1 || findExerBdry != 0)
    error = NewMidAtEval(
        EvalParms,                           /* Convolution parameters */
        tNow, eod, ycName, CalReq, LGMtsPtr, /* Information about today */
        GetVol, GetBeta,                     /* swaption/cap vols */
        dealtype, dealPtr,                   /* the deal */
        &LGMVal,                             /* value of the deal */
        &nExBdry, &tExBdry,
        (findExerBdry ? &xExBdry : NULL)); /* array of exercise points (in x) */

  if (error) {
    LGMFreeSwptns(&ExerIntoPtr);
    if (skipCalib != 1) {
      if (CalReq->LGMOneTwoFactor == 1)
        LGMFreeLGM_TS(LGMtsPtrPtr);
      else if (CalReq->LGMOneTwoFactor == 2)
        LGMFreeLGM2F_TS(LGMtsPtrPtr);
      LGMFreeRefSwptnData(lgmRefDeals);
    }
    if (convertTS != 0) {
      if (CalReq->LGMOneTwoFactor == 1) {
        LGMFreeSigKapTS(&SigKaptsPtr);
        srt_free(xExBdry);
        srt_free(tExBdry);
      }
    }
    return (error);
  }

  /* Step 6: Find exer boundary & store as set of swaptions in ExerBdryPtr */
  /* On error        , returns a NULL pointer. We're treating this as non-fatal
   */
  if (skipEval != 1 && findExerBdry != 0 && nExBdry > 0 && tExBdry != NULL &&
      xExBdry != NULL) {
    ExerBdryPtr = FindBdry(CalReq->LGMOneTwoFactor, tNow, ycName, LGMtsPtr,
                           nExBdry, tExBdry, xExBdry, LasttPay);
  }

  if (skipEval != 1 || findExerBdry != 0) {
    srt_free(xExBdry);
    srt_free(tExBdry); /* free unneeded structure */
  }
  /* Set outputs and return */
  if (skipEval != 1) {
    *LGMValPtr = LGMVal;
    *intrinValPtr = intrinVal;
  }
  if (skipCalib != 1) {
    *ExerIntoPtrPtr = ExerIntoPtr;
    *LGMtsPtrPtr = LGMtsPtr;
  } else
    LGMFreeSwptns(&ExerIntoPtr);
  if (skipEval != 1 && findExerBdry != 0)
    *ExerBdryPtrPtr = ExerBdryPtr;
  if (skipEval != 1 && convertTS != 0)
    *SigKaptsPtrPtr = SigKaptsPtr;

  return (error);
}

/********************************************************************************/
SrtSimMidAtPtr LGMCreateMidAtFromAmer(
    Date tfirst,         /* first day exercise is physically possible */
    double rfltCur,      /* floating rate for current period */
    String ycName,       /* yield curve (by name) */
    LGMMarkConv *conv,   /* market conventions */
    int dmos,            /* months between exercise dates */
    SrtSimAmer *AmerPtr) /* the American deal */
{
  SrtSimMidAt *MidAtPtr;
  SrtCurvePtr yldcrv;
  double rflt, cvg, dfa, dfb, dfSt, reduction;
  Date tfirstEx, tfirstStart, tlastEx, tlastStart;
  Date tNow, tLastDate, ta, tb, t1;
  long i, j, k, ifix, nfix, nflt, n, nlast, nfirst;
  long adjusted;
  int lag, calbus;
  BusDayConv bdconv;

  /* MidAt components */
  long nPay, nEx;
  Date *tPay, *tEx, *tStart;
  double *Payment, *Strike, *RedFirstPay;
  long *FirstPay;

  yldcrv = lookup_curve(ycName);
  tNow = get_clcndate_from_yldcrv(yldcrv); /* real eval date irrelevent */
  nfix = AmerPtr->nfix;
  nflt = AmerPtr->nflt;
  lag = AmerPtr->lagExerSettle;
  calbus = AmerPtr->CalBusLag;
  bdconv = AmerPtr->convSettle;

  /* find first exercise and start dates */
  tfirstStart = AmerPtr->tfixStart[0];
  tfirstEx = ExFromSt(tfirstStart, lag, calbus,
                      bdconv); /* no reason to exercise earlier */
  tfirstEx =
      max(tfirstEx, AmerPtr->tFirstExer); /* cannot exercise before allowed */
  tfirstEx = max(tfirstEx, tfirst);       /* cannot exercise before today+eod */
  tfirstStart = StFromEx(tfirstEx, lag, calbus, bdconv);
  if (tfirstStart < tfirst)
    tfirstStart = tfirst; /* should only occur if deal is mis-entered */

  /* find last exercise and start dates */
  if (tfirstStart >= AmerPtr->tfixEnd[nfix - 1]) /* can't make sense of deal */
    return (NULL);
  else if (tfirstStart >= AmerPtr->tfixEnd[nfix - 1] - 30 || nflt == 1) {
    tlastEx = tfirstEx;
    tlastStart = tfirstStart;
  } else {
    tlastEx = AmerPtr->tfltFixing[nflt - 1];
    tlastStart =
        StFromEx(tlastEx, lag, calbus,
                 bdconv); /* no value to exercising after last fixing */
    tlastStart = min(tlastStart, add_unit(AmerPtr->tfixEnd[nfix - 1], -1,
                                          SRT_MONTH, SUCCEEDING));
    tlastEx = ExFromSt(tlastStart, lag, calbus, bdconv);
    if (tlastEx <= tfirstEx ||
        tlastStart <= tfirstStart) /* effectively a vanilla swaption */
    {
      tlastEx = tfirstEx;
      tlastStart = tfirstStart;
    }
  }
  /* set fixed leg */
  for (ifix = 0; AmerPtr->tfixEnd[ifix] <= tfirstStart; ifix++)
    ;
  /* fixed leg runs from i=ifix        , ...        , nfix-1 */
  nPay = nfix - ifix; /* number of fixed periods */

  /* find start (upon exercise) dates */
  tLastDate = AmerPtr->tfixPay[nfix - 1];
  n = 1;
  while (n <= 601 &&
         tlastStart <= add_unit(tLastDate, -n * dmos, SRT_MONTH, conv->sbdconv))
    n++;
  nlast = n - 1;
  while (n <= 601 &&
         tfirstStart < add_unit(tLastDate, -n * dmos, SRT_MONTH, conv->sbdconv))
    n++;
  nfirst = n;
  nEx = nfirst - nlast + 1; /* this number may be reduced later */

  /* start dates are tLastDate - n*dmos        , n=nfirst        ,...        ,
   * nlast (plus adjustment for start & end) */
  MidAtPtr = LGMCreateSimMidAt(nEx, nPay);
  if (MidAtPtr == NULL)
    return (NULL);
  Payment = MidAtPtr->Payment;
  tPay = MidAtPtr->tPay;
  tEx = MidAtPtr->tEx;
  tStart = MidAtPtr->tStart;
  Strike = MidAtPtr->Strike;
  FirstPay = MidAtPtr->FirstPay;
  RedFirstPay = MidAtPtr->RedFirstPay;

  /* set fixed leg */
  for (i = ifix; i < nfix; i++) {
    tPay[i - ifix] = AmerPtr->tfixPay[i];
    Payment[i - ifix] = AmerPtr->fixCoupon[i];
  }
  Payment[nPay - 1] = Payment[nPay - 1] + 1.0;
  MidAtPtr->nPay = nPay;

  /* find start (upon exercise) dates */
  for (n = nfirst; n >= nlast; n--)
    tStart[nfirst - n] =
        add_unit(tLastDate, -n * dmos, SRT_MONTH, conv->sbdconv);

  /* adjust to match schedule (if possible) */
  j = 0;
  k = 0;
  for (i = 0; i < nEx; i++) {
    adjusted = 0;
    for (; j < nfix && AmerPtr->tfixPay[j] < tStart[i] - 10; j++)
      ;
    if (labs(tStart[i] - AmerPtr->tfixPay[j]) < 10 && adjusted == 0) {
      tStart[i] = AmerPtr->tfixPay[j];
      adjusted = 1;
    }

    for (; k < nflt && AmerPtr->tfltPay[k] < tStart[i] - 10; k++)
      ;
    if (labs(tStart[i] - AmerPtr->tfltPay[k]) < 10 && adjusted == 0) {
      tStart[i] = AmerPtr->tfltPay[k];
      adjusted = 1;
    }
  }

  /* adjust to match firstStart and lastStart exactly */
  tStart[0] = tfirstStart; /* eliminate dates too close together */
  if (nEx > 1) {
    tStart[nEx - 1] = tlastStart;
    if (tStart[1] <= tStart[0] + 10) {
      for (i = 2; i < nEx; i++)
        tStart[i - 1] = tStart[i];
      nEx--;
    }
    if (nEx > 2 && tStart[nEx - 2] >= tStart[nEx - 1] - 10) {
      tStart[nEx - 2] = tStart[nEx - 1];
      nEx--;
    }
    if (tStart[nEx - 1] <= tStart[0] + 10)
      nEx = 1;
  }

  MidAtPtr->nEx = nEx;
  MidAtPtr->PayRec = AmerPtr->PayRec;
  MidAtPtr->FirstExer = 0; /* dummy entry ... not used by code */

  /* exercise dates */
  for (i = 0; i < nEx; i++)
    tEx[i] = ExFromSt(tStart[i], lag, calbus, bdconv);
  tEx[0] = tfirstEx;
  if (nEx > 1)
    tEx[nEx - 1] = tlastEx;

  /* first payment        , strike */
  j = 0;
  k = 0;
  for (i = 0; i < nEx; i++) {
    for (; AmerPtr->tfixEnd[j] <= tStart[i]; j++)
      ;
    FirstPay[i] = j - ifix; /* settlement is in fixed period j */
    Strike[i] = 1.0 + AmerPtr->ExtraPrem[j];
    RedFirstPay[i] = 0.;
    reduction =
        AmerPtr->fixCoupon[j] /* compute accrued to settlement date */
        * coverage(AmerPtr->tfixStart[j], tStart[i], AmerPtr->fixBasis) /
        coverage(AmerPtr->tfixStart[j], AmerPtr->tfixEnd[j], AmerPtr->fixBasis);

    if (AmerPtr->EarlyFlagFix == 0)
      RedFirstPay[i] = reduction;
    else
      Strike[i] = Strike[i] + reduction;

    for (; AmerPtr->tfltPay[k] <= tStart[i]; k++)
      ; /* settlement is in floating period k */
    if (k > 0)
      ta = AmerPtr->tfltPay[k - 1];
    else
      ta = AmerPtr->tfixStart[0];
    tb = AmerPtr->tfltPay[k];

    dfb = swp_f_df(tNow, tb, ycName);
    dfSt = swp_f_df(tNow, tStart[i], ycName);

    if (dfb == SRT_DF_ERROR || dfSt == SRT_DF_ERROR) {
      LGMFreeSimMidAt(&MidAtPtr);
      return (NULL);
    }
    if (ta <= tfirst && rfltCur > 0)
      rflt = rfltCur;
    else {
      t1 = max(tfirst, ta); /* estimate the floating rate from today's curve */
      dfa = swp_f_df(tNow, t1, ycName);
      if (dfa == SRT_DF_ERROR) {
        LGMFreeSimMidAt(&MidAtPtr);
        return (NULL);
      }
      cvg = coverage(t1, tb, AmerPtr->fltBasis);
      if (cvg > 0)
        rflt = (dfa - dfb) / (dfb * cvg);
      else
        rflt = 0;
    }
    if (AmerPtr->EarlyFlagFlt == 0 && AmerPtr->ResetFlt == 0 && rflt > 0)
      Strike[i] = Strike[i] - 1.0 +
                  (1.0 + rflt * coverage(tStart[i], tb, AmerPtr->fltBasis)) *
                      dfb / dfSt;
    else if (AmerPtr->EarlyFlagFlt != 0 && rflt > 0)
      Strike[i] =
          Strike[i] - 1.0 - rflt * coverage(ta, tStart[i], AmerPtr->fltBasis) +
          (1.0 + rflt * coverage(ta, tb, AmerPtr->fltBasis)) * dfb / dfSt;
  }
  return (MidAtPtr);
}

/*********************************************************************/
/* start to exercise & exercise to start */
/*********************************************************************/

static Date StFromEx(Date tEx, int lag, int calbus, BusDayConv bdr) {
  Date tStart;
  if (calbus == 1)
    tStart = add_unit(tEx, lag, SRT_BDAY, bdr);
  else
    tStart = add_unit(tEx, lag, SRT_DAY, bdr);
  return (tStart);
}

static Date ExFromSt(Date tStart, int lag, int calbus, BusDayConv bdr) {
  Date tEx;
  if (calbus == 1)
    tEx = add_unit(tStart, -lag, SRT_BDAY, NO_BUSDAY_CONVENTION);
  else {
    tEx = add_unit(tStart, -lag, SRT_DAY, NO_BUSDAY_CONVENTION);
    while (add_unit(tEx, lag, SRT_DAY, bdr) > tStart)
      tEx--;
  }
  return (tEx);
}

/********************************************************************************/
/* Get info for calibration:
The effective number of exercise dates and the exercise dates - nEx        ,
TauArr[] The PV of fixed leg/PV of strike at each exercise date - FVArr[] The
last pay date of the deal - tlast

Also        , compute the intrinsic value of the deal        , and store the
underlying exercise dates        , enddates        , and the fixed rate Rfix of
the vanilla swaption that comes closest to describing the underlying in ExerInto
structure        , and return these to the calling routine as "ExerInto"
*/
LGMErr LGMExtractInfoFromDeal(
    /* information about today's market */
    Date tNow,     /* evaluation date */
    Date tfirst,   /* first possible exercise date (eval date + eod) */
    String ycName, /* yield curve for discount factors */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* function to provide exponent (beta) to interpret vols */
                   /* information about the deal */
    LGMDealType dealtype,    /* type of deal */
    void *dealPtr,           /* pointer to the deal */
                             /* output */
    long *EffnEx,            /* Effective number of exercise dates */
    Date *LasttPay,          /* Last pay date of deal */
    Date **TauArrPtr,        /* Array of effective exercise dates */
    double **FVArrPtr,       /* PV of fixed leg/PV of strike for these dates */
    double *intValPtr,       /* intrinsic value of the deal */
    LGMSwptns **ExerIntoPtr) /* closest European swaptions underlying deal */

{
  SrtGenEurPtr EurDeal;
  SrtSimMidAtPtr MidAtDeal;
  SrtGenMidAtPtr GenDeal;
  SrtCallInvFltPtr BerInvFltDeal;
  SrtCallCapFltPtr BerCapFltDeal;
  SrtCallTimeSwapPtr CallTimeSwapDeal;

  /* temporary variables */
  long i, j, n, nEx, nPay, check;
  Date nextdate, *tEx, *tStart, *tPay, *tCpn;
  double *Strike, *Payment, *redpay;
  SrtReceiverType PayRec;
  double ratio, intVal, tempIntVal, Rfix;
  LGMErr error = NULL;

  /* output variables */
  Date tlast, *TauArr; /* Output variables */
  double *FVArr; /* tlast        , nEx        , TauArr and FVArr go on to the
                    calibrator */
  LGMSwptns
      *ExerInto; /* Underlying swaptions        , returned to calling code */

  if (dealtype >= LGMLastDealType)
    return ("unknown deal type");

  /* First allocate enough space */
  switch (dealtype) {
  case GenEur:
    EurDeal = (SrtGenEur *)dealPtr;
    n = 1;
    if (EurDeal->tEx < tfirst)
      n = 0;
    break;

  case SimpleMidAt:
    MidAtDeal = (SrtSimMidAt *)dealPtr;
    n = MidAtDeal->nEx;
    if (MidAtDeal->tEx[n - 1] < tfirst)
      n = 0;
    break;

  case GenMidAt:
    GenDeal = (SrtGenMidAt *)dealPtr;
    n = GenDeal->nEx;
    if (GenDeal->tEx[n - 1] < tfirst)
      n = 0;
    break;

  case CallTimeSwap:
    CallTimeSwapDeal = (SrtCallTimeSwap *)dealPtr;
    n = CallTimeSwapDeal->nEx;
    if (CallTimeSwapDeal->tEx[n - 1] < tfirst)
      n = 0;
    break;

  case CallInvFlt:
    BerInvFltDeal = (SrtCallInvFlt *)dealPtr;
    n = BerInvFltDeal->nEx;
    if (BerInvFltDeal->tEx[n - 1] < tfirst)
      n = 0;
    break;

  case CallCapFlt:
    BerCapFltDeal = (SrtCallCapFlt *)dealPtr;
    n = BerCapFltDeal->nEx;
    if (BerCapFltDeal->tEx[n - 1] < tfirst)
      n = 0;
    break;

  default:
    return ("unknown deal type");
  }
  if (n <= 0)
    return ("no exercise dates left");

  /* allocate space */
  TauArr = (Date *)srt_calloc(n + 1, sizeof(Date));
  FVArr = (double *)srt_calloc(n + 1, sizeof(double));
  ExerInto = LGMCreateSwptns(n);

  if (TauArr == NULL || FVArr == NULL || ExerInto == NULL) {
    srt_free(TauArr);
    srt_free(FVArr);
    LGMFreeSwptns(&ExerInto);
    return ("alloc failed in CalibFromDeal");
  }

  switch (dealtype) {
  case GenEur: /* Compute intrinsic value and PV ratio */
    nPay = EurDeal->nPay;
    error = CompValGen(tNow, ycName, &ratio, &intVal, nPay, EurDeal->tPay,
                       EurDeal->Payment);
    if (error != NULL) {
      srt_free(TauArr);
      srt_free(FVArr);
      LGMFreeSwptns(&ExerInto);
      return (error);
    }

    nEx = 1; /* Set outputs */
    tlast = EurDeal->tPay[nPay - 1];
    TauArr[1] = EurDeal->tEx;
    FVArr[1] = ratio;
    ExerInto->n = 1;
    ExerInto->tEx[0] = TauArr[1];
    ExerInto->tEnd[0] = tlast;
    break;

  case SimpleMidAt:
    nPay = MidAtDeal->nPay;
    tEx = MidAtDeal->tEx; /* Rename to clarify code! */
    tStart = MidAtDeal->tStart;
    tPay = MidAtDeal->tPay;
    Strike = MidAtDeal->Strike;
    Payment = MidAtDeal->Payment;
    redpay = MidAtDeal->RedFirstPay;
    PayRec = MidAtDeal->PayRec;
    tlast = tPay[nPay - 1];

    intVal = INTRINSICDEFAULT;
    nEx = 0;
    nextdate = tfirst;
    for (j = 0; j < n; j++) /* find all exer dates on or after tfirst */
    {
      if (tEx[j] >= nextdate) /* which are at least 20 days apart */
      {
        nEx++;
        nextdate = tEx[j] + 20;

        TauArr[nEx] = tEx[j]; /* Record these dates */
        ExerInto->tEx[nEx - 1] = TauArr[nEx];
        ExerInto->tEnd[nEx - 1] = tlast;

        i = MidAtDeal->FirstPay[j]; /* Compute PV ratio and intrin val */
        error = CompValOpt(tNow, ycName, &ratio, &tempIntVal, PayRec, tStart[j],
                           Strike[j], redpay[j], nPay - i, &(tPay[i]),
                           &(Payment[i]));

        if (error != NULL) {
          srt_free(TauArr);
          srt_free(FVArr);
          LGMFreeSwptns(&ExerInto);
          return (error);
        }

        FVArr[nEx] = ratio; /* Record the PV ratio for this exercise date */
        if (tempIntVal >
            intVal) /* If this intrinsic value is larger        , replace */
          intVal = tempIntVal; /* the intrinsic value of the deal */
      }
    }
    ExerInto->n = nEx;
    break;

  case GenMidAt:
    tEx = GenDeal->tEx; /* Rename to clarify code! */
    intVal = INTRINSICDEFAULT;
    tlast = 0;

    nEx = 0;
    nextdate = tfirst;
    for (j = 0; j < n; j++) /* find all the exercise dates on or after tfirst */
    {
      if (tEx[j] >= nextdate) /* which are at least 20 days apart */
      {
        nEx++;
        nextdate = tEx[j] + 20;

        nPay = GenDeal->nPay[j];       /* Rename to clarify code! */
        tPay = GenDeal->tPay[j];       /* Rename to clarify code! */
        Payment = GenDeal->Payment[j]; /* Rename to clarify code! */

        TauArr[nEx] = tEx[j]; /* Record these dates */
        ExerInto->tEx[nEx - 1] = tEx[j];
        ExerInto->tEnd[nEx - 1] = tPay[nPay - 1];
        if (tPay[nPay - 1] > tlast)
          tlast = tPay[nPay - 1];

        error = CompValGen(tNow, ycName, &ratio,
                           &tempIntVal, /* Compute PV ratio and intrinsic val */
                           nPay, tPay, Payment); /* of this exer date */

        if (error != NULL) {
          srt_free(TauArr);
          srt_free(FVArr);
          LGMFreeSwptns(&ExerInto);
          return (error);
        }

        FVArr[nEx] = ratio; /* Record the PV ratio for this exercise date */
        if (tempIntVal >
            intVal) /* If this intrinsic value is larger        , replace */
          intVal = tempIntVal; /* the intrinsic value of the deal */
      }
    }
    ExerInto->n = nEx;
    break;

  case CallInvFlt:

    tlast = BerInvFltDeal->tCpnPay[BerInvFltDeal->nCpn - 1];

    intVal = INTRINSICDEFAULT;
    nEx = 0;
    nextdate = tfirst;

    for (j = 0; j < n; j++) /* find all exer dates on or after tfirst */
    {                       /* which are at least 20 days apart */
      if (BerInvFltDeal->tEx[j] >= nextdate) {
        nEx++;
        nextdate = BerInvFltDeal->tEx[j] + 20; /* Record these dates */
        TauArr[nEx] = BerInvFltDeal->tEx[j];

        i = BerInvFltDeal->iSet[j]; /* Compute PV ratio and intrin val */

        error =
            CompValOpt(tNow, ycName, &ratio, &tempIntVal, BerInvFltDeal->PayRec,
                       BerInvFltDeal->tCpnStart[i], BerInvFltDeal->strike[j],
                       0.0, BerInvFltDeal->nCpn - i,
                       &(BerInvFltDeal->tCpnPay[i]), &(BerInvFltDeal->a[i]));

        if (error != NULL) {
          srt_free(TauArr);
          srt_free(FVArr);
          LGMFreeSwptns(&ExerInto);
          return error;
        }

        FVArr[nEx] = ratio; /* Record the PV ratio for this exercise date */
        if (tempIntVal >
            intVal) /* If this intrinsic value is larger        , replace */
          intVal = tempIntVal; /* the intrinsic value of the deal */
      }
    }
    break;

  case CallTimeSwap:

    tlast = CallTimeSwapDeal->tCpnPay[CallTimeSwapDeal->nCpn - 1];
    intVal = INTRINSICDEFAULT;
    nEx = n;
    nextdate = tfirst;
    for (j = 0; j < nEx; j++) {
      TauArr[j + 1] = CallTimeSwapDeal->tEx[j];
    }

    error = CompVal_ForwardsTimeSwaps(CallTimeSwapDeal, tNow, ycName, FVArr,
                                      &intVal, CallTimeSwapDeal->tForwardTS[0],
                                      &(CallTimeSwapDeal->pv_fixedleg), nEx, 0,
                                      0, 0);

    if (error != NULL) {
      srt_free(TauArr);
      srt_free(FVArr);
      LGMFreeSwptns(&ExerInto);
      return error;
    }

    /*
                    tlast = CallTimeSwapDeal->tCpnPay[CallTimeSwapDeal->nCpn-1];

                    intVal = INTRINSICDEFAULT;
                    nEx = 0;
                    nextdate = tfirst;

                    for (j=0; j<n; j++)
                    {
                            if (CallTimeSwapDeal->tEx[j] >= nextdate)
                            {
                                    nEx++;
                                    nextdate = CallTimeSwapDeal->tEx[j] + 20;
                                    TauArr[nEx] = CallTimeSwapDeal->tEx[j];
                                    i = CallTimeSwapDeal->iSet[j];
                                    error =  CompValTimeSwap(
                                                                                            tNow        ,
                                                                                            ycName        ,
                                                                                            &ratio        ,
                                                                                            j        ,
                                                                                            &tempIntVal        ,
                                                                                            CallTimeSwapDeal->strike[j]        ,
                                                                                            CallTimeSwapDeal        ,
                                                                                            1        ,
                                                                                            0        ,
                                                                                            0        ,
                                                                                            0);


                                    CallTimeSwapDeal->tForwardTS[0][j] =
       tempIntVal;


                                    tempIntVal = max(0.        ,tempIntVal);

                                    if (error!=NULL)
                                    {
                                            srt_free (TauArr);
                                            srt_free (FVArr);
                                            LGMFreeSwptns (&ExerInto);
                                            return error;
                                    }
                                    if (1==nEx)
                                    {
                                            CallTimeSwapDeal->pv_fixedleg =
       tempIntVal;
                                    }

                                    FVArr[nEx] = ratio;
                                    if (tempIntVal>intVal)
                                            intVal = tempIntVal;
                            }
                    }

      */

    break;

  case CallCapFlt:
    nPay = BerCapFltDeal->nCpn;
    tEx = BerCapFltDeal->tEx; /* Rename to clarify code! */
    tCpn = BerCapFltDeal->tCpn;
    Strike = BerCapFltDeal->strike;
    PayRec = BerCapFltDeal->PayRec;

    tlast = tCpn[nPay];

    intVal = INTRINSICDEFAULT;
    nEx = 0;
    nextdate = tfirst;
    for (j = 0; j < n; j++) /* find all exer dates on or after tfirst */
    {
      if (tEx[j] >= nextdate) /* which are at least 20 days apart */
      {
        nEx++;
        nextdate = tEx[j] + 20;

        /* Record these dates */
        TauArr[nEx] = tEx[j];
        ExerInto->tEx[nEx - 1] = TauArr[nEx];
        ExerInto->tEnd[nEx - 1] = tlast;

        /* Compute PV ratio and intrin val */
        i = BerCapFltDeal->iSet[j]; /* settlement is on tCpn[i] */
        error = CompValCapFlt(
            tNow, ycName, GetVol, GetBeta, /* today's market info */
            &ratio, &tempIntVal,           /* return */
            PayRec, tCpn[i],
            Strike[j], /* pay/rec        , settlement date        , strike */
            nPay - i,
            &(tCpn[i]), /* number of periods left        , cpn dates */
            &(BerCapFltDeal->amax[i]), /* handles for cap max */
            &(BerCapFltDeal->amin[i]), /* handles for cap min */
            &(BerCapFltDeal->marg[i]), /* margins */
            BerCapFltDeal->lvg,        /* leverage */
            BerCapFltDeal->cpnBasis,   /* cpn payment basis */
            BerCapFltDeal->aBasis);    /* basis for rates */

        if (error != NULL) {
          srt_free(TauArr);
          srt_free(FVArr);
          LGMFreeSwptns(&ExerInto);
          return (error);
        }

        FVArr[nEx] = ratio; /* Record the PV ratio for this exercise date */
        if (tempIntVal >
            intVal) /* If this intrinsic value is larger        , replace */
          intVal = tempIntVal; /* the intrinsic value of the deal */
      }
    }
    ExerInto->n = nEx;
    break;
  }

  /* Check to see if output makes sense */
  check = 0;
  for (j = 1; j <= nEx; j++) {
    if (FVArr[j] < 0.2 || FVArr[j] > 5.0)
      check = 1;
  }
  if (check != 0) {
    for (j = 1; j <= nEx; j++)
      FVArr[j] = 1.0;
  }

  /* The swaption array ExerInto has n        , the exercise dates tEx[] , and
  the end dates tEnd[] set. We need to determine the fixed rates Rfix[]
  from the forward value array */

  /* Not used for CIF */
  /* Should actually be completely removed later */

  if (dealtype != CallInvFlt && dealtype != CallTimeSwap) {
    for (j = 0; j < nEx; j++) {
      error = LGMGetFixRate(ExerInto->tEx[j], ExerInto->tEnd[j], FVArr[j + 1],
                            ycName, &Rfix);
      if (error != NULL) {
        srt_free(TauArr);
        srt_free(FVArr);
        LGMFreeSwptns(&ExerInto);
        return (error);
      }
      ExerInto->Rfix[j] = Rfix;
    }
  }

  /* End of what should be removed */

  /* Output for calibration */
  *EffnEx = nEx;
  *LasttPay = tlast;
  *TauArrPtr = TauArr;
  *FVArrPtr = FVArr;
  *intValPtr = intVal;
  *ExerIntoPtr = ExerInto;
  return (error);
}

/********************************************************************************/
/*	Returns intrinsic value and the ratio of the PV of fixed leg to the PV
of the strike (strike is assumed to be the first payment). If deal is strictly
positive or strictly negative        , ratio is set to 0. Note that tNow is
needed only for computing the intrinsic value */
static LGMErr CompValGen(Date tNow, String ycName, double *ratioPtr,
                         double *intValPtr, long nPay, Date *tPay,
                         double *Payment) {
  long j;
  double Strike, df, intVal, ratio;

  intVal = INTRINSICDEFAULT;
  ratio = 0.;

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
  intVal = intVal + Strike;     /* intrinsic val */

  if (ratio > 5.0 || ratio < 0.2) /* clearly        , one sign deal */
    ratio = 0.;

  *ratioPtr = ratio;
  *intValPtr = intVal;
  return (NULL);
}

/********************************************************************************/
/*	Returns intrinsic value and the ratio of the PV of fixed leg to the PV
of the strike. If deal is strictly positive or strictly negative        , ratio
is set to 0. Note that tNow is needed only for computing the intrinsic value */
static LGMErr CompValOpt(Date tNow, String ycName, double *ratioPtr,
                         double *intValPtr, SrtReceiverType PayRec, Date tStart,
                         double Strike, double redfirstPay, long nPay,
                         Date *tPay, double *Payment) {
  long j;
  double sum, df, ratio, intVal;
  double signleg;

  signleg = 1.; /* intrin val = signleg*(fix leg - strike) */
  if (PayRec == SRT_PAYER)
    signleg = -1.;

  intVal = INTRINSICDEFAULT;
  ratio = 0.;

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

  if (ratio > 5.0 || ratio < 0.2) /* clearly        , one sign deal */
    ratio = 0.;

  *ratioPtr = ratio;
  *intValPtr = intVal;
  return (NULL);
}

/********************************************************************************/
/*	Returns intrinsic value and the ratio PV of coupon leg to PV of the
strike for an inverse floater. If deal is strictly positive or strictly negative
      , ratio is set to 0. */
static LGMErr CompValInvFlt(
    Date tNow, String ycName, /* yield curve */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double *), /* cap vols */
    LGMErr (*GetBeta)(Date, Date, double *), /* cap exponent (beta) */
    double *ratioPtr,                        /* ratio */
    double *intValPtr,                       /* intrinsic value */
    SrtReceiverType PayRec,                  /* pay/rec */
    Date tStart, double Strike,              /* settlement        , strike */
    long n, Date *t, /* num or periods        , cpn dates */
    double *a, double *marg,
    double lvg,            /* handle        , margin        , leverage */
    SrtBasisCode cpnbasis, /* coupon basis */
    SrtBasisCode basis)    /* basis for rates */
{
  SrtCurvePtr yldcrv = NULL; /* yield curve pointer */
  LGMMarkConv conventions;   /* standard market conventions */
  LGMErr error = NULL;
  Date tEx;
  long i, lag;
  double CEVvol, beta, Vfloor, pvCpn;
  double df, dfSt, dfEnd, cvgCap, cvgCpn;
  double sum, ratio, intVal;
  double Forward, Level, Yrtoexp;

  yldcrv = lookup_curve(ycName);
  error = LGMCcyDefaults(yldcrv, &conventions); /* get currency conventions */
  if (error != NULL)
    return (error);
  lag = conventions.lag; /* get lag */
  if (n < 1)
    return ("no coupons");

  sum = 0.;
  for (i = 1; i <= n; i++) /* find PV of coupon leg */
  {
    dfEnd = swp_f_df(tNow, t[i], ycName); /* calculate value of floor */
    dfSt = swp_f_df(tNow, t[i - 1], ycName);
    if (dfSt == SRT_DF_ERROR || dfEnd == SRT_DF_ERROR)
      return ("no discount factor");
    cvgCap = coverage(t[i - 1], t[i], basis);
    tEx = add_unit(t[i - 1], -lag, SRT_BDAY, SUCCEEDING);

    error = GetVol(t[i - 1], t[i], a[i], SRT_TRUE, &CEVvol); /* Get CEVvol */
    if (error)
      return (error);
    error = GetBeta(t[i - 1], t[i], &beta); /* Get CEV exponent */
    if (error)
      return (error);

    Yrtoexp = (tEx - tNow) * YEARS_IN_DAY;
    Level = cvgCap * dfEnd;
    Forward = (dfSt - dfEnd) / (Level);

    /* Vfloor = LGMCEVCapletGreek(SRT_RECEIVER        ,Forward        ,a[i]
     * ,Level ,Yrtoexp        ,CEVvol        , beta        ,PREMIUM); */

    Vfloor = LGMCEVCapletPrice(tNow, tEx, SRT_RECEIVER, a[i], cvgCap, dfSt,
                               dfEnd, CEVvol, beta);
    if (Vfloor < 0.0)
      return ("can't get floor price");

    cvgCpn = coverage(t[i - 1], t[i], cpnbasis);
    pvCpn = cvgCpn * (marg[i] * dfEnd + lvg * Vfloor / cvgCap);
    sum = sum + pvCpn;
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

  if (ratio > 5.0 || ratio < 0.2) /* clearly        , one sign deal */
    ratio = 0.;

  *ratioPtr = ratio;
  *intValPtr = intVal;
  return (NULL);
}

/********************************************************************************/
/*	Returns intrinsic value and the ratio PV of coupon leg to PV of the
strike for a cap floater. If deal is strictly positive or strictly negative ,
ratio is set to 0. */
static LGMErr CompValCapFlt(
    Date tNow, String ycName, /* yield curve */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double *), /* cap vols */
    LGMErr (*GetBeta)(Date, Date, double *), /* cap exponent (beta) */
    double *ratioPtr,                        /* ratio */
    double *intValPtr,                       /* intrinsic value */
    SrtReceiverType PayRec,                  /* pay/rec */
    Date tStart, double Strike,              /* settlement        , strike */
    long n, Date *t, /* num or periods        , cpn dates */
    double *amax, double *amin,
    double *marg,          /* handle max and min        , margin*/
    double lvg,            /* leverage */
    SrtBasisCode cpnbasis, /* coupon basis */
    SrtBasisCode basis)    /* basis for rates */
{
  SrtCurvePtr yldcrv = NULL; /* yield curve pointer */
  LGMMarkConv conventions;   /* standard market conventions */
  LGMErr error = NULL;
  Date tEx;
  long i, lag;
  double CEVvolmax, CEVvolmin, beta, Vcapmax, Vcapmin, pvCpn;
  double df, dfSt, dfEnd, cvgCap, cvgCpn;
  double Level, Forward, Yrtoexp;
  double sum, ratio, intVal;

  yldcrv = lookup_curve(ycName);
  error = LGMCcyDefaults(yldcrv, &conventions); /* get currency conventions */
  if (error != NULL)
    return (error);
  lag = conventions.lag; /* get lag */
  if (n < 1)
    return ("no coupons");

  sum = 0.;
  for (i = 1; i <= n; i++) /* find PV of coupon leg */
  {
    dfEnd = swp_f_df(tNow, t[i], ycName); /* calculate value of floor */
    dfSt = swp_f_df(tNow, t[i - 1], ycName);
    if (dfSt == SRT_DF_ERROR || dfEnd == SRT_DF_ERROR)
      return ("no discount factor");
    cvgCap = coverage(t[i - 1], t[i], basis);
    tEx = add_unit(t[i - 1], -lag, SRT_BDAY, SUCCEEDING);

    error = GetVol(t[i - 1], t[i], amax[i], SRT_TRUE,
                   &CEVvolmax); /* Get CEVvol for cap max*/
    if (error)
      return (error);
    error = GetVol(t[i - 1], t[i], amin[i], SRT_TRUE,
                   &CEVvolmin); /* Get CEVvol for cap min */
    if (error)
      return (error);
    error = GetBeta(t[i - 1], t[i], &beta); /* Get CEV exponent */
    if (error)
      return (error);

    Yrtoexp = (tEx - tNow) * YEARS_IN_DAY;
    Level = cvgCap * dfEnd;
    Forward = (dfSt - dfEnd) / (Level);

    Vcapmax = LGMCEVCapletPrice(tNow, tEx, SRT_PAYER, amax[i], cvgCap, dfSt,
                                dfEnd, CEVvolmax, beta);

    if (Vcapmax < 0.0)
      return ("can't get cap max price");
    Vcapmin = LGMCEVCapletPrice(tNow, tEx, SRT_PAYER, amin[i], cvgCap, dfSt,
                                dfEnd, CEVvolmin, beta);

    cvgCpn = coverage(t[i - 1], t[i], cpnbasis);
    pvCpn = cvgCpn * (marg[i] * dfEnd + lvg * (Vcapmin - Vcapmax) / cvgCap);
    sum = sum + pvCpn;
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

  if (ratio > 5.0 || ratio < 0.2) /* clearly        , one sign deal */
    ratio = 0.;

  *ratioPtr = ratio;
  *intValPtr = intVal;
  return (NULL);
}

/************************************************************************/
/* Convert xExBdry[i] = x(t) at tExBdry[i] into equivalent swaptions which
would be ATM at the exercise boundary        , for i = 0        ,... ,nEx-1 */
/* On error        , return NULL */
static LGMSwptnsPtr FindBdry(int iLGMOneTwoFactor,
                             Date tNow,       /* Evaluation date */
                             String ycName,   /* information about today */
                             LGM_TS *tsPtr,   /* term structure */
                             long nEx,        /* number of boundary points*/
                             Date *tExBdry,   /* exercise boundary dates */
                             double *xExBdry, /* exercise boundary (in x) */
                             Date tEnd)       /* last pay date */

{
  LGMErr error = NULL;
  LGMSwptns *ExerBdryPtr;
  SrtCurvePtr yldcrv = NULL;
  LGMMarkConv conv;
  Date tEx, tStart, tPay;
  double Rfix;
  double FloatRate;
  long j, k;

  yldcrv = lookup_curve(ycName);
  error = LGMCcyDefaults(yldcrv, &conv);
  if (error)
    return (NULL);

  /* allocate space */
  ExerBdryPtr = LGMCreateSwptns(nEx);
  if (ExerBdryPtr == NULL)
    return (NULL);

  k = -1;
  for (j = 0; j < nEx; j++) {
    tEx = tExBdry[j];
    if (tEx >= tNow && tEx < tEnd - 30) {
      if (iLGMOneTwoFactor == 1) {
        // error = LGMGetRfixFromX(tEx        , tEnd        , xExBdry[j] , tNow
        // , ycName
        //       , tsPtr        , &Rfix); if(error) return NULL;

        tStart = add_unit(tEx, conv.lag, SRT_BDAY, SUCCEEDING);
        tPay =
            add_unit(tStart, 12 / (int)conv.cxxfreq, SRT_MONTH, conv.cxxbdconv);

        error = LGMGetFloatRateFromX(tStart, tPay, xExBdry[j], ycName, tsPtr,
                                     &FloatRate);
        if (error)
          return NULL;

      } else
        Rfix = xExBdry[j];
      if (error == NULL) {
        k++;
        ExerBdryPtr->tEx[k] = tEx;
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

#undef INTRINSICDEFAULT
