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
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"
#include "swp_h_vol.h"
/****************************************************************************/
/* Prototypes of private functions */
/****************************************************************************/
static void TestOverride(LGMCalParm *CalReqPtr);
static LGMErr copyExerBdry(LGMSwptns *ExerBdryPtr,
                           SrtLgmExerBdryData *lgmExerBdryData);
static LGMErr copyTSData(SigKapTS *SKtsPtr, SrtLgmTSData *TSData);

/* Free arrays that will be used for output  , if already allocated */
static void free_outputs(
    int *convertTS,    /* 1=compute new sigs and taus; 0=don't bother */
    int *findExerBdry, /* 1=find swap rates at exercise boundary; 0=don't bother
                        */
    SrtLgmExerBdryData *lgmExerBdryData, /* ptr to exercise boundary data
                                            structure (NULL => not req'd) */
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
    SrtLgmTSData *lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr *atcTSData);   /* ptr to zeta/G data (NULL => not req'd) */

/****************************************************************************/
/**** Comments ****/
/* NewLGMautocalCaller is a temporizing step  , which does not use the
flexibility of the new program. It should be replaced with programs which call
the new LGMautocal directly from MAD and Westminster  , which can expose the new
functionality  , ASAP */

/* In creating the call to the new LGMautocal  , this routine ASSUMES that
the deal is exactly as follows. Let
        tfirst = tNow if end_of_day is off
                   = tNow + 1 business day if end of day is on

The option holder can exercise only once  , and can exercise on any of the dates
tEx[0]  , tEx[1]  , ...  , tEx[nEx-1] which are on or after tfirst.

Suppose the option holder exercises at tEx[j]  , and let tPay[i0] be the first
fixed leg paydate which is strictly after tStart[j]. Then the Receiver gets
the fixed payment stream
        Payment[i0] - RedFirstPay[j] paid at tPay[i0]  ,
        Payment[i] paid at tPay[i]  , i=i0+1  ,...  ,nPay-1.
In return the Receiver pays
        Strike[j] paid at tStart[j].
Note that the last fixed leg payment Payment[nPay-1] should include the notional
and the value of the floating leg should be in Strike[j]. Accruals for
inter-period starts can be handled by adjusting Strike[j] and/or RedFirstPay[j]
, depending on the deal  , and adjustments for the basis of the floating leg can
be handled by adjusting Strike[j] and Payment[i] */

/* Output
The outputs are all of the form *variable = answer for numbers and **variable =
*Array for arrays. The calling program must allocate space which is being
pointed to by *variable. Specifically  , *LGMPtr  , *nExerBdry  , *NsigPtr  ,
and *NtauPtr must point to space that has been allocated for double  , long  ,
long  , and long which can receive the outputs of the program;

*rArrExBdryPtr  , *tArrExBdryPtr  , *sigDateArrPtr  , *sigValArrPtr  ,
*tauDateArrPtr  , and *tauValArrPtr must point to space the has been allocated
for pointers to doubles and Dates; the program will free any space which is
being pointed to and reallocated it. Then *variable will be set to the address
of first element in the array of answers.

All output arrays which are not NULL are freed with srt_free and then
reallocated with srt_calloc. This will result in TERRIBLE errors if they have
not been allocated with calloc  , malloc  , srt_calloc  , or srt_malloc.
Specifically  , if the arrays *rArrExBdryPtr  , *tArrExBdryPtr  , *sigDateArrPtr
, *sigValArrPtr  , *tauDateArrPtr  , and *tauValArrPtr are not NULL  , they need
to have been allocated with calloc or malloc.

If findExerBdry is on  , the array *rArrExBdryPtr[i] gives the swap rate at the
exercise boundary for a swap that has exercise date *tArrExBdryPtr[i] and ends
at the final paydate tPay[nPay-1] of the underlying deal; all other parameters
(spot_lag  , daycountconv  , etc.) come from the currency code.

If convertTS is on  , the program finds the sigma-tau term structure equivalent
to the calibrated zeta-G term structure  , and outputs it in the arrays
*sigDateArrPtr and *sigValArrPtr (of length *NsigPtr) and *tauDateArrPtr and
*tauValArrPtr (of length *NtauPtr).	*/

/**** HOLIDAYS ***/
/* Wherever we need to use the add_unit function  ,
the currency code ccy is available */

/************************************************************************************/
/* This call is identical to the old LGMamerican call EXCEPT that the argument
rflt_current (rate for the current floating period) has been added (the program
ignores this argument if it is zero or negative) as well as the calibration
method arguments at the end. */

LGMErr TestLGMAmerCaller(
    /* definition of American deal */
    long nfix, /* number of fixed periods */
    Date
        tfixStart[], /* [0  ,1  ,...  ,nfix-1]  start dates for fixed coupons */
    Date tfixEnd[],  /* [0  ,1  ,...  ,nfix-1]  end dates for fixed coupons */
    Date tfixPay[],  /* [0  ,1  ,...  ,nfix-1]  pay dates for fixed coupons */
    double fixFullPayment[], /* [0  ,1  ,...  ,nfix-1]  fixed coupons (with
                                notional at end) */
    double Strike[],  /* [0  ,1  ,...  ,nfix-1]  premium for exer. at j (with
                         notional) */
    char *fixBasis,   /* basis for fixed coupons */
    int fixEarlyFlag, /* 0=subtract accrual from frst pymnt; 1=add accrual to
                         fee */
    long nfltdates,   /* number of floating dates */
    Date tflt[], /* [0  ,1  ,...  ,nfltdates-1] all flting dates (first start
                    date & all pay dates) */
    double rflt_current, /* current floating rate (if fixed); ignored if
                            non-positive */
    char *fltBasis,      /* basis for floating coupons */
    int fltEarlyFlag, /* 0 = subtrct accrual from frst pymnt; 1 = add accrual to
                         prem */
    Date tFirstExer,  /* first exercise date */
    int lagExerStart, /* days between exercise and start */
    int cal_or_bus,   /* 0 = cal. days  , 1 = bus. days for lag_exer_start */
    BusDayConv convStart, /* business day convention for start (typically none
                             or suceeding) */
    char *PayRecStr,      /* RECEIVER or PAYER */
                          /* info about today's marketplace */
    String ycName,        /* Yield curve name */
    int endofdayflag, /* 1=too late to exercise deals today  , 0=not too late */
    LGMErr (*GetVol)(long, long, double, SRT_Boolean,
                     double *), /* volatility function for reference swptns */
    char *char_vol_type,        /* normal or log normal swaption vols */
                                /* calibration method to use */
    int LGMOneTwoFactor, int usefixtau, /* 1=calibrate with fixed tau */
    int usecaps, /* 1=use caplets for calibration  , 0=use only swaptions */
    double tau,  /* if fixed tau  , use this value for tau (in years) */
    double alpha, double gamma, double rho,
    int calibrationmeth, /* 1 = fixexp  , 2=backboot  , 3=fixed sigma */
    int strikechoice,    /* 1 = IRR  , 2 = dIRR   , 3 = EMK1  , 4 = EMK2 */
                         /* task list */
    int skipEval,        /* 1=calibrate only  , 0=calibrate & evaluate deal */
    int convertTS, /* 1=cmpute new sigmas  , taus & store in ts; 0=don't bother
                    */
    int findExerBdry,  /* 1=find swap rates at exercise boundary; 0=don't bother
                        */
                       /* outputs */
    String outfile,    /* output log file (unused) */
    double *LGMValPtr, /* answer */
    SrtLgmExerBdryData *lgmExerBdryData, /* ptr to exercise boundary data
                                            structure (NULL => not req'd) */
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
                                         /* calibrated term structure */
    SrtLgmTSData *lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
                             /* miscellaneous */
    double *intrinsicValPtr) /* ptr to intrinsic value (NULL => not req'd) */
{
  return TestLGMAmerCallerRR(
      nfix,      /* number of fixed periods */
      tfixStart, /* [0  ,1  ,...  ,nfix-1]  start dates for fixed coupons */
      tfixEnd,   /* [0  ,1  ,...  ,nfix-1]  end dates for fixed coupons */
      tfixPay,   /* [0  ,1  ,...  ,nfix-1]  pay dates for fixed coupons */
      fixFullPayment, /* [0  ,1  ,...  ,nfix-1]  fixed coupons (with notional at
                         end) */
      Strike, /* [0  ,1  ,...  ,nfix-1]  premium for exer. at j (with notional)
               */
      fixBasis,     /* basis for fixed coupons */
      fixEarlyFlag, /* 0=subtract accrual from frst pymnt; 1=add accrual to fee
                     */
      nfltdates,    /* number of floating dates */
      tflt, /* [0  ,1  ,...  ,nfltdates-1] all flting dates (first start date &
               all pay dates) */
      rflt_current, /* current floating rate (if fixed); ignored if non-positive
                     */
      fltBasis,     /* basis for floating coupons */
      fltEarlyFlag, /* 0 = subtrct accrual from frst pymnt; 1 = add accrual to
                       prem */
      tFirstExer,   /* first exercise date */
      lagExerStart, /* days between exercise and start */
      cal_or_bus,   /* 0 = cal. days  , 1 = bus. days for lag_exer_start */
      convStart,    /* business day convention for start (typically none or
                       suceeding) */
      PayRecStr,    /* RECEIVER or PAYER */
      ycName,       /* Yield curve name */
      endofdayflag, /* 1=too late to exercise deals today  , 0=not too late */
      GetVol,       /* volatility function for reference swptns */
      char_vol_type,              /* normal or log normal swaption vols */
      LGMOneTwoFactor, usefixtau, /* 1=calibrate with fixed tau */
      usecaps, /* 1=use caplets for calibration  , 0=use only swaptions */
      tau,     /* if fixed tau  , use this value for tau (in years) */
      alpha, gamma, rho,
      calibrationmeth, /* 1 = fixexp  , 2=backboot  , 3=fixed sigma */
      strikechoice,    /* 1 = IRR  , 2 = dIRR   , 3 = EMK1  , 4 = EMK2 */
      skipEval,        /* 1=calibrate only  , 0=calibrate & evaluate deal */
      convertTS, /* 1=cmpute new sigmas  , taus & store in ts; 0=don't bother */
      findExerBdry, /* 1=find swap rates at exercise boundary; 0=don't bother */
      outfile,      /* output log file (unused) */
      LGMValPtr,    /* answer */
      lgmExerBdryData, /* ptr to exercise boundary data structure (NULL => not
                          req'd) */
      lgmRefSwptnData, /* ptr to reference swaption data structure (NULL => not
                          req'd) */
      lgmTSData,       /* ptr to tau/sigma data (NULL => not req'd) */
      intrinsicValPtr, /* ptr to intrinsic value (NULL => not req'd) */
      0);              /* optional refrate included for MUNI pricing */
}

/************************************************************************************/
/************************************************************************************/
/* Routine that takes in refrate for pricing muni americans */
LGMErr TestLGMAmerCallerRR(
    /* definition of American deal */
    long nfix, /* number of fixed periods */
    Date
        tfixStart[], /* [0  ,1  ,...  ,nfix-1]  start dates for fixed coupons */
    Date tfixEnd[],  /* [0  ,1  ,...  ,nfix-1]  end dates for fixed coupons */
    Date tfixPay[],  /* [0  ,1  ,...  ,nfix-1]  pay dates for fixed coupons */
    double fixFullPayment[], /* [0  ,1  ,...  ,nfix-1]  fixed coupons (with
                                notional at end) */
    double Strike[],  /* [0  ,1  ,...  ,nfix-1]  premium for exer. at j (with
                         notional) */
    char *fixBasis,   /* basis for fixed coupons */
    int fixEarlyFlag, /* 0=subtract accrual from frst pymnt; 1=add accrual to
                         fee */
    long nfltdates,   /* number of floating dates */
    Date tflt[], /* [0  ,1  ,...  ,nfltdates-1] all flting dates (first start
                    date & all pay dates) */
    double rflt_current, /* current floating rate (if fixed); ignored if
                            non-positive */
    char *fltBasis,      /* basis for floating coupons */
    int fltEarlyFlag, /* 0 = subtrct accrual from frst pymnt; 1 = add accrual to
                         prem */
    Date tFirstExer,  /* first exercise date */
    int lagExerStart, /* days between exercise and start */
    int cal_or_bus,   /* 0 = cal. days  , 1 = bus. days for lag_exer_start */
    BusDayConv convStart, /* business day convention for start (typically none
                             or suceeding) */
    char *PayRecStr,      /* RECEIVER or PAYER */
                          /* info about today's marketplace */
    String ycName,        /* Yield curve name */
    int endofdayflag, /* 1=too late to exercise deals today  , 0=not too late */
    LGMErr (*GetVol)(long, long, double, SRT_Boolean,
                     double *), /* volatility function for reference swptns */
    char *char_vol_type,        /* normal or log normal swaption vols */
                                /* calibration method to use */
    int LGMOneTwoFactor, int usefixtau, /* 1=calibrate with fixed tau */
    int usecaps, /* 1=use caplets for calibration  , 0=use only swaptions */
    double tau,  /* if fixed tau  , use this value for tau (in years) */
    double alpha, double gamma, double rho,
    int calibrationmeth, /* 1 = fixexp  , 2=backboot  , 3=fixed sigma */
    int strikechoice,    /* 1 = IRR  , 2 = dIRR   , 3 = EMK1  , 4 = EMK2 */
                         /* task list */
    int skipEval,        /* 1=calibrate only  , 0=calibrate & evaluate deal */
    int convertTS, /* 1=cmpute new sigmas  , taus & store in ts; 0=don't bother
                    */
    int findExerBdry,  /* 1=find swap rates at exercise boundary; 0=don't bother
                        */
                       /* outputs */
    String outfile,    /* output log file (unused) */
    double *LGMValPtr, /* answer */
    SrtLgmExerBdryData *lgmExerBdryData, /* ptr to exercise boundary data
                                            structure (NULL => not req'd) */
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
                                         /* calibrated term structure */
    SrtLgmTSData *lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
                             /* miscellaneous */
    double *intrinsicValPtr, /* ptr to intrinsic value (NULL => not req'd) */
    char *szRefRate)         /* optional refrate included for MUNI pricing */
{
  LGMAmerType ItsAnAmer;
  SrtSimAmer *AmerPtr = NULL;
  SrtReceiverType srt_p_r_flag;
  SrtDiffusionType srt_vol_type;
  SrtBasisCode fixBasisCode, fltBasisCode;
  LGMErr (*GetBeta)(Date, Date, double *);
  LGMErr (*RecVol)(Date, Date, double);
  LGMMarkConv conv;
  LGMCalParm *CalReqPtr = NULL;
  ConvParams *EvalParamsPtr = NULL;
  double intrinsicVal = 0;
  double LGMVal = 0;
  LGMSwptns *ExerIntoPtr = NULL;
  LGMSwptns *ExerBdryPtr = NULL;
  LGM_TS *LGMtsPtr = NULL;
  SigKapTS *SigKaptsPtr = NULL;
  Date tNow, tfirst;
  SrtCurvePtr yldcrv = NULL;
  int standardLag;

  LGMErr error = NULL;
  int skipCal;

  /* free arrays that will be used for output  , if already allocated */
  if (findExerBdry != 0 && lgmExerBdryData != NULL) {
    lgmExerBdryData->NexerBdry = 0;
    if (lgmExerBdryData->exerBdryArr != NULL)
      srt_free(lgmExerBdryData->exerBdryArr);
  } else
    findExerBdry = 0;

  if (convertTS != 0 && lgmTSData != NULL) {
    lgmTSData->NTS = 0;
    if (lgmTSData->TSArr != NULL)
      srt_free(lgmTSData->TSArr);
  } else
    convertTS = 0;

  /* Set operations to be done by LGMnewautocal */
  /* convertTS  , findExerBdry available from input */
  skipCal = 0; /* calibrate & evaluate */

  /* Set information about today  , and today's market */
  yldcrv = lookup_curve(ycName);
  tNow = get_clcndate_from_yldcrv(yldcrv); /* get calculation date */

  /* Modified code to take account of MJNi Case */
  if ((szRefRate != 0) && (strcmp(szRefRate, "UTFBM3") == 0))
    error = LGMMuniDefaults(szRefRate, &conv);
  else
    error = LGMCcyDefaults(yldcrv, &conv); /* get currency conventions  , lag */
                                           /* end modified code */

  if (error != NULL)
    return (error);
  standardLag = conv.lag;

  /* endofday and GetVol available from input */

  /* get the exponent beta to use to interpret market vols */
  /* transform char_vol_type to SrtDiffusionType */
  error = interp_diffusion_type(char_vol_type, &srt_vol_type);
  if (error != NULL)
    return (error);

  if (srt_vol_type == SRT_LOGNORMAL)
    GetBeta = LGMReturnExp1; /* This will return beta=1 for all swaptions */
  else
    GetBeta = LGMReturnExp0; /* This will return beta=0 for all swaptions */

  /* determine the first possible day of exercise */
  tfirst = tNow;
  if (endofdayflag == 1)
    tfirst = add_unit(tNow, 1, SRT_BDAY, SUCCEEDING);

  /* interpret pay-receive and the day count bases for the fixed and floating
   * legs */
  error = interp_rec_pay(PayRecStr, &srt_p_r_flag);
  if (error != NULL)
    return (error);

  error = interp_basis(fixBasis, &fixBasisCode);
  if (error)
    return (error);
  error = interp_basis(fltBasis, &fltBasisCode);
  if (error)
    return (error);

  /* Set deal type */
  ItsAnAmer = SimAmer;

  /* create simple American stucture  , check input data  , & copy input data
   * into structure */
  error =
      LGMFillSimAmer(tfirst, standardLag, &AmerPtr, nfix, tfixStart,
                     tfixEnd, /* fixed leg info */
                     tfixPay, fixFullPayment, fixBasisCode, fixEarlyFlag,
                     Strike, /* float leg info */
                     nfltdates, tflt, fltBasisCode, fltEarlyFlag, tFirstExer,
                     lagExerStart, cal_or_bus, convStart, /* exercise info */
                     srt_p_r_flag);
  if (error != NULL) {
    LGMFreeSimAmer(&AmerPtr);
    return (error);
  }

  /* Set calibration method */
  CalReqPtr = LGMSetCalibMeth(LGMOneTwoFactor, usefixtau, usecaps, tau, alpha,
                              gamma, rho, 2.0, NULL, NULL, NULL, NULL, NULL);

  if (CalReqPtr == NULL) {
    LGMFreeSimAmer(&AmerPtr);
    return ("alloc failed in LGMcaller");
  }
  if (LGMOneTwoFactor == 1) {
    if (calibrationmeth == 1 && usefixtau != 1)
      CalReqPtr->calmeth = FixExp;
    else if (calibrationmeth == 2 && usefixtau != 1)
      CalReqPtr->calmeth = TenorAndDiag;
    else if (calibrationmeth == 3 && usefixtau != 1)
      CalReqPtr->calmeth = FixSigma;
  }
  switch (strikechoice) {
  case 1:
    CalReqPtr->Rmeth = IRR;
    break;

  case 2:
    CalReqPtr->Rmeth = dIRR;
    break;

  case 3:
    CalReqPtr->Rmeth = EMK1;
    break;

  case 4:
    CalReqPtr->Rmeth = EMK2;
    break;
  }

  /*	TestOverride(calibrationmeth  , strikechoice  , CalReqPtr); */

  /* Set numerical parameters for convolution to their default values */
  EvalParamsPtr = LGMSetDefaultEvalParms();
  if (EvalParamsPtr == NULL) {
    LGMFreeSimAmer(&AmerPtr);
    LGMFreeCalParm(&CalReqPtr);
    return ("alloc failed in LGMcaller");
  }

  /* Record Vol */
  RecVol = LGMRecVolDummy;

  /* Call LGMautocal */
  error = LGMNewAmerican(
      skipEval, skipCal, convertTS, findExerBdry, /* task list */
      tNow, endofdayflag, ycName, /* info about today's market */
      rflt_current, /* current floating rate (if set) ignored if non-positive */
      GetVol, GetBeta, RecVol,    /* swaption & caplet price info */
      ItsAnAmer, (void *)AmerPtr, /* the deal */
      CalReqPtr, EvalParamsPtr,   /* calib & eval methods */
      &LGMVal, &intrinsicVal,     /* value & intrinsic value of deal */
      &ExerIntoPtr,
      &ExerBdryPtr, /* swaptions representing the underlying & exer boundary */
      lgmRefSwptnData,         /* reference swaptions */
      &LGMtsPtr, &SigKaptsPtr, /* calibrated term structures */
      szRefRate);

  /* free unneeded structures */ /* these are not needed to get output */
  LGMFreeSimAmer(&AmerPtr);      /* structure containing the deal */
  LGMFreeCalParm(
      &CalReqPtr);         /* structure containing the requested cal method */
  srt_free(EvalParamsPtr); /* structure containing the convolution parameters */
  LGMFreeSwptns(
      &ExerIntoPtr); /* structure containing underlying swaptions of MidAt */
  if (LGMOneTwoFactor == 1)
    LGMFreeLGM_TS(&LGMtsPtr);
  if (LGMOneTwoFactor == 2)
    LGMFreeLGM2F_TS(&LGMtsPtr);

  if (error != NULL) {
    LGMFreeSwptns(&ExerBdryPtr);
    LGMFreeSigKapTS(&SigKaptsPtr);
    return (error);
  }

  /* Unpack output */
  *LGMValPtr = LGMVal;
  if (intrinsicValPtr != NULL)
    *intrinsicValPtr = intrinsicVal;

  /* If exercise boundary has been found  , copy it into tExBdry & rExBdry */
  if (skipEval != 1 && findExerBdry != 0 && ExerBdryPtr != NULL &&
      ExerBdryPtr->n > 0 && lgmExerBdryData != NULL)
    error = copyExerBdry(ExerBdryPtr, lgmExerBdryData);
  if (error != NULL) {
    LGMFreeSwptns(&ExerBdryPtr);
    LGMFreeSigKapTS(&SigKaptsPtr);
    return (error);
  }
  LGMFreeSwptns(&ExerBdryPtr); /* no longer needed */

  /* If sigma-kappa term structure has been found  , copy it into output arrays
   */
  if (skipEval != 1 && convertTS != 0 && SigKaptsPtr != NULL &&
      SigKaptsPtr->numS > 0 && SigKaptsPtr->numK > 0 &&
      SigKaptsPtr->sdate != NULL && SigKaptsPtr->sig != NULL &&
      SigKaptsPtr->kdate != NULL && SigKaptsPtr->kap != NULL &&
      lgmTSData != NULL) {
    error = copyTSData(SigKaptsPtr, lgmTSData);
    LGMFreeSigKapTS(&SigKaptsPtr);
    if (error != NULL)
      return (error);
  }
  LGMFreeSigKapTS(&SigKaptsPtr); /* no longer needed */

  return (error);
}

/************************************************************************************/
/************************************************************************************/
/* Program to calibrate the term structure and compute value of a MidAtlantic
deal. Program has the identical call as the old LGMautocal. It uses its inputs
to create the call to the LGMnewautocal program. This file now calls a new
routine with an additional argument */
LGMErr TestLGMautocalCaller(
    long nEx,       /* nEx is number of exercises */
    Date *tEx,      /* notification (exercise) dates  , [0  ,1  ,...  ,nEx-1] */
    Date *tStart,   /* start dates for each exercise  , [0  ,1  ,...  ,nEx-1] */
    double *Strike, /* total value paid at tStart[j] for fixed leg  , [0  ,1
                       ,...  ,nEx-1] */
    long nPay,      /* nPay is number of fixed leg coupons */
    Date *tPay,     /* pay dates for period i  , [0  ,1  , ...  ,nPay-1] */
    double *Payment, /* total fixed leg payment(last includes notional)  , [0
                        ,...  ,nPay-1] */
    double *RedFirstPay, /* reduction in 1rst payment after exercise  , [0  ,...
                            ,nEx-1] */
    char *PayRecStr,     /* RECEIVER or PAYER */
                         /* information about today */
    String ycName,       /* pointer to market structures */
    int endofday, /* 1=too late to exercise deals today  , 0=not too late */
                  /* today's volatilities */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to get swaption vols */
    char *char_vol_type, /* determines whether vol is normal or log normal */
                         /* calibration method to use */
    int LGMOneTwoFactor, int usefixtau, /* 1=calibrate with fixed tau */
    int usecaps, /* 1=use caplets for calibration  , 0=use only swaptions */
    double tau,  /* if fixed tau  , use this value for tau (in years) */
    double alpha, double gamma, double rho,
    int calibrationmeth, /* 1 = fixexp  , 2=backboot  , 3=fixed sigma */
    int strikechoice,    /* 1 = IRR  , 2 = dIRR */
    double maxstd,       /* Maximum number of std between forward and strike */
                         /* requested operation */
    int skipEval,        /* 1=calibrate only  , 0=calibrate & evaluate deal */
    int convertTS,       /* 1=compute new sigs and taus; 0=don't bother */
    int findExerBdry, /* 1=find swap rates at exercise boundary; 0=don't bother
                       */
    long *Zeta1Dates, double *StartZeta1s, long *TauDates, double *StartTaus,
    double **HybridShortInstrsIndex,

    /* outputs */
    String outfile,    /* output file name for log file (unused) */
    double *LGMvalPtr, /* LGM value of mid-atlantic */
    SrtLgmExerBdryData *lgmExerBdryData, /* ptr to exercise boundary data
                                            structure (NULL => not req'd) */
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
                                         /* calibrated term structure */
    SrtLgmTSData *lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr *atcTSData,    /* ptr to zeta/G data (NULL => not req'd) */
                             /* miscellaneous */
    double *intrinsicValPtr) /* ptr to intrinsic value (NULL => not req'd) */
{
  return TestLGMautocalCallerRR(
      nEx,     /* nEx is number of exercises */
      tEx,     /* notification (exercise) dates  , [0  ,1  ,...  ,nEx-1] */
      tStart,  /* start dates for each exercise  , [0  ,1  ,...  ,nEx-1] */
      Strike,  /* total value paid at tStart[j] for fixed leg  , [0  ,1  ,...
                  ,nEx-1] */
      nPay,    /* nPay is number of fixed leg coupons */
      tPay,    /* pay dates for period i  , [0  ,1  , ...  ,nPay-1] */
      Payment, /* total fixed leg payment(last includes notional)  , [0  ,...
                  ,nPay-1] */
      RedFirstPay,   /* reduction in 1rst payment after exercise  , [0  ,...
                        ,nEx-1] */
      PayRecStr,     /* RECEIVER or PAYER */
      ycName,        /* pointer to market structures */
      endofday,      /* 1=too late to exercise deals today  , 0=not too late */
      GetVol,        /* function to get swaption vols */
      char_vol_type, /* determines whether vol is normal or log normal */
      LGMOneTwoFactor, usefixtau, /* 1=calibrate with fixed tau */
      usecaps, /* 1=use caplets for calibration  , 0=use only swaptions */
      tau,     /* if fixed tau  , use this value for tau (in years) */
      alpha, gamma, rho,
      calibrationmeth, /* 1 = fixexp  , 2=backboot  , 3=fixed sigma */
      strikechoice,    /* 1 = IRR  , 2 = dIRR */
      maxstd,          /* Maximum number of std between forward and strike */
      skipEval,        /* 1=calibrate only  , 0=calibrate & evaluate deal */
      convertTS,       /* 1=compute new sigs and taus; 0=don't bother */
      findExerBdry, /* 1=find swap rates at exercise boundary; 0=don't bother */
      Zeta1Dates, StartZeta1s, TauDates, StartTaus, HybridShortInstrsIndex,
      outfile,         /* output file name for log file (unused) */
      LGMvalPtr,       /* LGM value of mid-atlantic */
      lgmExerBdryData, /* ptr to exercise boundary data structure (NULL => not
                          req'd) */
      lgmRefSwptnData, /* ptr to reference swaption data structure (NULL => not
                          req'd) */
      lgmTSData,       /* ptr to tau/sigma data (NULL => not req'd) */
      atcTSData,       /* ptr to zeta/G data (NULL => not req'd) */
      intrinsicValPtr, /* ptr to intrinsic value (NULL => not req'd) */
      0);
}

LGMErr TestLGMautocalCallerRR(
    long nEx,       /* nEx is number of exercises */
    Date *tEx,      /* notification (exercise) dates  , [0  ,1  ,...  ,nEx-1] */
    Date *tStart,   /* start dates for each exercise  , [0  ,1  ,...  ,nEx-1] */
    double *Strike, /* total value paid at tStart[j] for fixed leg  , [0  ,1
                       ,...  ,nEx-1] */
    long nPay,      /* nPay is number of fixed leg coupons */
    Date *tPay,     /* pay dates for period i  , [0  ,1  , ...  ,nPay-1] */
    double *Payment, /* total fixed leg payment(last includes notional)  , [0
                        ,...  ,nPay-1] */
    double *RedFirstPay, /* reduction in 1rst payment after exercise  , [0  ,...
                            ,nEx-1] */
    char *PayRecStr,     /* RECEIVER or PAYER */
                         /* information about today */
    String ycName,       /* pointer to market structures */
    int endofday, /* 1=too late to exercise deals today  , 0=not too late */
                  /* today's volatilities */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to get swaption vols */
    char *char_vol_type, /* determines whether vol is normal or log normal */
                         /* calibration method to use */
    int LGMOneTwoFactor, int usefixtau, /* 1=calibrate with fixed tau */
    int usecaps, /* 1=use caplets for calibration  , 0=use only swaptions */
    double tau,  /* if fixed tau  , use this value for tau (in years) */
    double alpha, double gamma, double rho,
    int calibrationmeth, /* 1 = fixexp  , 2=backboot  , 3=fixed sigma */
    int strikechoice,    /* 1 = IRR  , 2 = dIRR */
    double maxstd,       /* Maximum number of std between forward and strike */
                         /* requested operation */
    int skipEval,        /* 1=calibrate only  , 0=calibrate & evaluate deal */
    int convertTS,       /* 1=compute new sigs and taus; 0=don't bother */
    int findExerBdry, /* 1=find swap rates at exercise boundary; 0=don't bother
                       */
    long *Zeta1Dates, double *StartZeta1s, long *TauDates, double *StartTaus,
    double **HybridShortInstrsIndex,

    /* outputs */
    String outfile,    /* output file name for log file (unused) */
    double *LGMvalPtr, /* LGM value of mid-atlantic */
    SrtLgmExerBdryData *lgmExerBdryData, /* ptr to exercise boundary data
                                            structure (NULL => not req'd) */
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
                                         /* calibrated term structure */
    SrtLgmTSData *lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr *atcTSData,    /* ptr to zeta/G data (NULL => not req'd) */
                             /* miscellaneous */
    double *intrinsicValPtr, /* ptr to intrinsic value (NULL => not req'd) */
    char *szRefRate) {
  LGMDealType ItsAMidAt;
  SrtSimMidAt *MidAtPtr = NULL;
  SrtReceiverType srt_p_r_flag;
  SrtDiffusionType srt_vol_type;
  LGMErr (*GetBeta)(Date, Date, double *);
  LGMErr (*RecVol)(Date, Date, double);
  LGMCalParm *CalReqPtr = NULL;
  ConvParams *EvalParamsPtr = NULL;
  double intrinsicVal = 0;
  double LGMVal = 0;
  LGMSwptns *ExerIntoPtr = NULL;
  LGMSwptns *ExerBdryPtr = NULL;
  LGM_TS *LGMtsPtr = NULL;
  SigKapTS *SigKaptsPtr = NULL;
  Date tNow, tfirst;
  SrtCurvePtr yldcrv = NULL;

  LGMErr error = NULL;
  long skipCal;

  /*	1.- Initialise generic data
          ---------------------------	*/

  /* Free arrays that will be used for output  , if already allocated */
  free_outputs(&convertTS, &findExerBdry, lgmExerBdryData, lgmRefSwptnData,
               lgmTSData, atcTSData);

  /* Set operations to be done by LGMnewautocal */
  skipCal = 0; /* always calibrate */

  /* Set information about today  , and today's market */
  yldcrv = lookup_curve(ycName);
  tNow = get_clcndate_from_yldcrv(yldcrv); /* get calculation date */

  /* Endofday and GetVol available from input */

  /* Get the exponent beta to use to interpret market vols */
  /* Transform char_vol_type to SrtDiffusionType */
  error = interp_diffusion_type(char_vol_type, &srt_vol_type);
  if (error)
    return error;

  if (srt_vol_type == SRT_LOGNORMAL) {
    GetBeta = LGMReturnExp1; /* This will return beta=1 for all swaptions */
  } else {
    GetBeta = LGMReturnExp0; /* This will return beta=0 for all swaptions */
  }

  /* Determine the first possible day that an exercise can take place */
  tfirst = tNow;
  if (endofday == 1) {
    tfirst = add_unit(tNow, 1, SRT_BDAY, SUCCEEDING);
  }

  /* Interpret pay-receive flag and basis */
  error = interp_rec_pay(PayRecStr, &srt_p_r_flag);
  if (error)
    return error;

  /* Set calibration method */
  CalReqPtr = LGMSetCalibMeth(LGMOneTwoFactor, usefixtau, usecaps, tau, alpha,
                              gamma, rho, maxstd, Zeta1Dates, StartZeta1s,
                              TauDates, StartTaus, HybridShortInstrsIndex);

  if (CalReqPtr == NULL) {
    return ("alloc failed in LGMcaller");
  }

  if (LGMOneTwoFactor == 2) {
    if (usefixtau == 0) {
      if ((calibrationmeth == 1) && (HybridShortInstrsIndex != NULL))
        CalReqPtr->calmeth = HybridShortAndDiag;
      else if ((calibrationmeth == 2) || (HybridShortInstrsIndex == NULL))
        CalReqPtr->calmeth = FullCapAndDiag;
      else {
        return serror("Calibration Method should be 1: HybridShortAndDiag or "
                      "2: FullCapAndDiag ");
        /* CalReqPtr->calmeth = TenorAndDiag; */
      }
    }
  }

  /* DEFINE INSTRUMENTS TO CALIBRATE ON */
  if (LGMOneTwoFactor == 1) {
    if (calibrationmeth == 1 && usefixtau != 1)
      /* CALIBRATION ON FIXED MATURITY SWAPTIONS + DIAG SWAPTIONS OR CAPLETS IF
       * usecaps==1 */
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
  if (!EvalParamsPtr) {
    LGMFreeCalParm(&CalReqPtr);
    return "alloc failed in LGMcaller";
  }

  /* Record Vol */
  RecVol = LGMRecVolDummy;

  /*	2.- Initialise specific Midat data
          ---------------------------------- */

  /* Set deal type */
  ItsAMidAt = SimpleMidAt;

  /* allocate MidAtlantic stucture  , and copy input data into it */
  /* This routine ASSUMES that the option holder will receive all payments
  whose pay dates tPay[i] are strictly after the settlement date. This shoud be
  replaced by user input */
  error = LGMFillSimMidAt(tfirst, &MidAtPtr, nPay, Payment, tPay, nEx, tEx,
                          tStart, Strike, RedFirstPay, srt_p_r_flag);

  if (error != NULL) {
    LGMFreeCalParm(&CalReqPtr);
    srt_free(EvalParamsPtr);
    LGMFreeSimMidAt(&MidAtPtr);
    return (error);
  }

  /*	3.- Call to LGMAutocal
          ---------------------- */

  error = LGMnewautocal(
      skipEval, skipCal, convertTS, findExerBdry, /* task list */
      tNow, endofday, ycName,      /* info about today's market */
      GetVol, GetBeta, RecVol,     /* swaption & caplet price info */
      ItsAMidAt, (void *)MidAtPtr, /* the deal */
      CalReqPtr, EvalParamsPtr,    /* calib & eval methods */
      &LGMVal, &intrinsicVal,      /* value & intrinsic value of deal */
      &ExerIntoPtr,
      &ExerBdryPtr, /* swaptions representing the underlying & exer boundary */
      lgmRefSwptnData,         /* reference swaptions */
      &LGMtsPtr, &SigKaptsPtr, /* calibrated term structures */
      szRefRate);              /* unused */

  /*	4.- Free and return
          ------------------- */

  /* free unneeded structures */ /* these are not needed to get output */
  LGMFreeSimMidAt(&MidAtPtr);    /* structure containing the deal */
  LGMFreeCalParm(
      &CalReqPtr);         /* structure containing the requested cal method */
  srt_free(EvalParamsPtr); /* structure containing the convolution parameters */
  LGMFreeSwptns(
      &ExerIntoPtr); /* structure containing underlying swaptions of MidAt */

  /* keep track of Autocal Term Struct if required  , otherwise free it */
  if (LGMOneTwoFactor == 1) {
    if (atcTSData) {
      *atcTSData = LGMtsPtr;
    } else {
      LGMFreeLGM_TS(&LGMtsPtr);
    }
  } else if (LGMOneTwoFactor == 2) {
    if (atcTSData) {
      *atcTSData = LGMtsPtr;
    } else if (LGMtsPtr)
      LGMFreeLGM2F_TS(&LGMtsPtr);
  }

  if (error != NULL) {
    LGMFreeSwptns(&ExerBdryPtr);
    LGMFreeSigKapTS(&SigKaptsPtr);
    if (LGMOneTwoFactor == 1 && atcTSData) {
      *atcTSData = NULL;
      LGMFreeLGM_TS(&LGMtsPtr);
    }
    return (error);
  }

  /* Unpack output */
  *LGMvalPtr = LGMVal;
  if (intrinsicValPtr != NULL)
    *intrinsicValPtr = intrinsicVal;

  /* If exercise boundary has been found  , copy it into tExBdry & rExBdry */
  if (skipEval != 1 && findExerBdry != 0 && ExerBdryPtr != NULL &&
      ExerBdryPtr->n > 0 && lgmExerBdryData != NULL)
    error = copyExerBdry(ExerBdryPtr, lgmExerBdryData);
  if (error != NULL) {
    LGMFreeSwptns(&ExerBdryPtr); /* no longer needed */
    LGMFreeSigKapTS(&SigKaptsPtr);
    if (LGMOneTwoFactor == 1 && atcTSData) {
      *atcTSData = NULL;
      LGMFreeLGM_TS(&LGMtsPtr);
    }
    if (lgmExerBdryData != NULL) {
      lgmExerBdryData->NexerBdry = 0;
      if (lgmExerBdryData->exerBdryArr != NULL)
        srt_free(lgmExerBdryData->exerBdryArr);
    }
    return (error);
  }
  LGMFreeSwptns(&ExerBdryPtr); /* no longer needed */

  /* If sigma-kappa term structure has been found  , copy it into output arrays
   */
  if (skipEval != 1 && convertTS != 0 && SigKaptsPtr != NULL &&
      SigKaptsPtr->numS > 0 && SigKaptsPtr->numK > 0 &&
      SigKaptsPtr->sdate != NULL && SigKaptsPtr->sig != NULL &&
      SigKaptsPtr->kdate != NULL && SigKaptsPtr->kap != NULL &&
      lgmTSData != NULL) {
    error = copyTSData(SigKaptsPtr, lgmTSData);
    if (error != NULL) {
      LGMFreeSigKapTS(&SigKaptsPtr);
      if (LGMOneTwoFactor == 1 && atcTSData) {
        *atcTSData = NULL;
        LGMFreeLGM_TS(&LGMtsPtr);
      }
      return (error);
    }
  }
  LGMFreeSigKapTS(&SigKaptsPtr); /* no longer needed */

  return error;
}

LGMErr LGMCallInvFloaterCaller(

    /*	Exercise	*/

    long nEx,     /* number of exercise dates		*/
    Date *tEx,    /* [0  ,1  ,...  ,nEx-1] notification dates:
                                  first coupon to be called is the first coupon
                     with a start date    on or after the notification date */
    Date *tStart, /* [0  ,1  ,...  ,nEx-1] settlement dates for each exercise
                                          i.e. dates on which the fee is paid */
    double
        *exerFee, /* [0  ,1  ,...  ,nEx-1] fee paid by option holder to exercise
                                          most frequently  , the bond is called
                     at par  , i.e. exerFee = 1.0 */
    char
        *PayRecStr, /* RECEIVER (call on the bond) or PAYER (put on the bond) */

    /*	Coupons: the structure must be decomposed into a swap that delivers:

- on the funding side  , Libor CASH FLAT

- on the exotic side  , a - gear * Libor CASH * libor_cvg

- the exotic side also includes a cap paying max (0  , Libor CASH - cap_str) *
libor_cvg */

    long nCpn,       /* number of coupon periods */
    Date *tCpnStart, /* [0  ,1  ,...  ,nCpn-1] coupon start date */
    Date *tCpnPay,   /* [0  ,1  ,...  ,nCpn-1] coupon pay date */

    double *a,         /* [0  ,1  ,...  ,nCpn-1] */
    double *gear,      /* [0  ,1  ,...  ,nCpn-1] */
    double *cvg,       /* [0  ,1  ,...  ,nCpn-1] */
    double *libor_cvg, /* [0  ,1  ,...  ,nCpn-1] */

    /*	Cap on cash libor	*/

    double *cap_str, /* [0  ,1  ,...  ,nCpn-1] */

    /*	Market & context	*/

    String ycName, /* yield curve name */
    int endofday,  /* 1 = too late to notify today  , 0 = not too late */
    LGMErr (*GetVol)(long, long, double, SRT_Boolean, double *),
    /* Autocal volatility function to get market cash vol */
    char *char_vol_type, /* normal or log normal swaption vols */

    /*	Calibration	*/

    int LGMOneTwoFactor, int usefixtau, /* 1 = calibrate with fixed tau */
    int usecaps,  /* 1 = use caplets for calibration  , 0=use only swaptions */
    double tau,   /* if fixed tau  , use this value for tau (in years) */
    double alpha, /* For the future 2F implementation  , usused so far */
    double gamma, /* For the future 2F implementation  , usused so far */
    double rho,   /* For the future 2F implementation  , usused so far */
    int calibrationmeth, /* Always choose 2 */
    int capletvolmethod, /* 1 = Model  , 2 = Sliding  , 3 = Converging  , 4 =
                            Calibrated */
    int strikechoice,    /* 1 = Not used  , 2 = MidatStrike/MidatStrike  ,
                                            3 = MidatStrike/std  , 4=
                            MidatStrike/ATM */
    double maxstd,       /* Maximum number of std between forward and strike */

    /*	Results	*/

    int convertTS,     /* 1 = cmpute new sigmas  , taus & store in ts; 0 = don't
                          bother */
    int findExerBdry,  /* 1 = find swap rates at exercise boundary; 0 = don't
                          bother */
    double *LGMValPtr, /* answer...value of the option */
    SrtLgmExerBdryData *lgmExerBdryData, /* ptr to exercise boundary data
                                            structure (NULL => not req'd) */
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
    SrtLgmTSData *lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr *atcTSData,    /* ptr to zeta/G data (NULL => not req'd) */
    FWD_VOL_STR *fwdVolStr)  /* ptr to fwd vol data (NULL => not req'd) */
{
  LGMDealType dealtype;
  SrtCurvePtr yldcrv = NULL;
  SrtReceiverType srt_p_r_flag;
  SrtDiffusionType srt_vol_type;
  SrtCallInvFlt *dealPtr;
  LGMCalParm *CalReqPtr = NULL;
  ConvParams *EvalParamsPtr = NULL;
  LGM_TS *LGMtsPtr = NULL;
  SigKapTS *SigKaptsPtr = NULL;
  LGMErr error = NULL;
  Date tNow, tfirst;
  long nExLeft;
  long skipCal;
  double LGMVal, intrinsicVal;
  LGMErr (*GetBeta)(Date, Date, double *);
  LGMErr (*RecVol)(Date, Date, double);
  LGMSwptns *ExerIntoPtr = NULL;
  LGMSwptns *ExerBdryPtr = NULL;
  int i, j;

  /*	1.- Initialise generic data
          ---------------------------	*/

  /* Free arrays that will be used for output  , if already allocated */
  free_outputs(&convertTS, &findExerBdry, lgmExerBdryData, lgmRefSwptnData,
               lgmTSData, atcTSData);

  /* Set operations to be done by LGMnewautocal */
  skipCal = 0; /* always calibrate */

  /* Set information about today  , and today's market */
  yldcrv = lookup_curve(ycName);
  tNow = get_clcndate_from_yldcrv(yldcrv); /* get calculation date */

  /* Endofday and GetVol available from input */

  /* Get the exponent beta to use to interpret market vols */
  /* Transform char_vol_type to SrtDiffusionType */
  error = interp_diffusion_type(char_vol_type, &srt_vol_type);
  if (error)
    return error;

  if (srt_vol_type == SRT_LOGNORMAL) {
    GetBeta = LGMReturnExp1; /* This will return beta=1 for all swaptions */
  } else {
    GetBeta = LGMReturnExp0; /* This will return beta=0 for all swaptions */
  }

  /* Determine the first possible day that an exercise can take place */
  tfirst = tNow;
  if (endofday == 1) {
    tfirst = add_unit(tNow, 1, SRT_BDAY, SUCCEEDING);
  }

  /* Interpret pay-receive flag and basis */
  error = interp_rec_pay(PayRecStr, &srt_p_r_flag);
  if (error)
    return error;

  /* Set calibration method */
  CalReqPtr = LGMSetCalibMeth(LGMOneTwoFactor, usefixtau, usecaps, tau, alpha,
                              gamma, rho, maxstd, NULL, NULL, NULL, NULL, NULL);

  if (CalReqPtr == NULL) {
    return ("alloc failed in LGMcaller");
  }

  if (LGMOneTwoFactor == 2) {
    return "Two factor model for CIF is not available yet";
  }

  /* DEFINE INSTRUMENTS TO CALIBRATE ON */
  if (calibrationmeth == 1 && usefixtau != 1)
    /* CALIBRATION ON FIXED MATURITY SWAPTIONS + DIAG SWAPTIONS OR CAPLETS IF
     * usecaps==1 */
    CalReqPtr->calmeth = FixExp;

  else if (calibrationmeth == 2 && usefixtau != 1)
    /* CALIBRATION ON LONG & SHORT SWAPTIONS  */
    CalReqPtr->calmeth = TenorAndDiag;

  else if (calibrationmeth == 3 && usefixtau != 1)
    CalReqPtr->calmeth = FixSigma;

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
  if (!EvalParamsPtr) {
    LGMFreeCalParm(&CalReqPtr);
    return "alloc failed in LGMcaller";
  }

  /* Record Vol */
  RecVol = LGMRecVolDummy;

  /*	2.- Initialise specific callable inverse floater data
          ----------------------------------------------------- */

  /* Set deal type */
  dealtype = CallInvFlt;

  /* Allocate and fill deal pointer  , also check that there is a least one
   * exerceise date left */
  error = LGMFillCallInvFloater(
      &dealPtr, &nExLeft, tfirst, ycName, nCpn, tCpnStart, tCpnPay, a, gear,
      cap_str, cvg, libor_cvg, nEx, tEx, tStart, exerFee, srt_p_r_flag);

  if (error) {
    LGMFreeCalParm(&CalReqPtr);
    srt_free(EvalParamsPtr);
    LGMFreeCallInvFlt(&dealPtr);
    return error;
  }

  if (nExLeft == 0) {
    *LGMValPtr = 0.0;
    LGMFreeCalParm(&CalReqPtr);
    srt_free(EvalParamsPtr);
    LGMFreeCallInvFlt(&dealPtr);
    return NULL;
  }

  /* DEFINE THE CAPLET VOLATILITY METHOD */
  if (capletvolmethod == 1)
    CalReqPtr->CapletVolMeth = MODEL;
  else if (capletvolmethod == 2)
    CalReqPtr->CapletVolMeth =
        SLIDING2; /* Just change the name SLIDING to SLIDING2 to avoid conflict
                     of definitions */
  else if (capletvolmethod == 3)
    CalReqPtr->CapletVolMeth = CONVERGING;
  else if (capletvolmethod == 4)
    CalReqPtr->CapletVolMeth = MARKET;

  /*	3.- Call to LGMAutocal
          ---------------------- */

  error = LGMnewautocal(
      0, skipCal, convertTS, findExerBdry, /* task list */
      tNow, endofday, ycName,              /* info about today's market */
      GetVol, GetBeta, RecVol,             /* swaption & caplet price info */
      dealtype, (void *)dealPtr,           /* the deal */
      CalReqPtr, EvalParamsPtr,            /* calib & eval methods */
      &LGMVal, &intrinsicVal,              /* value & intrinsic value of deal */
      &ExerIntoPtr,
      &ExerBdryPtr, /* swaptions representing the underlying & exer boundary */
      lgmRefSwptnData,         /* reference swaptions */
      &LGMtsPtr, &SigKaptsPtr, /* calibrated term structures */
      NULL);                   /* unused */

  /*	4.- Free and return
          ------------------- */

  /* Free unneeded structures */ /* these are not needed to get output */
  LGMFreeCalParm(&CalReqPtr);    /* contains the calibration request */
  srt_free(EvalParamsPtr);       /* contains the convolution parameters */
  LGMFreeSwptns(
      &ExerIntoPtr); /* structure containing underlying swaptions of MidAt */
  if (atcTSData) {
    *atcTSData = LGMtsPtr;
  } else {
    LGMFreeLGM_TS(&LGMtsPtr);
  }

  if (error) {
    LGMFreeCallInvFlt(&dealPtr);
    LGMFreeSwptns(&ExerBdryPtr);
    LGMFreeSigKapTS(&SigKaptsPtr);
    if (atcTSData) {
      *atcTSData = NULL;
      LGMFreeLGM_TS(&LGMtsPtr);
    }
    return error;
  }

  /* Unpack output */
  *LGMValPtr = LGMVal;

  /* If exercise boundary has been found  , copy it into tExBdry & rExBdry */
  if (findExerBdry != 0 && ExerBdryPtr != NULL && ExerBdryPtr->n > 0 &&
      lgmExerBdryData != NULL) {
    error = copyExerBdry(ExerBdryPtr, lgmExerBdryData);
  }

  if (error) {
    LGMFreeCallInvFlt(&dealPtr);
    LGMFreeSwptns(&ExerBdryPtr);
    LGMFreeSigKapTS(&SigKaptsPtr);
    if (atcTSData) {
      *atcTSData = NULL;
      LGMFreeLGM_TS(&LGMtsPtr);
    }
    if (lgmExerBdryData) {
      lgmExerBdryData->NexerBdry = 0;
      if (lgmExerBdryData->exerBdryArr) {
        srt_free(lgmExerBdryData->exerBdryArr);
      }
    }
    return error;
  }
  LGMFreeSwptns(&ExerBdryPtr); /* no longer needed */

  /* If sigma-kappa term structure has been found  , copy it into output arrays
   */
  if (convertTS != 0 && SigKaptsPtr != NULL && SigKaptsPtr->numS > 0 &&
      SigKaptsPtr->numK > 0 && SigKaptsPtr->sdate != NULL &&
      SigKaptsPtr->sig != NULL && SigKaptsPtr->kdate != NULL &&
      SigKaptsPtr->kap != NULL && lgmTSData != NULL) {
    error = copyTSData(SigKaptsPtr, lgmTSData);
    if (error) {
      LGMFreeCallInvFlt(&dealPtr);
      LGMFreeSigKapTS(&SigKaptsPtr);
      if (atcTSData) {
        *atcTSData = NULL;
        LGMFreeLGM_TS(&LGMtsPtr);
      }
      return error;
    }
  }
  LGMFreeSigKapTS(&SigKaptsPtr); /* no longer needed */

  /* Copy forward vol */
  if (fwdVolStr) {
    fwdVolStr->fwdVolCpn = dealPtr->fwdVolCpn;
    fwdVolStr->fwdVolEx = dealPtr->fwdVolEx;
    fwdVolStr->fwdVolMat = dmatrix(0, dealPtr->fwdVolCpn, 0, dealPtr->fwdVolEx);
    if (!fwdVolStr->fwdVolMat) {
      LGMFreeCallInvFlt(&dealPtr);
      return "Allocation of forward volatility matrix falied in "
             "LGMCallInvFloaterCaller";
    }

    for (i = 0; i <= dealPtr->fwdVolCpn; i++) {
      for (j = 0; j <= dealPtr->fwdVolEx; j++) {
        fwdVolStr->fwdVolMat[i][j] = dealPtr->fwdVolMat[i][j];
      }
    }
  }

  LGMFreeCallInvFlt(&dealPtr);

  return error;
}

/**********************************************************************************************/
/* Prices a Bermudan option on a cap floater
   Model is on a cash basis  , with all basis spreads being accounted for
before reaching this code. See below
   The underlying has two legs
Coupon leg:
        C[i] = cvg(t[i-1]  ,t[i]) * (margin[i] + lvg*max{rate[i]-amin[i]  ,0} -
lvg*max{rate[i]-amax[i]  ,0}) received at t[i] for i = 1  ,2  ,...  ,n where
rate[i] is the true (ie  , cash) floating rate for the period t[i-1] to t[i]

Funding leg:
        cvg(tau[k-1]  , tau[k]) * rate2[k]	for k = 1  ,2  , ...  , m
                where rate2[k] is the true (ie  , cash) floating rate for the
period tau[k-1] to tau[k]


Exercise (receiver):
        At any notification date tEx[j]  , j=0  ,1  , ...  , nEx-1  , the holder
can exercise. If exercised at tEx[j] then the settlement date is t[i0-1]  , the
next coupon date on or after tEx[j]. The option owner pays any exercise fee f[j]
at t[i0-1] receives all remaining coupons C[i]  , i = i0  , i0+1  , ...  , n
        receives notional ($1) at end date t[n] (this simplifies the funding
leg) pays the funding leg starting at t[i0-1]. If the start-upon-exercise date
t[i0-1] is a funding date tau[k-1]  , then the option owner receives the funding
payments at tau[k]  , tau[k+1]  , ...  , tau[m] receives the notional ($1) at
the end date The funding leg thus revalues to $1

If the start-upon-exercise date t[i0-1] is not a funding date  , so it occurs
somewhere in the middle of a funding period  , then there are small adjustments
that depend on whether the Early Flag and the ResetRate flags are on.

If the Early Flag is on (this represents an option that is cancelling an
existing CapFloater  , the option owner receives the accrued (funding leg)
interest cvg(tau[k-1]  ,t[i0-1]) * rate2[k] at t[i0-1]  , pays the entire first
funding payment cvg(tau[k-1]  ,tau[k]) * rate2[k] at tau[k] continues to make
the funding payments at tau[k+1]  , tau[k+2]  , ...tau[m]

If the start of this funding period occurs before spot-of-today (which happens
in re-valuing deals) then the floating rate for this period has already been
set; in this case the code uses the input rfund_cur (if it is positive); if
rfund_cur<=0 then the code will guess the correct rate from the discount curve.

If the Early Flag is off  , then the option owner
        pays remaining part of the first funding payment cvg(t[i0-1]  ,tau[k]) *
rate2[k] at tau[k] continues to make the funding payments at tau[k+1]  ,
tau[k+2]  , ...tau[m] If the Reset Rate flag is on  , the rate is assumed to
have been reset for the first (stub) period. Otherwise the code will use
rfund_cur (or guess the rate) as before.

Exercise (payer):
        At any notification date tEx[j]  , j=0  ,1  , ...  , nEx-1  , the holder
can exercise. If exercised  , the option owner pays an exercise fee f[j] pays
the coupon leg (reverse of above) receives the funding leg (reverse of above)

Notes:
        All margins and "basis spreads" must be moved from the funding leg to
the coupon leg prior to calling the code. The internal model is on a CASH basis
, with the FRA rates rate for t1 to t2 = [df(t1) - df(t2)]/cvg*df(t2) being used
in place of the floating rates. Before the code is called  , all basis spreads
must be accounted for by
(i) adjusting the funding leg margin (before moving the margin to the coupon
leg) (ii) adjusting the effective floor value a[i] The start dates t[0] and
tau[0] should be the same The end dates t[n] and tau[m] should be the same There
is no provision in the code for handling overlapping start/end/pay dates on
either leg Exercise cannot result in starting at any date other than one of the
coupon dates

Sign of the exercise fee: f[j]>0 means the exerciser pays the fee  , regardless
of whether its a payer or receiver. THIS IS UNLIKE LGMautocal  , where the
floating leg always paid the fee

Code automatically disregards any exercise dates that occur before today */
/* tCpn[i-1]  , tCpn[i]  , and tCpn[i] are the the start  , end  , and pay date
   for period i  , i=1  ,...  ,nCpn amin[0]  , amax[0] margin[0] are undefined:
   amin[*  ,1  ,2  ,...  ,nCpn]  , amax[*  ,1  ,2  ,...  ,nCpn]  , margin[*  ,1
   ,2  ,...  ,nCpn]
        tflt[j-1]  , tflt[j]  , tflt[j] are the start  , end  , and pay date for
   floating period j  , j=1  ,...  ,nflt */

LGMErr LGMCallCapFloaterCaller(
    /* definition of coupon leg */
    long nCpn,    /* number of coupon periods */
    Date *tCpn,   /* t[0  ,1  ,...  ,nCpn]...cpn for t[i-1] to t[i] paid at t[i]
                     , i=1  ,...  ,nCpn	*/
    double *amax, /* cpn[i] is
                   */
    double *amin, /*	cvg(t[i-1]  ,t[i]) * (margin[i] +
                   */
    double *margin, /*		lvg*max{rate[i]-amin[i]  ,0} - lvg*max{rate[i]-amax[i]
                       ,0})		*/
    double lvg,     /*			 paid at t[i] for i=1  ,2  , ...  , n
                     */
    char *cpnRateBasisName, /* day count basis for rate[i] */
    char *cpnBasisName,     /* day count basis for coupon leg */
                            /* definition of floating leg */
    long nflt,              /* number of funding periods */
    Date *tflt, /* tau[0  ,1  ,...  ,nflt]...float periofs are tau[j-1] to
                   tau[j]  , j=1  ,...  ,nflt */
    char *fltBasisName, /* day count basis for the floating leg */
                        /* exercise information */
    long nEx,           /* number of exercise dates */
    Date *tEx, /* [0  ,1  ,...  ,nEx-1] notification dates; deal starts on next
                  coupon date */
    double *exerFee,  /* [0  ,1  ,...  ,nEx-1] fee paid by option holder to
                         exercise */
    char *PayRecStr,  /* RECEIVER or PAYER */
    int earlyFlag,    /* if 1  , deal is a cancelation of an existing inverse
                         floater */
    int resetFlt,     /* float for stub period is reset upon exercise */
                      /* information about today and today's market place*/
    String ycName,    /* yield curve name */
    double rfund_cur, /* current (true) funding rate (if fixed); ignored if
                         non-positive */
    int endofday,     /* 1=too late to notify today  , 0=not too late */
    LGMErr (*GetVol)(long, long, double, SRT_Boolean,
                     double *), /* volatility function for reference swptns */
    char *char_vol_type,        /* normal or log normal swaption vols */
                                /* calibration method to use */
    int LGMOneTwoFactor, int usefixtau, /* 1=calibrate with fixed tau */
    int usecaps, /* 1=use caplets for calibration  , 0=use only swaptions */
    double tau,  /* if fixed tau  , use this value for tau (in years) */
    double alpha, double gamma, double rho,
    int calibrationmeth, /* 1 = fixexp  , 2=backboot  , 3=fixed sigma */
    int strikechoice,    /* 1 = IRR  , 2 = dIRR */
                         /* task list */
    int skipEval,        /* 1=calibrate only  , 0=calibrate & evaluate deal */
    int convertTS,       /* 1=compute new sigs and taus; 0=don't bother */
    int findExerBdry,  /* 1=find swap rates at exercise boundary; 0=don't bother
                        */
                       /* outputs */
    String outfile,    /* output log file (unused) */
    double *LGMValPtr, /* answer...value of the option */
    double *intrinValPtr, /* value if exercise date tNot[j] had to be chosen
                             today */
    SrtLgmExerBdryData *lgmExerBdryData, /* ptr to exercise boundary data
                                            structure (NULL => not req'd) */
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
                                         /* calibrated term structure */
    SrtLgmTSData *lgmTSData) /* ptr to tau/sigma data (NULL => not req'd) */
{
  LGMDealType dealtype;
  SrtCurvePtr yldcrv = NULL;
  SrtReceiverType srt_p_r_flag;
  SrtDiffusionType srt_vol_type;
  SrtBasisCode cpnBasis, aBasis, fltBasis;
  SrtCallCapFlt *dealPtr;
  LGMCalParm *CalReqPtr = NULL;
  ConvParams *EvalParamsPtr = NULL;
  LGMErr (*GetBeta)(Date, Date, double *);
  LGMErr (*RecVol)(Date, Date, double);
  LGMSwptns *ExerIntoPtr = NULL;
  LGMSwptns *ExerBdryPtr = NULL;
  LGM_TS *LGMtsPtr = NULL;
  SigKapTS *SigKaptsPtr = NULL;
  LGMErr error = NULL;
  Date tNow, tfirst;
  long nExLeft, nExEff;
  long skipCal;
  double LGMVal, intrinsicVal;

  /* free arrays that will be used for output  , if already allocated */
  if (findExerBdry != 0 && lgmExerBdryData != NULL) {
    lgmExerBdryData->NexerBdry = 0;
    if (lgmExerBdryData->exerBdryArr != NULL)
      srt_free(lgmExerBdryData->exerBdryArr);
  } else
    findExerBdry = 0;

  if (convertTS != 0 && lgmTSData != NULL) {
    lgmTSData->NTS = 0;
    if (lgmTSData->TSArr != NULL)
      srt_free(lgmTSData->TSArr);
  } else
    convertTS = 0;

  /* Set operations to be done by LGMnewautocal */
  *LGMValPtr = 0.0;
  skipCal = 0; /* always calibrate */

  /* Set information about today  , and today's market */
  yldcrv = lookup_curve(ycName);
  tNow = get_clcndate_from_yldcrv(yldcrv); /* get calculation date */

  /* endofday and GetVol available from input */

  /* get the exponent beta to use to interpret market vols */
  /* transform char_vol_type to SrtDiffusionType */
  error = interp_diffusion_type(char_vol_type, &srt_vol_type);
  if (error)
    return (error);

  if (srt_vol_type == SRT_LOGNORMAL)
    GetBeta = LGMReturnExp1; /* This will return beta=1 for all swaptions */
  else
    GetBeta = LGMReturnExp0; /* This will return beta=0 for all swaptions */

  /* determine the first possible day that an exercise can take place */
  tfirst = tNow;
  if (endofday == 1)
    tfirst = add_unit(tNow, 1, SRT_BDAY, SUCCEEDING);

  /* interpret pay-receive flag and basis */
  error = interp_rec_pay(PayRecStr, &srt_p_r_flag);
  if (error)
    return (error);
  error = interp_basis(cpnBasisName, &cpnBasis);
  if (error)
    return (error);
  error = interp_basis(cpnRateBasisName, &aBasis);
  if (error)
    return (error);
  error = interp_basis(fltBasisName, &fltBasis);
  if (error)
    return (error);

  /* Set deal type */
  dealtype = CallCapFlt;

  /* the next routine checks out the deal  , allocates the Bermudan on a inverse
  floater stucture  , and copies the input data into it */
  error = LGMFillCallCapFloater(&nExLeft, tfirst, ycName, rfund_cur, &dealPtr,
                                nCpn, tCpn, amax, amin, margin, lvg, cpnBasis,
                                aBasis, nflt, tflt, fltBasis, nEx, tEx, exerFee,
                                srt_p_r_flag, earlyFlag, resetFlt);

  if (error) {
    LGMFreeCallCapFlt(&dealPtr);
    return (error);
  }
  if (nExLeft == 0) /* no exercise dates left */
  {
    *LGMValPtr = 0.0;    /* deal is worth zero */
    *intrinValPtr = 0.0; /* free and return */
    LGMFreeCallCapFlt(&dealPtr);
    return (NULL);
  }
  error = LGMValidCallCapFlt(dealPtr, tfirst, &nExEff);
  if (error) {
    LGMFreeCallCapFlt(&dealPtr);
    return (error);
  }

  /* Set calibration method */
  CalReqPtr = LGMSetCalibMeth(LGMOneTwoFactor, usefixtau, usecaps, tau, alpha,
                              gamma, rho, 2.0, NULL, NULL, NULL, NULL, NULL);
  if (CalReqPtr == NULL) {
    LGMFreeCallCapFlt(&dealPtr);
    return ("alloc failed in LGMcaller");
  }
  if (calibrationmeth == 1 && usefixtau != 1)
    CalReqPtr->calmeth = FixExp;
  else if (calibrationmeth == 2 && usefixtau != 1)
    CalReqPtr->calmeth = TenorAndDiag;
  else if (calibrationmeth == 3 && usefixtau != 1)
    CalReqPtr->calmeth = FixSigma;
  if (strikechoice == 1)
    CalReqPtr->Rmeth = IRR;
  else
    CalReqPtr->Rmeth = dIRR;

  /* Set numerical parameters for convolution to their default values */
  EvalParamsPtr = LGMSetDefaultEvalParms();
  if (EvalParamsPtr == NULL) {
    LGMFreeCallCapFlt(&dealPtr);
    LGMFreeCalParm(&CalReqPtr);
    return ("alloc failed in LGMcaller");
  }

  /* Record Vol */
  RecVol = LGMRecVolDummy;

  /* Call LGMautocal */
  error = LGMnewautocal(
      skipEval, skipCal, convertTS, findExerBdry, /* task list */
      tNow, endofday, ycName,    /* info about today's market */
      GetVol, GetBeta, RecVol,   /* swaption & caplet price info */
      dealtype, (void *)dealPtr, /* the deal */
      CalReqPtr, EvalParamsPtr,  /* calib & eval methods */
      &LGMVal, &intrinsicVal,    /* value & intrinsic value of deal */
      &ExerIntoPtr,
      &ExerBdryPtr, /* swaptions representing the underlying & exer boundary */
      lgmRefSwptnData,         /* reference swaptions */
      &LGMtsPtr, &SigKaptsPtr, /* calibrated term structures */
      outfile);                /* unused */

  /* free unneeded structures */ /* these are not needed to get output */
  LGMFreeCallCapFlt(&dealPtr);   /* contains the deal */
  LGMFreeCalParm(&CalReqPtr);    /* contains the calibration request */
  srt_free(EvalParamsPtr);       /* contains the convolution parameters */
  LGMFreeLGM_TS(&LGMtsPtr);      /* contains the calibrated zeta_G term struc */

  if (error) {
    LGMFreeSwptns(&ExerBdryPtr);
    LGMFreeSigKapTS(&SigKaptsPtr);
    return (error);
  }

  /* Unpack output */
  *LGMValPtr = LGMVal;
  if (intrinValPtr != NULL)
    *intrinValPtr = intrinsicVal;

  /* If exercise boundary has been found  , copy it into tExBdry & rExBdry */
  if (skipEval != 1 && findExerBdry != 0 && ExerBdryPtr != NULL &&
      ExerBdryPtr->n > 0 && lgmExerBdryData != NULL)
    error = copyExerBdry(ExerBdryPtr, lgmExerBdryData);
  if (error) {
    LGMFreeSwptns(&ExerBdryPtr); /* no longer needed */
    LGMFreeSigKapTS(&SigKaptsPtr);
    if (lgmExerBdryData != NULL) {
      lgmExerBdryData->NexerBdry = 0;
      if (lgmExerBdryData->exerBdryArr != NULL)
        srt_free(lgmExerBdryData->exerBdryArr);
    }
    return (error);
  }
  LGMFreeSwptns(&ExerBdryPtr); /* no longer needed */

  /* If sigma-kappa term structure has been found  , copy it into output arrays
   */
  if (skipEval != 1 && convertTS != 0 && SigKaptsPtr != NULL &&
      SigKaptsPtr->numS > 0 && SigKaptsPtr->numK > 0 &&
      SigKaptsPtr->sdate != NULL && SigKaptsPtr->sig != NULL &&
      SigKaptsPtr->kdate != NULL && SigKaptsPtr->kap != NULL &&
      lgmTSData != NULL) {
    error = copyTSData(SigKaptsPtr, lgmTSData);
    if (error != NULL) {
      LGMFreeSigKapTS(&SigKaptsPtr);
      return (error);
    }
  }
  LGMFreeSigKapTS(&SigKaptsPtr); /* no longer needed */

  return (error);
}

/***************************************************************************************/
/* Routines that help make the conversion from old LGM call to new LGM call */
/* This routine creates a simple American deal and copies the input data into
the simple American structure */
LGMErr LGMFillSimAmer(Date tfirst, int standardLag, SrtSimAmer **AmerPtrPtr,
                      long nfix, Date *tfixStart, Date *tfixEnd, Date *tfixPay,
                      double *fixFullPayment, SrtBasisCode fixBasisCode,
                      int fixEarlyFlag, double *Strike, long nfltdates,
                      Date *tflt, SrtBasisCode fltBasisCode, int fltEarlyFlag,
                      Date tFirstExer, int lagExerStart, int calbus,
                      BusDayConv convStart, SrtReceiverType payrec) {
  long i, j, nflt;
  Date tfirstStart, LastFixing;
  SrtSimAmer *AmerPtr = NULL;

  /* check out deal */
  nflt = nfltdates - 1; /* number of floating periods */
  if (nfix < 1 || nfix > 601 || nflt < 1 || nflt > 601)
    return ("too few or too many dates");

  tFirstExer = max(tfirst, tFirstExer);
  LastFixing =
      add_unit(tflt[nflt - 1], -standardLag, SRT_BDAY, NO_BUSDAY_CONVENTION);
  if (tFirstExer >= LastFixing)
    return ("deal is a European swaption  , a swap  , or worthless");

  if (calbus == 1)
    tfirstStart = add_unit(tFirstExer, lagExerStart, SRT_BDAY, convStart);
  else
    tfirstStart = add_unit(tFirstExer, lagExerStart, SRT_DAY, convStart);
  if (tfirstStart >= tfixEnd[nfix - 1] || tfirstStart >= tflt[nfltdates - 1])
    return ("no underlying");

  for (i = 0; i < nfix - 1 && tfixStart[i] < tfixEnd[i] &&
              tfixStart[i] < tfixStart[i + 1] && tfixEnd[i] < tfixEnd[i + 1] &&
              tfixPay[i] <= tfixPay[i + 1] && tfixStart[i] < tfixPay[i];
       i++)
    ;
  if (tfixStart[nfix - 1] >= tfixEnd[nfix - 1] ||
      tfixStart[nfix - 1] >= tfixPay[nfix - 1])
    i = nfix - 2;

  for (j = 0; j < nfltdates - 1 && tflt[j] < tflt[j + 1]; j++)
    ;

  if (tflt[nfltdates - 1] <= tfirstStart)
    j = -1;

  if (i < (nfix - 1) || j < (nfltdates - 1))
    return ("dates are wrong");

  /* Deal looks OK ... create simple American deal structure */
  if (*AmerPtrPtr != NULL) /* free deal structure if allocated */
    LGMFreeSimAmer(AmerPtrPtr);
  AmerPtr = LGMCreateSimAmer(nfix, nflt);
  if (AmerPtr == NULL)
    return ("allocation of deal structure failed in LGMautocal");
  *AmerPtrPtr = AmerPtr; /* output is AmerPtrPtr */

  /* Fill in structure */
  /* fill in fixed leg */
  AmerPtr->nfix = nfix;
  AmerPtr->fixBasis = fixBasisCode;
  AmerPtr->EarlyFlagFix = fixEarlyFlag;

  for (i = 0; i < nfix; i++) {
    AmerPtr->tfixStart[i] = tfixStart[i];
    AmerPtr->tfixEnd[i] = tfixEnd[i];
    AmerPtr->tfixPay[i] = tfixPay[i];
    AmerPtr->fixCoupon[i] = fixFullPayment[i];
    AmerPtr->ExtraPrem[i] =
        Strike[i] - 1; /* remove notional from floating leg */
  }
  AmerPtr->fixCoupon[nfix - 1] =
      AmerPtr->fixCoupon[nfix - 1] - 1; /* remove notional from fixed leg*/

  /* fill in floating leg */
  AmerPtr->nflt = nflt;
  AmerPtr->fltBasis = fltBasisCode;
  AmerPtr->EarlyFlagFlt = fltEarlyFlag;
  AmerPtr->ResetFlt =
      0; /* do not reset floating rate for mid period exercise */

  for (i = 0; i < nflt; i++) {
    AmerPtr->tfltFixing[i] =
        add_unit(tflt[i], -standardLag, SRT_BDAY, NO_BUSDAY_CONVENTION);
    AmerPtr->tfltPay[i] = tflt[i + 1];
  }

  /* fill in exercise information */
  AmerPtr->PayRec = payrec;
  AmerPtr->tFirstExer = tFirstExer;
  AmerPtr->lagExerSettle = lagExerStart;
  AmerPtr->CalBusLag = calbus;
  AmerPtr->convSettle = convStart;

  return (NULL);
}

/*********************************************************************/
/* This routine creates a simple midatlantic deal and copies the input data into
the simple MidAtlantic structure */
/* This routine ASSUMES that if the option holder exercises on exer date tEx[j]
, then he will receive fixed leg payments on all pay dates tPay[i] which are
strictly after the settlement date Start[j]. Except for the first pay date after
settlement  , the amount of the payment is Payment[i] for all i. The first
payment is reduced by RedFirstPay[j]. */

LGMErr LGMFillSimMidAt(Date tfirst, SrtSimMidAt **MidAtPtrPtr, long nPay,
                       double *Payment, Date *tPay, long nEx, Date *tEx,
                       Date *tStart, double *Strike, double *RedFirstPay,
                       SrtReceiverType payrec) {
  long i, j;
  SrtSimMidAt *MidAtPtr = NULL;

  /* Check out deal */
  if (nPay < 1 || nEx < 1)
    return ("no pay date or no exer date");

  if (nEx > 600 || nPay > 600)
    return ("too many exercise or pay dates");

  for (i = 1; i < nPay; i++) {
    if (tPay[i] < tPay[i - 1])
      return ("Pay dates out of order");
  }

  for (j = 1; j < nEx; j++) {
    if (tEx[j] <= tEx[j - 1] || tStart[j] <= tStart[j - 1])
      return ("exer or settle dates out of order");

    if (tEx[j] > tStart[j] || tStart[j] >= tPay[nPay - 1])
      return ("settlement date too early or late");

    if (Strike[j] < 0)
      return ("strike is negative");
  }

  if (tEx[0] > tStart[0] || tStart[0] >= tPay[nPay - 1])
    return ("settlement date too early or late");

  if (Strike[0] < 0)
    return ("strike is negative");

  /* Deal looks OK ... create simple MidAtlantic deal structure */
  MidAtPtr = *MidAtPtrPtr; /* make code more transparent */
  if (MidAtPtr != NULL)    /* free deal structure if allocated */
    LGMFreeSimMidAt(&MidAtPtr);
  MidAtPtr = LGMCreateSimMidAt(nEx, nPay);
  if (MidAtPtr == NULL)
    return ("allocation of deal structure failed in LGMautocal");
  *MidAtPtrPtr = MidAtPtr; /* output is MidAtPtrPtr */

  /* Fill in structure */
  MidAtPtr->nPay = nPay; /* fill in fixed leg pay dates and payments */
  for (i = 0; i < nPay; i++) {
    MidAtPtr->Payment[i] = Payment[i];
    MidAtPtr->tPay[i] = tPay[i];
  }
  MidAtPtr->PayRec = payrec; /* fill in pay or receive */
  MidAtPtr->nEx = nEx;
  for (j = 0; j < nEx; j++) { /* fill in exercise information */
    MidAtPtr->tEx[j] = tEx[j];
    MidAtPtr->tStart[j] = tStart[j];
    MidAtPtr->Strike[j] = Strike[j];
    MidAtPtr->RedFirstPay[j] = RedFirstPay[j];
  }

  /* for each exercise  , find first payment received if option is exercised */
  for (j = 0; j < nEx; j++) {
    i = 0;
    while (i < nPay && tPay[i] <= tStart[j])
      i++;
    MidAtPtr->FirstPay[j] = i;
  }

  /* find first exercise date on or after tfirst */
  if (tEx[nEx - 1] < tfirst)
    MidAtPtr->FirstExer = nEx; /* No exercise dates left */
  else {
    for (j = 0; j < nEx && tEx[j] < tfirst;
         j++) /* first exer date after today+endofday */
      ;
    MidAtPtr->FirstExer = j;
  }

  return (NULL);
}

/* Create a Bermudan inverse floater structure  , and copy input data into it */
LGMErr
LGMFillCallInvFloater(SrtCallInvFlt **dealPtrPtr, /* output: the deal */
                      long *nExleftPtr, /* output: effective number of exer */
                      Date tfirst, char *ycName,        /* info about today */
                      long n, Date *tStart, Date *tPay, /* coupon leg */
                      double *a, double *gear, double *cap_str, double *cvg,
                      double *lcvg, /* more coupon leg */
                      long nEx, Date *tEx, Date *tExStart,
                      double *exFee,          /* exercise info */
                      SrtReceiverType payrec) /* PAYER or RECEIVER */
{
  LGMErr error = NULL;
  LGMMarkConv conventions; /* standard swaption conventions */
  Date tNow, tSpot;
  long i, j, iFirst;
  long jEff, nExEff, iPrev;
  SrtCallInvFlt *ptr = NULL;
  SrtCurvePtr yldcrv = NULL;

  *dealPtrPtr = NULL;

  /* get yield curve */
  yldcrv = lookup_curve(ycName);
  tNow = get_clcndate_from_yldcrv(yldcrv);      /* get calculation date */
  error = LGMCcyDefaults(yldcrv, &conventions); /* get currency conventions */
  if (error != NULL)
    return (error);
  tSpot = add_unit(tNow, conventions.lag, SRT_BDAY, SUCCEEDING);

  /* Step 1: Check out deal */

  if (n < 1 || nEx < 1) {
    return ("no pay date or no exer date");
  }

  if (n > 600 || nEx > 600) {
    return ("too many exercise or pay dates");
  }

  for (i = 1; i < n; i++) {
    if (tStart[i] < tStart[i - 1] || tPay[i] < tPay[i - 1]) {
      return ("coupon dates out of order");
    }
  }

  for (j = 1; j < nEx; j++) {
    if (tEx[j] <= tEx[j - 1] || tExStart[j] <= tExStart[j - 1]) {
      return ("exer or settle dates out of order");
    }

    if (tEx[j] > tExStart[j] || tEx[j] > tStart[n - 1]) {
      return ("settlement date too late");
    }

    if (exFee[j] < 0) {
      return ("strike is negative");
    }
  }

  if (tEx[0] > tExStart[0] || tExStart[0] >= tPay[n - 1]) {
    return ("settlement date too early or late");
  }

  if (exFee[0] < 0) {
    return ("strike is negative");
  }

  /* Step 2: Determine the relevant exercise and coupon dates */
  /* Figure out the effective number of exercise dates */

  nExEff = 0;
  iPrev = n;
  i = n - 1;
  for (j = nEx - 1; j >= 0 && tEx[j] >= tfirst; j--) {
    while (i >= 0 && tStart[i] >= tEx[j])
      i--;
    i++; /* i is first cpn date on or after tEx[j] */

    if (i < iPrev) {
      iPrev = i; /* exercise has value */
      nExEff++;
    }
  }
  iFirst = i;

  *nExleftPtr = nExEff;
  if (nExEff == 0) {
    return NULL; /* no exercise dates left; deal is worth zero */
  }

  /* Step 3: Allocate Bermudan inverse floater with nExEff exer dates and n-i
   * cpn periods */
  ptr = LGMCreateCallInvFlt(nExEff, n - iFirst);

  /* Step 4: Fill the Bermudan inverse floater structure with the input data */
  for (i = iFirst; i < n; i++) {
    ptr->tCpnStart[i - iFirst] = tStart[i];
    ptr->tCpnPay[i - iFirst] = tPay[i];
    ptr->a[i - iFirst] = a[i];
    ptr->gear[i - iFirst] = gear[i];
    ptr->cvg[i - iFirst] = cvg[i];
    ptr->lcvg[i - iFirst] = lcvg[i];
    ptr->cap_str[i - iFirst] = cap_str[i];
  }

  /* fill in the relevent exercise dates */
  iPrev = n;
  i = n - 1;
  jEff = nExEff;

  for (j = nEx - 1; j >= 0 && tEx[j] >= tfirst; j--) {
    while (i >= 0 && tStart[i] >= tEx[j])
      i--;
    i++; /* i is first cpn date on or after tEx[j] */
    if (i < iPrev) {
      iPrev = i;
      jEff--;
      ptr->tEx[jEff] = tEx[j];
      ptr->tSet[jEff] = tExStart[j];
      ptr->iSet[jEff] = i - iFirst;
      ptr->strike[jEff] = exFee[j];
    }
  }
  if (jEff != 0) {
    LGMFreeCallInvFlt(&ptr);
    return "Serious error";
  }

  ptr->PayRec = payrec;

  *dealPtrPtr = ptr;
  return NULL;
}

/***************************************************************************************/
/* Routines that help make the conversion from old LGM call to new LGM call */
/* This routine verifies the input data  , creates a Bermudan cap floater
structure and copies the input data into the simple cap structure */
LGMErr LGMFillCallCapFloater(
    long *nExleftPtr,                            /* effective number of exer */
    Date tfirst, char *ycName, double rfund_cur, /* info about today */
    SrtCallCapFlt **dealPtrPtr,                  /* output: the deal */
    long n, Date *t, double *amax, double *amin, double *marg, /* coupon leg */
    double lvg, SrtBasisCode cpnBasis,
    SrtBasisCode aBasis,                       /* more coupon leg */
    long m, Date *tflt, SrtBasisCode fltBasis, /* floating leg */
    long nEx, Date *tEx, double *exFee,        /* exercise info */
    SrtReceiverType payrec,                    /* PAYER or RECEIVER */
    int earlyFlag, int resetFlt)               /* early and reset flags */
{
  LGMErr error = NULL;
  LGMMarkConv conventions; /* standard swaption conventions */
  Date tNow, tSpot;
  long i, j, k, iFirst;
  long jEff, nExEff, iPrev;
  double isign, df1, df2, dfSet, cvg, parcvg, rflt, acc;
  SrtCallCapFlt *ptr = NULL;
  SrtCurvePtr yldcrv = NULL;

  *dealPtrPtr = NULL;

  /* get yield curve */
  yldcrv = lookup_curve(ycName);
  tNow = get_clcndate_from_yldcrv(yldcrv);      /* get calculation date */
  error = LGMCcyDefaults(yldcrv, &conventions); /* get currency conventions */
  if (error != NULL)
    return (error);
  tSpot = add_unit(tNow, conventions.lag, SRT_BDAY, SUCCEEDING);

  /* Step 1: Ensure the deal makes sense */
  if (cpnBasis < 0 || cpnBasis >= LASTBASISCODE || aBasis < 0 ||
      aBasis >= LASTBASISCODE || fltBasis < 0 || fltBasis >= LASTBASISCODE)
    return ("unknown basis");

  if (n < 1 || n > 600 || m < 1 || m > 600 || nEx < 1 || nEx > 600)
    return ("too few or too many dates");

  /* Coupon leg dates in increasing order and coupons between -100% and 100% */
  for (i = 1;
       i <= n &&
       t[i - 1] < t[i]
       /*					&& amax[i]>0.0 && amax[i]<1.0
                                               && amin[i]>0.0 && amin[i]<1.0
       */
       && marg[i] > -1.0 && marg[i] < 1.0;
       i++)
    ;
  /* Floating leg dates in increasing order */
  for (k = 1; k <= m && tflt[k - 1] < tflt[k]; k++)
    ;
  /* Exercise dates in increasing order and fees between -100% and 100% */
  for (j = 1;
       j < nEx && tEx[j - 1] < tEx[j] && exFee[j] > -1.0 && exFee[j] < 1.0; j++)
    ;
  /* Check there was no error before and gearing between 0 and 100 */
  if (i <= n || k <= m || j < nEx || exFee[0] <= -1.0 || exFee[0] >= 1.0 ||
      lvg <= 0.0 || lvg >= 100.0)
    return ("deal is mis-entered");

  if (t[n] != tflt[m])
    return ("end dates on legs don't match");

  /* Step 2: Determine the relevant exercise and coupon dates */
  /* Figure out the effective number of exercise dates */
  nExEff = 0;
  iPrev = n;
  i = n;
  for (j = nEx - 1; j >= 0 && tEx[j] >= tfirst; j--) {
    while (i >= 0 && t[i] >= tEx[j])
      i--;
    i++; /* i is first cpn date on or after tEx[j] */
    if (i < iPrev) {
      iPrev = i; /* exercise has value */
      nExEff++;
    }
  }
  *nExleftPtr = nExEff;
  if (nExEff == 0)
    return (NULL); /* no exercise dates left; deal is worth zero */

  /* there are nExEff exercise dates with value */
  /* the relevant coupon dates are all coupon dates on and after t[i] */
  iFirst = i;
  if (tflt[0] > t[iFirst])
    return ("funding dates don't start early enough");

  /* Step 3: Create Bermudan inverse floater with nExEff exer dates and n-i cpn
   * periods */
  ptr = LGMCreateCallCapFlt(nExEff, n - iFirst);

  /* Step 4: Fill the Bermudan inverse floater structure with the input data */
  for (i = iFirst; i <= n; i++) {
    ptr->tCpn[i - iFirst] = t[i];
    ptr->amax[i - iFirst] = amax[i];
    ptr->amin[i - iFirst] = amin[i];
    ptr->marg[i - iFirst] = marg[i];
  }
  ptr->amax[0] = 0.0; /* these values are irrelevent */
  ptr->amin[0] = 0.0; /* these values are irrelevent */
  ptr->marg[0] = 0.0; /* these values are irrelevent */

  ptr->nCpn = n - iFirst;
  ptr->lvg = lvg;
  ptr->cpnBasis = cpnBasis;
  ptr->aBasis = aBasis;

  /* fill in the relevent exercise dates */
  isign = 1.0;
  if (payrec == SRT_PAYER)
    isign = -1.0;

  iPrev = n;
  i = n;
  jEff = nExEff;
  for (j = nEx - 1; j >= 0 && tEx[j] >= tfirst; j--) {
    while (i >= 0 && t[i] >= tEx[j])
      i--;
    i++; /* i is first cpn date on or after tEx[j] */
    if (i < iPrev) {
      iPrev = i; /* exercise has value */
      jEff--;
      ptr->tEx[jEff] = tEx[j];
      ptr->iSet[jEff] = i - iFirst;
      ptr->strike[jEff] = isign * exFee[j];
    }
  }
  if (jEff != 0) {
    LGMFreeCallCapFlt(&ptr);
    return ("serous error");
  }

  ptr->nEx = nExEff;
  ptr->PayRec = payrec;

  /* Step 5: Add funding leg to strike  , accounting for early and reset flags
   */
  for (j = 0; j < ptr->nEx; j++) /* for each exercise j */
  {
    i = ptr->iSet[j]; /* find the floating period */
    for (k = 0; k < m && tflt[k + 1] <= ptr->tCpn[i]; k++)
      ; /* tflt[k]<=tSettle<=tflt[k+1] */

    /* get rflt for period */
    if (tflt[k] >= tSpot)
      df1 = swp_f_df(tNow, tflt[k], ycName);
    else
      df1 = 1.0;
    df2 = swp_f_df(tNow, tflt[k + 1], ycName);
    dfSet = swp_f_df(tNow, ptr->tCpn[i], ycName);
    cvg = coverage(tflt[k], tflt[k + 1], fltBasis);
    if (df1 == SRT_DF_ERROR || df2 == SRT_DF_ERROR || dfSet == SRT_DF_ERROR ||
        cvg <= 0.0) {
      LGMFreeCallCapFlt(&ptr);
      return ("can't get a discount factor");
    }

    /* if rate has been fixed  , use rfund_cur (if positive) */
    if (tflt[k] <= tSpot && rfund_cur > 0.0)
      rflt = rfund_cur;

    /* if rate has been fixed but not entered  , estimate from yield curve */
    else if (tflt[k] < tSpot) {
      parcvg = coverage(tNow, tflt[k + 1], fltBasis);
      if (parcvg <= 0.0) {
        LGMFreeCallCapFlt(&ptr);
        return ("can't get floating rate");
      }
      rflt = (1.0 - df2) / (df2 * parcvg);
    }

    /* if rate has not been fixed  , estimate rate from yield curve */
    else
      rflt = (df1 - df2) / (df2 * cvg);

    /* compute strike - equivalent to the PV of the funding leg at the
     * settlement date */
    if (earlyFlag != 0) {
      acc = coverage(tflt[k], ptr->tCpn[i], fltBasis) * rflt;
      ptr->strike[j] = ptr->strike[j] + (1.0 + cvg * rflt) * df2 / dfSet - acc;
    } else if (resetFlt != 0) {
      parcvg = coverage(ptr->tCpn[i], tflt[k + 1], fltBasis);
      ptr->strike[j] = ptr->strike[j] + (1.0 + parcvg * rflt) * df2 / dfSet;
    } else
      ptr->strike[j] = ptr->strike[j] + 1.0;
  }
  /* check it out */
  error = LGMValidCallCapFlt(ptr, tfirst, nExleftPtr);

  *dealPtrPtr = ptr;
  return (NULL);
}

/*********************************************************************/
/* Code to provide beta for lognormal (beta=1) or normal (beta=0) */
LGMErr LGMReturnExp1(Date swapstart, Date swapend, double *beta) {
  *beta = 1.0;
  return (NULL);
}

/*********************************************************************/
LGMErr LGMReturnExp0(Date swapstart, Date swapend, double *beta) {
  *beta = 0.0;
  return (NULL);
}

/*********************************************************************/
/* This is a dummy routine. It is called by LGMautocal whenever LGMautocal uses
a swaption or caplet volatility in pricing. The vol used is the vol for a
vanilla option where the underlying has start date swapstart  , end date swapend
, and fixed rate Rfix. It should be replaced by a routine which records the
volatilities used by the options. Note that LGM calls the vols for reasons not
directly related to pricing  , meaning these vols will have zero vega. Using
sticky vols to record the vols will thus be more inefficient than using these
vols.
*/
LGMErr LGMRecVolDummy(Date swapstart, Date swapend, double Rfix) {
  return (NULL);
}

/*********************************************************************/
/* This code sets the calibration methodology to the default method */
LGMCalParmPtr LGMSetCalibMeth(int iLGMOneTwoFactor, int FixedTauFlag,
                              int CapsFlag, double Tau, double Alpha,
                              double Gamma, double Rho, double MaxStd,
                              long *Zeta1Dates, double *StartZeta1s,
                              long *TauDates, double *StartTaus,
                              double **HybridShortInstrsIndex) {
  LGMCalParm *CalReqPtr = NULL;

  /* CREATE A SET OF DEFAULT CALIBRATION PARAMETERS */

  CalReqPtr = LGMCreateCalParm(0, 0, 0, 0);
  if (CalReqPtr == NULL)
    return (NULL);

  CalReqPtr->LGMOneTwoFactor = iLGMOneTwoFactor;

  if (iLGMOneTwoFactor == 1) {
    if (FixedTauFlag != 1) /*FIXED TAU = NO (ONLY IN ONE FACT)*/
      CalReqPtr->calmeth = FixSigma;
    else if (FixedTauFlag == 1) {
      CalReqPtr->calmeth = FixKappa;
      CalReqPtr->kap = 1.0 / Tau;
    }

    CalReqPtr->usecaps = CapsFlag;
  }
  if (iLGMOneTwoFactor == 2) {
    CalReqPtr->usecaps = CapsFlag;
    CalReqPtr->alpha = Alpha;
    CalReqPtr->gamma = Gamma;
    CalReqPtr->rho = Rho;
    CalReqPtr->Zeta1Dates = Zeta1Dates;
    CalReqPtr->StartZeta1s = StartZeta1s;
    CalReqPtr->TauDates = TauDates;
    CalReqPtr->StartTaus = StartTaus;
    CalReqPtr->HybridShortInstrsIndex = HybridShortInstrsIndex;

    if (FixedTauFlag == 1) {
      CalReqPtr->calmeth = FixKappa;
      CalReqPtr->kap = 1.0 / Tau;
    }
  }

  CalReqPtr->maxstd = MaxStd;

  return (CalReqPtr);
}

/*********************************************************************/
/* This code sets the Convolution parameters to their default values */
ConvParamsPtr LGMSetDefaultEvalParms() {
  ConvParams *EvalParamPtr = NULL;

  EvalParamPtr = (ConvParams *)srt_calloc(1, sizeof(ConvParams));
  if (EvalParamPtr == NULL)
    return (NULL);

  EvalParamPtr->gridwidth = 6.0;
  EvalParamPtr->nx = 192;
  EvalParamPtr->stencil = 6.0;
  EvalParamPtr->h = 0.0625;
  EvalParamPtr->killkinks = 1;
  EvalParamPtr->ywidth = 3.0;
  EvalParamPtr->ny = 30;
  EvalParamPtr->yGauss = 2.5;
  EvalParamPtr->hy = 0.2;
  return (EvalParamPtr);
}

/********************/
/* copying routines */
/********************/
static LGMErr copyExerBdry(LGMSwptns *ExerBdryPtr,
                           SrtLgmExerBdryData *lgmExerBdryData) {
  long i;

  if (lgmExerBdryData == NULL)
    return ("no lgmExerBdryData structure");

  lgmExerBdryData->NexerBdry = ExerBdryPtr->n;

  if (lgmExerBdryData->exerBdryArr != NULL)
    srt_free(lgmExerBdryData->exerBdryArr);
  lgmExerBdryData->exerBdryArr = (SrtLgmExerBdry *)srt_calloc(
      lgmExerBdryData->NexerBdry, sizeof(SrtLgmExerBdry));

  if (lgmExerBdryData->exerBdryArr == NULL)
    return ("alloc failed in LGM  , exerbdry");

  for (i = 0; i < lgmExerBdryData->NexerBdry; i++) {
    lgmExerBdryData->exerBdryArr[i].parRate = ExerBdryPtr->Rfix[i];
    lgmExerBdryData->exerBdryArr[i].bgnDate = ExerBdryPtr->tEx[i];
  }
  return (NULL);
}

static LGMErr copyTSData(SigKapTS *SKtsPtr, SrtLgmTSData *TSData) {
  long j;

  if (TSData == NULL)
    return (" no TSData structure");
  TSData->NTS = max(SKtsPtr->numK, SKtsPtr->numS);
  TSData->TSArr = (SrtLgmTS *)srt_calloc(TSData->NTS, sizeof(SrtLgmTS));
  if (TSData->TSArr == NULL)
    return ("alloc failed in LGM  , termstruct");

  for (j = 0; j < TSData->NTS; j++) {
    if (j < SKtsPtr->numS) {
      TSData->TSArr[j].sigma = SKtsPtr->sig[j];
      TSData->TSArr[j].sigmaDt = SKtsPtr->sdate[j];
    } else {
      TSData->TSArr[j].sigma = 0;
      TSData->TSArr[j].sigmaDt = 0;
    }

    if (j < SKtsPtr->numK) {
      TSData->TSArr[j].tauDt = SKtsPtr->kdate[j];
      TSData->TSArr[j].tau = SKtsPtr->kap[j];
      if (fabs(TSData->TSArr[j].tau) < 0.000001)
        TSData->TSArr[j].tau = 1000000.;
      else
        TSData->TSArr[j].tau = 1 / TSData->TSArr[j].tau;
    } else {
      TSData->TSArr[j].tau = 0;
      TSData->TSArr[j].tauDt = 0;
    }
  }
  return (NULL);
}

/* Free arrays that will be used for output  , if already allocated */
static void free_outputs(
    int *convertTS,    /* 1=compute new sigs and taus; 0=don't bother */
    int *findExerBdry, /* 1=find swap rates at exercise boundary; 0=don't bother
                        */
    SrtLgmExerBdryData *lgmExerBdryData, /* ptr to exercise boundary data
                                            structure (NULL => not req'd) */
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
    SrtLgmTSData *lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr *atcTSData)    /* ptr to zeta/G data (NULL => not req'd) */
{
  if (*findExerBdry != 0 && lgmExerBdryData != NULL) {
    lgmExerBdryData->NexerBdry = 0;
    if (lgmExerBdryData->exerBdryArr != NULL) {
      srt_free(lgmExerBdryData->exerBdryArr);
    }
  } else {
    *findExerBdry = 0;
  }

  if (*convertTS != 0 && lgmTSData != NULL) {
    lgmTSData->NTS = 0;
    if (lgmTSData->TSArr != NULL) {
      srt_free(lgmTSData->TSArr);
    }
  } else
    *convertTS = 0;

  if (*convertTS != 0 && atcTSData != NULL && *atcTSData != NULL) {
    if ((*atcTSData)->zdate)
      srt_free((*atcTSData)->zdate);
    if ((*atcTSData)->zeta)
      srt_free((*atcTSData)->zeta);
    if ((*atcTSData)->Gdate)
      srt_free((*atcTSData)->Gdate);
    if ((*atcTSData)->G)
      srt_free((*atcTSData)->G);
    (*atcTSData)->zdate = NULL;
    (*atcTSData)->zeta = NULL;
    (*atcTSData)->Gdate = NULL;
    (*atcTSData)->G = NULL;
  }
}