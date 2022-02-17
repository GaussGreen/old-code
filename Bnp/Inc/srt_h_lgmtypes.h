
#ifndef _SRT_H_LGMTYPES_H_
#define _SRT_H_LGMTYPES_H_

#include "opsabr.h"
/**************************************************************************/
/********* typedefs (enums and structures) for New LGM autocal*************/
/**************************************************************************/
/* Error handling */
typedef char *LGMErr;

/**************************************************************************/
/* Reference Swaptions */
/**************************************************************************/
/* This holds data for a single LGM reference swaption */
typedef struct srtlgmrefswptn {
  Date bgnDt;
  Date endDt;
  double strike;
  double bsPv;
  double bsVol;

} SrtLgmRefSwptn;

/* This holds all LGM reference swaption data */
typedef struct srtlgmrefswptndata {
  /* Arrays of short      , long      , caplets      , fixed term swaption
   * details */
  SrtLgmRefSwptn *refShortSwptnArr;
  SrtLgmRefSwptn *refLongSwptnArr;
  SrtLgmRefSwptn *refCapSwptnArr;
  SrtLgmRefSwptn *refFixSwptnArr;

  int NrefShortSwptn; /* The no. of short term swaptions */
  int NrefLongSwptn;  /* The no. of long term swaptions */
  int NrefCapSwptn;   /* The no. of caplets */
  int NrefFixSwptn;   /* The no. of fixed swaptions */
  int isSemi;         /* Semi-Annual      , as opposed to Annual */

} SrtLgmRefSwptnData;

/**************************************************************************/
/* This holds data for a single exercise boundary swap */
/**************************************************************************/
typedef struct srtlgmexerbdry {
  Date bgnDate;
  double parRate;

} SrtLgmExerBdry;

/* This holds all exercise boundary swap data */
typedef struct srtlgmexerbdrydata {
  SrtLgmExerBdry *exerBdryArr;
  int NexerBdry;

} SrtLgmExerBdryData;

/**************************************************************************/
/* Term structures for LGM model */
/**************************************************************************/

/* Structure for external call */
/* This holds LGM tau/sigma data for a single point */
typedef struct srtlgmts {
  Date tauDt;
  double tau;
  Date sigmaDt;
  double sigma;

} SrtLgmTS;

/* This holds all LGM tau/sigma data */
typedef struct srtlgmtsdata {
  SrtLgmTS *TSArr;
  int NTS;

} SrtLgmTSData;

typedef struct LGM_TS_ {
  long numZ;   /* number of zeta values */
  Date *zdate; /* zdate[0      ,1      , ...      , numZ-1] are dates for zeta
                  values */
  long numG;   /* number of G values */
  Date *Gdate; /* Gdate[0      ,1      , ...      , numG-1] are dates for G
                  values */

  /* ONE FACTOR CASE */
  double *zeta; /* zeta[0      ,1      , ...      , numZ-1] are zeta values at
                   above dates */
  double *
      G; /* G[0      ,1      , ...      , numG-1] are G values at above dates */

  /* TWO FACTOR CASE */
  double *Zeta1;  /* zeta1[0      ,1      , ...      , numZ-1] are zeta1 values
                     at above dates */
  double *Zeta2;  /* zeta2[0      ,1      , ...      , numZ-1] are zeta2 values
                     at above dates */
  double *Zeta12; /* zeta12[0      ,1      , ...      , numZ-1] are zeta12
                     values at above dates */
  double *H1; /* H1[0      ,1      , ...      , numG-1] are H1 values at above
                 dates */
  double *H2; /* H2[0      ,1      , ...      , numG-1] are H2 values at above
                 dates */

  double **Zeta1Scenari; /* Zeta1Scenari[0      ,...      , 2*numZ][0      ,...
                            ,numZ-1] */
  double **Zeta2Scenari;
  double **Zeta12Scenari;
  double **H1Scenari; /* H1Scenari[0      ,...      , 2*numZ][0      ,...
                         ,numG-1] */
  double **H2Scenari;

  double alpha;
  double gamma;
  double rho;
  double one_kappa;

} LGM_TS, *LGM_TSPtr;

/* Warning: numZ and numG MUST be at least 2 ****/
/* Interpolation: Use linear (linear at ends) to evaluate between dates */
/* Always use zeta=zeta(t)-zeta(tNow) in case "today" has changed since
 * calibration */

typedef struct SigKapTS_ {
  long numS;   /* number of sigma values */
  Date *sdate; /* sdate[0      ,1      , ...      , numS-1] are dates for sigma
                  values */
  double *sig; /* sig[0      ,1      , ...      , numS-1] are sigma values
                  between dates */

  long numK;   /* number of kap values */
  Date *kdate; /* kdate[0      ,1      , ...      , numK-1] are dates for kappa
                  values */
  double *kap; /* kap[0      ,1      , ...      , numK-1] are kappa values
                  between dates  */
} SigKapTS, *SigKapTSPtr;

/* Note that
        sigma(t) = sig[i] for sdate[i-1] < t <= sdate[i]
        kappa(t) = kap[i] for kdate[i-1] < t <= kdate[i]       ,
with sigma and kappa constant before the first and after the last date */

/* Forward vol data for callable inverse floater */
typedef struct {
  long fwdVolCpn;
  long fwdVolEx;
  double **fwdVolMat;
} FWD_VOL_STR;

/**************************************************************************/
/* Deal types */
/**************************************************************************/
/* Defined for conveniently passing exercise boundary info */
typedef struct LGMSwptns_ {
  long n;
  Date *tEx; /* exercise date of each swaption */
  Date *tTheoEnd;
  Date *tEnd; /* real end date of each swaption */
  double
      *Rfix; /* fixed rate (assumes standard dcb      , spot lag      , bdr) */
} LGMSwptns, *LGMSwptnsPtr;

/* Deal types understood by calibration module and evaluator **************/
typedef enum LGMDealType_ {
  GenEur = 0,  /* Options with a single exercise date and arbitrary payoff */
  SimpleMidAt, /* Standard Mid-Atlantics */
  GenMidAt,    /* Options with multiple exercise dates      , which can only be
              exercised once (like Mid-atlantics) */
  CallInvFlt,  /* Callable inverse floater */
  CallCapFlt,  /* Callable Cap Floater */
  CallTimeSwap,
  LGMLastDealType
} LGMDealType;

/* Structure for general European option */
typedef struct SrtGenEur_ {
  Date tEx;        /* notification (exercise) date */
  long nPay;       /* number of payments on fixed leg */
  double *Payment; /* [0      ,1      , ...      , nPay-1]      , amount paid on
                      each date */
  Date *tPay;      /* [0      ,...      , nPay-1]      , date of each payment */
} SrtGenEur, *SrtGenEurPtr;

/* If exercised at tEx      , the holder receives all positive payments      ,
   and pays all negative payments. For receivers      , all the payments should
   be positive except for the first payment (a negative payment of the strike)
   made at the settlement date. For payers      , all payments should be
   negative except for the strike. */

/* Structure for standard midatlantic deals */
typedef struct SrtSimMidAt_ {
  /* Fixed leg */
  long nPay;       /* number of payments on fixed leg */
  double *Payment; /* [0      ,1      , ...      , nPay-1]      , amount paid on
                      each date */
  Date *tPay;      /* [0      ,...      , nPay-1]      , date of each payment */
                   /* Exercise information */
  SrtReceiverType PayRec; /* Option to pay (SRT_PAYER) or receive (SRT_RECEIVER)
                             fixed leg */
  long nEx;               /* number of exercise dates */
  long FirstExer;         /* FOR USE BY LGMAUTOCAL */
  Date StartDate;
  Date
      *tEx; /* [0      ,...      ,nEx-1]      , notification (exercise) dates */
  Date *tStart; /* [0      ,...      ,nEx-1]      , settlement (start) dates for
                   each exercise */
  double *Strike; /* [0      ,...      ,nEx-1]      , amount paid on settle for
                     fixed leg at each exercise */
  long *FirstPay; /* [0      ,...      ,nEx-1]      , payment number of first
                     payment received for each exercise */
  double *RedFirstPay; /* [0      ,...      ,nEx-1]      , reduction in amount
                          of first payment received after exercise */

  Date tLast;
  double *MidAtCvgPay;
  double *MidAtCvgFirst;
  Date *HDates;

  double *MidAtStrike;
  double *StDLong;

} SrtSimMidAt, *SrtSimMidAtPtr;

/* The SrtSimMidAt structure assumes that a deal can only be exercised once , on
one of the dates tEx[0      ,1      ,...      ,nEx-1]. Suppose the deal is
exercised at , say      , tEx[j]      , and let	FirstPay[j] be i0. Then the
receiver will receive the series of payments

                Payment[i0] - RedFirstPay[j] paid on tPay[i0]
                Payment[i] paid on tPay[i] for i = i0+1      , i0+2      , ... ,
nPay-1.

In return      , the receiver will pay
                Strike[j] paid at tStart[j]

If PayRec = SRT_RECEIVER      , then the receiver of the fixed leg has the
option to exercise. If PayRec = SRT_PAYER      , then the payer of the fixed leg
has the option to exercise.

WARNING: The first payment is payment 0      , so for vanilla deals      , the
FirstPay array would start FirstPay[0]=0      , meaning that if the deal is
exercised at tEx[0]      , then it will receive payments at tPay[0]      ,
tPay[1]      , ...

Note: LGM autocal uses "FirstExer" in places to store the first exercise date
occurring after tNow. It should not be used to store anything permanent


/* Structure for general midatlantic */
typedef struct SrtGenMidAt_ {
  /* Exercise information */
  long nEx;       /* number of exercise dates */
  long FirstExer; /* only consider exercises k >= FirstExer */
  Date *tEx; /* [0      ,...      ,nEx-1]      , exercise (notification date) */
  long *nPay;       /* [0      ,...      ,nEx-1]      , number of payments if
                       exercised at each       date */
  double **Payment; /* Payment[j][i]      , j=0      ,1      ,...      ,nEx-1 ,
                       i=0      ,1      ,... ,nPay[j]-1 at */
  Date **tPay; /* tPay[j][i]      , j=0      ,1      ,...      ,nEx-1      , i=0
                  ,1      ,... ,nPay[j]-1 are the */
  /* the payments received if exercised at tEx[j] */
} SrtGenMidAt, *SrtGenMidAtPtr;

/* This structure assumes that a deal can only be exercised once      , on any
of the dates tEx[0      ,1      ,...      ,nEx-1]. If the deal is exercised at
, say      , tEx[j] , then the option exerciser will receive the series of
payments Payment[j][i] paid at tPay[j][i]      , i=0      ,1      ,...
,nPay[j]-1. If a Payment is negative      , then it represents an amount the
exerciser has to pay on tPay. To change between payers & receivers      , we
just need to flip the signs of the payments.

Note: LGM autocal uses "FirstExer" in places to store the first exercise date
occurring after tNow. It should not be used to store anything permanent */

/***********************************************************************/
/* American deal types understood by calibration module and evaluator **/
typedef enum LGMAmerType_ {
  SimAmer = 0, /* Options with a single exercise date and arbitrary payoff */
  LGMLastAmerType
} LGMAmerType;

/* Structure for simple American swaptions */
typedef struct SrtSimAmer_ {
  /* Fixed leg */
  long nfix;       /* number of payments on fixed leg */
  Date *tfixStart; /* [0      ,1      , ...      , nfix-1]      , start date of
                      each fixed period */
  Date *tfixEnd; /* [0      ,1      , ...      , nfix-1]      , end date of each
                    fixed period */
  Date *tfixPay; /* [0      ,1      , ...      , nfix-1]      , pay date for
                    each fixed period */
  double *fixCoupon; /* [0      ,1      , ...      , nfix-1]      , fixed leg
                        interest payments (no notional in final payment) */
  SrtBasisCode fixBasis; /* day count basis for fixed leg payments */
  int EarlyFlagFix; /* 0=subtract accruals from first payment; 1=pays accrual at
                       start */
                    /* Floating leg */
  long nflt;        /* number of floating payments !! */
  Date *tfltFixing; /* [0      ,1      , ...      , nflt-1]      , fixing dates
                       for each floating period */
  Date *tfltPay;    /* [0      ,1      , ...      , nflt-1]      , pay dates for
                       each floating    period */
  SrtBasisCode
      fltBasis; /* day count basis for floating leg payments (usually ACT360) */
  int EarlyFlagFlt; /* 0=subtract accruals from first payment; 1=pays accrual at
                       start */
  int ResetFlt;     /* 1=floating rate for first payment (stub) is reset upon
                       exercise */
                    /* Exercise information */
  SrtReceiverType PayRec; /* option to pay (SRT_PAYER) or receive (SRT_RECEIVER)
                             fixed leg */
  Date tFirstExer;        /* first allowed exercise date for deal */
  int lagExerSettle;      /* number of days between exercise (notification) and
                             settlement */
  int CalBusLag;          /* 0=lag is calendar days; 1=lag is business days */
  BusDayConv convSettle;  /* business day convention for the settlement date (if
                             calendar) */
  double *ExtraPrem; /* [0      ,1      , ...      , nfix-1]      , extra amount
                        paid on settlement */
} SrtSimAmer, *SrtSimAmerPtr;

/*
CAUTION: SrtSrtSimAmer swaption structure assumes that the notional is 1.
CAUTION: The final fixPayment should include only the interest rate coupon , NOT
the notional. CAUTION: Any floating leg spread (the difference between the
floating leg and cash) must be applied BEFORE the code is called. CAUTION: Sign
of the exercise premium: ExtraPrem>0 increases the value of the floating leg and
ExtraPrem<0 decreases the value of the floating leg      , regardless of whether
the deal is a payer of receiver. If the deal has an "exercise premium" which is
paid by the option owner to exercise      , then the sign must be changed if it
is a PAYER.

Deals which represent the cancellation feature of cancellable swaps should have
EarlyFlags of 1      , and typically ResetFlt 0. Deals which represent the rarer
"exercise into" swaptions should have EarlyFlags 0      , and ResetFlt may or
may not be 1. See below.

Fixed leg:
        The interest rate coupon for the period tfixStart[i] to tfixEnd[i] is
fixCoupon[i] and is paid on tfixPay[i]      , i=0      ,1      , ...      ,
nfix-1. The notional of 1 is also paid on tfixPay[nfix-1] Floating leg: The
interest rate for the period tfltPay[i-1] to tfltPay[i] is set (fixed) on
tfltFixing[i] and paid on tfltPay[i]      , i=0      ,1      , ...      ,
nflt-1. The notional is also paid on tfltPay[nflt-1]. The code assumes that the
start date of the first period is tfixStart[0]      , which causes no
appreciable errors.

Exercise:
        If PayRec = SRT_RECEIVER      , then the receiver of the fixed leg has
the option to exercise. If PayRec = SRT_PAYER      , then the payer of the fixed
leg has the option to exercise. This structure assumes the deal can be exercised
on any day that is on or after both tFirstExer and tFirst      , where tFirst is
either the evaluation date (if end-of-day is off) or the evaluation date plus 1
business day (if end-of-day is on). Neither the evaluation date nor the
end-of-day flag are part of the deal structure; they are supplied to the
valuation code independently.

If the deal is exercised on tEx      , the theoretical settlement date is
tSettle_theor      , which is lagExerSettle (calendar or business) days after
tEx. If business days      , then the actual settlement date tSettle is the same
as the theoretical settlement date. If calendar days      , then the business
day rule convSettle is applied to get the actual settlement date from the
theoretical date.

Suppose that the deal is exercised on tEx      , with actual settlement date
tSettle , and suppose that the settlement date falls in the k-th fixed period ,
        tfixEnd[k-1] <= tSettle < tfixEnd[k]
and the jth floating period:
        tfltPay[j-1] <= tSettle < tfltPay[j]

If EarlyFlagFix is 0      , then the first fixed payment the receiver gets is
reduced by the interest that has already accrued in the then current period , k.
So the receiver gets the payments fixCoupon[k](1 - coverage(tfixStart[k]      ,
tSettle)/coverage(tfixStart[k]      , tfixEnd[k])) paid at tfixPay[k]      ,
        fixCoupon[i] paid at tfixPay[i]      , i=k+1      , k+2      , ... ,
nfix-2 1 + fixCoupon[nfix-1] paid at tfixPay[nfix-1]

If EarlyFlagFix is 1      , then the receiver must pay the accrued interest
        fixCoupon[k] * coverage(tfixStart[k]      ,
tSettle)/coverage(tfixStart[k] , tfixEnd[k])) on the settlement date tSettle ,
but then receives the entire amount of all remaining coupons: fixCoupon[i] paid
at tfixPay[i]      , i=k      , k+1 , ...      , nfix-2 1 + fixCoupon[nfix-1]
paid at tfixPay[nfix-1]

In return for receiving the fixed leg payments      , the receiver makes
floating leg payments. Apart from some corrections (see below)      , these
floating rate payments are assumed to be worth the notional (of 1) plus any
extra exercise premium paid at settlement. ANY SPREAD due to the difference
between cash and the floating leg has to be accounted for BEFORE calling
LGMAmerican.

To be more precise      , suppose that the settlement date occurs in the k-th
fixed period: tfixEnd[k-1] <= tSettle < tfixEnd[k]. Then the floating leg is
assumed to be worth 1 + ExtraPremium[k] + "corrections" paid at tSettle.

Note that ExtraPrem[k]>0 INCREASES the value of the floating leg. If ExtraPrem
is an exercise premium      , paid by the exerciser      , then the sign of
ExtraPrem must be adjusted depending on whether the deal is a Payer or Receiver.

The corrections are:
        If EarlyFlagFlt is 0      , then the first floating payment received is
the interest for the period tSettle to tfltPay[k]. If ResetFlt is 0      , then
the rate for this first period is the previously fixed LIBOR rate which was
fixed at par on tfltFixing[k] for the period tfltPay[k-1] to tfltPay[k]); this
causes the floating leg to not be worth the same as 1 paid on tSettle. If
ResetFlt is 1      , then the rate for the first period is fixed on or after the
exercise date for the period tSettle to tfltPay[k]      , and the floating leg
is worth exactly 1 paid on tSettle. If EarlyFlagFlt is 1      , then the first
payment received is the floating interest rate applied to the whole period
tfltPay[k-1] to tfltPay[k] paid on tfltPay[k]      , but the floating leg has to
pay on the floating interest that has accrued from tfltPay[k-1] to tSettle on
tSettle. This is also slightly different than 1 being paid at tSettle.

If ResetFlt is 1 (usual) and if we are already beyond the fixing date for the
floating period that includes tFirstExer      , then accurate valuation requires
the already-fixed current floating rate. This is supplied to the code outside
the deal structure
*/

/**********************************************************************************/
/* If we ever get a significant book of deals in which there is a pattern
between exercise dates and the payments (as there is for MidAts) we can speed up
the re-val by developing special structures for the deal      , and developing
special payoff functions which take advantage of the relationships */
/**********************************************************************************/

/* Structure for Bermudan inverse floater */
typedef struct SrtCallInvFlt_ {

  /*	Exercise	*/

  long nEx;       /* number of exercise dates	*/
  Date *tEx;      /* [0      ,1      ,...      ,nEx-1] notification dates */
  Date *tSet;     /* [0      ,1      ,...      ,nEx-1] settlement dates */
  long *iSet;     /* [0      ,1      ,...      ,nEx-1] start coupon idx */
  double *strike; /* [0      ,1      ,...      ,nEx-1] fee paid by option holder
                     to exercise */
  SrtReceiverType PayRec; /* RECEIVER or PAYER */
  long FirstEx;           /* Required for evaluation */

  /*	Coupons: each coupon = fix - gear * libor (start      , end) * cvg
   */

  long nCpn;       /* number of coupon periods */
  Date *tCpnStart; /* [0      ,1      ,...      ,nCpn-1] coupon start date */
  Date *tCpnPay;   /* [0      ,1      ,...      ,nCpn-1] coupon pay date */

  double *a; /* [0      ,1      ,...      ,nCpn-1] fixed coupon      , total
              * coupon is a - gear cash libor * cvg */
  double *gear;
  double *cvg;
  double *lcvg;

  /*	Cap on cash libor	*/

  double *cap_str; /* [0      ,1      ,...      ,nCpn-1] strikes of cash libor
                      cap */

  /*	Record of fwd vol */
  long fwdVolCpn;
  long fwdVolEx;
  double **fwdVolMat;

} SrtCallInvFlt, *SrtCallInvFltPtr;

/**********************************************************************************/
/* Structure for Bermudan Cap Floater */
typedef struct SrtCallCapFlt_ {
  /* Coupon leg */
  long nCpn;    /* number of coupon periods */
  Date *tCpn;   /* [0      ,1      , ...      , nCpn]      , all coupon dates */
  double *amax; /* [*      ,1      , ...      , nCpn]      , effective max
                   handle for each coupon */
  double *amin; /* [*      ,1      , ...      , nCpn]      , effective min
                   handle for each coupon */
  double *marg; /* [*      ,1      , ...      , nCpn]      , effect margin */
  SrtBasisCode aBasis;   /* day count basis for rate on coupon leg */
  double lvg;            /* leverage */
  SrtBasisCode cpnBasis; /* day count basis for fixed leg payments */
                         /* Exercise information */
  long nEx;              /* number of exercise dates */
  long FirstEx;          /* Reserved for LGMautocal */
  Date *tEx;  /* [0      ,1      ,...      ,nEx-1] notification dates; deal
                 starts on next  coupon date */
  long *iSet; /* [0      ,1      ,...      ,nEx-1] on exer      , settlement is
                 tCpn[iSet[j]] */
  double *strike; /* [0      ,1      ,...      ,nEx-1] fee paid by option holder
                     to exercise */
  SrtReceiverType PayRec; /* RECEIVER or PAYER */
} SrtCallCapFlt, *SrtCallCapFltPtr;

/* Cash model for a Bermudan option on a Cap floater
Coupon leg:
        C[i] = cvg(t[i-1]      ,t[i]) * (margin[i] + lvg*max{rate[i]-amin[i] ,0}
- lvg*max{rate[i]-amax[i]      ,0}) received at t[i] for i = 1      ,2      ,...
,n where rate[i] is the true (ie      , cash) floating rate for the period
t[i-1] to t[i]

Exercise:
        At any notification date tEx[j]      , j=0      ,1      , ...      ,
nEx-1      , the holder can exercise. If exercised at tEx[j] the settlement date
is tCpn[iSet[j]]      , the next coupon date on or after tEx[j].

The receiver
        receives all remaining coupons C[i]      , i = iSet[j]+1      ,
iSet[j]+2      ,... , n receives notional ($1) at the end date tCpn[n] pays
strike[j] on the settlement date tCpn[i0-1]

In the coupon      , rate[i] is ASSUMED to be the cash FRA rate for the period
t[i-1] to t[i]; if the rate has a basis spread      , then this must be
incorporated into the handle a[i]. Any margin or basis spread on the floating
leg must be moved into the margin on the coupon leg */

/**********typedefs for the calibration module ****************************/

/* Market conventions */
typedef struct LGMMarkConv_ {
  CcyCode ccy;           /* currency (for holidays) */
  int lag;               /* standard spot lag (business days)*/
  SrtCompounding cfreq;  /* standard frequency for caps */
  SrtBasisCode cbasis;   /* standard day count basis for caplets */
  SrtBusDayConv cbdconv; /* standard business day convention for caplets */
  SrtCompounding sfreq;  /* standard frequency for swaps */
  SrtBasisCode sbasis;   /* standard day count basis for swaps */
  SrtBusDayConv sbdconv; /* standard business day convention for swaps */
  int minswap; /* number of fixed leg periods in minimum length swap */
} LGMMarkConv;

/* Calibration method enum */
typedef enum LGMCalMeth_ {
  FixKappa = 0, /* Use constant kappa (=1/tau) to get G(t)      , calib on
               swaptions or caplets for zeta(t) */
  GivenG,       /* G(t) given      , calib on swaptions or caplets for zeta */
  FixSigma,  /* Use constant sigma & calib on swaptions or caplets for G(t) */
  GivenZeta, /* zeta(t) given      , calib on swaptions or caplets for G(t) */
  FixExp, /* Calib on 1 into k swaptions to get G(t)      , calib on swaptions
         or caplets for zeta */
  TenorAndDiag, /* Calib on column (fixed tenor swaptions or caplets) & diagonal
               swaptions (old method) */
  FullCapAndDiag,
  HybridShortAndDiag,
  LGMLastCalMeth /* No more calibration methods */

} LGMCalMeth;

typedef enum LGMCalType_ {
  FIXEDCAL = 0,
  GLOBALCAL

} LGMCalType;

typedef enum LGMCAPVOLMeth_ {
  MODEL = 0,
  CONVERGING,
  SLIDING2,
  MARKET

} LGMCAPVOLMeth;

/* Strike method */
typedef enum LGMRMeth_ {
  IRR = 0,   /* strikes by proportional IRR method */
  dIRR,      /* long strikes by prop IRR      , short by delta IRR */
  addshift,  /* add Rdata to ATM swap rates */
  propshift, /* multiply ATM swap rates by (1+Rdata) */
  fracstd,   /* add Rdata*(std dev) to ATM swap rate */
  givenR,    /* use Rdata as strikes */
  EMK1,
  EMK2,
  LastRMeth /* No more methods */
} LGMRMeth;

/* Determines precise calibration to be done */
typedef struct LGMCalParm_ {
  LGMCalMeth calmeth; /*  method used */
  LGMRMeth Rmeth;
  LGMCAPVOLMeth CapletVolMeth;
  int respectExer; /* 1 = respect exer dates of deal in calib      , 0 = use
                      regular set of dates */
  int usecaps;     /* 1 = use caplets for one set of ref instruments      , 0 =
                      swaptions only */
  long MinMonToEx; /* time to earliest allowed exercise date for 1 into k
                      swaptions */
  int keep1intok;  /* 1 = give preference to zeta from 1 into k      , if
                      conflict  exists */
  /* 0 = give preference to other instruments */
  /* -1 = never use zeta from 1 into k */
  /* kappa information for fixed kappa method */
  double kap;    /* kappa to use if its a fixed kappa calculation */
  int usestarts; /* 1=G on start dates      , 0=G on standard pay dates (when
                    possible) */
                 /* G information for given G method */
  long numG;     /* number of G values */
  Date *Gdate;   /* G values for given G(t) method */
  double *G;     /* G values for given G(t) method */
                 /* zeta information for given zeta method */
  long numZ;     /* number of zeta values */
  Date *zdate;   /* dates of zeta values for given zeta(t) method */
  double *zeta;  /* zeta values for given zeta(t) method */
  /* strike information for addshift      , propshift      , fracstd      , abs
   * strike methods */
  long numR1;     /* number of values for long swaptions */
  Date *Rdate1;   /* dates of the strike data for long swaptions */
  double *Rdata1; /* values of the strike data for long swaptions */

  long numR2;     /* number of values for short swaptions */
  Date *Rdate2;   /* dates of strike data for short swaptions or caplets */
  double *Rdata2; /* values of the strike data for short swaptions */

  /* Autocal 2F */
  int LGMOneTwoFactor;
  double alpha;
  double gamma;
  double rho;

  double *Zeta1;
  double maxstd;

  /* Starting points for the vega && delta report */
  long *Zeta1Dates;
  double *StartZeta1s;
  long *TauDates;
  double *StartTaus;
  double **HybridShortInstrsIndex;

} LGMCalParm, *LGMCalParmPtr;

/* Contains information needed to construct all sets of calibration instruments
 */

typedef struct LGMCalSet_ {
  /* fixed leg pay dates for swaptions */
  long n;     /* tPay[n] is last date we need G(t) for */
  long nPay;  /* number of pay dates */
  Date *tPay; /* tPay[*      ,1      ,...      ,nPay] dates for swaption fixed
                 legs */
  double *cvgpay; /* cvgpay[*      ,1      ,...      ,nPay] cvg(tPay[i-1]
                     ,tPay[i]) */
  double *Dpay;   /* Dpay[*      ,1      ,...      ,nPay] discount factor to
                     tPay[i] */

  /* common exercise and start dates */
  Date tNow;    /* evaluation date */
  long nEx;     /* number of exercise dates */
  Date *tEx;    /* tEx[*      ,1      ,...      ,nEx] exercise (fixing) dates */
  Date *tStart; /* tStart[*      ,1      ,...      ,nEx] start (settlement)
                   dates */
  double *DStart; /* DStart[*      ,1      ,...      ,nEx] discount factor to
                     tStart[j] */

  /* long swaptions */
  int longflag; /* 1 means long swaptions have been created */
  long *ifirst; /* ifirst[*      ,1      ,...      ,nEx] fixed leg is i =
                   ifirst[j]      , ... , nlong[j] */
  long *nlong;
  double *cvgfirst; /* cvgfirst[*      ,1      ,...      ,nEx] cvg for first
                       payment of swaption j */
  double *FrLong;
  double *Rflong;   /* Rflong[*      ,1      ,...      ,nEx] fixed rate of long
                       swaption j */
  double *Vlong;    /* Vlong[*      ,1      ,...      ,nEx] market value of long
                       swaption j */
  double *VegaLong; /* VegaLong[*      ,1      ,...      ,nEx] vega of long
                       swaption j */
  double
      *StDLong; /* StDLong[*      ,1      ,...      ,nEx] Stdev of swaption j */
  double *BetaLong;
  double *CEVLong;

  /* short swaptions */
  int shortflag;   /* 1 means short swaptions and pairs have been created */
  long *nshort;    /* fixed leg is i = ifirst[j]      , ...      , nshort[j] */
  double *Rfshort; /* Rfshort[*      ,1      ,...      ,nEx] fixed rate of short
                      swaption j */
  double *Vshort;  /* Vsh[*      ,1      ,...      ,nEx] market value of short
                      swaption j */

  /* caplets */
  int capflag; /* 1 means caps have been created */
  Date *tEnd;  /* tEnd[*      ,1      ,...      ,nEx] fixing tEx[j]      , start
                  tStart[j]      , end  tEnd[j] */
  double *
      Dcap; /* Dcap[*      ,1      ,...      ,nEx] discount factor at tEnd[j] */
  double *cvgcap; /* cvgcap[*      ,1      ,...      ,nEx] cvg from tStart[j] to
                     tEnd[j] */
  double *Rfcap; /* Rfcap[*      ,1      ,...      ,nEx] fixed rate for each cap
                    j */
  double
      *Vcap; /* Vcap[*      ,1      ,...      ,nEx] market value of caplet j */
  double *VegaCap; /* VegaCap[*      ,1      ,...      ,nEx] vega of caplet j */

  /* 1 into k swaptions */
  int fixflag;  /* 1 means 1 into k have been created */
  Date ftEx;    /* exercise dates are all ftEx */
  Date ftStart; /* start dates are all ftStart */
  long ffirst;  /* pay dates of swap. j are tPay[ffirst      , ...      ,j] */
  long nfirst; /* the swaptions are j = nfirst      , nfirst+1      , ...      ,
                  nlast */
  long nlast;
  double
      cvgfix; /* cvgfix is the cvg for first period ftStart to tPay[ffirst] */
  double DfixStart; /* Dfixst is discount factor at ftStart[i] */
  double *Rfix;     /* Rfix[*      ,...      ,nlast] Rfix[j] is fixed rate for
                       swaption j */
  double *Vfix;     /* Vfix[*      ,...      ,nlast] Vfix[j] is fixed rate for
                       swaption j */
  long ExIndex;
  long MAXCpn;

} LGMCalSet, *LGMCalSetPtr;

/**********typedefs for the evalutation module ****************************/
/* Convolution parameters */
/* These parameters control the numerics used in evalutation */
typedef struct ConvParams_ {
  double gridwidth; /* (default 6.0) grid runs from -gridwidth*stddev < x <
                       +gridwidth*stddev */
  long nx;          /* (default 96) number of x gridpoints */
  double
      stencil;   /* (default 6.0) width of the Gaussian distribution included */
  double h;      /* (default is 0.0625) delta z of Gaussian */
  int killkinks; /* (default is 1) correct for kinks unless zero */

  /* the following variables are used for the 2D convoler */

  double ywidth; /* (default 3.0) grid runs from -ywidth*stddev < y <
                    +ywidth*stddev */
  long ny;       /* (default 30) number of y gridpoints */
  double yGauss; /* (default 2.5) width of the Gaussian distribution included */
  double hy;     /* (default is 0.2) delta v of Gaussian */

} ConvParams, *ConvParamsPtr;

typedef struct HermiteParms_ {
  long NbrHermitePoints;
  double *HermiteWeights;
  double *HermitePoints;

} HermiteParms;

/*-----------------2D functions ---------------------*/

LGMErr
GenEurpayoff2D(double **payoff, long nx, double *x, double *reductionx, long ny,
               double *y, double *reductiony, LGMDealType dealtype,
               void *dealPtr, Date tNow, long lj0, String ycName, LGM_TS *tsPtr,
               LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                                double *), /* swaption/cap vols */
               LGMErr (*GetBeta)(Date, Date,
                                 double *)); /* swaption/cap exponents (beta) */

LGMErr Convolver2D(
    /* info about convolutions */
    long nEx,                  /* number of exercises */
    Date *tEx, double *zeta11, /* [0      , 1      , ...      , nEx-1]      ,
                          values of zeta11 at the exercise dates */
    double *zeta12, /* [0      , 1      , ...      , nEx-1]      , values of
                   zeta12 at the exercise dates */
    double *zeta22, /* [0      , 1      , ...      , nEx-1]      , values of
                   zeta22 at the exercise dates */
    double *G1, /* [0      , 1      , ...      , nEx-1]      , values of G1 at
                   the exercise dates */
    double *G2, /* [0      , 1      , ...      , nEx-1]      , values of G2 at
                   the exercise dates */
    ConvParams *parms, /* convolution numerical constants */

    /* info about today's discount curve and swaption/cap vols */
    Date EvalDate, /* tNow for deal evaluation */
    String ycName, /* yield curve name */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *),              /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double *), /* swaption/cap exponents (beta) */

    /* information about the deal */
    LGMDealType dealtype, void *dealPtr, LGMErr (*payofffunc)(),

    /* output */
    double *answer,    /* value of the deal */
    double **xExBdry); /* array of exercise points (in x) */

#endif
