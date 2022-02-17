#ifndef _SRT_H_LGM_US_PROTOS_H
#define _SRT_H_LGM_US_PROTOS_H

#include "num_h_maxmin.h"
#include "srt_h_lgmUStypes.h"

/******************************************************************/
/******************************************************************/
/* LGMCaller.c routines */
/******************************************************************/
/******************************************************************/
/* These programs have the identical calls as the old LGMautocal  , except
that input parameters to choose the calibration method have been added.
In addition  , in the american call  , rflt_current is a new input parameter
which is ignored if it is zero or negative.
It uses its inputs to create the call to the LGMnewautocal program. */
/* NOTE: this program is a temporizing step  , which does not use the
flexibility of the new program  , and it should be replaced with programs which
call the new LGMautocal directly from MAD and Westminster  , which can expose
the new functionality  , ASAP */

// -------------------------------------------------------------------------------------------------------------------------------------------------
// //
//
// LGMCallableTimeSwapCaller
//
LGMErr LGMCallableTimeSwapCaller(
    //	Exercise
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
    const char
        *PayRecStr, /* RECEIVER (call on the bond) or PAYER (put on the bond) */

    /*	Coupons: the structure must be decomposed into a swap that delivers:
- on the funding side  , Libor CASH FLAT
- on the exotic side  , (cpn + gear*accrued_index)*accrued_ratio*coverage */
    double *cpn,                /* [0  ,1  ,...  ,nCpn-1]  coupons*/
    double *tgear_floatdigital, /*[0  ,1  ,...  ,nCpn-1] gearings for floating
                                   digitals */
    double *cvgCpn,             /* [0  ,1  ,...  ,nCpn-1]  coverages  */
    long nCpn,                  /* number of coupon periods */
    Date *tCpnStart,            /* [0  ,1  ,...  ,nCpn-1] coupon start date */
    Date *tCpnPay,              /* [0  ,1  ,...  ,nCpn-1] coupon pay date */
    double *tFundingPayment, /* [0  ,1  ,...  ,nCpn-1] payment of the floating
                                leg */

    // Spread:  not currently supported
    double *dvSpread,

    // floating digital
    SrtCompounding freqAccrued, double CorrelStart, double CorrelEnd,

    // Barriers details
    double call_spread_up, double call_spread_low, double **barriers,
    int iAccrueOnBarrier,

    // Index details
    char *szIndexRefRate, SrtCompounding freqIndex,

    //	Market details
    String ycName, /* yield curve name */
    String vcName, /* Vol curve name */

    // Stuff needed for the generic autocal code
    LGMErr (*GetVol)(long, long, double, SRT_Boolean, double *),
    /* Autocal volatility function to get market cash vol */
    char *char_vol_type, /* normal or log normal swaption vols */

    //	Calibration
    double tau, /* if fixed tau  , use this value for tau (in years) */
    int calibrationmethod, /* 1 = Swaptions ATM  , 2 = Swaptions Strike eq  , 3
                              = Caplets ATM  */
    int timeswapvolmethod, /* 1 = Model  , 2 = Sliding  , 3 = Converging  */
    int atmvolmethod,      /* 1 = Lognormal  , 2 = Normal  , 3 = SigmaBeta */
    long observation_freq, long num_subdiscret_obsfreq, double shiftvol,
    long nx,       /* number of state variable steps in the Pat's integral */
    double maxstd, /* Maximum number of std between forward and strike */
    int nofunding,
    double dStrikeVolCap, /* Maximum lognormal std dev for smile around ATM
                             (i.e.  , vol is flat outside) */

    // Exercise Stuff
    int isExercised,  // 0 not exercised; 1 exercised
    long lExDate,     // date when exercised
    int eod_fix_flag, /*	EOD Fixing Flag 0: I  , 1: E */
    int eod_pay_flag, /*	EOD Payment Flag 0: I  , 1: E */
    int eod_ex_flag,  /*	EOD Exercise Flag 0: I  , 1: E */

    // For period started in the past
    int nDaysInFirstCoupon, long *lvFirstCouponFixingDates,
    double *dvFirstCouponAccrualWeights, double *dvFirstCouponFixingOrSpread,
    double dFloatCouponFixing, long lPayDate,

    /* Accrual Period Day Rule -- not currently used:  1: only use business days
       , 3: count holidays using past setting */
    int iAccrualPeriodDayRule,

    /* Coupon Index Lag -- not currently used:  number of business days before
       start that each index sets */
    int iCouponIndexLag,

    /* End Lag -- not currently used:  number of business days before the end of
       the period that we stop checking the index */
    int iEndLag,

    /* calibrated term structure */
    SrtLgmTSData *lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr *atcTSData,    /* ptr to zeta/G data (NULL => not req'd) */

    /*	Results	*/
    double *Option, double *FixedLeg,
    //	double  *dFloatLeg  ,
    double **forwardTimeSwap, /*[0..1]*[0..nEx-1]  the first is the forward
                                 computed */
    double *dvStrikes         /* the strikes used in the calibration. */
);
//
// LGMCallableTimeSwapCaller
//
// -------------------------------------------------------------------------------------------------------------------------------------------------
// //

// -------------------------------------------------------------------------------------------------------------------------------------------------
// //
//
// calcRangeAccrual
//
LGMErr calcRangeAccrual(long lPayDate, double dCoupon, int iEODFixFlag,
                        int iAccrueOnBarrier, int nDaysInCoupon,
                        long *lvFixingDates, double *dvWeights,
                        double *dvCouponFixingOrSpread, double dUpperBarrier,
                        double dLowerBarrier, double dUpperSpread,
                        double dLowerSpread, char *szRefRate, double dCvg,
                        double *out_pdFixed);

//
// calcRangeAccrual
//
// -------------------------------------------------------------------------------------------------------------------------------------------------
// //

/* Create a Callable time swap structure  , and copy input data into it */
LGMErr LGMFillCallableTimeSwap(
    SrtCallTimeSwap **dealPtrPtr,   /* output: option deal */
    SrtCallTimeSwap **dealTSPtrPtr, /* output: time swap deal */
    long *nExleftPtr,               /* output: effective number of exer */
    Date tfirst, char *ycName,      /* info about today */
    long n, Date *tStartCpn, Date *tPayCpn, /* coupon leg */
    double *tCpn, double *tgear_floatdigital, double *tCvgCpn,
    double *tFundingPayment, long nEx, Date *tEx, Date *tExStart,
    double *exFee, /* exercise info */
    double call_spread_up, double call_spread_low, long observation_freq,
    double **barriers, SrtCompounding freqIndex, SrtCompounding freqAccrued,
    double CorrelStart, double CorrelEnd,
    SrtReceiverType payrec, /* PAYER or RECEIVER */
    int timeswapvolmethod, int atmvolmethod, int calibrationmethod,
    long num_subdiscret_obsfreq, double shiftvol, long nx, int nofunding,
    char *vol_id,
    LGMErr (*GetSABRParam)(long, long, char *, SABR_VOL_TYPE, double *),
    double exFeeTimeSwap, int eod_pay_flag, /*	0: I  , 1: E */
    int eod_fix_flag);                      /*	0: I  , 1: E */

LGMErr CompValTimeSwap(
    Date tNow,                        /* Value date of the deal */
    String ycName,                    /* yield curve */
    double *ratioPtr,                 /* ratio */
    long indexexo, double *intValPtr, /* intrinsic value */
    double Strike, /* fee to pay by the folder of the option */
    SrtCallTimeSwapPtr deal, int bullet_or_option, /* 0=bullet  ,1=option */
    int eod_pay_flag,
    int noDaysAccrued, /* for period which has already started */
    double floatingFixing /* fixing at the beginning of the period */);

LGMErr CompVal_ForwardsTimeSwaps(SrtCallTimeSwapPtr deal, Date tNow,
                                 String ycName,
                                 double *ratios, /* first index is 1 !!! */
                                 double *intVal, double *forwardsTimeSwaps,
                                 double *pvfixedleg, int nEx, int eod_pay_flag,
                                 int noDaysAccrued, double floatingfixing);

LGMErr CompVal_TimeSwap(SrtCallTimeSwapPtr deal, Date tNow, String ycName,
                        double *pvfixedleg, int eod_pay_flag, int noDaysAccrued,
                        double floatingfixing);

LGMErr GetVol_FudgeTimeSwap(long tNow, long value_date, long start, long end,
                            double strike, double fra, SrtCallTimeSwapPtr deal,
                            double *result);

/* Computes a fra in LGM autocal */
LGMErr fra_in_LGMAutocal(String ycName, LGM_TS *tsPtr, long value_date,
                         long tNow, long start, long end,
                         double xred, /* state variable in LGM autocal */
                         double *result);

void NewConvolverTimeSwap(long nEx, /* number of exercises */
                          long nx, /* number of values for the state variable */
                          long n0, long nz, double dx, double h, long m,
                          long killkinks, double *reduction, /*[0..nEx-1] */
                          double *weights,                   /* [0..nx] */
                          double *varr,                      /*[0..nEx-1] */
                          double **payofftable, double *answer);

LGMErr ComputeBermudanAndEuropeanTimeSwap(
    /* info about convolutions */
    long nEx,     /* number of exercises */
    double *zeta, /* [0  , 1  , ...  , nEx-1]  , values of zeta at the exercise
                     dates */
    ConvParams *parms, /* convolution numerical constants */

    /* info about today's discount curve */
    Date EvalDate, /* tNow ... the evaluation date */
    String ycname, /* yield curve name for discount factors */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *),              /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double *), /* swaption/cap exponents (beta) */
    LGM_TS *tsPtr,

    /* information aobut the deal */
    LGMDealType dealtype, void *dealptr, LGMErr (*payofffunc)(),

    LGMCalParm *CalReq,

    /* output */
    double *answer,    /* value of the deal */
    double **xExBdry); /* array of exercise points (in x) */

/* Calculates the payoff at each exercise date for Bermudan time swap */
LGMErr BerTimeSwapPayoff(
    double **payoff, long nx, double *x, double *reduc, LGMDealType dealtype,
    void *dealPtr, LGMCalParm *CalReq, Date tNow, String ycName, LGM_TS *tsPtr,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date,
                      double *)); /* swaption/cap exponents (beta) */

/******************************************************************/
/******************************************************************/
/* LGMUtil.c routines */
/******************************************************************/
/******* Utilities for Deal structures ****************************/

/* Create a callable time swap structure */
SrtCallTimeSwapPtr LGMCreateCallableTimeSwap(long nEx, long n,
                                             long observation_freq);

/* Frees a Callable Time Swap structure */
void LGMFreeCallableTimeSwap(SrtCallTimeSwap **ptr);

/* Check validity of a Callable Time Swap */
LGMErr LGMValidCallTimeSwap(SrtCallTimeSwap *ptr, Date tfirst, long *EffnEx);

/* Sort function called by the wrapper for the autocal time swap */
LGMErr LGMCallableTimeSwapCaller2(
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
    const char
        *PayRecStr, /* RECEIVER (call on the bond) or PAYER (put on the bond) */

    /*	Coupons: the structure must be decomposed into a swap that delivers:
- on the funding side  , Libor CASH FLAT
- on the exotic side  , (cpn + gear*accrued_index)*accrued_ratio*coverage */
    double *cpn, /* [0  ,1  ,...  ,nCpn-1]  coupons*/
    double *cpn2, double *tgear_floatdigital, /*[0  ,1  ,...  ,nCpn-1] gearings
                                                 for floating digitals */
    double *cvgCpn,          /* [0  ,1  ,...  ,nCpn-1]  coverages  */
    long nCpn,               /* number of coupon periods */
    Date *tCpnStart,         /* [0  ,1  ,...  ,nCpn-1] coupon start date */
    Date *tCpnPay,           /* [0  ,1  ,...  ,nCpn-1] coupon pay date */
    double *tFundingPayment, /* [0  ,1  ,...  ,nCpn-1] payment of the floating
                                leg */

    /* floating digital */
    SrtCompounding freqAccrued, double CorrelStart, double CorrelEnd,

    /* Barriers details */
    double call_spread_up, double call_spread_low, double **barriers,

    /* Index details */
    SrtCompounding freqIndex,

    /*	Market details	*/
    String ycName, /* yield curve name */
    String vcName, /* Vol curve name */

    /* Stuff needed for the generic autocal code */
    LGMErr (*GetVol)(long, long, double, SRT_Boolean, double *),
    /* Autocal volatility function to get market cash vol */
    char *char_vol_type, /* normal or log normal swaption vols */

    /*	Calibration	*/
    double tau, /* if fixed tau  , use this value for tau (in years) */
    int calibrationmethod, /* 1 = Swaptions ATM  , 2 = Swaptions Strike eq  , 3
                              = Caplets ATM  */
    int timeswapvolmethod, /* 1 = Model  , 2 = Sliding  , 3 = Converging  */
    int atmvolmethod,      /* 1 = Lognormal  , 2 = Normal  , 3 = SigmaBeta */
    long observation_freq, long num_subdiscret_obsfreq, double shiftvol,
    long nx,       /* number of state variable steps in the Pat's integral */
    double maxstd, /* Maximum number of std between forward and strike */
    int nofunding,

    /* For period started in the past */
    int noDaysAccrued, int floatingFixing,

    /*	EOD Fixing Flag */
    int eod_fix_flag, /*	0: I  , 1: E */
    /*	EOD Payment Flag */
    int eod_pay_flag, /*	0: I  , 1: E */
    /*	EOD Exercise Flag */
    int eod_ex_flag, /*	0: I  , 1: E */

    /* Calibrated term structure */
    SrtLgmTSData *lgmTSData, /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr *atcTSData,    /* ptr to zeta/G data (NULL => not req'd) */

    /*	Results	*/
    double *price, double *intrinsic_value,
    double *
        *forwardTimeSwap, /*[0..1]*[0..nEx-1]  the first is the forward computed
                                                  today; the second is with the
                             call put parity */
    double *dvStrikes     /* [0..nEx-1]  the strikes used in the calibration */
    /*	double	*reval_times	*/
);

/* Create a Callable time swap structure  , and copy input data into it */
LGMErr LGMFillCallableTimeSwap2(
    SrtCallTimeSwap **dealPtrPtr,   /* output: option deal */
    SrtCallTimeSwap **dealTSPtrPtr, /* output: time swap deal */
    long *nExleftPtr,               /* output: effective number of exer */
    Date tfirst, char *ycName,      /* info about today */
    long n, Date *tStartCpn, Date *tPayCpn, /* coupon leg */
    double *tCpn, double *tCpn2, /* should be 0 unless it is a double range */
    double *tgear_floatdigital, double *tCvgCpn, double *tFundingPayment,
    long nEx, Date *tEx, Date *tExStart, double *exFee, /* exercise info */
    double call_spread_up, double call_spread_low, long observation_freq,
    double **barriers, SrtCompounding freqIndex, SrtCompounding freqAccrued,
    double CorrelStart, double CorrelEnd,
    SrtReceiverType payrec, /* PAYER or RECEIVER */
    int timeswapvolmethod, int atmvolmethod, int calibrationmethod,
    long num_subdiscret_obsfreq, double shiftvol, long nx, int nofunding,
    char *vol_id,
    LGMErr (*GetSABRParam)(long, long, char *, SABR_VOL_TYPE, double *),
    double exFeeTimeSwap, int eod_pay_flag, /*	0: I  , 1: E */
    int eod_fix_flag);                      /*	0: I  , 1: E */

/* ------------------------------------------------------------------------------------------------------------------
 */
/*
        New midat interface with fixing date in vol function
*/
/* ------------------------------------------------------------------------------------------------------------------
 */

LGMErr TestLGMautocalCallerFix(
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
    LGMErr (*GetVol)(Date, Date, Date, double, SRT_Boolean,
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
    char *szRefRate);

#endif
