#ifndef _SRT_H_LGMPROTOS_H
#define _SRT_H_LGMPROTOS_H

#include "num_h_maxmin.h"

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

/* calibrate and evaluate an American swaption */
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
    int strikechoice,    /* 1 = IRR  , 2 = dIRR */
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
    SrtLgmTSData *lgmTSData,  /* ptr to tau/sigma data (NULL => not req'd) */
                              /* miscellaneous */
    double *intrinsicValPtr); /* ptr to intrinsic value (NULL => not req'd) */
/* calibrate and evaluate an American swaption */

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
    int strikechoice,    /* 1 = IRR  , 2 = dIRR */
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
    char *szRefRate);        /* optional refrate included for MUNI pricing */
/* calibrate and evaluate mid-atlantic swaption */

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
    double **ForwardVols,

    /* outputs */
    String outfile,    /* output file name for log file (unused) */
    double *LGMvalPtr, /* LGM value of mid-atlantic */
    SrtLgmExerBdryData *lgmExerBdryData, /* ptr to exercise boundary data
                                            structure (NULL => not req'd) */
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
                                         /* calibrated term structure */
    SrtLgmTSData *lgmTSData,  /* ptr to tau/sigma data (NULL => not req'd) */
    LGM_TSPtr *atcTSData,     /* ptr to zeta/G data (NULL => not req'd) */
                              /* miscellaneous */
    double *intrinsicValPtr); /* ptr to intrinsic value (NULL => not req'd) */

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
    double **ForwardVols,

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

- the exotic side also includes a cap paying gear * max (0  , Libor CASH -
cap_str) * libor_cvg */

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
    FWD_VOL_STR *fwdVolStr); /* ptr to fwd vol data (NULL => not req'd) */
/* calibrate and evaluate an Bermudan Cap floater */

/* calibrate and evaluate an Bermudan Cap floater */
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
    SrtLgmTSData *lgmTSData); /* ptr to tau/sigma data (NULL => not req'd) */

/******************************************************************/
/* Routines that help make the conversion from old LGM call to new LGM call */
/* Copies the input data into the Simple MidAtlantic Structure pointed to by
 * MidAtPtr */
/* This routine ASSUMES that the option holder will receive all payments
whose pay dates tPay[i] are strictly after the settlement date */

LGMErr LGMFillSimMidAt(Date tfirst, SrtSimMidAt **MidAtPtrPtr, long nPay,
                       double *Payment, Date *tPay, long nEx, Date *tEx,
                       Date *tStart, double *Strike, double *RedFirstPay,
                       SrtReceiverType payrec);

/* Create a simple American deal structure  , and copy input data into it */
LGMErr LGMFillSimAmer(Date tfirst, int standardLag, SrtSimAmer **AmerPtrPtr,
                      long nfix, Date *tfixStart, Date *tfixEnd, Date *tfixPay,
                      double *fixFullPayment, SrtBasisCode fixBasisCode,
                      int fixEarlyFlag, double *Strike, long nfltdates,
                      Date *tflt, SrtBasisCode fltBasisCode, int fltEarlyFlag,
                      Date tFirstExer, int lagExerStart, int calbus,
                      BusDayConv convStart, SrtReceiverType payrec);

/* Create a Bermudan inverse floater structure  , and copy input data into it */
LGMErr
LGMFillCallInvFloater(SrtCallInvFlt **dealPtrPtr, /* output: the deal */
                      long *nExleftPtr, /* output: effective number of exer */
                      Date tfirst, char *ycName,        /* info about today */
                      long n, Date *tStart, Date *tPay, /* coupon leg */
                      double *a, double *gear, double *cap_str, double *cvg,
                      double *lcvg, /* more coupon leg */
                      long nEx, Date *tEx, Date *tExStart,
                      double *exFee,           /* exercise info */
                      SrtReceiverType payrec); /* PAYER or RECEIVER */

/* Create a Bermudan cap floater structure  , and copy input data into it */
LGMErr LGMFillCallCapFloater(
    long *nExleftPtr,                            /* effective number of exer */
    Date tfirst, char *ycName, double rfund_cur, /* info about today */
    SrtCallCapFlt **dealPtrPtr,                  /* output: the deal */
    long n, Date *t, double *amax, double *amin, /* coupon leg */
    double *marg, double lvg, SrtBasisCode cpnBasis, /* more coupon leg */
    SrtBasisCode aBasis,                             /* more more coupon leg */
    long m, Date *tflt, SrtBasisCode fltBasis,       /* floating leg */
    long nEx, Date *tEx, double *exFee,              /* exercise info */
    SrtReceiverType payrec,                          /* PAYER or RECEIVER */
    int earlyFlag, int resetFlt);                    /* early and reset flags */

/*******************************************************************/
/* Codes to provide beta for lognormal (beta=1) or normal (beta=0) */
/*	This code always returns the value 1 */
LGMErr LGMReturnExp1(Date swapstart, Date swapend, double *beta);

/*******************************************************************/
/* Codes to provide beta for lognormal (beta=1) or normal (beta=0) */
/*	This code always returns the value 0 */
LGMErr LGMReturnExp0(Date swapstart, Date swapend, double *beta);

/******************************************************************/
/* This is a dummy code that does nothing.
It is called when and only when the code is going to use
the swaption or caplet for calibration. It should be replaced by a
code that records the swaptions used for calibration so we can
compute vega risks efficiently */
LGMErr LGMRecVolDummy(Date swapstart, Date swapend, double Rfix);

/* This code sets the calibration methodology to the default method */
LGMCalParmPtr LGMSetCalibMeth(int iLGMOneTwoFactor, int FixedTauFlag,
                              int CapsFlag, double Tau, double Alpha,
                              double Gamma, double Rho, double MaxStd,
                              long *Zeta1Dates, double *StartZeta1s,
                              long *TauDates, double *StartTaus,
                              double **HybridShortInstrsIndex);

/* This code sets the Convolution parameters to their default values */
ConvParamsPtr LGMSetDefaultEvalParms(void);

/******************************************************************/
/******************************************************************/
/* LGMNewAutocal.c routines */
/******************************************************************/
/******************************************************************/

/* MAIN PROGRAM to calibrate zeta(t) & G(t) and evaluate an American */
LGMErr LGMNewAmerican(
    /* task list */
    int skipEval,     /* 0=calibrate & value deal  , 1=calibrate only */
    int skipCalib,    /* 0=calibrate & value deal  , 1=value deal only */
    int convertTS,    /* 1=compute & output sig-kappa ts equivalent to LGM ts */
    int findExerBdry, /* 1=find swap rates at exer boundary  */
                      /* information about today & today's yield curve */
    Date tNow,        /* eval as if today is tNow */
    int eod,          /* 0=can exercise on tNow  , 1=cannot exercise on tNow */
    String ycName,    /* market pointer for discount factors */
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
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
    LGM_TS **LGMtsPtrPtr,                /* calibrated zeta-G term structure */
    SigKapTS **SigKaptsPtrPtr, /* an equivalent sigma-kappa term structure */
    String outfile);           /* output file name for log file (unused) */

/* MAIN PROGRAM to calibrate zeta(t) & G(t) to the market and evaluate a
mid_atlantic */
LGMErr LGMnewautocal(
    /* operations to be done */
    int skip_deal_eval,   /* 0=value deal  , 1=calibrate only */
    int skip_calibration, /* 0=calibrate & value deal  , 1=value deal only */
    int convert_ts_flag,  /* 1=compute new sigs and taus & store; 0=don't bother
                           */
    int find_exer_bdry, /* 1=find swap rates at exer boundary; 0=don't bother */

    /* information about today & today's market prices */
    Date tNow,     /* eval as if today is tNow */
    int eod,       /* 0=can exercise on tNow  , 1=cannot exercise on tNow */
    String ycname, /* yield curve name for discount factors */
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
    void *dealptr,        /* pointer to the deal */

    /* calibration and eval parameters */
    LGMCalParm *CalReq, /* structure defining calibration method to be used */
    ConvParams
        *EvalParms, /* structure containing numerical evaluation parameters */

    /* outputs */
    double *LGMvalue,  /* value of the deal */
    double *intrinval, /* intrinsic value of the deal */

    LGMSwptns *
        *ExerIntoPtrPtr, /* European swaptions most similar to underlying */
    LGMSwptns **
        ExerBdryPtrPtr, /* European swaptions struck at the exercise boundary */
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
    LGM_TS **LGMtsPtrPtr,                /* calibrated zeta-G term structure */
    SigKapTS **SigKaptsPtr, /* an equivalent sigma-kappa term structure */

    String outfile); /* output file name for log file (unused) */

/******************************************************************/
/* This routine creates a Mid Atlantic with k months between exercises out of
a Simple American deal */
SrtSimMidAtPtr LGMCreateMidAtFromAmer(Date tfirst, double rfltCur,
                                      String ycName, LGMMarkConv *conv,
                                      int dmos, SrtSimAmer *AmerPtr);

/******************************************************************/
/* Get info for calibration:
The effective number of exercise dates and the exercise dates - nEx  , TauArr[]
The PV of fixed leg/PV of strike at each exercise date - FVArr[]
The last pay date of the deal - tlast

Compute the intrinsic value of the deal  , and store the underlying
exercise dates  , enddates  , and the fixed rate Rfix of the vanilla swaption
that comes closest to describing the underlying in UnderSwptns structure  ,
and return these to the calling routine as "ExerInto"
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
    LGMDealType dealtype,     /* type of deal */
    void *dealPtr,            /* pointer to the deal */
                              /* output */
    long *EffnEx,             /* Effective number of exercise dates */
    Date *LasttPay,           /* Last pay date of deal */
    Date **TauArrPtr,         /* Array of effective exercise dates */
    double **FVArrPtr,        /* PV of fixed leg/PV of strike for these dates */
    double *intValPtr,        /* intrinsic value of the deal */
    LGMSwptns **ExerIntoPtr); /* closest European swaptions underlying deal */

/******************************************************************/
/******************************************************************/
/* LGMCalib.c routines */
/******************************************************************/
/******************************************************************/
/* Main LGM calibration routine */
LGMErr LGMCalibration(
    /* info about today's market */
    Date tNow,     /* first possible exercise date (today+eod) */
    String ycname, /* yield curve name for discount factors */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to provide swaption/cap vols */
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* function to provide exponent (beta) to interpret vols */
    LGMErr (*RecVol)(
        Date, Date,
        double), /* function called whenever LGM uses a vol for pricing */
                 /* info about the deal */
    LGMDealType dealType, void *dealptr,
    long nEx,      /* number of exercise dates to calibrate to */
    Date tlast,    /* last date that needs to be calibrated */
    Date *TauArr,  /* calibrate on these exercise dates */
    double *FVArr, /* use swaptions with these forward vals of fixed legs */
                   /* calibration methodology */
    LGMCalParm *CalReq,   /* structure defining calibration method to be used */
                          /* calibrated term structure */
    LGM_TS **LGMtsPtrPtr, /* calibrated zeta-G term structure */
    LGMCalSet **RefDealPtrPtr, /* Reference instruments used for calibration */
    SrtLgmRefSwptnData *lgmRefSwptnData); /* ptr to reference swaption data
                                             structure (NULL => not req'd) */

/******************************************************************/
/* Main routine to construct the reference instruments */
LGMErr LGMMakeRefDeals(
    LGMDealType DealType, void *dealptr, Date tNow, /* EVALUTAION DATE */
    Date tLast,   /* LAST PAY DATE OF THE REFERENCE INSTRUMENTS */
    long nArr,    /* NUMBER OF EXERCISE DATES IN DEAL*/
    Date *TauArr, /* [*  ,1  ,2  ,...  ,nEx] ARRAY OF EXERCISE DATES AFTER TODAY
                   */
    double *RatArr, /* [*  ,1  ,2  ,...  ,nEx] RATIO OF PV OF FIXED LEG TO PV OF
                       STRIKE */
    String ycname,  /* DISCOUNT CURVE */
    LGMMarkConv *conv,  /* MARKET PLACE CONVENTIONS */
    LGMCalParm *CalReq, /* METHODS AND DATA FOR CHOOSING STRIKES */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* SWAPTION CAP/VOLS*/
    LGMErr (*GetBeta)(
        Date, Date,
        double *), /* EXPONENT BETA TO INTERPRET VOLto interpret vols */
    LGMCalSet **RefDealsPtrPtr); /* REFERENCE INSTRUMENT (OUPUT) */

/******************************************************************/
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
    SrtLgmRefSwptnData *lgmRefSwptnData); /* ptr to reference swaption data
                                             structure (NULL => not req'd) */

/******************************************************************/
/* Returns the CEV model's price of a swaption */
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
    double beta);

/******************************************************************/
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
    double beta);  /* CEV exponent */

/******************************************************************/
/* Creates a term structure using a constant kappa  , and calibrate it to
the set of caplets (if usecaps=1) or long swaptions (if usecaps!=1) */
LGMErr
LGMCalFixKap(LGM_TS **LGMtsPtrPtr, /* Return: Calibrated term structure */
             double kap,           /* Fixed kappa to use to construct G */
             int usestarts,        /* construct G on start dates or pay dates */
             long nDate,           /* number of dates for G if usestarts=0 */
             Date *GDate,          /* dates for G if usestarts=0 */
             int usecaps,          /* calibrate on caplets or swaptions */
             LGMCalSet *CSPtr);    /* Reference caplets */

/******************************************************************/
/* Creates a term structure using a given G(t) data  , and calibrate zeta
to the set of caplets (if usecaps=1) or long swaptions (if usecaps!=1) */
LGMErr
LGMCalGivenG(LGM_TS **LGMtsPtrPtr, /* Return: Calibrated term structure */
             long nDate,           /* number of dates for G if usestarts=0 */
             Date *GDate,          /* dates for G(t) */
             double *G,            /* values for G(t) */
             int usecaps,          /* calibrate on caplets or swaptions */
             LGMCalSet *CSPtr);

/******************************************************************/
/* Creates a term structure using a constant sigma  , and calibrates G(t)
to the set of caplets if (usecaps=1) or long swaptions (if usecaps!=1) */
LGMErr LGMCalFixSig(
    LGM_TS **LGMtsPtrPtr, /* Return: Calibrated term structure */
    int usecaps,       /* calibrate on caplets or long k into n-k swaptions */
    LGMCalSet *CSPtr); /* Reference instruments */

/******************************************************************/
/* Creates a term structure using a given zeta(t)  , and calibrates G(t)
to the set of caplets (if usecaps=1) or long swaptions (if usecaps!=1) */
LGMErr LGMCalGivenZeta(
    LGM_TS **LGMtsPtrPtr, /* Return: Calibrated term structure */
    long numZ,            /* number of zeta(t) values */
    Date *Zdate,          /* dates for zeta(t) */
    double *Zval,         /* zeta(t) values */
    int usecaps,       /* calibrate on caplets or long k into n-k swaptions */
    LGMCalSet *CSPtr); /* Reference caplets */

/******************************************************************/
/* Creates a term structure  , and finds G(t) by calibrating on the 1 into k
swaptions  , and finds zeta(t) by calibrating to the caplets (if usecaps=1)
or long swaptions (if usecaps!=1) */
LGMErr LGMCalFixExp(
    LGM_TS **LGMtsPtrPtr, /* Return: Calibrated term structure */
    int usecaps,       /* calibrate on caplets or long "k into n-k" swaptions */
    int keep,          /* keep 1 into k zeta */
    LGMCalSet *CSPtr); /* Reference caplets */

/******************************************************************/
/* Creates a term structure  , and finds zeta(t) and G(t) by simultaneously
calibrating on the long (k into n-k) swaptions and either the caplets (if
usecaps==1) or the short (k into 1) swaptions (if usecaps!=1) */
LGMErr LGMCalTenorDiag(
    LGM_TS **LGMtsPtrPtr, /* Return: Calibrated term structure */
    int usecaps,      /* calibrate on caplets or long "k into n-k" swaptions */
    LGMCalSet *CSPtr, /* Reference caplets */
    long *noPairsFlag); /* signals whether there any pairs to work with */

/******************************************************************/
/******************************************************************/
/* LGM1Deval.c routines */
/******************************************************************/
/******************************************************************/
/* Main LGM deal evaluator
Evaluates a deal according to the calibrated LGM model termstruct*/
LGMErr NewMidAtEval(
    ConvParams *parms, /* convolution numerical constants */
                       /* info about today */
    Date tNow,         /* calculation date */
    int endofday,      /* 1 = cannot exercise today  , 0 = can exercise today */
    String ycname,     /* yield curve name for discount factors */
    LGMCalParm *CalReq, LGM_TS *tsPtr, /* calibrated term structure */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *),              /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double *), /* swaption/cap exponents (beta) */
                                             /* information aobut the deal */
    LGMDealType dealtype,                    /* type of deal */
    void *dealPtr,                           /* definition of the deal */
                                             /* output */
    double *LGMvalue,                        /* value of the deal */
    long *nExBdry,     /* number of points in the exercise boundary */
    Date **tExBdry,    /* number of dates in the exercise boundary */
    double **xExBdry); /* array of exercise points (in x) */

/* Computes a fra in LGM autocal */
LGMErr fra_in_LGMAutocal(String ycName, LGM_TS *tsPtr, long value_date,
                         long tNow, long start, long end,
                         double xred, /* state variable in LGM autocal */
                         double *result);

/******************************************************************/
/* LGM replacement for tree */
LGMErr Convolver(
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

/******************************************************************/
/******************************************************************/
/* LGMUtil.c routines */
/******************************************************************/
/******* Utilities for Deal structures ****************************/

/* Create an array of (exer date  , end date  , fixed rate) for
conveniently passing exer boundary info */
LGMSwptnsPtr LGMCreateSwptns(long nopt);

/* Free a LGMSwptn structure */
void LGMFreeSwptns(LGMSwptns **ptr);

/* Free contents of SrtLgmRefSwptnData structure */
void LGMFreeRefSwptnData(SrtLgmRefSwptnData *ptr);

/* Create a General European stucture with nPay payments */
SrtGenEurPtr LGMCreateGenEur(long nPay);

/* Frees the general European structure pointed to by GenEurPtr */
void LGMFreeGenEur(SrtGenEur **ptr);

/* Check validity of general European deal */
LGMErr LGMValidGenEur(SrtGenEur *ptr, Date tfirst, long *nExleft);

/* Create a simple MidAtlantic stucture with nEx exercise dates and nPay
 * payments */
SrtSimMidAtPtr LGMCreateSimMidAt(long nEx, long nPay);

/* Frees Simple MidAtlantic structure */
void LGMFreeSimMidAt(SrtSimMidAt **ptr);

/* Check validity of Simple MidAtlantic deal */
LGMErr LGMValidSimMidAt(SrtSimMidAt *ptr, Date tfirst, long *nExleft);

/* Create a general MidAtlantic stucture with nEx exercise */
SrtGenMidAtPtr LGMCreateGenMidAt(long nEx, long *nPay);

/* Frees the General MidAtlantic structure */
void LGMFreeGenMidAt(SrtGenMidAt **ptr);

/* Check validity of General MidAtlantic deal */
LGMErr LGMValidGenMidAt(SrtGenMidAt *ptr, Date tfirst, long *nExleft);

/* Create simple American with nfix fixed and nflt floating periods */
SrtSimAmerPtr LGMCreateSimAmer(long nfix, long nflt);

/* Frees Simple American structure */
void LGMFreeSimAmer(SrtSimAmer **ptr);

/* Check validity of Simple American deal */
LGMErr LGMValidSimAmer(SrtSimAmer *ptr, Date tfirst);

/* Create a Bermudan inverse floater structure */
SrtCallInvFltPtr LGMCreateCallInvFlt(long nEx, long n);

/* Frees a Bermudan inverse floater structure */
void LGMFreeCallInvFlt(SrtCallInvFlt **ptr);

/* Check validity of a Bermudan inverse floater */
LGMErr LGMValidCallInvFlt(SrtCallInvFlt *ptr, Date tfirst, long *ExFirst);

/* Create a Bermudan cap floater structure */
SrtCallCapFltPtr LGMCreateCallCapFlt(long nEx, long n);

/* Frees a Bermudan inverse floater structure */
void LGMFreeCallCapFlt(SrtCallCapFlt **ptr);

/* Check validity of a Bermudan inverse floater */
LGMErr LGMValidCallCapFlt(SrtCallCapFlt *ptr, Date tfirst, long *ExFirst);

/************* Utilities for term structures
 * ***********************************/
/*Get zeta at thedate from term structure */
double LGMZetaFromTS(Date theDate, Date tNow, LGM_TS *tsPtr);

/* Get G at thedate from term structure */
double LGMGFromTS(Date theDate, LGM_TS *tsPtr);

/* Get sigma at thedate from term structure */
double LGMSigFromTS(Date theDate, SigKapTS *tsPtr);

/* Get kappa at the date from term structure */
double LGMKapFromTS(Date theDate, SigKapTS *tsPtr);

/* Routine to validate and re-scale zetaG term structure */
LGMErr LGMVerifyLGM_TS(LGM_TS *tsPtr);

LGMErr LGMVerifyLGM_TS_For_DiagSwpts(LGMCalSet *CSPtr, LGM_TS *tsPtr);

/* Routine to verify that sigma-kappa term structure is sensible */
LGMErr LGMVerifySigKapTS(SigKapTS *tsPtr);

/*  Create LGM zeta-G term structure */

LGM_TSPtr LGMCreateLGM_TS(long numZ, long NumG);
LGM_TSPtr LGMCreateLGM2F_TS(long numZ, long NumG);

void LGMFreeLGM_TS(LGM_TS **LGMtsPtr);
void LGMFreeLGM2F_TS(LGM_TS **LGMtsPtr);

/* Create SigmaKappa term structure */
SigKapTSPtr LGMCreateSigKapTS(long NumS, long NumK);

/* free sigma kappa term struture */
void LGMFreeSigKapTS(SigKapTS **ptr);

/* This routine creates a LGM zeta G term structure equivalent
to the sigma kappa term structure */
LGM_TSPtr LGMConvertSigKaptoZG(Date tNow, SigKapTS *SKtsPtr);

/* This routine creates a sigma kappa term structure equivalent
to the zeta G term structure */
SigKapTSPtr LGMConvertZGtoSigKap(Date tNow, LGM_TS *LGMtsPtr);

/******* Create & Free LGM Calibration  and Calibration Parameters ****/
/* Create Calibration Set for reference deals */
LGMCalSetPtr LGMCreateCalSet(long nPay, long nEx, long nLast);

/* Create Calibration Set for reference deals */
void LGMFreeCalSet(LGMCalSet **ptr);

/* Create Calibration Set for reference deals */
LGMCalParmPtr LGMCreateCalParm(long numG, long numZ, long numR1, long numR2);

/* Create Calibration Set for reference deals */
void LGMFreeCalParm(LGMCalParm **ptr);

/******* Swap/Cap Utilities***************************************/
/* Get market conventions from currency */
LGMErr LGMCcyDefaults(SrtCurvePtr yldcrv, LGMMarkConv *conventions);
LGMErr LGMMuniDefaults(char *szRefRate, LGMMarkConv *conventions);

/* Construct the pay dates for a standard fixed leg with exercise
date tEx and last pay date tEnd  , using standard market conventions */
Date *LGMFixLegSched(Date tEx, Date tEnd, long extraperiods, LGMMarkConv *conv,
                     long *nPayPtr);

/* Construct the coverage array for a standard fixed leg with dates in tPay[] */
double *LGMFixLegCvg(long nPay, Date *tPay, LGMMarkConv *conv);

/* Construct the discount factor array for a fixed leg  with dates in tPay[] */
double *LGMFixLegDF(long nPay, Date *tPay, Date tNow, String ycName);

/* Construct the G(tPay) array for a fixed leg with dates in tPay */
double *LGMFixLegGs(long nPay, Date *tPay, LGM_TS *tsPtr);

/* Converts forward value of fixed leg into equivalent fixed rate Rfix */
LGMErr LGMGetFixRate(Date tEx, Date tEnd,
                     double FwdVal, /* definition of swaption */
                     String ycName, /* yield curve name */
                     double *Rfix);

/* Get fixed swap rate equivalent to state variable x */
LGMErr LGMGetRfixFromX(Date tEx, Date tEnd, double x, Date tNow, String ycname,
                       LGM_TS *tsPtr, double *Rfix);

LGMErr LGMGetFloatRateFromX(Date tStart, Date tPay, double x, String ycName,
                            LGM_TS *tsPtr, double *FloatRate);

/* LGM value of a receiver swaption */
double LGMRecVal(double *ystar, long n, double *a, double DStart, double sqzeta,
                 double *Gpay, double Gst);

/* finds ystar  , the state variable at which a swaption is ATM */
double LGMFindystar(long n, double *a, double DStart, double zeta, double *Gpay,
                    double Gst);

/******* Math Utilities ******************************************/
/* calculates Gaussian probability density  , protected against underflows */
double LGMsafeGauss(double z);

/* calculates cumulative normal distribution  , protected against underflows */
double LGMsafeNorm(double z);

/* Gets the equivalent normal vol from CEV vol */
double LGMNormalVolFromCEV(double fwd, double strike, double mat, double CEVvol,
                           double beta);

/* Returns the price of a European option for a normal model */
double LGMNormOptPrice(double fwd, double strike, double NormVol, double mat,
                       double dftoStart, SrtCallPutType callput);

/* Computes the implied Normal Vol for a European option */
double LGMNormImpVol(double price, double fwd, double strike, double mat,
                     double dftoStart, SrtCallPutType callput);

/* Computes the price of a European option from Black's formula	*/
double LGMBlackOptPrice(double fwd, double strike, double BlackVol, double mat,
                        double dftoStart, SrtCallPutType callput);

/* Computes the implied Black Vol for a European option */
double LGMBlackImpVol(double price, double fwd, double strike, double mat,
                      double dftoStart, SrtCallPutType callput);

int BasicDichotomie(Date *t, long Indexmin, long Indexmax, Date tRef);

Err update_autocal_midat_struct(SrtSimMidAt *MidAt, char *YcName,
                                Err (*GetVol)(Date, Date, double, SRT_Boolean,
                                              double *),
                                Err (*GetBeta)(Date, Date, double *));

#endif
