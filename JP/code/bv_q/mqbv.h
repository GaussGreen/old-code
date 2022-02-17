/******************************************************************************
 * Module:      
 * Submodule:
 * File:      mqbv.h    
 * Function:    
 * Author:        Interest Rates DR
 * Revision:  $Header$
 *****************************************************************************/
#ifndef  Q3BV_H
#define  Q3BV_H

#include "q3.h"

/*-----------------------------------------------------------------------------
 * Q3MQBivarPricer
 *
 * Multi-q bivariate option pricer.
 *
 */
int Q3MQBivarPricer(
    long    today,           /*  1 (I) Vol base date                        */
    long    valueDate,       /*  2 (I) Value date                           */
    long    numIndxPts1,     /*  3 (I) 1st index curve num points           */
    long   *indxDates1,      /*  4 (I) 1st index curve dates                */
    double *indxRates1,      /*  5 (I) 1st index curve rates                */
    long    rateFreq1,       /*  6 (I) 1st index curve frequency            */
    char   *rateDCC1,        /*  7 (I) 1st index curve day count convention */
    long    expiryDate1,     /*  8 (I) 1st index reset date                 */
    long    startDate1,      /*  9 (I) 1st index start date                 */
    long    matDate1,        /* 10 (I) 1st index rate end date              */
    double *smile1,          /* 11 (I) 1st index smile parameters           */
    double  sigATM1,         /* 12 (I) 1st index ATM vol                    */
    long    numVNFMParams1,  /* 13 (I) 1st index num VNFM parameters        */
    double *vnfmParams1,     /* 14 (I) 1st index VNFM parameters            */
    long    numIndxPts2,     /* 15 (I) 2nd index curve num points           */
    long   *indxDates2,      /* 16 (I) 2nd index curve dates                */
    double *indxRates2,      /* 17 (I) 2nd index curve rates                */
    long    rateFreq2,       /* 18 (I) 2nd index curve frequency            */
    char   *rateDCC2,        /* 19 (I) 2nd index curve day count convention */
    long    expiryDate2,     /* 20 (I) 2nd index reset date                 */
    long    startDate2,      /* 21 (I) 2nd index start date                 */
    long    matDate2,        /* 22 (I) 2nd index rate end date              */
    double *smile2,          /* 23 (I) 2nd index smile parameters           */
    double  sigATM2,         /* 24 (I) 2nd index ATM vol                    */
    long    numVNFMParams2,  /* 25 (I) 2nd index num VNFM parameters        */
    double *vnfmParams2,     /* 26 (I) 2nd index VNFM parameters            */
    long    numDiscPts,      /* 27 (I) Discount curve num points            */
    long   *discDates,       /* 28 (I) Discount curve dates                 */
    double *discRates,       /* 29 (I) Discount curve rates                 */
    double  corr,            /* 30 (I) Index correlation data               */
    double *fwdVolInput,     /* 31 (I) Forward vol betw rate start dates    */
    long    optType,         /* 32 (I) Option type                          */
    long    payDate,         /* 33 (I) Option payment date                  */
    long    setlType,        /* 34 (I) Cash or physical settlement          */
    long    numPayoffParams, /* 35 (I) Num payment parameters               */
    double *payoffParams,    /* 36 (I) Payment parameters                   */
    char   *holidayfile,     /* 37 (I) holiday file name                    */
    char   *BusVolDCC,       /* 38 (I) BUS/251F or BUS/BUS                  */
    double *results          /* 39 (O) Pricing results.                     */ 
    );


/*-----------------------------------------------------------------------------
 * Q3MQBivarPricerSimple
 * 
 * Multi-q bivariate option pricer with simplified interface
 *
 */
int Q3MQBivarPricerSimple(
    double  expiry,          /*  1 (I) expiration in years                 */
    long    rateType1,       /*  2 (I) rate or spread 1                    */
    double  fwdRate1,        /*  3 (I) fwd 1                               */
    double *smile1,          /*  4 (I) smile 1                             */
    double  sigATM1,         /*  5 (I) % vol or bp vol for fwd 1           */
    long    rateType2,       /*  6 (I) rate or spread 2                    */
    double  fwdRate2,        /*  7 (I) fwd 2                               */
    double *smile2,          /*  8 (I) smile 2                             */
    double  sigATM2,         /*  9 (I) % vol or bp vol for fwd 2           */
    double  corr,            /* 10 (I) correlation between rate drivers    */
    long    optType,         /* 11 (I) option type                         */
    long    numPayoffParams, /* 12 (I) number of payoff parameters         */
    double *payoffParams,    /* 13 (I) payoff parameters                   */
    double *results          /* 14 (O) pricing results                     */
    );

/*-----------------------------------------------------------------------------
 * Q3MQPSAPricer
 *
 * Multi-q bivariate PSA option pricer.
 *
 */
int Q3MQPSAPricer(
    long    today,           /*  1 (I) Vol base date                        */
    long    valueDate,       /*  2 (I) Value date                           */
    long    bsFixLegFreq,    /*  3 (I) BMA and CMS fixed leg frequency      */
    char   *bsFixLegDCC,     /*  4 (I) BMA and CMS fixed Leg DCC            */
    long    bsFltLegFreq,    /*  5 (I) BMA and CMS float leg frequency      */
    char   *bsFltLegDCC,     /*  6 (I) BMA and CMS float leg DCC            */
    long    numLiborZeroPts, /*  7 (I) LIBOR zero curve num points          */
    long   *liborZeroDates,  /*  8 (I) LIBOR zero curve dates               */
    double *liborZeroRates,  /*  9 (I) LIBOR zero curve rates               */
    long    numBsZeroPts,    /* 10 (I) basis zero curve num points          */
    long   *bsZeroDates,     /* 11 (I) basis zero curve dates               */
    double *bsZeroRates,     /* 12 (I) basis zero curve rates               */
    long    expiryDate,      /* 13 (I) reset date                           */
    long    startDate,       /* 14 (I) start date                           */
    long    matDate,         /* 15 (I) end date                             */
    double *cmsSmile,        /* 16 (I) CMS smile parameters                 */
    double  cmsSigATM,       /* 17 (I) CMS ATM vol                          */
    long    numVNFMParams,   /* 18 (I) 1st index num VNFM parameters        */
    double *vnfmParams,      /* 19 (I) 1st index VNFM parameters            */
    double *bsSpreadSmile,   /* 20 (I) spread smile parameters              */
    double  bsSpreadSigATM,  /* 21 (I) spread ATM vol                       */
    long    numDiscPts,      /* 22 (I) Discount curve num points            */
    long   *discDates,       /* 23 (I) Discount curve dates                 */
    double *discRates,       /* 24 (I) Discount curve rates                 */
    double  corr,            /* 25 (I) Index correlation data               */
    long    optType,         /* 26 (I) Option type                          */
    long    payDate,         /* 27 (I) Option payment date                  */
    long    setlType,        /* 28 (I) Cash or physical settlement          */
    long    numPayoffParams, /* 29 (I) Num payment parameters               */
    double *payoffParams,    /* 28 (I) Payment parameters                   */
    char   *holidayfile,     /* 29 (I) holiday file name                    */
    char   *BusVolDCC,       /* 30 (I) BUS/251F or BUS/BUS                  */
    double *results          /* 31 (O) Pricing results.                     */ 
    );


#define CORR_CALIB_RESN     1.0e-14  /*df/d\rho                            */
#define CORR_CALIB_TOL      1.0e-6   /*calibration relative error tolerance*/
#define CORR_CALIB_MAX_ITER 15       /*Max # iteration                     */
#define CORR_DELTA          1.0e-04  /*step size of FD to estimate gradient*/
/*-----------------------------------------------------------------------------
 * Q3SmileBivarSimpleCorrCalib
 * 
 * Multi-q bivariate correlation calibration with simplified interface
 *
 */
int Q3SmileBivarSimpleCorrCalib(
    double  premium,         /*  1 (I) premium                             */
    double  expiry,          /*  2 (I) expiration in years                 */
    double  fwdRate1,        /*  3 (I) fwd 1                               */
    double *smile1,          /*  4 (I) smile 1                             */
    double  sigATM1,         /*  5 (I) % vol or bp vol for fwd 1           */
    double  fwdRate2,        /*  6 (I) fwd 2                               */
    double *smile2,          /*  7 (I) smile 2                             */
    double  sigATM2,         /*  8 (I) % vol or bp vol for fwd 2           */
    long    optType,         /*  9 (I) option type                         */
    double  strike,          /* 10 (I) strike                              */
    double *payoffParams,    /* 11 (I) payoff parameters                   */
    double  initialGuess,    /* 12 (I) initial guess of corr               */
    double *results          /* 13 (O) pricing results                     */
    );

/*-----------------------------------------------------------------------------
 * Q3MQBivarCorr
 *
 * Multi-q bivariate correlation calculator
 *
 */
int Q3MQBivarCorr(
    long    today,           /*  1 (I) Vol base date                        */
    long    valueDate,       /*  2 (I) Value date                           */
    long    numIndxPts1,     /*  3 (I) 1st index curve num points           */
    long   *indxDates1,      /*  4 (I) 1st index curve dates                */
    double *indxRates1,      /*  5 (I) 1st index curve rates                */
    long    rateFreq1,       /*  6 (I) 1st index curve frequency            */
    char   *rateDCC1,        /*  7 (I) 1st index curve day count convention */
    long    expiryDate1,     /*  8 (I) 1st index reset date                 */
    long    startDate1,      /*  9 (I) 1st index start date                 */
    long    matDate1,        /* 10 (I) 1st index rate end date              */
    double *smile1,          /* 11 (I) 1st index smile parameters           */
    double  sigATM1,         /* 12 (I) 1st index ATM vol                    */
    long    numVNFMParams1,  /* 13 (I) 1st index num VNFM parameters        */
    double *vnfmParams1,     /* 14 (I) 1st index VNFM parameters            */
    long    numIndxPts2,     /* 15 (I) 2nd index curve num points           */
    long   *indxDates2,      /* 16 (I) 2nd index curve dates                */
    double *indxRates2,      /* 17 (I) 2nd index curve rates                */
    long    rateFreq2,       /* 18 (I) 2nd index curve frequency            */
    char   *rateDCC2,        /* 19 (I) 2nd index curve day count convention */
    long    expiryDate2,     /* 20 (I) 2nd index reset date                 */
    long    startDate2,      /* 21 (I) 2nd index start date                 */
    long    matDate2,        /* 22 (I) 2nd index rate end date              */
    double *smile2,          /* 23 (I) 2nd index smile parameters           */
    double  sigATM2,         /* 24 (I) 2nd index ATM vol                    */
    long    numVNFMParams2,  /* 25 (I) 2nd index num VNFM parameters        */
    double *vnfmParams2,     /* 26 (I) 2nd index VNFM parameters            */
    long    numDiscPts,      /* 27 (I) Discount curve num points            */
    long   *discDates,       /* 28 (I) Discount curve dates                 */
    double *discRates,       /* 29 (I) Discount curve rates                 */
    double *corr,            /* 30 (I) Index correlation data               */
    double *fwdVolInput,     /* 31 (I) Forward vol betw rate start dates    */
    long    payDate,         /* 32 (I) Option payment date                  */
    long    setlType,        /* 33 (I) Cash or physical settlement          */
    char   *holidayfile,     /* 34 (I) holiday file name                    */
    char   *BusVolDCC,       /* 35 (I) BUS/251F or BUS/BUS                  */
    double *results          /* 36 (O) Correlation                          */ 
    );

/*-----------------------------------------------------------------------------
 * Q3SwapRateCalc2
 *
 * Swap rate, annuity, etc calculator from curve (for quasi-vanilla 1d)
 *
 */
int Q3SwapRateCalc2 (
    long              VolBaseDate,       /* 1  (I) ATM vol base date       */
    long              SmileBaseDate,     /* 2  (I) Smile base date         */
    long              CurveBaseDate,     /* 3  (I) Index curve base date   */
    long              numIndxPts,        /* 4  (I) Nb of index points      */
    long              *indxDates,        /* 5  (I) Index dates             */
    double            *indxRates,        /* 6  (I) Act/365F index rates    */
    long              numDiscPts,        /* 7  (I) Nb of discount points   */ 
    long              *discDates,        /* 8  (I) Discount dates          */
    double            *discRates,        /* 9  (I) Act/365F discount rates */
    long              rateFreq,          /* 10 (I) rate frequency          */
    char              *rateDCC,          /* 11 (I) rate day count conv     */
    long              expiryDate,        /* 12 (I) rate reset date         */
    long              startDate,         /* 13 (I) rate effective date     */
    long              matDate,           /* 14 (I) rate maturity date      */
    long              payDate,           /* 15 (I) option payment date     */
    char              *HolidayVol,       /* 16 (I) holiday file name       */
    char              *BusVolDCC,        /* 17 (I) bus day disc convention */
    double            *output            /* 18 (O) fwd, annuity etc        */
    );

/*-----------------------------------------------------------------------------
 * Q3BasisLegCalc
 *
 * Basis leg calculator from curve (for quasi-vanilla 1d)
 *
 */
int Q3BasisLegCalc (
    long              CurveBaseDate,     /* 1  (I) Index curve base date   */
    long              numIndxPts,        /* 2  (I) Nb of index points      */
    long              *indxDates,        /* 3  (I) Basis index dates       */
    double            *indxRates,        /* 4  (I) Act/365F index rates    */
    long              numDiscPts,        /* 5  (I) Nb of discount points   */ 
    long              *discDates,        /* 6  (I) Discount dates          */
    double            *discRates,        /* 7  (I) Act/365F discount rates */
    long              fixLegFreq,        /* 8  (I) BMA and CMS fixed leg frequency      */
    char              *fixLegDCC,        /* 9  (I) BMA and CMS fixed Leg DCC            */
    long              rateFreq,          /* 10 (I) rate frequency          */
    char              *rateDCC,          /* 11 (I) rate day count conv     */
    long              expiryDate,        /* 12 (I) rate reset date         */
    long              startDate,         /* 13 (I) rate effective date     */
    long              matDate,           /* 14 (I) rate maturity date      */
    long              payDate,           /* 15 (I) rate maturity date      */
    char              *Holiday,          /* 16 (I) holiday file name       */
    double            *output            /* 17 (O) fwd, annuity etc        */
    );

/*f----------------------------------------------------------------------------
 * Q3SwapRateCalc3
 *
 * Swap rate, annuity, etc calculator from curve (for mid-curve products)
 */
int Q3SwapRateCalc3 (
    long              VolBaseDate,       /* 1  (I) ATM vol base date       */
    long              SmileBaseDate,     /* 2  (I) Smile base date         */
    long              CurveBaseDate,     /* 3  (I) Index curve base date   */
    long              numIndxPts,        /* 4  (I) Nb of index points      */
    long              *indxDates,        /* 5  (I) Index dates             */
    double            *indxRates,        /* 6  (I) Act/365F index rates    */
    long              numDiscPts,        /* 7  (I) Nb of discount points   */
    long              *discDates,        /* 8  (I) Discount dates          */
    double            *discRates,        /* 9  (I) Act/365F discount rates */
    long              rateFreq,          /* 10 (I) rate frequency          */
    char              *rateDCC,          /* 11 (I) rate day count conv     */
    char              *stubType,         /* 12 (I) stub type               */
    char              *stubAtEnd,        /* 13 (I) (F) or (B)              */
    long              expiryDate,        /* 14 (I) rate reset date         */
    long              startDate,         /* 15 (I) rate effective date     */
    long              matDate,           /* 16 (I) rate maturity date      */
    long              payDate,           /* 17 (I) option payment date     */
    double            *output            /* 18 (O) fwd, annuity etc        */
    );

/*-----------------------------------------------------------------------------
 * Q3MQBivarInit           
 *
 * Initialize payoff parameters and select payoff function
 *
 */
int Q3MQBivarInit(
    long      optType,
    long      numPayoffParams,
    double   *payoffParams,
    double   *corr,
    MQDATA   **mq,
    PAYOFF   *pf, 
    FPAYOFF **payFunc
    );


/*-----------------------------------------------------------------------------
 * Q3MQCalibJointFwd
 *
 *
 */
int Q3MQCalibJointFwd(
    double   target,
    double  *corr,
    long     index,
    MQDATA **mq
    );


/*-----------------------------------------------------------------------------
 * Q3MQCalibJoint
 *
 *
 */
int Q3MQCalibJoint(
    double  *corr,
    long     index,
    MQDATA **mq
    );


/*----------------------------------------------------------------------------
 * Q3QuasiErrBuffInit
 */
void Q3QuasiErrBuffInit(void);

int    Q3QuasiErrBuffRetrieve        (long    outBufLen,   
                                      char    *outBuffer);

int    Q3QuasiErrBuffLength          ();

/* copy from vanl.h */
int Q3AdjustFreqDCCATMVol_BVQ(long              CurveBaseDate,
                             long              numDiscPts,
                             long              *discDates,
                             double            *discRates,
                             long              startDate,
                             long              matDate,
                             long              fixedFreq1,
                             long              fixedFreq2,
                             char              *fixedType,
                             char              *fixedDCC1,
                             char              *fixedDCC2,
                             char              *stubType,
                             char              *stubAtEnd,
                             double            *output);

#endif /* define Q3BV_H */
