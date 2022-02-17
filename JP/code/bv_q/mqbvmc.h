/******************************************************************************
 * Module:      
 * Submodule:
 * File:      mqbvmc.h  
 * Function:    
 * Author:        Interest Rates DR
 * Revision:  $Header$
 *****************************************************************************/
#ifndef  Q3BVMC_H
#define  Q3BVMC_H

#include "q3.h"

/*-----------------------------------------------------------------------------
 * Q3MQBivarMCPricer
 *
 * Multi-q Bivariate Option Pricer Add-in
 *
 */
int Q3MQBivarMCPricer(
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
    double *results          /* 37 (O) Pricing results.                     */ 
    );   


/*-----------------------------------------------------------------------------
 * Q3MQBivarMCInit
 *
 * Initialize payoff parameters and select payoff function
 * 
 */
int Q3MQBivarMCInit(
    long      optType,         /* (I)                           */
    long      numPayPrms,      /* (I)   number of payoff params */
    double   *payPrms,         /* (I)   payoff coefficients     */
    double   *corr,            /* (I)   gaussian correlation    */
    MQDATA  **mq,              /* (I/O) measure data            */
    PAYOFF   *pf,              /* (O)   option payoff structure */
    FPAYOFF **payFunc          /* (O)   payoff function         */
    );


/*-----------------------------------------------------------------------------
 * Q3MQBivarMCPricerSimple
 * 
 * Multi-q bivariate option pricer with simplified interface
 *
 */
int Q3MQBivarMCPricerSimple(
    double  expiry,          /*  1 (I) exipation in years                  */
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
 * Q3MQBivarMCMidCurvePricer
 *
 * Multi-q Bivariate MidCurve Swaption Pricer Add-in
 *
 */
int Q3MQBivarMCMidCurvePricer(
    long    today,           /*  (I) Vol base date                        */
    long    valueDate,       /*  (I) Value date                           */
    long    numIndxPts,      /*  (I) index curve num points               */
    long   *indxDates,       /*  (I) index curve dates                    */
    double *indxRates,       /*  (I) index curve rates                    */
    long    rateFreq0,       /*  (I) mid curve frequency                  */
    char   *rateDCC0,        /*  (I) mid curve day count convention       */
    long    rateFreq1,       /*  (I) 1st index curve frequency            */
    char   *rateDCC1,        /*  (I) 1st index curve day count convention */
    long    expiryDate1,     /*  (I) 1st index reset date                 */
    long    startDate1,      /*  (I) 1st index start date                 */
    long    matDate1,        /*  (I) 1st index rate end date              */
    double *smile1,          /*  (I) 1st index smile parameters           */
    double  sigATM1,         /*  (I) 1st index ATM vol                    */
    long    rateFreq2,       /*  (I) 2nd index curve frequency            */
    char   *rateDCC2,        /*  (I) 2nd index curve day count convention */
    long    expiryDate2,     /*  (I) 2nd index reset date                 */
    long    startDate2,      /*  (I) 2nd index start date                 */
    long    matDate2,        /*  (I) 2nd index rate end date              */
    double *smile2,          /*  (I) 2nd index smile parameters           */
    double  sigATM2,         /*  (I) 2nd index ATM vol                    */
    char   *fixedType,       /*  (I) Type of payment(N,E,I,J)             */
    char   *stubType,        /*  (I) Stub type(S,B,D,N)                   */
    char   *stubAtEnd,       /*  (I) (F)ront or (B)ack                    */
    long    numDiscPts,      /*  (I) Discount curve num points            */
    long   *discDates,       /*  (I) Discount curve dates                 */
    double *discRates,       /*  (I) Discount curve rates                 */
    long    numVNFMParams,   /*  (I) index num VNFM parameters            */
    double *vnfmParams,      /*  (I) index VNFM parameters                */
    double  corr,            /*  (I) Index correlation data               */
    long    optType,         /*  (I) Option type                          */
    long    payDate,         /*  (I) Option payment date                  */
    long    setlType,        /*  (I) Cash or physical settlement          */
    double  strike,          /*  (I) Option strike                        */
    double *results          /*  (O) Pricing results.                     */
    );

#endif /* define Q3BVMC_H */
