/******************************************************************************
 * Module:       
 * Submodule:
 * File:      marginalDist.h   
 * Function:  Addin interface.  
 * Author:    Interest Rates DR
 * Revision:  $Header$
 *****************************************************************************/
#ifndef MQBVDIST_H
#define MQBVDIST_H

#include "gtomat.h"
#include "q3.h"

/*-----------------------------------------------------------------------------
 * Q3MQFAMarginalDistL
 *
 * Multi-q marginal distributions in FA.
 *
 */
int Q3MQFAMarginalDistL(
    long   *todayL,           /*  1 (I) Vol base date                     */
    long   *valueDateL,       /*  2 (I) Value date                        */
    long   *indxDates1L,      /*  3 (I) 1st index curve date              */
    double *indxRates1L,      /*  4 (I) 1st index curve rate              */
    long   *rateFreq1L,       /*  5 (I) 1st index frequency               */
    char   *rateDCC1L,        /*  6 (I) 1st index day count convention    */
    long   *rateRstDts1L,     /*  7 {I} 1st index expiry,start,mat. dates */
    double *smile1L,          /*  8 (I) 1st index smile parameters.       */
    double *sigATM1L,         /*  9 (I) 1st index ATM vol.                */
    double *vnfmParams1L,     /* 10 (I) 1st index VNFM parameters.        */
    long   *indxDates2L,      /* 11 (I) 2nd index curve date              */
    double *indxRates2L,      /* 12 (I) 2nd index curve rate              */
    long   *rateFreq2L,       /* 13 (I) 2nd index frequency               */
    char   *rateDCC2L,        /* 14 (I) 2nd index day count convention    */
    long   *rateRstDts2L,     /* 15 {I} 2nd index expiry,start,mat. dates */
    double *smile2L,          /* 16 (I) 2nd index smile parameters.       */
    double *sigATM2L,         /* 17 (I) 2nd index ATM vol.                */
    double *vnfmParams2L,     /* 18 (I) 2nd index VNFM parameters.        */
    long   *discDatesL,       /* 19 (I) Discount curve date               */
    double *discRatesL,       /* 20 (I) Discount curve rate               */
    double *corrL,            /* 21 (I) Index correlation data.           */
    long   *optTypeL,         /* 22 (I) Option type.                      */
    long   *payDateL,         /* 23 (I) Option payment date               */
    long   *setlTypeL,        /* 24 (I) Cash or physical settlement.      */
    double *payoffParamsL,    /* 25 (I) Payment parameters                */
    long   *traceL,           /* 26 (I) Trace                             */
    TMatrix2D **xloutputsL    /* 27 (O) Maginal distributions             */
    );


/*-----------------------------------------------------------------------------
 * Q3MQMarginalDistL
 * 
 * Multi-q marginal distributions
 *
 */
int Q3MQMarginalDistL(
    double *expiryL,          /*  1 (I) Exipation in years               */
    long   *rateType1L,       /*  2 (I) Rate or spread 1                 */
    double *fwdRate1L,        /*  3 (I) Fwd 1                            */
    double *smile1L,          /*  4 (I) Smile 1                          */
    double *sigATM1L,         /*  5 (I) % vol or bp vol for fwd 1        */
    long   *rateType2L,       /*  6 (I) Rate or spread 2                 */
    double *fwdRate2L,        /*  7 (I) Fwd 2                            */
    double *smile2L,          /*  8 (I) Smile 2                          */
    double *sigATM2L,         /*  9 (I) % vol or bp vol for fwd 2        */
    double *corrL,            /* 10 (I) Correlation between rate drivers */
    long   *optTypeL,         /* 11 (I) Option type                      */
    double *payoffParamsL,    /* 12 (I) Payoff parameters                */
    long   *traceL,           /* 13 (I) Trace                            */
    TMatrix2D  **xloutputsL   /* 14 (O) Marginal distributions           */
    );

/*-----------------------------------------------------------------------------
 * Q3MQFAMarginalDist
 *
 * Multi-q mariginal distributions in FA.
 *
 */
int Q3MQFAMarginalDist(
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
    double *results          /* 37 (O) Marignal distributions               */ 
    );

/*-----------------------------------------------------------------------------
 * Q3MQMarginalDist
 * 
 * Multi-q bivariate option pricer with simplified interface
 *
 */
int Q3MQMarginalDist(
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
    double *results          /* 14 (O) Marignal distributions              */ 
    );

#endif