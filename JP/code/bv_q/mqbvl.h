/******************************************************************************
 * Module:       
 * Submodule:
 * File:      mqbvl.h   
 * Function:  Addin interface.  
 * Author:    Interest Rates DR
 * Revision:  $Header$
 *****************************************************************************/
#ifndef MQBVL_H
#define MQBVL_H


/*-----------------------------------------------------------------------------
 * Q3MQBivarPricerL
 *
 * Multi-q bivariate option pricer add-in.
 *
 */
int Q3MQBivarPricerL(
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
    double *outputsL          /* 27 (O) Premium                           */
    );

/*-----------------------------------------------------------------------------
 * Q3MQBivarPricerBusL
 *
 * Multi-q bivariate option pricer add-in.
 * Use business day count fraction.
 *
 */
int Q3MQBivarBusPricerL(
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
    char   *holidayFileL,     /* 26 (I) Holiday file name and DCC         */
    long   *traceL,           /* 27 (I) Trace                             */
    double *outputsL          /* 28 (O) Premium                           */
    );

/*-----------------------------------------------------------------------------
 * Q3SmileBivarPricerL
 *
 * Smile bivariate option pricer add-in.
 *
 */
int Q3SmileBivarPricerL(
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
    double *outputsL          /* 27 (O) Premium                           */
    );

/*-----------------------------------------------------------------------------
 * Q3SmileBivarBusPricerL
 *
 * Smile bivariate option pricer add-in.
 * Use business day count fraction.
 *
 */
int Q3SmileBivarBusPricerL(
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
    char   *holidayFileL,     /* 26 (I) Holiday file name and DCC         */
    long   *traceL,           /* 27 (I) Trace                             */
    double *outputsL          /* 28 (O) Premium                           */
    );

/*-----------------------------------------------------------------------------
 * Q3MQBivarPricerSimpleL
 * 
 * Multi-q bivariate option pricer with simplified interface
 *
 */
int Q3MQBivarPricerSimpleL(
    double *expiryL,          /*  1 (I) Expiration in years              */
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
    double *outputsL          /* 14 (O) Pricing results                  */
    );

/*-----------------------------------------------------------------------------
 * Q3SmileBivarPricerSimpleL
 * 
 * Smile bivariate option pricer with simplified interface
 *
 */
int Q3SmileBivarPricerSimpleL(
    double *expiryL,          /*  1 (I) Expiration in years              */
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
    double *outputsL          /* 14 (O) Pricing results                  */
    );

/*-----------------------------------------------------------------------------
 * Q3SmileBivarPricerCorrCalibL
 * 
 * Smile bivariate correlation calibration with simplified interface
 *
 */
int Q3SmileBivarSimpleCorrCalibL(
    double *premiumL,         /*  1 (I) Premium                          */
    double *expiryL,          /*  2 (I) Expiration in years              */
    double *fwdRate1L,        /*  3 (I) Fwd 1                            */
    double *smile1L,          /*  4 (I) Smile 1                          */
    double *sigATM1L,         /*  5 (I) % vol or bp vol for fwd 1        */
    double *fwdRate2L,        /*  6 (I) Fwd 2                            */
    double *smile2L,          /*  7 (I) Smile 2                          */
    double *sigATM2L,         /*  8 (I) % vol or bp vol for fwd 2        */
    long   *optTypeL,         /*  9 (I) Option type CALL=1, PUT = 2      */
    double *strikeL,          /* 10 (I) strike                           */
    double *payoffParamsL,    /* 11 (I) payoff parameters                */
    double *initialGuessL,    /* 12 (I) initial guess of corr            */
    long   *traceL,           /* 13 (I) Trace                            */
    double *outputsL          /* 14 (O) Pricing results                  */
    );

/*-----------------------------------------------------------------------------
 * Q3MQBivarCorrL
 *
 * Multi-q bivariate correlation calculator
 *
 */
int Q3MQBivarCorrL(
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
    long   *payDateL,         /* 23 (I) Option payment date               */
    long   *setlTypeL,        /* 24 (I) Cash or physical settlement.      */
    long   *traceL,           /* 25 (I) Trace                             */
    double *outputsL          /* 26 (O) Correlation                       */
    );

/*-----------------------------------------------------------------------------
 * Q3MQBivarBusCorrL
 *
 * Multi-q bivariate correlation calculator
 * Use business day count fraction.
 *
 */
int Q3MQBivarBusCorrL(
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
    long   *payDateL,         /* 22 (I) Option payment date               */
    long   *setlTypeL,        /* 23 (I) Cash or physical settlement.      */
    char   *holidayFileL,     /* 24 (I) Holiday file name and DCC         */
    long   *traceL,           /* 25 (I) Trace                             */
    double *outputsL          /* 26 (O) Correlation                       */
    );

/*-----------------------------------------------------------------------------
 * Q3SmileBivarCorrL
 *
 * Smile bivariate correlation calculator
 *
 */
int Q3SmileBivarCorrL(
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
    long   *payDateL,         /* 23 (I) Option payment date               */
    long   *setlTypeL,        /* 24 (I) Cash or physical settlement.      */
    long   *traceL,           /* 25 (I) Trace                             */
    double *outputsL          /* 26 (O) Correlation                       */
    );

/*-----------------------------------------------------------------------------
 * Q3SmileBivarBusCorrL
 *
 * Smile bivariate correlation calculator
 * Use business day count fraction.
 */
int Q3SmileBivarBusCorrL(
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
    long   *payDateL,         /* 23 (I) Option payment date               */
    long   *setlTypeL,        /* 24 (I) Cash or physical settlement.      */
    char   *holidayFileL,     /* 25 (I) Holiday file name and DCC         */
    long   *traceL,           /* 26 (I) Trace                             */
    double *outputsL          /* 27 (O) Correlation                       */
    );

/*-----------------------------------------------------------------------------
 * Q3BivariateVersionL
 *
 * Add-in version number
 */

int Q3BivarVersionL(
    char* versionL
    );


/*----------------------------------------------------------------------------
 * Q3QuasiErrBuffInit
 *
 */
void Q3BivarErrBuffInit(
    void
    );


/*----------------------------------------------------------------------------
 * Q3QuasiErrBuffRetrieve  
 *
 */
int Q3BivarErrBuffRetrieve(
    long  outBufLen,   
    char *outBuffer
    );


/*----------------------------------------------------------------------------
 * Q3QuasiErrBuffLength
 *
 */
int Q3BivarErrBuffLength();


/*-----------------------------------------------------------------------------
 * load ALIB holidays
 */

int    Q3BivarHolidayLoadFromDateList   (char    *holName,
                                         long    *datesL,
                                         long    satIsAlwaysHoliday,
                                         long    sunIsAlwaysHoliday);

#endif /* MQBVL_H */
