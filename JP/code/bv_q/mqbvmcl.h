/******************************************************************************
 * Module:      
 * Submodule:
 * File:      mqbvmcl.h 
 * Function:  Addin interface.
 * Author:        Interest Rates DR
 * Revision:  $Header$
 *****************************************************************************/
#ifndef Q3BVMCL_H
#define Q3BVMCL_H


/*-----------------------------------------------------------------------------
 * Q3MQBivarMCPricerL 
 *
 * Multi-q bivariate option pricer add-in.
 *
 */
int Q3MQBivarMCPricerL(
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
    double *payoffParamsL,    /* 25 (I) payment parameters                */
    long   *traceL,           /* 26 (I) trace                             */
    double *outputsL          /* 27 (O) premium                           */
);

/*-----------------------------------------------------------------------------
 * Q3SmileBivarMCPricerL 
 *
 * Smile bivariate option pricer add-in.
 *
 */
int Q3SmileBivarMCPricerL(
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
    double *payoffParamsL,    /* 25 (I) payment parameters                */
    long   *traceL,           /* 26 (I) trace                             */
    double *outputsL          /* 27 (O) premium                           */
);

/*-----------------------------------------------------------------------------
 * Q3MQBivarMCPricerSimple
 * 
 * Multi-q bivariate option pricer with simplified interface
 *
 */
int Q3MQBivarMCPricerSimpleL(
    double *expiryL,          /*  1 (I) exipation in years               */
    long   *rateType1L,       /*  2 (I) rate or spread 1                 */
    double *fwdRate1L,        /*  3 (I) fwd 1                            */
    double *smile1L,          /*  4 (I) smile 1                          */
    double *sigATM1L,         /*  5 (I) % vol or bp vol for fwd 1        */
    long   *rateType2L,       /*  6 (I) rate or spread 2                 */
    double *fwdRate2L,        /*  7 (I) fwd 2                            */
    double *smile2L,          /*  8 (I) smile 2                          */
    double *sigATM2L,         /*  9 (I) % vol or bp vol for fwd 2        */
    double *corrL,            /* 10 (I) correlation between rate drivers */
    long   *optTypeL,         /* 11 (I) option type                      */
    double *payoffParamsL,    /* 12 (I) payoff parameters                */
    long   *traceL,           /* 13 (I) trace                            */
    double *outputsL          /* 14 (O) pricing results                  */
    );

/*-----------------------------------------------------------------------------
 * Q3MQBivarMCMidCurvePricerL
 *
 * Multi-q bivariate midcurve swaption pricer add-in.
 *
 */
int Q3MQBivarMCMidCurvePricerL( 
    long   *todayL,           /*  1 (I) Vol base date                     */
    long   *valueDateL,       /*  2 (I) Value date                        */
    long   *indxDatesL,       /*  3 (I) index curve date                  */
    double *indxRatesL,       /*  4 (I) index curve rate                  */
    long   *rateFreq0L,       /*  5 (I) midcurve freq                     */
    char   *rateDCC0L,        /*  6 (I) midcurve DCC                      */
    long   *rateFreq1L,       /*  7 (I) 1st index freq                    */
    char   *rateDCC1L,        /*  8 (I) 1st index DCC                     */
    long   *rateRstDts1L,     /*  9 {I} 1st index expiry,start,mat. dates */
    double *smile1L,          /* 10 (I) 1st index smile parameters.       */
    double *sigATM1L,         /* 11 (I) 1st index ATM vol.                */
    long   *rateFreq2L,       /* 12 (I) 2nd index freq                    */
    char   *rateDCC2L,        /* 13 (I) 2nd index DCC                     */
    long   *rateRstDts2L,     /* 14 (I) 2nd index expiry,start,mat. dates */
    double *smile2L,          /* 15 (I) 2nd index smile parameters.       */
    double *sigATM2L,         /* 16 (I) 2nd index ATM vol.                */
    char   *fixedTypeL,       /* 17 (I) Type of payment(N,E,I,J) */
    char   *stubTypeL,        /* 18 (I) Stub type (S,B,D,N)               */
    char   *stubAtEndL,       /* 19 (I) (F)ront or (B)ack                 */
    long   *discDatesL,       /* 20 (I) Discount curve date               */
    double *discRatesL,       /* 21 (I) Discount curve rate               */
    double *vnfmParamsL,      /* 22 (I) index VNFM parameters             */
    double *corrL,            /* 23 (I) Index correlation data.           */
    long   *optTypeL,         /* 24 (I) Option type.                      */
    long   *payDateL,         /* 25 (I) payment date                      */
    long   *setlTypeL,        /* 26 (I) Cash or physical settlement.      */
    double *strikeL,          /* 27 (I) payment parameters                */
    long   *traceL,           /* 28 (I) trace                             */
    double *outputsL          /* 29 (O) premium, etc                      */
    );

/*-----------------------------------------------------------------------------
 * Q3SmileBivarMCMidCurvePricerL
 *
 * Multi-q bivariate midcurve swaption pricer add-in.
 *
 */
int Q3SmileBivarMCMidCurvePricerL(
    long   *todayL,           /*  1 (I) Vol base date                     */
    long   *valueDateL,       /*  2 (I) Value date                        */
    long   *indxDatesL,       /*  3 (I) index curve date                  */
    double *indxRatesL,       /*  4 (I) index curve rate                  */
    long   *rateFreq0L,       /*  5 (I) midcurve freq                     */
    char   *rateDCC0L,        /*  6 (I) midcurve DCC                      */
    long   *rateFreq1L,       /*  7 (I) 1st index freq                    */
    char   *rateDCC1L,        /*  8 (I) 1st index DCC                     */
    long   *rateRstDts1L,     /*  9 {I} 1st index expiry,start,mat. dates */
    double *smile1L,          /* 10 (I) 1st index smile parameters.       */
    double *sigATM1L,         /* 11 (I) 1st index ATM vol.                */
    long   *rateFreq2L,       /* 12 (I) 2nd index freq                    */
    char   *rateDCC2L,        /* 13 (I) 2nd index DCC                     */
    long   *rateRstDts2L,     /* 14 (I) 2nd index expiry,start,mat. dates */
    double *smile2L,          /* 15 (I) 2nd index smile parameters.       */
    double *sigATM2L,         /* 16 (I) 2nd index ATM vol.                */
    char   *fixedTypeL,       /* 17 (I) Type of payment(N,E,I,J)          */
    char   *stubTypeL,        /* 18 (I) Stub type (S,B,D,N)               */
    char   *stubAtEndL,       /* 19 (I) (F)ront or (B)ack                 */
    long   *discDatesL,       /* 20 (I) Discount curve date               */
    double *discRatesL,       /* 21 (I) Discount curve rate               */
    double *vnfmParamsL,      /* 22 (I) index VNFM parameters             */
    double *corrL,            /* 23 (I) Index correlation data.           */
    long   *optTypeL,         /* 24 (I) Option type.                      */
    long   *payDateL,         /* 25 (I) payment date                      */
    long   *setlTypeL,        /* 26 (I) Cash or physical settlement.      */
    double *strikeL,          /* 27 (I) payment parameters                */
    long   *traceL,           /* 28 (I) trace                             */
    double *outputsL          /* 29 (O) premium, etc                      */
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


#endif /* Q3BVMCL_H */
