/******************************************************************************
 * Module:
 * Submodule:
 * File:        mqbvmc.c
 * Function:
 * Author:      Interest Rates DR
 * Revision:    $Header$
 *****************************************************************************/
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include "aliberr.h"
#include "aliblog.h"

#include "mqbvmcl.h"
#include "mqbvmc.h"
#include "mqbv.h"

#define Q3_BV_NCK  2000000769 /* Calibration control data */
#define Q3_BV_RATE   0
#define Q3_BV_SPREAD 1

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
    )
{
    static char routine[] = "Q3MQBivarMCPricerL";
    int status = FAILURE;
    double  smileSmile1L[20], smileSmile2L[20];

    memcpy(smileSmile1L+2, smile1L+1, ((long) smile1L[0])*sizeof(double));
    smileSmile1L[0] = smile1L[0] + 1.;
    smileSmile1L[1] = 10.;

    memcpy(smileSmile2L+2, smile2L+1, ((long) smile2L[0])*sizeof(double));
    smileSmile2L[0] = smile2L[0] + 1.;
    smileSmile2L[1] = 10.;

    if(Q3SmileBivarMCPricerL(
                             todayL,
                             valueDateL,
                             indxDates1L,
                             indxRates1L,
                             rateFreq1L,
                             rateDCC1L,
                             rateRstDts1L,
                             smileSmile1L,
                             sigATM1L,
                             vnfmParams1L,
                             indxDates2L,
                             indxRates2L,
                             rateFreq2L,
                             rateDCC2L,
                             rateRstDts2L,
                             smileSmile2L,
                             sigATM2L,
                             vnfmParams2L,
                             discDatesL,
                             discRatesL,
                             corrL,
                             optTypeL,
                             payDateL,
                             setlTypeL,
                             payoffParamsL,
                             traceL,
                             outputsL) == FAILURE)
         goto RETURN;

    status = SUCCESS;
    
RETURN:

    if (status == FAILURE) {
        Q3ErrMsg("%s: Failed\n", routine);
    }

    return status;

} /* Q3MQBivarMCPricerL */

/*-----------------------------------------------------------------------------
 * Q3SmileBivarMCPricerL
 * 
 * Multi-q bivariate option pricer add-in.
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
    )
{
    static char routine[] = "Q3SmileBivarMCPricerL";

    int status = FAILURE;

    long expiryDate1L[2] = {1, 0};
    long expiryDate2L[2] = {1, 0};
    long startDate1L[2]  = {1, 0};
    long startDate2L[2]  = {1, 0};
    long matDate1L[2]    = {1, 0};
    long matDate2L[2]    = {1, 0};

    double smile1Input[13], smile2Input[13];
    int    i;
    
    /* initialize alib callbacks */
    Q3ToAlibErrInit();
    
    /* initialize global error buffer */
    Q3BivarErrBuffInit();

    /* Test input */
    GTO_WRAP_CHECK_SCALAR(todayL,               "today");
    GTO_WRAP_CHECK_SCALAR(valueDateL,           "value date");
    GTO_WRAP_CHECK_SCALAR(rateFreq1L,           "rate frequency");
    GTO_WRAP_CHECK_SCALAR(rateDCC1L,            "rate day count convention");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(rateRstDts1L,"rate reset dates", 3);
    GTO_WRAP_CHECK_SCALAR(sigATM1L,             "ATM volatility");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile1L,     "smile", 4);
    GTO_WRAP_CHECK_SCALAR(rateFreq2L,           "rate frequency");
    GTO_WRAP_CHECK_SCALAR(rateDCC2L,            "rate day count convention");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(rateRstDts2L,"rate reset dates", 3);
    GTO_WRAP_CHECK_SCALAR(sigATM2L,             "ATM volatility");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile2L,     "smile", 4);
    GTO_WRAP_CHECK_SCALAR(optTypeL,             "option type");
    GTO_WRAP_CHECK_SCALAR(payDateL,             "payment date");
    GTO_WRAP_CHECK_SCALAR(setlTypeL,            "settlement type");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(corrL,       "corr", 1);
    GTO_WRAP_CHECK_VECTOR_LEN(outputsL,         "outputs", 3);
    
    /* unpack rate reset dates */
    expiryDate1L[1] = rateRstDts1L[1];
    startDate1L[1]  = rateRstDts1L[2];
    matDate1L[1]    = rateRstDts1L[3];
    expiryDate2L[1] = rateRstDts2L[1];
    startDate2L[1]  = rateRstDts2L[2];
    matDate2L[1]    = rateRstDts2L[3];
    
    /*  Alib RegTest logging if tracing is on */
    if (traceL[1] > 0) 
    {
        Q3LilVectLoggingFile(
            "q3.log",
            (traceL[1] == 1 ? "w" : "a"),
            routine,
            Q3_TDATE_L,      todayL,        "today",
            Q3_TDATE_L,      valueDateL,    "valueDate",
            Q3_TDATE_L,      indxDates1L,   "indxDates1",
            Q3_DOUBLE_L,     indxRates1L,   "indxRates1",
            Q3_LONG_L,       rateFreq1L,    "rateFreq1",
            Q3_CHAR_BLOCK_L, rateDCC1L,     "rateDCC1",     
            Q3_TDATE_L,      expiryDate1L,  "expiryDate1",
            Q3_TDATE_L,      startDate1L,   "startDate1",
            Q3_TDATE_L,      matDate1L,     "matDate1",
            Q3_DOUBLE_L,     smile1L,       "smile1",
            Q3_DOUBLE_L,     sigATM1L,      "sigATM1",
            Q3_DOUBLE_L,     vnfmParams1L,  "vnfmParams1",
            Q3_TDATE_L,      indxDates2L,   "indxDates2",
            Q3_DOUBLE_L,     indxRates2L,   "indxRates2",
            Q3_LONG_L,       rateFreq2L,    "rateFreq2",
            Q3_CHAR_BLOCK_L, rateDCC2L,     "rateDCC2",     
            Q3_TDATE_L,      expiryDate2L,  "expiryDate2",
            Q3_TDATE_L,      startDate2L,   "startDate2",
            Q3_TDATE_L,      matDate2L,     "matDate2",
            Q3_DOUBLE_L,     smile2L,       "smile2",
            Q3_DOUBLE_L,     sigATM2L,      "sigATM2",
            Q3_DOUBLE_L,     vnfmParams2L,  "vnfmParams2",
            Q3_TDATE_L,      discDatesL,    "discDates",
            Q3_DOUBLE_L,     discRatesL,    "discRates",
            Q3_DOUBLE_L,     corrL,         "corr",
            Q3_LONG_L,       optTypeL,      "optType",
            Q3_TDATE_L,      payDateL,      "payDate",
            Q3_LONG_L,       setlTypeL,     "setlType",
            Q3_DOUBLE_L,     payoffParamsL, "payoffParams",        
            Q3_LONG_L,       traceL,        "trace",
            0);
    }  

    memcpy(smile1Input, smile1L+1, ((int)smile1L[0])*sizeof(double));

    for (i=(int)smile1L[0] ; i<13 ; i++)
        smile1Input[i] = 0.;

    memcpy(smile2Input, smile2L+1, ((int)smile2L[0])*sizeof(double));

    for (i=(int)smile2L[0] ; i<13 ; i++)
        smile2Input[i] = 0.;

    if (Q3MQBivarMCPricer(
        todayL[1],         
        valueDateL[1],     
        indxDates1L[0],   
        indxDates1L+1,    
        indxRates1L+1,     
        rateFreq1L[1],      
        rateDCC1L+1,       
        expiryDate1L[1], 
        startDate1L[1], 
        matDate1L[1],   
        smile1Input,  
        sigATM1L[1],    
        (long) vnfmParams1L[0],
        vnfmParams1L+1,   
        indxDates2L[0],
        indxDates2L+1,     
        indxRates2L+1,     
        rateFreq2L[1],      
        rateDCC2L+1,       
        expiryDate2L[1], 
        startDate2L[1], 
        matDate2L[1],   
        smile2Input,  
        sigATM2L[1],   
        (long) vnfmParams2L[0],
        vnfmParams2L+1,   
        discDatesL[0],    
        discDatesL+1,     
        discRatesL+1,     
        corrL[1],          
        (corrL[0]>1 ? corrL+2 : NULL),
        optTypeL[1],       
        payDateL[1],       
        setlTypeL[1],     
        (long) payoffParamsL[0],
        payoffParamsL+1,
        outputsL+1) == FAILURE) goto RETURN;  
  
    status = SUCCESS;

 done:   /* for alib */
 RETURN:

    if (status == FAILURE){
        Q3ErrMsg("%s: failed.\n", routine);
    }
  
    return status;

} /* Q3SmileBivarMCPricerL */

/*-----------------------------------------------------------------------------
 * Q3MQBivarMCPricerSimpleL
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
    )
{
    static char routine[] = "Q3MQBivarMCPricerSimpleL";

    int status = FAILURE;
    int i;

    double smile1Input[13], smile2Input[13];
    
    double  smileSmile1L[20], smileSmile2L[20];

    memcpy(smileSmile1L+2, smile1L+1, ((long) smile1L[0])*sizeof(double));
    smileSmile1L[0] = smile1L[0] + 1.;
    smileSmile1L[1] = 10.;

    memcpy(smileSmile2L+2, smile2L+1, ((long) smile2L[0])*sizeof(double));
    smileSmile2L[0] = smile2L[0] + 1.;
    smileSmile2L[1] = 10.;

    /* initialize alib callbacks */
    Q3ToAlibErrInit();
    
    /* initialize global error buffer */
    Q3BivarErrBuffInit();

    /* Test input */
    GTO_WRAP_CHECK_SCALAR(expiryL,              "expiry");
    GTO_WRAP_CHECK_SCALAR(rateType1L,           "rateType1");
    GTO_WRAP_CHECK_SCALAR(fwdRate1L,            "fwdRate1");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile1L,     "smile1", 10);
    GTO_WRAP_CHECK_SCALAR(sigATM1L,             "sigATM1");
    GTO_WRAP_CHECK_SCALAR(rateType2L,           "rateType2");
    GTO_WRAP_CHECK_SCALAR(fwdRate2L,            "fwdRate2");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile2L,     "smile2", 10);
    GTO_WRAP_CHECK_SCALAR(sigATM2L,             "sigATM2");
    GTO_WRAP_CHECK_SCALAR(corrL,                "corr");
    GTO_WRAP_CHECK_SCALAR(optTypeL,             "optType");
    GTO_WRAP_CHECK_VECTOR_LEN(outputsL,         "outputs", 3);
    
    /*  Alib RegTest logging if tracing is on */
    if (traceL[1] > 0) 
    {
        Q3LilVectLoggingFile(
            "q3.log",
            (traceL[1] == 1 ? "w" : "a"),
            routine,
            Q3_DOUBLE_L,     expiryL,       "expiry",
            Q3_LONG_L,       rateType1L,    "rateType1",
            Q3_DOUBLE_L,     fwdRate1L,     "fwdRate1",
            Q3_DOUBLE_L,     smile1L,       "smile1",
            Q3_DOUBLE_L,     sigATM1L,      "sigATM1",
            Q3_LONG_L,       rateType2L,    "rateType2",
            Q3_DOUBLE_L,     fwdRate2L,     "fwdRate2",
            Q3_DOUBLE_L,     smile2L,       "smile2",
            Q3_DOUBLE_L,     sigATM2L,      "sigATM2",
            Q3_DOUBLE_L,     corrL,         "corr",
            Q3_LONG_L,       optTypeL,      "optType",
            Q3_DOUBLE_L,     payoffParamsL, "payoffParams",        
            Q3_LONG_L,       traceL,        "trace",
            0);
    }  

    memcpy(smile1Input, smileSmile1L+1, ((int)smileSmile1L[0])*sizeof(double));

    for (i=(int)smileSmile1L[0] ; i<13 ; i++)
        smile1Input[i] = 0.;

    memcpy(smile2Input, smileSmile2L+1, ((int)smileSmile2L[0])*sizeof(double));

    for (i=(int)smileSmile2L[0] ; i<13 ; i++)
        smile2Input[i] = 0.;

    
    if (Q3MQBivarMCPricerSimple(
        expiryL[1], 
        rateType1L[1],
        fwdRate1L[1],
        smile1Input,  
        sigATM1L[1],    
        rateType2L[1],
        fwdRate2L[1],
        smile2Input,  
        sigATM2L[1],   
        corrL[1],          
        optTypeL[1],       
        (long) payoffParamsL[0],
        payoffParamsL+1,
        outputsL+1) == FAILURE) goto RETURN;  

    status = SUCCESS;

 done: /* for alib */
 RETURN:

    if (status == FAILURE)
    {
        Q3ErrMsg("%s: failed.\n", routine);
    }
  
    return status;

} /* Q3MQBivarPricerSimpleL */

/*-----------------------------------------------------------------------------
 * Q3MQBivarMCMidCurvePricerL
 *
 * Multi-q bivariate midcurve swaption pricer add-in.
 *
 */
int Q3MQBivarMCMidCurvePricerL(
    long   *todayL,           /*  1 (I) Vol base date                     */
    long   *valueDateL,       /*  2 (I) Value date                        */
    long   *indxDatesL,       /*  3 (I) Index curve date                  */
    double *indxRatesL,       /*  4 (I) Index curve rate                  */
    long   *rateFreq0L,       /*  5 (I) Midcurve freq                     */
    char   *rateDCC0L,        /*  6 (I) Midcurve DCC                      */
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
    double *vnfmParamsL,      /* 22 (I) Index VNFM parameters             */
    double *corrL,            /* 23 (I) Index correlation                 */
    long   *optTypeL,         /* 24 (I) Option type.                      */
    long   *payDateL,         /* 25 (I) Option payment date               */
    long   *setlTypeL,        /* 26 (I) Cash or physical settlement.      */
    double *strikeL,          /* 27 (I) Strike                            */
    long   *traceL,           /* 28 (I) Trace                             */
    double *outputsL          /* 29 (O) Option premium, etc               */
    )
{
    static char routine[] = "Q3MQBivarMCMidCurvePricerL";
    int status = FAILURE;
    double smileSmile1L[20], smileSmile2L[20];

    memcpy(smileSmile1L+2, smile1L+1, ((long) smile1L[0])*sizeof(double));
    smileSmile1L[0] = smile1L[0] + 1.;
    smileSmile1L[1] = 10.;

    memcpy(smileSmile2L+2, smile2L+1, ((long) smile2L[0])*sizeof(double));
    smileSmile2L[0] = smile2L[0] + 1.;
    smileSmile2L[1] = 10.;

    if (Q3SmileBivarMCMidCurvePricerL(
        todayL,
        valueDateL,
        indxDatesL,
        indxRatesL,
        rateFreq0L,
        rateDCC0L,
        rateFreq1L,
        rateDCC1L,
        rateRstDts1L,
        smileSmile1L,
        sigATM1L,
        rateFreq2L,
        rateDCC2L,
        rateRstDts2L,
        smileSmile2L,
        sigATM2L,
        fixedTypeL,
        stubTypeL,
        stubAtEndL,
        discDatesL,
        discRatesL,
        vnfmParamsL,
        corrL,
        optTypeL,
        payDateL,
        setlTypeL,
        strikeL,
        traceL,
        outputsL) == FAILURE) goto RETURN;

    status = SUCCESS;

RETURN:

    if (status == FAILURE) {
        Q3ErrMsg("%s: Failed\n", routine);
    }

    return status;
} /* Q3MQBivarMCMidCurvePricerL */

/*-----------------------------------------------------------------------------
 * Q3SmileBivarMCMidCurvePricerL
 *
 * Multi-q bivariate midcurve MC swaption pricer add-in.
 *
 */ 
int Q3SmileBivarMCMidCurvePricerL(
    long   *todayL,           /*  1 (I) Vol base date                     */
    long   *valueDateL,       /*  2 (I) Value date                        */
    long   *indxDatesL,       /*  3 (I) Index curve date                  */
    double *indxRatesL,       /*  4 (I) Index curve rate                  */
    long   *rateFreq0L,       /*  5 (I) Midcurve freq                     */
    char   *rateDCC0L,        /*  6 (I) Midcurve DCC                      */
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
    double *vnfmParamsL,      /* 22 (I) Index VNFM parameters             */
    double *corrL,            /* 23 (I) Index correlation                 */
    long   *optTypeL,         /* 24 (I) Option type.                      */
    long   *payDateL,         /* 25 (I) Option payment date               */
    long   *setlTypeL,        /* 26 (I) Cash or physical settlement.      */
    double *strikeL,          /* 27 (I) Strike                            */
    long   *traceL,           /* 28 (I) Trace                             */
    double *outputsL          /* 29 (O) Option premium, etc               */
    )
{
    static char routine[] = "Q3SmileBivarMCMidCurvePricerL";
    int status = FAILURE;

    long expiryDate1L[2] = {1, 0};
    long expiryDate2L[2] = {1, 0};
    long startDate1L[2]  = {1, 0};
    long startDate2L[2]  = {1, 0};
    long matDate1L[2]    = {1, 0};
    long matDate2L[2]    = {1, 0};

    double smile1Input[13], smile2Input[13];
    int    i;
 
    /* initialize alib callbacks */
    Q3ToAlibErrInit();

    /* initialize global error buffer */
    Q3BivarErrBuffInit();

    /* Test input */
    GTO_WRAP_CHECK_SCALAR(todayL,               "today");
    GTO_WRAP_CHECK_SCALAR(valueDateL,           "value date");
    GTO_WRAP_CHECK_SCALAR(rateFreq0L,           "rate frequency");
    GTO_WRAP_CHECK_SCALAR(rateDCC0L,            "rate day count convention");
    GTO_WRAP_CHECK_SCALAR(rateFreq1L,           "rate frequency");
    GTO_WRAP_CHECK_SCALAR(rateDCC1L,            "rate day count convention");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(rateRstDts1L,"rate reset dates", 3);
    GTO_WRAP_CHECK_SCALAR(sigATM1L,             "ATM volatility");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile1L,     "smile", 4);
    GTO_WRAP_CHECK_SCALAR(rateFreq2L,           "rate frequency");
    GTO_WRAP_CHECK_SCALAR(rateDCC2L,            "rate day count convention");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(rateRstDts2L,"rate reset dates", 3);
    GTO_WRAP_CHECK_SCALAR(sigATM2L,             "ATM volatility");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile2L,     "smile", 4);
    GTO_WRAP_CHECK_SCALAR(fixedTypeL,           "fixed period type");
    GTO_WRAP_CHECK_SCALAR(stubTypeL,            "stub type");
    GTO_WRAP_CHECK_SCALAR(stubAtEndL,           "stub at end flag");
    GTO_WRAP_CHECK_SCALAR(optTypeL,             "option type");
    GTO_WRAP_CHECK_SCALAR(payDateL,             "payment date");
    GTO_WRAP_CHECK_SCALAR(setlTypeL,            "settlement type");
    GTO_WRAP_CHECK_SCALAR(corrL,                "corr");
    GTO_WRAP_CHECK_VECTOR_LEN(outputsL,         "outputs", 3);

    /* unpack rate reset dates */
    expiryDate1L[1] = rateRstDts1L[1];
    startDate1L[1]  = rateRstDts1L[2];
    matDate1L[1]    = rateRstDts1L[3];
    expiryDate2L[1] = rateRstDts2L[1];
    startDate2L[1]  = rateRstDts2L[2];
    matDate2L[1]    = rateRstDts2L[3];

    /*  Alib RegTest logging if tracing is on */
    if (traceL[1] > 0)
    {
        Q3LilVectLoggingFile(
            "q3.log",
            (traceL[1] == 1 ? "w" : "a"),
            routine,
            Q3_TDATE_L,      todayL,        "today",
            Q3_TDATE_L,      valueDateL,    "valueDate",
            Q3_TDATE_L,      indxDatesL,    "indxDates",
            Q3_DOUBLE_L,     indxRatesL,    "indxRates",
            Q3_LONG_L,       rateFreq0L,    "rateFreq0",
            Q3_CHAR_BLOCK_L, rateDCC0L,     "rateDCC0",
            Q3_LONG_L,       rateFreq1L,    "rateFreq1",
            Q3_CHAR_BLOCK_L, rateDCC1L,     "rateDCC1",   
            Q3_TDATE_L,      expiryDate1L,  "expiryDate1",
            Q3_TDATE_L,      startDate1L,   "startDate1",
            Q3_TDATE_L,      matDate1L,     "matDate1",
            Q3_DOUBLE_L,     smile1L,       "smile1",
            Q3_DOUBLE_L,     sigATM1L,      "sigATM1",
            Q3_LONG_L,       rateFreq2L,    "rateFreq2",
            Q3_CHAR_BLOCK_L, rateDCC2L,     "rateDCC2",   
            Q3_TDATE_L,      expiryDate2L,  "expiryDate2",
            Q3_TDATE_L,      startDate2L,   "startDate2",
            Q3_TDATE_L,      matDate2L,     "matDate2",
            Q3_DOUBLE_L,     smile2L,       "smile2",
            Q3_DOUBLE_L,     sigATM2L,      "sigATM2",
            Q3_CHAR_BLOCK_L, fixedTypeL,    "fixedType",
            Q3_CHAR_BLOCK_L, stubTypeL,     "stubType",
            Q3_CHAR_BLOCK_L, stubAtEndL,    "stubAtEnd",
            Q3_TDATE_L,      discDatesL,    "discDates",
            Q3_DOUBLE_L,     discRatesL,    "discRates",
            Q3_DOUBLE_L,     vnfmParamsL,   "vnfmParams",
            Q3_DOUBLE_L,     corrL,         "corr",
            Q3_LONG_L,       optTypeL,      "optType",
            Q3_TDATE_L,      payDateL,      "payDate",
            Q3_LONG_L,       setlTypeL,     "setlType",
            Q3_DOUBLE_L,     strikeL,       "strike",
            Q3_LONG_L,       traceL,        "trace",
            0);
    }

    memcpy(smile1Input, smile1L+1, ((int)smile1L[0])*sizeof(double));

    for (i=(int)smile1L[0] ; i<13 ; i++)
        smile1Input[i] = 0.;

    memcpy(smile2Input, smile2L+1, ((int)smile2L[0])*sizeof(double));

    for (i=(int)smile2L[0] ; i<13 ; i++)
        smile2Input[i] = 0.;

    if (Q3MQBivarMCMidCurvePricer(
        todayL[1],
        valueDateL[1],
        indxDatesL[0],
        indxDatesL+1,
        indxRatesL+1,
        rateFreq0L[1],
        rateDCC0L+1,
        rateFreq1L[1],
        rateDCC1L+1,
        expiryDate1L[1],
        startDate1L[1],
        matDate1L[1],
        smile1Input,
        sigATM1L[1],
        rateFreq2L[1],
        rateDCC2L+1,
        expiryDate2L[1],
        startDate2L[1],
        matDate2L[1],
        smile2Input,
        sigATM2L[1],
        fixedTypeL    + WRAP_STR_IDX(1),
        stubTypeL     + WRAP_STR_IDX(1),
        stubAtEndL    + WRAP_STR_IDX(1),
        discDatesL[0],
        discDatesL+1,
        discRatesL+1,
        (long) vnfmParamsL[0],
        vnfmParamsL+1,
        corrL[1],
        optTypeL[1],
        payDateL[1],
        setlTypeL[1],
        strikeL[1],
        outputsL+1) == FAILURE) goto RETURN;

    status = SUCCESS;

 done:   /* for alib */
 RETURN:

    if (status == FAILURE){
        Q3ErrMsg("%s: failed.\n", routine);
    }

    return status;

} /* Q3SmileBivarMCMidCurvePricerL */

/*-----------------------------------------------------------------------------
 * Q3MQBivarMCPricer 
 *
 * Multi-q bivariate option pricer.
 *
 */
int Q3MQBivarMCPricer(
    long    today,           /* (I) Vol base date                        */
    long    valueDate,       /* (I) Value date                           */
    long    numIndxPts1,     /* (I) 1st index curve num points           */
    long   *indxDates1,      /* (I) 1st index curve dates                */
    double *indxRates1,      /* (I) 1st index curve rates                */
    long    rateFreq1,       /* (I) 1st index curve frequency            */
    char   *rateDCC1,        /* (I) 1st index curve day count convention */
    long    expiryDate1,     /* (I) 1st index reset date                 */
    long    startDate1,      /* (I) 1st index start date                 */
    long    matDate1,        /* (I) 1st index rate end date              */
    double *smile1,          /* (I) 1st index smile parameters           */
    double  sigATM1,         /* (I) 1st index ATM vol                    */
    long    numVNFMParams1,  /* (I) 1st index num VNFM parameters        */
    double *vnfmParams1,     /* (I) 1st index VNFM parameters            */
    long    numIndxPts2,     /* (I) 2nd index curve num points           */
    long   *indxDates2,      /* (I) 2nd index curve dates                */
    double *indxRates2,      /* (I) 2nd index curve rates                */
    long    rateFreq2,       /* (I) 2nd index curve frequency            */
    char   *rateDCC2,        /* (I) 2nd index curve day count convention */
    long    expiryDate2,     /* (I) 2nd index reset date                 */
    long    startDate2,      /* (I) 2nd index start date                 */
    long    matDate2,        /* (I) 2nd index rate end date              */
    double *smile2,          /* (I) 2nd index smile parameters           */
    double  sigATM2,         /* (I) 2nd index ATM vol                    */
    long    numVNFMParams2,  /* (I) 2nd index num VNFM parameters        */
    double *vnfmParams2,     /* (I) 2nd index VNFM parameters            */
    long    numDiscPts,      /* (I) Discount curve num points            */
    long   *discDates,       /* (I) Discount curve dates                 */
    double *discRates,       /* (I) Discount curve rates                 */
    double  corr,            /* (I) Index correlation data               */
    double *fwdVolInput,     /* (I) Forward vol betw rate start dates    */
    long    optType,         /* (I) Option type                          */
    long    payDate,         /* (I) Option payment date                  */
    long    setlType,        /* (I) Cash or physical settlement          */
    long    numPayoffParams, /* (I) Num payment parameters               */
    double *payoffParams,    /* (I) Payment parameters                   */
    double *results          /* (O) Pricing results.                     */ 
    )   
{
    static char routine[] = "Q3MQBivarMCPricer";
    int         status    = FAILURE;

    union{
        struct{
            double fwdRate;
            double fwdRate30360;
            double fwdAnn30360;
            double zeroRateSwap;
            double zeroRatePay;
            double expiry;
            double expiryVolvol;
            double start;
            double swapMat;
            double payDelay;
            double freq;
        } s;
        double a[11];
    } rate[2];

    FPAYOFF  *payFunc       = NULL;
    double   *smile         = NULL;
    double   *indxRates     = NULL;
    double   *vnfmParams    = NULL;
    long     *indxDates     = NULL;
    long      numVNFMParams = 0;
    char     *rateDCC       = NULL;

    PAYOFF    pf;
    FADATA    fa[2];
    MQDATA    mq[2];
    MQDATA    pa[2];
    MQDATA   *ppa[2];
    double    sigATM;
    double    price;
    long      numIndxPts;
    long      rateFreq;
    long      startDate, matDate, expiryDate;
    int       i;

    /* calibrate PA measure for each index */
    for (i = 0; i < 2; i++)
    {
        
        /* select index */
        if (i == 0)
        {
            sigATM          = sigATM1;
            smile           = smile1;
            numIndxPts      = numIndxPts1;
            indxDates       = indxDates1;
            indxRates       = indxRates1;
            rateFreq        = rateFreq1;
            rateDCC         = rateDCC1;
            expiryDate      = expiryDate1;
            startDate       = startDate1;
            matDate         = matDate1;
            numVNFMParams   = numVNFMParams1;
            vnfmParams      = vnfmParams1;
        }
        else
        {
            sigATM          = sigATM2;
            smile           = smile2;
            numIndxPts      = numIndxPts2;
            indxDates       = indxDates2;
            indxRates       = indxRates2;
            rateFreq        = rateFreq2;
            rateDCC         = rateDCC2;
            expiryDate      = expiryDate2;
            startDate       = startDate2;
            matDate         = matDate2;
            numVNFMParams   = numVNFMParams2;
            vnfmParams      = vnfmParams2;
        }

        /* Fwd dates and rates for index */
        if (Q3SwapRateCalc2(
            today,
            today,
            valueDate,
            numIndxPts,
            indxDates,
            indxRates,
            numDiscPts,
            discDates,
            discRates,
            rateFreq,
            rateDCC,
            expiryDate,
            startDate,
            matDate,
            payDate,
            NULL,
            NULL,
            (rate[i]).a) == FAILURE) goto RETURN;

        if (sigATM * sqrt((rate[i].s.expiry)) < Q3_MIN_VOL_CALIB) corr = 0.;

        /* smile: initialize */
        if (Q3SmileInit(
                        (rate[i]).s.fwdRate,
                        sigATM,
                        0.,
                        (rate[i]).s.expiry,
                        (rate[i]).s.expiry,
                        1.,
                        smile,
                        &(mq[i])) == FAILURE) 
        {
            goto RETURN;
        }

        /* Calibrate MQ distribution to SV option prices. */
        if (Q3MQCalib(&(mq[i])) == FAILURE) goto RETURN;

        /* Initialize FA parameters. */
        (fa[i]).mq = &(mq[i]);
        if (Q3FASmileInit(
            (rate[i]).s.expiry,
            sigATM,
            (rate[i]).s.start,
            (long) (rate[i]).s.freq,
            (rate[i]).s.swapMat,
            (rate[i]).s.fwdRate30360,
            (rate[i]).s.fwdAnn30360,
            (rate[i]).s.zeroRateSwap,
            (rate[i]).s.payDelay,
            (rate[i]).s.zeroRatePay,
            numVNFMParams,
            vnfmParams,
            setlType,
            setlType,
            &(fa[i])) == FAILURE) goto RETURN;
     
        /* initialize PA MQ */
        if (Q3MQCopySmileFromMQ(
            &(mq[i]),
            &(pa[i])) == FAILURE) goto RETURN;

        /* set NCK parameters */
        if (Q3DecodeNCK(
            Q3_BV_NCK, 
            &(pa[i])) == FAILURE) goto RETURN;

        /* compute calibration targets */
        if (Q3MQTargetFA(
            &(fa[i]),
            &(pa[i])) == FAILURE) goto RETURN;

        /* bootstrap PA multi-q measure */
        if (Q3MQBootstrapQ(&(pa[i])) == FAILURE) goto RETURN;

        /* calibrate PA */
        if (Q3MQCalib(&(pa[i])) == FAILURE) goto RETURN;

    } /* i */

    /* adjust correlation for rate start date difference */
    if (startDate1 != startDate2)
    {
        double vol, fwdVol, fwdTenor, swapVol, swapFwdVol;
        double zeroVol, zeroSwapCorr, adj;
        int    iLate;

        fwdTenor = (rate[1]).s.start - (rate[0]).s.start;
        if (fwdTenor > 0) 
        {
            iLate = 1;
            vol   = sigATM2;
        }
        else
        {
            iLate    = 0;
            vol      = sigATM1;
            fwdTenor = fabs(fwdTenor);
        }

        if (fwdVolInput == NULL)
        {
            /* estimate forward vol */
                
            if (vol < TINY) goto RETURN;

            if (Q3VNFMZero2Swap(
                rate[iLate].s.expiry,
                rate[iLate].s.expiry,
                rate[iLate].s.start,
                (long) rate[iLate].s.freq,
                rate[iLate].s.swapMat,
                rate[iLate].s.fwdRate30360,
                rate[iLate].s.fwdAnn30360,
                rate[iLate].s.swapMat,
                rate[iLate].s.zeroRateSwap,
                numVNFMParams,
                vnfmParams,
                &swapVol,
                &zeroVol,
                &zeroSwapCorr) == FAILURE) goto RETURN;       

            if (Q3VNFMZero2Swap(
                rate[iLate].s.expiry,
                rate[iLate].s.expiry - fwdTenor,
                rate[iLate].s.start,
                (long) rate[iLate].s.freq,
                rate[iLate].s.swapMat,
                rate[iLate].s.fwdRate30360,
                rate[iLate].s.fwdAnn30360,
                rate[iLate].s.swapMat,
                rate[iLate].s.zeroRateSwap,
                numVNFMParams,
                vnfmParams,
                &swapFwdVol,
                &zeroVol,
                &zeroSwapCorr) == FAILURE) goto RETURN;       

            fwdVol = vol * swapFwdVol / swapVol;
        }
        else
        {
            fwdVol = fwdVolInput[0];
        }

        /* correlation adjustment */
        if (rate[iLate].s.start < TINY) goto RETURN;
        adj = (fwdVol*fwdVol*fwdTenor)/(vol*vol*rate[iLate].s.start);
        if (adj > 1.) goto RETURN;
        corr *= sqrt(1. - adj);
    }

    /* populate payoff structure, select payoff function      */
    /* note: ordering of ppa elements may be changed by call! */
    ppa[0] = &(pa[0]);
    ppa[1] = &(pa[1]);
    if (Q3MQBivarMCInit(
        optType,
        numPayoffParams,
        payoffParams,
        &corr,
        &(ppa[0]),
        &pf,
        &payFunc) == FAILURE) goto RETURN;

    /* Integrate payoff over density. */
    if (Q3MCPricer(
        2,
        ppa,
        &corr,
        &pf,
        payFunc,
        Q3_BV_RAN2,
        Q3_MC_NUM_PTS,
        &price) == FAILURE) goto RETURN;

    results[0] = price;
    results[1] = pa[0].fwdRate;
    results[2] = pa[1].fwdRate;

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) {
        Q3ErrMsg("%s: Failed\n", routine);  
    }
  
    return status;

} /* Q3MQBivarMCPricer */

/*-----------------------------------------------------------------------------
 * Q3MQBivarMCPricerSimple
 *
 * Multi-q bivariate option pricer.
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
    )
{
    static char routine[] = "Q3MQBivarPricerSimple";
    int         status    = FAILURE;

    FPAYOFF  *payFunc       = NULL;
    double   *smile         = NULL;
    double   *indxRates     = NULL;
    double   *vnfmParams    = NULL;
    long     *indxDates     = NULL;
    long      numVNFMParams = 0;
    char     *rateDCC       = NULL;

    PAYOFF    pf;
    MQDATA    mq[2];
    MQDATA   *ppa[2];
    double    vol;
    double    fwd;
    double    price;
    long      type;
    int       i;

    /* calibrate PA measure for each index */
    for (i = 0; i < 2; i++)
    {
        
        /* select index                         */
        /* note: ordering of rates is reversed! */
        if (i == 0)
        {
            type  = rateType1;
            vol   = sigATM1;
            smile = smile1;
            fwd   = fwdRate1;
        }
        else
        {
            type  = rateType2;
            vol   = sigATM2;
            smile = smile2;
            fwd   = fwdRate2;
        }

        /* different calibrations for rate and spread */
        switch (type)
        {
        case Q3_BV_RATE:
            if ( vol * sqrt(expiry) < Q3_MIN_VOL_CALIB) corr = 0.;

            /* smile: initialize */
            if (Q3SmileInit(
                            fwd,
                            vol,
                            0.,
                            expiry,
                            expiry,
                            1.,
                            smile,
                            &(mq[i])) == FAILURE) 
            {
                goto RETURN;
            }

            /* Calibrate MQ distribution to SV option prices. */
            if (Q3MQCalib(&(mq[i])) == FAILURE) goto RETURN;

            break;

        case Q3_BV_SPREAD:

            if ( vol * sqrt(expiry) < Q3_MIN_VOL_NORMAL) corr = 0.;

            if (Q3MQCalibSpread(
                fwd,
                vol,
                expiry,
                smile,
                &(mq[i])) == FAILURE) goto RETURN;

            break;

        default:

            Q3ErrMsg("%s: Illegal type.\n", routine);
            goto RETURN;

        }

    } /* i */
   
    /* populate payoff structure, select payoff function      */
    /* rates are being passed in this order: rate2, rate1     */
    /* note: ordering of ppa elements may be changed by call! */
    ppa[0] = &(mq[0]);
    ppa[1] = &(mq[1]);

    if (Q3MQBivarMCInit(
        optType,
        numPayoffParams,
        payoffParams,
        &corr,
        &(ppa[0]),
        &pf,
        &payFunc) == FAILURE) goto RETURN;

    /* Integrate payoff over density. */
    if (Q3MCPricer(
        2,
        ppa,
        &corr,
        &pf,
        payFunc,
        Q3_BV_RAN2,
        Q3_MC_NUM_PTS,
        &price) == FAILURE) goto RETURN;


    results[0] = price;
    results[1] = mq[0].fwdRate; /* index 1 */
    results[2] = mq[1].fwdRate; /* index 2 */

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: failed\n", routine);  
    }
  
    return status;
}

/*-----------------------------------------------------------------------------
 * Q3MQBivarMCMidCurvePricer
 *  
 * Multi-q Bivariate MidCurve Swaption Pricer
 *
 */
int Q3MQBivarMCMidCurvePricer(
    long    today,           /*  (I) Vol base date                        */
    long    valueDate,       /*  (I) Value date                           */
    long    numIndxPts,      /*  (I) Index curve num points               */
    long   *indxDates,       /*  (I) Index curve dates                    */
    double *indxRates,       /*  (I) Index curve rates                    */
    long    rateFreq0,       /*  (I) Midcurve frequency                   */
    char   *rateDCC0,        /*  (I) Midcurve DCC                         */
    long    rateFreq1,       /*  (I) 1st index curve frequency            */
    char   *rateDCC1,        /*  (I) 1st index curve DCC                  */
    long    expiryDate1,     /*  (I) 1st index reset date                 */
    long    startDate1,      /*  (I) 1st index start date                 */
    long    matDate1,        /*  (I) 1st index rate end date              */
    double *smile1,          /*  (I) 1st index smile parameters           */
    double  sigATM1,         /*  (I) 1st index ATM vol                    */
    long    rateFreq2,       /*  (I) 2nd index curve frequency            */
    char   *rateDCC2,        /*  (I) 2nd index curve DCC                  */
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
    )
{
    static char routine[] = "Q3MQBivarMCMidCurvePricer";
    int         status    = FAILURE;

    union{
        struct{
            double fwdRate;
            double fwdRate30360;
            double fwdAnn30360;
            double zeroRateSwap;
            double zeroRatePay;
            double expiry;
            double expiryVolvol;
            double start;
            double swapMat;
            double payDelay;
            double fwdAnn;    /* forward annuity to paydate                 */
            double fwdZero;   /* forward zero bond to paydate               */
            double shortTenor;/* flag: positive if matMonths < 12/rateFreq0 */
        } s;
        double a[13];
    } rate[2]; /* different structure compared with other bv-q pricers */

    FPAYOFF  *payFunc = NULL;
    double   *smile   = NULL;

    PAYOFF    pf;
    FADATA    fa[2];
    MQDATA    mq[2];
    MQDATA    pa[2];
    MQDATA   *ppa[2];
    double    sigATM;
    long      startDate, matDate, expiryDate;
    long      rateFreq;
    char     *rateDCC;
    int       i;

    double    adjustFreqDCCFactor;
    double    payPrms[18];
    int       numPayPrms = 18;
    double    outputs[2];

    /* Consistency check */
    if (expiryDate1 != expiryDate2 ||
        startDate1  != startDate2)
    {
        Q3ErrMsg("%s: failed: The midcurve pricer requires same "
                 "expiry and start dates.\n", routine);
        goto RETURN;
    }

    if (startDate1 >= matDate1 || startDate2 >= matDate2)
    {
        Q3ErrMsg("%s: failed: StartDate must be prior to matDate.\n", routine);
        goto RETURN;
    }

    if (optType == Q3_MIDCURVE_YIELD ||
        optType == Q3_MIDCURVE_RECEIVER ||
        optType == Q3_MIDCURVE_PAYER)
    {
        /* swaption payoffs: payDate is set as effective date */
        payDate = startDate1;
    }

    /* calibrate each benchmark rate */
    for (i = 0; i < 2; i++)
    {
        /* select rate */
        if (i == 0)
        {
            sigATM          = sigATM1;
            smile           = smile1;
            expiryDate      = expiryDate1;
            startDate       = startDate1;
            matDate         = matDate1;
            rateFreq        = rateFreq1;
            rateDCC         = rateDCC1;
        }
        else
        {
            sigATM          = sigATM2;
            smile           = smile2;
            expiryDate      = expiryDate2;
            startDate       = startDate2;
            matDate         = matDate2;
            rateFreq        = rateFreq2;
            rateDCC         = rateDCC2;
        }

        /* fwd dates and rates */
        if (Q3SwapRateCalc3(
            today,
            today,
            valueDate,
            numIndxPts,
            indxDates,
            indxRates,
            numDiscPts,
            discDates,
            discRates,
            rateFreq0,
            rateDCC0,
            stubType,
            stubAtEnd,
            expiryDate,
            startDate,
            matDate,
            payDate,
            (rate[i]).a) == FAILURE) goto RETURN;

        /* set correlation = 0 if total vol is too small */
        if (sigATM * sqrt((rate[i].s.expiry)) < Q3_MIN_VOL_CALIB) corr = 0.;

        /* adjust ATMVOL if freqs are different */
        if (Q3AdjustFreqDCCATMVol_BVQ(
            today,
            numDiscPts,
            discDates,
            discRates,
            startDate,
            matDate,
            rateFreq0,
            rateFreq,
            fixedType,
            rateDCC0,
            rateDCC,
            stubType,
            stubAtEnd,
            &adjustFreqDCCFactor) == FAILURE) goto RETURN;

        /* smile: initialize */
        if (Q3SmileInit(
                        (rate[i]).s.fwdRate,
                        sigATM,
                        0.,
                        (rate[i]).s.expiry,
                        (rate[i]).s.expiry,
                        adjustFreqDCCFactor,
                        smile,
                        &(mq[i])) == FAILURE) 
        {
            goto RETURN;
        }

        /* Calibrate MQ distribution to SV option prices. */
        if (Q3MQCalib(&(mq[i])) == FAILURE) goto RETURN;

        /* Initialize FA parameters. */
        (fa[i]).mq = &(mq[i]);
        if (Q3FASmileInit_Stub(
            (rate[i]).s.expiry,
            sigATM*adjustFreqDCCFactor,
            (rate[i]).s.start,
            rateFreq0,
            (rate[i]).s.swapMat,
            (rate[i]).s.fwdRate30360,
            (rate[i]).s.fwdAnn30360,
            (rate[i]).s.zeroRateSwap,
            (rate[i]).s.payDelay,
            (rate[i]).s.zeroRatePay,
            numVNFMParams,
            vnfmParams,
            (rate[i].s.shortTenor>0) ? Q3_CASH_SETL : setlType,
            /* no VNFM if tenor is too short */
            setlType,
            &(fa[i])) == FAILURE) goto RETURN;

        /* initialize PA MQ */
        if (Q3MQCopySmileFromMQ(
            &(mq[i]),
            &(pa[i])) == FAILURE) goto RETURN;

        /* set NCK parameters */
        if (Q3DecodeNCK(
            Q3_BV_NCK,
            &(pa[i])) == FAILURE) goto RETURN;

        /* compute calibration targets */
        if (Q3MQTargetFA(
            &(fa[i]),
            &(pa[i])) == FAILURE) goto RETURN;

        /* bootstrap PA multi-q measure */
        if (Q3MQBootstrapQ(&(pa[i])) == FAILURE) goto RETURN;

        /* calibrate PA */
        if (Q3MQCalib(&(pa[i])) == FAILURE) goto RETURN;

    } /* i */

    /* set parameters for MidCurve Pricer */

    payPrms[0]  = strike;
    payPrms[1]  = rate[0].s.swapMat;
    payPrms[2]  = rate[1].s.swapMat;
    payPrms[3]  = rate[0].s.payDelay;
    payPrms[4]  = rateFreq0;
    payPrms[5]  = (fa[0]).alphaAnn;
    payPrms[6]  = (fa[1]).alphaAnn;
    payPrms[7]  = (fa[0]).powerAnn;
    payPrms[8]  = (fa[1]).powerAnn;
    payPrms[9]  = (fa[0]).alphaDel;
    payPrms[10] = (fa[1]).alphaDel;
    payPrms[11] = (fa[0]).powerDel;
    payPrms[12] = (fa[1]).powerDel;
    payPrms[13] = rate[0].s.fwdZero;
    payPrms[14] = rate[1].s.fwdZero;
    payPrms[15] = rate[0].s.fwdAnn;
    payPrms[16] = rate[1].s.fwdAnn;
    payPrms[17] = ((setlType == Q3_CASH_SETL) ? 1 : 0);

    /* populate payoff structure, select payoff function      */
    ppa[0] = &(pa[0]);
    ppa[1] = &(pa[1]);

    if (Q3MQBivarMCInit(
        optType,
        numPayPrms,
        payPrms,
        &corr,
        &(ppa[0]),
        &pf,
        &payFunc) == FAILURE) goto RETURN;

    if (Q3MCPricer_MidCurve(
        2,
        ppa,
        &corr,
        &pf,
        payFunc,
        Q3_BV_RAN2,
        Q3_MC_NUM_PTS,
        &outputs[0]) == FAILURE) goto RETURN;

    results[0] = outputs[0]; /* forward premium */
    results[1] = pa[0].fwdRate; /* adjusted forward swap rate 1 */
    results[2] = pa[1].fwdRate; /* adjusted forward swap rate 2 */

   /* if optType is Q3_MIDCURVE_YIELD, the yield is returned */
    if (optType == Q3_MIDCURVE_YIELD)
        results[0] /= fabs(rate[0].s.fwdAnn - rate[1].s.fwdAnn);

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) Q3ErrMsg("%s: Failed\n", routine);

    return status;

} /* Q3MQBivarMCMidCurvePricer */

/*f----------------------------------------------------------------------------
 * Q3MQBivarMCInit
 *
 * Select payoff function and populate payoff parameters. 
 *
 */
int Q3MQBivarMCInit(
    long      optType,    /* (I) option type            */
    long      numPayPrms, /* (I) number of payoff params */
    double   *payPrms,    /* (I) payoff coefficients     */
    double   *corr,       /* (I) correlation             */
    MQDATA  **mq,         /* (I/O) measure data          */
    PAYOFF   *pf,         /* (O) option payoff structure */
    FPAYOFF **payFunc     /* (O) payoff function         */
    )
{
    static char routine[] = "Q3MQBivarMCInit";
    int         status   = FAILURE;
    long        calJoint = FALSE;
    long        calIndex = 0; /* TEST 1; */
    MQDATA     *mqTmp;
    long        i;
    double      RIB1, RIB2;
        
    /* number of payoff params */
    if (numPayPrms > Q3_MAX_PAY_PARAMS){
        Q3ErrMsg("%s: Number of payoff params (%d) exceeds max = %d.\n", 
                 routine, numPayPrms,Q3_MAX_PAY_PARAMS);
        goto RETURN;
    }

    switch (optType)
    {   

    case Q3_ADJ_FWD:

        *payFunc    = Q3Pay2D_NullYld;
        break;

    case Q3_JOINT_FWD:
        
        *payFunc = Q3Pay2D_YldYld;
        break;

    case Q3_VNL + Q3_CALL:

        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_VnlNull;
        pf->cop  = (Q3_COP_TYPE(optType) == Q3_CALL ? 1 : -1 );
        break;

    case Q3_VNL + Q3_PUT:

        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_VnlNull;
        pf->cop  = (Q3_COP_TYPE(optType) == Q3_CALL ? 1 : -1 );
        break;

    case Q3_CALL_SUM:           

        if (numPayPrms < 2)
        {
            Q3ErrMsg("%s: Payoff requires 2 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_Sum;
        pf->cop  = 1;
        break;

    case Q3_PUT_SUM:           

        if (numPayPrms < 2)
        {
            Q3ErrMsg("%s: Payoff requires 2 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_Sum;
        pf->cop  = -1;
        break;

    case Q3_CALL_PROD:           

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_Prod;
        pf->cop  = 1;
        break;

    case Q3_PUT_PROD:           

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_Prod;
        pf->cop  = -1;
        break;
    
    case Q3_CALL_PERC:           

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        pf->params[0] = payPrms[0];
        pf->params[1] = 0.; /* zero spread */
        *payFunc = Q3Pay2D_Perc;
        pf->cop  = 1;
        break;
    
    case Q3_PUT_PERC:           

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        pf->params[0] = payPrms[0];
        pf->params[1] = 0.; /* zero spread */
        *payFunc = Q3Pay2D_Perc;
        pf->cop  = -1;
        break;

    case Q3_CALL_PERC_WGT:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_PercWgt;
        pf->cop  = 1;
        break;
    
    case Q3_PUT_PERC_WGT:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_PercWgt;
        pf->cop  = -1;
        break;

    case Q3_FLR_W_FLR:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_FlrFlrOrCapCap;
        pf->cop  = 1;
        break;

    case Q3_CAP_W_CAP:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_FlrFlrOrCapCap;
        pf->cop  = -1;
        break;

    case Q3_FLR_W_CAP:           
        
        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_FlrCapOrCapFlr;
        pf->cop  = 1;
        break;
    
    case Q3_CAP_W_FLR:           
        
        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_FlrCapOrCapFlr;
        pf->cop  = -1;
        break;    
    
    case Q3_FLR_W_FLR_EMBED:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_FlrFlrOrCapCapEmbedFlt;
        pf->cop  = 1;
        break;
    
    case Q3_CAP_W_CAP_EMBED:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_FlrFlrOrCapCapEmbedFlt;
        pf->cop  = -1;
        break;

    case Q3_FLR_W_CAP_EMBED:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_FlrCapOrCapFlrEmbedFlt;
        pf->cop  = 1;
        break;
    
    case Q3_CAP_W_FLR_EMBED:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params, %d supplied.\n", 
                     routine, numPayPrms);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_FlrCapOrCapFlrEmbedFlt;
        pf->cop  = -1;
        break;

    case Q3_BS_PERC_RATE_CALL:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_Prod;
        pf->cop  = 1;
        calJoint = TRUE;
        break;

    case Q3_BS_PERC_RATE_PUT:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_Prod;
        pf->cop  = -1;
        calJoint = TRUE;
        break;

    case Q3_BS_PERC_CALL:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires at least 1 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        if ( numPayPrms == 1)
        {
            pf->params[1] = 0.;
        }
        *payFunc = Q3Pay2D_Perc;
        pf->cop  = 1;
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        calJoint = TRUE;
        calIndex = 1;
        break;

    case Q3_BS_PERC_PUT:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires at least 1 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        if ( numPayPrms == 1)
        {
            pf->params[1] = 0.;
        }
        *payFunc = Q3Pay2D_Perc;
        pf->cop  = -1;
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        calJoint = TRUE;
        calIndex = 1;
        break;

    case Q3_BS_SPRD_RATE_CALL:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->params[1] = 1.;         /* leverage */
        *payFunc      = Q3Pay2D_Sum;
        pf->cop       = 1;
        break;
        
    case Q3_BS_SPRD_RATE_PUT:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->params[1] = 1.;         /* leverage */
        *payFunc      = Q3Pay2D_Sum;
        pf->cop       = -1;
        break;

    case Q3_BS_PERC_SPRD_CALL:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_VnlNull; /* Vanilla on 1st variable */
        pf->cop  = 1.;
        calJoint = TRUE;
        calIndex = 0;
        break;

    case Q3_BS_PERC_SPRD_PUT:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_VnlNull; /* Vanilla on 1st variable */
        pf->cop  = -1.;
        calJoint = TRUE;
        calIndex = 0;
        break;

    case Q3_BS_SPRD_DIG_RIBIN:
        if (numPayPrms < 2)
        {
            Q3ErrMsg("%s: Payoff requires 2 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        /* payFunc expects leps, heps, weight1, weight2 */
        pf->params[2] = 0.;     
        pf->params[3] = 0.;
        pf->params[4] = 1.;
        pf->params[5] = 1.;
        /* payoff function expects order: rate1, rate2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
        *payFunc      = Q3Pay2D_DigSpdInRIB_EPS;
        break;

    case Q3_BS_SPRD_DIG_RIBOUT:
        if (numPayPrms < 2)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        /* payFunc expects leps, heps, weight1, weight2 */
        pf->params[2] = 0.;     
        pf->params[3] = 0.;
        pf->params[4] = 1.;
        pf->params[5] = 1.;
        /* payoff function expects order: rate1, rate2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
        *payFunc      = Q3Pay2D_DigSpdOutRIB_EPS;
        break;

    case Q3_BS_SPRD_DIG_RIBIN_EPS:
        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        /* payFunc expects weight1, weight2 */
        pf->params[4] = 1.;
        pf->params[5] = 1.;
        /* payoff function expects order: rate1, rate2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
        *payFunc      = Q3Pay2D_DigSpdInRIB_EPS;

        RIB1 = pf->params[1] - pf->params[0];
        if(RIB1 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[2] < 0 && -pf->params[2] > RIB1) ||
            (pf->params[3] < 0 && -pf->params[3] > RIB1))
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "leps (or heps) < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        break;

    case Q3_BS_SPRD_DIG_RIBOUT_EPS:
        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 1 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        /* payFunc expects weight1, weight2 */
        pf->params[4] = 1.;
        pf->params[5] = 1.;
        /* payoff function expects order: rate1, rate2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
        *payFunc      = Q3Pay2D_DigSpdOutRIB_EPS;
        
        RIB1 = pf->params[1] - pf->params[0];
        if(RIB1 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[2] < 0 && -pf->params[2] > RIB1) ||
            (pf->params[3] < 0 && -pf->params[3] > RIB1))
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "leps (or heps) < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        break;

    case Q3_BS_PERC_DIG_RIBIN:
        if (numPayPrms < 2)
        {
            Q3ErrMsg("%s: Payoff requires 2 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        /* payFunc expects leps, heps */
        pf->params[2] = 0.;     
        pf->params[3] = 0.;
        *payFunc      = Q3Pay2D_DigPercInRIB_EPS;
        /* payoff function expects order: rate1, rate2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
        calJoint = TRUE;
        calIndex = 0;
        break;

    case Q3_BS_PERC_DIG_RIBOUT:
        if (numPayPrms < 2)
        {
            Q3ErrMsg("%s: Payoff requires 2 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        /* payFunc expects leps, heps */
        pf->params[2] = 0.;     
        pf->params[3] = 0.;
        *payFunc      = Q3Pay1D_DigPercOutRIB_EPS;
        /* payoff function expects order: rate1, rate2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
        calJoint = TRUE;
        calIndex = 0;
        break;

    case Q3_BS_PERC_DIG_RIBIN_EPS:
        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc      = Q3Pay2D_DigPercInRIB_EPS;
        /* payoff function expects order: rate1, rate2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
        calJoint = TRUE;
        calIndex = 0;

        RIB1 = pf->params[1] - pf->params[0];
        if(RIB1 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[2] < 0 && -pf->params[2] > RIB1) ||
            (pf->params[3] < 0 && -pf->params[3] > RIB1))
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "leps (or heps) < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        break;

    case Q3_BS_PERC_DIG_RIBOUT_EPS:
        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc      = Q3Pay2D_DigPercOutRIB_EPS;
        /* payoff function expects order: rate1, rate2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
        calJoint = TRUE;
        calIndex = 0;

        RIB1 = pf->params[1] - pf->params[0];
        if(RIB1 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[2] < 0 && -pf->params[2] > RIB1) ||
            (pf->params[3] < 0 && -pf->params[3] > RIB1))
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "leps (or heps) < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        break;

    case Q3_IN_BIRIB:

        if (numPayPrms < 6)
        {
            Q3ErrMsg("%s: Payoff requires 6 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->params[6] = 0.; /* payFunc expects epsilon params */
        pf->params[7] = 0.; 
        *payFunc =  Q3Pay2D_MinMaxIn;
        break;

    case Q3_OUT_BIRIB:
        
        if (numPayPrms < 6)
        {
            Q3ErrMsg("%s: Payoff requires 6 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->params[6] = 0.; /* payFunc expects epsilon params */
        pf->params[7] = 0.; 
        *payFunc = Q3Pay2D_MinMaxOut;
        break;

    case Q3_IN_BIRIB_EPS:
        
        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        if( RIB1 < 0 )
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB1) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        *payFunc = Q3Pay2D_MinMaxIn;
        break;
            
    case Q3_OUT_BIRIB_EPS:
        
        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        if( RIB1 < 0 )
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB1) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        *payFunc = Q3Pay2D_MinMaxOut;
        break;

    case Q3_IN_SPDRIB:
        
        if (numPayPrms < 6)
        {
            Q3ErrMsg("%s: Payoff requires 6 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->params[6] = 0.; /* payFunc expects epsilon params */
        pf->params[7] = 0.; 
        *payFunc = Q3Pay2D_MinMaxSpdIn;
        break;

    case Q3_IN_SPDRIB_EPS:
        
        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        if( RIB1 < 0 )
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB1) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        *payFunc = Q3Pay2D_MinMaxSpdIn;
        break;

    case Q3_OUT_SPDRIB:
        
        if (numPayPrms < 6)
        {
            Q3ErrMsg("%s: Payoff requires 6 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->params[6] = 0.; /* payFunc expects epsilon params */
        pf->params[7] = 0.; 
        *payFunc = Q3Pay2D_MinMaxSpdOut;
        break;

    case Q3_OUT_SPDRIB_EPS:
        
        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        if( RIB1 < 0 )
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB1) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        *payFunc = Q3Pay2D_MinMaxSpdOut;
        break;
    case Q3_IN_AND_IN:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay2D_InAndIn;
        break;

    case Q3_IN_AND_IN_EPS:

        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        RIB2 = pf->params[3] - pf->params[2];
        if( RIB1 < 0 || RIB2 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        RIB2 /= 2.0;
        if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
            (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
            (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB2) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
            *payFunc = Q3Pay2D_InAndIn;
            break;

    case Q3_IN_OR_IN:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay2D_InOrIn;
        break;

    case Q3_IN_OR_IN_EPS:

        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        RIB2 = pf->params[3] - pf->params[2];
        if( RIB1 < 0 || RIB2 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        RIB2 /= 2.0;
        if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
            (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
            (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB2) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        *payFunc = Q3Pay2D_InOrIn;
        break;

    case Q3_IN_AND_OUT:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay2D_InAndOut;
        break;

    case Q3_IN_AND_OUT_EPS:

        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        RIB2 = pf->params[3] - pf->params[2];
        if( RIB1 < 0 || RIB2 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        RIB2 /= 2.0;
        if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
            (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
            (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB2) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        /* Ignore small epsilons */
        *payFunc = Q3Pay2D_InAndOut;
        break;

    case Q3_IN_OR_OUT:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay2D_InOrOut;
        break;

    case Q3_IN_OR_OUT_EPS:

        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        RIB2 = pf->params[3] - pf->params[2];
        if( RIB1 < 0 || RIB2 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        RIB2 /= 2.0;
        if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
            (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
            (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB2) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        *payFunc = Q3Pay2D_InOrOut;
        break;

    case Q3_OUT_AND_IN:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay2D_OutAndIn;
        break;

    case Q3_OUT_AND_IN_EPS:

        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        RIB2 = pf->params[3] - pf->params[2];
        if( RIB1 < 0 || RIB2 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        RIB2 /= 2.0;
        if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
            (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
            (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB2) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        *payFunc = Q3Pay2D_OutAndIn;
        break;

    case Q3_OUT_OR_IN:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay2D_OutOrIn;
        break;

    case Q3_OUT_OR_IN_EPS:

        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        RIB2 = pf->params[3] - pf->params[2];
        if( RIB1 < 0 || RIB2 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        RIB2 /= 2.0;
        if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
            (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
            (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB2) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        *payFunc = Q3Pay2D_OutOrIn;
        break;

    case Q3_OUT_AND_OUT:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay2D_OutAndOut;
        break;

    case Q3_OUT_AND_OUT_EPS:

        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        RIB2 = pf->params[3] - pf->params[2];
        if( RIB1 < 0 || RIB2 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        RIB2 /= 2.0;
        if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
            (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
            (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB2) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        *payFunc = Q3Pay2D_OutAndOut;
        break;

    case Q3_OUT_OR_OUT:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay2D_OutOrOut;
        break;

    case Q3_OUT_OR_OUT_EPS:

        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        RIB2 = pf->params[3] - pf->params[2];
        if( RIB1 < 0 || RIB2 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        RIB2 /= 2.0;
        if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
            (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
            (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB2) )
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "epsilon < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        *payFunc = Q3Pay2D_OutOrOut;
        break;

    case Q3_COMPLEX_SPD:

        if (numPayPrms < 19)
        {
            Q3ErrMsg("%s: Payoff requires 19 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay2D_ComplexSPD;
        break;

    case Q3_MIDCURVE_YIELD:

        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->cop  = 0;
        *payFunc = Q3Pay2D_MidCurve;
        break;

    case Q3_MIDCURVE_PAYER:

        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->cop  = -1;
        *payFunc = Q3Pay2D_MidCurve;
        break;

    case Q3_MIDCURVE_RECEIVER:

        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->cop  = 1;
        *payFunc = Q3Pay2D_MidCurve;
        break;

    case Q3_MIDCURVE_CMS_YIELD:

        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->cop  = 0;
        *payFunc = Q3Pay2D_MidCurve_CMS;
        break;

    case Q3_MIDCURVE_CMS_CAP:

        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->cop  = 1;
        *payFunc = Q3Pay2D_MidCurve_CMS;
        break;

    case Q3_MIDCURVE_CMS_FLOOR: 

        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->cop  = -1;
        *payFunc = Q3Pay2D_MidCurve_CMS;
        break;

    default: 
        
        Q3ErrMsg("%s: Unrecognized option type: %d\n", routine, optType);
        goto RETURN;

    }

    /* joint calibration */
    if (calJoint == TRUE)
    { 
        /* if no user calibration parameter, fwd cal. only */
        if (numPayPrms > 2 && fabs(payPrms[2] - Q3_BS_CAL_VOL)< TINY)
        {
            if(Q3MQCalibJoint(
                corr, 
                calIndex, 
                mq) == FAILURE) goto RETURN;
        }
        else
        {
            if(Q3MQCalibJointFwd(
                (mq[0])->fwdRate * (mq[1])->fwdRate,
                corr, 
                calIndex, 
                mq) == FAILURE) goto RETURN;
        }
    }

    /* fill in payoff parameters */
    pf->optType         = optType;
    pf->corr            = corr[0];
    (pf->mq)[0]         = mq[0];
    (pf->mq)[1]         = mq[1];
        
    status = SUCCESS;

 RETURN:

    if (status == FAILURE) {
        Q3ErrMsg("%s: Failed\n", routine);
    }
 
    return status;
} /* Q3MQBivarMCInit */
