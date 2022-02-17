/******************************************************************************
 * Module:      
 * Submodule:
 * File:      mqbv.c
 * Function:    
 * Author:        Interest Rates DR
 * Revision:  $Header$
 *****************************************************************************/
#include <math.h>
#include <ctype.h>                      
#include <stdio.h>

#include "aliberr.h"
#include "aliblog.h"
#include "mqbvl.h"
#include "mqbv.h"


/* Version Info:  U_VERSION is expected to be defined when compiling - see bv_q_defs.mk */
#define _STRINGIFY(x) #x
#define _Q3_BV_VERS(v) "Q3BV VERSION " _STRINGIFY(v) " COMPILED " __DATE__ " " __TIME__
#define Q3_BV_VERS _Q3_BV_VERS(U_VERSION)
 
#define Q3_BV_NCK  2000000769 /* Calibration control data */

#define Q3_BV_RATE   0
#define Q3_BV_SPREAD 1

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
    double *payoffParamsL,    /* 25 (I) payment parameters                */
    long   *traceL,           /* 26 (I) trace                             */
    double *outputsL          /* 27 (O) premium                           */
    )
{
    static char routine[] = "Q3MQBivarPricerL";
    int status = FAILURE;

    double  smileSmile1L[20], smileSmile2L[20];

    memcpy(smileSmile1L+2, smile1L+1, ((long) smile1L[0])*sizeof(double));
    smileSmile1L[0] = smile1L[0] + 1.;
    smileSmile1L[1] = 10.;

    memcpy(smileSmile2L+2, smile2L+1, ((long) smile2L[0])*sizeof(double));
    smileSmile2L[0] = smile2L[0] + 1.;
    smileSmile2L[1] = 10.;

    if(Q3SmileBivarPricerL(
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

} /* Q3MQBivarPricerL */

/*-----------------------------------------------------------------------------
 * Q3MQBivarBusPricerL
 *
 * Multi-q bivariate option pricer add-in.
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
    double *payoffParamsL,    /* 25 (I) payment parameters                */
    char   *holidayFileL,     /* 26 (I) Holiday file name         */
    long   *traceL,           /* 27 (I) trace                             */
    double *outputsL          /* 29 (O) premium                           */
    )
{
    static char routine[] = "Q3MQBivarPricerL";
    int status = FAILURE;

    double  smileSmile1L[20], smileSmile2L[20];

    memcpy(smileSmile1L+2, smile1L+1, ((long) smile1L[0])*sizeof(double));
    smileSmile1L[0] = smile1L[0] + 1.;
    smileSmile1L[1] = 10.;

    memcpy(smileSmile2L+2, smile2L+1, ((long) smile2L[0])*sizeof(double));
    smileSmile2L[0] = smile2L[0] + 1.;
    smileSmile2L[1] = 10.;

    if(Q3SmileBivarBusPricerL(
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
                           holidayFileL,
                           traceL,
                           outputsL) == FAILURE)
         goto RETURN;

    status = SUCCESS;
    
  RETURN:

    if (status == FAILURE) {
        Q3ErrMsg("%s: Failed\n", routine);
    }
  
    return status;

} /* Q3MQBivarBusPricerL */

/*-----------------------------------------------------------------------------
 * Q3SmileBivarPricerL
 *
 * Multi-q bivariate option pricer add-in.
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
    double *payoffParamsL,    /* 25 (I) payment parameters                */
    long   *traceL,           /* 26 (I) trace                             */
    double *outputsL          /* 27 (O) premium                           */
    )
{
    static char routine[] = "Q3SmileBivarPricerL";

    int status = FAILURE;

    long expiryDate1L[2] = {1, 0};
    long expiryDate2L[2] = {1, 0};
    long startDate1L[2]  = {1, 0};
    long startDate2L[2]  = {1, 0};
    long matDate1L[2]    = {1, 0};
    long matDate2L[2]    = {1, 0};
    int  i;

    double smile1Input[13], smile2Input[13];

    char   *HolidayVol  = NULL;
    char   *BusVolDCC   = NULL;

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

    memcpy(smile1Input, smile1L+1, ((long) smile1L[0])*sizeof(double));

    for (i = (long) smile1L[0]; i<13; i++)
        smile1Input[i] = 0.;

    memcpy(smile2Input, smile2L+1, ((long) smile2L[0])*sizeof(double));

    for (i = (long) smile2L[0]; i<13; i++)
        smile2Input[i] = 0.;

    if (Q3MQBivarPricer(
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
        HolidayVol,
        BusVolDCC,
        outputsL+1) == FAILURE) goto RETURN;  

    status = SUCCESS;

 done: /* for alib */
 RETURN:

    if (status == FAILURE)
    {
        Q3ErrMsg("%s: failed.\n", routine);
    }
  
    return status;

} /* Q3SmileBivarPricerL */

/*-----------------------------------------------------------------------------
 * Q3SmileBivarBusPricerL
 *
 * Multi-q bivariate option pricer add-in.
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
    double *payoffParamsL,    /* 25 (I) payment parameters                */
    char   *holidayFileL,     /* 26 (I) Holiday file name and DCC         */
    long   *traceL,           /* 27 (I) trace                             */
    double *outputsL          /* 28 (O) premium                           */
    )
{
    static char routine[] = "Q3SmileBivarPricerL";

    int status = FAILURE;

    long expiryDate1L[2] = {1, 0};
    long expiryDate2L[2] = {1, 0};
    long startDate1L[2]  = {1, 0};
    long startDate2L[2]  = {1, 0};
    long matDate1L[2]    = {1, 0};
    long matDate2L[2]    = {1, 0};
    int  i;

    double smile1Input[13], smile2Input[13];

    char   *BusVolDCC  = NULL;
    char   *HolidayVol = NULL;

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
            Q3_CHAR_BLOCK_L, holidayFileL,  "holidayFile",       
            Q3_LONG_L,       traceL,        "trace",
            0);
    }  

    memcpy(smile1Input, smile1L+1, ((long) smile1L[0])*sizeof(double));

    for (i = (long) smile1L[0]; i<13; i++)
        smile1Input[i] = 0.;

    memcpy(smile2Input, smile2L+1, ((long) smile2L[0])*sizeof(double));

    for (i = (long) smile2L[0]; i<13; i++)
        smile2Input[i] = 0.;

    switch (GTO_ARG_SIZE(holidayFileL)){
    case 0:
        break;
    case 1:
        HolidayVol = holidayFileL + WRAP_STR_IDX(1);
        break;
    case 2:
        HolidayVol = holidayFileL + WRAP_STR_IDX(1);
        BusVolDCC = holidayFileL + WRAP_STR_IDX(2);
        break;
    default:
        Q3ErrMsg("%s: Too many inputs '%s'. \n", 
                routine, 
                holidayFileL + WRAP_STR_IDX(3));
        goto RETURN;
    }

    if (Q3MQBivarPricer(
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
        HolidayVol,
        BusVolDCC,
        outputsL+1) == FAILURE) goto RETURN;  

    status = SUCCESS;

 done: /* for alib */
 RETURN:

    if (status == FAILURE)
    {
        Q3ErrMsg("%s: failed.\n", routine);
    }
  
    return status;

} /* Q3SmileBivarBusPricerL */

/*-----------------------------------------------------------------------------
 * Q3MQBivarPricerSimpleL
 * 
 * Multi-q bivariate option pricer with simplified interface
 *
 */
int Q3MQBivarPricerSimpleL(
    double *expiryL,          /*  1 (I) exirpation in years              */
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
    static char routine[] = "Q3MQBivarPricerSimpleL";
    int status = FAILURE;

    double  smileSmile1L[20], smileSmile2L[20];

    memcpy(smileSmile1L+2, smile1L+1, ((long) smile1L[0])*sizeof(double));
    smileSmile1L[0] = smile1L[0] + 1.;
    smileSmile1L[1] = 10.;

    memcpy(smileSmile2L+2, smile2L+1, ((long) smile2L[0])*sizeof(double));
    smileSmile2L[0] = smile2L[0] + 1.;
    smileSmile2L[1] = 10.;

    if(Q3SmileBivarPricerSimpleL(
                      expiryL,
                      rateType1L,
                      fwdRate1L,
                      smileSmile1L,
                      sigATM1L,
                      rateType2L,
                      fwdRate2L,
                      smileSmile2L,
                      sigATM2L,
                      corrL,
                      optTypeL,
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

} /* Q3MQBivarPricerSimpleL */


/*-----------------------------------------------------------------------------
 * Q3SmileBivarPricerSimpleL
 * 
 * Multi-q bivariate option pricer with simplified interface
 *
 */
int Q3SmileBivarPricerSimpleL(
    double *expiryL,          /*  1 (I) exirpation in years              */
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
    static char routine[] = "Q3SmileBivarPricerSimpleL";

    int status = FAILURE;
    int i;

    double smile1Input[13], smile2Input[13];
    
    /* initialize alib callbacks */
    Q3ToAlibErrInit();
    
    /* initialize global error buffer */
    Q3BivarErrBuffInit();

    /* Test input */
    GTO_WRAP_CHECK_SCALAR(expiryL,              "expiry");
    GTO_WRAP_CHECK_SCALAR(rateType1L,           "rateType1");
    GTO_WRAP_CHECK_SCALAR(fwdRate1L,            "fwdRate1");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile1L,     "smile1", 4);
    GTO_WRAP_CHECK_SCALAR(sigATM1L,             "sigATM1");
    GTO_WRAP_CHECK_SCALAR(rateType2L,           "rateType2");
    GTO_WRAP_CHECK_SCALAR(fwdRate2L,            "fwdRate2");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile2L,     "smile2", 4);
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

    memcpy(smile1Input, smile1L+1, ((long) smile1L[0])*sizeof(double));

    for (i = (long) smile1L[0]; i<13; i++)
        smile1Input[i] = 0.;

    memcpy(smile2Input, smile2L+1, ((long) smile2L[0])*sizeof(double));

    for (i = (long) smile2L[0]; i<13; i++)
        smile2Input[i] = 0.;

    if (Q3MQBivarPricerSimple(
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

} /* Q3SmileBivarPricerSimpleL */

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
    long   *optTypeL,         /*  9 (I) Option type                      */
    double *strikeL,          /* 10 (I) strike                           */
    double *payoffParamsL,    /* 11 (I) payoff Parameters                */
    double *initialGuessL,    /* 12 (I) initial guess of corr            */
    long   *traceL,           /* 13 (I) Trace                            */
    double *outputsL          /* 14 (O) Pricing results                  */
    )
{
    static char routine[] = "Q3SmileBivarSimpleCorrCalibL";

    int status = FAILURE;
    int i;

    double smileSmile1L[20], smileSmile2L[20];
    double smile1Input[13], smile2Input[13];
    
    /* initialize alib callbacks */
    Q3ToAlibErrInit();
    
    /* initialize global error buffer */
    Q3BivarErrBuffInit();

    /* Test input */
    GTO_WRAP_CHECK_SCALAR(premiumL,             "premium");
    GTO_WRAP_CHECK_SCALAR(expiryL,              "expiry");
    GTO_WRAP_CHECK_SCALAR(fwdRate1L,            "fwdRate1");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile1L,     "smile1", 4);
    GTO_WRAP_CHECK_SCALAR(sigATM1L,             "sigATM1");
    GTO_WRAP_CHECK_SCALAR(fwdRate2L,            "fwdRate2");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile2L,     "smile2", 4);
    GTO_WRAP_CHECK_SCALAR(sigATM2L,             "sigATM2");
    GTO_WRAP_CHECK_SCALAR(optTypeL,             "optType");
    GTO_WRAP_CHECK_SCALAR(strikeL,              "strike");
    GTO_WRAP_CHECK_VECTOR_LEN(payoffParamsL,    "payoffParams",2);
    GTO_WRAP_CHECK_SCALAR(initialGuessL,        "initialGuess");
    GTO_WRAP_CHECK_VECTOR_LEN(outputsL,         "outputs", 3);
    
    /*  Alib RegTest logging if tracing is on */
    if (traceL[1] > 0) 
    {
        Q3LilVectLoggingFile(
            "q3.log",
            (traceL[1] == 1 ? "w" : "a"),
            routine,
            Q3_DOUBLE_L,     premiumL,      "premium",
            Q3_DOUBLE_L,     expiryL,       "expiry",
            Q3_DOUBLE_L,     fwdRate1L,     "fwdRate1",
            Q3_DOUBLE_L,     smile1L,       "smile1",
            Q3_DOUBLE_L,     sigATM1L,      "sigATM1",
            Q3_DOUBLE_L,     fwdRate2L,     "fwdRate2",
            Q3_DOUBLE_L,     smile2L,       "smile2",
            Q3_DOUBLE_L,     sigATM2L,      "sigATM2",
            Q3_LONG_L,       optTypeL,      "optType",
            Q3_DOUBLE_L,     strikeL,       "strike",
            Q3_DOUBLE_L,     payoffParamsL, "payoffParams",
            Q3_DOUBLE_L,     initialGuessL, "initialGuess",        
            Q3_LONG_L,       traceL,        "trace",
            0);
    }  

    memcpy(smileSmile1L+2, smile1L+1, ((long) smile1L[0])*sizeof(double));
    smileSmile1L[0] = smile1L[0] + 1.;
    smileSmile1L[1] = 10.;

    memcpy(smileSmile2L+2, smile2L+1, ((long) smile2L[0])*sizeof(double));
    smileSmile2L[0] = smile2L[0] + 1.;
    smileSmile2L[1] = 10.;

    memcpy(smile1Input, smileSmile1L+1, ((long) smileSmile1L[0])*sizeof(double));

    for (i = (long) smileSmile1L[0]; i<13; i++)
        smile1Input[i] = 0.;

    memcpy(smile2Input, smileSmile2L+1, ((long) smileSmile2L[0])*sizeof(double));

    for (i = (long) smileSmile2L[0]; i<13; i++)
        smile2Input[i] = 0.;

    if (Q3SmileBivarSimpleCorrCalib(
        premiumL[1],
        expiryL[1], 
        fwdRate1L[1],
        smile1Input,  
        sigATM1L[1],    
        fwdRate2L[1],
        smile2Input,  
        sigATM2L[1],   
        optTypeL[1],       
        strikeL[1],
        payoffParamsL+1,
        initialGuessL[1],
        outputsL+1) == FAILURE) goto RETURN;  
    status = SUCCESS;

 done: /* for alib */
 RETURN:

    if (status == FAILURE)
    {
        Q3ErrMsg("%s: failed.\n", routine);
    }
  
    return status;
}/*Q3SmileBivarSimpleCorrCalibL*/

/*-----------------------------------------------------------------------------
 * Q3MQBivarCorrL
 *
 * Multi-q bivariate correlation calculator.
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
    long   *payDateL,         /* 22 (I) Option payment date               */
    long   *setlTypeL,        /* 23 (I) Cash or physical settlement.      */
    long   *traceL,           /* 24 (I) trace                             */
    double *outputsL          /* 25 (O) premium                           */
    )
{
    static char routine[] = "Q3MQBivarCorrL";
    int status = FAILURE;

    double  smileSmile1L[20], smileSmile2L[20];

    memcpy(smileSmile1L+2, smile1L+1, ((long) smile1L[0])*sizeof(double));
    smileSmile1L[0] = smile1L[0] + 1.;
    smileSmile1L[1] = 10.;

    memcpy(smileSmile2L+2, smile2L+1, ((long) smile2L[0])*sizeof(double));
    smileSmile2L[0] = smile2L[0] + 1.;
    smileSmile2L[1] = 10.;

    if(Q3SmileBivarCorrL(
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
                         payDateL,
                         setlTypeL,
                         traceL,
                         outputsL) == FAILURE)
         goto RETURN;

    status = SUCCESS;
    
  RETURN:

    if (status == FAILURE) {
        Q3ErrMsg("%s: Failed\n", routine);
    }
  

    return status;

} /* Q3MQBivarCorrL */

/*-----------------------------------------------------------------------------
 * Q3MQBivarBusCorrL
 *
 * Multi-q bivariate correlation calculator.
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
    long   *traceL,           /* 25 (I) trace                             */
    double *outputsL          /* 26 (O) premium                           */
    )
{
    static char routine[] = "Q3MQBivarCorrL";
    int status = FAILURE;

    double  smileSmile1L[20], smileSmile2L[20];

    memcpy(smileSmile1L+2, smile1L+1, ((long) smile1L[0])*sizeof(double));
    smileSmile1L[0] = smile1L[0] + 1.;
    smileSmile1L[1] = 10.;

    memcpy(smileSmile2L+2, smile2L+1, ((long) smile2L[0])*sizeof(double));
    smileSmile2L[0] = smile2L[0] + 1.;
    smileSmile2L[1] = 10.;

    if(Q3SmileBivarBusCorrL(
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
                         payDateL,
                         setlTypeL,
                         holidayFileL,
                         traceL,
                         outputsL) == FAILURE)
         goto RETURN;

    status = SUCCESS;
    
  RETURN:

    if (status == FAILURE) {
        Q3ErrMsg("%s: Failed\n", routine);
    }
  

    return status;

} /* Q3MQBivarBusCorrL */

/*-----------------------------------------------------------------------------
 * Q3SmileBivarCorrL
 *
 * Multi-q bivariate correlation calculator.
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
    long   *payDateL,         /* 22 (I) Option payment date               */
    long   *setlTypeL,        /* 23 (I) Cash or physical settlement.      */
    long   *traceL,           /* 24 (I) trace                             */
    double *outputsL          /* 25 (O) premium                           */
    )
{
    static char routine[] = "Q3SmileBivarCorrL";

    int status = FAILURE;
    int i;

    long expiryDate1L[2] = {1, 0};
    long expiryDate2L[2] = {1, 0};
    long startDate1L[2]  = {1, 0};
    long startDate2L[2]  = {1, 0};
    long matDate1L[2]    = {1, 0};
    long matDate2L[2]    = {1, 0};
    
    double smile1Input[13], smile2Input[13];

    char   *BusVolDCC  = NULL;
    char   *HolidayVol = NULL;

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
    GTO_WRAP_CHECK_SCALAR(payDateL,             "payment date");
    GTO_WRAP_CHECK_SCALAR(setlTypeL,            "settlement type");
    GTO_WRAP_CHECK_VECTOR_LEN(outputsL,         "outputs", 1);
    
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
            Q3_TDATE_L,      payDateL,      "payDate",
            Q3_LONG_L,       setlTypeL,     "setlType",
            Q3_LONG_L,       traceL,        "trace",
            0);
    }  

    memcpy(smile1Input, smile1L+1, ((long) smile1L[0])*sizeof(double));

    for (i = (long) smile1L[0]; i<13; i++)
        smile1Input[i] = 0.;

    memcpy(smile2Input, smile2L+1, ((long) smile2L[0])*sizeof(double));

    for (i = (long) smile2L[0]; i<13; i++)
        smile2Input[i] = 0.;

    if (Q3MQBivarCorr(
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
        corrL+1,          
        (corrL[0]>1 ? corrL+2 : NULL),
        payDateL[1],       
        setlTypeL[1],     
        HolidayVol,
        BusVolDCC,
        outputsL+1) == FAILURE) goto RETURN;  

    status = SUCCESS;

 done: /* for alib */
 RETURN:

    if (status == FAILURE)
    {
        Q3ErrMsg("%s: failed.\n", routine);
    }
  
    return status;

} /* Q3SmileBivarCorrL */

/*-----------------------------------------------------------------------------
 * Q3SmileBivarBusCorrL
 *
 * Multi-q bivariate correlation calculator.
 *
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
    long   *payDateL,         /* 22 (I) Option payment date               */
    long   *setlTypeL,        /* 23 (I) Cash or physical settlement.      */
    char   *holidayFileL,     /* 24 (I) Holiday file name and DCC         */
    long   *traceL,           /* 25 (I) trace                             */
    double *outputsL          /* 26 (O) premium                           */
    )
{
    static char routine[] = "Q3SmileBivarBusCorrL";

    int status = FAILURE;
    int i;

    long expiryDate1L[2] = {1, 0};
    long expiryDate2L[2] = {1, 0};
    long startDate1L[2]  = {1, 0};
    long startDate2L[2]  = {1, 0};
    long matDate1L[2]    = {1, 0};
    long matDate2L[2]    = {1, 0};
    
    double smile1Input[13], smile2Input[13];

    char   *BusVolDCC  = NULL;
    char   *HolidayVol = NULL;

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
    GTO_WRAP_CHECK_SCALAR(payDateL,             "payment date");
    GTO_WRAP_CHECK_SCALAR(setlTypeL,            "settlement type");
    GTO_WRAP_CHECK_VECTOR_LEN(outputsL,         "outputs", 1);
    
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
            Q3_TDATE_L,      payDateL,      "payDate",
            Q3_LONG_L,       setlTypeL,     "setlType",
            Q3_CHAR_BLOCK_L, holidayFileL,  "holidayFile",       
            Q3_LONG_L,       traceL,        "trace",
            0);
    }  

    memcpy(smile1Input, smile1L+1, ((long) smile1L[0])*sizeof(double));

    for (i = (long) smile1L[0]; i<13; i++)
        smile1Input[i] = 0.;

    memcpy(smile2Input, smile2L+1, ((long) smile2L[0])*sizeof(double));

    for (i = (long) smile2L[0]; i<13; i++)
        smile2Input[i] = 0.;

    switch (GTO_ARG_SIZE(holidayFileL)){
    case 0:
        break;
    case 1:
        HolidayVol = holidayFileL + WRAP_STR_IDX(1);
        break;
    case 2:
        HolidayVol = holidayFileL + WRAP_STR_IDX(1);
        BusVolDCC = holidayFileL + WRAP_STR_IDX(2);
        break;
    default:
        Q3ErrMsg("%s: Too many inputs '%s'. \n", 
                routine, 
                holidayFileL + WRAP_STR_IDX(3));
        goto RETURN;
    }

    if (Q3MQBivarCorr(
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
        corrL+1,          
        (corrL[0]>1 ? corrL+2 : NULL),
        payDateL[1],       
        setlTypeL[1],     
        HolidayVol,
        BusVolDCC,
        outputsL+1) == FAILURE) goto RETURN;  

    status = SUCCESS;

 done: /* for alib */
 RETURN:

    if (status == FAILURE)
    {
        Q3ErrMsg("%s: failed.\n", routine);
    }
  
    return status;

} /* Q3SmileBivarBusCorrL */

/*-----------------------------------------------------------------------------
 * Q3MQBivarPricer
 * 
 * Multi-q bivariate option pricer.
 *
 */
int Q3MQBivarPricer(
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
    char   *holidayfile,     /* (I) holiday file name                    */
    char   *BusVolDCC,       /* (I) BUS/251F or BUS/BUS                  */
    double *results          /* (O) Pricing results.                     */
    )
{
    static char routine[] = "Q3MQBivarPricer";
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

    /* CMS BMA option */
    if (optType == Q3_BS_PERC_CALL ||
        optType == Q3_BS_PERC_PUT  ||
        optType == Q3_BS_PERC_RATE_CALL ||
        optType == Q3_BS_PERC_RATE_PUT ||
        optType == Q3_BS_JOINT_FWD)
    {
        /* Consistency check */
        if (expiryDate1 != expiryDate2 ||
            startDate1  != startDate2 ||
            matDate1    != matDate2) 
        {
            Q3ErrMsg("%s: failed: The percentage spread and the rate require "
                     "same expiry, start & maturity dates.\n", routine);  
            goto RETURN;
        }


        if (Q3MQPSAPricer(
            today,           /*  1 (I) Vol base date                        */
            valueDate,       /*  2 (I) Value date                           */
            rateFreq2,       /*  3 (I) BMA and CMS fixed leg frequency      */
            rateDCC2,        /*  4 (I) BMA and CMS fixed Leg DCC            */
            rateFreq1,       /*  5 (I) BMA float leg frequency              */
            rateDCC1,        /*  6 (I) BMA float leg DCC                    */
            numIndxPts2,     /*  7 (I) LIBOR zero curve num points          */
            indxDates2,      /*  8 (I) LIBOR zero curve dates               */
            indxRates2,      /*  9 (I) LIBOR zero curve rates               */
            numIndxPts1,     /* 10 (I) basis zero curve num points          */
            indxDates1,      /* 11 (I) basis zero curve dates               */
            indxRates1,      /* 12 (I) basis zero curve rates               */
            expiryDate1,     /* 13 (I) reset date                           */
            startDate1,      /* 14 (I) start date                           */
            matDate1,        /* 15 (I) end date                             */
            smile2,          /* 16 (I) CMS smile parameters                 */
            sigATM2,         /* 17 (I) CMS ATM vol                          */
            numVNFMParams2,  /* 18 (I) 1st index num VNFM parameters        */
            vnfmParams2,     /* 19 (I) 1st index VNFM parameters            */
            smile1,          /* 20 (I) spread smile parameters              */
            sigATM1,         /* 21 (I) spread ATM vol                       */
            numDiscPts,      /* 22 (I) Discount curve num points            */
            discDates,       /* 23 (I) Discount curve dates                 */
            discRates,       /* 24 (I) Discount curve rates                 */
            corr,            /* 25 (I) Index correlation data               */
            optType,         /* 26 (I) Option type                          */
            payDate,         /* 27 (I) Option payment date                  */
            setlType,        /* 28 (I) Cash or physical settlement          */
            numPayoffParams, /* 29 (I) Number of payment parameters         */
            payoffParams,    /* 28 (I) Payment parameters                   */
            holidayfile,     /* 29 (I) holiday file name                    */
            BusVolDCC,       /* 30 (I) BUS/251F or BUS/BUS                  */
            results          /* 31 (O) Pricing results.                     */ 
            ) == FAILURE) goto RETURN;

        status = SUCCESS;
        return status;
    }


    /* calibrate PA measure for each index */
    for (i = 1; i >= 0; i--)
    {
        
        /* select index, note ordering is reversed */
        if (i == 1)
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

        /* Fwd dates and rates for index 1. */
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
            holidayfile,
            BusVolDCC,
            (rate[i]).a) == FAILURE) goto RETURN;
        
        if (sigATM * sqrt((rate[i].s.expiry)) < Q3_MIN_VOL_CALIB) corr = 0.;

        /* smile: initialize */
        if (Q3SmileInit(
                        (rate[i]).s.fwdRate,
                        sigATM,
                        0.,
                        (rate[i]).s.expiry,
                        (rate[i]).s.expiryVolvol,
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
            mq[i].sigATM,
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
    if (Q3MQBivarInit(
        optType,
        numPayoffParams,
        payoffParams,
        &corr,
        &(ppa[0]),
        &pf,
        &payFunc) == FAILURE) goto RETURN;
    /* integrate payoff over density. */
    if (Q3MQGridPricer(
        &pf,
        payFunc,
        ppa[0],
        &price) == FAILURE) goto RETURN;


    results[0] = price;
    results[1] = pa[1].fwdRate; /* index 1 */
    results[2] = pa[0].fwdRate; /* index 2 */

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: failed\n", routine);  
    }
  
    return status;

} /* Q3MQBivarPricer */


/*-----------------------------------------------------------------------------
 * Q3MQBivarPricerSimple
 * 
 * Multi-q bivariate option pricer with simplified interface
 *
 */
int Q3MQBivarPricerSimple(
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

    PAYOFF    pf;
    MQDATA    mq[2];
    MQDATA   *ppa[2];
    double    vol;
    double    fwd;
    double    price;
    long      type;
    int       i;

    /* calibrate PA measure for each index */
    for (i = 1; i >= 0; i--)
    {
        
        /* select index                         */
        /* note: ordering of rates is reversed! */
        if (i == 1)
        {
            type      = rateType1;
            vol       = sigATM1;
            smile     = smile1;
            fwd       = fwdRate1;
        }
        else
        {
            type      = rateType2;
            vol       = sigATM2;
            smile     = smile2;
            fwd       = fwdRate2;
        }

        /* different calibrations for rate and spread */
        switch (type)
        {
        case Q3_BV_RATE:

            if ( vol * sqrt(expiry) < Q3_MIN_VOL_CALIB) corr = 0.;

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
    if (Q3MQBivarInit(
        optType,
        numPayoffParams,
        payoffParams,
        &corr,
        &(ppa[0]),
        &pf,
        &payFunc) == FAILURE) goto RETURN;
    
    /* integrate payoff over density. */
    if (Q3MQGridPricer(
        &pf,
        payFunc,
        ppa[0],
        &price) == FAILURE) goto RETURN;

    results[0] = price;
    results[1] = mq[1].fwdRate; /* index 1 */
    results[2] = mq[0].fwdRate; /* index 2 */

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: failed\n", routine);  
    }
  
    return status;

} /* Q3MQBivarPricerSimple */


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
    )
{
    static char routine[] = "Q3MQPSAPricer";
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
    long     *indxDates     = NULL;
    char     *rateDCC       = NULL;

    PAYOFF    pf;
    FADATA    fa;       /* CMS Yield distribution in FA measure*/
    MQDATA    mq[2];
    MQDATA    pa;       /* re-parameterize CMS Yield distribution in FA */
    MQDATA   *ppa[2];
    double    sigATM;
    double    price;
    int       i;

    double    basisParRate, basisLegFv;
    double    payPrms[8];
    double    AnnFv, zeroDelay;
    double    FADensNorm, AAtoFAScaleAdj = 1.;


    /* Check optType validity */ 
    if ( optType != Q3_BS_PERC_CALL &&
         optType != Q3_BS_PERC_PUT &&
         optType != Q3_BS_PERC_RATE_CALL &&
         optType != Q3_BS_PERC_RATE_PUT &&
         optType != Q3_BS_JOINT_FWD)
    {
        Q3ErrMsg("%s failed: Invalid option type \n", routine);  
        goto RETURN;
    }


    /* Fwd dates and rates for CMS. */
    if (Q3SwapRateCalc2(
            today,
            today,
            valueDate,
            numLiborZeroPts,
            liborZeroDates,
            liborZeroRates,
            numDiscPts,
            discDates,
            discRates,
            bsFixLegFreq,
            bsFixLegDCC,
            expiryDate,
            startDate,
            matDate,
            payDate,
            holidayfile,
            BusVolDCC,
            (rate[0]).a) == FAILURE) goto RETURN;

    /* Basis par swap yield, basis leg PV*/
    if (Q3BasisLegCalc(
            valueDate,
            numBsZeroPts,
            bsZeroDates,
            bsZeroRates,
            numDiscPts,
            discDates,
            discRates,
            (long)rate[0].s.freq,
            bsFixLegDCC,
            bsFltLegFreq,
            bsFltLegDCC,
            expiryDate,
            startDate,
            matDate,
            payDate,
            holidayfile,
            (rate[1].a)) == FAILURE) goto RETURN;


    if (cmsSigATM * sqrt((rate[0].s.expiry)) < Q3_MIN_VOL_CALIB ||
        bsSpreadSigATM * sqrt(rate[0].s.expiry) < Q3_MIN_VOL_CALIB) corr = 0.;

    basisParRate = rate[1].a[0];
    basisLegFv   = rate[1].a[1];
    AnnFv        = rate[1].a[2];
    zeroDelay      = rate[1].a[3];

    /* Assign spread expiry, expiryVolvol, etc.*/
    rate[1].s.expiry = rate[0].s.expiry;
    rate[1].s.expiryVolvol = rate[0].s.expiryVolvol;

    /* spread fwd = basisSwapYld / CMSYld*/
    rate[1].s.fwdRate = basisParRate / rate[0].s.fwdRate;

    for (i = 0; i <= 1; i++)
    {
        
        /* select index, note ordering is reversed */
        if (i == 0)
        {
            sigATM          = cmsSigATM;
            smile           = cmsSmile;
        }
        else
        {
            sigATM          = bsSpreadSigATM;
            smile           = bsSpreadSmile;
        }

         /* smile: initialize */
         if (Q3SmileInit(
                         (rate[i]).s.fwdRate,
                         sigATM,
                         0.,
                         (rate[i]).s.expiry,
                         (rate[i]).s.expiryVolvol,
                         1.,
                         smile,
                         &(mq[i])) == FAILURE)
         {
             goto RETURN;
         }

        /* Calibrate MQ distribution to SV option prices. */
        if (Q3MQCalib(&(mq[i])) == FAILURE) goto RETURN;
    }

    /* Calibrate CMS yield distribution in FA measure    */
    fa.mq = &(mq[0]);
    if (Q3FASmileInit(
            (rate[0]).s.expiry,
            mq[0].sigATM,
            (rate[0]).s.start,
            (long) (rate[0]).s.freq,
            (rate[0]).s.swapMat,
            (rate[0]).s.fwdRate30360,
            (rate[0]).s.fwdAnn30360,
            (rate[0]).s.zeroRateSwap,
            (rate[0]).s.payDelay, 
            (rate[0]).s.zeroRatePay,
            numVNFMParams,
            vnfmParams,
            setlType,
            setlType,
            &fa) == FAILURE) goto RETURN;
     
    /* initialize PA MQ */
    if (Q3MQCopySmileFromMQ(
            &(mq[0]),
            &pa) == FAILURE) goto RETURN;
       
    /* set NCK parameters */
    if (Q3DecodeNCK(
            Q3_BV_NCK, 
            &pa) == FAILURE) goto RETURN;

    /* compute calibration targets */
    if (Q3MQTargetFA( 
            &fa,
            &pa) == FAILURE) goto RETURN;

    /* bootstrap PA multi-q measure */
    if (Q3MQBootstrapQ(&pa) == FAILURE) goto RETURN;

    /* calibrate PA */
    if (Q3MQCalib(&pa) == FAILURE) goto RETURN;

    if (Q3FADensNorm(&fa, &FADensNorm) == FAILURE) goto RETURN;
    
    AAtoFAScaleAdj = basisParRate / ( FADensNorm * basisLegFv );
   
    /* Adjust spread forward:
     * Compute E[spread * (1 - Z(T, T+M)/Z(T, Pay)] 
     */
    ppa[0] = &pa;
    ppa[1] = &(mq[1]);
    payPrms[0] = fa.alphaAnn;
    payPrms[1] = fa.powerAnn;
    payPrms[2] = fa.freqAnn;
    payPrms[3] = fa.matAnn;
    payPrms[4] = fa.alphaDel;
    payPrms[5] = fa.powerDel;
    payPrms[6] = fa.freqDel;
    payPrms[7] = fa.matDel;

    if (Q3MQBivarInit(
        Q3_BS_PERC_YLD_ANN,
        8,
        payPrms,
        &corr,
        &(ppa[0]),
        &pf,
        &payFunc) == FAILURE) goto RETURN;
    /* integrate payoff over density. */
    if (Q3MQGridPricer(
        &pf,
        payFunc,
        ppa[0],
        &price) == FAILURE) goto RETURN;
    
    /* Adjust spread forward :
     * The basis float leg FV is scaled first due to
     * numerical scaling from AA to FA.
     * See details in Changhong He's Doc.
     */
    (mq[1]).fwdRate *= AAtoFAScaleAdj * basisLegFv / price;

    /* Price under the payment measure */
    ppa[0] = &pa;
    ppa[1] = &(mq[1]);      
    if (Q3MQBivarInit(
        optType,
        numPayoffParams,
        payoffParams,
        &corr,
        &(ppa[0]),
        &pf,
        &payFunc) == FAILURE) goto RETURN;
    /* integrate payoff over density. */
    if (Q3MQGridPricer(
        &pf,
        payFunc,
        ppa[0],
        &price) == FAILURE) goto RETURN;

    results[0] = price;
    results[1] = mq[1].fwdRate;     /* spread forward */
    results[2] = pa.fwdRate;        /* CMS forward    */

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: failed\n", routine);  
    }
  
    return status;
}/* Q3MQPSAPricer */


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
    )
{
    static char routine[] = "Q3SmileBivarSimpleCorrCalib";
    int         status    = FAILURE;
    
    double      smile1_sav[13];
    double      smile2_sav[13];
    double      payoffParams_sav[2];

    double      lbound, ubound;
    double      corr, calibErr;
    double      bvprice, tol, delta, grad, step, stepSize;
    double      alpha;  /*restrict step size in line search */
    int         iter, i;
    int         mfound = FALSE;
    double      bvResults[3];
    double      MTMLbound, MTMUbound;   
    double      leverage = 1.;
    double      tmp;

    /* Save smiles */
    for (i = 0; i < 13; i++)
    {
        smile1_sav[i] = smile1[i];
        smile2_sav[i] = smile2[i];
    }

    /* Save payoff parameters*/
    payoffParams_sav[0] = payoffParams[0];
    payoffParams_sav[1] = payoffParams[1];

    /* Check payoff type */
    if ( optType == Q3_CALL)
        optType = Q3_CALL_SUM;
    else if (optType == Q3_PUT)
        optType = Q3_PUT_SUM;
    else
    {
        Q3ErrMsg("%s failed: option type %ld must be %ld or %ld \n", 
                 routine, optType, Q3_CALL, Q3_PUT);  
        goto RETURN;

    }

    /* Check payoff Params*/
    if ( fabs(payoffParams[0]) < TINY ||
         fabs(payoffParams[1]) < TINY)
    {
        Q3ErrMsg("%s failed: payoff parameters cannot be zero.\n", 
                 routine);  
        goto RETURN;

    }

    /* If both payoff parameters are negative, change payoff type.*/
    if ( payoffParams[0] < -TINY &&
         payoffParams[1] < -TINY)
    {
        optType = (optType == Q3_CALL_SUM)? Q3_PUT_SUM: Q3_CALL_SUM;
        payoffParams[0] *= -1.;
        payoffParams[1] *= -1.;
        strike          *= -1;
    }

    /* Setup payoff parameters*/
    if ( payoffParams[1] > TINY) 
    {
        /* The weight of rate 2 > 0*/
        leverage = payoffParams[1];
        payoffParams[1] = payoffParams[0] / leverage;
        payoffParams[0] = strike / leverage;
    }
    else
    {
        /* The weight of rate 2 < 0 and
         * the weight of rate 1 > 0.
         * Switch rate1 and rate2
         */
        tmp      = fwdRate1;
        fwdRate1 = fwdRate2;
        fwdRate2 = tmp;

        tmp      = sigATM1;
        sigATM1  = sigATM2;
        sigATM2  = tmp;

        for (i = 0; i < 13; i++)
        {
            smile1[i] = smile2_sav[i];
            smile2[i] = smile1_sav[i];
        }

        leverage = payoffParams[0];
        payoffParams[0] = strike / leverage;
        payoffParams[1] = payoffParams[1] / leverage;
    }
    
    
    /* Set MTM bounds */
    if (Q3MQBivarPricerSimple(
            expiry, 
            0,          /*rate type*/
            fwdRate1,
            smile1,  
            sigATM1,    
            0,          /*rate type*/
            fwdRate2,
            smile2,  
            sigATM2,   
            Q3_MAX_CORR,          
            optType,       
            2,
            payoffParams,
            bvResults) == FAILURE) goto RETURN;  
    MTMLbound = leverage * bvResults[0];

    if (Q3MQBivarPricerSimple(
            expiry, 
            0,
            fwdRate1,
            smile1,  
            sigATM1,    
            0,
            fwdRate2,
            smile2,  
            sigATM2,   
            Q3_MIN_CORR,          
            optType,       
            2,
            payoffParams,
            bvResults) == FAILURE) goto RETURN;  
    MTMUbound = leverage * bvResults[0];

    if ( (MTMLbound - MTMUbound) > TINY )
    {
        tmp = MTMLbound;
        MTMLbound = MTMUbound;
        MTMUbound = tmp;
    }

    if ( (premium - MTMLbound) < -TINY ||
         (premium - MTMUbound) > TINY )
    {
        Q3ErrMsg("%s failed: premium %f is out of the bound [%f, %f]\n", 
                routine,
                premium,
                MTMLbound,
                MTMUbound);  
        goto RETURN;
    }

    lbound   = Q3_MIN_CORR;
    ubound   = Q3_MAX_CORR;
    corr     = (fabs(initialGuess) > Q3_MAX_CORR) ? 0. : initialGuess;      
    calibErr = 0.;
    iter     = 0;
    alpha    = 0.5;
    tol      = premium * CORR_CALIB_TOL;

    /* Newton method  
     * Objective function : BvSimple(corr) - premium
     * Constraint:  lbound <= corr <= ubound
     */
    do 
    {
        /* Evaluate objective value */
        if (Q3MQBivarPricerSimple(
            expiry, 
            0,
            fwdRate1,
            smile1,  
            sigATM1,    
            0,
            fwdRate2,
            smile2,  
            sigATM2,   
            corr,          
            optType,       
            2,
            payoffParams,
            bvResults) == FAILURE) goto RETURN;  

        bvprice  = leverage * bvResults[0];
        calibErr = bvprice - premium;
        
        if ( fabs(calibErr) < tol )
        {
            mfound = TRUE;
            break;
        }

        /* Compute gradient by finite difference */
        delta = (corr + CORR_DELTA > ubound)? -CORR_DELTA : CORR_DELTA;

        if (Q3MQBivarPricerSimple(
            expiry, 
            0,
            fwdRate1,
            smile1,  
            sigATM1,    
            0,
            fwdRate2,
            smile2,  
            sigATM2,   
            corr+delta,          
            optType,       
            2,
            payoffParams,
            bvResults) == FAILURE) goto RETURN;  

        grad = (leverage * bvResults[0] - bvprice) / delta;
 
        if ( fabs(grad) < CORR_CALIB_RESN ) 
        {
            corr += delta;
        }
        else
        {
            /* Newton-Raphson step */
            step     = - calibErr / grad;

            /* Restrict step size */
            stepSize = MIN( fabs((corr - lbound) / step), 
                            fabs((ubound - corr) / step));
            stepSize = MIN(alpha * stepSize, 1.0);
            corr    += stepSize * step;
        }

        iter++;
    } while( fabs(calibErr) > tol && 
             iter < CORR_CALIB_MAX_ITER);

    if (mfound == FALSE)
    {
        if (Q3MQBivarPricerSimple(
            expiry, 
            0,
            fwdRate1,
            smile1,  
            sigATM1,    
            0,
            fwdRate2,
            smile2,  
            sigATM2,   
            corr,          
            optType,       
            2,
            payoffParams,
            bvResults) == FAILURE) goto RETURN;  

        calibErr = leverage * bvResults[1] - premium;
        mfound = (fabs(calibErr) < tol) ? TRUE : FALSE;
    }

    results[0] = corr;
    results[1] = calibErr;
    results[2] = (mfound == TRUE) ? 1.0 : 0.0;

    /* Recover smiles */
    for (i = 0; i < 13; i++)
    {
        smile1[i] = smile1_sav[i];
        smile2[i] = smile2_sav[i];
    }

    payoffParams[0] = payoffParams_sav[0];
    payoffParams[1] = payoffParams_sav[1];

    status = SUCCESS;

 RETURN:

    /* Recover smiles */
    for (i = 0; i < 13; i++)
    {
        smile1[i] = smile1_sav[i];
        smile2[i] = smile2_sav[i];
    }

    payoffParams[0] = payoffParams_sav[0];
    payoffParams[1] = payoffParams_sav[1];

    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: failed\n", routine);  
    }
  
    return status;

}/*Q3SmileBivarSimpleCorrCalib*/



/*-----------------------------------------------------------------------------
 * Q3MQBivarCorr
 * 
 * Multi-q bivariate correlation calculator.
 *
 */
int Q3MQBivarCorr(
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
    double *corr,            /* (I) Index correlation data               */
    double *fwdVolInput,     /* (I) Forward vol betw rate start dates    */
    long    payDate,         /* (I) Option payment date                  */
    long    setlType,        /* (I) Cash or physical settlement          */
    char   *holidayfile,     /* (I) holiday file name                    */
    char   *BusVolDCC,       /* (I) BUS/251F or BUS/BUS                  */
    double *results          /* (O) Pricing results.                     */
    )
{
    static char routine[] = "Q3MQBivarCorr";
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
    double    corrTmp, corrRate;
    double    e1, esq1, e2, esq2, e12, var1, var2, denom;
    long      numIndxPts;
    long      rateFreq;
    long      startDate, matDate, expiryDate;
    int       i;

    if (corr == NULL)
    {
        Q3ErrMsg("%s: Correlation input required.\n", routine);
        goto RETURN;
    }

    /* calibrate PA measure for each index */
    for (i = 1; i >= 0; i--)
    {
        
        /* select index, note ordering is reversed */
        if (i == 1)
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
        if (sigATM * (rate[i].s.expiry) < TINY) *corr = 0.;

        /* Fwd dates and rates for index 1. */
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
            holidayfile,
            BusVolDCC,
            (rate[i]).a) == FAILURE) goto RETURN;

        /* smile: initialize */
        if (Q3SmileInit(
                        (rate[i]).s.fwdRate,
                        sigATM,
                        0.,
                        (rate[i]).s.expiry,
                        (rate[i]).s.expiryVolvol,
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
            mq[i].sigATM,
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
        corr[0] *= sqrt(1. - adj);
    }
    
    /* compute joint expectation */
    ppa[0] = &(pa[0]);
    ppa[1] = &(pa[1]);
    if (Q3MQBivarInit(
        Q3_JOINT_FWD,
        0,    /* number of payment parameters */
        NULL, /* payment parameters           */
        corr,
        &(ppa[0]),
        &pf,
        &payFunc) == FAILURE) goto RETURN;

    if (Q3MQGridPricer(
        &pf,
        payFunc,
        ppa[0],
        &e12) == FAILURE) goto RETURN;

    /* compute expectation of squares */
    /* set correlation to 0.          */
    corrTmp = 1.;
    ppa[0]  = &(pa[0]);
    ppa[1]  = &(pa[0]);
    if (Q3MQBivarInit(
        Q3_JOINT_FWD,
        0,
        NULL,
        &corrTmp,
        &(ppa[0]),
        &pf,
        &payFunc) == FAILURE) goto RETURN;

    if (Q3MQGridPricer(
        &pf,
        payFunc,
        ppa[0],
        &esq2) == FAILURE) goto RETURN;
    
    ppa[0] = &(pa[1]);
    ppa[1] = &(pa[1]);
    if (Q3MQBivarInit(
        Q3_JOINT_FWD,
        0,
        NULL,
        &corrTmp,
        &(ppa[0]),
        &pf,
        &payFunc) == FAILURE) goto RETURN;

    if (Q3MQGridPricer(
        &pf,
        payFunc,
        ppa[0],
        &esq1) == FAILURE) goto RETURN;

    /* independent expectation from PA calibration */
    e1 = pa[1].fwdRate;
    e2 = pa[0].fwdRate;

    /* compute rate correlation */
    var1 = esq1 - e1 * e1;
    var2 = esq2 - e2 * e2;
    denom = sqrt(var1 * var2);
    if (fabs(denom) < TINY)
    {
        Q3ErrMsg("%s: Variance too small.\n", routine);
        goto RETURN;
    }

    corrRate = (e12 - e1 * e2) / denom;

    results[0] = corrRate;

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: failed\n", routine);  
    }
  
    return status;

} /* Q3MQBivarCorr */


/*-----------------------------------------------------------------------------
 * Q3MQBivarInit
 *
 * Select payoff function and populate payoff parameters. 
 * Routine may reorder measure data elements.
 *
 */
int Q3MQBivarInit(
    long      optType,    /* (I)   option type             */
    long      numPayPrms, /* (I)   number of payoff params */
    double   *payPrms,    /* (I)   payoff coefficients     */
    double   *corr,       /* (I)   index correlation       */
    MQDATA  **mq,         /* (I/O) measure data            */
    PAYOFF   *pf,         /* (O)   option payoff structure */
    FPAYOFF **payFunc     /* (O)   payoff function         */
    )
{
    static char routine[]   = "Q3MQBivarInit";
    int         status      = FAILURE;
    long        calJoint    = FALSE;
    long        calIndex = 0;
    MQDATA     *mqTmp;
    long        i;
    double      RIB1, RIB2;
        
    /* number of payoff params */
    if (numPayPrms > Q3_MAX_PAY_PARAMS)
    {
        Q3ErrMsg("%s: Number of payoff params exceeds maximum.\n", routine);
        goto RETURN;
    }

    switch (optType)
    {   

    case Q3_ADJ_FWD:

        *payFunc = Q3Pay1D_YldNull;
        break;

    case Q3_JOINT_FWD:
    case Q3_BS_JOINT_FWD:
        
        *payFunc = Q3Pay1D_YldYld;
        break;

    case Q3_VNL + Q3_CALL:

        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_VnlNull;
        pf->cop  = 1;
        break;

    case Q3_VNL + Q3_PUT:

        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_VnlNull;
        pf->cop  = -1;
        break;

    case Q3_CALL_SUM:           

        if (numPayPrms < 2)
        {
            Q3ErrMsg("%s: Payoff requires 2 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_Sum;
        pf->cop  = 1;
        break;

    case Q3_PUT_SUM:           

        if (numPayPrms < 2)
        {
            Q3ErrMsg("%s: Payoff requires 2 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_Sum;
        pf->cop  = -1;
        break;

    case Q3_CALL_PROD:           

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_Prod;
        pf->cop  = 1;
        break;

    case Q3_PUT_PROD:           

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_Prod;
        pf->cop  = -1;
        break;
    
    case Q3_CALL_PERC:           

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params.\n", routine);
            goto RETURN;
        }
        pf->params[0] = payPrms[0];
        pf->params[1] = 0.;  /* zero spread */
        *payFunc = Q3Pay1D_Perc;
        pf->cop  = 1;
        break;
    
    case Q3_PUT_PERC:           

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params.\n", routine);
            goto RETURN;
        }
        pf->params[0] = payPrms[0];
        pf->params[1] = 0.; /* zero spread */
        *payFunc = Q3Pay1D_Perc;
        pf->cop  = -1;
        break;

    case Q3_CALL_PERC_WGT:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_PercWgt;
        pf->cop  = 1;
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;
    
    case Q3_PUT_PERC_WGT:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_PercWgt;
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        pf->cop  = -1;
        break;



    case Q3_FLR_W_FLR:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_FlrFlrOrCapCap;
        pf->cop  = 1;
        break;

    case Q3_CAP_W_CAP:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_FlrFlrOrCapCap;
        pf->cop  = -1;
        break;

    case Q3_FLR_W_CAP:           
        
        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_FlrCapOrCapFlr;
        pf->cop  = 1;
        break;
    
    case Q3_CAP_W_FLR:           
        
        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_FlrCapOrCapFlr;
        pf->cop  = -1;
        break;    
    
    case Q3_FLR_W_FLR_EMBED:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_FlrFlrOrCapCapEmbedFlt;
        pf->cop  = 1;
        break;
    
    case Q3_CAP_W_CAP_EMBED:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_FlrFlrOrCapCapEmbedFlt;
        pf->cop  = -1;
        break;

    case Q3_FLR_W_CAP_EMBED:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_FlrCapOrCapFlrEmbedFlt;
        pf->cop  = 1;
        break;
    
    case Q3_CAP_W_FLR_EMBED:           

        if (numPayPrms < 3)
        {
            Q3ErrMsg("%s: Payoff requires 3 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_FlrCapOrCapFlrEmbedFlt;
        pf->cop  = -1;
        break;

    case Q3_FLR_CAP_SUM:
        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_FlrCapSum;
        pf->cop  = 1;
        break;

    case Q3_BS_PERC_RATE_CALL:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_Prod;
        pf->cop  = 1;
        calJoint = TRUE;
        calIndex = 1;
        break;

    case Q3_BS_PERC_RATE_PUT:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 params.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_Prod;
        pf->cop  = -1;
        calJoint = TRUE;
        calIndex = 1;
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
            /* Backward compatible */
            pf->params[1] = 0.;
        }
        *payFunc = Q3Pay1D_Perc;
        pf->cop  = 1;
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        calIndex = 0;
        calJoint = TRUE;
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
            /* Backward compatible */
            pf->params[1] = 0.;
        }
        *payFunc = Q3Pay1D_Perc;
        pf->cop  = -1;
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        calIndex = 0;
        calJoint = TRUE;
        break;

    case Q3_BS_SPRD_RATE_CALL:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->params[1] = 1.;         /* payFunc expects leverage */
        *payFunc      = Q3Pay1D_Sum;
        pf->cop       = 1;
        break;
        
    case Q3_BS_SPRD_RATE_PUT:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        pf->params[1] = 1.;         /* payFunc expects leverage */
        *payFunc      = Q3Pay1D_Sum;
        pf->cop       = -1;
        break;

    case Q3_BS_PERC_SPRD_CALL:

        if (numPayPrms < 1)
        {
            Q3ErrMsg("%s: Payoff requires 1 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_VnlNull; /* Vanilla on 1st variable */
        pf->cop  = 1.;
        /* payoff function expects order: rate1, rate2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
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
        *payFunc = Q3Pay1D_VnlNull;
        /* payoff function expects order: rate1, rate2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
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
        *payFunc      = Q3Pay1D_DigSpdInRIB_EPS;
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
        *payFunc      = Q3Pay1D_DigSpdOutRIB_EPS;
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
        *payFunc      = Q3Pay1D_DigSpdInRIB_EPS;

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
        *payFunc      = Q3Pay1D_DigSpdOutRIB_EPS;
        
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
        *payFunc      = Q3Pay1D_DigPercInRIB_EPS;
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
        *payFunc      = Q3Pay1D_DigPercInRIB_EPS;
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
        *payFunc      = Q3Pay1D_DigPercOutRIB_EPS;
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
        *payFunc = Q3Pay1D_InRIB;
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_IN_BIRIB_EPS:
        
        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        if(RIB1 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB1))
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "leps (or heps) < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        /* Ignore small epsilons */
        if ((fabs(payPrms[6]) < Q3_MIN_EPS) &&
            (fabs(payPrms[7]) < Q3_MIN_EPS))
        {
            *payFunc = Q3Pay1D_InRIB;
        }
        else
        {
            *payFunc = Q3Pay1D_InRIB_EPS;
        }
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;
            
    case Q3_OUT_BIRIB:

        if (numPayPrms < 6)
        {
            Q3ErrMsg("%s: Payoff requires 6 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_OutRIB;
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_OUT_BIRIB_EPS:
        
        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        if(RIB1 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB1))
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "leps (or heps) < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
        /* Ignore small epsilons */
        if ((fabs(payPrms[6]) < Q3_MIN_EPS) &&
            (fabs(payPrms[7]) < Q3_MIN_EPS))
        {
            *payFunc = Q3Pay1D_OutRIB;
        }
        else
        {
            *payFunc = Q3Pay1D_OutRIB_EPS;
        }
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_IN_SPDRIB:

        if (numPayPrms < 6)
        {
            Q3ErrMsg("%s: Payoff requires 6 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_InSpdRIB;
        /* swap measures to condition on 1st index in payoff*/
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
        break;

    case Q3_IN_SPDRIB_EPS:

        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        if(RIB1 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB1))
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "leps (or heps) < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
    /* Ignore small epsilons */
    if ((fabs(payPrms[6]) < Q3_MIN_EPS) &&
        (fabs(payPrms[7]) < Q3_MIN_EPS))
    {
        *payFunc = Q3Pay1D_InSpdRIB;
    }
    else
    {
        *payFunc = Q3Pay1D_InSpdRIB_EPS;
    }
    /* swap measures to condition on 1st index in payoff*/
    mqTmp = mq[0];
    mq[0] = mq[1];
    mq[1] = mqTmp;  
        break;

    case Q3_OUT_SPDRIB:

        if (numPayPrms < 6)
        {
            Q3ErrMsg("%s: Payoff requires 6 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_OutSpdRIB;
        /* swap measures to condition on 1st index in payoff*/
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_OUT_SPDRIB_EPS:

        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 8 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        RIB1 = pf->params[1] - pf->params[0];
        if(RIB1 < 0)
        {
            Q3ErrMsg("%s: Requires HB > LB.\n", routine);
            goto RETURN;
        }
        RIB1 /= 2.0;
        if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
            (pf->params[7] < 0 && -pf->params[7] > RIB1))
        {
            Q3ErrMsg("%s: When leps < 0 or heps < 0,",
                "leps (or heps) < 0.5 *(HB - LB). \n", routine);
            goto RETURN;
        }
    /* Ignore small epsilons */
    if ((fabs(payPrms[6]) < Q3_MIN_EPS) &&
        (fabs(payPrms[7]) < Q3_MIN_EPS))
    {
        *payFunc = Q3Pay1D_OutSpdRIB;
    }
    else
    {
        *payFunc = Q3Pay1D_OutSpdRIB_EPS;
    }
        /* swap measures to condition on 1st index in payoff*/
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;  
        break;

     case Q3_IN_AND_IN:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay1D_InAndIn;
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
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
        *payFunc = Q3Pay1D_InAndIn;

        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_IN_OR_IN:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay1D_InOrIn;
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
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
        *payFunc = Q3Pay1D_InOrIn;
        
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_IN_AND_OUT:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay1D_InAndOut;
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
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
        *payFunc = Q3Pay1D_InAndOut;
        
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_IN_OR_OUT:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay1D_InOrOut;
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
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
        *payFunc = Q3Pay1D_InOrOut;
        
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_OUT_AND_IN:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay1D_OutAndIn;
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
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
        *payFunc = Q3Pay1D_OutAndIn;

        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_OUT_OR_IN:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay1D_OutOrIn;
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
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
        *payFunc = Q3Pay1D_OutOrIn;
        
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_OUT_AND_OUT:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay1D_OutAndOut;
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
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
        *payFunc = Q3Pay1D_OutAndOut;
        
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_OUT_OR_OUT:

        if (numPayPrms < 4)
        {
            Q3ErrMsg("%s: Payoff requires 4 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        for (i=4; i<8; i++) pf->params[i] = 0.0;
        *payFunc = Q3Pay1D_OutOrOut;
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
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
        *payFunc = Q3Pay1D_OutOrOut;
        
        /* payoff function expects order: rate1, rate 2 */
        mqTmp = mq[0];
        mq[0] = mq[1];
        mq[1] = mqTmp;
        break;

    case Q3_COMPLEX_SPD:

        if (numPayPrms < 19)
        {
            Q3ErrMsg("%s: Payoff requires 19 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_ComplexSPD;
        break;

    case Q3_BS_PERC_YLD_ANN:
        
        if (numPayPrms < 8)
        {
            Q3ErrMsg("%s: Payoff requires 1 param.\n", routine);
            goto RETURN;
        }
        for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
        *payFunc = Q3Pay1D_SpdYldAnn;
        break;


    default: 
        
        Q3ErrMsg("%s: Unrecognized option type.\n", routine);
        goto RETURN;

    }
    
    /* joint calibration */
    if (calJoint == TRUE)
    { 
        /* if no user calibration parameter, fwd cal. only */
        if (numPayPrms > 2 && 
            fabs(payPrms[2] - Q3_BS_CAL_VOL) < TINY)
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

    /* fill in common payoff parameters */
    pf->optType = optType;
    pf->corr    = corr[0];
    (pf->mq)[0] = mq[0];
    (pf->mq)[1] = mq[1];

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: failed\n", routine);
    }
 
    return status;

} /* Q3MQBivarInit */


/*-----------------------------------------------------------------------------
 * Q3MQCalibJointFwd
 *
 * Calibrate joint forward to target by adjusting forward of selected index.
 *
 */
int Q3MQCalibJointFwd(
    double   target, /* (I)   calibration target     */
    double  *corr,   /* (I)   correlation of indices */
    long     index,  /* (I)   index to adjust (0,1)  */
    MQDATA **mq      /* (I/O) index measure data     */
    )
{
    const char routine[] = "Q3MQCalibJointFwd";
    int status = FAILURE;

    PAYOFF   pf;
    FPAYOFF *payFunc;    
    double   jointFwd;

    /* calculate joint forward */
    if (Q3MQBivarInit(
        Q3_JOINT_FWD,
        0,              /* num pay params */
        NULL,       
        corr,
        mq,
        &pf,
        &payFunc) == FAILURE) goto RETURN;

    if (Q3MQGridPricer(
        &pf,
        payFunc,
        mq[0],
        &jointFwd) == FAILURE) goto RETURN;

    /* fail on zero forward */
    if (jointFwd < TINY) goto RETURN;

    /* adjust selected forward */
    (mq[index])->fwdRate *= target / jointFwd;

    status = SUCCESS;

 RETURN: 

    if (status == FAILURE) Q3ErrMsg("%s: failed\n", routine);
 
    return status;

} /* Q3MQCalibJointFwd */


/*-----------------------------------------------------------------------------
 * Q3MQCalibJoint
 *
 *
 */
#define Q3_S_DELTA1 10*Q3_S_DELTA

int Q3MQCalibJoint(
    double  *corr,
    long     index,
    MQDATA **mq
    )
{
    const char routine[] = "Q3MQCalibJoint";
    int status = FAILURE;

    double targetVolPrice, targetFwd, compositeVol;
    double sigATM0, sigATM1, expiry;
    double volPrice, volPriceTwk, volPriceDelta, diff;
    long   i;

    PAYOFF   pf;
    FPAYOFF *payFunc;    

    /* local variables */
    sigATM0 = (mq[0])->sigATM; sigATM1 = (mq[1])->sigATM;
    expiry  = (mq[0])->expiry;

    /* compute forward target */
    targetFwd = (mq[0])->fwdRate * (mq[1])->fwdRate;

    /* composite volatility */
    compositeVol = 
        sigATM0*sigATM0 + sigATM1*sigATM1 + 2.*corr[0]*sigATM0*sigATM1;
    compositeVol = sqrt(compositeVol);

    /* compute volatility target */
    if (BSQPricer(
        targetFwd,    /* fwd       */
        targetFwd,    /* strike    */
        expiry,
        compositeVol,
        1,            /* lognormal */
        Q3_CALL,
        &targetVolPrice) == FAILURE) goto RETURN;    

    /* initialization for bivariate pricing */
    if (Q3MQBivarInit(
        Q3_CALL_PROD,
        1,
        &targetFwd,
        corr,
        mq,
        &pf,
        &payFunc) == FAILURE) goto RETURN;

    /* Newton-Raphson to match targets */
    for (i = 0; i < Q3_S_STEPS; i++)
    {
        /* calibrate forward */
        if (Q3MQCalibJointFwd(
            targetFwd,
            corr,
            index,
            mq) == FAILURE) goto RETURN;

        /* price vol */
        if (Q3MQGridPricer(
            &pf,
            payFunc,
            mq[0],
            &volPrice) == FAILURE) goto RETURN;
    
        diff = targetVolPrice - volPrice;

        if (fabs(diff) > .01 * targetVolPrice)
        {
            /* N-R */
            (mq[index])->sigMQ += Q3_S_DELTA1;

            if (Q3MQCalibJointFwd(
                targetFwd,
                corr,
                index,
                mq) == FAILURE) goto RETURN;

            if (Q3MQGridPricer(
                &pf,
                payFunc,
                mq[0],
                &volPriceTwk) == FAILURE) goto RETURN;
            
            (mq[index])->sigMQ -= Q3_S_DELTA1;
            volPriceDelta = volPrice - volPriceTwk;
            if (fabs(volPriceDelta) < Q3_MQ_RESN * Q3_S_DELTA1)
            {
                Q3ErrMsg("%s: df/ds = 0 in Newton-Raphson.\n", routine);
                goto RETURN;
            }
            
            (mq[index])->sigMQ -= diff*Q3_S_DELTA1/volPriceDelta;
        }
        else
        {
            if (Q3MQCalibJointFwd(
                targetFwd,
                corr,
                index,
                mq) == FAILURE) goto RETURN;
            status = SUCCESS;
            goto RETURN;
        }
    } /* N-R loop */
    

 RETURN:

    if (status == FAILURE) Q3ErrMsg("%s: failed\n", routine);

    return status;

} /* Q3MQCalibJoint */


/*f----------------------------------------------------------------------------
 * Q3BivarVersionL
 *
 * Add-in version number
 */

int Q3BivarVersionL(char* versionL) 
{
    versionL[0] = 1; 
    /* U_VERSION is expected to be defined when compiling - see bv_q_defs.mk */
    strcpy(versionL+1, Q3_BV_VERS);  
    return SUCCESS; 
}

