/******************************************************************************
 * Module:      
 * Submodule:
 * File:      marginalDist.c
 * Function:    
 * Author:    Interest Rates DR
 * Revision:  $Header$
 *****************************************************************************/
#include <math.h>
#include <ctype.h>                      
#include <stdio.h>

#include "aliberr.h"
#include "aliblog.h"
#include "mqbvl.h"
#include "mqbv.h"
#include "marginalDist.h"

#define Q3_BV_NCK  2000000769 /* Calibration control data */

#define Q3_BV_RATE   0
#define Q3_BV_SPREAD 1

/*-----------------------------------------------------------------------------
 * Q3MQFAMarginalDistL
 *
 * Multi-q bivariate option pricer add-in.
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
    double *payoffParamsL,    /* 25 (I) payment parameters                */
    long   *traceL,           /* 26 (I) trace                             */
    TMatrix2D **xloutputsL    /* 27 (O) Marginal distributions            */
    )
{
    static char routine[] = "Q3MQFAMarginalDistL";

    int status = FAILURE;


    long expiryDate1L[2] = {1, 0};
    long expiryDate2L[2] = {1, 0};
    long startDate1L[2]  = {1, 0};
    long startDate2L[2]  = {1, 0};
    long matDate1L[2]    = {1, 0};
    long matDate2L[2]    = {1, 0};


    double smile1Input[13], smile2Input[13];

    int     i;
    double  *outputsL;

    double  smileSmile1L[20], smileSmile2L[20];

    memcpy(smileSmile1L+2, smile1L+1, ((long) smile1L[0])*sizeof(double));
    smileSmile1L[0] = smile1L[0] + 1.;
    smileSmile1L[1] = 10.;

    memcpy(smileSmile2L+2, smile2L+1, ((long) smile2L[0])*sizeof(double));
    smileSmile2L[0] = smile2L[0] + 1.;
    smileSmile2L[1] = 10.;

    outputsL   = (double *) Q3_DR_Array (DOUBLE, 0, Q3_1D_INT_PTS * 17);

    *xloutputsL = GtoMatrixNewEmpty( Q3_1D_INT_PTS, 17);
    if (*xloutputsL == NULL) goto RETURN;

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
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile1L,     "smile", 10);
    GTO_WRAP_CHECK_SCALAR(rateFreq2L,           "rate frequency");
    GTO_WRAP_CHECK_SCALAR(rateDCC2L,            "rate day count convention");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(rateRstDts2L,"rate reset dates", 3);
    GTO_WRAP_CHECK_SCALAR(sigATM2L,             "ATM volatility");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(smile2L,     "smile", 10);
    GTO_WRAP_CHECK_SCALAR(optTypeL,             "option type");
    GTO_WRAP_CHECK_SCALAR(payDateL,             "payment date");
    GTO_WRAP_CHECK_VECTOR_BIG_ENUF(corrL,       "corr", 1);
    GTO_WRAP_CHECK_SCALAR(setlTypeL,            "settlement type");
    
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

    memcpy(smile1Input, smileSmile1L+1, ((long) smileSmile1L[0])*sizeof(double));

    for (i = (long) smileSmile1L[0]; i<13; i++)
        smile1Input[i] = 0.;

    memcpy(smile2Input, smileSmile2L+1, ((long) smileSmile2L[0])*sizeof(double));

    for (i = (long) smileSmile2L[0]; i<13; i++)
        smile2Input[i] = 0.;


    if (Q3MQFAMarginalDist(
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
        outputsL) == FAILURE) goto RETURN;  

    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        (*xloutputsL)->data[i][0]  = outputsL[1+17*i];
        (*xloutputsL)->data[i][1]  = outputsL[2+17*i];
        (*xloutputsL)->data[i][2]  = outputsL[3+17*i];
        (*xloutputsL)->data[i][3]  = outputsL[4+17*i];
        (*xloutputsL)->data[i][4]  = outputsL[5+17*i];
        (*xloutputsL)->data[i][5]  = outputsL[6+17*i];
        (*xloutputsL)->data[i][6]  = outputsL[7+17*i];
        (*xloutputsL)->data[i][7]  = outputsL[8+17*i];
        (*xloutputsL)->data[i][8]  = outputsL[9+17*i];
        (*xloutputsL)->data[i][9]  = outputsL[10+17*i];
        (*xloutputsL)->data[i][10] = outputsL[11+17*i];
        (*xloutputsL)->data[i][11] = outputsL[12+17*i];
        (*xloutputsL)->data[i][12] = outputsL[13+17*i];
        (*xloutputsL)->data[i][13] = outputsL[14+17*i];
        (*xloutputsL)->data[i][14] = outputsL[15+17*i];
        (*xloutputsL)->data[i][15] = outputsL[16+17*i];
    }

    status = SUCCESS;

 done: /* for alib */
 RETURN:

    if (outputsL != NULL)
    {
        Q3_Free_DR_Array(outputsL, DOUBLE, 0, Q3_1D_INT_PTS * 9);
    }

    if (status == FAILURE)
    {
        Q3ErrMsg("%s: failed.\n", routine);
    }
  
    return status;

} /* Q3MQFAMarginalDistL */

/*-----------------------------------------------------------------------------
 * Q3MQMarginalDistL
 * 
 * Multi-q marginal distribution with simplified interface
 *
 */
int Q3MQMarginalDistL(
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
    TMatrix2D **xloutputsL    /* 14 (O) Marginal distributions           */
    )
{
    static char routine[] = "Q3MQMarginalDistL";

    int status = FAILURE;

    double smile1Input[13], smile2Input[13];

    int     i;
    double  *outputsL;

    double  smileSmile1L[20], smileSmile2L[20];

    memcpy(smileSmile1L+2, smile1L+1, ((long) smile1L[0])*sizeof(double));
    smileSmile1L[0] = smile1L[0] + 1.;
    smileSmile1L[1] = 10.;

    memcpy(smileSmile2L+2, smile2L+1, ((long) smile2L[0])*sizeof(double));
    smileSmile2L[0] = smile2L[0] + 1.;
    smileSmile2L[1] = 10.;

    outputsL   = (double *) Q3_DR_Array (DOUBLE, 0, Q3_1D_INT_PTS * 9);

    *xloutputsL = GtoMatrixNewEmpty( Q3_1D_INT_PTS, 9);
    if (*xloutputsL == NULL) goto RETURN;
    
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

    memcpy(smile1Input, smileSmile1L+1, ((long) smileSmile1L[0])*sizeof(double));

    for (i = (long) smileSmile1L[0]; i<13; i++)
        smile1Input[i] = 0.;

    memcpy(smile2Input, smileSmile2L+1, ((long) smileSmile2L[0])*sizeof(double));

    for (i = (long) smileSmile2L[0]; i<13; i++)
        smile2Input[i] = 0.;

    if (Q3MQMarginalDist(
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
        outputsL) == FAILURE) goto RETURN;  

    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        (*xloutputsL)->data[i][0] = outputsL[1+9*i];
        (*xloutputsL)->data[i][1] = outputsL[2+9*i];
        (*xloutputsL)->data[i][2] = outputsL[3+9*i];
        (*xloutputsL)->data[i][3] = outputsL[4+9*i];
        (*xloutputsL)->data[i][4] = outputsL[5+9*i];
        (*xloutputsL)->data[i][5] = outputsL[6+9*i];
        (*xloutputsL)->data[i][6] = outputsL[7+9*i];
        (*xloutputsL)->data[i][7] = outputsL[8+9*i];
    }
    status = SUCCESS;

 done: /* for alib */
 RETURN:
    if (outputsL != NULL)
    {
        Q3_Free_DR_Array(outputsL, DOUBLE, 0, Q3_1D_INT_PTS * 9);
    }

    if (status == FAILURE)
    {
        Q3ErrMsg("%s: failed.\n", routine);
    }
  
    return status;

} /* Q3MQMarginalDistL */

/*-----------------------------------------------------------------------------
 * Q3MQFAMarginalDist
 * 
 * Multi-q bivariate option pricer.
 *
 */
int Q3MQFAMarginalDist(
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
    static char routine[] = "Q3MQFAMarginalDist";
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

    SVDATA    sv[2];
    FADATA    fa[2];
    MQDATA    mq[2];
    MQDATA    pa[2];
    double    sigATM;
    long      numIndxPts, smileType;
    long      rateFreq;
    long      startDate, matDate, expiryDate;
    int       i;

    double grid[Q3_1D_INT_PTS];
    double yield[Q3_1D_INT_PTS];
    double dens[Q3_1D_INT_PTS];
    double densNorm;
    FILE  *fstream;

    /* calibrate PA measure for each index */
    for (i = 0; i <=1 ; i++)
    {
        
        /* select index, note ordering is reversed */
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
            smileType       = (long)(smile1[0]+0.5) / 10;
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
            smileType       = (long)(smile2[0]+0.5) / 10;
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
            NULL,
            NULL,
            (rate[i]).a) == FAILURE) goto RETURN;

        
        if (sigATM * sqrt((rate[i].s.expiry)) < Q3_MIN_VOL_CALIB) corr = 0.;

        /* initialize MQ smile. */
        if (smileType == Q3_Q2_SMILE)
        {
             if (Q3MQSmileInitFromQ2(
                          (rate[i]).s.fwdRate,
                          sigATM,
                          0.,
                          (rate[i]).s.expiry,
                          smile,
                          &(mq[i]))) goto RETURN;

        }
        else if (smileType == Q3_SV_SMILE)
        {
            if (Q3MQSmileInitFromSV(
                (rate[i]).s.fwdRate,                 
                sigATM,
                0.,
                (rate[i]).s.expiry,
                (rate[i]).s.expiryVolvol,
                1.,
                smile,
                &(sv[i]),
                &(mq[i])) == FAILURE) goto RETURN;

        }
        else if (smileType == Q3_MQ_SMILE)
        {
            /* Initialize MQ smile. */
            if (Q3MQSmileInitFromMQ(
                (rate[i]).s.fwdRate,                 
                sigATM,
                0.,
                (rate[i]).s.expiry,
                (rate[i]).s.expiry,
                1.,
                smile,
                &(sv[i]),
                &(mq[i])) == FAILURE) goto RETURN;
        }
        else
        {
            goto RETURN;
        }

            
        /* Calibrate MQ distribution to SV option prices. */
        if (Q3MQCalib(&(mq[i])) == FAILURE) goto RETURN;

        /*test ATM */
        {
            double ATMcloseForm, ATMSimpson;
            PAYOFF pf;
            if (Q3MQPricer(&(mq[i]),
                            Q3_CALL,
                            mq[i].fwdRate,
                            &ATMcloseForm) == FAILURE) goto RETURN;
            if (Q3MQDens(            
                &(mq[i]),
                -Q3_1D_NUM_STDEV,   /* start */
                Q3_1D_NUM_STDEV,    /* end   */
                Q3_1D_INT_PTS,
                grid,
                yield,
                dens,
                &densNorm) == FAILURE) goto RETURN;
            pf.strike = mq[i].fwdRate; pf.cop = 1;
            if (Q3SimpsonPricer1D(
                &pf,
                Q3Pay1D_Vnl,
                grid,
                yield,
                dens,
                densNorm,
                &ATMSimpson) == FAILURE) goto RETURN;

            if (Q3SimpsonPricer1D(
                &pf,
                Q3Pay1D_Yield,
                grid,
                yield,
                dens,
                densNorm,
                &ATMSimpson) == FAILURE) goto RETURN;
                          
        }
        
        /* Initialize FA parameters. */
        (fa[i]).mq = &(mq[i]);
        if (Q3FASmileInit(
            (rate[i]).s.expiry,
            sv[i].sigATM,
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
    
    /* Density function and yields */        
    if (mq[0].sigMQ > Q3_MIN_VOL)
    {
        if (Q3MQDens(            
            &(mq[0]),
            -Q3_1D_NUM_STDEV,   /* start */
            Q3_1D_NUM_STDEV,    /* end   */
            Q3_1D_INT_PTS,
            grid,
            yield,
            dens,
            &densNorm) == FAILURE) goto RETURN;
    }
    else
    {
        densNorm = 1.0;
        grid[0]  = 0.0;
        dens[0]  = 1.0;
        yield[0] = pa[0].fwdRate;
        if( mq[0].calcFwd == TRUE)
        {
            if(Q3MQMap(
               &(mq[0]),
               mq[0].muMQ,
               &(yield[0])) == FAILURE) goto RETURN;
        }
    }


    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        results[1+17*i] = grid[i];
        results[2+17*i] = yield[i];
        results[3+17*i] = dens[i];
        results[4+17*i] = densNorm;
    }

    if (pa[0].sigMQ > Q3_MIN_VOL)
    {
        if (Q3MQDens(            
            &(pa[0]),
            -Q3_1D_NUM_STDEV,   /* start */
            Q3_1D_NUM_STDEV,    /* end   */
            Q3_1D_INT_PTS,
            grid,
            yield,
            dens,
            &densNorm) == FAILURE) goto RETURN;
    }
    else
    {
        densNorm = 1.0;
        grid[0]  = 0.0;
        dens[0]  = 1.0;
        yield[0] = pa[0].fwdRate;
        if( pa[0].calcFwd == TRUE)
        {
            if(Q3MQMap(
               &(pa[0]),
               pa[0].muMQ,
               &(yield[0])) == FAILURE) goto RETURN;
        }
    }


    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        results[5+17*i] = grid[i];
        results[6+17*i] = yield[i];
        results[7+17*i] = dens[i];
        results[8+17*i] = densNorm;
    }

    if (mq[1].sigMQ > Q3_MIN_VOL)
    {
        if (Q3MQDens(            
            &(mq[1]),
            -Q3_1D_NUM_STDEV,   /* start */
            Q3_1D_NUM_STDEV,    /* end   */
            Q3_1D_INT_PTS,
            grid,
            yield,
            dens,
            &densNorm) == FAILURE) goto RETURN;
    }
    else
    {
        densNorm = 1.0;
        grid[0]  = 0.0;
        dens[0]  = 1.0;
        yield[0] = pa[1].fwdRate;

        if( mq[1].calcFwd == TRUE)
        {
            if(Q3MQMap(
               &(mq[1]),
               mq[1].muMQ,
               &(yield[0])) == FAILURE) goto RETURN;
        }
    }


    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        results[9+17*i]  = grid[i];
        results[10+17*i] = yield[i];
        results[11+17*i] = dens[i];
        results[12+17*i] = densNorm;
    }

    if (pa[1].sigMQ > Q3_MIN_VOL)
    {
        if (Q3MQDens(            
            &(pa[1]),
            -Q3_1D_NUM_STDEV,   /* start */
            Q3_1D_NUM_STDEV,    /* end   */
            Q3_1D_INT_PTS,
            grid,
            yield,
            dens,
            &densNorm) == FAILURE) goto RETURN;
    }
    else
    {
        densNorm = 1.0;
        grid[0]  = 0.0;
        dens[0]  = 1.0;
        yield[0] = pa[1].fwdRate;

        if( pa[1].calcFwd == TRUE)
        {
            if(Q3MQMap(
               &(pa[1]),
               pa[1].muMQ,
               &(yield[0])) == FAILURE) goto RETURN;
        }
    }


    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        results[13+17*i] = grid[i];
        results[14+17*i] = yield[i];
        results[15+17*i] = dens[i];
        results[16+17*i] = densNorm;
    }

    fstream = fopen("marginal.dat", "w");
    if (fstream == NULL)
        goto RETURN;

    fprintf(fstream, "rate1: Ann measure\n");
    fprintf(fstream, "grid\tyield\tdensity\tdensnorm\n");
    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        fprintf(fstream, "%e\t%e\t%e\t%e\n",
                results[1+17*i],
                results[2+17*i],
                results[3+17*i],
                results[4+17*i]);
    }
    fprintf(fstream, "rate1: Fwd measure\n");
    fprintf(fstream, "grid\tyield\tdensity\tdensnorm\n");
    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        fprintf(fstream, "%e\t%e\t%e\t%e\n",
                results[5+17*i],
                results[6+17*i],
                results[7+17*i],
                results[8+17*i]);
    }
    fprintf(fstream, "rate2 Ann measure\n");
    fprintf(fstream, "grid\tyield\tdensity\tdensnorm\n");
    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        fprintf(fstream, "%e\t%e\t%e\t%e\n",
                results[9+17*i],
                results[10+17*i],
                results[11+17*i],
                results[12+17*i]);
    }
    fprintf(fstream, "rate2 Fwd measure\n");
    fprintf(fstream, "grid\tyield\tdensity\tdensnorm\n");
    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        fprintf(fstream, "%e\t%e\t%e\t%e\n",
                results[13+17*i],
                results[14+17*i],
                results[15+17*i],
                results[16+17*i]);
    }
    
    fclose(fstream);


    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: failed\n", routine);  
    }
  
    return status;

} /* Q3MQFAMarginalDist */

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
    double *results
    )
{
    static char routine[] = "Q3MQMarginalDist";
    int         status    = FAILURE;

    FPAYOFF  *payFunc       = NULL;
    double   *smile         = NULL;
    double   *indxRates     = NULL;
    double   *vnfmParams    = NULL;
    long     *indxDates     = NULL;
    long      numVNFMParams = 0;
    char     *rateDCC       = NULL;

    SVDATA    sv[2];
    MQDATA    mq[2];
    double    vol;
    double    fwd;
    long      type, smileType;
    int       i;

    double grid[Q3_1D_INT_PTS];
    double yield[Q3_1D_INT_PTS];
    double dens[Q3_1D_INT_PTS];
    double densNorm;

    FILE  *fstream;

    /* calibrate PA measure for each index */
    for (i = 0; i <= 1; i++)
    {
        
        /* select index                         */
        /* note: ordering of rates is reversed! */
        if (i == 0)
        {
            type  = rateType1;
            vol   = sigATM1;
            smile = smile1;
            fwd   = fwdRate1;
            smileType = (long)(smile1[0]+0.5) / 10;
        }
        else
        {
            type  = rateType2;
            vol   = sigATM2;
            smile = smile2;
            fwd   = fwdRate2;
            smileType = (long)(smile1[0]+0.5) / 10;
        }

        /* different calibrations for rate and spread */
        switch (type)
        {
        case Q3_BV_RATE:

            if ( vol * sqrt(expiry) < Q3_MIN_VOL_CALIB) corr = 0.;

            if (smileType == Q3_Q2_SMILE)
            {
                 if (Q3MQSmileInitFromQ2(
                             fwd,
                             vol,
                             0.,
                             expiry,
                             smile,
                             &(mq[i])) == FAILURE) goto RETURN;
            }
            else if (smileType == Q3_SV_SMILE)
            {
                  if (Q3MQSmileInitFromSV(
                             fwd,
                             vol,
                             0.,
                             expiry,
                             expiry,
                             1.,
                             smile,
                             &(sv[i]),
                             &(mq[i])) == FAILURE) goto RETURN;
            }
            else if (smileType == Q3_MQ_SMILE)
            {
                if (Q3MQSmileInitFromMQ(
                        fwd,
                        vol,
                        0.,
                        expiry,
                        expiry,
                        1.,
                        smile,
                        &(sv[i]),
                        &(mq[i])) == FAILURE) goto RETURN;
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
   
    /* Density function and yields */        
    if (mq[0].sigMQ > Q3_MIN_VOL)
    {
        if (Q3MQDens(            
            &(mq[0]),
            -Q3_1D_NUM_STDEV,   /* start */
            Q3_1D_NUM_STDEV,    /* end   */
            Q3_1D_INT_PTS,
            grid,
            yield,
            dens,
            &densNorm) == FAILURE) goto RETURN;
    }
    else
    {
        densNorm = 1.0;
        grid[0]  = 0.0;
        dens[0]  = 1.0;
        yield[0] = mq[0].fwdRate;
        if( mq[0].calcFwd == TRUE)
        {
            if(Q3MQMap(
               &(mq[0]),
               mq[0].muMQ,
               &(yield[0])) == FAILURE) goto RETURN;
        }
    }


    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        results[1+9*i] = grid[i];
        results[2+9*i] = yield[i];
        results[3+9*i] = dens[i];
        results[4+9*i] = densNorm;
    }

    if (mq[1].sigMQ > Q3_MIN_VOL)
    {
        if (Q3MQDens(            
            &(mq[1]),
            -Q3_1D_NUM_STDEV,   /* start */
            Q3_1D_NUM_STDEV,    /* end   */
            Q3_1D_INT_PTS,
            grid,
            yield,
            dens,
            &densNorm) == FAILURE) goto RETURN;
    }
    else
    {
        densNorm = 1.0;
        grid[0]  = 0.0;
        dens[0]  = 1.0;
        yield[0] = mq[1].fwdRate;

        if( mq[1].calcFwd == TRUE)
        {
            if(Q3MQMap(
               &(mq[1]),
               mq[1].muMQ,
               &(yield[0])) == FAILURE) goto RETURN;
        }
    }


    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        results[5+9*i] = grid[i];
        results[6+9*i] = yield[i];
        results[7+9*i] = dens[i];
        results[8+9*i] = densNorm;
    }

    fstream = fopen("marginal.dat", "w");
    if (fstream == NULL)
        goto RETURN;

    fprintf(fstream, "rate1\n");
    fprintf(fstream, "grid\tyield\tdensity\tdensnorm\n");
    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        fprintf(fstream, "%e\t%e\t%e\t%e\n",
                results[1+9*i],
                results[2+9*i],
                results[3+9*i],
                results[4+9*i]);
    }
    fprintf(fstream, "rate2\n");
    fprintf(fstream, "grid\tyield\tdensity\tdensnorm\n");
    for (i = 0; i < Q3_1D_INT_PTS; i++)
    {
        fprintf(fstream, "%e\t%e\t%e\t%e\n",
                results[5+9*i],
                results[6+9*i],
                results[7+9*i],
                results[8+9*i]);
    }
    
    fclose(fstream);
    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3ErrMsg("%s: failed\n", routine);  
    }
  
    return status;

} /* Q3MQMarginalDist */
