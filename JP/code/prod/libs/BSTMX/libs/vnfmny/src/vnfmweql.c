#include "drlstd.h"
#include "drlvtype.h"                   /* DrlLilStructGet() */
#include "drlproc.h"                    /* DRL_PROCFLAG_ID_ADDIN_LOG */

#include "ldate.h"
#include "macros.h"
#include "tcurve.h"
#include "yearfrac.h"
#include "vnfmeqma.h"	/* prototype consistency */
#include "swaptn2l.h"                   /* Prototype consistency */

/*f-------------------------------------------------------------------------
 * <b> Add-in Function Name:</b> DR_WRAPPER\_EQMAT
 * 
 * <br><br>
 * VnfmDRWrapperMREquivalentMaturity  - SPECIAL WRAPPER FOR DR WRAPPER<br>
 *
 * This routine is a wrapper for  VnfmMREquivalentMaturity, which takes all 
 * the specific deal inforation from the same inputs as SWAPTION\_FLOWS2.
 * (for description of the routine see VnfmMREquivalentMaturity).\\
 * 
 * <br>
 * <br>[1]  Value date
 * <br>[2]  Zero curve dates
 * <br>[3]  Zero curve rates
 * <br>[4]  Array of dates on which equivalent maturities are needed.
 * <br>[5]  Array of coupon dates.
 * <br>[6]  Array of coupon payments.
 * <br>[7]  Array of principal payment dates.
 * <br>[8]  Array of principal payments.
 * <br>[9]  Floating accrue start dates 
 * <br>[10] Floating pay dates ( accrue end is assumed to be pay date).
 * <br>[11] Floating notionals. 
 * <br>[12] Index weights (the last element may contain the spread over the 
 *           floating leg). 
 * <br>[13] IndexFreqs. If the number of index frequencies is less than the 
 *           number of index weights the last element of index weights is the 
 *           spread over the floating leg.
 * <br>[14] UnderConvs: containes the daycount convention for the floating leg.
 * <br>[15] Exercise dates: containes the first exercise date.
 * <br>[16] Strike amounts: containes the first strike amount.
 * <br>[17] Coupons per year in the benchmark swaps. 
 * <br>[18] FDaycount convention for vanilla benchmark swaps.
 * <br>[19] (1) first mean reversion.
 * 	  (2) second mean reversion (>first mean reversion. NOT for a 2nd fact).
 * 	  (3) backboneq, the power of the interest rate that is taken before
 *                 multiplying by the volatility to obtain the basis point
 *                 volatility. (0 for lognormal, 0.5 for normal, between 0 and 0.5
 *                 cevPower = 1-2*q for other CEV processes).
 * <br>[20] (1) Minimum maturity.  The routine will output this maturity if 
 * 	      the equivalent maturity is shorter.
 * 	  (2) Maximum maturity. The routine will output this maturity if the 
 *                         equivalent maturity is larger than it. 
 *                         Note that it is important to enter here a relevant no.
 *                         since maybe no swap gives similar behaviour and no 
 *                         maturity will be found. In this case one should use 
 *                         the longest maturity available.
 * <br>[21] An array of equivalent maturities is the output of this routine.
 * <br> 
 * 
 */

DLL_EXPORT(int)
VnfmDRWrapperMREquivalentMaturityL(
        long    *valueDateL,	       /*  1 'D' (I) */
        long    *zeroDatesL,	       /*  2 'D' (I) */
        double  *zeroRatesL,	       /*  3 'F' (I) */
    	long    *resetDatesL,	       /*  4 'D' (I) */
	long    *couponDatesL,	       /*  5 'D' (I) */
        double  *couponPaymentsL,      /*  6 'F' (I) */
        long    *principalDatesL,      /*  7 'D' (I) */
        double  *principalPaymentsL,   /*  8 'F' (I) */
        long    *floatAccrueDatesL,    /*  9 'D' (I) */
        long    *floatPayDatesL,       /* 10 'D' (I) */
        double  *floatNotionalsL,      /* 11 'F' (I) */
        double  *indexWeightsL,        /* 12 'F' (I) */
                                       /* last element may be a spread */
        long    *indexFreqsL,          /* 13 'L' (I) */
                                       /* to determine if a spread is given */
        char    *underConvsL,          /* 14 'C' (I) */
        long    *exerDatesL,           /* 15 'D' (I) */
        double  *exerAmountsL,         /* 16 'F' (I) */
        long    *swapPayFreqL,	       /* 17 'L' (I) */
        char    *swapDayCountConvL,    /* 18 'C' (I) */
                                       /*  ACT/365   */
                                       /*  ACT/365F  */
                                       /*  ACT/360   */
                                       /*  30/360    */
                                       /*  30E/360   */
                                       /*  ACT/365J  */
                                       /*  30E/360I  */
        double  *mrL,                  /* 19 'F' (I) */
        double  *minmaxMaturityL,      /* 20 'F' (I) */
     	double  *maturitiesL)	       /* 21 'F' (O) */
{
    TCurve   *zeroCurve = NULL;   
    long numZeroDates;
    long numCoupons;
    long numPrincipals;
    long swapDayCountConv;
    long numFloatPayments;
    double spread;
    long floatPayDayCountConv;
    long couponIdx;
    long floatIdx;
    double notional;
    double sum;
    double yearFrac;
    double discount;
    double *bondPrincipalPayments = NULL;
    long   principalIdx;
    double floatPrincipalPay;
    double strikeAdjustment;
    TDate  *couponDates;   
    double *couponPayments; 
    TDate  *principalDates; 
    double *principalPayments;
    TDate  *floatAccrueDates; 
    TDate  *floatPayDates;    
    double *floatNotionals;   

    double mr1, mr2, backboneq, minMat, maxMat;

    long   status = FAILURE;
    static char routine[] = "VnfmDRWrapperMREquivalentMaturityL";


    /* Check length of arrays coming from Lotus
     */
    if (NOT (
         IS_VECTOR(floatAccrueDatesL)   &&
         IS_VECTOR(floatPayDatesL)      &&
         IS_VECTOR(indexWeightsL)       &&
         IS_VECTOR(indexFreqsL)       ))
    {
        GtoWrapVectorErrMsg(routine);
        return(FAILURE);
    }

    
   if (!IS_SCALAR(valueDateL))
   {
       GtoErrMsg("%s: value date can't be a range.\n", routine);
       goto error;
   }

   numCoupons = (long)couponDatesL[0];
   if (numCoupons != (long)couponPaymentsL[0])
   {
    GtoErrMsg("%s: inconsistent number of coupon dates and payments"
    ".\n",routine);
       goto error;
   }

   numPrincipals = (long)principalDatesL[0];
   if (numPrincipals != (long)principalPaymentsL[0])
   {
    GtoErrMsg("%s: inconsistent number of principal dates and payments"
    ".\n",routine);
       goto error;
   }

   if (ARGSIZE(maturitiesL) < ARGSIZE(resetDatesL)) 
   {
       GtoErrMsg("%s: inconsistent number of reset dates and maturities.\n",
                  routine);
       goto error;
   }

   if ( exerDatesL[0] < 1 )
   {
       GtoErrMsg("%s: no exercise dates were given\n",
                  routine);
       goto error;
   }

   if ( exerAmountsL[0] < 1 )
   {
       GtoErrMsg("%s: no exercise amounts were given\n",
                  routine);
       goto error;
   }

   if (!IS_SCALAR(swapPayFreqL))
   {
       GtoErrMsg("%s: swapPayFreq can't be a range.\n",routine);
       goto error;
   }

   if (!IS_SCALAR(swapDayCountConvL))
   {
       GtoErrMsg("%s: swapDaycountConv can't be a range.\n",routine);
       goto error;
   }

   if (mrL[0] != 3)
   {
       GtoErrMsg("%s: mr must contains mr1, mr2 and backbone q.\n",routine);
       goto error;
   }

   mr1 = mrL[1];
   mr2 = mrL[2];
   backboneq = mrL[3];

   if ( mr1 < 0.)
   {
       GtoErrMsg("%s: mr1=%12.6f  can't be a negative.\n",routine,mrL[1]);
       goto error;
   }

   if ( mr2 < 0.)
   {
       GtoErrMsg("%s: mr2=%12.6f  can't be a negative.\n",routine,mrL[2]);
       goto error;
   }

   if ( mr2 <= mr1)
   {
       GtoErrMsg("%s: mr2=%12.6f  has to be larger than mr1=%12.6f\n",routine,
                 mr2,mr1);
       goto error;
   }

   if (minmaxMaturityL[0] != 2)
   {
       GtoErrMsg("%s: minmaxMaturity must have two arguments, minimumMaturity"
		 " and maximumMaturity.\n",routine);
       goto error;
   }

   minMat = minmaxMaturityL[1];
   maxMat = minmaxMaturityL[2];

   /* Log wrapper inputs  */
   if (GtoLoggingGet() > 0) {
	DrlLilVectLoggingFile("wrapper.log", "w", "DR_WRAPPER_EQMAT",
                DRL_TDATE_L,  valueDateL,         "VALUE_DATE",
                DRL_TDATE_L,  zeroDatesL,         "ZERO_DATES",
                DRL_FLOAT_L,  zeroRatesL,         "ZERO_RATES",
 
                DRL_TDATE_L,  resetDatesL,         "RESET_DATES",
                DRL_TDATE_L,  couponDatesL,        "COUPON_DATES",
                DRL_FLOAT_L,  couponPaymentsL,     "COUPON_PAYMENTS",
                DRL_TDATE_L,  principalDatesL,     "PRINCIPAL_DATES",
                DRL_FLOAT_L,  principalPaymentsL,  "PRINCIPAL_PAYMENTS",
                DRL_TDATE_L,  floatAccrueDatesL,   "FLOAT_ACCRUE_DATES",
                DRL_TDATE_L,  floatPayDatesL,      "FLOAT_PAY_DATES",
                DRL_FLOAT_L,  floatNotionalsL,     "FLOAT_NOTIONALS",
                DRL_FLOAT_L,  indexWeightsL,       "INDEX_WEIGHTS",
 
                DRL_LONG_L,   indexFreqsL,         "INDEX_FREQS",

                DRL_CHAR_BLOCK_L, underConvsL,     "UNDER_CONVS",
                DRL_TDATE_L,  exerDatesL,          "EXER_DATES",
                DRL_FLOAT_L,  exerAmountsL,        "EXER_AMOUNTS",
                DRL_LONG_L,   swapPayFreqL,        "SWAP_PAY_FREQ",
                DRL_CHAR_BLOCK_L, swapDayCountConvL, "SWAP_DAY_COUNT_CONV",
                DRL_FLOAT_L,  mrL,                 "MR",
 
                DRL_FLOAT_L,  minmaxMaturityL,     "MINMAX_MATURITY",
                0);
   }

    if (GtoConvertTCurveWrap(
                   valueDateL,
                   zeroDatesL,
                   zeroRatesL,
	            	   1L,
                   GTO_ACT_365F,
		               routine,
                   &zeroCurve) != SUCCESS)
				goto error;


   if (zeroCurve == NULL)
       {
           GtoErrMsg("%s: Error setting up zero TCurve structure.\n", routine);
           goto error;
       }

   if ( GtoStringToDayCountConv(&swapDayCountConvL[1],&swapDayCountConv) 
                       != SUCCESS)
       {
           GtoErrMsg("%s: conversion of daycount string to long failed.\n", 
                                               routine);
           goto error;
       }

      numFloatPayments = floatPayDatesL[0];

      if (numFloatPayments != floatNotionalsL[0] ||
          numFloatPayments != floatAccrueDatesL[0])
      {
        GtoErrMsg("%s: number of floating payments accrual starts and"
                     " notionals is not consistent\n",routine);
        return(FAILURE);
      }     

      couponDates       = couponDatesL+1;
      couponPayments    = couponPaymentsL+1;
      principalDates    = principalDatesL+1;
      principalPayments = principalPaymentsL+1;
      floatAccrueDates  = floatAccrueDatesL+1;
      floatPayDates     = floatPayDatesL+1;
      floatNotionals    = floatNotionalsL+1;

    if ( numFloatPayments > 0 ) 
   {
       /* Floating leg is modeled. Since eqmat assumes that it is given a 
          bond the principal payments have to be adjusted and the coupons 
          should be adjusted if there is a spread.                       */
    
    if ( indexFreqsL[0]+1 == ((long) indexWeightsL[0]))
    {
      spread = indexWeightsL[((long)indexWeightsL[0])];

      if ((long)underConvsL[0] ISNT GTO_SW2_NUM_UNDER_CONVS)
      {
        GtoErrMsg("%s: Received %ld underlying conventions, instead of %d.\n",
                  routine, (long)underConvsL[0], GTO_SW2_NUM_UNDER_CONVS);
        return(FAILURE);
      }

      if (GtoStringToDayCountConv
        (&underConvsL[WRAP_STR_IDX(GTO_SW2_FL_PAY_DCC_IDX)],
         &floatPayDayCountConv) IS FAILURE)
        goto error;



   
      for (couponIdx = 0; couponIdx < numCoupons &&
                          couponDates[couponIdx] < floatPayDates[0];
                                                               couponIdx++);
      if ( couponIdx >= numCoupons)
      {
         GtoErrMsg("%s: first floating payment is past last fixed payment\n", 
                     routine);
           goto error;
       }        

      for ( floatIdx = 0; 
            couponIdx < numCoupons && floatIdx < numFloatPayments;
                                couponIdx++)
      {
       notional = floatNotionals[floatIdx]; 
       sum = 0.;
       while ( floatIdx < numFloatPayments &&
               floatPayDates[floatIdx] <= couponDates[couponIdx] )
       {
         if ( GtoDayCountFraction( floatAccrueDates[floatIdx],
                                   floatPayDates[floatIdx],
                                   floatPayDayCountConv,
                                  &yearFrac) ISNT SUCCESS) goto error; 

         if (GtoDiscountDate( floatPayDates[floatIdx],
                            zeroCurve,
                            GTO_LINEAR_INTERP,
                           &discount     ) ISNT SUCCESS ) goto error;

 
         sum += yearFrac * discount;  
         
         if ( floatNotionals[floatIdx] != notional )
         {
          GtoErrMsg("%s: notional cannot change between fixed coupon dates \n",                     routine);
           goto error;
         } 
  
         floatIdx++;
       }
       if (floatPayDates[floatIdx-1] != couponDates[couponIdx] && notional > 0)
       {
         GtoErrMsg("%s: coupon date no. %d is not a floating date\n", 
                     routine,couponIdx);
           goto error;
       } 
       couponPayments[couponIdx] -= notional * spread * (sum / discount);   
 
      }
    }

    bondPrincipalPayments = NEW_ARRAY(double, numCoupons);
    if ( bondPrincipalPayments IS NULL) goto error;    
    
    couponIdx = 0;
    for (principalIdx = 0; principalIdx < numPrincipals; principalIdx++)
    {
     if (principalPayments[principalIdx] != 0.)
     {
      while (couponIdx < numCoupons 
             && couponDates[couponIdx] < principalDates[principalIdx] )        
                                                             couponIdx++;
      if (  couponDates[couponIdx] != principalDates[principalIdx])
      {   
        GtoErrMsg("%s: principal dates must be fixed coupon dates \n",      
               routine);
           goto error;
      }

      bondPrincipalPayments[couponIdx] += principalPayments[principalIdx];
     }
    }

    if ( couponDates[0] > floatAccrueDates[0])
      {   
        GtoErrMsg("%s: first floating acrrue start must be on or after"
                  "first coupon date (this first coupon date is actually"
                  " the accrual start for the fixed leg).",      
               routine);
           goto error;
      }

    for (couponIdx = 0; couponIdx < numCoupons 
             && couponDates[couponIdx] < floatAccrueDates[0];
                                                         couponIdx++);
    if ( exerDatesL[1] <  couponDates[couponIdx])
    {
      bondPrincipalPayments[couponIdx] -= floatNotionals[0];
      strikeAdjustment = 0;
    }
    else
    {
     for (floatIdx = 0; floatIdx < numFloatPayments &&
                      floatPayDates[floatIdx] < exerDatesL[1] ; floatIdx++);

     strikeAdjustment = floatNotionals[floatIdx];
    }


    couponIdx = 0;
    for (floatIdx = 0; floatIdx < numFloatPayments; floatIdx++)
    {
     if (floatIdx == numFloatPayments-1)
     {
      floatPrincipalPay =  floatNotionals[floatIdx];
     }
     else
     {
      floatPrincipalPay = floatNotionals[floatIdx]- floatNotionals[floatIdx+1];
     }
     if ( floatPrincipalPay !=0 )
     {
      while (couponIdx < numCoupons 
             && couponDates[couponIdx] < floatPayDates[floatIdx] )        
                                                             couponIdx++;
      if (  couponDates[couponIdx] != floatPayDates[floatIdx])
      {   
        GtoErrMsg("%s: float amortization dates must be fixed coupon dates \n",      
               routine);
           goto error;
      }

      bondPrincipalPayments[couponIdx] += floatPrincipalPay;
     }
    }
     principalPayments = bondPrincipalPayments;
   }


   if ( VnfmMREquivalentMaturity( zeroCurve,
                                   (long) resetDatesL[0],
                                  &resetDatesL[1],
                                   numCoupons,
                                   couponDates,
                                   couponPayments,
                                   numCoupons,
                                   couponDates,
                                   principalPayments,
                                   exerDatesL[1],
                                   exerAmountsL[1]+strikeAdjustment,
                                   swapPayFreqL[1],
                                   swapDayCountConv,
                                   mr1,
                                   mr2,
                                   backboneq,
                                   minMat,
                                   maxMat,
                                  &maturitiesL[1]
                                          ) != SUCCESS ) goto error;
   maturitiesL[0] = resetDatesL[0];

   status = SUCCESS;

   error:
 
   GtoFreeTCurve(zeroCurve);

   if ( bondPrincipalPayments != NULL) FREE(bondPrincipalPayments);

   return status;
}
