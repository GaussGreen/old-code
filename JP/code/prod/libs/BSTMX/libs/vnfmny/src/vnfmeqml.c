#include "drlstd.h"
#include "ldate.h"
#include "macros.h"
#include "tcurve.h"
#include "yearfrac.h"
#include "vnfmeqma.h"	/* prototype consistency */

/*f-------------------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_EQMAT
 * 
 * <br><br>
 * This routine is a wrapper for  VnfmMREquivalentMaturity. 
 *  <br>
 * The inputs to the routine are:
 * <br>
 * <br>[1]  Value Date
 * <br>[2]  Zero curve Dates
 * <br>[3]  Zero curve Rates
 * <br>[4]  Array of dates on which equivalent maturities are needed.
 * <br>[5]  Array of coupon dates.
 * <br>[6]  Array of coupon payments.
 * <br>[7]  Array of principal payment dates.
 * <br>[8]  Array of principal payments.
 * <br>[9]  First Exercise date.
 * <br>[10] First Strike amount.
 * <br>[11] Coupons per year in the benchmark swaps.
 * <br>[12] Daycount convention for benchmark swaps. 
 * <br>[13] First mean reversion.
 * <br>[14] Second mean reversion.( > than first mean reversion this is NOT for 
 *           a second factor).
 * <br>[15] backboneq, the power of the interest rate that is taken before
 *                 multiplying by the volatility to obtain the basis point
 *                 volatility. (0 for lognormal, 0.5 for normal, between 0 and 0.5
 *                 cevPower = 1-2*q for other CEV processes).
 * <br>[16] Minimum maturity. The routine will output this maturity if the 
 *           equivalent maturity is shorter.
 * <br>[17] Maximum maturity. The routine will output this maturity if the 
 *           equivalent maturity is larger than it. Note that it is important to 
 *           enter here a relevant no. since maybe no swap gives similar 
 *           behaviour and no maturity will be found. In this case one should 
 *           use the longest maturity available.
 * <br>[18] An array of equivalent maturities is the output of this routine.
 * <br> 
 */

DLL_EXPORT(int)
VnfmMREquivalentMaturityL(
	long    *valueDateL,	       /*  1 'D' (I) */
	long    *zeroDatesL,	       /*  2 'D' (I) */
	double  *zeroRatesL,	       /*  3 'F' (I) */
	long    *resetDatesL,	       /*  4 'D' (I) */
	long    *couponDatesL,	       /*  5 'D' (I) */
	double  *couponPaymentsL,      /*  6 'F' (I) */
	long    *principalDatesL,      /*  7 'D' (I) */
	double  *principalPaymentsL,   /*  8 'F' (I) */
	long    *firstExerDateL,       /*  9 'D' (I) */
	double  *exerAmountL,          /* 10 'F' (I) */
	long    *swapPayFreqL,	       /* 11 'L' (I) */
	char    *swapDayCountConvL,    /* 12 'C' (I) ACT/365, etc. */
	double  *mr1L,                 /* 13 'F' (I) */
	double  *mr2L,                 /* 14 'F' (I) */
	double  *backboneqL,           /* 15 'F' (I) */
	double  *minimumMaturityL,     /* 16 'F' (I) */
	double  *maximumMaturityL,     /* 17 'F' (I) */
	double  *maturitiesL)	       /* 18 'F' (O) */
{
    TCurve   *zeroCurve = NULL;   
    long numZeroDates;
    long numCoupons;
    long numPrincipals;
    long swapDayCountConv;

    long   status = FAILURE;
    static char routine[] = "VnfmMREquivalentMaturityL";

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

   if (!IS_SCALAR(firstExerDateL))
   {
       GtoErrMsg("%s: firstExerDate can't be a range.\n", routine);
       goto error;
   }
   if (!IS_SCALAR(exerAmountL))
   {
       GtoErrMsg("%s: exerAmount can't be a range.\n", routine);
       goto error;
   }

   if (ARGSIZE(maturitiesL) < ARGSIZE(resetDatesL)) 
   {
       GtoErrMsg("%s: inconsistent number of reset dates and maturities.\n",
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

   if (!IS_SCALAR(mr1L))
   {
       GtoErrMsg("%s: mr1 can't be a range.\n",routine);
       goto error;
   }
   if ( mr1L[1] < 0.)
   {
       GtoErrMsg("%s: mr1=%12.6f  can't be a negative.\n",routine,mr1L[1]);
       goto error;
   }

   if (!IS_SCALAR(mr2L))
   {
       GtoErrMsg("%s: mr2 can't be a range.\n",routine);
       goto error;
   }
   if ( mr2L[1] < 0.)
   {
       GtoErrMsg("%s: mr2=%12.6f  can't be a negative.\n",routine,mr2L[1]);
       goto error;
   }

   if ( mr2L[1] <= mr1L[1])
   {
       GtoErrMsg("%s: mr2=%12.6f  has to be larger than mr1=%12.6f\n",routine,
                 mr2L[1],mr1L[1]);
       goto error;
   }
   if (!IS_SCALAR(backboneqL))
   {
       GtoErrMsg("%s: backboneq can't be a range.\n",routine);
       goto error;
   }
   if (!IS_SCALAR(minimumMaturityL))
   {
       GtoErrMsg("%s: minimumMaturity can't be a range.\n",routine);
       goto error;
   }
   if (!IS_SCALAR(maximumMaturityL))
   {
       GtoErrMsg("%s: maximumMaturity  can't be a range.\n",routine);
       goto error;
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

   if ( VnfmMREquivalentMaturity( zeroCurve,
                                   (long) resetDatesL[0],
                                  &resetDatesL[1],
                                   numCoupons,
                                  &couponDatesL[1],
                                  &couponPaymentsL[1],
                                   numPrincipals,
                                  &principalDatesL[1],
                                  &principalPaymentsL[1],
                                   firstExerDateL[1],
                                   exerAmountL[1],
                                   swapPayFreqL[1],
                                   swapDayCountConv,
                                   mr1L[1],
                                   mr2L[1],
                                   backboneqL[1],
                                   minimumMaturityL[1],
                                   maximumMaturityL[1],
                                  &maturitiesL[1]
                                          ) != SUCCESS ) goto error;
   maturitiesL[0] = resetDatesL[0];

   status = SUCCESS;

   error:
 
   GtoFreeTCurve(zeroCurve);

   return status;
}
