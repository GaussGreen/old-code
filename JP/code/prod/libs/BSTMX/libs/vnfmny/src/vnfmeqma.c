#include <math.h>
#include "drlstd.h"
#include "ldate.h"
#include "cashflow.h"
#include "stub.h"
#include "convert.h"
#include "dateconv.h"
#include "macros.h"
#include "busday.h"
#include "tcurve.h"
#include "interp.h"
#include "vnfmeqma.h"		/* prototype consistency */

long vnfmBetaParBond (long    initialize,
                      TCashFlowList *cfl,
                      TCurve *zeroCurve,
                      TDate   resetDate,
                      long    lastCouponIdx,
                      long    frequency,
                      double  beta1,
                      double  beta2,
                      double  backboneq,
                      double *value1,
                      double *value2);

long  vnfmBetaBondGeneral (TCurve *zeroCurve,
                       TDate  resetDate,
               	       long numCouponDates,	/* (I) */
	                     TDate *couponDates,	/* (I) */
                       double *couponPayments,           /*(I)  */
                       long    numPrincipalDates,
                       TDate  *principalDates,
                       double *principalPayments,
                       double  beta1,
                       double  beta2,
                       double  backboneq,
                       double *value1,
                       double *value2);

/*f--------------------------------------------------------------------------
 * Equivalent maturity calculation.
 *                                                             
 * <br><br>
 *  Given a generalised bond (or equivalently a generalized swap where the 
 *  floating leg always values to par), we compute the maturity of a vanilla swap
 *  that has the same sensitivity to mean reversion as the generalised bond 
 *  (or swap). 
 *  We do this for each of the given reset dates, i.e. we compute a maturity 
 *  that corresponds to part of the bond that remains past the reset date. 
 *  These maurities will be the ones one should calibrate to, on the 
 *  corresponding reset dates.
 *   E.g if we are given a vanilla par swap the equivalent maturity for a reset 
 *  date of k years from today would be M-k, where M is the final maturity of the
 *  deal in years from today (this is a diagonal on the swaption matrix). If we  
 *  had a zero coupon bond with final maturity M we would get a maturity which is
 *  larger than M-k for reset in K years. If we had an amortizing bond we would 
 *  get a maturity which is smaller than M-k for reset in year k.\\
 * 
 *  The choice of an equivalent maturity depends on the mean reversion and on the
 *  specification of the benchmark swaps, i.e their coupon frequency and daycount
 *  convention. This routine takes two mean reversions as an input. It actually 
 *  finds an equivalent maturity swap that satisfies the following:\\
 *   Given a variance for the price of an equivalent swap and given a mean 
 *  reversion, a constant spot volatility can be computed. Thus if we keep the 
 *  variance of the swap constant we will have a spot vol for each choice of 
 *  mean reversion. 
 *   Given a spot volatility and mean reversion we can approximate the variance 
 *  of the price of the generalized bond. Thus using two mean reversions we get 
 *  two variances for the generalised bond. We pick an equivalent maturity that 
 *  will cause these two variances to coincide. What this will mean is that if 
 *  mean reversion is in the vicinity of the interval 
 *  (mean reversion 1, mean reversion 2), and if this interval is not large, 
 *  a change in volatility of the generalised bond, should be observed by a 
 *  change of the volatility of the equivalent swap. The relative change in vol 
 *  for the bench mark swap will approximately be the relative change in the 
 *  volatility for the generalized bond. Thus mean reversions 1 and 2 should 
 *  bracket the estimation of the one factor mean reversion and be relatively 
 *  close to it. \\ 
 * 
 * the inputs to the routine are:
 * <br>
 * <br>[1]  Zero curve
 * <br>[2]  Number of reset dates: 
 *           the number of dates on which equivalent maturities are needed.
 * <br>[3]  Array of dates on which equivalent maturities are needed.
 * <br>[4]  Number of coupon dates.
 * <br>[5]  Array of coupon dates.
 * <br>[6]  Array of coupon payments.
 * <br>[7]  Number of principal Payments.
 * <br>[8]  Array of principal payment dates.
 * <br>[9]  Array of principal payments.
 * <br>[10] First Exercise date.
 * <br>[11] Strike amount.
 * <br>[12] Coupons per year in the benchmark swaps. 
 * <br>[13] First mean reversion.
 * <br>[14] Second mean reversion.
 *           ( > than first mean reversion. This is NOT for a second factor).
 * <br>[15] backboneq, the power of the interest rate that is taken before 
 *                 multiplying by the volatility to obtain the basis point 
 *                 volatility. (0 for lognormal, 0.5 for normal, between 0 and 0.5 
 *                 cevPower = 1-2*q for other CEV processes).
 * <br>[16] Minimum maturity.  The routine will output this maturity 
 *                          if the equivalent maturity is shorter.
 * <br>[17] Maximum maturity. The routine will output this maturity if the 
 *                         equivalent maturity is larger than it. 
 *                         Note that it is important to enter here a relevant no.
 *                         since maybe no swap gives similar behaviour and no 
 *                         maturity will be found. In this case one should use 
 *                         the longest maturity available.
 * <br>[18] An array of equivalent maturities is the output of this routine.
 * <br> 
 */

DLL_EXPORT(long)
VnfmMREquivalentMaturity(
  TCurve *zeroCurve,		/* (I) */
  long    numResetDates,	/* (I) */
  TDate  *resetDates,	        /* (I) */
  long    numCouponDates,	/* (I) */
  TDate  *couponDates,	      /* (I) */
  double *couponPayments,     /* (I) */
  long    numPrincipalDates,  /* (I) */
  TDate  *principalDates,     /* (I) */
  double *principalPayments,  /* (I) */
  TDate   firstExerDate,      /* (I) */
  double  exerAmount,         /* (I) */
  long    swapPayFreq,	      /* (I) */
  long    swapDayCountConv,   /* (I) */
                              /* GTO_ACT_365    1  Actual/365 */
                              /* GTO_ACT_365F   2  Actual/365 Fixed */
                              /* GTO_ACT_360    3  Actual/360 */
                              /* GTO_B30_360    4  30/360 */
                              /* GTO_B30E_360   5  30E/360 */
                              /* GTO_ACT_365FJ  6     */             
                              /* GTO_B30E_360I  7 */
                              /* GTO_B30_360_FIXED 8 For bond coupon payments*/
  double  mr1,                /* (I) mr1<mr2  */
  double  mr2,                /* (I) */
  double  backboneq,          /* (I) */
  double  minimumMaturity,    /* (I) */
  double  maximumMaturity,    /* (I) */
  double *maturities)	        /* (O) */
{
 TDateInterval matMonths;
 long resetIdx;
 long cflIdx;
 long principalIdx;
 double value1;
 double value2;
 double  ratioGeneral;
 double ratioRegular;
 double oldRatio;
 double maturity;
 TCashFlowList *regularCfl = NULL;
 TDateInterval betCoupons;
 TDate maturityDate;
 double lastCouponIdx;
 TDate *ammendedPrincipalDates = NULL;
 double *ammendedPrincipalPayments = NULL;
 long numAmmendedPrincipalDates;
 TDate *ammendedCouponDates = NULL;
 double *ammendedCouponPayments = NULL;
 long numAmmendedCouponDates;
 long idx;
 long couponIdx;

 static   char routine[] = "MREquivalentMaturity";
 int status = FAILURE;

/*
  printf(" exer date %d  strike %12.0f\n\n",firstExerDate,exerAmount);
               

 for (idx=0; idx<numPrincipalDates; idx++)
 {
  printf(" %d principal date %d  payment %12.0f \n",idx,principalDates[idx],
                principalPayments[idx]);
 }
  printf(" \n\n\n");

 for (idx=0; idx<numCouponDates; idx++)
 {
  printf(" %d coupon date %d  payment %12.0f \n",idx,couponDates[idx],
                couponPayments[idx]);
 }
 fflush(stdout);
*/

 if (mr1 >= mr2)
 {
    GtoErrMsg("\n%s: mr1 >= mr2\n", routine);
    goto done;
 }

 if (backboneq < 0. || backboneq > 0.5)
 {
    GtoErrMsg("\n%s: backboneq = %f isn't in the interval [0,0.5]", 
                                                           routine,backboneq);
    goto done;
 }


 if (GtoMakeDateInterval(12*maximumMaturity, 'M', 
                               &matMonths) ISNT SUCCESS)
 {
         GtoErrMsg("\n%s: Cannot compute matMonths\n", routine);
         goto done;
 }

 if (GtoMakeDateInterval(12/swapPayFreq, 'M', 
                               &betCoupons) ISNT SUCCESS)
 {
         GtoErrMsg("\n%s: Cannot compute couponInterval\n", routine);
         goto done;
 }


/* if the first exercise date is past the first date on which a maturity is 
   requested a principal payment has to be added on that date and thus the 
   list might have to be extended */

 ammendedPrincipalPayments = NEW_ARRAY(double, numPrincipalDates + 1);
     if ( ammendedPrincipalPayments IS NULL) goto done; 

 ammendedPrincipalDates = NEW_ARRAY(TDate, numPrincipalDates + 1);
     if ( ammendedPrincipalDates IS NULL) goto done;

 ammendedCouponPayments = NEW_ARRAY(double, numCouponDates);
     if ( ammendedCouponPayments IS NULL) goto done; 

 ammendedCouponDates = NEW_ARRAY(TDate, numCouponDates);
     if ( ammendedCouponDates IS NULL) goto done;


 /* find the first coupon past first exercise date */
 for (couponIdx = 0; couponIdx < numCouponDates 
                        && couponDates[couponIdx] <= firstExerDate; 
                                                               couponIdx++);

     /* we are interested only on payments past first exercise date. 
         Previous coupon date is important for stub calculations.
         We delete all cash flows on or before first exercise date since we 
         want to recieve the same maturities for resets before first exercise 
         date, as we have obtained if we modelled only the forward starting 
	 swap */

 if (couponIdx > 0) /* First coupon date is on or before first exercise date */
 {
  ammendedCouponDates[0] =  couponDates[couponIdx-1];
  ammendedCouponPayments[0] =  0.;
  numAmmendedCouponDates = numCouponDates - couponIdx + 1;
  idx = 1;
 }
 else 
 {
  numAmmendedCouponDates = numCouponDates - couponIdx;
  idx = 0;
 }

 for (; couponIdx < numCouponDates; couponIdx++ , idx++)
 {
  ammendedCouponDates[idx] =  couponDates[couponIdx];
  ammendedCouponPayments[idx] =  couponPayments[couponIdx];
 }

 for (principalIdx = 0; principalIdx < numPrincipalDates 
                        && principalDates[principalIdx] < firstExerDate; 
                                                               principalIdx++);

 if ( principalDates[principalIdx] < firstExerDate)
 {
         GtoErrMsg("\n%s: firstExerDate cannot be after last principal date\n",
                                                                  routine);
         goto done;
 }     


 if ( principalDates[principalIdx] == firstExerDate )
 {
     /* first exercise date falls on a principal date. The effective principal
        payment is thus the payment minus the strike price */

    numAmmendedPrincipalDates = numPrincipalDates - principalIdx;

    for (idx=0; principalIdx < numPrincipalDates;  principalIdx++, idx++)
    {
     ammendedPrincipalPayments[idx] = principalPayments[principalIdx];
     ammendedPrincipalDates[idx]    = principalDates[principalIdx];
    }
    ammendedPrincipalPayments[0] -= exerAmount;
 }
 else
 {
     /* first exercise date does not fall on a principal date. We add an 
        effective payment of minus the strike price on this date */
     ammendedPrincipalPayments[0] = - exerAmount;
     ammendedPrincipalDates[0]    = firstExerDate;
     numAmmendedPrincipalDates = numPrincipalDates - principalIdx + 1;    

     for (idx=1; principalIdx < numPrincipalDates;  principalIdx++, idx++)
     {
     ammendedPrincipalPayments[idx] = principalPayments[principalIdx];
     ammendedPrincipalDates[idx]    = principalDates[principalIdx];
     }

 }

  
for (resetIdx = 0; resetIdx < numResetDates; resetIdx++)
  {
    /* compute the standard deviation of the generalized bond devided by the
       spot volatility for the two mean reversions */ 

   if( vnfmBetaBondGeneral( zeroCurve,
                        resetDates[resetIdx],
                        numAmmendedCouponDates,
                        ammendedCouponDates,
                        ammendedCouponPayments,
                        numAmmendedPrincipalDates,
                        ammendedPrincipalDates,
                        ammendedPrincipalPayments,
                        mr1,
                        mr2,
                        backboneq,
                       &value1,
                       &value2) ISNT SUCCESS ) goto done;

   if ( value1 == 0e0 ) 
   {
    maturities[resetIdx] = minimumMaturity;
    continue;
   }

   ratioGeneral = value2/value1;

   ratioRegular = 1.;
   lastCouponIdx = -1;


    if( GtoDtFwdAny (resetDates[resetIdx],
                     &matMonths,
                     &maturityDate) ISNT SUCCESS)
    {
         GtoErrMsg("\n%s: Cannot compute maturity date\n", routine);
         goto done;
    }
   
   /* create the coupon dates and daycount fractions for the longest allowed 
      equivalent regular swap */

   if( (regularCfl = GtoMakeCFL( 1,
                       resetDates[resetIdx],
                      &betCoupons,
                       maturityDate,
                       swapDayCountConv,
                       GTO_STUB_BOND,
                       TRUE,
                       0,
                       GTO_BAD_DAY_NONE,
                       GTO_BAD_DAY_NONE,
                       "No_Weekends"))  IS NULL )
   {
    GtoErrMsg("\n%s: Cannot compute CFL \n", routine);
                goto done;
   }

   /* initialize vnfmBetaParBond */
   if( vnfmBetaParBond(1,
                     regularCfl,
                     zeroCurve,
                     resetDates[resetIdx],
                     lastCouponIdx,
                     swapPayFreq,
                     mr1,
                     mr2,
                     backboneq,
                    &value1,
                    &value2) ISNT SUCCESS) goto done;

   /* find the two consecutive (with maturities differing by one coupon period)
      regular bonds such that the ration of the standard deviations of each of
      them with respect to the chage in mean reversion, brackets that of the 
      general bond */
 
   while (ratioRegular > ratioGeneral && lastCouponIdx < regularCfl->fNumItems)
   {
    lastCouponIdx++; 
    oldRatio = ratioRegular;
    
    if( vnfmBetaParBond(0,
                     regularCfl,
                     zeroCurve,
                     resetDates[resetIdx],
                     lastCouponIdx,
                     swapPayFreq,
                     mr1,
                     mr2,
                     backboneq,
                    &value1,
                    &value2) ISNT SUCCESS) goto done;

    ratioRegular = value2/value1;
   }

   if ( lastCouponIdx >= regularCfl->fNumItems )
   {
     maturities[resetIdx] = maximumMaturity;
   }
   else
   {
     maturities[resetIdx] = MAX( (((double) lastCouponIdx+1)/swapPayFreq)
        - (ratioRegular - ratioGeneral)/(ratioRegular - oldRatio)/swapPayFreq,
                                                minimumMaturity);
   }

  GtoFreeCFL(regularCfl);

  }
  status = SUCCESS;

  done:

  if( ammendedPrincipalDates != NULL) FREE(ammendedPrincipalDates);
  if( ammendedPrincipalPayments != NULL) FREE(ammendedPrincipalPayments);
  if( ammendedCouponDates != NULL) FREE(ammendedCouponDates);
  if( ammendedCouponPayments != NULL) FREE(ammendedCouponPayments);

  if (status ISNT SUCCESS)
  {     
    GtoErrMsg("%s: FAILED\n",routine);
    GtoFreeCFL(regularCfl);
  }
  return (status);

  }  



/* routine to compute the standard deviation of the price of a regular bond, 
  devided by the  spot volatility for the two mean reversions */ 

   
long vnfmBetaParBond (long    initialize,
                      TCashFlowList *cfl,
                      TCurve *zeroCurve,
                      TDate   resetDate,
                      long    lastCouponIdx,
                      long    frequency,
                      double  beta1,
                      double  beta2,
                      double  backboneq,
                      double *value1,
                      double *value2)
{
 static   char routine[] = "vnfmBetaParBond";
 int status = FAILURE;
 double parCoupon;
 double annuityInc;
 static double annuity;
 double discount;
 static double resetDiscount;
 long i;
 double avgRate;
 double mRMultip1;
 double mRMultip2;
 double rateToCev;
 long couponIdx;
 static double intValue1,intValue2;
 static long prevCouponIdx;
 /* This routine is build in such a way as to copute only the incremental 
    computations from the previous bond to the next */
 if (initialize == 1)
 {
  annuity = 0.;
  intValue1 = 0.;
  intValue2 = 0.;
  prevCouponIdx =-1;

  if ( GtoDiscountDate( resetDate,
                       zeroCurve,
                       /*  GTO_PARABOLIC_FORWARDS */ GTO_LINEAR_INTERP,
                      &resetDiscount) ISNT SUCCESS) goto done; 
  return SUCCESS;
 }


 /* add the contribution to variance of coupons which were not included 
    in the previous call to the routine */
 for (couponIdx=lastCouponIdx; couponIdx>prevCouponIdx; couponIdx--)
 {
  if ( GtoDiscountDate( cfl->fArray[couponIdx].fDate,
                        zeroCurve,
                         /*  GTO_PARABOLIC_FORWARDS */ GTO_LINEAR_INTERP,
                       &discount) ISNT SUCCESS) goto done;

  annuityInc = discount * cfl->fArray[couponIdx].fAmount;

  annuity += annuityInc;
  
  avgRate = log (resetDiscount/discount) 
                     /( (cfl->fArray[couponIdx].fDate - resetDate)/365.);

  if ( fabs(beta1)<1e-4 )
  {
    mRMultip1 = (cfl->fArray[couponIdx].fDate - resetDate)/365.;
  }
  else
  {
    mRMultip1 = (1-exp(-beta1*(cfl->fArray[couponIdx].fDate - resetDate)/365.))
                     /beta1;
    
  }
  if ( fabs(beta2)<1e-4 )
  {
    mRMultip2 = (cfl->fArray[couponIdx].fDate - resetDate)/365.;
  }
  else
  {
    mRMultip2 = (1-exp(-beta2*(cfl->fArray[couponIdx].fDate - resetDate)/365.))
                     /beta2;
    
  }

   rateToCev = exp(log(avgRate)*(1.-2.*backboneq));
   intValue1 += annuityInc * rateToCev * mRMultip1;
   intValue2 += annuityInc * rateToCev  * mRMultip2;
 }


 parCoupon = (resetDiscount - discount)/annuity;

 *value1 = intValue1 * parCoupon + discount * rateToCev * mRMultip1;
 *value2 = intValue2 * parCoupon + discount * rateToCev * mRMultip2;

 prevCouponIdx = lastCouponIdx;

 status = SUCCESS;

 done:

 if ( status ISNT SUCCESS)
                GtoErrMsg("%s: FAILED\n",routine);

 return status;
}

/* Routine to compute the standard deviation of the price of a generalized 
   bond devided by the  spot volatility for the two mean reversions */ 

long  vnfmBetaBondGeneral (TCurve *zeroCurve,
                       TDate  resetDate,
               	       long numCouponDates,
	                     TDate *couponDates,	
                       double *couponPayments,
                       long    numPrincipalDates,
                       TDate  *principalDates,
                       double *principalPayments,
                       double  beta1,
                       double  beta2,
                       double  backboneq,
                       double *value1,
                       double *value2)
{
 static   char routine[] = "vnfmBetaBondGeneral";
 int status = FAILURE;

 double resetZero;
 double couponZero;
 double principalZero;
 long couponIdx;
 long principalIdx;
 double avgRate;
 double mRMultip1;
 double mRMultip2;
 double daycountFraction;


   if (GtoDiscountDate(resetDate,
                       zeroCurve,
                       /*  GTO_PARABOLIC_FORWARDS */ GTO_LINEAR_INTERP,
                       &resetZero) ISNT SUCCESS) goto done;
   
   *value1 = 0.;
   *value2 = 0.; 


 /* Compute contribution to variace of the coupon payments */
   for ( couponIdx = numCouponDates-1; couponIdx >= 0 &&
               couponDates[couponIdx] > resetDate; couponIdx--)
   {

    if (GtoDiscountDate(couponDates[couponIdx],
                        zeroCurve,
                        /*  GTO_PARABOLIC_FORWARDS */ GTO_LINEAR_INTERP,
                       &couponZero) ISNT SUCCESS) goto done;

    if ( couponIdx == 0 )
    {
/*
     if (couponPayments[0] !=0)
     {
       GtoErrMsg("%s: The first coupon payment appears after reset date"
                " and is not zero\n",routine);
       goto done;
     }
*/
     daycountFraction = 1;
    }
    else
    {
     daycountFraction = MIN( (double)(couponDates[couponIdx] - resetDate) 
               /(double)( couponDates[couponIdx] - couponDates[couponIdx-1]),
                                                                        1. );
    }

    avgRate = log (resetZero/couponZero) 
                     / ((couponDates[couponIdx] - resetDate)/365.);

    if ( fabs(beta1)<1e-4 )
    {
      mRMultip1 = (couponDates[couponIdx] - resetDate)/365.;
    }
    else
    {
      mRMultip1 = (1-exp(-beta1*(couponDates[couponIdx] - resetDate)/365.))
                   /beta1;
    }

    *value1 += couponPayments[couponIdx] * daycountFraction
                      * exp(log(avgRate) * (1.-2.*backboneq)) * mRMultip1 * couponZero;
   
    if ( fabs(beta2)<1e-4 )
    {
      mRMultip2 = (couponDates[couponIdx] - resetDate)/365.;
    }
    else
    {
      mRMultip2 = (1-exp(-beta2*(couponDates[couponIdx] - resetDate)/365.))
                   /beta2;
    }

    *value2 += couponPayments[couponIdx] * daycountFraction
                       * exp(log(avgRate) * (1.-2.*backboneq)) * mRMultip2 * couponZero;
   }


 /* Compute contribution to variace of the principal payments */
   for ( principalIdx = numPrincipalDates-1; principalIdx >=0 &&
               principalDates[principalIdx] > resetDate; principalIdx--)
   {

    if (GtoDiscountDate(principalDates[principalIdx],
                        zeroCurve,
                        /*  GTO_PARABOLIC_FORWARDS */ GTO_LINEAR_INTERP,
                       &principalZero) ISNT SUCCESS) goto done;

    avgRate = log (resetZero/principalZero) 
                     / ((principalDates[principalIdx] - resetDate)/365.);

    if ( fabs(beta1)<1e-4 )
    {
      mRMultip1 = (principalDates[principalIdx] - resetDate)/365.;
    }
    else
    {
      mRMultip1=(1-exp(-beta1*(principalDates[principalIdx] - resetDate)/365.))
                       /beta1;
    }
    *value1 += principalPayments[principalIdx] 
                   * exp(log(avgRate) * (1.-2.*backboneq)) * mRMultip1 * principalZero;


    if ( fabs(beta2)<1e-4 )
    {
      mRMultip2 = (principalDates[principalIdx] - resetDate)/365.;
    }
    else
    {
      mRMultip2=(1-exp(-beta2*(principalDates[principalIdx] - resetDate)/365.))
                       /beta2;
    }
    *value2 += principalPayments[principalIdx] 
                    * exp(log(avgRate) * (1.-2.*backboneq)) * mRMultip2 * principalZero;
   }   

   status = SUCCESS;

  done:
   
  if (status !=SUCCESS)

    GtoErrMsg("%s: FAILED\n",routine);

  return (status);

}







