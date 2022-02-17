
/*  annuityEquivalentMaturity

    created by Arnon Levy 11/11/96

   Given a swap with non constant notional we find the maturity of a fixed 
   notional swap, with a possible stub at the end and swapPayFreq payment
   frequency that has the same annuity as the input swap. If this swap
   has a larger maturity than MAX_COUPONS /swapPayFreq, 
   MAX_COUPONS /swapPayFreq   is returned .
   NOTE: the first coupon date is the accrual start of the swap, the notional
         given for that date is immaterial. notional[i] corresponds to the 
         payment on couponDate[i].
 $Header$
*/
#include "drlstd.h"	/* platform compatibility */

#include "convert.h"
#include "macros.h"

#include "drlaneqm.h"		/* prototype consistency */


/*f---------------------------------------------------------------------
 * Given a swap with non constant notional we find the maturity of a fixed 
 * notional swap, with a possible stub at the end and swapPayFreq payment
 * frequency that has the same annuity as the input swap. If this swap
 * has a larger maturity than MAX\_COUPONS/swapPayFreq, 
 * MAX\_COUPONS /swapPayFreq  is returned.\\
 * <b> NOTE:</b>
 * the first coupon date is the accrual start of the swap, the notional
 * given for that date is immaterial. <i> notional[i]</i> corresponds to the 
 * payment on <i> couponDate[i]</i>.\\
 * <br>
 * <br>[notionals] notionals on which coupons
 * are paid on coupon dates.
 * <br>[swapPayFreq] We map varying notional
 * swaps to ones with fixed notionals and swapPayFreq payment frequency.
 * <br>[maturities] For each reset date we find the maturity of a fixed 
 * notional swap, with a possible stub at the end, and swapPayFreq payment
 * frequency that has the same annuity as the input swap.
 * <br>
 * Returns SUCCESS/FAILURE.
 */


DLL_EXPORT(int)
DrlAnnuityEquivalentMaturity(
	long numResetDates,	/* (I) number of reset dates */
	TDate *resetDates,	/* (I) reset dates [0..numResetDate-1] */
	TCurve *zeroCurve,	/* (I) zero coupon curve */
	long numCouponDates,	/* (I) number of coupon dates */
	TDate *couponDates,	/* (I) coupon dates */
	double *notionals,	/* (I) see below */
	long swapPayFreq,	/* (I) see below */
	double *maturities)	/* (O) see below */
{
 TDateInterval matInterval;
 TDate maxMat;
 long couponIdx;
 double continuingAnnuity;
 double annuity;
 long resetIdx;
 double discount;
 double swapAnnuity;
 double prevSwapAnnuity;
 TDate  prevDate;
 TDate  couponDate;
 TDateInterval  couponInterval;
 long cflIdx;
 double couponPeriod;

 static   char routine[] = "DrlAnnuityEquivalentMaturity";
 int status = FAILURE;

  couponPeriod = 1./ (double) swapPayFreq;
  couponIdx = numCouponDates-1;
  continuingAnnuity =0.;

  /*  go from last reset to first reset */

  for (resetIdx = numResetDates-1; resetIdx > -1; resetIdx--)
  {

   annuity = continuingAnnuity;

   /* go from last coupon to reset date accumulating the annuity */
   for (;  couponIdx >0 &&
               couponDates[couponIdx] > resetDates[resetIdx]; couponIdx--)
   {

    if (GtoDiscountDate(couponDates[couponIdx], 
                        zeroCurve, 
                        GTO_LINEAR_INTERP, 
                        &discount)
        ISNT SUCCESS) goto failed;

    annuity += discount * notionals[couponIdx] * 
                      (couponDates[couponIdx] - couponDates[couponIdx-1])/365.;
    
   }


   /* Since we are going to modify the annuity due to a stub period
      and relevant notional we keep the accumulated annuity for the
      computation of the annuity for the next (really previous) reset date
   */
 
   continuingAnnuity = annuity;

   /* reduce annuity if there is a stub period note: last coupon to be added
      was couponIdx+1 */
   if ( couponDates[couponIdx] < resetDates[resetIdx])
         annuity -= discount * notionals[couponIdx+1] * 
                      (resetDates[resetIdx] - couponDates[couponIdx])/365.;

   /* normalize annuity to reflect notional on reset date. NOTE that the
      annuity is not a forward annuity. We do not need the forward annuity
      to do the interpolation, since the discount factor to reset date is
      a factor that would multiply all annuities we'll calculate */


     if (IS_ALMOST_ZERO(notionals[couponIdx+1])) {
         annuity = 0e0;
     } else {
         annuity /= notionals[couponIdx+1];
     }


   /*
     Add one coupon at a time to a swap that starts on reset date, 
     untill the annuity exceeds that of the inputed swap, on the reset date
   */

   if ( annuity < 1.e-10 )
   {
    maturities[resetIdx] = 0.;
   }
   else 
   {
    prevDate = resetDates[resetIdx];
    swapAnnuity = 0.;

    for(cflIdx=0; annuity > swapAnnuity && cflIdx <= MAX_COUPONS;   cflIdx++)
    {
     prevSwapAnnuity = swapAnnuity;

     if ( GtoMakeDateInterval( (int) (couponPeriod * (cflIdx+1) * 12 +1.e-6) 
                               ,'M', &couponInterval)
                           ISNT SUCCESS) goto failed;

     if ( GtoDtFwdAny(resetDates[resetIdx],&couponInterval,&couponDate)
                           ISNT SUCCESS) goto failed;


     if (GtoDiscountDate(couponDate, 
                        zeroCurve, 
                        GTO_LINEAR_INTERP, 
                        &discount)
        ISNT SUCCESS) goto failed; 



     swapAnnuity += (couponDate - prevDate)/365. * discount;

     prevDate = couponDate;
    
    }

    if ( cflIdx > MAX_COUPONS)
    {
     maturities[resetIdx] = (double) MAX_COUPONS * couponPeriod; 
     
     GtoErrMsg("\n%s: WARNING: need longer maturity than provided by ",
                   routine);
     GtoErrMsg("MAX_COUPONS.\n resetIdx = %d annuity = %16.4f\n", 
                                               resetIdx, annuity);
    }  
    else
    {
    /* equivalent maturity computation */

    maturities[resetIdx] = ( (double)(cflIdx - 1) 
             + (annuity - prevSwapAnnuity)/ ( swapAnnuity - prevSwapAnnuity))
                                                * couponPeriod;
    }
   }
  }

  status = SUCCESS;

  failed:

  return status;
}
  




















