/*  annuityEquivalentMaturity
 *
 *  created by Arnon Levy 11/11/96
 *
 * Given a swap with non constant notional we find the maturity of a fixed 
 * notional swap, with a possible stub at the end and swapPayFreq payment
 * frequency that has the same annuity as the input swap. If this swap
 * has a larger maturity than MAX_COUPONS /swapPayFreq, 
 * MAX_COUPONS /swapPayFreq   is returned 
 * NOTE: the first coupon date is the accrual start of the swap, the notional
 *       given for that date is immaterial. notional[i] corresponds to the 
 *       payment on couponDate[i].
 *
 * RCS: $Header$
 */
#include "drlstd.h"	/* platform compatibility */
#include "ldate.h"
#include "macros.h"
#include "tcurve.h"

#include "drlaneqm.h"	/* prototype consistency */

/*f---------------------------------------------------------------------
 * Wrapper for <i> DrlAnnuityEquivalentMaturity</i>.
 */

DLL_EXPORT(int)
DrlAnnuityEquivalentMaturityL(
	long *valueDate,	/*  1 'D' (I) value date */
	long *zeroDates,	/*  2 'D' (I) zero curve dates */
	double *zeroRates,	/*  3 'F' (I) zero curve rates */
	long *resetDates,	/*  4 'D' (I) */
	long *couponDates,	/*  5 'D' (I) */
	double *notionals,	/*  6 'F' (I) */
	long *swapPayFreq,	/*  7 'L' (I) */
	double *maturities)	/*  8 'F' (O) */
{
    TCurve   *zeroCurve = NULL;   
    long numZeroDates;
    long numCoupons;

    long   status = FAILURE;
    static char routine[] = "DrlAnnuityEquivalentMaturityL";

   if (!IS_SCALAR(valueDate))
   {
       GtoErrMsg("%s: value date can't be a range.\n", routine);
       goto error;
   }
#ifdef	_SKIP
   numZeroDates = (long)zeroDates[0];
   if (numZeroDates != (long)zeroRates[0])
   {
       GtoErrMsg("%s: inconsistent number of zero rates and dates.\n",
                  routine);
       goto error;
   }
#endif

   numCoupons = (long)couponDates[0];
   if (numCoupons != (long)notionals[0])
   {
       GtoErrMsg("%s: inconsistent number of coupon notionals and dates.\n",
                  routine);
       goto error;
   }

   if (ARGSIZE(maturities) < ARGSIZE(resetDates)) 
   {
       GtoErrMsg("%s: inconsistent number of reset dates and maturities.\n",
                  routine);
       goto error;
   }



   if (!IS_SCALAR(swapPayFreq))
   {
       GtoErrMsg("%s: swapPayFreq can't be a range.\n",routine);
       goto error;
   }


   /*zeroCurve = GtoMakeTCurve(
                   valueDate[1],
                   &zeroDates[1],
                   &zeroRates[1],
                   numZeroDates,
                   1,
                   GTO_ACT_365F);*/

    if (GtoConvertTCurveWrap(
                   valueDate,
                   zeroDates,
                   zeroRates,
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

   if ( DrlAnnuityEquivalentMaturity(
                                   resetDates[0],
                                  &resetDates[1],
                                   zeroCurve, 
                                   numCoupons,
                                  &couponDates[1],
                                  &notionals[1],
                                   swapPayFreq[1], 
                                  &maturities[1]
                                          ) ISNT SUCCESS ) goto error;
   maturities[0] = resetDates[0];

   status = SUCCESS;

   error:
 
   GtoFreeTCurve(zeroCurve);

   return status;
}
