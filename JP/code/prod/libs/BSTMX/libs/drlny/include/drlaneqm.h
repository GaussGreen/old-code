/*  annuityEquivalentMaturity

    created by Arnon Levy 11/11/96

   Given a swap with non constant notional we find the maturity of a fixed 
   notional swap, with a possible stub at the end and swapPayFreq payment
   frequency that has the same annuity as the input swap. If this swap
   has a larger maturity than MAX_COUPONS /swapPayFreq, 
   MAX_COUPONS /swapPayFreq   is returned 
   NOTE: the first coupon date is the accrual start of the swap, the notional
         given for that date is immaterial. notional[i] corresponds to the 
         payment on couponDate[i].

   RCS: $Header$
*/
#ifndef _drlaneqm_H
#define	_drlaneqm_H
#include "drlstd.h"	/* DLL_EXPORT() */


#include "bastypes.h"
#include "dateconv.h"
#include "tcurve.h"

#define MAX_COUPONS 1000

extern	DLL_EXPORT(int)
DrlAnnuityEquivalentMaturity(
	long numResetDates,	/* (I) */
	TDate *resetDates,	/* (I) */
	TCurve *zeroCurve,	/* (I) */
	long numCouponDates,	/* (I) */
	TDate *couponDates,	/* (I) */
	double *notionals,	/* (I) notionals on which coupons
				       are paid on coupon dates */
	long swapPayFreq,	/* (I) We map varying notional
				       swaps to ones with fixed
				       notionals and swapPayFreq
				       payment frequency */
	double *maturities);	/* (O) For each reset date we find
					the maturity of a fixed 
					notional swap, with a 
					possible stub at the end,
					and swapPayFreq payment
					frequency that has the same
					annuity as the input swap */


/* wrapper for annuityEquivalentMaturity */

extern	DLL_EXPORT(int)
DrlAnnuityEquivalentMaturityL(
	long *valueDate,		/*  1 'D' (I) */
	long *zeroDates,		/*  2 'D' (I) */
	double *zeroRates,		/*  3 'F' (I) */
	long *resetDates,		/*  4 'D' (I) */
	long *couponDates,		/*  5 'D' (I) */
	double *notionals,		/*  6 'F' (I) */
	long *swapPayFreq,		/*  7 'L' (I) */
	double *maturities);		/*  8 'F' (O) */


#endif /*_drlaneqm_H*/







