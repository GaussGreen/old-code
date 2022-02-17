#ifndef _vnfmeqma_H
#define	_vnfmeqma_H
#include "drlstd.h"
#include "bastypes.h"
#include "tcurve.h"

/*---------------------------------------------------------------------------
 *    VnfmMREquivalentMaturity 
 * 
 *    Given a generalised bond (or equivalently a generalized swap where the 
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
 * \begin{itemize}
 * \item[1]  Zero curve
 * \item[2]  Number of reset dates: 
 *           the number of dates on which equivalent maturities are needed.
 * \item[3]  Array of dates on which equivalent maturities are needed.
 * \item[4]  Number of coupon dates.
 * \item[5]  Array of coupon dates.
 * \item[6]  Array of coupon payments.
 * \item[7]  Number of principal Payments.
 * \item[8]  Array of principal payment dates.
 * \item[9]  Array of principal payments.
 * \item[10] First Exercise date.
 * \item[11] Strike amount.
 * \item[12] Coupons per year in the benchmark swaps. 
 * \item[13] First mean reversion.
 * \item[14] Second mean reversion.
 *           ( > than first mean reversion. This is NOT for a second factor).
 * \item[15] CevPower, the power of the interest rate that is taken before 
 *                 multiplying by the volatility to obtain the basis point 
 *                 volatility. (0 for normal 1 for log-normal between 0 and 1 
 *                 for other CEV processes).
 * \item[16] Minimum maturity.  The routine will output this maturity 
 *                          if the equivalent maturity is shorter.
 * \item[17] Maximum maturity. The routine will output this maturity if the 
 *                         equivalent maturity is larger than it. 
 *                         Note that it is important to enter here a relevant no.
 *                         since maybe no swap gives similar behaviour and no 
 *                         maturity will be found. In this case one should use 
 *                         the longest maturity available.
 * \item[18] An array of equivalent maturities is the output of this routine.
 * \end{itemize} 
 */


DLL_EXPORT(long)
VnfmMREquivalentMaturity(
	TCurve *zeroCurve,	        /* (I) */
	long    numResetDates,	    /* (I) */
	TDate  *resetDates,	        /* (I) */
	long    numCouponDates,	    /* (I) */
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
	double *maturities);	        /* (O) */





/*-----------------------------------------------------------------------
 * Wrapper for "VnfmMREquivalentMaturityL".
 * Addin function name  EQMAT 
 */

DLL_EXPORT(int)
   VnfmMREquivalentMaturityL(
	    long    *valueDateL,	        /*  1 'D' (I) */
      long    *zeroDatesL,	        /*  2 'D' (I) */
    	double  *zeroRatesL,	        /*  3 'F' (I) */
    	long    *resetDatesL,	        /*  4 'D' (I) */
	    long    *couponDatesL,	      /*  5 'D' (I) */
      double  *couponPaymentsL,     /*  6 'F' (I) */
      long    *principalDatesL,     /*  7 'D' (I) */
      double  *principalPaymentsL,  /*  8 'F' (I) */
      long    *firstExerDateL,      /*  9 'D' (I) */
      double  *exerAmountL,         /* 10 'F' (I) */
	    long    *swapPayFreqL,	      /* 11 'L' (I) */
      char    *swapDayCountConvL,   /* 12 'C' (I) */
                                     /*  ACT/365   */
                                     /*  ACT/365F  */
                                     /*  ACT/360   */
                                     /*  30/360    */
                                     /*  30E/360   */
                                     /*  ACT/365J  */
                                     /*  30E/360I  */
      double  *mr1L,                /* 13 'F' (I) */
      double  *mr2L,                /* 14 'F' (I) */
      double  *backboneqL,           /* 15 'F' (I) */
      double  *minimumMaturityL,    /* 16 'F' (I) */
      double  *maximumMaturityL,    /* 17 'F' (I) */
     	double  *maturitiesL);        /* 18 'F' (O) */



/*--------------------------------------------------------------------------
 *    VnfmDRWrapperMREquivalentMaturityL  - SPECIAL WRAPPER FOR DR WRAPPER\\
 * 
 *    Addin function name DR\_WRAPPER\_EQMAT\\
 * 
 *   This routine is a wrapper for  VnfmMREquivalentMaturity, which takes all 
 *   the specific deal inforation from the same inputs as SWAPTION_FLOWS2.
 *   (for description of the routine see VnfmMREquivalentMaturity).
 * 
 * \begin{itemize}
 * \item[1]  Value date
 * \item[2]  Zero curve dates
 * \item[3]  Zero curve rates
 * \item[4]  Array of dates on which equivalent maturities are needed.
 * \item[5]  Array of coupon dates.
 * \item[6]  Array of coupon payments.
 * \item[7]  Array of principal payment dates.
 * \item[8]  Array of principal payments.
 * \item[9]  Floating accrue start dates 
 * \item[10] Floating pay dates ( accrue end is assumed to be pay date).
 * \item[11] Floating notionals. 
 * \item[12] Index weights (the last element may contain the spread over the 
 *           floating leg). 
 * \item[13] IndexFreqs. If the number of index frequencies is less than the 
 *           number of index weights the last element of index weights is the 
 *           spread over the floating leg.
 * \item[14] UnderConvs: containes the daycount convention for the floating leg.
 * \item[15] Exercise dates: containes the first exercise date.
 * \item[16] Strike amounts: containes the first strike amount.
 * \item[17] Coupons per year in the benchmark swaps. 
 * \item[18] FDaycount convention for vanilla benchmark swaps.
 * \item[19] (1) first mean reversion.
 *           (2) second mean reversion (>first mean reversion. 
 *	         NOT for a 2nd fact).
 *           (3) backboneq, the power of the interest rate that is taken before
 *               multiplying by the volatility to obtain the basis point
 *               volatility. (0 for lognormal, 0.5 for normal, between 0 and 0.5
 *               cevPower = 1-2*q for other CEV processes).
 * \item[20] (1) Minimum maturity.  The routine will output this maturity if
 *               the equivalent maturity is shorter.
 *           (2) Maximum maturity. The routine will output this maturity if the
 *                      equivalent maturity is larger than it.
 *                       Note that it is important to enter here a relevant no.
 *                       since maybe no swap gives similar behaviour and no
 *                       maturity will be found. In this case one should use
 *                       the longest maturity available.
 * \item[21] An array of equivalent maturities is the output of this routine.
 * \end{itemize} 
 * 
 */

DLL_EXPORT(int)
   VnfmDRWrapperMREquivalentMaturityL(
	    long    *valueDateL,	         /*  1 'D' (I) */
      long    *zeroDatesL,	         /*  2 'D' (I) */
    	double  *zeroRatesL,	         /*  3 'F' (I) */
    	long    *resetDatesL,	         /*  4 'D' (I) */
	    long    *couponDatesL,	       /*  5 'D' (I) */
      double  *couponPaymentsL,      /*  6 'F' (I) */
      long    *principalDatesL,      /*  7 'D' (I) */
      double  *principalPayments,    /*  8 'F' (I) */
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
     double  *mrL,                   /* 19 'F' (I) */
     double  *minmaxMaturityL,       /* 20 'F' (I) */
     double  *maturitiesL);	     /* 21 'F' (O) */


#endif


