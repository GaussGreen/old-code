/****************************************************************
 * Module:	VNFM
 * Submodule:	WRAP
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <math.h>

#include "convert.h"
#include "macros.h"
#include "ldate.h"
#include "yearfrac.h"
#include "date_sup.h"

#include "drlmem.h"
#include "drlvtype.h"			/* DrlLilStructGet() */
#include "drlsmat.h"
#include "drlts.h"			/* DrlTCurveWrap() */
#include "drlio.h"			/* DrlFPrintf() */

#include "swaptn2l.h"                   /* Prototype consistency */

#define	_vnfm_SOURCE
#include "vnfmanly.h"
#include "vnfmcali.h"
#include "vnfmwrap.h"

#if defined(_WINDLL) 
# undef __DEBUG__
#endif



/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CAL2F1SPVOL_GEN.
 *                                                             
 * <br><br>
 * Wrapper for routine <i> VnfmCalib1V2FGeneralL</i>
 * that simultaneously performs spot volatility calibration 
 * of arbitrary volatity points
 * and generate implied base volatility curve and model swaption matrix.
 * <br>
 * <br>[refDateL] zero curve value date.
 * <br>[zcDateL] zero curve coupon dates.
 * <br>[zcRateL] zero curve coupon rates.
 * <br>[floatScalarsL] numerical  scalars:\\
 *	(1) back bone $q$ (0 for log-normal, 0.5 for normal),\\
 *	(2) generate base vol (0 or 1),\\
 *	(3) generate swaption matrix(0 or 1), \\
 *	(4) minimum volatility rate maturity($&gt;= 0$).
 * <br>[nfParamsL] LIL array containg the mean reversion, weight and
 * correlation coefficients.
 * Tree formats are available: an array format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n}),$$
 * that has ${n(n+3)/ 2}$ elements,
 * or a ``short'' array format where the 1st weight,
 * assumed to be 1.0, is omitted 
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_2,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n})$$
 * that has ${n(n+3)/ 2}-1$ elements,
 * or a matrix format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,1},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n,n})$$
 * that has $n(n+2)$ elements,
 * {\bf Remark: because of possible ambiguity in the format
 * (e.g. 8 elements can correspond to $n=3$ of the short array format
 * or $n=2$ of the matrix form) the format is checked in 
 * the previous order).}
 * <br>[nfDatesL] array of dates used for the spot volatility
 * bootstrapping (the spot volatility is assumed to be constant between
 * two consecutive dates).\\
 * {\bf WARNING: Unlike the ALIB routines, the convention is that
 * the spot volatility $i$ applies between date $i$ and date $i+1$.
 * The first date in the array MUST be the today date}.
 * <br>[rateMatL] array of volatility rate maturity (in yrs).
 * The element of index $i$ in the array corresponds
 * to an option expiring at date $i$ of array <i> nfDatesL</i>.
 * This array must have same length as the array <i> nfDatesL</i>.
 * {\bf
 * <br>
 *  <br> The volatility 1 (corresponding to the date 1, i.e. today)
 *     is never calibrated.
 * <br> The calibration stops at the first negative maturity encountered
 *     according to the following rules:
 *     <br>
 *     <br> If the maturity 1 (corresponding to the date 1, i.e. today)
 *         is negative, an error is returned.
 *     <br> If the maturity 1 is positive,
 *         but maturity 2 is negative, then volatility 2 is
 *         nevertheless calibrated assuming a maturity 0.25.
 *         No other volatility point is calibrated.
 *     <br> Otherwise, the volatility points are calibrated until a negative
 *         maturity is encountered.
 * <br>
 * <br> Any maturity smaller than floatScalarsL[4], but positive,
 *     will be rounded to floatScalarsL[4] (and the corresponding volatility 
 *     will be calibrated as such).
 * <br>
 * }
 * <br>[rateFreqL] array of volatility frequency (0,1,2,4,12).
 * {\it The element of index $i$ in the array corresponds
 * to an option expiring at date $i$ of array <i> nfDatesL</i>}.
 * <br>[rateVolL] array of volatilities.
 * {\it The element of index $i$ in the array corresponds
 * to an option expiring at date $i$ of array <i> nfDatesL</i>}.
 * <br>[bvMatL] base volatility underlying rate  maturity (in yrs)
 * for the output <i> bvRatesL</i>.
 * <br>[bvDatesL] base volatility expiration dates
 * for the output <i> bvRatesL</i>.
 * <br>[swTypeL] array of 2 elements:\\
 * 	(1) swaption matrix output type (0 for vertical, 1 for diagonal)
 *		for the output <i> SwVolL</i>.\\
 * 	(2) swaption matrix volatility frequency (0,1,2,4,12)
 *		for the output <i> SwVolL</i>.
 * <br>[swMatL] array of maturity intervals (in yrs)
 *		for the output <i> SwVolL</i>.
 * <br>[swExpL] array of expration intervals (in yrs)
 *		for the output <i> SwVolL</i>.
 * <br>[spotVolsL] output spot volatility corresponding to
 * array of dates <i> nfDatesL</i>.
 * {\bf WARNING: Unlike the ALIB routines, the convention is that
 * the spot volatility $i$ applies between date $i$ and date $i+1$.
 * Therefore, the volatility dates sould be offset by one when
 * passing the curve to an ALIB tree routine}.
 * <br>[bvRatesL] output base volatility curve
 * specified by arguments <i> bvMatL</i> and <i> bvDatesL</i>
 * <br>[swVolL] output swaption matrix
 * specified by arguments <i> swTypeL</i>, <i> swMatL</i> and <i> swExpL</i>.
 * <br>
 */

DLL_EXPORT(int)
VnfmEqmatL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */
	
	FloatL *floatScalarsL,	/*  4 'F' (I) num scalars:
				 *        [1] back bone q */
	FloatL *nfParamsL,	/*  5 'F' (I) N-fact params */
	TDateL *nfDatesL,	/*  6 'D' (I) array of dates */

	FloatL *cpnL,

	FloatL *teqMatL)	/* 17 'F' (O) swaption vol matrix */
{
static	char		routine[] = "VnfmEqmatL";
	int		status = FAILURE;

	VnfmData	*vnfmData = NULL;
	TCurve		*zcCurve = NULL;
	int		nDates;

	int		idxStart = 0,
			idxEnd;
	double		tExp,
			tMatMin,
			tMatMax;
	int		freq = 2;
	double		B[VNFM_NDIMMAX],
			B1[VNFM_NDIMMAX];	/* bullet V-coefficients */
	double		tEqMat;

	double		tMat = 10e0;


	VNFM_LOGOPEN

	/* check arg len */
	WRAP_CHECK_SCALAR(refDateL);
	WRAP_CHECK_VECTOR(zcDateL);
	WRAP_CHECK_VECTOR(zcRateL);
	WRAP_CHECK_VECTOR_LEN(floatScalarsL, 1);
	WRAP_CHECK_VECTOR(nfParamsL);
	WRAP_CHECK_VECTOR(nfDatesL);

	WRAP_CHECK_VECTOR_LEN(cpnL, 1);


	WRAP_CHECK_VECTOR(teqMatL);


	/* get zero curve */
	IF_FAILED_DONE( DrlTCurveWrapRead(
		&zcCurve,
		refDateL,
		zcDateL,
		zcRateL));

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve,\
		vnfmFpLog, DRL_TCURVE_FMT_PERCENT));
#endif

	nDates = ARGSIZE(nfDatesL);

	/* Create Vnfm structure */
	IF_FAILED_DONE( VnfmWrapReadSimple(
		&vnfmData,
		floatScalarsL,
		nfParamsL,
		nfDatesL,
		zcCurve));

	IF_FAILED_DONE( VnfmComputeCoeff(
		vnfmData));



#ifndef	NO_LOGGING
	GTO_IF_LOGGING( VnfmFpWrite(vnfmData,\
		vnfmFpLog));
#endif


	/* */
	tExp = 0e0;
	tMat = 10e0;
	tMatMin = 0.25e0;
	tMatMax = 30e0;
	freq = 2;


	IF_FAILED_DONE( VnfmBulletB1(
		vnfmData,
		tExp,
		cpnL[1],
		tMat,
		freq,
		B,
		B1));


	IF_FAILED_DONE( VnfmSolveBulletEqmat(
		vnfmData,
		tExp,
		B,
		B1,
		tMatMin,
		tMatMax,
		freq,
		&tEqMat));


#ifndef	NO_LOGGING
	GTO_IF_LOGGING( DrlFPrintf(vnfmFpLog, \
		"%s: cpn=%8.4f%% tMat=%12.8f tEqMat=%12.8f.\n",
			routine, cpnL[1]*1e2, tMat, tEqMat));
#endif






	/* made it through OK */
	status = SUCCESS;
done:
	VnfmFree(vnfmData);
	GtoFreeTCurve(zcCurve);

	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}
	VNFM_LOGCLOSE
	return(status);

}





/*f--------------------------------------------------------------
 * Wrapper for SWAPTION_FLOWS.
 *                                                             
 * <br><br>
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
 *            multiplying by the volatility to obtain the basis point
 *            volatility. (0 for lognormal, 0.5 for normal, between 0 and 0.5
 *            cevPower = 1-2*q for other CEV processes).
 * <br>[20] (1) Minimum maturity.  The routine will output this maturity if 
 * 	      the equivalent maturity is shorter.
 * 	  (2) Maximum maturity. The routine will output this maturity if the 
 *            equivalent maturity is larger than it. 
 *            Note that it is important to enter here a relevant no.
 *            since maybe no swap gives similar behaviour and no 
 *            maturity will be found. In this case one should use 
 *            the longest maturity available.
 * <br>[21] An array of equivalent maturities is the output of this routine.
 * <br> 
 * 
*/

DLL_EXPORT(int)
VnfmSW2EqmatL(
	long    *valueDateL,	       /*  1 'D' (I) */
	long    *zeroDatesL,	       /*  2 'D' (I) */
	double  *zeroRatesL,	       /*  3 'F' (I) */
	long    *resetDatesL,	       /*  4 'D' (I) */
	long    *coupDatesL,	       /*  5 'D' (I) */
	double  *coupPaymentsL,		/*  6 'F' (I) */
	long    *prinDatesL,		/*  7 'D' (I) */
	double  *prinPaymentsL,		/*  8 'F' (I) */
	long    *floatAccDatesL,	/*  9 'D' (I) */
	long    *floatPayDatesL,       /* 10 'D' (I) */
	double  *floatNotionalsL,      /* 11 'F' (I) */
	double  *indexWeightsL,        /* 12 'F' (I) last ele may be a spd */
	long    *indexFreqsL,          /* 13 'L' (I) determine if a spd given */
	char    *underConvsL,          /* 14 'C' (I) */
	long    *exerDatesL,           /* 15 'D' (I) */
	double  *exerAmountsL,         /* 16 'F' (I) */
	long    *swapPayFreqL,	       /* 17 'L' (I) */
	char    *swapDccL,		/* 18 'C' (I) */
	double  *mrL,                  /* 19 'F' (I) */
	double  *minmaxMaturityL,      /* 20 'F' (I) */
	double  *maturitiesL)	       /* 21 'F' (O) */
{
static	char	routine[] = "VnfmSW2EqmatL";
	int	status = FAILURE;


	TCurve   *zcCurve = NULL;
	long	numZeroDates;
	long	swapDcc;
	double	spread;
	long	floatPayDcc;
	long	cpIdx;
	long	flIdx;
	double	notional;
	double	sum;
	double	yearFrac;
	double	discount;
	double	*bondPrincipalPayments = NULL;
	long	prinIdx;

	double	floatPrincipalPay;
	double	strikeAdjustment;

	long	numCoup;
	TDate	*coupDates;
	double	*coupPayments;

	long	numPrin;
	TDate	*prinDates;
	double	*prinPayments;

	long	numFloatPayments;
	TDate	*floatAccDates;
	TDate	*floatPayDates;
	double	*floatNotionals;

	double mr1, mr2, backboneq, minMat, maxMat;


	int	newMethod = FALSE;

	/*
	 * Log inputs
	 */
	if (GtoLoggingGet() > 0) {
		DrlLilVectLoggingFile("wrapper.log", "w", "DR_WRAPPER_EQMAT",
		    DRL_TDATE_L,  valueDateL,      "ZC_VALUEDATE",
		    DRL_TDATE_L,  zeroDatesL,      "ZC_DATES",
		    DRL_FLOAT_L,  zeroRatesL,      "ZC_RATES",

		    DRL_TDATE_L,  resetDatesL,     "RESET_DATES",
		    DRL_TDATE_L,  coupDatesL,      "COUPON_DATES",
		    DRL_FLOAT_L,  coupPaymentsL,   "COUPON_PAYMENTS",
		    DRL_TDATE_L,  prinDatesL,      "PRINCIPAL_DATES",
		    DRL_FLOAT_L,  prinPaymentsL,   "PRINCIPAL_PAYMENTS",
		    DRL_TDATE_L,  floatAccDatesL,  "FLOATACC_DATES",
		    DRL_TDATE_L,  floatPayDatesL,  "FLOATPAY_DATES",
		    DRL_FLOAT_L,  floatNotionalsL, "FLOAT_NOTIONALS",
		    DRL_FLOAT_L,  indexWeightsL,   "INDEX_WEIGHTS",

		    DRL_LONG_L,   indexFreqsL,     "INDEX_FREQS",

		    DRL_CHAR_BLOCK_L, underConvsL, "UNDER_CONVS",
		    DRL_TDATE_L,  exerDatesL,      "EXER_DATES",
		    DRL_FLOAT_L,  exerAmountsL,    "EXER_AMNT",
		    DRL_LONG_L,   swapPayFreqL,    "SWAP_FREQ",
		    DRL_CHAR_BLOCK_L, swapDccL,    "SWAP_DCC",
		    DRL_FLOAT_L,  mrL,             "MOD_PAR",

		    DRL_FLOAT_L,  minmaxMaturityL, "MINMAX_MATS",
		    0);
	}


	/* Check length of arrays coming from Lotus
	     */
	if (NOT (
	    IS_VECTOR(floatAccDatesL)   &&
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
		goto done;
	}

	numCoup = (long)coupDatesL[0];
	if (numCoup != (long)coupPaymentsL[0])
	{
		GtoErrMsg("%s: inconsistent number of coup dates "
			"and payments.\n",routine);
		goto done;
	}

	numPrin = (long)prinDatesL[0];
	if (numPrin != (long)prinPaymentsL[0])
	{
		GtoErrMsg("%s: inconsistent number of prin dates "
			"and payments.\n",routine);
		goto done;
	}

	if (ARGSIZE(maturitiesL) < ARGSIZE(resetDatesL))
	{
		GtoErrMsg("%s: inconsistent number of reset dates "
			"and maturities.\n", routine);
		goto done;
	}

	if ( exerDatesL[0] < 1 )
	{
		GtoErrMsg("%s: no exercise dates were given\n",
			routine);
		goto done;
	}

	if ( exerAmountsL[0] < 1 )
	{
		GtoErrMsg("%s: no exercise amounts were given\n",
		    routine);
		goto done;
	}

	if (!IS_SCALAR(swapPayFreqL))
	{
		GtoErrMsg("%s: swapPayFreq can't be a range.\n",routine);
		goto done;
	}

	if (!IS_SCALAR(swapDccL))
	{
		GtoErrMsg("%s: swapDaycountConv can't be a range.\n",routine);
		goto done;
	}

	/* Check mean reversion and bb input */
	if (mrL[0] != 3)
	{
		GtoErrMsg("%s: mr must contains mr1, mr2 and "
			"backbone q.\n",routine);
		goto done;
	}

	mr1 = mrL[1];
	mr2 = mrL[2];
	backboneq = mrL[3];

	if (mr1 < 0e0)
	{
		GtoErrMsg("%s: mr1=%12.6f  can't be a negative.\n",
			routine,mrL[1]);
		goto done;
	}

	if (mr2 < 0e0)
	{
		GtoErrMsg("%s: mr2=%12.6f  can't be a negative.\n",
			routine,mrL[2]);
		goto done;
	}

	if (mr2 <= mr1)
	{
		GtoErrMsg("%s: mr2=%12.6f  has to be larger "
			"than mr1=%12.6f\n",routine, mr2, mr1);
		goto done;
	}



	if (minmaxMaturityL[0] != 2)
	{
		GtoErrMsg("%s: minmaxMaturity must have two arguments,"
			"minimumMaturity and maximumMaturity.\n",routine);
		goto done;
	}

	minMat = minmaxMaturityL[1];
	maxMat = minmaxMaturityL[2];

	/* Convert curve */
	IF_FAILED_DONE( GtoConvertTCurveWrap(
	    valueDateL,
	    zeroDatesL,
	    zeroRatesL,
	    1L,
	    GTO_ACT_365F,
	    routine,
	    &zcCurve));

	ASSERT_OR_DONE(zcCurve != NULL);


	IF_FAILED_DONE( GtoStringToDayCountConv(
		&swapDccL[1],
		&swapDcc));




	numFloatPayments = floatPayDatesL[0];

	if (numFloatPayments != floatNotionalsL[0] ||
	    numFloatPayments != floatAccDatesL[0])
	{
		GtoErrMsg("%s: number of floating payments accrual starts and"
		    " notionals is not consistent\n",routine);
		return(FAILURE);
	}

	coupDates       = coupDatesL+1;
	coupPayments    = coupPaymentsL+1;
	prinDates       = prinDatesL+1;
	prinPayments    = prinPaymentsL+1;
	floatAccDates   = floatAccDatesL+1;
	floatPayDates   = floatPayDatesL+1;
	floatNotionals  = floatNotionalsL+1;



	/****************************************************
	 ** Convert floating leg to fixed payments
	 ****************************************************/

	if ( numFloatPayments > 0 )
	{
	    /* Floating leg is modeled. Since eqmat assumes that it is given a 
	     * bond the prin payments have to be adjusted and the coups 
	     * should be adjusted if there is a spread.
	     */

	    if ( indexFreqsL[0]+1 == ((long) indexWeightsL[0]))
	    {
		/* Get spread and floating dcc
		 */
		spread = indexWeightsL[((long)indexWeightsL[0])];
		WRAP_CHECK_VECTOR_LEN(underConvsL, GTO_SW2_NUM_UNDER_CONVS);
		/*if ((long)underConvsL[0] ISNT GTO_SW2_NUM_UNDER_CONVS)
		{
			GtoErrMsg("%s: Received %ld underlying conventions, "
				"instead of %d.\n", routine,
				(long)underConvsL[0],
				GTO_SW2_NUM_UNDER_CONVS);
			goto done;
		}*/
		IF_FAILED_DONE( GtoStringToDayCountConv(
			&underConvsL[WRAP_STR_IDX(GTO_SW2_FL_PAY_DCC_IDX)],
			&floatPayDcc));


		/* Skip coupon payments before floating
		 */
		for (cpIdx = 0;
		    cpIdx < numCoup && coupDates[cpIdx] < floatPayDates[0];
		    cpIdx++);

		if (cpIdx >= numCoup)
		{
			GtoErrMsg("%s: first floating payment is past "
				"last fixed payment\n", routine);
			goto done;
		}

		for (flIdx = 0; 
		     cpIdx < numCoup && flIdx < numFloatPayments;
		     cpIdx++)
		{
		    notional = floatNotionals[flIdx];
		    sum = 0.;
		    while ( flIdx < numFloatPayments &&
			    floatPayDates[flIdx] <= coupDates[cpIdx] )
		    {
			IF_FAILED_DONE( GtoDayCountFraction(
				floatAccDates[flIdx],
				floatPayDates[flIdx],
				floatPayDcc,
				&yearFrac));

			IF_FAILED_DONE( GtoDiscountDate(
				floatPayDates[flIdx],
				zcCurve,
				GTO_LINEAR_INTERP,
				&discount));


			sum += yearFrac * discount;

			if ( floatNotionals[flIdx] != notional )
			{
				GtoErrMsg("%s: notional cannot change between"
					" fixed coupon dates \n", routine);
				goto done;
			}

			flIdx++;
		    }

		    if (floatPayDates[flIdx-1] != coupDates[cpIdx] &&
			notional > 0)
		    {
			GtoErrMsg("%s: coup date no. %d is not a "
				"floating date\n", routine,cpIdx);
			goto done;
		    }

		    coupPayments[cpIdx] -= notional * spread * (sum / discount);

		}
	    }

	    bondPrincipalPayments = NEW_ARRAY(double, numCoup);
	    if ( bondPrincipalPayments IS NULL) goto done;

	    cpIdx = 0;
	    for (prinIdx = 0; prinIdx < numPrin; prinIdx++)
	    {
		if (prinPayments[prinIdx] != 0.)
		{
			while (cpIdx < numCoup 
			    && coupDates[cpIdx] < prinDates[prinIdx] )
				cpIdx++;
			if (  coupDates[cpIdx] != prinDates[prinIdx])
			{
				GtoErrMsg("%s: prin dates must be fixed coup dates \n",      
				    routine);
				goto done;
			}

			bondPrincipalPayments[cpIdx] += prinPayments[prinIdx];
		}
	    }

	    if ( coupDates[0] > floatAccDates[0])
	    {
		GtoErrMsg("%s: first floating acrrue start must be on or after"
		    "first coup date (this first coup date is actually"
		    " the accrual start for the fixed leg).",      
		    routine);
		goto done;
	    }

	    for (cpIdx = 0; cpIdx < numCoup 
	    	&& coupDates[cpIdx] < floatAccDates[0];
		    cpIdx++);
	    if ( exerDatesL[1] <  coupDates[cpIdx])
	    {
		bondPrincipalPayments[cpIdx] -= floatNotionals[0];
		strikeAdjustment = 0;
	    } else {
		for (flIdx = 0; flIdx < numFloatPayments &&
		    floatPayDates[flIdx] < exerDatesL[1] ; flIdx++);

		strikeAdjustment = floatNotionals[flIdx];
	    }


	    cpIdx = 0;
	    for (flIdx = 0; flIdx < numFloatPayments; flIdx++)
	    {
		if (flIdx == numFloatPayments-1)
		{
			floatPrincipalPay =  floatNotionals[flIdx];
		} else {
			floatPrincipalPay = floatNotionals[flIdx]- floatNotionals[flIdx+1];
		}
		if ( floatPrincipalPay !=0 )
		{
			while (cpIdx < numCoup 
			    && coupDates[cpIdx] < floatPayDates[flIdx] )
				cpIdx++;
			if (  coupDates[cpIdx] != floatPayDates[flIdx])
			{
				GtoErrMsg("%s: float amortization dates must be fixed coup dates \n",      
				    routine);
				goto done;
			}

			bondPrincipalPayments[cpIdx] += floatPrincipalPay;
		}
	    }
	    prinPayments = bondPrincipalPayments;
	}


	/*
	 * Call Analytivs
	 */
	if ( VnfmMREquivalentMaturity(
		zcCurve,
		(long) resetDatesL[0],
		&resetDatesL[1],
		numCoup,
		coupDates,
		coupPayments,
		numCoup,
		coupDates,
		prinPayments,
		exerDatesL[1],
		exerAmountsL[1]+strikeAdjustment,
		swapPayFreqL[1],
		swapDcc,
		mr1,
		mr2,
		backboneq,
		minMat,
		maxMat,
		&maturitiesL[1]
		) != SUCCESS ) goto done;


	maturitiesL[0] = resetDatesL[0];


#define	_NEW
#ifdef	_NEW
	{ 
	double		B[VNFM_NDIMMAX],
			B1[VNFM_NDIMMAX];	/* bullet V-coefficients */
	VnfmData	*vnfmData = NULL;
	double		backboneqL[2];	/* back bone Q */
	double		paramsL[2];	/* array of mr, weights, corrs */
	int		swapFreq = swapPayFreqL[1];
	double		tExp, tMat[1024];

	double		tEqmat[1024];

	int		numResetDates = ARGSIZE(resetDatesL),
			idxR;
	TDate		*resetDates = &resetDatesL[1];
	TDate		vnfmDatesL[3];


	/* Create Vnfm structure */
	backboneqL[0] = 1e0;
	backboneqL[1] = backboneq;
	paramsL[0] = 1e0;
	paramsL[1] = mr1;
	/* Set up a flat vol vnfm for coeff calculation */
	vnfmDatesL[0] = 2L;
	vnfmDatesL[1] = zcCurve->fBaseDate;
	vnfmDatesL[2] = zcCurve->fBaseDate + 10000L;
	IF_FAILED_DONE( VnfmWrapReadSimple(
		&vnfmData,
		backboneqL,
		paramsL,
		vnfmDatesL,
		zcCurve));

	IF_FAILED_DONE( VnfmComputeCoeff(
		vnfmData));


	for (idxR=0; idxR<numResetDates; idxR++) {

		/* COmpute V-coefficients of adjusted cash-flows
		 */
		IF_FAILED_DONE( VnfmFlowsB1(
			vnfmData,
			resetDates[idxR],
			numCoup,
			coupDates,
			coupPayments,
			numCoup,
			coupDates,
			prinPayments,

			(int) ARGSIZE(exerDatesL),
			&exerDatesL[1],
			&exerAmountsL[1],

			TRUE,	/* FALSE=euro */
			B,
			B1));

		IF_FAILED_DONE( GtoDayCountFraction(
			resetDates[idxR],
			coupDates[numCoup-1],
			GTO_ACT_365F,
			&tMat[idxR]));

		/** Check any coupons */
		if (IS_ALMOST_ZERO(B[0])) {
			tEqmat[idxR] = minMat;
			continue;
		}




		/* COmpute time to expiration */
		IF_FAILED_DONE( GtoDayCountFraction(
			vnfmData->fDate[0],
			resetDates[idxR],
			GTO_ACT_365F,
			&tExp));



		IF_FAILED_DONE( VnfmSolveBulletEqmat(
			vnfmData,
			tExp,
			B,
			B1,
			minMat,
			maxMat,
			swapFreq,
			&tEqmat[idxR]));
	}

	GtoErrMsg("%s: \t----------  ---TMAT-----  "
			"---OLD------  ---NEW------\n",
			routine);
	for (idxR=0; idxR<numResetDates; idxR++) {
		GtoErrMsg("%s: \t%10s  %12.8f  %12.8f  %12.8f\n",
			routine,
			DrlTDatePrint(NULL, resetDates[idxR]),
			tMat[idxR],
			maturitiesL[idxR+1],
			tEqmat[idxR]);
	}

	VnfmFree(vnfmData);

	}
#endif






















	/* OK */
	status = SUCCESS;
done:
	GtoFreeTCurve(zcCurve);
	FREE(bondPrincipalPayments);

	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}
	return (status);
}




