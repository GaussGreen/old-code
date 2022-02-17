/************************************************************************
 * Module:      Brute Force Mean-reverting Volatility Calibration
 * File:
 * Function:    
 * Author:      C. Daher
 * Revision:    $Header$
 ************************************************************************/
#include "drlstd.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "cgeneral.h"
#include "bastypes.h"
#include "macros.h"             /* MAX */
#include "ldate.h"              /* GtoDayCountFraction */
#include "stub.h"               /* GTO_STUB_NONE */
#include "datelist.h"           /* GtoNewDateList */
#include "tcurve.h"             /* GtoDiscountDate */
#include "cerror.h"             /* GtoErrMsg */
#include "date_sup.h"
#include "convert.h"
#include "yearfrac.h"
#include "zr2coup.h"
#include "zr2simp.h"
#include "duration.h"		/* GtoBondModDuration */
#include "convex.h"		/* GtoBondConvexity */

#include "cashflow.h"
#include "barrier.h"

#include "drlsmat.h"
#include "drltime.h"
#include "drlio.h"
#include "drloptio.h"
#include "drlstr.h" 		/* DrlStrSubsEnv */
#include "dritkwrp.h"

#include "drirlopt.h"		/* Prototype consistency */

#if defined(VOLCONV_LN)
#include "voladj.h"		/* GtoSwaptionVolNumToBlack  */
#endif

#define	__DEBUG__
#undef	__DEBUG__

#define	ISSTR(s1, s2)	(!strcmp((s1),(s2)))

#undef	SQR
#define	SQR(x)	((x)*(x))

#define	NO_ADJUSTMENTS



/*----------------------------------------------------------------------
 * Convenience routine to compute par sensitivities:
 * duration, convexity and cubicity (-d^3B / dY^3).
 */

static	int
bndParSens(
	double mat,			/* (I) maturity in years */
	int freq,			/* (I) bond frequency */
	double ytm,			/* (I) yield to maturity */
	double *dur,			/* (I) par IRR duration */
	double *cvx,			/* (I) par IRR convexity */
	double *cub)			/* (I) par IRR cubicity */
{
static	char	routine[] = "bndParSens";
	int	status = FAILURE;
	double	dfreq = (double) freq;


	*dur = (1e0 - pow(1 + ytm / dfreq, -mat*dfreq)) / ytm;

	*cvx = (2e0 / ytm)
		* (*dur - mat * pow(1 + ytm / dfreq, -mat*dfreq - 1e0 ));

	*cub = (3e0 / ytm)
		* (-0.5*(*cvx) - mat * (mat + 1e0/dfreq) * 
			pow(1 + ytm / dfreq, -mat*dfreq - 2e0 ));


	/* OK */
	status = SUCCESS;
/*done:*/
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}



/*----------------------------------------------------------------------
 * Convenience routine to compute a forward rate.
 */

static	int
computeFwdRate(
	TCurve *zcCurve,	/* (I) zero curve */
	TDateInterval rateMat,	/* (I) rate maturity interval */
	int rateFreq,		/* (I) rate frequency */
	TDayCount rateDcc,	/* (I) rate day count convention */
	TDate rateEffDate,	/* (I) rate effective date */
	double *fwdRate)	/* (O) forward rate */
{
static	char	routine[] = "computeFwdRate";
	TDateInterval	payInterval;
	TDate		maturityDate;
	int		status = 0;


	if (GtoDtFwdAny(rateEffDate, &rateMat, &maturityDate)
		!= SUCCESS) goto done;

	switch (rateFreq) {
	case 1:
	case 2:
	case 4:
	case 12:
		if (GtoFreq2TDateInterval((long) rateFreq, &payInterval)
			!= SUCCESS) goto done;

		if (GtoZerosToCouponsPoint(zcCurve,
			GTO_LINEAR_INTERP,
			rateEffDate,
			&payInterval,
			maturityDate,
			rateDcc,
			GTO_STUB_BOND,
			FALSE,
			fwdRate) != SUCCESS) goto done;
		break;

	case 0:
		if (GtoZerosToSimplePoint(zcCurve,
			GTO_LINEAR_INTERP,
			rateEffDate,
			maturityDate,
			rateDcc,
			fwdRate) != SUCCESS) goto done;
		break;
	default:
		GtoErrMsg("%s: bad frequency %d.\n", routine, rateFreq);
		goto done;
	}

	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}

/*----------------------------------------------------------------------
 * Convenience routine to compute a forward rate, its ATM and smile volatility.
 * possibly adjusted for smile and compounding frequency.
 */

static	int
computeFwdRateVol(
	TDate todayDate,		/* (I) today's date */
	TCurve *zcCurve,		/* (I) zero curve */
	TSwaptionMatrix2D *swMat,	/* (I) swaption matrix */
	TSmile3Data *smlData,		/* (I) smile data (or NULL) */

	TDateInterval rateMat,		/* (I) rate maturity interval */
	int rateFreq,			/* (I) rate frequency */
	TDayCount rateDcc,		/* (I) rate day count convention */
	TDate rateEffDate,		/* (I) rate start date */

	TDate rateObsDate,		/* (I) rate observation date */
	double strikeRate,		/* (I) strike rate for options */

	double *rateFwd,		/* (O) rate forward */
	double *rateVolAtm,		/* (O) rate ATM volatility */
	double *rateVolSmile)		/* (O) rate smile volatility */
{
static	char	routine[] = "computeFwdRate";
	int		status = 0;
	double		rateMatYrs,
			rateVolSmileAdj;
	int		volFreq = swMat->swapPayFreq;


	IF_FAILED_DONE( GtoDateIntervalToYears(&rateMat, &rateMatYrs));


	/* Compute forward rate */

	IF_FAILED_DONE( computeFwdRate(
		zcCurve,
		rateMat,
		rateFreq,
		rateDcc,
		rateEffDate,
		rateFwd));



	/* Interpolate vol curve from swaption matrix */

	IF_FAILED_DONE(DrlTSwaptionMatrix2DDateInterpValue(
		swMat,		/* swaption matrix */
		todayDate,
		rateObsDate,
		FALSE, 		/* TRUE=final, FALSE=cms */
		(TDate)0L,	/* Not used */
		rateMat,	/* maturity interval */
		rateVolAtm,
		FALSE)); 	/* TRUE=rebucket, FALSE=interp */


	/* Perform frequency adjustment for volatility
	 * if the swaption matrix frequency is not equal
	 * to the rate.
	 */
/*
	!!! TO BE TESTED !!!!

	if (rateFreq != volFreq) {
		double	r1 = *rateFwd,
			f1 = (double) rateFreq,
			v1,
			r2,
			f2 = (double) swMat->swapPayFreq,
			v2 = *rateVolAtm;

		r2 = f2 * (pow(1e0 + r1/f1,f1/f2) - 1e0);
		v1 = v2 * (r2 * (1e0 + r1/f1)) / 
		          (r1 * (1e0 + r2/f2));

		DrlFPrintf(NULL, "VOLADJ: "
			" r1=%7.4f f1=%3.1f v1=%7.4f"
			" r2=%7.4f f2=%3.1f v2=%7.4f\n",
			r1*1e2, f1, v1*1e2, r2*1e2, f2, v2*1e2);

		*rateVolAtm = v1;
	}
*/


	/* If smile ON, compute smile adjustment factor */
	if (smlData != NULL) {
		IF_FAILED_DONE(DrlTSmile3DataInterp(
			zcCurve,
			smlData,
			todayDate,
			strikeRate,
			rateFreq,
			rateDcc,	
			rateObsDate,
			rateMatYrs,
			(int)swMat->swapPayFreq,
			*rateVolAtm,
			&rateVolSmileAdj));

		*rateVolSmile = *rateVolAtm + rateVolSmileAdj;
	} else {
		*rateVolSmile = *rateVolAtm + rateVolSmileAdj;
	}


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}


/*f-------------------------------------------------------------
 * Convenience routine to compute
 * the convexity/delay adjustement of a forward rate
 * paid over an IRR annuity of a different rate.
 */


static int
computeAdjustmentRateIRRPay(
	TDateInterval rateMat,		/* (I) rate maturity interval */
	int rateFreq,			/* (I) rate frequency */
	TDayCount rateDcc,		/* (I) rate day count convention */
	double rateFwd,			/* (I) rate forward (unadjusted) */
	double rateVolAtm,		/* (I) rate ATM volatility */

	TDateInterval payMat,		/* (I) payment maturity interval */
	int payFreq,			/* (I) payment frequency */
	double payRateFwd,		/* (I) payment rate forward (unadjusted) */
	double payRateVolAtm,		/* (I) payment rate ATM volatility */
	double corr,			/* (I) corr of index with pay rate */

	double tExp,			/* (I) time top expiration */

	double *rateFwdAdj)		/* (O) adjusted forward rate */
{
static	char	routine[] = "computeAdjustmentRateIRRPay";
	int	status = FAILURE;
	double	rateMatYrs,
		payMatYrs;

	/* Convert maturity to years */
	IF_FAILED_DONE(GtoDateIntervalToYears(&rateMat, &rateMatYrs));
	IF_FAILED_DONE(GtoDateIntervalToYears(&payMat,  &payMatYrs));

	/* Compute forward rate */
	*rateFwdAdj = rateFwd;


	/*
	 * Perform convexity adjustment.
	 */
	if (payFreq != rateFreq) {
		GtoErrMsg("%s: cannot handle rateFreq (%d) != payFreq (%d), sorry.\n",
			routine, rateFreq, payFreq);
		goto done;
	}

#ifdef	_OLD
	if (GtoIntervalEqualsInterval(&rateMat, &payMat) == FALSE)
	{
		double	rateDur, rateCvx,
			payDur, payCvx,
			rateAdj;


		IF_FAILED_DONE(GtoBondModDuration(
			rateFwd, rateFwd, (long)rateFreq, rateMatYrs,
			GTO_STUB_SIMPLE, &rateDur));
		IF_FAILED_DONE(GtoBondConvexity(
			rateFwd, rateFwd, (long)rateFreq, rateMatYrs,
			GTO_STUB_SIMPLE, &rateCvx));

		IF_FAILED_DONE(GtoBondModDuration(
			payRateFwd, payRateFwd, (long) payFreq, payMatYrs,
			GTO_STUB_SIMPLE, &payDur));
		IF_FAILED_DONE(GtoBondConvexity(
			payRateFwd, payRateFwd, (long) payFreq, payMatYrs,
			GTO_STUB_SIMPLE, &payCvx));

		rateAdj = rateFwd * (exp(
			      0.5e0 * (rateCvx/rateDur) * SQR(rateVolAtm)
				* rateFwd * tExp
			    - 0.5e0 * (payCvx  / payDur) * rateVolAtm * payRateVolAtm
				* payRateFwd * tExp)
			- 1e0);

		DrlFPrintStruct(NULL, NULL, '\t',
		    DRL_CVAR_T, "rateAdj",      DRL_PERCENT_T, (void*) &rateAdj,
		    DRL_NULL_T);

		*rateFwdAdj += rateAdj;
	}
#endif
	{
		double	rateDur, rateCvx,
		payDur, payCvx,
		rateAdj;


		IF_FAILED_DONE(GtoBondModDuration(
			rateFwd, rateFwd, (long)rateFreq, rateMatYrs,
			GTO_STUB_SIMPLE, &rateDur));
		IF_FAILED_DONE(GtoBondConvexity(
			rateFwd, rateFwd, (long)rateFreq, rateMatYrs,
			GTO_STUB_SIMPLE, &rateCvx));

		IF_FAILED_DONE(GtoBondModDuration(
			payRateFwd, payRateFwd, (long) payFreq, payMatYrs,
			GTO_STUB_SIMPLE, &payDur));
		IF_FAILED_DONE(GtoBondConvexity(
			payRateFwd, payRateFwd, (long) payFreq, payMatYrs,
			GTO_STUB_SIMPLE, &payCvx));

		rateAdj = rateFwd * (exp(
			      0.5e0 * (rateCvx/rateDur) * SQR(rateVolAtm)
				* rateFwd * tExp
			    - 0.5e0 * (payCvx  / payDur) * rateVolAtm * payRateVolAtm
				* corr * payRateFwd * tExp)
			- 1e0);

		DrlFPrintStruct(NULL, NULL, '\t',
		    DRL_CVAR_T, "rateAdj",      DRL_PERCENT_T, (void*) &rateAdj,
		    DRL_NULL_T);

		*rateFwdAdj += rateAdj;
	}




	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*----------------------------------------------------------------------
 *
 */

int
DriRLockOption(
	TDate rateObsDate,		/* (I) notification date */
	TDate rateEffDate,		/* (I) start date */

	TDateInterval rateMat1,	/* (I) rate maturity interval */
	int rateFreq1,			/* (I) rate frequency */
	TDayCount rateDcc1,		/* (I) day count convention */
	double rateWeight1,		/* (I) rate weight */

	TDateInterval rateMat2,	/* (I) rate maturity interval */
	int rateFreq2,			/* (I) rate frequency */
	TDayCount rateDcc2,		/* (I) day count convention */
	double rateWeight2,		/* (I) rate weight */

	TDate payDate,			/* (I) payment date */

	TDateInterval payMat,		/* (I) payment maturity interval */
	int payFreq,			/* (I) payment frequency */

	double strikeRate,		/* (I) strike rate for options */
	char *callPutType,		/* (I) "C"all, "P"ut, "F"orward */

	TDate todayDate,		/* (I) todayDates's date */
	TCurve *discZcCurve,		/* (I) discount zero curve */
	TCurve *indZcCurve1,		/* (I) index zero curve rate 1 */
	TCurve *indZcCurve2,		/* (I) index zero curve rate 2 */
	TSwaptionMatrix2D *swMat,	/* (I) CMS swaption matrix */
	TSmile3Data *smlData,		/* (I) smile data (or NULL) */
	VnfmData *vnfmData,		/* (I) calibration info */
	int noAdj,			/* (I) disable adjustment */

	double *pv)			/* (O) present value */
{
static	char	routine[] = "DriRTLockOption";
	int	status = FAILURE;

	double		rateMatYrs1,
			rateMatYrs2,
			rateFwd1,
			rateFwd2,
			rateFwdAdj1,
			rateFwdAdj2,
			rateVolAtm1,
			rateVolAtm2,
			rateVolSmile1,
			rateVolSmile2,
			corr12,
			corr01,
			corr02,

			tExp,		/* time to expiration */

			payRateFwd,
			payRateFwdAdj,
			payRateVolAtm,
			payRateVolSmile,

			payMatYrs,
			annRateFwd,
			annRateFwdAdj,

			optFv,
			discFact,
			annuityAdj,		/* IRR convexity adjustment */
			annuity;

	char		optType[4];


#ifndef	NO_ADJUSTMENTS
	double		pD, pC, pQ,		/* par sensitivities of pay rate */
			r1D, r1C, r1Q,		/* par sensitivities of rate 1 */
			r2D, r2C, r2Q;		/* par sensitivities of rate 2 */
#endif


	/**
	 ** CHECK INPUTS
	 **/
	if (rateObsDate > rateEffDate) {
		GtoErrMsg("%s: rateObsDate (%s) > rateEffDate (%s).\n",
			routine,
			DrlTDatePrint(NULL, rateObsDate),
			DrlTDatePrint(NULL, rateEffDate));
		goto done;
	}
	if (rateObsDate > payDate) {
		GtoErrMsg("%s: rateObsDate (%s) > payDate (%s).\n",
			routine,
			DrlTDatePrint(NULL, rateObsDate),
			DrlTDatePrint(NULL, payDate));
		goto done;
	}

	if (todayDate >= rateObsDate) {
		GtoErrMsg("%s: todayDate (%s) >= rateObsDate (%s).\n",
			routine,
			DrlTDatePrint(NULL, todayDate),
			DrlTDatePrint(NULL, rateObsDate));
		goto done;
	}

	/* Call/price = put/yield, etc. */
	if (ISSTR(callPutType, "C")) {
		strcpy(optType, "P");
	} else if (ISSTR(callPutType, "P")) {
		strcpy(optType, "C");
	} else {
		GtoErrMsg("%s: bad callPutType `%s' (only C or P).\n",
			routine, callPutType);
		goto done;
	}

	/* time to expiration */
	IF_FAILED_DONE( GtoDayCountFraction(
		todayDate,
		rateObsDate,
		GTO_ACT_365F,
		&tExp));

	/**
	 ** COMPUTE UNADJUSTED PAID RATE AS WEIGHTED AVERAGE OF TWO RATES
	 **/

	IF_FAILED_DONE( GtoDateIntervalToYears(&rateMat1,  &rateMatYrs1));
	IF_FAILED_DONE( GtoDateIntervalToYears(&rateMat2,  &rateMatYrs2));

	IF_FAILED_DONE( computeFwdRateVol(
		todayDate,
		indZcCurve1,
		swMat,
		smlData,
		
		rateMat1,
		rateFreq1,
		rateDcc1,
		rateEffDate,

		rateObsDate,
		strikeRate,

		&rateFwd1,
		&rateVolAtm1,
		&rateVolSmile1));


	IF_FAILED_DONE( computeFwdRateVol(
		todayDate,
		indZcCurve2,
		swMat,
		smlData,
		
		rateMat2,
		rateFreq2,
		rateDcc2,
		rateEffDate,

		rateObsDate,
		strikeRate,

		&rateFwd2,
		&rateVolAtm2,
		&rateVolSmile2));

	payRateFwd = rateWeight1 * rateFwd1 + rateWeight2 * rateFwd2;

	/* Compute volatility and correlations of paid rate */
	IF_FAILED_DONE( VnfmAvgQBCorr(
		vnfmData,
		tExp,
		tExp, rateMatYrs1, rateFreq1,
		tExp, rateMatYrs2, rateFreq2,
		&corr12));


	payRateVolSmile = sqrt(
		  SQR(rateWeight1 * rateFwd1 * rateVolSmile1)
		+ SQR(rateWeight2 * rateFwd2 * rateVolSmile2)
		+ 2e0 * corr12
		      * rateWeight1 * rateFwd1 * rateVolSmile1
		      * rateWeight2 * rateFwd2 * rateVolSmile2) /
			payRateFwd;

	payRateVolAtm = sqrt(
		  SQR(rateWeight1 * rateFwd1 * rateVolAtm1)
		+ SQR(rateWeight2 * rateFwd2 * rateVolAtm2)
		+ 2e0 * corr12
		      * rateWeight1 * rateFwd1 * rateVolAtm1
		      * rateWeight2 * rateFwd2 * rateVolAtm2) /
			payRateFwd;



	corr01 =  (  rateWeight1 * SQR(rateFwd1 * rateVolAtm1)
		   + rateWeight2 * corr12 
		      * rateFwd1 * rateVolAtm1
		      * rateFwd2 * rateVolAtm2)
		/ (rateFwd1 * rateVolAtm1 * payRateFwd * payRateVolAtm);

	corr02 =  (  rateWeight2 * SQR(rateFwd2 * rateVolAtm2)
		   + rateWeight1 * corr12 
		      * rateFwd1 * rateVolAtm1
		      * rateFwd2 * rateVolAtm2)
		/ (rateFwd2 * rateVolAtm2 * payRateFwd * payRateVolAtm);


	/**
	 ** COMPUTE FORWARD RATES ADJUSTED FOR CONVEXITY/DELAY
	 **/

	IF_FAILED_DONE( computeAdjustmentRateIRRPay(
		rateMat1,
		rateFreq1,
		rateDcc1,
		rateFwd1,
		rateVolAtm1,
		payMat,
		payFreq,
		payRateFwd,
		payRateVolAtm,
		corr01,
		tExp,
		&rateFwdAdj1));

	IF_FAILED_DONE( computeAdjustmentRateIRRPay(
		rateMat2,
		rateFreq2,
		rateDcc2,
		rateFwd2,
		rateVolAtm2,
		payMat,
		payFreq,
		payRateFwd,
		payRateVolAtm,
		corr02,
		tExp,
		&rateFwdAdj2));


	payRateFwdAdj = rateWeight1 * rateFwdAdj1 + rateWeight2 * rateFwdAdj2;



#ifdef	__DEBUG__
#endif
	{
	double	cvxAdj;
	DrlFPrintf(NULL, "%s: input data:\n", routine);
	cvxAdj = rateFwdAdj1 - rateFwd1;
	DrlFPrintStruct(NULL, NULL, '\t',
	    DRL_CVAR_T, "rateMat1",      DRL_TDATEINTERVAL_T, (void*) &rateMat1,
	    DRL_CVAR_T, "rateDcc1",      DRL_TDAYCOUNT_T, (void*) &rateDcc1,
	    DRL_CVAR_T, "rateFreq1",     DRL_INT_T,       (void*) &rateFreq1,
	    DRL_CVAR_T, "rateWeight1",   DRL_DOUBLE_T,    (void*) &rateWeight1,
	    DRL_CVAR_T, "rateFwd1",      DRL_PERCENT_T,   (void*) &rateFwd1,
	    DRL_CVAR_T, "cvxAdj",        DRL_PERCENT_T,   (void*) &cvxAdj,
	    DRL_CVAR_T, "rateFwdAdj1",   DRL_PERCENT_T,   (void*) &rateFwdAdj1,
	    DRL_CVAR_T, "rateVolAtm1",   DRL_PERCENT_T,   (void*) &rateVolAtm1,
	    DRL_CVAR_T, "rateVolSmile1", DRL_PERCENT_T,   (void*) &rateVolSmile1,
	    DRL_NULL_T);
	cvxAdj = rateFwdAdj2 - rateFwd2;
	DrlFPrintStruct(NULL, NULL, '\t',
	    DRL_CVAR_T, "rateMat2",      DRL_TDATEINTERVAL_T, (void*) &rateMat2,
	    DRL_CVAR_T, "rateDcc2",      DRL_TDAYCOUNT_T, (void*) &rateDcc2,
	    DRL_CVAR_T, "rateFreq2",     DRL_INT_T,       (void*) &rateFreq2,
	    DRL_CVAR_T, "rateWeight2",   DRL_DOUBLE_T,    (void*) &rateWeight2,
	    DRL_CVAR_T, "rateFwd2",      DRL_PERCENT_T,   (void*) &rateFwd2,
	    DRL_CVAR_T, "cvxAdj",        DRL_PERCENT_T,   (void*) &cvxAdj,
	    DRL_CVAR_T, "rateFwdAdj2",   DRL_PERCENT_T,   (void*) &rateFwdAdj2,
	    DRL_CVAR_T, "rateVolAtm2",   DRL_PERCENT_T,   (void*) &rateVolAtm2,
	    DRL_CVAR_T, "rateVolSmile2", DRL_PERCENT_T,   (void*) &rateVolSmile2,
	    DRL_NULL_T);

	DrlFPrintStruct(NULL, NULL, '\t',
	    DRL_CVAR_T, "corr12",        DRL_DOUBLE_T,    (void*) &corr12,
	    DRL_CVAR_T, "corr01",        DRL_DOUBLE_T,    (void*) &corr01,
	    DRL_CVAR_T, "corr02",        DRL_DOUBLE_T,    (void*) &corr02,
	    DRL_NULL_T);

	DrlFPrintStruct(NULL, NULL, '\t',
	    DRL_CVAR_T, "payMat",          DRL_TDATEINTERVAL_T, (void*) &payMat,
	    DRL_CVAR_T, "payFreq",         DRL_INT_T,       (void*) &payFreq,
	    DRL_CVAR_T, "payRateFwd",      DRL_PERCENT_T,   (void*) &payRateFwd,
	    DRL_CVAR_T, "payRateFwdAdj",   DRL_PERCENT_T,   (void*) &payRateFwdAdj,
	    DRL_CVAR_T, "payRateVolAtm",   DRL_PERCENT_T,   (void*) &payRateVolAtm,
	    DRL_CVAR_T, "payRateVolSmile", DRL_PERCENT_T,   (void*) &payRateVolSmile,
	    DRL_NULL_T);

	}


	/**
	 ** COMPUTE ANNUITY RATE & ANNUITY WITH POSSIBLE ADJUSTMENT 
	 **/

	IF_FAILED_DONE(GtoDateIntervalToYears(&payMat,    &payMatYrs));
	annRateFwd    = rateWeight1 * rateFwd1    + rateWeight2 * rateFwd2;

#ifndef	NO_ADJUSTMENTS
	/* Compute payment sensitivites */
	IF_FAILED_DONE(bndParSens(
		payMatYrs, payFreq, annRateFwd,
		&pD, &pC, &pQ));

	annuity = pD;

	DrlFPrintf(NULL, "annuity unadjusted = %12.8f\n", annuity);


	/* Adjust it for IRR convexity.  */
	if (IS_ALMOST_ZERO(rateWeight2)) {
		/* Compute bnd sensitivities */
		IF_FAILED_DONE(bndParSens(
			rateMatYrs1, rateFreq1, rateFwdAdj1,
			&r1D, &r1C, &r1Q));

		annuityAdj = 0.5 * SQR(payRateFwdAdj*payVolAtm) * tExp *
			(pQ / 3e0 - 0.5 * pC * r1C / r1D);

	} else if (IS_ALMOST_ZERO(rateWeight1)) {
		/* Compute bnd sensitivities */
		IF_FAILED_DONE(bndParSens(
			rateMatYrs2, rateFreq2, rateFwdAdj2,
			&r2D, &r2C, &r2Q));

		annuityAdj = 0.5 * SQR(payRateFwdAdj*payVolAtm) * tExp *
			(pQ / 3e0 - 0.5 * pC * r2C / r2D);
	} else {
		/* in this case, we simply use the IRR sensitivities
		 * It is not really correct..., but there is not
		 * much else to do
		 */

		annuityAdj = 0.5 * SQR(payRateFwdAdj*payVolAtm) * tExp *
			(pQ / 3e0 - 0.5 * SQR(pD) / pD);
	}
	DrlFPrintf(NULL, "annuity adjustment = %12.8f\n", annuityAdj);
	annuity += annuityAdj;
#else
	/* No adjustments */
	annRateFwdAdj = annRateFwd;
	annuity = (1 - pow(1e0 + annRateFwdAdj / payFreq, - payMatYrs*payFreq))
				/ annRateFwdAdj;
#endif



	/* compute discount rate */
	IF_FAILED_DONE( GtoDiscountDate(
		payDate,
		discZcCurve,
		GTO_LINEAR_INTERP,
		&discFact));

	/*
	 * Black Analytics (log normal)
	 */
	IF_FAILED_DONE( DrlBlack(
		tExp,
		payRateFwdAdj,
		payRateVolSmile,
		strikeRate,
		optType,
		"P",
		&optFv));


	/* Compute PV */
	*pv = optFv * annuity * discFact;



#ifdef	__DEBUG__
#endif
	DrlFPrintStruct(NULL, NULL, '\t',
	    DRL_CVAR_T, "strikeRate",  DRL_PERCENT_T, (void*) &strikeRate,
	    DRL_CVAR_T, "tExp",        DRL_DOUBLE_T,  (void*) &tExp,
	    DRL_CVAR_T, "optType",     DRL_CHAR_ARRAY_T, (void*) optType,

	    DRL_CVAR_T, "optFv",       DRL_PERCENT_T, (void*) &optFv,
	    DRL_CVAR_T, "annRateFwdAdj",DRL_DOUBLE_T,  (void*) &annRateFwdAdj,
	    DRL_CVAR_T, "annuity",     DRL_DOUBLE_T,  (void*) &annuity,
	    DRL_CVAR_T, "discFact",    DRL_DOUBLE_T,  (void*) &discFact,
	    DRL_CVAR_T, "optFv",       DRL_PERCENT_T, (void*) &optFv,
	    DRL_NULL_T);

	DrlFPrintf(NULL, "%s: MTM = %14.10f \n", routine, *pv);


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}






/*----------------------------------------------------------------------
 */

int
DriRLockOptionW(char *dataFnam)
{
static	char	routine[] = "DriRLockOption";
	int	status = FAILURE;


	TDate		rateEffDate;
	TDate		rateObsDate;

	TDateInterval	rateMat1;
	int		rateFreq1;
	TDayCount	rateDcc1;
	double		rateWeight1;

	TDateInterval	rateMat2;
	int		rateFreq2;
	TDayCount	rateDcc2;
	double		rateWeight2;

	TDate		payDate;
	TDateInterval	payMat;
	int		payFreq;

	double		strikeRate;
	char		callPutType[64],
			dummy[256];

	char		pricingOptions[1024];

	char		smileCurCode[1024];
	char		smileFeedDir[1024];
	int		noAdj;


	FILE		*fp = NULL;
static	char		defDataFnam[] = "kotlockw.dat";
	TDrWrapperData	*drWrap = NULL;
	TSmile3Data	*smlData = NULL;
	VnfmData	*vnfmData = NULL;


	double		pv;



	/* read deal data */
#define	READ_DATA(type,ptr,str)	\
		{ if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
		    { GtoErrMsg("%s: can't read %s.\n", routine, str); \
		    goto done;}}


	if (dataFnam == NULL) dataFnam = defDataFnam;

	if ((fp = fopen(dataFnam, "r")) == NULL) {
		GtoErrMsg("%s: can't open `%s' (%s).\n",
                        routine, dataFnam, strerror(errno));
		goto done;
	}


	READ_DATA(DRL_TDATE_T,         &rateObsDate,	"rate observation date");
	READ_DATA(DRL_TDATE_T,         &rateEffDate,	"rate effective date");
	READ_DATA(DRL_TDATE_T,         &payDate,	"payment date");

	READ_DATA(DRL_TDATEINTERVAL_T, &rateMat1,	"rate 1 maturity interval");
	READ_DATA(DRL_INT_T,           &rateFreq1,	"rate 1 frequency");
	READ_DATA(DRL_TDAYCOUNT_T,     &rateDcc1,	"rate 1 day count conv");
	READ_DATA(DRL_DOUBLE_T,        &rateWeight1,	"rate 1 weight ");

	READ_DATA(DRL_TDATEINTERVAL_T, &rateMat2,	"rate 2 maturity interval");
	READ_DATA(DRL_INT_T,           &rateFreq2,	"rate 2 frequency");
	READ_DATA(DRL_TDAYCOUNT_T,     &rateDcc2,	"rate 2 day count conv");
	READ_DATA(DRL_DOUBLE_T,        &rateWeight2,	"rate 2 weight ");

	READ_DATA(DRL_PERCENT_T,       &strikeRate,	"strike rate");
	READ_DATA(DRL_CHAR_ARRAY_T,    callPutType,	"call/put type");

	READ_DATA(DRL_TDATEINTERVAL_T, &payMat,		"annuity maturity interval");
	READ_DATA(DRL_INT_T,           &payFreq,	"annuity frequency");

	READ_DATA(DRL_CHAR_ARRAY_T,    dummy,		"first calibration index");
	READ_DATA(DRL_CHAR_ARRAY_T,    dummy,		"second calibration index");


	READ_DATA(DRL_CHAR_ARRAY_T,    pricingOptions,	"PricingOptions");

	/*
	READ_DATA(DRL_CHAR_ARRAY_T, smileCurCode,	"Smile Currency Code");
	READ_DATA(DRL_CHAR_ARRAY_T, smileFeedDir,	"Smile Feed Directory");
	*/


	/* read market data */
	if (DriTDrWrapperDataGetFull(
		NULL,
		DRI_DRW_TYPE2_3CURVES,
		&drWrap) != SUCCESS)
			goto done;


	/*
	 * Parse pricing options
	 */
	strcpy(smileCurCode, "nil");
	strcpy(smileFeedDir, "nil");
	noAdj = FALSE;

	IF_FAILED_DONE( DrlStrParseOptions(
		pricingOptions,
		DRL_CHAR_ARRAY_T, (void*) smileCurCode, "smileCurCode",
		DRL_CHAR_ARRAY_T, (void*) smileFeedDir, "smileFeedDir",
		DRL_INT_T,        (void*) &noAdj,       "noAdj",
		DRL_NULL_T));




	/* Read smile data from directory */
	IF_FAILED_DONE(DrlStrSubsEnv(smileFeedDir));

	if ((smlData = DrlTSmile3DataNewEmpty()) == NULL)
		goto done;

	if (DrlTSmile3DataReadFromExport(
		smlData,
		smileCurCode,
		smileFeedDir,
		drWrap->fToday) != SUCCESS)
			goto done;


	/* Get calibration info */
	vnfmData = DriTDrWrapperDataGetVnfmData(
        	drWrap,
		3,			/* number of factors */
		VNFM_LOGNORMAL_DIST,	/* back bone q */
		0);			/* 0=flat, 1=bv dates, 2=swvol dates */
	if (vnfmData == NULL) goto done;


	/* call pricing routine */
	if (DriRLockOption(
		rateObsDate,
		rateEffDate,

		rateMat1,
		rateFreq1,
		rateDcc1,
		rateWeight1,
		rateMat2,
		rateFreq2,
		rateDcc2,
		rateWeight2,

		payDate,
		payMat,
		payFreq,

		strikeRate,
		callPutType,

		drWrap->fToday,
		drWrap->fDiscZcCurve,
		drWrap->fZcCurve,		/* zero.dat is 1st index curve */
		drWrap->fRiskZcCurve,		/* riskzero.dat is 2nd index curve */
		drWrap->fCmsSwMat,
		smlData,
		vnfmData,
		noAdj,
		&pv) != SUCCESS)
			goto done;

	/* */
	if (DriTDrWrapperDataPutPrice(pv) != SUCCESS)
		goto done;


	/* */
	status = SUCCESS;
done:
	if (fp) fclose(fp);
	DrlTSmile3DataDelete(smlData);
	DriTDrWrapperDataFree(drWrap);
	VnfmFree(vnfmData);
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}



