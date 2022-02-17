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
#include <ctype.h>

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

#include "check.h"
#include "zr2coup.h"		/* Analytics C Library */
#include "zr2simp.h"		/* Analytics C Library */
#include "duration.h"
#include "convex.h"
#include "swapadj.h"

#include "cashflow.h"
#include "barrier.h"
#include "barbin.h"

#include "drlio.h"
#include "drlstr.h"
#include "drlsmat.h"

#include "dritkwrp.h"		/* TDrWrapperData routines */

#include "driwbrt.h"		/* Prototype consistency */

#define	__DEBUG__
#undef	__DEBUG__

#define	ISSTR(s1, s2)	(!strcmp((s1),(s2)))



/*----------------------------------------------------------------------
 * Convenience routine to compute a forward rate.
 */

static	int
computeFwdRate(
	TCurve *zcCurve,	/* (I) zero curve */

	TDateInterval rateMat,	/* (I) rate maturity interval */
	int rateFreq,		/* (I) rate frequency */
	TDayCount rateDcc,	/* (I) rate day count convention */

	TDate resetDate,	/* (I) */

	double *fwdRate)	/* (O) */
{
static	char	routine[] = "computeFwdRate";
	TDateInterval	payInterval;
	TDate		maturityDate;
	int		status = 0;


	if (GtoDtFwdAny(resetDate, &rateMat, &maturityDate)
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
			resetDate,
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
			resetDate,
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



/*f-------------------------------------------------------------
 * COnvenience routine to compute
 * the convexity adjustement of a forward rate
 * on a zero curve.
 */


static int
computeFwdAdjRate(
	TCurve *zcCurve,	/* (I) zero curve (base date only used) */

	TDateInterval rateMat,	/* (I) rate maturity interval */
	int rateFreq,		/* (I) rate frequency */
	TDayCount rateDcc,	/* (I) rate day count convention */

	TDate resetDate,	/* (I) reset date */

	double yldVol,		/* (I) yield volatility */
	double *fwdAdjRate)	/* (O) cvx adjustment */
{
static	char	routine[] = "computeFwdAdjRate";
	int	status = FAILURE;
	TDate	maturityDate;
	long	numDays;
	long	den;
	double	fwdRate, tExp, tMat, dur, cvx;

	/* Compute forward rate */
	if (computeFwdRate(
		zcCurve,
		rateMat,
		rateFreq,
		rateDcc,
		resetDate,
		&fwdRate) != SUCCESS)
			goto done;


	/* Compute time to reset */
	if (GtoDayCountFraction(
		zcCurve->fBaseDate,
		resetDate,
		GTO_ACT_365F, &tExp) != SUCCESS)
			goto done;

	/* compute rate maturity date */
	if (GtoDtFwdAny(	
		resetDate,
		&rateMat,	
		&maturityDate)
		!= SUCCESS)
			goto done;

	if (GtoDateIntervalToYears(&rateMat, &tMat) != SUCCESS)
		goto done;


	/* compute dur and cvx */
	switch (rateFreq) {
	case 0:
		/* Simple Rate */
		if (GtoDaysInYearFromDayCountConv(rateDcc, &den)
			!= SUCCESS) goto done;

		if (GtoDaysDiff(resetDate, maturityDate, rateDcc,
			&numDays) != SUCCESS) goto done;

		if (GtoMoneyMarketModDuration(
			fwdRate,
			numDays,
			(long) den,
			&dur) != SUCCESS) goto done;

		if (GtoMoneyMarketConvexity(
			fwdRate,
			numDays,
			(long) den,
			&cvx) != SUCCESS) goto done;
		break;
	case 1:
	case 2:
	case 4:
	case 12:
		/* Compounded Rate */
		if (GtoBondModDuration(
			fwdRate,
			fwdRate,
			(long) rateFreq,
			tMat,
			GTO_STUB_BOND,
			&dur) != SUCCESS) goto done;

		if (GtoBondConvexity(
			fwdRate,
			fwdRate,
			(long) rateFreq,
			tMat,
			GTO_STUB_BOND,
			&cvx) != SUCCESS) goto done;
		break;
	default:
		GtoErrMsg("%s: bad frequency %d.\n", routine, rateFreq);
		goto done;
	}


	/* compute convexity adjustment */
	if (GtoSwapConvexityAdj(
		fwdRate,
		dur,
		cvx,
		yldVol,
		tExp,
		fwdAdjRate) != SUCCESS)
			goto done;


	(*fwdAdjRate) += fwdRate;

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}




/*f---------------------------------------------------------------------
 * Pricing routine for American binary on a spread.
 * Returns SUCCESS/FAILURE.
 */

int
DriRateBarrierBinary(
	TDate expDate,		/* (I) expiration date */
	TDate resetDate,	/* (I) rate reset date */

	TDateInterval rateMat,	/* (I) rate maturity interval */
	int rateFreq,		/* (I) rate frequency */
	TDayCount rateDcc,	/* (I) rate day count convention */

	double strike,		/* (I) used if asset */
	double barrier,		/* (I) barrier level */

	char *callPutS,		/* (I) "C"all, "P"ut */
	char *upDownS,		/* (I) "U"p, "D"own */
	char *inOutS,		/* (I) "I"n, "O"ut */
	char *assetCashS,	/* (I) "A"sset, "C"ash */
	char *hitExpS,		/* (I) "H"it, "E"xp */

	double spdVol,		/* (I) spread volatility in bp/day */
	char *distTypeS,	/* (I) "L"ognormal, "N"ormal */

	TDate today,		/* (I) todays's date */
	TCurve *discZcCurve,	/* (I) discount zero curve */
	TCurve *indxZcCurve,	/* (I) cmt zero curve */
	TSwaptionMatrix2D *swMat,/* (I) CMS swaption matrix */

	double *pv)		/* (O) present value */
{
static	char	routine[] = "DrlRateBarrierBinary";
	int	status = FAILURE;

	long   callPutType;      /* (I) GTO_OPTION_CALL or GTO_OPTION_PUT */
	long   upDownType;       /* (I) GTO_BBIN_UP     or GTO_BBIN_DOWN */
	long   inOutType;        /* (I) GTO_BBIN_IN     or GTO_BBIN_OUT */
	long   assetCashType;    /* (I) GTO_BBIN_ASSET  or GTO_BBIN_CASH */
	long   atExpHitType;     /* (I) GTO_BBIN_AT_EXP or GTO_BBIN_AT_HIT */

	TDate		maturityDate;
	double		spotYIndx, fwdYIndx,
			spotY1, fwdY1,
			yldVol,
			discRate, discFact;
	double		correlation = 1e0,	/*  underlying 1 and 2 */
			volatilityY1LN,
			yuShift,
			volatilityY1,
			tExp, tPay, tMat,
			price,
			rebateProb;
	char		distType;

	int		optionType;	/* (I) GTO_BAR_UP_CALL,
						GTO_BAR_DOWN_PUT,etc */
	int		inOut;		/* (I) GTO_BAR_IN, GTO_BAR_OUT */ 


	/* */
	if (DrlStrLongValueScan(callPutS, "call/put", &callPutType,
		"C", (long) GTO_OPTION_CALL,
		"P", (long) GTO_OPTION_PUT,
		NULL) != SUCCESS) goto done;

	if (DrlStrLongValueScan(upDownS, "up/down", &upDownType,
		"U", (long) GTO_BBIN_UP,
		"D", (long) GTO_BBIN_DOWN,
		NULL) != SUCCESS) goto done;

	if (DrlStrLongValueScan(inOutS, "in/out", &inOutType,
		"I", (long) GTO_BBIN_IN,
		"O", (long) GTO_BBIN_OUT,
		NULL) != SUCCESS) goto done;

	if (DrlStrLongValueScan(assetCashS, "asset/cash", &assetCashType,
		"A", (long) GTO_BBIN_ASSET,
		"C", (long) GTO_BBIN_CASH,
		NULL) != SUCCESS) goto done;

	if (DrlStrLongValueScan(hitExpS, "asset/cash", &atExpHitType,
		"E", (long) GTO_BBIN_AT_EXP,
		"H", (long) GTO_BBIN_AT_HIT,
		NULL) != SUCCESS) goto done;

	distType = toupper(distTypeS[0]);
	switch (distType) {
	case 'N':
		if (assetCashType == GTO_BBIN_ASSET) {
		    GtoErrMsg("%s: N distType not supported with Asset "
			"(use Cash).\n", routine);
		    goto done;
		}
		break;
	case 'L':
		break;
	default:
		GtoErrMsg("%s: unknown dist type `%s' (L or N).\n",
			routine, distTypeS);
		goto done;
	}



	/* interp volatility */
	if (GtoDateIntervalToYears(&rateMat, &tMat) != SUCCESS)
		goto done;

	if (DrlTSwaptionMatrix2DInterpDate(
		swMat,
		&yldVol,
		today,
		expDate,
		tMat,
		FALSE)	/* (I) TRUE=adjoint, FALSE=direct*/
			!= SUCCESS)
				goto done;


	/* compute adjusted rates */

	if (computeFwdRate(
		indxZcCurve,
		rateMat,
		rateFreq,
		rateDcc,
		discZcCurve->fBaseDate,
		&spotYIndx) != SUCCESS)
			goto done;


	if (computeFwdAdjRate(
		indxZcCurve,
		rateMat,
		rateFreq,
		rateDcc,
		resetDate,
		yldVol,		/* (I) yield volatility */
		&fwdYIndx) != SUCCESS)
			goto done;



	spotY1 = spotYIndx;
	fwdY1  = fwdYIndx;


	/* compute discount rate */
	if (GtoDiscountDate(
			resetDate,
			discZcCurve,
			GTO_LINEAR_INTERP,
			&discFact) != SUCCESS)
				goto done;

	/* compute discount rate */
	if (GtoInterpRate(
			resetDate,
			discZcCurve,
			GTO_LINEAR_INTERP,
			&discRate) != SUCCESS)
				goto done;

	/* compute time to expiration and payment */
	if (GtoDayCountFraction(
		today,
		expDate,
		GTO_ACT_365F,
		&tExp) != SUCCESS)
			goto done;

	if (GtoDayCountFraction(
		discZcCurve->fBaseDate,
		resetDate,
		GTO_ACT_365F,
		&tPay) != SUCCESS)
			goto done;


	/* Convert volatility */
	volatilityY1LN = spdVol * 1e-4 * 15.8745e0 / fwdY1;
	volatilityY1 = volatilityY1LN;

	/* Perform the Yu-Transform */
	if (distType == 'L') {
		yuShift = 0e0;
	} else if (distType == 'N') {
		yuShift = 100e0;
		volatilityY1 *= fwdY1 / (fwdY1 + yuShift);
		fwdY1 += yuShift;
		spotY1 += yuShift;
		barrier += yuShift;
	} else {
		goto done;
	}




#ifdef	__DEBUG__
#endif
	DrlFPrintf(NULL, "%s: input data:\n", routine);

	DrlFPrintStruct(NULL, NULL, '\t',
	    DRL_CVAR_T, "expDate",     DRL_TDATE_T, (void*) &expDate,
	    DRL_CVAR_T, "resetDate",   DRL_TDATE_T, (void*) &resetDate,
	    DRL_CVAR_T, "rateMat",     DRL_TDATEINTERVAL_T, (void*) &rateMat,
	    DRL_CVAR_T, "rateFreq",    DRL_INT_T, (void*) &rateFreq,
	    DRL_CVAR_T, "rateDcc",     DRL_TDAYCOUNT_T, (void*) &rateDcc,
	    DRL_NULL_T);

	DrlFPrintStruct(NULL, NULL, '\t',
	    DRL_CVAR_T, "fwdYIndx",     DRL_PERCENT_T, (void*) &fwdYIndx,
	    DRL_CVAR_T, "spotYIndx",    DRL_PERCENT_T, (void*) &spotYIndx,
	    DRL_CVAR_T, "fwdY1",        DRL_PERCENT_T, (void*) &fwdY1,
	    DRL_CVAR_T, "spotY1",       DRL_PERCENT_T, (void*) &spotY1,
	    DRL_CVAR_T, "volY1",        DRL_PERCENT_T, (void*) &volatilityY1,
	    DRL_CVAR_T, "distType",     DRL_CHAR_T,    (void*) &distType,
	    DRL_CVAR_T, "yuShift",      DRL_PERCENT_T, (void*) &yuShift,
	    DRL_CVAR_T, "strike",       DRL_PERCENT_T, (void*) &strike,
	    DRL_CVAR_T, "barrier",      DRL_PERCENT_T, (void*) &barrier,
	    DRL_CVAR_T, "correlation",  DRL_DOUBLE_T,  (void*) &correlation,
	    DRL_CVAR_T, "tExp",	        DRL_DOUBLE_T,  (void*) &tExp,
	    DRL_CVAR_T, "tPay",	        DRL_DOUBLE_T,  (void*) &tPay,
	    DRL_CVAR_T, "discRate",     DRL_DOUBLE_T,  (void*) &discFact,
	    DRL_NULL_T);









	/*
	 * The main routine for computing price of outside barrier21
	 * options with binary payoffs.
	 */
	if (GtoOutsideBarrierBinary(
		callPutType,	/* GTO_OPTION_CALL or GTO_OPTION_PUT */
		upDownType,	/* GTO_BBIN_UP     or GTO_BBIN_DOWN */
		inOutType,	/* GTO_BBIN_IN     or GTO_BBIN_OUT */
		assetCashType,	/* GTO_BBIN_ASSET  or GTO_BBIN_CASH */
		atExpHitType,	/* GTO_BBIN_AT_EXP or GTO_BBIN_AT_HIT */
		spotY1,		/* Volatility of asset 1 */
		fwdY1,		/* Fwd at Texp of asset used in payoff*/
		volatilityY1, 
		spotY1,		/* Price at T=0 of asset used w/ barrier */
		fwdY1,		/* Fwd at Texp of asset with barrier */
		volatilityY1,
		correlation, 	/* Correlation of asset 1 and 2 */
		strike,		/* Strike for asset 1 payoff function */
		barrier,	/* Barrier for asset 2 */
		tExp,		/* Num years to expiration */
		tPay,		/* Num years to payment */
		discRate,	/* Discount rate (ann. compound) */ 
		&price) != SUCCESS)
			goto done;



	/* */
	*pv = price;

	DrlFPrintf(NULL, "pv=%0.8f\n", *pv);


	/* */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}



/*f---------------------------------------------------------------------
 * Dr-Wrapper for {\tt DrlRateBarrierBinary}.
 * The argument {\tt dataFnam} specifies the name of the file
 * containing the data (if it is NULL, the default name
 * "wbinrate.dat" is used.
 * Returns SUCCESS/FAILURE.
 */

int
DriRateBarrierBinaryW(char *dataFnam)
{
static	char	routine[] = "DriRateBarrierBinaryW";
	int	status = FAILURE;


	TDate		expDate;
	TDate		resetDate;
	TDateInterval	rateMat;
	int		rateFreq;
	TDayCount	rateDcc;

	double		strike, barrier;

	char		callPutS[64];
	char		upDownS[64];
	char		inOutS[64];
	char		assetCashS[64];
	char		hitExpS[64];
	char		distTypeS[64];

	double		spdVol;

	FILE		*fp = NULL;
static	char		defDataFnam[] = "wbinrate.dat";
	TDrWrapperData	*drWrap = NULL;
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


	READ_DATA(DRL_TDATE_T, &expDate,	"expiration date");
	READ_DATA(DRL_TDATE_T, &resetDate,	"reset date");
	READ_DATA(DRL_TDATEINTERVAL_T, &rateMat,"rate maturity");
	READ_DATA(DRL_INT_T, &rateFreq,		"rate frequency");
	READ_DATA(DRL_TDAYCOUNT_T, &rateDcc,	"rate day count conv");

	READ_DATA(DRL_PERCENT_T, &strike,	"strike rate");
	READ_DATA(DRL_PERCENT_T, &barrier,	"barrier rate");


	READ_DATA(DRL_CHAR_ARRAY_T, callPutS,	"call/put type");
	READ_DATA(DRL_CHAR_ARRAY_T, upDownS,	"up/down type");
	READ_DATA(DRL_CHAR_ARRAY_T, inOutS,	"in/out type");
	READ_DATA(DRL_CHAR_ARRAY_T, assetCashS,	"asset/cash type");
	READ_DATA(DRL_CHAR_ARRAY_T, hitExpS,	"hit/exp type");

	READ_DATA(DRL_DOUBLE_T, &spdVol,	"spread volatility");
	READ_DATA(DRL_CHAR_ARRAY_T, distTypeS,	"distribution type (L,N)");


	/* read market data */
	if (DriTDrWrapperDataGet(NULL, &drWrap) != SUCCESS)
		goto done;



	/* Call pricing routine */
	if (DriRateBarrierBinary(
		expDate,
		resetDate,

		rateMat,
		rateFreq,
		rateDcc,

		strike,
		barrier,

		callPutS,
		upDownS,
		inOutS,	
		assetCashS,
		hitExpS,

		spdVol,
		distTypeS,

		drWrap->fToday,
		drWrap->fDiscZcCurve,
		drWrap->fZcCurve,
		drWrap->fCmsSwMat,
		&pv) != SUCCESS)
			goto done;


	/* */
	if (DriTDrWrapperDataPutPrice(pv) != SUCCESS)
		goto done;


	/* */
	status = SUCCESS;
done:
	if (fp) fclose(fp);
	DriTDrWrapperDataFree(drWrap);
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}



