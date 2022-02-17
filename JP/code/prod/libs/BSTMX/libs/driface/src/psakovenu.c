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

#include "swapkol.h"

#include "drlsmat.h"
#include "drltime.h"
#include "drlio.h"
#include "drloptio.h"
#include "drlstr.h" 		/* DrlStrSubsEnv */
#include "dritkwrp.h"

#include "drirlopt.h"		/* Prototype consistency */


#define	__DEBUG__
#undef	__DEBUG__



#undef	ARGSIZE
#define ARGSIZE(arg)	arg[0]








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
 */

int
DriPsaSwapKoVenuW(char *dataFnam)
{
static	char	routine[] = "DriPsaSwapKoVenuW";
	int	status = FAILURE;

	FILE		*fp = NULL;
static	char		defDataFnam[] = "psakovenu.dat";
static	char		logfile[] = "psakovenu_t.log";
	TDrWrapperData	*drWrap = NULL;
	double		pv;




#define	MAX_CF	360
	long     fixedCashDates[MAX_CF];
	double   fixedCashAmounts[MAX_CF];

	long     floatResetDates[MAX_CF];
	long     floatAccrueDates[MAX_CF];
	long     floatPayDates[MAX_CF];
	double   floatNotionals[MAX_CF];

	long     barrierDates[MAX_CF];
	double   loBarrierRates[MAX_CF];
	double   hiBarrierRates[MAX_CF];
	double   koCurveRebates[MAX_CF];

	double	numericScalars[GTO_SWKO_NUM_NUM_SCALARS+1];
	char	stringScalars[(GTO_SWKO_NUM_STR_SCALARS+1)*(GTO_MAX_STR_LEN+1)];

#define	NUM_RATES	1
	long	indexFreqs[NUM_RATES+1];
	char	indexMaturities[(NUM_RATES+1)*(GTO_MAX_STR_LEN+1)];
	char	indexDayCountStrs[(NUM_RATES+1)*(GTO_MAX_STR_LEN+1)];
	long	indexCurves[NUM_RATES+1];
	double	floatWeights[NUM_RATES+1];
	double	knockWeights[NUM_RATES+1];


	long     baseDates[3];
	long     discZeroDates[MAX_CF];
	double   discZeroRates[MAX_CF];
	long     indexZeroDates[MAX_CF];
	double   indexZeroRates[MAX_CF];

	long     volDates[MAX_CF];
	double   vols[MAX_CF];


	double   calibParams[MAX_CF];
	double   output[MAX_CF];




	int		idx, numFl, numItems;
	double		strike, dcf;
	TDayCount	basPayDcc,
			fltPayDcc;

	TDateInterval	basMat, fltMat;
	int		basFreq, fltFreq;
	TDayCount	basDcc, fltDcc;
	double		basFwd, fltFwd;
	TCurve		*basCurve, *fltCurve;







	/*
	 *
	 */

	/*
	 * Read market data
	 */
	if (DriTDrWrapperDataGetFull(
		NULL,
		DRI_DRW_TYPE2_2CURVES,
		&drWrap) != SUCCESS)
			goto done;

	basCurve = drWrap->fZcCurve;
	fltCurve = drWrap->fDiscZcCurve;

	/*
	 * Read deal data
	 */
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




	/*
	 * Read flows
	 */
	READ_DATA(DRL_INT_T, &numFl, "num flows");
	ASSERT_OR_DONE(numFl < MAX_CF);


	ARGSIZE(fixedCashDates) = numFl;
	ARGSIZE(fixedCashAmounts) = numFl;

	ARGSIZE(floatResetDates) = numFl;
	ARGSIZE(floatAccrueDates) = numFl;
	ARGSIZE(floatPayDates) = numFl;
	ARGSIZE(floatNotionals) = numFl;

	ARGSIZE(barrierDates) = numFl;
	ARGSIZE(loBarrierRates) = numFl;
	ARGSIZE(hiBarrierRates) = numFl;
	ARGSIZE(koCurveRebates) = numFl;

	for (idx=1; idx<=numFl; idx++) {

		READ_DATA(DRL_TDATE_T,
			(void*) &floatResetDates[idx], "floatResetDates");
		READ_DATA(DRL_TDATE_T,
			(void*) &floatAccrueDates[idx], "floatAccrueDates");
		READ_DATA(DRL_TDATE_T,
			(void*) &floatPayDates[idx], "floatPayDates");
		READ_DATA(DRL_DOUBLE_T,
			(void*) &strike, "strike");
		READ_DATA(DRL_DOUBLE_T,
			(void*) &floatNotionals[idx], "floatNotionals");

		barrierDates[idx] = floatResetDates[idx];

		READ_DATA(DRL_DOUBLE_T,
			(void*) &loBarrierRates[idx], "loBarrierRates");
		READ_DATA(DRL_DOUBLE_T,
			(void*) &hiBarrierRates[idx], "hiBarrierRates");
		READ_DATA(DRL_DOUBLE_T,
			(void*) &koCurveRebates[idx], "koCurveRebates");
	}

	/* Read pay dcc */

	READ_DATA(DRL_CHAR_ARRAY_T,
		(void*)&stringScalars[WRAP_STR_IDX(6)],
		"basis payment dcc");
	IF_FAILED_DONE( GtoStringToDayCountConv(
		(char*)&stringScalars[WRAP_STR_IDX(6)],
		&basPayDcc));

	READ_DATA(DRL_TDAYCOUNT_T,
		(void*)&fltPayDcc,
		"flt payment dcc");




	READ_DATA(DRL_CHAR_ARRAY_T,
		(void*)&stringScalars[WRAP_STR_IDX(7)],
		"fixed/float stub conv");


	/* Read rates */
	ARGSIZE(indexFreqs) = NUM_RATES;
	ARGSIZE(indexMaturities) = NUM_RATES;
	ARGSIZE(indexDayCountStrs) = NUM_RATES;
	ARGSIZE(indexCurves) = NUM_RATES;
	ARGSIZE(floatWeights) = NUM_RATES;
	ARGSIZE(knockWeights) = NUM_RATES;

	READ_DATA(DRL_CHAR_ARRAY_T,
		(void*)&indexMaturities[WRAP_STR_IDX(1)],
		"bas rate mat");
	READ_DATA(DRL_LONG_T,
		(void*)&indexFreqs[1],
		"bas rate freq");
	READ_DATA(DRL_CHAR_ARRAY_T,
		(void*)&indexDayCountStrs[WRAP_STR_IDX(1)],
		"bas rate dcc");
	indexCurves[1] = 2;
	floatWeights[1] = 1.0;
	knockWeights[1] = 1.0;



	IF_FAILED_DONE( GtoStringToDateInterval(
		(void*)&indexMaturities[WRAP_STR_IDX(1)],
		routine,
		&basMat));
	basFreq = indexFreqs[1];
	IF_FAILED_DONE( GtoStringToDayCountConv(
		(void*)&indexDayCountStrs[WRAP_STR_IDX(1)],
		&basDcc));

	READ_DATA(DRL_TDATEINTERVAL_T,
		(void*)&fltMat,
		"flt rate mat");
	READ_DATA(DRL_INT_T,
		(void*)&fltFreq,
		"flt rate freq");
	READ_DATA(DRL_TDAYCOUNT_T,
		(void*)&fltDcc,
		"flt rate dcc");



	READ_DATA(DRL_CHAR_ARRAY_T,
		(void*)&stringScalars[WRAP_STR_IDX(8)],
		"float stub rate type");


	DrlFilePrintf(logfile, "FORWARDS:\n"
		"  RESET     ACCSTART    PAY         BASFWD    FLTFWD     \n");



	/* Create fixed leg and V-method ! */
	for (idx=1; idx<=numFl; idx++) {

		fixedCashDates[idx] = floatPayDates[idx];

		IF_FAILED_DONE( GtoDayCountFraction(
			floatAccrueDates[idx],
			floatPayDates[idx],
			fltPayDcc,
			&dcf));

		fixedCashAmounts[idx] = strike * dcf;


		IF_FAILED_DONE( computeFwdRate(
			basCurve,
			basMat,
			basFreq,
			basDcc,
			floatResetDates[idx],
			&basFwd));

		IF_FAILED_DONE( computeFwdRate(
			fltCurve,
			fltMat,
			fltFreq,
			fltDcc,
			floatResetDates[idx],
			&fltFwd));

		fixedCashAmounts[idx] *= fltFwd;
		fixedCashAmounts[idx] *= floatNotionals[idx];


		loBarrierRates[idx] *= fltFwd;
		hiBarrierRates[idx] *= fltFwd;




		DrlFilePrintf(logfile,
			" %10s %10s %10s %10.6f %10.6f \n",
			DrlTDatePrint(NULL, floatResetDates[idx]),
			DrlTDatePrint(NULL, floatAccrueDates[idx]),
			DrlTDatePrint(NULL, floatPayDates[idx]),
			basFwd*1e2,
			fltFwd*1e2);

	}


	/* Numeric Scalars */
	ARGSIZE(numericScalars) = GTO_SWKO_NUM_NUM_SCALARS;

	READ_DATA(DRL_DOUBLE_T, (void*)&numericScalars[1],
			"numPpy");

	numericScalars[2] = (long) (drWrap->fBvCurve->fBasis);
	numericScalars[3] = 2;	/* BV interp type */

	READ_DATA(DRL_DOUBLE_T, (void*)&numericScalars[4],
			"firstFloatRate");




	/* String Scalars*/
	ARGSIZE(stringScalars) = GTO_SWKO_NUM_STR_SCALARS;

	READ_DATA(DRL_CHAR_ARRAY_T,
		(void*)&stringScalars[WRAP_STR_IDX(1)],
		"instrument type");

	READ_DATA(DRL_CHAR_ARRAY_T,
		(void*)&stringScalars[WRAP_STR_IDX(2)],
		"pay in/out");

	READ_DATA(DRL_CHAR_ARRAY_T,
		(void*)&stringScalars[WRAP_STR_IDX(3)],
		"barrier type");

	READ_DATA(DRL_CHAR_ARRAY_T,
		(void*)&stringScalars[WRAP_STR_IDX(4)],
		"date/interval");

	READ_DATA(DRL_CHAR_ARRAY_T,
		(void*)&stringScalars[WRAP_STR_IDX(5)],
		"knock in/out interval");




	/* Parameters */
	numItems = 6;
	ARGSIZE(calibParams) = 6;
	for (idx=1; idx<=numItems; idx++) {
		READ_DATA(DRL_DOUBLE_T,
			(void*)&calibParams[idx],
			"calibParams");
	}




	/*
	 * Market environment
	 */

	ARGSIZE(baseDates) = 1;
	baseDates[1] = basCurve->fBaseDate;

	ASSERT_OR_DONE(fltCurve->fNumItems < MAX_CF);
	ARGSIZE(discZeroDates) = fltCurve->fNumItems;
	ARGSIZE(discZeroRates) = fltCurve->fNumItems;
	for (idx=0; idx< fltCurve->fNumItems; idx++) {
		discZeroDates[idx+1] = fltCurve->fArray[idx].fDate;
		discZeroRates[idx+1] = fltCurve->fArray[idx].fRate;
	}


	ASSERT_OR_DONE(basCurve->fNumItems < MAX_CF);
	ARGSIZE(indexZeroDates) = basCurve->fNumItems;
	ARGSIZE(indexZeroRates) = basCurve->fNumItems;
	for (idx=0; idx< basCurve->fNumItems; idx++) {
		indexZeroDates[idx+1] = basCurve->fArray[idx].fDate;
		indexZeroRates[idx+1] = basCurve->fArray[idx].fRate;
	}


	/* Vol */
	READ_DATA(DRL_INT_T, &numItems, "num vol dates");
	ASSERT_OR_DONE(numItems < MAX_CF);
	ARGSIZE(volDates) = numItems;
	ARGSIZE(vols) = numItems;
	for (idx=1; idx<=numItems; idx++) {
		READ_DATA(DRL_TDATE_T,
			(void*) &volDates[idx], "volDates");
		READ_DATA(DRL_DOUBLE_T,
			(void*) &vols[idx], "vols");
	}






	/*
	 * LIL call
	 */
	ARGSIZE(output) = 1;


	IF_FAILED_DONE( DrlLilVectLoggingFile(
		logfile, "a", "SWAPKO",
		DRL_TDATE_L, fixedCashDates, "fixedCashDates",
		DRL_DOUBLE_L, fixedCashAmounts, "fixedCashAmounts",
		DRL_TDATE_L, floatResetDates, "floatResetDates",
		DRL_TDATE_L, floatAccrueDates, "floatAccrueDates",
		DRL_TDATE_L, floatPayDates, "floatPayDates",
		DRL_DOUBLE_L, floatNotionals, "floatNotionals",
		DRL_TDATE_L, barrierDates, "barrierDates",
		DRL_DOUBLE_L, loBarrierRates, "loBarrierRates",
		DRL_DOUBLE_L, hiBarrierRates, "hiBarrierRates",
		DRL_DOUBLE_L, koCurveRebates, "koCurveRebates",
		DRL_DOUBLE_L, numericScalars, "numericScalars",
		DRL_CHAR_BLOCK_L, stringScalars, "stringScalars",
		DRL_LONG_L, indexFreqs, "indexFreqs",
		DRL_CHAR_BLOCK_L, indexMaturities, "indexMaturities",
		DRL_CHAR_BLOCK_L, indexDayCountStrs, "indexDayCountStrs",
		DRL_LONG_L, indexCurves, "indexCurves",
		DRL_DOUBLE_L, floatWeights, "floatWeights",
		DRL_DOUBLE_L, knockWeights, "knockWeights",
		DRL_TDATE_L, baseDates, "baseDates",
		DRL_TDATE_L, discZeroDates, "discZeroDates",
		DRL_DOUBLE_L, discZeroRates, "discZeroRates",
		DRL_TDATE_L, indexZeroDates, "indexZeroDates",
		DRL_DOUBLE_L, indexZeroRates, "indexZeroRates",
		DRL_TDATE_L, volDates, "volDates",
		DRL_DOUBLE_L, vols, "vols",
		DRL_DOUBLE_L, calibParams, "calibParams",
		DRL_NULL_T));

	DrlFilePrintf(logfile, "FLTRATE: %3s %d %10s (pay %10s)\n",
			DrlTDateIntervalPrint(NULL, fltMat),
			fltFreq,
			DrlTDayCountPrint(NULL, fltDcc),
			DrlTDayCountPrint(NULL, fltPayDcc));
	DrlFilePrintf(logfile, "BASRATE: %3s %d %10s (pay %10s)\n",
			DrlTDateIntervalPrint(NULL, basMat),
			basFreq,
			DrlTDayCountPrint(NULL, basDcc),
			DrlTDayCountPrint(NULL, basPayDcc));





	IF_FAILED_DONE( GtoSwapKOModelL(
		fixedCashDates,        /* (I) Fixed cashflow dates */
		fixedCashAmounts,      /* (I) Fixed cashflow amounts */

		floatResetDates,       /* (I) Reset dates of floating side */
		floatAccrueDates,      /* (I) Accrue start dates of floating */
		floatPayDates,         /* (I) Payment dates of floating */
		floatNotionals,        /* (I) Floating notionals */

		barrierDates,          /* (I) Barrier/rebate dates      */
		loBarrierRates,        /* (I) Low barrier rates         */
		hiBarrierRates,        /* (I) High barrier rates        */
		koCurveRebates,        /* (I) Rebates                   */

		numericScalars,        /* (I) See #defines above.       */
		stringScalars,         /* (I) See #defines above.       */

		indexFreqs,            /* (I) # payments/year for indices */
		indexMaturities,       /* (I) Time (as TDateIntervals) till */
                                      /*     maturity for indices */
		indexDayCountStrs,     /* (I) Day count conventions */
		indexCurves,           /* (I) Which curve index is from */
		floatWeights,          /* (I) Weights for floating side  */
		knockWeights,          /* (I) Weights for knockout index */

		baseDates,
		discZeroDates,
		discZeroRates,
		indexZeroDates,
		indexZeroRates,

		volDates,
		vols,

		calibParams,
		output));


	pv = output[1];


	GtoErrMsg("%s: NPV= %12.8f\n", routine, pv);



	/* */
	if (DriTDrWrapperDataPutPrice(pv) != SUCCESS)
		goto done;


	/* */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed (compiled %s %s).\n",
			routine, __DATE__, __TIME__);
	return(status);
}



