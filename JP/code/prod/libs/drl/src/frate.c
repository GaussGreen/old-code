/****************************************************************
 * Module:	DRL
 * Submodule:	FRATE
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

/* C analytics library */
#ifdef	DRL_CLIB
#include "date_sup.h"
#include "zr2simp.h"
#include "zr2coup.h"
#include "duration.h"
#include "convex.h"
#include "swapadj.h"
#include "convert.h"		/* DrlDIntervalSet() */
#include "macros.h"		/* WRAP_STR_IDX */
#endif

#include "drltime.h"
#include "drlstr.h"
#include "drlio.h"
#include "drlts.h"

#include "drlfrate.h"		/* prototype consistency */

/*f-------------------------------------------------------------
 * DFloatRate : computes a forward rate on a zero curve.
 */

int
DrlDFloatRateForward(
	DFloatRate *that,	/* (I) input rate */
        DDate startDate,	/* (I) reset date */
        DCurve* zcCurve,	/* (I) zero curve */
	double *rate)		/* (O) forward rate */
{
static	char	routine[] = "DrlDFloatRateForward";
	int	status = FAILURE;
	int	freq;

	IF_FAILED_DONE( DrlDIntervalToFreq(
		&that->payInterval,
		&freq));
	IF_FAILED_DONE( DrlDCurveForwardRate2(
		zcCurve,
		startDate,
		that->matInterval,
		freq,
		that->dayCountConv,
		rate));

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*--------------------------------------------------------------
 * 
 */

int
DrlDFloatRateIsValid(DFloatRate *that)
{
	return(TRUE);
}


/*f-------------------------------------------------------------
 * Reads a <i> DFloatRate</i> in a string <i> s</i> and puts the result
 * (if successful) in <i> that</i>.
 * The recognized string format is 
 * \begin{verbatim}
 * FLRATE(Mat=<value>,Freq=<value>,DayCount=<value>,
 *        Spread=<value>,Weight=<value>)
 * \end{verbatim}
 * Returns 0 iff successful.
 */


int
DrlDFloatRateScan(char *s, DFloatRate *that)
{
/*static	char	routine[] = "DrlDFloatRateScan";*/
	int	status = FAILURE;
	double	spreadDef = 0e0,
		weightDef = 1e0;

	if (DrlStrParseFuncScan(s, "FLRATE", 5,
		"Mat",  DRL_TDATEINTERVAL_T, (void*) &that->matInterval,
			0, (void*) NULL,
		"Freq", DRL_TDATEINTERVAL_T, (void*) &that->payInterval,
			0, (void*) NULL,
		"DayCount", DRL_TDAYCOUNT_T, (void*) &that->dayCountConv,
			0, (void*) NULL,
		"Spread", DRL_DOUBLE_T, (void*) &that->spread,
			1, (void*) &spreadDef,
		"Weight", DRL_DOUBLE_T, (void*) &that->weight,
			1, (void*) &weightDef
		) != SUCCESS)
			goto done;
#ifdef	DRL_CLIB
	GTO_SET_ADJ_INTERVAL_DAYS(that->spotOffset, 0);
#endif

	status = SUCCESS;
done:
	return(status);
}

/*f-------------------------------------------------------------
 * Prints a <i> DFloatRate</i> <i> that</i> in a string <i> s</i>.
 * The output format is the same as in <i> DrlDFloatRateScan</i>.
 * Returns <i> s</i> (or a private copy if <i> s</i> is NULL).
 */

char*
DrlDFloatRatePrint(char *s, DFloatRate *that)
{
static	char	buf[256];
	s =  (s == NULL ? buf : s);


	DrlStrParseFuncPrint(s, "RATE", 5,
		"Mat",  DRL_TDATEINTERVAL_T, (void*) &that->matInterval,
			0, (void*) NULL,
		"Freq", DRL_TDATEINTERVAL_T, (void*) &that->payInterval,
			0, (void*) NULL,
		"DayCount", DRL_TDAYCOUNT_T, (void*) &that->dayCountConv,
			0, (void*) NULL,
		"Spread", DRL_DOUBLE_T, (void*) &that->spread,
			0, (void*) NULL,
		"Weight", DRL_DOUBLE_T, (void*) &that->weight,
			0, (void*) NULL
		);

	return(s);
}

#ifdef	_SKIP
/*---------------------------------------------------------------
 * Prints a <i> DFloatRateArray</i> <i> that</i> in a string <i> s</i>.
 * The output format is the same as in <i> DrlDFloatRateScan</i>.
 * Returns <i> s</i> (or a private copy if <i> s</i> is NULL).
 */

int
DrlDFloatRateArrayFpWrite(DFloatRateArray *that, FILE *fp)
{
	int	idx;

	DrlFPrintf(fp, "# NUM_RATES\n%d\n", that->numRates);
	DrlFPrintf(fp, "# IDX  CURVE_IDX    TFLOATRATE\n");
	for (idx=0; idx<= that->numRates-1; idx++) {
	    DrlFPrintf(fp, "  %2d   %2d    %s\n",
		idx,
		that->curveIndices[idx],
		DrlDFloatRatePrint(NULL, &that->defs[idx])
	    );
	}
	return(SUCCESS);
}
#endif


/*f--------------------------------------------------------------
 * Scans a swaption volatility point in a char string "s" and puts
 * the result in "that". The recognized format for the string
 * is
 * \begin{verbatim}
 * FLRATE_VOL(Start=<interval>,Exp=<interval>,Reset=<interval>,
 *            Mat=<interval>,Freq=<int>)
 * \end{verbatim}
 * Returns 0 iff successful.
 */

int
DrlDFloatRateVolPointScan(char *s, DFloatRateVolPoint *that)
{
/*static	char	routine[] = "DrlDFloatRateVolPointScan";*/
	int	freq, freqDef = 2,
		status = FAILURE;
	DInterval	startDef;

	if (DrlDIntervalSet(0, 'D', &startDef) != SUCCESS)
		goto done;

	status = DrlStrParseFuncScan(s, "FLRATE_VOL", 5,
	    "Exp",  DRL_TDATEINTERVAL_T, (void*) &that->fExp,
			0, (void*) NULL,
	    "Reset",DRL_TDATEINTERVAL_T, (void*) &that->fReset,
			0, (void*) NULL,
	    "Mat",  DRL_TDATEINTERVAL_T, (void*) &that->fRate.matInterval,
			0, (void*) NULL,
	    "Start",DRL_TDATEINTERVAL_T, (void*) &that->fStart,
			1, (void*) &startDef,
	    "Freq", DRL_INT_T,           (void*) &freq,
			1, (void*) &freqDef);
	if (status != SUCCESS) goto done;

	if (freq == 0) {
	    that->fRate.payInterval = that->fRate.matInterval;
	} else {
	    if (DrlFreqToDInterval((long)freq, &that->fRate.payInterval)
		!= SUCCESS) goto done;
	}
	that->fRate.spread = 0e0;
	that->fRate.weight = 0e0;
#ifdef	DRL_CLIB
	GTO_SET_ADJ_INTERVAL_DAYS(that->fRate.spotOffset, 0);
#endif

	/* set proper day count */
	that->fRate.dayCountConv = (freq > 0 ? DRL_B30_360 : DRL_ACT_360);

	/* made it through */
	status = SUCCESS;
done:
	/*if (status != SUCCESS) {
	    DrlErrMsg("%s: failed (can't read string `%s')\n", routine, s);
	}*/
	return(status);
}


/*f--------------------------------------------------------------
 * Prints a swaption volatility point "that" in a char string "s".
 * The output format is the same as in <i> DrlDFloatRateVolPointScan</i>.
 */

char*
DrlDFloatRateVolPointPrint(char *s, DFloatRateVolPoint *that)
{
static	char	buf[256];
	int	freq;
	s =  (s == NULL ? buf : s);

	if (DrlDIntervalToFreq(&that->fRate.payInterval, &freq)
	    != SUCCESS) {
		sprintf(s, "FLRATE_CORR=ERROR");
		return(s);
	}

	DrlStrParseFuncPrint(s, "FLRATE_VOL", 5,
	    "Start",DRL_TDATEINTERVAL_T, (void*) &that->fStart,
			0, (void*) NULL,
	    "Exp",  DRL_TDATEINTERVAL_T, (void*) &that->fExp,
			0, (void*) NULL,
	    "Reset",DRL_TDATEINTERVAL_T, (void*) &that->fReset,
			0, (void*) NULL,
	    "Mat",  DRL_TDATEINTERVAL_T, (void*) &that->fRate.matInterval,
			0, (void*) NULL,
	    "Freq", DRL_INT_T,           (void*) &freq,
			0, (void*) NULL);

	return(s);
}



/*f--------------------------------------------------------------
 * Scans a correlation point in a char string "s" and puts
 * the result in "that". The recognized format for the string
 * is
 * \begin{verbatim}
 * FLRATE_CORR(Exp=<interval>,Mat1=<interval>,Mat2=<interval>,
 *            Freq1=<int>,Freq2=<int)
 * \end{verbatim}
 * Returns 0 if successful.
 */

int
DrlDFloatRateCorrPointScan(char *s, DFloatRateCorrPoint *that)
{
/*static	char	routine[] = "DrlDFloatRateCorrPointScan";*/
	int	freq1, freq2, freqDef = 2;
	int	status = FAILURE;

	status = DrlStrParseFuncScan(s, "FLRATE_CORR", 7,
	    "Exp",    DRL_TDATEINTERVAL_T, &that->fExp, 0, NULL,
	    "Reset1", DRL_TDATEINTERVAL_T, &that->fReset1, 0, NULL,
	    "Reset2", DRL_TDATEINTERVAL_T, &that->fReset2, 0, NULL,
	    "Mat1",   DRL_TDATEINTERVAL_T, &that->fRate1.matInterval, 0, NULL,
	    "Mat2",   DRL_TDATEINTERVAL_T, &that->fRate2.matInterval, 0, NULL,
	    "Freq1",  DRL_INT_T,           &freq1, 1, &freqDef,
	    "Freq2",  DRL_INT_T,           &freq2, 1, &freqDef
	);
	if (status != SUCCESS) goto done;

	if (freq1 == 0) {
	    that->fRate1.payInterval = that->fRate1.matInterval;
	} else {
	    if (DrlFreqToDInterval((long)freq1, &that->fRate1.payInterval)
		!= SUCCESS) goto done;
	}
	that->fRate1.spread = 0e0;
	that->fRate1.weight = 0e0;
#ifdef	DRL_CLIB
	GTO_SET_ADJ_INTERVAL_DAYS(that->fRate1.spotOffset, 0);
#endif

	if (freq2 == 0) {
	    that->fRate2.payInterval = that->fRate2.matInterval;
	} else {
	    if (DrlFreqToDInterval((long)freq2, &that->fRate2.payInterval)
		!= SUCCESS) goto done;
	}
	that->fRate2.spread = 0e0;
	that->fRate2.weight = 0e0;
#ifdef	DRL_CLIB
	GTO_SET_ADJ_INTERVAL_DAYS(that->fRate2.spotOffset, 0);
#endif

	/* set proper day count */
	that->fRate1.dayCountConv = (freq1 > 0 ? DRL_B30_360 : DRL_ACT_360);
	that->fRate2.dayCountConv = (freq2 > 0 ? DRL_B30_360 : DRL_ACT_360);


	/* made it through */
	status = SUCCESS;
done:
	/*if (status != SUCCESS) {
	    DrlErrMsg("%s: failed (can't read string `%s')\n", routine, s);
	}*/
	return(status);
}


/*f--------------------------------------------------------------
 * Prints a correlation point "that" in a char string "s".
 * The output format is the same as in <i> DrlDFloatRateCorrPointScan</i>.
 */

char*
DrlDFloatRateCorrPointPrint(char *s, DFloatRateCorrPoint *that)
{
static	char	buf[256];
	int	freq1, freq2;
	s =  (s == NULL ? buf : s);

	if ((DrlDIntervalToFreq(&that->fRate1.payInterval, &freq1)
		!= SUCCESS) || 
	    (DrlDIntervalToFreq(&that->fRate2.payInterval, &freq2)
		!= SUCCESS)) {
		sprintf(s, "FLRATE_CORR=ERROR");
		return(s);
	}
	
	DrlStrParseFuncPrint(s, "FLRATE_CORR", 7,
	    "Exp",    DRL_TDATEINTERVAL_T, &that->fExp,  0, NULL,
	    "Reset1", DRL_TDATEINTERVAL_T, &that->fReset1,  0, NULL,
	    "Reset2", DRL_TDATEINTERVAL_T, &that->fReset2,  0, NULL,
	    "Mat1",   DRL_TDATEINTERVAL_T, &that->fRate1.matInterval, 0, NULL,
	    "Mat2",   DRL_TDATEINTERVAL_T, &that->fRate2.matInterval, 0, NULL,
	    "Freq1",  DRL_INT_T,           &freq1,       0, NULL,
	    "Freq2",  DRL_INT_T,           &freq2,       0, NULL
	);

	return(s);
}


/*f-------------------------------------------------------------
 * Convenience routine to construct set a DVolBenchmark
 * to be a cap/swaption volatility point.
 * Returns 0 iff successful.
 */

int
DrlDVolBenchmarkSetVol(
	DVolBenchmark *that,		/* (O) benchmark set */
	DInterval rateReset,	/* (I) rate reset */
	DInterval rateMat,		/* (I) rate maturity */
	int rateFreq,			/* (I) rate frequency */
	DDayCount rateDcc)		/* (I) rate dcc (-1L for default) */
{
	int	status = FAILURE;

	DFloatRateVolPoint	*vpt;
	DInterval		startDef;

	that->fType = DRL_TVOLBENCHMARK_VOL;
	vpt  = &(that->fU.fVol);

	IF_FAILED_DONE( DrlDIntervalSet(0, 'D', &startDef));

	vpt->fExp = rateReset;
	vpt->fReset = rateReset;
	vpt->fRate.matInterval = rateMat;
	vpt->fStart = startDef;
	if (rateFreq == 0) {
	    vpt->fRate.payInterval = vpt->fRate.matInterval;
	} else {
	    IF_FAILED_DONE( DrlFreqToDInterval(
		(long)rateFreq,
		&vpt->fRate.payInterval));
	}

	vpt->fRate.spread = 0e0;
	vpt->fRate.weight = 0e0;
#ifdef	DRL_CLIB
	GTO_SET_ADJ_INTERVAL_DAYS(vpt->fRate.spotOffset, 0);
#endif

	vpt->fRate.dayCountConv = (rateDcc == -1L ? 
			(rateFreq > 0 ? DRL_B30_360 : DRL_ACT_360) :
			rateDcc);

	/* made it through */
	status = SUCCESS;
done:
	/*if (status != SUCCESS) {
	    DrlErrMsg("%s: failed.\n", routine, s);
	}*/
	return(status);
}


/*f-------------------------------------------------------------
 * Convenience routine to construct set a DVolBenchmark
 * to be a corrleation point.
 * (the default day count is assumed to be ACT/360 for simple rates
 * and 30/360 for coupon rates).
 * Returns 0 iff successful.
 */

int
DrlDVolBenchmarkSetCorr(
	DVolBenchmark *that,		/* (O) benchmark set */
	DInterval rateReset,	/* (I) rate reset */
	DInterval rateMat1,		/* (I) rate maturity */
	int rateFreq1,			/* (I) rate frequency */
	DDayCount rateDcc1,		/* (I) rate dcc (-1L for default) */
	DInterval rateMat2,		/* (I) rate maturity */
	int rateFreq2,			/* (I) rate frequency */
	DDayCount rateDcc2)		/* (I) rate dcc (-1L for default) */
{
	int	status = FAILURE;

	DFloatRateCorrPoint	*cpt;
	DInterval		startDef;

	that->fType = DRL_TVOLBENCHMARK_CORR;
	cpt  = &(that->fU.fCorr);

	IF_FAILED_DONE( DrlDIntervalSet(0, 'D', &startDef));

	cpt->fExp = rateReset;
	cpt->fReset1 = rateReset;
	cpt->fReset2 = rateReset;
	cpt->fRate1.matInterval = rateMat1;
	cpt->fRate2.matInterval = rateMat2;


	if (rateFreq1 == 0) {
	    cpt->fRate1.payInterval = cpt->fRate1.matInterval;
	} else {
	    IF_FAILED_DONE( DrlFreqToDInterval(
		(long)rateFreq1,
		&cpt->fRate1.payInterval));
	}
	cpt->fRate1.spread = 0e0;
	cpt->fRate1.weight = 0e0;
#ifdef	DRL_CLIB
	GTO_SET_ADJ_INTERVAL_DAYS(cpt->fRate1.spotOffset, 0);
#endif

	if (rateFreq2 == 0) {
	    cpt->fRate2.payInterval = cpt->fRate2.matInterval;
	} else {
	    IF_FAILED_DONE( DrlFreqToDInterval(
		(long)rateFreq2,
		&cpt->fRate2.payInterval));
	}
	cpt->fRate2.spread = 0e0;
	cpt->fRate2.weight = 0e0;
#ifdef	DRL_CLIB
	GTO_SET_ADJ_INTERVAL_DAYS(cpt->fRate2.spotOffset, 0);
#endif

	/* set proper day count */
	cpt->fRate1.dayCountConv = (rateDcc1 == -1L ? 
			(rateFreq1 > 0 ? DRL_B30_360 : DRL_ACT_360) :
			rateDcc1);
	cpt->fRate2.dayCountConv = (rateDcc2 == -1L ? 
			(rateFreq2 > 0 ? DRL_B30_360 : DRL_ACT_360) :
			rateDcc2);

	/* made it through */
	status = SUCCESS;
done:
	/*if (status != SUCCESS) {
	    DrlErrMsg("%s: failed.\n", routine, s);
	}*/
	return(status);
}


/*f-------------------------------------------------------------
 * Returns TRUE/FALSE.
 */

int
DrlDVolBenchmarkIsValid(DVolBenchmark *pt)
{
static 	char	routine[] = "DrlDVolBenchmarkIsValid";
	switch(pt->fType) {
	case DRL_TVOLBENCHMARK_VOL:
	case DRL_TVOLBENCHMARK_CORR:
		return TRUE;
	default:
		DrlErrMsg("%s: bad type `%d'\n", routine, pt->fType);
		return FALSE;
	}
}


/*f-------------------------------------------------------------
 * Read a volatility benchmark
 * in a string "s" and places the result in "pt".\\
 * Currently recognized formats are
 * <br>
 * <br> <i> Swaption Volatility Points.</i>
 * \begin{verbatim}
 * FLRATE_VOL(Exp=3M,Reset=2Y,Mat=5Y,Freq=2)
 * \end{verbatim}
 * <br> <i> Average Correlation Points.</i>
 * \begin{verbatim}
 * FLRATE_CORR(Exp=3M,Reset1=6M,Mat1=1Q,Reset2=9M,Mat2=10Y)
 * \end{verbatim}
 * <br>
 * Returns 0 if scan successful.
 */

int
DrlDVolBenchmarkScan(char *s, DVolBenchmark *pt)
{
	int	errCode;

	errCode = DrlDFloatRateVolPointScan(s, &pt->fU.fVol);
	if (errCode == 0) {
		pt->fType = DRL_TVOLBENCHMARK_VOL;
		return(0);
	}
	if (errCode != -1) goto done;


	errCode = DrlDFloatRateCorrPointScan(s, &pt->fU.fCorr);
	if (errCode == 0) {
		pt->fType = DRL_TVOLBENCHMARK_CORR;
		return(0);
	}
	if (errCode != -1) goto done;


	errCode = 2;
done:
	if (errCode != 0) {
	    DrlErrMsg("DrlDVolBenchmarkScan: can't parse `%s'\n", s);
	}

	return(errCode);
}


/*f-------------------------------------------------------------
 * Prints a volatility benchmark "pt"
 * in the string "s"
 * (same format as in <i> DrlDVolBenchmarkScan</i>).
 * Returns a pointer to "s".
 */

char*
DrlDVolBenchmarkPrint(char *s, DVolBenchmark *pt)
{
static 	char	routine[] = "DrlDVolBenchmarkPrint";
	switch(pt->fType) {
	case DRL_TVOLBENCHMARK_VOL:
		return DrlDFloatRateVolPointPrint(s, &pt->fU.fVol);
	case DRL_TVOLBENCHMARK_CORR:
		return DrlDFloatRateCorrPointPrint(s, &pt->fU.fCorr);
	default:
		DrlErrMsg("%s: bad type `%c'\n", routine, pt->fType);
		return("ERROR()");
	}
}


/*--------------------------------------------------------------
 * Stream I/O for an array of DVolBenchmark
 */

static	int
_vDVolBenchmarkFpRead(DVolBenchmark *pt, int *nPt, FILE *fp,
		va_list arg)
{
#ifdef	_SKIP
	int	line=0, i, j, jmax;
	double	*x[32];
	char	buf[256], *s0, *s1, sTmp0[256], *sTmp1;
static	char	toks[] = " \t;\n";


	for(j=0; (x[j] = (double*) va_arg(arg, double*)) != NULL; j++);
	jmax = j;

	if (GetNextVariableFp(DRL_INT_T, (void*) nPt, fp, &line) != 0)
	return(line);

	for (i=0; i<=*nPt-1; i++) {
		if (FGetLine(buf, sizeof(buf), fp, &line) == NULL)
			return(line);
		s0 = buf;

		if ((s1 = StrToken(s0, toks, sTmp0, &sTmp1)) == NULL)
			return(line);
		s0 = NULL;
		if (DVolBenchmarkScan(s1, pt+i) != 0) return(line);

		for(j=0; j<=jmax-1; j++) {
			if ((s1 = StrToken(s0, toks, sTmp0, &sTmp1)) == NULL)
				return(line);
			if (Sscanf(s1, "%lf", &x[j][i]) != 1)
				return(line);
		}
	}
	return(0);
#endif
	DrlErrMsg("_vDVolBenchmarkFpRead: not implemented\n");
	return(1);
}

static	int
_vDVolBenchmarkFpWrite(DVolBenchmark *pt, int nPt, FILE *fp,
		va_list arg)
{
	int	i, j, jmax;
	double	*x[32];

	for(j=0; (x[j] = (double*) va_arg(arg, double*)) != NULL; j++);
	jmax = j;

	DrlFPrintf(fp, "#\n# Vtfm Voltility Benchmarks\n#\n");
	DrlFPrintf(fp, "%d\n#\n", nPt);

	for (i=0; i<=nPt-1; i++) {
		DrlFPrintf(fp, "\t%s", DrlDVolBenchmarkPrint(NULL, pt+i));

		for(j=0; j<=jmax-1; j++) {
			DrlFPrintf(fp, "\t%lf", x[j][i]);
		}
		DrlFPrintf(fp, "\n");
	}

	DrlFPrintf(fp, "#\n");

	return(0);
}

/*--------------------------------------------------------------
 * Stream I/O for an array of DVolBenchmark
 */

int
DrlDVolBenchmarkFpRead(DVolBenchmark *pt, int *nPt, FILE *fp, ...)
{
	int	errCode;
	va_list	ap;

	va_start(ap, fp);
	errCode = _vDVolBenchmarkFpRead(pt, nPt, fp, ap);
	va_end(ap);
	return(errCode);
}


int
DrlDVolBenchmarkFpWrite(DVolBenchmark *pt, int nPt, FILE *fp, ...)
{
	int	errCode;
	va_list	ap;

	va_start(ap, fp);
	errCode = _vDVolBenchmarkFpWrite(pt, nPt, fp, ap);
	va_end(ap);
	return(errCode);
}


/*--------------------------------------------------------------
 * File I/O for an array of DVolBenchmark
 */

int
DrlDVolBenchmarkFileRead(DVolBenchmark *pt, int *nPt, char *fnam, ...)
{
	FILE	*fp ;
	int	errCode;
	va_list	ap;

	va_start(ap, fnam);
	if ((fp = fopen(fnam, "r")) == NULL) {
		errCode = -1;
	} else {
		errCode = _vDVolBenchmarkFpRead(pt, nPt, fp, ap);
		fclose(fp) ;
	}
	va_end(ap);
	return(errCode);
}


int
DrlDVolBenchmarkFileWrite(DVolBenchmark *pt, int nPt, char *fnam, ...)
{
	FILE	*fp ;
	int	errCode ;
	va_list	ap;

	va_start(ap, fnam);
	if ((fp = fopen(fnam, "w")) == NULL) {
		errCode = -1;
	} else {
		errCode = _vDVolBenchmarkFpWrite(pt, nPt, fp, ap);
		fclose(fp) ;
	}
	va_end(ap);
	return(errCode);
}


/*f-------------------------------------------------------------
 * Reads an array of volatility benchmarks from a wrapper interface.
 * On entry,
 * "nPtL" is an array of int length 1 containing the number of benchmarks
 * to be used,
 * "ptL" is an array of chararcter strings containing the description
 * of the benchmarks (same as in <i> DrlDVolBenchmarkScan</i>).
 * On exit, "nPt" is the number of items and "pt" is the array
 * of read benchmarks (should be long enough).
 * Returns 0 iff OK.
 */


int
DrlDVolBenchmarkWrapRead(
	DVolBenchmark *pt,	/* (O) array of benchmarks */
	int *nPt,		/* (O) actual # of benchmarks */
	long *nPtL,		/* (I) LIL vector (or NULL) */
	char *ptL)		/* (I) LIL vector */
{
static	char	routine[] = "DrlDVolBenchmarkWrapRead";
	int	i, status = FAILURE;

	/* get # of arguments */
	DRL_WRAP_CHECK_VECTOR(ptL);
	*nPt = DRL_ARGSIZE(ptL);

	if (nPtL != NULL) {
	    DRL_WRAP_CHECK_VECTOR(nPtL);
	    *nPt = MIN(*nPt, (int) nPtL[1]);
	}

	/* scan strings */
	for (i=0; i<=*nPt-1; i++) {
	    if (DrlDVolBenchmarkScan(&ptL[DRL_WRAP_STR_IDX(i+1)],
		&pt[i]) != 0) {
		goto done;
	    }
	    if (!DrlDVolBenchmarkIsValid(&pt[i])) {
		goto done;
	    }
	}

	/* made it through OK */
#ifdef	__DEBUG__
	DrlDVolBenchmarkFpWrite(pt, *nPt, stdout, NULL);
#endif
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed\n", routine);
	}
	return(status);
}


