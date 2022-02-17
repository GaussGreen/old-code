/****************************************************************
 * Module:	DRL
 * Submodule:	FRATE - Rate Data Structure
 * File:	drlfrate.h
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_drlfrate_H
#define	_drlfrate_H

#include "drlstd.h"
#include "drlts.h"

#include <stdio.h>


#if !defined(DRL_CLIB)
/*t-@CDOC(idxn=DFloatRate,catn=structdef)
 * DFloatRate defines a floating rate. For example, a 5 year semi-annual
 * swap rate, LIBOR or a compounding rate.
 * <br><br>
 * Note that payInterval is ONLY used if rateType = GTO_SIMPLE_BASIS. 
 * This setting is used to compute LIBOR and par swap rates.
 * In order to compute compounding rates, set rateType = GTO_ANNUAL_BASIS (1)
 * or 2 for semi-annual compounding, etc.
 * See the documentation for GtoDayCountFraction for a list of possible 
 * day count convention constants. Note that there is no difference
 * between a DFloatRate which has a payInterval *equal* to the 
 * matInterval (for the same instrument) one which has a payInterval
 * which is *greater* than the matInterval, and one which has payInterval = 0.
 * This is because the instruments *always* have a payment at maturity, 
 * and also payments at intervals of payInterval, for dates on or 
 * before the maturity date. 
 * Simple LIBOR, etc (a zero coupon rate) is expressed by making
 * the payInterval = matInterval. 
 * Note that the spread is added to the rate *AFTER* multiplying by
 * the weight. 
 * Note that the holiday file and bad day convention in the
 * DDateAdjIntvl are used to adjust dates used to compute the rate. 
 */
typedef struct _DFloatRate
{
    DInterval   matInterval;        /* Time to maturity of rate */
    DInterval   payInterval;        /* Time between payments for rate */
    long        dayCountConv;       /* Day count convention of rate */
    DInterval   spotOffset;         /* From reset to rate effective date */
    double      spread;             /* Added to the rate  */
    long        rateType;           /* GTO_SIMPLE_BASIS, GTO_ANNUAL_BASIS*/
    double      weight;             /* Multiplied by rate */
} DFloatRate;                           
/*e*/
#else
# include "fltrate.h"
typedef	TFloatRate	DFloatRate;
#endif	/*DRL_CLIB*/

/*
 * Routines
 */

extern	int	DrlDFloatRateScan(char *s, DFloatRate *aRate);
extern	char*	DrlDFloatRatePrint(char *s, DFloatRate *aRate);

extern	int	DrlDFloatRateForward(
			DFloatRate *that,	/* (I) input rate */
        		DDate startDate,	/* (I) reset date */
        		DCurve* zcCurve,	/* (I) zero curve */
			double *rate);		/* (O) forward rate */


/*t-@CDOC(idxn=DFloatRateVolPoint,catn=structdef)
 * Structure describing a volatility point for a forward rate.
 */
typedef	struct {
	DFloatRate	fRate;		/* forward rate */
	DInterval	fReset,		/* rate reset */
			fExp,		/* option expiration */
			fStart;		/* option start */
	double		fValue;		/* volatility value */
} DFloatRateVolPoint;

/*e*/

extern	int	DrlDFloatRateVolPointScan(char *s, DFloatRateVolPoint *that);
extern	char*	DrlDFloatRateVolPointPrint(char *s, DFloatRateVolPoint *that);


/*t-@CDOC(idxn=DFloatRateCorrPoint,catn=structdef)
 * Correlation of forward rates point data structure
 */
typedef	struct {
	DFloatRate	fRate1,		/* rate 1 */
			fRate2;		/* rate 2 */
	DInterval	fExp,		/* option expiration */
			fReset1,	/* reset rate 1 */
			fReset2;	/* reset rate 2 */
	double		fValue;
} DFloatRateCorrPoint;

/*e*/

extern	int	DrlDFloatRateCorrPointScan(char *s, DFloatRateCorrPoint *that);
extern	char*	DrlDFloatRateCorrPointPrint(char *s, DFloatRateCorrPoint *that);



/*--------------------------------------------------------------
 * Two-Factor Model Benchmark Routines
 */

#define	DRL_TVOLBENCHMARK_VOL	((int) 0x01)
#define	DRL_TVOLBENCHMARK_CORR	((int) 0x02)


/*t-@CDOC(idxn=DVolBenchmark,catn=structdef)
 * Structure describing a volatility benchmark.
 */

typedef	struct	{
	int	fType;			/* type  */
	union {
	    DFloatRateVolPoint	fVol;
	    DFloatRateCorrPoint	fCorr;
	}	fU;
} DVolBenchmark;

/*e*/


extern	int	DrlDVolBenchmarkIsValid(DVolBenchmark *pt);

extern	int	DrlDVolBenchmarkSetVol(
	DVolBenchmark *that,		/* (O) benchmark set */
	DInterval rateReset,	/* (I) rate reset */
	DInterval rateMat,		/* (I) rate maturity */
	int rateFreq,			/* (I) rate frequency */
	DDayCount rateDcc);		/* (I) rate dcc (-1L for default) */

extern	int	DrlDVolBenchmarkSetCorr(
	DVolBenchmark *that,		/* (O) benchmark set */
	DInterval rateReset,	/* (I) rate reset */
	DInterval rateMat1,		/* (I) rate maturity */
	int rateFreq1,			/* (I) rate frequency */
	DDayCount rateDcc1,		/* (I) rate dcc (-1L for default) */
	DInterval rateMat2,		/* (I) rate maturity */
	int rateFreq2,			/* (I) rate frequency */
	DDayCount rateDcc2);		/* (I) rate dcc (-1L for default) */

extern	int	DrlDVolBenchmarkScan(char *s, DVolBenchmark *pt);
extern	char*	DrlDVolBenchmarkPrint(char *s, DVolBenchmark *pt);
extern	int	DrlDVolBenchmarkFpRead(DVolBenchmark *pt, int *nPt,
			FILE *fp, ...);
extern	int	DrlDVolBenchmarkFpWrite(DVolBenchmark *pt, int nPt,
			FILE *fp, ...);
extern	int	DrlDVolBenchmarkFileRead(DVolBenchmark *pt, int *nPt,
			char *fnam, ...);
extern	int	DrlDVolBenchmarkFileWrite(DVolBenchmark *pt, int nPt,
			char *fnam, ...);
extern	int	DrlDVolBenchmarkWrapRead(DVolBenchmark *pt, int *nPt,
			long *nPtL, char *ptL);
extern	int	DrlDVolBenchmarkWrapWrite(DVolBenchmark *pt, int *nPt,
			long *nPtL, char *ptL);



#endif /* _drlfrate_H */


