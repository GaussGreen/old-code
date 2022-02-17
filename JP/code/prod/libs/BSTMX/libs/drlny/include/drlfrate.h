/****************************************************************
 * Module:	DRL
 * Submodule:	FRATE - Rate Data Structure
 * File:	drlfrate.h
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_drlfrate_H
#define	_drlfrate_H
# include "drlstd.h"


#if !defined(CLIB)
/*t-@CDOC(idxn=TFloatRate,catn=structdef)
 * TFloatRate defines a floating rate. For example, a 5 year semi-annual
 * swap rate, LIBOR or a compounding rate.
 * <br><br>
 * Note that payInterval is ONLY used if rateType = GTO_SIMPLE_BASIS. 
 * This setting is used to compute LIBOR and par swap rates.
 * In order to compute compounding rates, set rateType = GTO_ANNUAL_BASIS (1)
 * or 2 for semi-annual compounding, etc.
 * See the documentation for GtoDayCountFraction for a list of possible 
 * day count convention constants. Note that there is no difference
 * between a TFloatRate which has a payInterval *equal* to the 
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
 * TDateAdjIntvl are used to adjust dates used to compute the rate. 
 */
typedef struct _TFloatRate
{
    TDateInterval   matInterval;        /* Time to maturity of rate */
    TDateInterval   payInterval;        /* Time between payments for rate */
    long            dayCountConv;       /* Day count convention of rate */
    TDateAdjIntvl   spotOffset;         /* From reset to rate effective date */
    double          spread;             /* Added to the rate  */
    long            rateType;           /* GTO_SIMPLE_BASIS, GTO_ANNUAL_BASIS*/
    double          weight;             /* Multiplied by rate */
} TFloatRate;                           
/*e*/
#else
# include "fltrate.h"
#endif	/*CLIB*/

/*
 * Routines
 */

extern	DLL_EXPORT(int)		DrlTFloatRateScan(char *s, TFloatRate *aRate);
extern	DLL_EXPORT(char*)	DrlTFloatRatePrint(char *s, TFloatRate *aRate);
extern	DLL_EXPORT(int)		DrlTFloatRateArrayFpWrite(TFloatRateArray *that,
					FILE *fp);


#ifdef	_SKIP
extern int GtoForwardRate(
    TCurve     *zeroCurve,     /* (I) Zero curve for generating forward rate */
    TFloatRate *fwdRateInfo,   /* (I) Information for the forward rate */
    TDate      startDate,      /* (I) Start date of forward */
    double     *forwardRate);  /* (O) Desired forward rate */
#endif




/*t-@CDOC(idxn=TFloatRateVolPoint,catn=structdef)
 * Structure describing a volatility point for a forward rate.
 */
typedef	struct {
	TFloatRate	fRate;		/* forward rate */
	TDateInterval	fReset,		/* rate reset */
			fExp,		/* option expiration */
			fStart;		/* option start */
	double		fValue;		/* volatility value */
} TFloatRateVolPoint;

/*e*/

extern	DLL_EXPORT(int)		DrlTFloatRateVolPointScan(char *s,
					TFloatRateVolPoint *that);
extern	DLL_EXPORT(char*)	DrlTFloatRateVolPointPrint(char *s,
					TFloatRateVolPoint *that);

/*t-@CDOC(idxn=TFloatRateCorrPoint,catn=structdef)
 * Correlation of forward rates point data structure
 */
typedef	struct {
	TFloatRate	fRate1,		/* rate 1 */
			fRate2;		/* rate 2 */
	TDateInterval	fExp,		/* option expiration */
			fReset1,	/* reset rate 1 */
			fReset2;	/* reset rate 2 */
	double		fValue;
} TFloatRateCorrPoint;

/*e*/

extern	DLL_EXPORT(int)		DrlTFloatRateCorrPointScan(char *s,
					TFloatRateCorrPoint *that);
extern	DLL_EXPORT(char*)	DrlTFloatRateCorrPointPrint(char *s,
					TFloatRateCorrPoint *that);




/*--------------------------------------------------------------
 * Two-Factor Model Benchmark Routines
 */

#define	DRL_TVOLBENCHMARK_VOL	((int) 0x01)
#define	DRL_TVOLBENCHMARK_CORR	((int) 0x02)


/*t-@CDOC(idxn=TVolBenchmark,catn=structdef)
 * Structure describing a volatility benchmark.
 */

typedef	struct	{
	int	fType;			/* type  */
	union {
	    TFloatRateVolPoint	fVol;
	    TFloatRateCorrPoint	fCorr;
	}	fU;
} TVolBenchmark;

/*e*/


extern	DLL_EXPORT(int)	DrlTVolBenchmarkIsValid(TVolBenchmark *pt);

extern	DLL_EXPORT(int)	DrlTVolBenchmarkSetVol(
	TVolBenchmark *that,		/* (O) benchmark set */
	TDateInterval rateReset,	/* (I) rate reset */
	TDateInterval rateMat,		/* (I) rate maturity */
	int rateFreq,			/* (I) rate frequency */
	TDayCount rateDcc);		/* (I) rate dcc (-1L for default) */

extern	DLL_EXPORT(int)	DrlTVolBenchmarkSetCorr(
	TVolBenchmark *that,		/* (O) benchmark set */
	TDateInterval rateReset,	/* (I) rate reset */
	TDateInterval rateMat1,		/* (I) rate maturity */
	int rateFreq1,			/* (I) rate frequency */
	TDayCount rateDcc1,		/* (I) rate dcc (-1L for default) */
	TDateInterval rateMat2,		/* (I) rate maturity */
	int rateFreq2,			/* (I) rate frequency */
	TDayCount rateDcc2);		/* (I) rate dcc (-1L for default) */

extern	DLL_EXPORT(int)	DrlTVolBenchmarkScan(char *s, TVolBenchmark *pt);
extern	DLL_EXPORT(char*) DrlTVolBenchmarkPrint(char *s, TVolBenchmark *pt);
extern	DLL_EXPORT(int)	DrlTVolBenchmarkFpRead(TVolBenchmark *pt, int *nPt,
			FILE *fp, ...);
extern	DLL_EXPORT(int)	DrlTVolBenchmarkFpWrite(TVolBenchmark *pt, int nPt,
			FILE *fp, ...);
extern	DLL_EXPORT(int)	DrlTVolBenchmarkFileRead(TVolBenchmark *pt, int *nPt,
			char *fnam, ...);
extern	DLL_EXPORT(int)	DrlTVolBenchmarkFileWrite(TVolBenchmark *pt, int nPt,
			char *fnam, ...);
extern	DLL_EXPORT(int)	DrlTVolBenchmarkWrapRead(TVolBenchmark *pt, int *nPt,
			IntL *nPtL, CharBlockL *ptL);
extern	DLL_EXPORT(int)	DrlTVolBenchmarkWrapWrite(TVolBenchmark *pt, int *nPt,
			IntL *nPtL, CharBlockL *ptL);



#endif /* _drlfrate_H */


