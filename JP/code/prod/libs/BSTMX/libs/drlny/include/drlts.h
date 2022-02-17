/****************************************************************
 * Module:	DRL
 * Submodule:	TS - Curve Data Structure
 * File:	drlts.h
 * Function:	Curve routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_drlts_H
#define	_drlts_H
#include "drlstd.h"

#include <stdio.h>

#include "drlvtype.h"	/* Wrapper type definitions (TDateL, FLoatL, etc..) */

/*
 * Term Structure (zero coupons, base vol, etc ...)
 */
#if !defined(CLIB)
 
typedef struct
{
	TDate	fDate;
	double	fRate;
} TRatePt;

typedef struct
{
	int	fNumItems;	/* Number of TRatePts in fArray */
	TRatePt	*fArray;	/* dates & rates */
	TDate	fBaseDate;	/* the discount date */
	double	fBasis;		/* number compounding periods/year */
	long	fDayCountConv;	/* year fraction computation method */
} TCurve;
#else
# include "tcurve.h"	/* C analytics library */
#endif



#define	DrlTCurveBaseDate(object)		((object)->fBaseDate)
#define	DrlTCurveLastDate(object)		((object)->fArray[(object)->fNumItems-1].fDate)
#define	DrlTCurveBasis(object)		((object)->fBasis)
#define	DrlTCurveDayCountConv(object)	((object)->fDayCountConv)
#define	DrlTCurveNumItems(object)		((object)->fNumItems)
#define	DrlTCurveDate(object, index)	((object)->fArray[(index)].fDate)
#define	DrlTCurveRate(object, index)	((object)->fArray[(index)].fRate)


/*
 *
 */
extern	DLL_EXPORT(TDate) DrlTCurveRefDate(TCurve *that);
extern	DLL_EXPORT(int)	DrlTCurveDatesArray(TCurve *that,
				int *nDates, TDate **dates);

/*
 * I/O
 */

#define	DRL_TCURVE_FMT_STD	(1L)
#define	DRL_TCURVE_FMT_CUR	(2L)
#define	DRL_TCURVE_FMT_PERCENT	(3L)
#define	DRL_TCURVE_FMT_LON	(4L)
#define	DRL_TCURVE_FMT_REGTEST	(5L)


extern	DLL_EXPORT(int)	DrlTCurveFileRead(TCurve **that, char *, long fmt);
extern	DLL_EXPORT(int)	DrlTCurveFileWrite(TCurve *that, char *, long fmt);
extern	DLL_EXPORT(int)	DrlTCurveFpRead(TCurve **that, FILE *fp, long fmt);
extern	DLL_EXPORT(int) DrlTCurveFpWrite(TCurve *that, FILE *fp, long fmt);

extern	DLL_EXPORT(int)	DrlTCurveLondonBaseVolFpRead(TCurve **that, FILE *fp);
extern	DLL_EXPORT(int)	DrlTCurveLondonBaseVolFpRead_Mod(TCurve ***that, FILE *fp, int nbCurves);
extern	DLL_EXPORT(int)	DrlTCurveLondonBaseVolFileRead(TCurve **that, char *fnam);
extern	DLL_EXPORT(int)	DrlTCurveLondonBaseVolFileRead_Mod(TCurve ***that, char *fnam, int nbCurves);
extern	DLL_EXPORT(char*) DrlTCurvePrintYields(TCurve *that, char *s);

extern	DLL_EXPORT(int)	DrlTCurveWrapRead(TCurve **that, TDateL *refDate,
				TDateL *zDate, FloatL *zRate);
extern	DLL_EXPORT(int)	DrlTCurveWrapWrite(TCurve *that, TDateL *refDate,
				TDateL *zDate, FloatL *zRate);

/*
 *
 */
extern	DLL_EXPORT(int) DrlTCurveFwdRate(
	TCurve* that,		/* (I) zero curve */
	double tExp,		/* (I) time to reset interval */
	double tMat,		/* (I) forward maturity interval */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	TDayCount dayCountConv,	/* (I) day count convention */
	double *yield);		/* (O) forward rate */


extern	DLL_EXPORT(int)	DrlTCurveFwdSwapsSens(
	TCurve* that,		/* (I) zero curve */
	double tExp,		/* (I) time to reset interval */
	double tMat,		/* (I) forward maturity interval */
	int freq,		/* (I) rate frequency (0=simple, 1,2,4,12) */
	TDayCount dayCountConv,	/* (I) day count convention */
	char *what,		/* (I) "Y"=yield, "D"=duration */
	double *retVal);	/* (O) output fwd sensitivity */

extern	DLL_EXPORT(int)	DrlTCurveForwardDateTCurve(
	TCurve *that,		/* (I) zero curve */
	TDate resetDate,	/* (I) reset date */
	TCurve **outCurve);	/* (O) output zero curve */

extern	DLL_EXPORT(int)	DrlTCurveShiftRates(
	TCurve *that,		/* (B) zero curve */
	double amount);		/* (I) amount */

extern	DLL_EXPORT(int)	DrlTCurveShiftDates(
	TCurve *that,		/* (B) zero curve */
	TDate newBaseDate);	/* (I) new base date */

extern	DLL_EXPORT(TCurve*)	DrlTCurveInterpTCurve(
	TCurve *templateCurve,		/* (I) template */
	TCurve *valueCurve);		/* (I) to be interpolated */


/*
 *
 */
extern	DLL_EXPORT(int)	DrlTCurveBondSens(
	TCurve *zcCurve,
	TDate expDate,		/* (I) time to reset */
	TDate matDate,		/* (I) fwd maturity */
	int freq,		/* (I) yield frequency (1,2,4,12) */
	double *dur,		/* (O) regular duration */
	double *zDur,		/* (O) zCurve duration */
	double *zCvx,		/* (O) zCurve convexity */
	double *zP);		/* (O) zCurve cubicity */



#endif	/* _drlts_H */

