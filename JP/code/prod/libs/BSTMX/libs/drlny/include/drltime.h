/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	TIME - TDate Routines
 * File:	drltime.h
 * Function:	Date and time routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef _drltime_H
#define _drltime_H
#include "drlstd.h"

#include <stdarg.h>
#include <limits.h>

#if defined(CLIB)
#include "ldate.h"
#include "convert.h"
#include "cerror.h"
#endif

/*---------------------------------------------------------------
 * Definition of the type TDate and basic functions
 */
#if !defined(CLIB)

#define GTO_ACT_365		1L		/* Actual/365 */
#define GTO_ACT_365F		2L		/* Actual/365 Fixed */
#define GTO_ACT_360		3L		/* Actual/360 */
#define GTO_B30_360		4L		/* 30/360 */
#define GTO_B30E_360		5L		/* 30E/360 */
#define GTO_ACT_365FJ		6L		/*  */
#define GTO_B30E_360I		7L
#define GTO_ACT_ACT		GTO_ACT_365
#define GTO_B30_360_FIXED	8L		/* For bond coupon */
#define GTO_B30EP_360		9L		/* 30E+/360 */
#define GTO_ACT_ACT_FRF		10L		/* For French TAM  */

#endif /*!CLIB*/


#define DRL_DATE_NULL	((TDate)0)
#define DRL_DATE_MIN	((TDate)LONG_MIN)
#define DRL_DATE_MAX	((TDate)LONG_MAX)

#define	GTO_DCNT_ERROR	(100L)

/*
 * Type for day count convention.
 */
/*    typedef	long		TDayCount;*/

/*
 * Useful macros
 */

#define DRL_YEAR_IS_LEAP(year)				 	\
		(((year) % 4 == 0 && (year) % 100 != 0) || ((year) % 400 == 0))

#define	DRL_DATE_IS_LEAP(dt)						\
		YEAR_IS_LEAP(YEAR(dt))


/*---------------------------------------------------------------
 * Routines
 */

extern	DLL_EXPORT(TDate)	DrlTDateMake(int m, int d, int y);
extern	DLL_EXPORT(int)		DrlDAY(TDate dt);
extern	DLL_EXPORT(int)		DrlMONTH(TDate dt);
extern	DLL_EXPORT(int)		DrlYEAR(TDate dt);
extern	DLL_EXPORT(int)		DrlTDateSplit(TDate aTDate,
					int *mm, int *dd, int *yy);
extern	DLL_EXPORT(int)		DrlTDateCheckValid(TDate aTDate);

extern	DLL_EXPORT(int)		DrlTDateScan(char *buf, TDate *aTDate);
extern	DLL_EXPORT(char*)	DrlTDatePrint(char *buf, TDate aTDate);
extern	DLL_EXPORT(int)		DrlTDateScanYMD(char *buf, TDate *aTDate);
extern	DLL_EXPORT(char*)	DrlTDatePrintYMD(char *buf, TDate aTDate);
extern	DLL_EXPORT(int)		DrlTDayCountScan(char *buf, TDayCount *dcb);
extern	DLL_EXPORT(char*)	DrlTDayCountPrint(char *buf, TDayCount dcb);
extern	DLL_EXPORT(TDate)	DrlTDateToday(void);


extern	DLL_EXPORT(int)		DrlWeekDay(TDate date);
extern	DLL_EXPORT(char*)	DrlWeekDayPrint(char *string, int wDay);


extern	DLL_EXPORT(double)	DrlTDate2Time(TDate date, TDayCount dcb);
extern	DLL_EXPORT(TDate)	DrlTime2TDate(double tm, TDayCount dcb);
extern	DLL_EXPORT(long)	DrlTDate2Lotus(TDate aTDate);
extern	DLL_EXPORT(TDate)	DrlLotus2TDate(long n123);
extern	DLL_EXPORT(long)	DrlTDate2Julday(TDate aDate);
extern	DLL_EXPORT(TDate)	DrlJulday2TDate(long julday);
extern	DLL_EXPORT(TDate)	DrlDrDate2TDate(long drDate);
extern	DLL_EXPORT(long)	DrlTDate2DrDate(TDate tDate);

extern	DLL_EXPORT(int)		DrlTDayCountFromNumDen(TDayCount *daycount,
				char num, int den);
extern	DLL_EXPORT(int)		DrlTDayCountToNumDen(TDayCount daycount,
				char *num, int *den);

/*extern	DLL_EXPORT(TDate)	DrlTDateAddTime(TDate startDate, double tt);*/
extern	DLL_EXPORT(int)		DrlTDateAddYears(TDate startDate, double tt,
					TDate *endDate);

/*---------------------------------------------------------------
 * Operations on date arrays
 */
#ifdef	_SKIP
extern	DLL_EXPORT(TDate*)	DrlTDateArrayAlloc(int size);
extern	DLL_EXPORT(void)	DrlTDateArrayFree(TDate *);
extern	DLL_EXPORT(int)		DrlTDateArrayNewSimpleUniform(TDate startDate,
					int freq, double maturity,
					TDate **dates, int *nDates);
extern	DLL_EXPORT(int)		DrlLocateTDateInArray(int nTDateArray,
					TDate *dateArray,
					int *j, TDate aTDate);
extern	DLL_EXPORT(int)		DrlLocateClosestTDateInArray(int n,
					TDate *dateArray,
					int *j, TDate aTDate);
extern	DLL_EXPORT(int)		DrlSLUBTDateInArray(int n, TDate *dateArray,
					int *j, TDate aTDate);
extern	DLL_EXPORT(void)	DrlSortTDateArray(int nTDateArray,
					TDate *dateArray);
extern	DLL_EXPORT(void)	DrlSortNTDateArray(int *n, TDate *dateArray);
extern	DLL_EXPORT(int)		DrlAddTDateToArray(int *nTDateArray,
					TDate *dateArray, TDate aTDate);
extern	DLL_EXPORT(int)		DrlTDateArrayGetVaList(TDate minTDate,
					TDate maxTDate,
					int *nTDate, TDate **date,
					va_list ap);

extern	DLL_EXPORT(TDateList*)	DrlTDateArray2TDateList(int nTDates,
					TDate *dates);
#endif


/*---------------------------------------------------------------
 * Structure TDateInterval
 */

extern	DLL_EXPORT(TDateInterval)	DrlTDateIntervalFromDMY(int nDays,
						int nMonths, int nYears);


extern	DLL_EXPORT(int)	DrlTDateAdvanceYears(TDate baseDate, double years,
				TDayCount dayCount, TDate *newDate);
extern	DLL_EXPORT(int)	DrlTDateIntervalFromYears(TDate baseDate, double years,
				TDayCount dayCount, TDateInterval *interval);
extern	DLL_EXPORT(int)	DrlTDateIntervalToYears(TDate baseDate,
				TDateInterval interval, TDayCount dayCount,
				double *years);


extern	DLL_EXPORT(TDateInterval)	DrlTDateIntervalNil(void);
extern	DLL_EXPORT(int)		DrlTDateIntervalIsEqual(TDate refTDate,
					TDateInterval t1, TDateInterval t2);
extern	DLL_EXPORT(int)		DrlTDateIntervalIsLowerOrEqual(TDate refTDate,
					TDateInterval t1, TDateInterval t2);
extern	DLL_EXPORT(int)		DrlTDateIntervalIsGreaterOrEqual(TDate refTDate,
					TDateInterval t1, TDateInterval t2);
extern	DLL_EXPORT(int)		DrlTDateIntervalIsLower(TDate refTDate,
					TDateInterval t1, TDateInterval t2);
extern	DLL_EXPORT(int)		DrlTDateIntervalIsGreater(TDate refTDate,
					TDateInterval t1, TDateInterval t2);



extern	DLL_EXPORT(TDateInterval*)	DrlTDateIntervalArrayAlloc(int size);
extern	DLL_EXPORT(int)	DrlTDateIntervalArrayFree(TDateInterval* ptr);

extern	DLL_EXPORT(int)	DrlTDateIntervalScan(char *s, TDateInterval *tl);
extern	DLL_EXPORT(char*)	DrlTDateIntervalPrint(char *string,
					TDateInterval st);

extern	DLL_EXPORT(int)	DrlFreq2TDateInterval(int freq,
					TDateInterval *interval);
extern	DLL_EXPORT(int)	DrlTDateInterval2Freq(TDateInterval laps,
					int *freq);
extern	DLL_EXPORT(int)	DrlCheckFreqValid(int freq);



extern	DLL_EXPORT(int)	DrlTDateAdvanceToRefDate(TDate oldDate,
				int numNotifDays, TDateInterval interval,
				TDate refDate, TDate *newDate);



#endif	/* _drltime_H */


