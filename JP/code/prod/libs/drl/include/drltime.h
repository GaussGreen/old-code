/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	TIME - DDate Routines
 * File:	drltime.h
 * Function:	Date and time routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef _drltime_H
#define _drltime_H
#include "drlstd.h"

#if defined(DRL_CLIB)
  /* Using ALIB */
# include "ldate.h"
# include "convert.h"
# include "cerror.h"
# include "date_sup.h"
#endif

/*
 * Routines
 */

extern	int	DrlDDateMake(int m, int d, int y, DDate *date);
extern	int	DrlDDateSplit(DDate aDDate, int *mm, int *dd, int *yy);
extern	int	DrlDDateCheckValid( DDate aDDate);
extern	int	DrlDDateScan(char *buf, DDate *aDDate);
extern	char*	DrlDDatePrint(char *buf, DDate aDDate);
extern	int	DrlDDateScanYMD(char *buf, DDate *aDDate);
extern	char*	DrlDDatePrintYMD(char *buf, DDate aDDate);
extern	int	DrlDDayCountScan(char *buf, DDayCount *dcb);
extern	char*	DrlDDayCountPrint(char *buf, DDayCount dcb);
extern	int	DrlDDateToday(DDate *today);
extern	int	DrlWeekDay(DDate date, int *dayno);
extern	char*	DrlWeekDayPrint(char *string, int wDay);

extern	int	DrlDDateToTime(
			DDate date, DDayCount dcb, double *yearfrac);
extern	int	DrlTimeToDDate(
			double tm, DDayCount dcb, DDate *date);
extern	int	DrlDDateToJulday(
			DDate aDate, long *julday);
extern	int	DrlJuldayToDDate(
			long julday, DDate *date);

extern	int	DrlAlibDateToDrDate(
			long alibDate, DDate *drDate);
extern	int	DrlAlibDateToDrDateArray(int numDates, long *alibDates,
			DDate **drDates);
extern	int	DrlDrDateToAlibDate(
			DDate drDate, long *alibDate);

extern	int	DrlDDayCountFromNumDen(
			DDayCount *daycount, char num, int den);
extern	int	DrlDDayCountToNumDen(
			DDayCount daycount, char *num, int *den);



#if !defined(DRL_CLIB) 

extern	int	DrlDIntervalSet(
			int numPeriods, char periodType,
			DInterval *interval);
extern	int	DrlDDateAddDays(
			DDate startDate, int numDays, DDate *endDate);
extern	int	DrlDDateDaysDiff(
				DDate startDate, DDate endDate);
extern	int	DrlDDateFwdAny(
			DDate startDate, DInterval *interval,
			DDate *endDate);
extern	int	DrlDayCountFract(
			DDate date1, DDate date2,
			DDayCount method, double *fract);
extern	int	DrlDDateAdvanceYears(
			DDate baseDate, double years, DDate *newDate)	;
extern	int	DrlDIntervalToYears(
			DInterval *interval, double *years);
extern	int	DrlYearsToDInterval(
			double years, DInterval *interval);
extern	int	DrlFreqToDInterval(
			int freq, DInterval *interval);
extern	int	DrlDIntervalToFreq(
			DInterval *interval, int *freq);

#else 

/*
 * If ALIB, simply use the ALIB routines
 */
#define	DrlDIntervalSet	GtoMakeDateInterval
#define	DrlDDateAddDays(startDate, numDays, endDate)	\
		 *(endDate) = (startDate) + (numDays)
#define	DrlDDateDaysDiff(startDate, endDate)	\
		 ((int)((endDate) - (startDate)))
#define	DrlDDateFwdAny		GtoDtFwdAny
#define	DrlDayCountFract	GtoDayCountFraction
#define	DrlDDateAdvanceYears	GtoTDateAdvanceYears
#define	DrlDIntervalToYears	GtoDateIntervalToYears
#define	DrlYearsToDInterval	GtoYearsToDateInterval
#define	DrlFreqToDInterval	GtoFreq2TDateInterval
#define	DrlDIntervalToFreq	GtoDateIntervalToFreq

#endif





/*---------------------------------------------------------------
 * Structure DInterval
 */

extern	int	DrlDIntervalScan(char *s, DInterval *tl);
extern	char*	DrlDIntervalPrint(char *s, DInterval *tl);
extern	int	DrlCheckFreqValid(int freq);
extern	DInterval*	DrlDIntervalArrayAlloc(int size);
extern	int	DrlDIntervalArrayFree(DInterval* ptr);




#endif	/* _drltime_H */


