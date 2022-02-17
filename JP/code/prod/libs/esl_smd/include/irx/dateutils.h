/*
***************************************************************************
** HEADER FILE: dateutils.h
**
** Various date utilities. This is in contradiction to date.c which
** represents core date functions.
**
** This module understands day count conventions, date intervals etc.
**
** $Header$
***************************************************************************
*/

#ifndef IRX_DATEUTILS_H
#define IRX_DATEUTILS_H

#include "irxutils.h"   /* basic data types */

/**
***************************************************************************
** Computes the day count fraction between two dates using a day count
** convention. 
***************************************************************************
*/
int irxDayCountFraction
(IrxTDate         startDate,
 IrxTDate         endDate,
 IrxTDayCountConv dcc,
 double          *time);

/**
***************************************************************************
** Computes the number of days (numerator in day count fraction) between
** two dates using a day count convention.
***************************************************************************
*/
int irxDayCountDays
(IrxTDate         startDate,
 IrxTDate         endDate,
 IrxTDayCountConv dcc,
 long            *days);

/**
***************************************************************************
** Validates a date interval object.
***************************************************************************
*/
int irxDateIntervalValidate
(const IrxTDateInterval *interval);   /* (I) Date interval */

/**
***************************************************************************
** Makes a date interval from a string.
***************************************************************************
*/
IrxTDateInterval* irxDateIntervalMakeFromString (const char* str);

/**
***************************************************************************
** Convert a string to a date interval.
***************************************************************************
*/
int irxStringToDateInterval (const char* str, IrxTDateInterval* interval);

/**
***************************************************************************
** Converts a date interval to a string.
**
** The resulting string must be freed by the caller.
***************************************************************************
*/
char* irxDateIntervalToString (const IrxTDateInterval* ivl);

/**
***************************************************************************
** Adds a single date interval to a date.
***************************************************************************
*/
int irxDateAddInterval
(IrxTDate         startDate,      /* (I) Start date */
 IrxTDateInterval interval,       /* (I) Interval to add */
 IrxTDate        *endDate);       /* (O) Output date */

/*
***************************************************************************
** Adds a number of date intervals to a date.
**
** Note that adding a number of date intervals can be different from adding
** the date interval a number of times due to the effect of the end of the
** month. For example, add one month to 31 May and you get 30 Jun. If you
** then add another month to 30 Jun then depending on your end of month
** policy you may get 31 Jul or 30 Jul. However adding 2 months to 31 May
** you would always get 31 Jul.
***************************************************************************
*/
int irxDateAddMultiInterval
(IrxTDate         startDate,      /* (I) Start date */
 int              multi,          /* (I) Number of times to add interval */
 IrxTDateInterval interval,       /* (I) Interval to add */
 IrxTDate        *endDate);       /* (O) Output date */



/**
 * Counts number of date intervals in a range of dates. Note that fromDate
 * can be either later than or earlier than toDate.
 */
extern int irxCountDateIntervals(
    IrxTDate          fromDate,        /* (I) Date to count from */
    IrxTDate          toDate,          /* (I) Date to count to */
    IrxTDateInterval  interval,        /* (I) Interval to count */
    int              *numIntervals,    /* (O) Answer (Quotient) */
    int              *extraDays);      /* (O) Days left over(remainder) */

/**
 * Converts a date interval to a frequency basis.
 */
extern int  irxDateIntervalToBasis
(IrxTDateInterval interval,        /* (I) */
 double *freq);                    /* (O) # times per year */

/**
 * Converts a date interval to a frequency in char
 */
extern int  irxDateIntervalToFreq
(IrxTDateInterval interval,        /* (I) */
 char *freq);                    /* (O) 'A', 'S', 'Q', 'M', 'W', 'D', etc */


/**
 * Converts a date interval to a number of years. Note that if the
 * date interval is a month type (A,S,Q), the routine uses 30/360 to
 * compute the year fraction. If it is a day type (D,W), it
 * uses Act/365F.
 */
extern int  irxDateIntervalToYears
(IrxTDateInterval interval,        /* (I) */
 double *years);                   /* (O) # Years */


/**
 * Is the stub location at the front?
 */
IrxTBool irxStubLocationAtFront(IrxTStubLocation stubLocation);

/**
 * Is the stub location at the back?
 */
IrxTBool irxStubLocationAtBack(IrxTStubLocation stubLocation);

#endif

