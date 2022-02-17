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

#ifndef CX_DATEUTILS_H
#define CX_DATEUTILS_H

#include "cxutils.h"   /* basic data types */

#ifdef __cplusplus
extern "C"
{
#endif

/*f
***************************************************************************
** Computes the day count fraction between two dates using a day count
** convention. 
**
** If the day count convention is invalid, then behaviour is undefined.
***************************************************************************
*/
double CxDayCountFraction
(TDate startDate,
 TDate endDate,
 CxTDayCountConv dcc);

/*f
***************************************************************************
** Computes the number of days (numerator in day count fraction) between
** two dates using a day count convention.
**
** If the day count convention is invalid, then behaviour is undefined.
***************************************************************************
*/
long CxDayCountDays
(TDate startDate,
 TDate endDate,
 CxTDayCountConv dcc);

/*f
***************************************************************************
** Adds a single date interval to a date.
***************************************************************************
*/
int CxDateAddInterval
(TDate          startDate,      /* (I) Start date */
 TDateInterval *interval,       /* (I) Interval to add */
 TDate         *endDate);       /* (O) Output date */

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
int CxDateAddMultiInterval
(TDate          startDate,      /* (I) Start date */
 int            multi,          /* (I) Number of times to add interval */
 TDateInterval *interval,       /* (I) Interval to add */
 TDate         *endDate);       /* (O) Output date */

#ifdef __cplusplus
}
#endif

#endif

