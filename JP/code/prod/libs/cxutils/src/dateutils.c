/*
***************************************************************************
** SOURCE FILE: dateutils.c
**
** Various date utilities. This is in contradiction to date.c which
** represents core date functions.
**
** This module understands day count conventions, date intervals etc.
**
** $Header$
***************************************************************************
*/

#include "dateutils.h"

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "date.h"
#include "cxmacros.h"
#include "strutils.h"

#include <alib/ldate.h>

#define NaN sqrt(-1.0)

static double CxDayCountFractionActActISDA(TDate date1, TDate date2);

/*f
***************************************************************************
** Computes the day count fraction between two dates using a day count
** convention. 
**
** If the day count convention is invalid, then behaviour is undefined.
***************************************************************************
*/
double CxDayCountFraction
(TDate         date1,
 TDate         date2,
 CxTDayCountConv dcc)
{
    static char routine[] = "CxDayCountFraction";

    switch (dcc)
    {
    case CX_ACT_ACT:
        return CxDayCountFractionActActISDA(date1, date2);
    case CX_ACT_365F:
        return (date2-date1)/365.0;
    case CX_ACT_360:
        return (date2-date1)/360.0;
    case CX_B30_360:
    case CX_B30E_360:
    {
        long days = CxDayCountDays(date1,date2,dcc);
        return days/360.0;
    }
    }

    GtoErrMsg ("%s: Bad day count convention %ld\n", routine, (long)dcc);
    return NaN;
}

/*
***************************************************************************
** ACT/ACT day count fraction using the ISDA convention.
**
** Each day in a leap year counts as 1/366.
** Each day in a non-leap year counts as 1/365.
***************************************************************************
*/
static double CxDayCountFractionActActISDA
(TDate date1,
 TDate date2)
{
    double         result;
    double         sign;
    long           leap_days;
    long           non_leap_days;
    CxTYearMonthDay ymd1;
    CxTYearMonthDay ymd2;
    long           year;
    long           days;

    if (date1 == date2) return 0.0;

    if (date1 > date2)
    {
        TDate tmp = date1;
        date1 = date2;
        date2 = tmp;
        sign  = -1.0;
    }
    else
    {
        sign = 1.0;
    }

    if (CxDateToYMD (date1, &ymd1) != SUCCESS) return NaN;
    if (CxDateToYMD (date2, &ymd2) != SUCCESS) return NaN;

    leap_days     = 0;
    non_leap_days = 0;
    /* handle full years */
    for (year = ymd1.year+1; year < ymd2.year; ++year)
    {
        if (CX_IS_LEAP(year))
            leap_days += 366;
        else
            non_leap_days += 365;
    }

    /* handle partial years */
    if (ymd1.year == ymd2.year)
    {
        days = date2-date1;
        if (CX_IS_LEAP(ymd1.year))
            leap_days += date2-date1;
        else
            non_leap_days += date2-date1;
    }
    else
    {
        days = CxDateToEndOfYear(date1) - date1 + 1;
        if (CX_IS_LEAP(ymd1.year))
            leap_days += days;
        else
            non_leap_days += days;

        days = date2 - CxDateToStartOfYear(date2) + 1;
        if (CX_IS_LEAP(ymd2.year))
            leap_days += days;
        else
            non_leap_days += days;
    }

    result =  leap_days/366.0 + non_leap_days/365.0;
    result *= sign;

    return result;
}


/*f
***************************************************************************
** Computes the number of days (numerator in day count fraction) between
** two dates using a day count convention.
**
** If the day count convention is invalid, then behaviour is undefined.
***************************************************************************
*/
long CxDayCountDays
(TDate         date1,
 TDate         date2,
 CxTDayCountConv dcc)
{
    CxTYearMonthDay ymd1;
    CxTYearMonthDay ymd2;
    long           Y1,Y2,M1,M2,D1,D2;
    long           sign;
    long           result = FAILURE;

    if (date1 == date2) return 0;

    if (date1 > date2)
    {
        TDate tmp = date1;
        date1 = date2;
        date2 = tmp;
        sign = -1;
    }
    else
    {
        sign = 1;
    }

    switch (dcc)
    {
    case CX_ACT_ACT:
    case CX_ACT_365F:
    case CX_ACT_360:
        result = date2-date1;
        break;
    case CX_B30_360:
        if (CxDateToYMD(date1, &ymd1) != SUCCESS) return LONG_MAX;
        if (CxDateToYMD(date2, &ymd2) != SUCCESS) return LONG_MAX;

        Y1 = ymd1.year;
        M1 = ymd1.month;
        D1 = ymd1.day;

        Y2 = ymd2.year;
        M2 = ymd2.month;
        D2 = ymd2.day;

        if (D1 == 31) D1 = 30;
        if (D2 == 31 && D1 == 30) D2 = 30;

        result = (Y2-Y1)*360 + (M2-M1)*30 + (D2-D1);
        break;

    case CX_B30E_360:
        if (CxDateToYMD(date1, &ymd1) != SUCCESS) return LONG_MAX;
        if (CxDateToYMD(date2, &ymd2) != SUCCESS) return LONG_MAX;

        Y1 = ymd1.year;
        M1 = ymd1.month;
        D1 = ymd1.day;

        Y2 = ymd2.year;
        M2 = ymd2.month;
        D2 = ymd2.day;

        if (D1 == 31) D1 = 30;
        if (D2 == 31) D2 = 30;

        result = (Y2-Y1)*360 + (M2-M1)*30 + (D2-D1);
        break;
    }

    result = result * sign;
    return result;
}

/*f
***************************************************************************
** Adds a single date interval to a date.
***************************************************************************
*/
int CxDateAddInterval
(TDate          startDate,      /* (I) Start date */
 TDateInterval *interval,       /* (I) Interval to add */
 TDate         *endDate)        /* (O) Output date */
{
    static char    routine[]="CxDateAddInterval";
    int            status = FAILURE; /* Until proven successful */

    REQUIRE (interval != NULL);
    REQUIRE (startDate > 0);
    REQUIRE (endDate != NULL);

    if (GtoDtFwdAny (startDate, interval, endDate) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;
    
 done:
    if (status != SUCCESS)
        GtoErrMsgFailure(routine);
    
    return(status);
}



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
 TDate         *endDate)        /* (O) Output date */
{
    static char routine[] = "CxDateAddMultiInterval";
    int         status    = FAILURE;

    TDateInterval multiInterval;

    REQUIRE (interval != NULL);
    REQUIRE (endDate != NULL);

    SET_TDATE_INTERVAL (multiInterval, 
                        interval->prd * multi, 
                        interval->prd_typ);
    if (CxDateAddInterval (startDate, &multiInterval, endDate) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:
    
    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}
                                              

