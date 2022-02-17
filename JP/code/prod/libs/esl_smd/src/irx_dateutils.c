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

#include "irx/dateutils.h"

#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "irx/convert.h"
#include "irx/macros.h"
#include "irx/strutils.h"

static int irxDayCountFractionActAct(IrxTDate date1, IrxTDate date2, double *time);


/*f
***************************************************************************
** Computes the day count fraction between two dates using a day count
** convention. 
**
** If the day count convention is invalid, then behaviour is undefined.
***************************************************************************
*/
int irxDayCountFraction
(IrxTDate          date1,
 IrxTDate          date2,
 IrxTDayCountConv  dcc,
 double           *time)
{
    static char routine[] = "irxDayCountFraction";

    switch (dcc)
    {
    case IRX_ACT_ACT:
        return irxDayCountFractionActAct(date1, date2, time);
    case IRX_ACT_365F:
        *time = (date2-date1)/365.0;
        return SUCCESS;
    case IRX_ACT_360:
        *time = (date2-date1)/360.0;
        return SUCCESS;
    case IRX_B30_360:
    case IRX_B30E_360:
    {
        long days;
        if (irxDayCountDays(date1,date2,dcc,&days) != SUCCESS)
            return irxErrorFailure(routine);
        *time = (double)days/360.0;
        return SUCCESS;
    }
    }

    irxError ("%s: Bad day count convention %ld\n", routine, (long)dcc);
    return FAILURE;
}

/*
***************************************************************************
** ACT/ACT day count fraction using swap conventions.
**
** Each day in a leap year counts as 1/366.
** Each day in a non-leap year counts as 1/365.
***************************************************************************
*/
static int irxDayCountFractionActAct
(IrxTDate date1,
 IrxTDate date2,
 double  *time)
{
    double         result;
    double         sign;
    long           leap_days;
    long           non_leap_days;
    IrxTYearMonthDay ymd1;
    IrxTYearMonthDay ymd2;
    long           year;
    long           days;

    if (date1 == date2) 
    {
        *time = 0.0;
        return SUCCESS;
    }

    if (date1 > date2)
    {
        IrxTDate tmp = date1;
        date1 = date2;
        date2 = tmp;
        sign  = -1.0;
    }
    else
    {
        sign = 1.0;
    }

    if (irxDateToYMD (date1, &ymd1) != SUCCESS) return FAILURE;
    if (irxDateToYMD (date2, &ymd2) != SUCCESS) return FAILURE;

    leap_days     = 0;
    non_leap_days = 0;
    /* handle full years */
    for (year = ymd1.year+1; year < ymd2.year; ++year)
    {
        if (IRX_IS_LEAP(year))
            leap_days += 366;
        else
            non_leap_days += 365;
    }

    /* handle partial years */
    if (ymd1.year == ymd2.year)
    {
        days = date2-date1;
        if (IRX_IS_LEAP(ymd1.year))
            leap_days += date2-date1;
        else
            non_leap_days += date2-date1;
    }
    else
    {
        days = irxDateToEndOfYear(date1) - date1 + 1;
        if (IRX_IS_LEAP(ymd1.year))
            leap_days += days;
        else
            non_leap_days += days;

        days = date2 - irxDateToStartOfYear(date2) + 1;
        if (IRX_IS_LEAP(ymd2.year))
            leap_days += days;
        else
            non_leap_days += days;
    }

    result =  leap_days/366.0 + non_leap_days/365.0;
    result *= sign;

    *time = result;
    return SUCCESS;
}


/*f
***************************************************************************
** Computes the number of days (numerator in day count fraction) between
** two dates using a day count convention.
**
** If the day count convention is invalid, then behaviour is undefined.
***************************************************************************
*/
int irxDayCountDays
(IrxTDate         date1,
 IrxTDate         date2,
 IrxTDayCountConv dcc,
 long            *days)
{
    static char routine[] = "irxDayCountDays";

    IrxTYearMonthDay ymd1;
    IrxTYearMonthDay ymd2;
    long           Y1,Y2,M1,M2,D1,D2;
    long           sign;
    long           result;

    if (date1 == date2) 
    {
        *days = 0;
        return SUCCESS;
    }

    if (date1 > date2)
    {
        IrxTDate tmp = date1;
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
    case IRX_ACT_ACT:
    case IRX_ACT_365F:
    case IRX_ACT_360:
        result = date2-date1;
        break;
    case IRX_B30_360:
        if (irxDateToYMD(date1, &ymd1) != SUCCESS) return LONG_MAX;
        if (irxDateToYMD(date2, &ymd2) != SUCCESS) return LONG_MAX;

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

    case IRX_B30E_360:
        if (irxDateToYMD(date1, &ymd1) != SUCCESS) return LONG_MAX;
        if (irxDateToYMD(date2, &ymd2) != SUCCESS) return LONG_MAX;

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
    default:
        irxError ("%s: Bad day count convention %ld\n", routine, (long)dcc);
        return FAILURE;
    }

    result = result * sign;
    *days  = result;
    return SUCCESS;
}


/*f
***************************************************************************
** Validates a date interval object.
**
** Returns SUCCESS or FAILURE.
***************************************************************************
*/
int irxDateIntervalValidate
(const IrxTDateInterval *interval)   /* (I) Date interval */
{
    static char routine[] = "irxDateIntervalValidate";
    int         status    = FAILURE;

    REQUIRE(interval != NULL);
    
    switch(interval->prd_typ)
    {
    case IRX_PRD_TYPE_ANNUAL:
    case IRX_PRD_TYPE_SEMI_ANNUAL:
    case IRX_PRD_TYPE_YEAR:
    case IRX_PRD_TYPE_MONTH:
    case IRX_PRD_TYPE_QUARTER:
        break; /* OK */
    case IRX_PRD_TYPE_WEEK:
    case IRX_PRD_TYPE_DAY:
    case IRX_PRD_TYPE_LUNAR_MONTH:
        if (interval->eom)
        {
            irxError ("%s: Should not specify end of month adjustment for "
                     "non-monthly adjustment %c\n", 
                     routine, interval->prd_typ);
            goto RETURN; /* failure */
        }
        break; /* OK */
    default:
        irxError ("%s: Invalid period type %c\n", routine, interval->prd_typ);
        goto RETURN; /* failure */
    }

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS) irxErrorFailure (routine);
    return status;
}



/*
***************************************************************************
** Makes a date interval from a string.
***************************************************************************
*/
int irxStringToDateInterval (const char* str, IrxTDateInterval* interval)
{
    static char routine[] = "irxStringToDateInterval";
    int         status    = FAILURE;

    int            numPeriods;
    char           *ep, ec;
    IrxTBool       eom;
    char           *uc = NULL;
    char           *eomStr;
    char           *p;

    REQUIRE (str      != NULL);
    REQUIRE (interval != NULL);

    uc = irxStringToUpper (str);
    if (irxStringParser (uc, ",", &eomStr) != SUCCESS)
        goto RETURN;
    
    if (eomStr == NULL)
    {
        eom = FALSE;
    }
    else if (strcmp (eomStr, "EOM") == 0)
    {
        eom = TRUE;
    }
    else
    {
        irxError ("%s: Bad format string %s - should be <number>"
                 "<interval_character>{,EOM} where the optional ,EOM "
                 "indicates end of month adjustment\n", routine, str);
        goto RETURN; /* failure */
    }
    

    /* Special case */
    if ((!strcmp(uc, "O/N")) ||
        (!strcmp(uc, "ON"))  ||
        (!strcmp(uc, "SN"))  ||
        (!strcmp(uc, "S/N")))
    {
        numPeriods    = 1;
        ec            = 'D';
    }

    /* Scan for IMMn intervals */
    else if ((p = strstr(uc, "IMM")) != NULL)
    {
        p +=3;
        if (sscanf(p, "%d", &numPeriods) != 1)
            goto RETURN;

        numPeriods    = 1;
        ec            = 'I';
    }
    else
    {
        /* format of the string is <number><interval_character> */
        /* hence strtol is ideal for getting the number and interval */

        numPeriods = (int) strtol(uc, &ep, 10);
        if (ep == uc)
        { 
            /* means that the user went straight to the interval 
             * type without number 
             */
            numPeriods = 1; /* first character is the interval */
        }

        ec = *ep;
    }

    /* replace '3y' by '3A' */
    if (ec == 'Y')
        ec = 'A';

    interval->prd             = numPeriods;
    interval->prd_typ         = ec;
    interval->eom             = eom;

    if (irxDateIntervalValidate(interval) != SUCCESS) goto RETURN; /* failure */

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS) irxErrorFailure (routine);
    FREE (uc);
    return status;
}



/*
***************************************************************************
** Makes a date interval from a string.
***************************************************************************
*/
IrxTDateInterval* irxDateIntervalMakeFromString (const char* str)
{
    static char routine[] = "irxDateIntervalMakeFromString";

    IrxTDateInterval *interval = NULL;

    interval = NEW(IrxTDateInterval);
    if (interval==NULL) goto RETURN;

    if (irxStringToDateInterval(str, interval) != SUCCESS) goto RETURN;

 RETURN:

    if (interval == NULL) irxErrorFailure(routine);
    return interval;
}



/*
***************************************************************************
** Converts a date interval to a string.
**
** The resulting string must be freed by the caller.
***************************************************************************
*/
char* irxDateIntervalToString (const IrxTDateInterval* ivl)
{
    static char routine[] = "irxDateIntervalToString";

    char *str = NULL;
    char  buf[128];

    REQUIRE (ivl != NULL);

    if (ivl->eom)
    {
        sprintf(buf, "%d%c,EOM", ivl->prd, ivl->prd_typ);
    }
    else
    {
        sprintf(buf, "%d%c", ivl->prd, ivl->prd_typ);
    }

    str = STRDUP(buf);

 RETURN:

    if (str == NULL) 
        irxErrorFailure(routine);

    return str;
}


/*f
***************************************************************************
** Adds a single date interval to a date.
***************************************************************************
*/
int irxDateAddInterval
(IrxTDate         startDate,      /* (I) Start date */
 IrxTDateInterval interval,       /* (I) Interval to add */
 IrxTDate        *endDate)        /* (O) Output date */
{
    static char    routine[]="irxDateAddInterval";
    int            status = FAILURE; /* Until proven successful */

    REQUIRE (irxDateIntervalValidate (&interval) == SUCCESS);
    REQUIRE (startDate > 0);
    REQUIRE (endDate != NULL);

    switch (interval.prd_typ)
    {
    case IRX_PRD_TYPE_MONTH:
    case IRX_PRD_TYPE_ANNUAL:
    case IRX_PRD_TYPE_YEAR:
    case IRX_PRD_TYPE_SEMI_ANNUAL:
    case IRX_PRD_TYPE_QUARTER:
    {
        long months;
        switch (interval.prd_typ)
        {
        case IRX_PRD_TYPE_MONTH:
            months = interval.prd;
            break;
        case IRX_PRD_TYPE_ANNUAL:
        case IRX_PRD_TYPE_YEAR:
            months = interval.prd * 12;
            break;
        case IRX_PRD_TYPE_SEMI_ANNUAL:
            months = interval.prd * 6;
            break;
        case IRX_PRD_TYPE_QUARTER:
            months = interval.prd * 3;
            break;
        default: 
            PROGRAM_BUG();
            goto RETURN; 
        }
            
        *endDate = irxDateAddMonths (startDate, months, interval.eom);
        break;
    }
    case IRX_PRD_TYPE_DAY:
        *endDate = startDate + interval.prd;
        break;

    case IRX_PRD_TYPE_WEEK:
        *endDate = startDate + 7 * interval.prd;
        break;                

    case IRX_PRD_TYPE_LUNAR_MONTH:
        *endDate = startDate + 28 * interval.prd;
        break;

    default:
        PROGRAM_BUG(); /* due to validation call earlier */
        goto RETURN;
    }
    
    status = SUCCESS;
    
 RETURN:
    if (status != SUCCESS)
        irxErrorFailure(routine);
    
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
int irxDateAddMultiInterval
(IrxTDate         startDate,      /* (I) Start date */
 int              multi,          /* (I) Number of times to add interval */
 IrxTDateInterval interval,       /* (I) Interval to add */
 IrxTDate        *endDate)        /* (O) Output date */
{
    static char routine[] = "irxDateAddMultiInterval";
    int         status    = FAILURE;

    IrxTDateInterval multiInterval;

    REQUIRE (endDate != NULL);

    multiInterval.prd     = interval.prd * multi;
    multiInterval.prd_typ = interval.prd_typ;
    multiInterval.eom     = interval.eom;
    if (irxDateAddInterval (startDate, multiInterval, endDate) != SUCCESS)
        goto RETURN; /* failure */

    status = SUCCESS;

 RETURN:
    
    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}


/**
 * Counts number of date intervals in a range of dates. Note that fromDate
 * can be either later than or earlier than toDate.
 *
 * Promises the following:
 *
 *      fromDate + numIntervals*interval
 *
 * will be in the range [fromDate,toDate] (remembering that toDate can be
 * less than fromDate).
 *
 *      extraDays = toDate - (fromDate - numIntervals*interval)
 *
 * This can be negative - but only when toDate < fromDate.
 */
int irxCountDateIntervals(
    IrxTDate          fromDate,        /* (I) Date to count from */
    IrxTDate          toDate,          /* (I) Date to count to */
    IrxTDateInterval  interval,        /* (I) Interval to count */
    int              *numIntervals,    /* (O) Answer (Quotient) */
    int              *extraDays)       /* (O) Days left over(remainder) */
{
    static char routine[]="irxCountDates";
    int status = FAILURE;     /* Until proven successful */
    double intervalYears;     /* TDateInterval expressed in years */
    int lowNumIntervals;      /* low estimation of *numIntervals */
    int index;                /* running approximation of *numIntervals */
    IrxTDate date;            /* date index times interval away from fromDate*/
    double approxIntervals;
    int increment;

    REQUIRE(interval.prd > 0);

    if (fromDate == toDate)
    {
        *numIntervals = 0;
        *extraDays    = 0;
        status = SUCCESS;
        goto RETURN; /* success */
    }

    if (irxDateIntervalToYears(interval, &intervalYears) != SUCCESS)
        goto RETURN; /* failure */
    
    approxIntervals = (double)(toDate-fromDate)/365.0 / intervalYears;

    if (approxIntervals > 0)
    {
        lowNumIntervals = (int)approxIntervals-2;
        if (lowNumIntervals < 0)
            lowNumIntervals = 0;
        increment = 1;
    }
    else
    {
        /* note that (int) always rounds towards 0 */
        lowNumIntervals = (int)approxIntervals+2;
        if (lowNumIntervals > 0)
            lowNumIntervals = 0;
        increment = -1;
    }

    /* Initialize the LOOP by index of intervals
     */
    index = lowNumIntervals;
    if (irxDateAddMultiInterval(fromDate,index,interval,&date) != SUCCESS)
        goto RETURN;                       /* Failed */

   /* Keep advancing date while date is between fromDate and toDate.
    */
   while ( IS_BETWEEN (date,fromDate,toDate) )
   {   
       index   += increment;
       if (irxDateAddMultiInterval(fromDate,index,interval,&date) != SUCCESS)
           goto RETURN; /* failure */
   }
   
   index -= increment;

   /* do a final check that all is good */
   if (irxDateAddMultiInterval (fromDate, index, interval, &date) != SUCCESS)
       goto RETURN; /* failure */

   if (IS_BETWEEN(date,fromDate,toDate))
   {
       *numIntervals = index;
       *extraDays    = (toDate - date);
   }
   else
   {
       PROGRAM_BUG();
       goto RETURN; /* failure */
   }
   
   status = SUCCESS;
   
 RETURN:

   if (status != SUCCESS)
       irxErrorFailure (routine);
   return status;
}
    

/*
** Converts a date interval to a basis in double
*/
int irxDateIntervalToBasis(IrxTDateInterval interval, /* (I) */
                           double *basis)             /* (O) # times per year */
{
    static char routine[] = "irxDateIntervalToBasis";
    int         status    = FAILURE;
    
    REQUIRE(interval.prd > 0);

    switch(toupper(interval.prd_typ))
    {
    case IRX_PRD_TYPE_ANNUAL:
    case IRX_PRD_TYPE_YEAR:
        *basis = 1.0/(double)interval.prd;
        break;
    case IRX_PRD_TYPE_SEMI_ANNUAL:
        *basis = 2.0/(double)(interval.prd);
        break;
    case IRX_PRD_TYPE_QUARTER:
        *basis = 4.0/(double)(interval.prd);
        break;
    case IRX_PRD_TYPE_MONTH:
        *basis = 12.0/(double)(interval.prd);
        break;
    case IRX_PRD_TYPE_LUNAR_MONTH:
        *basis = 365.0/(28.0*(double)interval.prd);
        break;
    case IRX_PRD_TYPE_WEEK:
        *basis = 365.0/(7.0*(double)interval.prd);
        break;
    case IRX_PRD_TYPE_DAY:
        *basis = 365.0/(double)interval.prd;
        break;
    default:
        irxError("%s: unknown interval type %c.\n",
                 routine, interval.prd_typ);
        goto RETURN; /* failure */
    }

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);
    return status;
}
    


/*
** Converts a date interval to a frequency in char
*/
int irxDateIntervalToFreq(IrxTDateInterval interval, /* (I) */
                          char *freq)                /* (O) # times per year */
{
    static char routine[] = "irxDateIntervalToFreq";
    int         status    = FAILURE;

    double basis;
    char   c;
    
    REQUIRE(interval.prd > 0);

    if (irxDateIntervalToBasis(interval, &basis) != SUCCESS) goto RETURN;

    switch((int)(basis)){
    case 1:
        c = 'A';
        break;
    case 2:
        c = 'S';
        break;
    case 4:
        c = 'Q';
        break;
    case 12:
        c = 'M';
        break;
    case 52:
        c = 'W';
        break;
    case 365:
        c = 'D';
        break;
    default:
        irxError("%s: invalid interval %s.\n",
                 routine, irxDateIntervalToString(&interval));
        goto RETURN; /* failure */
    }

    freq[0] = c;

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);
    return status;
}
    

/**
 * Converts a date interval to a number of years. Note that if the
 * date interval is a month type (A,S,Q), the routine uses 30/360 to
 * compute the year fraction. If it is a day type (D,W), it
 * uses Act/365F.
 */
extern int irxDateIntervalToYears
(IrxTDateInterval interval,        /* (I) */
 double          *years)           /* (O) # Years */
{
    static char routine[] = "irxDateIntervalToYears";
    int         status    = FAILURE;
    switch(toupper(interval.prd_typ))
    {
    case IRX_PRD_TYPE_ANNUAL:
    case IRX_PRD_TYPE_YEAR:
        *years = (double)interval.prd;
        break;
    case IRX_PRD_TYPE_SEMI_ANNUAL:
        *years = (double)(interval.prd)/2.0;
        break;
    case IRX_PRD_TYPE_QUARTER:
        *years = (double)(interval.prd)/4.0;
        break;
    case IRX_PRD_TYPE_MONTH:
        *years = (double)(interval.prd)/12.0;
        break;
    case IRX_PRD_TYPE_LUNAR_MONTH:
        *years = (double)interval.prd * 28.0 / 365.0;
        break;
    case IRX_PRD_TYPE_WEEK:
        *years = (double)interval.prd * 7.0 / 365.0;
        break;
    case IRX_PRD_TYPE_DAY:
        *years = (double)interval.prd/365.0;
        break;
    default:
        irxError("%s: unknown interval type %c.\n", routine, interval.prd_typ);
        goto RETURN; /* failure */
    }

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}





/**
 * Is the stub location at the front?
 */
IrxTBool irxStubLocationAtFront(IrxTStubLocation stubLocation)
{
    switch (stubLocation)
    {
    case IRX_SHORT_FRONT_STUB:
    case IRX_LONG_FRONT_STUB:
        return TRUE;
    case IRX_SHORT_BACK_STUB:
    case IRX_LONG_BACK_STUB:
        return FALSE;
    default:
        PROGRAM_BUG();
        return FALSE;
    }
}


/**
 * Is the stub location at the back?
 */
IrxTBool irxStubLocationAtBack(IrxTStubLocation stubLocation)
{
    switch (stubLocation)
    {
    case IRX_SHORT_FRONT_STUB:
    case IRX_SHORT_BACK_STUB:
        return TRUE;
    case IRX_LONG_FRONT_STUB:
    case IRX_LONG_BACK_STUB:
        return FALSE;
    default:
        PROGRAM_BUG();
        return FALSE;
    }
}
