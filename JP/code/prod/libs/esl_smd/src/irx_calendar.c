/*
***************************************************************************
** SOURCE FILE: calendar.c
**
** Implements the calendar object, the calendar cache (so we can access
** calendar objects by name) and basic business day functions.
**
** The authors of this module include Doug Gallager, Simon Meldrum, Peter
** Taylor, Jim Mitsiopoulos.
**
** $Header$
***************************************************************************
*/

#include "irx/calendar.h"

#include <assert.h>
#include <ctype.h>
#include <stdio.h>

#include "irx/bsearch.h"             /* irxBinarySearchLong */
#include "irx/datelist.h"            /* irxDateList... */
#include "irx/error.h"               /* irxError */
#include "irx/hash.h"                /* IRX_HASH_... */
#include "irx/macros.h" 
#include "irx/strutils.h"

/*
***************************************************************************
** Private function declarations
***************************************************************************
*/
static int getNextBusDateMulti
(IrxTDate         startDate,     /* (I) starting date. */
 long           direction,     /* (I) +1=forwards, -1=backwards */
 int            numCalendars,  /* (I) Size of calendars array. */
 IrxTCalendar const  **calendars,      /* (I) calendars */
 IrxTDate         *nextDate);    /* (O) next business day */

static int   getNextBusDate 
(IrxTDate      startDate,    /* (I) starting date. */
 long        direction,    /* (I) +1=forwards, -1=backwards */
 const IrxTCalendar *cal,          /* (I) calendar */
 IrxTDate     *nextDate);    /* (O) next business day */

static int   forwardNonStandardWeekends
(IrxTDate         fromDate,       /* (I) start date */
 long           numBusDaysLeft, /* (I) abs. num. bus. days */
 long           intervalSign,   /* (I) +1=forwards, -1=backwards */
 long           weekends,       /* (I) Weekends flags */
 long           busDaysPerWeek, /* (I) num. bus. days per week */
 const IrxTDateList    *hlDatelist,     /* (I) holiday dateList */
 IrxTDate        *resultDate);    /* (O) forwarded date. */

static int   holAdjustStandardWeekends 
(IrxTDate          fromDate,       /* (I) start date */
 const IrxTDateList     *hlDatelist,     /* (I) holiday dateList */
 IrxTDate         *resultDate);    /* (I/O) forwarded date. */

static long  calcNumWeekdayHolidays 
(IrxTDate             toDate,     /* (I) End date.   */
 long               startHolIdx,/* (I) Starting holiday index. */
 long               direction,  /* (I) 1=forwards, 2=backwards*/
 long               weekends,   /* (I) Weekends flags */
 const IrxTDateList        *hlDatelist, /* (I) Holiday datelist. */
 long              *endHolIdx); /* (O) Idx where hol=toDate, +1 */

static int   findFirstHolidayIdx 
(IrxTDate         date,           /* (I) Date to search for. */
 const IrxTDateList    *hlDatelist,     /* (I) Holiday date list. */
 long           direction,      /* (I) +1=forwards, -1=backwards. */
 long          *holIdx,         /* (O) Index of first holiday. */
 IrxTBool        *doneSearching); /* (O) TRUE=date is beyond list bounds. */

static IrxTHashTable* calendarCacheGet(void);
static char* cacheNameGet(const char *name);

/*---------------------------------------------------------------------------
 *                              Static Data.
 *---------------------------------------------------------------------------
 */

static IrxTHashTable* g_calendarCache = NULL;


/*
 * Table containing the offset to add to a date to skip a given number of
 * business days.  Only valid for less than 5 business days, and when saturday
 * and sunday are both holidays.
 *
 * So:
 * g_offsetTable[dayOfWeek(date)][nDays] is the number of days to add to date
 * to get a new business day.
 *
 * Note: The indexes are set up as IrxTDate mod 7, which is slightly different
 * than GtoDayOfWeek(), which adds one to IrxTDate.
 */
static long  g_offsetTable[7][5] =
{/*n:  0  1  2  3  4 */
    
    {  0, 1, 2, 3, 4 },  /* Monday */
    {  0, 1, 2, 3, 6 },  /* Tuesday */
    {  0, 1, 2, 5, 6 },  /* Wednesday */
    {  0, 1, 4, 5, 6 },  /* Thursday */
    {  0, 3, 4, 5, 6 },  /* Friday */
    { -1, 2, 3, 4, 5 },  /* Saturday */
    { -2, 1, 2, 3, 4 }   /* Sunday */
};
/*---------------------------------------------------------------------------
 *                              End Static Data.
 *---------------------------------------------------------------------------
 */

/*f
***************************************************************************
** Converts weekends to string format.
**
** Returns a static string which should not be freed - as a result this
** routine is not thread safe since the string is different for different
** values of weekends.
**
** The string returned consists of all the weekends seperated by commas.
** The names are the first three letters of the weekend name.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxWeekendsToString (long weekends)
{
    static char str[32];

    strcpy(str, "");

    if (weekends < IRX_WEEKEND_MIN_VALUE || weekends > IRX_WEEKEND_MAX_VALUE)
    {
        strcat (str, "***ERROR***");
    }
    else
    {
        if (weekends & IRX_WEEKEND_MONDAY)    strcat(str, "MON,");
        if (weekends & IRX_WEEKEND_TUESDAY)   strcat(str, "TUE,");
        if (weekends & IRX_WEEKEND_WEDNESDAY) strcat(str, "WED,");
        if (weekends & IRX_WEEKEND_THURSDAY)  strcat(str, "THU,");
        if (weekends & IRX_WEEKEND_FRIDAY)    strcat(str, "FRI,");
        if (weekends & IRX_WEEKEND_SATURDAY)  strcat(str, "SAT,");
        if (weekends & IRX_WEEKEND_SUNDAY)    strcat(str, "SUN,");
        if (strlen(str) > 0) str[strlen(str)-1] = '\0';
    }
    return str;
}

/*f
***************************************************************************
** Converts from string to weekends.
**
** The format is that the weekends are separated by commas. The names are
** the first three letters of the weekend name.
**
** Returns SUCCESS or FAILURE.
***************************************************************************
*/
int irxWeekendsFromString (char *str, long *weekends)
{
    static char routine[] = "irxWeekendsFromString";
    int         status    = FAILURE;

    size_t  i;
    size_t  numItems;
    char  **split = NULL;
    long    myWeekends = 0;

    REQUIRE(str != NULL);
    REQUIRE(weekends != NULL);

    if (strlen(str) == 0)
        goto success;

    if (irxStringSplit (str, ',', &numItems, &split) != SUCCESS)
        goto RETURN; /* failure */

    for (i = 0; i < numItems; ++i)
    {
        if (irxStrcmpi (split[i], "MON") == 0)
            myWeekends |= IRX_WEEKEND_MONDAY;
        else if (irxStrcmpi (split[i], "TUE") == 0)
            myWeekends |= IRX_WEEKEND_TUESDAY;
        else if (irxStrcmpi (split[i], "WED") == 0)
            myWeekends |= IRX_WEEKEND_WEDNESDAY;
        else if (irxStrcmpi (split[i], "THU") == 0)
            myWeekends |= IRX_WEEKEND_THURSDAY;
        else if (irxStrcmpi (split[i], "FRI") == 0)
            myWeekends |= IRX_WEEKEND_FRIDAY;
        else if (irxStrcmpi (split[i], "SAT") == 0)
            myWeekends |= IRX_WEEKEND_SATURDAY;
        else if (irxStrcmpi (split[i], "SUN") == 0)
            myWeekends |= IRX_WEEKEND_SUNDAY;
        else
        {
            irxError ("%s: Invalid date string %s\n", routine, split[i]);
            goto RETURN; /* failure */
        }
    }

 success:
    *weekends = myWeekends;
    status = SUCCESS;

 RETURN:

    FREE(split);
    if (status != SUCCESS)
        irxErrorFailure(routine);

    return status;
}


/*f
***************************************************************************
** Creates a new calendar structure.
**
** The date list can be NULL, in which case the resulting date list in the
** holiday structure will be a date list with no dates, e.g.
**      cal->dateList->fNumItems = 0;
***************************************************************************
*/
IrxTCalendar* irxCalendarNew
(const IrxTDateList *dateList,              /* (I) Date list to use */
 IrxTBool      satIsAlwaysHoliday,    /* (I) TRUE if saturday always holiday*/
 IrxTBool      sunIsAlwaysHoliday     /* (I) TRUE if sunday always holiday */
)
{
    static char   routine[] = "irxCalendarNew";

    IrxTCalendar *cal = NULL;
    long        weekends;

    weekends = IRX_WEEKEND_NO_WEEKENDS;
    if (satIsAlwaysHoliday)
        weekends |= IRX_WEEKEND_SATURDAY;

    if (sunIsAlwaysHoliday)
        weekends |= IRX_WEEKEND_SUNDAY;

    cal = irxCalendarMake (dateList, weekends);
    if (cal == NULL) irxErrorFailure(routine);

    return cal;
}

/*f
***************************************************************************
** Validate a newly created calendar structure.
**
** This will do the following:
** 1. Sort the calendar date list and remove weekends.
** 2. Ensure that we don't have all days as weekends.
***************************************************************************
*/
int irxCalendarValidate
(IrxTCalendar *cal)     /* (I/O) Calendar to be adjusted */
{
    static char   routine[] = "irxCalendarValidate";
    int           status = FAILURE;

    long        i1;
    long        i2;

    REQUIRE (cal != NULL);
    REQUIRE (cal->weekends >= IRX_WEEKEND_MIN_VALUE);
    REQUIRE (cal->weekends <= IRX_WEEKEND_MAX_VALUE);

    if (cal->weekends == IRX_WEEKEND_MAX_VALUE)
    {
        irxError ("%s: All weekdays are defined as weekends\n", routine);
        goto RETURN; /* failure */
    }

    if (cal->dateList == NULL) goto success;

    /*
    ** Remove weekends from the date list
    */
    i1 = 0;
    i2 = 0;
    while (i1 < cal->dateList->fNumItems)
    {
        IrxTDate thisDate = cal->dateList->fArray[i1];
        if (IRX_IS_WEEKEND(thisDate, cal->weekends))
        {
            i1++;
        }
        else
        {
            cal->dateList->fArray[i2] = thisDate;
            i1++;
            i2++;
        }
    }
    cal->dateList->fNumItems = i2;

    /* irxDateListSort also removes duplicates */
    irxDateListSort(cal->dateList);

 success:
    status = SUCCESS;

 RETURN:
    
    if (status != SUCCESS) irxErrorFailure(routine);
    return status;
}


/*
***************************************************************************
** Indicates whether a date is a business day according to a calendar.
***************************************************************************
*/
int irxIsBusinessDay
(IrxTDate      date,               /* (I) Input Date */
 const IrxTCalendar *cal,                /* (I) Holiday list structure */
 IrxTBool     *isBusinessDay)      /* (O) TRUE if a business day */
{
    static char routine[] = "irxIsBusinessDay";
    int         status = FAILURE;

    IrxTBool      aHoliday;

    REQUIRE (cal != NULL);
    REQUIRE (isBusinessDay != NULL);

    if (IRX_IS_WEEKEND (date, cal->weekends))
    {
        *isBusinessDay = FALSE;
        return SUCCESS;
    }

    if (irxIsHoliday (date, cal, &aHoliday) != SUCCESS)
        goto RETURN;

    *isBusinessDay = !aHoliday;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure(routine);
        
    return status;
}

/*
***************************************************************************
** Checks to see if a date is a holiday. Does not take week-ends into
** account.
***************************************************************************
*/
int irxIsHoliday
(IrxTDate      date,      /* (I) Arbitrary date      */
 const IrxTCalendar *cal,       /* (I) Calendar structure  */
 IrxTBool     *isHoliday) /* (O) 0 = not a holiday, 1 = is a holiday */
{
    static char routine[] = "irxIsHoliday";
    int         status = FAILURE;

    long        exact;

    REQUIRE (cal != NULL);
    REQUIRE (cal->dateList != NULL);

    if (cal->dateList->fNumItems == 0)
    {
        *isHoliday = FALSE;
        return SUCCESS;
    }

    /*
    ** Do a binary search of the date list to see if date is in the list.
    */
    if (irxBinarySearchLong (date,
                             cal->dateList->fArray,
                             sizeof(IrxTDate),
                             cal->dateList->fNumItems,
                             &exact,
                             NULL,
                             NULL) != SUCCESS)
        goto RETURN; /* failure */

    *isHoliday = (exact >= 0);

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS) irxErrorFailure(routine);
    return status;
}

/*
***************************************************************************
** Converts date to a valid business date using the calendar and bad day
** convention.
**
** Returns SUCCESS if convention OK, FAILURE otherwise.
***************************************************************************
*/
int irxBusinessDay
(IrxTDate        date,       /* (I) Arbitrary date. */
 long          badDayConv, /* (I) Bad day convention for adjusting
                              non-business days. Use one of the following:
                              IRX_BAD_DAY_FOLLOW
                              IRX_BAD_DAY_PREVIOUS
                              IRX_BAD_DAY_MODIFIED
                              IRX_BAD_DAY_NONE */
 const IrxTCalendar  *cal,         /* (I) Calendar. */
 IrxTDate      *outDate)     /* (O) Valid business day. */

{
    static char routine[] = "irxBusinessDay";
    int         status;

    status = irxBusinessDayMultiCalendar(date, badDayConv, 1, &cal, outDate);

    if (status != SUCCESS) irxErrorFailure(routine);
    return status;
}

/*
***************************************************************************
** Converts date to a valid business date using all the given calendars and
** bad day convention.
**
** Returns SUCCESS if convention OK, FAILURE otherwise.
***************************************************************************
*/
int irxBusinessDayMultiCalendar
(IrxTDate        date,          /* (I) Arbitrary date. */
 long          badDayConv,    /* (I) Bad day convention for adjusting
                                 non-business days. Use one of the following:
                                 IRX_BAD_DAY_FOLLOW
                                 IRX_BAD_DAY_PREVIOUS
                                 IRX_BAD_DAY_MODIFIED
                                 IRX_BAD_DAY_NONE */
 int           numCalendars,  /* (I) Number of calendars. */
 IrxTCalendar const  **calendars,     /* (I) [numCalendars] Calendars. */
 IrxTDate       *outDate)       /* (O) Valid business day. */

{
    static char          routine[] = "GtoMultiHolidayListBusinessDay";
    int                  status = FAILURE;

    IrxTYearMonthDay        ymd1, ymd2;
    IrxTDate                nextDate = 0L;

    REQUIRE (numCalendars > 0);
    REQUIRE (calendars != NULL);

    switch (badDayConv)
    {
    case IRX_BAD_DAY_NONE:
        nextDate = date;
        break;

    case IRX_BAD_DAY_FOLLOW:
        if (getNextBusDateMulti (date, +1, numCalendars, calendars,
                                 &nextDate) != SUCCESS)
            goto RETURN; /* failure */
        break;
        
    case IRX_BAD_DAY_PREVIOUS:
        if (getNextBusDateMulti (date, -1, numCalendars, calendars,
                                 &nextDate) != SUCCESS)
            goto RETURN; /* failure */
        break;

    case IRX_BAD_DAY_MODIFIED:
        /*
        ** Go forwards first. If you wind up in a different
        ** month, then go backwards.
        */
        if (getNextBusDateMulti (date, +1, numCalendars, calendars,
                                 &nextDate) != SUCCESS)
            goto RETURN; /* failure */

        if (irxDateToYMD (nextDate, &ymd1) != SUCCESS)
            goto RETURN;       /* failed */

        if (irxDateToYMD (date, &ymd2) != SUCCESS)
            goto RETURN;       /* failed */

        if (ymd1.month != ymd2.month)
        {
            if (getNextBusDateMulti (date, -1, numCalendars, calendars, 
                                     &nextDate) != SUCCESS)
                goto RETURN; /* failure */
        }
        break;
        
    default:
        irxError ("%s: Unrecognized bad day convention = %ld.\n",
                 routine, badDayConv);
        goto RETURN; /* failure */
    }
    
    *outDate = nextDate;
    status = SUCCESS;

 RETURN:

    if (status != SUCCESS) irxErrorFailure(routine);
    return status;
}

/*
***************************************************************************
** Calculates the number of business days between two dates (FROM & TO).
**
** Algorithm:
** 1. if FROM = TO, the result is 0
** 2. if FROM < TO, the result is the number of business days in the CLOSED
**    interval of [FROM+1,TO]
** 3. if FROM > TO, the result is negated number of business days in the
**    CLOSED interval of [TO,FROM-1]
****************************************************************************/
int irxBusinessDaysDiff
(IrxTDate      fromDate,     /* (I) Earlier date            */
 IrxTDate      toDate,       /* (I) Later date              */
 const IrxTCalendar *cal,          /* (I) Holiday calendar        */
 long       *result)       /* (O) Number of business days */
{
    static char routine[] = "irxBusinessDaysDiff";
    int         status    = FAILURE;

    int         busDaysPerWeek = -1;

    long        nrExtraDays;
    long        numWeeks = 0L;

    int         signum;
    long        numHolidays = 0L;
    long        numBusinessDays;
    IrxTDate      curDate = fromDate;

/*
 * The following tables are used to calculate the number of business days
 * between weekdays.
 * If (toDate > fromDate) use fwdDiffTable.
 * eg NumDays = fwdDiffTable[ fromDate % 7 ][ toDate % 7 ]
 *
 * If (fromDate > toDate) use -bwdDiffTable.
 * eg NumDays = bwdDiffTable[ toDate % 7 ][ fromDate % 7 ]
 *
 * The values of the table were chosen so that the current ALIB behaviour
 * was left unchanged. Note that because fwdDiffTable is not identical
 * to (-1 x bwdDiffTable), swapping fromDate and toDate in the call to this
 * function does not simply reverse the sign of the result.
 */
    static long fwdDiffTable[7][7] =
        
    {/*   Mo Tu We Th Fr Sa Su*/
        
        {  0, 1, 2, 3, 4, 4, 4 },  /* Monday */
        {  4, 0, 1, 2, 3, 3, 3 },  /* Tuesday */
        {  3, 4, 0, 1, 2, 2, 2 },  /* Wednesday */
        {  2, 3, 4, 0, 1, 1, 1 },  /* Thursday */
        {  1, 2, 3, 4, 0, 0, 0 },  /* Friday */
        {  1, 2, 3, 4, 5, 0, 0 },  /* Saturday */
        {  1, 2, 3, 4, 5, 5, 0 }   /* Sunday */
    };

    static long bwdDiffTable[7][7] =
        
    {/*   Mo Tu We Th Fr Sa Su*/
        
        {  0,-1,-2,-3,-4,-5,-5 },  /* Monday */
        { -4, 0,-1,-2,-3,-4,-4 },  /* Tuesday */
        { -3,-4, 0,-1,-2,-3,-3 },  /* Wednesday */
        { -2,-3,-4, 0,-1,-2,-2 },  /* Thursday */
        { -1,-2,-3,-4, 0,-1,-1 },  /* Friday */
        { -0,-1,-2,-3,-4, 0, 0 },  /* Saturday */
        { -0,-1,-2,-3,-4,-5, 0 }   /* Sunday */
    };
    
    REQUIRE (cal != NULL);

    if (fromDate == toDate)
    {
        numBusinessDays = 0;
        goto success;
    }

    /* Set the direction. */
    if (toDate < fromDate)
        signum = -1;
    else
        signum = +1;

    /*
    ** First, compute the date difference adjusting for weekends only.
    */
    switch (cal->weekends)
    {
    case IRX_WEEKEND_STANDARD:
        busDaysPerWeek = 5;
        numWeeks = (toDate - fromDate) / 7;
        curDate += 7 * numWeeks;
        if (curDate > toDate)  /* backwards */
        {
            nrExtraDays = bwdDiffTable [ toDate %7 ][ fromDate % 7 ];
        }
        else
        {
            nrExtraDays = fwdDiffTable [ fromDate %7 ][ toDate % 7 ];
        }
        numBusinessDays = numWeeks * busDaysPerWeek + nrExtraDays;
        break;
    case IRX_WEEKEND_NO_WEEKENDS:
        numBusinessDays = toDate - fromDate;
        break;
    default:
        busDaysPerWeek = 7;
        if (cal->weekends & IRX_WEEKEND_SUNDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_MONDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_TUESDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_WEDNESDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_THURSDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_FRIDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_SATURDAY) busDaysPerWeek--;
        numWeeks    = (toDate - fromDate) / 7;
        curDate    += 7 * numWeeks;
        nrExtraDays = 0;

        while (curDate != toDate)
        {
            curDate += signum;
            if (IRX_IS_WEEKDAY( curDate, cal->weekends ))
            {
                nrExtraDays++;
            }
        }

        numBusinessDays = ((numWeeks * busDaysPerWeek) + nrExtraDays) * signum;
        break;
    }

    /*
    ** Now count the number of weekday holidays
    ** and adjust the date difference by the result. The
    ** call to calcNumWeekdayHolidays searches the
    ** holiday list at most once.
    */

    if (cal->dateList->fNumItems > 0)
    {
        IrxTBool RETURN = FALSE;
        long   holIdx = -1;

        if (findFirstHolidayIdx (fromDate + signum, cal->dateList, signum,
                                 &holIdx, &RETURN) != SUCCESS)
            goto RETURN; /* failure */

        if (!RETURN)
        {
            numHolidays = calcNumWeekdayHolidays (toDate, holIdx,
                                                  signum, cal->weekends,
                                                  cal->dateList, &holIdx);
        }
    }
    
    numBusinessDays -= numHolidays * signum;

 success:
    *result = numBusinessDays;
    status  = SUCCESS;

 RETURN:

    if (status != SUCCESS) irxErrorFailure (routine);
    return status;
}


/*
***************************************************************************
** Calculates a business date being at offset business days
** from the original date
***************************************************************************
*/
int  irxDateAddBusinessDays
(IrxTDate      fromDate,       /* (I) Start date                 */
 long        offset,         /* (I) Number of business days    */
 const IrxTCalendar *cal,            /* (I) Holiday calendar           */
 IrxTDate     *result)         /* (O) resulting business date    */
{

    static char routine[] = "irxDateAddBusinessDays";
    int         status = FAILURE;

    long        intervalSign;
    long        numBusDays;
    long        busDaysPerWeek;
    IrxTDate      toDate;

    if (offset >= 0)
    {
        intervalSign = +1;
        numBusDays   = offset;
    }
    else
    {
        intervalSign = -1;
        numBusDays   = -offset;
    }

    REQUIRE (cal != NULL);

    /*
    ** Get the number of business days per week. In-line for speed.
    */
    switch (cal->weekends)
    {
    case IRX_WEEKEND_STANDARD:
        busDaysPerWeek = 5;
        break;
    case IRX_WEEKEND_NO_WEEKENDS:
        busDaysPerWeek = 7;
        break;
    default:
        busDaysPerWeek = 7;
        if (cal->weekends & IRX_WEEKEND_SUNDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_MONDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_TUESDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_WEDNESDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_THURSDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_FRIDAY) busDaysPerWeek--;
        if (cal->weekends & IRX_WEEKEND_SATURDAY) busDaysPerWeek--;
        break;
    }

    if (cal->dateList->fNumItems == 0 &&
        cal->weekends == IRX_WEEKEND_NO_WEEKENDS)
    {
        /*
        ** No adjustments at all.
        */
        toDate = fromDate + (intervalSign * numBusDays);
    }
    else if (intervalSign > 0 && cal->weekends == IRX_WEEKEND_STANDARD)
    {
        /* First, check for special case. The code for this case
         * was provided by Tom Robbins-Milne. It is much faster
         * than the general algorithm.
         */

        if (numBusDays == 0)
        {
            toDate = fromDate;
        }
        else
        {
           /*
            * Calculate result if no holidays are involved.
            *
            * Move forward over a week for every 5 business days.
            * Use a table for moving 0..4 bus days from each day of the week.
            */
            toDate = fromDate +  7 * (numBusDays / 5) +
                g_offsetTable[ fromDate % 7 ][ numBusDays % 5 ];
            
            /*
             * Now adjust for any holidays
             */
            if (cal->dateList->fNumItems > 0)
            {
                if (holAdjustStandardWeekends (fromDate,
                                               cal->dateList,
                                               &toDate) != SUCCESS)
                    goto RETURN;
            }
        }
    }
    else
    {
        if (forwardNonStandardWeekends (fromDate,
                                        numBusDays,
                                        intervalSign,
                                        cal->weekends,
                                        busDaysPerWeek,
                                        cal->dateList,
                                        &toDate) != SUCCESS)
            goto RETURN;
    }

    *result = toDate;
    status  = SUCCESS;
    
 RETURN:
    
    if (status != SUCCESS) irxErrorFailure(routine);
    return status;
}


/*
***************************************************************************
** Computes the last business day of the month.
***************************************************************************
*/
int irxDateToBusinessEOM
(IrxTDate      inDate,  /* (I) Input date */
 const IrxTCalendar *cal,     /* (I) Holiday calendar */
 IrxTDate     *result   /* (O) Last business day of the month */
)
{
    /*
    ** Calculate the last day of the month.
    ** Adjust backwards for holidays.
    */

    static char routine[] = "irxDateToBusinessEOM";
    int         status    = FAILURE;

    IrxTDate      outDate;

    REQUIRE (cal != NULL);
    REQUIRE (result != NULL);

    outDate = irxDateToEndOfMonth (inDate);
    if (irxBusinessDay (outDate,
                       IRX_BAD_DAY_PREVIOUS,
                       cal,
                       &outDate) != SUCCESS)
        goto RETURN; /* failure */

    *result = outDate;
    status  = SUCCESS;

 RETURN:

    if (status != SUCCESS) irxErrorFailure (routine);
    return status;
}

/*f
***************************************************************************
** Returns the calendar from a name using the calendar cache.
**
** You can use name = "NONE" or name = "NO_WEEKENDS" and the resulting
** calendar list will be defined appropriately.
**
** Returns NULL on failure, a valid IRXCalendar pointer on success. Does
** not report errors.
***************************************************************************
*/
IrxTCalendar* irxCalendarCacheSearch
(const char  *name          /* (I) Name of calendar */
)
{
    IrxTCalendar      *cal   = NULL;

    if (name == NULL) return NULL;
    IRX_HASH_SEARCH(calendarCacheGet(), cacheNameGet(name), &cal);
    return cal;
}

/*f
***************************************************************************
** Adds a calendar list to the calendar cache. If the entry already exists
** in the cache, then the old version will be deleted.
***************************************************************************
*/
int irxCalendarCacheAdd
(char       *name,    /* (I) Name to associate calendars with */
 IrxTCalendar *calendar /* (I) Adds shallow copy */
)
{
    static char routine[] = "irxCalendarCacheAdd";
    int         status    = FAILURE;

    IrxTHashTable* cache = calendarCacheGet();
    IrxTCalendar*  old;
    char* cacheName;

    REQUIRE (name != NULL);
    REQUIRE (calendar != NULL);
    REQUIRE (cache != NULL);

    cacheName = cacheNameGet(name);

    while (IRX_HASH_SEARCH(cache, cacheName, &old) == SUCCESS)
        IRX_HASH_DELETE_ENTRY(cache, cacheName);

    status = IRX_HASH_INSERT(cache, cacheName, calendar);

 RETURN:
    
    if (status != SUCCESS) irxErrorFailure(routine);

    return status;
}


/*f
***************************************************************************
** Reads a holiday file into memory as a calendar.
**
** Structure of holiday file, ascii text file of lines like:
**     #            - commment (blank lines are also skipped)
**     19631026     - means 10/26/1963 is a holiday
**     # SATURDAY_NOT_ALWAYS_HOLIDAY  - sets "saturday isn't always a holiday"
**     # SUNDAY_NOT_ALWAYS_HOLIDAY    - sets "sunday isn't always a holiday"
**     # MONDAY_ALWAYS_HOLIDAY        - sets "monday as always a holiday"
**     # TUESDAY_ALWAYS_HOLIDAY       - sets "tuesday as always a holiday"
**     # WEDNESDAY_ALWAYS_HOLIDAY     - sets "wednesday as always a holiday"
**     # THURDSAY_ALWAYS_HOLIDAY      - sets "thursday as always a holiday"
**     # FRIDAY_ALWAYS_HOLIDAY        - sets "friday as always a holiday"
**
** Dates must be in increasing order.
***************************************************************************
*/
IrxTCalendar* irxCalendarRead
(const char *fileName   /* (I) Name of file to read (may differ) */
)
{
    static char routine[] = "irxCalendarRead";

    char        buffer[512];         /* line of text from file */
    FILE       *fp = NULL;           /* file handle */
    IrxTCalendar *cal = NULL;          /* to be returned */
    IrxTDateList *dl = NULL;           /* date list read in */
    long        weekends;

    IrxTDate      dates[32];
    int         numDates = 0;
    int         maxDates = sizeof(dates)/sizeof(IrxTDate);

    fp = fopen(fileName, "r");
    if (fp == NULL)
    {
        irxError ("%s: Couldn't open file %s.\n", routine, fileName);
        goto RETURN; /* failure */
    }

    weekends = IRX_WEEKEND_STANDARD; /* Includes SAT and SUN */

    while (fgets (buffer, sizeof(buffer), fp) != NULL)
    {
        if (buffer[0] == '#')
        {
            static char monString[] = "# MONDAY_ALWAYS_HOLIDAY";
            static char tueString[] = "# TUESDAY_ALWAYS_HOLIDAY";
            static char wedString[] = "# WEDNESDAY_ALWAYS_HOLIDAY";
            static char thuString[] = "# THURSDAY_ALWAYS_HOLIDAY";
            static char friString[] = "# FRIDAY_ALWAYS_HOLIDAY";
            static char satString[] = "# SATURDAY_NOT_ALWAYS_HOLIDAY";
            static char sunString[] = "# SUNDAY_NOT_ALWAYS_HOLIDAY";

            if (strncmp (buffer, monString, sizeof(monString)-1) == 0)
                weekends |= IRX_WEEKEND_MONDAY;

            if (strncmp (buffer, tueString, sizeof(tueString)-1) == 0)
                weekends |= IRX_WEEKEND_TUESDAY;

            if (strncmp (buffer, wedString, sizeof(wedString)-1) == 0)
                weekends |= IRX_WEEKEND_WEDNESDAY;

            if (strncmp (buffer, thuString, sizeof(thuString)-1) == 0)
                weekends |= IRX_WEEKEND_THURSDAY;

            if (strncmp (buffer, friString, sizeof(friString)-1) == 0)
                weekends |= IRX_WEEKEND_FRIDAY;

            if (strncmp (buffer, satString, sizeof(satString)-1) == 0)
                weekends &= ~IRX_WEEKEND_SATURDAY;

            if (strncmp (buffer, sunString, sizeof(sunString)-1) == 0)
                weekends &= ~IRX_WEEKEND_SUNDAY;
        }
        else
        {
            long dateVal = atol (buffer);
            if (dateVal > 16010101)
            {
                long year   = dateVal / 10000;
                long month  = (dateVal % 10000) / 100;
                long day    = dateVal % 100;
                IrxTDate date = irxDate (year, month, day);
                
                if (date == FAILURE)
                {
                    irxError ("%s: invalid date: %s", routine, buffer);
                    goto RETURN; /* failure */
                }

                dates[numDates] = date;
                ++numDates;

                assert (numDates <= maxDates);
                if (numDates == maxDates)
                {
                    IrxTDateList* tmp = irxDateListAddDates (dl, numDates, dates);
                    if (tmp == NULL) goto RETURN; /* failure */
                    irxDateListFree (dl);
                    dl = tmp;
                    numDates = 0;
                }
            }
        }
    }

    if (numDates > 0)
    {
        IrxTDateList* tmp = irxDateListAddDates (dl, numDates, dates);
        if (tmp == NULL) goto RETURN; /* failure */
        irxDateListFree (dl);
        dl = tmp;
    }

    if (dl == NULL && weekends == IRX_WEEKEND_STANDARD)
    {
        irxError ("%s: No holiday information found in %s.\n",
                 routine, fileName);
        irxError ("   Either week-end information or dates must be "
                 "provided.\n");
        goto RETURN; /* failure */
    }

    cal = irxCalendarMake (dl, weekends);

RETURN:

    if (cal == NULL) irxErrorFailure (routine);

    if (fp != NULL) fclose(fp);
    irxDateListFree (dl);
    return cal;
}



/*f
***************************************************************************
** Deletes all cache entries in hash table.  Useful for being fastidious
** about memory leaks.
***************************************************************************
*/
void irxCalendarCacheClear (void)
{
    if (g_calendarCache != NULL)
    {
        irxHashTableFree(g_calendarCache);
        g_calendarCache = NULL;
    }
}

/*
***************************************************************************
** Static functions
***************************************************************************
*/

/*
** Gets the calendar cache hash table - creating it if necessary
*/
static IrxTHashTable* calendarCacheGet(void)
{
    if (g_calendarCache == NULL)
    {
        g_calendarCache = irxHashTableMakeStringKey (
            101, FALSE,
            (IrxTHashDataFreeFunc*)irxCalendarFree);
        IRX_HASH_INSERT(g_calendarCache, 
                       "NONE",
                       irxCalendarMake(NULL, IRX_WEEKEND_STANDARD));
        IRX_HASH_INSERT(g_calendarCache,
                       "NO_WEEKENDS",
                       irxCalendarMake(NULL, IRX_WEEKEND_NO_WEEKENDS));
    }
    return g_calendarCache;
}

/*
** Returns the name entry as used within the cache
*/
static char* cacheNameGet(const char *name)
{
    static char buf[256];
    size_t len = strlen(name);
    size_t i;
    if (len >= sizeof(buf)) 
        len = sizeof(buf)-1;

    for (i = 0; i < len; ++i) buf[i] = toupper(name[i]);
    buf[len] = '\0';
    return buf;
}


/*
***************************************************************************
** Compute the date of next business day, given a start date and
** an array of calendars.
***************************************************************************
*/
static int getNextBusDateMulti
(IrxTDate         startDate,        /* (I) starting date. */
 long           direction,        /* (I) +1=forwards, -1=backwards */
 int            numCalendars,     /* (I) Size of calendars array. */
 IrxTCalendar const   **calendars,        /* (I) Calendars. */
 IrxTDate        *nextDate          /* (O) next business day */
)
{
    static char          routine[] = "getNextBusDateMulti";
    int                  status = FAILURE;
    int                  n;
    IrxTDate               adjDate;

    *nextDate = startDate;

    /* this loop executes getNextBusDate for each market in sequence; 
       whenever weekends/holidays force an adjustment, the process recommences
       at the first market; eventually, a date will be reached which is a 
       business day in all markets, whereupon the loop will execute to 
       completion

       we can be certain that all of the skipped dates were weekends/holidays
       in at least one of the markets, i.e. output will be the _first_ 
       suitable day */

    n = 0;
    while (n < numCalendars)
    {
        if ( getNextBusDate(*nextDate, direction, calendars[n], 
                            &adjDate) != SUCCESS )
        {
            goto RETURN;
        }
        
        if ( adjDate != *nextDate )
        {
            *nextDate = adjDate;

            /* date has changed to acommodate this market; now
               have to recommence testing for all other markets */

            /* in fact, this date does not need retesting for the current
               market, and we exploit this in the special case below, which
               avoids a redundant duplicate test in "single holiday list"
               usage */

            n = ( n == 0 ) ? 1 : 0;
            
            /* ^-- if adjusted for first market, now move on to second; in all
               other cases, recommence testing from scratch */
        }
        else
        {
            n++; /* next market */
        }
    }

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS) irxErrorFailure(routine);
    return status;
}

/*
***************************************************************************
** Compute the date of next business day, given a start date and
** a calendar.
***************************************************************************
*/
static int getNextBusDate
(IrxTDate      startDate,    /* (I) starting date. */
 long        direction,    /* (I) +1=forwards, -1=backwards */
 const IrxTCalendar *cal,          /* (I) calendar */
 IrxTDate     *nextDate      /* (O) next business day */
)
{
    static char routine[] = "getNextBusDate";
    int         status    = FAILURE;

    IrxTBool               RETURN;
    IrxTBool               aHoliday;
    long                 holIdx;

    long                 numHols = cal->dateList->fNumItems;
    IrxTDate              *holArray = cal->dateList->fArray;
    IrxTDate               curDate = startDate;

    if (numHols == 0)
    {
        RETURN = TRUE;
    }
    else
    {
        if (findFirstHolidayIdx (startDate, cal->dateList, direction,
                                 &holIdx, &RETURN) != SUCCESS)
            goto RETURN; /* failure */
    }

    /*
    ** holIdx should now point to the first holiday on or before
    ** (after) the start date. Now loop though each day until
    ** you reach the first business day.
    */

    aHoliday = TRUE;

    do
    {   /* Check if it's a holiday first. */
        if (! RETURN  &&
            curDate == holArray[ holIdx ])
        {
            holIdx  += direction;
            curDate += direction;
            RETURN =
                (IrxTBool)(holIdx < 0 || holIdx >= numHols);
        }
        else
        {
            /* Day is not a holiday. Check if it's a week-end day. */
            if (IRX_IS_WEEKEND (curDate, cal->weekends))
            {
                curDate += direction;
            }
            else
            {
                aHoliday = FALSE;
            }
        }
    }
    while (aHoliday);
    
    *nextDate = curDate;
    status    = SUCCESS;

 RETURN:

    if (status != SUCCESS) irxErrorFailure(routine);
    return status;
}

/*
***************************************************************************
** This function does the following:
**
** If going forward in time, it finds the index of the 1st holiday date
** which is >= the input date.
**
** If going backwards in time, it finds the index of the first of holiday
** which is <= the input date.
**
** The function calls GtoBSearchLong to do the search. It returns SUCCESS
** or FAILURE.
**
** Note that the function DOES NOT check if the input holiday list is
** empty or NULL.
***************************************************************************
*/
static int  findFirstHolidayIdx (
    IrxTDate          date,         /* (I) Date to search for. */
    const IrxTDateList     *calDateList,  /* (I) Holiday date list. */
    long            direction,    /* (I) +1=forwards, -1=backwards. */
    long           *holIdx,       /* (O) Index of first holiday. */
    IrxTBool         *doneSearching /* (O) TRUE=date is beyond list bounds. */
)
{
    static char routine[] = "findFirstHolidayIdx";
    int         status    = FAILURE;

    long      idxExact;
    long      idxLo;
    long      idxHi;

    /*
    ** Do a binary search of the date list to find the
    ** bracketing low and high index values.
    */
    if (irxBinarySearchLong (date,
                            calDateList->fArray,
                            sizeof(IrxTDate),
                            calDateList->fNumItems,
                            &idxExact,
                            &idxLo,
                            &idxHi) != SUCCESS)
        goto RETURN; /* failure */
    
    if (direction > 0)
    {
        /* going forward in time */
        /* return index of first holiday which is >= input date */
        if (idxExact >= 0)
        {
            *holIdx = idxExact;
            *doneSearching = FALSE;
        }
        else
        {
            *holIdx = idxHi;
            *doneSearching = idxHi >= calDateList->fNumItems;
        }
    }
    else
    {
        /* going backward in time */
        /* return index of first holiday which is <= input date */
        if (idxExact >= 0)
        {
            *holIdx = idxExact;
            *doneSearching = FALSE;
        }
        else
        {
            *holIdx = idxLo;
            *doneSearching = idxLo < 0;
        }
    }

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS) irxErrorFailure(routine);
    return status;
}

/*f
***************************************************************************
** Handles the special case where the number of business days to forward
** by is very large.
***************************************************************************
*/
static int forwardNonStandardWeekends (
    IrxTDate         fromDate,       /* (I) start date */
    long           numBusDays,     /* (I) abs. num. bus. days */
    long           direction,      /* (I) +1=forwards, -1=backwards */
    long           weekends,       /* (I) Weekends flags */
    long           busDaysPerWeek, /* (I) num. bus. days per week */
    const IrxTDateList    *calDateList,    /* (I) holiday dateList */
    IrxTDate        *resultDate      /* (O) forwarded date. */
)
{
    static char routine[] = "forwardNonStandardWeekends";
    int         status = FAILURE;

    long        numWeeks;
    long        numHolidays = 0L;
    IrxTDate      curDate = fromDate;
    long        holIdx = -1;
    IrxTBool      doneSearchingList = FALSE;

    IrxTDate     *holArray = calDateList->fArray;
    long        numHols  = calDateList->fNumItems;

    /*
    ** First, adjust for weekends only, pretending there are no holidays.
    */
    numWeeks  = MAX(((numBusDays / busDaysPerWeek) - 1), 0);
    curDate  += (7 * numWeeks * direction);

    if (numHols > 0)
    {
        /*
        ** Search the holiday list for the first holiday
        ** strictly after (if going forward in time) or strictly
        ** before (if going backward in time) the start date. Note
        ** that the holiday list is assumed to be sorted in increasing
        ** order.
        */
        if (findFirstHolidayIdx (fromDate + direction,
                                 calDateList,
                                 direction,
                                 &holIdx,
                                 &doneSearchingList) != SUCCESS)
        {
            goto RETURN; /* failure */
        }
    }
    else
    {
        doneSearchingList = TRUE; /* No holidays - no searching */
    }

    /*
    ** Count the number of weekday holidays starting from the
    ** current holiday index and going up to the curDate.
    ** Count weekday holidays only because holidays occurring
    ** on week-end days have been skipped by the previous calculation.
    */
    if (! doneSearchingList)
    {
        numHolidays = calcNumWeekdayHolidays (curDate,
                                              holIdx,
                                              direction,
                                              weekends,
                                              calDateList,
                                              &holIdx);
        doneSearchingList = (holIdx  <  0  ||  holIdx  >=  numHols);
    }
    else
    {
        numHolidays = 0;
    }
    numBusDays -= ((busDaysPerWeek * numWeeks) - numHolidays);

    /*
    ** Now search one day at a time starting one day beyond
    ** the current date.
    */
    while (numBusDays > 0)
    {
        curDate += direction;

        /*
        ** Check if the day is a holiday first. If it is, don't decrement
        ** numBusDays.
        */
        if ( ! doneSearchingList  &&
             curDate  ==  holArray[ holIdx ] )
        {
            holIdx  += direction;
            doneSearchingList = (holIdx  <  0  ||  holIdx  >=  numHols);
        }
        else  /* Not a holiday. */
        {
            /*
            ** If the day is a weekday, then decrement numBusDays,
            ** otherwise continue looping.
            */
            if (IRX_IS_WEEKDAY( curDate, weekends ))
            {
                numBusDays--;
            }
        }
    }

    *resultDate = curDate;
    status = SUCCESS;

 RETURN:
    
    if (status != SUCCESS) irxErrorFailure(routine);
    return status;
}


/*
***************************************************************************
** This function adjusts a forwarded date to take account of any
** holidays which might have occured during the forward calculation.
** The algorithm assumes the following conditions:
**
** (1) Saturdays and Sundays are both holidays.
** (2) The direction is forwards in time.
**
** The function is optimized for speed.
***************************************************************************
*/
static int   holAdjustStandardWeekends (
    IrxTDate           fromDate,       /* (I) start date */
    const IrxTDateList     * calDateList,    /* (I) holiday dateList */
    IrxTDate         * resultDate      /* (I/O) forwarded date. */
)
{
    static char routine[] = "holAdjustStandardWeekends";
    int         status    = FAILURE;

    IrxTDate      dateEnd;      /* Date to return */
    int         i;

    long        idxExact;
    long        idxLo;
    long        idxHi;

    int         numHols  = calDateList->fNumItems;
    IrxTDate     *holDates = calDateList->fArray;

    dateEnd = *resultDate;

    /*
    ** If there are any holidays, then we want to find the first holiday
    ** strictly after fromDate. 
    **
    ** As a special case, if there are no holidays or fromDate >= last
    ** holiday, then no further adjustments are necessary.
    */
    if (numHols <= 0 ||
        fromDate >= holDates[numHols-1])
    {
        /* Special case as above */
    }
    else
    {
        /*
         * Find date in date list using a binary search - this routine
         * returns idxHi as the first date strictly after fromDate
         */
        if (irxBinarySearchLong (fromDate,
                                holDates,
                                sizeof(IrxTDate),
                                numHols,
                                &idxExact,
                                &idxLo,
                                &idxHi) != SUCCESS)
            goto RETURN; /* failure */
        
        /* idxHi will be the index of the first date > fromDate */
        /* if there is no such date, then idxHi = numHols */
        i = idxHi;

        while (i < numHols && holDates[i] <= dateEnd)
        {
            dateEnd += g_offsetTable[ dateEnd % 7 ][ 1 ];
            i++;
        }
    }

    *resultDate = dateEnd;

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS) irxErrorFailure(routine);
    return status;
}


/*
***************************************************************************
** Compute the number of week-day holidays starting at the input
** holiday index and continuing until the holiday list is at
** or beyond toDate.
**
** Week-day holidays are holidays which occur on a non-weekend day.
**
** Since when we construct the calendar we remove weekdays, we do not
** actually test for weekends, but assume that everything is a weekday
** already.
***************************************************************************
*/
static long  calcNumWeekdayHolidays (
    IrxTDate       toDate,      /* (I) End date.   */
    long         startHolIdx, /* (I) Starting holiday index. */
    long         direction,   /* (I) 1=forwards, 2=backwards*/
    long         weekends,    /* (I) Weekend flags */
    const IrxTDateList  *calDateList, /* (I) Holiday datelist. */
    long        *endHolIdx)    /* (O) Idx where hol=toDate, +1 */
{
    long       k;
    
    long       numHols   = 0L;
    IrxTDate    *dateArray = calDateList->fArray;
    long       numItems  = calDateList->fNumItems;

    if (direction == 1)     /* forward */
    {
        for (k = startHolIdx;
             (k < numItems && dateArray[k] <= toDate);
             k++)
        {
            numHols++;
        }
        *endHolIdx = k;
        if (k < numItems   &&
            dateArray[k] == toDate)
        {
            (*endHolIdx)++;
        }
    }
    else                   /* backward */
    {
        for (k = startHolIdx;
             (k >= 0 && dateArray[k] >= toDate);
             k--)
        {
            numHols++;
        }
        *endHolIdx = k;
        if (k >= 0  &&
            dateArray[k] == toDate)
        {
            (*endHolIdx)--;
        }
    }
    return numHols;
}

