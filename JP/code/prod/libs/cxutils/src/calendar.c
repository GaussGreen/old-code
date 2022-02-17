/*
***************************************************************************
** SOURCE FILE: calendar.c
**
** Implements the calendar object, the calendar cache (so we can access
** calendar objects by name) and basic business day functions.
**
** Adapts the ALIB holiday functions and improves the interface.
**
** $Header$
***************************************************************************
*/

#include "calendar.h"

#include <assert.h>
#include <ctype.h>
#include <stdio.h>

#include "bsearch.h"             /* CxBinarySearchLong */
#include "date.h"
#include "datelist.h"            /* CxDateList... */
#include "hash.h"                /* CX_HASH_... */
#include "cxmacros.h" 
#include "strutils.h"

#include <alib/busday.h>

#define CX_WEEKEND_MIN_VALUE   0x0000
#define CX_WEEKEND_MAX_VALUE   0x007f

#define DEFAULT_CALENDAR CxCalendarCacheSearch("NONE")

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
char* CxWeekendsToString (long weekends)
{
    static char str[32];

    strcpy(str, "");

    if (weekends < CX_WEEKEND_MIN_VALUE || weekends > CX_WEEKEND_MAX_VALUE)
    {
        strcat (str, "***ERROR***");
    }
    else
    {
        if (weekends & GTO_WEEKEND_MONDAY)    strcat(str, "MON,");
        if (weekends & GTO_WEEKEND_TUESDAY)   strcat(str, "TUE,");
        if (weekends & GTO_WEEKEND_WEDNESDAY) strcat(str, "WED,");
        if (weekends & GTO_WEEKEND_THURSDAY)  strcat(str, "THU,");
        if (weekends & GTO_WEEKEND_FRIDAY)    strcat(str, "FRI,");
        if (weekends & GTO_WEEKEND_SATURDAY)  strcat(str, "SAT,");
        if (weekends & GTO_WEEKEND_SUNDAY)    strcat(str, "SUN,");
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
int CxWeekendsFromString (char *str, long *weekends)
{
    static char routine[] = "CxWeekendsFromString";
    int         status    = FAILURE;

    size_t  i;
    size_t  numItems;
    char  **split = NULL;
    long    myWeekends = 0;

    REQUIRE(str != NULL);
    REQUIRE(weekends != NULL);

    if (strlen(str) == 0)
        goto success;

    if (CxStringSplit (str, ',', &numItems, &split) != SUCCESS)
        goto done; /* failure */

    for (i = 0; i < numItems; ++i)
    {
        if (CxStrcmpi (split[i], "MON") == 0)
            myWeekends |= GTO_WEEKEND_MONDAY;
        else if (CxStrcmpi (split[i], "TUE") == 0)
            myWeekends |= GTO_WEEKEND_TUESDAY;
        else if (CxStrcmpi (split[i], "WED") == 0)
            myWeekends |= GTO_WEEKEND_WEDNESDAY;
        else if (CxStrcmpi (split[i], "THU") == 0)
            myWeekends |= GTO_WEEKEND_THURSDAY;
        else if (CxStrcmpi (split[i], "FRI") == 0)
            myWeekends |= GTO_WEEKEND_FRIDAY;
        else if (CxStrcmpi (split[i], "SAT") == 0)
            myWeekends |= GTO_WEEKEND_SATURDAY;
        else if (CxStrcmpi (split[i], "SUN") == 0)
            myWeekends |= GTO_WEEKEND_SUNDAY;
        else
        {
            GtoErrMsg ("%s: Invalid date string %s\n", routine, split[i]);
            goto done; /* failure */
        }
    }

 success:
    *weekends = myWeekends;
    status = SUCCESS;

 done:

    FREE(split);
    if (status != SUCCESS)
        GtoErrMsgFailure(routine);

    return status;
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
int CxCalendarValidate
(CxTCalendar *cal)     /* (I/O) Calendar to be adjusted */
{
    static char   routine[] = "CxCalendarValidate";
    int           status = FAILURE;

    long          i1;
    long          i2;

    REQUIRE (cal != NULL);
    REQUIRE (cal->weekends >= CX_WEEKEND_MIN_VALUE);
    REQUIRE (cal->weekends <= CX_WEEKEND_MAX_VALUE);

    if (cal->weekends == CX_WEEKEND_MAX_VALUE)
    {
        GtoErrMsg ("%s: All weekdays are defined as weekends\n", routine);
        goto done; /* failure */
    }

    if (cal->dateList == NULL)
    {
        cal->dateList = GtoNewDateListFromDates (NULL, 0);
        if (cal->dateList == NULL) 
            goto done; /* failure */
    }

    /*
    ** Remove weekends from the date list
    */
    i1 = 0;
    i2 = 0;
    while (i1 < cal->dateList->fNumItems)
    {
        TDate thisDate = cal->dateList->fArray[i1];
        if (GTO_IS_WEEKEND(thisDate, cal->weekends))
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

    /* CxDateListSort also removes duplicates */
    CxDateListSort(cal->dateList);
    status = SUCCESS;

 done:
    
    if (status != SUCCESS) GtoErrMsgFailure(routine);
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
int CxBusinessDay
(TDate        date,       /* (I) Arbitrary date. */
 CxTBadDayConv badDayConv, /* (I) Bad day convention. */
 CxTCalendar  *cal,        /* (I) Calendar. */
 TDate       *outDate)    /* (O) Valid business day. */

{
    static char routine[] = "CxBusinessDay";
    int         status = FAILURE;

    if (cal == NULL) cal = DEFAULT_CALENDAR;

    REQUIRE (cal != NULL);
    REQUIRE (outDate != NULL);

    if (GtoHolidayListBusinessDay (date, 
                                   badDayConv,
                                   (THolidayList*)cal,
                                   outDate) != SUCCESS)
        goto done; /* failure */
    
    status = SUCCESS;

 done:

    if (status != SUCCESS) GtoErrMsgFailure(routine);
    return status;
}


/*f
***************************************************************************
** Returns the calendar from a name using the calendar cache.
**
** You can use name = "NONE" or name = "NO_WEEKENDS" and the resulting
** calendar list will be defined appropriately.
**
** Returns NULL on failure, a valid CxTCalendar pointer on success. Does
** not report errors.
***************************************************************************
*/
CxTCalendar* CxCalendarCacheSearch
(char  *name          /* (I) Name of calendar */
)
{
    CxTCalendar *cal   = NULL;

    if (GtoIsHolidayLoaded (name))
    {
        cal = (CxTCalendar*) GtoHolidayListFromCache(name);
    }
    else
    {
        cal = NULL;
    }
    return cal;
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
CxTCalendar* CxCalendarRead
(char *fileName   /* (I) Name of file to read (may differ) */
)
{
    static char routine[] = "CxCalendarRead";

    CxTCalendar *cal = NULL;          /* to be returned */

    REQUIRE (fileName != NULL);

    cal = (CxTCalendar*)GtoHolidayListRead (fileName);

done:

    if (cal == NULL) GtoErrMsgFailure (routine);
    return cal;
}


/*
** This section includes functions from ALIB which have been wrapped but
** are not currently in use. Hence this section is not being compiled, but
** kept to hand in case we need it later in a hurry!
*/

#if 0
/*
***************************************************************************
** Indicates whether a date is a business day according to a calendar.
***************************************************************************
*/
int CxIsBusinessDay
(TDate       date,               /* (I) Input Date */
 CxTCalendar *cal,                /* (I) Holiday list structure */
 TBoolean   *isBusinessDay)      /* (O) TRUE if a business day */
{
    static char routine[] = "CxIsBusinessDay";
    int         status = FAILURE;

    if (cal == NULL) cal = DEFAULT_CALENDAR;

    REQUIRE (cal != NULL);
    REQUIRE (isBusinessDay != NULL);

    if (GtoHolidayListIsBusinessDay (date,
                                     (THolidayList*)cal,
                                     isBusinessDay) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure(routine);
        
    return status;
}

/*
***************************************************************************
** Checks to see if a date is a holiday. Does not take week-ends into
** account.
***************************************************************************
*/
int CxIsHoliday
(TDate       date,      /* (I) Arbitrary date      */
 CxTCalendar *cal,       /* (I) Calendar structure  */
 TBoolean   *isHoliday) /* (O) 0 = not a holiday, 1 = is a holiday */
{
    static char routine[] = "CxIsHoliday";
    int         status = FAILURE;

    if (cal == NULL) cal = DEFAULT_CALENDAR;

    REQUIRE (cal != NULL);
    REQUIRE (cal->dateList != NULL);

    if (GtoHolidayListIsHoliday (date,
                                 (THolidayList*)cal,
                                 isHoliday) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS) GtoErrMsgFailure(routine);
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
int CxBusinessDayMultiCalendar
(TDate        date,          /* (I) Arbitrary date. */
 CxTBadDayConv badDayConv,    /* (I) Bad day convention. */
 int          numCalendars,  /* (I) Number of calendars. */
 CxTCalendar **calendars,     /* (I) [numCalendars] Calendars. */
 TDate       *outDate)       /* (O) Valid business day. */

{
    static char routine[] = "CxBusinessDayMultiCalendar";
    int         status = FAILURE;

    REQUIRE (numCalendars > 0);
    REQUIRE (calendars != NULL);

    if (GtoMultiHolidayListBusinessDay (date,
                                        badDayConv,
                                        numCalendars,
                                        (THolidayList**)calendars,
                                        outDate) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS) GtoErrMsgFailure(routine);
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
int CxBusinessDaysDiff
(TDate       fromDate,     /* (I) Earlier date            */
 TDate       toDate,       /* (I) Later date              */
 CxTCalendar *cal,          /* (I) Holiday calendar        */
 long       *result)       /* (O) Number of business days */
{
    static char routine[] = "CxBusinessDaysDiff";
    int         status    = FAILURE;

    if (cal == NULL) cal = DEFAULT_CALENDAR;

    REQUIRE (cal != NULL);
    REQUIRE (result != NULL);

    if (GtoHolidayListBusinessDaysDiff (fromDate,
                                        toDate,
                                        (THolidayList*)cal,
                                        result) != SUCCESS)
        goto done; /* failure */

    status  = SUCCESS;

 done:

    if (status != SUCCESS) GtoErrMsgFailure (routine);
    return status;
}

/*
***************************************************************************
** Calculates a business date being at offset business days
** from the original date
***************************************************************************
*/
int CxDateAddBusinessDays
(TDate       fromDate,       /* (I) Start date                 */
 long        offset,         /* (I) Number of business days    */
 CxTCalendar *cal,            /* (I) Holiday calendar           */
 TDate      *result)         /* (O) resulting business date    */
{
    static char routine[] = "CxDateAddBusinessDays";
    int         status = FAILURE;

    if (cal == NULL) cal = DEFAULT_CALENDAR;

    REQUIRE (cal != NULL);
    REQUIRE (result != NULL);

    if (GtoHolidayListAddBusinessDays (fromDate,
                                       offset,
                                       (THolidayList*)cal,
                                       result) != SUCCESS)
        goto done; /* failure */

    status  = SUCCESS;
    
 done:
    
    if (status != SUCCESS) GtoErrMsgFailure(routine);
    return status;
}

/*
***************************************************************************
** Computes the last business day of the month.
***************************************************************************
*/
int CxDateToBusinessEOM
(TDate       inDate,  /* (I) Input date */
 CxTCalendar *cal,     /* (I) Holiday calendar */
 TDate      *result   /* (O) Last business day of the month */
)
{
    /*
    ** Calculate the last day of the month.
    ** Adjust backwards for holidays.
    */
    static char routine[] = "CxDateToBusinessEOM";
    int         status    = FAILURE;

    if (cal == NULL) cal = DEFAULT_CALENDAR;

    REQUIRE (cal != NULL);
    REQUIRE (result != NULL);

    if (GtoHolidayListDateToBusinessEOM (inDate,
                                         (THolidayList*)cal,
                                         result) != SUCCESS)
        goto done; /* failure */

    status  = SUCCESS;

 done:

    if (status != SUCCESS) GtoErrMsgFailure (routine);
    return status;
}

/*f
***************************************************************************
** Adds a calendar list to the calendar cache. If the entry already exists
** in the cache, then the old version will be deleted.
***************************************************************************
*/
int CxCalendarCacheAdd
(char       *name,    /* (I) Name to associate calendars with */
 CxTCalendar *calendar /* (I) Adds shallow copy */
)
{
    static char routine[] = "CxCalendarCacheAdd";
    int         status    = FAILURE;

    REQUIRE (name != NULL);
    REQUIRE (calendar != NULL);

    if (GtoHolidayListAddToCache (name, (THolidayList*)calendar) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:
    
    if (status != SUCCESS) GtoErrMsgFailure(routine);

    return status;
}

/*f
***************************************************************************
** Deletes all cache entries in hash table.  Useful for being fastidious
** about memory leaks.
***************************************************************************
*/
void CxCalendarCacheClear (void)
{
    GtoHolidayEmptyCache();
}


#endif
