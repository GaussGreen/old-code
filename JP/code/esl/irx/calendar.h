/*
***************************************************************************
** HEADER FILE: calendar.h
**
** Defines the calendar object, the calendar cache (so we can access
** calendar objects by name) and basic business day functions.
**
** $Header$
***************************************************************************
*/

#ifndef IRX_CALENDAR_H
#define IRX_CALENDAR_H

#include "irxutils.h"    /* basic data types */

#ifdef __cplusplus
extern "C"
{
#endif
/*
** These constants define flags for determining whether a day is a weekend
** or not. This is to cope with countries for which weekends fall on Fridays
** and Saturdays etc.
**
** Note the holiday file "NONE" corresponds to an empty date list and
** IRX_WEEKEND_STANDARD, whereas the holiday file "NO_WEEKENDS" corresponds
** to an empty date list and IRX_WEEKEND_NO_WEEKENDS.
**
** By default when reading holidays from file, the weekend status will be
** IRX_WEEKEND_STANDARD unless specifically over-ridden.
*/
#define IRX_WEEKEND_SUNDAY      0x0001
#define IRX_WEEKEND_MONDAY      0x0002
#define IRX_WEEKEND_TUESDAY     0x0004
#define IRX_WEEKEND_WEDNESDAY   0x0008
#define IRX_WEEKEND_THURSDAY    0x0010
#define IRX_WEEKEND_FRIDAY      0x0020
#define IRX_WEEKEND_SATURDAY    0x0040
#define IRX_WEEKEND_NO_WEEKENDS 0x0000
#define IRX_WEEKEND_STANDARD    (IRX_WEEKEND_SATURDAY | IRX_WEEKEND_SUNDAY)
#define IRX_WEEKEND_MIN_VALUE   0x0000
#define IRX_WEEKEND_MAX_VALUE   0x007f
/* note that if weekends is the MAX_VALUE then every day is a weekend */
/* nice but probably not very useful! */

#define IRX_IS_WEEKEND(date, weekends) ((1 << irxDayOfWeek(date)) & (weekends))
#define IRX_IS_WEEKDAY(date, weekends) (! IRX_IS_WEEKEND(date, weekends))

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
char* irxWeekendsToString (long weekends);

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
int irxWeekendsFromString (char *str, long *weekends);

/*f
***************************************************************************
** Creates a new holiday structure.
**
** The date list can be NULL, in which case the resulting date list in the
** holiday structure will be a date list with no dates, e.g.
**      hl->dateList->fNumItems = 0;
**
** Note that this routine strips dates from the date list where the date
** is a weekend.
***************************************************************************
*/
IrxTCalendar* irxCalendarNew
(const IrxTDateList *dateList,              /* (I) Date list to use */
 IrxTBool      satIsAlwaysHoliday,    /* (I) TRUE if saturday always holiday*/
 IrxTBool      sunIsAlwaysHoliday);   /* (I) TRUE if sunday always holiday */

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
(IrxTCalendar *cal);    /* (I/O) Calendar to be adjusted */

/*f
***************************************************************************
** Indicates whether a date is a business day according to a calendar.
***************************************************************************
*/
int irxIsBusinessDay
(IrxTDate      date,              /* (I) Input Date */
 const IrxTCalendar *cal,               /* (I) Calendar structure */
 IrxTBool     *isBusinessDay);    /* (O) TRUE if a business day */

/*f
***************************************************************************
** Checks to see if a date is a holiday. Does not take week-ends into
** account.
***************************************************************************
*/
int irxIsHoliday
(IrxTDate      date,       /* (I) Arbitrary date          */
 const IrxTCalendar *cal,        /* (I) Calendar structure  */
 IrxTBool     *isHoliday); /* (O) 0 = not a holiday, 1 = is a holiday */

/*f
***************************************************************************
** Converts date to a valid business date using the calendar and bad day
** convention.
**
** Returns SUCCESS if convention OK, FAILURE otherwise.
***************************************************************************
*/
int irxBusinessDay
(IrxTDate        date,       /* (I) Arbitrary date. */
 long            badDayConv, /* (I) Bad day convention for adjusting
                              non-business days. Use one of the following:
                              IRX_BAD_DAY_FOLLOW
                              IRX_BAD_DAY_PREVIOUS
                              IRX_BAD_DAY_MODIFIED
                              IRX_BAD_DAY_NONE */
 const IrxTCalendar  *cal,         /* (I) Calendar. */
 IrxTDate      *outDate);    /* (O) Valid business day. */

/*f
***************************************************************************
** Converts date to a valid business date using all the given calendars and
** bad day convention.
**
** Returns SUCCESS if convention OK, FAILURE otherwise.
***************************************************************************
*/
int irxBusinessDayMultiCalendar
(IrxTDate       date,             /* (I) Arbitrary date. */
 long         badDayConv,       /* (I) Bad day convention */
 int          numCalendars,     /* (I) Number of calendars */
 IrxTCalendar const **calendars,        /* (I) [numCalendars] Calendars. */
 IrxTDate      *outDate);         /* (O) Valid business day. */

/*f
***************************************************************************
** Calculates the number of business days between two dates (FROM & TO).
**
** Algorithm:
**   1. if FROM = TO, the result is 0
**   2. if FROM < TO, the result is
**      the number of business days in the CLOSED interval
**      of [FROM+1,TO]
**   3. if FROM > TO, the result is
**      negated number of business days in the CLOSED interval
**      of [TO,FROM-1]
***************************************************************************
*/
int irxBusinessDaysDiff
(IrxTDate      from,     /* (I) Earlier date            */
 IrxTDate      to,       /* (I) Later date              */
 const IrxTCalendar *calendar, /* (I) Holiday calendar        */
 long       *result);  /* (O) Number of business days */

/*f
***************************************************************************
** Calculates a business date being at offset business days
** from the original date
***************************************************************************
*/
int  irxDateAddBusinessDays
(IrxTDate      fromDate,       /* (I) Start date                 */
 long        offset,         /* (I) Number of business days    */
 const IrxTCalendar *calendar,       /* (I) Holiday calendar           */
 IrxTDate     *result);        /* (O) resulting business date    */

/*f
***************************************************************************
** Computes the last business day of the month.
***************************************************************************
*/
int irxDateToBusinessEOM
(IrxTDate      inDate,    /* (I) Input date */
 const IrxTCalendar *calendar,  /* (I) Holiday calendar */
 IrxTDate     *result     /* (O) Last business day of the month */
);

/*f
***************************************************************************
** Returns the calendar from a name using the calendar cache.
**
** You can use name = "NONE" or name = "NO_WEEKENDS" and the resulting
** calendar list will be defined appropriately.
**
** Returns NULL on failure, a valid IrxTCalendar pointer on success. Does
** not report errors.
***************************************************************************
*/
IrxTCalendar* irxCalendarCacheSearch
(const char  *name          /* (I) Name of calendar */
);

/*f
***************************************************************************
** Adds a calendar to the calendar cache. If the entry already exists
** in the cache, then the old version will be deleted.
**
** On SUCCESS, the calendar is owned by the calendar cache.
***************************************************************************
*/
int irxCalendarCacheAdd
(char       *name,    /* (I) Name to associate calendars with */
 IrxTCalendar *calendar /* (I) Adds shallow copy */
);

/*f
***************************************************************************
** Deletes all cache entries in hash table.  Useful for being fastidious
** about memory leaks.
***************************************************************************
*/
void irxCalendarCacheClear (void);

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
);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif







