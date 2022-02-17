/*
***************************************************************************
** HEADER FILE: calendar.h
**
** Defines the calendar object, the calendar cache (so we can access
** calendar objects by name) and basic business day functions.
**
** Calendar functions are all implemented by calling the equivalent ALIB
** function via the THolidayList interface. Only those functions currently
** used by CX have been wrapped in this manner.
***************************************************************************
*/

#ifndef CX_CALENDAR_H
#define CX_CALENDAR_H

#include "cxutils.h"    /* basic data types */

#include <alib/busday.h>

#ifdef __cplusplus
extern "C"
{
#endif

/**
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
char* CxWeekendsToString (
/** Sum of ALIB constants GTO_WEEKEND_MONDAY, GTO_WEEKEND_TUESDAY etc as
    defined in <alib/buscache.h>. Therefore a typical value might be 
    GTO_WEEKEND_SATURDAY + GTO_WEEKEND_SUNDAY which is also defined directly
    as GTO_WEEKEND_STANDARD. */
    long weekends);

/**
***************************************************************************
** Converts from string to weekends.
**
** The format is that the weekends are separated by commas. The names are
** the first three letters of the weekend name.
**
** Returns SUCCESS or FAILURE.
***************************************************************************
*/
int CxWeekendsFromString (
/** Comma separated string indicating the weekend, e.g. "SAT,SUN" would be
    the standard weekend. */
    char *str, 
/** Returns the value corresponding to the string as sum of ALIB constants
    GTO_WEEKEND_MONDAY, GTO_WEEKEND_TUESDAY etc. */
    long *weekends);

/**
***************************************************************************
** Validate a newly created calendar structure.
**
** This will do the following:
** 1. Sort the calendar date list and remove weekends.
** 2. Ensure that we don't have all days as weekends.
***************************************************************************
*/
int CxCalendarValidate
(
/** Calendar to be amended in place. */
    CxTCalendar *cal);

/**
***************************************************************************
** Converts date to a valid business date using the calendar and bad day
** convention.
**
** Returns SUCCESS if convention OK, FAILURE otherwise.
**
** Calls ALIB function GtoHolidayListBusinessDay.
***************************************************************************
*/
int CxBusinessDay(
/** Initial date. */
    TDate        date,
/** Bad day convention. */
    CxTBadDayConv badDayConv,
/** Calendar structure. Can be NULL in which case a calendar with just
    standard weekends will be used. */
    CxTCalendar  *cal,
/** Output. Date adjusted according to badDayConv and cal. This is not
    necessarily a business day (if badDayConv == CX_BAD_DAY_CONV_NONE then
    there is no adjustment. */
    TDate       *outDate);

/**
***************************************************************************
** Returns the calendar from a name using the calendar cache.
**
** You can use name = "NONE" or name = "NO_WEEKENDS" and the resulting
** calendar list will be defined appropriately.
**
** Returns NULL on failure, a valid CxTCalendar pointer on success. Does
** not report errors.
**
** Calls GtoIsHolidayLoaded and GtoHolidayListFromCache.
***************************************************************************
*/
CxTCalendar* CxCalendarCacheSearch(
/** Name of calendar in the cache */
    char  *name
);

/**
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
CxTCalendar* CxCalendarRead(
/** Name of file to read */
    char *fileName
);

#ifdef __cplusplus
}
#endif

#endif







