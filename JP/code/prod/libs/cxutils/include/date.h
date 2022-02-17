/*
***************************************************************************
** HEADER FILE: date.h
** 
** Defines the date types used in credit exotics.
** If these seem remarkably familiar then that is not surprising.
**
** $Header$
***************************************************************************
*/

#ifndef CX_DATE_H
#define CX_DATE_H

#include <alib/cdate.h>

/* Some useful macros */
#define CX_IS_LEAP(year) ( (((year)%4 == 0) && ((year)%100 != 0)) || \
                            ((year)%400 == 0) )

#ifdef __cplusplus
extern "C"
{
#endif

/*t
***************************************************************************
** Defines the date type used in credit exotics.
**
** This is an offset from year 1601. However note that a pedant might
** point out (with reason) that the date system changed from Julian to
** Gregorian since that date, so that early dates are meaningless.
**
** Any negative date is a bad date.
**
** A date of zero should often be used to indicate that the date is
** undefined. There is very little justification for using a date of zero
** as a real date (for the reasons given above regarding the switch from
** Julian to Gregorian calendar).
***************************************************************************
*/
#define CX_INVALID_DATE(date) ((date) < 0)

/* Some constants for days of the week */
#define CX_SUNDAY     0 
#define CX_MONDAY     1
#define CX_TUESDAY    2
#define CX_WEDNESDAY  3
#define CX_THURSDAY   4
#define CX_FRIDAY     5
#define CX_SATURDAY   6

#define CX_DAY_OF_WEEK(date) ( ((date)+1) % 7)

/*t
***************************************************************************
** Breakdown of a date into constituent parts year,month,day.
***************************************************************************
*/
typedef struct _CxTYearMonthDay
{
    long year;      /* In range [1600-] */
    long month;     /* In range [1,12] */
    long day;       /* In range [1-31] */
} CxTYearMonthDay;


/*f
***************************************************************************
** Returns today's date
***************************************************************************
*/
TDate CxToday(void);
                           
/*f
***************************************************************************
** Converts year, month, day to a date and returns the date. This is useful
** for initializing dates within test routines. Returns -1 for bad dates.
***************************************************************************
*/
TDate CxDate
(long year,  /* (I) Year */
 long month, /* (I) Month */
 long day    /* (I) Day */
);

/*f
***************************************************************************
** Converts TDate to Year, Month, Day.
***************************************************************************
*/
int CxDateToYMD
(TDate          date,           /* (I) TDate format */
 CxTYearMonthDay *ymd);           /* (O) Date in yyyy-mm-dd format */
                 
/*f
***************************************************************************
** Converts Year, Month, Day to TDate.
***************************************************************************
*/
int CxYMDToDate
(CxTYearMonthDay *ymd,            /* (I) Date in yyyy-mm-dd format */
 TDate         *date);          /* (O) TDate format */

/*f
***************************************************************************
** Normalizes a month/day/year. If month is out of range, it is brought
** into range, and the years are incremented or decremented as appropriate.
** If day belongs to a month/year combination which does not exist (such as
** April 31), the day is reduced so that it becomes valid (to April 30).
***************************************************************************
*/
int CxNormalizeYMD
(CxTYearMonthDay *ymd);          /* (I/O) */

/*f
***************************************************************************
** Adds a number of months to a date.
**
** EOM means the following:
**
** If date is at the end of the month, then the resulting date is at the
** end of the month.
**
** If adding a number of months to date takes you to a date beyond the
** end of the month, then this is adjusted back to the end of the month.
***************************************************************************
*/
TDate CxDateAddMonths
(TDate date,    /* (I) Initial date */
 long   months,  /* (I) Number of months to add */
 TBoolean eom      /* (I) Perform end of month adjustment */
);

/*f
***************************************************************************
** Goes to the end of the current month.
**
** Returns -1 on FAILURE.
***************************************************************************
*/
TDate CxDateToEndOfMonth
(TDate date    /* (I) */
);

/*f
***************************************************************************
** Goes to the start of the current month.
**
** Returns -1 on FAILURE.
***************************************************************************
*/
TDate CxDateToStartOfMonth
(TDate date   /* (I) */
);

/*f
***************************************************************************
** Goes to the end of the current year.
**
** Returns -1 on FAILURE.
***************************************************************************
*/
TDate CxDateToEndOfYear
(TDate date    /* (I) */
);

/*f
***************************************************************************
** Goes to the start of the current year.
**
** Returns -1 on FAILURE.
***************************************************************************
*/
TDate CxDateToStartOfYear
(TDate date   /* (I) */
);

/*f
***************************************************************************
** Returns the weekday (see CX_SUNDAY etc.)
***************************************************************************
*/
long CxDayOfWeek
(TDate date   /* (I) */
);

/*f
***************************************************************************
** Returns the number of days in the given month.
***************************************************************************
*/
long CxDaysInMonth
(long year,   /* (I) */
 long month   /* (I) Range 1-12 */
);

/*f
***************************************************************************
** Computes the Nth weekday of a month, e.g. the 3rd Wednesday of September
**
** Returns -1 for a bad request
***************************************************************************
*/
TDate CxDateNthWeekDay
(long year,    /* (I) Year */
 long month,   /* (I) Month */
 long weekday, /* (I) Weekday */
 long whichOne /* (I) n as in nth */
);

/*
***************************************************************************
** Formats the day of the week. You get the long form of the name.
***************************************************************************
*/
char* CxFormatDayOfWeek(long dayOfWeek); /* (I) */

/*f
***************************************************************************
** Formats a date according to a format string. If the format string
** is invalid, then a NULL string is returned. If the caller provides an
** output buffer, then this will be populated and returned. Otherwise the
** returned string needs to be deleted by the caller.
**
** Format as follows:
**
** D indicates day
** Q indicates quarter in the year
** M indicates month in the year
** Y indicates year
** U indicates upper case
** L indicates lower case
** X indicates mixed case - mixed case is the default
**
** All case indicators apply from there-on.
**
** Any other letters are invalid.
** Any numbers are invalid.
** Any other characters are appended verbatim.
**
** Number of letters indicates further format information
**
** D: 1 - date as number - no padding
**    2 - date as number - padded with zero
**    3 - date as weekday - abbreviated (3 letters)
**    4 - date as weekday - full
**
** Q: 1 - quarter as number
**    2 - quarter as number - padded with Q
**
** M: 1 - month as number - no padding
**    2 - month as number - padded with zero
**    3 - month as text - abbreviated (3 letters)
**    4 - month as text - full
**
** Y: 2 - year as number - last two digits only
**    4 - year as number - full
**
** Any other number of letters is an error.
*/
char* CxDateFormat
(TDate      date,    /* (I) Date to be formatted */
 const char *format,  /* (I) Method of formatting the date */
 char       *output   /* (O) Outputs to this buffer (if provided) */
);

/*f
***************************************************************************
** Converts a string in various formats into a date.
**
** The format is the same as used in CxDateFormat, but only a limited
** subset of the formats are accepted.
**
** In fact at present, the only format we accept is "YYYYMMDD"
***************************************************************************
*/
int CxStringToDate
(const char* str,
 const char* format,
 TDate*     result);

#ifdef __cplusplus
}
#endif

#endif




