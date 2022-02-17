/*
***************************************************************************
** HEADER FILE: date.h
** 
** Defines the date types used in credit & rates exotics.
** If these seem remarkably familiar then that is not surprising.
***************************************************************************
*/

#ifndef IRX_DATE_H
#define IRX_DATE_H

#include "cgeneral.h"

/* Some useful macros */
#define IRX_IS_LEAP(year) ( (((year)%4 == 0) && ((year)%100 != 0)) || \
                            ((year)%400 == 0) )

#ifdef __cplusplus
extern "C"
{
#endif

/**
***************************************************************************
** Defines the date type used in credit & rate exotics.
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
**
** This date is identical to the dates used in ALIB and QLIB.
***************************************************************************
*/
typedef long IrxTDate;

#define IRX_INVALID_DATE(date) ((date) < 0)

#define IRX_DAY_OF_WEEK(date) ( ((date)+1) % 7)

/*t
***************************************************************************
** Breakdown of a date into constituent parts year,month,day.
***************************************************************************
*/
typedef struct _IrxTYearMonthDay
{
    long year;      /* In range [1600-] */
    long month;     /* In range [1,12] */
    long day;       /* In range [1-31] */
} IrxTYearMonthDay;


/**
***************************************************************************
** Returns today's date
***************************************************************************
*/
IrxTDate irxToday(void);
                           
/**
***************************************************************************
** Converts year, month, day to a date and returns the date. This is useful
** for initializing dates within test routines. Returns -1 for bad dates.
***************************************************************************
*/
IrxTDate irxDate
(long year,  /* (I) Year */
 long month, /* (I) Month */
 long day    /* (I) Day */
);

/**
***************************************************************************
** Converts IrxTDate to Year, Month, Day.
***************************************************************************
*/
int irxDateToYMD
(IrxTDate          date,           /* (I) IrxTDate format */
 IrxTYearMonthDay *ymd);           /* (O) Date in yyyy-mm-dd format */
                 
/**
***************************************************************************
** Converts Year, Month, Day to IrxTDate.
***************************************************************************
*/
int irxYMDToDate
(IrxTYearMonthDay *ymd,            /* (I) Date in yyyy-mm-dd format */
 IrxTDate         *date);          /* (O) IrxTDate format */

/**
***************************************************************************
** Normalizes a month/day/year. If month is out of range, it is brought
** into range, and the years are incremented or decremented as appropriate.
** If day belongs to a month/year combination which does not exist (such as
** April 31), the day is reduced so that it becomes valid (to April 30).
***************************************************************************
*/
int irxNormalizeYMD
(IrxTYearMonthDay *ymd);          /* (I/O) */

/**
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
IrxTDate irxDateAddMonths
(IrxTDate date,    /* (I) Initial date */
 long     months,  /* (I) Number of months to add */
 IrxTBool eom      /* (I) Perform end of month adjustment */
);

/**
***************************************************************************
** Goes to the end of the current month.
**
** Returns -1 on FAILURE.
***************************************************************************
*/
IrxTDate irxDateToEndOfMonth
(IrxTDate date    /* (I) */
);

/**
***************************************************************************
** Goes to the start of the current month.
**
** Returns -1 on FAILURE.
***************************************************************************
*/
IrxTDate irxDateToStartOfMonth
(IrxTDate date   /* (I) */
);

/**
***************************************************************************
** Goes to the end of the current year.
**
** Returns -1 on FAILURE.
***************************************************************************
*/
IrxTDate irxDateToEndOfYear
(IrxTDate date    /* (I) */
);

/**
***************************************************************************
** Goes to the start of the current year.
**
** Returns -1 on FAILURE.
***************************************************************************
*/
IrxTDate irxDateToStartOfYear
(IrxTDate date   /* (I) */
);

/**
***************************************************************************
** Returns the weekday (see IRX_SUNDAY etc.)
***************************************************************************
*/
long irxDayOfWeek
(IrxTDate date   /* (I) */
);

/**
***************************************************************************
** Returns the number of days in the given month.
***************************************************************************
*/
long irxDaysInMonth
(long year,   /* (I) */
 long month   /* (I) Range 1-12 */
);

/**
***************************************************************************
** Computes the Nth weekday of a month, e.g. the 3rd Wednesday of September
**
** Returns -1 for a bad request
***************************************************************************
*/
IrxTDate irxDateNthWeekDay
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
const char* irxFormatDayOfWeek(long dayOfWeek); /* (I) */

/**
***************************************************************************
** Formats a date according to a format string. If the format string
** is invalid, then a NULL string is returned. If the caller provides an
** output buffer, then this will be populated and returned. Otherwise the
** returned string needs to be deleted by the caller.
**
** If the format is not provided (NULL or empty string), then uses the
** format YYYYMMDD.
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
char* irxDateFormat
(IrxTDate    date,    /* (I) Date to be formatted */
 const char *format,  /* (I) Method of formatting the date */
 char       *output   /* (O) Outputs to this buffer (if provided) */
);

/**
***************************************************************************
** Converts a string in various formats into a date.
**
** The formats we accept at present are
** YYYYMMDD,  MM/DD/YYYY, MM-DD-YYYY, Mth-DD-YYYY, DD-Mth-YYYY
** 20050201,  02/01/2005, 02-01-2005, Feb-01-2005, 01-Feb-2005, FEB-01-2005
***************************************************************************
*/
int irxStringToDate
(const char* str,
 IrxTDate*   result);

#ifdef __cplusplus
}
#endif

#endif

