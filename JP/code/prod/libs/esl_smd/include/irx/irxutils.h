/*
***************************************************************************
** HEADER FILE: irxutils.h
**
** Defines data structures and enums used in the irxutils library.
***************************************************************************
*/

#ifndef _IRX_IRXUTILS_H
#define _IRX_IRXUTILS_H

#include "date.h"
#include "cgeneral.h"
#include "matrix2d.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** Bad day conventions used for converting bad days (weekends and holidays)
    into good days (working days where the market is open for business).
   
    Defined as an enumerated data type. Only the first letter is significant
    in the text format.
   
    Note that the values used are the same as in the ALIB. */
typedef enum
{
/** No adjustment is made for bad days. */
    IRX_BAD_DAY_NONE = 78,                   /* None */
/** Bad days are adjusted to the first good day after the bad day. */
    IRX_BAD_DAY_FOLLOW = 70,                 /* Following */
/** Bad days are adjusted to the first good day before the bad day. */
    IRX_BAD_DAY_PREVIOUS = 80,               /* Previous */
/** Bad days are adjusted to the first good day after the bad day, unless
    this takes you into a new month. In the latter case, the adjustment is
    then made to the first good day before the bad day. */
    IRX_BAD_DAY_MODIFIED = 77                /* Modified Following */
} IrxTBadDayConv;

/** Day count conventions are used to calculate the time elapsed between
    two dates. These are often used for the calculation of accrued interest
    payments and coupon payments.
   
    Designed to be cast as long into equivalent ALIB type. */
typedef enum
{
/** The ISDA ACT/ACT convention. Not to be confused
    with the Bond ACT/ACT convention which is different! */
    IRX_ACT_ACT = 1,                         /* ACT/ACT */
/** Day count fraction calculated as
    \begin{equation*}
    \frac{date2-date1}{365}
    \end{equation*} */
    IRX_ACT_365F = 2,                        /* ACT/365F */
/** Day count fraction calculated as
    \begin{equation*}
    \frac{date2-date1}{360}
    \end{equation*}
    
    This is a very common convention for CDS coupon calculation. */
    IRX_ACT_360 = 3,                         /* ACT/360 */
/** Day count fraction calculated as
    \begin{equation*}
    \frac{DaysDiff(date1,date2)}{360}
    \end{equation*}
    where the days difference is calculated as
    $360.(year2-year1) + 30.(month2-month1) + (day2 - day1)$
    where we convert day1 from 31 to 30 and we convert day2 from 31 to 30
    but only if day1 is 30 (or 31).
    
    This is a very common convention for swap coupon calculation. */
    IRX_B30_360 = 4,                         /* 30/360 */
/** Day count fraction calculated as
    \begin{equation*}
    \frac{DaysDiff(date1,date2)}{360}
    \end{equation*}
    where the days difference is calculated as
    $360.(year2-year1) + 30.(month2-month1) + (day2 - day1)$
    where we convert day1 from 31 to 30 and we convert day2 from 31 to 30
    in all cases. */
    IRX_B30E_360 = 5                         /* 30E/360 */
} IrxTDayCountConv;

typedef enum
{
    IRX_SHORT_FRONT_STUB,                    /* SF */
    IRX_SHORT_BACK_STUB,                     /* SB */
    IRX_LONG_FRONT_STUB,                     /* LF */
    IRX_LONG_BACK_STUB                       /* LB */
} IrxTStubLocation;

#define IRX_PRD_TYPE_ANNUAL  (char)'A'              /* Annual */
#define IRX_PRD_TYPE_SEMI_ANNUAL (char)'S'              /* Semi-Annual */
#define IRX_PRD_TYPE_YEAR    (char)'Y'              /* Year */
#define IRX_PRD_TYPE_MONTH   (char)'M'              /* Month */
#define IRX_PRD_TYPE_QUARTER (char)'Q'              /* Quarter */
#define IRX_PRD_TYPE_WEEK    (char)'W'              /* Week */
#define IRX_PRD_TYPE_DAY     (char)'D'              /* Day */
#define IRX_PRD_TYPE_LUNAR_MONTH (char)'L'              /* Lunar Month */

#define IRX_SUNDAY           (long)0                /* Sunday */
#define IRX_MONDAY           (long)1                /* Monday */
#define IRX_TUESDAY          (long)2                /* Tuesday */
#define IRX_WEDNESDAY        (long)3                /* Wednesday */
#define IRX_THURSDAY         (long)4                /* Thursday */
#define IRX_FRIDAY           (long)5                /* Friday */
#define IRX_SATURDAY         (long)6                /* Saturday */

/** Date list structure - a simple list of dates. */
typedef struct _IrxTDateList
{
    int             fNumItems;
    /** Array of size fNumItems. */
    IrxTDate*       fArray;
} IrxTDateList;

/** Define the calendar structure. The constructor is clever - it removes
    weekends from the list of dates.
   
    Note that the structure is identical to the one used within the ALIB.
   
    Note that a NULL calendar will typically be equivalent to no holidays
    and standard weekends. */
typedef struct _IrxTCalendar
{
    IrxTDateList*   dateList;
    long            weekends;
} IrxTCalendar;

/** The date interval is slightly different from the ALIB equivalent in that
    we have added the end of month adjustment flag to allow for monthly
    type intervals to specify end of month adjustment as well. */
typedef struct _IrxTDateInterval
{
    int             prd;
    char            prd_typ;
    IrxTBool        eom;
} IrxTDateInterval;

/** Defines the current state of a multi-dimensional minimization. Designed so
    that when this structure is returned from the multi-dimensional solver
    that it can be used as the initial state for a subsequent call. */
typedef struct _IrxTMultiDimMinState
{
    /** Number of dimensions */
    int             n;
    /** Array of size n.
        Initial starting point (on initialization).
        Solution (on exit). */
    double*         x;
    /** Initial directions for search (on initialization). Final directions for
        search (on exit). */
    IrxTMatrix2D*   direction;
    /** Number of iterations. */
    int             iter;
    /** Minimal value. */
    double          value;
    /** Last step-size in the algorithm. */
    double          vdiff;
} IrxTMultiDimMinState;

/**
***************************************************************************
** Converts BadDayConv to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxBadDayConvToString (IrxTBadDayConv value);

/**
***************************************************************************
** Converts BadDayConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxBadDayConvFromString (const char* str, IrxTBadDayConv *val);

/**
***************************************************************************
** Converts DayCountConv to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxDayCountConvToString (IrxTDayCountConv value);

/**
***************************************************************************
** Converts DayCountConv from string format.
**
** Support following formats (case insensitive):
** ALib styles:  30/360, 30E/360, ACT/360, ACT/365F, ACT/ACT
** London style: '5', '3', '0', 'A'
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxDayCountConvFromString (const char* str, IrxTDayCountConv *val);

/**
***************************************************************************
** Converts StubLocation to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxStubLocationToString (IrxTStubLocation value);

/**
***************************************************************************
** Converts StubLocation from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxStubLocationFromString (const char* str, IrxTStubLocation *val);

/**
***************************************************************************
** Converts PeriodType to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxPeriodTypeToString (char value);

/**
***************************************************************************
** Converts PeriodType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxPeriodTypeFromString (const char* str, char *val);

/**
***************************************************************************
** Converts WeekDay to string format.
**
** Returns a static string which should not be freed (hence this can be
** used within printf(...) style statements.
**
** If the input is not valid, then returns "***ERROR***"
***************************************************************************
*/
char* irxWeekDayToString (long value);

/**
***************************************************************************
** Converts WeekDay from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxWeekDayFromString (const char* str, long *val);

/**
***************************************************************************
** Constructor for IrxTDateList
***************************************************************************
*/
IrxTDateList* irxDateListMake(
int             fNumItems,
/** Array of size fNumItems. */
IrxTDate const* fArray
);

/**
***************************************************************************
** Memory allocator for IrxTDateList
***************************************************************************
*/
IrxTDateList* irxDateListMakeEmpty(
int             fNumItems
);

/**
***************************************************************************
** Copy constructor for IrxTDateList
***************************************************************************
*/
IrxTDateList* irxDateListCopy(IrxTDateList const* src);

/**
***************************************************************************
** Destructor for IrxTDateList
***************************************************************************
*/
void irxDateListFree(IrxTDateList *p);

/**
***************************************************************************
** Constructor for IrxTCalendar
***************************************************************************
*/
IrxTCalendar* irxCalendarMake(
IrxTDateList const* dateList,
long            weekends
);

/**
***************************************************************************
** Memory allocator for IrxTCalendar
***************************************************************************
*/
IrxTCalendar* irxCalendarMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for IrxTCalendar
***************************************************************************
*/
IrxTCalendar* irxCalendarCopy(IrxTCalendar const* src);

/**
***************************************************************************
** Destructor for IrxTCalendar
***************************************************************************
*/
void irxCalendarFree(IrxTCalendar *p);

/**
***************************************************************************
** Constructor for IrxTDateInterval
***************************************************************************
*/
IrxTDateInterval* irxDateIntervalMake(
int             prd,
char            prd_typ,
IrxTBool        eom
);

/**
***************************************************************************
** Memory allocator for IrxTDateInterval
***************************************************************************
*/
IrxTDateInterval* irxDateIntervalMakeEmpty(void);

/**
***************************************************************************
** Copy constructor for IrxTDateInterval
***************************************************************************
*/
IrxTDateInterval* irxDateIntervalCopy(IrxTDateInterval const* src);

/**
***************************************************************************
** Destructor for IrxTDateInterval
***************************************************************************
*/
void irxDateIntervalFree(IrxTDateInterval *p);

/**
***************************************************************************
** Constructor for IrxTMultiDimMinState
***************************************************************************
*/
IrxTMultiDimMinState* irxMultiDimMinStateMake(
/** Number of dimensions */
int             n,
/** Array of size n.
    Initial starting point (on initialization).
    Solution (on exit). */
double const*   x,
/** Initial directions for search (on initialization). Final directions for
    search (on exit). */
IrxTMatrix2D const* direction,
/** Number of iterations. */
int             iter,
/** Minimal value. */
double          value,
/** Last step-size in the algorithm. */
double          vdiff
);

/**
***************************************************************************
** Memory allocator for IrxTMultiDimMinState
***************************************************************************
*/
IrxTMultiDimMinState* irxMultiDimMinStateMakeEmpty(
/** Number of dimensions */
int             n
);

/**
***************************************************************************
** Copy constructor for IrxTMultiDimMinState
***************************************************************************
*/
IrxTMultiDimMinState* irxMultiDimMinStateCopy(IrxTMultiDimMinState const* src);

/**
***************************************************************************
** Destructor for IrxTMultiDimMinState
***************************************************************************
*/
void irxMultiDimMinStateFree(IrxTMultiDimMinState *p);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _IRX_IRXUTILS_H */
