/*
***************************************************************************
** SOURCE FILE: irxutils.c
**
** Defines data structures and enums used in the irxutils library.
***************************************************************************
*/

#include "irx/irxutils.h"
#include "irx/macros.h"
#include "irx/calendar.h"
#include "irx/dateutils.h"
#include "irx/mdmin.h"

#include <ctype.h>

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
char* irxBadDayConvToString (IrxTBadDayConv value)
{
    static char routine[] = "irxBadDayConvToString";

    switch (value)
    {
    case IRX_BAD_DAY_NONE:              return "None";
    case IRX_BAD_DAY_FOLLOW:            return "Following";
    case IRX_BAD_DAY_PREVIOUS:          return "Previous";
    case IRX_BAD_DAY_MODIFIED:          return "Modified Following";
    }

    irxError("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts BadDayConv from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxBadDayConvFromString (const char* str, IrxTBadDayConv *val)
{
    static char routine[] = "irxBadDayConvFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        irxError ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'F':
        if (strcmp (buf, "F") == 0) *val = IRX_BAD_DAY_FOLLOW;
        else goto RETURN;
        break;
    case 'M':
        if (strcmp (buf, "M") == 0) *val = IRX_BAD_DAY_MODIFIED;
        else goto RETURN;
        break;
    case 'N':
        if (strcmp (buf, "N") == 0) *val = IRX_BAD_DAY_NONE;
        else goto RETURN;
        break;
    case 'P':
        if (strcmp (buf, "P") == 0) *val = IRX_BAD_DAY_PREVIOUS;
        else goto RETURN;
        break;
    default:
        goto RETURN;
    }
    status = SUCCESS;

  RETURN:

    if (status != SUCCESS && str != NULL && val != NULL)
        irxError("%s: Unknown value %s\n", routine, str);
    return status;
}

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
char* irxDayCountConvToString (IrxTDayCountConv value)
{
    static char routine[] = "irxDayCountConvToString";

    switch (value)
    {
    case IRX_ACT_ACT:                   return "ACT/ACT";
    case IRX_ACT_365F:                  return "ACT/365F";
    case IRX_ACT_360:                   return "ACT/360";
    case IRX_B30_360:                   return "30/360";
    case IRX_B30E_360:                  return "30E/360";
    }

    irxError("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}


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
int irxDayCountConvFromString (const char* str, IrxTDayCountConv *val)
{
    static char routine[] = "irxDayCountConvFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[9];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        irxError ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case '3':
        if (strcmp (buf, "30/360") == 0) *val = IRX_B30_360;
        else if (strcmp (buf, "30E/360") == 0) *val = IRX_B30E_360;
        else if (strcmp (buf, "3") == 0) *val = IRX_B30_360;   /*Fix3 style */
        else goto RETURN;
        break;
    case 'A':
        if (strcmp (buf, "ACT/360") == 0) *val = IRX_ACT_360;
        else if (strcmp (buf, "ACT/365F") == 0) *val = IRX_ACT_365F;
        else if (strcmp (buf, "ACT/ACT") == 0) *val = IRX_ACT_ACT;
        else if (strcmp (buf, "A") == 0) *val = IRX_ACT_ACT;   /*Fix3 style */
        else goto RETURN;
        break;
    case '5':   /*Fix3 style */
        if (strcmp (buf, "5") == 0) *val = IRX_ACT_365F;
        else goto RETURN;
    case '0':   /*Fix3 style */
        if (strcmp (buf, "0") == 0) *val = IRX_ACT_360;
        else goto RETURN;
    default:
        goto RETURN;
    }
    status = SUCCESS;

  RETURN:

    if (status != SUCCESS && str != NULL && val != NULL)
        irxError("%s: Unknown value %s\n", routine, str);
    return status;
}

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
char* irxStubLocationToString (IrxTStubLocation value)
{
    static char routine[] = "irxStubLocationToString";

    switch (value)
    {
    case IRX_SHORT_FRONT_STUB:          return "SF";
    case IRX_SHORT_BACK_STUB:           return "SB";
    case IRX_LONG_FRONT_STUB:           return "LF";
    case IRX_LONG_BACK_STUB:            return "LB";
    }

    irxError("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts StubLocation from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxStubLocationFromString (const char* str, IrxTStubLocation *val)
{
    static char routine[] = "irxStubLocationFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[3];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        irxError ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'L':
        if (strcmp (buf, "LB") == 0) *val = IRX_LONG_BACK_STUB;
        else if (strcmp (buf, "LF") == 0) *val = IRX_LONG_FRONT_STUB;
        else goto RETURN;
        break;
    case 'S':
        if (strcmp (buf, "SB") == 0) *val = IRX_SHORT_BACK_STUB;
        else if (strcmp (buf, "SF") == 0) *val = IRX_SHORT_FRONT_STUB;
        else goto RETURN;
        break;
    default:
        goto RETURN;
    }
    status = SUCCESS;

  RETURN:

    if (status != SUCCESS && str != NULL && val != NULL)
        irxError("%s: Unknown value %s\n", routine, str);
    return status;
}

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
char* irxPeriodTypeToString (char value)
{
    static char routine[] = "irxPeriodTypeToString";

    switch (value)
    {
    case IRX_PRD_TYPE_ANNUAL:           return "Annual";
    case IRX_PRD_TYPE_SEMI_ANNUAL:      return "Semi-Annual";
    case IRX_PRD_TYPE_YEAR:             return "Year";
    case IRX_PRD_TYPE_MONTH:            return "Month";
    case IRX_PRD_TYPE_QUARTER:          return "Quarter";
    case IRX_PRD_TYPE_WEEK:             return "Week";
    case IRX_PRD_TYPE_DAY:              return "Day";
    case IRX_PRD_TYPE_LUNAR_MONTH:      return "Lunar Month";
    }

    irxError("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts PeriodType from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxPeriodTypeFromString (const char* str, char *val)
{
    static char routine[] = "irxPeriodTypeFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[2];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        irxError ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'A':
        if (strcmp (buf, "A") == 0) *val = IRX_PRD_TYPE_ANNUAL;
        else goto RETURN;
        break;
    case 'D':
        if (strcmp (buf, "D") == 0) *val = IRX_PRD_TYPE_DAY;
        else goto RETURN;
        break;
    case 'L':
        if (strcmp (buf, "L") == 0) *val = IRX_PRD_TYPE_LUNAR_MONTH;
        else goto RETURN;
        break;
    case 'M':
        if (strcmp (buf, "M") == 0) *val = IRX_PRD_TYPE_MONTH;
        else goto RETURN;
        break;
    case 'Q':
        if (strcmp (buf, "Q") == 0) *val = IRX_PRD_TYPE_QUARTER;
        else goto RETURN;
        break;
    case 'S':
        if (strcmp (buf, "S") == 0) *val = IRX_PRD_TYPE_SEMI_ANNUAL;
        else goto RETURN;
        break;
    case 'W':
        if (strcmp (buf, "W") == 0) *val = IRX_PRD_TYPE_WEEK;
        else goto RETURN;
        break;
    case 'Y':
        if (strcmp (buf, "Y") == 0) *val = IRX_PRD_TYPE_YEAR;
        else goto RETURN;
        break;
    default:
        goto RETURN;
    }
    status = SUCCESS;

  RETURN:

    if (status != SUCCESS && str != NULL && val != NULL)
        irxError("%s: Unknown value %s\n", routine, str);
    return status;
}

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
char* irxWeekDayToString (long value)
{
    static char routine[] = "irxWeekDayToString";

    switch (value)
    {
    case IRX_SUNDAY:                    return "Sunday";
    case IRX_MONDAY:                    return "Monday";
    case IRX_TUESDAY:                   return "Tuesday";
    case IRX_WEDNESDAY:                 return "Wednesday";
    case IRX_THURSDAY:                  return "Thursday";
    case IRX_FRIDAY:                    return "Friday";
    case IRX_SATURDAY:                  return "Saturday";
    }

    irxError("%s: Unknown value %ld\n", routine, (long)value);
    return ("***ERROR***");
}

/**
***************************************************************************
** Converts WeekDay from string format.
**
** Returns SUCCESS or FAILURE depending on whether the string is valid
** or not.
***************************************************************************
*/
int irxWeekDayFromString (const char* str, long *val)
{
    static char routine[] = "irxWeekDayFromString";
    int         status    = FAILURE;

    size_t i;
    size_t n;
    char   buf[4];

    REQUIRE(val != NULL);
    if (str == NULL || *str == '\0')
    {
        irxError ("%s: Undefined string input\n", routine);
        return FAILURE;
    }
    
    n = strlen(str);
    if (n >= sizeof(buf)) n = sizeof(buf)-1;
    for (i = 0; i < n; ++i) buf[i] = (char) toupper(str[i]);
    buf[n] = '\0';
    
    switch (buf[0])
    {
    case 'F':
        if (strcmp (buf, "FRI") == 0) *val = IRX_FRIDAY;
        else goto RETURN;
        break;
    case 'M':
        if (strcmp (buf, "MON") == 0) *val = IRX_MONDAY;
        else goto RETURN;
        break;
    case 'S':
        if (strcmp (buf, "SAT") == 0) *val = IRX_SATURDAY;
        else if (strcmp (buf, "SUN") == 0) *val = IRX_SUNDAY;
        else goto RETURN;
        break;
    case 'T':
        if (strcmp (buf, "THU") == 0) *val = IRX_THURSDAY;
        else if (strcmp (buf, "TUE") == 0) *val = IRX_TUESDAY;
        else goto RETURN;
        break;
    case 'W':
        if (strcmp (buf, "WED") == 0) *val = IRX_WEDNESDAY;
        else goto RETURN;
        break;
    default:
        goto RETURN;
    }
    status = SUCCESS;

  RETURN:

    if (status != SUCCESS && str != NULL && val != NULL)
        irxError("%s: Unknown value %s\n", routine, str);
    return status;
}

/**
***************************************************************************
** Constructor for IrxTDateList
***************************************************************************
*/
IrxTDateList* irxDateListMake(
int             fNumItems,           /* (I) */
IrxTDate const* fArray               /* (I) [fNumItems] */
)
{
    static char routine[] = "irxDateListMake";
    int status = FAILURE;

    IrxTDateList* p = NULL;


    p = NEW(IrxTDateList);
    if (p==NULL) goto RETURN;

    p->fNumItems       = fNumItems;
    if (p->fNumItems > 0 && fArray != NULL)
    {
        p->fArray = NEW_ARRAY(IrxTDate, p->fNumItems);
        if (p->fArray == NULL) goto RETURN;
        COPY_ARRAY (p->fArray, fArray, IrxTDate, p->fNumItems);
    }
    else
    {
        p->fArray = NULL;
    }


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxDateListFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTDateList
***************************************************************************
*/
IrxTDateList* irxDateListMakeEmpty(
int             fNumItems            /* (I) */
)
{
    static char routine[] = "irxDateListMakeEmpty";
    int status = FAILURE;

    IrxTDateList* p = NULL;


    p = NEW(IrxTDateList);
    if (p==NULL) goto RETURN;

    p->fNumItems       = fNumItems;

    if (p->fNumItems > 0)
    {
        p->fArray = NEW_ARRAY(IrxTDate, p->fNumItems);
        if (p->fArray == NULL) goto RETURN;
    }
    else
    {
        p->fArray = NULL;
    }


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxDateListFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTDateList
***************************************************************************
*/
IrxTDateList* irxDateListCopy(IrxTDateList const* src)
{
    IrxTDateList* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxDateListMake(src->fNumItems,
                          src->fArray);

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTDateList
***************************************************************************
*/
void irxDateListFree(IrxTDateList *p)
{
    if (p != NULL)
    {
        FREE(p->fArray);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for IrxTCalendar
***************************************************************************
*/
IrxTCalendar* irxCalendarMake(
IrxTDateList const* dateList,            /* (I) */
long            weekends             /* (I) */
)
{
    static char routine[] = "irxCalendarMake";
    int status = FAILURE;

    IrxTCalendar* p = NULL;


    p = NEW(IrxTCalendar);
    if (p==NULL) goto RETURN;

    if (dateList != NULL)
    {
        p->dateList = irxDateListCopy(dateList);
        if (p->dateList == NULL) goto RETURN;
    }
    else
    {
        p->dateList = NULL;
    }

    p->weekends        = weekends;

    if (irxCalendarValidate(p) != SUCCESS) goto RETURN; /* failure */

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxCalendarFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTCalendar
***************************************************************************
*/
IrxTCalendar* irxCalendarMakeEmpty(void)
{
    static char routine[] = "irxCalendarMakeEmpty";
    int status = FAILURE;

    IrxTCalendar* p = NULL;


    p = NEW(IrxTCalendar);
    if (p==NULL) goto RETURN;


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxCalendarFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTCalendar
***************************************************************************
*/
IrxTCalendar* irxCalendarCopy(IrxTCalendar const* src)
{
    IrxTCalendar* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxCalendarMake(src->dateList,
                          src->weekends);

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTCalendar
***************************************************************************
*/
void irxCalendarFree(IrxTCalendar *p)
{
    if (p != NULL)
    {
        irxDateListFree(p->dateList);
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for IrxTDateInterval
***************************************************************************
*/
IrxTDateInterval* irxDateIntervalMake(
int             prd,                 /* (I) */
char            prd_typ,             /* (I) */
IrxTBool        eom                  /* (I) */
)
{
    static char routine[] = "irxDateIntervalMake";
    int status = FAILURE;

    IrxTDateInterval* p = NULL;


    p = NEW(IrxTDateInterval);
    if (p==NULL) goto RETURN;

    p->prd             = prd;
    p->prd_typ         = prd_typ;
    p->eom             = eom;

    if (irxDateIntervalValidate(p) != SUCCESS) goto RETURN; /* failure */

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxDateIntervalFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTDateInterval
***************************************************************************
*/
IrxTDateInterval* irxDateIntervalMakeEmpty(void)
{
    static char routine[] = "irxDateIntervalMakeEmpty";
    int status = FAILURE;

    IrxTDateInterval* p = NULL;


    p = NEW(IrxTDateInterval);
    if (p==NULL) goto RETURN;


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxDateIntervalFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTDateInterval
***************************************************************************
*/
IrxTDateInterval* irxDateIntervalCopy(IrxTDateInterval const* src)
{
    IrxTDateInterval* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxDateIntervalMake(src->prd,
                              src->prd_typ,
                              src->eom);

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTDateInterval
***************************************************************************
*/
void irxDateIntervalFree(IrxTDateInterval *p)
{
    if (p != NULL)
    {
        FREE(p);
    }
}

/**
***************************************************************************
** Constructor for IrxTMultiDimMinState
***************************************************************************
*/
IrxTMultiDimMinState* irxMultiDimMinStateMake(
int             n,                   /* (I) */
double const*   x,                   /* (I) [n] */
IrxTMatrix2D const* direction,           /* (I) */
int             iter,                /* (I) */
double          value,               /* (I) */
double          vdiff                /* (I) */
)
{
    static char routine[] = "irxMultiDimMinStateMake";
    int status = FAILURE;

    IrxTMultiDimMinState* p = NULL;

    REQUIRE(n > 0);
    REQUIRE(x != NULL);

    p = NEW(IrxTMultiDimMinState);
    if (p==NULL) goto RETURN;

    p->n               = n;
    p->x = NEW_ARRAY(double, p->n);
    if (p->x == NULL) goto RETURN;
    COPY_ARRAY (p->x, x, double, p->n);

    if (direction != NULL)
    {
        p->direction = irxMatrixCopy(direction);
        if (p->direction == NULL) goto RETURN;
    }
    else
    {
        p->direction = NULL;
    }

    p->iter            = iter;
    p->value           = value;
    p->vdiff           = vdiff;

    if (irxMultiDimMinStateValidate(p) != SUCCESS) goto RETURN; /* failure */

    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxMultiDimMinStateFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for IrxTMultiDimMinState
***************************************************************************
*/
IrxTMultiDimMinState* irxMultiDimMinStateMakeEmpty(
int             n                    /* (I) */
)
{
    static char routine[] = "irxMultiDimMinStateMakeEmpty";
    int status = FAILURE;

    IrxTMultiDimMinState* p = NULL;

    REQUIRE(n > 0);

    p = NEW(IrxTMultiDimMinState);
    if (p==NULL) goto RETURN;

    p->n               = n;

    p->x = NEW_ARRAY(double, p->n);
    if (p->x == NULL) goto RETURN;


    status = SUCCESS;

RETURN:
    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        irxMultiDimMinStateFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for IrxTMultiDimMinState
***************************************************************************
*/
IrxTMultiDimMinState* irxMultiDimMinStateCopy(IrxTMultiDimMinState const* src)
{
    IrxTMultiDimMinState* dst = NULL;
    if (src==NULL) return NULL;

    dst = irxMultiDimMinStateMake(src->n,
                                  src->x,
                                  src->direction,
                                  src->iter,
                                  src->value,
                                  src->vdiff);

    return dst;
}

/**
***************************************************************************
** Destructor for IrxTMultiDimMinState
***************************************************************************
*/
void irxMultiDimMinStateFree(IrxTMultiDimMinState *p)
{
    if (p != NULL)
    {
        FREE(p->x);
        irxMatrixFree(p->direction);
        FREE(p);
    }
}

