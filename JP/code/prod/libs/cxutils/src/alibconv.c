#include "alibconv.h"

#include <math.h>

#include <alib/buscache.h>
#include <alib/cdate.h>
#include <alib/cerror.h>
#include <alib/convert.h>
#include <alib/cversion.h>
#include <alib/ldate.h>
#include <alib/tcurve.h>
#if GTO_VERSION_NUM >= 1300
#include <alib/zccustom.h>
#endif
#include <alib/zcurve.h>

#include "cxmacros.h"
#include "zerocurve.h"

/*f
 * Converts a string representation ("2Y", "3M", "1S", etc) to a TDateInterval
 */
int CxStringToDateInterval
(char          *input,      /* (I) String w/ 1A, 3M, 4D, etc */
 TDateInterval *interval)   /* (O) Value read from file */
{
    return GtoStringToDateInterval (input, "DateInterval", interval);
}

/*f
 * Formats a TDateInterval. Can be called twice from the same print
 * statement, but not more.
 */
char* CxFormatDateInterval
(TDateInterval interval)
{
    return GtoFormatDateInterval(&interval);
}

/*
** calendar coercion uses the string version
*/
CxTCalendar* CxCalendarFromAlibHolidayName(char* holidayName)
{
    static char routine[] = "CxCalendarFromAlibHolidayName";

    THolidayList* hl;
    CxTCalendar*   cal = NULL;

    REQUIRE(holidayName != NULL);

    hl = GtoHolidayListFromCache(holidayName);
    if (hl == NULL) goto done; /* failure */

    /* THolidayList and CxTCalendar are identical structures */
    cal = CxCalendarCopy ((CxTCalendar*)(hl));
    if (cal == NULL) goto done; /* failure */

 done:

    if (cal == NULL) GtoErrMsgFailure(routine);
    return cal;
}

#if GTO_VERSION_NUM >= 1300
/*
***************************************************************************
** The ALIB custom zero curve interface is designed to allow a 3rd party
** curve to be easily represented via the encapsulated zero curve object.
**
** All you need to do is define the relevant C-functions in a custom zero
** curve type, and then create the zero curve object.
***************************************************************************
*/
static int CxCustomZeroCurveDiscount
(TCurve *zc,
 TDate   date,
 double *discount)
{
    *discount = CxZeroPrice (zc, date);
    return SUCCESS;
}

static TDate CxCustomZeroCurveBaseDate
(TCurve *zc)
{
    return zc->fBaseDate;
}

static TCustomZCType CxCustomZCType =
{
    sizeof(TCustomZCType),
    "CX_FLAT_FORWARDS",
    (TCustomZCDiscountFunc*) CxCustomZeroCurveDiscount,
    (TCustomZCBaseDateFunc*) CxCustomZeroCurveBaseDate,
    (TCustomZCCopyFunc*)     GtoCopyCurve,
    (TCustomZCFreeFunc*)     GtoFreeTCurve,
    NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL,
    NULL
};
#endif

/*f
***************************************************************************
** Converts a TCurve into an ALIB custom curve.
**
** The resulting TCurve will have fBaseDate = baseDate.
**
** This TCurve will have the correct interpolation if you use it using
** standard ALIB functions, but will not necessarily have good values in
** the dates and rates parameters.
**
** This curve will only be good for a single session since it does not
** have serialisation built in at present.
***************************************************************************
*/
TCurve* CxZeroCurveToAlibCustomCurve (TCurve *in)
{
    static char routine[] = "CxZeroCurveToAlibCustomCurve";

    TZeroCurve *zc;
    TCurve     *tc = NULL;

    REQUIRE (in != NULL);

#if GTO_VERSION_NUM > 1300
    zc = GtoCustomZCNew (&CxCustomZCType, in);
    if (zc == NULL) goto done; /* failure */
    tc = zc->curve;
    GtoFreeSafe (zc);
#else
    GtoErrMsg ("%s: Custom curve not supported in %s\n", routine, GTO_VERSION);
#endif

 done:

    if (tc == NULL) GtoErrMsgFailure (routine);
    return tc;
}


