/*
***************************************************************************
** FILE NAME: timeline.c
**
** Some internal routines for generating timelines.
**
** $Header$
***************************************************************************
*/

#include "timeline.h"

#include "creditcurve.h"

#include <cxutils/include/datelist.h>
#include <cxutils/include/cxmacros.h>

#include <alib/datelist.h>
#include <alib/tcurve.h>

/*f
** Returns a timeline for use with risky integrations assuming flat
** forward curves.
**
** The timeline will contain the following dates:
**
** startDate
** endDate
** all the points in the discount curve moved back by delay
** all the points in the risky curve
** all the points in the loss curve (can be NULL)
** nothing before startDate
** nothing after endDate
*/
TDateList* CxRiskyTimeLine
(TDate             startDate,
 TDate             endDate,
 TCurve*           discCurve,
 CxTCreditCurve*   spreadCurve,
 CxTRecoveryCurve* recoveryCurve,
 long              delay)
{
    static char routine[] = "CxRiskyTimeLine";

    int i;
    TDateList *tl  = NULL;
    TDate     *dates = NULL;

    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (delay >= 0);
    REQUIRE (endDate > startDate);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (spreadCurve->type == CX_CURVE_TYPE_FLOW); /* TBD */

    /*
    ** Timeline is points on the spreadCurve between startDate and endDate,
    ** combined with delay adjusted points from the discCurve, plus
    ** the startDate and endDate, plus the critical dates.
    */
    tl = GtoNewDateListFromTCurve (discCurve);
    if (tl == NULL) goto done;
    if (delay > 0)
        for (i = 0; i < tl->fNumItems; ++i)
            tl->fArray[i] -= delay;

    dates = GtoDatesFromCurve (spreadCurve->tc);
    tl = CxDateListAddDatesFreeOld (tl, spreadCurve->tc->fNumItems, dates);
    if (tl == NULL) goto done;
    if (recoveryCurve != NULL)
    {
        tl = CxDateListAddDatesFreeOld (tl, recoveryCurve->numItems,
                                        recoveryCurve->dates);
        if (tl == NULL) goto done;
    }
    tl = CxDateListAddDatesFreeOld (tl, 1, &startDate);
    if (tl == NULL) goto done;
    tl = CxDateListAddDatesFreeOld (tl, 1, &endDate);
    if (tl == NULL) goto done;

    /* remove dates strictly before startDate and strictly after endDate */
    tl = CxDateListTruncate (tl, startDate, TRUE, TRUE, TRUE);
    tl = CxDateListTruncate (tl, endDate, TRUE, FALSE, TRUE);
    
 done:

    if (tl == NULL)
        GtoErrMsgFailure (routine);

    FREE (dates);

    return tl;
}

/*f
** Truncate timeline.
**
** Truncates a timeline so that it will contain the following:
**
** startDate
** endDate
** criticalDates (can be NULL)
** nothing before startDate
** nothing after endDate
*/
TDateList* CxTruncateTimeLine
(TDateList* criticalDates,
 TDate      startDate,
 TDate      endDate)
{
    static char routine[] = "CxTruncateTimeLine";
    
    TDateList *tl = NULL;
    TDate      startEndDate[2];

    REQUIRE (endDate > startDate);

    startEndDate[0] = startDate;
    startEndDate[1] = endDate;

    tl = CxDateListAddDates (criticalDates, 2, startEndDate);
    if (tl == NULL) goto done; /* failure */

    /* remove dates strictly before startDate and strictly after endDate */
    tl = CxDateListTruncate (tl, startDate, TRUE, TRUE, TRUE);
    tl = CxDateListTruncate (tl, endDate, TRUE, FALSE, TRUE);
    
 done:

    if (tl == NULL)
        GtoErrMsgFailure (routine);
    
    return tl;
}
