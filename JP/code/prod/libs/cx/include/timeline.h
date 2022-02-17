/*
***************************************************************************
** FILE NAME: timeline.h
**
** Some internal routines for generating timelines.
**
** $Header$
***************************************************************************
*/

#ifndef CX_TIMELINE_H
#define CX_TIMELINE_H

#include "cx.h"

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
 CxTCreditCurve*   riskyCurve,
 CxTRecoveryCurve* recoveryCurve,
 long              delay);

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
 TDate      endDate);

#endif



