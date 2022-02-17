/*
***************************************************************************
** FILE NAME: cdsw.h
**
** CDSW functions. These are functions which emulate the CDSW page on
** Bloomberg.
**
** $Header$
***************************************************************************
*/

#ifndef CX_CDSW_H
#define CX_CDSW_H

#include "cx.h"

/*
 * Computes the upfront charge for a flat spread par curve.
 */
int CxCdswUpfrontCharge
(TDate           today,
 TDate           valueDate,
 TDate           startDate,
 TDate           endDate,
 double          couponRate,
 TBoolean        payAccruedOnDefault,
 TDateInterval  *dateInterval,
 CxTStubType     stubType,
 CxTDayCountConv accrueDCC,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar,
 TCurve         *discCurve,
 double          bbgSpread,
 double          recoveryRate,
 TBoolean        inclusiveProtection,
 TBoolean        payAccruedAtStart,
 double         *upfrontCharge);

/*
 * Computes the flat spread required to match the upfront charge.
 */
int CxCdswSpread
(TDate           today,
 TDate           valueDate,
 TDate           startDate,
 TDate           endDate,
 double          couponRate,
 TBoolean        payAccruedOnDefault,
 TDateInterval  *dateInterval,
 CxTStubType     stubType,
 CxTDayCountConv accrueDCC,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar,
 TCurve         *discCurve,
 double          upfrontCharge,
 double          recoveryRate,
 TBoolean        inclusiveProtection,
 TBoolean        payAccruedAtStart,
 double         *bbgSpread);

#endif

