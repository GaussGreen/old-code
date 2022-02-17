/*
***************************************************************
**FILE NAME : cdscarry.h
**
**CDS Carry functions.
**
** $Header$
**************************************************************
*/

#ifndef CX_CDS_CARRY_H
#define CX_CDS_CARRY_H

#include "cx.h"


/*Calculates the carry of a CDS up to horizon date.
 * Returns a CdsCarryOutputs object.
 */

CxTCdsCarryOutputs* CxCdsCarry(
 TDate           today,
 TDate           valueDate,
 TDate           startDate,
 TDate           endDate,
 TDate            horizonDate,
 double          couponRate,
 TDateInterval  *dateInterval,
 CxTStubType      stubType,
 CxTDayCountConv  paymentDcc,
 CxTBadDayConv    badDayConv,
 CxTCalendar     *calendar,
 TCurve         *discCurve,
 TBoolean        protectStart
);

int CxCdsSlideBreakeven(
		      TDate horizonDate,
		      TCurve *discCurve, 
		      TDate startDate,
		      TDate maturityDate1,
		      TDate maturityDate2,
		      double signedNotional1,
		      double signedNotional2,
		      double spread1, 
		      double fee1,
		      double fee2,
		      CxTRecoveryCurve *recoveryCurve,
		      CxTCdsConventions   *cdsConventions, 
		      CxTCalendar *calendar,
		      double target,
		      double *spread2);



#endif
