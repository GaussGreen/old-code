/*
***************************************************************************
** FILE NAME: feeLeg.h
**
** Analytics for a fee leg.
**
** $Header$
***************************************************************************
*/

#ifndef CX_FEE_LEG_H
#define CX_FEE_LEG_H

#include "cx.h"

/**
** Calculates the PV of a fee leg with fixed fee payments.
*/
int CxFeeLegPV
(CxTFeeLeg       *fl,
 TDate            today,
 TDate            valueDate,
 TCurve          *discCurve,
 CxTCreditCurve  *spreadCurve,
 TBoolean         payAccruedAtStart,
 double          *pv);

/**
** Calculates the PV of the accruals which occur on default with delay.
**
** We assume that the accruals are linear in time to make the integration
** easier.
*/
int CxAccrualOnDefaultPV
(TDate            today,
 TDate            startDate,
 TDate            endDate,
 long             delay,
 double           amount,
 TCurve          *discCurve,
 CxTCreditCurve  *spreadCurve,
 TBoolean         obsStartOfDay,
 double          *pv);

/**
 * Calculates the PV of a single fee payment.
 *
 * Calling this function repeatedly with sensible inputs and summing the
 * result would give the same answer as using CxFeeLegPV directly
 * with one difference - the pv returned by this function is for today,
 * whereas the pv returned by CxFeeLegPV is for valueDate.
 *
 * The conversion is a matter of dividing by the discount factor between
 * today and valueDate using the risk-free curve.
 */
int CxFeePaymentPV
(CxTAccrualPayConv accrualPayConv,
 TDate            today,
 TDate            accStartDate,
 TDate            accEndDate,
 TDate            payDate,
 CxTDayCountConv  accrueDCC,
 double           notional,
 double           couponRate,
 TCurve          *discCurve,
 CxTCreditCurve  *spreadCurve,
 TBoolean         obsStartOfDay,
 double          *pv);

/*f
** Calculates the PV of the accruals which occur on default with delay.
*/
int AccrualOnDefaultPVWithTimeLine
(TDate           today,
 TDate           startDate,
 TDate           endDate,
 long            delay,
 double          amount,
 TCurve         *discCurve,
 CxTCreditCurve *spreadCurve,
 TDateList      *tl,
 double         *pv);


/**
 * Computes the non-contingent cash flows for a fee leg. These are the
 * cash flows you will receive if there is no default.
 */
TCashFlowList* CxFeeLegFlows
(CxTFeeLeg      *fl);


/**
 * Computes the expected cash flows for a fee leg. This returns the cash
 * flows on the fee leg payment dates taking into account the survival
 * probability implied by the credit curve.
 *
 * The calculation involves computing the PV of each fee period, and then
 * forward valuing that value to the cash flow payment date. Although this
 * is not very sensitive to the discount curve, we therefore need a
 * discount curve.
 *
 * The guarantee is that if you subsequently PV these flows using the
 * discount curve, then you should get the same result as CxFeeLegPV.
 */
TCashFlowList* CxFeeLegExpectedFlows
(CxTFeeLeg      *fl,
 TDate           today,
 TCurve         *discCurve,
 CxTCreditCurve *creditCurve);


#endif

