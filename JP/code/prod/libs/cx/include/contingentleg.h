/*
***************************************************************************
** FILE NAME: contingentLeg.h
**
** Analytics for a contingent leg.
**
** $Header$
***************************************************************************
*/

#ifndef CX_CONTINGENT_LEG_H
#define CX_CONTINGENT_LEG_H

#include "cx.h"

/**
 * Computes the PV of a contingent leg as a whole
 *
 * For each payment period this is the integral of LGD(t) . Z(t) . dS/dt dt
 * where S is the survival function and LGD is the loss given default
 * function and Z is the discount function. Discounting is calculated at the
 * payment date and not at the observation date.
 *
 * If the recoveryCurve is provided, then this gives the recovery rate 
 * (given default) at different times.
 *
 * Otherwise a constant recovery rate of 0 is assumed - you can then multiply
 * the answer by (1-R) externally.
 */
int CxContingentLegPV
(CxTContingentLeg *cl,               /* (I) Contingent leg                  */
 TDate             today,            /* (I) No observations before today    */
 TDate             valueDate,        /* (I) Value date for discounting      */
 TCurve           *discountCurve,    /* (I) Risk-free curve                 */
 CxTCreditCurve   *spreadCurve,      /* (I) Spread curve                    */
 CxTRecoveryCurve *recoveryCurve,    /* (I) Recovery curve - can be NULL    */
 double           *pv);              /* (O) Present value of contingent leg */

/**
 * Computes the PV of a single payment of a contingent leg.
 *
 * Calling this function repeatedly with the correct sensible inputs would
 * give the same result as CxContingentLegPV on the whole contingent leg,
 * with the difference that this function gives the PV at today, whereas
 * CxContingentLegPV gives the PV at the value date (using the risk-free
 * rate to adjust from the PV at today).
 *
 * This is the integral of LGD(t) . Z(t) . dS/dt dt where S is the survival
 * function and LGD is the loss given default function and Z is the discount
 * function. Discounting is calculated at the payment date and not at the
 * observation date.
 *
 * If the recoveryCurve is provided, then this gives the recovery rate 
 * (given default) at different times.
 *
 * Otherwise a constant recovery rate of 0 is assumed - you can then multiply
 * the answer by (1-R) externally.
 */
int CxContingentPaymentPV
(CxTProtPayConv    payType,        /* (I) Contingent leg delay type       */
 TDate             today,          /* (I) PV to today                     */
 double            notional,       /* (I) Notional value                  */
 TDate             startDate,      /* (I) Observation start date          */
 TDate             endDate,        /* (I) Observation end date            */
 long              payDelay,       /* (I) Delay in payment                */
 TCurve           *discountCurve,  /* (I) Risk-free curve                 */
 CxTCreditCurve   *spreadCurve,    /* (I) Spread curve                    */
 CxTRecoveryCurve *recoveryCurve,  /* (I) Recovery curve - can be NULL    */
 double           *pv);            /* (O) Present value of contingent leg */

/**
 * Computes the expected cash flows for a contingent leg. This returns the
 * cash flows on the corresponding payment dates taking into account
 * the survival probability implied by the credit curve.
 *
 * The calculation involves computing the PV of each fee period, and then
 * forward valuing that value to the fee cash flow payment date. Although
 * this is not very sensitive to the discount curve, we therefore need
 * a discount curve. We also need a fee leg to provide the dates.
 *
 * The guarantee is that if you subsequently PV these flows using the
 * discount curve, then you should get the same results as
 * CxContingentLegPV.
 */
TCashFlowList* CxContingentLegExpectedFlows
(CxTContingentLeg *cl,
 TDateList        *flowDates,
 TDate             today,
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 CxTRecoveryCurve *recoveryCurve);

#endif

