#ifndef IRX_CFL_H
#define IRX_CFL_H

#include "irxflow.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
***************************************************************************
** Merges two cash flow lists without scaling.
***************************************************************************
*/
IrxTCashFlowList* irxCashFlowListMerge
(IrxTCashFlowList const *cfl1,
 IrxTCashFlowList const *cfl2);

/**
***************************************************************************
** Merges an array of cash flow lists, simultaneously scaling each cash
** flow list by a factor.
**
** The cash flow lists in the array can be NULL. If all the provided
** cash flow lists are NULL, then the output will also be NULL, but
** there will be no message in the error log.
**
** Failure is really only expected if one runs out of memory!
***************************************************************************
*/
IrxTCashFlowList* irxCashFlowListMergeAndScale
(long               arraySize,   /* (I) Size of arrays */
 IrxTCashFlowList const **cflArray,    /* (I) [arraySize] */
 double            *factors);    /* (I) [arraySize] */


/**
 * Copies and sorts the cash flow list. Repeating dates are removed from
 * the list and the amount summed. Zero amounts are not removed.
 */
IrxTCashFlowList *irxCashFlowListCopySort(
    const IrxTCashFlowList *cfl);

/**
 * Copies and scales the cash flow list.
 */
IrxTCashFlowList *irxCashFlowListCopyScale(
    const IrxTCashFlowList *cfl,
    double                  scale);

/**
 * Scales a cash flow list in place.
 */
void irxCashFlowListScale
(IrxTCashFlowList    *cfl,
 double               scale);

/**
 * Copies and strips out cash flows before a particular date.
 */
IrxTCashFlowList *irxCashFlowListCopyRemoveHistory(
    const IrxTCashFlowList *cfl,
    IrxTDate                minDate);

/**
 * Strip out cash flows before a particular date.
 */
void irxCashFlowListRemoveHistory
(IrxTCashFlowList *cfl,
 IrxTDate          minDate);

/**
 * Buckets cash flows for a particular range of dates from a cash flow list.
 * You will get cash flows for which startDate <= date <= endDate.
 *
 * To get cash flows on a particular date, set endDate=startDate.
 * If you put startDate=0 then you get all cash flows before endDate. 
 * If you put endDate=0, then you get all cash flows after startDate.
 *
 * If both are zero, then you get all cash flows.
 */
int irxCashFlowListBucket
(const IrxTCashFlowList *cfl,
 IrxTDate                startDate,
 IrxTDate                endDate,
 double                 *amount
);


/**
 * Splits a cash flow lists into cash flows on or before a particular date
 * and cash flows after that date. May return cash flow lists with zero
 * elements, but will always attempt to return something.
 */
int irxCashFlowListSplit
(IrxTCashFlowList const *cfl,
 IrxTDate                splitDate,
 IrxTCashFlowList      **before,
 IrxTCashFlowList      **after);


/**
 * Computes the PV of a cash flow list - this is calculated to the base date
 * of the zero curve and all cash flows are included whether in the past
 * or not.
 */
int irxCashFlowListPV
(const IrxTCashFlowList *cfl,
 const IrxTZeroCurve    *zc,
 double                 *pv);

/**
 * Computes the FV of a cash flow list - this is calculated to the given
 * value date, and all cash flows are included whether in the past or not.
 */
int irxCashFlowListFV
(const IrxTCashFlowList *cfl,
 const IrxTZeroCurve    *zc,
 IrxTDate                valueDate,
 double                *fv);

#ifdef __cplusplus
}
#endif

#endif
