#include "irx/cfl.h"

#include <assert.h>
#include <limits.h>

#include <irx/macros.h>
#include <irx/sortutil.h>

#include "irx/zerocurve.h"


/**
***************************************************************************
** Merges two cash flow lists without scaling.
***************************************************************************
*/
IrxTCashFlowList* irxCashFlowListMerge
(IrxTCashFlowList const *cfl1,
 IrxTCashFlowList const *cfl2)
{
    static char routine[] = "irxCashFlowListMerge";
    
    IrxTCashFlowList *cfl = NULL;

    IrxTCashFlowList const* array[2];

    REQUIRE (cfl1 != NULL);
    REQUIRE (cfl2 != NULL);

    array[0] = cfl1;
    array[1] = cfl2;

    cfl = irxCashFlowListMergeAndScale (2, array, NULL);

 RETURN:

    if (cfl == NULL)
        irxErrorFailure (routine);

    return cfl;
}


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
(long                     arraySize,   /* (I) Size of arrays */
 IrxTCashFlowList const **cflArray,    /* (I) [arraySize] */
 double                  *factors)     /* (I) [arraySize] */
{
    static char routine[] = "irxCashFlowListMergeAndScale";
    int         status    = FAILURE;
    
    IrxTCashFlowList *cfl = NULL;

    long i;

    REQUIRE (cflArray != NULL);

    for (i = 0; i < arraySize; ++i)
    {
        double            factor;
        const IrxTCashFlowList *add = cflArray[i];
        IrxTCashFlowList *tmp; /* Temporary result of adding this cfl */
        int               tmpSize;
        int               cflSize;
        int               addSize;
        int               iTmp, iCfl, iAdd; /* iterators */

        if (add == NULL)
            continue;

        if (factors == NULL)
            factor = 1.0;
        else
            factor = factors[i];

        if (IS_ZERO(factor))
            continue;
        
        if (cfl != NULL)
            cflSize = cfl->numItems;
        else
            cflSize = 0;
        
        addSize = add->numItems;
        tmpSize = addSize + cflSize;
        
        if (tmpSize == 0)
            continue;
        
        tmp = irxCashFlowListMakeEmpty (tmpSize);
        if (tmp == NULL)
            goto RETURN; /* failure */
        
        iTmp = iCfl = iAdd = 0;
        
        while (iCfl < cflSize && iAdd < addSize)
        {
            IrxTDate cflDate = cfl->dates[iCfl];
            IrxTDate addDate = add->dates[iAdd];
            double   cflAmount = cfl->amounts[iCfl];
            double   addAmount = add->amounts[iAdd];
                
            if (cflDate < addDate)
            {
                tmp->dates[iTmp]   = cflDate;
                tmp->amounts[iTmp] = cflAmount;
                ++iCfl;
            }
            else if (cflDate > addDate)
            {
                tmp->dates[iTmp]   = addDate;
                tmp->amounts[iTmp] = factor * addAmount;
                ++iAdd;
            }
            else
            {
                assert (cflDate == addDate);
                tmp->dates[iTmp]   = cflDate;
                tmp->amounts[iTmp] = cflAmount + factor * addAmount;
                ++iAdd;
                ++iCfl;
            }
            ++iTmp;
        }
        
        while (iCfl < cflSize)
        {
            tmp->dates[iTmp]   = cfl->dates[iCfl];
            tmp->amounts[iTmp] = cfl->amounts[iCfl];
            ++iCfl;
            ++iTmp;
        }
        
        while (iAdd < addSize)
        {
            tmp->dates[iTmp]   = add->dates[iAdd];
            tmp->amounts[iTmp] = factor * add->amounts[iAdd];
            ++iAdd;
            ++iTmp;
        }
        
        tmp->numItems = iTmp;
        
        irxCashFlowListFree (cfl);
        cfl = tmp;
    }

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
    {
        irxCashFlowListFree(cfl);
        cfl = NULL;
        irxErrorFailure(routine);
    }
    return cfl;
}


/*
 * Comparison routine required for irxGetSortOrder
 */
static int compareDatesInList (const IrxTDate *dates, size_t i, size_t j)
{
    if (dates[i] == dates[j]) return 0;
    if (dates[i] < dates[j]) return -1;
    return 1;
}

/*
 * Copies and sorts the cash flow list. Repeating dates are removed from
 * the list and the amount summed. Zero amounts are not removed.
 */
IrxTCashFlowList *irxCashFlowListCopySort(
    const IrxTCashFlowList *cfl)
{
    static char routine[] = "irxCashFlowListCopySort";
    int         status    = FAILURE;
        
    /* I have made it difficult for myself by the structure of
       the IrxTCashFlowList - i.e. dates and amounts in separate
       arrays. Hence I cannot use qsort directly. However there is
       a little known ALIB routine which handles this sort of thing
       which I have adapted.
    */
    size_t *sortOrder = NULL;
    int i,j,k;

    IrxTCashFlowList *copy = NULL;

    REQUIRE (cfl != NULL);

    copy = irxCashFlowListMakeEmpty (cfl->numItems);
    if (copy == NULL)
        goto RETURN; /* failure */

    if (cfl->numItems > 0)
    {
        sortOrder = irxGetSortOrder (cfl->dates,
                                     (IrxTSortOrderCompFunc)compareDatesInList,
                                     cfl->numItems);
        if (sortOrder == NULL)
            goto RETURN; /* failure */
        
        /* sortOrder gives the sort order of the dates in cfl->dates.
           However we also want to sum amounts for the same date.
           This means that we need three iterators.
           
           i: Iterator through the sortOrder
           j: Corresponding value from cfl extracted from sortOrder
           k: Insertion position back into copy.
        */
        k = 0;
        for (i = 0; i < cfl->numItems; ++i)
        {
            j = sortOrder[i];
            if (k > 0 && (cfl->dates[j] == copy->dates[k-1]))
            {
                copy->amounts[k-1] += cfl->amounts[j];
            }
            else
            {
                copy->dates[k]   = cfl->dates[j];
                copy->amounts[k] = cfl->amounts[j];
                ++k;
            }
        }
        copy->numItems = k;
    }

    status = SUCCESS;

 RETURN:

    FREE(sortOrder);
    if (status != SUCCESS)
    {
        irxCashFlowListFree (copy);
        copy = NULL;
        irxErrorFailure (routine);
    }
    return copy;
}

/**
 * Copies and scales the cash flow list.
 */
IrxTCashFlowList *irxCashFlowListCopyScale(
    const IrxTCashFlowList *cfl,
    double                  scale)
{
    IrxTCashFlowList *copy = irxCashFlowListCopy(cfl);
    if (copy != NULL)
        irxCashFlowListScale(copy, scale);

    return copy;
}


/**
 * Scales a cash flow list in place.
 */
void irxCashFlowListScale
(IrxTCashFlowList    *cfl,
 double               factor)
{
    if (cfl != NULL)
    {
        int i;
        for (i = 0; i < cfl->numItems; ++i)
            cfl->amounts[i] *= factor;
    }
}



/**
 * Copies and strips out cash flows before a particular date.
 */
IrxTCashFlowList *irxCashFlowListCopyRemoveHistory(
    const IrxTCashFlowList *cfl,
    IrxTDate                minDate)
{
    IrxTCashFlowList *copy = irxCashFlowListCopy(cfl);
    if (copy != NULL)
        irxCashFlowListRemoveHistory(copy, minDate);

    return copy;
}

/*
***************************************************************************
** Buckets cash flows for a particular range of dates from a cash flow list.
** You will get cash flows for which startDate <= date <= endDate.
**
** To get cash flows on a particular date, set endDate=startDate.
** If you put startDate=0 then you get all cash flows before endDate. 
** If you put endDate=0, then you get all cash flows after startDate.
**
** If both are zero, then you get all cash flows period.
***************************************************************************
*/
int irxCashFlowListBucket
(const IrxTCashFlowList *cfl,
 IrxTDate                startDate,
 IrxTDate                endDate,
 double                 *amount
)
{
    static char routine[] = "irxCashFlowListBucket";
    int         status    = FAILURE;

    int idx;

    REQUIRE (cfl != NULL);
    REQUIRE (amount != NULL);

    *amount = 0.0;

    if (endDate == 0)
        endDate = LONG_MAX;

    for (idx = 0; idx < cfl->numItems; ++idx)
    {
        if (cfl->dates[idx] >= startDate &&
            cfl->dates[idx] <= endDate)
        {
            *amount += cfl->amounts[idx];
        }
    }

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}




/**
 * Splits a cash flow lists into cash flows on or before a particular date
 * and cash flows after that date. May return cash flow lists with zero
 * elements, but will always attempt to return something.
 */
int irxCashFlowListSplit
(IrxTCashFlowList const *cfl,
 IrxTDate                splitDate,
 IrxTCashFlowList      **before,
 IrxTCashFlowList      **after)
{
    static char routine[] = "irxCashFlowListSplit";
    int         status    = FAILURE;

    int i,j,k;
    
    IrxTCashFlowList *myBefore = NULL;
    IrxTCashFlowList *myAfter  = NULL;

    REQUIRE (cfl != NULL);
    REQUIRE (before != NULL);
    REQUIRE (after != NULL);

    myBefore = irxCashFlowListMakeEmpty (cfl->numItems);
    myAfter  = irxCashFlowListMakeEmpty (cfl->numItems);
    if (myBefore == NULL || myAfter == NULL)
        goto RETURN; /* failure */

    j = 0;
    k = 0;
    for (i = 0; i < cfl->numItems; ++i)
    {
        if (cfl->dates[i] <= splitDate)
        {
            myBefore->dates[j]   = cfl->dates[i];
            myBefore->amounts[j] = cfl->amounts[i];
            ++j;
        }
        else
        {
            myAfter->dates[k]   = cfl->dates[i];
            myAfter->amounts[k] = cfl->amounts[i];
            ++k;
        }
    }
    myBefore->numItems = j;
    myAfter->numItems  = k;

    *before = myBefore;
    *after  = myAfter;

    myBefore = NULL;
    myAfter  = NULL;
    
    status = SUCCESS;

 RETURN:

    irxCashFlowListFree (myAfter);
    irxCashFlowListFree (myBefore);

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}






/**
 * Strip out cash flows before a particular date.
 */
void irxCashFlowListRemoveHistory
(IrxTCashFlowList *cfl,
 IrxTDate          minDate)
{
    if (cfl != NULL)
    {
        int i, j;

        j = 0;
        for (i = 0; i < cfl->numItems; ++i)
        {
            if (cfl->dates[i] >= minDate)
            {
                if (j < i)
                {
                    cfl->dates[j]   = cfl->dates[i];
                    cfl->amounts[j] = cfl->amounts[i];
                }
                ++j;
            }
        }
        cfl->numItems = j;
    }
}




/**
 * Computes the PV of a cash flow list - this is calculated to the base date
 * of the zero curve and all cash flows are included whether in the past
 * or not.
 */
int irxCashFlowListPV
(const IrxTCashFlowList *cfl,
 const IrxTZeroCurve    *zc,
 double                 *pv)
{
    static char routine[] = "irxCashFlowListPV";
    int         status    = FAILURE;

    int i;
    double myPv;

    REQUIRE (cfl != NULL);
    REQUIRE (zc != NULL);
    REQUIRE (pv != NULL);

    myPv = 0.0;
    for (i = 0; i < cfl->numItems; ++i)
    {
        double thisPv;

        if (irxZeroPrice(zc, cfl->dates[i], &thisPv) != SUCCESS)
            goto RETURN; /* failure */

        thisPv *= cfl->amounts[i];
        myPv   += thisPv;
    }

    *pv = myPv;
    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure(routine);

    return status;
}



/**
 * Computes the FV of a cash flow list - this is calculated to the given
 * value date, and all cash flows are included whether in the past or not.
 */
int irxCashFlowListFV
(const IrxTCashFlowList *cfl,
 const IrxTZeroCurve    *zc,
 IrxTDate                valueDate,
 double                 *fv)
{
    static char routine[] = "irxCashFlowListFV";
    int         status    = FAILURE;

    int i;
    double myFv;

    REQUIRE (cfl != NULL);
    REQUIRE (zc != NULL);
    REQUIRE (fv != NULL);

    myFv = 0.0;
    for (i = 0; i < cfl->numItems; ++i)
    {
        double thisFv;

        if (irxFwdZeroPrice(zc, valueDate, cfl->dates[i], &thisFv) != SUCCESS)
            goto RETURN; /* failure */

        thisFv *= cfl->amounts[i];
        myFv   += thisFv;
    }

    *fv = myFv;
    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure(routine);

    return status;
}










