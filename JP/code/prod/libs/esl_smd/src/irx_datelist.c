/*
***************************************************************************
** HEADER FILE: datelist.c
**
** Basic date and date list functions.
**
** A lot of the original code was imported from the ALIB - such code is of
** course beyond criticism.
**
** $Header$
***************************************************************************
*/

#include "irx/datelist.h"

#include <assert.h>

#include "irx/bsearch.h"
#include "irx/dateutils.h"
#include "irx/macros.h"

static int qsortDateCompare(const void *date1, const void *date2)
{
    IrxTDate d1 = *(IrxTDate*)date1;
    IrxTDate d2 = *(IrxTDate*)date2;

    if (d1 == d2) return 0;
    if (d1 < d2)  return -1;
    return 1;
}

/*f
***************************************************************************
** Sorts a IrxTDateList in place.
***************************************************************************
*/
void irxDateListSort(IrxTDateList* dl)
{
    if (dl != NULL && dl->fNumItems > 0 && dl->fArray != NULL)
    {
        qsort(dl->fArray,
              (size_t)dl->fNumItems,
              sizeof(IrxTDate),
              qsortDateCompare);
    }
}

/*f
***************************************************************************
** Removes duplicates from a IrxTDateList. This returns a new date list
** which will be sorted in ascending order.
***************************************************************************
*/
IrxTDateList* irxDateListCopyUnique(IrxTDateList* dl)
{
    static char routine[] = "irxDateListCopyUnique";
    int         status    = FAILURE;
    
    int i, j;
    IrxTDateList* dst = NULL;
    
    REQUIRE(dl != NULL);
    
    dst = irxDateListCopy(dl);
    if (dst == NULL) goto RETURN; /* failure */
    
    irxDateListSort(dst);
    j = 1; /* index in original list */
    i = 1;
    while (i < dst->fNumItems)
    {
        if (dst->fArray[i] == dst->fArray[j-1])
        {
            i += 1;
        }
        else
        {
            dst->fArray[j] = dst->fArray[i];
            j += 1;
            i += 1;
        }
    }
    /* j is one higher than the last element of dst that we used */
    dst->fNumItems = j;
    
    status = SUCCESS;
    
 RETURN:

    if (status != SUCCESS)
    {
        irxErrorFailure (routine);
        irxDateListFree (dst);
        dst = NULL;
    }
    return dst;
}
    
/*f
***************************************************************************
** Adds dates to a IrxTDateList. 
**
** If the original date list and date list to be added are sorted, then
** the resulting date list will be sorted and duplicate dates will be
** removed. Sorting assumes ascending order ([0] < [1] etc).
**
** If either of the inputs are not sorted, then the resulting date list
** will not be sorted, and some duplicates may remain.
**
** For efficiency, we do not automatically try to sort the resulting
** date list for unsorted inputs. Sorting the date list each time appears
** to be a huge performance issue in some algorithms (where the input
** dates would all be sorted anyway).
**
** Note that if dl is NULL, then this will create a new date list from
** the given dates.
**
** Note that if numItems=0, this will copy the given date list.
***************************************************************************
*/
IrxTDateList* irxDateListAddDates
(IrxTDateList *dl,          /* (I) Initial date list            */
 int         numItems,    /* (I) Number of dates to add       */
 IrxTDate      *array)      /* (I) [numItems] Dates to be added */
{
    static char routine[] = "irxDateListAddDates";
    int         status    = FAILURE;

    IrxTDateList tmp = {0, NULL};
    IrxTDateList *result = NULL;

    REQUIRE (numItems >= 0);
    REQUIRE (dl == NULL || dl->fNumItems >= 0);

    if (dl == NULL)
    {
        result = irxDateListMake (numItems, array);
    }
    else if (numItems <= 0)
    {
        result = irxDateListCopy (dl);
    }
    else if (dl->fNumItems == 0 && numItems == 0)
    {
        result = irxDateListMake (0, NULL);
    }
    else
    {
        int totalItems = dl->fNumItems + numItems;
        int i = 0;
        int j = 0;
        int k = 0;

        result = irxDateListMakeEmpty (totalItems);
        if (result == NULL) goto RETURN; /* failure */

        while (i < dl->fNumItems && j < numItems)
        {
            if (dl->fArray[i] == array[j])
            {
                /* exclude duplicates */
                ++j;
                --totalItems;
            }
            else if (dl->fArray[i] < array[j])
            {
                result->fArray[k] = dl->fArray[i];
                ++i;
                ++k;
            }
            else
            {
                assert (dl->fArray[i] > array[j]);
                result->fArray[k] = array[j];
                ++j;
                ++k;
            }
        }

        if (i < dl->fNumItems)
        {
            int n = dl->fNumItems - i;
            COPY_ARRAY (result->fArray+k, dl->fArray+i, IrxTDate, n);
            k += n;
        }

        if (j < numItems)
        {
            int n = numItems - j;
            COPY_ARRAY (result->fArray+k, array+j, IrxTDate, n);
            k += n;
        }

        assert (k == totalItems);
        result->fNumItems = totalItems;
    }
    if (result == NULL) goto RETURN; /* failure */

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
    {
        irxErrorFailure (routine);
        irxDateListFree (result);
        result = NULL;
    }

    IRX_FREE(tmp.fArray);
    
    return result;
}

/*f
***************************************************************************
** Adds dates to a IrxTDateList and frees the input date list.
**
** If the original date list and date list to be added are sorted, then
** the resulting date list will be sorted and duplicate dates will be
** removed. Sorting assumes ascending order ([0] < [1] etc).
**
** If either of the inputs are not sorted, then the resulting date list
** will not be sorted, and some duplicates may remain.
**
** For efficiency, we do not automatically try to sort the resulting
** date list for unsorted inputs. Sorting the date list each time appears
** to be a huge performance issue in some algorithms (where the input
** dates would all be sorted anyway).
**
** Note that if dl is NULL, then this will create a new date list from
** the given dates.
**
** Note that if numItems=0, this will copy the given date list.
**
** The input date list is FREE'd by this routine. Thus if you have an
** algorithm which involves building up a datelist gradually, you can
** do something like this:
**
**  IrxTDateList* dl = NULL;
**  ...
**  dl = irxDateListAddDatesFreeOld (dl, numItems, array);
**  if (dl == NULL) goto RETURN;
**  ..
**  dl = irxDateListAddDatesFreeOld (dl, numItems, array);
**  if (dl == NULL) goto RETURN;
**  ..
**  etc.
**
** with the point being that you don't have to worry about the original
** date list at each step since this routine frees it for you.
***************************************************************************
*/
IrxTDateList* irxDateListAddDatesFreeOld
(IrxTDateList *dl,         /* (I/O) Initial date list - gets freed */
 int         numItems,   /* (I) Number of dates to add           */
 IrxTDate     *array)      /* (I) [numItems] Dates to be added     */
{
    static char routine[] = "irxDateListAddDatesFreeOld";
    
    IrxTDateList *output;

    output = irxDateListAddDates (dl, numItems, array);
    irxDateListFree (dl);
    
    if (output == NULL) irxErrorFailure (routine);
    return output;
}

/*
***************************************************************************
** Searches a date list for the previous date and the next date.
** By definition:
**    previousDate <= date, maximum such date in the date list
**    nextDate     >  date, minimum such date in the date list
**
** If there is no such date, then returns zero for that date.
**
** You can state that the date list is sorted - in which case a binary
** search will be performed, or unsorted - in which case a linear search
** will be performed.
**
** You can pass NULL for the output pointers if you only care about
** either the previous date or the next date and not both.
***************************************************************************
*/
int irxDateListSearch(IrxTDateList* dl, IrxTBool sorted, IrxTDate date,
                     IrxTDate* prevDate, IrxTDate* nextDate)
{
    static char routine[] = "irxDateListSearch";
    int         status    = FAILURE;

    IrxTDate myPrevDate;
    IrxTDate myNextDate;
    int   i;

    if (sorted)
    {
        /* do the fastest possible binary search */
        long exact;
        long lo;
        long hi;

        if (irxBinarySearchLong (date,
                                 dl->fArray,
                                 sizeof(IrxTDate),
                                 dl->fNumItems,
                                 &exact,
                                 &lo,
                                 &hi) != SUCCESS)
        {
            goto RETURN; /* failure */
        }

        if (exact >= 0)
            myPrevDate = dl->fArray[exact];
        else if (lo >= 0)
            myPrevDate = dl->fArray[lo];
        else
            myPrevDate = 0;

        if (hi >= 0)
            myNextDate = dl->fArray[hi];
        else
            myNextDate = 0;
    }
    else
    {
        myPrevDate = 0;
        myNextDate = 0;
                
        for (i = 0; i < dl->fNumItems; ++i)
        {
            IrxTDate thisDate = dl->fArray[i];

            if (thisDate <= date)
            {
                if (thisDate > myPrevDate)
                    myPrevDate = thisDate;
            }
            else
            {
                if (thisDate > myNextDate && myNextDate != 0)
                    myNextDate = thisDate;
            }
        }
    }

    if (prevDate != NULL) *prevDate = myPrevDate;
    if (nextDate != NULL) *nextDate = myNextDate;

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}


/*f
***************************************************************************
** Truncates a date list at the specified date. The resulting date list
** will contain all dates previous to (or following) the specified date.
** Dates in the datelist which match the specified date may be optionally
** included.
**
** The datelist may optionally be modified in place or a new copy is
** returned.
**
** The input date list must be sorted.
***************************************************************************
*/
IrxTDateList* irxDateListTruncate
(IrxTDateList *dateList,       /* (I/O) Date list to be modified in place */
 IrxTDate      truncationDate, /* (I) Date on which to perform trunctation */
 IrxTBool      inclusive,      /* (I) TRUE=include truncation date if in list*/
 IrxTBool      excludeBefore,  /* (I) TRUE=exclude dates before truncation date*/
 IrxTBool      inPlace         /* (I) TRUE=modify date list in place */
)
{
    static char  routine[] = "irxDateListTruncate";
    IrxTDateList  *truncated = NULL;
    int          numItems;
    int          size;
    int          offset;
    int          truncatePt;
    int          i;

    REQUIRE (dateList != NULL);

    /*  First we find the truncation point in the date list and the size
     *  of the new date list
     */

    /* perhaps we could do this more efficiently with a binary search -
       but this is the code which was imported and therefore should work */
    numItems = dateList->fNumItems;
    if (excludeBefore)
    {
        truncatePt = 0;
        for (i = 0; i < numItems; i++)
        {
            if (dateList->fArray[i] > truncationDate)
            {
                truncatePt = i;
                break;
            }
            if (inclusive && dateList->fArray[i] == truncationDate)
            {
                truncatePt = i;
                break;
            }
        }
        size = numItems - truncatePt;
        offset = truncatePt;
    }
    else
    {
        truncatePt = numItems - 1;
        for (i = numItems - 1; i > 0; i--)
        {
            if (dateList->fArray[i] < truncationDate)
            {
                truncatePt = i;
                break;
            }
            if (inclusive && dateList->fArray[i] == truncationDate)
            {
                truncatePt = i;
                break;
            }
        }
        size = truncatePt + 1;
        offset = 0;
    }

    /* Next we get a pointer to the datelist where we store the result */
    if (inPlace)
    {
        truncated = dateList;
    }
    else
    {
        truncated = irxDateListMakeEmpty (size);
        if (truncated == NULL)
            goto RETURN;
    }

    /* finally we populate the result */
    if (!inPlace || offset != 0)
    {
        for (i = 0; i < size; i++)
        {
            truncated->fArray[i] = dateList->fArray[i + offset];
        }
    }
    truncated->fNumItems = size;

 RETURN:

    if (truncated == NULL)
        irxErrorFailure(routine);

    return(truncated);
}





/*f
***************************************************************************
** Makes a date list from a given start date to a given end date with dates
** seperated by a given interval.
**
** Use the stub parameter to determine whether the stub appears at the
** start or the end of the date list, and whether the stub is long or
** short.
**
** The start date and end date are always both in the date list.
** The end date must be strictly after the start date.
** The date interval must be positive.
***************************************************************************
*/
IrxTDateList* irxDateListMakeRegular
(IrxTDate         startDate,  /* (I) Start date */
 IrxTDate         endDate,    /* (I) End date */
 IrxTDateInterval interval,   /* (I) Date interval */
 IrxTStubLocation stubLocation)   /* (I) Stub type */
{
    static char routine[] = "irxDateListMakeRegular";

    IrxTDateList *dl = NULL;
    IrxTDate      dates[100];
    int         i;
    int         numIntervals;
    int         numDates = sizeof(dates) / sizeof(IrxTDate);
    IrxTDate      date;

    REQUIRE (interval.prd > 0);
    REQUIRE (endDate > startDate);

    /* we calculate dates in blocks of 100 and add to the datelist */
    switch (stubLocation)
    {
    case IRX_SHORT_FRONT_STUB:
    case IRX_LONG_FRONT_STUB:
        /* front stub - so we start at endDate and work backwards */
        numIntervals = 0;
        i            = numDates;
        date         = endDate;
        while (date > startDate)
        {
            if (i == 0)
            {
                dl = irxDateListAddDatesFreeOld (dl, numDates, dates);
                if (dl == NULL) goto RETURN; /* failure */
                i = numDates;
            }
            --i;
            --numIntervals;
            assert (i >= 0);
            dates[i] = date;
            if (irxDateAddMultiInterval (endDate,
                                         numIntervals,
                                         interval,
                                         &date) != SUCCESS)
                goto RETURN; /* failure */
        }
        assert (date <= startDate);
        if (date == startDate || stubLocation == IRX_SHORT_FRONT_STUB)
        {
            /* don't change existing dates[] but need to add startDate */
            if (i == 0)
            {
                dl = irxDateListAddDatesFreeOld (dl, numDates, dates);
                if (dl == NULL) goto RETURN; /* failure */
                i = numDates;
            }
            --i;
            dates[i] = startDate;
        }
        else
        {
            assert (stubLocation == IRX_LONG_FRONT_STUB);
            assert (date < startDate);
            /* the existing date in dates[] should be changed to be
               the start date */
            dates[i] = startDate;
        }
        /* now add from dates[i] to dates[numDates-1] to the date list */
        dl = irxDateListAddDatesFreeOld (dl, numDates-i, dates+i);
        if (dl == NULL) goto RETURN; /* failure */
        break;
    case IRX_SHORT_BACK_STUB:
    case IRX_LONG_BACK_STUB:
        /* back stub - so we start at startDate and work forwards */
        numIntervals = 0;
        i            = -1;
        date         = startDate;
        while (date < endDate)
        {
            ++i;
            if (i == numDates)
            {
                dl = irxDateListAddDatesFreeOld (dl, numDates, dates);
                if (dl == NULL) goto RETURN; /* failure */
                i = 0;
            }
            ++numIntervals;
            assert (i < numDates);
            dates[i] = date;
            if (irxDateAddMultiInterval (startDate,
                                         numIntervals,
                                         interval,
                                         &date) != SUCCESS)
                goto RETURN; /* failure */
        }
        assert (date >= endDate);
        if (date == endDate || stubLocation == IRX_SHORT_BACK_STUB)
        {
            /* don't change existing dates[] but need to add endDate */
            ++i;
            if (i == numDates)
            {
                dl = irxDateListAddDatesFreeOld (dl, numDates, dates);
                if (dl == NULL) goto RETURN; /* failure */
                i = 0;
            }
            dates[i] = endDate;
        }
        else
        {
            assert (stubLocation == IRX_LONG_BACK_STUB);
            assert (date > endDate);
            /* the existing date in dates[] should be changed to be
               the start date */
            dates[i] = endDate;
        }
        /* now add from dates[0] to dates[i] to the date list */
        dl = irxDateListAddDatesFreeOld (dl, i+1, dates);
        if (dl == NULL) goto RETURN; /* failure */
        break;
    }

 RETURN:

    if (dl == NULL)
        irxErrorFailure (routine);

    return dl;
}





/*f
***************************************************************************
** Makes a date list from a given start date to a given end date with dates
** seperated by a given interval.
**
** All dates are on cycle with the given roll date.
**
** The startDate and endDate are not necessarily in the date list.
** The first date will be <= startDate.
** The last date will be >= endDate.
** The end date must be strictly after the start date.
** The date interval must be positive.
***************************************************************************
*/
IrxTDateList* irxDateListMakeWithRoll
(IrxTDate         startDate,     /* (I) Start date */
 IrxTDate         endDate,       /* (I) End date */
 IrxTDate         rollDate,      /* (I) Roll date */
 IrxTDateInterval interval)      /* (I) Date interval */
{
    static char routine[] = "irxDateListMakeWithRoll";

    IrxTDateList *dl = NULL;

    IrxTDate      firstDate;
    int           i;
    int           numIntervals;
    int           extraDays;
    int           numDates;
    IrxTDate      date;

    REQUIRE (interval.prd > 0);
    REQUIRE (endDate > startDate);
    REQUIRE (rollDate != 0);
    REQUIRE (startDate != 0);
    REQUIRE (endDate != 0);
    REQUIRE (rollDate <= endDate);

    /* firstDate is the first date that is on cycle with rollDate and
       is on or before startDate - calculate this by counting the
       number of intervals from the rollDate - note that this routine
       is capable of counting backwards and guarantees the following:

       rollDate + numIntervals*interval + extraDays = startDate

       Hence startDate - extraDays is on cycle with rollDate and is
       either one date after or one date before startDate.
    */
    if (irxCountDateIntervals (rollDate,
                               startDate,
                               interval,
                               &numIntervals,
                               &extraDays) != SUCCESS)
        goto RETURN; /* failure */

    firstDate = startDate - extraDays;
    /* we expect this loop to be short (no steps or one step) */
    while (firstDate > startDate)
    {
        --numIntervals;
        if (irxDateAddMultiInterval (rollDate,
                                     numIntervals,
                                     interval,
                                     &firstDate) != SUCCESS)
            goto RETURN; /* failure */
    }

    /* Now count the number of intervals from firstDate to endDate */
    if (irxCountDateIntervals (firstDate,
                               endDate,
                               interval,
                               &numIntervals,
                               &extraDays) != SUCCESS)
        goto RETURN; /* failure */

    if (extraDays > 0)
        ++numIntervals;

    /* if we add numIntervals to firstDate then we are at a date on or
       after endDate - therefore the number of dates in our date list
       should be numIntervals+1 */

    numDates = numIntervals + 1;
    dl = irxDateListMakeEmpty(numDates);
    if (dl == NULL)
        goto RETURN; /* failure */
 
    dl->fArray[0] = firstDate;
    for (i = 1; i < numDates; ++i)
    {
        if (irxDateAddMultiInterval (firstDate,i,interval,&date) != SUCCESS)
            goto RETURN; /* failure */
        dl->fArray[i] = date;
    }

    if (dl->fArray[numDates-1] < endDate)
        PROGRAM_BUG();
    if (dl->fArray[numDates-2] >= endDate)
        PROGRAM_BUG();

 RETURN:

    if (dl == NULL)
        irxErrorFailure (routine);

    return dl;
}
