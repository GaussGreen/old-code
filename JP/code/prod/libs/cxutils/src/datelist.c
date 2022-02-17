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

#include "datelist.h"

#include <assert.h>

#include "bsearch.h"
#include "dateutils.h"
#include "cxmacros.h"

static int qsortDateCompare(const void *date1, const void *date2)
{
    TDate d1 = *(TDate*)date1;
    TDate d2 = *(TDate*)date2;

    if (d1 == d2) return 0;
    if (d1 < d2)  return -1;
    return 1;
}

/*f
***************************************************************************
** Sorts a TDateList in place.
***************************************************************************
*/
void CxDateListSort(TDateList* dl)
{
    if (dl != NULL && dl->fNumItems > 0 && dl->fArray != NULL)
    {
        qsort(dl->fArray,
              (size_t)dl->fNumItems,
              sizeof(TDate),
              qsortDateCompare);
    }
}

/*f
***************************************************************************
** Removes duplicates from a TDateList. This returns a new date list
** which will be sorted in ascending order.
***************************************************************************
*/
TDateList* CxDateListCopyUnique(TDateList* dl)
{
    static char routine[] = "CxDateListCopyUnique";
    int         status    = FAILURE;
    
    int i, j;
    TDateList* dst = NULL;
    
    REQUIRE(dl != NULL);
    
    dst = GtoCopyDateList(dl);
    if (dst == NULL) goto done; /* failure */
    
    CxDateListSort(dst);
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
    
 done:

    if (status != SUCCESS)
    {
        GtoErrMsgFailure (routine);
        GtoFreeDateList (dst);
        dst = NULL;
    }
    return dst;
}
    
/*f
***************************************************************************
** Adds dates to a TDateList. 
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
TDateList* CxDateListAddDates
(TDateList *dl,          /* (I) Initial date list            */
 int         numItems,    /* (I) Number of dates to add       */
 TDate      *array)      /* (I) [numItems] Dates to be added */
{
    static char routine[] = "CxDateListAddDates";
    int         status    = FAILURE;

    TDateList tmp = {0, NULL};
    TDateList *result = NULL;

    REQUIRE (numItems >= 0);
    REQUIRE (dl == NULL || dl->fNumItems >= 0);

    if (dl == NULL)
    {
        result = GtoNewDateListFromDates (array, numItems);
    }
    else if (numItems <= 0)
    {
        result = GtoCopyDateList (dl);
    }
    else if (dl->fNumItems == 0 && numItems == 0)
    {
        result = GtoNewDateListFromDates (NULL, 0);
    }
    else
    {
        int totalItems = dl->fNumItems + numItems;
        int i = 0;
        int j = 0;
        int k = 0;

        result = GtoNewEmptyDateList (totalItems);
        if (result == NULL) goto done; /* failure */

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
            COPY_ARRAY (result->fArray+k, dl->fArray+i, TDate, n);
            k += n;
        }

        if (j < numItems)
        {
            int n = numItems - j;
            COPY_ARRAY (result->fArray+k, array+j, TDate, n);
            k += n;
        }

        assert (k == totalItems);
        result->fNumItems = totalItems;
    }
    if (result == NULL) goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        GtoErrMsgFailure (routine);
        GtoFreeDateList (result);
        result = NULL;
    }

    CX_FREE(tmp.fArray);
    
    return result;
}

/*f
***************************************************************************
** Adds dates to a TDateList and frees the input date list.
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
**  TDateList* dl = NULL;
**  ...
**  dl = CxDateListAddDatesFreeOld (dl, numItems, array);
**  if (dl == NULL) goto done;
**  ..
**  dl = CxDateListAddDatesFreeOld (dl, numItems, array);
**  if (dl == NULL) goto done;
**  ..
**  etc.
**
** with the point being that you don't have to worry about the original
** date list at each step since this routine frees it for you.
***************************************************************************
*/
TDateList* CxDateListAddDatesFreeOld
(TDateList *dl,         /* (I/O) Initial date list - gets freed */
 int         numItems,   /* (I) Number of dates to add           */
 TDate     *array)      /* (I) [numItems] Dates to be added     */
{
    static char routine[] = "CxDateListAddDatesFreeOld";
    
    TDateList *output;

    output = CxDateListAddDates (dl, numItems, array);
    GtoFreeDateList (dl);
    
    if (output == NULL) GtoErrMsgFailure (routine);
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
int CxDateListSearch(TDateList* dl, TBoolean sorted, TDate date,
                     TDate* prevDate, TDate* nextDate)
{
    static char routine[] = "CxDateListSearch";
    int         status    = FAILURE;

    TDate myPrevDate;
    TDate myNextDate;
    int   i;

    if (sorted)
    {
        /* do the fastest possible binary search */
        long exact;
        long lo;
        long hi;

        if (CxBinarySearchLong (date,
                                dl->fArray,
                                sizeof(TDate),
                                dl->fNumItems,
                                &exact,
                                &lo,
                                &hi) != SUCCESS)
        {
            goto done; /* failure */
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
            TDate thisDate = dl->fArray[i];

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

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

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
TDateList* CxDateListTruncate
(TDateList *dateList,       /* (I/O) Date list to be modified in place */
 TDate      truncationDate, /* (I) Date on which to perform trunctation */
 TBoolean      inclusive,      /* (I) TRUE=include truncation date if in list*/
 TBoolean      excludeBefore,  /* (I) TRUE=exclude dates before truncation date*/
 TBoolean      inPlace         /* (I) TRUE=modify date list in place */
)
{
    static char  routine[] = "CxDateListTruncate";
    TDateList  *truncated = NULL;
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
        truncated = GtoNewEmptyDateList (size);
        if (truncated == NULL)
            goto done;
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

 done:

    if (truncated == NULL)
        GtoErrMsgFailure(routine);

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
TDateList* CxDateListMakeRegular
(TDate          startDate,  /* (I) Start date */
 TDate          endDate,    /* (I) End date */
 TDateInterval *interval,   /* (I) Date interval */
 CxTStubType    stubType)   /* (I) Stub type */
{
    static char routine[] = "CxDateListMakeRegular";
    int         status = FAILURE;

    TDateList *dl = NULL;
    TDate      tmpDates[100];
    int        i;
    int        numIntervals;
    int        numTmpDates = sizeof(tmpDates) / sizeof(TDate);
    int        totalDates = 0;
    TDate      date;

    REQUIRE (interval != NULL);
    REQUIRE (interval->prd > 0);
    REQUIRE (endDate > startDate);

    /* we calculate tmpDates in blocks of 100 and add to the datelist */
    switch (stubType)
    {
    case CX_SHORT_FRONT_STUB:
    case CX_LONG_FRONT_STUB:
        /* front stub - so we start at endDate and work backwards */
        numIntervals = 0;
        i            = numTmpDates;
        date         = endDate;
        while (date > startDate)
        {
            if (i == 0)
            {
                dl = CxDateListAddDatesFreeOld (dl, numTmpDates, tmpDates);
                if (dl == NULL) goto done; /* failure */
                i = numTmpDates;
            }
            --i;
            --numIntervals;
            ++totalDates;
            assert (i >= 0);
            tmpDates[i] = date;
            if (CxDateAddMultiInterval (endDate,
                                        numIntervals,
                                        interval,
                                        &date) != SUCCESS)
                goto done; /* failure */
        }
        assert (totalDates > 0);
        assert (date <= startDate);
        if (date == startDate ||
            totalDates == 1 ||
            stubType == CX_SHORT_FRONT_STUB)
        {
            /* don't change existing tmpDates[] but need to add startDate */
            if (i == 0)
            {
                dl = CxDateListAddDatesFreeOld (dl, numTmpDates, tmpDates);
                if (dl == NULL) goto done; /* failure */
                i = numTmpDates;
            }
            --i;
            ++totalDates;
            tmpDates[i] = startDate;
        }
        else
        {
            assert (stubType == CX_LONG_FRONT_STUB);
            assert (date < startDate);
            /* the existing date in tmpDates[] should be changed to be
               the start date */
            tmpDates[i] = startDate;
        }
        /* now add from tmpDates[i] to tmpDates[numTmpDates-1] to date list */
        dl = CxDateListAddDatesFreeOld (dl, numTmpDates-i, tmpDates+i);
        if (dl == NULL) goto done; /* failure */
        break;
    case CX_SHORT_BACK_STUB:
    case CX_LONG_BACK_STUB:
        /* back stub - so we start at startDate and work forwards */
        numIntervals = 0;
        i            = -1;
        date         = startDate;
        while (date < endDate)
        {
            ++i;
            ++totalDates;
            if (i == numTmpDates)
            {
                dl = CxDateListAddDatesFreeOld (dl, numTmpDates, tmpDates);
                if (dl == NULL) goto done; /* failure */
                i = 0;
            }
            ++numIntervals;
            assert (i < numTmpDates);
            tmpDates[i] = date;
            if (CxDateAddMultiInterval (startDate,
                                        numIntervals,
                                        interval,
                                        &date) != SUCCESS)
                goto done; /* failure */
        }
        assert (totalDates > 0);
        assert (date >= endDate);
        if (date == endDate ||
            totalDates == 1 ||
            stubType == CX_SHORT_BACK_STUB)
        {
            /* don't change existing tmpDates[] but need to add endDate */
            ++i;
            ++totalDates;
            if (i == numTmpDates)
            {
                dl = CxDateListAddDatesFreeOld (dl, numTmpDates, tmpDates);
                if (dl == NULL) goto done; /* failure */
                i = 0;
            }
            tmpDates[i] = endDate;
        }
        else
        {
            assert (stubType == CX_LONG_BACK_STUB);
            assert (date > endDate);
            /* the existing date in tmpDates[] should be changed to be
               the end date */
            tmpDates[i] = endDate;
        }
        /* now add from tmpDates[0] to tmpDates[i] to the date list */
        dl = CxDateListAddDatesFreeOld (dl, i+1, tmpDates);
        if (dl == NULL) goto done; /* failure */
        break;
    }
    ASSERT (totalDates >= 2);
    ASSERT (dl->fNumItems == totalDates);
    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        GtoFreeDateList (dl);
        dl = NULL;
        GtoErrMsgFailure (routine);
    }

    return dl;
}

