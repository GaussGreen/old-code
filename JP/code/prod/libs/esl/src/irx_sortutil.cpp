/*
***************************************************************************
** SOURCE FILE: sortutil.c
**
** CREATED BY:  Peter Taylor (13 July 2000)
**
** Various sort utilities.
***************************************************************************
*/

#include "irx/sortutil.h"

#include <stdlib.h>     /* qsort */

#include "irx/convert.h"
#include "irx/error.h"      /* irxError */
#include "irx/macros.h"     /* NEW_ARRAY etc */

static int compare(const void* p1, const void* p2);

typedef struct _SORT_INDEX
{
    size_t                idx;
    IrxTSortOrderCompFunc compFunc;
    void                 *data;
} SORT_INDEX;

/*
***************************************************************************
** Computes the sort order for a list of items in the case that the
** list to be sorted cannot be conveniently changed in place (and
** hence not suitable for direct use of qsort).
**
** For example we have a list of dates and another list of corresponding
** items. Sorting the dates on their own is useless without also putting
** the corresponding items in the same order.
**
** Returns the sortOrder as an array of size_t. This consists of an array
** of index values into the structure of size arraySize, containing values
** [0,...,arraySize-1].
**
** Involves calls to \texttt{compFunc} using the \texttt{sortable} input
** to this routine as the first argument to \texttt{compFunc}. The
** comparison function should be capable of comparing elements of the
** \texttt{sortable} given two positions within the list.
***************************************************************************
*/
size_t* irxGetSortOrder
(void                  *data,
 IrxTSortOrderCompFunc  compFunc,
 size_t                 arraySize
)
{
    static char routine[] = "irxGetSortOrder";
    int         status    = FAILURE;

    size_t     *sortOrder = NULL; /* to be returned */

    size_t      idx;
    SORT_INDEX *sortIndex = NULL;

    sortIndex = NEW_ARRAY (SORT_INDEX, arraySize);
    if (sortIndex == NULL)
        goto RETURN; /* failure */

    for (idx = 0; idx < arraySize; ++idx)
    {
        sortIndex[idx].idx      = idx;
        sortIndex[idx].compFunc = compFunc;
        sortIndex[idx].data     = data;
    }

    qsort (sortIndex, arraySize, sizeof(SORT_INDEX), compare);

    sortOrder = EXTRACT_ARRAY(size_t, SORT_INDEX, sortIndex, idx, arraySize);

    if (sortOrder == NULL)
        goto RETURN; /* failure */

    status = SUCCESS;

 RETURN:

    FREE_ARRAY (sortIndex);
    if (status != SUCCESS)
    {
        irxErrorFailure (routine);
        FREE (sortOrder);
        sortOrder = NULL;
    }

    return sortOrder;
}


/*
***************************************************************************
** FUNCTION: compare
** AUTHOR:   Peter Taylor (13 July 2000)
**
** qsort comparison function.
***************************************************************************
*/
static int compare(const void* p1, const void* p2)
{
    SORT_INDEX *ix1 = (SORT_INDEX*)p1;
    SORT_INDEX *ix2 = (SORT_INDEX*)p2;

    return (ix1->compFunc) (ix1->data, ix1->idx, ix2->idx);
}
