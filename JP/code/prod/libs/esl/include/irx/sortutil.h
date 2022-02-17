/*
***************************************************************************
** HEADER FILE: sortutil.h
** CREATED BY:  Peter Taylor (13 July 2000)
**
** Various sort utilities.
***************************************************************************
*/

#ifndef IRX_SORTUTIL_H
#define IRX_SORTUTIL_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif

/**
***************************************************************************
** Sort order comparison function. This function will take the list of
** items to be sorted as its first argument, plus two positions within
** that list.
**
** If the arguments are called sortable, pos1 and pos2, then this function
** should return 0 if sortable[pos1] = sortable[pos2], should return 1 if
** sortable[pos1] > sortable[pos2] and should return -1 if 
** sortable[pos1] < sortable[pos2].
***************************************************************************
*/
typedef int (*IrxTSortOrderCompFunc) (void*, size_t, size_t);

/**
***************************************************************************
** Computes the sort order for a list of items in the case that the
** list to be sorted cannot be conveniently changed in place (and
** hence not suitable for direct use of qsort).
**
** For example we have a list of dates and another list of corresponding
** items. Sorting the dates on their own is useless without also putting
** the corresponding items in the same order.
**
** Returns an allocated array of size_t. This consists of an array
** of index values into the structure of size \texttt{arraySize},
** containing values [0,...,arraySize-1].
**
** Involves calls to \texttt{compFunc} using the \texttt{sortable} input
** to this routine as the first argument to \texttt{compFunc}. The
** comparison function should be capable of comparing elements of the
** \texttt{sortable} given two positions within the list.
***************************************************************************
*/
extern size_t* irxGetSortOrder
(void                  *sortable,
 IrxTSortOrderCompFunc  compFunc,
 size_t                 arraySize);

#ifdef __cplusplus
}
#endif

#endif    /* IRX_SORTUTIL_H */


