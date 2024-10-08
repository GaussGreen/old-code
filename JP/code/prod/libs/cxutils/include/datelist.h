/*
***************************************************************************
** HEADER FILE: datelist.h
** 
** Defines the date list type used in credit exotics.
**
** $Header$
***************************************************************************
*/

#ifndef CX_DATELIST_H
#define CX_DATELIST_H

#include "cxutils.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*f
***************************************************************************
** Sorts a TDateList in place.
***************************************************************************
*/
void CxDateListSort(TDateList* dl /* (I/O) */);

/*f
***************************************************************************
** Copies a date list removing duplicate dates and sorting the resulting
** date list in ascending order.
***************************************************************************
*/
TDateList* CxDateListCopyUnique(TDateList* dl /* (I) */);

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
 TDate      *array);     /* (I) [numItems] Dates to be added */

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
 TDate     *array);     /* (I) [numItems] Dates to be added     */

/*f
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
int CxDateListSearch
(TDateList* dl,       /* (I) */
 TBoolean      sorted,   /* (I) */
 TDate      date,     /* (I) */
 TDate*     prevDate, /* (O) */ 
 TDate*     nextDate  /* (O) */
);

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
 TBoolean   inclusive,      /* (I) TRUE=include truncation date if in list*/
 TBoolean   excludeBefore,  /* (I) TRUE=exclude dates before truncation date*/
 TBoolean   inPlace         /* (I) TRUE=modify date list in place */
);

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
 CxTStubType    stubType);  /* (I) Stub type */

#ifdef __cplusplus
}
#endif

#endif

