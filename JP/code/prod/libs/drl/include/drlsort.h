/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	SORT - Sorting Routines
 * File:	drlsort.h
 * Function:	Indexing and sorting routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlsort_H
#define	_drlsort_H
#include "drlstd.h"

/*
 * Searching ordered arrays
 */

extern	int	DrlDoubleArrayFloorIdx(double *y, int n,
				double x, int *j);
extern	int	DrlDDateArrayFloorIdx(DDate *y, int n,
				DDate x, int *j);
extern	int	DrlDoubleArrayCeilIdx(double *y, int n,
				double x, int *j);
extern	int	DrlDDateArrayCeilIdx(DDate *y, int n,
				DDate x, int *j);
extern	int	DrlDoubleArrayClosestIdx(double *y, int n,
				double x, int *j);


/*
 * Sorting
 */

extern	void DrlQSort(
				void *base,
				size_t nmemb,
				size_t size,
        			int (*compar)(const void *, const void *));
extern	int	DrlDoubleVectRevert(double *x, int n);
extern	int	DrlDoubleVectSort(double *x, int n);
extern	int	DrlDoubleVectLexSort(double tolerance,
				int removeDoubleItems, int *numItems,
				/* double*, ..., double*, NULL */ ...);



#endif	/* _drlsort_H */
