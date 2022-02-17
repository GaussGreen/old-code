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

extern	DLL_EXPORT(int)	DrlDoubleArrayFloorIdx(double *y, int n,
				double x, int *j);
extern	DLL_EXPORT(int)	DrlTDateArrayFloorIdx(TDate *y, int n,
				TDate x, int *j);
extern	DLL_EXPORT(int)	DrlDoubleArrayCeilIdx(double *y, int n,
				double x, int *j);
extern	DLL_EXPORT(int)	DrlTDateArrayCeilIdx(TDate *y, int n,
				TDate x, int *j);
extern	DLL_EXPORT(int)	DrlDoubleArrayClosestIdx(double *y, int n,
				double x, int *j);
extern  DLL_EXPORT(int) DrlTDateArrayCeilIdx_Exclude(TDate *y, int n, 
                TDate x, int *j);


/*
 * Sorting
 */

extern	DLL_EXPORT(void) DrlQSort(
				void *base,
				size_t nmemb,
				size_t size,
        			int (*compar)(const void *, const void *));
extern	DLL_EXPORT(int)	DrlDoubleVectRevert(double *x, int n);
extern	DLL_EXPORT(int)	DrlDoubleVectSort(double *x, int n);
extern	DLL_EXPORT(int)	DrlDoubleVectLexSort(double tolerance,
				int removeDoubleItems, int *numItems,
				/* double*, ..., double*, NULL */ ...);



#endif	/* _drlsort_H */
