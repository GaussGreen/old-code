/****************************************************************************/
/*      Memory allocation routines for C++.                                 */
/****************************************************************************/
/*      alloc_plus.cpp                                                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C"
{
#include "fix123head.h"
}
#include "fix123plus.h"



double *Fix3_Alloc_Slice_CPP(
	FIX3_TREE_DATA			*tree_data)			// (I) Tree data
//////////////////////////////////////////////////////////////////////
// Allocate a slice, but throw an exception on failure
{
	double				*slice=NULL;		// Result slice

	// Allocate
	if ((slice = Fix3_Alloc_Slice(tree_data)) == NULL)
		throw "Fix3_Alloc_Slice_CPP: Failed to allocate slice.";

	// Return the slice
	return slice;
}


void *Fix3_DR_Array_CPP(
	int					type,				// (I) Type
    int					nl,					// (I) Lower bound
    int					nh)					// (I) Higher bound
//////////////////////////////////////////////////////////////////////
// Allocate an array, but throw an exception on failure
{
	void				*array=NULL;		// Result array

	// Allocate
	if ((array = DR_Array(type, nl, nh)) == NULL)
		throw "Fix3_DR_Array_CPP: Failed to allocate array.";

	// Return the array
	return array;
}

