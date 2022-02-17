/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	dvtspar.h
 * Function:	
 * Author:	
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_dvtspar_H
#define	_dvtspar_H
#include "kvtree.h"


/**
 * Virtual Tree dynamic time slice formula parser tool.
 * WARNING: since the routine uses static variables, the code
 * is NOT reentrant.
 * @param vt the virtual tree,
 * @param tpIdx timepoint index at which evaluated,
 * @param nx the number of X-variables in the formula.
 * @param x array of length "nx" of the timeslices to be substituted
 *          to x-variables.
 * @param formula a formula of variables labelled
 *   "x0", "x1",\dots,"xNX-1","y0",\dots,"yNY-2" that contains<br>
 *    numerical constants, such as 1.0 or 4.56e-2,<br>
 *    basic arithmetic operations "+", "-", "*" and "/"
 *    MIN(E1,E2) and MAX(E1,E2)
 *    IFPOS(E0,E1,E2) statement which returns
 *           the expression E1 (resp. E2) if expression E0 is
 *           positive (resp.negative).
 * @param constTable The table of constants as a map
 *           between string and double values
 * @param retVal on successful exit, contains the timeslice obtained
 *    by evaluating the formula by replacing the formal arguments
 *    "x0", "x1",\dots,"xNX-1", "y0", "y1",\dots,"yNY-1",
 *     by the values given in the input array
 *   "x".
 */
void	KTSliceParEval(
	KVTree &vt,		// (I) virtual tree 
	int nx,			// (I) # of X args
	KTSlice *x,		// (I) arrays of X args timeslices
	const char *formula,	// (I) formula
	KMap(String,double)& constTable,	// (I) constants table
	KTSlice *retVal);	// (O) return value timeslice

/*
 * Variable argument call of
 * @see KTSliceParEval.
 */

void	KTSliceParEvalV(
	KVTree &vt,		// (I) virtual tree 
	const char *formula,	// (I) formula
	KMap(String,double)& constTable,	// (I) constants table
	KTSlice *retVal,	// (O) return value timeslice
	int nx,			// (I) # of X args
	// KTSlice *x1		// (I) X1  arg timeslice
	// ..
	// KTSlice *xNx		// (I) Xnx arg timeslice
	...);





#endif /* _dvtspar_H */
