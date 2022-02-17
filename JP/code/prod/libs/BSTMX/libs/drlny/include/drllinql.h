/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	LINEQ - Wrapper
 * File:	drllinql.h
 * Function:	header
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drllinql_H
#define	_drllinql_H

#include "drlstd.h"

extern	DLL_EXPORT(int)	DrlMatrixRealEigenVectL(
	double *aL,		/*  1 F (I) matrix */
	double *vapL,		/*  2 F (O) array of eigen values */
	double *vepL);		/*  3 F (O) array of eigen vectors */

extern	DLL_EXPORT(int)	DrlMatrixInvertL(
	double *aL,		/*  1 F (I) input matrix */
	double *ainvL);		/*  2 F (O) output inverse matrix */


DLL_EXPORT(int)		DrlVectSortL(
	double *invecL,		/*  1 F (I) input vector */
	double *outvecL);	/*  2 F (O) output vector */


#endif	/* _drllinql_H */
