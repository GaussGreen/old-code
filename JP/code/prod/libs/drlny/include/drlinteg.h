/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	INTEG - Integration of functions
 * File:	drlinteg.h
 * Function:	header
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlinteg_H
#define	_drlinteg_H

#include "drlstd.h"


extern	DLL_EXPORT(int)	DrlQSimp(
	double (*func)(double),		/* (I) function to integrate */
	double a,			/* (I) lower bound */
	double b,			/* (I) upper bound */
	double *retVal);		/* (I) integral */


#endif	/* _drlinteg_H */
