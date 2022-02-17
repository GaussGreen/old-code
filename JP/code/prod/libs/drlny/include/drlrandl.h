/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	RAND - Simulation Routines
 * File:	drlrandl.h
 * Function:	Random number generators.
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlrandl_H
#define	_drlrandl_H
#include "drlstd.h"

/* Monte-Carlo on a parsed formula (The Parser) - Wrapper */


extern	DLL_EXPORT(int)	DrlFormMonteCarloL(
	long *nDimL,		/* 0 'L' (I) number of dimensions */
	char *whatL,		/* 1 'C' (I) 'n' normal, 'l' log-normal */
	double *pL,		/* 2 'F' (I) expectation vector[0..nDim-1] */
	double *texpL,		/* 3 'F' (I) time factor */
	double *volL,		/* 4 'F' (I) stddev vector[0..nDim-1] */
	double *corrMatL,	/* 5 'F' (I) corr matrix[0..nDim-1][0..nDim-1]*/
	char *formL,		/* 6 'C' (I) formula */
	long *nSampleL,		/* 7 'L' (I) number of deviates */
	long *doStatFlagL,	/* 8 'L' (I) if TRUE, performs statistics */
	long *traceL,           /* 9 'L' (I) if TRUE, output to log file */
	double *retValL);	/*   'F' (O) */


#endif	/* _drlrandl_H */


