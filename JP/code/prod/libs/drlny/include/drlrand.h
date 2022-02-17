/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	RAND - Simulation Routines
 * File:	drlrand.h
 * Function:	Random number generators.
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlrand_H
#define	_drlrand_H
#include "drlstd.h"

#include <stdlib.h>


/* Returns a random number between a and b */

extern	DLL_EXPORT(double)	DrlDoubleRand(long *idum,
					double xLow, double xHigh);
extern	DLL_EXPORT(double)	DrlGaussSimul(long *idum);
extern	DLL_EXPORT(int)		DrlIntRand(long *idum, int iLow, int iHigh);


/* Generation of multidimensional normal deviates */

extern	DLL_EXPORT(int)	DrlMultiNormSimulInit(long seed, int nDim,
				double **corrMat);
extern	DLL_EXPORT(int)	DrlMultiNormSimulGet(double *deviate);
extern	DLL_EXPORT(int)	DrlMultiNormSimulDispose(void);

/* Monte-Carlo on a parsed formula (The Parser) */

extern	DLL_EXPORT(int)	DrlFormMonteCarlo(
	int	nDim,	 	/* (I) number of dimensions */
	char	what,		/* (I) 'n' form normal, 'l' for log-normal */
	char	mcType,		/* (I) 'm' for MC, 's' for Sobol */
	double	*p,		/* (I) expectation vector[0..nDim-1] */
	double	texp,		/* (I) time factor */
	double	*vol,		/* (I) stddev vector[0..nDim-1] */
	double	**corrMat,	/* (I) corr matrix[0..nDim-1][0..nDim-1]*/
	char	*form,		/* (I) formula */
	long	nSample,	/* (I) number of deviates */
	int	doStatFlag,	/* (I) if TRUE, performs statistics */
	double	*retVal);	/* (O) */


/* Quasi Monte Carlo */

extern	DLL_EXPORT(int)	DrlSobolSequence(int *n, double *x);


#endif	/* _drlrand_H */


