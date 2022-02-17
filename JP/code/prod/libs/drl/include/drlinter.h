/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	INTER - Interpolation and Extrapolation
 * File:	drlinter.h
 * Function:	header
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef _drlinter_H
#define _drlinter_H

#include "drlstd.h"

/*  Linear interpolation */
extern	int	DrlLinearInterp1d(double *xa, double *ya, int m,
				double x, double *y);
extern	int	DrlLinearInterp1dWeights(
	double *xa,	/* (I) array of X(i) (i=0,..,m-1) */
	int m,		/* (I) arrays size */
	double x,	/* (I) point to intepolate */
	int *jlo,	/* (O) lo index */
	int *jhi,	/* (O) hi index */
	double *wlo,	/* (O) lo weight */
	double *whi);	/* (O) hi weight */

extern	int	DrlLinearInterp2d(double *x1a, double *x2a, double **ya,
				int n1, int n2, double x1,
				double x2, double *y);
extern	int	DrlStepInterp1d(
	double *xa,	/* (I) array of X(i) (i=0,..,n-1) */
	double *ya,	/* (I) array of Y(i) (i=0,..,n-1) */
	int n,		/* (I) # of points */
	double x,	/* (I) point to interpolate */
	double *y);	/* (O) interpolated value */

/*  Date interpolation */
extern	int	DrlDDateLinearInterp1d(
	DDate *xa,	/* (I) array of X(i) (i=0,..,m-1) */
	double *ya,	/* (I) array of Y(i) (i=0,..,m-1) */
	int m,		/* (I) arrays size */
	DDate x,	/* (I) point to intepolate */
	double *y);	/* (O) interpolated value */

extern	int	DrlDDateStepInterp1d(DDate *xa, double *ya, int n,
				DDate x, double *y);


/* Spline interpolation */
extern	int	DrlSplineInterp1dInit(
	double *x,		/* (I) array of X(i) (i=0,..,n-1) */
	double *y,		/* (I) array of X(i) (i=0,..,n-1) */
	int n,			/* (I) # of points */
	double yp1,		/* (I) dx/dy at X(0) */
	double ypn,		/* (I) dx/dy at X(n-1) */
	double *y2);		/* (O) should be allocated on entry */
extern	int	DrlSplineInterp1dInterp(
	double *xa,		/* (I) array of X(i) (i=0,..,n-1) */
	double *ya,		/* (I) array of X(i) (i=0,..,n-1) */
	double *y2a,		/* (I) from SplineInit */
	int n,			/* (I) # of points */
	double x,		/* (I) point to intepolate */
	double *y);		/* (O) interpolated value */

extern	int	DrlSplineInterp1d(
	double *xa,		/* (I) array of X(i) (i=0,..,n-1) */
	double *ya,		/* (I) array of X(i) (i=0,..,n-1) */
	int n,			/* (I) # of points */
	double yp1,		/* (I) dx/dy at X(0) */
	double ypn,		/* (I) dx/dy at X(n-1) */
	double x,		/* (I) point to intepolate */
	double *y);		/* (O) interpolated value */



#endif	/*_drlinter_H*/
