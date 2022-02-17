/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	MATH - Mathematical Functions
 * File:	drlmath.h
 * Function:	header
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlmath_H
#define	_drlmath_H
#include "drlstd.h"


/* Useful Gaussian Integrals */

extern	DLL_EXPORT(double)	DrlDenNorm(double x);
extern	DLL_EXPORT(double)	DrlCumNorm(double x);
extern	DLL_EXPORT(double)	DrlCumNormInv(double p, double *x);
extern	DLL_EXPORT(double)	DrlCumBiNorm(double a, double b, double p);
extern	DLL_EXPORT(double)	DrlCumMultiNorm(int dim, double *a, double *r);

/* Generic numerical procedures - JC */

extern	DLL_EXPORT(double)	DrlSimpsIntegral(double step, long numPoints,
						 double *fValues);       
extern	DLL_EXPORT(double)	DrlSmoothStepFcn(double x, double step);
extern	DLL_EXPORT(double)	DrlSmoothMAX(double x, double step);

#endif	/* _drlmath_H */
