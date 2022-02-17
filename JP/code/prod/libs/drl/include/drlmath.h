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

extern	double	DrlDenNorm(double x);
extern	double	DrlCumNorm(double x);
extern	double	DrlCumNormInv(double p, double *x);
extern	double	DrlCumBiNorm(double a, double b, double p);
extern	double	DrlCumMultiNorm(int dim, double *a, double *r);

/* Generic numerical procedures - JC */

extern	double	DrlSimpsIntegral(double step, long numPoints,
						 double *fValues);       
extern	double	DrlSmoothStepFcn(double x, double step);
extern	double	DrlSmoothMAX(double x, double step);

#endif	/* _drlmath_H */
