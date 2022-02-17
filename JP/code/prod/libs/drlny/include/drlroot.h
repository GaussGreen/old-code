/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	ROOT - Root Finding
 * File:	drlroot.h
 * Function:	Root finding.
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlroot_H
#define	_drlroot_H
#include "drlstd.h"


/* constants used in optional argumsnts */

#define	DRL_ROOT_MAXITER			((int) 0x01)
#define	DRL_ROOT_FPLOG			((int) 0x02)
#define	DRL_ROOT_TOLF			((int) 0x03)
#define	DRL_ROOT_TOLX			((int) 0x04)
#define	DRL_ROOT_MEM			((int) 0x05)
#define	DRL_ROOT_NOERRMSG		((int) 0x06)


#define	DRL_ROOT_NEWDIC_DFX		((int) 0x21)

#define	DRL_ROOT_BROYDEN_JACOBIAN	((int) 0x31)


/* 1-D Brent */

extern	DLL_EXPORT(int)	DrlRootFinderBrent(
	double (*func)(double),	/* (I) function to solve */
	double x1,		/* (I) low bound */
	double x2,		/* (I) high bound */
	double tol,		/* (I) tolerance */
	int itmax,		/* (I) max iterations */
	double *xSol);		/* (I) solution */



DLL_EXPORT(int)	DrlRootFinderBrentRP(
	double *X,	/* (I) root estimate */
	double Y,	/* (I) function value: WARNING: ON FIRST CALL,
			 *     MUST SET TO f(XMIN) */
	double XMIN,	/* (I) low vbound */
	double XMAX,	/* (I) hghh bound */
	int *iter,	/* (I) iteration counter */
	int *indic,	/* (O) TRUE if root found */
	/* DRL_ROOT_TOLX, double tolx,		// tolerance 
	 * DRL_ROOT_MEM, double c[20],		// working memory 
	 * DRL_ROOT_MAXITER, int itermax,	// max # iterations
	 * 0) (last arg MUST be 0) 
	 */
	...);



/* 1-D Newton/Dicho */
extern	DLL_EXPORT(int)	DrlRootFinderNewtonDicho(
	double *x,	/* (I) current estimate of the root */
	double fx,
	double xmin,	/* (I) lower bound */
	double xmax,	/* (I) upper bound */
	int *count,
	int *ind,
	/* DRL_ROOT_TOLX, double tolx,		// tolerance 
	 * DRL_ROOT_MEM, double c[20],		// working memory 
	 * 0) (last arg MUST be 0) 
	 */
	...);



/* Multi-D Broyden */

/*t-@CDOC(idxn=TMDFunc,catn=structdef)
 * Type definition for a n-dimensional function f
 * of an n-dimensional vector x=(x_0,...,x_{n-1})
 * <blockquote>
 * x = (x_0,...,x_{n-1}) to
 *     f(x)=(f_0(x),...,f_{n-1}(x)).
 * </blockquote>
 * The routine return value should be SUCCESS/FAILURE.
 */
typedef	int (*TMDFunc)(
	int n,			/* (I) # dim */
	double *v,		/* (I) point where evaluated [0..n-1] */
	double *f);		/* (O) vector values [0..n-1] */
/*e*/


extern	DLL_EXPORT(int)	DrlRootFindBroyd(
	double *x,		/* (I/O) start/end point */
	int n,			/* (I) number of dim */
	int *check,		/* (I) */
	TMDFunc vecfunc,	/* (I) */
	/* DRL_ROOT_BROYDEN_JACOBIAN, double **jac,
	 * DRL_ROOT_MAXITER, int iterMax,
	 * DRL_ROOT_FPLOG, FILE *fpLog,
	 * DRL_ROOT_TOLF, double tolf,
	 * 0,
	 */
	...);


extern	DLL_EXPORT(int)	DrlRootNonlinearSystem(
	double *x,		/* (I/O) start/end point */
	int n,			/* (I) number of dim */
	TMDFunc vecfunc,	/* (I) user defined function to solve */
	/* DRL_ROOT_BROYDEN_JACOBIAN, double **jac,
	 * DRL_ROOT_MAXITER, int iterMax,
	 * DRL_ROOT_FPLOG, FILE *fpLog,
	 * 0,
	 */
	...);


extern DLL_EXPORT(int) DrlNRPoly(
	double  guess,		/* (I) First guess for the root       */
	double  *a,		/* (I) Coefficients of the polynomial */
	int n,			/* (I) Degree of the polynomial       */
	double *retVal);	/* (O) Root (when successful)         */

#endif	/* _drlroot_H */
