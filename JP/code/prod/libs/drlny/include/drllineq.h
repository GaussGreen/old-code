/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	LINEQ - System of Linear Equations - Matrix Algebra Routines
 * File:	drllineq.h
 * Function:	header
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drllineq_H
#define	_drllineq_H

#include "drlstd.h"
#include <stdio.h>


/* Vector operations */
extern	DLL_EXPORT(int)		DrlVectScalarOper(double *v, int n,
					char *operation, double scalarArg);
extern	DLL_EXPORT(int)		DrlVectUnaryOper(double *v, int n,
					char *operation, double *vectArg);

extern	DLL_EXPORT(int)		DrlVectNorm(double *v, int nv, char *normType,
					double *retVal);


/* Matrix Algebra */

extern	DLL_EXPORT(double**)	DrlMatrixNew(int nl, int nc);
extern	DLL_EXPORT(int)		DrlMatrixFree(double **p, int nl, int nc);
extern	DLL_EXPORT(double**)	DrlMatrixNewDiag(int nl, double *w);
extern	DLL_EXPORT(int)		DrlMatrixCopy(double ***toMat,
					double **fromMat,
					int nl, int nc);
extern	DLL_EXPORT(int)		DrlMatrixPrint(FILE *fp, double **p,
					int nl, int nc);
extern	DLL_EXPORT(int)		DrlMatrixUnaryMatrixOper(
					double **aMat, int aNl, int aNc,
					double **bMat, int bNl, int bNc,
					char *oper);
extern	DLL_EXPORT(int)		DrlMatrixTranspose(double ***mat, int nl,
					int nc);
extern	DLL_EXPORT(int)		DrlMatrixProduct(double ***prodMat,
					double **aMat, int aNl, int aNc,
					double **bMat, int bNl, int bNc);
extern	DLL_EXPORT(int)		DrlMatrixInvert(double ***invMat,
					double **aMat,
					int nl);

/* Matrix norms */

extern	DLL_EXPORT(double)	DrlMatrixNorm(double **mat, int nl, int nc,
					int p);
extern	DLL_EXPORT(double)	DrlMatrixL2Product(double **aMat,
					double **bMat,
					int nl, int nc);
extern	DLL_EXPORT(int)		DrlMatrixCond(double **mat, int nl, int p,
					double *cond);

/* Decomposition routines */

extern	DLL_EXPORT(int)		DrlMatrixCholesky(double ***choMat,
					double **aMat, int nl);
extern	DLL_EXPORT(int)		DrlMatrixCheckPositiveDefinite(double **aMat,
					int nl);

/* Eigen value problems */

extern	DLL_EXPORT(int)		DrlMatrixRealEigenVect(int n, double **a,
					double *vap, double **vep);

/* Misc */

extern	DLL_EXPORT(int)		DrlMatrixSquareAnalyze(FILE *fp,
					double **p, int nl);


/*  System of real linear equations */

extern	DLL_EXPORT(int)		DrlRealLinSys(
	int n,			/* (I) # of dim */
	double **a,		/* (I) system coeffs [0..n-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b);		/* (I) inhomogeneous term [0..n-1] */

extern	DLL_EXPORT(int)		DrlRealLinSysCheck(
	int m,			/* (I) # of equations */
	int n,			/* (I) # of unknowns */
	double **a,		/* (I) system coeffs [0..m-1][0..n-1] */
	double *x,		/* (I) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	double *axmb);		/* (O) A*X-B [0..m-1] */

extern	DLL_EXPORT(int)		DrlRealLinSysPivot(
	int m,			/* (I) # of equations */
	int n,			/* (I) # of unknowns */
	double **a,		/* (I) system coeffs [0..m-1][0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int npiv,		/* (I) # of pivots */
	int *ipiv,		/* (I) array of line to pivot [0..npiv] */
	int *jpiva,		/* (I) array col to pivot (NULL if implicit) */
	double ***apo,		/* (I) piv system [0..m-1-npiv][0..n-1-npiv] */
	double **bpo,		/* (I) piv inhomogeneous term [0..m-1-npiv] */
	double ***cpo,		/* (I) sol transf matr [0..n-1][0..n-1-npiv] */
	double **dpo);		/* (I) sol transf 2ndt [0..n-1] */

/* Singular Value Decomposition */

extern	DLL_EXPORT(int)	DrlRealLinSysLog(
	int m,			/* (I) # of equations */
	int n,			/* (I) # of unknowns */
	double **a,		/* (I) system coeffs [0..m-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	FILE *fp);		/* (I) */

extern	DLL_EXPORT(int)	DrlMatrixSvdDecomp(
	int m,			/* (I) # of lines */
	int n,			/* (I) # of columns */
	double **a,		/* (I) input matrix m x n [0..m-1][0..n-1] */
	double ***u, 		/* (O) output [0..m-1][0..n-1] */
	double **w, 		/* (O) output vap [0..n-1] */
	double ***v); 		/* (O) output [0..n-1][0..n-1] */


extern	DLL_EXPORT(int)	DrlMatrixSvdPseudoInverse(
	int m,			/* (I) # of lines */
	int n,			/* (I) # of columns */
	double **u, 		/* (I) from SVD [0..m-1][0..n-1] */
	double *lambda,		/* (I) psudo inverse EV [0..n-1] */
	double **v, 		/* (I) from SVD [0..n-1][0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	double *x);		/* (O) solution [0..n-1] */


extern	DLL_EXPORT(int)	DrlMatrixSvdRegularizeEV(
	int n,			/* (I) # of columns */
	double *w, 		/* (I) EV [0..n-1] */
	double *lambda,		/* (O) pseudo inverse EV [0..n-1] */
	int regType,		/* (I) see below */
	double alpha,		/* (I) see below */
	double param);		/* (I) see below */

extern	DLL_EXPORT(int)	DrlRealLinSysSvd(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int regType,		/* (I) see below */
	double alpha,		/* (I) see below */
	double param);		/* (I) see below */


/* Singular systems */

extern	DLL_EXPORT(int)	DrlRealLinSysSvdSelectNz(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int numNz,		/* (I) max non zero terms required */
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double alpha,		/* (I) see DrlMatrixSvdRegularizeEV */
	double param);		/* (I) see DrlMatrixSvdRegularizeEV */

extern	DLL_EXPORT(int)	DrlRealLinSysSvdPivot(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *xo,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int npiv,		/* (I) # of pivots */
	int *ipiv,		/* (I) array of line to pivot [0..npiv] */
	int *jpiva,		/* (I) array of col to pivot [0..npiv] */
	double *axmb,		/* (O) residual ax-b [0..m-1] (or NULL) */
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double alpha,		/* (I) see DrlMatrixSvdRegularizeEV */
	double regParam);	/* (I) see DrlMatrixSvdRegularizeEV */

extern	DLL_EXPORT(int)	DrlRealLinSysSvdSelectNzPivot(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *xo,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int npiv,		/* (I) # of pivots */
	int *ipiv,		/* (I) array of line to pivot [0..npiv] */
	int *jpiv,		/* (I) array col to pivot (NULL if implicit) */
	double *axmb,		/* (O) residual ax-b [0..m-1] (or NULL) */
	int numNz,		/* (I) max non zero terms required */
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double alpha,		/* (I) see DrlMatrixSvdRegularizeEV */
	double param);		/* (I) see DrlMatrixSvdRegularizeEV */


extern	DLL_EXPORT(int)	DrlRealLinSysSvdSelectNzPivotOptim(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *xo,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int npiv,		/* (I) # of pivots */
	int *ipiv,		/* (I) array of line to pivot [0..npiv] */
	int *jpiv,		/* (I) array col to pivot (NULL if implicit) */
	double *axmb,		/* (O) residual ax-b [0..m-1] (or NULL) */
	int numNzMax,		/* (I) max non zero terms required */
	double rLinfMax,	/* (I) max residual linf norm */
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double regParam,	/* (I) see DrlMatrixSvdRegularizeEV */
	double *alphaO);	/* (O) optimal */


/* utilities */

extern	DLL_EXPORT(int)	drlInitPermutation(int numNz, int *idx, int n);
extern	DLL_EXPORT(int)	drlAdvancePermutation(int numNz, int *idx, int n);
extern	DLL_EXPORT(int)	drlComplPermutation(int numNz, int *idx, int n,
				int *idxComp);





#endif	/* _drllineq_H */
