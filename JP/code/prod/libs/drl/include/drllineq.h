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
extern	int		DrlVectScalarOper(double *v, int n,
					char *operation, double scalarArg);
extern	int		DrlVectUnaryOper(double *v, int n,
					char *operation, double *vectArg);

extern	int		DrlVectNorm(double *v, int nv, char *normType,
					double *retVal);


/* Matrix Algebra */

extern	double**	DrlMatrixNew(int nl, int nc);
extern	int		DrlMatrixFree(double **p, int nl, int nc);
extern	double**	DrlMatrixNewDiag(int nl, double *w);
extern	int		DrlMatrixCopy(double ***toMat,
					double **fromMat,
					int nl, int nc);
extern	int		DrlMatrixPrint(FILE *fp, double **p,
					int nl, int nc);
extern	int		DrlMatrixUnaryMatrixOper(
					double **aMat, int aNl, int aNc,
					double **bMat, int bNl, int bNc,
					char *oper);
extern	int		DrlMatrixTranspose(double ***mat, int nl,
					int nc);
extern	int		DrlMatrixProduct(double ***prodMat,
					double **aMat, int aNl, int aNc,
					double **bMat, int bNl, int bNc);
extern	int		DrlMatrixInvert(double ***invMat,
					double **aMat,
					int nl);

/* Matrix norms */

extern	double	DrlMatrixNorm(double **mat, int nl, int nc,
					int p);
extern	double	DrlMatrixL2Product(double **aMat,
					double **bMat,
					int nl, int nc);
extern	int		DrlMatrixCond(double **mat, int nl, int p,
					double *cond);

/* Decomposition routines */

extern	int		DrlMatrixCholesky(double ***choMat,
					double **aMat, int nl);
extern	int		DrlMatrixCheckPositiveDefinite(double **aMat,
					int nl);

/* Eigen value problems */

extern	int		DrlMatrixRealEigenVect(int n, double **a,
					double *vap, double **vep);

/* Misc */

extern	int		DrlMatrixSquareAnalyze(FILE *fp,
					double **p, int nl);


/*  System of real linear equations */

extern	int		DrlRealLinSys(
	int n,			/* (I) # of dim */
	double **a,		/* (I) system coeffs [0..n-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b);		/* (I) inhomogeneous term [0..n-1] */

extern	int		DrlRealLinSysCheck(
	int m,			/* (I) # of equations */
	int n,			/* (I) # of unknowns */
	double **a,		/* (I) system coeffs [0..m-1][0..n-1] */
	double *x,		/* (I) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	double *axmb);		/* (O) A*X-B [0..m-1] */

extern	int		DrlRealLinSysPivot(
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

extern	int	DrlRealLinSysLog(
	int m,			/* (I) # of equations */
	int n,			/* (I) # of unknowns */
	double **a,		/* (I) system coeffs [0..m-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	FILE *fp);		/* (I) */

extern	int	DrlMatrixSvdDecomp(
	int m,			/* (I) # of lines */
	int n,			/* (I) # of columns */
	double **a,		/* (I) input matrix m x n [0..m-1][0..n-1] */
	double ***u, 		/* (O) output [0..m-1][0..n-1] */
	double **w, 		/* (O) output vap [0..n-1] */
	double ***v); 		/* (O) output [0..n-1][0..n-1] */


extern	int	DrlMatrixSvdPseudoInverse(
	int m,			/* (I) # of lines */
	int n,			/* (I) # of columns */
	double **u, 		/* (I) from SVD [0..m-1][0..n-1] */
	double *lambda,		/* (I) psudo inverse EV [0..n-1] */
	double **v, 		/* (I) from SVD [0..n-1][0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	double *x);		/* (O) solution [0..n-1] */


extern	int	DrlMatrixSvdRegularizeEV(
	int n,			/* (I) # of columns */
	double *w, 		/* (I) EV [0..n-1] */
	double *lambda,		/* (O) pseudo inverse EV [0..n-1] */
	int regType,		/* (I) see below */
	double alpha,		/* (I) see below */
	double param);		/* (I) see below */

extern	int	DrlRealLinSysSvd(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int regType,		/* (I) see below */
	double alpha,		/* (I) see below */
	double param);		/* (I) see below */


/* Singular systems */

extern	int	DrlRealLinSysSvdSelectNz(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int numNz,		/* (I) max non zero terms required */
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double alpha,		/* (I) see DrlMatrixSvdRegularizeEV */
	double param);		/* (I) see DrlMatrixSvdRegularizeEV */

extern	int	DrlRealLinSysSvdPivot(
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

extern	int	DrlRealLinSysSvdSelectNzPivot(
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


extern	int	DrlRealLinSysSvdSelectNzPivotOptim(
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

extern	int	drlInitPermutation(int numNz, int *idx, int n);
extern	int	drlAdvancePermutation(int numNz, int *idx, int n);
extern	int	drlComplPermutation(int numNz, int *idx, int n,
				int *idxComp);


/* eigen systems */


extern	int	DrlMatrixRealSymEigen(
	int n,		/* (I) dimension */
	double **a,	/* (I) matrix */
	double *vap,	/* (O) array of eigen values */
	double **vep);	/* (O) array of eigen vectors */



#endif	/* _drllineq_H */
