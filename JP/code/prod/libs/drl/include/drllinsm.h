/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	Linear algebra routine with static memory
 * File:	drllineq.h
 * Function:	header
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drllinsm_H
#define	_drllinsm_H

#include "drllineq.h"


int	 DrlMatrixCopySM(
	double ***toMat,	/* (O) new copy matrix */
	double **fromMat,	/* (I) input matrix */
	int nl,			/* (I) # lines */
	int nc);			/* (I) # columns */

int
DrlMatrixTransposeSM(
	double ***mat,		/* (I/O) matrix to transpose */
	int nl,			/* (I) # lines */
	int nc);		/* (I) # columns */

int
DrlMatrixProductSM(
	double ***prodMat,	/* (O) new product matrix */
	double **aMat,		/* (I) 1st matrix */
	int aNl,		/* (I) # of lines 1st matrix */
	int aNc,		/* (I) # of columns 1st matrix */
	double **bMat,		/* (I) 2nd matrix */
	int bNl,		/* (I) # of lines 2nd matrix */
	int bNc);		/* (I) # of columns 2nd matrix */

int
DrlMatrixSvdDecompSM(
	int m,			/* (I) # of lines */
	int n,			/* (I) # of columns */
	double **a,		/* (I) input matrix m x n [0..m-1][0..n-1] */
	double ***u, 		/* (O) output [0..m-1][0..n-1] */
	double **w, 		/* (O) output vap [0..n-1] */
	double ***v); 		/* (O) output [0..n-1][0..n-1] */

int
DrlMatrixPseudoInverseSM(
	int m,			       /* (I) # of lines */
	int n,			       /* (I) # of columns */
	double **u, 		   /* (I) from SVD [0..m-1][0..n-1] */
	double *lambda,		   /* (I) psudo inverse EV [0..n-1] */
	double **v, 		   /* (I) from SVD [0..n-1][0..n-1] */
	double ***outputMatrix /* Inverse matrix as output */
);

int
DrlRealLinSysSvdSM(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double alpha,		/* (I) see DrlMatrixSvdRegularizeEV */
	double param);		/* (I) see DrlMatrixSvdRegularizeEV */

int
DrlPseudoInverseSM(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double ***outputMatrix, /* (I) outputMatrix */ 
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double alpha,		/* (I) see DrlMatrixSvdRegularizeEV */
	double param);		/* (I) see DrlMatrixSvdRegularizeEV */

#endif	/*_drllinsm_H */
