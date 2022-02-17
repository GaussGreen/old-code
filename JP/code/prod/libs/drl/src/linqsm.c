/************************************************************************
 * Module:	DRL - LINEQ
 * Function:	Static memory version of some routines
 * Author:	
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <math.h>
#include <float.h>
#include "drlmem.h"

#include "drllineq.h"		/* prototype consistency */


int	drlSvdcmp(double **a, int m, int n, double w[], double **v);

/*CHANGES IN LINQMTR.C*/


/*f---------------------------------------------------------------------
 * Matrix copy.
 *
 * No allocation
 */

int
DrlMatrixCopySM(
	double ***toMat,	/* (O) new copy matrix */
	double **fromMat,	/* (I) input matrix */
	int nl,			/* (I) # lines */
	int nc)			/* (I) # columns */
{
static	char	routine[] = "DrlMatrixCopy";
	int	i, j;

	for(i=0; i<=nl-1; i++)
	for(j=0; j<=nc-1; j++) {
		(*toMat)[i][j] = fromMat[i][j];
	}

	return(0);
}



/*f---------------------------------------------------------------------
 * Matrix transpose.
 *                        
 * Max size reallocation
 */
int
DrlMatrixTransposeSM(
	double ***mat,		/* (I/O) matrix to transpose */
	int nl,			/* (I) # lines */
	int nc)			/* (I) # columns */
{
static	char	routine[] = "DrlMatrixTranspose";
	register int	i, j;
	double	**tMat;

    if ((tMat = DrlMatrixNew(MAX(nc,nl), MAX(nl,nc))) == NULL) {
		DrlErrMsg("%s: matrix allocation failed\n", routine);
		return(-4);
	}

	for(i=0; i<=nl-1; i++)
	for(j=0; j<=nc-1; j++) {
		tMat[j][i] = (*mat)[i][j];
	}

    DrlMatrixFree(*mat, MAX(nl,nc), MAX(nc,nl));
	*mat = tMat;

	return(0);
}

/*f---------------------------------------------------------------------
 * Matrix product.
 * No allocation
 *
 */

int
DrlMatrixProductSM(
	double ***prodMat,	/* (O) new product matrix */
	double **aMat,		/* (I) 1st matrix */
	int aNl,		/* (I) # of lines 1st matrix */
	int aNc,		/* (I) # of columns 1st matrix */
	double **bMat,		/* (I) 2nd matrix */
	int bNl,		/* (I) # of lines 2nd matrix */
	int bNc)		/* (I) # of columns 2nd matrix */
{
static	char	routine[] = "DrlMatrixProduct";
	int	i, j, k;
	double	x;


	if (aNc != bNl) {
		DrlErrMsg("%s: can't multiply (aNc = %d, bNl = %d)\n",
			routine, aNc, bNl);
		*prodMat = NULL;
		return(1);
	}

	for(i=0; i<=aNl-1; i++)
	for(j=0; j<=bNc-1; j++) {
		x = 0e0;
		for(k=0; k<=aNc-1; k++) {
			x += aMat[i][k]*bMat[k][j];
		}
		(*prodMat)[i][j] = x;
	}

	return(0);
}



/*END CHANGES IN LINQMTR.C*/

/*CHANGES IN LINQSVD*/

/*f---------------------------------------------------------------------
 * Matrix singular value decomposition.
 * No allocation
 */

int
DrlMatrixSvdDecompSM(
	int m,			/* (I) # of lines */
	int n,			/* (I) # of columns */
	double **a,		/* (I) input matrix m x n [0..m-1][0..n-1] */
	double ***u, 		/* (O) output [0..m-1][0..n-1] */
	double **w, 		/* (O) output vap [0..n-1] */
	double ***v) 		/* (O) output [0..n-1][0..n-1] */
{
static	char	routine[] = "DrlMatrixSvdDecomp";
	int	status = FAILURE;


	if (DrlMatrixCopy(u, a, m, n) != SUCCESS)
		goto done;

	if ((*w = DrlDoubleVectAlloc(0, n-1)) == NULL)
		goto done;


	DrlDoubleMatrCToNr(&a, 1, m, 1, n);
	DrlDoubleMatrCToNr(u, 1, m, 1, n);
	(*w) -= 1;
	DrlDoubleMatrCToNr(v, 1, n, 1, n);

	/*svdcmp(*u, m, n, *w, *v);*/
	drlSvdcmp(*u, m, n, *w, *v);

	DrlDoubleMatrNrToC(&a, 1, m, 1, n);
	DrlDoubleMatrNrToC(u, 1, m, 1, n);
	(*w) += 1;
	DrlDoubleMatrNrToC(v, 1, n, 1, n);


	status = SUCCESS;
	/* made it through OK */
done:
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed\n", routine);
	}
	return(status);
}

/*f---------------------------------------------------------------------
 * Matrix pseudo-inversion.
 *
 * Only inverses the matrix
 * No allocation
 *
 * Returns 0 iff successful.
 */

int
DrlMatrixPseudoInverseSM(
	int m,			       /* (I) # of lines */
	int n,			       /* (I) # of columns */
	double **u, 		   /* (I) from SVD [0..m-1][0..n-1] */
	double *lambda,		   /* (I) psudo inverse EV [0..n-1] */
	double **v, 		   /* (I) from SVD [0..n-1][0..n-1] */
	double ***outputMatrix /* Inverse matrix as output */
)		
{
	register int	i, j, k;
	
	/* compute X = V W^-1 tU b */
		for (k=0; k<=n-1; k++) 
		{
		    for (j=0; j<=m-1; j++)
			{
			(*outputMatrix)[k][j] = 0e0;	
			for (i=0; i<=n-1; i++)
			(*outputMatrix)[k][j] += v[k][i] * lambda[i] * u[j][i];
			}
		}

	return(SUCCESS);
}


/*f---------------------------------------------------------------------
 * Linear system solving with eigenvalue regularisation.
 *
 * <br><br>
 * Solves a possibly ill conditioned linear system
 * <blockquote>
 * sum <sub>j=1,...,n</sub> a<sub>ij</sub> x<sub>j</sub> = b<sub>i</sub>, for i=1,...,m,
 * </blockquote>
 * by singular value decomposition and regularization.
 * See <i> DrlMatrixSvdRegularizeEV</i> for a description
 * of the regularization parameters.
 */
/********************/
int
DrlRealLinSysSvdSM(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double alpha,		/* (I) see DrlMatrixSvdRegularizeEV */
	double param)		/* (I) see DrlMatrixSvdRegularizeEV */
{
static	char	routine[] = "DrlRealLinSysSvd";
	int	status = FAILURE;

	double	**u = NULL,	/* [0..m-1][0..n-1] */
		*w = NULL,	/* [0..n-1] */
		**v = NULL;	/* [0..n-1][0..n-1] */

	//CGE
	if ((v = DrlMatrixNew(n, n)) == NULL)
		goto done;

	if ((u = DrlMatrixNew(m, n)) == NULL)
		goto done;
	//CGE

	if (DrlMatrixSvdDecomp(m, n, a, &u, &w, &v) != SUCCESS)
		goto done;

	/* compute pseudo inverse of w */
	if (DrlMatrixSvdRegularizeEV(n, w, w, regType, alpha, param) != SUCCESS)
			goto done;

	/* compute X = V W^-1 tU b */
	if (DrlMatrixSvdPseudoInverse(m, n, u, w, v, b, x) != SUCCESS)
		goto done;

	status = SUCCESS;
	/* made it through OK */
done:
	DrlMatrixFree(u, m, n);
	DrlDoubleVectFree(w, 0, n-1);
	DrlMatrixFree(v, n, n);

	if (status != SUCCESS) {
		DrlErrMsg("%s: failed\n", routine);
	}
	return(status);
}

/*f---------------------------------------------------------------------
 * Pseudo-Inverse Matrix .
 * Do not reallocate the pseudo-inverse matrix.
 *
 */

int
DrlPseudoInverseSM(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double ***outputMatrix, /* (I) outputMatrix */ 
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double alpha,		/* (I) see DrlMatrixSvdRegularizeEV */
	double param)		/* (I) see DrlMatrixSvdRegularizeEV */
{
static	char	routine[] = "DrlPseudoInverse";
	int	status = FAILURE;

	double	**u = NULL,	/* [0..m-1][0..n-1] */
		*w = NULL,	/* [0..n-1] */
		**v = NULL;	/* [0..n-1][0..n-1] */


	if ((v = DrlMatrixNew(n, n)) == NULL)
		goto done;

	if ((u = DrlMatrixNew(m, n)) == NULL)
		goto done;

	if (DrlMatrixSvdDecomp(m, n, a, &u, &w, &v) != SUCCESS)
		goto done;

	/* compute pseudo inverse of w */
	if (DrlMatrixSvdRegularizeEV(n, w, w, regType, alpha, param) != SUCCESS)
			goto done;

	/* compute X = V W^-1 tU b */
	if (DrlMatrixPseudoInverseSM(m, n, u, w, v, outputMatrix) != SUCCESS)
		goto done;

	status = SUCCESS;
	/* made it through OK */
done:
	DrlMatrixFree(u, m, n);
	DrlDoubleVectFree(w, 0, n-1);
	DrlMatrixFree(v, n, n);

	if (status != SUCCESS) {
		DrlErrMsg("%s: failed\n", routine);
	}
	return(status);
}

/*END CHANGES IN LINQSVD*/
