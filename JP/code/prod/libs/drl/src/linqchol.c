/************************************************************************
 * Module:	DRL - LINEQ
 * Function:	Linear Equations and Matrix Algebra
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <math.h>

#include "drllineq.h"		/* prototype consistency */


/*f---------------------------------------------------------------------
 * Matrix Cholesky decomposition.
 *
 * <br><br>
 * Performs the Cholesky decomposition of the real symmetric positive definite
 * matrix <i> aMat</i> of size $nx n$.
 * On exit, <i> choMat</i> is a pointer
 * to an new allocated <i> double**</i> matrix
 * containing the decomposition. <br>
 * Returns 0 iff OK.
 */


int
DrlMatrixCholesky(
	double ***choMat,	/* (O) Cholesky decomposition */
	double **aMat,		/* (I) input square matrix */
	int nl)			/* (I) # lines */
{
static	char	routine[] = "DrlMatrixCholesky";
	int	i, j, k, errCode=0;
	double	sum;

	*choMat = NULL;

	if ((*choMat = DrlMatrixNew(nl, nl)) == NULL) {
		DrlErrMsg("%s: matrix allocation failed.\n", routine);
		return(-4);
	}

	for (i=0;i<=nl-1;i++) {
	    for (j=i;j<=nl-1;j++) {
		for (sum=aMat[i][j],k=i-1;k>=0;k--) {
			sum -= aMat[i][k]*aMat[j][k];
		}
		if (i == j) {
		    if (sum <= 0.0) {
			DrlErrMsg("%s: Cholesky decomposition failed "
				"(i=%d, sum=%g).\n",
				routine, i, sum);
			DrlMatrixFree(*choMat, nl, nl);
			(*choMat) = NULL;
			errCode = 1;
			goto done;
		    }
		    (*choMat)[i][i] = sqrt(sum);
		} else {
		    aMat[j][i] = sum / ((*choMat)[i][i]);
		}
	    }
	}

	/* lower triangle of aMat is the Cholesky decomposition */
	for (i=0;i<=nl-1;i++)
	for (j=i+1;j<=nl-1;j++) {
		(*choMat)[j][i] = aMat[j][i];
		(*choMat)[i][j] = 0e0;
	}

done:
	/* restore aMat in its original form */
	for (i=0;i<=nl-1;i++)
	for (j=i+1;j<=nl-1;j++)
		aMat[j][i] = aMat[i][j];


	return(errCode);
}




/*f---------------------------------------------------------------------
 * Matrix positive definite test.
 *
 * <br><br>
 * Checks if a  matrix <i> aMat</i> of size $nx n$ is positive definite.
 * Returns SUCCESS if it is, FAILURE otherwise.
 */


int
DrlMatrixCheckPositiveDefinite(
	double **aMat,		/* (I) input square matrix */
	int nl)			/* (I) # lines (and columns) */
{
static	char	routine[] = "DrlMatrixCheckPositiveDefinite";
	int	i, j, k,
		errCode=FAILURE;
	double	sum,
		**choMat = NULL;

	if ((choMat = DrlMatrixNew(nl, nl)) == NULL) {
		DrlErrMsg("%s: matrix allocation failed\n", routine);
		return(FAILURE);
	}

	for (i=0;i<=nl-1;i++) {
		for (j=i;j<=nl-1;j++) {
			for (sum=aMat[i][j],k=i-1;k>=0;k--) {
				sum -= aMat[i][k]*aMat[j][k];
			}
			if (i == j) {
				if (sum <= 0.0) {
					/* Cholesky decomposition failed */
					errCode = SUCCESS;
					goto done;
				}
				choMat[i][i] = sqrt(sum);
			} else {
				aMat[j][i] = sum / (choMat[i][i]);
			}
		}
	}
	/* lower triangle of aMat is the Cholesky decomposition */
	for (i=0;i<=nl-1;i++)
	for (j=i+1;j<=nl-1;j++) {
		choMat[j][i] = aMat[j][i];
		choMat[i][j] = 0e0;
	}

	/* made it through OK */
	errCode = SUCCESS;
done:
	if (choMat != NULL) DrlMatrixFree(choMat, nl, nl);
	/* restore aMat in its original form */
	for (i=0;i<=nl-1;i++)
	for (j=i+1;j<=nl-1;j++)
		aMat[j][i] = aMat[i][j];

	return(errCode);
}







