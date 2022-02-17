/************************************************************************
 * Module:	DRL - LINQ
 * Function:	Wrapper
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <math.h>
#include "drlmem.h"
#include "drlsort.h"		/* DoubleVectSort() */

#include "drllineq.h"


/*----------------------------------------------------------------------
 * 
 */

static	double**
drlMatrixGetL(double *matL, int nr, int nc)
{
static	char	routine[] = "drlMatrixGetL";
	int	status = FAILURE;


	/* made ith through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s:failed.\n", routine);
	}
	return(status);
}


/*----------------------------------------------------------------------
 * Wrapper for <i> DrlMatrixRealEigenVect</i>.
 */

DLL_EXPORT(int)
DrlMatrixRealEigenVectL(
	double *aL,	/*  1 F (I) input matrix */
	double *vapL,	/*  2 F (O) array of eigen values */
	double *vepL)	/*  3 F (O) array of eigen vectors */
{
static	char	routine[] = "DrlMatrixRealEigenVectL";
	int	status = FAILURE;
	int	i, j, n;
	double	**a = NULL,
		*vap = NULL,
		**vep = NULL;

	WRAP_CHECK_VECTOR(aL);
	WRAP_CHECK_VECTOR(vapL);
	WRAP_CHECK_VECTOR(vepL);

	n = (int) sqrt((double) ARGSIZE(aL));
	ASSERT_OR_DONE(ARGSIZE(aL) == (unsigned) (n*n));
	ASSERT_OR_DONE(ARGSIZE(vapL) == (unsigned) n);
	ASSERT_OR_DONE(ARGSIZE(vepL) == (unsigned) (n*n));

	if ((a = DrlDoubleMatrAlloc(0, n-1, 0, n-1)) == NULL)
		goto done;
	if ((vep = DrlDoubleMatrAlloc(0, n-1, 0, n-1)) == NULL)
		goto done;
	if ((vap = DrlDoubleVectAlloc(0, n-1)) == NULL)
		goto done;


	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++)
		a[i][j] = aL[WRAP_MATR_IDX(n, n, i, j)];


	if (DrlMatrixRealEigenVect(n, a, vap, vep) != SUCCESS)
		goto done;


	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++)
		vepL[WRAP_MATR_IDX(n, n, i, j)] = vep[i][j];
	for (i=0; i<=n-1; i++)
		vapL[i+1] = vap[i];


	/* made ith through */
	status = SUCCESS;
done:
	DrlDoubleMatrFree(a, 0, n-1, 0, n-1);
	DrlDoubleMatrFree(vep, 0, n-1, 0, n-1);
	DrlDoubleVectFree(vap, 0, n-1);

	if (status != SUCCESS) {
		GtoErrMsg("%s:failed.\n", routine);
	}
	return(status);
}



/*----------------------------------------------------------------------
 * Wrapper for <i> DrlMatrixInvert</i>.
 */

DLL_EXPORT(int)
DrlMatrixInvertL(
	double *aL,	/*  1 F (I) input matrix */
	double *ainvL)	/*  2 F (O) output inverse matrix */
{
static	char	routine[] = "DrlMatrixInvertL";
	int	status = FAILURE;
	int	i, j, n;
	double	**a = NULL,
		**ainv = NULL;

	WRAP_CHECK_VECTOR(aL);
	WRAP_CHECK_VECTOR(ainvL);

	n = (int) sqrt((double) ARGSIZE(aL));
	ASSERT_OR_DONE(ARGSIZE(aL) == (unsigned) (n*n));
	ASSERT_OR_DONE(ARGSIZE(ainvL) == (unsigned) (n*n));

	if ((a = DrlDoubleMatrAlloc(0, n-1, 0, n-1)) == NULL)
		goto done;

	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++)
		a[i][j] = aL[WRAP_MATR_IDX(n, n, i, j)];


	if (DrlMatrixInvert(&ainv, a, n) != SUCCESS)
		goto done;


	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++)
		ainvL[WRAP_MATR_IDX(n, n, i, j)] = ainv[i][j];


	/* made ith through */
	status = SUCCESS;
done:
	DrlDoubleMatrFree(a, 0, n-1, 0, n-1);
	DrlDoubleMatrFree(ainv, 0, n-1, 0, n-1);

	if (status != SUCCESS) {
		GtoErrMsg("%s:failed.\n", routine);
	}
	return(status);
}




/*----------------------------------------------------------------------
 * Wrapper for <i> DrlDoubleVectSort</i>.
 */

DLL_EXPORT(int)
DrlVectSortL(
	double *invecL,		/*  1 F (I) input vector */
	double *outvecL)	/*  2 F (O) output vector */
{
static	char	routine[] = "DrlVectSortL";
	int	status = FAILURE;
	int	n, i;
	double	x;

	WRAP_CHECK_VECTOR(invecL);
	WRAP_CHECK_VECTOR(outvecL);

	n = (int) ARGSIZE(invecL);
	ASSERT_OR_DONE(ARGSIZE(invecL) == ARGSIZE(outvecL));
	for (i=0; i<=n-1; i++) outvecL[i+1] = invecL[i+1];

	if (DrlDoubleVectSort(&outvecL[1], n) != SUCCESS)
		goto done;

	if (DrlDoubleVectRevert(&outvecL[1], n) != SUCCESS)
		goto done;


	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s:failed.\n", routine);
	}
	return(status);
}

