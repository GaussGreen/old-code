/************************************************************************
 * Module:	DRL
 * Submodule:	RAND
 * Function:	Random Number Generation and Simulation
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */
#include <math.h>

#include "drlmem.h"
#include "drlmath.h"
#include "drllineq.h"

#include "drlrand.h"

#define	DIM_MAX	20

static	double	**mAMatat,
		mSavDeviate[DIM_MAX];
static	int	mNDim,
		mAntiFlag=0;
static	int	mIsSobol;
static	long	mSeedSave;

/*#define	__DEBUG__*/

/*f---------------------------------------------------------------------
 * Initializes a multivariate normal simulation. Should be
 * used in conjunction with  <i> DrlMultiNormSimulGet</i> and
 * <i> DrlMultiNormSimulDispose</i>.
 * <br>
 * <br>[seed] if nonnegative, the seed of the random number generator.
 * If -1, initializes a Sobol sequence.
 * <br>[nDim] the number of dimensions of the simulation.
 * <br>[corrMat] the correlation matrix (size <i> nDim</i>\*<i> nDim</i>.
 * <br>
 * Returns 0 iff successful.
 */

int
DrlMultiNormSimulInit(long seed, int nDim, double **corrMat)
{
static	char	routine[] = "DrlMultiNormSimulInit";
	int	status = FAILURE;
	int	n;

	if ((nDim < 1) || (nDim >= DIM_MAX)) {
		DrlErrMsg("%s: bad dimension %d\n", routine, nDim);
		goto done;
	}

	mNDim = nDim;

	if (DrlMatrixCholesky(&mAMatat, corrMat, mNDim) != 0) {
		DrlErrMsg("%s: correlation matrix decomposition "
			"failed\n", routine);
		goto done;
	}

	mAntiFlag = 0;

	if (seed == -1L) {
	    if (nDim > 6) {
		DrlErrMsg("%s: dimension %d too high for Sobol.\n",
			routine, nDim);
		goto done;
	    }
	    mIsSobol = TRUE;
	    n = (-1);
	    DrlSobolSequence(&n, NULL);
	} else {
	    mIsSobol = FALSE;
	    mSeedSave = seed;
	    mSeedSave = 0L;
	}


	status = SUCCESS;
done:
	if (status != SUCCESS)
		DrlErrMsg("%s: Failed.\n", routine);
	return(status);
}

/*f---------------------------------------------------------------------
 * Returns, at each call, a new  multidimensional normal deviate.
 * To be called after <i> DrlMultiNormSimulInit</i>:
 * the array "deviate" should have length of "nDim".
 */

int
DrlMultiNormSimulGet(double *deviate)
{
static	char	routine[] = "DrlMultiNormSimulGet";
	int	i, j;
	double	x;

	if (mIsSobol != FALSE) {
		/* QMC */
		DrlSobolSequence(&mNDim, mSavDeviate);
		for (i=0; i<=mNDim-1; i++) {
			DrlCumNormInv(mSavDeviate[i], mSavDeviate+i);
		}
	} else {
		/* Regular MC */
		if (mAntiFlag == 0) {
			for (i=0; i<=mNDim-1; i++) {
				mSavDeviate[i] = DrlGaussSimul(&mSeedSave);
			}
			mAntiFlag = 1;
		} else {
			for (i=0; i<=mNDim-1; i++) {
				mSavDeviate[i] = (-mSavDeviate[i]);
			}
			mAntiFlag = 0;
		}
	}

#ifdef	__DEBUG__
	DrlErrMsg("DEVIATES(%3s):",
		(mIsSobol != FALSE ? "QMC" : "MC"));
	for (i=0; i<=mNDim-1; i++)
		DrlErrMsg(" %12.4f", mSavDeviate[i]);
	DrlErrMsg("\n");
#endif

	for (i=0; i<=mNDim-1; i++) {
		x = 0e0;
		for (j=0; j<=i; j++)
			x += mAMatat[i][j] * mSavDeviate[j];
		deviate[i] = x;
	}

	return(SUCCESS);
}

/*f---------------------------------------------------------------------
 * Frees the memory allocated in <i> DrlMultiNormSimulInit</i>.
 */

int
DrlMultiNormSimulDispose(void)
{
	DrlMatrixFree(mAMatat, mNDim, mNDim);
	mAMatat = NULL;

	return(SUCCESS);
}



