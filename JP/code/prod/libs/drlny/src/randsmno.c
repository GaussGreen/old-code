/************************************************************************
 * Module:	DRL
 * Submodule:	RAND
 * Function:	Random Number Generation and Simulation
 * Author:	C. Daher
 * Revision:	$Id: randsmno.c 55172 2006-03-27 20:33:17Z shuxiang.a.liu $
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include <math.h>

#include "drlmem.h"
#include "drllineq.h"

#ifdef CLIB
# include "gtobf.h"
# define CumNormInv	GtoNormalCumInv
#endif

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

DLL_EXPORT(int)
DrlMultiNormSimulInit(long seed, int nDim, double **corrMat)
{
static	char	routine[] = "DrlMultiNormSimulInit";
	int	status = FAILURE;
	int	n, i, j;

	if ((nDim < 1) || (nDim >= DIM_MAX)) {
		GtoErrMsg("%s: bad dimension %d\n", routine, nDim);
		goto done;
	}

        /* test the symmetry of corrMat */
        for (i = 0; i < nDim; i++) {
            if (!ARE_ALMOST_EQUAL(corrMat[i][i], 1.0)) {
                GtoErrMsg("%s: Bad correlation matrix.\n", routine);
                goto done;
            }
            for (j = i+1; j < nDim; j++)
                if (!ARE_ALMOST_EQUAL(corrMat[i][j],corrMat[j][i])) {
                    GtoErrMsg("%s: Bad correlation matrix.\n", routine);
                    goto done;
                }
        }

	mNDim = nDim;

        /* positive definite of corrMat is tested inside */
	if (DrlMatrixCholesky(&mAMatat, corrMat, mNDim) != 0) {
		GtoErrMsg("%s: correlation matrix decomposition "
			"failed\n", routine);
		goto done;
	}

	mAntiFlag = 0;

	if (seed == -1L) {
	    if (nDim > 6) {
		GtoErrMsg("%s: dimension %d too high for Sobol.\n",
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
		GtoErrMsg("%s: Failed.\n", routine);
	return(status);
}

/*f---------------------------------------------------------------------
 * Returns, at each call, a new  multidimensional normal deviate.
 * To be called after <i> DrlMultiNormSimulInit</i>:
 * the array "deviate" should have length of "nDim".
 */

DLL_EXPORT(int)
DrlMultiNormSimulGet(double *deviate)
{
static	char	routine[] = "DrlMultiNormSimulGet";
	int	i, j;
	double	x;

	if (mIsSobol != FALSE) {
		/* QMC */
		DrlSobolSequence(&mNDim, mSavDeviate);
		for (i=0; i<=mNDim-1; i++) {
			CumNormInv(mSavDeviate[i], mSavDeviate+i);
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
	GtoErrMsg("DEVIATES(%3s):",
		(mIsSobol != FALSE ? "QMC" : "MC"));
	for (i=0; i<=mNDim-1; i++)
		GtoErrMsg(" %12.4f", mSavDeviate[i]);
	GtoErrMsg("\n");
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

DLL_EXPORT(int)
DrlMultiNormSimulDispose(void)
{
	DrlMatrixFree(mAMatat, mNDim, mNDim);
	mAMatat = NULL;

	return(SUCCESS);
}



