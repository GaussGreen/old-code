/************************************************************************
 * Module:	DRL - LINEQ
 * Function:	Linear Equations and Matrix Algebra
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <math.h>
#include <stdio.h>
#include "drlmem.h"
#include "drlio.h"

#include "drllineq.h"		/* prototype consistency */

#if defined(USEIMSL)
#include <imsl.h>
#endif

/*f---------------------------------------------------------------------
 * Matrix inversion.
 *
 * <br><br>
 * Inverts a square
 * <i> double **</i> matrix <i> mat</i> of size * <i> nl</i>$x$<i> nl</i>.
 * and puts the result in <i> invMat</i>.
 * Returns 0 iff successful.<br>
 * <b> Remark: uses the IMSL library. </b>
 */

DLL_EXPORT(int)
DrlMatrixInvert(double ***invMat, double **aMat, int nl)
{
static	char	routine[] = "DrlMatrixInvert";

#if defined(USEIMSL)
	/*
	 * Interface for C-Math Imsl Library
	 */
	int	status = FAILURE;
	int	i, j,
		n = nl;
	long	imslErrCode;
	double	*b1 = NULL,	/* array containing inverse matrix */
		*a1 = NULL;	/* array of containing matrix */ 
	double	**iMat = NULL;

	/**
	 **
	 **/
	if (nl == 1) {
		if (IS_ALMOST_ZERO(aMat[0][0])) {
			GtoErrMsg("%s: singular matrix.\n", routine);
			goto done;
		}
		if ((*invMat = DrlMatrixNew(nl, nl)) == NULL)
			goto done;
		(*invMat)[0][0] = 1e0 / aMat[0][0];
		status = SUCCESS;
		goto done;
	}

	/**
	 **
	 **/
	*invMat = NULL;
	if ((a1 = DrlDoubleVectAlloc(0, (n+1)*(n+1)-1)) == NULL)
		goto done;

	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++)
		a1[i*n+j] = aMat[i][j];

	/*
	 * IMSL error setting
	 */
	imsl_error_options(
		IMSL_SET_STOP, IMSL_FATAL, 0,
		IMSL_SET_STOP, IMSL_TERMINAL, 0,
		IMSL_SET_ERROR_FILE, stdout,
		0);

	/*
	 * Call to the imsl routine
	 */
	imsl_d_lin_sol_gen(n, a1, NULL,
		IMSL_INVERSE, &b1,
		IMSL_INVERSE_ONLY,
		0);
	imslErrCode = imsl_error_code();


	if (imslErrCode != 0L) {
	    switch (imslErrCode) {
	    case IMSL_ILL_CONDITIONED:
		GtoErrMsg("%s: [imsl] %s\n", routine,
			"input matrix ill-conditioned");
#ifdef	__DEBUG__
		fprintf(stderr, "%s: [imsl] %s\n", routine,
			"input matrix ill-conditioned");
		FPrintfDoubleVector(stderr, NULL, a1, n*n);
		FPrintfDoubleDrlMatrix(stderr, NULL, aMat, n, n);
#endif
		break;

	    case IMSL_SINGULAR_MATRIX:
		GtoErrMsg("%s: [imsl] %s\n", routine,
			"input matrix singular");
#ifdef	__DEBUG__
		fprintf(stderr, "%s: [imsl] %s\n", routine,
			"input matrix singular");
		FPrintfDoubleVector(stderr, NULL, a1, n*n);
		FPrintfDoubleDrlMatrix(stderr, NULL, aMat, n, n);
#endif
		break;

	    default:
		GtoErrMsg("%s: [imsl] %s\n", routine, "unknown error");
		break;
	    }
	    GtoErrMsg("%s: IMSL failed (error code %d)\n",
			routine, imslErrCode);

	    goto done;
	}

	/* inversion OK: copy inverse matrix  */
	if ((iMat = DrlMatrixNew(nl, nl)) == NULL)
		goto done;

	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++)
		iMat[i][j] = b1[i*n+j];

	*invMat = iMat;


	/* made it through OK */
	status = SUCCESS;
done:
	DrlDoubleVectFree(a1, 0, 0);
	if (b1 != NULL) free((char*) b1);
	if (status != SUCCESS) {
		GtoErrMsg("%s: input matrix to be inverted:\n", routine);
		DrlFPrintDoubleMatr(NULL, NULL, aMat, nl, nl);
		GtoErrMsg("%s: failed.\n", routine);
		*invMat = NULL;
	}
	return(status);
#else
	GtoErrMsg("%s: routine not available\n", routine);
	return(FAILURE);
#endif
}



