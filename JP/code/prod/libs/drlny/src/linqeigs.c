/************************************************************************
 * Module:	DCU - LINQ
 * Function:	Eigen stuff
 * Author:	C. Daher
 * REMARK:	<<< call to NAG Fortran library >>
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <math.h>

#include "drlmem.h"

#include "drllineq.h"		/* prototype consistency */

#if defined(USEIMSL)
# include "imsl.h"
#endif

#if defined(USENAG)
extern	void	c_f04axf(int *, ...) ;
#endif

/*f---------------------------------------------------------------------
 * Matrix eigenvalues and vector.
 * 
 * <br><br>
 * This routines computes the eigen values and vectors
 * of the square matrix <i> a</i> of size n.
 * The <i>n</i> eigen values are placed in the vector <i> vap</i> and
 * the eigen vectors in <i> vep</i>  (<i> vep[i]</i> is the vector 
 * for eigen value <i> vap[i]</i>).
 * Both <i> vap</i> and <i> vep</i> should be allocated <i> before</i> the 
 * routine call.
 * The routine returns a failure if one of the eigenvalues
 * has a non-zero imaginary part.
 * Returns 0 on success.<br>
 * <b> Remark: Must be linked with IMSL</b>.
 */

DLL_EXPORT(int)
DrlMatrixRealEigenVect(
	int n,		/* (I) dimension */
	double **a,	/* (I) matrix */
	double *vap,	/* (O) array of eigen values */
	double **vep)	/* (O) array of eigen vectors */
{
static	char	routine[] = "DrlMatrixRealEigenVect";
#if defined(USENAG)
	int	i, j, ifail;
	double	*vfar,	/* contains Re a[i][j] */
		*vfai,	/* Im a[i][j] */
		*vfvr,	/* Re v_i[j] (i-th eigen vec) */
		*vfvi ;	/* Im v_i[j] (i-th eigen vec) */

	ifail = 0 ;
	vfar = NULL;
	vfai = NULL;
	vfvr = NULL;
	vfvi = NULL;

	/*
	 * Allocate memory
	 */
	if ( ((vfar = DoubleVectAlloc(0, n*n-1)) == NULL)
	  || ((vfai = DoubleVectAlloc(0, n*n-1)) == NULL)
	  || ((vfvr = DoubleVectAlloc(0, n*n-1)) == NULL)
	  || ((vfvi = DoubleVectAlloc(0, n*n-1)) == NULL)) {
		ifail = -4 ;
		goto done;
	}


	/*
	 * transfer matrix
	 */
	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++) {
		vfar[i*n+j] = a[i][j] ; /* Re a[i][j] */
		vfai[i*n+j] = 0.0 ;	/* Im a[i][j] */
	}


	/*
	 * Appel routine fortran
	 */
	c_f04axf(&n, vfar, vfai, vap, vfvr, vfvi, &ifail) ;


	/*
	 * Transfer vectors
	 */
	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++) {
		vep[i][j] = vfvr[i*n+j] ;
	}


	/*
	 * FREE memory
	 */

done:
	if (vfar != NULL) DoubleVectFree(0, vfar-1);
	if (vfai != NULL) DoubleVectFree(0, vfai-1);
	if (vfvr != NULL) DoubleVectFree(0, vfvr-1);
	if (vfvi != NULL) DoubleVectFree(0, vfvi-1);

	if (ifail != 0) {
		GtoErrMsg("%s: Failed\n", routine);
	}

	return(ifail);
#elif defined(USEIMSL)
	int		status = FAILURE;
	int		i, j;
	double		*avec = NULL;
	d_complex	*evalu,
			*evecu,
			*eval;

	if ((avec = DrlDoubleVectAlloc(0, n*n-1)) == NULL)
		goto done;
	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++)
		avec[n*i+j] = a[i][j];
		
	if ((evalu = NEW_ARRAY(d_complex, n)) == NULL)
		goto done;
	if ((evecu = NEW_ARRAY(d_complex, n*n)) == NULL)
		goto done;


	eval = imsl_d_eig_gen(n, avec,
		IMSL_VECTORS_USER, evecu,
		IMSL_RETURN_USER, evalu,
		0);
	if (eval == NULL) {
		GtoErrMsg("%s: IMSL failure.\n", routine);
		goto done;
	}


	for (i=0; i<=n-1; i++) {
	    if (!IS_ALMOST_ZERO(evalu[i].im)) {
		GtoErrMsg("%s: eigen val #%d (%lf,%lf) has im part.\n",
			routine, i,
	    		evalu[i].re,
	    		evalu[i].im);
		goto done;
	    }
	    vap[i] = evalu[i].re;
	}


	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++) {
		vep[i][j] = evecu[n*i+j].re;
	}


	status = SUCCESS;
done:
	if (evalu != NULL) FREE(evalu);
	if (evecu != NULL) FREE(evecu);
	DrlDoubleVectFree(avec, 0, n*n-1);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
#else

	GtoErrMsg("%s: not implemented\n", routine);
	return(-1);
#endif
}
