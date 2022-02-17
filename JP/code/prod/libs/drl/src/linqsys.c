/************************************************************************
 * Module:	DRL - LINEQ
 * Function:	Linear Equations and Matrix Algebra
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <math.h>
#include <float.h>
#include "drlmem.h"
#include "drlio.h"		/* DrlFPrintXXXX */
#include "drlstr.h"		/* DrlFloatPrint() */

#include "drllineq.h"		/* prototype consistency */

#if defined(USENAG)
extern	void	c_f04atf(double *a, double *b, int *n,
			double *c, int *ifail);
#elif defined(USEIMSL) 
#include <imsl.h>
#endif

/*#define	__DEBUG__*/
#undef	__DEBUG__

/*f---------------------------------------------------------------------
 * Linear system resolution (standard).
 *
 * <br><br>
 * Solves a linear system 
 * <blockquote>
 * sum<sub>j=1,...,n</sub> a<sub>ij</sub> x<sub>j</sub> = b<sub>i</sub>
 * </blockquote>
 * where <i>a</i> is <i>a</i> nxn-matrix and <i>b</i> an n-vector.
 * Returns 0 iff successful. <br>
 */

int
DrlRealLinSys(
	int n,			/* (I) # of dim */
	double **a,		/* (I) system coeffs [0..n-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b)		/* (I) inhomogeneous term [0..n-1] */
{
static	char	routine[] = "DrlRealLinSys";
	int	status = FAILURE;

#if defined(USENAG)

	int	i, j, ifail = 0;
	double	*a1 ;

	/*
	 * memory allocation
	 */
	if ((a1 = DrlDoubleVectAlloc((n+1)*(n+1))) == NULL) 
		goto done;

	/*
	 * Fortran Interface
	 */
	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++)
		a1[i*n+j] = a[i][j] ;


	c_f04atf(a1, b, &n, x, &ifail);
	if (ifail != 0) goto done;


	status = SUCCESS;
	/* made it through OK */
done:
	DrlDoubleVectFree(a1, 0, (n+1)*(n+1)-1);
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed\n", routine);
	}
	return(status);

#elif defined(USEIMSL) 
	/*****************************************
	 * Interface for C-Math Imsl Library
	 *****************************************/

	int	i, j;
	long	imslErrCode;
	double	*xSol = NULL,	/* final point */
		*a1 = NULL;	/* array of containing matrix */ 


	/**
	 ** Simple closed form cases
	 **/
	if (n == 1) {
		if (IS_ALMOST_ZERO(a[0][0])) {
			DrlErrMsg("%s: singular 1x1 matrix.\n", routine);
			goto done;
		}
		x[0] = b[0] / a[0][0];
		status = SUCCESS;
		goto done;
	}



	/**
	 ** General case
	 **/
	if ((a1 = DrlDoubleVectAlloc(0, (n+1)*(n+1)-1)) == NULL) 
		goto done;

	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++)
		a1[i*n+j] = a[i][j] ;

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
	xSol = imsl_d_lin_sol_gen(n, a1, b,
		0);
	imslErrCode = imsl_error_code();


	if ((xSol == NULL) || (imslErrCode != 0L)) {
	    switch (imslErrCode) {
	    case IMSL_ILL_CONDITIONED:
		DrlErrMsg("%s: [imsl] %s\n", routine,
			"input matrix ill-conditioned");
		DrlFPrintDoubleVect(NULL, NULL, a1, n*n);
		DrlFPrintDoubleMatr(NULL, NULL, a, n, n);

		break;
	    case IMSL_SINGULAR_MATRIX:
		DrlErrMsg("%s: [imsl] %s\n", routine,
			"input matrix singular");
		DrlFPrintDoubleVect(NULL, NULL, a1, n*n);
		DrlFPrintDoubleMatr(NULL, NULL, a, n, n);

		break;
	    default:
		DrlErrMsg("%s: [imsl] %s\n", routine, "unknown error");
		break;
	    }
	    DrlErrMsg("%s: [imsl] solving linear system "
			"failed (code %d)\n",
			routine, imslErrCode);

	    status = ((int) imslErrCode == 0L ? SUCCESS : FAILURE);
	    for (i=0; i<=n-1; i++) x[i] = 0.0;
	    goto done;
	} else {
	    for (i=0; i<=n-1; i++) x[i] = xSol[i];
	}



	status = SUCCESS;
	/* made it through OK */
done:
	if (xSol != NULL) free((char*) xSol);
	DrlDoubleVectFree(a1, 0, (n+1)*(n+1)-1);
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed\n", routine);
	}
	return(status);

#elif defined(_DRL_USENRC_)

	ludcmp(a, n, &indx, d);

#else
	DrlErrMsg("%s: routine not available\n", routine);
	return(-1);
#endif
}


/*f---------------------------------------------------------------------
 * Linear system solution checking.
 *
 * <br><br>
 * Computes
 * <blockquote>
 * a*x - b
 * </blockquote>
 * where a is a mxn-matrix,
 * x an n-vector and b an m-vector.
 * <br>
 * Returns 0 iff successful.
 */

int
DrlRealLinSysCheck(
	int m,			/* (I) # of equations */
	int n,			/* (I) # of unknowns */
	double **a,		/* (I) system coeffs [0..m-1][0..n-1] */
	double *x,		/* (I) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	double *axmb)		/* (O) A*X-B [0..m-1] */
{
/*static	char	routine[] = "DrlRealLinSysCheck";*/
	int	status = FAILURE;
	int	i, j;

	for (i=0; i<=m-1; i++) {
		axmb[i] = 0e0;
		for (j=0; j<=n-1; j++)
			axmb[i] += a[i][j] * x[j];
	}
	if (b != NULL) {
		for (i=0; i<=m-1; i++) 
			axmb[i] -= b[i];
	}

	status = SUCCESS;
	return(status);
}


/*f---------------------------------------------------------------------
 * Linear system pivoting.
 *
 * <br><br>
 * Performs a pivot of $p$ equations of the system
 * <blockquote>
 * sum<sub>j=1,...,n</sub> a<sub>ij</sub> x<sub>j</sub> - b<sub>i</sub> = 0, i=0,...m-1
 * </blockquote>
 * where a is a mxn-matrix and b an m-vector,
 * and creates the pivoted (m-p)x(n-p)-matrix <i> ap</i>
 * and vectors <i>bp</i> (of length p).
 * It also computes a nx(n-p)-matrix (c<sub>ij<sub>)
 * and a vector (d<sub>i</sub>) of length p  such that
 * if y = (y<sub>j</sub>)<sub>1&lt;= j &lt;= n-p</sub> is a solution of the
 * reduced system, then
 * x = cy-d is a solution of 
 * the original system.
 * <br>
 * Returns SUCCESS/FAILURE<br>
 */

int
DrlRealLinSysPivot(
	int m,			/* (I) # of equations */
	int n,			/* (I) # of unknowns */
	double **a,		/* (I) system coeffs [0..m-1][0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int npiv,		/* (I) # of pivots */
	int *ipiv,		/* (I) array of line to pivot [0..npiv] */
	int *jpiva,		/* (I) array col to pivot (NULL if implicit) */
	double ***apo,		/* (O) piv system [0..m-1-npiv][0..n-1-npiv] */
	double **bpo,		/* (O) piv inhomogeneous term [0..m-1-npiv] */
	double ***cpo,		/* (O) sol transf matr [0..n-1][0..n-1-npiv] */
	double **dpo)		/* (O) sol transf 2ndt [0..n-1] */
{
static	char	routine[] = "DrlRealLinSysPivot";
	int	status = FAILURE;

	double	**at = NULL,
		**ut = NULL,
		val,
		cond, condmax;
	int	i, j, k, l,
		mp, np,		/* size after pivot */
		*jpiv = NULL,
		*jpivt = NULL,
		*iciv = NULL,
		*jciv = NULL;

	/*
	 *
	 */
#ifdef	__DEBUG__
	DrlFPrintf(stdout, "%s: Input System\n", routine);
	DrlFPrintf(stdout, "a:\n");
	DrlMatrixPrint(stdout, a, m, n);
	DrlFPrintf(stdout, "b:\n");
	DrlFPrintDoubleVect(stdout, NULL, b, m);
#endif

	*apo = NULL;
	*bpo = NULL;
	*cpo = NULL;
	*dpo = NULL;

	/*
	 * 
	 */
	if (npiv <= 0) {
		mp = m;
		np = n;
		/* Allocate output matrices */
		if ((*apo = DrlDoubleMatrAlloc(0, m-1, 0, n-1)) == NULL)
			goto done;
		if ((*bpo = DrlDoubleVectAlloc(0, m-1)) == NULL)
			goto done;
		if ((*cpo = DrlDoubleMatrAlloc(0, n-1, 0, n-1)) == NULL)
			goto done;
		if ((*dpo = DrlDoubleVectAlloc(0, n-1)) == NULL)
			goto done;

		for (i=0; i<=m-1; i++)
		for (j=0; j<=n-1; j++) 
			(*apo)[i][j] = a[i][j];

		for (i=0; i<=m-1; i++)
			(*bpo)[i] = b[i];

		/* solution transform */
		for (i=0; i<=n-1; i++)
		for (j=0; j<=n-1; j++)
			(*cpo)[i][j] = (i == j ? 1e0 : 0e0);
		for (i=0; i<=n-1; i++)
			(*dpo)[i] = 0e0;

		status = SUCCESS;
		goto done;
	}

	/*
	 *
	 */
	if ((m <= 1) || (n <= 1) || (npiv >= m) || (npiv >= n)) {
		DrlErrMsg("%s: can't %d-pivot system %d x %d.\n",
			routine, npiv, m, n);
		goto done;
	}
	for (i=0; i<=npiv-1; i++) {
	    if ((ipiv[i] < 0) || (ipiv[i] >= m)) {
		DrlErrMsg("%s: wrong i-pivot #%d (%d) for system %d x %d.\n",
			routine, i, ipiv[i], m, n);
		goto done;
	    }
	}



	/* Allocate tmp working space */
	if ((at = DrlDoubleMatrAlloc(0, npiv-1, 0, npiv-1)) == NULL)
		goto done;
	if ((jpiv = DrlIntVectAlloc(0, n-1)) == NULL)
		goto done;

	/* Get j-indices of pivot elements */
	if (jpiva != NULL) {
	    /* j-indices of pivot elements given */
	    for (j=0; j<=npiv-1; j++)
		jpiv[j] = jpiva[j];

	} else {
	    /* Look for best pivot: check all combinations
	     * of indices to get lowest conditionant
	     */
	    if ((jpivt = DrlIntVectAlloc(0, n-1)) == NULL)
		goto done;


	    if (drlInitPermutation(npiv, jpivt, n) != SUCCESS)
		goto done;

	    condmax = 1e32;
	    do {
		for (i=0; i<=npiv-1; i++)
		for (j=0; j<=npiv-1; j++)
			at[i][j] = a[ipiv[i]][jpivt[j]];

		if (DrlMatrixCond(at, npiv, npiv, &cond) != SUCCESS)
			cond = 1e32;

#ifdef	__DEBUG__
		DrlFPrintf(stdout, "%s: ======== Checking \n", routine);
		DrlFPrintf(stdout, "ipiv:\n");
		DrlFPrintStruct(stdout, NULL, '\t',
			DRL_CVECTOR_T, "ipiv", npiv, DRL_INT_T, (void*)ipiv,
			DRL_NULL_T);
		DrlFPrintf(stdout, "jpiv:\n");
		DrlFPrintStruct(stdout, NULL, '\t',
			DRL_CVECTOR_T, "jpiv", npiv, DRL_INT_T, (void*)jpivt,
			DRL_NULL_T);
		DrlFPrintf(stdout, "\t COND=%12.8g CONDMAX=%12.8g\n",
			cond, condmax);
#endif



		if (cond <= condmax) {
		    for (j=0; j<=npiv-1; j++)
			jpiv[j] = jpivt[j];
		    condmax = cond;
		}

	    } while (drlAdvancePermutation(npiv, jpivt, n) != FALSE);
	}


	/* Check pivot indices */
	for (j=0; j<=npiv-1; j++) {
	      if ((jpiv[j] < 0) || (jpiv[j] >= n)) {
		DrlErrMsg("%s: wrong i-pivot #%d (%d) for system %d x %d.\n",
			routine, j, jpiv[i], m, n);
		goto done;
	      }
	}



	/* Arrays ipv and jpiv contain the indices of pivot elements
	 * compute map bet new and old indices */
	mp = m - npiv;
	np = n - npiv;
	if ((iciv = DrlIntVectAlloc(0, mp-1)) == NULL)
		goto done;
	if ((jciv = DrlIntVectAlloc(0, np-1)) == NULL)
		goto done;

	if (drlComplPermutation(npiv, ipiv, m, iciv) != SUCCESS)
		goto done;
	if (drlComplPermutation(npiv, jpiv, n, jciv) != SUCCESS)
		goto done;

#ifdef	__DEBUG__
	DrlFPrintf(stdout, "%s: ====== Final Pivot Selection:\n", routine);
	DrlFPrintf(stdout, "ipiv:\n");
	DrlFPrintStruct(stdout, NULL, '\t',
	    DRL_CVECTOR_T, "ipiv", npiv, DRL_INT_T, (void*)ipiv, DRL_NULL_T);
	DrlFPrintf(stdout, "jpiv:\n");
	DrlFPrintStruct(stdout, NULL, '\t',
	    DRL_CVECTOR_T, "jpiv", npiv, DRL_INT_T, (void*)jpiv, DRL_NULL_T);
	DrlFPrintf(stdout, "iciv:\n");
	DrlFPrintStruct(stdout, NULL, '\t',
	    DRL_CVECTOR_T, "iciv", mp, DRL_INT_T, (void*)iciv, DRL_NULL_T);
	DrlFPrintf(stdout, "jciv:\n");
	DrlFPrintStruct(stdout, NULL, '\t',
	    DRL_CVECTOR_T, "jciv", np, DRL_INT_T, (void*)jciv, DRL_NULL_T);

#endif

	/* invert pivot matrix */
	for (i=0; i<=npiv-1; i++)
	for (j=0; j<=npiv-1; j++)
		at[i][j] = a[ipiv[i]][jpiv[j]];
	if (DrlMatrixInvert(&ut, at, npiv) != SUCCESS) {
		PROGRAM_BUG();
		goto done;
	}

#ifdef	__DEBUG__
	DrlFPrintf(stdout, "at:\n");
	DrlMatrixPrint(stdout, at, npiv, npiv);
	DrlFPrintf(stdout, "ut:\n");
	DrlMatrixPrint(stdout, ut, npiv, npiv);
#endif

	/* Allocate output matrices */
	if ((*apo = DrlDoubleMatrAlloc(0, mp-1, 0, np-1)) == NULL)
		goto done;
	if ((*bpo = DrlDoubleVectAlloc(0, mp-1)) == NULL)
		goto done;
	if ((*cpo = DrlDoubleMatrAlloc(0, n-1, 0, np-1)) == NULL)
		goto done;
	if ((*dpo = DrlDoubleVectAlloc(0, n-1)) == NULL)
		goto done;


	/* reduced system */
	for (i=0; i<=mp-1; i++)
	for (j=0; j<=np-1; j++) {
		(*apo)[i][j] = a[iciv[i]][jciv[j]];
		val = 0e0;
		for (k=0; k<=npiv-1; k++)
		for (l=0; l<=npiv-1; l++)
			val += ut[k][l] * a[ipiv[l]][jciv[j]]
				* a[iciv[i]][jpiv[k]];
		(*apo)[i][j] -= val;
	}

	/* reduced second term */
	for (i=0; i<=mp-1; i++) {
		(*bpo)[i] = b[iciv[i]];
		val = 0e0;
		for (k=0; k<=npiv-1; k++)
		for (l=0; l<=npiv-1; l++)
			val += ut[k][l] * b[ipiv[l]] * a[iciv[i]][jpiv[k]];
		(*bpo)[i] -= val;
	}


	/* solution transform */
	for (i=0; i<=n-1; i++)
	for (j=0; j<=np-1; j++)
		(*cpo)[i][j] = 0e0;
	for (i=0; i<=n-1; i++)
		(*dpo)[i] = 0e0;

	for (i=0; i<=np-1; i++) 
		(*cpo)[jciv[i]][i] = 1e0;

#ifdef	__DEBUG__
	DrlFPrintf(stdout, "%s, %d\n", __FILE__,__LINE__);
	DrlMatrixPrint(stdout, *cpo, n, np);
#endif

	for (k=0; k<=npiv-1; k++) {
		i = jpiv[k];


		for (j=0; j<=np-1; j++) {
			val = 0e0;
			for (l=0; l<=npiv-1; l++)
				val += ut[k][l] * a[ipiv[l]][jciv[j]];
			(*cpo)[i][j] = (-val);
		}

		val = 0e0;
		for (l=0; l<=npiv-1; l++)
			val += ut[k][l] * b[ipiv[l]];
		(*dpo)[i] = (-val);
	}

#ifdef	__DEBUG__
	DrlFPrintf(stdout, "%s, %d\n", __FILE__,__LINE__);
	DrlMatrixPrint(stdout, *cpo, n, np);
#endif






	/* made it through OK */
	status = SUCCESS;
done:
	DrlDoubleMatrFree(at, 0, npiv-1, 0, npiv-1);
	DrlIntVectFree(jpiv, 0, n-1);
	DrlIntVectFree(jpivt, 0, n-1);
	DrlIntVectFree(iciv, 0, mp-1);
	DrlIntVectFree(jciv, 0, np-1);
	DrlDoubleMatrFree(ut, 0, npiv-1, 0, npiv-1);

	if (status != SUCCESS) {
		DrlDoubleMatrFree(*apo, 0, mp-1, 0, np-1);
		DrlDoubleVectFree(*bpo, 0, mp-1);
		DrlDoubleMatrFree(*cpo, 0, n-1, 0, np-1);
		DrlDoubleVectFree(*dpo, 0, n-1);

		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


/*f---------------------------------------------------------------------
 * Linear system evaluate.
 *
 * <br><br>
 * Prints a report on the mxn linear system 
 * <blockquote>
 * a*x-b=0
 * </blockquote>
 * on the file pointer <i> fp</i>.
 */

int
DrlRealLinSysLog(
	int m,			/* (I) # of equations */
	int n,			/* (I) # of unknowns */
	double **a,		/* (I) system coeffs [0..m-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	FILE *fp)		/* (I) */
{
static	char	routine[] = "DrlRealLinErrMinCSolve";
	int	status = FAILURE;

	int	i, j;
	double	*axmb = NULL;
	double	normx, normaxmb;

	/*
	 *
	 */
	if ((axmb = DrlDoubleVectAlloc(0, m-1)) == NULL) goto done;

	DrlFPrintf(fp, "%s: AX-B:\n", routine);
	for (i=0; i<=m-1; i++) {
		axmb[i] = 0e0;
		for (j=0; j<=n-1; j++)
			axmb[i] += x[j] * a[i][j];
		axmb[i] -= b[i];
		DrlFPrintf(fp, "\t[%3d] %8.4f (%6.4f%%)\n", i,
			axmb[i], 1e2*axmb[i] / b[i]);
	}

	/* norms */
	DrlVectNorm(x,    n,   "L1", &normx);
	DrlVectNorm(axmb, m,   "L1", &normaxmb);
	DrlFPrintf(fp, "%s: |X|L1   = %10.6g   |AX-B|L1   = %10.6g\n",
		routine, normx, normaxmb);

	DrlVectNorm(x,    n,   "L2", &normx);
	DrlVectNorm(axmb, m,   "L2", &normaxmb);
	DrlFPrintf(fp, "%s: |X|L2   = %10.6g   |AX-B|L2   = %10.6g\n",
		routine, normx, normaxmb);

	DrlVectNorm(x,    n, "LINF", &normx);
	DrlVectNorm(axmb, m, "LINF", &normaxmb);
	DrlFPrintf(fp, "%s: |X|LINF = %10.6g   |AX-B|LINF = %10.6g\n",
		routine, normx, normaxmb);

	DrlVectNorm(x,    n, "LMIN", &normx);
	DrlVectNorm(axmb, m, "LMIN", &normaxmb);
	DrlFPrintf(fp, "%s: |X|LMIN = %10.6g   |AX-B|LMIN = %10.6g\n",
		routine, normx, normaxmb);

	DrlVectNorm(x,    n, "LMAX", &normx);
	DrlVectNorm(axmb, m, "LMAX", &normaxmb);
	DrlFPrintf(fp, "%s: |X|LMAX = %10.6g   |AX-B|LMAX = %10.6g\n",
		routine, normx, normaxmb);



	status = SUCCESS;
	/* made it through OK */
done:
	DrlDoubleVectFree(axmb, 0, m-1);
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed\n", routine);
	}
	return(status);
}

