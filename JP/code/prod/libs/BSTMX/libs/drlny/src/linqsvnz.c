/************************************************************************
 * Module:	DRL - LINEQ
 * Function:	Linear Equations and Matrix Algebra
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <math.h>
#include <float.h>
#include "drlmem.h"
#include "drlio.h"		/* DrlFPrintXXXX */
#include "drllno.h"		/* DrlNLinProg() */
#include "drlstr.h"		/* DrlFloatPrint() */

#include "drllineq.h"		/* prototype consistency */

#define	__DEBUG__
#undef	__DEBUG__


/*----------------------------------------------------------------------
 * Initializes a multi-index <i> idx</i> (vector of length <i> numNz</i>
 * containing indices (i_0,...\,i_<i> numNz</i>-1)
 * for an iteration over all permutations of <i> numNz</i> elements
 * on a total of <i> n</i>.
 */

DLL_EXPORT(int)
drlInitPermutation(int numNz, int *idx, int n)
{
static	char	routine[] = "drlInitPermutation";
	int	inz;

	if (numNz > n) {
		GtoErrMsg("%s: combination size (%d) > num. elements (%d).\n",
			routine, numNz, n);
		return(FAILURE);
	}

	for (inz=0; inz<=numNz-1; inz++) {
		idx[inz] = inz;
	}
	return(SUCCESS);
}

/*----------------------------------------------------------------------
 * Advance a multi-index <i> idx</i> (vector of length <i> numNz</i>)
 * to the next permutation.
 * Returns TRUE if advance successful, FALSE if terminated.
 */

DLL_EXPORT(int)
drlAdvancePermutation(int numNz, int *idx, int n)
{
	int	inz, j;

	for (inz=numNz-1; inz>=0; inz--) {
		/* increment */
		idx[inz]++;
		/* check past indices */
		for (j=inz+1; j<=numNz-1; j++) {
			idx[j] = idx[j-1]+1;
			if (idx[j] >= n) {
				/* need to increment previous */
				goto do_next;
			}
		}

		if (idx[inz] <= n-1) {
			return(TRUE);
		}

		/* need to increment previous */
do_next:
		;
	}
	/* finished */
	return(FALSE);
}


/*----------------------------------------------------------------------
 * Given a permutation 
 * 1&lt;= i_1 &lt;= ...&lt;= i_{NZ} &lt;= n,
 * computes the complement indices 
 * 1&lt;= j_1 &lt;= ...&lt;= j_{n-\mbox{NZ}} &lt;= n
 * such that i_k \neq j_l for all 
 * 1&lt;= k \mbox{NZ} and 1&lt;= k &lt;= n-\mbox{NZ}.
 */

DLL_EXPORT(int)
drlComplPermutation(int numNz, int *idx, int n, int *idxComp)
{
	int	i, k, ic, isinidx;
	ic = 0;
	for (i=0; i<=n-1; i++) {

		/* check if in idx */
		isinidx = 0;
		for (k=0; k<=numNz-1; k++)
			if (idx[k] == i) { isinidx = 1; break; }
		if (isinidx == 0) {
			idxComp[ic] = i;
			ic++;
		}
	}
	if (ic != n - numNz) {PROGRAM_BUG(); return(FAILURE);}
	return(SUCCESS);
}


/*f---------------------------------------------------------------------
 * Linear system solving with non-zero elements and regularisation.
 *
 * <br><br>
 * Solves a possibly ill conditioned linear system by
 * singular value decomposition and selecting at most
 * <i> numNz</i> non zero terms in the solution.
 * The best solution is obtained by checking all combination
 * and keeping the solution <i>x</i> that minimizes Ax-b$ in L<sup>2</sup> sense.
 * Each subsystem (of size $n x N_z$) is solved 
 * with SVD and regarisation parameters given by arguments
 * <i> regType</i>, <i> alpha</i> and <i> param</i>
 * (see <i> DrlMatrixSvdRegularizeEV</i>).
 */

DLL_EXPORT(int)
DrlRealLinSysSvdSelectNz(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *x,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int numNz,		/* (I) max non zero terms required */
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double alpha,		/* (I) see DrlMatrixSvdRegularizeEV */
	double param)		/* (I) see DrlMatrixSvdRegularizeEV */
{
static	char	routine[] = "DrlRealLinSysSvdSelectNz";
	int	status = FAILURE;


	double	**aR = NULL,		/* shrunk matrix */
		*xR = NULL,		/* shrunk solution */
		ax,			/* tmp storage product Ax */
		l2norm,
		optl2norm;
	int	i, j, *idx = NULL;		/* indices storage */

	/*
	 *
	 */
	if ((aR = DrlMatrixNew(m, numNz)) == NULL)
		goto done;
	if ((xR = DrlDoubleVectAlloc(0, numNz-1)) == NULL)
		goto done;
	if ((idx = DrlIntVectAlloc(0, numNz-1)) == NULL)
		goto done;


	for (j=0; j<=n-1; j++)
		x[j] = 0e0;

	if (numNz == 0) {
		status = SUCCESS;
		goto done;
	}

	optl2norm = 1e32;


	if (drlInitPermutation(numNz, idx, n) != SUCCESS)
		goto done;

	do {
		/* shrink matrix */
		for (i=0; i<=m-1; i++)
		for (j=0; j<=numNz-1; j++)
			aR[i][j] = a[i][idx[j]];

		/* do SVD */
		if (DrlRealLinSysSvd(
			m,
			numNz,
			aR,
			xR,
			b,
			regType,
			alpha,
			param) != SUCCESS)
				goto done;
		/* compute l2 norm */
		l2norm = 0e0;
		for (i=0; i<=m-1; i++) {
			ax = 0e0;
			for (j=0; j<=numNz-1; j++)
				ax +=aR[i][j]*xR[j];
			ax -= b[i];
			l2norm += ax*ax;
		}
		l2norm = sqrt(l2norm);

		/* check l2 norm optimal */
		if (l2norm <= optl2norm) {
			for (j=0; j<=n-1; j++)
				x[j] = 0e0;
			for (j=0; j<=numNz-1; j++)
				x[idx[j]] = xR[j];
			optl2norm = l2norm;

		}


	} while (drlAdvancePermutation(numNz, idx, n) != FALSE);


	/* made it through OK */
	status = SUCCESS;
done:
	DrlMatrixFree(aR, m, numNz);
	DrlDoubleVectFree(xR, 0, numNz-1);
	DrlIntVectFree(idx, 0, numNz-1);


	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*f---------------------------------------------------------------------
 * Linear system solving with pivoting and regularisation.
 *
 * <br><br>
 * Solves a possibly ill conditioned linear system by
 * singular value decomposition.
 * The <i>mx n</i> system is first reduced to a
 * (m-n<sub>piv</sub>)x(n-n<sub>piv</sub>) system
 * using a pivot given by the arguments <i> npiv</i>, <i> ipiv</i> and <i> jpiv</i>
 * (see <i> DrlRealLinSysPivot</i>).
 * The best solution is obtained by checking all combination
 * of the subsystem
 * and keeping the solution x (in the original system)
 * that minimizes Ax-b in L<sup>2</sup> sense.
 * Each subsystem (of size (m-n<sub>piv</sub>)x(N_z-n<sub>piv</sub>)
 * is solved  with SVD and regarisation parameters given by arguments
 * <i> regType</i>, <i> alpha</i> and <i> param</i>
 * (see <i> DrlMatrixSvdRegularizeEV</i>).
 */

DLL_EXPORT(int)
DrlRealLinSysSvdPivot(
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
	double regParam)	/* (I) see DrlMatrixSvdRegularizeEV */
{
static	char	routine[] = "DrlRealLinSysSvdPivot";
	int	status = FAILURE;

	int	nR, mR;
	double	*x = NULL,
		*r = NULL,
		**aR = NULL,
		*bR = NULL,
		**cR = NULL,
		*dR = NULL,
		*xR = NULL;
	int	j, jpivmin, jpivmax, jpiv;
	double	rLinf,
		rLinfO;

#ifdef	__DEBUG__
	DrlFPrintf(stdout, "%s: START\n", routine);
	DrlFPrintf(stdout, "npiv=%d, ipiv[0]=%d jpiv[0]=%d\n",
		npiv, (npiv > 0 ? ipiv[0] : -1),
		((jpiva == NULL) || (jpiva[0] == -1) ? -1 : jpiva[0]));
	DrlFPrintf(stdout, "regType=%d, alpha=%g, regParam=%g\n", 
		regType, alpha, regParam);
	DrlFPrintf(stdout, "a:\n");
        DrlFPrintDoubleMatr(stdout, NULL, a, m, n);
	DrlFPrintf(stdout, "b:\n");
        DrlFPrintDoubleVect(stdout, NULL, b, m);
#endif





	/* Allocate working space */
	if ((r = DrlDoubleVectAlloc(0, m-1)) == NULL)
		goto done;
	if ((x = DrlDoubleVectAlloc(0, n-1)) == NULL)
		goto done;
	if ((xR = DrlDoubleVectAlloc(0, n-1)) == NULL)
		goto done;

	/* Set optimal values */
	rLinfO = 1e32;

	/*  pivot system */
	if (npiv > 1) {
		GtoErrMsg("%s: more than 1 pivot not implemented.\n",
			routine);
		goto done;
	}

	/* if one pivot and j-pivot not specified, optimize on pivot */
	if ((jpiva == NULL) || (jpiva[0] == -1)) {
		jpivmin = 0;
		jpivmax = n-1;
	} else {
		jpivmin = jpiva[0];
		jpivmax = jpiva[0];
	}


	/* Check all possible j-pivot */
	for (jpiv=jpivmin; jpiv<=jpivmax; jpiv++)
	{


	    /* Check pivoting can be done */
	    if (IS_ALMOST_ZERO(a[ipiv[0]][jpiv])) {
		continue;
	    }

	    /* perform pivoting */
	    if (DrlRealLinSysPivot(
		m, n, a, b,
		npiv,
		ipiv,
		&jpiv,
		&aR, &bR, &cR, &dR) != SUCCESS)
			goto done;

	    mR = m - npiv;
	    nR = n - npiv;

	    /*
	     * (2) solving reduced mR x nR system aR*xR-bR = xR
	     */
	    if (DrlRealLinSysSvd(
			mR,
			nR,
			aR,
			xR,
			bR,
			regType,
			alpha,
			regParam) != SUCCESS)
				goto done;

	    /* reconstruct solution x of original system */
	    if (DrlRealLinSysCheck(n, nR, cR, xR, dR, x) != SUCCESS)
			goto done;
	
	    /* compute original system residual vector */
	    if (DrlRealLinSysCheck(m, n, a, x, b, r) != SUCCESS)
			goto done;


	    /* compute L1-norm of x and Linf-norm of r  */
	    if (DrlVectNorm(r, m, "L2", &rLinf) != SUCCESS)
				goto done;

	    /* Chheck optimal */
	    if (rLinf <= rLinfO) {
		for (j=0; j<=n-1; j++) xo[j] = x[j];
		rLinfO = rLinf;
	    }

	    DrlDoubleMatrFree(aR, 0, mR-1, 0, nR-1);
	    DrlDoubleVectFree(bR, 0, mR-1);
	    DrlDoubleMatrFree(cR, 0, n-1, 0, nR-1);
	    DrlDoubleVectFree(dR, 0, n-1);

	    aR = NULL;
	    bR = NULL;
	    cR = NULL;
	    dR = NULL;

	} /* for(jpiv */


	/* compute system residual */
	if (axmb != NULL) {
		if (DrlRealLinSysCheck(m, n, a, x, b, axmb) != SUCCESS)
			goto done;
	}


#ifdef	__DEBUG__
	DrlFPrintf(stdout, "x solution:\n");
        DrlFPrintDoubleVect(stdout, NULL, xo, n);
	DrlFPrintf(stdout, "%s: END\n", routine);
#endif


	/* made it through OK */
	status = SUCCESS;
done:
	DrlDoubleMatrFree(aR, 0, m-npiv-1, 0, n-npiv-1);
	DrlDoubleVectFree(bR, 0, m-npiv-1);
	DrlDoubleMatrFree(cR, 0, n-1, 0, n-npiv-1);
	DrlDoubleVectFree(dR, 0, n-1);
	DrlDoubleVectFree(xR, 0, n-1);


	DrlDoubleVectFree(r, 0, m-1);
	DrlDoubleVectFree(x, 0, n-1);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}



/*f---------------------------------------------------------------------
 * Linear system solving with pivoting, non-zero elements and regularisation.
 *
 * <br><br>
 * Solves a possibly ill conditioned linear system by
 * singular value decomposition and selecting at most
 * <i> numNz</i> non zero terms in the solution.
 * The <i>mx n</i> system is first reduced to a
 * (m-n<sub>piv</sub>)x(n-n<sub>piv</sub>) system
 * using a pivot given by the arguments <i> npiv</i>, <i> ipiv</i> and <i> jpiv</i>
 * (see <i> DrlRealLinSysPivot</i>).
 * The best solution is obtained by checking all combination
 * of the subsystem
 * and keeping the solution x (in the original system)
 * that minimizes <i>Ax-b<i> in L<sup>2</sup> sense.
 * Each subsystem (of size (m-n<sub>piv</sub>)x(N_z-n<sub>piv</sub>)
 * is solved  with SVD and regarisation parameters given by arguments
 * <i> regType</i>, <i> alpha</i> and <i>param</i>
 * (see <i> DrlMatrixSvdRegularizeEV</i>).
 */

DLL_EXPORT(int)
DrlRealLinSysSvdSelectNzPivot(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *xo,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int npiv,		/* (I) # of pivots */
	int *ipiv,		/* (I) array of line to pivot [0..npiv] */
	int *jpiva,		/* (I) array col to pivot (NULL if implicit) */
	double *axmb,		/* (O) residual ax-b [0..m-1] (or NULL) */
	int numNzMax,		/* (I) max non zero terms required */
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double alpha,		/* (I) see DrlMatrixSvdRegularizeEV */
	double regParam)	/* (I) see DrlMatrixSvdRegularizeEV */
{
static	char	routine[] = "DrlRealLinSysSvdSelectNzPivot";
	int	status = FAILURE;

	int	nR, mR;
	double	*x = NULL,
		*r = NULL,
		**aR = NULL,
		*bR = NULL,
		**cR = NULL,
		*dR = NULL,
		*xR = NULL,
		**aNz = NULL,
		*xNz = NULL;
	int	i, j, jpivmin, jpivmax, jpiv,
		numNz,
		*idxNz = NULL;
	double	rLinf,
		rLinfO;


#ifdef	__DEBUG__
	DrlFPrintf(stdout, "%s: START\n", routine);
	DrlFPrintf(stdout, "numNzMax=%d, npiv=%d, ipiv[0]=%d jpiv[0]=%d\n",
		numNzMax, npiv, (npiv > 0 ? ipiv[0] : -1),
		((jpiva == NULL) || (jpiva[0] == -1) ? -1 : jpiva[0]));
	DrlFPrintf(stdout, "regType=%d, alpha=%g, regParam=%g\n", 
		regType, alpha, regParam);
	DrlFPrintf(stdout, "a:\n");
        DrlFPrintDoubleMatr(stdout, NULL, a, m, n);
	DrlFPrintf(stdout, "b:\n");
        DrlFPrintDoubleVect(stdout, NULL, b, m);
#endif





	/* Allocate working space */
	if ((r = DrlDoubleVectAlloc(0, m-1)) == NULL)
		goto done;
	if ((x = DrlDoubleVectAlloc(0, n-1)) == NULL)
		goto done;
	if ((xR = DrlDoubleVectAlloc(0, n-1)) == NULL)
		goto done;



	if ((aNz = DrlMatrixNew(m, numNzMax)) == NULL)
		goto done;
	if ((xNz = DrlDoubleVectAlloc(0, numNzMax-1)) == NULL)
		goto done;
	if ((idxNz = DrlIntVectAlloc(0, numNzMax-1)) == NULL)
		goto done;


	/* Set optimal values */
	rLinfO = 1e32;

	/*  pivot system */
	if (npiv > 1) {
		GtoErrMsg("%s: more than 1 pivot not implemented.\n",
			routine);
		goto done;
	}
	/* if one pivot and j-pivot not specified, optimize on pivot */
	if ((jpiva == NULL) || (jpiva[0] == -1)) {
		jpivmin = 0;
		jpivmax = n-1;
	} else {
		jpivmin = jpiva[0];
		jpivmax = jpiva[0];
	}


	/* Check all possible j-pivot */
	for (jpiv=jpivmin; jpiv<=jpivmax; jpiv++)
	{

	    /* Check pivoting can be done */
	    if (IS_ALMOST_ZERO(a[ipiv[0]][jpiv]))
			continue;

	    /* perform pivoting */
	    if (DrlRealLinSysPivot(
		m, n, a, b,
		npiv,
		ipiv,
		&jpiv,
		&aR, &bR, &cR, &dR) != SUCCESS)
			goto done;

	    mR = m - npiv;
	    nR = n - npiv;

	    /*
	     * (2) solving reduced system aR*xR-bR = xR
	     */

	    /* shrink system to only Nz variables */
	    for (numNz = 0; numNz <= numNzMax-npiv; numNz++) {

	        /* loop over combinations of numNz indices */
	        if (drlInitPermutation(numNz, idxNz, nR) != SUCCESS) goto done;
	        do {
		    if (numNz == 0) {
		        for (j=0; j<=nR-1; j++) xR[j] = 0e0;
		    } else {
		        /* shrink matrix */
		        for (i=0; i<=mR-1; i++)
		        for (j=0; j<=numNz-1; j++)
			    aNz[i][j] = aR[i][idxNz[j]];

		        /* perform SVD */
		        if (DrlRealLinSysSvd(
				mR,
				numNz,
				aNz,
				xNz,
				bR,
				regType,
				alpha,
				regParam) != SUCCESS)
					goto done;

		        /* copy solution back to R system */
		        for (j=0; j<=nR-1; j++) xR[j] = 0e0;
		        for (j=0; j<=numNz-1; j++) xR[idxNz[j]] = xNz[j];
		    }

		    /* reconstruct solution x of original system */
		    if (DrlRealLinSysCheck(n, nR, cR, xR, dR, x) != SUCCESS)
				goto done;
	
		    /* compute original system residual vector */
		    if (DrlRealLinSysCheck(m, n, a, x, b, r) != SUCCESS)
				goto done;


		    /* compute L1-norm of x and Linf-norm of r  */
		    /*if (DrlVectNorm(x, n, "L1", &xL1) != SUCCESS)
					goto done;*/
		    if (DrlVectNorm(r, m, "L2", &rLinf) != SUCCESS)
				goto done;

		    /* Chheck optimal */
		    if (rLinf <= rLinfO) {
			for (j=0; j<=n-1; j++) xo[j] = x[j];
			rLinfO = rLinf;
		    }

	        } while (drlAdvancePermutation(numNz, idxNz, nR) != FALSE);



	    } /* for(numNz */

	    DrlDoubleMatrFree(aR, 0, mR-1, 0, nR-1);
	    DrlDoubleVectFree(bR, 0, mR-1);
	    DrlDoubleMatrFree(cR, 0, n-1, 0, nR-1);
	    DrlDoubleVectFree(dR, 0, n-1);

	    aR = NULL;
	    bR = NULL;
	    cR = NULL;
	    dR = NULL;

	} /* for(jpiv */


	/* compute system residual */
	if (axmb != NULL) {
		if (DrlRealLinSysCheck(m, n, a, x, b, axmb) != SUCCESS)
			goto done;
	}


#ifdef	__DEBUG__
	DrlFPrintf(stdout, "x solution:\n");
        DrlFPrintDoubleVect(stdout, NULL, xo, n);
	DrlFPrintf(stdout, "%s: END\n", routine);
#endif


	/* made it through OK */
	status = SUCCESS;
done:
	DrlDoubleMatrFree(aR, 0, m-npiv-1, 0, n-npiv-1);
	DrlDoubleVectFree(bR, 0, m-npiv-1);
	DrlDoubleMatrFree(cR, 0, n-1, 0, n-npiv-1);
	DrlDoubleVectFree(dR, 0, n-1);
	DrlDoubleVectFree(xR, 0, n-1);


	DrlDoubleVectFree(r, 0, m-1);
	DrlDoubleVectFree(x, 0, n-1);

	DrlMatrixFree(aNz, m, numNzMax);
	DrlDoubleVectFree(xNz, 0, numNzMax-1);
	DrlIntVectFree(idxNz, 0, numNzMax-1);


	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}



/*----------------------------------------------------------------------
 * NOT USED (experimental)
 */


DLL_EXPORT(int)
DrlRealLinSysSvdSelectNzPivotOptim(
	int m,			/* (I) number of equations */
	int n,			/* (I) number of unknowns  */
	double **a,		/* (I) input matrix [0..m-1][0..n-1] */
	double *xo,		/* (O) solution [0..n-1] */
	double *b,		/* (I) inhomogeneous term [0..m-1] */
	int npiv,		/* (I) # of pivots */
	int *ipiv,		/* (I) array of line to pivot [0..npiv] */
	int *jpiva,		/* (I) array col to pivot (NULL if implicit) */
	double *axmb,		/* (O) residual ax-b [0..m-1] (or NULL) */
	int numNzMax,		/* (I) max non zero terms required */
	double rLinfMax,	/* (I) max residual linf norm */
	int regType,		/* (I) see DrlMatrixSvdRegularizeEV */
	double regParam,	/* (I) see DrlMatrixSvdRegularizeEV */
	double *alphaO)		/* (O) optimal */
{
static	char	routine[] = "DrlRealLinSysSvdSelectNzPivotOptim";
	int	status = FAILURE;

	int	i, j, nR, mR;
	double	*r = NULL,
		*x = NULL,
		alphaL, alphaU,
		alpha,
		xL1,
		xL1Min,
		rLinf,
		rLinfO;
	double	**aR = NULL,
		*bR = NULL,
		**cR = NULL,
		*dR = NULL,
		*xR = NULL;

	double	**aNz = NULL,
		*xNz = NULL;
	int	jpivmin, jpivmax, jpiv,
		numNz, numNzO,
		*idxNz = NULL;

	/* Allocate working space */
	if ((r = DrlDoubleVectAlloc(0, m-1)) == NULL)
		goto done;
	if ((x = DrlDoubleVectAlloc(0, n-1)) == NULL)
		goto done;

	/* Set optimal values */
	*alphaO = -1e0;
	xL1Min = 1e32;
	rLinfO = 1e32;
	numNzO = -1;


	/* 
	 * (1) pivot system
	 */
	if (npiv > 1) {
		GtoErrMsg("%s: more than 1 pivot not implemented.\n",
			routine);
		goto done;
	}
	/* if one pivot and j-pivot not specified, optimize on pivot */
	if ((jpiva == NULL) || (jpiva[0] == -1)) {
		jpivmin = 0;
		jpivmax = n-1;
	} else {
		jpivmin = jpiva[0];
		jpivmax = jpiva[0];
	}


	/* Check all possible j-pivot */
	for (jpiv=jpivmin; jpiv<=jpivmax; jpiv++)
	{

	    /* Check pivoting can be done */
	    if (IS_ALMOST_ZERO(a[ipiv[0]][jpiv]))
			continue;

	    /* perform pivoting */
	    if (DrlRealLinSysPivot(
		m, n, a, b,
		npiv,
		ipiv,
		&jpiv,
		&aR, &bR, &cR, &dR) != SUCCESS)
			goto done;

	    mR = m - npiv;
	    nR = n - npiv;

	    if ((xR = DrlDoubleVectAlloc(0, nR-1)) == NULL)
		goto done;

	    /*
	     * (2) solving reduced system aR*xR-bR = xR
	     */


	    if ((aNz = DrlMatrixNew(m, numNzMax)) == NULL)
		goto done;
	    if ((xNz = DrlDoubleVectAlloc(0, numNzMax-1)) == NULL)
		goto done;
	    if ((idxNz = DrlIntVectAlloc(0, numNzMax-1)) == NULL)
		goto done;



	    /* shrink system to only Nz variables */
	    for (numNz = 0; numNz <= numNzMax-npiv; numNz++) {
	        /* loop over combinations of numNz indices */
	        if (drlInitPermutation(numNz, idxNz, nR) != SUCCESS) goto done;
	        do {
		    /* shrink matrix */
		    for (i=0; i<=mR-1; i++)
		    for (j=0; j<=numNz-1; j++)
			aNz[i][j] = aR[i][idxNz[j]];

		    /* Check */
		    alphaL = 1e-6;
		    alphaU = 1e+6;
		    alpha = 1e0;
		    do {

			if (DrlRealLinSysSvd(
				mR,
				numNz,
				aNz,
				xNz,
				bR,
				regType,
				alpha,
				regParam) != SUCCESS)
					goto done;

			/* copy solution back to R system */
			for (j=0; j<=nR-1; j++) xR[j] = 0e0;
			for (j=0; j<=numNz-1; j++) xR[idxNz[j]] = xNz[j];

			/* reconstruct solution x of original system */
			if (DrlRealLinSysCheck(n, nR, cR, xR, dR, x) != SUCCESS)
				goto done;
	
			/* compute original system residual vector */
			if (DrlRealLinSysCheck(m, n, a, x, b, r) != SUCCESS)
				goto done;
	
			/* compute L1-norm of x and Linf-norm of r  */
			if (DrlVectNorm(x, n, "L1", &xL1) != SUCCESS)
					goto done;
			if (DrlVectNorm(r, m, "LINF", &rLinf) != SUCCESS)
				goto done;


#ifdef	__DEBUG__
			fprintf(stdout,
				"NZ=%d PER=%d alpha(%8.2e,%8.2e,%8.2e)|"
				" rLinf: %8.4g %8.4g %8.4g"
				" xL1=%8.4g  xL1Min=%8.4g\n",
				numNz, idxNz[0],
				alphaL, alpha, alphaU,
				rLinf, rLinfMax, rLinfO,
				xL1, xL1Min);
#endif

			/* stop optimization */
			if (rLinf <= rLinfMax) {
				/* rLinf <= rLinfMax */

				if (xL1 < xL1Min) {
					for (j=0; j<=n-1; j++) xo[j] = x[j];
					numNzO = numNz;
					*alphaO = alpha;
					rLinfO = rLinf;
					xL1Min = xL1;
				}

				/* accurate enough, INCREASE reg */
				alphaL = alpha;
				alpha = sqrt(alpha*alphaU);

			} else {
				/* rLinf >  rLinfMax */

				if ((rLinf <= rLinfO) && (rLinfO > rLinfMax)) {
					for (j=0; j<=n-1; j++) xo[j] = x[j];
					numNzO = numNz;
					*alphaO = alpha;
					rLinfO = rLinf;
					xL1Min = xL1;
				}

				/* NOT accurate enough, DECREASE reg */
				alphaU = alpha;
				alpha = sqrt(alpha*alphaL);
			}



		    } while (log(alphaU/alphaL) > 0.1);

#ifdef	__DEBUG__
		    fprintf(stdout, "optimal alpha %8.4g\n", *alphaO);
#endif

	        } while (drlAdvancePermutation(numNz, idxNz, nR) != FALSE);
	    }
	}

#ifdef	__DEBUG__
	fprintf(stdout, "optimal: ");
	fprintf(stdout, "numNz=%d  rLinf=%lf xL1=%lf\n",
			numNzO, rLinfO, xL1Min);
#endif



	/* compute system residual */
	if (axmb != NULL) {
		if (DrlRealLinSysCheck(m, n, a, x, b, axmb) != SUCCESS)
			goto done;
	}



	/* made it through OK */
	status = SUCCESS;
done:

	DrlDoubleVectFree(r, 0, m-1);
	DrlDoubleVectFree(x, 0, n-1);
	DrlDoubleVectFree(xR, 0, nR-1);

	DrlDoubleMatrFree(aR, 0, mR-1, 0, nR-1);
	DrlDoubleVectFree(bR, 0, mR-1);
	DrlDoubleMatrFree(cR, 0, n-1, 0, nR-1);
	DrlDoubleVectFree(dR, 0, n-1);


	DrlMatrixFree(aNz, m, numNzMax);
	DrlDoubleVectFree(xNz, 0, numNzMax-1);
	DrlIntVectFree(idxNz, 0, numNzMax-1);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}




