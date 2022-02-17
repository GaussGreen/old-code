/****************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	
 * Function:	
 * Author:	Christian Daher, David Liu
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include "date_sup.h"

#include "drlio.h"		/* FScanStruct */
#include "drlstr.h"
#include "drltime.h"		/* TDateAdvanceToRefDate */
#include "drlmem.h"
#include "drlsmat.h"
#include "drlsort.h"
#include "drlts.h"		/* DrlTCurveWrap() */
#include "drllineq.h"		/* DrlMatrixPrint  */


#define USE_OOQP		/* USE_OOQP, USE_IMSL, USE_BOOTSTRP */


#if !defined(NO_OOQP) 		/* if you don't have OOQP! (CD) */
#include "cQpGenSparse.h"	/* OOQP  */
#endif

#include "imsl.h"		/* IMSL  */

#define	_vnfm_SOURCE
#include "vnfmcali.h"
#include "vnfmopca.h"
#include "vnfmwrap.h"	/* Prototype Consistency */



extern  int		_VnfmCalibSpotVolEqualFailedExpIdx;
extern	double		_VnfmCalibSpotVolEqualFailedExp,
			_VnfmCalibSpotVolEqualFailedMat,
			_VnfmCalibSpotVolEqualAdjFailedVol;

#undef	__DEBUG__
/*#define	__DEBUG__*/

typedef struct  {
	VnfmData	*tfData;

	double		*volMat;
	int		*volFreq;
	double		*variances;

	double		spotVolMin;
	int   		volBoundType;
	int   		spotVolRatioFlag;
	double 		spotVolRatio;

} TOptimData;  
 

static	TOptimData	optimData;


static int
ComputeIMSLQPCoeff(
	double	*a, 	/* (O) linear constraint coeffs */
	double	*b,	/* (O) linear constraint consts */
	double	*g, 	/* (O) QP linear coeffs */
	double	*h);	/* (O) QP Hessian matrix */


static int
ComputeOOQPCoeff(
			/* Indices of QP coefficients */
	int	*irowQ,		/* row index for Q matrix */
	int	*jcolQ,		/* column index for Q matrix */
	int	*irowA,  	/* row index for A matrix */
	int	*jcolA,		/* column index for A matrix */
	int	*irowC,  	/* row index for C matrix */
	int	*jcolC,		/* column index for C matrix */
				/* Coefficients of QP specifications */
	double	*c,		/* linear vector in obj function */
	double	*dQ,		/* Q matrix */
	double	*dA,		/* A matrix */
	double	*dC,		/* C matrix */
			/* Constraint boundaries */
	char	*ixlow,		/* index of non-zero lower bound of X */
	double	*xlow,		/* lower bound of X */
	char	*ixupp,		/* index of non-zero upper bound of X */
	double	*xupp,		/* upper bound of X */
	char	*iclow,		/* index of lower bound of ineq cons */
	double	*clow,		/* lower bound of ineq cons */
	char	*icupp,		/* index of upper bound of ineq cons */
	double	*cupp,		/* upper bound of ineq cons */
	double	*bA);		/* right-hand side vector for eq-cons */
	


static void
ComputeMatrixC(double	**C);	/* (O) C matrix */


static int _degenerateAndSelectMaxArray(
	double	tolerance,	/* (I) tolerance (positive)	   */
	int	*numItem,       /* (B) # of elements in each array */
        double  *arrayA,        /* (B) array A                     */
        double  *arrayB,        /* (B) array B                     */
        double  *arrayC);




/*----------------------------------------------------------------------
 *
 *      NEW VERSION   Convert to QP problem and use optimization 
 *      routine OOQP to solve. 
 *
 */

static int
ComputeAdjustment_OOQP(
	VnfmData *that,			/* (I) vnfm data */
	TSwaptionMatrix2D *swoptMat,	/* (I) input matrix */

	TDate	*volDates,		/* (I) Benchmark vol dates */
	double	*volMat,		/* (I) Benchmark vol maturities */ 
	int	*volFreq,		/* (I) Benchmark vol frequencies */ 
	double	*volRates,		/* (I) Benchmark vols */

	double spotVolMin,		/* (I) minimum spot volatility */
	double tMatMin,			/* (I) Minimum muturity */

	int   volBoundType,		/* (I) vol adj type: 1, -1, 0  */
	int   spotVolRatioFlag,		/* (I) spot vol ratio check flag  */
	double spotVolRatio,		/* (I) Max spot vol ratio  */

	int *nFailed,			/* (O) # of failures (or NULL) */
	double *tExpFailed,		/* (O) failed exp times (or NULL) */
	double *tMatFailed,		/* (O) failed mat times (or NULL) */
	double *volMinAdj)		/* (O) needed min vol corrections 
					 *     for failed points */
{
static	char		routine[] = "ComputeAdjustment_OOQP";
	int		status = FAILURE;

#if !defined(NO_OOQP) 		/* if you don't have OOQP! (CD) */

#define	SP_VOL_MIN	0.25e-2	     /* minimum spot vol */
#define	TOL_VOL_ADJ	1.0e-5	     /* minimum mkt vol adjustment */
#define	TOL_VOL_RATIO	5.0e-2	     /* 1% extra margin for vol ratio check 
				      *	so ratio of 5 is actually 4.95 */

	/* For convenient use of some macros */

	int		idx, vIdx;
	int		nDim;

	double		impVol, volDiff;
	double		**jt = NULL;


	double		*X 	= NULL; /* optimized variables, spot vols */ 

	int		nVar; 		/* total number of variables */

	int		nnzQ,	/* number of non-zeros in Q matrix */
			nnzA,	/* number of non-zeros in A (equality) matrix */
			nnzC;	/* number of non-zeros in C (inequal matrix */
	
	/* Indices of QP coefficients */
	int		*irowQ = NULL,	/* row index for Q matrix */
			*jcolQ = NULL,	/* column index for Q matrix */
			*irowA = NULL,  /* row index for A matrix */
			*jcolA = NULL,	/* column index for A matrix */
			*irowC = NULL,  /* row index for C matrix */
			*jcolC = NULL;	/* column index for C matrix */

	/* Coefficients of QP specifications */
	double		*c     = NULL,	/* linear vector in obj function */
			*dQ    = NULL,	/* Q matrix */
			*dA    = NULL,	/* A matrix */
			*dC    = NULL;	/* C matrix */

	/* Constraint boundaries */
	char		*ixlow = NULL;	/* index of non-zero lower bound of X */
	double		*xlow  = NULL;	/* lower bound of X */

	char		*ixupp = NULL;	/* index of non-zero upper bound of X */
	double		*xupp  = NULL;	/* upper bound of X */

	int		mz;		/* number of ineqality constraints */
	char		*iclow = NULL;	/* index of lower bound of ineq cons */
	double		*clow  = NULL;	/* lower bound of ineq cons */
	char		*icupp = NULL;	/* index of upper bound of ineq cons */
	double		*cupp  = NULL;	/* upper bound of ineq cons */

	double		*bA    = NULL;	/* right-hand side vector for eq-cons */
	int		my;		/* size of bA */
	

	/* Output */
	double		*gamma = NULL,	/* Lagrangian multipliers l-bound */	
			*phi   = NULL,	/* Lagrangian multipliers u-bound */ 
			*y     = NULL,	/* Lagrangian multipliers eq cons */
			*lambda= NULL,	/* Lagrangian multipliers lb ineq cons*/
			*pi    = NULL,	/* Lagrangian multipliers ub ineq cons*/
			*z     = NULL;	/* z = lambda - pi */

	int		print_level = 0,/* info print level. 0=silent */
			ierr;		/* status: 0 = success */ 

	char		status_ooqp[256];



	/* Check for failed vol calibration */

	if (VnfmVolCalib1VArbitrary(
		that,
		0,
		that->fNDates,
		volMat-1,
		volFreq-1,
		volRates-1,
		LOGVOL,
		FALSE,
		NULL) != SUCCESS ||
	   VnfmVolCalib1VArbitraryMinSpotVol(
		that,
		0,
		that->fNDates,
		volMat-1,
		volFreq-1,
		volRates-1, 
		FALSE,
		spotVolRatioFlag,
		spotVolRatio,
		spotVolMin) != SUCCESS) 
	{

	    /* number of factors */
	    nDim = that->fNf;

	    /* number of spot vols or optimized variables */
	    nVar = NDATES - 1;

	    /* 
	     * Shift by dummy base date at beginning 
	     * so the valid benchmark vols start from 1, ... NDATES-1.
	     */
	    optimData.volMat   = volMat   - 1;
	    optimData.volFreq  = volFreq  - 1;

	    /* Convert the vol to variance = vol^2 * T  */	
	    ASSERT_OR_DONE((optimData.variances = NEW_ARRAY(double, nVar+1)) 
			!= NULL);

	    for (vIdx=1; vIdx<=nVar; vIdx++)
		optimData.variances[vIdx] = 
				volRates[vIdx-1]*volRates[vIdx-1]*TT[vIdx];


	    /*
	     * Spot vol constraints
	     */
	    optimData.spotVolMin       = MAX(spotVolMin, SP_VOL_MIN);
	    optimData.spotVolRatioFlag = spotVolRatioFlag;
            optimData.spotVolRatio     = spotVolRatio;


	    /*
	     * Allocate memory for vol calculation
	     */
	    ASSERT_OR_DONE((jt = DrlDoubleMatrAlloc(
						0,
						nDim*(nDim+1)/2-1,
						0,
						nVar)) != NULL);



	    /** 
	     ** OPTIMIZATION:
	     ** Adjust spot vols subject to the contraints
	     **/
	

	    /* 
	     * Q is a symmetric full matrix. Only the lower triangle 
	     * is specified.
	     */
	    nnzQ = nVar * (nVar + 1) / 2;

	
	    /* 
	     * Number of equality constraints 
	     */
	    my   = 0;
	    nnzA = 0;

	    /* 
	     * Number of inequality constraints: 
	     */
	    mz   = 0;	/* initialize */
	    nnzC = 0;   /* initialize */

	    switch (optimData.volBoundType)
	    {
	    case 1: 	/* adjust vol >= market vol */
	    case -1: 	/* adjust vol <= market vol */
		mz   += nVar;
		nnzC += nVar * (nVar + 1) / 2;
		break;
	    case 0: 	/* no constraints */
		break;
	    default:
		GtoErrMsg("%s: invalid vol bound type (%d).\n",
                        routine, optimData.volBoundType);
		goto done;
	    }

	    /*
	     * 2*(nVar-1) ratio bounds
	     */
	    if (spotVolRatioFlag)	
	    {
		mz   += 2*(nVar-1);
		nnzC += 4*(nVar-1);
	    }


	    /* 
	     * Allocate memory for solver inputs
	     */
	    newQpGenSparse(&c, nVar,
			   &irowQ, nnzQ, &jcolQ, &dQ,
			   &xlow,	&ixlow,
			   &xupp,	&ixupp,
			   &irowA, nnzA, &jcolA, &dA,
			   &bA,   my,
			   &irowC, nnzC, &jcolC, &dC,
			   &clow,  mz,   &iclow,
			   &cupp,        &icupp,
			   &ierr);
			  

	    if (ierr != 0) {
		GtoErrMsg("%s:, Failed to allocate memory for OOQP solver.\n",
			  routine);
		goto done;
	    }


	    /*
	     * Allocate memory for solver outputs
	     * X need an extra element for the vol calculation
	     * in VnfmAvgQBVolSquare, though the last one is obsolete.  
	     */
	    ASSERT_OR_DONE((X     = NEW_ARRAY(double, nVar+1)) != NULL);
	    ASSERT_OR_DONE((gamma = NEW_ARRAY(double, nVar)) != NULL);
	    ASSERT_OR_DONE((phi   = NEW_ARRAY(double, nVar)) != NULL);

	    if (my>0) ASSERT_OR_DONE((y  = NEW_ARRAY(double, my)) != NULL);
	    if (mz>0) {
		ASSERT_OR_DONE((z       = NEW_ARRAY(double, mz)) != NULL);
		ASSERT_OR_DONE((lambda  = NEW_ARRAY(double, mz)) != NULL);
		ASSERT_OR_DONE((pi      = NEW_ARRAY(double, mz)) != NULL);
	    }
	


	    /*
	     * Compute the quadratic coefficients of objective function
	     * and the constraint coefficients
	     */
	    IF_FAILED_DONE(ComputeOOQPCoeff(
				irowQ,
				jcolQ,	
				irowA, 
				jcolA,
				irowC, 
				jcolC,
				c,	
				dQ,
				dA,
				dC,
				ixlow,	
				xlow,
				ixupp,	
				xupp,	
				iclow,	
				clow,	
				icupp,	
				cupp,
				bA));
	
	    /*
	     * Call OOQP solver
	     */
	    qpsolvesp(c, nVar,
		      irowQ, nnzQ, jcolQ, dQ,
		      xlow,	ixlow,
		      xupp,	ixupp,
		      irowA, nnzA, jcolA, dA,
		      bA,     my,
		      irowC, nnzC, jcolC, dC,
		      clow,  mz,   iclow,
		      cupp,        icupp,
		      X,    gamma,   phi,
		      y,
		      z,     lambda,  pi,
		      print_level,	/* print level. No print */
		      &ierr);
			  
	    if (ierr != 0L) {
		switch (ierr) {
		case 1:
			strcpy(status_ooqp, "NOT_FINISHED");
			break;
		case 2:
			strcpy(status_ooqp, "MAX_ITS_EXCEEDED");
			break;
		case 3:
			strcpy(status_ooqp, "INFEASIBLE");
			break;
		case 4:
			strcpy(status_ooqp, "UNKNOWN");
			break;
	        default:
			GtoErrMsg("%s: OOQP failed with invalid "
				  "termination status(%d).\n",
                        	  routine, ierr);
			goto done;
		}

		GtoErrMsg("%s: [OOQP] solving QP system failed (status: %s).\n",
			routine, status_ooqp);
		goto done;
	    }



	   /*
	    * Optimization success, solutions are found,
	    * Now need to compute the implied vol 
	    */


  
	   /* 
	    * Scale back the spot vol^2
	    */
	    for(idx=1; idx<=NDATES-1;idx++)
		    X[idx-1] /= 1.0e2;	


	    /*
	     * Compute the implied market vols 
	     */
	    
	    /* 1. Compute the pseudo-J */
	    if (VnfmComputeJSquare(
	        	that,
	        	X,
	        	jt) != SUCCESS)
		goto done;


	    /* 
	     * Compute the changes if there any and record them
	     */
	    for(idx=1; idx<=NDATES-1;idx++) {
	        /* Compute the benchmark vols */
	        if (VnfmAvgQBVolSquare(
			that,
			0.,
			TT[idx],
			TT[idx],
			optimData.volMat[idx],
			(int)optimData.volFreq[idx],
			X, 
			jt,
			&impVol) != SUCCESS)
		    goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "Vol[%d] = %12.8f%%\n", idx, impVol*1e2);
	)
#endif

		/* Vol difference */
	        volDiff = impVol - volRates[idx-1];

		if (fabs(volDiff) > TOL_VOL_ADJ)
		{
	    	    tExpFailed[*nFailed] = swoptMat->table->dim1Values[idx-1];
		    tMatFailed[*nFailed] = optimData.volMat[idx];
		    volMinAdj[*nFailed]  = volDiff;
		    (*nFailed)++;
		}
	    }


	} /* end of vol calib check */



	/* OK */
	status = SUCCESS;

done:
	if (optimData.variances) FREE(optimData.variances);
	optimData.variances = NULL;

	if (jt)       		DrlDoubleMatrFree(jt,
						  0,
                                                  nDim*(nDim+1)/2-1,
                                                  0,
                                                  nVar);

	freeQpGenSparse(&c,
			&irowQ, &jcolQ, &dQ,
			&xlow,	&ixlow,
			&xupp,	&ixupp,
			&irowA, &jcolA, &dA,
			&bA, 
			&irowC, &jcolC, &dC,
			&clow,  &iclow,
			&cupp,  &icupp);
			  
	if (X      != NULL)	FREE(X);
	if (gamma  != NULL) 	FREE(gamma);
	if (phi    != NULL) 	FREE(phi);
	if (y      != NULL) 	FREE(y);
	if (lambda != NULL) 	FREE(lambda);
	if (pi     != NULL) 	FREE(pi);
	if (z      != NULL) 	FREE(z);

	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);

#else /* defined(USE_OOQP) */
	GtoErrMsg("%s: OOQP not available.\n", routine);
	return (FAILURE);
#endif
}



/*----------------------------------------------------------------------
 *
 *      NEW VERSION   Convert to QP problem and use optimization 
 *      routine IMSL to solve. 
 *
 */
static int
ComputeAdjustment_IMSL(
	VnfmData	*that,		/* (I) vnfm data */
	TSwaptionMatrix2D *swoptMat,	/* (I) input matrix */

	TDate	*volDates,		/* (I) Benchmark vol dates */
	double	*volMat,		/* (I) Benchmark vol maturities */ 
	int	*volFreq,		/* (I) Benchmark vol frequencies */ 
	double	*volRates,		/* (I) Benchmark vols */

	double spotVolMin,		/* (I) minimum spot volatility */
	double tMatMin,			/* (I) Minimum muturity */

	int   volBoundType,		/* (I) vol adj type: 1, -1, 0  */
	int   spotVolRatioFlag,		/* (I) spot vol ratio check flag  */
	double spotVolRatio,		/* (I) Max spot vol ratio  */

	int *nFailed,			/* (O) # of failures (or NULL) */
	double *tExpFailed,		/* (O) failed exp times (or NULL) */
	double *tMatFailed,		/* (O) failed mat times (or NULL) */
	double *volMinAdj)		/* (O) needed min vol corrections */
{
static	char		routine[] = "ComputeAdjustment_IMSL";
	int		status = FAILURE;

#define	SP_VOL_MIN	0.25e-2	     /* minimum spot vol */
#define	TOL_VOL_ADJ	1.0e-5	     /* minimum mkt vol adjustment */
#define	TOL_VOL_RATIO	5.0e-2	     /* 1% extra margin for vol ratio check 
				      *	so ratio of 5 is actually 4.95 */

	/* For convenient use of some macros */
	int		idx, vIdx;
	int		nDim;


	double		impVol, volDiff;
	double		**jt = NULL;

	double		*X 	= NULL; /* optimized variables, spot vols */ 
	double		*sqspv 	= NULL; /* spot vols sqaure with n+1 size */ 

	int		nVar; 		/* total number of variables */

	int		nCon, 	/* total number of constrains */
			nEq;	/* number of equality constrains */

	double		*g	= NULL, /* linear coefficients of obj func */
			*h	= NULL; /* Hessian matrix of obj func */
 
	double		*a	= NULL, /* coeffs of linear constraints */
			*b	= NULL; /* right hand side constraint consts */
 
	long    	imslErrCode;  	/* IMSL error code */




	/* Check for failed vol calibration */

	if (VnfmVolCalib1VArbitrary(
		that,
		0,
		that->fNDates,
		volMat-1,
		volFreq-1,
		volRates-1,
		LOGVOL,
		FALSE,
		NULL) != SUCCESS ||
	   VnfmVolCalib1VArbitraryMinSpotVol(
		that,
		0,
		that->fNDates,
		volMat-1,
		volFreq-1,
		volRates-1, 
		FALSE,
		spotVolRatioFlag,
		spotVolRatio,
		spotVolMin) != SUCCESS) 
	{

	    /* number of factors */
	    nDim = that->fNf;

	    /* number of spot vols or optimized variables */
	    nVar = NDATES - 1;


	    /* 
	     * Shift by dummy base date at beginning 
	     * so the valid benchmark vols start from 1, ... NDATES-1.
	     */
	    optimData.volMat   = volMat   - 1;
	    optimData.volFreq  = volFreq  - 1;

	    /* Convert the vol to variance = vol^2 * T  */	
	    ASSERT_OR_DONE((optimData.variances = NEW_ARRAY(double, nVar+1)) 
			!= NULL);

	    for (vIdx=1; vIdx<=nVar; vIdx++)
		optimData.variances[vIdx] = 
				volRates[vIdx-1]*volRates[vIdx-1]*TT[vIdx];


	    /*
	     * Spot vol constraints
	     */
	    optimData.spotVolMin       = MAX(spotVolMin, SP_VOL_MIN);
	    optimData.spotVolRatioFlag = spotVolRatioFlag;
            optimData.spotVolRatio     = spotVolRatio;


	    /*
	     * Allocate memory for vol calculation
	     */
	    ASSERT_OR_DONE((jt = DrlDoubleMatrAlloc(
						0,
						nDim*(nDim+1)/2-1,
						0,
						nVar)) != NULL);



	    /** 
	     ** OPTIMIZATION:
	     ** Adjust spot vols subject to the contraints
	     **/
	

	    /* 
	     * Number of equality constraints 
	     */
	    nEq = 0;

	    /* 
	     * Total number of constraints: 
	     */

	    /* Lower bounds of variables (spot vols^2) */
	    nCon = nVar;

	    switch (optimData.volBoundType)
	    {
	    case 1: 	/* adjust vol >= market vol */
	    case -1: 	/* adjust vol <= market vol */
		nCon += nVar;
		break;
	    case 0: 	/* no constraints */
		nCon += 0;
		break;
	    default:
		GtoErrMsg("%s: invalid vol bound type (%d).\n",
                        routine, optimData.volBoundType);
		goto done;
	    }

	    /*
	     * 2*(nVar-1) ratio bounds
	     */
	    if (spotVolRatioFlag)	
		nCon += 2*(nVar-1);

	    /*
	     * Allocate memory for constraints
	     */
	    ASSERT_OR_DONE((sqspv  = NEW_ARRAY(double, nVar+1)) != NULL);
	    ASSERT_OR_DONE((g  = NEW_ARRAY(double, nVar)) != NULL);
	    ASSERT_OR_DONE((h  = NEW_ARRAY(double, nVar*nVar)) != NULL);

	    ASSERT_OR_DONE((a = NEW_ARRAY(double, nVar*nCon)) != NULL);
	    ASSERT_OR_DONE((b = NEW_ARRAY(double, nCon)) != NULL);


	    /*
	     * Compute the quadratic coefficients of objective function
	     * and the constraint coefficients
	     */
	    IF_FAILED_DONE(ComputeIMSLQPCoeff(a, b, g, h));



	    /*
	     * IMSL error setting
	     */
            imsl_error_options(
                IMSL_SET_STOP, IMSL_FATAL, 0,
                IMSL_SET_STOP, IMSL_TERMINAL, 0,
                IMSL_SET_ERROR_FILE, stdout,
                0);


	    /*
	     * Call to IMSL routine
	     */
	    X = imsl_d_quadratic_prog(
		nCon, nVar, nEq, a, b, g, h,
		0);

	    imslErrCode = imsl_error_code();


	    if ((X == NULL) || (imslErrCode != 0L)) {
		GtoErrMsg("%s: [imsl] solving QP system failed (code %d)\n",
			routine, imslErrCode);
		goto done;
	    }



	   /*
	    * Optimization success, solutions are found,
	    * Now need to compute the implied vol 
	    */

  
	   /* 
	    * Scale back the spot vol^2
	    *
	    * Memory of X is allocated by IMSL of size nVar,
	    * need nVar + 1 for VnfmAvgQBVolSquare. 
	    */
	    for(idx=1; idx<=NDATES-1;idx++)
	    {
	    	    sqspv[idx-1] = X[idx-1]/1.0e2; 
	    }



	    /*
	     * Compute the implied market vols 
	     */
	    
	    /* 1. Compute the pseudo-J */
	    if (VnfmComputeJSquare(
	        	that,
	        	sqspv,
	        	jt) != SUCCESS)
		goto done;


	    /* 
	     * Compute the changes if there any and record them
	     */
	    for(idx=1; idx<=NDATES-1;idx++) {
	        /* Compute the bechmark vols */
	        if (VnfmAvgQBVolSquare(
			that,
			0.,
			TT[idx],
			TT[idx],
			optimData.volMat[idx],
			(int)optimData.volFreq[idx],
			sqspv, 
			jt,
			&impVol) != SUCCESS)
		    goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "Vol[%d] = %12.8f%%\n", idx, impVol*1e2);
	)
#endif

		/* Vol difference */
	        volDiff = impVol - volRates[idx-1];

		if (fabs(volDiff) > TOL_VOL_ADJ)
		{
	    	    tExpFailed[*nFailed] = swoptMat->table->dim1Values[idx-1];
		    tMatFailed[*nFailed] = optimData.volMat[idx];
		    volMinAdj[*nFailed]  = volDiff;
		    (*nFailed)++;
		}
	    }


	} /* end of vol calib check */



	/* OK */
	status = SUCCESS;
done:
	if (optimData.variances) FREE(optimData.variances);
	optimData.variances = NULL;

	if (jt)       		DrlDoubleMatrFree(jt,
						  0,
                                                  nDim*(nDim+1)/2-1,
                                                  0,
                                                  nVar);

	if (X != NULL)		free((char*) X);
	if (sqspv)		FREE(sqspv);
	if (a)  		FREE(a);
	if (b)  		FREE(b);
	if (g) 		 	FREE(g);
	if (h) 			FREE(h);


	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);


}



/*----------------------------------------------------------------------
 *	THE BOOTSTRAP METHOD DOES NOT INCLUDE ALGORITHM FOR SPOT VOL RATIO
 *	CHECK, THEREFORE OFTEN FAILED UNDER STRESS CASES.  
 */

static int
ComputeAdjustment_BootStrp(
	VnfmData	*tfData,	/* (I) vnfm data */
	TSwaptionMatrix2D *swoptMat,	/* (I) input matrix */

	TDate	*volDates,		/* (I) Benchmark vol dates */
	double	*volMat,		/* (I) Benchmark vol maturities */ 
	int	*volFreq,		/* (I) Benchmark vol frequencies */ 
	double	*volRates,		/* (I) Benchmark vols */

	double spotVolMin,		/* (I) minimum spot volatility */
	double tMatMin,			/* (I) Minimum muturity */

	int   volBoundType,		/* (I) vol adj type: 1, -1, 0  */
	int   spotVolRatioFlag,		/* (I) spot vol ratio check flag  */
	double spotVolRatio,		/* (I) Max spot vol ratio  */

	int *nFailed,			/* (O) # of failures (or NULL) */
	double *tExpFailed,		/* (O) failed exp times (or NULL) */
	double *tMatFailed,		/* (O) failed mat times (or NULL) */
	double *volMinAdj)		/* (O) needed min vol corrections */
{
static	char		routine[] = "ComputeAdjustment_BootStrp";
	int		status = FAILURE;


	/* Call bootstraping routine.
	 * We need to offset by 1 the input vectors */
	if (VnfmVolCalib1VArbitrary(
		tfData,
		0,
		tfData->fNDates,
		volMat-1,
		volFreq-1,
		volRates-1,
		LOGVOL,
		FALSE,
		NULL) != SUCCESS) {

		tExpFailed[*nFailed] = 
	swoptMat->table->dim1Values[_VnfmCalibSpotVolEqualFailedExpIdx];
		tMatFailed[*nFailed] = _VnfmCalibSpotVolEqualFailedMat;
		volMinAdj[*nFailed] = _VnfmCalibSpotVolEqualAdjFailedVol;
		(*nFailed)++;

#ifndef	NO_LOGGING
		GTO_IF_LOGGING(\
	   	    DrlFPrintf(vnfmFpLog, \
			"\t failed at %7.4f yrs into %7.4f yrs.\n" \
			"\t Minimum vol increase is  %7.4f.\n", \
			_VnfmCalibSpotVolEqualFailedExp, \
			_VnfmCalibSpotVolEqualFailedMat, \
			_VnfmCalibSpotVolEqualAdjFailedVol));
#endif

	}

	/*  Check minimum spot volatility */
	if (spotVolMin > 0e0 || spotVolRatioFlag) {
		if ( VnfmVolCalib1VArbitraryMinSpotVol(
			tfData,
			0,
			tfData->fNDates,
			volMat-1,
			volFreq-1,
			volRates-1, 
			FALSE,
			spotVolRatioFlag,
			spotVolRatio,
			spotVolMin) != SUCCESS) {

			tExpFailed[*nFailed] = 
	swoptMat->table->dim1Values[_VnfmCalibSpotVolEqualFailedExpIdx];
			tMatFailed[*nFailed] = _VnfmCalibSpotVolEqualFailedMat;
		        volMinAdj[*nFailed] = 0e0;
			(*nFailed)++;

#ifndef	NO_LOGGING
			GTO_IF_LOGGING(\
			    DrlFPrintf(vnfmFpLog, \
				"Min spot volatility failed " \
				"at %7.4f yrs into %7.4f yrs\n", \
				_VnfmCalibSpotVolEqualFailedExp, \
				_VnfmCalibSpotVolEqualFailedMat));
#endif

		}
	}




	/* OK */
	status = SUCCESS;
done:
 
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);

}




/*----------------------------------------------------------------------
 *
 *      NEW VERSION   Convert to QP problem and use optimization 
 *      routine OOQP or IMSL to solve. 
 *	
 *
 */
static int
VnfmCheckSwaptionCalibrationOne(
	TCurve *zcCurve,		/* (I) zero curve */
	TSwaptionMatrix2D *swoptMat,	/* (I) input matrix */

	double *nfParamsL,		/* (I) mr, alphas, etc. */
	double backboneQ,		/* (I) 0.5=N, 0=LN, etc. */
	double volTwkSize,		/* (I) vol tweak size */
	double spotVolMin,		/* (I) minimum spot volatility */
	double tMatMin,			/* (I) Minimum muturity */

	int   volBoundType,		/* (I) vol adj type: 1, -1, 0  */
	int   spotVolRatioFlag,		/* (I) spot vol ratio check flag  */
	double spotVolRatio,		/* (I) Max spot vol ratio  */

        int   finalFlag,		/* (I) TRUE=final, FALSE=cms  */
        TDate maturityDate,		/* (I) final mat date (if final) */
        TDateInterval matInterval,	/* (I) fwd mat (if cms)  */

	int nFailedMax,			/* (I) max # failures */
	int *nFailed,			/* (O) # of failures (or NULL) */
	double *tExpFailed,		/* (O) failed exp times (or NULL) */
	double *tMatFailed,		/* (O) failed mat times (or NULL) */
	double *volMinAdj)		/* (O) needed min vol corrections 
					 *     for failed points */
{
static	char		routine[] = "VnfmCheckSwaptionCalibrationOne";
	int		status = FAILURE;

#define	SP_VOL_MIN	0.25e-2	     /* minimum spot vol */
#define	TOL_VOL_ADJ	1.0e-5	     /* minimum mkt vol adjustment */
#define	TOL_VOL_RATIO	5.0e-2	     /* 1% extra margin for vol ratio check 
				      *	so ratio of 5 is actually 4.95 */

	/* For convenient use of some macros */
	VnfmData	*that = NULL;

	int		vIdx;

	int		numVolDates = 0;
	TDate		*volDates = NULL;
	double		*volMat = NULL;
	int		*volFreq = NULL;
	double		*volRates = NULL;

	TDate		*datesL = NULL;
	double		backboneqL[2];



	/* Input logging */

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "%s: Checking %s %s\n", \
		routine, \
		(finalFlag ? "Final" : "CMS"), \
		(finalFlag ? GtoFormatDate(maturityDate) : \
			GtoFormatDateInterval(&matInterval))); \
	    DrlFPrintf(vnfmFpLog, "\nInput matrix:\n", routine); \
            DrlTSwaptionMatrix2DFpWrite(swoptMat, \
                        vnfmFpLog, TSWAPTION_MATRIX_FMT_STD); \
	);
#endif


	
	/*
	 * Initialized market vol data
	 */
	optimData.volBoundType = volBoundType;
	optimData.tfData    = NULL;
	optimData.volMat    = NULL;
	optimData.volFreq   = NULL;
	optimData.variances = NULL;

	/*  interpolate vol curve from swaption matrix */
	IF_FAILED_DONE( DrlTSwaptionMatrix2DInterpVolCurve(
		swoptMat,
		zcCurve->fBaseDate,
		finalFlag,
		maturityDate,
		matInterval,
		tMatMin,
		FALSE,		/*  TRUE=rebucket, FALSE=interp  */
		&numVolDates,
		&volDates,
		NULL,
		&volMat,
		&volFreq,
		&volRates,
		NULL));


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "%s: Interpolated volatility:\n", routine);\
	    for (vIdx=0; vIdx<=numVolDates-1; vIdx++) { \
		DrlFPrintf(vnfmFpLog, "%3d/%3d  %10s %7.4f  %d   %12.8f%%\n",\
		vIdx, numVolDates, \
		GtoFormatDate(volDates[vIdx]), \
		volMat[vIdx], \
		volFreq[vIdx], \
		volRates[vIdx]*1e2);})
#endif



	/*  Copy dates in LIL array */


	ASSERT_OR_DONE((datesL = NEW_ARRAY(TDate, numVolDates+2)) != NULL);
	datesL[1] = zcCurve->fBaseDate;
	for (vIdx=1; vIdx<=numVolDates; vIdx++)
		datesL[vIdx+1] = volDates[vIdx-1];
	datesL[0] = (long)numVolDates + 1;

	/*  Convert process pwr to bbq */
	backboneqL[0] = 1e0;
	backboneqL[1] = backboneQ;


	/*  Create vnfm */
	IF_FAILED_DONE( VnfmWrapReadSimple(
		&that,		/* optimData.tfData */
		backboneqL,
		nfParamsL,
		datesL,
		zcCurve));

	optimData.tfData = that;

	IF_FAILED_DONE( VnfmComputeCoeff(that));


	/*  logging */
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "%s: VNFM data:\n", routine); \
		VnfmFpWrite(that, vnfmFpLog));
#endif


/* =====================================================================
 *
 *	OOQP IS A METHOD THAT SEEMS TO WORK THE BEST AMONG THE THREE.
 *
 *	THE IMSL QP ROUTINE IS NOT VERY ROBUST WHEN NUMBER
 *	OF OPTIMIZING VARIABLES AND CONSTRAINTS ARE LARGE.
 *
 *	THE BOOTSTRAP METHOD DOES NOT INCLUDE ALGORITHM FOR SPOT VOL RATIO
 *	CHECK, THEREFORE OFTEN FAILED UNDER STRESS CASES.  
 *
 * ====================================================================*/

	/* Check for failed vol calibration in base case */

#if defined(USE_OOQP)
	IF_FAILED_DONE(ComputeAdjustment_OOQP(
				that,	
				swoptMat,
				volDates,	
				volMat,
				volFreq,	
				volRates,
				spotVolMin,	
				tMatMin,	
				volBoundType,	
				spotVolRatioFlag,	
				spotVolRatio,
				nFailed,	
				tExpFailed,
				tMatFailed,
				volMinAdj));

#elif defined (USE_IMSL)	
	IF_FAILED_DONE(ComputeAdjustment_IMSL(
				that,	
				swoptMat,
				volDates,	
				volMat,
				volFreq,	
				volRates,
				spotVolMin,	
				tMatMin,	
				volBoundType,	
				spotVolRatioFlag,	
				spotVolRatio,
				nFailed,	
				tExpFailed,
				tMatFailed,
				volMinAdj));
#else	/* USE_BOOTSTRAP */
	IF_FAILED_DONE(ComputeAdjustment_BootStrp(
				that,	
				swoptMat,
				volDates,	
				volMat,
				volFreq,	
				volRates,
				spotVolMin,	
				tMatMin,	
				volBoundType,	
				spotVolRatioFlag,	
				spotVolRatio,
				nFailed,	
				tExpFailed,
				tMatFailed,
				volMinAdj));
#endif


	/*  We are done */
	if (fabs(volTwkSize)<=TOL_VOL_ADJ) {
		status = SUCCESS;
		goto done;
	}

	/**
	 ** Tweak the vol
	 **/
	for (vIdx=0; vIdx<=numVolDates-1; vIdx++) {

	    /* tweak vol */
	    volRates[vIdx] += volTwkSize;

#if defined(USE_OOQP)
	    IF_FAILED_DONE(ComputeAdjustment_OOQP(
				that,	
				swoptMat,
				volDates,	
				volMat,
				volFreq,	
				volRates,
				spotVolMin,	
				tMatMin,	
				volBoundType,	
				spotVolRatioFlag,	
				spotVolRatio,
				nFailed,	
				tExpFailed,
				tMatFailed,
				volMinAdj));
#elif defined (USE_IMSL)	
	IF_FAILED_DONE(ComputeAdjustment_IMSL(
				that,	
				swoptMat,
				volDates,	
				volMat,
				volFreq,	
				volRates,
				spotVolMin,	
				tMatMin,	
				volBoundType,	
				spotVolRatioFlag,	
				spotVolRatio,
				nFailed,	
				tExpFailed,
				tMatFailed,
				volMinAdj));
#else	/* USE_BOOTSTRAP */
	IF_FAILED_DONE(ComputeAdjustment_BootStrp(
				that,	
				swoptMat,
				volDates,	
				volMat,
				volFreq,	
				volRates,
				spotVolMin,	
				tMatMin,	
				volBoundType,	
				spotVolRatioFlag,	
				spotVolRatio,
				nFailed,	
				tExpFailed,
				tMatFailed,
				volMinAdj));
#endif

	    /* Restore the original vol */
	    volRates[vIdx] -= volTwkSize;
	}


	/*  Sort and remove double items */
	IF_FAILED_DONE( DrlDoubleVectLexSort(
		0.0025,
		TRUE,
		nFailed,
		tExpFailed,
		tMatFailed,
		volMinAdj,
		NULL));


	/*  Select only the max value of volMinAdj in degenerate
	 *  set of tExpFailed and tMatFailed
	 */
	IF_FAILED_DONE( _degenerateAndSelectMaxArray(
		0.0025,   	/* less than ONE day */
        	nFailed,
        	tExpFailed,  
		tMatFailed,
        	volMinAdj));



	/* OK */
	status = SUCCESS;

done:
	if (datesL) 		FREE(datesL);
	if (volDates) 		FREE(volDates);
	if (volMat)   		FREE(volMat);
	if (volFreq)  		FREE(volFreq);
	if (volRates) 		FREE(volRates);

	VnfmFree(optimData.tfData);

	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);


}




/*f--------------------------------------------------------------
 * Swaption volatility check for calibration failures.
 *                                                             
 * <br><br>
 * Performs a test of the swaption matrix for non negative
 * culative variances (i.e. that the spot volatility bootstrapping
 * does not fail) for a selected array of maturities.
 */

DLL_EXPORT(int)
VnfmCheckSwaptionCalibration(
	TCurve *zcCurve,		/*  (I) zero curve */
	TSwaptionMatrix2D *swoptMat,	/*  (I) input matrix */

	double *nfParamsL,		/*  (I) mr, alphas, etc. */
	double backboneQ,		/*  (I) 0.5=N, 0=LN, etc. */
	double volTwkSize,		/* (I) vol tweak size */
	double spotVolMin,		/*  (I) minimum spot volatility */
	double tMatMin,			/*  (I) Minimum muturity */

	int     volBoundType,		/*  (I) vol adj type: 1, -1, 0  */
	int	spotVolRatioFlag,	/*  (I) spot vol ratio check flag  */
	double  spotVolRatio,		/*  (I) Max spot vol ratio  */

	int 	nMat,			/*  (I) # of maturities  */
	double *tMat,			/*  (I) array of maturities */
	long   *finalFlag,		/*  (I) array of final flags */

	double *tExpFailed,		/*  (O) failed exp times (or NULL) */
	double *tMatFailed,		/*  (O) failed mat times (or NULL) */
	double *volMinAdj,		/*  (O) vol adjustment needed  */
	int    *nFailed)		/*  (O) # of failures (or NULL) */
{
static	char	routine[] = "VnfmCheckSwaptionCalibration";
	int	status = FAILURE;

	int		n;
	TDate		maturityDate = -1L;
	TDateInterval	matInterval;
static	double		tMatDef[] = {
			0.083333333333e0, 0.25e0, 0.5e0, 1e0,
			2e0, 3e0, 4e0, 5e0, 6e0, 7e0, 8e0, 9e0, 10e0,
			11e0, 12e0, 13e0, 14e0, 15e0,
			18e0, 20e0, 25e0, 30e0, 35e0, 40e0};
static	long		finalFlagDef[] = {
			0, 0, 0, 0,
			1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1};
static	int		nMatDef = 24;
	int		nFailedMax;

	nFailedMax = nMat * swoptMat->table->matrix->numDim1;

	/*  for PURIFY: wont actually be used  */
	IF_FAILED_DONE( GtoMakeDateInterval(1, 'Q', &matInterval));

	/*  If no maturity specified, use defaults. */
	if ((nMat < 0) || (tMat == NULL) || (finalFlag == NULL)) {
		nMat = nMatDef;
		tMat = tMatDef;
		finalFlag = finalFlagDef;
	}

	*nFailed = 0;


	for (n=0; n<=nMat-1; n++) {

	    if ((int)finalFlag[n] == TRUE) {
		IF_FAILED_DONE( GtoTDateAdvanceYears(
			zcCurve->fBaseDate, tMat[n], &maturityDate));

#ifndef	NO_LOGGING
		GTO_IF_LOGGING(\
		    DrlFPrintf(vnfmFpLog,
			"\nChecking Final %7.4f yrs (%s->%s)\n", \
			tMat[n],  GtoFormatDate(zcCurve->fBaseDate), \
			GtoFormatDate(maturityDate)));
#endif

	    } else {
		IF_FAILED_DONE( GtoYearsToDateInterval(
			tMat[n], &matInterval));

#ifndef	NO_LOGGING
		GTO_IF_LOGGING(\
		    DrlFPrintf(vnfmFpLog, \
			"\nChecking Const %7.4f yrs (%s)\n", \
			tMat[n], \
			GtoFormatDateInterval(&matInterval)));
#endif
	    }



	    if (VnfmCheckSwaptionCalibrationOne(
		zcCurve,
		swoptMat,

		nfParamsL,
		backboneQ,
		volTwkSize,
		spotVolMin,
		tMatMin,

		volBoundType,	
		spotVolRatioFlag,	
		spotVolRatio,

        	(int)finalFlag[n],
        	maturityDate,
        	matInterval,
		nFailedMax,
		nFailed,
		tExpFailed,
		tMatFailed,
		volMinAdj) != SUCCESS)
			goto done;
	}


	
	/*  Sort and remove double items */
	IF_FAILED_DONE( DrlDoubleVectLexSort(
		0.0025,
		TRUE,
		nFailed,
		tExpFailed,
		tMatFailed,
		volMinAdj,
		NULL));


	/*  Select only the max value of volMinAdj in degenerate
	 *  set of tExpFailed and tMatFailed
	 */
	IF_FAILED_DONE( _degenerateAndSelectMaxArray(
		0.0025,   	/* less than ONE day */
        	nFailed,
        	tExpFailed,  
		tMatFailed,
        	volMinAdj));


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "Total Errors: %d.\n", *nFailed));
#endif

	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed\n", routine);
	return(status);
}




/*f--------------------------------------------------------------
 * Swaption volatility smoothing for calibration failures.
 *                                                             
 * <br><br>
 * Performs a swaption vol adjustment if the calibration failed
 * for the spot volatility bootstrapping. For each failed point,
 * increase the vol by minimum required amount until bootstrapping 
 * succeed. Repeat the process until it succeeds for all the maturities.
 */

DLL_EXPORT(int)
VnfmCheckAdjustSwaptionCalibration(
	TCurve *zcCurve,		/*  (I) zero curve */
	TSwaptionMatrix2D *swoptMat,	/*  (B) input matrix */

	double *nfParamsL,		/*  (I) mr, alphas, etc. */
	double backboneQ,		/*  (I) 0.5=N, 0=LN, etc. */
	double volTwkSize,		/* (I) vol tweak size */
	double spotVolMin,		/*  (I) minimum spot volatility */
	double tMatMin,			/*  (I) Minimum muturity */

	int	volAdjFlag,		/*  (I) vol adjustment flag   */
	double	volAdjAmount,		/*  (I) vol adjustment amount */

	int     volBoundType,		/*  (I) vol adj type: 1, -1, 0  */
	int	spotVolRatioFlag,	/*  (I) spot vol ratio check flag  */
	double  spotVolRatio,		/*  (I) Max spot vol ratio  */

	int nMat,			/*  (I) # of maturities  */
	double *tMat,			/*  (I) array of maturities */
	long   *finalFlag,		/*  (I) array of final flags
					 * 	TRUE=final, FALSE=fwd mat */

	double **tExpFailedTotal,	/*  (O) failed exp times (or NULL) */
	double **tMatFailedTotal,	/*  (O) failed mat times (or NULL) */
	int    *nFailedTotal)		/*  (O) Total # of failures/adjust */
{
static	char	routine[] = "VnfmCheckAdjustSwaptionCalibration";
	int	status = FAILURE;

	int             nFailed, i, idx;
 
	int		nAdj, j, jExp, jMat;

        int             MAXNUMCHK = 500;  	/* max number of loops for
                                                vol calibration check  */
	int		nFailedMax;

	double		*tExpFailed = NULL;
	double		*tMatFailed = NULL;
	double		*volMinAdj  = NULL;

	double		*expIdx = NULL;
	double		*matIdx = NULL;
	double		*volAdjIdx = NULL;

	int     k, ie[4], im[4];
        double  w[4], wSQTotal, wVol;

	double	volAdj;

	/* Maximum number of failures: 
	 * every checked vol points + all the tweaked vol scenarios
	 */
	nFailedMax = nMat * swoptMat->table->matrix->numDim1
		     + (swoptMat->table->matrix->numDim1)
			*(swoptMat->table->matrix->numDim1);

	if ((*tExpFailedTotal = NEW_ARRAY(double, nFailedMax*2)) == NULL)
		goto done;
	if ((*tMatFailedTotal = NEW_ARRAY(double, nFailedMax*2)) == NULL)
		goto done;

	if ((tExpFailed = NEW_ARRAY(double, nFailedMax)) == NULL)
		goto done;
	if ((tMatFailed = NEW_ARRAY(double, nFailedMax)) == NULL)
		goto done;
	if ((volMinAdj = NEW_ARRAY(double, nFailedMax)) == NULL)
		goto done;

    	if ((expIdx = NEW_ARRAY(double, nFailedMax*4)) == NULL)
		goto done;
	if ((matIdx = NEW_ARRAY(double, nFailedMax*4)) == NULL)
		goto done;
	if ((volAdjIdx = NEW_ARRAY(double, nFailedMax*4)) == NULL)
		goto done;

        *nFailedTotal = 0;
 
	/* 
	 * Repeat the calibration process, if vol adjustment is required,
         * until either no failed point is found or the max number of
         * iterations is reached.
         */
        for (i=0; i<=MAXNUMCHK-1; i++)
        {
                /* routine call */
                if (VnfmCheckSwaptionCalibration(
                        zcCurve,
                        swoptMat,
 
                        nfParamsL,         /* mr, alphas, etc. */
                        backboneQ,  	   /* process power: 0=N, 1=LN, etc. */
			volTwkSize,        /* vol tweak size */
                        spotVolMin,        /* minimum spot volatility */
                        tMatMin,           /* minimum maturity */
 
			volBoundType,
			spotVolRatioFlag,  /* spot vol ratio check flag  */
			spotVolRatio, 

                        nMat,
                        tMat,
                        finalFlag,         /* final flag */

                        tExpFailed,
                        tMatFailed,
			volMinAdj,
                        &nFailed) != SUCCESS)
                                goto done;
 
                /* Aggregate failed points */
                for (idx=0; idx<=nFailed-1; idx++) 
		{
                    (*tExpFailedTotal)[*nFailedTotal+idx] = tExpFailed[idx];
                    (*tMatFailedTotal)[*nFailedTotal+idx] = tMatFailed[idx];
                }
 
                *nFailedTotal += nFailed;
 
                if (!volAdjFlag  || nFailed == 0 )
                        break;
		else
		{	
#ifndef NO_LOGGING
		    GTO_IF_LOGGING(\
			 GtoErrMsg("\n");
			 GtoErrMsg("tExpFailed\ttMatFailed\tvolMinAdj\n");
			);
#endif

		    nAdj = 0;

	    	    /* Adjust the swaption vols */
		    for (idx=0; idx<=nFailed-1; idx++) 
		    {

		    	/* get interpolation coefficients */
        	    	if (DrlGetInterpCoeffs(swoptMat, 
					       tExpFailed[idx], 
					       tMatFailed[idx], 
					       ie, im, w) != SUCCESS)
                		goto done;

			wSQTotal = 0e0;
			for (k=0; k<=3; k++)
			{
			    wSQTotal += w[k]*w[k];
			}



		    	/* The rebucketing of the adjustment is proportional
			 * to the weight in four nearest neighbors.
		     	 */
		    	for (k=0; k<=3; k++) 
			{
			    if (IS_ALMOST_ZERO(w[k]))
				continue;

			    /* 
			     * Weight in vol increase is proportional to 
			     * w[k]
			     */
			    wVol = w[k] / wSQTotal;

#if (defined(USE_OOQP) || defined(USE_IMSL))
			    volAdj =  volMinAdj[idx];
#else
			    volAdj =  volMinAdj[idx] + volAdjAmount;
#endif


			    expIdx[nAdj] = (double)(ie[k]);
			    matIdx[nAdj] = (double)(im[k]);
			    volAdjIdx[nAdj] = volAdj * wVol;
			    nAdj++;
		    	}
		    } 
		
		    /* sort, remove duplicate items and select max volAdjIdx */
		    IF_FAILED_DONE( DrlDoubleVectLexSort(
					1e-2,
					TRUE,
					&nAdj,
					expIdx,
					matIdx,
					volAdjIdx,
					NULL));


		    IF_FAILED_DONE( _degenerateAndSelectMaxArray(
					0.0025,   /* less than ONE day */
        			   	&nAdj,
					expIdx,
					matIdx,
					volAdjIdx));
		   

		    for (j = 0; j <= nAdj-1; j++)
		    {	
			jExp = (int)(expIdx[j]);
			jMat = (int)(matIdx[j]);
			swoptMat->table->matrix->data[jExp][jMat] +=
			              volAdjIdx[j];	
		    }
			
		    if (i == MAXNUMCHK -1) 
		    {
			GtoErrMsg("%s: Too many vol adjustments to report."
			      " Number of iterations for calibration is %d.\n",
				routine, i+1);	

#ifdef USE_BOOTSTRAP
			if (spotVolRatioFlag && volAdjFlag)
			      GtoErrMsg("\tOld implementation does NOT "
					"include volatility adjustment to "
					"ensure the success of spot volatility"
					" ratio check.\n");
#endif

			goto done;
		    }
    		}


		/*  Sort and remove dublicate items */
		IF_FAILED_DONE( DrlDoubleVectLexSort(
					     0.0025, /* < ONE day */
                                             TRUE,
                                             nFailedTotal,
                                             *tExpFailedTotal,
                                             *tMatFailedTotal,
                                             NULL));

	}	/* for i<=MAXNUMCHK-1 loop */


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "# of adjusted vols: %d,"
				      " # of iterations: %d. \n", 
				*nFailedTotal, i+1));
#endif


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed\n", routine);

	if (expIdx)     FREE(expIdx);
	if (matIdx)     FREE(matIdx);
	if (volAdjIdx)  FREE(volAdjIdx);
	if (tExpFailed) FREE(tExpFailed);
	if (tMatFailed) FREE(tMatFailed);
	if (volMinAdj)  FREE(volMinAdj);

	return(status);
}





static void
ComputeMatrixC(double	**C)	/* (O) result */
{
	/* For convenient use of some macros */
	VnfmData	*that = optimData.tfData;

	double  QB[VNFM_NDIMMAX];

	int	iRow;		/* ith benchmart vol, iRow = 1, ... NDATES-1 */
	int	jCol;		/* jth spot vol , jCol = 1, ... NDATES-1 */

	double	fwdRate, dt, lambda,
		j0;

	int	i, j, rIdx,
		nDim = NDIM;


	/* the first date - base date, is ignored */
	for (iRow=1; iRow<=NDATES-1; iRow++)
	{
	    /* compute the QB coefficients */
	    if (optimData.volMat[iRow] > 0.019230769e0) {
	        if (optimData.volFreq[iRow] > 0) {
		    	VnfmB(that, TT[iRow], optimData.volMat[iRow], 
		          	optimData.volFreq[iRow], &fwdRate, QB);
		} else {
		    	VnfmQ(that, TT[iRow], TT[iRow]+optimData.volMat[iRow], 
			      &fwdRate, QB);
		}
	    } else {
		    VnfmQ(that, TT[iRow], TT[iRow]+optimData.volMat[iRow], 
		          &fwdRate, QB);
	    }


	    /*
	     * Convert to the percentage B/Q coefficients
	     */
	    for (i=0; i<=nDim-1;i++)
		QB[i] /= fwdRate;


	    for (jCol=1; jCol<=NDATES-1; jCol++)
	    {
		/* The matrix is a lower left triangular matrix */
	   	C[iRow-1][jCol-1] = 0.0;

		if (jCol <= iRow ) 
		{	
	    	    dt = TT[jCol] - TT[jCol-1];

		    for (i=0; i<=nDim-1; i++)
		    for (j=i; j<=nDim-1; j++) 
		    {
		    	rIdx = RHOIDX(i, j);
		    	lambda = BETA[i] + BETA[j];
		    	j0 = ALPHA[i] * ALPHA[j] *
			     (i == j ? 1e0 : RHO[rIdx][jCol-1]) *
			     exp(-(TT[iRow]-TT[jCol])*lambda)*L(dt, lambda);

		    	C[iRow-1][jCol-1] += (i == j ? 1e0 : 2e0) 
					     * QB[i] * QB[j] * j0; 
		    }
		}
	    }  	/* for jCol */
	}	/* for iRow */

	return;

}



static int
ComputeIMSLQPCoeff(
	double	*a, 	/* (O) linear constraint coeffs */
	double	*b,	/* (O) linear constraint consts */
	double	*g, 	/* (O) QP linear coeffs */
	double	*h)	/* (O) QP Hessian matrix */
{
static  char    routine[] = "ComputeIMSLQPCoeff";
	int	status = FAILURE;
	
	/* For convenient use of some macros */
	VnfmData	*that = optimData.tfData;

	int		iRow, jCol;

	int		gIdx, hIdx,
			aIdx, bIdx;

	double		norm,
			ratioSquare;

	int		i, j, rIdx,
			nDim = NDIM;

	double		**C = NULL;	/* C coeff matrix */

	int		nMinVolCon, nMktVolCon, nVolRatioCon;

	/* Allocate memory for C matrix: (NDATES-1)*(NDATES-1) */
	ASSERT_OR_DONE((C = DrlDoubleMatrAlloc(
	    				0, 
					NDATES-2, 
					0,
					NDATES-2)) != NULL);

	/* Compute the C coefficients */
	ComputeMatrixC(C);
	

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "\n%s: IMSL QP Coefficients:\n", routine); \
	    DrlFPrintf(vnfmFpLog, "C Matrix:\n"); \
	    DrlMatrixPrint(vnfmFpLog, C, NDATES-1, NDATES-1);\
	);
#endif




	/* Initialize indices */
	aIdx = 0;	/* initialze the index of a */
	bIdx = 0;	/* initialze the index of b */
	gIdx = 0;	/* initialze the index of g */
	hIdx = 0;	/* initialze the index of h */


	/* 
	 * Matrix a is (NDATES-1) * nCon and allocated
	 * in 1-D array. 
	 */

	/* 
	 * 1. Bounds of variables (spot vols^2)
	 */ 
	for(iRow=1; iRow<=NDATES-1; iRow++, bIdx++)
	{
	    for(jCol=1; jCol<=NDATES-1; jCol++, aIdx++)
	    	a[aIdx] = (iRow == jCol ? 1e0 : 0e0);

	    /*
	     * Scale the spot volatility by weighting factors 
	     * and correlation
	     */
	    norm = 0e0;
	    for (i=0; i<=nDim-1; i++)
	    for (j=i; j<=nDim-1; j++) 
	    {
	    	rIdx = RHOIDX(i, j);
	    	norm += ALPHA[i] * ALPHA[j] *
			(i == j ? 1e0 : RHO[rIdx][iRow-1]);
	    }
	    b[bIdx] = optimData.spotVolMin*optimData.spotVolMin / norm * 1.0e2;
	}

	nMinVolCon = bIdx - 1;

	/* 
	 * 2. Mkt Vol bound will have (NDATES - 1) x (NDATES - 1) 
	 * coeffs in matrix a.
	 */
	if(optimData.volBoundType)
	{

	    /* the first date - base date, is ignored */
	    for(iRow=1; iRow<=NDATES-1; iRow++, bIdx++)
	    {
		/* a) Compute the matrix "a" coefficients: nVar x nVar 
		 * Assign the coefficients and 
		 * reverse the sign if adjust vol <= market vol 
		 */
		for (jCol=1; jCol<=NDATES-1; jCol++, aIdx++)
		    a[aIdx] = C[iRow-1][jCol-1] * optimData.volBoundType * 1.0e2;


		/* b) Compute the constraint constants "b", total variance. 
		 * Reverse sign if adjust vol <= market vol 
		 */
		b[bIdx] = optimData.variances[iRow] * optimData.volBoundType * 1.0e4;
	    }
	}

	nMktVolCon = bIdx - 1;

	/* 
	 * 3. Spot vol ratio constraints: 2*(NDATES-1)*(NDATES-2) 
	 */
	if (optimData.spotVolRatioFlag)	
	{
	    ratioSquare = optimData.spotVolRatio * optimData.spotVolRatio
			* (1e0 - TOL_VOL_RATIO) * (1e0 - TOL_VOL_RATIO);

	    /* check ratio between Sigma[jCol]/Sigma[0] */
	    for (iRow=2; iRow<=NDATES-1; iRow++)
	    {
		/* Upper limit */
	    	for (jCol=1; jCol<=NDATES-1; jCol++, aIdx++)
	    	{
		    /* only coeff of sigma[0] and sigma[jCol] are non-zero */
		    a[aIdx] = 0.0;
		    if (jCol==1)    a[aIdx] = ratioSquare;
		    if (jCol==iRow) a[aIdx] = -1.0;
	    	}	

		/* bIdx has been increased by 1 from previous loop */
		b[bIdx++] = 0.0;

		/* Lower limit */
	    	for (jCol=1; jCol<=NDATES-1; jCol++, aIdx++)
	    	{
		    /* only coeff of sigma[0] and sigma[jCol] are non-zero */
		    a[aIdx] = 0.0;
		    if (jCol==1)    a[aIdx] = -1.0;
		    if (jCol==iRow) a[aIdx] = ratioSquare;
	    	}	

		b[bIdx++] = 0.0;
	    }
	}

	nVolRatioCon = bIdx - 1;



#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    aIdx = 0;
	    DrlFPrintf(vnfmFpLog, "\nQP Linear Constraints:\n");\
	    for(iRow=0; iRow<=bIdx-1; iRow++) { \
	        for(jCol=1; jCol<=NDATES-1; jCol++, aIdx++) \
		    DrlFPrintf(vnfmFpLog, "\t%e", a[aIdx]); \
		DrlFPrintf(vnfmFpLog, "   |   "); \
		DrlFPrintf(vnfmFpLog, "%e\n", b[iRow]); \
		/* Separate constraint types */ \
		if (iRow == nMinVolCon || \
		    iRow == nMktVolCon )  \
			DrlFPrintf(vnfmFpLog, "\n"); \
	    } \
	);
#endif



	/*
	 * 4. Hessian matrix and linear coeffs of QP objective function
	 */
	for (iRow=1; iRow<=NDATES-1; iRow++, gIdx++)
	{

	    /* a) QP Hessian matrix */
	    for (jCol=1; jCol<=NDATES-1; jCol++, hIdx++)
	    {
		h[hIdx] = 0.0;	/* initialize */
		
		for (i=MAX(iRow, jCol); i<=NDATES-1; i++)
		    h[hIdx] += C[i-1][iRow-1] * C[i-1][jCol-1] / (TT[i]*TT[i]);

		h[hIdx] *= 2e4;		/* 1/2 xHx, IMSL convention */
	    }

	    /* b) QP linear coeffs */
	    g[gIdx] = 0.0;	/* initialize */
	    for (i=iRow; i<=NDATES-1; i++)
		g[gIdx] += C[i-1][iRow-1] * optimData.variances[i] / (TT[i]*TT[i]);

	    g[gIdx] *= -2.0e6;
	}



#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    gIdx = hIdx = 0; \
	    DrlFPrintf(vnfmFpLog, "\nQP Hessian Matrix:\n"); \
	    for (iRow=1; iRow<=NDATES-1; iRow++) { \
	        for (jCol=1; jCol<=NDATES-1; jCol++, hIdx++) \
		    DrlFPrintf(vnfmFpLog, "\t%e", h[hIdx]);\
		DrlFPrintf(vnfmFpLog, "\n"); \
	    } \
	    DrlFPrintf(vnfmFpLog, "\nQP Linear Coeffs:\n"); \
	    for (iRow=1; iRow<=NDATES-1; iRow++, gIdx++) { \
		DrlFPrintf(vnfmFpLog, "\t%e\n", g[gIdx]); \
	    } \
	);
#endif




	status = SUCCESS;

done:
	if (status != SUCCESS)
                GtoErrMsg("%s: failed\n", routine);

	if (C)  DrlDoubleMatrFree(
				C,
	    			0, 
				NDATES-2, 
				0,
				NDATES-2);

        return(status);
}





static int
ComputeOOQPCoeff(
			/* Indices of QP coefficients */
	int	*irowQ,		/* row index for Q matrix */
	int	*jcolQ,		/* column index for Q matrix */
	int	*irowA,  	/* row index for A matrix */
	int	*jcolA,		/* column index for A matrix */
	int	*irowC,  	/* row index for C matrix */
	int	*jcolC,		/* column index for C matrix */
				/* Coefficients of QP specifications */
	double	*c,		/* linear vector in obj function */
	double	*dQ,		/* Q matrix */
	double	*dA,		/* A matrix */
	double	*dC,		/* C matrix */
			/* Constraint boundaries */
	char	*ixlow,		/* index of non-zero lower bound of X */
	double	*xlow,		/* lower bound of X */
	char	*ixupp,		/* index of non-zero upper bound of X */
	double	*xupp,		/* upper bound of X */
	char	*iclow,		/* index of lower bound of ineq cons */
	double	*clow,		/* lower bound of ineq cons */
	char	*icupp,		/* index of upper bound of ineq cons */
	double	*cupp,		/* upper bound of ineq cons */
	double	*bA)		/* right-hand side vector for eq-cons */
{
static  char    routine[] = "ComputeOOQPCoeff";
	int	status = FAILURE;
	
	/* For convenient use of some macros */
	VnfmData	*that = optimData.tfData;

	int		iRow, jCol;

	int		QIdx, CIdx, rowIdx, iSt, iFlag;

	double		norm,
			ratioSquare;

	int		i, j, rIdx,
			nDim = NDIM;

	double		**C = NULL;	/* C coeff matrix */


	/* Allocate memory for C matrix: (NDATES-1)*(NDATES-1) */
	ASSERT_OR_DONE((C = DrlDoubleMatrAlloc(
	    				0, 
					NDATES-2, 
					0,
					NDATES-2)) != NULL);

	/* 
	 * Compute the C coefficients 
	 * NOTE: this is not the same "C" in OOQP constraint spec 
	 */
	ComputeMatrixC(C);
	

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "\n%s: OOQP Coefficients:\n", routine); \
	    DrlFPrintf(vnfmFpLog, "C Matrix:\n"); \
	    DrlMatrixPrint(vnfmFpLog, C, NDATES-1, NDATES-1);\
	);
#endif

	/* Initialize indices */
	QIdx = 0;	/* initialze the index of Q matrix  */
	CIdx = 0;	/* initialze the index of C matrix  */
	rowIdx = 0;	/* initialze the index of lower bound of C*/

	/* 
	 * 1. Bounds of variables (spot vols^2)
	 */ 
	for(iRow=1; iRow<=NDATES-1; iRow++)
	{
	    /*
	     * Scale the spot volatility by weighting factors 
	     * and correlation
	     */
	    norm = 0e0;
	    for (i=0; i<=nDim-1; i++)
	    for (j=i; j<=nDim-1; j++) 
	    {
	    	rIdx = RHOIDX(i, j);
	    	norm += ALPHA[i] * ALPHA[j] *
			(i == j ? 1e0 : RHO[rIdx][iRow-1]);
	    }

	    /* lower limit */
	    ixlow[iRow-1] = 1;
	    xlow[iRow-1]  = optimData.spotVolMin*optimData.spotVolMin 
			  / norm * 1.0e2;

	    /* no upper limit */
	    ixupp[iRow-1] = 0;
	    xupp[iRow-1] = 0e0;
	}


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "\nOOQP Variable Bounds:\n");\
	    for(iRow=0; iRow<NDATES-1; iRow++) { \
		DrlFPrintf(vnfmFpLog, "x[%d] > %e\n", iRow, xlow[iRow]); \
	    } \
	);
#endif

	/* 
	 * 2. Mkt Vol bound will have nVar x (nVar+1)/2 lower triangle
	 * coeffs in matrix C.
	 */
	if(optimData.volBoundType)
	{

	    /* the first date - base date, is ignored */
	    for(i=1; i<=NDATES-1; i++, rowIdx++)
	    {
		/* a) Compute the matrix "C" coefficients: nVar x nVar 
		 * Assign the coefficients and 
		 * reverse the sign if adjust vol <= market vol 
		 */
		for (j=1; j<=i; j++, CIdx++)
		{
		    irowC[CIdx] = rowIdx;
		    jcolC[CIdx] = j - 1;
		    dC[CIdx]    = C[i-1][j-1] * optimData.volBoundType * 1.0e2;
		}

		/* b) Compute the constraint constants "clow", total variance. 
		 * Reverse sign if adjust vol <= market vol 
		 */
		iclow[rowIdx] = 1;
		clow[rowIdx] = optimData.variances[i] * optimData.volBoundType * 1.0e4;
		/* no upper limit */
		icupp[rowIdx] = 0;
		cupp[rowIdx] = 0e0;
	    }
	}


	/* 
	 * 3. Spot vol ratio constraints: 4*(nVar-1) 
	 */
	if (optimData.spotVolRatioFlag)	
	{
	    ratioSquare = optimData.spotVolRatio * optimData.spotVolRatio
			* (1e0 - TOL_VOL_RATIO) * (1e0 - TOL_VOL_RATIO);

	    /* check ratio between Sigma[jCol]/Sigma[0] */
	    for (i=2; i<=NDATES-1; i++)
	    {
		/* 
		 * Only coeff of sigma[0] and sigma[i-1] are non-zero 
		 */

		/* Upper ratio limit */
		irowC[CIdx] = rowIdx;
		jcolC[CIdx] = 0;
		dC[CIdx]    = ratioSquare;

		CIdx++;

		irowC[CIdx] = rowIdx;
		jcolC[CIdx] = i - 1;
		dC[CIdx]    = -1.0;

		CIdx++;

		/* lower bound only */
		iclow[rowIdx] = 1;
		clow[rowIdx]  = 0.0;
		icupp[rowIdx] = 0;		/* no upper limit */
		cupp[rowIdx] = 0e0;		/* no upper limit */

		rowIdx++;

		/* Lower ratio limit */
		irowC[CIdx] = rowIdx;
		jcolC[CIdx] = 0;
		dC[CIdx]    = -1.0;

		CIdx++;

		irowC[CIdx] = rowIdx;
		jcolC[CIdx] = i - 1;
		dC[CIdx]    = ratioSquare;

		CIdx++;

		/* lower bound only */
		iclow[rowIdx] = 1;
		clow[rowIdx] = 0.0;
		icupp[rowIdx] = 0;  	/* no upper limit */
		cupp[rowIdx] = 0e0;  	/* no upper limit */

		rowIdx++;
	    }
	}


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "\nOOQP Linear Constraints: (%d)\n", rowIdx);\
	    iSt = 0; \
	    iFlag = FALSE; \
	    for(iRow=0; iRow<=rowIdx-1; iRow++) { \
	        for(jCol=0; jCol<NDATES-1; jCol++) { \
		    for (i=iSt; i<= CIdx-1; i++) { \
		    	if (irowC[i] == iRow && jcolC[i] == jCol) { \
			    iFlag = TRUE; \
			    iSt++; \
		    	    DrlFPrintf(vnfmFpLog, "%13.6e  ",dC[i]); \
			    break; \
			} \
		    } \
		    \
		    if (!iFlag) \
		    	DrlFPrintf(vnfmFpLog, "       0       "); \
		    else \
			iFlag = FALSE; \
		} \
		DrlFPrintf(vnfmFpLog, " >   "); \
		DrlFPrintf(vnfmFpLog, "%13.6e\n", clow[iRow]); \
	    } \
	);
#endif


	/*
	 * 4. Hessian matrix and linear coeffs of QP objective function
	 */
	for (iRow=1; iRow<=NDATES-1; iRow++)
	{

	    /* a) QP Hessian matrix Q, lower triangle only */
	    for (jCol=1; jCol<=iRow; jCol++, QIdx++)
	    {
		dQ[QIdx] = 0.0;	/* initialize */
		
		for (i=iRow; i<=NDATES-1; i++)
		    dQ[QIdx] += C[i-1][iRow-1] * C[i-1][jCol-1] / (TT[i]*TT[i]);

		dQ[QIdx] *= 2e4;		/* 1/2 x*H*x */

		irowQ[QIdx] = iRow - 1;
		jcolQ[QIdx] = jCol - 1;

	    }

	    /* b) QP linear coeffs c */
	    c[iRow-1] = 0.0;	/* initialize */
	    for (i=iRow; i<=NDATES-1; i++)
		c[iRow-1] += C[i-1][iRow-1] * optimData.variances[i] 
			     / (TT[i]*TT[i]);

	    c[iRow-1] *= -2.0e6;
	}

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "\nOOQP Hessian Matrix (Lower triangular): %d x %d\n", NDATES-1, NDATES-1); \
	    iSt = 0;
	    iFlag = FALSE; \
	    for(iRow=0; iRow<NDATES-1; iRow++) { \
	        for(jCol=0; jCol<NDATES-1; jCol++) { \
		    for (i=iSt; i<= QIdx-1; i++) { \
		    	if (irowQ[i] == iRow && jcolQ[i] == jCol) { \
			    iSt++; \
			    iFlag = TRUE; \
		    	    DrlFPrintf(vnfmFpLog, "%13.6e  ",dQ[i]); \
			    break; \
			} \
		    } \
		    \
		    if (!iFlag) \
		    	DrlFPrintf(vnfmFpLog, "               "); \
		    else \
			iFlag = FALSE; \
		} \
		DrlFPrintf(vnfmFpLog, "\n"); \
	    } \
	    DrlFPrintf(vnfmFpLog, "\nOOQP Linear Coeffs:\n"); \
	    for (iRow=1; iRow<=NDATES-1; iRow++) { \
		DrlFPrintf(vnfmFpLog, "\t%13.6e\n", c[iRow-1]); \
	    } \
	);
#endif



	status = SUCCESS;

done:
	if (status != SUCCESS)
                GtoErrMsg("%s: failed\n", routine);

	if (C)  DrlDoubleMatrFree(
				C,
	    			0, 
				NDATES-2, 
				0,
				NDATES-2);

        return(status);
}






/*  
 *  Remove the duplicate items in arrayA and arrayB, and select the 
 *  maximum elements in arrayC.
 */   
static int 
_degenerateAndSelectMaxArray(
	double	tolerance,	/* (I) tolerance (Positive)	   */
	int	*numItem,	/* (B) # of elements in each array */
	double	*arrayA,	/* (B) array A (to be degenerated) */
	double	*arrayB,	/* (B) array B (to be degenerated) */
	double	*arrayC)	/* (B) array C (select the max)	   */
{
static	char	routine[] = "_degenerateAndSelectMaxArray";
	int	status = FAILURE;

	int	i, j;
	double	tol;

	j = 0;

	tol = (tolerance > DBL_EPSILON ? tolerance : DBL_EPSILON*1e1);

	if (*numItem <= 1)
	{
	   	status = SUCCESS;
		goto done;
	}

	for (i=0; i<=(*numItem)-2; i++)
	{
	    if( fabs(arrayA[j] - arrayA[i+1]) < tol &&
	        fabs(arrayB[j] - arrayB[i+1]) < tol ) 
	    {
		if (arrayC[j] >= 0e0 && arrayC[i+1] >= 0e0) 
			arrayC[j] = MAX(arrayC[j], arrayC[i+1]);
		else if (arrayC[j] <= 0e0 && arrayC[i+1] <= 0e0)
			arrayC[j] = MIN(arrayC[j], arrayC[i+1]);
		else	/* average (a fudge, not gaurantee to work!!!) */
			arrayC[j] = arrayC[j] + arrayC[i+1];
	    }
	    else
	    {
		j++;
		arrayA[j] = arrayA[i+1];
		arrayB[j] = arrayB[i+1];
		arrayC[j] = arrayC[i+1];
	    }
	}
		
	*numItem = j + 1;
	
	status = SUCCESS;

done:
	if (status != SUCCESS)
                GtoErrMsg("%s: failed\n", routine);
        return(status);
}
