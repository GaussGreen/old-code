/****************************************************************
 * Module:	VNFM
 * Submodule:	CALI
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "drlmem.h"
#include "drllno.h"
#include "drlts.h"			/* DrlTCurveXXX */
#include "drlio.h"			/* DrlFPrintf */

#define	_vnfm_SOURCE

#include "vnfmopca.h"

#undef __DEBUG__

/* This flag enables the positive correlation matrix constraint
 */
#define	POSITIVE_CORR_CONSTRAINT



/*
 * Static data passed to optimization and constraint functions.
 */

typedef	struct	{
	VnfmData	*nfData;	/* model parameters */

	char		optimMethod;
	char		*optimType;

	int		nOptBench;
	TVolBenchmark	*optBench;
	double		*optMid;
	double		*optWei;

	int		nConstBench;
	TVolBenchmark	*constBench;
	double		*constMin;
	double		*constMax;
} TUserData;
static	TUserData	userData;


static	int		nIterF,
			nIterC;


static	int		VnfmFromVect(double *x, VnfmData *m);
static	void		VnfmToVect(double *x, VnfmData *m);

static	int		objectiveFunc(int *MODE, int *N, double *X,
				void *userData,
				double *OBJF, double *OBJGRD);
static	int		constraintFunc(int *NCNLN, int *N, int *MODE, int *NEEDC,
				double *X,
				void *userData,
				double *C, double *cjac1);

#define	NCMAX		32		/* max # of correlation constraints */

static	double		*blcnln1,
			*bucnln1;




			/* vector containng the 2-fact model
			  parameters during optimization.
			  the i-th  optmization variable corresponds
			  to the the k-th variable of the 2-fact model */
static	double		Y[VNFM_NMAX_PARAM];
			/* i-th optim variable -> k-th model variable */
static	int		i2k[VNFM_NMAX_PARAM],
			/* k-th model variable -> i-th optim variable */
			k2i[VNFM_NMAX_PARAM];


/*f-------------------------------------------------------------
 * Parameter optimization using constant spot volatility.
 *                                                             
 * <br><br>
 * Performs the calibration of the two-factor parameters over
 * a series of short-term volatility benchmarks.
 * <br>
 * <br>[that] On entry, parameter data structure with an already
 * initialized timeline. On exit, contains the optimal parameters.
 * <br>[optimType] Optimization type:
 *      "F" for the optimization of the weighted quadratic difference 
 *      between market and benchmark volatilities,
 *      "R" for the optimization of the weighted quadratic difference 
 *      between the market and benchmark ratio of volatilities of two
 *	consecutive benchmarks (the number of bencharks must be even).
 * <br>[nOptBench] number of optimized benchmarks,
 * <br>[optBench] array of optimized benchmarks,
 * <br>[optMid] market value of benchmarks,
 * <br>[optMod] On exit, model value of benchmarks,
 * <br>[nConstBench] number of constrained benchmarks,
 * <br>[constBench] array of constrained benchmarks,
 * <br>[constMin] array of low  constraints,
 * <br>[constMax] array of high constraints,
 * <br>[constMod] On exit, array of model values of constraints,
 * <br>[paOptFlag] array of flags (TRUE=optimize),
 * <br>[paOpt] array of optimal paramters (O),
 * <br>[paMin] array of low  constraints on parameters,
 * <br>[paMax] array high constraints on parameters,
 * <br>[paMid] array of initial guess (or value if no optimization),
 * <br>
 * Returns 0 iff successful.
 */

int
VnfmCalibParamShortTerm(
	VnfmData *that,		/* (I/O) model parameters */
	char *optimType,	/* (I) "F"lat, "R"atio */

	int nOptBench,		/* (I) # of optimized benchmarks */
	TVolBenchmark *optBench,/* (I) array of optimized benchmarks */
	double *optWei,		/* (O) weight of benchmarks in optim */
	double *optMid,		/* (I) market value of benchmarks */
	double *optMod,		/* (O) model  value of benchmarks */

	int nConstBench,	/* (I) # of constrained benchmarks */
	TVolBenchmark *constBench,/* (I) array of constrained benchmarks */
	double *constMin,	/* (I) array of low  constraints */
	double *constMax,	/* (I) array of high constraints */
	double *constMod,	/* (O) model values of constraints */
				/* Optim Params array:
				   beta1,..,betaN,alpha1,..,alphaN,rho11,..*/
	long *paOptFlag,	/* (I) array of flags (TRUE=optimize) */
	double *paOpt,		/* (O) optimal paramters */
	double *paMin,		/* (I) low  constraints on parameters */
	double *paMax,		/* (I) high constraints on parameters */
	double *paMid)		/* (I) initial guess (or value if no opt)*/
{
static	char	routine[] = "VnfmCalibParamShortTerm";
	int	status = FAILURE;
	int	NVAR ;			/* number of variables */
	double	*blsc=NULL,
		*busc=NULL;		/* constraints on variables */

	int	NCLIN;			/* # linear constraints */
	double	*blclin=NULL,
		*buclin=NULL,
		**a=NULL;

	int	NCNLN;
	double	*blcnln=NULL,
		*bucnln=NULL;
	int	ITERMAX;
	double	*C=NULL;
	double	OBJF;
	double	*X=NULL;

	int	nDim = NDIM,
		nPa,
		idxF, idxT,
		i, k;


	long		optimMethod;	/* optimizatiom method */
	TLNOParams	optimParams;	/* optimization parameters */



#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: start\n", routine));
#endif
	/*
	 * Copy pointers in static variables (visible outside the routine)
	 */
	userData.nfData = that;
	userData.optimType = optimType;

	userData.nOptBench = nOptBench;
	userData.optBench = optBench;
	userData.optWei = optWei;
	userData.optMid = optMid;
	userData.nConstBench = nConstBench;
	userData.constBench = constBench;
	userData.constMax = constMax;
	userData.constMax = constMax;

	/*
	 * Initialization of the structure
	 */
	nPa = that->fNf * (that->fNf+3) / 2;

	/* set spot vols to 1 (weights are optimized) */
	for (idxF=0; idxF<=that->fNf-1; idxF++) {
		for (idxT=0; idxT<=that->fNDates-1; idxT++) {
			that->fSigma[idxF][idxT] = 1e0;
		}
	}

	/* Fill param struct with mid values and compute coeff
	 */
	IF_FAILED_DONE( VnfmFromVect(paMid, that));

	/* compute forward rates */
	if (VnfmComputeCoeff(that) != 0) goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlTVolBenchmarkFpWrite(optBench, nOptBench, vnfmFpLog, \
		optMid, NULL);\
	DrlTVolBenchmarkFpWrite(constBench, nConstBench, vnfmFpLog,\
		constMin, constMax, NULL);\
	VnfmFpWrite(that, vnfmFpLog);\
	DrlTCurveFpWrite(that->fZcCurve, vnfmFpLog, DRL_TCURVE_FMT_PERCENT));
#endif


	/*
	 * Check if anything to do
	 */
	for (k=0, NVAR = 0; k<=nPa-1; k++)
			if (paOptFlag[k] != 0L) NVAR++;

	/* no optimization */
	if (NVAR == 0) {
		/*
		 * Nothing to optimize
		 */
		for (k=0; k<=nPa-1; k++) Y[k] = paMid[k];
		VnfmToVect(Y, that);
		VnfmToVect(paOpt, that);
		VnfmComputeJ(that);

		/* Compute benchmark values */
		VnfmComputeCoeff(that);
		for(i=0; i<=nOptBench-1;i++) {
		    if (VnfmTVolBenchmarkValue(that, &optBench[i],
		    	&optMod[i]) != SUCCESS) goto done;
		}
		for(i=0; i<=nConstBench-1;i++) {
		    if (VnfmTVolBenchmarkValue(that, &constBench[i],
		    	&optMod[i]) != SUCCESS) goto done;
		}
		return(0);
	}


	/*-----------------------------------------------------
	 * Fill the default values with mid values
	 *-----------------------------------------------------*/

	for (k=0; k<=nPa-1; k++) Y[k] = paMid[k];
	VnfmToVect(Y, that);


	/*-----------------------------------------------------
	 * Bounds on the state variables
	 *-----------------------------------------------------*/

	/* compute number of optimized variables [0..n-1] */
	for (k=0, NVAR = 0; k<=nPa-1; k++)
		if (paOptFlag[k] != 0) NVAR++;

	/* allocate memory for bounds */
	if ((X    = DrlDoubleVectAlloc(0, NVAR-1)) == NULL) goto done;
	if ((blsc = DrlDoubleVectAlloc(0, NVAR-1)) == NULL) goto done;
	if ((busc = DrlDoubleVectAlloc(0, NVAR-1)) == NULL) goto done;


	/* fill the arrays */
	i = 0;
	for (k=0; k<=nPa-1; k++) {
		k2i[k] = -1;

		if (paOptFlag[k] != 0L) {
			/*
			 * the i-th  optmization variable corresponds
			 * to the the k-th variable of the 2-fact model
			 */
			i2k[i] = k;
			k2i[k] = i;

			/* set bounds */
			blsc[i] = paMin[k];
			busc[i] = paMax[k];

			/* one more opt variable */
			i++;
		}
	}

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "\tNVAR=%d\n", NVAR);\
	for (i=0; i<=NVAR-1; i++)\
		DrlFPrintf(vnfmFpLog, "\t\ti2k[%2d] = %d\n", i, i2k[i]);\
	DrlFPrintf(vnfmFpLog, "\tNParam=%d\n", nPa);\
	for(k=0; k<=nPa-1; k++)\
		DrlFPrintf(vnfmFpLog, "\t\tk2i[%2d] = %d\n", k, k2i[k]));
#endif



	/*-----------------------------------------------------
	 * LINEAR CONSTRAINTS
	 *-----------------------------------------------------*/

	/* check if beta1 and beta2 are BOTH optimized */
	NCLIN = 0;
#ifdef	_SKIP
	if ((k2i[0] >= 0) && (k2i[1] >= 0)) {
		NCLIN = 1 ;
		if ((blclin = DrlDoubleVectAlloc(0, NCLIN-1)) == NULL)
			goto done;
		if ((buclin = DrlDoubleVectAlloc(0, NCLIN-1)) == NULL)
			goto done;
		if ((a = DrlDoubleMatrAlloc(0, NCLIN-1, 0, NVAR)) == NULL)
			goto done;

		/* set the linear constraints */
		/* beta1 > beta2 */
		blclin[0] = 0e0;
		buclin[0] = 1e12;
		for (i=0; i<=NVAR-1; i++) a[0][i] = 0. ;
		a[0][k2i[0]] =  1. ;
		a[0][k2i[1]] = -1. ;
	}
#endif

	/*-----------------------------------------------------
	 * NONLINEAR CONSTRAINTS
	 *-----------------------------------------------------*/

	/*
	 * Nonlinear Constraints
	 */
#ifdef	POSITIVE_CORR_CONSTRAINT
	NCNLN = nConstBench + MAX(nDim-2, 0);
#else
	NCNLN = nConstBench;
#endif
	if (NCNLN >= 1) {
		if ((blcnln = DrlDoubleVectAlloc(0, NCNLN-1)) == NULL)
			goto done;
		if ((bucnln = DrlDoubleVectAlloc(0, NCNLN-1)) == NULL)
			goto done;
		if ((C =      DrlDoubleVectAlloc(0, NCNLN-1)) == NULL)
			goto done;

		/* set the nonlinear constraints */
		for(k=0; k<=nConstBench-1; k++) {
			blcnln[k] = constMin[k] ;
			bucnln[k] = constMax[k] ;
		}
#ifdef	POSITIVE_CORR_CONSTRAINT
		for(k=nConstBench; k<nConstBench+MAX(nDim-2,0); k++) {
		    blcnln[k] = 1e-1;
		    bucnln[k] = 1e12;
		}
#endif
	}

	/*
	 * Initial Value
	 */
	OBJF = 0e0;

	VnfmToVect(Y, that) ;
	for (i=0; i<=NVAR-1; i++)
		X[i] = Y[i2k[i]];


	/*
	 * static variables used in func1
	 */
	nIterF = nIterC = 0 ;
	blcnln1 = blcnln;
	bucnln1 = bucnln;

	/*
	 * Do The Minimization
	 */
	TLNOParamsSetDefault(&optimParams);
	optimMethod = DRL_LNO_COMPBOX;
	optimMethod = DRL_LNO_IMSL_NLN;
	optimParams.ftol = 1e-6;
	optimParams.xtol = 1e-6;

	userData.optimMethod = optimMethod;

	ITERMAX = 1000;


	if (DrlNLinProg(
		NVAR, blsc, busc,
		NCLIN, blclin, buclin, a,
		NCNLN, blcnln, bucnln, constraintFunc,
		objectiveFunc, ITERMAX, C, &OBJF, X,
		(void*)NULL,		/* user data */
		optimMethod,
		&optimParams)
	    != 0) {
		GtoErrMsg("%s: optimization failed.\n", routine);
		goto done;
	}


	/*
	 * store the optimal parameters
	 */
	
	for (i=0; i<=NVAR-1; i++) Y[i2k[i]] = X[i];
	IF_FAILED_DONE(VnfmFromVect(Y, that));

	/* Rescale alphas */
	if (optimType[0] == 'R') {
		IF_FAILED_DONE( VnfmTVolBenchmarkRescaleAlpha(
        		that,
			nOptBench,
			optWei,
			optMid,
			optBench));
	}



	VnfmToVect(paOpt, that);

	/*
	 * Compute benchmark values
	 */
	VnfmComputeCoeff(that);
	for(i=0; i<=nOptBench-1;i++) {
	    if (VnfmTVolBenchmarkValue(that, &optBench[i],
	    	&optMod[i]) != SUCCESS) goto done;
	}
	for(i=0; i<=nConstBench-1;i++) {
	    if (VnfmTVolBenchmarkValue(that, &constBench[i],
	    	&constMod[i]) != SUCCESS) goto done;
	}


	/*
	 * exit the routine
	 */

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: nIterF=%4d \n",\
		routine, ++nIterF));
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "Optimal:\n");\
		for (k=0; k<=NVAR-1; k++)
			DrlFPrintf(vnfmFpLog, "\t%lf", X[k]);\
		DrlFPrintf(vnfmFpLog, "\n"));
	GTO_IF_LOGGING(VnfmFpWrite(that, vnfmFpLog));
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
		"\tbenchmark volatilities\n\tmkt\tmod\tdiff\n");\
		for(i=0; i<=nOptBench-1;i++) {\
			DrlFPrintf(vnfmFpLog, "\t%.2f\t%.2f\t%.2f\n",\
			optMid[i]*1e2, optMod[i]*1e2,\
			(optMid[i]-optMod[i])*1e2);\
		}\
	);
#endif


	/* made it through OK */
	status = SUCCESS;
done:
	/* free alloctaed memory */
	DrlDoubleVectFree(X,    0, NVAR-1);
	DrlDoubleVectFree(blsc, 0, NVAR-1);
	DrlDoubleVectFree(busc, 0, NVAR-1);

	DrlDoubleVectFree(blclin, 0, NCLIN-1);
	DrlDoubleVectFree(buclin, 0, NCLIN-1);
	DrlDoubleMatrFree(a, 0, NCLIN-1, 0, NVAR);

	DrlDoubleVectFree(blcnln, 0, NCNLN-1);
	DrlDoubleVectFree(bucnln, 0, NCNLN-1);
	DrlDoubleVectFree(C,      0, NCNLN-1);

	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}

/*--------------------------------------------------------------
 */

static	int
VnfmFromVect(double *y, VnfmData *m)
{
	int	status = FAILURE;
	int	idxF1, idxF2, idxT, idxR,
		nDim = m->fNf;

	for (idxF1=0; idxF1<=nDim-1; idxF1++) {
		m->fBeta[idxF1]  = y[idxF1];
		m->fAlpha[idxF1] = y[idxF1+nDim];
	}


	for (idxF1=0;       idxF1<=nDim-1; idxF1++)
	for (idxF2=idxF1+1; idxF2<=nDim-1; idxF2++) {
		idxR = RHOIDX(idxF1, idxF2);
		for (idxT=0; idxT<=m->fNDates-1; idxT++) {
			m->fRho[idxR][idxT] = y[idxR+nDim+nDim];
		}
	}

	status = SUCCESS;
	return(status);
}

static	void
VnfmToVect(double *y, VnfmData *m)
{
	int	idxF1, idxF2, idxR,
		nDim = m->fNf;

	for (idxF1=0; idxF1<=nDim-1; idxF1++) {
		y[idxF1]      = m->fBeta[idxF1];
		y[idxF1+nDim] = m->fAlpha[idxF1];
	}


	for (idxF1=0;       idxF1<=nDim-1; idxF1++)
	for (idxF2=idxF1+1; idxF2<=nDim-1; idxF2++) {
		idxR = RHOIDX(idxF1, idxF2);
		y[idxR+nDim+nDim] =  m->fRho[idxR][0];
	}
}


/*ARGSUSED*/
/*--------------------------------------------------------------
 * Function:	Objective function for the minimization routine
 */

static	int
objectiveFunc(
	int *MODE,
	int *N,
	double *X,
	void *data,
	double *OBJF,
	double *OBJGRD)
{
	int	status = FAILURE;

	int	i, numRatios;
	double	vMod0, vMkt0,
		vMod1, vMkt1,
		objVal;

	/*
	 * 1: new values of the model parameters
	 */
	for (i=0; i<=*N-1; i++) Y[i2k[i]] = X[i];
	VnfmFromVect(Y, userData.nfData);
	VnfmComputeJ(userData.nfData);


#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING(VnfmFpWrite(userData.nfData, vnfmFpLog));
#endif

	/*
	 * 2: compute the objective function
	 */
	switch (userData.optimType[0]) {
	case 'F':
	    /* Difference between actual and model vol squared
	     */
	    objVal = 0e0;
	    for (i=0; i<=userData.nOptBench-1; i++) {
	        if (!IS_ALMOST_ZERO(userData.optWei[i])) {
		    IF_FAILED_DONE( VnfmTVolBenchmarkValue(
			    userData.nfData,
			    &userData.optBench[i],
			    &vMod0));
		    vMkt0 = userData.optMid[i];
		    objVal += userData.optWei[i]*(vMod0 - vMkt0)
			* (userData.optimMethod == DRL_LNO_COMPBOX ?
				1e0 : vMod0 - vMkt0);
	        }
	    }
	    break;

	case 'R':
	    /* Ratio of two consecutive volatilities
	     */
	    numRatios = (int) userData.nOptBench / 2;
	    objVal = 0e0;
	    for (i=0; i<numRatios; i++) {
	        if (!IS_ALMOST_ZERO(userData.optWei[2*i])) {
		    IF_FAILED_DONE( VnfmTVolBenchmarkValue(
			    userData.nfData,
			    &userData.optBench[2*i],
			    &vMod0));
		    vMkt0 = userData.optMid[2*i];

		    IF_FAILED_DONE( VnfmTVolBenchmarkValue(
			    userData.nfData,
			    &userData.optBench[2*i+1],
			    &vMod1));
		    vMkt1 = userData.optMid[2*i+1];


		    objVal += userData.optWei[2*i]*
			((vMod1 / vMod0) / (vMkt1 / vMkt0) - 1e0) *
			((vMod1 / vMod0) / (vMkt1 / vMkt0) - 1e0);
	        }
	    }

	    break;
	default:
	    goto done;
	}

	++nIterF;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "O[%3d] ", nIterF);\
		for (i=0; i<=*N-1; i++) \
			DrlFPrintf(vnfmFpLog, "%6.4f ", X[i]);\
		DrlFPrintf(vnfmFpLog, "o=%.4e\n", objVal);\
	);
#endif

	*OBJF = objVal;

	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("objectiveFunc: failed.\n");
		*OBJF = objVal;
	}
	return(status);
}


/*ARGSUSED*/
/*--------------------------------------------------------------
 * Function:	Nonlinear constraints
 */

static	int
constraintFunc(
	int *NCNLN,
	int *N,
	int *MODE,
	int *NEEDC,
	double *X,
	void *data,
	double *C,
	double *cjac1)
{
	int	k,		/* constraint index */
		k0,		/* offset for 1st set of constraints */
		idx;
	int	i;


	k0 = userData.nConstBench;


	for(k=0; k<=*NCNLN-1; k++) {
	    if (NEEDC[k] > 0) {
		/* constraint need to be computed */

		/* new values for model parameters */
		for (i=0; i<=*N-1; i++) Y[i2k[i]] = X[i];
		VnfmFromVect(Y, userData.nfData);
		VnfmComputeJ(userData.nfData);

#ifdef	POSITIVE_CORR_CONSTRAINT
		if (k <= k0-1) {
			idx = k;
			if (VnfmTVolBenchmarkValue(
				userData.nfData,
				&userData.constBench[k],
				&C[k]) != SUCCESS)
				    C[k] = 1e32;

		} else {
			/* Corr mat positivity constraint
			 */
			idx = k - k0;

			if (VnfmCorrMinorDeterm(
				userData.nfData, 0,
				idx+2, 			/* order */
				&C[k]) != SUCCESS) 	/* determinant */
			    C[k] = -1e32;
		}

#else
		/* Compute model correlation */
		if (VnfmTVolBenchmarkValue(
			userData.nfData,
			&userData.constBench[k],
			&C[k]) != SUCCESS) {
			C[k] = 1e32;
		}
#endif



#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "C[%3d] ", nIterC++);\
		for (i=0; i<=*N-1; i++) \
			DrlFPrintf(vnfmFpLog, "%6.4f ", X[i]);\
		DrlFPrintf(vnfmFpLog, "C[%1d]=%8.4f %c\n",\
			k, C[k],\
			(C[k] < blcnln1[k] - 1e-2 ? '-' :\
				(C[k] >= bucnln1[k] + 1e-2 ? '+' : 's')));\
	);
#endif

	    }
	}
	return(SUCCESS);
}





