/************************************************************************
 * Module:	DRL
 * Submodule:	RAND
 * Function:	Random Number Generation and Simulation
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <string.h>
#include <math.h>

#include "drlio.h"
#include "drlstr.h"
#include "drllineq.h"
#include "drlmem.h"
#include "drlrand.h"
#include "drlgpars.h"
#include "drlvtype.h"
#include "drlproc.h"

#include "drlrand.h"
#include "drlrandl.h"


/*#define	__DEBUG__*/

#undef	SQR
#define	SQR(x)	((x)*(x))

/*----------------------------------------------------------------------
 * Simulation of multivariate normal distribution
 * Parsing formula
 */

DLL_EXPORT(int)
DrlFormMonteCarlo(
	int	nDim,	 	/* (I) number of dimensions */
	char	what,		/* (I) 'n' form normal, 'l' for log-normal */
	char	mcType,		/* (I) 'm' for MC, 's' for Sobol */
	double	*p,		/* (I) expectation vector[0..nDim-1] */
	double	texp,		/* (I) time factor */
	double	*vol,		/* (I) stddev vector[0..nDim-1] */
	double	**corrMat,	/* (I) corr matrix[0..nDim-1][0..nDim-1]*/
	char	*form,		/* (I) formula */
	long	nSample,	/* (I) number of deviates */
	int	doStatFlag,	/* (I) if TRUE, performs statistics */
	double	*retVal)	/* (O) retVal[0] = expected, retVal[1] = sdev */
{
static	char	routine[] = "DrlFormMonteCarlo";
	int	status = FAILURE;

#define	MAXDIM		20		/* maximum dimension */
#define	SUB_SAMPLE	4		/* num of subsets */

	long	seed, k, kMax;
	double	gausDev[MAXDIM],
		y[MAXDIM],
		val,
		sVal,			/* expect. */
		subVal[SUB_SAMPLE+1],	/* sub sample expect. */
		vVal,			/* variance */
#ifdef	__DEBUG__
		testVal,
#endif
		eVec[MAXDIM],
		*cMat[MAXDIM],
		cMat1[MAXDIM][MAXDIM];
	int	i, j,
		cntVal,
		subIdxMax=SUB_SAMPLE,
		subIdx;

	/*
	 * Init
	 */
	if (nDim > MAXDIM-1) {
		GtoErrMsg("%s: dim too large %d\n", routine, nDim);
		return(2);
	}


	/* select bet MC and QMC */
	switch (mcType) {
	case 's':
		seed = -1L;
		break;
	default:
		seed = 0L;
		break;
	}

	/*
	 * Initialize the Random Number Generator
	 */
	if (DrlMultiNormSimulInit(seed, nDim, corrMat) != SUCCESS) {
		GtoErrMsg("%s: Init multi-normal simulation failed.\n",
			routine);
		GtoErrMsg("%s: Input correlation matrix:\n", routine);
		DrlFPrintDoubleMatr(NULL, NULL, corrMat, nDim, nDim);
		goto done;
	}

	if (doStatFlag != 0) {
		/*cMat = MatrixAlloc(nDim, nDim);
		eVec = DoubleVectAlloc(nDim);*/
		for(i=0; i<=nDim-1; i++)
			cMat[i] = &cMat1[i][0];

		for(i=0; i<=nDim-1; i++)
			eVec[i] = 0e0;
		for(i=0; i<=nDim-1; i++)
		for(j=0; j<=nDim-1; j++)
			cMat[i][j] = 0e0;
	}

	/*
	 * MC loop
	 */

	/* sub sampling for stddev estimation */
	kMax =  (nSample / subIdxMax);
	nSample = kMax * subIdxMax;

	sVal = 0e0;
	vVal = 0e0;
	cntVal = 0;
#ifdef	__DEBUG__
	testVal = 0e0;
#endif
	for (subIdx=0; subIdx<=subIdxMax-1; subIdx++) {
		subVal[subIdx] = 0.;
	}

	/*
	 * MC loop
	 */
	for (subIdx=0; subIdx<=subIdxMax-1; subIdx++) {
	for (k=0; k<=kMax-1; k++) {

		/* get the deviates */
		DrlMultiNormSimulGet(gausDev);

		if (doStatFlag != 0) {
			for(i=0; i<=nDim-1; i++) {
				eVec[i] += gausDev[i];
			}
			for(i=0; i<=nDim-1; i++)
			for(j=0; j<=nDim-1; j++) {
				cMat[i][j] += gausDev[i]*gausDev[j];
			}
		}

		/* generate the deviates */
		switch (what) {
		case 'l':
		case 'L':
			/* lognormal deviates */
			for(i=0; i<=nDim-1; i++) {
				y[i] = p[i]*exp(-0.5*vol[i]*vol[i]*texp
					+ vol[i]*sqrt(texp)*gausDev[i]);
			}
			break;
		case 'n':
		case 'N':
			/* normal deviates */
			for(i=0; i<=nDim-1; i++) {
				y[i] = p[i] + vol[i]*sqrt(texp)*gausDev[i];
			}
			break;
		default:
			GtoErrMsg("%s: bad type `%c'\n", routine, what);
			return(2);
		}




		/* compute the payoff */
		if (DrlGParEval(nDim, y, form, &val) != SUCCESS) {
			goto done;
		}

		cntVal++;
		sVal += val;
		vVal += SQR(val);
		subVal[subIdx] += val;


#ifdef	__DEBUG__
		testVal +=  MAX(y[0]-((y[1]>2.8) && (y[1]<3.2) ? 3.8 : 4.2),0);
#endif



	}
	}

	/*
	 *
	 */
	DrlMultiNormSimulDispose();

	sVal /= (double) nSample;
	vVal /= (double) nSample;
	for (subIdx=0; subIdx<=subIdxMax-1; subIdx++) {
		subVal[subIdx] /= ((double) kMax);
	}

	if (cntVal != nSample) {
		GtoErrMsg("%s: cntVal (%ld) != nSample (%ld).\n", 
			routine, ((long) cntVal), ((long) nSample));
		PROGRAM_BUG();
		goto done;
	}


#ifdef	__DEBUG__
	testVal /= (double) nSample;
#endif


#ifdef	__DEBUG__
	fprintf(stdout, "M1 = %lf\n", sVal);
	fprintf(stdout, "M2 = %lf\n", vVal);
	for (subIdx=0; subIdx<=subIdxMax-1; subIdx++) {
		fprintf(stdout, "\tM1[%d] = %lf\n", subIdx, subVal[subIdx]);
	}
	fflush(stdout);
#endif

	/* statistics */
	if (doStatFlag != 0) {
		FILE	*fpLog = NULL;
		for(i=0; i<=nDim-1; i++)
			eVec[i] /= (double) nSample;
		for(i=0; i<=nDim-1; i++)
		for(j=0; j<=nDim-1; j++)
			cMat[i][j] /= (double) nSample;

#define	FMT	"\t%10.6f"
		DrlFPrintf(fpLog,"----- Sample Expectation -----\n");
		DrlFPrintDoubleVect(fpLog, FMT, eVec, nDim);
		DrlFPrintf(fpLog, "----- Input Covariance Matrix -----\n");
		DrlFPrintDoubleMatr(fpLog, FMT, corrMat, nDim, nDim);
		DrlFPrintf(fpLog, "----- Sample Covariance Matrix -----\n");
		DrlFPrintDoubleMatr(fpLog, FMT, cMat, nDim, nDim);
		for(i=0; i<=nDim-1; i++)
		for(j=0; j<=nDim-1; j++)
			cMat[i][j] /= sqrt(cMat[i][i] * cMat[j][j]);
		DrlFPrintf(fpLog, "----- Sample Correlation Matrix -----\n");
		DrlFPrintDoubleMatr(fpLog, FMT, cMat, nDim, nDim);

		DrlFPrintf(fpLog, "----- Sample ------------------------\n");
		DrlFPrintf(fpLog, "NDev = %d\n", nSample);
		DrlFPrintf(fpLog, "Form = `%s'\n", form);
		DrlFPrintf(fpLog, "Avg  = %10.8f\n", sVal);
		DrlFPrintf(fpLog, "Std  = %10.8f\n",
				sqrt((vVal-sVal*sVal)/(double) nSample));
	}

#ifdef	__DEBUG__
	fprintf(stdout, "testVal = %lf\n", testVal);
#endif

	/*
	 * End
	 */
	status = SUCCESS;
done:
	if (doStatFlag != 0) {
		/*MatrixFree(cMat, nDim, nDim);
		DoubleVectorFree(eVec);*/
	}

	retVal[0] = sVal;

	/* stddev estimation */
	/*retVal[1] = sqrt((vVal-sVal*sVal)/(double) nSample);*/
	retVal[1] = 0.;
	for (subIdx=0; subIdx<=subIdxMax-1; subIdx++) {
		retVal[1] += SQR(subVal[subIdx] - sVal);
	}
	retVal[1] = sqrt(retVal[1]/((double) subIdxMax));

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}

	return(status);
}




/*f---------------------------------------------------------------------
 * Wrapper for the formula parser.
 */

DLL_EXPORT(int)
DrlFormMonteCarloL(
	long *nDimL,		/* 0 'L' (I) number of dimensions */
	char *whatL,		/* 1 'C' (I) 'n' normal, 'l' log-normal */
	double *pL,		/* 2 'F' (I) expectation vector[0..nDim-1] */
	double *texpL,		/* 3 'F' (I) time factor */
	double *volL,		/* 4 'F' (I) stddev vector[0..nDim-1] */
	double *corrMatL,	/* 5 'F' (I) corr matrix[0..nDim-1][0..nDim-1]*/
	char *formL,		/* 6 'C' (I) formula */
	long *nSampleL,		/* 7 'L' (I) number of deviates */
	long *doStatFlagL,	/* 8 'L' (I) if TRUE, performs statistics */
	double *retValL)	/*   'F' (O) [1]=value, [2]=stdev. */
{
static	char	routine[] = "DrlFormMonteCarloL";
	int	status = FAILURE;
#define	MAXDIM	20
	int	nDim,
		n1, n2, n3;
	char	what,
		mcType,
		buf[64];
	double	texp,		/* time to expiration */
		*p=NULL,	/* vector of exp values */
		*vol=NULL,	/* vector of vlatilities */
		**corrMat=NULL;	/* correlation matrix */
	char	form[128];	/* formula */
	long	nSample;
	int	doStatFlag;
	double	retVal[5];

	/*
	 *
	 */
	WRAP_CHECK_VECTOR(nDimL);
	WRAP_CHECK_VECTOR(whatL);
	WRAP_CHECK_VECTOR(pL);
	WRAP_CHECK_VECTOR(texpL);
	WRAP_CHECK_VECTOR(volL);
	WRAP_CHECK_VECTOR(corrMatL);
	WRAP_CHECK_VECTOR(formL);
	WRAP_CHECK_VECTOR(nSampleL);
	WRAP_CHECK_VECTOR(doStatFlagL);
	WRAP_CHECK_VECTOR(retValL);

	/*
	 *
	 */
#ifdef	__DEBUG__
	fprintf(stdout, "#\n# MC_FORMULA\n#\n"); fflush(stdout);
#endif /* __DEBUG__ */


	if (DrlLilStructGet(0,
	    DRL_LILVAR_L,
		"nDim", DRL_INT_L, (void*) nDimL, (void*) &nDim,
	    0) != SUCCESS) goto done;

	if (DrlLilStructGet(0,
	    DRL_LILVAR_L,
		"what", DRL_CHAR_L, (void*) whatL, (void*) &what,
	    DRL_LILVAR_L,
		"texp", DRL_FLOAT_L, (void*) texpL, (void*) &texp,
	    DRL_LILVAR_L,
		"nSample", DRL_INT_L, (void*) nSampleL, (void*) &nSample,
	    DRL_LILVAR_L,
		"doStatFlag", DRL_INT_L,(void*)doStatFlagL, (void*)&doStatFlag,
	    DRL_LILVECT_L, &n1, 1, MAXDIM,
		"p", DRL_FLOAT_L, TRUE, (void*) pL, (void**) &p,
	    DRL_LILVECT_L, &n2, 1, MAXDIM,
		"vol", DRL_FLOAT_L, TRUE, (void*) volL, (void**) &vol,
	    DRL_LILMATR_L, nDim, nDim,
		"corrMat", DRL_FLOAT_L, TRUE,(void*)corrMatL, (void*)&corrMat,
	    0) != SUCCESS) goto done;


	strcpy(form, &formL[WRAP_STR_IDX(1)]);
#ifdef	__DEBUG__
	fprintf(stdout, "Form: `%s'\n", form);
#endif /* __DEBUG__ */

	if ((n1 != nDim) || (n2 != nDim)) {
	    GtoErrMsg("%s: array or matrix size incompatible with dimension\n",
			routine);
	    goto done;
	}


	/*
	 *
	 */
	if (nSample < 0) {
		mcType = 's';
		nSample = -nSample;
	} else {
		mcType = 'm';
	}

	if (DrlFormMonteCarlo(
		nDim,
		what,
		mcType,
		p,
		texp,
		vol,
		corrMat,
		form,
		nSample,
		doStatFlag,
		retVal) != SUCCESS) goto done;

	/*
	 *
	 */


	switch (ARGSIZE(retValL)) {
	case 1:
		retValL[1] = retVal[0];
		break;
	default:
		retValL[1] = retVal[0];
		retValL[2] = retVal[1];
		break;
	}


	/* made it through OK */
	status = SUCCESS;
done:
	/* free mem */
	DrlVTypeVectFree((void*) p, 0, DRL_DOUBLE_T);
	DrlVTypeVectFree((void*) vol, 0, DRL_DOUBLE_T);
	DrlVTypeMatrFree((void*) corrMat, nDim, nDim, DRL_DOUBLE_T);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed \n", routine);
	}
	return(status);
}



