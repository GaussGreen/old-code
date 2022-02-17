/************************************************************************
 * Module:	DRL
 * Submodule:	RAND
 * Function:	Random Number Generation and Simulation
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

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

#undef	SQR
#define	SQR(x)	((x)*(x))

/*----------------------------------------------------------------------
 * Simulation of multivariate normal distribution
 * Parsing formula
 */

int
DrlFormMonteCarlo(
	int	nDim,	 	/* (I) number of dimensions */
	char	what,		/* (I) 'n' normal, 'l' log-normal */
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
		DrlErrMsg("%s: dim too large %d\n", routine, nDim);
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
		DrlErrMsg("%s: Init multi-normal simulation failed.\n",
			routine);
		DrlErrMsg("%s: Input correlation matrix:\n", routine);
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
			DrlErrMsg("%s: bad type `%c'\n", routine, what);
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
		DrlErrMsg("%s: cntVal (%ld) != nSample (%ld).\n", 
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
		DrlErrMsg("%s: failed.\n", routine);
	}

	return(status);
}


