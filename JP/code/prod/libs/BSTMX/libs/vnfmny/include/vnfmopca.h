/****************************************************************
 * Module:	VNFM
 * Submodule:	CALI
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vnfmopca_H
#define	_vnfmopca_H
#include "drlstd.h"		/* DLL_EXPORT() */

#include <stdio.h>
#include <math.h>

#include "bastypes.h"		/* TFloatRate */
#include "drlsmat.h"
#include "drlfrate.h"		/* TVolBenchmark, etc.. */

#include "vnfmanly.h"



extern	int	VnfmTVolBenchmarkValue(
	VnfmData *that,		/* (I) model parameters */
	TVolBenchmark *pt,	/* (I) benchmark */
	double *retVal);	/* (O) value */

extern	int	VnfmTVolBenchmarkRescaleAlpha(
	VnfmData *that,		/* (B) model parameters */
	int nOptBench,		/* (I) # of benchmarks */
	double *optWei,		/* (I) weight of benchmarks in average */
	double *optMid,		/* (I) market value of benchmarks */
	TVolBenchmark *optBench);/* (I) array of benchmarks */





/* Parameter calibration (optimization) routine */
extern	int	VnfmCalibParamShortTerm(
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
				   0=beta1, 1=beta2, 3=sigma1, 4=sigma2,
				   5=rho */
	long *paOptFlag,	/* (I) array of flags (TRUE=optimize) */
	double *paOpt,		/* (O) optimal paramters */
	double *paMin,		/* (I) low  constraints on parameters */
	double *paMax,		/* (I) high constraints on parameters */
	double *paMid);		/* (I) initial guess (or value if no opt)*/


#define	VNFM_NPARAM	16	/* number of possibly optimized parameters */
#define	VNFM_NMAX_PARAM	32
#define	VNFM_NMAX_BENCH	32
#define	VNFM_NMAX_CORR	32


extern int
VnfmComputeJSquare(
        VnfmData *that,         /* (I) model parameters */
        double *sqspv,          /* (I) array of square spot vols */
        double **jt);           /* (O) pseudo-J coeffs */

extern int
VnfmAvgQBVolSquare(
        VnfmData *that,
        double S1,              /* (I) start */
        double S2,              /* (I) expiration */
        double tReset,          /* (I) rate reset */
        double tMat,            /* (I) rate maturity */
        int freq,               /* (I) rate frequency */
        double *sqspv,          /* (I) array of square spot vols */
        double **jt,            /* (I) pseudo-J coeffs */
        double *retVal);        /* (O) volatility */


extern	int	VnfmCalibParamSquareVol(
	VnfmData *that,		/* (I/O) model parameters */

	int nOptBench,		/* (I) # of optimized benchmarks */
	TVolBenchmark *optBench,/* (I) array of optimized benchmarks */
	double *optWeight,	/* (I) benchmarks weight */
	double *optMid,		/* (I) market value of benchmarks */
	double *optMod,		/* (O) model  value of benchmarks */

	int nConstBench,	/* (I) # of constrained benchmarks */
	TVolBenchmark *constBench,/* (I) array of constrained benchmarks */
	double *constMin,	/* (I) array of low  constraints */
	double *constMax,	/* (I) array of high constraints */
	double *constMod,	/* (O) model values of constraints */

	double sqspvMin,	/* (I) min spot volatility */
	double *cRateMat,	/* (I) cal swaption mat [0..idxEnd] */
	int *cRateFreq,		/* (I) cal swaption freq [0..idxEnd] */
	double *cRateVol,	/* (I) cal swaption vol [0..idxEnd] */

				/* Optim Params array: see below  */
	int nPa,		/* (I) len of arrays */
	long *paOptFlag,	/* (I) array of flags (TRUE=optimize) */
	double *paOpt,		/* (O) optimal paramters */
	double *paMin,		/* (I) low  constraints on parameters */
	double *paMax,		/* (I) high constraints on parameters */
	double *paMid);		/* (I) initial guess (or value if no opt)*/


/* */
extern	DLL_EXPORT(int)	VnfmSmoothSwaptionMatrix(
	TCurve *zcCurve,		/* (I) zero curve */
	TSwaptionMatrix2D *midMkt,	/* (I) */
	TSwaptionMatrix2D *bidToMid,	/* (I) */
	TSwaptionMatrix2D *matWeigth,	/* (I) */
	TSwaptionMatrix2D *expWeigth,	/* (I) */
	int optimType,			/* (I) */
	int normType,			/* (I) */
	double smoothParam,		/* (I) */
	double *nfParamsIn,		/* (I) N-fact params in wrapper form */
	TSwaptionMatrix2D *origMkt,	/* (I) used for timeline only */
	TSwaptionMatrix2D **outputMat,	/* (O) output matrix */
	TSwaptionMatrix2D **spvolMat,	/* (O) output spot vol */
	double *nfParamsOut);		/* (O) */



#endif	/* _vnfmopca_H */

