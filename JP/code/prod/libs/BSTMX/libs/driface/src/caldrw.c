/************************************************************************
 * Module:	DRIIFACE
 * Function:	DRWrapper Utilities
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <errno.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>

#include "date_sup.h"

#include "drlio.h"
#include "drlmem.h"
#include "drlsmat.h"
#include "drlts.h"
#include "drlstr.h"
#include "drlvtype.h"
#include "drltime.h"		/* DrlTDateScanYMD */

#define	 _vnfm_SOURCE		/* to use macros (RHOIDX) */
#include "vnfmanly.h"
#include "vnfmopca.h"

#include "dritkwrp.h"
#include "drical.h"	/* prototype consistency */

#define	__DEBUG__
#undef	__DEBUG__

static	int
drWrapperDataInterpVolD(
	TDrWrapperData *drwData,	/* (I) wrapper data */
	TDate expDate,			/* (I) rate reset date */
	TDateInterval matInt,		/* (I) rate maturity */
	int calibFinal,			/* (I) TRUE=final, FALSE=cms */
	double *cvolMat,		/* (O) calib rate maturity (fwd) */
        int *cvolFreq,			/* (O) calib rate frequency */
	double *cvolVol);		/* (O) calib rate volatility */

static	int
drWrapperDataInterpVol(
	TDrWrapperData *drWrapper,	/* (I) wrapper data */
        TDateInterval rateReset,	/* (I) rate reset */
        TDateInterval rateMat,		/* (I) rate maturity */
	int calibFinal,			/* (I) TRUE=final, FALSE=cms */
	double *cvolMat,		/* (O) calib rate maturity (fwd) */
        int *cvolFreq,			/* (O) calib rate frequency */
	double *cvolVol);		/* (O) calib rate volatility */


static	int
driCreateReport(
	VnfmData *vnfmData,		/* (I) model parameters */
	int nBenchMat,			/* (I) number of benchmarks */
	double *benchMat,		/* (I) array of benchmarks mat */
	int *benchFreq,			/* (I) array of benchmarks freq */
	FILE *fpOut);


#define	NBENCH		13
static		double	benchMat[] = {
	1e0/365e0,
	0.0833333333333e0, 0.25e0, 0.5e0, 1e0, 2e0,
	3e0, 5e0, 7e0, 10e0, 15e0,
	20e0, 30e0};
static		int	nBenchMat = NBENCH;
static		int	benchFreq[NBENCH];



/*---------------------------------------------------------------
 */

DLL_EXPORT(int)
DriCalibrateVnfmParams(
	const char *logFnam,		/* (I) log file (or NULL)*/
	const char *logLink,		/* (I) log file link (or NULL)*/
	const char *curCode,		/* (I) currency code (or NULL)*/
	const char *tagName,		/* (I) calibration name (or NULL) */
	FILE *fpIn,			/* (I) input */
	TDrWrapperData *drwData)	/* (I) mkt env */
{
static	char	routine[] = "DriCalibrateVnfmParams";
	int	status = FAILURE;

#define	NFMAX	8
#define	NBMAX	32

	int	nDim;			/* (I) number of factors */
	TDate		baseDate = drwData->fToday;

	char		calName[256];

	int		numVols;
	TDate		*volDates = NULL;
	double		*volMat = NULL;
	int		*volFreq = NULL;
	double		*volRates = NULL;

	int		numCvol;
	double		*cvolMat = NULL;
	int		*cvolFreq = NULL;
	double		*cvolRates = NULL;

	double		rateVol;
	char		buf[256];

	int		nPa;
	long		*paOptFlag = NULL;
	double		*paOpt = NULL;
	double		*paMin = NULL;
	double		*paMax = NULL;
	double		*paMid = NULL;

	int		idxB, idxF1, idxF2, idxR, idxE, idxM,
			nDim3;
	VnfmData	*vnfmData = NULL;

 	double	backboneq;		/* (I) backbone */

	long betaOpt[NFMAX];		/* (O) optimal */
	double betaMin[NFMAX];		/* (I) low  constr */
	double betaMax[NFMAX];		/* (I) high constr */
	double betaMid[NFMAX];		/* (I) guess (or value if no opt)*/

	long alphaOpt[NFMAX];		/* (O) optimal */
	double alphaMin[NFMAX];		/* (I) low  constr */
	double alphaMax[NFMAX];		/* (I) high constr */
	double alphaMid[NFMAX];		/* (I) guess (or value if no opt)*/

	long rhoOpt[NFMAX*NFMAX];	/* (O) optimal */
	double rhoMin[NFMAX*NFMAX];	/* (I) low  constr */
	double rhoMax[NFMAX*NFMAX];	/* (I) high constr */
	double rhoMid[NFMAX*NFMAX];	/* (I) guess (or value if no opt)*/

	int optimType;			/* (I) 0=const spot vol, 1=spsq */
	TDateInterval bcalMat;		/* (I) bootstrapping maturity */
	TDate	bcalDate;
	int bcalFinal;			/* (I) TRUE/FALSE */
	double sqspvMin;		/* (I) min spot volatility */

	int nOptBench;			/* (I) # of optimized benchmarks */
	TDateInterval optExpInt[NBMAX];	/* optimized benchmark exp (intvl) */
	TDateInterval optMatInt[NBMAX];	/* optimized benchmark mat (intvl) */
	double	optMatYrs[NBMAX];	/* optimized benchmark mat (yrs) */
	int optFreq[NBMAX];		/* optimized benchmark freq */
	TVolBenchmark optBench[NBMAX];
	double optWei[NBMAX];		/* (I) benchmarks weight */
	double optMod[NBMAX];		/* (O) model  value of benchmarks */
	double optMid[NBMAX];

	int nConstBench;		/* (I) # of constrained benchmarks */
	TDateInterval constExpInt[NBMAX];/* constimized benchmark exp (intvl) */
	TDateInterval constMatInt1[NBMAX];/* constimized benchmark mat  */
	TDateInterval constMatInt2[NBMAX];/* constimized benchmark mat  */
	double	constMatYrs1[NBMAX];	/* constimized benchmark mat (yrs) */
	double	constMatYrs2[NBMAX];	/* constimized benchmark mat (yrs) */
	int constFreq1[NBMAX];		/* constimized benchmark freq */
	int constFreq2[NBMAX];		/* constimized benchmark freq */
	TVolBenchmark constBench[NBMAX];
	double constMin[NBMAX];		/* (I) array of low  constraints */
	double constMax[NBMAX];		/* (I) array of high constraints */
	double constMod[NBMAX];		/* (O) model values of constraints */

	double	optL2 = 0e0;		/* Optimal quadratic distance */

	FILE	*fpLog = stdout;	/* (I) detailed report */
	FILE	*fpOut;			/* (I) summary 1-line output */

	/**
	 ** Read input data
	 **/

#define	READ_DATA(type, ptr, str) \
	{ if (DrlFScanVType(fpIn, type, (void*) ptr) != SUCCESS) \
	{ GtoErrMsg("%s: can't read %s.\n", routine, str); goto done;}}

	/* If tagName not provided, read from file */
	if (tagName == NULL) {
		READ_DATA(DRL_CHAR_ARRAY_T, calName, "boostrapping type");
	} else {
		strcpy(calName, tagName);
	}


	/* Open log file */
	if (logFnam && logFnam[0]) {
	    if ((fpLog = fopen(logFnam, "w")) == NULL) {
                GtoErrMsg("%s: can't open `%s' (%s).\n",
                        routine, logFnam, strerror(errno));
                goto done;
            }
	}

	/*
	 * Print 1-line summary 
	 */ 
	fpOut = fpLog;
	DrlFPrintf(fpOut, "\n");
	DrlFPrintf(fpOut, "SUMMARY:");
	DrlFPrintf(fpOut, " %-4s", (tagName ? curCode : "----"));
	DrlFPrintf(fpOut, "  %04d-%02d-%02d",
		DrlYEAR(baseDate), DrlMONTH(baseDate), DrlDAY(baseDate));
	DrlFPrintf(fpOut, "  %-30s", 
			(tagName ? tagName : "CALIBRATION"));

	if (logLink && logLink[0]) {
		DrlFPrintf(fpOut, " <a href=\"%s\">REPORT</a>", logLink);
	} else {
		DrlFPrintf(fpOut, "       ");
	}




	READ_DATA(DRL_INT_T, &nDim,"nDim");
	ASSERT_OR_DONE(nDim < NFMAX);

	READ_DATA(DRL_DOUBLE_T, &backboneq, "backboneq");


	/* read params spec */
	for (idxF1=0; idxF1<=nDim-1; idxF1++) {
	    READ_DATA(DRL_LONG_T,   &betaOpt[idxF1], "betaOpt");
	    READ_DATA(DRL_DOUBLE_T, &betaMin[idxF1], "betaMin");
	    READ_DATA(DRL_DOUBLE_T, &betaMax[idxF1], "betaMax");
	    READ_DATA(DRL_DOUBLE_T, &betaMid[idxF1], "betaMid");
	}
	for (idxF1=0; idxF1<=nDim-1; idxF1++) {
	    READ_DATA(DRL_LONG_T,   &alphaOpt[idxF1], "alphaOpt")
	    READ_DATA(DRL_DOUBLE_T, &alphaMin[idxF1], "alphaMin")
	    READ_DATA(DRL_DOUBLE_T, &alphaMax[idxF1], "alphaMax")
	    READ_DATA(DRL_DOUBLE_T, &alphaMid[idxF1], "alphaMid")
	}
 
	for (idxF1=0;       idxF1<=nDim-1; idxF1++)
	for (idxF2=idxF1+1; idxF2<=nDim-1; idxF2++) {
	    idxR = RHOIDX(idxF1, idxF2);
	    READ_DATA(DRL_LONG_T,   &rhoOpt[idxR], "rhoOpt");
	    READ_DATA(DRL_DOUBLE_T, &rhoMin[idxR], "rhoMin");
	    READ_DATA(DRL_DOUBLE_T, &rhoMax[idxR], "rhoMax");
	    READ_DATA(DRL_DOUBLE_T, &rhoMid[idxR], "rhoMid");
	}
 
	/* Read optim type */
	READ_DATA(DRL_INT_T, &optimType,  "optim type");

	/* Read bootstrapping params */
	READ_DATA(DRL_TDATEINTERVAL_T, &bcalMat,
		"bootstrapping maturiy");
	READ_DATA(DRL_CHAR_ARRAY_T, buf,
		"boostrapping type");
	switch (buf[0]) {
	case 'F':
	case 'f':
		bcalFinal = TRUE;
		break;
	case 'C':
	case 'c':
		bcalFinal = FALSE;
		break;
	default:
		GtoErrMsg("%s: bad calib type %s.\n", routine, buf);
		goto done;
	}




	/* Read optimize benchmarks */
	READ_DATA(DRL_INT_T, &nOptBench,"nOptBench");
	ASSERT_OR_DONE(nOptBench < NBMAX);

	for (idxB=0; idxB<nOptBench; idxB++) {
	    READ_DATA(DRL_TDATEINTERVAL_T, &optExpInt[idxB],  "optim exp");
	    READ_DATA(DRL_TDATEINTERVAL_T, &optMatInt[idxB], "optim mat");
	    READ_DATA(DRL_DOUBLE_T, &optWei[idxB], "optim weight");

	    IF_FAILED_DONE( drWrapperDataInterpVol(
		drwData,
		optExpInt[idxB],
		optMatInt[idxB],
		FALSE,
		&optMatYrs[idxB],
		&optFreq[idxB],
		&optMid[idxB]));


	    IF_FAILED_DONE( DrlTVolBenchmarkSetVol(
		&optBench[idxB],
		optExpInt[idxB],
		optMatInt[idxB],
        	optFreq[idxB],
        	-1L));	/* default */

	}

	/* Read constraint benchmarks */
	READ_DATA(DRL_INT_T, &nConstBench,"nConstBench");
	ASSERT_OR_DONE(nConstBench < NBMAX);

	for (idxB=0; idxB<nConstBench; idxB++) {
	    READ_DATA(DRL_TDATEINTERVAL_T, &constExpInt[idxB],  "const exp");
	    READ_DATA(DRL_TDATEINTERVAL_T, &constMatInt1[idxB], "const mat1");
	    READ_DATA(DRL_TDATEINTERVAL_T, &constMatInt2[idxB], "const mat2");
	    READ_DATA(DRL_DOUBLE_T, &constMin[idxB], "const min");
	    READ_DATA(DRL_DOUBLE_T, &constMax[idxB], "const max");


	    IF_FAILED_DONE( drWrapperDataInterpVol(
		drwData,
		constExpInt[idxB],
		constMatInt1[idxB],
		FALSE,
		&constMatYrs1[idxB],
		&constFreq1[idxB],
		&rateVol));

	    IF_FAILED_DONE( drWrapperDataInterpVol(
		drwData,
		constExpInt[idxB],
		constMatInt2[idxB],
		FALSE,
		&constMatYrs2[idxB],
		&constFreq2[idxB],
		&rateVol));


	    IF_FAILED_DONE( DrlTVolBenchmarkSetCorr(
		&constBench[idxB],
		constExpInt[idxB],
		constMatInt1[idxB],
        	constFreq1[idxB],
        	-1L,	/* default */
		constMatInt2[idxB],
		constFreq2[idxB],
		-1L));	/* default */
	}




	/*
	 * Interp volatility curve
	 */
	IF_FAILED_DONE( GtoDtFwdAny(
		baseDate,
		&bcalMat,
		&bcalDate));

	IF_FAILED_DONE( DriTDrWrapperDataGetInterpVol(
		drwData,
		bcalFinal,
		bcalDate,
		bcalMat,
		&numVols,
		&volDates,
		&volMat,
		&volFreq,
		&volRates));

	/*
	 * Set up vnfm
	 */
	vnfmData = VnfmNewTimeLine(
		nDim,
 		backboneq,
		baseDate,
		numVols,
		volDates,
		NULL,
		drwData->fDiscZcCurve);
	ASSERT_OR_DONE(vnfmData != NULL);
	nDim = nDim;
	
	IF_FAILED_DONE( VnfmComputeCoeff(vnfmData));



	/*
	 * Interp calibrated volatility
	 */
	numCvol = vnfmData->fNDates;
	cvolMat   = NEW_ARRAY(double, numCvol);
	cvolFreq  = NEW_ARRAY(int,    numCvol);
	cvolRates = NEW_ARRAY(double, numCvol);
	ASSERT_OR_DONE(cvolMat   != NULL);
	ASSERT_OR_DONE(cvolFreq  != NULL);
	ASSERT_OR_DONE(cvolRates != NULL);
	for (idxB=0; idxB<numCvol; idxB++) {

	    IF_FAILED_DONE( drWrapperDataInterpVolD(
		drwData,
		vnfmData->fDate[idxB],
		bcalMat,
		bcalFinal,
		&cvolMat[idxB],
		&cvolFreq[idxB],
		&cvolRates[idxB]));
	}





	/*
	 * Log
	 */
	GTO_IF_LOGGING( \
	DrlFPrintf(NULL, "Optim type: %d\n", optimType); \
	DrlFPrintf(NULL, "Bootstrap type: %s %s (%s->%s)\n", \
		(bcalFinal  ? "Fix" : "Cms"), \
		DrlTDateIntervalPrint(NULL, bcalMat), \
		DrlTDatePrint(NULL, baseDate), \
		DrlTDatePrint(NULL, bcalDate)); \
	for (idxB=0; idxB<numCvol; idxB++) { \
		DrlFPrintf(NULL, " %2d/%2d  %10s %8.4f %2d %12.8f\n", \
			idxB, numCvol, \
			DrlTDatePrint(NULL, vnfmData->fDate[idxB]), \
			cvolMat[idxB], \
			cvolFreq[idxB], \
			cvolRates[idxB]); \
	} \
	DrlFPrintf(NULL, "Optimize:\n"); \
	for (idxB=0; idxB<nOptBench; idxB++) { \
		DrlFPrintf(NULL, "%18s    %5.2f   %12.8f.\n", \
			DrlTVolBenchmarkPrint(NULL, &optBench[idxB]), \
			optWei[idxB], \
			optMid[idxB]); \
	} \
	for (idxB=0; idxB<nConstBench; idxB++) { \
		DrlFPrintf(NULL, "%18s    %12.8f  %12.8f.\n", \
			DrlTVolBenchmarkPrint(NULL, &constBench[idxB]), \
			constMin[idxB], \
			constMax[idxB]); \
	} \
	VnfmFpWrite(vnfmData, NULL); \
	);


	/*
	 * Calibrate
	 */
	if (optimType == 2) {

		/* Number of parameters to optimize:
		 * 1: mr and others
		 * mr:          (b1,...,bn)           n
		 * weights:     (a1,...,bn)           n
		 * corr:        (r11,r12,...,rn-1n)   n(n-1)/2
		 * Total:                             2n + n(n-1)/2
		 */
		nPa = 2*nDim + nDim*(nDim-1)/2;

		paOptFlag = NEW_ARRAY(long,   nPa);
		paOpt = NEW_ARRAY(double, nPa);
		paMin = NEW_ARRAY(double, nPa);
		paMax = NEW_ARRAY(double, nPa);
		paMid = NEW_ARRAY(double, nPa);


		/* fill with inputs */
		for (idxF1=0; idxF1<nDim; idxF1++) {
			paOptFlag[idxF1] = betaOpt[idxF1];
			paMin[idxF1] = MAX(betaMin[idxF1], 1e-8);
			paMax[idxF1] = betaMax[idxF1];
			paMid[idxF1] = MAX(betaMid[idxF1], 1e-8);
		}
		/* Rescale weights so that alpha1 = 1 */
		for (idxF1=0; idxF1<nDim; idxF1++) {
			paOptFlag[idxF1+nDim] = alphaOpt[idxF1];
			paMin[idxF1+nDim] = alphaMin[idxF1] ;
			paMax[idxF1+nDim] = alphaMax[idxF1] ;
			paMid[idxF1+nDim] = alphaMid[idxF1] ;
		}
 
		for (idxF1=0;       idxF1<nDim; idxF1++)
		for (idxF2=idxF1+1; idxF2<nDim; idxF2++) {
			idxR = RHOIDX(idxF1, idxF2);
			paOptFlag[idxR+nDim+nDim] = rhoOpt[idxR];
			paMin[idxR+nDim+nDim] = rhoMin[idxR];
			paMax[idxR+nDim+nDim] = rhoMax[idxR];
			paMid[idxR+nDim+nDim] = rhoMid[idxR];
		}
 
		GTO_IF_LOGGING( \
		DrlFPrintf(fpLog, "Pa:\n"); \
		for (idxB=0; idxB<nPa; idxB++) { \
			DrlFPrintf(fpLog, \
				"%2d/%2d  %ld  %12.8g  %12.8g %12.8g\n", \
				idxB, nPa, \
				paOptFlag[idxB], \
				paMin[idxB], \
				paMax[idxB], \
				paMid[idxB]); \
		});

		/* Min spot volatility */
		sqspvMin = 1e-8;


		IF_FAILED_DONE( VnfmCalibParamSquareVol(
			vnfmData,
			nOptBench,
			optBench,
			optWei,
			optMid,
			optMod,
			nConstBench,
			constBench,
			constMin,
			constMax,
			constMod,
			sqspvMin,
			cvolMat,
			cvolFreq,
			cvolRates,
			nPa,
			paOptFlag,
			paOpt,
			paMin,
			paMax,
			paMid));

		/* fill with ouputs */
		for (idxF1=0; idxF1<=nDim-1; idxF1++) {
			betaMid[idxF1] = paOpt[idxF1];
		}
		for (idxF1=0; idxF1<=nDim-1; idxF1++) {
			alphaMid[idxF1] = paOpt[idxF1+nDim];
		}
		for (idxF1=0;       idxF1<=nDim-1; idxF1++)
		for (idxF2=idxF1+1; idxF2<=nDim-1; idxF2++) {
			idxR = RHOIDX(idxF1, idxF2);
			rhoMid[idxR] = paOpt[idxR+nDim+nDim];
		}

		optL2 = 0e0;
		for (idxB=0; idxB<nOptBench; idxB++) {
			optL2 += optWei[idxB] * (optMod[idxB] - optMid[idxB])
				* optWei[idxB] * (optMod[idxB] - optMid[idxB]);
		}
		optL2 = sqrt(optL2/nOptBench);


 
 
	} else if ((optimType == 1)  || (optimType == 3)) {
		/* Number of parameters to optimize:
		 * 1: mr and others
		 * mr:          (b1,...,bn)           n
		 * weights:     (a1,...,bn)           n
		 * corr:        (r11,r12,...,rn-1n)   n(n-1)/2
		 * Total:                             2n + n(n-1)/2
		 */
		nPa = 2*nDim + nDim*(nDim-1)/2;

		paOptFlag = NEW_ARRAY(long,   nPa);
		paOpt = NEW_ARRAY(double, nPa);
		paMin = NEW_ARRAY(double, nPa);
		paMax = NEW_ARRAY(double, nPa);
		paMid = NEW_ARRAY(double, nPa);


		/* fill with inputs */
		for (idxF1=0; idxF1<nDim; idxF1++) {
			paOptFlag[idxF1] = betaOpt[idxF1];
			paMin[idxF1] = MAX(betaMin[idxF1], 1e-8);
			paMax[idxF1] = betaMax[idxF1];
			paMid[idxF1] = MAX(betaMid[idxF1], 1e-8);
		}
		/* Rescale weights so that alpha1 = 1 */
		for (idxF1=0; idxF1<nDim; idxF1++) {
			paOptFlag[idxF1+nDim] = alphaOpt[idxF1];
			paMin[idxF1+nDim] = alphaMin[idxF1];
			paMax[idxF1+nDim] = alphaMax[idxF1];
			paMid[idxF1+nDim] = alphaMid[idxF1];
		}
 
		for (idxF1=0;       idxF1<nDim; idxF1++)
		for (idxF2=idxF1+1; idxF2<nDim; idxF2++) {
			idxR = RHOIDX(idxF1, idxF2);
			paOptFlag[idxR+nDim+nDim] = rhoOpt[idxR];
			paMin[idxR+nDim+nDim] = rhoMin[idxR];
			paMax[idxR+nDim+nDim] = rhoMax[idxR];
			paMid[idxR+nDim+nDim] = rhoMid[idxR];
		}
 
		GTO_IF_LOGGING( \
		DrlFPrintf(fpLog, "Pa:\n"); \
		for (idxB=0; idxB<nPa; idxB++) { \
			DrlFPrintf(fpLog, \
				"%2d/%2d  %ld  %12.8g  %12.8g %12.8g\n", \
				idxB, nPa, \
				paOptFlag[idxB], \
				paMin[idxB], \
				paMax[idxB], \
				paMid[idxB]); \
		});





		IF_FAILED_DONE( VnfmCalibParamShortTerm(
			vnfmData,
			(optimType == 1 ? "F" : "R"),
			nOptBench,
			optBench,
			optWei,
			optMid,
			optMod,
			nConstBench,
			constBench,
			constMin,
			constMax,
			constMod,
			paOptFlag,
			paOpt,
			paMin,
			paMax,
			paMid));

		/* fill with ouputs */
		for (idxF1=0; idxF1<=nDim-1; idxF1++) {
			betaMid[idxF1] = paOpt[idxF1];
		}
		for (idxF1=0; idxF1<=nDim-1; idxF1++) {
			alphaMid[idxF1] = paOpt[idxF1+nDim];
		}
		for (idxF1=0;       idxF1<=nDim-1; idxF1++)
		for (idxF2=idxF1+1; idxF2<=nDim-1; idxF2++) {
			idxR = RHOIDX(idxF1, idxF2);
			rhoMid[idxR] = paOpt[idxR+nDim+nDim];
		}

		optL2 = 0e0;
		for (idxB=0; idxB<nOptBench; idxB++) {
			optL2 += optWei[idxB] * (optMod[idxB] - optMid[idxB])
				* optWei[idxB] * (optMod[idxB] - optMid[idxB]);
		}
		optL2 = sqrt(optL2/nOptBench);




		/* Perform bootstrapping for report */
			IF_FAILED_DONE( VnfmVolCalib1VArbitrary(
				vnfmData,
				0,
				vnfmData->fNDates-1,
				cvolMat,
				cvolFreq,
				cvolRates,
				LOGVOL,
				TRUE,
				NULL));
			IF_FAILED_DONE( VnfmComputeCoeff(vnfmData));

 
	} else if (optimType == 0) {
		/* No optimization */
		for (idxF1=0; idxF1<nDim; idxF1++) {
			vnfmData->fBeta[idxF1] = MAX(betaMid[idxF1], 1e-8);
			vnfmData->fAlpha[idxF1] = MAX(alphaMid[idxF1], 1e-8);
			for (idxB=0; idxB<numCvol; idxB++) {
				vnfmData->fSigma[idxF1][idxB] = 1e0;
			}
		}
		for (idxF1=0;       idxF1<nDim; idxF1++)
		for (idxF2=idxF1+1; idxF2<nDim; idxF2++) {
			idxR = RHOIDX(idxF1, idxF2);
			for (idxB=0; idxB<numCvol; idxB++) {
				vnfmData->fRho[idxR][idxB] = rhoMid[idxR];
			}
		}

		IF_FAILED_DONE( VnfmCheckValid(vnfmData));

		IF_FAILED_DONE( VnfmComputeCoeff(vnfmData));

		IF_FAILED_DONE( VnfmVolCalib1VArbitrary(
			vnfmData,
			0,
			INT_MAX,
			cvolMat,
			cvolFreq,
			cvolRates,
			LOGVOL,
			TRUE,
			NULL));

		GTO_IF_LOGGING(VnfmFpWrite(vnfmData, NULL));


		for (idxB=0; idxB<nOptBench; idxB++) {
			IF_FAILED_DONE( VnfmTVolBenchmarkValue(
				vnfmData,
				&optBench[idxB],
				&optMod[idxB]));
		}

		for (idxB=0; idxB<nConstBench; idxB++) {
			IF_FAILED_DONE( VnfmTVolBenchmarkValue(
				vnfmData,
				&constBench[idxB],
				&constMod[idxB]));
		}




	} else {
		GtoErrMsg("%s: bad optim type %d (0, 1, 2).\n",
			routine, optimType);
		goto done;
	}

	/*
	 * 1-line report
	 */
	nDim3 = MAX(nDim, 3);

	/*DrlFPrintf(fpOut, "%10s", DrlTDatePrint(NULL, baseDate));
	DrlFPrintf(fpOut, " %-30s", calName);
	DrlFPrintf(fpOut, " %-40s", calName); */

	DrlFPrintf(fpOut, " %d ", nDim);

	for (idxF1=0; idxF1<nDim3; idxF1++) {
	    DrlFPrintf(fpOut, "%c %8s ", (idxF1?' ':'|'),
		(idxF1 < nDim ? (sprintf(buf,"%8.4f", betaMid[idxF1]), buf) :
		"-"));
	}
	for (idxF1=0; idxF1<nDim3; idxF1++) {
	    DrlFPrintf(fpOut, "%c %8s ", (idxF1?' ':'|'),
		(idxF1 < nDim ? (sprintf(buf,"%8.4f", alphaMid[idxF1]), buf) :
		"-"));
	}
	for (idxF1=0;       idxF1<nDim3; idxF1++)
	for (idxF2=idxF1+1; idxF2<nDim3; idxF2++) {
		idxR = RHOIDX(idxF1, idxF2);
		DrlFPrintf(fpOut, "%c %8s ", (idxF1==0&&idxF2==1?'|':' '),
		    (idxR <= RHOIDX(nDim-2,nDim-1) ?
			(sprintf(buf,"%8.4f", rhoMid[idxR]), buf) : "-"));
	}
	/*DrlFPrintf(fpOut, " |  (%8.4f)", optL2*1e2);*/
	DrlFPrintf(fpOut, "\n");
	DrlFPrintf(fpOut, "\n");


	/*
	 * Create report
	 */
	if (TRUE) {
	    double		volMod, volMkt, expYrs,
				matYrs, matYrs2,
				volMat, volMat2;
	    int			volFreq, volFreq2;
	    TDate		expDate;
	    TDateInterval	expInt, matInt, matInt2;

	    /* Print calibration config and results
	     */
	    DrlFPrintf(fpLog, "DATE: %10s\n", DrlTDatePrint(NULL, baseDate));
	    DrlFPrintf(fpLog, "CALIBRATION NAME: %s\n", calName);

	    DrlFPrintf(fpLog, "NUM FACT: %d\n", nDim);
	    DrlFPrintf(fpLog, "BACKBONE: %s\n", 
		(IS_ALMOST_ZERO(vnfmData->fBackBoneQ-1e0) ? "LOGNORMAL" :
		(IS_ALMOST_ZERO(vnfmData->fBackBoneQ-0e0) ? "NORMAL" :
		(sprintf(buf, "%lf", vnfmData->fBackBoneQ), buf))));

	    DrlFPrintf(fpLog, "BETA:    ");
	    for (idxF1=0; idxF1<=nDim-1; idxF1++) {
		DrlFPrintf(fpLog, "  %8.4f", betaMid[idxF1]);
	    }
	    DrlFPrintf(fpLog, "\n");
	    DrlFPrintf(fpLog, "ALPHA:   ");
	    for (idxF1=0; idxF1<=nDim-1; idxF1++) {
		DrlFPrintf(fpLog, "  %8.4f", alphaMid[idxF1]);
	    }
	    DrlFPrintf(fpLog, "\n");
	    DrlFPrintf(fpLog, "RHO:     ");
	    for (idxF1=0;       idxF1<=nDim-1; idxF1++)
	    for (idxF2=idxF1+1; idxF2<=nDim-1; idxF2++) {
		idxR = RHOIDX(idxF1, idxF2);
		DrlFPrintf(fpLog, "  %8.4f", rhoMid[idxR]);
	    }
	    DrlFPrintf(fpLog, "\n");



	    switch (optimType) {
	    case 0:
	        DrlFPrintf(fpLog, "NO OPTIMIZATION:\n");
	    	DrlFPrintf(fpLog, "    BOOSTRAP-SPVOL %s %s\n",
			DrlTDateIntervalPrint(NULL, bcalMat),
			(bcalFinal ? "FIX" : "CMS"));
		break;
	    case 1:
	        DrlFPrintf(fpLog, "OPTIMIZATION TYPE:\n");
	    	DrlFPrintf(fpLog, "    CST-SPVOL (REPORT %s %s)\n",
			DrlTDateIntervalPrint(NULL, bcalMat),
			(bcalFinal ? "FIX" : "CMS"));
		break;
	    case 2:
	        DrlFPrintf(fpLog, "OPTIMIZATION TYPE:\n");
	    	DrlFPrintf(fpLog, "    BOOSTRAP-SPVOL %s %s\n",
			DrlTDateIntervalPrint(NULL, bcalMat),
			(bcalFinal ? "FIX" : "CMS"));
		break;
	    case 3:
	        DrlFPrintf(fpLog, "OPTIMIZATION TYPE:\n");
	    	DrlFPrintf(fpLog, "    CST-SPVOL VOL RATIO (REPORT %s %s)\n",
			DrlTDateIntervalPrint(NULL, bcalMat),
			(bcalFinal ? "FIX" : "CMS"));
		break;
	    default:
	        DrlFPrintf(fpLog, "OPTIMIZATION TYPE:\n");
	    	DrlFPrintf(fpLog, "    ??\n");
	    }
	    DrlFPrintf(fpLog, "\n");

	    DrlFPrintf(fpLog, "OPTIMIZED BENCHMARKS: %d\n", nOptBench); 
	    if (nOptBench > 0) DrlFPrintf(fpLog,
		"   EXP x  MAT    WEIGHT   MARKET     MODEL      DIFF    \n");
	    for (idxB=0; idxB<nOptBench; idxB++) {
		DrlFPrintf(fpLog,
			"  %4s x %4s    %5.2f   %8.4f%%  %8.4f%%  %8.4f%%\n",
			DrlTDateIntervalPrint(NULL, optExpInt[idxB]),
			DrlTDateIntervalPrint(NULL, optMatInt[idxB]),
			optWei[idxB],
			optMid[idxB]*1e2,
			optMod[idxB]*1e2,
			(optMod[idxB]-optMid[idxB])*1e2);

		/*DrlFPrintf(fpLog, "    %-75s   %8.4f%%    %5.2f\n",
			DrlTVolBenchmarkPrint(NULL, &optBench[idxB]),
			optMid[idxB]*1e2,
			optWei[idxB]);*/
	    }
	    DrlFPrintf(fpLog, "OPTIMAL QUADRATIC SPREAD: %8.4f%%\n", optL2*1e2);
	    DrlFPrintf(fpLog, "\n");


	    DrlFPrintf(fpLog, "CONSTRAINED BENCHMARKS: %d\n", nConstBench); 
	    if (nConstBench > 0) DrlFPrintf(fpLog,
		"   EXP CORR MAT1 x MAT2     MIN       MAX        MODEL  \n");
	    for (idxB=0; idxB<nConstBench; idxB++) {
		DrlFPrintf(fpLog,
			"  %4s      %4s x %4s   %8.4f  %8.4f   %8.4f.\n",
			DrlTDateIntervalPrint(NULL, constExpInt[idxB]),
			DrlTDateIntervalPrint(NULL, constMatInt1[idxB]),
			DrlTDateIntervalPrint(NULL, constMatInt2[idxB]),
			constMin[idxB],
			constMax[idxB],
			constMod[idxB]);

		/*DrlFPrintf(fpLog, "    %-75s   %12.8f  %12.8f.\n",
			DrlTVolBenchmarkPrint(NULL, &constBench[idxB]),
			constMin[idxB],
			constMax[idxB]);*/
	    }
	    DrlFPrintf(fpLog, "\n");



	    /* Base/swaption volatilities
	     */
	    DrlFPrintf(fpLog, "SWAPTION VOLATILITIES:\n");
	    DrlFPrintf(fpLog, "MOD-MKT");
	    for(idxM=0; idxM<nBenchMat; idxM++)
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxM]);
	    DrlFPrintf(fpLog, "\n");
	    for(idxE=0; idxE<nBenchMat; idxE++) {
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxE]);
	        for(idxM=0; idxM<nBenchMat; idxM++) {
			expYrs = benchMat[idxE];
			matYrs = benchMat[idxM];
			IF_FAILED_DONE( GtoYearsToDateInterval(
				expYrs, &expInt));
			IF_FAILED_DONE( GtoYearsToDateInterval(
				matYrs, &matInt));


			IF_FAILED_DONE( drWrapperDataInterpVol(
				drwData,
				expInt,
				matInt,
				FALSE,
				&volMat,
				&volFreq,
				&volMkt));
			benchFreq[idxM] = volFreq;
			IF_FAILED_DONE( VnfmAvgQBVol(
				vnfmData, 0e0, expYrs,
				expYrs, volMat, volFreq, LOGVOL, &volMod));


			DrlFPrintf(fpLog, " %7.2f", (volMod-volMkt)*1e2);
		}
		DrlFPrintf(fpLog, "\n");
	    }


	    /* Fwd volatilities
	     */
	    DrlFPrintf(fpLog, "FORWARD 1M-SPOT VOLATILITIES:\n");
	    DrlFPrintf(fpLog, "MOD    ");
	    for(idxM=0; idxM<nBenchMat; idxM++)
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxM]);
	    DrlFPrintf(fpLog, "\n");
	    for(idxE=0; idxE<nBenchMat; idxE++) {
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxE]);
	        for(idxM=0; idxM<nBenchMat; idxM++) {
			expYrs = benchMat[idxE];
			matYrs = benchMat[idxM];
			IF_FAILED_DONE( GtoYearsToDateInterval(
				expYrs, &expInt));
			IF_FAILED_DONE( GtoYearsToDateInterval(
				matYrs, &matInt));


			IF_FAILED_DONE( drWrapperDataInterpVol(
				drwData,
				expInt,
				matInt,
				FALSE,
				&volMat,
				&volFreq,
				&volMkt));

			IF_FAILED_DONE( VnfmAvgQBVol(
				vnfmData,
				expYrs, expYrs+1e0/12e0, expYrs+1e0/12e0,
				volMat, volFreq, LOGVOL, &volMod));


			DrlFPrintf(fpLog, " %7.2f", volMod*1e2);
		}
		DrlFPrintf(fpLog, "\n");
	    }


	    /* Correlations
	     */
	    DrlFPrintf(fpLog, "3M-SPOT CORRELATIONS:\n");
	    DrlFPrintf(fpLog, "MOD    ");
	    for(idxM=0; idxM<nBenchMat; idxM++)
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxM]);
	    DrlFPrintf(fpLog, "\n");
	    for(idxE=0; idxE<nBenchMat; idxE++) {
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxE]);
	        for(idxM=0; idxM<nBenchMat; idxM++) {
			expYrs = 0.25e0;
			matYrs  = benchMat[idxE];
			matYrs2 = benchMat[idxM];
			IF_FAILED_DONE( GtoYearsToDateInterval(
				expYrs, &expInt));
			IF_FAILED_DONE( GtoYearsToDateInterval(
				matYrs, &matInt));
			IF_FAILED_DONE( GtoYearsToDateInterval(
				matYrs2, &matInt2));

			IF_FAILED_DONE( drWrapperDataInterpVol(
				drwData,
				expInt,
				matInt,
				FALSE,
				&volMat,
				&volFreq,
				&volMkt));

			IF_FAILED_DONE( drWrapperDataInterpVol(
				drwData,
				expInt,
				matInt2,
				FALSE,
				&volMat2,
				&volFreq2,
				&volMkt));


			IF_FAILED_DONE( VnfmAvgQBCorr(
				vnfmData,
				expYrs,
				expYrs, matYrs,  volFreq,
				expYrs, matYrs2, volFreq2,
				&volMod));


			DrlFPrintf(fpLog, " %7.2f", volMod*1e2);
		}
		DrlFPrintf(fpLog, "\n");
	    }


	    /* Base/swaption volatilities
	     */
	    DrlFPrintf(fpLog, "SWAPTION VOLATILITIES:\n");
	    DrlFPrintf(fpLog, "MKT    ");
	    for(idxM=0; idxM<nBenchMat; idxM++)
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxM]);
	    DrlFPrintf(fpLog, "\n");
	    for(idxE=0; idxE<nBenchMat; idxE++) {
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxE]);
	        for(idxM=0; idxM<nBenchMat; idxM++) {
			expYrs = benchMat[idxE];
			matYrs = benchMat[idxM];
			IF_FAILED_DONE( GtoYearsToDateInterval(
				expYrs, &expInt));
			IF_FAILED_DONE( GtoYearsToDateInterval(
				matYrs, &matInt));


			IF_FAILED_DONE( drWrapperDataInterpVol(
				drwData,
				expInt,
				matInt,
				FALSE,
				&volMat,
				&volFreq,
				&volMkt));
			IF_FAILED_DONE( VnfmAvgQBVol(
				vnfmData, 0e0, expYrs,
				expYrs, volMat, volFreq, LOGVOL, &volMod));


			DrlFPrintf(fpLog, " %7.2f", volMkt*1e2);
		}
		DrlFPrintf(fpLog, "\n");
	    }

	    /* Base/swaption volatilities
	     */
	    DrlFPrintf(fpLog, "SWAPTION VOLATILITIES:\n");
	    DrlFPrintf(fpLog, "MOD    ");
	    for(idxM=0; idxM<nBenchMat; idxM++)
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxM]);
	    DrlFPrintf(fpLog, "\n");
	    for(idxE=0; idxE<nBenchMat; idxE++) {
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxE]);
	        for(idxM=0; idxM<nBenchMat; idxM++) {
			expYrs = benchMat[idxE];
			matYrs = benchMat[idxM];
			IF_FAILED_DONE( GtoYearsToDateInterval(
				expYrs, &expInt));
			IF_FAILED_DONE( GtoYearsToDateInterval(
				matYrs, &matInt));


			IF_FAILED_DONE( drWrapperDataInterpVol(
				drwData,
				expInt,
				matInt,
				FALSE,
				&volMat,
				&volFreq,
				&volMkt));
			IF_FAILED_DONE( VnfmAvgQBVol(
				vnfmData, 0e0, expYrs,
				expYrs, volMat, volFreq, LOGVOL, &volMod));


			DrlFPrintf(fpLog, " %7.2f", volMod*1e2);
		}
		DrlFPrintf(fpLog, "\n");
	    }


	    IF_FAILED_DONE( driCreateReport(
			vnfmData,
			nBenchMat,
			benchMat,
			benchFreq,
			fpLog));


	}



	/* OK */
	status = SUCCESS;
done:
	FREE(paOptFlag);
	FREE(paOpt);
	FREE(paMin);
	FREE(paMax);
	FREE(paMid);

	FREE(volDates);
	FREE(volMat);
	FREE(volFreq);
	FREE(volRates);

	FREE(cvolMat);
	FREE(cvolFreq);
	FREE(cvolRates);

	VnfmFree(vnfmData);

	if (status != SUCCESS) {
		/* Print failure report */
		DrlFPrintf(fpOut, " %d ", nDim);
		DrlFPrintf(fpOut, "| FAILED\n");

		GtoErrMsg("%s: failed\n", routine);
	}

	if (logFnam && logFnam[0]) {
		if (fpLog) fclose(fpLog);
	}

	return(status);
}




/*---------------------------------------------------------------
 * Convenience routine to return the interpolated volatility
 * and rate frequency (i.e. base or swaption) from
 * a market environment.
 */

static	int
drWrapperDataInterpVolD(
	TDrWrapperData *drwData,	/* (I) wrapper data */
	TDate expDate,			/* (I) rate reset date */
	TDateInterval matInt,		/* (I) rate maturity */
	int calibFinal,			/* (I) TRUE=final, FALSE=cms */
	double *cvolMat,		/* (O) calib rate maturity (fwd) */
        int *cvolFreq,			/* (O) calib rate frequency */
	double *cvolVol)		/* (O) calib rate volatility */
{
static	char	routine[] = "DriTDrWrapperDataGetInterpVol";
	int	status = FAILURE;

	TDate		baseDate = drwData->fToday;
	TDate		matDate;
	double		bvMat,
			matYrs,
			matYrsMin,
			matYrsMax,
			expYrs,
			minMat = 0.25;


	/* base vol maturity */
	bvMat = 1e0 / (double) drwData->fBvCurve->fBasis;
	minMat = bvMat;

	if (calibFinal) {
		IF_FAILED_DONE( GtoDtFwdAny(
			baseDate,
			&matInt,
			&matDate));
		IF_FAILED_DONE( GtoDayCountFraction(
			expDate,
			matDate,
			GTO_B30_360,
			&matYrs));

	} else {
		IF_FAILED_DONE( GtoDateIntervalToYears(
			&matInt,
			&matYrs));
	}


	/*
	 * If final calibration or if CMS calib with
	 * mat larger than base vol mat, use swaptions
	 * Otherwise use base vol.
	 */
	if ((calibFinal) ||
	    (matYrs > bvMat + 1e-4)) {
		/*
		 * Interp from swaption matrix
		 */
		matYrsMin = 1e0/drwData->fCmsSwMat->swapPayFreq;
		matYrsMax = TSWAPTION_MATRIX2D_MAT(drwData->fCmsSwMat,
				TSWAPTION_MATRIX2D_NMAT(drwData->fCmsSwMat)-1);
	
		matYrs = MAX(matYrs, matYrsMin);
		matYrs = MIN(matYrs, matYrsMax);

		IF_FAILED_DONE( GtoDayCountFraction(
			baseDate,
			expDate,
			GTO_B30_360,
			&expYrs));

		/* Perform interp
		 */
		IF_FAILED_DONE( DrlTSwaptionMatrix2DInterpExpMat(
			drwData->fCmsSwMat,
			cvolVol,
			expYrs,
			matYrs,
			FALSE));	/* (I) TRUE=adjoint, FALSE=direct*/

		*cvolFreq = drwData->fCmsSwMat->swapPayFreq;
		*cvolMat = matYrs;

	} else {
		/*
		 * Interp from base volatility
		 */

		IF_FAILED_DONE( GtoInterpRate(
			expDate,
			drwData->fBvCurve,
                        GTO_LINEAR_INTERP,
			cvolVol));

		*cvolFreq = (int) (1e0 / bvMat);
		*cvolMat = bvMat;
	}
	
	
	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}






static	int
drWrapperDataInterpVol(
	TDrWrapperData *drwData,	/* (I) wrapper data */
	TDateInterval expInt,		/* (I) rate reset */
	TDateInterval matInt,		/* (I) rate maturity */
	int calibFinal,			/* (I) TRUE=final, FALSE=cms */
	double *cvolMat,		/* (O) calib rate maturity (fwd) */
        int *cvolFreq,			/* (O) calib rate frequency */
	double *cvolVol)		/* (O) calib rate volatility */
{
	int	status = FAILURE;
	TDate	expDate;
	TDate	baseDate = drwData->fToday;

	IF_FAILED_DONE( GtoDtFwdAny(
		baseDate,
		&expInt,
		&expDate));

	IF_FAILED_DONE( drWrapperDataInterpVolD(
		drwData,
		expDate,
		matInt,
		calibFinal,
		cvolMat,
        	cvolFreq,
		cvolVol));

	/* OK */
	status = SUCCESS;
done:
	return(status);
}



/*---------------------------------------------------------------
 *
 */
static	int
driCreateReport(
	VnfmData *vnfmData,		/* (I) model parameters */
	int nBenchMat,			/* (I) number of benchmarks */
	double *benchMat,		/* (I) array of benchmarks mat */
	int *benchFreq,			/* (I) array of benchmarks freq */
	FILE *fpLog)
{
static	char	routine[] = "DriGetFromLiverate";
	int	status = FAILURE;

	double	re[NBENCH],	/* auxiliary */
		rm[NBENCH],
		t0, t1;
	int	rf[NBENCH];
	double  *eigVal = NULL, eigValNorm;
	double  **eigVect = NULL;
	int	nDim = vnfmData->fNf,
		idxM, idxE;


	eigVal = DrlDoubleVectAlloc(0, nBenchMat-1);
	ASSERT_OR_DONE(eigVal != NULL);
	eigVect = DrlDoubleMatrAlloc(0, nBenchMat-1, 0, nBenchMat-1);
	ASSERT_OR_DONE(eigVect != NULL);
 


	IF_FAILED_DONE( VnfmComputeCoeff(vnfmData));


	/* Spot Yields Factors
	 */
	t0 = 0e0;
	t1 = 0e0;
	for(idxM=0; idxM<nBenchMat; idxM++) {
		re[idxM] = t1;
		rf[idxM] = benchFreq[idxM];
		rm[idxM] = benchMat[idxM];
	}

	IF_FAILED_DONE( VnfmAvgQBOrthFactors(
		vnfmData,
		t0, t1,
		nBenchMat,
		re,
		rf,
		rm,
		eigVect,
		eigVal));

	for(idxM=0, eigValNorm=0e0 ; idxM<nBenchMat; idxM++)
		eigValNorm += eigVal[idxM];

	DrlFPrintf(fpLog, "YIELD FACTORS:\n");
	DrlFPrintf(fpLog, "NO EXPL");
	for(idxM=0; idxM<nBenchMat; idxM++)
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxM]);
	DrlFPrintf(fpLog, "\n");
	for(idxE=0; idxE<MIN(nBenchMat,nDim); idxE++) {
		/*DrlFPrintf(fpLog, " %7.2f", eigVal[idxE]*1e2);*/
		DrlFPrintf(fpLog, " %d %4.0f%%", idxE+1,
			eigVal[idxE]/eigValNorm*1e2);
	        for(idxM=0; idxM<nBenchMat; idxM++) {
			DrlFPrintf(fpLog, " %7.2f", 
				eigVect[idxM][idxE]
				*sqrt(fabs(eigVal[idxE])) *1e2);
		}
		DrlFPrintf(fpLog, "\n");
	}



	/* Forward 3M Factors
	 */
	t0 = 0e0;
	t1 = 0e0;
	for(idxM=0; idxM<nBenchMat; idxM++) {
		re[idxM] = benchMat[idxM];
		rf[idxM] = 0;
		rm[idxM] = 0.25e0;
	}

	IF_FAILED_DONE( VnfmAvgQBOrthFactors(
		vnfmData,
		t0, t1,
		nBenchMat,
		re,
		rf,
		rm,
		eigVect,
		eigVal));

	for(idxM=0, eigValNorm=0e0 ; idxM<nBenchMat; idxM++)
		eigValNorm += eigVal[idxM];

	DrlFPrintf(fpLog, "3M-FRA FACTORS:\n");
	DrlFPrintf(fpLog, "NO EXPL");
	for(idxM=0; idxM<nBenchMat; idxM++)
		DrlFPrintf(fpLog, " %7.2f", benchMat[idxM]);
	DrlFPrintf(fpLog, "\n");
	for(idxE=0; idxE<MIN(nBenchMat,nDim); idxE++) {
		DrlFPrintf(fpLog, " %d %4.0f%%", idxE+1,
			eigVal[idxE]/eigValNorm*1e2);
	        for(idxM=0; idxM<nBenchMat; idxM++) {
			DrlFPrintf(fpLog, " %7.2f", 
				eigVect[idxM][idxE]
				*sqrt(fabs(eigVal[idxE])) *1e2);
		}
		DrlFPrintf(fpLog, "\n");
	}




	/* OK */
	status = SUCCESS;
done:
	DrlDoubleVectFree(eigVal, 0, nBenchMat-1);
	DrlDoubleMatrFree(eigVect, 0, nBenchMat-1, 0, nBenchMat-1);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}


