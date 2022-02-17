/****************************************************************
 * Module:	VNFM
 * Submodule:	CALI
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include "macros.h"
#include "convert.h"
#include "date_sup.h"

#include "drlfrate.h"			/* TVolBenchmark */
#include "drlio.h"			/* DrlFPrintf */

#define	_vnfm_SOURCE
#include "vnfmopca.h"

#undef __DEBUG__
#if defined(_WINDLL) || !defined(TESTLIB)
# undef __DEBUG__
#endif


/*f-------------------------------------------------------------
 * Compute volatility of an arbitrary benchmark.
 *                                                             
 * <br><br>
 * Computes the value of a volatility benchmark "pt".
 */

int
VnfmTVolBenchmarkValue(
	VnfmData *that,		/* (I) model parameters */
	TVolBenchmark *pt,	/* (I) benchmark */
	double *retVal)		/* (O) value */
{
static  char    routine[] = "VnfmTVolBenchmarkValue";
	int     status = FAILURE;

	TDate		expDate;
        double		freq1, freq2;
	double		U1=0e0, U2, tReset1, tMat1, tReset2, tMat2;

	TFloatRateVolPoint      *vpt;
	TFloatRateCorrPoint     *cpt;

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: %s\n",\
		routine, DrlTVolBenchmarkPrint(NULL, pt)))
#endif

	switch(pt->fType) {
	case DRL_TVOLBENCHMARK_VOL:
	    /* volatility point */
	    vpt = &pt->fU.fVol;

	    if (GtoDateIntervalToFreq(&vpt->fRate.payInterval, &freq1)
		!= SUCCESS) return(FAILURE);

	    if (GtoDtFwdAny(REFDATE, &vpt->fStart, &expDate)
		!= SUCCESS) return(FAILURE);
	    if (GtoDayCountFraction(REFDATE, expDate, GTO_ACT_365F, &U1)
		!= SUCCESS) return(FAILURE);

	    if (GtoDtFwdAny(REFDATE, &vpt->fExp, &expDate)
		!= SUCCESS) return(FAILURE);
	    if (GtoDayCountFraction(REFDATE, expDate, GTO_ACT_365F, &U2)
		!= SUCCESS) return(FAILURE);

	    if (GtoDtFwdAny(REFDATE, &vpt->fReset, &expDate)
		!= SUCCESS) return(FAILURE);
	    if (GtoDayCountFraction(REFDATE, expDate, GTO_ACT_365F, &tReset1)
		!= SUCCESS) return(FAILURE);

	    if (GtoDateIntervalToYears(&vpt->fRate.matInterval, &tMat1)
		!= SUCCESS) return(FAILURE);

	    if (VnfmAvgQBVol(that, U1, U2, tReset1,
			tMat1, (int) freq1, LOGVOL, retVal)
		!= SUCCESS) return(FAILURE);



#if !defined(NO_LOGGING) && defined(__DEBUG__)
	   GTO_IF_LOGGING(\
	   DrlFPrintf(vnfmFpLog, "%s: U1=%8.4f U2=%8.4f tReset=%8.4f "\
		"tMat=%8.4f freq=%d vol=%12.8f\n",\
		routine,\
		U1, U2, tReset1, tMat1, (int) freq1, *retVal))
#endif

		break;

	case DRL_TVOLBENCHMARK_CORR:
	    /* average correlation point */
	    cpt = &pt->fU.fCorr;


	    if (GtoDateIntervalToFreq(&cpt->fRate1.payInterval, &freq1)
		!= SUCCESS) goto done;
	    if (GtoDateIntervalToFreq(&cpt->fRate2.payInterval, &freq2)
		!= SUCCESS) goto done;

	    if (GtoDtFwdAny(REFDATE, &cpt->fExp, &expDate)
		!= SUCCESS) goto done;
	    if (GtoDayCountFraction(REFDATE, expDate, GTO_ACT_365F, &U2)
		!= SUCCESS) goto done;

	    if (GtoDtFwdAny(REFDATE, &cpt->fReset1, &expDate)
		!= SUCCESS) goto done;
	    if (GtoDayCountFraction(REFDATE, expDate, GTO_ACT_365F, &tReset1)
		!= SUCCESS) goto done;
	    if (GtoDtFwdAny(REFDATE, &cpt->fReset2, &expDate)
		!= SUCCESS) goto done;
	    if (GtoDayCountFraction(REFDATE, expDate, GTO_ACT_365F, &tReset2)
		!= SUCCESS) goto done;

	    if (GtoDateIntervalToYears(&cpt->fRate1.matInterval, &tMat1)
		!= SUCCESS) goto done;
	    if (GtoDateIntervalToYears(&cpt->fRate2.matInterval, &tMat2)
		!= SUCCESS) goto done;

	    if (VnfmAvgQBCorr(that,
		U2,
		tReset1, tMat1, (int) freq1,
		tReset2, tMat2, (int) freq2,
		retVal)
	    != SUCCESS) goto done;



#if !defined(NO_LOGGING) && defined(__DEBUG__)
	    GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "%s: U1=%8.4f U2=%8.4f "\
		"tReset1=%lf, tMat1=%lf, freq1=%d, "\
		"tReset2=%lf, tMat2=%lf, freq2=%d, "\
		"retVal=%lf\n",\
		routine, U1, U2,\
		tReset1, tMat1, (int) freq1,\
		tReset2, tMat2, (int) freq2,\
		*retVal))
#endif
 	    break;
	default:
	    GtoErrMsg("%s: bad type %d.\n", routine, pt->fType);
	    goto done;
	}


	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}

/*f-------------------------------------------------------------
 * Rescale factor weights to best fit a set of benchmarks.
 *                                                             
 * <br><br>
 * Rescales the alpha parameters to fit the average volatility
 * of a set of benchmarks.
 */

int
VnfmTVolBenchmarkRescaleAlpha(
	VnfmData *that,		/* (B) model parameters */
	int nOptBench,		/* (I) # of benchmarks */
	double *optWei,		/* (I) weight of benchmarks in average */
	double *optMid,		/* (I) market value of benchmarks */
	TVolBenchmark *optBench)/* (I) array of benchmarks */
{
static	char	routine[] = "VnfmRescaleParams";
	int	status = FAILURE;
	int	idxB, idxF;
	double	vMod,
		vMkt,
		vModAvg = 0e0,
		vMktAvg = 0e0;

	/* Compute avg model and market values
	 */
	for (idxB=0; idxB<nOptBench; idxB++) {
	    if (!IS_ALMOST_ZERO(optWei[idxB])) {
		IF_FAILED_DONE( VnfmTVolBenchmarkValue(
			    that,
			    &optBench[idxB],
			    &vMod));
		vMkt = optMid[idxB];

		vModAvg += vMod * optWei[idxB];
		vMktAvg += vMkt * optWei[idxB];
	    }
	}
	if (IS_ALMOST_ZERO(vModAvg)) {
		GtoErrMsg("%s: model volatility is zero.\n", routine);
		goto done;
	}

	GTO_IF_LOGGING(\
	    VnfmFpWrite(that, vnfmFpLog);\
	    DrlFPrintf(vnfmFpLog, "%s: vModAvg=%8.4f vMktAvg=%8.4f\n",\
		routine, vModAvg, vMktAvg));

	/* Rescale alphas
	 */
	for (idxF=0; idxF<that->fNf; idxF++) {
		that->fAlpha[idxF] *= (vMktAvg / vModAvg);
	}

	GTO_IF_LOGGING(\
	    VnfmFpWrite(that, vnfmFpLog));


	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}








