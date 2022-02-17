/****************************************************************
 * Module:	VNFM
 * Submodule:	CALI
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include "cerror.h"
#include "macros.h"
#include "ldate.h"
#include "convert.h"
#include "date_sup.h"

#include "drlsort.h"
#include "drlio.h"			/* DrlFPrintf */

#define	_vnfm_SOURCE
#include "vnfmcali.h"

# undef __DEBUG__
#if defined(_WINDLL)
# undef __DEBUG__
#endif


/*f-------------------------------------------------------------
 * Bootstrap 2D volatility curve (TCurve input, arbitrary start/end).
 *                                                             
 * <br><br>
 *                                                             
 * <br><br>
 * Calibration of the two spot volatility curves.
 */

DLL_EXPORT(int)
VnfmCalib2VVolCurves(
	VnfmData *that,		/* (I/O) model parameters */

	int calType1,		/* (I) 0=sw cms, 1=sw final */
	double volMat1,		/* (I) swaption maturity */
	int volFreq1,		/* (I) 0=simple, 1, 2, 4, 12 */
	TCurve *volCurve1,	/* (I) input base volatility curve */

	int calType2,		/* (I) 0=sw cms, 1=sw final */
	double volMat2,		/* (I) swaption maturity */
	int volFreq2,		/* (I) 0=bvol, 1=swaption */
	TCurve *volCurve2,	/* (I) input base volatility curve */

	double tStart1,		/* (I) time to start 1-vol calib */
	double tStart2,		/* (I) time to start 2-vol calib */
	double tEnd2,		/* (I) time to   end 2-vol calib */
	double tEnd1,		/* (I) time to   end 1-vol calib */

	int errMsgFlag,		/* (I) reports bootstrap failure iff TRUE */
	double *resValue)	/* (O) <0 if calib OK, >0 otherwise */
{
static	char	routine[] = "VnfmCalib2VVolCurves";
	int	status = FAILURE, errCode;
	double	*tMat1 = NULL,
		*tMat2 = NULL,
		*vol1 = NULL,
		*vol2 = NULL,
		S, k0;
	int	*freq1 = NULL,
		*freq2 = NULL;
	int	i,
		idxStart1, idxStart2, idxEnd2, idxEnd1;


	/*
	 *
	 */
	/* check two-factor */
	if (that->fNf != 2) {
	    GtoErrMsg("%s: must be 2-factor (got %d factors)\n",
		routine, that->fNf);
	    goto done;
	}

	/* check that refdate of zc and bv curve agree */
	if (REFDATE != volCurve1->fBaseDate) {
	    GtoErrMsg("%s: zero ref date %10s != vol ref date %10s\n",
		routine, GtoFormatDate(REFDATE),
		GtoFormatDate(volCurve1->fBaseDate));
	    goto done;
	}
	if (VnfmComputeCoeff(that) != 0) goto done;


	/* get benchmarks vols */
	tMat1 = NEW_ARRAY(double, that->fNDates);
	tMat2 = NEW_ARRAY(double, that->fNDates);
	vol1  = NEW_ARRAY(double, that->fNDates);
	vol2  = NEW_ARRAY(double, that->fNDates);
	freq1 = NEW_ARRAY(int, that->fNDates);
	freq2 = NEW_ARRAY(int, that->fNDates);
	if ((tMat1 == NULL) || (vol1 == NULL) || (freq1 == NULL) ||
	    (tMat2 == NULL) || (vol2 == NULL) || (freq2 == NULL)) {
		GtoErrMsg("%s: malloc failed\n", routine);
		goto done;
	}




	/*
	 *
	 */
	for (i=0; i<=that->fNDates-1; i++) {
	    tMat1[i] = volMat1;
	    freq1[i] = volFreq1;
	    GtoInterpRate(that->fDate[i], volCurve1,
			GTO_LINEAR_INTERP, &vol1[i]);
	    GtoDayCountFraction(REFDATE, that->fDate[i],
		GTO_ACT_365F, &S);
	    tMat1[i] = volMat1 - (calType1 != 0 ? S : 0e0);
	}

	for (i=0; i<=that->fNDates-1; i++) {
	    tMat2[i] = volMat2;
	    freq2[i] = volFreq2;
	    GtoInterpRate(that->fDate[i], volCurve2,
			GTO_LINEAR_INTERP, &vol2[i]);
	    GtoDayCountFraction(REFDATE, that->fDate[i],
		GTO_ACT_365F, &S);
	    tMat2[i] = volMat2 - (calType2 != 0 ? S : 0e0);
	}


	/*
	 *
	 */
	VNFMIDX(tStart1, &idxStart1);
	VNFMIDX(tStart2, &idxStart2);
	VNFMIDX(tEnd2,   &idxEnd2);
	VNFMIDX(tEnd1,   &idxEnd1);

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s:\n", routine);\
	DrlFPrintf(vnfmFpLog, "\ttStart1=%lf idx=%d\n", tStart1, idxStart1);\
	DrlFPrintf(vnfmFpLog, "\ttStart2=%lf idx=%d\n", tStart2, idxStart2);\
	DrlFPrintf(vnfmFpLog, "\t  tEnd2=%lf idx=%d\n",   tEnd2,   idxEnd2);\
	DrlFPrintf(vnfmFpLog, "\t  tEnd1=%lf idx=%d\n",   tEnd1,   idxEnd1);)
	GTO_IF_LOGGING(\
	for (i=0; i<=that->fNDates-1; i++) {\
	    DrlFPrintf(vnfmFpLog, "\t  [%3d] tMat1=%8.4f freq1=%d vol1=%8.4f%%"\
		" tMat1=%8.4f freq1=%d vol1=%8.4f%%\n",\
		i, tMat1[i], freq1[i], vol1[i]*1e2,\
		tMat2[i], freq2[i], vol2[i]*1e2);\
	})
#endif

	/* (1) 1-vol calib */
	if (VnfmVolCalib1VArbitrary(
			that,	
			idxStart1, idxStart2-1,
			tMat1, freq1, vol1,
			LOGVOL,
			errMsgFlag,	
			resValue) != 0) goto done;



	/* (2) 2-vol calib */
	if (VnfmVolCalib2VArbitrary(
			that,	
			idxStart2,
			idxEnd2-1,
			that->fDate, tMat1, freq1, vol1,
			that->fDate, tMat2, freq2, vol2,
			errMsgFlag,	
			resValue) != 0) goto done;

	/* rescale spot vol2 to have equal last 
	 * calibrated spot vols */
	k0 = SIGMA[0][idxEnd2] / SIGMA[1][idxEnd2];
	for (i=0; i<=that->fNDates-1; i++) {
		SIGMA[1][i] *= k0;
	}
	ALPHA[1] /= k0;

	errCode = VnfmVolCalib1VArbitrary(
			that,	
			idxEnd2, idxEnd1,
			tMat1, freq1, vol1,
			LOGVOL,
			errMsgFlag,	
			resValue);

	/* rescale back */
	for (i=0; i<=that->fNDates-1; i++) {
	    	SIGMA[1][i] /= k0;
	}
	ALPHA[1] *= k0;

	if (errCode != 0) goto done;



	/* made it through OK */
	status = 0;
done:
	/* tesing of accuracy */
	VnfmComputeCoeff(that);

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: final parameters\n", routine);\
	VnfmFpWrite(that, vnfmFpLog);\
	VnfmPrintCoeff(that, vnfmFpLog);\
	DrlFPrintf(vnfmFpLog, "%s: done.", routine));
#endif

	/* free allocated memory */
	if (tMat1 != NULL) FREE((void*) tMat1);
	if (vol1  != NULL) FREE((void*) vol1);
	if (freq1 != NULL) FREE((void*) freq1);
	if (tMat2 != NULL) FREE((void*) tMat2);
	if (vol2  != NULL) FREE((void*) vol2);
	if (freq2 != NULL) FREE((void*) freq2);

	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}

