/****************************************************************
 * Module:	VNFM
 * Submodule:	WRAP
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <math.h>

#include "convert.h"			/* GtoFormatDate() */
#include "drlvtype.h"
#include "drlts.h"			/* DrlTCurveWrap() */
#include "drlsort.h"			/* DrlTDateArrayFloorIdx() */

#define	_vnfm_SOURCE
#include "vnfmcali.h"

#include "vnfmwrap.h"	/* Prototype Consistency */

#if defined(_WINDLL)
# undef __DEBUG__
#endif


/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CAL2V_GEN.
 *                                                             
 * <br><br>
 * Wrapper for double spot volatility calibration (2 factor only).
 * Calls <i> VnfmVolCalib1VArbitrary</i> and <i> VnfmVolCalib2VArbitrary</i>.
 */

DLL_EXPORT(int)
VnfmCalib2VVolCurvesL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDatesL,	/*  2 'D' (I) zero coupon dates */
	FloatL *zcRatesL,	/*  3 'F' (I) zero coupon rates */
				/*        Input Model: */
 	FloatL *backboneqL,	/*  4 'F' (I) 0=lognormal, 0.5=normal */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array correlations */
				/*        Input base vol curve: */
	TDate *rateReset1L,	/* 10 'D' (I) rate 1 reset [0..nDates-1] */
	FloatL *rateMat1L,	/* 11 'F' (I) rate 1 mat [0..nDates-1] */
	IntL *rateFreq1L,	/* 12 'L' (I) rate 1 freq [0..nDates-1] */
	FloatL *rateVol1L,	/* 13 'F' (I) rate 1 vol [0..nDates-1] */
	TDate *rateReset2L,	/* 14 'D' (I) rate 2 reset [0..nDates-1] */
	FloatL *rateMat2L,	/* 15 'F' (I) rate 2 mat [0..nDates-1] */
	IntL *rateFreq2L,	/* 16 'L' (I) rate 2 freq [0..nDates-1] */
	FloatL *rateVol2L,	/* 17 'F' (I) rate 2 vol [0..nDates-1] */

	TDateL *calDatesL,	/* 18 'D' (I) array of calib start/end dates*/
				/*        [0] = time start 1-vol calib */
				/*        [1] = time start 2-vol calib */
				/*        [2] = time   end 2-vol calib */
				/*        [3] = time   end 1-vol calib */
				/*        Output Calibrated Model: */

	FloatL *outSigma1L,	/* 19 'F' (O) spot vol 1 array */
	FloatL *outSigma2L)	/* 20 'F' (O) spot vol 2 array */
{
static	char		routine[] = "VnfmCalib2VVolCurvesL";
	int		status = FAILURE;
	VnfmData	*tfcalData = NULL;
	TCurve		*zcCurve = NULL;

	int		*freq1 = NULL,
			*freq2 = NULL;
	int		i, errCode,
			idxStart1, idxStart2, idxEnd2, idxEnd1;
	double		ratioVol;



	/* log inputs */
	if (GtoLoggingGet() > 0) {
	    DrlLilVectLoggingFile("wrapper.log", "w", "TF_SMOOTH_SWMAT",
		DRL_TDATE_L,  refDateL, "ZC_REFDATE",
		DRL_TDATE_L,  zcDatesL, "ZC_DATES",
		DRL_FLOAT_L,  zcRatesL, "ZC_RATES",
 		DRL_FLOAT_L, backboneqL,"BACKBONE_q",
		DRL_FLOAT_L, betaL,	"IN_BETA",
		DRL_FLOAT_L, alphaL, "IN_ALPHA",
		DRL_TDATE_L, dateL, "IN_DATES",
		DRL_FLOAT_L, sigmaL, "IN_SIGMA",
		DRL_FLOAT_L, rhoL, "IN_RHO",
		DRL_TDATE_L, rateReset1L, "RATE1_RESET",
		DRL_FLOAT_L, rateMat1L, "RATE1_MAT",
		DRL_LONG_L,  rateFreq1L, "RATE1_FREQ",
		DRL_FLOAT_L, rateVol1L, "RATE1_VOL",
		DRL_TDATE_L, rateReset2L, "RATE2_RESET",
		DRL_FLOAT_L, rateMat2L, "RATE2_MAT",
		DRL_LONG_L,  rateFreq2L, "RATE2_FREQ",
		DRL_FLOAT_L, rateVol2L, "RATE2_VOL",
		DRL_TDATE_L, calDatesL, "CAL_DATES",
		DRL_FLOAT_L, outSigma1L, "OUT_SIGMA1",
		DRL_FLOAT_L, outSigma2L, "OUT_SIGMA2",
		0);
	}


	VNFM_LOGOPEN

	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDatesL, zcRatesL)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlTCurveFpWrite(zcCurve, vnfmFpLog, DRL_TCURVE_FMT_PERCENT));
#endif

	/* get model parameyters */
	if (VnfmWrapRead(&tfcalData,
			backboneqL, betaL, alphaL,
			dateL, sigmaL, rhoL,
			zcCurve)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(tfcalData, vnfmFpLog));
#endif
	/* check two-factor */
	if (tfcalData->fNf != 2) {
	    GtoErrMsg("%s: must be 2-factor (got %d factors)\n",
		routine, tfcalData->fNf);
	    goto done;
	}
	if (VnfmComputeCoeff(tfcalData) != 0) goto done;



	/*
	 *
	 */
	WRAP_CHECK_VECTOR_LEN(rateReset1L,tfcalData->fNDates);
	WRAP_CHECK_VECTOR_LEN(rateMat1L,  tfcalData->fNDates);
	WRAP_CHECK_VECTOR_LEN(rateFreq1L, tfcalData->fNDates);
	WRAP_CHECK_VECTOR_LEN(rateVol1L,  tfcalData->fNDates);
	WRAP_CHECK_VECTOR_LEN(rateReset2L,tfcalData->fNDates);
	WRAP_CHECK_VECTOR_LEN(rateMat2L,  tfcalData->fNDates);
	WRAP_CHECK_VECTOR_LEN(rateFreq2L, tfcalData->fNDates);
	WRAP_CHECK_VECTOR_LEN(rateVol2L,  tfcalData->fNDates);

	WRAP_CHECK_VECTOR_LEN(outSigma1L,  tfcalData->fNDates);
	WRAP_CHECK_VECTOR_LEN(outSigma2L,  tfcalData->fNDates);

	/* convert vol frequencies from long to int */
	ASSERT_OR_DONE((freq1 = NEW_ARRAY(int, tfcalData->fNDates)) != NULL);
	ASSERT_OR_DONE((freq2 = NEW_ARRAY(int, tfcalData->fNDates)) != NULL);
	for (i=0; i<=tfcalData->fNDates-1; i++) {
		freq1[i] = (int) rateFreq1L[i+1];
		freq2[i] = (int) rateFreq2L[i+1];
	}


	/*
	 *
	 */
	WRAP_CHECK_VECTOR_LEN(calDatesL, 4);
	DrlTDateArrayFloorIdx(tfcalData->fDate, tfcalData->fNDates,
		calDatesL[1], &idxStart1);
	DrlTDateArrayFloorIdx(tfcalData->fDate, tfcalData->fNDates,
		calDatesL[2], &idxStart2);
	DrlTDateArrayFloorIdx(tfcalData->fDate, tfcalData->fNDates,
		calDatesL[3],   &idxEnd2);
	DrlTDateArrayFloorIdx(tfcalData->fDate, tfcalData->fNDates,
		calDatesL[4],   &idxEnd1);

	
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s:\n", routine);\
	DrlFPrintf(vnfmFpLog, "\tcalDatesL[1]=%s idx=%2d\n",\
			GtoFormatDate(calDatesL[1]), idxStart1);\
	DrlFPrintf(vnfmFpLog, "\tcalDatesL[2]=%s idx=%2d\n",\
			GtoFormatDate(calDatesL[2]), idxStart2);\
	DrlFPrintf(vnfmFpLog, "\tcalDatesL[3]=%s idx=%2d\n",\
			GtoFormatDate(calDatesL[3]),   idxEnd2);
	DrlFPrintf(vnfmFpLog, "\tcalDatesL[4]=%s idx=%2d\n",\
			GtoFormatDate(calDatesL[4]),   idxEnd1);)
	GTO_IF_LOGGING(\
	for (i=0; i<=tfcalData->fNDates-1; i++) {\
	    DrlFPrintf(vnfmFpLog, "\t[%3d]"\
		" rate1: reset=%10s mat=%8.4f freq=%d vol=%8.4f%%|"\
		" rate2: reset=%10s mat=%8.4f freq=%d vol=%8.4f%%\n",\
		i,\
		GtoFormatDate(rateReset1L[i+1]),\
		rateMat1L[i+1], freq1[i], rateVol1L[i+1]*1e2,\
		GtoFormatDate(rateReset2L[i+1]),\
		rateMat2L[i+1], freq2[i], rateVol2L[i+1]*1e2);\
	})
#endif

	/*
	 * (1) 1-vol calib
	 */
	if (VnfmVolCalib1VArbitrary(
		tfcalData,	
		idxStart1, idxStart2-1,
		&rateMat1L[1], freq1, &rateVol1L[1],
		LOGVOL,
		TRUE, NULL) != SUCCESS)
			goto done;



	/*
	 * (2) 2-vol calib
	 */
	if (VnfmVolCalib2VArbitrary(
		tfcalData,	
		idxStart2,
		idxEnd2,
		&rateReset1L[1], &rateMat1L[1], freq1, &rateVol1L[1],
		&rateReset2L[1], &rateMat2L[1], freq2, &rateVol2L[1],
		TRUE, NULL) != SUCCESS)
			goto done;

	/* rescale spot vol2 to have equal last 
	 * calibrated spot vols */
	ratioVol = tfcalData->fSigma[0][idxEnd2] /
			tfcalData->fSigma[1][idxEnd2];
	for (i=0; i<=tfcalData->fNDates-1; i++) {
		tfcalData->fSigma[1][i] *= ratioVol;
	}
	tfcalData->fAlpha[1] /= ratioVol;


	/*
	 * (3) 1-spot vol calibration
	 */
	errCode = VnfmVolCalib1VArbitrary(
		tfcalData,	
		idxEnd2+1,
		idxEnd1,
		&rateMat1L[1], freq1, &rateVol1L[1],
		LOGVOL,
		TRUE, NULL);

	/* rescale back */
	for (i=0; i<=tfcalData->fNDates-1; i++) {
	    	tfcalData->fSigma[1][i] /= ratioVol;
	}
	tfcalData->fAlpha[1] *= ratioVol;
	if (errCode != 0) goto done;


	/*
	 * wrap out the calibrated spot volatilities
	 */

	for (i=0; i<=tfcalData->fNDates-1; i++) {
		outSigma1L[i+1] = tfcalData->fSigma[0][i];
		outSigma2L[i+1] = tfcalData->fSigma[1][i];
	}


	/* log outputs */
	if (GtoLoggingGet() > 0) {
	    DrlLilVectLoggingFile("wrapper.log", "a", "TF_SMOOTH_SWMAT: Output",
		DRL_FLOAT_L, outSigma1L, "OUT_SIGMA1",
		DRL_FLOAT_L, outSigma2L, "OUT_SIGMA2",
		0);
	}


	/* made it through OK */
	status = SUCCESS;
done:

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: final parameters\n", routine);\
	VnfmFpWrite(tfcalData, vnfmFpLog));
#endif

	/* free allocated memory */
	GtoFreeTCurve(zcCurve);
	VnfmFree(tfcalData);
	if (freq1 != NULL) FREE((void*) freq1);
	if (freq2 != NULL) FREE((void*) freq2);


	VNFM_LOGCLOSE

	return (status);
}





