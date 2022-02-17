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
#include <ctype.h>

#include "convert.h"
#include "yearfrac.h"
#include "date_sup.h"

#include "drlmem.h"
#include "drlvtype.h"			/* DrlLilStructGet() */
#include "drlsmat.h"
#include "drlts.h"			/* DrlTCurveWrap() */
#include "drlio.h"			/* DrlFPrintf() */
#include "gtomat.h"

#define	_vnfm_SOURCE
#include "vnfmanly.h"
#include "vnfmcali.h"
#include "vnfmwrap.h"

#if defined(_WINDLL) 
# undef __DEBUG__
#endif


/* 
 * Convert vol type from string to KVolType
 */
int
ConvertVolType(
	char *volTypeStr,
	KVolType *volType);




/*---------------------------------------------------------------
 * <b> Add-in Function Name:</b>
 *                                                             
 * <br><br>
 * This is an old function for backward compatible assuming volatilities 
 * are all lognormal. Call the new function with lognormal vol input.
 */

DLL_EXPORT(int)
VnfmCalib1V2FGeneralL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */
	
	FloatL *floatScalarsL,	/*  4 'F' (I) num scalars:
				 *        [1] back bone q
				 *        [2] generate base vol 
				 *        [3] generate swaption matrix 
				 *        [4] minimum volatility rate mat */
	FloatL *nfParamsL,	/*  5 'F' (I) N-fact params */
	TDateL *nfDatesL,	/*  6 'D' (I) array of dates */

	FloatL *rateMatL,	/*  7 'F' (I) rate mat [0..nDates-1] */
	IntL *rateFreqL,	/*  8 'L' (I) rate freq [0..nDates-1] */
	FloatL *rateVolL,	/*  9 'F' (I) rate vol [0..nDates-1] */

	TDateIntervalL *bvMatL,	/* 10 'F' (I) vol maturity */
	TDateL *bvDatesL,	/* 11 'D' (I) vol dates */
	IntL *swTypeL,		/* 12 'L' (I) [0]=type, [1]=freq */
	TDateIntervalL* swMatL,	/* 13 'F' (I) array of mat intervals */
	TDateIntervalL* swExpL,	/* 14 'F' (I) array of exp intervals */

	FloatL *spotVolsL,	/* 15 'F' (O) spot vol */
	FloatL *bvRatesL,	/* 16 'F' (O) base vol */
	FloatL *swVolL)		/* 17 'F' (O) swaption vol matrix */
{
static	char		routine[] = "VnfmCalib1V2FGeneralL";
	int		status = FAILURE;

	char	inVolTypeL[WRAP_STR_IDX(2)], outVolTypeL[WRAP_STR_IDX(2)];
	double	numScalarsL[8];

	/*
	 * Lognormal volatility
	 */
	inVolTypeL[0] = (char)1;
	inVolTypeL[WRAP_STR_IDX(1)] = 'L';
	inVolTypeL[WRAP_STR_IDX(1)+1] = '\0';
	outVolTypeL[0] = (char)1;
	outVolTypeL[WRAP_STR_IDX(1)] = 'L';
	outVolTypeL[WRAP_STR_IDX(1)+1] = '\0';

	/*
	 * Lognormal distribution
	 */
	numScalarsL[0] = 7;			/* size */
	numScalarsL[1] = floatScalarsL[1];	/* backboneq */
	numScalarsL[2] = 0.0;			/* qL */
	numScalarsL[3] = 0.0;			/* qR */
	numScalarsL[4] = 0.0;			/* Fsh */
	numScalarsL[5] = floatScalarsL[2];	/* base vol flag */
	numScalarsL[6] = floatScalarsL[3];	/* swaption vol flag */
	numScalarsL[7] = floatScalarsL[4];	/* min maturity */
	
	status =  VnfmCalib1V2FGeneralNewL(
			refDateL,
			zcDateL,	
			zcRateL,
			numScalarsL,
			nfParamsL,
			nfDatesL,
			inVolTypeL, 
			rateMatL,	
			rateFreqL,
			rateVolL,

			outVolTypeL, 
			bvMatL,	
			bvDatesL,
			swTypeL,
			swMatL,	
			swExpL,	

			spotVolsL,
			bvRatesL,
			swVolL);

	return status;
}



/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CAL2F1SPVOL_GEN.
 *
 * <br><br>
 * Wrapper for routine <i> VnfmCalib1V2FGeneralL</i>
 * that simultaneously performs spot volatility calibration 
 * of arbitrary volatity points
 * and generate implied base volatility curve and model swaption matrix.
 * <br>
 * <br>[refDateL] zero curve value date.
 * <br>[zcDateL] zero curve coupon dates.
 * <br>[zcRateL] zero curve coupon rates.
 * <br>[floatScalarsL] numerical  scalars:\\
 *	(1) back bone $q$ (0 for log-normal, 0.5 for normal),\\
 *	(2) smile parameter 1: $q_L$  (0 for log-normal, 1.0 for normal),\\
 *	(3) smile parameter 2: $q_R$  (0 for log-normal, 1.0 for normal),\\
 *	(4) smile parameter 3: $Fsh$  forward shift,\\
 *	(5) generate base vol (0 or 1),\\
 *	(6) generate swaption matrix(0 or 1), \\
 *	(7) minimum volatility rate maturity($&gt;= 0$).
 * <br>[nfParamsL] LIL array containg the mean reversion, weight and
 * correlation coefficients.
 * Tree formats are available: an array format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n}),$$
 * that has ${n(n+3)/ 2}$ elements,
 * or a ``short'' array format where the 1st weight,
 * assumed to be 1.0, is omitted 
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_2,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n})$$
 * that has ${n(n+3)/ 2}-1$ elements,
 * or a matrix format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,1},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n,n})$$
 * that has $n(n+2)$ elements,
 * {\bf Remark: because of possible ambiguity in the format
 * (e.g. 8 elements can correspond to $n=3$ of the short array format
 * or $n=2$ of the matrix form) the format is checked in 
 * the previous order).}
 * <br>[nfDatesL] array of dates used for the spot volatility
 * bootstrapping (the spot volatility is assumed to be constant between
 * two consecutive dates).\\
 * {\bf WARNING: Unlike the ALIB routines, the convention is that
 * the spot volatility $i$ applies between date $i$ and date $i+1$.
 * The first date in the array MUST be the today date}.
 * <br>[rateMatL] array of volatility rate maturity (in yrs).
 * The element of index $i$ in the array corresponds
 * to an option expiring at date $i$ of array <i> nfDatesL</i>.
 * This array must have same length as the array <i> nfDatesL</i>.
 * {\bf
 * <br>
 *  <br> The volatility 1 (corresponding to the date 1, i.e. today)
 *     is never calibrated.
 * <br> The calibration stops at the first negative maturity encountered
 *     according to the following rules:
 *     <br>
 *     <br> If the maturity 1 (corresponding to the date 1, i.e. today)
 *         is negative, an error is returned.
 *     <br> If the maturity 1 is positive,
 *         but maturity 2 is negative, then volatility 2 is
 *         nevertheless calibrated assuming a maturity 0.25.
 *         No other volatility point is calibrated.
 *     <br> Otherwise, the volatility points are calibrated until a negative
 *         maturity is encountered.
 * <br>
 * <br> Any maturity smaller than floatScalarsL[4], but positive,
 *     will be rounded to floatScalarsL[4] (and the corresponding volatility 
 *     will be calibrated as such).
 * <br>
 * }
 * <br>[inVTypeL] input volatility type (Lognormal, Normal).
 * <br>[rateFreqL] array of volatility frequency (0,1,2,4,12).
 * {\it The element of index $i$ in the array corresponds
 * to an option expiring at date $i$ of array <i> nfDatesL</i>}.
 * <br>[rateVolL] array of volatilities.
 * {\it The element of index $i$ in the array corresponds
 * to an option expiring at date $i$ of array <i> nfDatesL</i>}.
 * <br>[outVTypeL] output volatility type (Lognormal, Normal).
 * <br>[bvMatL] base volatility underlying rate  maturity (in yrs)
 * for the output <i> bvRatesL</i>.
 * <br>[bvDatesL] base volatility expiration dates
 * for the output <i> bvRatesL</i>.
 * <br>[swTypeL] array of 2 elements:\\
 * 	(1) swaption matrix output type (0 for vertical, 1 for diagonal)
 *		for the output <i> SwVolL</i>.\\
 * 	(2) swaption matrix volatility frequency (0,1,2,4,12)
 *		for the output <i> SwVolL</i>.
 * <br>[swMatL] array of maturity intervals (in yrs)
 *		for the output <i> SwVolL</i>.
 * <br>[swExpL] array of expration intervals (in yrs)
 *		for the output <i> SwVolL</i>.
 * <br>[spotVolsL] output spot volatility corresponding to
 * array of dates <i> nfDatesL</i>.
 * {\bf WARNING: Unlike the ALIB routines, the convention is that
 * the spot volatility $i$ applies between date $i$ and date $i+1$.
 * Therefore, the volatility dates sould be offset by one when
 * passing the curve to an ALIB tree routine}.
 * <br>[bvRatesL] output base volatility curve
 * specified by arguments <i> bvMatL</i> and <i> bvDatesL</i>
 * <br>[swVolL] output swaption matrix
 * specified by arguments <i> swTypeL</i>, <i> swMatL</i> and <i> swExpL</i>.
 * <br>
 */

DLL_EXPORT(int)
VnfmCalib1V2FGeneralNewL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */
	
	FloatL *floatScalarsL,	/*  4 'F' (I) num scalars:
				 *        [1] back bone q
				 *        [2] qL (1=normal,0=lognormal) 
				 *        [3] qR (1=normal,0=lognormal) 
				 *        [4] Fsh 
				 *        [5] generate base vol 
				 *        [6] generate swaption matrix 
				 *        [7] minimum volatility rate mat */
	FloatL *nfParamsL,	/*  5 'F' (I) N-fact params */
	TDateL *nfDatesL,	/*  6 'D' (I) array of dates */


	CharBlockL *inVTypeL, 	/*  7 'L' (I) input vol type:LOG, NORM */
	FloatL *rateMatL,	/*  8 'F' (I) rate mat [0..nDates-1] */
	IntL *rateFreqL,	/*  9 'L' (I) rate freq [0..nDates-1] */
	FloatL *rateVolL,	/* 10 'F' (I) rate vol [0..nDates-1] */

	CharBlockL *outVTypeL, 	/* 11 'L' (I) output vol type:LOG, NORM */
	TDateIntervalL *bvMatL,	/* 12 'F' (I) vol maturity */
	TDateL *bvDatesL,	/* 13 'D' (I) vol dates */
	IntL *swTypeL,		/* 14 'L' (I) [0]=type, [1]=freq */
	TDateIntervalL* swMatL,	/* 15 'F' (I) array of mat intervals */
	TDateIntervalL* swExpL,	/* 16 'F' (I) array of exp intervals */

	FloatL *spotVolsL,	/* 17 'F' (O) spot vol */
	FloatL *bvRatesL,	/* 18 'F' (O) base vol */
	FloatL *swVolL)		/* 19 'F' (O) swaption vol matrix */
{
static	char		routine[] = "VnfmCalib1V2FGeneralNewL";
	int		status = FAILURE;

/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
/*	NEW VERSION						*/
/*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

	VnfmData	*vnfmData = NULL;
	TCurve		*zcCurve = NULL;
	TDateInterval	bvMat;
	int		nDates, i;

	int		idxStart = 0,
			idxEnd;
	double		*rateMat = &rateMatL[1],
			*rateVol = &rateVolL[1]; 
	int		*rateFreq = NULL;

	int		swNExp, swNMat, idxE, idxM;
	TDateInterval	*swLExp = NULL, *swLMat = NULL;
	double		**swVol = NULL;

	KVolType	inVType, outVType;


	VNFM_LOGOPEN

	/* check arg len */
	WRAP_CHECK_SCALAR(refDateL);
	WRAP_CHECK_VECTOR(zcDateL);
	WRAP_CHECK_VECTOR(zcRateL);
	WRAP_CHECK_VECTOR_LEN(floatScalarsL, 7);
	WRAP_CHECK_VECTOR(nfParamsL);
	WRAP_CHECK_VECTOR(nfDatesL);

	WRAP_CHECK_SCALAR(inVTypeL);
	WRAP_CHECK_VECTOR(rateMatL);
	WRAP_CHECK_VECTOR(rateFreqL);
	WRAP_CHECK_VECTOR(rateVolL);

	WRAP_CHECK_SCALAR(outVTypeL);
	WRAP_CHECK_VECTOR_LEN(bvMatL, 1);
	WRAP_CHECK_VECTOR(bvDatesL);
	WRAP_CHECK_VECTOR_LEN(swTypeL, 2);
	WRAP_CHECK_VECTOR(swMatL);
	WRAP_CHECK_VECTOR(swExpL);

	WRAP_CHECK_VECTOR(spotVolsL);
	WRAP_CHECK_VECTOR(bvRatesL);
	WRAP_CHECK_VECTOR(swVolL);


	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDateL, zcRateL) != SUCCESS)
		goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve,\
		vnfmFpLog, DRL_TCURVE_FMT_PERCENT));
#endif

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "bvMat=%s\n",\
		GtoFormatDateInterval(&bvMat)));
#endif

	nDates = ARGSIZE(nfDatesL);
	ASSERT_OR_DONE(ARGSIZE(rateMatL) == ARGSIZE(nfDatesL));
	ASSERT_OR_DONE(ARGSIZE(rateFreqL) == ARGSIZE(nfDatesL));
	ASSERT_OR_DONE(ARGSIZE(rateVolL) == ARGSIZE(nfDatesL));
	if ((rateFreq = NEW_ARRAY(int, ARGSIZE(rateFreqL))) == NULL)
		goto done;
	for (i=0; i<=ARGSIZE(rateFreqL)-1; i++)
		rateFreq[i] = (int) rateFreqL[i+1];



	/**
	 ** Create Vnfm structure
	 **/
	if (VnfmWrapReadSimple(
		&vnfmData,
		floatScalarsL,
		nfParamsL,
		nfDatesL,
		zcCurve) != SUCCESS)
			goto done;

	
	/**
	 ** Set smile
	 **/
	if (VnfmSetSmile(
		vnfmData,
		1.0 - 2.0 * floatScalarsL[1],	/* bbq: 1=lognormal, 0=normal */
		1.0 - floatScalarsL[2],		/* qL:  1=lognormal, 0=normal */
		1.0 - floatScalarsL[3],		/* qR:  1=lognormal, 0=normal */
		floatScalarsL[4]) != SUCCESS)
			goto done;

	/* 
	 * To avoid conflict of spot vol calibration combined with
	 * other routine, we default the backbone reference vols
	 * to be 1.
	 */
	vnfmData->fVolNorm = 1e0;
	vnfmData->fVolLogn = 1e0;
	 
	/** 
	 ** Check the cut-off rate maturity 
	 **/
	if (floatScalarsL[7] < 0e0) {
	    GtoErrMsg("%s: the cut-off vol rate maturity is negative(%lf).\n",
                routine, floatScalarsL[7]);
            goto done;
        }

	/**
	 ** Adjust array maturity 
	 **/
	/* Set last calibrated point to be the 1st negative maturity
	 * (or the 1st less than 1w that is not the 1st point).
	 */
	if (rateMat[0] < 0e0) {
	    GtoErrMsg("%s: first maturity found in calibrated "
		"rate array is negative (%lf).\n", routine, rateMat[0]);
	    goto done;
	}

	/* treat the special case of maturity between today (date idx 0)
	 * and next vol point (date idx 1)
	 */
	if (rateMat[1] < 0e0) {
		idxEnd = 1;
	} else {
		for (idxEnd=1; idxEnd<=nDates-1; idxEnd++) {
		    if (rateMat[idxEnd] < 0e0) 
			break;
		}
		idxEnd--;
	}

	/* Cutoff for small maturities */
	for (i=0; i<=idxEnd; i++) {
		rateMat[i] = MAX(rateMat[i], floatScalarsL[7]);
	}


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: idxEnd=%d "\
		"(last idx calibrated).\n", routine, idxEnd))
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "\tInput vols:\n"\
		"\tidx         dates     rateMat rateFreq  rateVol\n"))
	GTO_IF_LOGGING(\
	for (i=0; i<=nDates-1; i++) {\
	    DrlFPrintf(vnfmFpLog, "\t[%2d/%2d]  %10s  %8.4f  %d  %8.4f%%\n",\
		i, nDates,
		GtoFormatDate(nfDatesL[i+1]),
		rateMat[i], rateFreq[i], rateVol[i]*1e2);\
	})
#endif

	if (idxEnd < 0) {
	    GtoErrMsg("%s: no positive maturity found in calibrated "
		"rate array.\n", routine);
	    goto done;
	}



	/* 
	 * Convert vol type from string to KVolType
	 */
	if (ConvertVolType(&inVTypeL[WRAP_STR_IDX(1)], &inVType)   != SUCCESS ||
	    ConvertVolType(&outVTypeL[WRAP_STR_IDX(1)], &outVType) != SUCCESS)
		goto done;


	/* perform the calibration */
	if (VnfmVolCalib1VArbitrary(
		vnfmData,
		idxStart,
		idxEnd,
		rateMat,
		rateFreq,
		rateVol,
		inVType,
		TRUE, NULL) != SUCCESS)
			goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: calibrated:\n", routine))
	GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog))
#endif

	/* output spot volatilities */
	if ((int) ARGSIZE(spotVolsL) < vnfmData->fNDates) {
		GtoErrMsg("%s: output array of spot vols not large enough "
			"(got %d, timeline length is %d).\n", routine,
	    		(int)ARGSIZE(spotVolsL), vnfmData->fNDates);
		goto done;
	}
	for(i=0; i<=vnfmData->fNDates-1; i++) {
		spotVolsL[i+1] = vnfmData->fSigma[0][i];
	}




	/**
	 ** Generate base vol curve
	 **/
	if ((long)floatScalarsL[5] != 0L) {
		int	n1, n2; 	/* tmp variables */

		/* get desired base vola maturity */
		if (GtoYearsToDateInterval(bvMatL[1], &bvMat) != SUCCESS)
			goto done;
#ifndef	NO_LOGGING
		GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "bvMat=%s\n",\
			GtoFormatDateInterval(&bvMat)));
#endif

		/* Check size consistency */
		n1 = DrlLilVectActiveSize(DRL_TDATE_L,  bvDatesL, "bvDatesL");
		n2 = DrlLilVectActiveSize(DRL_DOUBLE_L, bvRatesL, "bvRatesL");
		if (n1 != n2) {
			GtoErrMsg("%s: num base vol dates (%d) !="
				" num base vol rates (%d).\n", routine, n1, n2);
			goto done;
		}

		/* Generate curve */
		if (VnfmVolCurve(
			vnfmData,
			bvMat,
			0,	/* simple rate */
			0,	/* base vol */
			vnfmData->fDate[0], /* bv ref date */
			ARGSIZE(bvDatesL),
			&bvDatesL[1],
			outVType,
			&bvRatesL[1]) != SUCCESS)
				goto done;

	}

	/**
	 ** Generate swaption volatilities
	 **/
	if ((long)floatScalarsL[6] != 0L) {


		WRAP_CHECK_VECTOR_LEN(swTypeL, 2);
		WRAP_CHECK_VECTOR(swExpL);
		WRAP_CHECK_VECTOR(swMatL);
		swNExp = (int) ARGSIZE(swExpL);
		swNMat = (int) ARGSIZE(swMatL);

		if (ARGSIZE(swVolL) != swNExp * swNMat) {
			GtoErrMsg("%s: swVolL size (%d) != "
				"swNExp (%d) * swNMat (%d).\n",
				routine, ARGSIZE(swVolL), swNExp, swNMat);
			goto done;
		}

		if ((swLExp = DrlTDateIntervalVectAlloc(0, swNExp-1))
			== NULL) goto done;
		for (idxE=0; idxE<=swNExp-1;idxE++)
		    if (GtoYearsToDateInterval(swExpL[idxE+1], &swLExp[idxE])
			!= SUCCESS) goto done;

		if ((swLMat = DrlTDateIntervalVectAlloc(0, swNMat-1))
			== NULL) goto done;
		for (idxM=0; idxM<=swNMat-1;idxM++)
		    if (GtoYearsToDateInterval(swMatL[idxM+1], &swLMat[idxM])
			!= SUCCESS) goto done;

		if ((swVol = DrlDoubleMatrAlloc(0, swNExp-1, 0, swNMat-1))
			== NULL) goto done;

		if (VnfmSwaptionVolMatrix(
			vnfmData,
			(int) swTypeL[1],
			(int) swTypeL[2],
			FALSE,	/* no adjustment */
			swNExp, swLExp,
			swNMat, swLMat,
			outVType,
			swVol) != SUCCESS)
				goto done;

		/* */
		for (idxM=0; idxM<=swNMat-1;idxM++)
		for (idxE=0; idxE<=swNExp-1;idxE++) {
			swVolL[1+idxM+idxE*swNMat] = swVol[idxE][idxM];
		}

	}



	/* made it through OK */
	status = SUCCESS;
done:

	DrlTDateIntervalVectFree(swLExp, 0, swNExp-1);
	DrlTDateIntervalVectFree(swLMat, 0, swNMat-1);
	DrlDoubleMatrFree(swVol, 0, swNExp-1, 0, swNMat-1);


	FREE(rateFreq);
	VnfmFree(vnfmData);
	GtoFreeTCurve(zcCurve);

	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}
	VNFM_LOGCLOSE
	return(status);

}




/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CAL1SPVOL.
 *                                                             
 * <br><br>
 * Wrapper for single spot volatility calibration (arbitrary number
 * of factors and arbitrary volatility specification).
 * <br>
 * <br>[refDateL] zero curve value date.
 * <br>[zcDateL] zero curve coupon dates.
 * <br>[zcRateL] zero curve coupon rates.
 * <br>[backboneqL] array of 1 double:\\
 *	(1) backboneq (0 for lognormal, 0.5 for normal).
 * <br>[betaL] array of mean-reversion coefficients (as may as factors).
 * <br>[alphaL] array of factor weights (as may as factors).
 * <br>[dateL] array of dates used for the spot volatility
 * (the spot volatility is assumed to be constant between
 * two consecutive dates).\\
 * {\bf WARNING: Unlike the ALIB routines, the convention is that
 * the spot volatilities and correlations of index $i$
 * apply between date $i$ and date $i+1$.
 * The first date in the array MUST be the volatility reference date (today)}.
 * <br>[sigmaL] range of spot volatilities, indexed by dates in row
 * and factors in column. {\bf These volatilities are
 * currently all overwritten, the argument being here for future use.
 * Pass a dummy input value of 0 for all elements in the array}.
 * <br>[rhoL] range of correlations, indexed by dates in in row
 * and by factor in columns in the following way:
 * the correlations $(\rho_{ij}(t))_{1&lt;= i<j &lt;= N}$
 * between factors $i$ and $j$ is stored in the
 * column of index $k$ (where $k=0,\dots\,N-1$) where
 * $k = {N(N-1)/ 2} - {(N-i)(N-i-1)/ 2} + (j-i-1)$
 * with  $0&lt;= k &lt;= {N(N-1)/ 2}-1$.
 * <br>[rateMatL] array of maturities (in years) of rates
 *      to be calibrated (corresponding to the dates
 *      in the timeline <i> dateL</i>).
 *      <b> The calibration stops at the first negative maturity encountered</b>.
 *      This array must have same length as the array {\t dateL}.
 * <br>[rateFreqL] array of rates frequencies (1,2,4,12).
 *      This array must have same length as the array {\t dateL}.
 * <br>[rateVolL] array of rates volatilities.
 *      This array must have same length as the array {\t dateL}.
 * <br>[tStartL] Not curently used. Pass a dummy calue of 0.
 * <br>[oSigmaL] Output aspot volatility array (should have same length
 *      as <i> dateL</i>.
 * <br>
 */

DLL_EXPORT(int)
VnfmCalib1VVolCurveNewL(
	TDateL *refDateL,	/*  1 'D' (I) zero coupon value date */
	TDateL *zcDatesL,	/*  2 'D' (I) array of zero coupon dates */
	FloatL *zcRatesL,	/*  3 'F' (I) array of zero coupon rates */
				/*        Input Model: */
 	FloatL *backboneqL,	/*  4 'F' (I) [1] 0=lognormal, 0.5=normal */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlations */
				/*        Input base vol curve: */
	double *rateMatL,	/* 10 'F' (I) array of vol  mat [0..nDates-1] */
	IntL *rateFreqL,	/* 11 'L' (I) array of vol freq [0..nDates-1] */
	FloatL *rateVolL,	/* 12 'F' (I) array of vol [0..nDates-1] */
	FloatL *tStartL,	/* 13 'F' (I) start/end time for calibration */
				/*        Output Calibrated Model: */
	FloatL *oSigmaL)	/* 14 'F' (O) volatility vector (1-D) */
{
static	char		routine[] = "VnfmCalib1VVolCurveNewL";
	int		status = FAILURE;
	VnfmData	*vnfmData = NULL;
	TCurve		*zcCurve = NULL;
	int		idx,
			idxStart = 0,
			idxEnd;
	int		*rateFreq = NULL;
	double		*rateMat = NULL,
			*rateVol = NULL;


	VNFM_LOGOPEN



	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDatesL, zcRatesL)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,
		DRL_TCURVE_FMT_PERCENT));
#endif

	/* get model parameters */
	if (VnfmWrapRead(&vnfmData,
			backboneqL, betaL, alphaL,
			dateL, sigmaL, rhoL,
			zcCurve)
		!= SUCCESS) goto done;

	/* compute coeffs */
	if (VnfmComputeCoeff(vnfmData) != SUCCESS)
		goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog);\
	VnfmPrintCoeff(vnfmData, vnfmFpLog));
#endif

	/* Check size of arrays is consistent */
	ASSERT_OR_DONE(ARGSIZE(rateMatL)   == vnfmData->fNDates);
	ASSERT_OR_DONE(ARGSIZE(rateFreqL)  == vnfmData->fNDates);
	ASSERT_OR_DONE(ARGSIZE(rateVolL)   == vnfmData->fNDates);
	ASSERT_OR_DONE(ARGSIZE(oSigmaL)    == vnfmData->fNDates);

	if ((rateMat  = DrlDoubleVectAlloc(0, vnfmData->fNDates-1)) == NULL)
		goto done;
	if ((rateVol  = DrlDoubleVectAlloc(0, vnfmData->fNDates-1)) == NULL)
		goto done;
	if ((rateFreq = DrlIntVectAlloc(0, vnfmData->fNDates-1)) == NULL)
		goto done;
	for (idx=0; idx<=vnfmData->fNDates-1; idx++) {
		rateMat[idx]  = rateMatL[idx+1];
		rateFreq[idx] = (int) rateFreqL[idx+1];
		rateVol[idx]  = rateVolL[idx+1];
	}


	/* find last cali time = 1st negative maturity */
	for (idxEnd=1; idxEnd<=vnfmData->fNDates-1; idxEnd++) {
	    if (rateMat[idxEnd] <= 0e0) {
		break;
	    }
	}
	idxEnd--;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: idxStart=%d idxEnd=%d.\n",\
		routine, idxStart, idxEnd))
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "\tInput vols:\n"\
		"\tidx\tmatt\tfreq\tvol\n"))
	GTO_IF_LOGGING(\
	for (idx=0; idx<=vnfmData->fNDates-1; idx++) {
	    DrlFPrintf(vnfmFpLog, "\t[%2d]\t%s\t%7.4f\t%d\t%8.4f%%\n",\
		idx, GtoFormatDate(vnfmData->fDate[idx]),\
		rateMat[idx], rateFreq[idx], rateVol[idx]*1e2);\
	})
#endif



	/* perform the calibration */
	if (VnfmVolCalib1VArbitrary(
		vnfmData,
		idxStart,
		idxEnd,
		rateMat,
		rateFreq,
		rateVol,
		LOGVOL,
		TRUE, NULL) != SUCCESS)
			goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: calibrated:\n", routine))
	GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog))
#endif

	/* write vol to output */
	for (idx=0; idx<=vnfmData->fNDates-1; idx++) {
		oSigmaL[idx+1] = vnfmData->fSigma[0][idx];
	}


	/* made it through OK */
	status = SUCCESS;
done:
	DrlDoubleVectFree(rateMat, 0, (vnfmData!=NULL ? vnfmData->fNDates-1:0));
	DrlDoubleVectFree(rateVol, 0, (vnfmData!=NULL ? vnfmData->fNDates-1:0));
	DrlIntVectFree(rateFreq, 0, (vnfmData!=NULL ? vnfmData->fNDates-1:0));
	GtoFreeTCurve(zcCurve);
	VnfmFree(vnfmData);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
}



/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CAL1SPVOL_ARB.
 *                                                             
 * <br><br>
 * Wrapper for single spot volatility calibration (arbitrary number
 * of factors, arbitrary volatility specification and arbitrary resets).
 * <br>
 * <br>[refDateL] zero curve value date.
 * <br>[zcDateL] zero curve coupon dates.
 * <br>[zcRateL] zero curve coupon rates.
 * <br>[backboneqL] array of 1 double:\\
 *	(1) backboneq (0 for lognormal, 0.5 for normal).
 * <br>[betaL] array of mean-reversion coefficients (as may as factors).
 * <br>[alphaL] array of factor weights (as may as factors).
 * <br>[dateL] array of dates used for the spot volatility
 * (the spot volatility is assumed to be constant between
 * two consecutive dates).\\
 * {\bf WARNING: Unlike the ALIB routines, the convention is that
 * the spot volatilities and correlations of index $i$
 * apply between date $i$ and date $i+1$.
 * The first date in the array MUST be the volatility reference date (today)}.
 * <br>[sigmaL] range of spot volatilities, indexed by dates in row
 * and factors in column. {\bf These volatilities are
 * currently all overwritten, the argument being here for future use.
 * Pass a dummy input value of 0 for all elements in the array}.
 * <br>[rhoL] range of correlations, indexed by dates in in row
 * and by factor in columns in the following way:
 * the correlations $(\rho_{ij}(t))_{1&lt;= i<j &lt;= N}$
 * between factors $i$ and $j$ is stored in the
 * column of index $k$ (where $k=0,\dots\,N-1$) where
 * $k = {N(N-1)/ 2} - {(N-i)(N-i-1)/ 2} + (j-i-1)$
 * with  $0&lt;= k &lt;= {N(N-1)/ 2}-1$.
 * <br>[rateResetL] array of rate reset dates.
 *      This array must have same length as the array {\t dateL}.
 * <br>[rateMatL] array of maturities (in years) of rates
 *      to be calibrated (corresponding to the dates
 *      in the timeline <i> dateL</i>).
 *      <b> The calibration stops at the first negative maturity encountered</b>.
 *      This array must have same length as the array {\t dateL}.
 * <br>[rateFreqL] array of rates frequencies (1,2,4,12).
 *      This array must have same length as the array {\t dateL}.
 * <br>[rateVolL] array of rates volatilities.
 *      This array must have same length as the array {\t dateL}.
 * <br>[tStartL] Not curently used. Pass a dummy calue of 0.
 * <br>[oSigmaL] Output aspot volatility array (should have same length
 *      as <i> dateL</i>.
 * <br>
 */

DLL_EXPORT(int)
VnfmCalib1VVolCurveArbL(
	TDateL *refDateL,	/*  1 'D' (I) zero coupon value date */
	TDateL *zcDatesL,	/*  2 'D' (I) array of zero coupon dates */
	FloatL *zcRatesL,	/*  3 'F' (I) array of zero coupon rates */
				/*        Input Model: */
 	FloatL *backboneqL,	/*  4 'F' (I) [1] 0=lognormal, 0.5=normal */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlations */
				/*        Input benchmark vol curve: */
	TDateL *rateResetL,	/* 10 'F' (I) array of rate reset[0..nDates-1]*/
	double *rateMatL,	/* 11 'F' (I) array of vol  mat [0..nDates-1] */
	IntL *rateFreqL,	/* 12 'L' (I) array of vol freq [0..nDates-1] */
	FloatL *rateVolL,	/* 13 'F' (I) array of vol [0..nDates-1] */
	FloatL *tStartL,	/* 14 'F' (I) start/end time for calibration */
				/*        Output Calibrated Model: */
	FloatL *oSigmaL)	/* 15 'F' (O) volatility vector (1-D) */
{
static	char		routine[] = "VnfmCalib1VVolCurveArbL";
	int		status = FAILURE;
	VnfmData	*vnfmData = NULL;
	TCurve		*zcCurve = NULL;
	int		idx,
			idxStart = 0,
			idxEnd;
	int		*rateFreq = NULL;
	double		*tReset = NULL,
			*rateMat = NULL,
			*rateVol = NULL;


	VNFM_LOGOPEN



	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDatesL, zcRatesL)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,
		DRL_TCURVE_FMT_PERCENT));
#endif

	/* get model parameters */
	if (VnfmWrapRead(&vnfmData,
			backboneqL, betaL, alphaL,
			dateL, sigmaL, rhoL,
			zcCurve)
		!= SUCCESS) goto done;

	/* compute coeffs */
	if (VnfmComputeCoeff(vnfmData) != SUCCESS)
		goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog);\
	VnfmPrintCoeff(vnfmData, vnfmFpLog));
#endif

	/* Check size of arrays is consistent */
	ASSERT_OR_DONE(ARGSIZE(rateResetL) == vnfmData->fNDates);
	ASSERT_OR_DONE(ARGSIZE(rateMatL)   == vnfmData->fNDates);
	ASSERT_OR_DONE(ARGSIZE(rateFreqL)  == vnfmData->fNDates);
	ASSERT_OR_DONE(ARGSIZE(rateVolL)   == vnfmData->fNDates);
	ASSERT_OR_DONE(ARGSIZE(oSigmaL)    == vnfmData->fNDates);

	if ((tReset  = DrlDoubleVectAlloc(0, vnfmData->fNDates-1)) == NULL)
		goto done;
	if ((rateMat  = DrlDoubleVectAlloc(0, vnfmData->fNDates-1)) == NULL)
		goto done;
	if ((rateVol  = DrlDoubleVectAlloc(0, vnfmData->fNDates-1)) == NULL)
		goto done;
	if ((rateFreq = DrlIntVectAlloc(0, vnfmData->fNDates-1)) == NULL)
		goto done;
	for (idx=0; idx<=vnfmData->fNDates-1; idx++) {
		if(GtoDayCountFraction(vnfmData->fDate[0], rateResetL[idx+1],
                        GTO_ACT_365F, &tReset[idx]) != SUCCESS)
			goto done;
		rateMat[idx]  = rateMatL[idx+1];
		rateFreq[idx] = (int) rateFreqL[idx+1];
		rateVol[idx]  = rateVolL[idx+1];
	}


	/* find last cali time = 1st negative maturity */
	for (idxEnd=1; idxEnd<=vnfmData->fNDates-1; idxEnd++) {
	    if (rateMat[idxEnd] <= 0e0) {
		break;
	    }
	}
	idxEnd--;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: idxStart=%d idxEnd=%d.\n",\
		routine, idxStart, idxEnd))
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "\tInput vols:\n"\
		"\tidx\tmatt\tfreq\tvol\n"))
	GTO_IF_LOGGING(\
	for (idx=0; idx<=vnfmData->fNDates-1; idx++) {
	    DrlFPrintf(vnfmFpLog, "\t[%2d]\t%s\t%7.4f\t%d\t%8.4f%%\n",\
		idx, GtoFormatDate(vnfmData->fDate[idx]),\
		rateMat[idx], rateFreq[idx], rateVol[idx]*1e2);\
	})
#endif



	/* perform the calibration */
	if (VnfmVolCalib1VArbitraryNew(
		vnfmData,
		idxStart,
		idxEnd,
		tReset,
		rateMat,
		rateFreq,
		rateVol,
		LOGVOL,
		TRUE, NULL) != SUCCESS)
			goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: calibrated:\n", routine))
	GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog))
#endif

	/* write vol to output */
	for (idx=0; idx<=vnfmData->fNDates-1; idx++) {
		oSigmaL[idx+1] = vnfmData->fSigma[0][idx];
	}


	/* made it through OK */
	status = SUCCESS;
done:
	DrlDoubleVectFree(tReset, 0, (vnfmData!=NULL ? vnfmData->fNDates-1:0));
	DrlDoubleVectFree(rateMat, 0, (vnfmData!=NULL ? vnfmData->fNDates-1:0));
	DrlDoubleVectFree(rateVol, 0, (vnfmData!=NULL ? vnfmData->fNDates-1:0));
	DrlIntVectFree(rateFreq, 0, (vnfmData!=NULL ? vnfmData->fNDates-1:0));
	GtoFreeTCurve(zcCurve);
	VnfmFree(vnfmData);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
}


/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CAL1SPVOL_OLD.
 *                                                             
 * <br><br>
 * Wrapper for single spot volatility calibration (arbitrary number
 * of factors).  Calls <i> VnfmCalib1VVolCurve</i>.
 */

DLL_EXPORT(int)
VnfmCalib1VVolCurveOldL(
	TDateL *refDateL,	/*  1 'D' (I) zero coupon value date */
	TDateL *zcDatesL,	/*  2 'D' (I) array of zero coupon dates */
	FloatL *zcRatesL,	/*  3 'F' (I) array of zero coupon rates */
				/*        Input Model: */
 	FloatL *backboneqL,	/*  4 'F' (I) [1] 0=lognormal, 0.5=normal */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlations */
				/*        Input base vol curve: */
	IntL *calTypeL,		/* 10 'L' (I) calibration type [0] */
	double *volMatL,	/* 11 'F' (I) vol maturity [0] */
	IntL *volFreqL,		/* 12 'L' (I) volatility frequency [0] */
	TDateL *volDatesL,	/* 13 'D' (I) array of vol dates */
	FloatL *volRates1L,	/* 14 'F' (I) array of vol values index # 1*/
	FloatL *tStartL,	/* 15 'F' (I) start/end time for calibration */
				/*        Output Calibrated Model: */
	FloatL *oSigmaL)	/* 16 'F' (O) volatility vector (1-D) */
{
static	char		routine[] = "VnfmCalib1VVolCurveOldL";
	int		status = FAILURE;
	VnfmData	*theModel = NULL;
	TCurve		*zcCurve = NULL,
			*volCurve1 =  NULL;
	int		calType[1];
	double		volMat[1];
	int		volFreq[1],
			n;
	double		tStart=0e0,
			tEnd=1e12,
			resValue;


	VNFM_LOGOPEN



	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDatesL, zcRatesL)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,
		DRL_TCURVE_FMT_PERCENT));
#endif

	/* get model parameyters */
	if (VnfmWrapRead(&theModel,
			backboneqL, betaL, alphaL,
			dateL, sigmaL, rhoL,
			zcCurve)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(theModel, vnfmFpLog));
#endif
	/* make sure size enough for output vol */
	ASSERT_OR_DONE(ARGSIZE(oSigmaL) == ARGSIZE(dateL));

	/* get benchmarks vols */
	if (DrlTCurveWrapRead(&volCurve1, refDateL, volDatesL, volRates1L)
		!= SUCCESS) goto done;


	if (DrlLilStructGet(0,
	    DRL_LILVECT_L, &n, 1, 1,
		"calType", DRL_INT_L,
			FALSE, (void*) calTypeL, (void*) calType,
	    DRL_LILVECT_L, &n, 1, 1,
		"volMat", DRL_FLOAT_L,
			FALSE, (void*) volMatL, (void*) volMat,
	    DRL_LILVECT_L, &n, 1, 1,
		"volFreq", DRL_INT_L,
			FALSE, (void*) volFreqL, (void*) volFreq,
	    0) != SUCCESS) goto done;

	volCurve1->fBasis = (volFreq[0] > 0 ? volFreq[0] :
			(int) (1e0/volMat[0]));

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: calType=%d volMat=%lf volFreq=%d\n",\
		routine, calType[0], volMat[0],\
		volFreq[0]);\
	DrlFPrintf(vnfmFpLog, "%s: input base volatility curve\n", routine);\
	DrlTCurveFpWrite(volCurve1, vnfmFpLog, DRL_TCURVE_FMT_PERCENT));
#endif

	/* start time */
	ASSERT_OR_DONE(ARGSIZE(tStartL) == 1);
	tStart = tStartL[1];
	tEnd   = 1e12;

	/* compute coeffs */
	if (VnfmComputeCoeff(theModel) != SUCCESS)
		goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(theModel, vnfmFpLog);\
	VnfmPrintCoeff(theModel, vnfmFpLog));
#endif


	/* perform the calibration */
	if (VnfmCalib1VVolCurve(
			theModel,
			calType[0], volMat[0], volFreq[0], volCurve1,
			tStart, tEnd,
			NULL,	/* reoprt all errors */
			&resValue)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(theModel, vnfmFpLog));
#endif

	/* write vol to output */
	for (n=0; n<=theModel->fNDates-1; n++) {
		oSigmaL[n+1] = theModel->fSigma[0][n];
	}


	/* made it through OK */
	status = SUCCESS;
done:
	GtoFreeTCurve(zcCurve);
	GtoFreeTCurve(volCurve1);
	VnfmFree(theModel);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
}


/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_FWDFUT.
 *                                                             
 * <br><br>
 * Wrapper for base volatility calibration (n-factor) and
 * forward/futures adjustment.
 * Calls <i> VnfmCalib1VVolCurve</i> and <i> VnfmFwdFutAdjustment</i>.
 */

DLL_EXPORT(int)
VnfmCalib1V2FVolCurveAndGenFutAdjL(
	TDateL *refDateL,		/*  1 'D' (I) reference date */
	TDateL *zcDateL,		/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,		/*  3 'F' (I) zero coupon rates */
	
	FloatL *nfParamsL,		/*  4 'F' (I) N-fact params */
	TDateL *nfDatesL,		/*  5 'D' (I) array of dates */

	FloatL *volMatL,		/*  6 'F' (I) vol maturity (cms) */
	TDateL *volDatesL,		/*  7 'D' (I) vol dates */
	FloatL *volRatesL,		/*  8 'F' (I) vol values */

	TDateIntervalL *rateMatL,	/*  9 'F' (I) rate maturity */
	IntL *rateFreqL,		/* 10 'L' (I) rate frequency */
	CharBlockL *rateDayCountL,	/* 11 'C' (I) rate day count */
	TDateL *resetDatesL,		/* 12 'D' (I) adj dates */
	FloatL *fwdRateL,		/* 13 'F' (O) forward rate */
	FloatL *cvxAdjL,		/* 14 'F' (O) cvx adjustment */
	FloatL *mtmAdjL)		/* 15 'F' (O) mtm adjustment */
{
static	char		routine[] = "VnfmCalib1V2FVolCurveAndGenFutAdjL";
	int		status = FAILURE;
	VnfmData	*theModel = NULL;
	TCurve		*zcCurve = NULL,
			*volCurve = NULL;
	int		rateFreq,
			volFreq,
			nDates,
			n;
	TDateInterval	rateMat;
	long		rateDayCount;
	double		volMat,
			resValue = 0e0;




	VNFM_LOGOPEN


	/* check arg len */
	WRAP_CHECK_SCALAR(refDateL);
	WRAP_CHECK_VECTOR(zcDateL);
	WRAP_CHECK_VECTOR(zcRateL);
	WRAP_CHECK_VECTOR_LEN(nfParamsL, 4);
	WRAP_CHECK_VECTOR(nfDatesL);

	WRAP_CHECK_VECTOR(volMatL);
	WRAP_CHECK_VECTOR(volDatesL);
	WRAP_CHECK_VECTOR(volRatesL);

	WRAP_CHECK_VECTOR(rateMatL);
	WRAP_CHECK_VECTOR(rateFreqL);
	WRAP_CHECK_VECTOR(rateDayCountL);
	WRAP_CHECK_VECTOR(resetDatesL);
	WRAP_CHECK_VECTOR(fwdRateL);
	WRAP_CHECK_VECTOR(cvxAdjL);
	WRAP_CHECK_VECTOR(mtmAdjL);


	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDateL, zcRateL) != SUCCESS)
		goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,\
		DRL_TCURVE_FMT_PERCENT));
#endif

	/* get vol curve */
	if (DrlTCurveWrapRead(&volCurve, refDateL, volDatesL, volRatesL)
		!= SUCCESS) goto done;
	volCurve->fBasis = 0;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(volCurve, vnfmFpLog,
		DRL_TCURVE_FMT_PERCENT));
#endif

	/* convert other arguments */
	if (GtoYearsToDateInterval(rateMatL[1], &rateMat) != SUCCESS)
		goto done;
	rateFreq = rateFreqL[1];
	if (GtoStringToDayCountConv(&rateDayCountL[WRAP_STR_IDX(1)],
		&rateDayCount) != SUCCESS) goto done;

	/* get number of nonzero dates */
	if ((nDates = DrlLilVectActiveSize(DRL_TDATE_L, resetDatesL,
		"resetDatesL")) <= 0)
			goto done;

	ASSERT_OR_DONE(ARGSIZE(fwdRateL) >= nDates);
	ASSERT_OR_DONE(ARGSIZE(cvxAdjL) >= nDates);
	ASSERT_OR_DONE(ARGSIZE(mtmAdjL) >= nDates);

	/*
	 * create a new TF data structure
	 */
	if ((theModel = VnfmNew2FactSimple(
		zcCurve->fBaseDate,
		(int) ARGSIZE(nfDatesL), &nfDatesL[1],
		(TCurve*)NULL,
		zcCurve,
		VNFM_LOGNORMAL_DIST,
		nfParamsL[1], nfParamsL[2],
		1e0, nfParamsL[3],
		2e-1, 2e-2,
		nfParamsL[4])) == NULL) goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: parameters:\n", routine);\
	VnfmFpWrite(theModel, vnfmFpLog);\
	DrlFPrintf(vnfmFpLog, "%s: volCurve:\n", routine);\
	DrlTCurveFpWrite(volCurve, vnfmFpLog, DRL_TCURVE_FMT_PERCENT));
#endif

	/*
	 * perform the calibration
	 */
	volFreq = 0;	/*  simple rate */
	volMat = volMatL[1];
	if (VnfmCalib1VVolCurve(
		theModel,
		0, /* CMS calibration */
		volMat,
		volFreq,
		volCurve,
		0e0, 1e2,
		NULL,	/* report all errors */
		&resValue) != 0) goto done;

	/*
	 * Compute adjustment
	 */
	if (VnfmFwdFutAdjustment(
		theModel,
		zcCurve,
		rateMat,
		rateFreq,
		rateDayCount,
		nDates,
		&resetDatesL[1],
		&fwdRateL[1],
		&cvxAdjL[1],
		&mtmAdjL[1]) != SUCCESS) goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "Rate Maturity: %s\n",\
		GtoFormatDateInterval(&rateMat));\
	DrlFPrintf(vnfmFpLog, "Dates:  %d\n", nDates);\
	for (n=0; n<=nDates-1; n++) {\
		DrlFPrintf(vnfmFpLog, "\t[%3d/%3d] %s\n",\
			n, nDates, GtoFormatDate(resetDatesL[n]));\
	});
#endif

	/* wrap out */

	/* made it through OK */
	status = SUCCESS;
done:
	GtoFreeTCurve(zcCurve);
	GtoFreeTCurve(volCurve);
	VnfmFree(theModel);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
}





/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_SWMAT_CHECK.
 *                                                             
 * <br><br>
 * Convenience routine to check that a selected array of maturities
 * in a given swaption matrix can be bootstrapped with a given
 * set of n-factor parameters.
 * <br>
 * <br>[refDateL] zero curve value date.
 * <br>[zcDateL] zero curve coupon dates.
 * <br>[zcRateL] zero curve coupon rates.
 * <br>[nfParamsL] array containg the mean reversion, weight and
 * correlation coefficients.
 * Tree formats are available: an array format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n}),$$
 * that has ${n(n+3)/ 2}$ elements,
 * or a ``short'' array format where the 1st weight,
 * assumed to be 1.0, is omitted 
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_2,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n})$$
 * that has ${n(n+3)/ 2}-1$ elements,
 * or a matrix format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,1},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n,n})$$
 * that has $n(n+2)$ elements,
 * {\bf Remark: because of possible ambiguity in the format
 * (e.g. 8 elements can correspond to $n=3$ of the short array format
 * or $n=2$ of the matrix form) the format is checked in 
 * the previous order).}
 * <br>[floatScalarsL] numerical scalars:\\
 *	(1) back bone $q$ (0 for log-normal, 0.5 for normal),\\
 *	(2) minimum admissible spot volatility ($&gt;= 0$), \\
 *	(3) minimum admissible calibration maturity.
 * <br>[swTypeL] array of 2 elements:\\
 * 	(1) swaption matrix type (0 for vertical, 1 for diagonal)
 *		for the output <i> SwVolL</i>.\\
 * 	(2) swaption matrix volatility frequency (0,1,2,4,12)
 *		for the output <i> SwVolL</i>.
 * <br>[swMatL] array of maturity intervals (in yrs)
 *		for the output <i> SwVolL</i>.
 * <br>[swExpL] array of expration intervals (in yrs)
 *		for the output <i> SwVolL</i>.
 * <br>[swVolL] input swaption matrix volatilities
 * specified by arguments <i> swTypeL</i>, <i> swMatL</i> and <i> swExpL</i>.
 * %
 * <br>[volCheckFlagL] vol check flags:\\
 *	(1) check spot vol ratio flag  (1=Yes, 0=No) \\
 *	(2) vol adjustment flag (1=Yes, 0=No)  \\
 *	(3) output vol type (1=adjusted vol, 0=vol corretion). \\
 *	(4) vol adjustment type (1=higher, -1=lower, 0=no restrictions).
 * <br>[numScalarsL] numerical scalars:\\
 *	(1) maximum admissible spot vol ratio ($>1.0$).\\
 *	(2) minimum volatility adjustment step ($>0$. Not used with optimization). 
 * <br>[tExpToCheckL] array of maturities to check.
 * <br>[finalMatL]   array of final flags. constant (0) or final (1) \\
 * <br>
 * The output consists in two arrays <i> tExpFailed</i> and <i> tMatFailed</i>
 * of same length that gives the expiration/maturity of points
 * in swaption matrix that caused failure. The output arrays contain
 * the number 0 whwn no failure is to be reported.
 */

DLL_EXPORT(int)
VnfmCheckSwaptionCalibrationL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates	*/
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */
	
	FloatL *nfParamsL,	/*  4 'F' (I) N-F params (b1,b2,r,a) */
	FloatL *floatScalarsL,	/*  5 'F' (I) floating scalars */
				/*	  [1] vol tweak size */
				/*	  [2] distribution type (0.5=N, 0=LN) */
				/*	  [3] min spot vol  */
				/*	  [4] min vol calibration maturity  */ 
	IntL *swTypeL,		/*  6 'L' (I) mattype, freq  */
	TDateIntervalL *swMatL,	/*  7 'F' (I) array of mat intervals */
	TDateIntervalL *swExpL,	/*  8 'F' (I) array of exp intervals */
	FloatL *swVolL,		/*  9 'F' (I) swaption volatilities */

	IntL   *volCheckFlagL,	/* 10 'L' (I) Vol check flags  */
				/*	  [1] check spot vol ratio flag */ 
				/*	      1=Yes, 0=No		*/
				/*	  [2] vol adjutment flag      */
				/*	      1=Yes, 0=No        */ 
				/*	  [3] 1 = output adjusted swap vol*/ 
				/*	      0 = output vol correction */ 
				/*	  [4] 1 = adj vol >= mkt vol */ 
				/*	     -1 = adj vol <= mkt vol */ 
				/*	      0 = none   */ 
	FloatL *numScalarsL,	/* 11 'L' (I) numerical scalars  */
				/*	  [1] max spot vol raio >1.0     */ 
				/*	  [2] min vol adjust amount (Not used 
					      with optimiztion */ 

	FloatL *tExpToCheckL,	/* 12 'F' (I) array of exp to check */
	IntL   *finalMatL,	/* 13 'L' (I) array of final flags  */
				/*            TRUE = final; FALSE = CMS  */


	FloatL *tExpFailedL,	/* 14 'F' (O) output failed exps */
	FloatL *tMatFailedL,	/* 15 'F' (O) output failed mats  */
	FloatL *outputSWVolL)	/* 16 'F' (O) output modified swaption vols  */
{
static	char	routine[] = "VnfmCheckSwaptionCalibrationL";
	int	status = FAILURE;

	TCurve			*zcCurve = NULL;
	TSwaptionMatrix2D	*swoptMat = NULL;

	double			spotVolMin = 0e0;

	int			volAdjFlag = FALSE;  /* default no vol adj */
	int			outputVolFlag = TRUE;/* adjusted vol matrix */
	int			spotVolRatioFlag = TRUE; /* chech vol ratio */
	int			volBoundFlag = 1;  /*  adj vol bound */

	double			volAdjAmount = 0.0025;  /* default 1/4 vega */
	double			spotVolRatio = 5.0;  /* default 5 times */
	double			tMatMin = 0.25;      /* default 3M  */

	double                  *tExpFailedTotal = NULL;
        double                  *tMatFailedTotal = NULL;
        int                     nFailedTotal, idx;

	VNFM_LOGOPEN

	/* Log wrapper inputs  */
	if (GtoLoggingGet() > 0) {
	    DrlLilVectLoggingFile("wrapper.log", "w", "TF_SWMAT_CHECK",
		DRL_TDATE_L,  refDateL,		"ZC_REFDATE",
		DRL_TDATE_L,  zcDateL,		"ZC_DATES",
		DRL_FLOAT_L,  zcRateL,		"ZC_RATES",

		DRL_FLOAT_L,  nfParamsL,	"NF_PARAMS",
		DRL_FLOAT_L,  floatScalarsL,	"FLOAT_SCALARS",

                DRL_LONG_L,   swTypeL,		"SW_TYPE",
                DRL_FLOAT_L,  swMatL,		"SW_MAT",
                DRL_FLOAT_L,  swExpL,		"SW_EXP",
                DRL_FLOAT_L,  swVolL,		"SW_VOLS",

                DRL_LONG_L,   volCheckFlagL,	"VOL_CHECK_FLAGS",
                DRL_FLOAT_L,  numScalarsL,	"NUM_SCALARS",

                DRL_FLOAT_L,  tExpToCheckL,	"EXPS_CHECK",
                DRL_LONG_L,   finalMatL,	"FINAL_FLAGS",

                DRL_FLOAT_L,  tExpFailedL,	"EXPS_FAILED_OUT",
                DRL_FLOAT_L,  tMatFailedL,	"MATS_FAILED_OUT",
                DRL_FLOAT_L,  outputSWVolL,	"SW_VOL_OUT",
                0);
	}


	/* Check arg length */
	/* */
	WRAP_CHECK_SCALAR(refDateL);
	WRAP_CHECK_VECTOR(zcDateL);
	WRAP_CHECK_VECTOR(zcRateL);
	/*WRAP_CHECK_VECTOR_LEN(nfParamsL, 5); */
	WRAP_CHECK_VECTOR(nfParamsL);
	WRAP_CHECK_VECTOR(swTypeL);
	WRAP_CHECK_VECTOR(swMatL);
	WRAP_CHECK_VECTOR(swExpL);
	WRAP_CHECK_VECTOR(swVolL);

	if (volCheckFlagL)
	{
	   if((long)volCheckFlagL[0] == 3)
	   {	
		spotVolRatioFlag = volCheckFlagL[1];
		volAdjFlag = volCheckFlagL[2];
		outputVolFlag = volCheckFlagL[3];
		volBoundFlag = 1;  /* default adj vol >= mkt vol */
	   }
	   else if ((long)volCheckFlagL[0] == 4)
	   {	
		spotVolRatioFlag = volCheckFlagL[1];
		volAdjFlag = volCheckFlagL[2];
		outputVolFlag = volCheckFlagL[3];
		volBoundFlag = volCheckFlagL[4]; 
	   }
	   else {
		GtoErrMsg("%s: invalid number of arguments (%ld) "
			  "in volCheckFlag.  Either 3 or 4.\n",
			 routine, volCheckFlagL[0]);
		goto done;	
	   }
	
	}
				
	if (numScalarsL)
	{
		spotVolRatio = numScalarsL[1];

	    if (ARGSIZE(numScalarsL)>=2) {
		volAdjAmount = numScalarsL[2]; 
		if(volAdjAmount<= 0.0e0){
			GtoErrMsg("%s: Minimum volatility adjustment (%lf),"
				  " in numScalars must be positive.\n",
				  routine, numScalarsL[2]);
			goto done;
		}
	    }
	}

	if (tExpToCheckL)
	{
		WRAP_CHECK_VECTOR(tExpToCheckL);
		ASSERT_OR_DONE(ARGSIZE(tExpToCheckL) == ARGSIZE(finalMatL));
	}	

	WRAP_CHECK_VECTOR(tExpFailedL);
	WRAP_CHECK_VECTOR(tMatFailedL);

	/* Clear the memory */
	for (idx=0; idx<=(int) ARGSIZE(tExpFailedL)-1; idx++)
		tExpFailedL[idx+1] = 0.0;

	for (idx=0; idx<=(int) ARGSIZE(tMatFailedL)-1; idx++)
		tMatFailedL[idx+1] = 0.0;


	/* Get arguments  */

	/* get zero curve  */
	IF_FAILED_DONE( DrlTCurveWrapRead(
		&zcCurve,
		refDateL,
		zcDateL,
		zcRateL));
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "%s: input zero curve:\n", routine); \
		DrlTCurveFpWrite(zcCurve, vnfmFpLog, DRL_TCURVE_FMT_PERCENT));
#endif

	/* get swaption matrix */
	IF_FAILED_DONE( DrlTSwaptionMatrix2DWrapRead(
		&swoptMat,
		swTypeL,
		swMatL,
		swExpL,
		swVolL));

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "%s: input matrix:\n", routine); \
		DrlTSwaptionMatrix2DFpWrite(swoptMat, \
			vnfmFpLog, TSWAPTION_MATRIX_FMT_STD));
#endif

	switch (ARGSIZE(floatScalarsL)) {
	case 4:
		tMatMin = floatScalarsL[4];	
		spotVolMin = floatScalarsL[3];
	case 3:
		spotVolMin = floatScalarsL[3];
	case 2:
		break;
	default:
		GtoErrMsg("%s: invalid floatScalarsL range length (%d).\n",
			routine, ARGSIZE(floatScalarsL));
		goto done;
	}


	switch ((int) floatScalarsL[2]) {
	case 0:
	case 1:
		break;
	default:
		GtoErrMsg("%s: invalid dist type %lf (0=N,1=LN).\n",
			routine, floatScalarsL[2]);
		goto done;
	}

	
	/* routine call */
	if (VnfmCheckAdjustSwaptionCalibration(
			zcCurve,
			swoptMat,

			nfParamsL,	   /* mr, alphas, etc. */
			floatScalarsL[2],  /* process power: 0=N, 1=LN, etc. */
			floatScalarsL[1],  /* vol tweak size */
			spotVolMin,	   /* minimum spot volatility */
			tMatMin,	   /* minimum calibration maturity */

			volAdjFlag,	   /* adjust vol flag */
			volAdjAmount,      
			volBoundFlag,
			spotVolRatioFlag,
			spotVolRatio,

			(int) tExpToCheckL[0],
			tExpToCheckL+1,
			finalMatL+1,	   /* final flag */

			&tExpFailedTotal,
			&tMatFailedTotal,
			&nFailedTotal) != SUCCESS)
				goto done;
	

	/* Output failed points */
        nFailedTotal = MIN(nFailedTotal, (int) ARGSIZE(tExpFailedL));
        nFailedTotal = MIN(nFailedTotal, (int) ARGSIZE(tMatFailedL));
        for (idx=0; idx<=nFailedTotal-1; idx++) {
                tExpFailedL[idx+1] = tExpFailedTotal[idx];
                tMatFailedL[idx+1] = tMatFailedTotal[idx];
        }

	/* wrap swaption vol output */
	if(outputSWVolL != NULL)
		if (DrlTSwaptionMatrix2DWrapWrite(swoptMat, NULL, NULL, NULL, 
					  outputSWVolL) 
					  != SUCCESS) goto done;
	if(!outputVolFlag){
		for (idx = 1; idx <= outputSWVolL[0]; idx++)
			outputSWVolL[idx] -= swVolL[idx];
	}

	if (GtoLoggingGet() > 0) {
	    DrlLilVectLoggingFile("wrapper.log", "a", "TF_SWMAT_CHECK",
                DRL_FLOAT_L,  tExpFailedL, "TEXP_FAILED",
                DRL_FLOAT_L,  tMatFailedL, "TMAT_FAILED",
		0);
	}


	/* made it */
	status = SUCCESS;
done:
	if (tExpFailedTotal) FREE(tExpFailedTotal);
	if (tMatFailedTotal) FREE(tMatFailedTotal);
	DrlTSwaptionMatrix2DFree(swoptMat);
	GtoFreeTCurve(zcCurve);


	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}
	VNFM_LOGCLOSE
	return(status);
}




/*---------------------------------------------------------------
 * For backwards compatibility
 */

DLL_EXPORT(int)
VnfmCheckSwaptionCalibrationOldL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates	*/
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */
	
	FloatL *nfParamsL,	/*  4 'F' (I) N-F params (b1,b2,r,a) */
	FloatL *floatScalarsL,	/*  5 'F' (I) numerical scalars */
				/*	  [1] distribution type (0=N, 1=LN) */
				/*	  [2] min spot vol  */
	IntL *swTypeL,		/*  6 'L' (I) mattype, freq  */
	TDateIntervalL* swMatL,	/*  7 'F' (I) array of mat intervals */
	TDateIntervalL* swExpL,	/*  8 'F' (I) array of exp intervals */
	FloatL *swVolL,		/*  9 'F' (I) swaption volatilities */

	FloatL *tExpFailedL,	/* 10 'F' (O) output swaption vols */
	FloatL *tMatFailedL)	/* 11 'F' (O) output swaption vols  */
{
	return VnfmCheckSwaptionCalibrationL(
		refDateL,
		zcDateL,
		zcRateL,
		nfParamsL,
		floatScalarsL,
		swTypeL,
		swMatL,
		swExpL,
		swVolL,
		NULL,
		NULL,
		NULL,
		NULL,
		tExpFailedL,
		tMatFailedL,
		NULL);
}




/*f-------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_AVGQB_VOL.
 *                                                             
 * <br><br>
 * For a set of parameters "that",
 * computes and returns the average volatility between
 * two times "S1" and "S2" of a forward rate
 * of forward maturity "tMat", time to reset "tReset" and
 * frequency "freq" (1,2,4,12).\\
 * <b> Warning.</b> The routine assumes that the VnfmData
 * structure has been updated.
 */
 
DLL_EXPORT(int)
VnfmAvgQBVolL(
        VnfmData *that,	  	  /* (I) VnfmData structure */
        double *S1L,              /* (I) start */
        double *S2L,              /* (I) expiration */
        double *tResetL,          /* (I) rate reset */
        double *tMatL,            /* (I) rate maturity */
        int    *freqL,		  /* (I) rate frequency */
        double *retValL)          /* (O) volatility */
{
static  char    routine[] = "VnfmAvgQBVolL";
	int	status = FAILURE;

	double	vol;

	if (VnfmCheckValid(that) != SUCCESS)
		goto done;

	/* Check arg length */
	WRAP_CHECK_SCALAR(S1L);
	WRAP_CHECK_SCALAR(S2L);
	WRAP_CHECK_SCALAR(tResetL);
	WRAP_CHECK_SCALAR(tMatL);
	WRAP_CHECK_SCALAR(freqL);
	WRAP_CHECK_SCALAR(retValL);

	/* routine call */
        if (VnfmAvgQBVol(that,
			 S1L[1],
			 S2L[1],
		      	 tResetL[1], 
		      	 tMatL[1],  
		      	 freqL[1], 
			 LOGVOL,
		      	 &vol) != SUCCESS)
		goto done;

	retValL[1] = vol;

	status = SUCCESS;

done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}

	return(status);
}




/*f-------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_AVGQB_CORR.
 *                                                             
 * <br><br>
 * For a set of parameters "that",
 * computes and returns the average correlation between
 * time $0$ and "S" of two rates
 * that have time to reset "tReset1" and "tReset2",
 * forward maturitites "tMat1" and "tMat2",
 * and frequencies "freq1" and "freq2" (0 for MM rate, 1,2,4,12).\\
 * <b> Warning.</b> The routine assumes that the VnfmData
 * structure has been updated.
 */
 
DLL_EXPORT(int)
VnfmAvgQBCorrL(
        VnfmData *that,	  	  /* (I) VnfmData structure */
        double *SL,               /* (I) Expiration */
        double *tReset1L,         /* (I) rate reset 1 */
        double *tMat1L,           /* (I) rate maturity 1 */
        int    *freq1L,		  /* (I) rate frequency 1 */
        double *tReset2L,         /* (I) rate reset 2 */
        double *tMat2L,           /* (I) rate maturity 2 */
        int    *freq2L,		  /* (I) rate frequency 2 */
        double *retValL)          /* (O) correlation */
{
static  char    routine[] = "VnfmAvgQBCorrL";
	int	status = FAILURE;

	double	corr;

	if (VnfmCheckValid(that) != SUCCESS)
		goto done;

	/* Check arg length */
	WRAP_CHECK_SCALAR(SL);
	WRAP_CHECK_SCALAR(tReset1L);
	WRAP_CHECK_SCALAR(tMat1L);
	WRAP_CHECK_SCALAR(freq1L);
	WRAP_CHECK_SCALAR(tReset2L);
	WRAP_CHECK_SCALAR(tMat2L);
	WRAP_CHECK_SCALAR(freq2L);
	WRAP_CHECK_SCALAR(retValL);

	/* routine call */
        if (VnfmAvgQBCorr(that,
		       	  SL[1],
		       	  tReset1L[1], 
		       	  tMat1L[1],  
		       	  freq1L[1], 
		       	  tReset2L[1], 
		       	  tMat2L[1],  
		       	  freq2L[1], 
		       	  &corr) != SUCCESS)
		goto done;

	retValL[1] = corr;

	status = SUCCESS;

done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}

	return(status);
}



/*f-------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_AVGQB_CORR2.
 *                                                             
 * <br><br>
 * For a set of parameters "that",
 * computes and returns the average correlation between
 * time $0$ and "S" of two rates
 * that have time to reset "tReset1" and "tReset2",
 * forward maturitites "tMat1" and "tMat2",
 * and frequencies "freq1" and "freq2" (0 for MM rate, 1,2,4,12).\\
 * <b> Warning.</b> The routine assumes that the VnfmData
 * structure has been updated.
 */
 
DLL_EXPORT(int)
VnfmAvgQBCorr2L(
        VnfmData *that,	  	  /* (I) VnfmData structure */
	double *ObsL,             /* (I) observation */
        double *SL,               /* (I) Expiration */
        double *tReset1L,         /* (I) rate reset 1 */
        double *tMat1L,           /* (I) rate maturity 1 */
        int    *freq1L,		  /* (I) rate frequency 1 */
        double *tReset2L,         /* (I) rate reset 2 */
        double *tMat2L,           /* (I) rate maturity 2 */
        int    *freq2L,		  /* (I) rate frequency 2 */
        double *retValL)          /* (O) correlation */
{
static  char    routine[] = "VnfmAvgQBCorrL";
	int	status = FAILURE;

	double	corr;

	if (VnfmCheckValid(that) != SUCCESS)
		goto done;

	/* Check arg length */
	WRAP_CHECK_SCALAR(ObsL);
	WRAP_CHECK_SCALAR(SL);
	WRAP_CHECK_SCALAR(tReset1L);
	WRAP_CHECK_SCALAR(tMat1L);
	WRAP_CHECK_SCALAR(freq1L);
	WRAP_CHECK_SCALAR(tReset2L);
	WRAP_CHECK_SCALAR(tMat2L);
	WRAP_CHECK_SCALAR(freq2L);
	WRAP_CHECK_SCALAR(retValL);

	/* routine call */
        if (VnfmAvgQBCorr2(that,
		       	   ObsL[1],
		       	   SL[1],
		       	   tReset1L[1], 
		       	   tMat1L[1],  
		       	   freq1L[1], 
		       	   tReset2L[1], 
		       	   tMat2L[1],  
		       	   freq2L[1], 
		       	   &corr) != SUCCESS)
		goto done;

	retValL[1] = corr;

	status = SUCCESS;

done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}

	return(status);
}



/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CALNF1SPVOL_AVGON.
 *                                                             
 * <br><br>
 * Wrapper for routine <i> VnfmCalib1VNFAverageONL</i>
 * that simultaneously performs spot volatility calibration 
 * of arbitrary volatity points
 * and generate implied volatility for a simple average rate of
 * over night rates.
 * <br>
 * <br>[refDateL] zero curve value date.
 * <br>[zcDateL] zero curve coupon dates.
 * <br>[zcRateL] zero curve coupon rates.
 * <br>[floatScalarsL] numerical  scalars:\\
 *	(1) back bone $q$ (0 for log-normal, 0.5 for normal),\\
 *	(2) minimum volatility rate maturity($&gt;= 0$).
 * <br>[nfParamsL] LIL array containg the mean reversion, weight and
 * correlation coefficients.
 * Tree formats are available: an array format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n}),$$
 * that has ${n(n+3)/ 2}$ elements,
 * or a ``short'' array format where the 1st weight,
 * assumed to be 1.0, is omitted 
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_2,\dots,\alpha_{n},
 *    \rho_{1,2},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n-1,n})$$
 * that has ${n(n+3)/ 2}-1$ elements,
 * or a matrix format
 * $$(\beta_1,\dots,\beta_{n},
 *    \alpha_1,\dots,\alpha_{n},
 *    \rho_{1,1},\dots,\rho_{1,n},\rho_{2,1},\dots,\rho_{n,n})$$
 * that has $n(n+2)$ elements,
 * {\bf Remark: because of posort array format
 * or $n=2$ of the matrix form) the format is checked in 
 * the previous order).}
 * <br>[rateDatesL] array of dates used for the spot volatility
 * bootstrapping (the spot volatility is assumed to be constant between
 * two consecutive dates).\\
 * {\bf WARNING: Unlike the ALIB routines, the convention is that
 * the spot volatility $i$ applies between date $i$ and date $i+1$.
 * The first date in the array MUST be the today date}.
 * <br>[rateMatL] array of volatility rate maturity (in yrs).
 * The element of index $i$ in the array corresponds
 * to an option expiring at date $i$ of array <i> rateDatesL</i>.
 * This array must have same length as the array <i> rateDatesL</i>.
 * {\bf
 * <br>
 *  <br> The volatility 1 (corresponding to the date 1, i.e. today)
 *     is never calibrated.
 * <br> The calibration stops at the first negative maturity encountered
 *     according to the following rules:
 *     <br>
 *     <br> If the maturity 1 (corresponding to the date 1, i.e. today)
 *         is negative, an error is returned.
 *     <br> If the maturity 1 is positive,
 *         but maturity 2 is negative, then volatility 2 is
 *         nevertheless calibrated assuming a maturity 0.25.
 *         No other volatility point is calibrated.
 *     <br> Otherwise, the volatility points are calibrated until a negative
 *         maturity is encountered.
 * <br>
 * <br> Any maturity smaller than floatScalarsL[4], but positive,
 *     will be rounded to floatScalarsL[4] (and the corresponding volatility 
 *     will be calibrated as such).
sible ambiguity in the format
 * (e.g. 8 elements can correspond to $n=3$ of the short array format
 * or $n=2$ of the matrix form) the format is checked in 
 * the previous order).}
 * <br>[rateDatesL] array of dates used for the spot volatility
 * bootstrapping (the spot volatility is assumed to be constant between
 * two consecutive dates).\\
 * {\bf WARNING: Unlike the ALIB routines, the convention is that
 * the spot volatility $i$ applies between date $i$ and date $i+1$.
 * The first date in the array MUST be the today date}.
 * <br>[rateMatL] array of volatility rate maturity (in yrs).
 * The element of index $i$ in the array corresponds
 * to an option expiring at date $i$ of array <i> rateDatesL</i>.
 * This array must have same length as the array <i> rateDatesL</i>.
 * {\bf
 * <br>
 *  <br> The volatility 1 (corresponding to the date 1, i.e. today)
 *     is never calibrated.
 * <br> The calibration stops at the first negative maturity encountered
 *     according to the following rules:
 *     <br>
 *     <br> If the maturity 1 (corresponding to the date 1, i.e. today)
 *         is negative, an error is returned.
 *     <br> If the maturity 1 is positive,
 *         but maturity 2 is negative, then volatility 2 is
 *         nevertheless calibrated assuming a maturity 0.25.
 *         No other volatility point is calibrated.
 *     <br> Otherwise, the volatility points are calibrated until a negative
 *         maturity is encountered.
 * <br>
 * <br> Any maturity smaller than floatScalarsL[4], but positive,
 *     will be rounded to floatScalarsL[4] (and the corresponding volatility 
 *     will be calibrated as such).
 * <br>
 * }
 * <br>[rateFreqL] array of volatility frequency (0,1,2,4,12).
 * {\it The element of index $i$ in the array corresponds
 * to an option expiring at date $i$ of array <i> rateDatesL</i>}.
 * <br>[rateVolL] array of volatilities.
 * {\it The element of index $i$ in the array corresponds
 * to an option expiring at date $i$ of array <i> rateDatesL</i>}.
 *
 * <br>[startDateL] starting date of averaging period.
 * <br>[endDateL] end date of averaging period.
 * <br>[resetDateL] reset date of option.
 *
 * <br>[avgONVolL] output volatility of average ON rate.
 * <br>
 */

DLL_EXPORT(int)
VnfmCalib1VNFAverageONL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */
	
	FloatL *floatScalarsL,	/*  4 'F' (I) num scalars:
				 *        [1] back bone q
				 *        [2] minimum volatility rate mat */
	FloatL *nfParamsL,	/*  5 'F' (I) N-fact params */

	TDateL *rateDatesL,	/*  6 'D' (I) array of dates */
	FloatL *rateMatL,	/*  7 'F' (I) rate mat [0..nDates-1] */
	IntL   *rateFreqL,	/*  8 'L' (I) rate freq [0..nDates-1] */
	FloatL *rateVolL,	/*  9 'F' (I) rate vol [0..nDates-1] */

	TDateL *startDateL,	/* 10 'D' (I) start date of averaging period */
	TDateL *endDateL,	/* 11 'D' (I) end date of averaging period */
	TDateL *resetDateL,	/* 12 'D' (I) reset date of option */

	FloatL *avgONVolL)	/* 13 'F' (O) vol of average ON rate */
{
static	char		routine[] = "VnfmCalib1VNFAverageONL";
	int		status = FAILURE;


	VnfmData	*vnfmData = NULL;
	TCurve		*zcCurve = NULL;
	int		nDates, i;

	int		idxStart = 0,
			idxEnd;
	double		*rateMat = &rateMatL[1],
			*rateVol = &rateVolL[1]; 
	int		*rateFreq = NULL;


	VNFM_LOGOPEN

	/* check arg len */
	WRAP_CHECK_SCALAR(refDateL);
	WRAP_CHECK_VECTOR(zcDateL);
	WRAP_CHECK_VECTOR(zcRateL);
	WRAP_CHECK_VECTOR_LEN(floatScalarsL, 2);
	WRAP_CHECK_VECTOR(nfParamsL);

	WRAP_CHECK_VECTOR(rateDatesL);
	WRAP_CHECK_VECTOR(rateMatL);
	WRAP_CHECK_VECTOR(rateFreqL);
	WRAP_CHECK_VECTOR(rateVolL);

	WRAP_CHECK_VECTOR_LEN(startDateL, 1);
	WRAP_CHECK_VECTOR_LEN(endDateL, 1);
	WRAP_CHECK_VECTOR_LEN(resetDateL, 1);

	WRAP_CHECK_VECTOR(avgONVolL);


	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDateL, zcRateL) != SUCCESS)
		goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve,\
		vnfmFpLog, DRL_TCURVE_FMT_PERCENT));
#endif


	nDates = ARGSIZE(rateDatesL);
	ASSERT_OR_DONE(ARGSIZE(rateMatL) == ARGSIZE(rateDatesL));
	ASSERT_OR_DONE(ARGSIZE(rateFreqL) == ARGSIZE(rateDatesL));
	ASSERT_OR_DONE(ARGSIZE(rateVolL) == ARGSIZE(rateDatesL));
	if ((rateFreq = NEW_ARRAY(int, ARGSIZE(rateFreqL))) == NULL)
		goto done;
	for (i=0; i<=ARGSIZE(rateFreqL)-1; i++)
		rateFreq[i] = (int) rateFreqL[i+1];



	/**	
	 * Check dates, minimum 1 day.
	 */
	ASSERT_OR_DONE(endDateL[1] > startDateL[1]);
	ASSERT_OR_DONE(endDateL[1] >= resetDateL[1]);

	/** Only normal vnfm allowed */
	if(!IS_ALMOST_ZERO(floatScalarsL[1]-0.5))
	{
		GtoErrMsg("%s: only normal backbone allowed(bbq=%lf !=0.5).\n",	
			 routine, floatScalarsL[1]);
		goto done;
	}

	/**
	 ** Create Vnfm structure
	 **/
	if (VnfmWrapReadSimple(
		&vnfmData,
		floatScalarsL,
		nfParamsL,
		rateDatesL,
		zcCurve) != SUCCESS)
			goto done;

	/** 
	 ** Check the cut-off rate maturity 
	 **/
	if (floatScalarsL[4] < 0e0) {
	    GtoErrMsg("%s: the cut-off vol rate maturity is negative(%lf).\n",
                routine, floatScalarsL[4]);
            goto done;
        }

	/**
	 ** Adjust array maturity 
	 **/
	/* Set last calibrated point to be the 1st negative maturity
	 * (or the 1st less than 1w that is not the 1st point).
	 */
	if (rateMat[0] < 0e0) {
	    GtoErrMsg("%s: first maturity found in calibrated "
		"rate array is negative (%lf).\n", routine, rateMat[0]);
	    goto done;
	}

	/* treat the special case of maturity between today (date idx 0)
	 * and next vol point (date idx 1)
	 */
	if (rateMat[1] < 0e0) {
		idxEnd = 1;
	} else {
		for (idxEnd=1; idxEnd<=nDates-1; idxEnd++) {
		    if (rateMat[idxEnd] < 0e0) 
			break;
		}
		idxEnd--;
	}

	/* Cutoff for small maturities */
	for (i=0; i<=idxEnd; i++) {
		rateMat[i] = MAX(rateMat[i], floatScalarsL[4]);
	}


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: idxEnd=%d "\
		"(last idx calibrated).\n", routine, idxEnd))
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "\tInput vols:\n"\
		"\tidx         dates     rateMat rateFreq  rateVol\n"))
	GTO_IF_LOGGING(\
	for (i=0; i<=nDates-1; i++) {\
	    DrlFPrintf(vnfmFpLog, "\t[%2d/%2d]  %10s  %8.4f  %d  %8.4f%%\n",\
		i, nDates,
		GtoFormatDate(rateDatesL[i+1]),
		rateMat[i], rateFreq[i], rateVol[i]*1e2);\
	})
#endif

	if (idxEnd < 0) {
	    GtoErrMsg("%s: no positive maturity found in calibrated "
		"rate array.\n", routine);
	    goto done;
	}



	/* perform the calibration */
	if (VnfmVolCalib1VArbitrary(
		vnfmData,
		idxStart,
		idxEnd,
		rateMat,
		rateFreq,
		rateVol,
		LOGVOL,
		TRUE, NULL) != SUCCESS)
			goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: calibrated:\n", routine))
	GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog))
#endif



	/**
	 ** Generate vol of average rate
	 **/
	if (VnfmAvgVolONRate(
			vnfmData,
			startDateL[1],
			endDateL[1],
			resetDateL[1],
			&avgONVolL[1]) != SUCCESS)
		goto done;


	/* made it through OK */
	status = SUCCESS;
done:


	FREE(rateFreq);
	VnfmFree(vnfmData);
	GtoFreeTCurve(zcCurve);

	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}
	VNFM_LOGCLOSE
	return(status);

}



/* 
 * Convert vol type from string to KVolType
 */
int
ConvertVolType(
	char *volTypeStr,
	KVolType *volType)
{ 
	int	status = SUCCESS;

	if (toupper(volTypeStr[0]) == 'L')
		*volType = LOGVOL;
	else if((char)toupper(volTypeStr[0]) == 'N')
		*volType = NORMVOL;
	else
	{
		GtoErrMsg("invalid input vol type (%s).\n",
			  volTypeStr);
		status = FAILURE;
	}

	return status;
}



/*f-------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_AVGQB_SPRD_VOL.
 *                                                             
 * <br><br>
 * For a set of parameters "that",
 * computes and returns the spread volatility between
 * time $0$ and "S" of two rates
 * that have time to reset "tReset1" and "tReset2",
 * forward maturitites "tMat1" and "tMat2",
 * and frequencies "freq1" and "freq2" (0 for MM rate, 1,2,4,12).\\
 * <b> Warning.</b> The routine assumes that the VnfmData
 * structure has been updated.
 */
 
DLL_EXPORT(int)
VnfmAvgQBSprdVolL(
        VnfmData *that,	  	  /* (I) VnfmData structure */
        double *SL,               /* (I) Expiration */
        double *tReset1L,         /* (I) rate reset 1 */
        double *tMat1L,           /* (I) rate maturity 1 */
        int    *freq1L,		  /* (I) rate frequency 1 */
        double *tReset2L,         /* (I) rate reset 2 */
        double *tMat2L,           /* (I) rate maturity 2 */
        int    *freq2L,		  /* (I) rate frequency 2 */
        double *retValL)          /* (O) spread vol */
{
static  char    routine[] = "VnfmAvgQBSpdVolL";
	int	status = FAILURE;

	double	corr;

	if (VnfmCheckValid(that) != SUCCESS)
		goto done;

	/* Check arg length */
	WRAP_CHECK_SCALAR(SL);
	WRAP_CHECK_SCALAR(tReset1L);
	WRAP_CHECK_SCALAR(tMat1L);
	WRAP_CHECK_SCALAR(freq1L);
	WRAP_CHECK_SCALAR(tReset2L);
	WRAP_CHECK_SCALAR(tMat2L);
	WRAP_CHECK_SCALAR(freq2L);
	WRAP_CHECK_SCALAR(retValL);

	/* routine call */
        if (VnfmAvgQBSprdVol(that,
		       	     SL[1],
		       	     tReset1L[1], 
		       	     tMat1L[1],  
		       	     freq1L[1], 
		       	     tReset2L[1], 
		       	     tMat2L[1],  
		       	     freq2L[1], 
		       	     &corr) != SUCCESS)
	         goto done;

	retValL[1] = corr;

	status = SUCCESS;

done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}

	return(status);
}


/*f-------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_AVGQB_ORTH_FACT.
 *                                                             
 * <br><br>
 * For a set of parameters "that",
 * computes and returns the factor shapes over an observation
 * period determined by start and end times of a number of rates
 * defined by their resets and maturities in years, and 
 * coupon frequencies (0 for MM rate, 1,2,4,12).\\
 * <b> Warning.</b> The routine assumes that the VnfmData
 * structure has been updated.
 */
 
DLL_EXPORT(int)
VnfmAvgQBOrthFactorsL(
        VnfmData *that,	  	  /* (I) VnfmData structure */
	double   *tStartL,        /* (I) Start of observation period */
	double   *tEndL,          /* (I) End of observation period */
	double   *rateExpYrsL,    /* (I) Rate expirations in years */
	long     *rateFreqL,      /* (I) Rate frequencies in years */
	double   *rateMatYrsL,    /* (I) Rate maturities in years */
	TMatrix2D **eigVectL,     /* (O) Array of eigenvectors */
	TMatrix2D **eigValL)      /* (O) Array of eigenvalues */
{
static  char    routine[] = "VnfmAvgQBOrthFactorsL";
	int	status = FAILURE;

	int    n         = (int) rateExpYrsL[0];
	double **eigVect = NULL;

	if (VnfmCheckValid(that) != SUCCESS)
		goto done;

	/* Check arg length */
	WRAP_CHECK_SCALAR(tStartL);
	WRAP_CHECK_SCALAR(tEndL);
	WRAP_CHECK_VECTOR_LEN (rateExpYrsL, n);
	WRAP_CHECK_VECTOR_LEN (rateFreqL,   n);
	WRAP_CHECK_VECTOR_LEN (rateMatYrsL, n);

	/* Create matrices */
	if ((*eigVectL = GtoMatrixNewEmpty(n,n)) == NULL)
                goto done;
	if ((*eigValL  = GtoMatrixNewEmpty(1,n)) == NULL)
                goto done;

	/* routine call */
        if (VnfmAvgQBOrthFactors(that,
				 tStartL[1],
				 tEndL[1],
				 n,
				 rateExpYrsL+1,
				 (int *) rateFreqL+1,
				 rateMatYrsL+1,
				 (*eigVectL)->data,
				 ((*eigValL)->data)[0]) != SUCCESS)
		goto done;

	status = SUCCESS;

done:

	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}

	return(status);
}
