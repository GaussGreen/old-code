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
#include "date_sup.h"
#include "convert.h"

#include "drlsmat.h"		/* TSwaptionMatrix2D */
#include "drlio.h"			/* DrlFPrintf */

#define	_vnfm_SOURCE
#include "vnfmcali.h"

#undef __DEBUG__

#if defined(_WINDLL) 
# undef __DEBUG__
#endif

/*f-------------------------------------------------------------
 * Bootstrap 1D spot volatility (2F convenience).
 *                                                             
 * <br><br>
 * Convenience routine to perform a single spot volatility
 * calibration of the 2F model on arbitrary rates and
 * generate a model base volatility and swaption matrix.
 * \vspace{2mm}\par\noindent<b> Argument Details:</b>
 * <br>
 * <br>[zcCurve] zero coupon term structure.
 * <br>[backboneq] backbone q (0 = lognormal, 0.5 = normal,
 * 			i.e. cevPower = 1-2q)
 * <br>[beta1] first mean-reversion coefficient.
 * <br>[beta2] second mean-reversion coefficient.
 * <br>[alpha] weight of the second factor.
 * <br>[rho] correlation between factors.
 * <br>[nDates] length of array <i> dates</i>.
 * <br>[dates] array of dates used to generate the 2F timeline 
 * and the output base volatility curve. Pass NULL to use
 * a default of 12 periods per year.
 * <b> WARNING: the first date in the array must be the today date.</b>
 * <br>[rateMat] array of maturities (in years) of rates
 * to be calibrated (corresponding to the dates in the time line).
 * This array must have same length as the array "dates".
 * {\bf
 * <br>
 * <br> The volatility 1 (corresponding to the date 1, i.e. today)
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
 * <br> Any maturity smaller than 0.25 years, but positive,
 *     will be rounded to 0.25 (and the corresponding volatility 
 *     will be calibrated as such).
 * <br>
 * }
 * <br>[rateFreq] array of rates frequencies (1,2,4,12).
 * This array must have same length as the array "dates".
 * <br>[rateVol] array of rates volatilities.
 * This array must have same length as the array "dates".
 * <br>[nSpotVols] size of output spot volatilities vector "spotVols".
 * <br>[spotVols] on exit, contains the spot vols (should be allocated
 * before entry).
 * <br>[nBvDates] number of base volatility dates for output
 * base volatility.
 * <br>[bvDates] array of base volatility dates.
 * <br>[bvMat]	output base volatility curve rate  maturity.
 * <br>[bvCurve] on exit, base volatility curve. To be freed
 * by <i> GtoFreeCurve</i>.
 * <br>[swMat] model swaption matrix to be generated. The matrix 
 * should already be initialized on entry. On exit, contains the
 * model swaption volatilities.
 * <br>
 */

DLL_EXPORT(int)
VnfmCalib1V2FGeneral(
	TCurve *zcCurve,	/* (I) zero coupon term structure */
 	double backboneq,	/* (I) backbone q */
	double beta1,		/* (I) 2F parameters */
	double beta2,		/* (I) idem  */
	double alpha,		/* (I) idem  */
	double rho,		/* (I) idem  */
	TDate refDate,		/* (I) today's date */
	int nDates,		/* (I) # dates */
	TDate *dates,		/* (I) timeline dates [0..nDates-1] */
	double *rateMat,	/* (I) arrays of rate mat [0..nDates-1] */
	int *rateFreq,		/* (I) arrays of rate freq [0..nDates-1] */
	double *rateVol,	/* (I) arrays of rate vol [0..nDates-1] */

	int nSpotVols,		/* (I) # of spot vols (can be 0) */
	double *spotVols,	/* (O) array of spot vols (can be NULL) */

	int nBvDates,		/* (I) # of base vol dates */
	TDate *bvDates,		/* (I) array of base vol dates */
	TDateInterval bvMat,	/* (I) base vol maturity */
	TCurve **bvCurve,	/* (O) model base volatility curve (or NULL) */
	TSwaptionMatrix2D *swMat)/* (O) model swaption matrix (or NULL) */
{
static	char		routine[] = "VnfmCalib1V2FGeneral";
	int		status = FAILURE;

	int		i;
	VnfmData	*that = NULL;

	int		nSwoptExp, nSwoptMat;
	TDateInterval	*swoptExp = NULL,
			*swoptMat = NULL;

	int		idxStart = 0,
			idxEnd;
	double		rateMatMin = 0.25;	/* min calibrated maturity */

	/* create a new TF data structure */
	if ((that = VnfmNew2FactSimple(refDate,
		nDates, dates,
		(TCurve*)NULL,
		zcCurve,
		backboneq,
		beta1, beta2,
		1e0, alpha,
		0e0,	/* default spot volatility */
		0e0,	/* default spot volatility */
		rho)) == NULL) goto done;


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

	/* For Kapital: cutoff for small maturities */
	for (i=0; i<=idxEnd; i++) {
		rateMat[i] = MAX(rateMat[i], rateMatMin);
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
		GtoFormatDate(dates[i]),
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
		that,
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
	GTO_IF_LOGGING(VnfmFpWrite(that, vnfmFpLog))
#endif

	/* generate base volatility curve */
	if (bvCurve != NULL) {
	    if ((nBvDates <= 0) || (bvDates == NULL)) {
		if (VnfmGenerateVolTCurve(that,
			bvMat, 0, 0,
			REFDATE,
			nDates, dates,
			bvCurve) != SUCCESS) goto done;
	    } else {
		if (VnfmGenerateVolTCurve(that,
			bvMat, 0, 0,
			REFDATE,
			nBvDates, bvDates,
			bvCurve) != SUCCESS) goto done;
	    }
	}


	/* generate swaption matrix */
	if (swMat != NULL) {
	    nSwoptExp = TSWAPTION_MATRIX2D_NEXP(swMat);
	    nSwoptMat = TSWAPTION_MATRIX2D_NMAT(swMat);
	    if ((swoptExp = NEW_ARRAY(TDateInterval, nSwoptExp)) == NULL)
		goto done;
	    if ((swoptMat = NEW_ARRAY(TDateInterval, nSwoptMat)) == NULL)
		goto done;
	    for (i=0; i<=nSwoptExp-1; i++) {
		if (GtoYearsToDateInterval(TSWAPTION_MATRIX2D_EXP(swMat, i),
			swoptExp+i) != SUCCESS) goto done;
	    }
	    for (i=0; i<=nSwoptMat-1; i++) {
		if (GtoYearsToDateInterval(TSWAPTION_MATRIX2D_MAT(swMat, i),
			swoptMat+i) != SUCCESS) goto done;
	    }

	    if (VnfmSwaptionVolMatrix(
		that,
		(swMat->diagonal ? 1 : 0),
		(int)swMat->swapPayFreq, 0,
		nSwoptExp, swoptExp,
		nSwoptMat, swoptMat,
		LOGVOL,
		swMat->table->matrix->data) != SUCCESS)
			goto done;
#ifndef	NO_LOGGING
	    GTO_IF_LOGGING(DrlTSwaptionMatrix2DFpWrite(\
		swMat, vnfmFpLog, TSWAPTION_MATRIX_FMT_STD))
#endif
	}

	/* output spot volatilities */
	if ((nSpotVols > 0) && (spotVols != NULL)) {
	    if (nSpotVols < that->fNDates) {
		GtoErrMsg("%s: output array of spot vols not large enough "
			"(got %d, timeline length is %d).\n", routine,
	    		nSpotVols, that->fNDates);
		goto done;
	    }
	    for(i=0; i<=that->fNDates-1; i++) {
		spotVols[i] = SIGMA[0][i];
	    }
	}

	/* made it through OK */
	status = SUCCESS;
done:
	FREE(swoptExp);
	FREE(swoptMat);
	VnfmFree(that);

	if (status != SUCCESS)
	{
/*
	    if ((bvCurve != NULL) && (*bvCurve != NULL))
		GtoFreeTCurve(*bvCurve);
*/
	    GtoErrMsg("%s: failed \n", routine);
	}
	return(status);
}




