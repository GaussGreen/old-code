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

#include "date_sup.h"
#include "convert.h"

#include "drlsort.h"			/* DrlDoubleArrayFloorIdx */
#include "drlio.h"			/* DrlFPrintf */

#define	_vnfm_SOURCE
#include "vnfmcali.h"

#undef  __DEBUG__
#if defined(_WINDLL) 
# undef __DEBUG__
#endif


/*f-------------------------------------------------------------
 * Bootstrap 1D volatility curve (TCurve input).
 *                                                             
 * <br><br>
 * Routine that performs a multi-factor/single spot volatility
 * calibration (i.e. spot volatilities of all factors are equal)
 * on a given volatility curve.
 * \vspace{2mm}\par\noindent<b> Argument Details:</b>
 * <br>
 * <br>[that] parameter data structure. Should be initialized
 *	on entry (time grid, MR, correlations,etc.). Only spot
 *	volatilities are modified on exit.
 * <br>[calType1] type of calibration. Pass 0 for CMS calibration,
 * 	1 for final.
 * <br>[volMat1] maturity (constant or final) in years of the
 * 	input volatilities.
 * <br>[volFreq1] frequency of the input volatility.
 * <br>[tStart] time to start the spot volatility calibration
 *	(use $0.0$ as a default value).
 * <br>[tEnd] time to end the spot volatility calibration
 *	(use $1e12$ as a default value).
 * <br>[bootstrapError] if NULL, will write all error messages to
 *	the error log if a failure occurs.
 *      If not NULL, will be set to TRUE only if a bootstrapping error
 * 	occurs (nothing reported on error log), FALSE otherwise.
 * <br>[volCurve] curve of volatilities corresponding
 *	to the desired calibrated maturity (CMS of final) as
 *	a function of expiration date. The <i> fBasis</i> argument
 *	of the <i> TCurve</i> should contain the frequency of the volatility.
 * <br>[resValue] returns a bootstrap resildual.
 * <br>
 */

DLL_EXPORT(int)
VnfmCalib1VVolCurve(
	VnfmData *that,		/* (I/O) model parameters */

	int calType1,		/* (I) 0=cms, 1=final */
	double volMat1,		/* (I) swaption maturity */
	int volFreq1,		/* (I) 0=simple, 1, 2, 4, 12 */
	TCurve *volCurve1,	/* (I) input base volatility curve */

	double tStart,		/* (I) time to start vol calib */
	double tEnd,		/* (I) time to end vol calib */

	int *bootstrapError,	/* (O) reports bootstrap failure iff NULL */
	double *resValue)	/* (O) <0 if calib OK, >0 otherwise */
{
static	char	routine[] = "VnfmCalib1VVolCurve";
	int	status = FAILURE;
	double	S,
		*tMat1 = NULL,
		*vol1 = NULL;
	int	*freq1 = NULL,
		i,
		idxStart, idxEnd;


	/* set to FALSE */
	if (bootstrapError != NULL) *bootstrapError = FALSE;

	/* check that refdate of zc and bv curve agree */
	if (REFDATE != volCurve1->fBaseDate) {
	    GtoErrMsg("%s: calib timeline ref date %s != "
		"vol curve base date %s\n",
		routine,
		GtoFormatDate(REFDATE),
		GtoFormatDate(volCurve1->fBaseDate));
	    goto done;
	}


	/* get benchmarks vols */
	tMat1 = NEW_ARRAY(double, that->fNDates);
	vol1  = NEW_ARRAY(double, that->fNDates);
	freq1 = NEW_ARRAY(int, that->fNDates);
	if ((tMat1 == NULL) || (vol1 == NULL) || (freq1 == NULL)) {
		GtoErrMsg("%s: malloc failed\n", routine);
		goto done;
	}


	if (calType1 == 0) {
	    /*
	     * CMS calibration
	     */
	    VNFMIDX(tStart, &idxStart);
	    VNFMIDX(tEnd, &idxEnd);


#ifndef	NO_LOGGING
	    GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
		"%s: calType1=%d volMat1=%lf volFreq1=%d\n"\
		"\ttStart=%12.8e tEnd=%12.8e idxStart=%d, idxEnd=%d\n",\
		routine, calType1, volMat1, volFreq1,\
		tStart, tEnd, idxStart, idxEnd))
#endif


	    /* interpolate vol */
	    for (i=idxStart; i<=idxEnd; i++) {
		tMat1[i] = volMat1;
		freq1[i] = volFreq1;
	    	GtoInterpRate(that->fDate[i], volCurve1,
			GTO_LINEAR_INTERP, vol1+i);

	    }

	    if (VnfmVolCalib1VArbitrary(
			that,
			idxStart, idxEnd,
			tMat1, freq1, vol1,
			LOGVOL,
			(bootstrapError != NULL ?  FALSE :  TRUE),
			resValue) != SUCCESS) {
				if (bootstrapError != NULL)
					*bootstrapError = TRUE;
				goto done;
	    }

	} else if (calType1 == 1) {
	    /*
	     * Final maturity calibration
	     */

	    VNFMIDX(tStart, &idxStart);
	    VNFMIDX(volMat1, &idxEnd);
	    VNFMIDX(tEnd, &i);

	    idxEnd = MIN(i, idxEnd);

	    for (i=idxStart; i<=idxEnd; i++) {
		GtoDayCountFraction(REFDATE, that->fDate[i], GTO_ACT_365F, &S);
		tMat1[i] = volMat1 - S;
		freq1[i] = volFreq1;
	    	GtoInterpRate(that->fDate[i], volCurve1,
			GTO_LINEAR_INTERP, vol1+i);
		/* check maturity long enough (1w+1d) */
		if (tMat1[i] <= 2.1970e-2) {
			idxEnd = i-1;
			break;
		}

	    }
	    if (VnfmVolCalib1VArbitrary(
			that,
			idxStart, idxEnd,
			tMat1, freq1, vol1,
			LOGVOL,
			(bootstrapError != NULL ?  FALSE :  TRUE),
			resValue) != SUCCESS) {
				if (bootstrapError != NULL)
					*bootstrapError = TRUE;
				goto done;
	    }



	} else {
	    GtoErrMsg("%s: bad calType (%d).\n", routine, calType1);
	    goto done;
	}



	/* made it through OK */
	status = SUCCESS;
done:
	/* tesing of accuracy */
	VnfmComputeCoeff(that);


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: final parameters\n", routine))
	GTO_IF_LOGGING(VnfmFpWrite(that, vnfmFpLog))
	GTO_IF_LOGGING(VnfmPrintCoeff(that, vnfmFpLog))
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "\n#\n# %s:END\n", routine))
#endif

	/* free allocated memory */
	if (tMat1 != NULL) FREE((void*) tMat1);
	if (vol1  != NULL) FREE((void*) vol1);
	if (freq1 != NULL) FREE((void*) freq1);

	if (status != SUCCESS) {
	    if ((bootstrapError == NULL) || (*bootstrapError == FALSE)) {
		GtoErrMsg("%s: failed\n", routine);
	    }
	}
	return(status);
}



/*f-------------------------------------------------------------
 * Bootstrap 1D volatility curve (TCurve input, convenience).
 *                                                             
 * <br><br>
 * Convenience routine that performs a multi-factor/single spot volatility
 * calibration (i.e. spot volatilities of all factors are equal)
 * on a given volatility curve.
 * <br>
 * <br>[that] parameter data structure. Should be initialized
 *  on entry (time grid, MR, correlations,etc.). Only spot
 *  volatilities are modified on exit.
 * <br>[idxStart] TL idx of expiry date to start spot vol calibration
 *  (use 1 as a default value).
 * <br>[idxEnd] TL idx of expiry date to start spot vol calibration
 *  (use fNDates-1 as a default value).
 * <br>[finalMaturity] type of calibration. Pass 0 for CMS calibration,
 *  1 for final.
 * <br>[volMatDate] maturity date (only for finalMaturity = TRUE)
 * <br>[volMat1] constant maturity (only for finalMaturity = FALSE)
 * <br>[volCurve] curve of volatilities corresponding
 *  to the desired calibrated maturity (CMS of final) as
 *  a function of expiration date. The <i> fBasis</i> argument
 *  of the <i> TCurve</i> should contain the frequency of the volatility.
 * <br>
 */

DLL_EXPORT(int)
VnfmCalib1VCMSFinalVolTCurve(
	VnfmData *that,		/* (I/O) Model parameters */
	int idxStart,		/* (I) TL idx of 1st expiry date calibrated */
	int idxEnd,		/* (I) TL idx of last expiry date calibrated */
	TBoolean finalMaturity,	/* (I) F=cms(use intval), T=final(use date) */
	TDate finalMatDate,	/* (I) Final mat date (used if finalMat=T)*/
	TDateInterval volMat1,	/* (I) Mat interval (used if finalMat=F)*/
	TCurve *volCurve1)	/* (I) Input base volatility curve */
{
static	char routine[] = "VnfmCalib1VCMSFinalVolTCurve";
	int	status = FAILURE;
	double	finalMatYears;
	double	cmsMatYears;
	double	*tMat1 = NULL;
	double	*vol1 = NULL;
	int	*freq1 = NULL;
	int	idx;
	int	numVols;
	int	maxIdx;				/* Max for idxEnd */
	int	volFreq1 = (int) volCurve1->fBasis;

	/* Check that refdate of zc and bv curve agree 
	 */
	if (REFDATE != volCurve1->fBaseDate) {
	        GtoErrMsg("%s: calib timeline ref date %s != "
			"vol curve base date %s\n",
			routine,
			GtoFormatDate(REFDATE),
			GtoFormatDate(volCurve1->fBaseDate));
		goto done;
	}

	/* compute coefficients */
	/* compute coefficients */
	if (VnfmComputeCoeff(that) != SUCCESS)
		goto done;


	if (finalMaturity) {
		/* Final maturity calibration */
		if (GtoDayCountFraction(REFDATE, finalMatDate,
				GTO_ACT_365F, &finalMatYears) IS FAILURE)
			goto done;

		DrlDoubleArrayFloorIdx(TT, NDATES, finalMatYears, &maxIdx);
		/*VnfmFloorIdx(that, finalMatYears, &maxIdx);*/

		idxEnd = MIN(maxIdx, idxEnd);
	}


	/* Get benchmarks vols 
	 */
	numVols = idxEnd-idxStart+1;
	if (numVols <= 0)
	{
	    GtoErrMsg("%s: Calib start index (%d) > calib end index (%d).\n",
		  routine, idxStart, idxEnd);
	    goto done;
	}

	ASSERT_OR_DONE((tMat1 = NEW_ARRAY(double, numVols)) != NULL);
	ASSERT_OR_DONE((vol1  = NEW_ARRAY(double, numVols)) != NULL);
	ASSERT_OR_DONE((freq1 = NEW_ARRAY(int,    numVols)) != NULL);


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: finalMaturity=%d volMat1=%f volFreq1=%d\n"\
		"\tidxStart=%d, idxEnd=%d\n",\
		routine, finalMaturity, volMat1, volFreq1,\
		idxStart, idxEnd));
#endif

    if (NOT finalMaturity)	/* CMS */
    {
	if (GtoDateIntervalToYears(&volMat1, &cmsMatYears) IS FAILURE)
	    goto done;

	/* Interpolate vols
	 */
	for (idx = 0; idx < numVols; idx++)
	{
	    tMat1[idx] = cmsMatYears;
	    freq1[idx] = volFreq1;
	    if (GtoInterpRate(that->fDate[idx+idxStart], volCurve1,
			      GTO_LINEAR_INTERP, vol1+idx) IS FAILURE)
		goto done;	/* Failed */
	}
    } 
    else			/* Final maturity (diagonal) */
    {
	double	S;
	for (idx = 0; idx < numVols; idx++)
	{
	    if (GtoDayCountFraction(REFDATE, that->fDate[idx+idxStart], 
				    GTO_ACT_365F, &S) IS FAILURE)
		goto done;	/* Failed */
	    tMat1[idx] = finalMatYears - S;
	    freq1[idx] = volFreq1;
	    if (GtoInterpRate(that->fDate[idx+idxStart], volCurve1,
			      GTO_LINEAR_INTERP, vol1+idx) IS FAILURE)
		goto done;	/* Failed */

	    /* check maturity long enough (1w+1d) */
	    if (tMat1[idx] <= 2.1970e-2) 
	    {
		idxEnd = idx-1+idxStart;
		break;
	    }
	} /* for */
    } 


	/* Now do the calibration.  */
	if (VnfmVolCalib1VArbitrary(
		that,
		idxStart,
		idxEnd,
		&tMat1[-idxStart], &freq1[-idxStart], &vol1[-idxStart],
		LOGVOL,
		TRUE,
		NULL) != SUCCESS)
			goto done;




	/* made it through OK */
	status = SUCCESS;
done:
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: final parameters\n", routine);\
	VnfmFpWrite(that, vnfmFpLog);\
	VnfmPrintCoeff(that, vnfmFpLog);\
	DrlFPrintf(vnfmFpLog, "\n#\n# %s:END\n", routine));
#endif

	/* free allocated memory */
	if (tMat1) FREE(tMat1);
	if (vol1) FREE(vol1);
	if (freq1) FREE(freq1);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed.\n", routine);
	}
	return(status);
}




