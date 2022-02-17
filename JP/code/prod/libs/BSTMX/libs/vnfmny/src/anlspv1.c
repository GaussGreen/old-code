/****************************************************************
 * Module:	VNFM
 * Submodule:	CALI
 * File:	
 * Function:	
 * Author:	Christian Daher, David Liu
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "macros.h"
#include "date_sup.h"
#include "convert.h"
#include "cerror.h"

#include "drloptio.h"			/* DrlBSTo2QVol() */
#include "drlmem.h"
#include "drlio.h"			/* DrlFPrintf */
#include "drlts.h"			/* DrlTCurveFwdRate */

#define	_vnfm_SOURCE
#include "vnfmanly.h"

#undef __DEBUG__

#if defined(_WINDLL) 
# undef __DEBUG__
#endif

	int	_VnfmCalibSpotVolEqualFailedExpIdx;
	double	_VnfmCalibSpotVolEqualFailedExp,
		_VnfmCalibSpotVolEqualFailedMat,
		_VnfmCalibSpotVolEqualAdjFailedVol;

/*f-------------------------------------------------------------
 * Bootstrap 1D spot volatility (arbitrary, simple version).
 *                                                             
 * <br><br>
 * Low-level routine for calibration of a single time-dependent
 * spot volatility in an N-factor model (spot volatilities
 * of factors are all equal).
 * On entry, "that" contains the model parameters.
 * The routine takes as input an array of swaption volatilities
 * of (increasing) expirations corresponding exactly to the timeline
 * of the structure "that" (argument "fDate" of the structure "VnfmData").
 * The calibration starts at timeline index "idxStart" and
 * end at "idxEnd" (both included).
 * The input arrays "tMat1", "freq1" and "vol1" should have length
 * of at least "idxEnd" and contain respectively the swaptions
 * maturities, frequencies (0,1,2,4,12) and volatilities
 * (the elements with index between 0 and "idxStart"$-1$ are
 * disregarded).

 * {\it Remark: the timeline indices "idxStart" and "idxEnd"
 * correspond to the expirations of the options to be calibrated.
 * For the option od timeline index "idx", the spot volatility
 * of index "idx-1" is fitted.}
 * The routine returns 0 if the bootstrapping can be done.
 * If a failure occurs (because of a too low input volatility),
 * the bootstrapping stops and the routine exits with an nonzero
 * error code. In this case, the global variables
 * \begin{verbatim}
 *  extern int	   _VnfmCalibSpotVolEqualFailedExpIdx;
 *  extern double  _VnfmCalibSpotVolEqualFailedExp,
 *                 _VnfmCalibSpotVolEqualFailedMat,
 *		   _VnfmCalibSpotVolEqualAdjFailedVol;
 * \end{verbatim}
 * contain the expiration and maturity (in years) of the 
 * points where the failure coccurs.
 * An error message is written to the error log only
 * if the flag "errMsgFlag" is set to TRUE.
 */

DLL_EXPORT(int)
VnfmVolCalib1VArbitrary(
	VnfmData *that,		/* (I/O) model parameters */
	int idxStart,		/* (I) first timeline idx calibrated */
	int idxEnd,		/* (I) last timeline idx calibrated */
	double *tMat1,		/* (I) array of swaption mat [0..idxEnd] */
	int *freq1,		/* (I) array of swaption freq [0..idxEnd] */
	double *vol1,		/* (I) array of swaption vol [0..idxEnd] */
	KVolType vType,		/* (I) LOGVOL, NORMVOL */
	int errMsgFlag,		/* (I) write to errlog iff TRUE */
	double *resValue)	/* (O) <0 if calib OK, >0 otherwise */
{
static	char		routine[] = "VnfmVolCalib1VArbitrary";
	int		status = FAILURE;
	int		nExp;
	double		*tReset = NULL;

	/*
	 *
	 */
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: start.\n", routine))
#endif
	idxStart = MAX(idxStart, 1);
	idxEnd   = MIN(idxEnd, that->fNDates-1);

	if ((tReset  = DrlDoubleVectAlloc(0, that->fNDates-1)) == NULL)
		goto done;

	/* Rate resets at option expiration time  */
	for (nExp=idxStart; nExp <= idxEnd; nExp++) 
	    tReset[nExp] = that->fTime[nExp];


	
	if (VnfmVolCalib1VArbitraryNew( that,
					idxStart,
					idxEnd,	
					tReset,
					tMat1,
					freq1,
					vol1,
					vType,
					errMsgFlag,	
					resValue) != SUCCESS)
		goto done;

	/* made it through OK */
	status = SUCCESS;
done:

	if ((errMsgFlag != FALSE) && (status != SUCCESS)) {
	    GtoErrMsg("%s: failed\n", routine);
	}

	DrlDoubleVectFree(tReset, 0, (that!=NULL ? that->fNDates-1:0));

	return(status);
}



/*f-------------------------------------------------------------
 * Bootstrap 1D spot volatility (arbitrary, general version).
 *                                                             
 * <br><br>
 * Low-level routine for calibration of a single time-dependent
 * spot volatility in an N-factor model (spot volatilities
 * of factors are all equal).
 * On entry, "that" contains the model parameters.
 * The routine takes as input an array of swaption volatilities
 * of (increasing) expirations corresponding exactly to the timeline
 * of the structure "that" (argument "fDate" of the structure "VnfmData").
 * The calibration starts at timeline index "idxStart" and
 * end at "idxEnd" (both included).
 * The input arrays "reset1", "tMat1", "freq1" and "vol1" should have 
 * length of at least "idxEnd" and contain respectively the swaptions
 * reset times (in years), maturities, frequencies (0,1,2,4,12) and 
 * volatilities (the elements with index between 0 and "idxStart"$-1$ are
 * disregarded).

 * {\it Remark: the timeline indices "idxStart" and "idxEnd"
 * correspond to the expirations of the options to be calibrated.
 * For the option of timeline index "idx", the spot volatility
 * of index "idx-1" is fitted.}
 * The routine returns 0 if the bootstrapping can be done.
 * If a failure occurs (because of a too low input volatility),
 * the bootstrapping stops and the routine exits with an nonzero
 * error code. In this case, the global variables
 * \begin{verbatim}
 *  extern int	   _VnfmCalibSpotVolEqualFailedExpIdx;
 *  extern double  _VnfmCalibSpotVolEqualFailedExp,
 *                 _VnfmCalibSpotVolEqualFailedMat,
 *		   _VnfmCalibSpotVolEqualAdjFailedVol;
 * \end{verbatim}
 * contain the expiration and maturity (in years) of the 
 * points where the failure coccurs.
 * An error message is written to the error log only
 * if the flag "errMsgFlag" is set to TRUE.
 */

DLL_EXPORT(int)
VnfmVolCalib1VArbitraryNew(
	VnfmData *that,		/* (I/O) model parameters */
	int idxStart,		/* (I) first timeline idx calibrated */
	int idxEnd,		/* (I) last timeline idx calibrated */
	double *reset1,		/* (I) array of swaption resets [0..idxEnd] */
	double *tMat1,		/* (I) array of swaption mat [0..idxEnd] */
	int *freq1,		/* (I) array of swaption freq [0..idxEnd] */
	double *vol1,		/* (I) array of swaption vol [0..idxEnd] */
	KVolType vType,		/* (I) LOGVOL, NORMVOL */
	int errMsgFlag,		/* (I) write to errlog iff TRUE */
	double *resValue)	/* (O) <0 if calib OK, >0 otherwise */
{
static	char		routine[] = "VnfmVolCalib1VArbitraryNew";
	int		status = FAILURE;
	int		i, j,
			nDim = NDIM,
			jIdx, rIdx,
			nExp, n;
	double		x1[VNFM_NDIMMAX],
			vol12q,		/* target vol 2q adjusted */
			v0, v1, K0, spv,
			S, dt,
			factor,
			fwdRate,
			lambda,
			dVol = 1e9;

	/*
	 *
	 */
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: start.\n", routine))
#endif
	idxStart = MAX(idxStart, 1);
	idxEnd   = MIN(idxEnd, that->fNDates-1);

	if (idxStart > idxEnd) {
	    /* nothing to do */
#ifndef	NO_LOGGING
	    GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: idxStart=%d "\
		"\tidxEnd=%d : nothing to do.\n", routine, idxStart, idxEnd))
#endif
	    status = SUCCESS; goto done;
	}


	VnfmComputeCoeff(that);
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "\tidxStart=%d (S=%lf)"\
		"\tidxEnd=%d (S=%lf)\n",\
		idxStart, TT[idxStart], idxEnd, TT[idxEnd]))
#endif

	/* main calibration loop */
	for (nExp=idxStart; nExp <= idxEnd; nExp++) {
	    /* nExp = timeline exp idx of base vol point to be fitted
	     * n = idx of base vol point to be modified
	     */
	    n = nExp-1;

	    /*this is fucked up ! if (nExp == that->fNDates-1) {
		spv = SIGMA[0][n];
		goto bootstrap_done;
	    }*/

	    S = that->fTime[nExp];
	    dt = S - that->fTime[n];

	    /* check maturity */
	    /*if (tMat1[nExp] < 1.92307692e-2) {*/
	    if (tMat1[nExp] < 0e0) {
		GtoErrMsg("%s: (%s) calibrated rate "
		    "has too low maturity (%lf yrs)\n",
		    routine, GtoFormatDate(that->fDate[nExp]), tMat1[nExp]);
		goto done;
	    }
	    /* rate must be at least one week */
	    tMat1[nExp] = MAX(tMat1[nExp], 1.92307692e-2);

	    /* check reset >= option expiration */
	    if (reset1[nExp] < S) {
		GtoErrMsg("%s: (%s) rate can not reset (%lf) before "
		    "option expiration (%lf) \n",
		    routine, reset1[nExp], S);
		goto done;
	    }


	    /* compute the Q coefficients */
	    if (freq1[nExp] > 0) {
		VnfmB(that, reset1[nExp], tMat1[nExp], freq1[nExp], &fwdRate, x1);
	    } else {
		VnfmQ(that, reset1[nExp], reset1[nExp]+tMat1[nExp], &fwdRate, x1);
	    }

	    /* compute coefficients */
	    v0 = 0e0;
	    K0 = 0e0;
	    for (i=0; i<=nDim-1; i++)
	    for (j=i; j<=nDim-1; j++) {
		jIdx = JIDX(i, j);
		rIdx = RHOIDX(i, j);
		lambda = BETA[i] + BETA[j];
		factor = exp(-(reset1[nExp] - S)*lambda);

		v0 += (i==j ? 1e0 : 2e0) *
			x1[i] * x1[j] *
			factor *
			exp(-dt*lambda) * JJ[jIdx][n];

		K0 += (i==j ? 1e0 : 2e0 * RHO[rIdx][n]) *
			x1[i] * x1[j] * ALPHA[i] * ALPHA[j] *
			factor *
	    		L(dt, lambda);
	    }

	    /*
	     * Adjust market volatility for 2q duration
	     */
	     IF_FAILED_DONE( DrlATMTo2QVol(
			S,
			fwdRate,
		    	vol1[nExp],	/* volatility */
			vType,
			fwdRate,
			that->fQLeft,
			that->fQRight,
			that->fQFShift,
			"C",
	    		&vol12q));

	    /* market volatility */
	    v1 = SQR(vol12q) * S;

#ifndef	NO_LOGGING
	    GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
		"[%3d] S=%lf dt=%lf tMat=%lf freq=%d vol=%10.6f\n",\
		nExp, S, dt, tMat1[nExp], freq1[nExp], vol12q))
	    GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "\tv0=%lf\tv1=%lf\tK0=%lf\n",\
		v0, v1, K0));
#endif

	    dVol = MIN(v1-v0, dVol);


	    /* Check if bootstrapping possible
	     */
	    if (v1 - v0 <= 0e0) {
		/* failed */
		if (errMsgFlag != FALSE) {
		    GtoErrMsg("%s: failure : ", routine);
		    GtoErrMsg(" tExp=%8.4f (%10s)"
		    	      " tMat=%8.4f ",
				S, GtoFormatDate(that->fDate[nExp]),
				tMat1[nExp]);

		    /* The vol12q and v0 are all bp vols */
		    if (vType == NORMVOL)
		    	GtoErrMsg(" mktVol=%.4f%% 2qVol=%.4f%% modVol=%.4f%% "
			      "diff=%.4f%%\n",
			vol1[nExp]*1e2, sqrt(v1/S)*1e2, sqrt(v0/S)*1e2,
				(sqrt(v1/S)-sqrt(v0/S))*1e2);
		    else	/* convert to % vol */
		    	GtoErrMsg(" mktVol=%.4f%% 2qVol=%.4f%% modVol=%.4f%% "
			      "diff=%.4f%%\n",
			vol1[nExp]*1e2, sqrt(v1/S)/fwdRate*1e2, 
			sqrt(v0/S)/fwdRate*1e2,
			(sqrt(v1/S)-sqrt(v0/S))/fwdRate*1e2);
			
		 }

		 _VnfmCalibSpotVolEqualFailedExpIdx = nExp - 1;
		 _VnfmCalibSpotVolEqualFailedExp = S;
		 _VnfmCalibSpotVolEqualFailedMat = tMat1[nExp];
		 _VnfmCalibSpotVolEqualAdjFailedVol = sqrt(v0 / S) - vol12q;
			
		 goto done;
		/*spv = 0e0;*/
	    } else {
		spv = sqrt((v1 - v0) / K0);
	    }

bootstrap_done:
	    /* update spot vola arrays */
	    for (i=0; i<=nDim-1; i++)
		SIGMA[i][n] = spv;

	    /* update coefficients J for next step */
	    n++;
	    dt = that->fTime[n] - that->fTime[n-1];
	    for (i=0; i<=nDim-1; i++)
	    for (j=i; j<=nDim-1; j++) {
		jIdx = JIDX(i, j);
		rIdx = RHOIDX(i, j);
		lambda = BETA[i] + BETA[j];
		JJ[jIdx][n] =
		    (i==j ? 1e0 : RHO[rIdx][n-1]) *
		    ALPHA[i] * ALPHA[j] *
		    SIGMA[i][n-1] * SIGMA[j][n-1] * L(dt, lambda) +
		    exp(-dt * lambda) * JJ[jIdx][n-1];
	    }

	}


	/* fill the remaining spot vol */
	for (; n<=that->fNDates-1; n++) {
	    for (i=0; i<=nDim-1; i++) {
		SIGMA[i][n] = SIGMA[i][n-1];
	    }
	}



	if (resValue != NULL) {
		*resValue = dVol;
	}

	/* tesing of accuracy */
	VnfmComputeCoeff(that);

	/* made it through OK */
	status = SUCCESS;
done:
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmPrintCoeff(that, vnfmFpLog))
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: done.\n", routine))
#endif

	if ((errMsgFlag != FALSE) && (status != SUCCESS)) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*f-------------------------------------------------------------
 * Bootstrap failure retrieve.
 *                                                             
 * <br><br>
 * Routine returns maturity and expiration, in years
 * of rate which failed to calibrate
 * during a call to <i> VnfmVolCalib1VArbitrary</i>
 * (must be called immediately after).
 */

DLL_EXPORT(int)
VnfmVolCalib1VFailureId(
	VnfmData *that,		/* (I) Model parameters */
	double *failedExp,	/* (O) Expiration which failed (in years) */
	double *failedMat)	/* (O) Fwd maturity which failed (in years) */
{
	*failedExp = _VnfmCalibSpotVolEqualFailedExp;
	*failedMat = _VnfmCalibSpotVolEqualFailedMat;
	return(SUCCESS);
}




/*f-------------------------------------------------------------
 * Bootstrap minimum spot volatility check.
 *                                                             
 * <br><br>
 * Routine check that the minimum spot volatility of the instataneous
 * continuously compounded rate is larger than 
 * a specifed <i> minSpotVol</i>. If <i> spotVolRatioFlag</i> is set,
 * it also checks that the ratio of spot volatilities is within the 
 * specified range.  If it is the case, it returns SUCCESS.
 * Otherwise it returns FAILURE and sets the global variables 
 * \begin{verbatim}
 *  extern int     _VnfmCalibSpotVolEqualFailedExpIndx;
 *  extern double  _VnfmCalibSpotVolEqualFailedExp,
 *                 _VnfmCalibSpotVolEqualFailedMat;
 * \end{verbatim}
 * at the point where the minimum spot volatility or ratio is attained.
 * This routine can be called immediately after a successful return
 * from <i> VnfmVolCalib1VArbitrary</i> to check that the spot volatility 
 * is large enough every where.
 */


DLL_EXPORT(int)
VnfmVolCalib1VArbitraryMinSpotVol(
	VnfmData *that,		/* (I) model parameters */
	int idxStart,		/* (I) first timeline idx calibrated */
	int idxEnd,		/* (I) last timeline idx calibrated */
	double *tMat1,		/* (I) array of swaption mat [0..idxEnd] */
	int *freq1,		/* (I) array of swaption freq [0..idxEnd] */
	double *vol1,		/* (I) array of swaption vol [0..idxEnd] */
	int errMsgFlag,		/* (I) write to errlog iff TRUE */
	int    spotVolRatioFlag,/* (I) spot vol ratio check flag  */
	double spotVolRatio,    /* (I) spot vol ratio check flag  */
	double minSpotVol)	/* (I) minmum spot volatility */
{
static	char		routine[] = "VnfmVolCalib1VArbitraryMinSpotVol";
	int		status = FAILURE;
	int		i, j,
			nDim = NDIM,
			rIdx,
			tpIdx, tpIdxMax = 0;
	double		S,
			v00, v0, v0Max;

	/*
	 *
	 */
	idxStart = MAX(idxStart, 1);
	idxEnd   = MIN(idxEnd, that->fNDates-1);

	if (idxStart > idxEnd) {
	    /* nothing to do */
	    status = SUCCESS; goto done;
	}

	v0Max = minSpotVol;

	for (tpIdx=idxStart; tpIdx <= idxEnd; tpIdx++) {
	    S = that->fTime[tpIdx];

	    /* compute spot volatility */
	    v0 = 0e0;
	    for (i=0; i<=nDim-1; i++)
	    for (j=i; j<=nDim-1; j++) {
		rIdx = RHOIDX(i, j);
		v0 += (i==j ? 1e0 : RHO[rIdx][tpIdx-1]) *
			ALPHA[i] * ALPHA[j] *
			SIGMA[i][tpIdx-1] * SIGMA[j][tpIdx-1];
	    }
	    v0 = sqrt(v0);
	
	    if (tpIdx == idxStart)
		v00 = v0;

	    if (v0 < v0Max) {
		v0Max = v0;
		tpIdxMax = tpIdx;
		_VnfmCalibSpotVolEqualFailedExpIdx = tpIdx - 1;
		_VnfmCalibSpotVolEqualFailedExp = S;
		_VnfmCalibSpotVolEqualFailedMat = tMat1[tpIdx];
		goto done;
	    }

	    if (spotVolRatioFlag){
		if ((v0/v00) > spotVolRatio ||
		    (v0/v00) < 1.0/spotVolRatio) {
			_VnfmCalibSpotVolEqualFailedExpIdx = tpIdx - 1;
			_VnfmCalibSpotVolEqualFailedExp = S;
			_VnfmCalibSpotVolEqualFailedMat = tMat1[tpIdx];
			goto done;
	    	}
	    }
	}

	/* made it through OK */
	status = SUCCESS;
done:
	if ((errMsgFlag != FALSE) && (status != SUCCESS)) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}



/*---------------------------------------------------------------
 * Bootstrap 1D spread spot volatility.
 *                                                             
 * <br><br>
 * Calibration of spread spot volatilities given
 * the base volatilities of spread.  
 * Analagous to VnfmVolCalib1VArbitrary, but with all the Q and B 
 * coeff set to backbone function.  The zero curve in "that" is obsolete 
 * in Q and B calculation, but used instead to store the forward spread curve 
 * for the calculation of 2q vol correction.
 */
DLL_EXPORT(int)
VnfmSpreadVolCalib1VArbitrary(
	VnfmData *that,		/* (I/O) model parameters */
	int idxStart,		/* (I) first timeline idx calibrated */
	int idxEnd,		/* (I) last timeline idx calibrated */
	double *vol1,		/* (I) array of spread base vol [0..idxEnd] */
	KVolType vType,         /* (I) LOGVOL, NORMVOL */
	int errMsgFlag,		/* (I) write to errlog iff TRUE */
	double *resValue)	/* (O) <0 if calib OK, >0 otherwise */
{
static	char		routine[] = "VnfmSpreadVolCalib1VArbitrary";
	int		status = FAILURE;
	int		i, j,
			nDim = NDIM,
			jIdx, rIdx,
			m, nExp, n;
	double		v0, v1, K0, spv,
			S, dt,
			lambda, dVol = 1e9;

	TDate		interpDate;	/* Interpolated date */
	double		spread;		/* Interpolated forward spread */

	double		vol12q;         /* target vol 2q adjusted */

	double		x, x2;		/* backbone */

	/*
	 *
	 */
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: start.\n", routine))
#endif
	idxStart = MAX(idxStart, 1);
	idxEnd   = MIN(idxEnd, that->fNDates-1);

	if (idxStart > idxEnd) {
	    /* nothing to do */
#ifndef	NO_LOGGING
	    GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: idxStart=%d "\
		"\tidxEnd=%d : nothing to do.\n", routine, idxStart, idxEnd))
#endif
	    return SUCCESS;
	}


	VnfmComputeCoeff(that);
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "\tidxStart=%d (S=%lf)"\
		"\tidxEnd=%d (S=%lf)\n",\
		idxStart, TT[idxStart], idxEnd, TT[idxEnd]))
#endif

	/* main calibration loop */
	for (nExp=idxStart; nExp <= idxEnd; nExp++) {
	    /* nExp = timeline exp idx of base vol point to be fitted
	     * n = idx of base vol point to be modified
	     */
	    n = nExp-1;

	    S = that->fTime[nExp];
	    dt = S - that->fTime[n];

	    /* 	
	     * Compute the spread at S 
	     */
	    IF_FAILED_DONE(GtoTDateAdvanceYears(
			   that->fZcCurve->fBaseDate, 
			   S, 
			   &interpDate));

	    /* Interpolate the spread from spread curve
	     */
	    IF_FAILED_DONE(GtoInterpRate(
			   interpDate,
			   that->fZcCurve,
                           GTO_LINEAR_INTERP,
                           &spread)); 

	    /*
	     * Check for zero spread
	     */
	    if (IS_ALMOST_ZERO(spread) &&
		IS_ALMOST_ZERO(that->fBackBoneQ - 1.))
		GtoErrMsg("%s: invalid case with lognormal backbone (q=1) "
			 "and zero initial spread on %s.\n",
			 routine,
			 GtoFormatDate(interpDate));

	    /*
 	     * compute the backbone 
	     */
	    x = FN(spread, that->fBackBoneQ);
	    x2 = x * x;

	    /* compute coefficients */
	    v0 = 0e0;
	    K0 = 0e0;
	    for (i=0; i<=nDim-1; i++)
	    for (j=i; j<=nDim-1; j++) {
		jIdx = JIDX(i, j);
		rIdx = RHOIDX(i, j);
		lambda = BETA[i] + BETA[j];

		v0 += (i==j ? 1e0 : 2e0) *
			x2 * exp(-dt*lambda) * JJ[jIdx][n];

		K0 += (i==j ? 1e0 : 2e0 * RHO[rIdx][n]) *
			x2 *
			ALPHA[i] * ALPHA[j] *
	    		L(dt, lambda);
	    }

	    /*
	     * Adjust market volatility for 2q duration
	     */

	    /* Campute the 2q vol adjustment
 	     */
	    IF_FAILED_DONE( DrlATMTo2QVol(
			S,
			spread,
		    	vol1[nExp],	/* volatility */
			vType,
			spread,
			that->fQLeft,
			that->fQRight,
			that->fQFShift,
			"C",
	    		&vol12q));

	    /* market volatility */
	    v1 = SQR(vol12q) * S;


	    dVol = MIN(v1-v0, dVol);

	    /* Check if bootstrapping possible
	     */
	    if (v1 - v0 <= 0e0) {
		/* failed */
		if (errMsgFlag != FALSE) {
		    GtoErrMsg("%s: failure : ", routine);
		    if (vType == NORMVOL)
		    	GtoErrMsg(" tExp=%8.4f (%10s)\n"
		    	      "mktVol=%.4f%% 2qVol=%.4f%% modVol=%.4f%%\n"
			      "diff=%.4f%%\n",
				S, 
				GtoFormatDate(that->fDate[nExp]),
				vol1[nExp]*1e2,
				sqrt(v1/S)*1e2, sqrt(v0/S)*1e2,
				(sqrt(v1/S)-sqrt(v0/S))*1e2);
		    else
		    	GtoErrMsg(" tExp=%8.4f (%10s)\n"
		    	      "mktVol=%.4f%% 2qVol=%.4f%% modVol=%.4f%%\n"
			      "diff=%.4f%%\n",
				S, 
				GtoFormatDate(that->fDate[nExp]),
				vol1[nExp]*1e2/spread,
				sqrt(v1/S)*1e2/spread, sqrt(v0/S)*1e2/spread,
				(sqrt(v1/S)-sqrt(v0/S))*1e2/spread);
		    goto done;
		 }
	    } else {
		spv = sqrt((v1 - v0) / K0);
	    }


	    /* update spot vola arrays */
	    for (i=0; i<=nDim-1; i++)
		SIGMA[i][n] = spv;

	    /* update coefficients J for next step */
	    n++;
	    dt = that->fTime[n] - that->fTime[n-1];
	    for (i=0; i<=nDim-1; i++)
	    for (j=i; j<=nDim-1; j++) {
		jIdx = JIDX(i, j);
		rIdx = RHOIDX(i, j);
		lambda = BETA[i] + BETA[j];
		JJ[jIdx][n] =
		    (i==j ? 1e0 : RHO[rIdx][n-1]) *
		    ALPHA[i] * ALPHA[j] *
		    SIGMA[i][n-1] * SIGMA[j][n-1] * L(dt, lambda) +
		    exp(-dt * lambda) * JJ[jIdx][n-1];
	    }

	}


	/* fill the remaining spot vol */
	for (; n<=that->fNDates-1; n++) {
	    for (i=0; i<=nDim-1; i++) {
		SIGMA[i][n] = SIGMA[i][n-1];
	    }
	}



	if (resValue != NULL) {
		*resValue = dVol;
	}

	/* tesing of accuracy */
	VnfmComputeCoeff(that);

	/* made it through OK */
	status = SUCCESS;


done:
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmPrintCoeff(that, vnfmFpLog))
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: done.\n", routine))
#endif

	if ((errMsgFlag != FALSE) && (status != SUCCESS)) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}
