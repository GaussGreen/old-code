/****************************************************************
 * Module:	VNFM
 * Submodule:	ANLY
 * File:	
 * Function:	Computes fwd/fut adjustment 
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <math.h>
#include <string.h>

#include "cerror.h"
#include "ldate.h"
#include "macros.h"
#include "convert.h"
#include "date_sup.h"
#include "fltrate.h"
#include "gtomat.h"

#include "drlsort.h"		/* DrlDoubleArrayFloorIdx() */
#include "drlio.h"		/* DrlFPrintf() */
#define	_MYMEM
#ifdef	_MYMEM
#include "drlmem.h"		/* DrlDoubleMatrAlloc() */
#endif

#define	_vnfm_SOURCE
#include "vnfmanly.h"

#undef	__DEBUG__

/*f-------------------------------------------------------------
 * Compute the forward/future adjustment curve.
 *                                                             
 * <br><br>
 * Computes the forward/future adjusments of a rate
 * of forward maturity "rateMat", frequency "rateFreq" (0,1,2,4,12),
 * day count convention "rateDayCount".
 * The adjustments are computed at every date of the "resetDates" array
 * The result (forward rate computed form the input zero curve,
 * convexity adjustment and MTM adjustment)
 * are put in the three arrays "fwdRate", "cvxAdj",
 * "mtmAdj" (all arrays are assumed to have length at least "nResetDates")
 */

int
VnfmFwdFutAdjustment(
	VnfmData *that,		/* (I) model parameters */
	TCurve *zcCurve,	/* (I) zero curve */
	TDateInterval rateMat,	/* (I) rate maturity */
	int rateFreq,		/* (I) rate frequency */
	long rateDayCount,	/* (I) day count convention */
	int nResetDates,	/* (I) # of desired dates */
	TDate *resetDates,	/* (I) array of desired dates [0..nResetD] */
	double *fwdRate,	/* (O) array of fwd rates (or NULL) */
	double *cvxAdj,		/* (O) array of cvx adj (or NULL) */
	double *mtmAdj)		/* (O) array of mtm adj (or NULL) */
{
static	char	routine[] = "VnfmFutAdjRateCompute";
	int	status = FAILURE;
	double	**faInteg = NULL,/* contains integrals to be interp */
		tMat,		/* rate maturity */
		tReset,		/* rate reset */
		fRate, fAdjustedRate, yieldVol,
		futAdjCoeff;
	int	i, nDim = NDIM;
	TFloatRate	floatRate;


	/* this routine only works in the lognormal case */
 	if (!(IS_ALMOST_ZERO(that->fBackBoneQ - 1e0))) {
 		GtoErrMsg("%s: only supports lognormal backbone (%lf != 1e0).\n",
 			  routine, that->fBackBoneQ);
		goto done;
	}


	/*
	 *
	 */
	if (VnfmComputeCoeff(that) != SUCCESS)
		goto done;
#ifdef	_MYMEM
	if ((faInteg = DrlDoubleMatrAlloc(0, nDim*nDim-1, 0, that->fNDates-1))
		== NULL) goto done;
#else
	faInteg = (double**) GtoArray2DNew(nDim*nDim,
			that->fNDates, (int) sizeof(double));
	if (faInteg == NULL) goto done;
#endif

	if (VnfmFutAdjIntegralsCompute(that, faInteg)
		!= SUCCESS) goto done;

	if (GtoDateIntervalToYears(&rateMat, &tMat) != SUCCESS)
		goto done;

	/* create floating rate structure */
	floatRate.matInterval = rateMat;
	if (rateFreq != 0) {
	    if (GtoFreq2TDateInterval((long) rateFreq, &floatRate.payInterval)
		!= SUCCESS) goto done;
	} else {
	    floatRate.payInterval = rateMat;
	}
	floatRate.dayCountConv = rateDayCount;
	GTO_SET_ADJ_INTERVAL_DAYS(floatRate.spotOffset, 0);
	floatRate.spread = 0e0;
	floatRate.weight= 0e0;



	for (i=0; i<=nResetDates-1; i++) {
		/* compute time to reset */
		if (GtoDayCountFraction(REFDATE, resetDates[i],
			GTO_ACT_365F, &tReset) != SUCCESS) goto done;

		/* compute fwd rate */
#ifdef	CLIB7X
		if (GtoForwardRate(
			zcCurve,
			&floatRate,
			resetDates[i],
			GTO_LINEAR_INTERP,
			&fRate) != SUCCESS)
				goto done;
#else
		if (GtoForwardRate(
			zcCurve,
			GTO_LINEAR_INTERP,
			&floatRate,
			resetDates[i],
			&fRate) != SUCCESS)
				goto done;
#endif

		if (fwdRate != NULL) fwdRate[i] = fRate;

		if (cvxAdj != NULL) {
		    /* compute rate vol */
		    if (VnfmAvgQBVol(that,
			0e0, tReset, tReset, tMat, rateFreq, LOGVOL,
			&yieldVol) != SUCCESS) goto done;

		    /* compute adj rate */
#ifdef	CLIB7X
		    if (GtoAdjustedFwdRate(
			zcCurve,
			resetDates[i],
			&floatRate,
			yieldVol,
			0e0,		/*(I) Fwd MM rate for delay payment */
			resetDates[i],	/*(I) Payment date */
			0e0,		/*(I) Vol for the fwd MM rate */
			GTO_ACT_360,	/*(I) Day count conv for MM rate */
			1e0,		/*(I) Corr fwdRate and fwdMMRate */
			zcCurve->fBaseDate,
			&fAdjustedRate) != SUCCESS)
				goto done;
#else
		    if (GtoForwardRateAdjusted(
			zcCurve,
			zcCurve,
			GTO_LINEAR_INTERP,
			resetDates[i],
			resetDates[i],	/*(I) Payment date */
			&floatRate,
			zcCurve->fBaseDate,
			yieldVol,
			0e0,		/* (I) Vol for the fwd MM pay rate */
			1e0,		/*(I) Corr fwdRate and fwdMMRate */
			&fAdjustedRate) != SUCCESS)
				goto done;
#endif

		    cvxAdj[i] = fAdjustedRate - fRate;

		}


		/* compute mtm adj */
		if (mtmAdj != NULL) {
		    if (VnfmFutAdjRateComputeAdjFactor(that, faInteg, tReset,
			tMat, rateFreq, &futAdjCoeff) != SUCCESS)
				goto done;
		    mtmAdj[i] = fRate * futAdjCoeff;
		}

	}

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s:\n", routine);\
	for (i=0; i<=nResetDates-1; i++) {\
	    DrlFPrintf(vnfmFpLog, "[%3d/%3d] %10s  f=%lf  c=%lf  m=%lf\n",\
		i, nResetDates, GtoFormatDate(resetDates[i]),\
		(fwdRate != NULL ? fwdRate[i] : -1e0),\
		(cvxAdj  != NULL ? cvxAdj[i]  : -1e0),\
		(mtmAdj  != NULL ? mtmAdj[i]  : -1e0));\
	});
#endif



	/* made it through OK */
	status = SUCCESS;
done:
#ifdef	_MYMEM
	DrlDoubleMatrFree(faInteg, 0, nDim*nDim-1, 0, that->fNDates-1);
#else
	GtoArray2DFree((void**) faInteg);
#endif
	if (status!= SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}





/*f-------------------------------------------------------------
 * Compute the forward/future adjustment at a single date.
 *                                                             
 * <br><br>
 * Computes the forward/futures adjustment
 * a forward rate
 * of forward maturity "tMat", time to reset "tReset" and
 * frequency "freq" (1,2,4,12).\\
 * <b> Warning.</b> The routine assumes that the internal
 * coefficients have been computed (using "VnfmComputeCoeff")
 * and also expects the array of futures adjustment volatility
 * integrals "faInteg" computed using <i> VnfmFutAdjIntegralsCompute</i>.
 */

int
VnfmFutAdjRateComputeAdjFactor(
	VnfmData *that,		/* (I) model parameters */
	double **faInteg,	/* (I) contains integrals to be interp */
	double tReset,		/* (I) rate reset */
	double tMat,		/* (I) rate maturity */
	int freq,		/* (I) rate frequency */
	double *retVal)		/* (O) adjustment */
{
/*static	char	routine[] = "VnfmFutAdjRateCompute";*/
	double	x[VNFM_NDIMMAX],
		faIntegInterp[VNFM_NDIMMAX*VNFM_NDIMMAX],
		adj = 0e0;
	int	i, j, nDim = NDIM, hIdx;

	double	yield;

	/* check */

	/* compute the Q coefficients */
	if (VnfmQBCoeff(that, tReset, tReset, tMat, freq, &yield, x) != SUCCESS)
		return(FAILURE);

	/* compute integrals */
	VnfmFutAdjIntegralsCompute(that, faInteg);

	/* interpolate integrals */
	VnfmFutAdjIntegralsInterp(that, faInteg, tReset, faIntegInterp);


	/* vol computation */
	adj = 0e0;
	for (i=0; i<=nDim-1; i++)
	for (j=0; j<=nDim-1; j++) {
	    hIdx = HIDX(i, j);
	    adj += x[i] * faIntegInterp[hIdx] / yield;
	}


	*retVal = adj;
	return(SUCCESS);
}



/*f-------------------------------------------------------------
 * Compute the H integrals.
 *                                                             
 * <br><br>
 * Computes the volatility integrals
 * <blockquote>
 * H_{ij}(T) = integral_0^S
 * 		alpha_i sigma_i(t) alpha_j sigma_j(t) rho_{ij}(t)
 *		e^{-beta_j(S-t)} A(t,S,beta_i) dt.
 * </blockquote>
 * used ind the forward/futures adjustment (0&lt;= i,j&lt;= Nf)
 * at each point of the time line T=t_0,t_1,....
 * and stores the values in the array <i> faInteg</i>.
 * Should not be called directly in applications
 * (use <i> VnfmFwdFutAdjustment</i> instead.
 * Returns 0 if OK.
 */

int
VnfmFutAdjIntegralsCompute(
	VnfmData *that,		/* (I) model parameters */
	double **faInteg)	/* (O) storage [0..nF*nF-1][0..nDates-1] */
{
static	char	routine[] = "VnfmFutAdjIntegralsCompute";
	int	i, j,
		nDim = NDIM,
		jIdx, rIdx, hIdx,
		n;
	double	fnm1, dt;


	/* This routine only works if fZTime == fTime */
	if (NDATES ISNT NZDATES) {
		GtoErrMsg("%s: # dates (%d) != # zero dates (%d).\n",
			routine, (int) NDATES, (int) NZDATES);
		GtoErrMsg("%s: routine N/A if zero curve dates "
			"not equal to volatility dates (sorry...).\n",
			routine);
		return (FAILURE);
	}
	for (i=0;i<NZDATES;i++) {
		if (TT[i] ISNT ZTT[i]) {
			GtoErrMsg("%s: timept  #%d (%lf) != zero "
				"timept %d (%lf).\n", routine,
				i, TT[i],
				i, ZTT[i]);
			GtoErrMsg("%s: routine N/A if zero curve dates "
				"not equal to volatility dates (sorry...).\n",
				routine);
			return (FAILURE);
		}
	}

	/* compute integrals at each time point */
	for (i=0; i<=nDim-1; i++)
	for (j=0; j<=nDim-1; j++) {
	    jIdx = JIDX(i, j);
	    rIdx = RHOIDX(i, j);
	    hIdx = HIDX(i, j);

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "VnfmComputeH: i=%d j=%d jIdx=%2d  rIdx=%2d\n",\
		i, j, jIdx, rIdx));
#endif
	    n = 0;
	    faInteg[hIdx][n] = 0e0;

	    for (n=1; n<=NDATES-1; n++) {
		dt = TT[n] - TT[n-1];
		fnm1 = ALPHA[i] * ALPHA[j] *
			SIGMA[i][n-1] * SIGMA[j][n-1] *
			(i == j ? 1e0 : RHO[rIdx][n-1]);

		faInteg[hIdx][n] =
		    faInteg[hIdx][n-1] * exp(-dt*BETA[j])
		  + RATE[n-1] * (
		      fnm1 * L(dt, BETA[j]) / (BETA[i] + BETA[j])
		    + (that->fJ[jIdx][n-1] - fnm1 / (BETA[i] + BETA[j]))
		    * exp(-dt*BETA[i]) * L(dt, BETA[i])
		  );
	    }
	}
	return(SUCCESS);
}


/*f-------------------------------------------------------------
 * Interpolate the H intergrals.
 *                                                             
 * <br><br>
 * Interpolates the volatility integrals
 * used ind the forward/futures adjustment.
 * Should not be called directly in applications
 * (use <i> VnfmFwdFutAdjustment</i> instead.
 * Returns 0 if OK.
 */

int
VnfmFutAdjIntegralsInterp(
	VnfmData *that,		/* (I) model parameters */
	double **faInteg,	/* (I) contains integrals to be interp */
	double T,		/* (I) time to perform the interp */
	double *faIntegInterp)	/* (O) should be allocated */
{
static	char	routine[] = "VnfmFutAdjIntegralsInterp";
	int	i, j,
		nDim = NDIM,
		jIdx, rIdx, hIdx,
		n;
	double	fn, dt;

	/* This routine only works if fZTime == fTime */
	if (NDATES ISNT NZDATES) {
		GtoErrMsg("%s: # dates (%d) != # zero dates (%d).\n",
			routine, (int) NDATES, (int) NZDATES);
		GtoErrMsg("%s: routine N/A if zero curve dates "
			"not equal to volatility dates (sorry...).\n",
			routine);
		return (FAILURE);
	}
	for (i=0;i<NZDATES;i++) {
		if (TT[i] ISNT ZTT[i]) {
			GtoErrMsg("%s: timept  #%d (%lf) != zero "
				"timept %d (%lf).\n", routine,
				i, TT[i],
				i, ZTT[i]);
			GtoErrMsg("%s: routine N/A if zero curve dates "
				"not equal to volatility dates (sorry...).\n",
				routine);
			return (FAILURE);
		}
	}


	DrlDoubleArrayFloorIdx(TT, NDATES, T, &n);

	for (i=0; i<=nDim-1; i++)
	for (j=0; j<=nDim-1; j++) {
	    jIdx = JIDX(i, j);
	    rIdx = RHOIDX(i, j);
	    hIdx = HIDX(i, j);

		dt = T - TT[n];
		fn = ALPHA[i] * ALPHA[j] *
			SIGMA[i][n] * SIGMA[j][n] *
			(i == j ? 1e0 : RHO[rIdx][n]);

		faIntegInterp[hIdx] =
		    faInteg[hIdx][n] * exp(-dt*BETA[j])
		  + RATE[n] * (
		      fn * L(dt, BETA[j]) / (BETA[i] + BETA[j])
		    + (that->fJ[jIdx][n] - fn / (BETA[i] + BETA[j]))
		    * exp(-dt*BETA[i]) * L(dt, BETA[i])
		  );
	}
	return(SUCCESS);
}







