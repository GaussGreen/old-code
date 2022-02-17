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
#include <float.h>

#include "macros.h"
#include "cerror.h"
#include "date_sup.h"
#include "convert.h"

#include "drllno.h"
#include "drlio.h"
#include "drlmem.h"
#include "drlsort.h"			/* DrlDoubleArrayFloorIdx() */
#include "drlts.h"			/* DrlTCurveFpWrite() */
#include "drlstr.h"			/* DrlFloatPrint() */

#define	_vnfm_SOURCE
#include "vnfmopca.h"


/* This flag enables the positive correlation matrix constraint
 */
#define	POSITIVE_CORR_CONSTRAINT




#define	__DEBUG__
#undef __DEBUG__

#if defined(_WINDLL) || !defined(TESTLIB)
# undef __DEBUG__
#endif

/*--------------------------------------------------------------
 */

int
VnfmCalibSquareSpotVolEqual(
	VnfmData *that,		/* (I/O) model parameters */
	int idxStart,		/* (I) first timeline idx calibrated */
	int idxEnd,		/* (I) last timeline idx calibrated */
	double *tMat1,		/* (I) array of swaption mat [0..idxEnd] */
	int *freq1,		/* (I) array of swaption freq [0..idxEnd] */
	double *vol1,		/* (I) array of swaption vol [0..idxEnd] */
	double *sqspv)		/* (O) array of square spot vols [0..idxEnd] */
{
static	char		routine[] = "VnfmCalibSquareSpotVolEqual";
	int		status = FAILURE;
	int		i, j,
			nDim = NDIM,
			jIdx, rIdx,
			nExp, n;
	double		x1[VNFM_NDIMMAX],
			v0, v1, K0,
			S, dt,
			lambda;
	double		**jt = NULL;


	double		fwdRate;

	/*
	 *
	 */
	idxStart = MAX(idxStart, 1);
	idxEnd   = MIN(idxEnd, that->fNDates-1);

	if (idxStart > idxEnd) {
	    /* nothing to do */
	    GtoErrMsg("%s: idxStart=%d > idxEnd=%d\n",
		routine, idxStart, idxEnd);
	    goto done;
	}


	VnfmComputeCoeff(that);
#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: start.\n", routine));
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
		"\tidxStart=%d (S=%lf)\tidxEnd=%d (S=%lf)\n",\
		idxStart, TT[idxStart], idxEnd, TT[idxEnd]));
#endif

	/* allocate memory for temporary storage */
	jt = DrlDoubleMatrAlloc(0, nDim*(nDim+1)/2-1, 0, that->fNDates+1);
	if (jt == NULL) goto done;

	/* copy J coeffs in temporary storage */
	VnfmComputeCoeff(that);
	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) {
	    for (n=0; n<=that->fNDates-1; n++) {
		jIdx = JIDX(i, j);
		jt[jIdx][n] = JJ[jIdx][n];
	    }
	}


	/* main calibration loop */
	for (nExp=idxStart; nExp <= idxEnd; nExp++) {
	    /* nExp = timeline exp idx of base vol point to be fitted
	     * n = idx of base vol point to be modified
	     */
	    n = nExp-1;
	    S = that->fTime[nExp];
	    dt = S - that->fTime[n];

	    /* check maturity */
	    if (tMat1[nExp] < 1.92307692e-2) {
		GtoErrMsg("%s: (%s) calibrated rate "
		    "has too low maturity (%lf yrs)\n",
		    routine, GtoFormatDate(that->fDate[nExp]), tMat1[nExp]);
		goto done;
	    }

	    /* compute the Q coefficients */
	    if (freq1[nExp] > 0) {
		VnfmB(that, S, tMat1[nExp], freq1[nExp], &fwdRate, x1);
	    } else {
		VnfmQ(that, S, S+tMat1[nExp], &fwdRate, x1);
	    }

	    /*
             * Convert to the old percentage B/Q coefficients
             */
            for (i=0; i<=nDim-1;i++)
                    x1[i] /= fwdRate;


	    /* compute coefficients: with sigma(S)^2 S = v0, the
	     * equation is
	     *	v0 + K0 * sigma^2 = v1
	     */
	    v0 = 0e0;
	    K0 = 0e0;
	    for (i=0; i<=nDim-1; i++)
	    for (j=i; j<=nDim-1; j++) {
		jIdx = JIDX(i, j);
		rIdx = RHOIDX(i, j);
		lambda = BETA[i] + BETA[j];

		v0 += (i==j ? 1e0 : 2e0) *
			x1[i] * x1[j] *
			exp(-dt*lambda) * jt[jIdx][n];

		K0 += (i==j ? 1e0 : 2e0 * RHO[rIdx][n]) *
			x1[i] * x1[j] * ALPHA[i] * ALPHA[j] *
	    		L(dt, lambda);
	    }
	    /* market volatility */
	    v1 = SQR(vol1[nExp]) * S;

	    /* compute square spot vol */
	    sqspv[n] = (v1 - v0) / K0;


#if !defined(NO_LOGGING) && defined(__DEBUG__)
	    GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
		"[%3d] S=%5.2f dt=%lf tMat=%7.4f freq=%d ",\
		nExp, S, dt, tMat1[nExp], freq1[nExp]));
	    GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
		"\tv0=%lf\tv1=%lf\tK0=%lf\tsqv=%lf\n",\
		v0, v1, K0, sqspv[n]));
#endif


	    /* update coefficients J for next step */
	    n++;
	    dt = that->fTime[n] - that->fTime[n-1];
	    for (i=0; i<=nDim-1; i++)
	    for (j=i; j<=nDim-1; j++) {
		jIdx = JIDX(i, j);
		rIdx = RHOIDX(i, j);
		lambda = BETA[i] + BETA[j];
		jt[jIdx][n] =
		    (i==j ? 1e0 : RHO[rIdx][n-1]) *
		    ALPHA[i] * ALPHA[j] *
		    sqspv[n-1] * L(dt, lambda) +
		    exp(-dt * lambda) * jt[jIdx][n-1];
	    }

	}


	/* made it through OK */
	status = SUCCESS;
done:
	/* tesing of accuracy */
#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING(VnfmFpWrite(that, vnfmFpLog));
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: END\n", routine));
#endif
	/* free memory */
	DrlDoubleMatrFree(jt, 0, nDim*(nDim+1)/2-1, 0, that->fNDates+1);

	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}

/*--------------------------------------------------------------
 * Computes and stores internally the pseudo-volatility integrals
 * <br>
 * \Jt_{ij}(S) &=&integral_0^S
 * 		\alpha_i\sigma_i(t) \alpha_j\sigma_j(t) \rho_{ij}(t)
 *		e^{-(\beta_i+\beta_j) (S-t)} dt.
 * <br>
 * Should not be called directly in applications
 * (use <i> VnfmComputeCoeff</i> instead.
 * Returns 0 if OK.
 */

int
VnfmComputeJSquare(
	VnfmData *that,		/* (I) model parameters */
	double *sqspv,		/* (I) array of square spot vols */
	double **jt)		/* (O) pseudo-J coeffs */
{
	int	i, j,
		nDim = NDIM,
		jIdx, rIdx,
		n;
	double	lambda, dt;

	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) {
	    jIdx = JIDX(i, j);
	    rIdx = RHOIDX(i, j);
	    lambda = BETA[i] + BETA[j];

	    n = 0;
	    jt[jIdx][n] = 0e0;

	    for (n=1; n<=NDATES-1; n++) {
		dt = TT[n] - TT[n-1];
		jt[jIdx][n] =
		    ALPHA[i] * ALPHA[j] *
			sqspv[n-1] *
			(i == j ? 1e0 : RHO[rIdx][n-1]) *
			L(dt, lambda)
		    + exp(-dt*lambda) * jt[jIdx][n-1];
	    }
	}

	return(SUCCESS);
}

/*--------------------------------------------------------------
 */

int
VnfmJSquare(
	VnfmData *that,		/* (I) model parameters */
	double S1,		/* (I) option start */
	double S2,		/* (I) option expiration */
	double T,		/* (I) rate reset */
	double *sqspv,		/* (I) array of square spot vols */
	double **jt,		/* (I) pseudo-J coeffs */
	double *J0)		/* (O) J[i][j] stored as vector (see above) */
{

	int	i, j,
		nDim = NDIM,
		jIdx, rIdx,
		m, n;
	double	lambda;
	/* find closest but lower */
	DrlDoubleArrayFloorIdx(TT, NDATES, S1, &m);
	DrlDoubleArrayFloorIdx(TT, NDATES, S2, &n);

	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) {
	    jIdx = JIDX(i, j);
	    rIdx = RHOIDX(i, j);
	    lambda = BETA[i] + BETA[j];

	    J0[jIdx] = 
		  exp(-(T - TT[n])*lambda) * jt[jIdx][n] 
	    	+ ALPHA[i] * ALPHA[j] *
		    sqspv[n] *
		    (i == j ? 1e0 : RHO[rIdx][n]) *
		    exp(-(T-S2)*lambda) * L(S2-TT[n], lambda)
		- exp(-(T-TT[m])*lambda) * jt[jIdx][m]
	    	- ALPHA[i] * ALPHA[j] *
		    sqspv[m] *
		    (i == j ? 1e0 : RHO[rIdx][m]) *
		    exp(-(T-S1)*lambda)* L(S1-TT[m], lambda);
	}


#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
		"VnfmJ(%3d): S1=%lf S2=%lf T=%lf m=%2d n=%2d\n",\
		__LINE__, S1, S2, T, m, n));
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "VnfmJ     :"));
	GTO_IF_LOGGING(\
	for (i=0; i<=nDim-1; i++)\
	for (j=i; j<=nDim-1; j++)\
		DrlFPrintf(vnfmFpLog, " J[%2d]=%lf ", JIDX(i,j), J0[JIDX(i,j)]);\
	DrlFPrintf(vnfmFpLog, "\n"));
#endif

	return(0);
}


/*--------------------------------------------------------------
 * Calculates the algebraic variance of a rate (can be negative).
 */

int
VnfmAvgQBVolSquare(
	VnfmData *that,
	double S1,		/* (I) start */
	double S2,		/* (I) expiration */
	double tReset,		/* (I) rate reset */
	double tMat,		/* (I) rate maturity */
	int freq,		/* (I) rate frequency */
	double *sqspv,		/* (I) array of square spot vols */
	double **jt,		/* (I) pseudo-J coeffs */
	double *retVal)		/* (O) volatility^2 */
{
static	char	routine[] = "VnfmAvgQBVolSquare";
	double	x[VNFM_NDIMMAX],
		j00[VNFM_NDIMMAX*(VNFM_NDIMMAX+1)/2],
		vol;
	int	i, j, nDim = NDIM, jIdx;
	double	fwdRate;

	/* check */
	if ((S2 > tReset) || (S1 > S2)) {
	    GtoErrMsg("%s: tExp=%lf > tReset=%lf or tStart=%lf > tExp\n",
		routine, S2, tReset, S1);
	    return(FAILURE);
	}

	/* compute the J integrals */
	VnfmJSquare(that, S1, S2, tReset, sqspv, jt, j00);

	/* compute the Q coefficients */
	if (tMat > 0.019230769e0) {
	    if (freq > 0) {
		VnfmB(that, tReset, tMat, freq, &fwdRate, x);
	    } else {
		VnfmQ(that, tReset, tReset+tMat, &fwdRate, x);
	    }
	} else {
	    VnfmQ(that, tReset, tReset+tMat, &fwdRate, x);
	}

	/*
	 * Convert to the old percentage B/Q coefficients
	 */
	for (i=0; i<=nDim-1;i++)
		x[i] /= fwdRate;

	/* vol computation */
	vol = 0e0;
	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) {
	    jIdx = JIDX(i, j);
	    vol += (i == j ? 1e0 : 2e0) * x[i] * x[j] * j00[jIdx];
	}
	vol = SQRT_SIGN(vol / (S2-S1));

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
		"VnfmAvgBVol(%d): S1=%lf S2=%lf T=%lf tMat=%lf "\
		"freq=%d vol=%lf\n",\
		__LINE__, S1, S2, tReset, tMat, freq, vol));
#endif

	*retVal = vol;
	return(SUCCESS);
}

/*--------------------------------------------------------------
 */

int
VnfmAvgQBCorrSquare(
	VnfmData *that,		/* (I) model parameters */
	double S1,		/* (I) vol start */
	double S2,		/* (I) vol expiration */
	double tReset1,		/* (I) reset rate 1 */
	double tMat1,		/* (I) maturity rate 1 */
	int freq1,		/* (I) freq rate 1 */
	double tReset2,		/* (I) reset rate 2 */
	double tMat2,		/* (I) maturity rate 2 */
	int freq2,		/* (I) freq rate 2 */
	double *sqspv,		/* (I) array of square spot vols */
	double **jt,		/* (I) pseudo-J coeffs */
	double *retVal)		/* (O) correlation */
{
static	char	routine[] = "VnfmAvgQBCorr";
	double	x[VNFM_NDIMMAX],
		y[VNFM_NDIMMAX],
		j00[VNFM_NDIMMAX*(VNFM_NDIMMAX+1)/2],
		vol1, vol2, cor;
	int	i, j, nDim = NDIM, jIdx;

	double	fwdRate1, fwdRate2;

	/* check */
	if ((S2 > tReset1) || (S2 > tReset2)) {
	    GtoErrMsg("%s: volExp (%lf) > tReset1 (%lf) or tReset2 (%lf).\n",
		routine, S2, tReset1, tReset2);
	    return(FAILURE);
	}
	if (S1 > S2) {
	    GtoErrMsg("%s: volStart (%lf) > volExp (%lf).\n", routine, S1, S2);
	    return(FAILURE);
	}

	/* compute vol coefficients */
	if (VnfmQBCoeff(that, S2, tReset1, tMat1, freq1, &fwdRate1, x) != SUCCESS)
		return(FAILURE);
	if (VnfmQBCoeff(that, S2, tReset2, tMat2, freq2, &fwdRate2, y) != SUCCESS)
		return(FAILURE);

	/* compute the J integrals */
	VnfmJSquare(that, S1, S2, S2, sqspv, jt, j00);

	/* vol computation */
	vol1 = vol2 = cor = 0e0;
	for (i=0; i<=nDim-1; i++)
	for (j=0; j<=nDim-1; j++) {
	    jIdx = (i <= j ? JIDX(i, j) : JIDX(j, i));
	    vol1 += x[i] * x[j] * j00[jIdx];
	    vol2 += y[i] * y[j] * j00[jIdx];
	    cor  += x[i] * y[j] * j00[jIdx];
	}
	cor = cor / (SQRT_SIGN(vol1) * SQRT_SIGN(vol2));

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: "\
	    "S1=%lf S2=%lf tReset1=%lf tReset2=%lf tMat1=%lf tMat2=%lf\n",\
	    routine, S1, S2, tReset1, tReset2, tMat1, tMat2));
	GTO_IF_LOGGING(\
	for (i=0; i<=nDim-1; i++)\
	    DrlFPrintf(vnfmFpLog, "\t[%2d]\tx=%lf\ty=%lf\n", i, x[i], y[i]));
#endif


	*retVal = cor;

	return(SUCCESS);
}



/*--------------------------------------------------------------
 * Computes the value of a volatility benchmark "pt".
 */

int
VnfmTVolBenchmarkValueSquare(
	VnfmData *that,		/* (I) model parameters */
	TVolBenchmark *pt,	/* (I) benchmark */
	double *sqspv,		/* (I) array of square spot vols */
	double **jt,		/* (I) pseudo-J coeffs */
	double *retVal)		/* (O) value */
{
static  char    routine[] = "VnfmTVolBenchmarkValueSquare";
	int     status = FAILURE;

	TDate		expDate;
        double		freq1, freq2;
	double		U1=0e0, U2, tReset1, tMat1, tReset2, tMat2;

	TFloatRateVolPoint      *vpt;
	TFloatRateCorrPoint     *cpt;

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
		"%s: %s\n", routine, DrlTVolBenchmarkPrint(NULL, pt)));
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


	    if (VnfmAvgQBVolSquare(
	    	that,
		U1, U2,
		tReset1, tMat1, (int) freq1,
		sqspv, jt,
		retVal) != SUCCESS)
			goto done;

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	   GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
		"%s: U1=%8.4f U2=%8.4f tReset=%8.4f "\
		"tMat=%8.4f freq=%d vol=%12.8f\n",\
		routine,\
		U1, U2, tReset1, tMat1, (int) freq1, *retVal));
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

	    if (VnfmAvgQBCorrSquare(that,
		U1, U2,
		tReset1, tMat1, (int) freq1,
		tReset2, tMat2, (int) freq2,
		sqspv, jt,
		retVal) != SUCCESS)
			goto done;

#if !defined(NO_LOGGING) && defined(__DEBUG__)
	    GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: U1=%8.4f U2=%8.4f "\
		"tReset1=%lf, tMat1=%lf, freq1=%d, "\
		"tReset2=%lf, tMat2=%lf, freq2=%d, "\
		"retVal=%lf\n",\
		routine, U1, U2,\
		tReset1, tMat1, (int) freq1,\
		tReset2, tMat2, (int) freq2,\
		*retVal));
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



/*--------------------------------------------------------------
 * Static variables used in optimisation
 */



static	int		nIterF,
			nIterDF,
			nIterC;

/*
 * User data for the optimisation
 */
typedef	struct	{
	VnfmData	*tfData;	/* model parameters */
	double		**jt;		/* pseudo-J coeffs */
	double		*sqspv;		/* array of square spot vol */

	int		nOptBench;
	TVolBenchmark	*optBench;
	double		*optMid;
	double		*optWeight;

	int		nConstBench;
	TVolBenchmark	*constBench;
	double		*constMin;
	double		*constMax;

	double		*crMat;		/* array of swaption mat [0..idxEnd] */
	int		*crFreq;	/* array of swaption freq [0..idxEnd] */
	double		*crVol;		/* array of swaption vol [0..idxEnd] */

} TUserData;


static	TUserData	*userData;	/* static data passed to optimization
					 * and constraint functions */


static	int		objectiveFunction(int *MODE, int *N, double *X,
				void *optData,
				double *OBJF, double *OBJGRD);
static	int		constraintFunction(int *NCNLN, int *N, int *MODE,
				int *NEEDC, double *X,
				void *optData,
				double *C, double *cjac1);
static	int		updateUserData(double *y, TUserData *usrData);


#define	NCMAX		10		/* max # of correlation constraints */
#define	NVMAX		30		/* max # of volatility points */


#undef	__DEBUG__
#define	__DEBUG__
#ifdef 	__DEBUG__
static	double	*blcnln1, *bucnln1;
#endif


/*f-------------------------------------------------------------
 * Parameter optimization using bootstrapping with constraints.
 *                                                             
 * <br><br>
 * Performs the calibration of the two-factor parameters over
 * a series of short-term volatility benchmarks.
 * <br>
 * <br>[that] On entry, parameter data structure with an already
 * initialized timeline. On exit, contains the optimal parameters.
 * <br>[nOptBench] number of optimized benchmarks,
 * <br>[optBench] array of optimized benchmarks,
 * <br>[optMid] market value of benchmarks,
 * <br>[optMod] On exit, model value of benchmarks,
 * <br>[nConstBench] number of constrained benchmarks,
 * <br>[constBench] array of constrained benchmarks,
 * <br>[constMin] array of low  constraints,
 * <br>[constMax] array of high constraints,
 * <br>[constMod] On exit, array of model values of constraints,
 * <br>[paOptFlag] array of flags (TRUE=optimize),
 * <br>[paOpt] array of optimal paramters (O),
 * <br>[paMin] array of low  constraints on parameters,
 * <br>[paMax] array high constraints on parameters,
 * <br>[paMid] array of initial guess (or value if no optimization),
 * <br>
 * Returns 0 iff successful.
 */

int
VnfmCalibParamSquareVol(
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
	double *paMid)		/* (I) initial guess (or value if no opt)*/
{
static	char	routine[] = "VnfmCalibParamSquareVol";
	int	status = FAILURE;
	int	NVAR ;			/* number of variables */
	double	*blsc=NULL,
		*busc=NULL;		/* constraints on variables */

	int	NCLIN;			/* # linear constraints */
	double	*blclin=NULL,
		*buclin=NULL,
		**a=NULL;

	int	NCNLN;
	double	*blcnln=NULL,
		*bucnln=NULL;
	int	ITERMAX;
	double	*C=NULL;
	double	OBJF;
	double	*X=NULL;

	int	nDim = NDIM,
		nDates = that->fNDates,
		i, k, idx,
		idxF1, idxF2, idxR, idxT;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s: start.\n", routine));
#endif
	/*
	 * Copy pointers in static variables (visible outside the routine)
	 */
	userData = NULL;
	if ((userData = NEW(TUserData)) == NULL) goto done;

	userData->tfData = that;
	userData->nOptBench = nOptBench;
	userData->optBench = optBench;
	userData->optMid = optMid;
	userData->optWeight = optWeight;
	userData->nConstBench = nConstBench;
	userData->constBench = constBench;
	userData->constMax = constMax;
	userData->constMax = constMax;


	/* Number of parameters passed:
	 * 1: mr and others
	 * mr:		(b1,...,bn)			n
	 * weights:	(a1,...,bn)			n
	 * corr:	(r11,r12,...,rn-1n)		n(n-1)/2
	 * Total:				nPa=	2n + n(n-1)/2
	 * Parameters to optimize:
	 * a1 is NOT optimized (rescaling)
	 * Total:				nPa-1
	 */
	if (nPa != 2*nDim + nDim*(nDim-1)/2) {
	    GtoErrMsg("%s: received optimisation arrays of length %d,"
		" expected %d (nDim=%d).\n", routine,
		nPa, 2*nDim + nDim*(nDim-1)/2, nDim);
	    goto done;
	}


	/* alloc memory */
	if ((userData->sqspv = NEW_ARRAY(double, that->fNDates)) == NULL)
		goto done;
	if ((userData->crMat = NEW_ARRAY(double, that->fNDates)) == NULL)
		goto done;
	if ((userData->crFreq = NEW_ARRAY(int, that->fNDates)) == NULL)
		goto done;
	if ((userData->crVol = NEW_ARRAY(double, that->fNDates)) == NULL)
		goto done;
	if ((userData->jt = DrlDoubleMatrAlloc(0, nDim*(nDim+1)/2-1, 0,
		that->fNDates)) == NULL) goto done;


	/*
	 * Store vol bootstrap data
	 */
	for (i=0; i<=that->fNDates-1; i++) {
		userData->crMat[i]  = cRateMat[i];
                userData->crFreq[i] = cRateFreq[i];
		userData->crVol[i]  = cRateVol[i];
                userData->sqspv[i] = 0e0;
        }


	/* fill data structure with input (mid) parameters */
	for (idxF1=0; idxF1<=nDim-1; idxF1++) {
		that->fBeta[idxF1]  = paMid[idxF1];
	}
	for (idxF1=0; idxF1<=nDim-1; idxF1++) {
		that->fAlpha[idxF1] = paMid[idxF1+nDim];
	}

	for (idxF1=0;       idxF1<=nDim-1; idxF1++) {
		for (idxT=0; idxT<=that->fNDates-1; idxT++) {
			that->fSigma[idxF1][idxT] = 1e0;
		}
	}

	for (idxF1=0;       idxF1<=nDim-1; idxF1++)
	for (idxF2=idxF1+1; idxF2<=nDim-1; idxF2++) {
		idxR = RHOIDX(idxF1, idxF2);
		for (idxT=0; idxT<=that->fNDates-1; idxT++) {
			that->fRho[idxR][idxT] = paMid[idxR+nDim+nDim];
		}
	}



	/* compute forward rates */
	if (VnfmComputeCoeff(that) != 0) goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: opt benchmarks:\n", routine);\
	DrlTVolBenchmarkFpWrite(optBench, nOptBench, vnfmFpLog, \
			optMid, NULL);\
	DrlFPrintf(vnfmFpLog, "%s: const benchmarks:\n", routine);\
	DrlTVolBenchmarkFpWrite(constBench, nConstBench, vnfmFpLog,\
			constMin, constMax, NULL));
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: initial model data:\n", routine);\
	VnfmFpWrite(that, vnfmFpLog);\
	DrlFPrintf(vnfmFpLog, "%s: zero curve:\n", routine);\
	DrlTCurveFpWrite(that->fZcCurve, vnfmFpLog, DRL_TCURVE_FMT_PERCENT));
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: calibrated volatilities:\n", routine);\
	DrlFPrintf(vnfmFpLog, "\tidx\tmat\tfreq\tvol\n");\
	for (i=0; i<=that->fNDates-1; i++) {\
	    DrlFPrintf(vnfmFpLog, "\t%d\t%lf\t%d\t%lf\n",\
		i, userData->crMat[i], userData->crFreq[i],\
		userData->crVol[i]);\
        });

#endif


	/*-----------------------------------------------------
	 * COUNT NUMBER OF VARIABLES TO OPTIMIZE
	 *-----------------------------------------------------*/
	NVAR = nPa-1;

	/*-----------------------------------------------------
	 * BOUNDS ON VARIABLES AND INITIAL VALUE
	 *-----------------------------------------------------*/

	/* allocate memory for bounds */
	if ((X    = DrlDoubleVectAlloc(0, NVAR-1)) == NULL) goto done;
	if ((blsc = DrlDoubleVectAlloc(0, NVAR-1)) == NULL) goto done;
	if ((busc = DrlDoubleVectAlloc(0, NVAR-1)) == NULL) goto done;

	/* fill variable bounds constraints */
	/* i-th optim variable -> k-th model variable */
	i = 0;
	k = 0;
	for (idxF1=0; idxF1<=nDim-1; idxF1++) {
	    blsc[i] = (paOptFlag[k] != 0L ? paMin[k] : paMid[k]);
	    busc[i] = (paOptFlag[k] != 0L ? paMax[k] : paMid[k]);
	    X[i]    = paMid[k];
	    i++; k++;
	}
	k++;	/* skip alpha 1 */
	for (idxF1=1; idxF1<=nDim-1; idxF1++) {
	    blsc[i] = (paOptFlag[k] != 0L ? paMin[k] : paMid[k]);
	    busc[i] = (paOptFlag[k] != 0L ? paMax[k] : paMid[k]);
	    X[i]    = paMid[k];
	    i++; k++;
	}


	for (idxF1=0;       idxF1<=nDim-1; idxF1++)
	for (idxF2=idxF1+1; idxF2<=nDim-1; idxF2++) {
	    blsc[i] = (paOptFlag[k] != 0L ? paMin[k] : paMid[k]);
	    busc[i] = (paOptFlag[k] != 0L ? paMax[k] : paMid[k]);
	    X[i]    = paMid[k];
	    i++; k++;
	}

	if ((i != NVAR) || (k != nPa)) {
	    PROGRAM_BUG(); goto done;
	}


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "\tNVAR=%d\n", NVAR);\
	for (i=0; i<=NVAR-1; i++)\
	    DrlFPrintf(vnfmFpLog, "\t\t[%2d] blsc=%lf\tbusc=%lf\n",\
			i, blsc[i], busc[i]);\
	);
#endif



	/*-----------------------------------------------------
	 * LINEAR CONSTRAINTS
	 *-----------------------------------------------------*/
	NCLIN = 0;
#ifdef	_SKIP
	/* Linear constraints: beta0 >= beta1 >= beta2 ..
	 * beta:					n-1
	 * Total:					n-1
	 */
	NCLIN = nDim-1;
	if ((blclin = DrlDoubleVectAlloc(0, NCLIN-1)) == NULL)
		goto done;
	if ((buclin = DrlDoubleVectAlloc(0, NCLIN-1)) == NULL)
		goto done;
	if ((a = DrlDoubleMatrAlloc(0, NCLIN-1, 0, NVAR)) == NULL)
		goto done;

	/* set the linear constraints */
	for (i=0; i<=nDim-2; i++) {
		/* beta[i] > beta[i+1] */
		blclin[i] = 0e0;
		buclin[i] = 1e6;
		for (k=0; k<=NVAR-1; k++) a[i][k] = 0e0;
		a[i][i]   =  1e0 ;
		a[i][i+1] = -1e0 ;
	}
#endif


	/*-----------------------------------------------------
	 * NONLINEAR CONSTRAINTS
	 *-----------------------------------------------------*/


	/* Nonlinear Constraints:
	 * sqspv >= sqspvMin			nDates-1
	 * additional bench			nConstBench
	 * positivity of corr matrix 		nDim-2
	 */
#ifdef	POSITIVE_CORR_CONSTRAINT
	NCNLN = (nDates-1) + nConstBench + MAX(nDim-2, 0);
#else
	NCNLN = (nDates-1) + nConstBench;
#endif
	if (NCNLN >= 1) {
		if ((blcnln = DrlDoubleVectAlloc(0, NCNLN-1)) == NULL)
			goto done;
		if ((bucnln = DrlDoubleVectAlloc(0, NCNLN-1)) == NULL)
			goto done;
		if ((C =      DrlDoubleVectAlloc(0, NCNLN-1)) == NULL)
			goto done;

		for(k=0; k<=nDates-2; k++) {
		    blcnln[k] = sqspvMin;
		    bucnln[k] = 1e32;
		}
		for(k=nDates-1; k<=(nDates-1)+nConstBench-1; k++) {
		    idx = k - (userData->tfData->fNDates-1);
		    blcnln[k] = constMin[idx];
		    bucnln[k] = constMax[idx];
		}
#ifdef	POSITIVE_CORR_CONSTRAINT
		for(k=(nDates-1)+nConstBench;
		    k<(nDates-1)+nConstBench+MAX(nDim-2,0); k++) {
		    blcnln[k] = 1e-1;
		    bucnln[k] = 1e12;
		}
#endif
	}

	/*
	 * static variables used in func1
	 */
	ITERMAX = 400;
	OBJF = 0e0;

	nIterF = nIterDF = nIterC = 0 ;
#ifdef 	__DEBUG__
	blcnln1 = blcnln;
	bucnln1 = bucnln;
#endif

	/*-----------------------------------------------------
	 * Do The Minimization
	 *-----------------------------------------------------*/
	if (DrlNLinProg(
		NVAR, blsc, busc,
		NCLIN, blclin, buclin, a,
		NCNLN, blcnln, bucnln, constraintFunction,
		objectiveFunction, ITERMAX, C, &OBJF, X,
		(void*)NULL, /* user data */
		DRL_LNO_IMSL_NLN,	/* method */
		0)
	    != 0) {
		GtoErrMsg("%s: optimization failed.\n", routine);
		goto done;
	}


	/*
	 * Store the optimal parameters
	 */
	i = 0;
	for (idxF1=0; idxF1<=nDim-1; idxF1++) {
		that->fBeta[idxF1]  = X[i];
		i++;
	}
	for (idxF1=1; idxF1<=nDim-1; idxF1++) {
		that->fAlpha[idxF1] = X[i];
		i++;
	}
	for (idxF1=0;       idxF1<=nDim-1; idxF1++)
	for (idxF2=idxF1+1; idxF2<=nDim-1; idxF2++) {
		idxR = RHOIDX(idxF1, idxF2);
		for (idxT=0; idxT<=that->fNDates-1; idxT++) {
			that->fRho[idxR][idxT] = X[i];
		}
		i++;
	}



	/* Rescale alpha on the avregage vol of all benchmarks */
	IF_FAILED_DONE( VnfmComputeCoeff(that));
	IF_FAILED_DONE( VnfmTVolBenchmarkRescaleAlpha(
		that,
		nOptBench,
		optWeight,
		optMid,
		optBench));
#ifdef	_SKIP
#endif



	/* Return parameters */
	i = 0;
	for (idxF1=0; idxF1<=nDim-1; idxF1++) {
		paOpt[i] = that->fBeta[idxF1];
		i++;
	}
	for (idxF1=0; idxF1<=nDim-1; idxF1++) {
		paOpt[i] = that->fAlpha[idxF1];
		i++;
	}
	for (idxF1=0;       idxF1<=nDim-1; idxF1++)
	for (idxF2=idxF1+1; idxF2<=nDim-1; idxF2++) {
		idxR = RHOIDX(idxF1, idxF2);
		paOpt[i] = that->fRho[idxR][0];
		i++;
	}

	
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: optimization successful.\n", routine));
#endif

	/*-----------------------------------------------------
	 * Optimization successful: finish
	 *-----------------------------------------------------*/


	if (VnfmCalibSquareSpotVolEqual(
		userData->tfData,
		0,
		userData->tfData->fNDates-1,
		userData->crMat,
		userData->crFreq,
		userData->crVol,
		userData->sqspv) != SUCCESS)
			goto done;

	/* compute pseudo-J coefficients */
	if (VnfmComputeJSquare(
		userData->tfData,
		userData->sqspv,
		userData->jt) != SUCCESS)
			goto done;

	/* Compute benchmark values */
	for(i=0; i<=nOptBench-1;i++) {
	    if (VnfmTVolBenchmarkValueSquare(
		userData->tfData,
		&optBench[i],
		userData->sqspv, userData->jt,
		&optMod[i]) != SUCCESS)
			goto done;
	}
	for(i=0; i<=nConstBench-1;i++) {
	    if (VnfmTVolBenchmarkValueSquare(
		userData->tfData,
		&constBench[i],
		userData->sqspv, userData->jt,
		&constMod[i]) != SUCCESS)
			goto done;
	}



	/* recalibrate with real volatilities */
	if (VnfmVolCalib1VArbitrary(
		userData->tfData,
		0,
		userData->tfData->fNDates-1,
		userData->crMat,
		userData->crFreq,
		userData->crVol,
		LOGVOL,
		TRUE,
		NULL) != SUCCESS)
			goto done;

	if (VnfmComputeCoeff(that) != SUCCESS)
		goto done;

	/* Compute benchmark values */
	for(i=0; i<=nOptBench-1;i++) {
	    if (VnfmTVolBenchmarkValue(that, &optBench[i],
	    	&optMod[i]) != SUCCESS) goto done;
	}
	for(i=0; i<=nConstBench-1;i++) {
	    if (VnfmTVolBenchmarkValue(that, &constBench[i],
	    	&constMod[i]) != SUCCESS) goto done;
	}


	/*
	 * exit the routine
	 */

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: done.\n", routine);\
	DrlFPrintf(vnfmFpLog, "nIterF=%4d nIterDF=%4d\n", ++nIterF, nIterDF);\
	VnfmFpWrite(that, vnfmFpLog));
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "benchmark volatilities\n\tmkt\tmod\tdiff\n");\
	for(i=0; i<=nOptBench-1;i++) {\
		DrlFPrintf(vnfmFpLog, "\t%.2f\t%.2f\t%.2f\n",\
		optMid[i]*1e2, optMod[i]*1e2, (optMid[i]-optMod[i])*1e2);\
	});
#endif


	/* made it through OK */
	status = SUCCESS;
done:
	/* free alloctaed memory */
	FREE(userData->crMat);
	FREE(userData->crFreq);
	FREE(userData->crVol);
	FREE(userData->sqspv);
	DrlDoubleMatrFree(userData->jt, 0, nDim*(nDim+1)/2-1, 0, that->fNDates);
	FREE(userData);

	DrlDoubleVectFree(X,    0, NVAR-1);
	DrlDoubleVectFree(blsc, 0, NVAR-1);
	DrlDoubleVectFree(busc, 0, NVAR-1);

	DrlDoubleVectFree(blclin, 0, NCLIN-1);
	DrlDoubleVectFree(buclin, 0, NCLIN-1);
	DrlDoubleMatrFree(a, 0, NCLIN-1, 0, NVAR);

	DrlDoubleVectFree(blcnln, 0, NCNLN-1);
	DrlDoubleVectFree(bucnln, 0, NCNLN-1);
	DrlDoubleVectFree(C,      0, NCNLN-1);

	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*--------------------------------------------------------------
 * y
 */

static	int
updateUserData(double *y, TUserData *usrData)
{
int	status = FAILURE;
	register int	i, j, n;
	VnfmData	*m = usrData->tfData;
	int		k,
			rIdx,
			nDates = m->fNDates,
			nDim = m->fNf;
	double		*rvect;


	k = 0;

	/* update parameters */
	for (i=0; i<=nDim-1; i++)
		m->fBeta[i]  = y[k++];
	m->fAlpha[0] = 1e0;
	for (i=1; i<=nDim-1; i++)
		m->fAlpha[i] = y[k++];

	for (i=0;   i<=nDim-1; i++) 
	for (j=i+1; j<=nDim-1; j++) {
		rIdx = RHOIDX(i,j);
		rvect = &m->fRho[rIdx][0];
		rvect[0] = y[k++];
		for (n=0; n<=nDates-1; n++) {
			rvect[n] = rvect[0];
		}
	}

	/* calibrate square spot vols */
	if (VnfmCalibSquareSpotVolEqual(
		usrData->tfData,
		0,
		usrData->tfData->fNDates-1,
		usrData->crMat,
		usrData->crFreq,
		usrData->crVol,
		usrData->sqspv) != SUCCESS)
			goto done;
#ifdef	_SKIP
	VnfmCalib1VArbitrary(
		usrData->tfData,
		0,
		usrData->tfData->fNDates-1,
		usrData->tMat,
		usrData->freq,
		usrData->vol,
		LOGVOL,
		TRUE,
		NULL);
	DrlFPrintf(vnfmFpLog, "userUpdate:\n");
	for (n=0; n<=nDates-1; n++) {
	    DrlFPrintf(vnfmFpLog, "%lf\t%lf\n",
		usrData->tfData->fSigma[0][n]
		* usrData->tfData->fSigma[0][n],
		usrData->sqspv[n]);
	}
#endif


	/* compute pseudo-J coefficients */
	if (VnfmComputeJSquare(
		usrData->tfData,
		usrData->sqspv,
		usrData->jt) != SUCCESS)
			goto done;


	status = SUCCESS;
done:
	return(status);
}


/*ARGSUSED*/
/*--------------------------------------------------------------
 * Function:	Objective function for the minimization routine
 */

static	int
objectiveFunction(
	int *MODE, int *N, double *X,
	void *optData,
	double *OBJF, double *OBJGRD)
{
/*static	char	routine[] = "objectiveFunction";*/
int	status = FAILURE;
	int	i ;
	double	v0, v1, s0;

	/*
	 * 1: new values of the model parameters
	 */
	if (updateUserData(X, userData) != SUCCESS)
		goto done;

	/*
	 * 2: compute the objective function
	 */
	s0 = 0e0 ;
	for (i=0; i<=userData->nOptBench-1; i++) {
	    if (!IS_ALMOST_ZERO(userData->optWeight[i])) {
		/* model volatility estimate */
		if (VnfmTVolBenchmarkValueSquare(
			userData->tfData,
			&userData->optBench[i],
			userData->sqspv, userData->jt,
			&v0) != SUCCESS) {
				goto done;
		}

		/* actual observed volatility */
		v1 = userData->optMid[i];

		/* diff squared */
		s0 += (v0 - v1)*(v0 - v1)*userData->optWeight[i];
	    }
	}

	++nIterF;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "XOPT:");\
	for(i=0; i<=*N-1; i++)\
		DrlFPrintf(vnfmFpLog, " %s", DrlFloatPrint(NULL, X[i], 8));\
	DrlFPrintf(vnfmFpLog, " o=%.4e\n", s0));
#endif

	*OBJF = s0;

	status = SUCCESS;
done:
	return(status);
}


/*ARGSUSED*/
/*--------------------------------------------------------------
 * Function:	Nonlinear constraints
 */

static	int
constraintFunction(
	int *NCNLN, int *N, int *MODE, int *NEEDC,
	double *X,
	void *optData,
	double *C, double *cjac1)
{
/*static	char	routine[] = "constraintFunction";*/
	int	status = FAILURE;
	int	k,		/* constraint index */
		k0,		/* offset for 1st set of constraints */
		k1,		/* offset for 2st set of constraints */
		idx,
		i;


	k0 = userData->tfData->fNDates-1;
	k1 = userData->tfData->fNDates-1 + userData->nConstBench;

	/*
	 * 1: new values of the model parameters
	 */
	if (updateUserData(X, userData) != SUCCESS)
		goto done;


	for(k=0; k<=*NCNLN-1; k++) {
	    /* constraint need to be computed */
	    if (NEEDC[k] > 0) {
		if (k <= k0-1) {
		    /* sq spot vol > 0 */
		    idx = k;
		    C[k] = userData->sqspv[k];

		} else if (k <= k1-1) {

		    idx = k - k0;


		    if (VnfmTVolBenchmarkValueSquare(
			userData->tfData,
			&userData->constBench[idx],
			userData->sqspv, userData->jt,
			&C[k]) != SUCCESS) {
			C[k] = sqrt(-1e0);
		    }
		} else {
			/* Corr mat positivity constraint
			 */
			idx = k - k1;

			if (VnfmCorrMinorDeterm(
				userData->tfData, 0,
				idx+2, 			/* order */
				&C[k]) != SUCCESS) 	/* determinant */
			    C[k] = sqrt(-1e0);
		}

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "XOPT:");\
	for(i=0; i<=*N-1; i++)\
		DrlFPrintf(vnfmFpLog, " %s", DrlFloatPrint(NULL, X[i], 8));\
	DrlFPrintf(vnfmFpLog, " C[%2d]=%8.4f %c\n",\
		k, C[k], (C[k] < blcnln1[k] - 1e-2 ? '-' :\
			(C[k] >= bucnln1[k] + 1e-2 ? '+' : 's'))));
#endif

	    }
	}

	status = SUCCESS;
done:
	return(status);
}




