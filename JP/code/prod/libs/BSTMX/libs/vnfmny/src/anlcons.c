/****************************************************************
 * Module:	VNFM
 * Submodule:	ANLY
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include "drlio.h"			/* platform compatibility */
#include <math.h>

#include "convert.h"

#define	_vnfm_SOURCE
#include "vnfmanly.h"


static	int	VnfmMemoryAlloc(VnfmData *that, int nDates,
				int nf, int nZDates);
static	int	VnfmMemoryFree(VnfmData *that);


/*f-------------------------------------------------------------
 * Allocate a new VnfmData structure.
 *                                                             
 * <br><br>
 * Returns a pointer to a new two-factor parameter data
 * structure.
 * Allocate memory for "nDim" dimensions and a timeline of size "nDates"
 * in the data structure "that" (for 
 * storage of parameters and internal coefficients).
 * Returns NULL if failed.
 */

VnfmData*
VnfmNew(int nDates, int nDim, int nZDates)
{
static	char		routine[] = "VnfmNew";
	VnfmData*	that = (VnfmData*)NULL;
	int		status = FAILURE;

	if ((that = NEW(VnfmData)) == NULL) {
		GtoErrMsg("%s: memory allocation failed\n", routine);
		goto done;
	}

	if (VnfmMemoryAlloc(that, nDates, nDim, nZDates) != 0) goto done;

	status = SUCCESS;
done:
	if (status != SUCCESS) {
		if (that != NULL) {
			FREE(that);
			that = NULL;
		}
		GtoErrMsg("%s: failed\n", routine);
	}
	return(that);
}


/*f-------------------------------------------------------------
 * Free a VnfmData structure.
 *                                                             
 * <br><br>
 * Frees a pointer to a two-factor parameter data structure
 * created by <i> VnfmNew</i>.
 */

int
VnfmFree(VnfmData *that)
{
	VnfmMemoryFree(that);
	if (that != NULL) FREE((char*) that);
	return(0);
}


/*f-------------------------------------------------------------
 * Allocate a new VnfmData structure with given input timeline.
 *                                                             
 * <br><br>
 * Creates and initializes the timeline of an N-factor parameter set.
 * The timeline is build according to the following rule:
 * <br>
 * <br> If "dates" is not NULL and "nDates" is positive,
 *	"refDate" and the array dates are used to build the timeline.
 * <br> If "dates" is NULL (or "nDates" is 0), and "volCurve" is not null,
 *	the dates of "volCurve" are used for the timeline ("refDate" is not
 *	used).
 * <br> If "dates" is NULL (or "nDates" is 0), and "volCurve" is also NULL,
 *	a default timeline of 360 months is used (starting from "refDate").
 * <br>
 * Returns NULL if failed.
 */

VnfmData*
VnfmNewTimeLine(
	int nDim,		/* (I) # of factors */
 	double backboneq,	/* (I) backbone WARNING: 0.5=normal */
	TDate refDate,		/* (I) reference date (today) */
	int nDates,		/* (I) # dates (can be -1) */
	TDate *dates,		/* (I) array of dates (can be NULL) */
	TCurve *volCurve,	/* (I) used for timeline only (can be NULL) */
	TCurve *zcCurve)	/* (I) zero curve */
{
static	char	routine[] = "VnfmNewTimeLine";
	int	status = FAILURE;
	int	i, k,
		nf = nDim,			/* # of factors */
		nZeroDates = zcCurve->fNumItems+1;
	VnfmData *that = NULL;

	/* */
	if ((nDates > 0) && (dates != (TDate*)NULL))  {
	    /* use input dates */
	    if (refDate == dates[0]) {
		if ((that = VnfmNew(nDates, nf, nZeroDates)) == NULL)
			goto done;

		/* set dates */
		for (i=0; i<=that->fNDates-1; i++) {
			that->fDate[i] = dates[i];
		}
	    } else {
		if ((that = VnfmNew(nDates+1, nf, nZeroDates)) == NULL)
			goto done;

		/* set dates */
		that->fDate[0] = refDate;
		for (i=1; i<=that->fNDates-1; i++) {
			that->fDate[i] = dates[i-1];
		}
	    }
	} else if (volCurve != NULL) {
	    /* use vol curve dates */
	    if (volCurve->fBaseDate == volCurve->fArray[0].fDate) {
		nDates = volCurve->fNumItems;
		if ((that = VnfmNew(nDates, nf, nZeroDates)) == NULL)
			goto done;

		/* set dates */
		for (i=0; i<=that->fNDates-1; i++) {
			that->fDate[i] = volCurve->fArray[i].fDate;
		}
	    } else {
		nDates = volCurve->fNumItems;
		if ((that = VnfmNew(nDates+1, nf, nZeroDates)) == NULL)
			goto done;

		/* set dates */
		that->fDate[0] = volCurve->fBaseDate;
		for (i=1; i<=that->fNDates-1; i++) {
			that->fDate[i] = volCurve->fArray[i-1].fDate;
		}
	    }
	} else {
	    /* used default timeline */
	    TDateInterval	interval;

	    nDates = 360;
	    if ((that = VnfmNew(nDates, nf, nZeroDates)) == NULL)
		goto done;
	    for (i=0; i<=that->fNDates-1; i++) {
	        if (GtoMakeDateInterval(i, 'M', &interval) != SUCCESS)
		    goto done;
		GtoDtFwdAny(refDate, &interval, &that->fDate[i]);
	    }
	}

	/*
	 * Set parameters
	 */
	/* Set distribution type */
	that->fBackBoneQ = 1e0 -2*backboneq;

	for (k=0; k<=nf-1; k++) {
	    that->fBeta[k] = 0e0;
	    that->fAlpha[k] = 0e0;
	}


	/* set time */
	for (i=0; i<=that->fNDates-1; i++) {
	    GtoDayCountFraction(REFDATE, that->fDate[i],
		GTO_ACT_365F, &that->fTime[i] );
	}
	for (i=0; i<=that->fNDates-2; i++) {
	    that->fDt[i] = that->fTime[i+1] - that->fTime[i];
	}
	that->fDt[that->fNDates-1] = 0.;


	/* Set zero curve */
	that->fZcCurve = zcCurve;

#ifndef	VNFM_V5X
	/* Set fZTime */
	that->fZTime[0] = 0.0;    /* equals REFDATE */
	for (i=1;i<=NZDATES-1;i++) {
		GtoDayCountFraction(REFDATE, zcCurve->fArray[i-1].fDate,
			GTO_ACT_365F, &that->fZTime[i]);
	}
#endif	/*VNFM_V5X*/

#ifdef	_SKIP
	/* check valid data */
	if (VnfmCheckValid(that) != 0) {
	    goto done;
	}
#endif

	/*
	 *
	 */
	status = SUCCESS;
done:
	if (status != 0) {
		GtoErrMsg("%s: failed\n", routine);
		VnfmFree(that);
		return(NULL);
	}

	return(that);
}




/*f-------------------------------------------------------------
 * Allocate a new VnfmData structure with 2F parameters.
 *                                                             
 * <br><br>
 * Creates and initializes
 * the timeline of the two-factor parameter set with given
 * arrays of spot volatilities and correlation.
 * Returns NULL if failed.
 */

VnfmData*
VnfmNew2SpotVols(
	TDate refDate,		/* (I) vol reference date */
 	double backboneq,	/* (I) backbone WARNING: 0.5=normal */
	double beta1,		/* (I) 2F parameter */
	double beta2,		/* (I) 2F parameter */
	double alpha1,		/* (I) 2F parameter */
	double alpha2,		/* (I) 2F parameter */
	int nDates,		/* (I) # dates */
	TDate *dates,		/* (I) array of dates [0..nDates-1] */
	double *sigma1,		/* (I) 2F parameter [0..nDates-1] */
	double *sigma2,		/* (I) 2F parameter [0..nDates-1] */
	double *rho,		/* (I) 2F parameter [0..nDates-1] */
	TCurve *zcCurve)	/* (I) zero curve */
{
static	char	routine[] = "VnfmNew2SpotVols";
	int	status = FAILURE;
	int	i,
		nf = 2;			/* # of factors */
	VnfmData *that = NULL;


	/* use input dates */
	if (refDate == dates[0]) {
		if ((that = VnfmNew(nDates, nf, zcCurve->fNumItems+1)) == NULL)
			goto done;

		/* set dates */
		for (i=0; i<=that->fNDates-1; i++) {
			that->fDate[i] = dates[i];
	    		that->fSigma[0][i] = sigma1[i];
	    		that->fSigma[1][i] = sigma2[i];
	    		that->fRho[0][i] = rho[i];
		}
	} else {
		if ((that = VnfmNew(nDates+1, nf,zcCurve->fNumItems+1)) == NULL)
			goto done;

		/* set dates */
		that->fDate[0] = refDate;
	    	that->fSigma[0][0] = sigma1[0];
	    	that->fSigma[1][0] = sigma2[0];
	    	that->fRho[0][0] = rho[0];
		for (i=1; i<=that->fNDates-1; i++) {
			that->fDate[i] = dates[i-1];
	    		that->fSigma[0][i] = sigma1[i];
	    		that->fSigma[1][i] = sigma2[i];
	    		that->fRho[0][i] = rho[i];
		}
	}

	/* Set distribution type */
	that->fBackBoneQ = 1e0 - 2e0 * backboneq;

	/*
	 * Set parameters
	 */
	that->fBeta[0] = (IS_ALMOST_ZERO(beta1) ? 1e-8 : beta1);
	that->fBeta[1] = (IS_ALMOST_ZERO(beta2) ? 1e-8 : beta2);
	that->fAlpha[0] = alpha1;
	that->fAlpha[1] = alpha2;

	/* set time */
	for (i=0; i<=that->fNDates-1; i++) {
	    GtoDayCountFraction(REFDATE, that->fDate[i],
		GTO_ACT_365F, &that->fTime[i] );
	}
	for (i=0; i<=that->fNDates-2; i++) {
	    that->fDt[i] = that->fTime[i+1] - that->fTime[i];
	}
	that->fDt[that->fNDates-1] = 0.;

	/* Set zero curve */
	that->fZcCurve = zcCurve;

#ifndef	VNFM_V5X
	/* Set fZTime */
	that->fZTime[0] = 0.0;    /* equals REFDATE */
	for (i=1;i<=NZDATES-1;i++) {
		GtoDayCountFraction(REFDATE, zcCurve->fArray[i-1].fDate,
			GTO_ACT_365F, &that->fZTime[i]);
	}    
#endif	/*VNFM_V5X*/

	/* check valid data */
	if (VnfmCheckValid(that) != 0) {
	    goto done;
	}

	/* compute coefficients */
	if (VnfmComputeCoeff(that) != SUCCESS)
		goto done;


	/*
	 *
	 */
	status = SUCCESS;
done:
	if (status != 0) {
		GtoErrMsg("%s: failed.\n", routine);
		VnfmFree(that);
		return(NULL);
	}
	return(that);
}


/*f-------------------------------------------------------------
 * Allocate a new VnfmData strtucture with 2F parameters.
 *                                                             
 * <br><br>
 * Creates and initializes the timeline of the two-factor parameter set.
 * The timeline is build according to the following:
 * <br>
 * <br> If "dates" is not and "nDates" is positive, "refDate" and the
 *	array dates are used to build the timeline.
 * <br> If "dates" is NULL (or "nDates" is 0), and "volCurve" is not null,
 *	the dates of "volCurve" are used for the timeline ("refDate" is not
 *	used).
 * <br> If "dates" is NULL (or "nDates" is 0), and "volCurve" is also NULL,
 *	a default timeline of 360 months is used (starting from "refDate").
 * <br>
 * The argument <i> backboneq</i> can be one either
 *      <i> VNFM\_LOGNORMAL\_DIST</i> or
 *	<i> VNFM\_NORMAL\_DIST</i> or
 *      any number between 0 and 0.5.
 * The arguments "beta1", "beta2", "alpha", "alpha2" define
 * the factor parameters
 * and "sigma" and "rho" are the diffusion parameters (the routine
 * initializes the time dependent $\sigma(t)$ and $\rho(t)$ to
 * constant values).
 * Returns NULL if failed.
 */

VnfmData*
VnfmNew2FactSimple(
	TDate refDate,		/* (I) reference date */
	int nDates,		/* (I) # dates (can be -1) */
	TDate *dates,		/* (I) array of dates (can be NULL) */
	TCurve *volCurve,	/* (I) used for timeline only (can be NULL) */
	TCurve *zcCurve,	/* (I) zero curve */
 	double backboneq,	/* (I) backbone WARNING: 0.5=normal */
	double beta1,		/* (I) 2F parameter */
	double beta2,		/* (I) 2F parameter */
	double alpha1,		/* (I) 2F parameter */
	double alpha2,		/* (I) 2F parameter */
	double sigma1,		/* (I) 2F parameter */
	double sigma2,		/* (I) 2F parameter */
	double rho)		/* (I) 2F parameter */
{
static	char	routine[] = "VnfmNew2FactSimple";
	int	status = FAILURE;
	int	i,
		nf = 2;			/* # of factors */
	VnfmData *that = NULL;


	/* */
	if ((nDates > 0) && (dates != (TDate*)NULL))  {
	    /* use input dates */
	    if (refDate == dates[0]) {
		if ((that=VnfmNew(nDates, nf,zcCurve->fNumItems+1)) == NULL)
			goto done;

		/* set dates */
		for (i=0; i<=that->fNDates-1; i++) {
			that->fDate[i] = dates[i];
		}
	    } else {
		if ((that=VnfmNew(nDates+1, nf, zcCurve->fNumItems+1)) == NULL)
		goto done;

		/* set dates */
		that->fDate[0] = refDate;
		for (i=1; i<=that->fNDates-1; i++) {
			that->fDate[i] = dates[i-1];
		}
	    }
	} else if (volCurve != NULL) {
	    /* use vol curve dates */
	    if (volCurve->fBaseDate == volCurve->fArray[0].fDate) {
		nDates = volCurve->fNumItems;
		if ((that = VnfmNew(nDates, nf,zcCurve->fNumItems+1)) == NULL)
			goto done;

		/* set dates */
		for (i=0; i<=that->fNDates-1; i++) {
			that->fDate[i] = volCurve->fArray[i].fDate;
		}
	    } else {
		nDates = volCurve->fNumItems;
		if ((that=VnfmNew(nDates+1, nf,zcCurve->fNumItems+1)) == NULL)
			goto done;

		/* set dates */
		that->fDate[0] = volCurve->fBaseDate;
		for (i=1; i<=that->fNDates-1; i++) {
			that->fDate[i] = volCurve->fArray[i-1].fDate;
		}
	    }
	} else {
	    /* used default timeline */
	    TDateInterval	interval;

	    nDates = 360;
	    if ((that = VnfmNew(nDates, nf,zcCurve->fNumItems+1)) == NULL)
			goto done;
	    for (i=0; i<=that->fNDates-1; i++) {
	        if (GtoMakeDateInterval(i, 'M', &interval) != SUCCESS)
		    goto done;
		GtoDtFwdAny(refDate, &interval, &that->fDate[i]);
	    }
	}

	/*
	 * Set parameters
	 */
	that->fBackBoneQ = 1e0 - 2e0 * backboneq;
  
	that->fBeta[0] = (IS_ALMOST_ZERO(beta1) ? 1e-8 : beta1);
	that->fBeta[1] = (IS_ALMOST_ZERO(beta2) ? 1e-8 : beta2);
	that->fAlpha[0] = alpha1;
	that->fAlpha[1] = alpha2;
	for (i=0; i<=that->fNDates-1; i++) {
	    that->fSigma[0][i] = sigma1;
	    that->fSigma[1][i] = sigma2;
	    that->fRho[0][i]   = rho;
	}

	/* set time */
	for (i=0; i<=that->fNDates-1; i++) {
	    GtoDayCountFraction(REFDATE, that->fDate[i],
		GTO_ACT_365F, &that->fTime[i] );
	}
	for (i=0; i<=that->fNDates-2; i++) {
	    that->fDt[i] = that->fTime[i+1] - that->fTime[i];
	}
	that->fDt[that->fNDates-1] = 0.;


	/* Set zero curve */
	that->fZcCurve = zcCurve;

#ifndef	VNFM_V5X
	/* Set fZTime */
	that->fZTime[0] = 0.0;    /* equals REFDATE */
	for (i=1;i<=NZDATES-1;i++) {
		GtoDayCountFraction(REFDATE, zcCurve->fArray[i-1].fDate,
			GTO_ACT_365F, &that->fZTime[i]);
	}
#endif	/*VNFM_V5X*/



	/* check valid data */
	if (VnfmCheckValid(that) != 0) {
	    GtoErrMsg("%s: invalid parameters.\n", routine);
	    goto done;
	}

	/* compute coefficients */
	if (VnfmComputeCoeff(that) != SUCCESS) {
		goto done;
	}

	/*
	 *
	 */
	status = SUCCESS;
done:
	if (status != 0) {
		GtoErrMsg("%s: failed\n", routine);
		VnfmFree(that);
		return(NULL);
	}

	return(that);
}



/*f-------------------------------------------------------------
 * Smile parameter set.
 *                                                             
 * <br><br>
 * Sets and verifies the smile in a vnfm data structure
 * (default is lognormal/lognormal backbone).
 * Must be called AFTER the mean reversion and betas
 * have been loaded in the object (because in internally sets
 * the reference normal and lognormal volatilities based
 * on the norm of the alpha vector).
 */

DLL_EXPORT(int)
VnfmSetSmile(
	VnfmData *that,		/* (B) vnfm data */
	double bbq,		/* (I) backbone q (0=normal, 1=lognormal) */
	double qLeft,		/* (I) mapping coefficient */
	double qRight,		/* (I) mapping coefficient */
	double qFShift)		/* (I) mapping coefficient */
{
static	char	routine[] = "VnfmSetSmile";
	int	status = FAILURE;
	double	alphaN;		/* alpha norm */
	int	k;

	/* Compute alpha norm */
	alphaN = 0e0;
	for (k=0; k<=that->fNf-1; k++) alphaN += that->fAlpha[k] * that->fAlpha[k];
	if (alphaN <= DBL_EPSILON) {
		GtoErrMsg("%s: norm of alphas (%lf) too small\n.",
			routine, alphaN);
		goto done;
	}
	alphaN = sqrt(alphaN);


	/* Compute reference volatilities */
	if (IS_ALMOST_ZERO(bbq - 0e0)) {
		that->fBackBoneQ = bbq;
		that->fVolLogn = 0e0;
		that->fVolNorm = alphaN;
	} else if (IS_ALMOST_ZERO(bbq - 1e0)) {
		that->fBackBoneQ = bbq;
		that->fVolLogn = alphaN;
		that->fVolNorm = 0e0;
	} else {
		GtoErrMsg("%s: back bone Q (%lf) "
			"must be either 0.0 (normal) and 1.0 (lognormal).\n",
			routine, bbq);
		goto done;
	}

	that->fQLeft   = qLeft;
	that->fQRight  = qRight;
	that->fQFShift = qFShift;


	status = SUCCESS;
done:
	if (status != 0) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*f-------------------------------------------------------------
 * Check valid data.
 *                                                             
 * <br><br>
 * Checks if a set of parameters is valid (valid and increasing dates,
 * positive MR, correlations between -1 and 1, etc\dots).
 * Returns 0 iff input is valid.\\
 * WARNING: does not check for nonnegativity of the correlation matrix
 * above dimension 3.
 */

DLL_EXPORT(int)
VnfmCheckValid(VnfmData *that)
{
static	char	routine[] = "VnfmCheckValid";
	int	status = FAILURE;
	int	i, k;
	double	res;
#undef	CHECK
#define	CHECK(cond, str)	if (!(cond)) {GtoErrMsg("%s: %s.\n",\
				 routine, (str)); goto done;}

	CHECK(that->fNDates > 0, "# of timeline points <= 0");


	for (i=0; i<=that->fNDates-2; i++) {
	    CHECK((that->fDate[i] > 1L), "invalid date");
	}

	for (i=0; i<=that->fNDates-2; i++) {
	    if (that->fDate[i] >= that->fDate[i+1]) {
		GtoErrMsg("%s: date # %d (%s) >= date # %d (%s).\n", 
			routine,
			i,   GtoFormatDate(that->fDate[i]),
			i+1, GtoFormatDate(that->fDate[i+1]));
		goto done;
	    }
	}

	/* beta > 0 */
	for (k=0; k<=that->fNf-1; k++) {
	    if (IS_ALMOST_ZERO(that->fBeta[k])) 
		that->fBeta[k] = 1e-8;
	    CHECK(that->fBeta[k] > 0e0, "beta[k] <= 0");
	}

	/* Check distribution type */

	/* volatilities */
	for (k=0; k<=that->fNf-1; k++) {
	    for (i=0; i<=that->fNDates-1; i++) {
		CHECK(that->fSigma[k][i] >= 0e0, "sigma < 0");
	    }
	}

	/* check correlation positivity */
	if (that->fNf >= 2) {
	    for (i=0; i<=that->fNDates-1; i++) {
		for (k=0; k<(that->fNf)*(that->fNf-1)/2; k++) {
		    if ((res = (1e0 - SQR(that->fRho[k][i]))) < 0e0) {
		    	GtoErrMsg("%s: corr matrix < 0 "
			    "(step %d, nDim=%d, res=%lf)\n",
			    routine, i, that->fNf, res);
		    	goto done;
		    }
		}
	    }
	}

	/* Check correlation positivity
	 * $$$ ONLY DOES IT FOR NDIM <= 3 
	 */
	if (that->fNf >= 3) {
	    for (i=0; i<=that->fNDates-1; i++) {
		res = 1e0 - SQR(that->fRho[0][i])
		          - SQR(that->fRho[1][i])
		          - SQR(that->fRho[2][i])
		          + 2e0 * that->fRho[0][i]
		                * that->fRho[1][i]
		                * that->fRho[2][i];

		if (res < 0e0) {
		    GtoErrMsg("%s: corr matrix < 0 "
			"(step %d, nDim=%d, res=%lf)\n",
			routine, i, that->fNf, res);
		    goto done;
		}
	    }
	}

/*
	if (that->fNf >= 4) {
		GtoErrMsg("%s: not implemented for dim > 3\n", routine);
		goto done;
	}
*/

	/* Check backbone Q */
	if ((!IS_ALMOST_ZERO(that->fBackBoneQ - 0e0)) &&
	    (!IS_ALMOST_ZERO(that->fBackBoneQ - 1e0))) {
		GtoErrMsg("%s: back bone Q (%lf) "
			"must be either 0.0 (normal) and 1.0 (lognormal).\n",
			routine, that->fBackBoneQ);
		goto done;
	}


	/* check valid curve */
	if (that->fZcCurve->fBaseDate < that->fDate[0]) {
	    GtoErrMsg("%s: zero ref date %s < vol ref date %s\n",
		routine,
		GtoFormatDate(that->fZcCurve->fBaseDate),
		GtoFormatDate(that->fDate[0]));
	    goto done;
	}
  

#ifndef	VNFM_V5X
	CHECK((that->fNZDates IS that->fZcCurve->fNumItems+1), 
		"Size of fZTime <> Number of Zero curve points + 1");
#endif	/*VNFM_V5X*/


	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
#ifndef NO_LOGGING
	    GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
		"%s: invalid model parameters:\n", routine);\
		VnfmFpWrite(that, vnfmFpLog));
#endif
	    GtoErrMsg("%s: invalid model parameters.\n", routine);
	}

	return(status);
#undef	CHECK
}



/*f-------------------------------------------------------------
 * Get correlation matrix determinant.
 *                                                             
 * <br><br>
 * Calculates and returns the minor determinant of order "order"
 * of the correlation matrix at time point index "tpIdx".
 * Returns SUCCESS or FAILURE.
 */

DLL_EXPORT(int)
VnfmCorrMinorDeterm(
	VnfmData *that,		/* (I) model parameters */
	int tpIdx,		/* (I) time point index */
	int order, 		/* (I) order */
	double *determ)		/* (O) determinant */
{
static	char	routine[] = "VnfmCorrMinorDeterm";
	int	i, k;
	double	res;

	if (order > that->fNf-1) {
		GtoErrMsg("%s: order %d > num fact %d.\n",
			routine, order, that->fNf);
		return(FAILURE);
	}

	switch (order) {
	case 0:
		*determ = 1e0;
		break;
	case 1:
		*determ = 1e0 - SQR(that->fRho[0][tpIdx]);
		break;
	case 2:
		*determ = 1e0 - SQR(that->fRho[0][tpIdx])
		          - SQR(that->fRho[1][tpIdx])
		          - SQR(that->fRho[2][tpIdx])
		          + 2e0 * that->fRho[0][tpIdx]
		                * that->fRho[1][tpIdx]
		                * that->fRho[2][tpIdx];
		break;
	default:
		GtoErrMsg("%s: can't compute order %d.\n",
			routine, order);
		return(FAILURE);
	}
	return(SUCCESS);
}




/*------------------------------------------------------
 * Allocate memory for "nf" dimensions and a timeline of size "nDates"
 * in the data structure "that" (for 
 * storage of parameters and internal coefficients).
 * <b> Remark.</b> This routine should be called internally only.
 * Returns 0 iff OK.
 */

static	int
VnfmMemoryAlloc(VnfmData *that, int nDates, int nf, int nZDates)
{
	int	errCode=0,
		  i;

#undef	CHECK
#define	CHECK(ptr)	if ((ptr) == NULL) {errCode = (-4); goto done;}

	if ((nDates <= 0) || (nf <= 0)) return(-1);
  
	CHECK(that->fDate = NEW_ARRAY(TDate,  nDates+1))
	CHECK(that->fTime = NEW_ARRAY(double, nDates+1))
	CHECK(that->fDt   = NEW_ARRAY(double, nDates+1))

#ifndef	VNFM_V5X
	if (nZDates > 0) {
		CHECK(that->fZTime = NEW_ARRAY(double, nZDates+1))
		CHECK(that->fRate =  NEW_ARRAY(double, nZDates+1))
		CHECK(that->fZero =  NEW_ARRAY(double, nZDates+1))
	}
	else {
		that->fZTime = NULL;
		that->fRate  = NULL;
		that->fZero  = NULL;
	}
#else
	CHECK(that->fRate =  NEW_ARRAY(double, nDates+1))
	CHECK(that->fZero =  NEW_ARRAY(double, nDates+1))
#endif	/*VNFM_V5X*/

	/*
	 *
	 */

	CHECK(that->fBeta  = NEW_ARRAY(double, nf))
	CHECK(that->fAlpha = NEW_ARRAY(double, nf))

	CHECK(that->fSigma = NEW_ARRAY(double*, nf))
	for (i=0; i<=nf-1; i++) {
		CHECK(that->fSigma[i] = NEW_ARRAY(double, nDates+1))
	}

	if (nf > 1) {
		CHECK(that->fRho   = NEW_ARRAY(double*, nf*(nf-1)/2))
		for (i=0; i<=nf*(nf-1)/2-1; i++) {
			CHECK(that->fRho[i] = NEW_ARRAY(double, nDates+1))
		}
	}

	CHECK(that->fJ     = NEW_ARRAY(double*, nf*(nf+1)/2))
	for (i=0; i<=nf*(nf+1)/2-1; i++) {
		CHECK(that->fJ[i] = NEW_ARRAY(double, nDates+1))
	}

	that->fNDates = nDates;
#ifndef	VNFM_V5X
	that->fNZDates = nZDates;
#endif	/*VNFM_V5X*/
	that->fNf = nf;

	/* Default backbone and smile parameters */
	that->fBackBoneQ     = 1e0;
	that->fVolNorm = 1e0;
	that->fVolLogn = 1e0;
	that->fQLeft   = 1e0;	/* lognormal mapping default */
	that->fQRight  = 1e0;
	that->fQFShift = 0e0;

done:
	if (errCode != 0) {
	     GtoErrMsg("VnfmMemoryAlloc: malloc failure "
		"(nDates=%d, nDim=%d, nZDates=%d)\n",
		nDates, nf, nZDates);
	}
	return(errCode);

#undef	CHECK
}


/*------------------------------------------------------
 * Frees memory allocated by <i> VnfmAlloc</i>.
 * <b> Remark.</b> This routine should be called internally only.
 * Returns 0 iff OK.
 */

static	int
VnfmMemoryFree(VnfmData *that)
{
	int	i;

	if ((that != NULL) && (that->fNDates > 0)) {

		FREE_ARRAY(that->fDate);
		FREE_ARRAY(that->fTime);
		FREE_ARRAY(that->fDt);
		FREE_ARRAY(that->fRate);
		FREE_ARRAY(that->fZero);
#ifndef	VNFM_V5X
		FREE_ARRAY(that->fZTime);
#endif	/*VNFM_V5X*/

		FREE_ARRAY(that->fBeta);
		FREE_ARRAY(that->fAlpha);

		for (i=0; i<=that->fNf-1; i++) {
			FREE_ARRAY(that->fSigma[i]);
		}
		for (i=0; i<=that->fNf*(that->fNf-1)/2-1; i++) {
			FREE_ARRAY(that->fRho[i]);
		}
		for (i=0; i<=that->fNf*(that->fNf+1)/2-1; i++) {
			FREE_ARRAY(that->fJ[i]);
		}
		FREE_ARRAY(that->fSigma);
		FREE_ARRAY(that->fJ);
		if (that->fNf > 1) 
			FREE_ARRAY(that->fRho);


		that->fNDates = -1;
		that->fNf = -1;
#ifndef	VNFM_V5X
		that->fNZDates = -1;
#endif	/*VNFM_V5X*/

	}
	return(SUCCESS);
}


