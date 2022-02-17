/****************************************************************
 * Module:	VNFM
 * Submodule:	ANLY
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <math.h>

#include "convert.h"
#include "date_sup.h"

#define	_vnfm_SOURCE
#include "vnfmanly.h"

#include "drlmem.h"			/* Drl memory management */
#include "drllineq.h"			/* DrlEigen */
#include "drlio.h"			/* DrlFPrintf */

#undef	VNFM_PRINT_LEVEL
#define	VNFM_PRINT_LEVEL	2	/* print level 2 for logging */

/*f-------------------------------------------------------------
 * Compute the average volatility of a forward rate.
 *                                                             
 * <br><br>
 * For a set of parameters "that",
 * computes and returns the average volatility between
 * two times "S1" and "S2" of a forward rate
 * of forward maturity "tMat", time to reset "tReset" and
 * frequency "freq" (1,2,4,12).\\
 * <b> Warning.</b> The routine assumes that the internal
 * coefficients have been computed (using "VnfmComputeCoeff").
 */

DLL_EXPORT(int)
VnfmAvgQBVol(
	VnfmData *that,
	double S1,		/* (I) start */
	double S2,		/* (I) expiration */
	double tReset,		/* (I) rate reset */
	double tMat,		/* (I) rate maturity */
	int freq,		/* (I) rate frequency */
	KVolType vType,		/* (I) LOGVOL, NORMVOL */
	double *retVal)		/* (O) volatility */
{
static	char	routine[] = "VnfmAvgQBVol";
	double	x[VNFM_NDIMMAX],
		j00[VNFM_NDIMMAX*(VNFM_NDIMMAX+1)/2],
		vol;
	int	i, j, nDim = NDIM, jIdx;
	double	yield;

	/* check */
	if ((S2 > tReset) || (S1 > S2)) {
	    GtoErrMsg("%s: tExp=%lf > tReset=%lf or tStart=%lf > tExp\n",
		routine, S2, tReset, S1);
	    return(FAILURE);
	}

    if (IS_ALMOST_ZERO(S2-S1))
    {
        *retVal = 0.0;
        return SUCCESS;
    }

	/* compute the J integrals */
	VnfmJ(that, S1, S2, tReset, j00);

	/* compute the Q coefficients */
	if (VnfmQBCoeff(that, tReset, tReset, tMat, freq, &yield, x) != SUCCESS)
		return(FAILURE);

	/* vol computation */
	vol = 0e0;
	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) {
	    jIdx = JIDX(i, j);
	    vol += (i == j ? 1e0 : 2e0) * x[i] * x[j] * j00[jIdx];
	}

	/* $$$ Return NEGATIVE volatility if the covariance
	 * matrix is negative.
	 */
	vol = SQRT_SIGN(-vol / (S2-S1));

	/* vol being computed is bp vol */
	if (vType == LOGVOL)
	{
		if (IS_ALMOST_ZERO(yield) ||
		    yield < 0e0)
		{
			GtoErrMsg("%s: forward rate (%lf) can not be <= 0 "
				  "while the output vol is lognormal.\n", 
				routine, yield);
			return (FAILURE);
		}
		else
			vol /= yield;
	}

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "VnfmAvgBVol(%d): "\
		"S1=%lf S2=%lf T=%lf tMat=%lf freq=%d vol=%lf\n",\
		__LINE__, S1, S2, tReset, tMat, freq, vol));
#endif

	*retVal = vol;
	return(SUCCESS);
}



/*f-------------------------------------------------------------
 * Compute the average correlation of two forward rates.
 *                                                             
 * <br><br>
 * For a set of parameters "that",
 * computes and returns the average correlation between 
 * time $0$ and "S" of two rates
 * that have time to reset "tReset1" and "tReset2",
 * forward maturitites "tMat1" and "tMat2",
 * and frequencies "freq1" and "freq2" (0 for MM rate, 1,2,4,12).\\
 * <b> Warning.</b> The routine assumes that the internal
 * coefficients have been computed \\ (using "VnfmComputeCoeff").
 */

DLL_EXPORT(int)
VnfmAvgQBCorr(
	VnfmData *that,		/* (I) model parameters */
	double S,		/* (I) expiration */
	double tReset1,		/* (I) reset rate 1 */
	double tMat1,		/* (I) maturity rate 1 */
	int freq1,		/* (I) freq rate 1 */
	double tReset2,		/* (I) reset rate 2 */
	double tMat2,		/* (I) maturity rate 2 */
	int freq2,		/* (I) freq rate 2 */
	double *retVal)		/* (O) correlation */
{
static	char	routine[] = "VnfmAvgQBCorr";
	int	status = FAILURE;

	double	x[VNFM_NDIMMAX],
		y[VNFM_NDIMMAX],
		j00[VNFM_NDIMMAX*(VNFM_NDIMMAX+1)/2],
		vol1, vol2, cor;
	int	i, j, nDim = NDIM, jIdx;

	double	fwdRate1, fwdRate2;

	/* check */
	if ((S > tReset1) || (S > tReset2)) {
	    GtoErrMsg("%s: tExp=%lf > tReset1=%lf or tReset2=%lf\n",
		routine, S, tReset1, tReset2);
	    return(FAILURE);
	}

    if (IS_ALMOST_ZERO(S))
    {
        *retVal = 0.0;
        return SUCCESS;
    }

	IF_FAILED_DONE(VnfmQBCoeff(that, tReset1, tReset1, tMat1, freq1, &fwdRate1, x));
	IF_FAILED_DONE(VnfmQBCoeff(that, tReset2, tReset2, tMat2, freq2, &fwdRate2, y));

	/* compute the J integrals */
	VnfmJ(that, 0e0, S, S, j00);

	/* vol computation */
	vol1 = vol2 = cor = 0e0;
	for (i=0; i<=nDim-1; i++)
	for (j=0; j<=nDim-1; j++) {
	    jIdx = (i <= j ? JIDX(i, j) : JIDX(j, i));

	    vol1 += x[i] * x[j] * j00[jIdx]
			* exp(-(BETA[i]+BETA[j])*(tReset1 - S));
	    vol2 += y[i] * y[j] * j00[jIdx]
			* exp(-(BETA[i]+BETA[j])*(tReset2 - S));
	    cor  += x[i] * y[j] * j00[jIdx]
			* exp(-BETA[i]*(tReset1 - S))
			* exp(-BETA[j]*(tReset2 - S));
	}

	cor = cor / (SQRT_SIGN(vol1) * SQRT_SIGN(vol2));

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(
	    DrlFPrintf(vnfmFpLog, "%s: "\
		"S=%lf tReset1=%lf tReset2=%lf tMat1=%lf tMat2=%lf\n",\
		routine, S, tReset1, tReset2, tMat1, tMat2);\
	    for (i=0; i<=nDim-1; i++)\
		DrlFPrintf(vnfmFpLog, "\t[%2d]\tx=%lf\ty=%lf\n",\
			i, x[i], y[i]));
#endif


	*retVal = cor;

	/* OK */
	status = SUCCESS;
done:
	return(status);
}


/*f-------------------------------------------------------------
 * Compute the average correlation of two forward rates (arbitrary).
 *                                                             
 * <br><br>
 * For a set of parameters "that",
 * computes and returns the average correlation between 
 * time $0$ and "S" of two rates
 * that have time to reset "tReset1" and "tReset2",
 * forward maturitites "tMat1" and "tMat2",
 * and frequencies "freq1" and "freq2" (0 for MM rate, 1,2,4,12).\\
 * <b> Warning.</b> The routine assumes that the internal
 * coefficients have been computed \\ (using "VnfmComputeCoeff").
 */

DLL_EXPORT(int)
VnfmAvgQBCorr2(
	VnfmData *that,		/* (I) model parameters */
	double Obs,		/* (I) observation */
	double S,		/* (I) expiration */
	double tReset1,		/* (I) reset rate 1 */
	double tMat1,		/* (I) maturity rate 1 */
	int freq1,		/* (I) freq rate 1 */
	double tReset2,		/* (I) reset rate 2 */
	double tMat2,		/* (I) maturity rate 2 */
	int freq2,		/* (I) freq rate 2 */
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
	if ((S > tReset1) || (S > tReset2)) {
	    GtoErrMsg("%s: tExp=%lf > tReset1=%lf or tReset2=%lf\n",
		routine, S, tReset1, tReset2);
	    return(FAILURE);
	}

	if (Obs > S)  {
	    GtoErrMsg("%s: tObs=%lf > tExp=%lf. \n",
		routine, Obs, S);
	    return(FAILURE);
	}

    if (IS_ALMOST_ZERO(S-Obs))
    {
        *retVal = 0.0;
        return SUCCESS;
    }


	/* compute vol coefficients */
	if (VnfmQBCoeff(that, tReset1, tReset1, tMat1, freq1, &fwdRate1, x) != SUCCESS)
		return(FAILURE);
	if (VnfmQBCoeff(that, tReset2, tReset2, tMat2, freq2, &fwdRate2, y) != SUCCESS)
		return(FAILURE);

	/* compute the J integrals */
	VnfmJ(that, Obs, S, S, j00);

	/* vol computation */
	vol1 = vol2 = cor = 0e0;
	for (i=0; i<=nDim-1; i++)
	for (j=0; j<=nDim-1; j++) {
	    jIdx = (i <= j ? JIDX(i, j) : JIDX(j, i));

	    vol1 += x[i] * x[j] * j00[jIdx]
			* exp(-(BETA[i]+BETA[j])*(tReset1 - S));
	    vol2 += y[i] * y[j] * j00[jIdx]
			* exp(-(BETA[i]+BETA[j])*(tReset2 - S));
	    cor  += x[i] * y[j] * j00[jIdx]
			* exp(-BETA[i]*(tReset1 - S))
			* exp(-BETA[j]*(tReset2 - S));
	}

	cor = cor / (sqrt(vol1) * sqrt(vol2));

	*retVal = cor;

	return(SUCCESS);
}



/*f-------------------------------------------------------------
 * For a set of parameters "that",
 * computes and returns the average normal spread vol between
 * time $0$ and "S" of two rates
 * that have time to reset "tReset1" and "tReset2",
 * forward maturitites "tMat1" and "tMat2",
 * and frequencies "freq1" and "freq2" (0 for MM rate, 1,2,4,12).\\
 * {\bf Warning.} The routine assumes that the internal
 * coefficients have been computed \\ (using "VnfmComputeCoeff").
 */

DLL_EXPORT(int)
VnfmAvgQBSprdVol(
    VnfmData *that,     /* (I) model parameters */
    double S,       /* (I) expiration */
    double tReset1,     /* (I) reset rate 1 */
    double tMat1,       /* (I) maturity rate 1 */
    int freq1,      /* (I) freq rate 1 */
    double tReset2,     /* (I) reset rate 2 */
    double tMat2,       /* (I) maturity rate 2 */
    int freq2,      /* (I) freq rate 2 */
    double *retVal)     /* (O) correlation */
{
static  char    routine[] = "VnfmAvgQBSprdVol";
    int status = FAILURE;

    double  x[VNFM_NDIMMAX],
        y[VNFM_NDIMMAX],
        j00[VNFM_NDIMMAX*(VNFM_NDIMMAX+1)/2],
        sprdVol;
    int i, j, nDim = NDIM, jIdx;

    double  fwdRate1, fwdRate2;

    /* check */
    if ((S > tReset1) || (S > tReset2)) {
        GtoErrMsg("%s: tExp=%lf > tReset1=%lf or tReset2=%lf\n",
        routine, S, tReset1, tReset2);
        return(FAILURE);
    }

    if (IS_ALMOST_ZERO(S))
    {
        *retVal = 0.0;
        return SUCCESS;
    }

    IF_FAILED_DONE(VnfmQBCoeff(that, tReset1, tReset1, tMat1, freq1, &fwdRate1,
x));
    IF_FAILED_DONE(VnfmQBCoeff(that, tReset2, tReset2, tMat2, freq2, &fwdRate2,
y));

    /* compute the J integrals */
    VnfmJ(that, 0e0, S, S, j00);

    /* vol computation */
    sprdVol = 0e0;
    for (i=0; i<=nDim-1; i++)
    for (j=0; j<=nDim-1; j++) {
        jIdx = (i <= j ? JIDX(i, j) : JIDX(j, i));

        sprdVol += (  x[i] * exp(-BETA[i]*(tReset1 - S))
                    - y[i] * exp(-BETA[i]*(tReset2 - S)))
                 * (  x[j] * exp(-BETA[j]*(tReset1 - S))
                    - y[j] * exp(-BETA[j]*(tReset2 - S)))
                 * j00[jIdx];

    }

    sprdVol = SQRT_SIGN(sprdVol / S);

#ifndef NO_LOGGING
    GTO_IF_LOGGING(
        DrlFPrintf(vnfmFpLog, "%s: "\
        "S=%lf tReset1=%lf tReset2=%lf tMat1=%lf tMat2=%lf\n",\
        routine, S, tReset1, tReset2, tMat1, tMat2);\
        for (i=0; i<=nDim-1; i++)\
        DrlFPrintf(vnfmFpLog, "\t[%2d]\tx=%lf\ty=%lf\n",\
            i, x[i], y[i]));
#endif


    *retVal = sprdVol;

    /* OK */
    status = SUCCESS;
done:
    return(status);
}


/*f-------------------------------------------------------------
 * Compute the average volatility of a forward O/N rate.
 *                                                             
 * <br><br>
 * For a set of parameters "that",
 * computes and returns the average volatility of a forward O/N rate
 * of reset time "S".\\
 * <b> Warning.</b> The routine assumes that the internal
 * coefficients have been computed (using "VnfmComputeCoeff").
 */

double
VnfmAvgONVol(
	VnfmData *that,		/* (I) model parameters */
	double S)		/* (I) time to expiration */
{
	double	j00[VNFM_NDIMMAX*(VNFM_NDIMMAX+1)/2],
		vol;
	int	i, j, nDim = NDIM, jIdx;

	/* compute the J integrals */
	VnfmJ(that, 0e0, S, S, j00);

	/* vol computation */
	vol = 0e0;
	for (i=0; i<=nDim-1; i++)
	for (j=i; j<=nDim-1; j++) {
	    jIdx = JIDX(i, j);
	    vol += (i == j ? 1e0 : 2e0) * j00[jIdx];
	}
	vol = sqrt(vol / S);


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "VnfmONVol(%d):  S=%lf vol=%lf\n",\
		__LINE__, S, vol));
#endif

	return(vol);
}



/*f--------------------------------------------------------------
 * Computes the average volatility curve of a given rate.
 *                                                             
 * <br><br>
 * Generates the model base or forward volatility curve for a given index.
 * On entry, "that" contains the model parameters.
 * "lMat" is the maturity of the rate, "freq" its 
 * frequency (0 for MM), * and
 * "dates" is an array of length "nDates" of dates at which
 * the base volatility is to be computed, and "vol" is an
 * array of length "nDates" that contains, on exit, the 
 * base volatilities.
 * Returns 0 iff OK.
 */

int
VnfmVolCurve(
	VnfmData *that,		/* (I) model parameters */
	TDateInterval lMat,	/* (I) maturity of the rate */
	int freq,		/* (I) frequency of rate */
	int volType,		/* (I) 0=basevol, 1=fwdvol */
	TDate bvRefDate,	/* (I) base vol reference date */
	int nBvDates,		/* (I) # base vol dates */
	TDate *bvDates,		/* (I) vol dates */
	KVolType vType,		/* (I) LOGVOL, NORMVOL */
	double *bvRates)	/* (O) vol curve */
{
static	char		routine[] = "VnfmVolCurve";
	double		S1, T, tMat;
	int		status = FAILURE,
			i,
			onFlag;		/* TRUE = generate O/N vol */

	/* compute coefficients (you may have forgotten it before ..) */
	if (VnfmComputeCoeff(that) != 0)
			goto done;

	if (GtoDateIntervalToYears(&lMat, &tMat) != SUCCESS) goto done;
	onFlag = IS_ALMOST_ZERO(tMat);


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "%s: tMat=%lf freq=%d "\
		"volType=%d bvRefDate=%10s\n",\
		routine, tMat, freq, volType, GtoFormatDate(bvRefDate)));
#endif

	/* compute volatilities */
	switch (volType) {
	case 0:
	    /* base volatility */
	    GtoDayCountFraction(REFDATE, bvRefDate, GTO_ACT_365F, &S1);

	    for (i=0; i<=nBvDates-1; i++) {
		GtoDayCountFraction(REFDATE, bvDates[i], GTO_ACT_365F, &T);
		T = MAX(T, 1e0/365e0); /* to avoid T=0*/
		if (T <= S1) {
		    GtoErrMsg("%s: date # %d (%s) <= bv ref date (%s)\n",
			routine, i,
			GtoFormatDate(bvDates[i]),
			GtoFormatDate(bvRefDate));
		    goto done;
		}
		/* compute vol */
		if (VnfmAvgQBVol(that, S1, T, T, tMat, freq, vType, &bvRates[i])
			!= SUCCESS) goto done;

#ifndef	NO_LOGGING
		GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
			 "\tS1=%lf\tT=%lf\tvol=%lf\n", S1, T, bvRates[i]));
#endif

	    }
	    break;

	case 1:
	    /* fwd volatility */
	    GtoDayCountFraction(REFDATE, bvRefDate, GTO_ACT_365F, &S1);

	    for (i=0; i<=nBvDates-1; i++) {
		GtoDayCountFraction(REFDATE, bvDates[i], GTO_ACT_365F, &T);
		T = MAX(T, 1e0/365); /* to avoid T=0*/
		if (T <= S1) {
		    GtoErrMsg("%s: date # %d (%s) <= bv ref date (%s)\n",
			routine, i,
			GtoFormatDate(bvRefDate),
			GtoFormatDate(bvDates[i]));
		    goto done;
		}
		if (VnfmAvgQBVol(that, S1, T, T, tMat, freq, vType, &bvRates[i])
			!= SUCCESS) goto done;
		S1 = T;
	    }
	    break;

	default:
	    GtoErrMsg("%s: bad volType %d\n", routine, volType);
	    goto done;
	}

	/*
	 * We are done
	 */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*f--------------------------------------------------------------
 * Compute the swaption volatility matrix.
 *                                                             
 * <br><br>
 * Compute swaption matrix for a set of parameters "that".
 * On entry, "lExp" is the array of length "nExp" of 
 * intervals to expiration,
 * "lMat" the array of length "nMat" of maturities,
 * "type" defines the type of the matrix (1 for upper-triangular,
 * 0 for vertical),
 * "freq" the frequency of the rate,
 * "adjFlag" is TRUE if the trinomial adjustment is to be made,
 * and "vol" is a matrix of size "nExp" $x$ "nMat".\\
 * On exit, "vol[i][j]" contains the volatility for
 * the $i$-th expiration and $j$-th maturity.\\
 * Returns 0 iff OK.
 */

int
VnfmSwaptionVolMatrix(
	VnfmData *that,
	int type,		/* (I) vertical(0), diagonal (1)  */
	int freq,		/* (I) frequency of rate */
	int adjFlag,		/* (I) no adj (0), trino (1) */
	int nExp,		/* (I) # of exp intervals */
	TDateInterval *lExp,	/* (I) arry of exp intervals */
	int nMat,		/* (I) # of mat intervals */
	TDateInterval *lMat,	/* (I) array of mat intervals */
	KVolType vType,		/* (I) LOGVOL, NORMVOL */
	double **vol)		/* (O) swaption volatilities */
{
static	char	routine[] = "VnfmSwaptionMatrix";
	int	status = FAILURE;
	int	i, j;
	TDate	expDate, matDate;
	double	S, T, tExp, tMat;
	/*
	 *
	 */

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "%s: refDate=%10s type=%d freq=%d\n",\
		routine, GtoFormatDate(REFDATE), type, freq));
#endif


	switch (type) {
	case 1:
	    /*
	     * DIAGONAL
	     */
	    for (j=0; j<=nExp-1;j++)
	    for (i=0; i<=nMat-1; i++) {

		if (GtoDateIntervalToYears(&lExp[j], &tExp) != SUCCESS)
				goto done;
		if (GtoDateIntervalToYears(&lMat[i], &tMat) != SUCCESS)
				goto done;

		if (tMat > tExp) {
		    GtoDtFwdAny(REFDATE, &lExp[j], &expDate);
		    GtoDtFwdAny(REFDATE, &lMat[i], &matDate);
		    GtoDayCountFraction(REFDATE, expDate, GTO_ACT_365F, &S);
		    if (IS_ALMOST_ZERO(S)) S = 1e0/365e0;
		    GtoDayCountFraction(REFDATE, matDate, GTO_ACT_365F, &T);
		    if (VnfmAvgQBVol(that, 0e0, S, S, T-S, freq, vType, &vol[j][i])
			!= SUCCESS) goto done;
		} else {
		    vol[j][i] = -1.0;
		}

#ifndef	NO_LOGGING
		GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
			"\t[nExp=%2d, %3s] [nMat=%2d, %3s]  vol = %lf\n",\
			j, GtoFormatDateInterval(&lExp[j]),\
			i, GtoFormatDateInterval(&lMat[i]),\
			vol[j][i]));
#endif


		/* Volatility adjustment */
		switch (adjFlag) {
		case 1:
			GtoErrMsg("%s: trino vol not available\n", routine);
			goto done;
		default:
			break;
		}

	    }
	    break;

	case 0:
	    /* vertical matrix */
	    for (j=0; j<=nExp-1;j++)
	    for (i=0; i<=nMat-1; i++) {
		GtoDtFwdAny(REFDATE, &lExp[j], &expDate);
		GtoDtFwdAny(expDate, &lMat[i], &matDate);
		GtoDayCountFraction(REFDATE, expDate, GTO_ACT_365F, &S);
		if (IS_ALMOST_ZERO(S)) S = 1e0/365e0;

		GtoDayCountFraction(expDate, matDate, GTO_ACT_365F, &tMat);

		if (VnfmAvgQBVol(that, 0e0, S, S, tMat, freq, vType, &vol[j][i])
		    != SUCCESS) goto done;


#ifndef	NO_LOGGING
		GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
			"\t[nExp=%2d, %3s] [nMat=%2d, %3s]  vol = %lf\n",\
			j, GtoFormatDateInterval(&lExp[j]),\
			i, GtoFormatDateInterval(&lMat[i]),\
			vol[j][i]));
#endif


		/* Volatility adjustment */
		switch (adjFlag) {
		case 1:
			GtoErrMsg("%s: trino vol not available\n", routine);
			goto done;
		default:
			break;
		}


	    }
	    break;
	default:
	    GtoErrMsg("%s: bad type %d\n", routine, type);
	    goto done;
	}

	/* made it through OK */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*f--------------------------------------------------------------
 * Compute the average correlation matrix.
 *                                                             
 * <br><br>
 * Compute average correlation matrix for given
 * two-factor parameters. On entry,
 * "that" contains the two-factor parameters,
 * "lExp" is the date interval
 * to option expiration,
 * and "lMat" is an array of size "nMat" containing the maturities
 * of the rate between which the corraelation is to be computed.
 * The matrix "corr" should be allocated before the function call.
 * On exit, it contains the cross correlations. 
 * The routine assumes that all internal intermediate coefficients
 * have been computed (calling <i> VnfmComputeCoeff</i>).\\
 * Returns 0 iff OK.
 */

int
VnfmAvgCorrMatrix(
	VnfmData *that,
	TDateInterval lExp,	/* (I) exp interval */
	int nMat,		/* (I) # of mat intervals */
	TDateInterval *lMat,	/* (I) array of mat intervals */
	double **corr)		/* (O) correlations (nMat x nMat) */
{
static	char	routine[] = "VnfmAvgCorrMatrix";
	int	status = FAILURE,
		i, j;
	TDate	expDate, matDate1, matDate2;
	double	S, T1, T2;
	/*
	 *
	 */

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "%s: REFDATE=%s exp=%s\n", routine,\
		GtoFormatDate(REFDATE),\
		GtoFormatDateInterval(&lExp)));
#endif



	/*
	 *
	 */
	if (GtoDtFwdAny(REFDATE, &lExp, &expDate) != SUCCESS)
		goto done;

	if (GtoDayCountFraction(REFDATE, expDate, GTO_ACT_365F, &S) != SUCCESS)
		goto done;


	for (i=0; i<=nMat-1; i++) 
	for (j=0; j<=i;      j++) {
		GtoDtFwdAny(expDate, &lMat[i], &matDate1);
		GtoDtFwdAny(expDate, &lMat[j], &matDate2);
		GtoDayCountFraction(expDate, matDate1, GTO_ACT_365F, &T1);
		GtoDayCountFraction(expDate, matDate2, GTO_ACT_365F, &T2);

		if (VnfmAvgQBCorr(that,
			S,
			S, T1, (T1 > .75 ? 2 : 0),
			S, T2, (T2 > .75 ? 2 : 0),
			&corr[i][j]) != SUCCESS) goto done;


		corr[j][i] = corr[i][j];
	}

	/* made it through */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}

/*f--------------------------------------------------------------
 * Compute average volatility curve of a rate (TCurve output).
 *                                                             
 * <br><br>
 * Generates the model base volatility curve.
 * \vspace{2mm}\par\noindent<b> Argument Details:</b>
 * <br>
 * <br>[that] parameter data structure.
 * <br>[lMat] the maturity date interval of the rate,
 * <br>[dates] is an array of length "nDates" of dates at which
 * the base volatility is to be computed
 * (if "dates" is <i> NULL</i>, the set of dates of the model
 * data structure is used).
 * <br>[volType] the type of volatility curve to be computed:
 * 0 for base volatility, 1 for forward volatility.
 * <br>[bvRefDate] the forward reference date for the volatility
 * curve to be computed. For today's base volatility curve,
 * pass the refernce date of the model parameters.
 * <br>
 * On exit, "bvCurve" contains the base volatility curve.
 * Returns 0 iff OK.
 */

int
VnfmGenerateVolTCurve(
	VnfmData *that,		/* (I) model parameters */
	TDateInterval lMat,	/* (I) maturity of rate */
	int freq,		/* (I) rate frequency */
	int volType,		/* (I) 0=base volatility, 1=forward */
	TDate bvRefDate,	/* (I) fwd vol ref date */
	int nDates,		/* (I) # base vol dates */
	TDate *dates,		/* (I) base vol dates (or NULL) */
	TCurve **bvCurve)	/* (O) output base vol curve */
{
static	char		routine[] = "VnfmGenerateBaseVolCurve";
	double		tMat,
			*vols = NULL;

	int		status = FAILURE,
			i,
			freeDatesFlag = FALSE,
			onFlag;		/* TRUE = generate O/N vol */

	*bvCurve = NULL;
	if (bvRefDate <= 0L) bvRefDate = REFDATE;



	/* compute coefficients (you may have forgotten it before ..) */
	if (VnfmComputeCoeff(that) != 0) goto done;



#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog, "%s:\n", routine);\
	    VnfmPrintCoeff(that, vnfmFpLog));
#endif


	if (GtoDateIntervalToYears(&lMat, &tMat) != SUCCESS)
		goto done;
	onFlag = (tMat <= 2.739726e-4 ? TRUE : FALSE);

	/*----------------------------------------------*
	 * Generate the base Vol dates
	 *----------------------------------------------*/

	if ((dates == (TDate*)NULL) || (nDates <= 0)) {
		if ((dates = NEW_ARRAY(TDate, that->fNDates)) == NULL)
			goto done;
		nDates  = that->fNDates;
		freeDatesFlag = TRUE;
		for (i=0; i<=nDates-1; i++) {
			dates[i] = that->fDate[i];
		}
		/* to avaoid spot date */
		dates[0] += 1L;

	}

	/* create the base vol curve */
	*bvCurve = GtoNewTCurve(
			bvRefDate,
			nDates,
			1e0,
			GTO_ACT_365F);
	if (*bvCurve == NULL) goto done;

	/* basis = frequency */
	if (onFlag) {
		(*bvCurve)->fBasis = 1e0;
	} else if (freq == 0) {
		(*bvCurve)->fBasis = (double) 1e0/tMat;
	} else {
		(*bvCurve)->fBasis = (double) freq;
	}


	/* compute base volatilities and put it in the curve */
	if ((vols = NEW_ARRAY(double, nDates)) == NULL) goto done;
	if (VnfmVolCurve(that,
			lMat, freq, volType,
			bvRefDate, nDates, dates, LOGVOL, vols) != SUCCESS)
		goto done;

	for (i=0; i<=nDates-1; i++) {
		(*bvCurve)->fArray[i].fDate = dates[i];
		(*bvCurve)->fArray[i].fRate = vols[i];
	}


	/*
	 * We are done
	 */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		if (*bvCurve != NULL) GtoFreeTCurve(*bvCurve);
		GtoErrMsg("%s: failed\n", routine);
	}
	if (vols != NULL) FREE((void*) vols);
	if (freeDatesFlag != FALSE) FREE(dates);
	return(status);
}


/*f--------------------------------------------------------------
 * Calculates orthogonal factors.
 *                                                             
 * <br><br>
 * This routine computes the sets of orthogonal factors of 
 * a series of rates in a given observation period.
 * <br>
 * <br>[that] parameter data structure.
 * <br>[tStart] the start ot the observation period,
 * <br>[tEnd] the end ot the observation period (can be equal to tStart for
 *	instantaneous factors,
 * <br>[numRates] number of rates in the factor matrix,
 * <br>[rateExpYrs] array of rate resets as offset in years from "tEnd"
 * 	(length "numRates"),
 * <br>[rateFreq] array of rate frequencies (0,1,2,4,12)
 * 	(length "numRates"),
 * <br>[rateMatYrs] array of rate underlying maturities in years
 * 	(length "numRates"),
 * <br>[eigVect] on exit, contains array of eigenvectors.
 *	All but the N first (where N is the number of factors of the
 *	vnfmData) should be zero.
 * 	(must be allcated before call of length "numRates"),
 * <br>[eigVal] on exit, contains array of normalized eigenvalues
 * 	(must be allcated before call of size "numRates"*"numRates"),
 * <br>
 * Returns 0 iff OK.
 */

int
VnfmAvgQBOrthFactors(
	VnfmData *that,		/* (I) model parameters */
	double tStart,		/* (I) start of observation period (yrs) */
	double tEnd,		/* (I) end of observation period (yrs) */

	int numRates,		/* (I) # of maturities  */
	double *rateExpYrs,	/* (I) array of expirations (yrs from tEnd) */
	int *rateFreq,		/* (I) array of rate frequencies */
	double *rateMatYrs,	/* (I) array of maturities (yrs) */

	double **eigVect,	/* (O) array of eigen vectors */
	double *eigVal)		/* (O) array of eigen values */
{
static	char	routine[] = "VnfmAvgQBOrthFactors";
	int	status = FAILURE;

	double	*rateVol = NULL;
	double	**rateCov = NULL;
	double	rateCorr;
	int	idxR1, idxR2;

	/* (1) compute covariance matrix
	 */
	if (IS_ALMOST_ZERO(tStart - tEnd)) {
		tEnd = tStart + 1e0/365e0;	/* 1Day vol */
		for (idxR1=0; idxR1<numRates; idxR1++) {
			rateExpYrs[idxR1] = MAX(rateExpYrs[idxR1], tEnd);
		}
	}


	rateVol = DrlDoubleVectAlloc(0, numRates-1);
	ASSERT_OR_DONE(rateVol != NULL);
	rateCov = DrlDoubleMatrAlloc(0, numRates-1, 0, numRates-1);
	ASSERT_OR_DONE(rateCov != NULL);

	for (idxR1=0; idxR1<numRates; idxR1++) {

		IF_FAILED_DONE(VnfmAvgQBVol(
			that,
			tStart,
			tEnd,
			tEnd + rateExpYrs[idxR1],
			rateMatYrs[idxR1],
			rateFreq[idxR1],
			LOGVOL,
			&rateVol[idxR1]));
	}
	for (idxR1=0; idxR1<numRates; idxR1++)
	for (idxR2=0; idxR2<=idxR1  ; idxR2++) {

		IF_FAILED_DONE(	VnfmAvgQBCorr2(
			that,
			tStart,
			tEnd,

			tEnd + rateExpYrs[idxR1],
			rateMatYrs[idxR1],
			rateFreq[idxR1],

			tEnd + rateExpYrs[idxR2],
			rateMatYrs[idxR2],
			rateFreq[idxR2],

			&rateCorr));

		rateCov[idxR1][idxR2] = 
			rateVol[idxR1] * rateVol[idxR2] * rateCorr;
		rateCov[idxR2][idxR1] = rateCov[idxR1][idxR2];
	}


	/* (2) Diagonalize matrix
	 */
	IF_FAILED_DONE(DrlMatrixRealEigenVect(
		numRates,
		rateCov,
		eigVal,
		eigVect));


	/* OK ! */
	status = SUCCESS;
done:
	DrlDoubleVectFree(rateVol, 0, numRates-1);
	DrlDoubleMatrFree(rateCov, 0, numRates-1, 0, numRates-1);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}








/*f--------------------------------------------------------------
 * Compute average volatility of O/N rate.
 *                                                             
 * <br><br>
 * Generates the vol of average overnigh rate over specified period.
 * On entry, "that" contains the model parameters.
 * "startDate" is the start date of the averaging period, 
 * "endDate" is the end date of the averaging period, 
 * "resetDate" is the reset date of the option, 
 * Returns 0 if OK.
 */

int
VnfmAvgVolONRate(
	VnfmData *that,		/* (I) model parameters */
	TDate	startDate,	/* (I) start date of period */
	TDate	endDate,	/* (I) end date of period */
	TDate	resetDate,	/* (I) reset date of option */
	double *avgVol)		/* (O) vol of average overnight rate */
{
static	char		routine[] = "VnfmAvgVolONRatee";
	double		S, 	/* reset time point */
			T1, 	/* start date time point */
			T2, 	/* end date time point */
			avgMat = 0.;

	int		status = FAILURE,
			i, jIdx, j, nDim = NDIM;

	double		vol, avgVol1, avgVol2;

	double		x,
			i00[VNFM_NDIMMAX*(VNFM_NDIMMAX+1)/2];

	/* compute coefficients (you may have forgotten it before ..) */
	if (VnfmComputeCoeff(that) != 0)
			goto done;



	/* Set time line */
	GtoDayCountFraction(REFDATE, resetDate, GTO_ACT_365F, &S);

	GtoDayCountFraction(REFDATE, startDate, GTO_ACT_365F, &T1);
	GtoDayCountFraction(REFDATE, endDate,   GTO_ACT_365F, &T2);

	S  = MAX(S, 1e0/365e0); /* to avoid S=0*/
	T2 = MAX(T2, 1e0/365e0); /* to avoid T2=0*/

	avgMat = T2 - T1;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	    DrlFPrintf(vnfmFpLog, "%s: avgMat=%lf\n",\
		routine, avgMat));
#endif

	if (T2 < S) {
	    GtoErrMsg("%s: reset date (%s) > end date of avg period (%s)\n",
		routine, 
		GtoFormatDate(resetDate),
		GtoFormatDate(endDate));
	    goto done;
	}

	/* Total vol consists of two parts:
	 * 1. vol between 0 and MIN(T1, S), where rates are
	 * all forward. 
	 */
	if (IS_ALMOST_ZERO(MIN(T1, S)))
	    avgVol1 = 0.0;
	else {
	    if (VnfmAvgQBVol(that, 0e0, MIN(T1, S), T1, avgMat, 0, LOGVOL, &vol)
			!= SUCCESS) goto done;
	    avgVol1 = vol*vol*MIN(T1, S);
	}

	/* 2. compute vol between T1 and S, where rate
	 * starting to reset.
	 */
	avgVol2 = 0e0;

	if (S > T1)
	{
	    /* compute the average I integrals */
	    VnfmI(that, T1, S, T2, i00);

	    /* compute the average G coefficients */
	    if (VnfmG(that, T1, T2, &x) != SUCCESS)
		return(FAILURE);


	    /* vol computation */
	    for (i=0; i<=nDim-1; i++)
	    for (j=i; j<=nDim-1; j++) {
	    	jIdx = JIDX(i, j);
	    	avgVol2 += (i == j ? 1e0 : 2e0) * x * x * i00[jIdx];
	    }
	}


	/* $$$ Return NEGATIVE volatility if the covariance
	 * matrix is negative.
	 */
	*avgVol = SQRT_SIGN( (avgVol1 + avgVol2) / S );


	/*
	 * We are done
	 */
	status = SUCCESS;
done:
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}
	return(status);
}

