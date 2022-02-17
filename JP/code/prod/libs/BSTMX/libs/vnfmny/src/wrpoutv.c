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

#include "date_sup.h"
#include "convert.h"
#include "yearfrac.h"

#include "drlmem.h"		/* DrlDoubleMatrAlloc() */
#include "drlvtype.h"		/* DrlLilStructGet() */
#include "drlts.h"		/* DrlTCurveWrapRead() */
#include "drlio.h"		/* DrlFPrintf() */


#define	_vnfm_SOURCE
#include "vnfmanly.h"
#include "vnfmwrap.h"

#if defined(_WINDLL) 
# undef __DEBUG__
#endif

/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_GENVOLFWD.
 *                                                             
 * <br><br>
 * Wrapper for arbitrary rate volatility computation
 * given input factor mean-reversions, weights and spot volatilities.
 * <br>
 * <br>[refDateL] zero coupon curve value date.
 * <br>[zcDateL] zero coupon dates.
 * <br>[zcRateL] zero coupon rates.
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
 * and factors in column.
 * <br>[rhoL] range of correlations, indexed by dates in in row
 * and by factor in columns in the following way:
 * the correlations $(\rho_{ij}(t))_{1&lt;= i<j &lt;= N}$
 * between factors $i$ and $j$ is stored in the
 * column of index $k$ (where $k=0,\dots\,N-1$) where
 * $k = {N(N-1)/ 2} - {(N-i)(N-i-1)/ 2} + (j-i-1)$
 * with  $0&lt;= k &lt;= {N(N-1)/ 2}-1$.
 * %
 * <br>[obsStartDatesL] array of rate observation start dates
 * (must be of the same length as <i> obsStartDatesL</i>).
 * <br>[obsEndDatesL] array of rate observation end dates
 * (must be of the same length as <i> obsStartDatesL</i>).
 * <br>[resetDatesL] array of rate reset dates
 * (must be of the same length as <i> obsStartDatesL</i>).
 * <br>[maturitiesL] array rate maturities in years
 * (must be of the same length as <i> obsStartDatesL</i>).
 * <br>[frequenciesL] array of rate frequencies (0,1,2,12)
 * (must be of the same length as <i> obsStartDatesL</i>).
 * <br>[volsL] output array of volatilities
 * (must be of the same length as <i> obsStartDatesL</i>).
 * <br>
 */

DLL_EXPORT(int)
VnfmGenerateFwdVolL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

 	FloatL *backboneqL,	/*  4 'F' (I) [1] backbone q (0=LN, 0.5=N) */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of spot volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlation */

	TDateL *obsStartDatesL,	/* 10 'D' (I) observation start dates */
	TDateL *obsEndDatesL,	/* 11 'D' (I) observation end dates */
	TDateL *resetDatesL,	/* 12 'D' (I) rate reset dates */
	double *maturitiesL,	/* 13 'F' (I) rate maturities */
	long *frequenciesL,	/* 14 'L' (I) rate maturities */

	FloatL *volsL)		/* 15 'F' (O) zero coupon rates */
{
static	char		routine[] = "VnfmGenerateFwdVolL";
	int		status = FAILURE;
	VnfmData	*that = NULL;
	TCurve		*zcCurve = NULL;
	int		idx, nDates;
	double		S1, S2, T;

	VNFM_LOGOPEN

	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDateL, zcRateL)
		!= SUCCESS) goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,\
		 DRL_TCURVE_FMT_PERCENT));
#endif

	/* get model parameyters */
	if (VnfmWrapRead(&that,
			backboneqL, betaL, alphaL,
			dateL, sigmaL, rhoL,
			zcCurve)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(that, vnfmFpLog));
#endif

	/* Check array size */
	nDates = DrlLilVectActiveSize(DRL_TDATE_L, obsStartDatesL,
		"obsStartDatesL");
	if (nDates <= 0) {
		GtoErrMsg("%s: length array obsStartDatesL (%d) < 1\n",
			routine, nDates);
		goto done;
	}
#undef	CHECKARG
#define	CHECKARG(arg)	\
		{if (ARGSIZE(arg) != nDates) {\
		GtoErrMsg( "%s: len array %s (%d) != len array %s (%d)\n",\
		routine, #arg, ARGSIZE(arg), "obsStartDatesL", nDates);\
		goto done;}}

	CHECKARG(obsStartDatesL);
	CHECKARG(obsEndDatesL);
	CHECKARG(resetDatesL);
	CHECKARG(maturitiesL);
	CHECKARG(frequenciesL);
	CHECKARG(volsL);

#undef	CHECKARG

	/* Compute volatilities */
#undef	DCC
#define	DCC(d1, d2, t) \
		{if (GtoDayCountFraction((d1), (d2), GTO_ACT_365F, (t)) \
		!= SUCCESS) goto done;}

	if (VnfmComputeCoeff(that) != SUCCESS)
		goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
		DrlFPrintf(vnfmFpLog, "%s: nDates: %d\n", routine, nDates))
#endif

	for (idx=0; idx<=nDates-1; idx++) {
		DCC(REFDATE, obsStartDatesL[idx+1], &S1);
		DCC(REFDATE, obsEndDatesL[idx+1],   &S2);
		DCC(REFDATE, resetDatesL[idx+1],    &T);


		if (VnfmAvgQBVol(
			that,
			S1,
			S2,
			T,
			maturitiesL[idx+1],
			(int) frequenciesL[idx+1],
			LOGVOL,
			&volsL[idx+1])  != SUCCESS)
				goto done;


#ifndef	NO_LOGGING
		GTO_IF_LOGGING(DrlFPrintf(vnfmFpLog,\
			"  %10s (%8.4f) -> %10s (%8.4f) -"  \
			" %10s (%8.4f) %8.4f %2d : %8.6f %%\n",  \
			GtoFormatDate(obsStartDatesL[idx+1]), S1, \
			GtoFormatDate(obsEndDatesL[idx+1]),   S2, \
			GtoFormatDate(resetDatesL[idx+1]),    T, \
			maturitiesL[idx+1], (int) frequenciesL[idx+1], \
			volsL[idx+1]))
#endif
	}


	/* made it through OK */
	status = SUCCESS;
done:
	GtoFreeTCurve(zcCurve);
	VnfmFree(that);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
}


/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_GENBVOL.
 *                                                             
 * <br><br>
 * Wrapper for arbitrary volatility curve generation given
 * input factor mean-reversions, weights and spot volatilities.
 * <br>
 * <br>[refDateL] zero coupon curve value date.
 * <br>[zcDateL] zero coupon dates.
 * <br>[zcRateL] zero coupon rates.
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
 * and factors in column.
 * <br>[rhoL] range of correlations, indexed by dates in in row
 * and by factor in columns in the following way:
 * the correlations $(\rho_{ij}(t))_{1&lt;= i<j &lt;= N}$
 * between factors $i$ and $j$ is stored in the
 * column of index $k$ (where $k=0,\dots\,N-1$) where
 * $k = {N(N-1)/ 2} - {(N-i)(N-i-1)/ 2} + (j-i-1)$
 * with  $0&lt;= k &lt;= {N(N-1)/ 2}-1$.
 * %
 * <br>[bvMatL] base volatility underlying rate  maturity (in yrs)
 * for the output <i> bvRatesL</i>.
 * <br>[volFreqL] volatility frequency (0,1,2,4,12).
 * <br>[volTypeL] volatility curve type (0 for base, 1 for forward).
 * <br>[bvRefDateL] starting reference date for the volatility curve to be
 * generated (if this date is not today, the a forward volatility curve
 * is generated).
 * <br>[bvDatesL] base volatility expiration dates
 * for the output <i> bvRatesL</i>.
 * <br>[bvRatesL] output base volatility curve
 * specified by previous arguments.
 * <br>
 */

DLL_EXPORT(int)
VnfmGenerateVolCurveL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

 	FloatL *backboneqL,	/*  4 'F' (I) [1] backbone q (0=LN, 0.5=N) */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of spot volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlation */

	TDateIntervalL *bvMatL,	/* 10 'F' (I) vol maturity */
	IntL *volFreqL,		/* 11 'L' (I) vol frequency */
	IntL *volTypeL,		/* 12 'L' (I) 0=base, 1=fwd */
	TDateL *bvRefDateL,	/* 13 'D' (I) vol ref date */
	TDateL *bvDatesL,	/* 14 'D' (I) array of vol dates */
	FloatL *bvRatesL)	/* 15 'F' (O) array of vol */
{
static	char		routine[] = "VnfmGenerateBaseVolL";
	int		status = FAILURE;
	VnfmData	*theModel = NULL;
	TCurve		*zcCurve = NULL,
			*bvCurve = NULL;
	int		n, volType, volFreq,
			nBvDates;
	TDateInterval	bvMat;
	TDate		bvRefDate, *bvDates;


	VNFM_LOGOPEN

	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDateL, zcRateL)
		!= SUCCESS) goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,\
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

	/* get basevol interval */
	if (DrlLilStructGet(0,
	    DRL_LILVAR_L, "bvRefDate", DRL_TDATE_L,
		(void*) bvRefDateL, (void*) &bvRefDate,
	    DRL_LILVAR_L, "bvMat", DRL_TDATEINTERVAL_L,
		(void*) bvMatL, (void*) &bvMat,
	    DRL_LILVAR_L, "volType", DRL_INT_L,
		(void*) volTypeL, (void*) &volType,
	    DRL_LILVAR_L, "volFreq", DRL_INT_L,
		(void*) volFreqL, (void*) &volFreq,
	    0) != SUCCESS) goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "Base Vol Maturity: %s\n",\
		GtoFormatDateInterval(&bvMat)));
#endif

	/* get number of nonzero dates */
	nBvDates = DrlLilVectActiveSize(DRL_TDATE_L, bvDatesL, "bvDates");
	if (nBvDates <= 0) goto done;
	bvDates = bvDatesL + 1;


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "nBvDates: %d\n", nBvDates);\
	DrlFPrintf(vnfmFpLog, "Base Vol TDates:  %d\n", nBvDates);\
	for (n=0; n<=nBvDates-1; n++) {\
		DrlFPrintf(vnfmFpLog, "\t[%3d/%3d] %s\n",\
			n, nBvDates, GtoFormatDate(bvDates[n]));\
	});
#endif


	/* Generate base voltility curve */
	if (VnfmGenerateVolTCurve(
			theModel,
			bvMat, volFreq, volType,
			bvRefDate, nBvDates, bvDates,
			&bvCurve)
		!= SUCCESS) goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlTCurveFpWrite(bvCurve, vnfmFpLog, DRL_TCURVE_FMT_PERCENT));
#endif

	if (DrlTCurveWrapWrite(bvCurve, NULL, NULL, bvRatesL)
		!= SUCCESS) goto done;


	/* made it through OK */
	status = SUCCESS;
done:
	GtoFreeTCurve(zcCurve);
	GtoFreeTCurve(bvCurve);
	VnfmFree(theModel);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
}


/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_GENSWMAT.
 *                                                             
 * <br><br>
 * Wrapper for swaption matrix generation.
 * Assumes all model parameters known (mean reversions, weights,
 * spot volatilties, correlations).
 * Wrapper for arbitrary volatility curve generation given
 * input factor mean-reversions, weights and spot volatilities.
 * <br>
 * <br>[refDateL] zero coupon curve value date.
 * <br>[zcDateL] zero coupon dates.
 * <br>[zcRateL] zero coupon rates.
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
 * and factors in column.
 * <br>[rhoL] range of correlations, indexed by dates in in row
 * and by factor in columns in the following way:
 * the correlations $(\rho_{ij}(t))_{1&lt;= i<j &lt;= N}$
 * between factors $i$ and $j$ is stored in the
 * column of index $k$ (where $k=0,\dots\,N-1$) where
 * $k = {N(N-1)/ 2} - {(N-i)(N-i-1)/ 2} + (j-i-1)$
 * with  $0&lt;= k &lt;= {N(N-1)/ 2}-1$.
 * %
 * <br>[swRefDateL] starting reference date for the volatility matrix to be
 *  generated (if this date is not today, the a forward volatility curve
 * <br>[swTypeL] array of 2 elements:\\
 * 	(1) swaption matrix output type (0 for vertical, 1 for diagonal)
 *		for the output <i> volL</i>.\\
 * 	(2) swaption matrix volatility frequency (0,1,2,4,12)
 *		for the output <i> volL</i>.
 * <br>[tMatL] array of maturity intervals (in yrs)
 *		for the output <i> volL</i>.
 * <br>[tExpL] array of expration intervals (in yrs)
 *		for the output <i> volL</i>.
 * <br>[swVolL] output swaption matrix
 * specified by arguments <i> swTypeL</i>, <i> tMatL</i> and <i> tExpL</i>.
 * <br>
 */

DLL_EXPORT(int)
VnfmGenerateSwaptionMatrixL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

 	FloatL *backboneqL,	/*  4 'F' (I) [1] backbone q (0=LN, 0.5=N) */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of spot volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlation */

	TDateL *swRefDateL,	/* 10 'D' (I) vol ref date */
	IntL *swTypeL,		/* 11 'L' (I) array of matrix param [0..2]: */
				/*        [0] matr type (0=vertical, 1=diag) */
				/*        [1] vol frequency (1,2,4,12) */
				/*        [2] NOT USED (pass 0) */
	TDateIntervalL* tMatL,	/* 12 'F' (I) array of mat intervals */
	TDateIntervalL* tExpL,	/* 13 'F' (I) array of exp intervals */
	FloatL *volL)		/* 14 'F' (O) matrix of swaption vols */
{
static	char		routine[] = "VnfmGenerateSwaptionMatrixL";
	VnfmData	*theModel = NULL;
	TCurve		*zcCurve = NULL;
	int		status = FAILURE;

#undef	NMAX
#define	NMAX		64

	int		nExp = 0, nMat = 0,
			n,
			*swType = NULL;
	TDate		swRefDate;
	TDateInterval	*lExp = NULL,
			*lMat = NULL;
	double		**vol = NULL;




	VNFM_LOGOPEN


	/*
	 * Wrap arguments
	 */

	/* read the zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDateL, zcRateL)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,\
		DRL_TCURVE_FMT_PERCENT));
#endif

	/* read the model parameters */
	if (VnfmWrapRead(&theModel,
			backboneqL, betaL, alphaL,
			dateL, sigmaL, rhoL,
			zcCurve)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(theModel, vnfmFpLog));
#endif

	/* read swaption maturities and expirations */
	if (DrlLilStructGet(0,
	    DRL_LILVAR_L,
	    "swaptionRefDate", DRL_TDATE_L,
			(void*) swRefDateL, (void*) &swRefDate,
	    DRL_LILVECT_L, &nExp, 1, NMAX,
	    "swaptionExp", DRL_TDATEINTERVAL_L, TRUE,
			(void*)tExpL, (void*)&lExp,
	    DRL_LILVECT_L, &nMat, 1, NMAX,
	    "swaptionMat", DRL_TDATEINTERVAL_L, TRUE,
			(void*)tMatL, (void*)&lMat,
	    DRL_LILVECT_L, &n, 2, 2,
	    "swaptionTyp", DRL_INT_L, TRUE, (void*)swTypeL, (void*)&swType,
	    0) != SUCCESS) goto done;

	if ((vol = DrlDoubleMatrAlloc(0, NMAX, 0, NMAX)) == NULL)
		goto done;



#ifndef	NO_LOGGING
	GTO_IF_LOGGING( {int	i;\
	for (i=0; i<=nExp-1; i++)\
		DrlFPrintf(vnfmFpLog, "\tnExp=%2d\t%s\n",\
			i, GtoFormatDateInterval(&lExp[i]));\
	for (i=0; i<=nMat-1; i++)\
		DrlFPrintf(vnfmFpLog, "\tnMat=%2d\t%s\n",\
			i, GtoFormatDateInterval(&lMat[i]));\
	});
#endif


	/*
	 * Generate the swaption matrix
	 */
	if (VnfmComputeCoeff(theModel)
		!= SUCCESS) goto done;

	if (VnfmSwaptionVolMatrix(
			theModel,
			swType[0], swType[1], FALSE,
			nExp, lExp,
			nMat, lMat,
			LOGVOL,
			vol)
		!= SUCCESS) goto done;

#ifdef	_SKIP
	GTO_IF_LOGGING(\
	DrlFPrintfDrlDoubleMatrix(vnfmFpLog, "\t%lf", vol, nExp, nMat));
#endif

	/*
	 *
	 */
	if (DrlLilStructPut(0,
	    DRL_LILMATR_L, nExp, nMat,
		"swaptionVol", DRL_FLOAT_L, (void*) volL, (void*) vol,
	    0) != SUCCESS) goto done;

	/* made it through OK */
	status = SUCCESS;
done:
	GtoFreeTCurve(zcCurve);
	VnfmFree(theModel);
	DrlVTypeVectFree(lMat, nMat, DRL_TDATEINTERVAL_T);
	DrlVTypeVectFree(lExp, nExp, DRL_TDATEINTERVAL_T);
	DrlVTypeVectFree(swType, 3, DRL_INT_T);
	DrlDoubleMatrFree(vol, 0, NMAX, 0, NMAX);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
#undef	NMAX
}


/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_GENACORR.
 *                                                             
 * <br><br>
 * Wrapper for average correlation matrix generation.
 * Assumes all model parameters known (mean reversions, weights,
 * spot volatilties, correlations).
 * <br>
 * <br>[refDateL] zero coupon curve value date.
 * <br>[zcDateL] zero coupon dates.
 * <br>[zcRateL] zero coupon rates.
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
 * and factors in column.
 * <br>[rhoL] range of correlations, indexed by dates in in row
 * and by factor in columns in the following way:
 * the correlations $(\rho_{ij}(t))_{1&lt;= i<j &lt;= N}$
 * between factors $i$ and $j$ is stored in the
 * column of index $k$ (where $k=0,\dots\,N-1$) where
 * $k = {N(N-1)/ 2} - {(N-i)(N-i-1)/ 2} + (j-i-1)$
 * with  $0&lt;= k &lt;= {N(N-1)/ 2}-1$.
 * %
 * <br>[tExpL] time to expiration of correlations to be computed.
 * <br>[tMatL] array of time to maturity of rates whose
 * correlation is to be computed.
 * <br>[corrL] output range (square matrix) of correlations.
 * <br>
 */

DLL_EXPORT(int)
VnfmGenerateCorrelationMatrixL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

 	FloatL *backboneqL,	/*  4 'F' (I) [1] backbone q (0=LN, 0.5=N) */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of spot volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlation */

	TDateIntervalL* tExpL,	/* 10 'F' (I) expiration time */
	TDateIntervalL* tMatL,	/* 11 'F' (I) array of mat times */
	FloatL *corrL)		/* 12 'F' (O) matrix of correlations */
{
static	char		routine[] = "VnfmGenerateCorrelationMatrixL";
	int		status = FAILURE;
	VnfmData	*theModel = NULL;
	TCurve		*zcCurve = NULL;

#undef	NMAX
#define	NMAX		64

	int		nMat;
	TDateInterval	lExp,
			*lMat = NULL;
	double		**corr = NULL;


	VNFM_LOGOPEN


	/*
	 * Wrap arguments
	 */

	/* read the zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDateL, zcRateL)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,\
		DRL_TCURVE_FMT_PERCENT));
#endif

	/* read the model parameters */
	if (VnfmWrapRead(&theModel,
			backboneqL, betaL, alphaL,
			dateL, sigmaL, rhoL,
			zcCurve)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(theModel, vnfmFpLog));
#endif

	/* read corr maturities and expirations */
	if (DrlLilStructGet(0,
	    DRL_LILVAR_L,
	    "corrExp", DRL_TDATEINTERVAL_L, (void*)tExpL, (void*)&lExp,
	    DRL_LILVECT_L, &nMat, 1, NMAX,
	    "corrMat", DRL_TDATEINTERVAL_L, TRUE, (void*)tMatL, (void*)&lMat,
	    0) != SUCCESS) goto done;

	if ((corr = DrlDoubleMatrAlloc(0, NMAX, 0, NMAX)) == NULL)
		goto done;


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	{int	i;\
	DrlFPrintf(vnfmFpLog, "\tnExp=%s\n", GtoFormatDateInterval(&lExp));\
	for (i=0; i<=nMat-1; i++)\
		DrlFPrintf(vnfmFpLog, "\tnMat=%2d\t%s\n",\
			i, GtoFormatDateInterval(&lMat[i]));
	});
#endif


	/* Generate corr matrix */
	if (VnfmComputeCoeff(theModel)
		!= SUCCESS) goto done;

	if (VnfmAvgCorrMatrix(
			theModel,
			lExp, nMat,
			lMat, corr)
		!= SUCCESS) goto done;

#ifdef	_SKIP
	GTO_IF_LOGGING(\
	drlFPrintDoubleMatr(vnfmFpLog, "\t%lf", corr, nMat, nMat));
#endif

	/* put corr in wrapper */
	if (DrlLilStructPut(0,
	    DRL_LILMATR_L, nMat, nMat,
		"corr", DRL_FLOAT_L, (void*) corrL, (void*) corr,
	    0) != SUCCESS) goto done;

	/* made it through OK */
	status = SUCCESS;
done:
	GtoFreeTCurve(zcCurve);
	VnfmFree(theModel);
	DrlVTypeVectFree(lMat, nMat, DRL_TDATEINTERVAL_T);
	DrlDoubleMatrFree(corr, 0, NMAX, 0, NMAX);

	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
#undef	NMAX
}


/*f--------------------------------------------------------------
 * Wrapper for forward/futures adjustment.
 *                                                             
 * <br><br>
 * Assumes all model parameters known (mean reversions, weights,
 * spot volatilties, correlations).
 * <br>
 * <br>[refDateL] zero coupon curve value date.
 * <br>[zcDateL] zero coupon dates.
 * <br>[zcRateL] zero coupon rates.
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
 * and factors in column.
 * <br>[rhoL] range of correlations, indexed by dates in in row
 * and by factor in columns in the following way:
 * the correlations $(\rho_{ij}(t))_{1&lt;= i<j &lt;= N}$
 * between factors $i$ and $j$ is stored in the
 * column of index $k$ (where $k=0,\dots\,N-1$) where
 * $k = {N(N-1)/ 2} - {(N-i)(N-i-1)/ 2} + (j-i-1)$
 * with  $0&lt;= k &lt;= {N(N-1)/ 2}-1$.
 * %
 * <br>[rateMatL] rate maturity (in yrs).
 * <br>[rateFreqL] rate frequency (0,1,2,4,12).
 * <br>[rateDayCountL]	rate day count (ACT/360, 30/360, etc).
 * <br>[adjTypeL] NOT USED (pass dummy variable).
 * <br>[resetDatesL] array of reset dates at which the adjustment
 * is to be computed.
 * <br>[outputL] output range index, indexed by dates in in row
 * and having 3 columns containing respecively the forward rate
 * (computed with no adjustment), the convexity adjustment
 * and the MTM adjustment.
 * <br>
 */

DLL_EXPORT(int)
VnfmGenerateFwdFutAdjustmentL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

 	FloatL *backboneqL,	/*  4 'F' (I) [1] backbone q (0=LN, 0.5=N) */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of spot volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlation */

	TDateIntervalL *rateMatL,/* 10 'F' (I) rate maturity */
	IntL *rateFreqL,	/* 11 'L' (I) rate frequency */
	CharBlockL *rateDayCountL,/* 12 'C' (I) rate day count */
	CharBlockL *adjTypeL,	/* 13 'C' (I) NOT USED */
	TDateL *resetDatesL,	/* 14 'D' (I) array of dates */
	FloatL *outputL)	/* 15 'F' (O) array of adjustment values */
{
static	char		routine[] = "VnfmGenerateFwdFutAdjusmentL";
	int		status = FAILURE;
	VnfmData	*theModel = NULL;
	TCurve		*zcCurve = NULL;
	int		rateFreq,
			nDates,
			n;
	TDateInterval	rateMat;
	long		rateDayCount;
	double		*fwdRate = NULL,
			*cvxAdj = NULL,
			*mtmAdj = NULL;

	/*
	 *
	 */

	WRAP_CHECK_VECTOR(refDateL);
	WRAP_CHECK_VECTOR(zcDateL);
	WRAP_CHECK_VECTOR(zcRateL);
	WRAP_CHECK_VECTOR(backboneqL);
	WRAP_CHECK_VECTOR(betaL);
	WRAP_CHECK_VECTOR(alphaL);
	WRAP_CHECK_VECTOR(dateL);
	WRAP_CHECK_VECTOR(sigmaL);
	WRAP_CHECK_VECTOR(rhoL);

	WRAP_CHECK_VECTOR(rateMatL);
	WRAP_CHECK_VECTOR(rateFreqL);
	WRAP_CHECK_VECTOR(rateDayCountL);
	WRAP_CHECK_VECTOR(resetDatesL);
	WRAP_CHECK_VECTOR(outputL);





	VNFM_LOGOPEN



	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDateL, zcRateL)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,\
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

	/* convert other arguments */
	if (GtoYearsToDateInterval(rateMatL[1], &rateMat) != SUCCESS)
		goto done;

	rateFreq = rateFreqL[1];

	if (GtoStringToDayCountConv(&rateDayCountL[WRAP_STR_IDX(1)],
		&rateDayCount) != SUCCESS) goto done;

	/* get number of nonzero dates */
	if ((nDates = DrlLilVectActiveSize(DRL_TDATE_L,
		resetDatesL, "resetDatesL")) <= 0)
			goto done;

	if (ARGSIZE(outputL) != 3*ARGSIZE(resetDatesL)) {
		GtoErrMsg("%s: output range size (%d) != "
		        "3 * reset dates range size (%d)\n",
			routine, ARGSIZE(outputL), ARGSIZE(resetDatesL));
		goto done;
	}
	if ((fwdRate = DrlDoubleVectAlloc(0, nDates-1)) == NULL)
		goto done;
	if ((cvxAdj  = DrlDoubleVectAlloc(0, nDates-1)) == NULL)
		goto done;
	if ((mtmAdj  = DrlDoubleVectAlloc(0, nDates-1)) == NULL)
		goto done;


#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlFPrintf(vnfmFpLog, "Rate Maturity: %4s\n",\
		GtoFormatDateInterval(&rateMat));\
	DrlFPrintf(vnfmFpLog, "Dates:  %d\n", nDates);\
	for (n=0; n<=nDates-1; n++) {\
		DrlFPrintf(vnfmFpLog, "\t[%3d/%3d] %10s\n",\
			n, nDates, GtoFormatDate(resetDatesL[n]));\
	});
#endif


	/* Compute adjustment */
	if (VnfmFwdFutAdjustment(
		theModel,
		zcCurve,
		rateMat,
		rateFreq,
		rateDayCount,
		nDates,
		&resetDatesL[1],
		fwdRate,
		cvxAdj,
		mtmAdj) != SUCCESS) goto done;


	/* wrap out */
	for (n=0; n<=nDates-1; n++) {
	    outputL[WRAP_MATR_IDX(ARGSIZE(resetDatesL), 3, n, 0)] =
			fwdRate[n];
	    outputL[WRAP_MATR_IDX(ARGSIZE(resetDatesL), 3, n, 1)] =
			cvxAdj[n];
	    outputL[WRAP_MATR_IDX(ARGSIZE(resetDatesL), 3, n, 2)] =
			mtmAdj[n];
	}


	/* made it through OK */
	status = SUCCESS;
done:
	DrlDoubleVectFree(fwdRate, 0, nDates-1);
	DrlDoubleVectFree(cvxAdj, 0, nDates-1);
	DrlDoubleVectFree(mtmAdj, 0, nDates-1);
	GtoFreeTCurve(zcCurve);
	VnfmFree(theModel);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return(status);
}




