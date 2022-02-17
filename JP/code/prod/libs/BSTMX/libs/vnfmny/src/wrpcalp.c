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

#include "drlts.h"
#include "drlvtype.h"

#define	_vnfm_SOURCE

#include "vnfmopca.h"
#include "vnfmwrap.h"

#if defined(_WINDLL) || !defined(TESTLIB)
# undef __DEBUG__
#endif


/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_GENBENCH.
 *                                                             
 * <br><br>
 * Wrapper for benchmarks volatilities computation.
 * Assumes all model parameters known (mean reversions, weights,
 * spot volatilties, correlations).
 * The argument <tt> optBenchL</tt> is an array of strings
 * defining volatility or correleation points.
 * The format recognized for the string is described below:
 * <ol>
 * <li> <b> Volatility Points</b> defined to be the volatility
 * between two dates (<tt> Start</tt> and <i> Exp</i>) of a forward rate
 * that resets at some time <tt> Reset</tt> in the future.
 * The format recognized is
 * <br><tt>
 * FLRATE_VOL(Start=interval,Exp=interval,Reset=interval,
 *            Mat=interval,Freq=int)
 * <br></tt>
 * where all intervals are from the base date (today's date).
 * For Example, to define the foward volatility between 1M and  3M
 * from now of a 10Y semiannual rate that reset in 1Y from now,
 * pass the string
 * <br><tt>
 * FLRATE_VOL(Start=1M,Exp=3M,Reset=1A,Mat=10A,Freq=2)
 * <br></tt>
 * <li> <b> Correlation Points</b> defined to be the correlation
 * of two rates between today and a date. The format is similar
 * that the volatility point:
 * all intervals are from the base date (today's date).
 * <br><tt>
 * FLRATE_CORR(Exp=interval,Mat1=interval,Mat2=interval,
 *            Freq1=int,Freq2=int)
 * <br></tt>
 * For Example, to define the 3M correlation between
 * the 3M and 10Y rates, pass the string
 * <br><tt>
 * FLRATE_CORR(Exp=3M,Mat1=3M,Freq1=0,Mat2=10A,Freq2=2)
 * <br></tt>
 * </ol>
 */

DLL_EXPORT(int)
VnfmGenerateBenchmarkVolL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

 	FloatL *backboneqL,	/*  4 'F' (I) 0=lognormal, 0.5=normal */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of spot volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlation */

	CharBlockL *optBenchL,	/* 10 'C' (I) string describing benchmarks */
	FloatL *optModL)	/* 11 'F' (O) benchmark values */
{
static	char		routine[] = "VnfmGenerateBenchmarkVolL";
	int		status = FAILURE;
	VnfmData	*vnfmData = NULL;
	TCurve		*zcCurve = NULL;
	int		nOptBench, k;
	TVolBenchmark	optBench[256];


	VNFM_LOGOPEN

	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDateL, zcRateL)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(\
	DrlTCurveFpWrite(zcCurve, vnfmFpLog, DRL_TCURVE_FMT_PERCENT));
#endif

	/* get model parameyters */
	if (VnfmWrapRead(&vnfmData,
			backboneqL, betaL, alphaL,
			dateL, sigmaL, rhoL,
			zcCurve)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog));
#endif


	/* get benchmarks */
	if (DrlTVolBenchmarkWrapRead(
		optBench, &nOptBench, NULL, optBenchL)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTVolBenchmarkFpWrite(optBench, nOptBench,\
		vnfmFpLog, NULL));
#endif

	ASSERT_OR_DONE(ARGSIZE(optModL) >= nOptBench);


	/* Generate corr matrix */
	if (VnfmComputeCoeff(vnfmData)
		!= SUCCESS) goto done;

	/* compute vols */
	for (k=0; k<=nOptBench-1; k++) {
	    if (VnfmTVolBenchmarkValue(vnfmData,
			&optBench[k], &optModL[k+1]) != SUCCESS)
		goto done;
	}


	/* made it through OK */
	status = SUCCESS;
done:
	GtoFreeTCurve(zcCurve);
	VnfmFree(vnfmData);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
}


/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CALPAR_SHORT.
 *                                                             
 * <br><br>
 * Wrapper for NF parameter calibration assuming constant spot volatilities.
 * Given benchmarks U_1,...,U_N and C_1,...,C_M
 * (a benchmark being as either a volatility or a correlation point,
 * see <tt> VnfmGenerateBenchmarkVolL</tt>), 
 * performs the optimization of the two factor parameters
 * (beta_i,alpha_i,rho_ij) (assuming a <i> constant</i>
 * spot volatility) to minimize the weighted sum of
 * deviation from the market values of the benchmarks U_i
 * <br><br> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
 * min_{beta_i,alpha_i,rho_ij} sum_{i=1,N} w_i
 *            (U_i(market) - U_i(model))^2
 * <br><br>
 * under the constraints that the model values of the benchmarks C_i
 * are within prescribed bounds.
 * <br><br> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
 * for all j, C_i(min) < C_i(model) < C_i(max)
 * <br><br>
 * where w_i, U_i(market),
 * C_i(min),  C_i(max) are specified by the user.
 * Calls <tt> VnfmCalibParamShortTerm</tt>.
 * <br><br><b>Arguments Details</b><br><br><ul>
 * <li> <b> refDateL :</b> zero curve value date.
 * <li> <b> zcDatesL :</b> zero curve dates.
 * <li> <b> zcRatesL :</b> zero curve rates.
 * <li> <b> datesL :</b> array of dates for calibration timeline.
 * <li> <b> paOptFlagL :</b> array of length $n(n+3)/2$ where $n$ is the number
 * of factors corresponding to the model parameters in the order
 * $$(\beta_1,...,\beta_n,\alpha_1,...,\alpha_n,
 *     \rho_{11},...,\rho_{1n},\rho_{21},...,\rho_{n-1,n})$$
 * containing the flag 1 if
 * the corresponding parameter is to be optimized, 0 if not.
 * <li> <b> paMinL :</b> array of lower bounds on the 2F parameters
 *               (same length as <tt> paOptFlagL</tt>).
 * <li> <b> paMaxL :</b> array of upper bounds on the 2F parameters
 *               (same length as <tt> paOptFlagL</tt>).
 * <li> <b> paMidL :</b> array of initial values of 2F parameters (or value used
 *      if the  parameter is not optimized
 *               (same length as <tt> paOptFlagL</tt>).
 * <li> <b> nOptBenchL :</b> number of optimized benchmarks.
 * <li> <b> optBenchL :</b> array of optimized benchmarks definition as a string
 *   See <tt> VnfmGenerateBenchmarkVolL</tt> for a description of the format.
 *   The length of the array can be more than <tt> nOptBenchL</tt>,
 *   but only the first <tt> nOptBenchL</tt> are used.
 * <li> <b> optWeiL :</b> array of optimized benchmarks weights in the optimization.
 *   (same length as <tt> optBenchL</tt>).
 * <li> <b> optMidL :</b> array of optimized benchmarks market values
 *   (same length as <tt> optBenchL</tt>).
 * <li> <b> nConstBenchL :</b> number of constraints benchmarks
 * <li> <b> constBenchL :</b> array of constraint benchmark 
 *   See <tt> VnfmGenerateBenchmarkVolL</tt> for a description of the format.
 *   The length of the array can be more than <tt> nConstBenchL</tt>,
 *   but only the first <tt> nConstBenchL</tt> are used.
 * <li> <b> constMinL :</b> array of constraint min value 
 *   (same length as <tt> constBenchL</tt>).
 * <li> <b> constMaxL :</b> array of constraint max value 
 *   (same length as <tt> constBenchL</tt>).
 * <li> <b> floatScalarsL :</b>	array of length 1.\\
 * Element 1: back bone $q$ (0 for log-normal, 0.5 for normal).
 * <li> <b> paOptL :</b> output optimal parameters,
 *   (same length as <tt> paOptFlagL</tt>).
 * <li> <b> optModL :</b> output optimized benchmarks model value,
 *   (same length as <tt> optBenchL</tt>).
 * <li> <b> constModL :</b> output constraint benchmarks model value,
 *   (same length as <tt> constBenchL</tt>).
 * </ul>
 *
 */

DLL_EXPORT(int)
VnfmCalibParamShortTermL(
	TDateL *refDateL,	/*  1 'D' (I) zero curve value date */
	TDateL *zcDatesL,	/*  2 'D' (I) zero curve dates */
	FloatL *zcRatesL,	/*  3 'F' (I) zero curve rates */

	TDateL *datesL,		/*  4 'D' (I) array of dates */
				/*        Input Model: */
	IntL *paOptFlagL,	/*  5 'L' (I) optimize flag */
	FloatL *paMinL,		/*  6 'F' (I) min value */
	FloatL *paMaxL,		/*  7 'F' (I) max value */
	FloatL *paMidL,		/*  8 'F' (I) initial value */
	
	IntL *nOptBenchL,	/*  9 'L' (I) # of optimized benchmarks */
	CharBlockL *optBenchL,	/* 10 'C' (I) optim. benchmarks */
	FloatL *optWeiL,	/* 11 'F' (I) optim. market weight */
	FloatL *optMidL,	/* 12 'F' (I) optim. market value */

	IntL *nConstBenchL,	/* 13 'L' (I) # of constraints benchmarks */
	CharBlockL *constBenchL,	/* 14 'C' (I) constr. benchmark */
	FloatL *constMinL,	/* 15 'F' (I) constr. min value */
	FloatL *constMaxL,	/* 16 'F' (I) constr. max value */

 	FloatL *floatScalarsL,	/* 17 'F' (I) [1] = backbone q (0=LN,0.5=N) */

	FloatL *paOptL,		/* 18 'F' (O) optimal parameters */
	FloatL *optModL,	/* 19 'F' (O) optim. model value */
	FloatL *constModL)	/* 20 'F' (O) constr. model value */

{
static	char		routine[] = "VnfmCalibParamShortTerm";
	int		status = FAILURE;
	VnfmData	*vnfmData = NULL;
	TCurve		*zcCurve = NULL;
	int		n;		/* dummy */
	int		nF;		/* number of factors*/
			/* optimization benchmarks */
	int		nOptBench;
	TVolBenchmark	optBench[VNFM_NMAX_BENCH];
			/* constraint benchmarks */
	int		nConstBench;
	TVolBenchmark	constBench[VNFM_NMAX_BENCH];



	/* Addin argumrnts Logging */
	if (GtoLoggingGet() > 0) {
	    DrlLilVectLoggingFile("wrapper.log", "w", "TF_CALPAR_SHORT",
		DRL_TDATE_L, (void*) refDateL,      "ZC_REFDATE",
		DRL_TDATE_L, (void*) zcDatesL,      "ZC_DATES",
		DRL_FLOAT_L, (void*) zcRatesL,      "ZC_RATES",

		DRL_TDATE_L, (void*) datesL,        "IN_DATES",

		DRL_LONG_L,  (void*) paOptFlagL,    "OPT_FLAG",
		DRL_FLOAT_L, (void*) paMinL,        "OPT_MIN",
		DRL_FLOAT_L, (void*) paMaxL,        "OPT_MAX",
		DRL_FLOAT_L, (void*) paMidL,        "OPT_MID",

		DRL_LONG_L,  (void*) nOptBenchL,    "OBJ_N",
		DRL_CHAR_BLOCK_L,(void*) optBenchL, "OBJ_BENCH",
		DRL_FLOAT_L, (void*) optWeiL,       "OBJ_WEI",
		DRL_FLOAT_L, (void*) optMidL,       "OBJ_MID",

		DRL_LONG_L,  (void*) nConstBenchL,  "CON_N",
		DRL_CHAR_BLOCK_L,(void*) constBenchL,"CON_BENCH",
		DRL_FLOAT_L, (void*) constMinL,      "CON_MIN",
		DRL_FLOAT_L, (void*) constMaxL,      "CON_MAX",

		DRL_FLOAT_L, (void*) floatScalarsL,  "FLOAT_SCALARS",

		DRL_FLOAT_L, (void*) paOptL,         "OPT_OPT",
		DRL_FLOAT_L, (void*) optModL,        "OBJ_MOD",
		DRL_FLOAT_L, (void*) constModL,      "CON_MOD",
		0);
	}


	/* open log */
	VNFM_LOGOPEN


	/* get zero curve */
	if (DrlTCurveWrapRead(&zcCurve, refDateL, zcDatesL, zcRatesL)
		!= SUCCESS) goto done;
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,\
		DRL_TCURVE_FMT_PERCENT));
#endif
	/* check array length */
	WRAP_CHECK_VECTOR_LEN(floatScalarsL, 1);


	/* get number of factors from length of arrays */
	nF = (int) floor(sqrt(2e0*ARGSIZE(paOptFlagL)+2.25e0) - 1.5 +0.5);
#undef	CHECK_LEN
#define	CHECK_LEN(arg)    {if (ARGSIZE(arg) != nF*(nF+3)/2) {\
		GtoErrMsg("%s: argument %s has length %d incompatible"\
		" with num fact %d.\n", routine,  #arg, (int) ARGSIZE(arg),\
		 nF); goto done;}}

	CHECK_LEN(paOptFlagL);
	CHECK_LEN(paMinL);
	CHECK_LEN(paMaxL);
	CHECK_LEN(paMidL);
	CHECK_LEN(paOptL);

#undef	CHECK_LEN

	/* create new data structure */
	vnfmData = VnfmNewTimeLine(
			nF,
			floatScalarsL[1], /* back bone */
			zcCurve->fBaseDate,
			ARGSIZE(datesL),
			datesL+1,
			NULL,	/* used for timeline only (can be NULL) */
			zcCurve);
	if (vnfmData == NULL) goto done;


#ifdef	_SKIP
	GTO_IF_LOGGING(\
	DrlFPrintStruct(vnfmFpLog, NULL, '\t',\
	    DRL_CARRAY_T, "paraopt", n, 5,\
		INT_T, (void*) (paOptFlagL+1),\
		DRL_DOUBLE_T, (void*) (paOptL+1),\
		DRL_DOUBLE_T, (void*) (paMinL+1),\
		DRL_DOUBLE_T, (void*) (paMaxL+1),\
		DRL_DOUBLE_T, (void*) (paMidL+1),\
	0));
	GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog));
#endif

	/* get benchmarks vols */
	/* input vol benchmarks */
	if (DrlTVolBenchmarkWrapRead(optBench, &nOptBench,
		nOptBenchL, optBenchL) != SUCCESS)
			goto done;

	if (DrlLilStructGet(0,
	    DRL_LILVECTARRAY_L, &n, nOptBench, 128, 3,
		"optWei", DRL_FLOAT_L, FALSE, (void*) optWeiL, NULL,
		"optMid", DRL_FLOAT_L, FALSE, (void*) optMidL, NULL,
		"optMod", DRL_FLOAT_L, FALSE, (void*) optModL, NULL,
	    0) != SUCCESS) goto done;


	ASSERT_OR_DONE(ARGSIZE(optWeiL) >= nOptBench);
	ASSERT_OR_DONE(ARGSIZE(optMidL) >= nOptBench);
	ASSERT_OR_DONE(ARGSIZE(optModL) >= nOptBench);

	/* correlation benchmarks */
	if (DrlTVolBenchmarkWrapRead( constBench, &nConstBench,
		nConstBenchL, constBenchL) != SUCCESS)
			goto done;

	if (DrlLilStructGet(0,
	    DRL_LILVECTARRAY_L, &n, nConstBench, 128, 3,
		"constMax", DRL_FLOAT_L, FALSE, (void*) constMaxL, NULL,
		"constMin", DRL_FLOAT_L, FALSE, (void*) constMinL, NULL,
		"constMod", DRL_FLOAT_L, FALSE, (void*) constModL, NULL,
	    0) != SUCCESS) goto done;

	ASSERT_OR_DONE(ARGSIZE(constMaxL) >= nConstBench);
	ASSERT_OR_DONE(ARGSIZE(constMinL) >= nConstBench);
	ASSERT_OR_DONE(ARGSIZE(constModL) >= nConstBench);


	/* Perform the calibration */
	if (VnfmCalibParamShortTerm(
		vnfmData,
		"F",	/* flat spot vol optimization */
		nOptBench, optBench, optWeiL+1, optMidL+1, optModL+1,
		nConstBench, constBench, constMinL+1, constMaxL+1, constModL+1,
		paOptFlagL+1,
		paOptL+1,
		paMinL+1,
		paMaxL+1,
		paMidL+1)
		!= SUCCESS) goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog));
#endif


	/* made it through OK */
	status = SUCCESS;
done:
	GtoFreeTCurve(zcCurve);
	VnfmFree(vnfmData);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
}




/*f--------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CALPAR_SQ1V.
 *                                                             
 * <br><br>
 * Wrapper for parameter calibration 
 * combining parameters with base volatility bootstrapping
 * (single time-dependent spot volatility).
 * Given benchmarks U_1,...,UN and C1,...,CM
 * (a benchmark being as either a volatility or a correlation point,
 * see <tt>VnfmGenerateBenchmarkVolL</tt>), 
 * performs the optimization of the two factor parameters \\
 * 
 * (beta_i,alpha_i,rho_ij), \\ \\
 *
 * with the time-dependent spot volatility t->sigma(t)
 * calibrated to the base volatility curve.
 * Minimizes the weighted sum of
 * deviation from the market values of the benchmarks Ui
 * <br><br> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
 * min (beta_i,alpha_2,..,alpha_n,rho_ij} Sum{i=1,...N}
 *		w_i * (U_i(market) - U_i(model))^2 <br><br>
 * <br><br>
 * under the constraints that the model values of the benchmarks C_i
 * are within prescribed bounds.
 * <br><br> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
 * for all j,    C_i(min) <= C_i(model) <= C_i(max)
 * <br><br>
 * where w_i,  U_i(market),
 * C_i(min), C_i(max) are specified by the user.<br>
 * Remarks: because of the scaling between the weights alpha and
 * spot volatilities, the parameters alpha1 is not optimized,
 * but instead rescaled in a way consistent with the average
 * volatility of the optimized benchmarks.
 *
 * <br><br>
 * <b>Arguments Details</b><br><br>
 * <blockquote>
 * <li> <b><tt> refDateL </tt></b>  zero curve value date.
 * <li> <b><tt> zcDatesL </tt></b>  zero curve dates.
 * <li> <b><tt> zcRatesL </tt></b>  zero curve rates.
 * <li> <b><tt> datesL </tt></b>  array of dates for calibration timeline.
 * <li> <b><tt> paOptFlagL </tt></b>  array of length N(N+3)/2-1 where N is the number
 * of factors corresponding to the model parameters in the order \\
 * (beta_1,...,beta_n,alpha_1,...,alpha_n,
 *     rho_1,...,rho_1,rho_21,..,rho_n-1,n)\\
 * containing the flag 1 if
 * the corresponding parameter is to be optimized, 0 if not.
 * <li> <b><tt> paMinL </tt></b>  array of lower bounds on the parameters
 *               (same length as <tt> paOptFlagL </tt>).
 * <li> <b><tt> paMaxL </tt></b>  array of upper bounds on the parameters
 *               (same length as <tt> paOptFlagL </tt>).
 * <li> <b><tt> paMidL </tt></b>  array of initial values of parameters
 *              (or value used if the  parameter is not optimized
 *               (same length as <tt> paOptFlagL </tt>).
 * <li> <b><tt> nOptBenchL </tt></b>  number of optimized benchmarks.
 * <li> <b><tt> optBenchL </tt></b>  array of optimized benchmarks definition as a string
 *   See <tt> VnfmGenerateBenchmarkVolL </tt> for a description of the format.
 *   The length of the array can be more than <tt> nOptBenchL </tt>,
 *   but only the first <tt> nOptBenchL </tt> are used.
 * <li> <b><tt> optWeiL </tt></b>  array of optimized benchmarks weights in the optimization.
 *   (same length as <tt> optBenchL </tt>).
 * <li> <b><tt> optMidL </tt></b>  array of optimized benchmarks market values
 *   (same length as <tt> optBenchL </tt>).
 * <li> <b><tt> nConstBenchL </tt></b>  number of constraints benchmarks
 * <li> <b><tt> constBenchL </tt></b>  array of constraint benchmark 
 *   See <tt> VnfmGenerateBenchmarkVolL </tt> for a description of the format.
 *   The length of the array can be more than <tt> nConstBenchL </tt>,
 *   but only the first <tt> nConstBenchL </tt> are used.
 * <li> <b><tt> constMinL </tt></b>  array of constraint min value 
 *   (same length as <tt> constBenchL </tt>).
 * <li> <b><tt> constMaxL </tt></b>  array of constraint max value 
 *   (same length as <tt> constBenchL </tt>).
 * <li> <b><tt> floatScalarsL </tt></b> 	array of length 1.\\
 * Element 1: constraint on minimum spot volatility during
 *   base volatility bootstrapping ( > 0).\\
 * Element 2: back bone q (0 for log-normal, 0.5 for normal).
 * <li> <b><tt> rateMatL </tt></b>  array of volatility rate maturity (in yrs).
 * The element of index i in the array corresponds
 * to an option expiring at date i of array <tt> datesL </tt>.
 * This array must have same length as the array <tt> nfDatesL </tt>.
 * <li> <b><tt> rateFreqL </tt></b>  array of volatility frequency (0,1,2,4,12).
 * {\it The element of index i in the array corresponds
 * to an option expiring at date i of array <tt> datesL </tt>.
 * <li> <b><tt> rateVolL </tt></b>  array of volatilities.
 * {\it The element of index i in the array corresponds
 * to an option expiring at date i of array <tt> datesL </tt>.
 * <li> <b><tt> paOptL </tt></b>  output optimal parameters,
 *   (same length as <tt> paOptFlagL </tt>).
 * <li> <b><tt> optModL </tt></b>  output optimized benchmarks model value,
 *   (same length as <tt> optBenchL </tt>).
 * <li> <b><tt> constModL </tt></b>  output constraint benchmarks model value,
 *   (same length as <tt> constBenchL </tt>).
 * </blockquote>
 */

DLL_EXPORT(int)
VnfmCalibParamSquareL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDatesL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRatesL,	/*  3 'F' (I) zero coupon rates */

	TDateL *datesL,		/*  4 'D' (I) array of dates */

	IntL *paOptFlagL,	/*  5 'L' (I) optimize param flags */
	FloatL *paMinL,		/*  6 'F' (I) optimize param min value */
	FloatL *paMaxL,		/*  7 'F' (I) optimize param max value */
	FloatL *paMidL,		/*  8 'F' (I) optimize param init value */
	
	IntL *nOptBenchL,	/*  9 'L' (I) # of optim. benchmarks */
	CharBlockL *optBenchL,	/* 10 'C' (I) optim. benchmarks */
	FloatL *optWeightL,	/* 11 'F' (I) optim. benchmarks weights */
	FloatL *optMidL,	/* 12 'F' (I) optim. benchmarks market value */

	IntL *nConstBenchL,	/* 13 'L' (I) # of constr. benchmarks */
	CharBlockL *constBenchL,/* 14 'C' (I) constr. benchmarks */
	FloatL *constMinL,	/* 15 'F' (I) constr. benchmarks min value */
	FloatL *constMaxL,	/* 16 'F' (I) constr. benchmarks max value */

	FloatL *floatScalarsL,	/* 17 'F' (I) [1] = min spot vol */
				/*            [2] = dist (0=N,1=LN) */

	FloatL *rateMatL,	/* 18 'F' (I) rate mat [0..nDates-1] */
	IntL *rateFreqL,	/* 19 'L' (I) rate freq [0..nDates-1] */
	FloatL *rateVolL,	/* 20 'F' (I) rate vol [0..nDates-1] */

	FloatL *paOptL,		/* 21 'F' (O) optimal param */
	FloatL *optModL,	/* 22 'F' (O) optim. benchmarks optim value */
	FloatL *constModL)	/* 23 'F' (O) constr. benchmarks optim value */
{
static	char		routine[] = "VnfmCalibParamSquareL";
	int		status = FAILURE;
	VnfmData	*vnfmData = NULL;
	TCurve		*zcCurve = NULL;
	int		nPa,		/* len. of param array */
			nF,		/* # of factors */
			n;		/* dummy */
			/* optimization benchmarks */
	int		nOptBench;
	TVolBenchmark	optBench[VNFM_NMAX_BENCH];
			/* constraint benchmarks */
	int		nConstBench;
	TVolBenchmark	constBench[VNFM_NMAX_BENCH];

	double		*rateMat,
			*rateVol;
	int		*rateFreq = NULL;
	int		nDates;



	/* Addin Loggingf
	 */
	if (GtoLoggingGet() > 0) {
	    DrlLilVectLoggingFile("wrapper.log", "w", "TF_CALPAR_SQ1V",
		DRL_TDATE_L, (void*) refDateL, "ZC_REFDATE",
		DRL_TDATE_L, (void*) zcDatesL, "ZC_DATES",
		DRL_FLOAT_L, (void*) zcRatesL, "ZC_RATES",
		DRL_TDATE_L, (void*) datesL, "IN_DATES",
		DRL_LONG_L,  (void*) paOptFlagL, "OPT_FLAG",
		DRL_FLOAT_L, (void*) paMinL, "OPT_MIN",
		DRL_FLOAT_L, (void*) paMaxL, "OPT_MAX",
		DRL_FLOAT_L, (void*) paMidL, "OPT_MID",
		DRL_LONG_L,  (void*) nOptBenchL, "OBJ_N",
		DRL_CHAR_BLOCK_L,(void*) optBenchL, "OBJ_BENCH",
		DRL_FLOAT_L, (void*) optWeightL, "OBJ_WEI",
		DRL_FLOAT_L, (void*) optMidL,      "OBJ_MID",
		DRL_LONG_L,  (void*) nConstBenchL, "CON_N",
		DRL_CHAR_BLOCK_L,(void*) constBenchL, "CON_BENCH",
		DRL_FLOAT_L, (void*) constMinL,	"CON_MIN",
		DRL_FLOAT_L, (void*) constMaxL,	"CON_MAX",
		DRL_FLOAT_L, (void*) floatScalarsL, "FLOAT_SCALARS",

		DRL_FLOAT_L, (void*) rateMatL,	"RATE_MAT",
		DRL_LONG_L,  (void*) rateFreqL,	"RATE_FREQ",
		DRL_FLOAT_L, (void*) rateVolL,	"RATE_VOL",

		DRL_FLOAT_L, (void*) paOptL,	"OPT_OPT",
		DRL_FLOAT_L, (void*) optModL,	"OBJ_MOD",
		DRL_FLOAT_L, (void*) constModL,	"CON_MOD",
		0);
	}


	VNFM_LOGOPEN


	/*
	 * (1) Get Parameterr
	 */

	/* get zero curve */
	IF_FAILED_DONE( DrlTCurveWrapRead(
		&zcCurve,
		refDateL,
		zcDatesL,
		zcRatesL));
#ifndef	NO_LOGGING
	GTO_IF_LOGGING(DrlTCurveFpWrite(zcCurve, vnfmFpLog,\
		DRL_TCURVE_FMT_PERCENT));
#endif

	/* check array length */
	WRAP_CHECK_VECTOR_LEN(floatScalarsL, 2);

	/* create new data structure */

	/* get number of factors from length of arrays */
	nF = (int) floor(sqrt(2e0*(ARGSIZE(paOptFlagL))+2.25e0) - 1.5 +0.5);
#undef	CHECK_LEN
#define	CHECK_LEN(arg)    {if (ARGSIZE(arg) != nF*(nF+3)/2) {\
		GtoErrMsg("%s: argument %s has length %d incompatible"\
		" with num fact %d.\n", routine,  #arg, (int) ARGSIZE(arg),\
		 nF); goto done;}}

	CHECK_LEN(paOptFlagL);
	CHECK_LEN(paMinL);
	CHECK_LEN(paMaxL);
	CHECK_LEN(paMidL);
	CHECK_LEN(paOptL);

#undef	CHECK_LEN
	nPa = nF*(nF+3)/2;

	/* create new data structure */
	vnfmData = VnfmNewTimeLine(
			nF,
			floatScalarsL[2], /* backboneq */
			zcCurve->fBaseDate,
			ARGSIZE(datesL),
			datesL+1,
			NULL,	/* used for timeline only (can be NULL) */
			zcCurve);
	if (vnfmData == NULL) goto done;



#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog));
#endif

#ifdef	_SKIP
	DrlFPrintStruct(vnfmFpLog, NULL, '\t',
	    DRL_CARRAY_T, "paraopt", nPa, 5,
		DRL_INT_T,    (void*) (paOptFlagL+1),
		DRL_DOUBLE_T, (void*) (paOptL+1),
		DRL_DOUBLE_T, (void*) (paMinL+1),
		DRL_DOUBLE_T, (void*) (paMaxL+1),
		DRL_DOUBLE_T, (void*) (paMidL+1),
	0);
#endif

	/* get benchmarks vols */
	/* input vol benchmarks */
	if (DrlTVolBenchmarkWrapRead(optBench, &nOptBench,
		nOptBenchL, optBenchL) != SUCCESS)
			goto done;

	if (DrlLilStructGet(0,
	    DRL_LILVECTARRAY_L, &n, nOptBench, 128, 3,
		"optWeight", DRL_FLOAT_L, FALSE, (void*) optWeightL, NULL,
		"optMid",    DRL_FLOAT_L, FALSE, (void*) optMidL,    NULL,
		"optMod",    DRL_FLOAT_L, FALSE, (void*) optModL,    NULL,
	    0) != SUCCESS) goto done;


	ASSERT_OR_DONE(ARGSIZE(optWeightL) >= nOptBench);
	ASSERT_OR_DONE(ARGSIZE(optMidL) >= nOptBench);
	ASSERT_OR_DONE(ARGSIZE(optModL) >= nOptBench);

	/* correlation benchmarks */
	if (DrlTVolBenchmarkWrapRead(constBench, &nConstBench,
		nConstBenchL, constBenchL) != SUCCESS)
			goto done;

	if (DrlLilStructGet(0,
	    DRL_LILVECTARRAY_L, &n, nConstBench, 128, 3,
		"constMax", DRL_FLOAT_L, FALSE, (void*) constMaxL, NULL,
		"constMin", DRL_FLOAT_L, FALSE, (void*) constMinL, NULL,
		"constMod", DRL_FLOAT_L, FALSE, (void*) constModL, NULL,
	    0) != SUCCESS) goto done;

	ASSERT_OR_DONE(ARGSIZE(constMaxL) >= nConstBench);
	ASSERT_OR_DONE(ARGSIZE(constMinL) >= nConstBench);
	ASSERT_OR_DONE(ARGSIZE(constModL) >= nConstBench);


	/*
	 * Get bootstrap vol data
	 */

	nDates = ARGSIZE(datesL);
	ASSERT_OR_DONE(ARGSIZE(rateMatL)  == nDates);
	ASSERT_OR_DONE(ARGSIZE(rateFreqL) == nDates);
	ASSERT_OR_DONE(ARGSIZE(rateVolL)  == nDates);
	rateMat = &rateMatL[1];
	rateVol = &rateVolL[1]; 
	rateFreq = NEW_ARRAY(int, nDates); ASSERT_OR_DONE(rateFreq != NULL);
	for (n=0; n<=ARGSIZE(rateFreqL)-1; n++)
		rateFreq[n] = (int) rateFreqL[n+1];


#ifndef	NO_LOGGING


#endif


	/*
	 * (2) Perform the calibration
	 */
	if (VnfmCalibParamSquareVol(
		vnfmData,
		nOptBench,
		optBench,
		optWeightL+1,
		optMidL+1,
		optModL+1,
		nConstBench, constBench, constMinL+1, constMaxL+1, constModL+1,
		floatScalarsL[1],

		rateMat,
		rateFreq,
		rateVol,

		nPa,
		paOptFlagL+1,
		paOptL+1,
		paMinL+1,
		paMaxL+1,
		paMidL+1) != SUCCESS)
			goto done;

#ifndef	NO_LOGGING
	GTO_IF_LOGGING(VnfmFpWrite(vnfmData, vnfmFpLog));
#endif


	/* made it through OK */
	status = SUCCESS;
done:
	GtoFreeTCurve(zcCurve);
	if (rateFreq) FREE(rateFreq);
	VnfmFree(vnfmData);
	if (status != SUCCESS) {
		GtoErrMsg("%s: failed\n", routine);
	}

	VNFM_LOGCLOSE

	return (status);
}




