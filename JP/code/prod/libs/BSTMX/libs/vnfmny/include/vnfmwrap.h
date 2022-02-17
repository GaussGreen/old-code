/****************************************************************
 * Module:	VNFM
 * Submodule:	WRAP
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vnfmwrap_H
#define	_vnfmwrap_H
#include "drlstd.h"		/* DLL_EXPORT() */

#ifdef	__DEBUG__
#include "drlio.h"		/* redirect stdout to log */
#endif	/*__DEBUG__*/

#ifdef	_vnfm_SOURCE
#include "drlproc.h"		/* DrlGlobFlag */
#endif	/*_vnfm_SOURCE*/

#include "vnfmanly.h"



/*
 * Routine for general output generation
 * (N-factor, takes spot volatilities as input)
 */

extern	DLL_EXPORT(int)	VnfmGenerateFwdVolL(
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

	FloatL *volsL);		/* 15 'F' (O) zero coupon rates */


extern	DLL_EXPORT(int)	VnfmGenerateVolCurveL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcTDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

 	FloatL *backboneqL,	/*  4 'F' (I) 0=lognormal, 0.5=normal */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) volatility arrays */
	FloatL *rhoL,		/*  9 'F' (I) correlation arrays */

	TDateIntervalL *bvMatL,	/* 10 'F' (I) vol maturity */
	IntL *volFreqL,		/* 11 'L' (I) vol frequency */
	IntL *volTypeL,		/* 12 'L' (I) 0=base, 1=fwd */
	TDateL *bvRefDateL,	/* 13 'D' (I) vol ref date */
	TDateL *bvDatesL,	/* 14 'D' (I) vol dates */
	FloatL *bvRatesL);	/* 15 'F' (O) vol values */

extern	DLL_EXPORT(int)	VnfmGenerateSwaptionMatrixL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

 	FloatL *backboneqL,	/*  4 'F' (I) 0=lognormal, 0.5=normal */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) volatility arrays */
	FloatL *rhoL,		/*  9 'F' (I) correlation arrays */

	TDateL *swRefDateL,	/* 10 'D' (I) base vol ref date */
	IntL *typeL,		/* 11 'L' (I) mattype, freq, adj */
	TDateIntervalL* tMatL,	/* 12 'F' (I) array of mat intervals */
	TDateIntervalL* tExpL,	/* 13 'F' (I) array of exp intervals */
	FloatL *volL);		/* 14 'F' (O) swaption volatilities */

extern	DLL_EXPORT(int)	VnfmGenerateCorrelationMatrixL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

 	FloatL *backboneqL,	/*  4 'F' (I) 0=lognormal, 0.5=normal */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) volatility arrays */
	FloatL *rhoL,		/*  9 'F' (I) correlation arrays */

	TDateIntervalL* tExpL,	/* 10 'F' (I) exp interval */
	TDateIntervalL* tMatL,	/* 11 'F' (I) array of mat intervals */
	FloatL *corrL);		/* 12 'F' (O) correlations */


extern	DLL_EXPORT(int)	VnfmGenerateBenchmarkVolL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcTDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

 	FloatL *backboneqL,	/*  4 'F' (I) 0=lognormal, 0.5=normal */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) volatility arrays */
	FloatL *rhoL,		/*  9 'F' (I) correlation arrays */
	CharBlockL *optBenchL,	/* 10 'C' (I) */
	FloatL *optModL);	/* 11 'F' (O) */

DLL_EXPORT(int)	VnfmGenerateFwdFutAdjustmentL(
	TDateL *refDateL,		/*  1 'D' (I) reference date */
	TDateL *zcDateL,		/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,		/*  3 'F' (I) zero coupon rates */

 	FloatL *backboneqL,	/*  4 'F' (I) 0=lognormal, 0.5=normal */
	FloatL *betaL,			/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,			/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,			/*  7 'D' (I) array of dates */
	FloatL *sigmaL,			/*  8 'F' (I) volatility arrays */
	FloatL *rhoL,			/*  9 'F' (I) correlation arrays */

	TDateIntervalL *rateMatL,	/* 10 'F' (I) rate maturity */
	IntL *rateFreqL,		/* 11 'L' (I) rate frequency */
	CharBlockL *rateDayCountL,		/* 12 'C' (I) rate day count */
	CharBlockL *adjTypeL,		/* 13 'C' (I) adjustment type */
	TDateL *resetDatesL,		/* 14 'D' (I) adj dates */
	FloatL *outputL);		/* 15 'F' (O) adj values */



/*
 * Convenience spot volatility calibration routines (2F form)
 */

extern	DLL_EXPORT(int)	VnfmCalib1V2FGeneralL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */
	
	FloatL *numScalarsL,	/*  4 'F' (I) num scalars:
				 *        [1] distType
				 *        [2] generate base vol 
				 *        [3] generate swaption matrix */
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

	FloatL *spotVolsL,	/* 15 'F' (O) spot vol values */
	FloatL *bvRatesL,	/* 16 'F' (O) base vol values */
	FloatL *swVolL);	/* 17 'F' (O) swaption vol matrix */


/*
 * Generalized from VnfmCalib1V2FGeneralL to allow both lognormal and normal
 * volatility. 
 */

extern	DLL_EXPORT(int)	VnfmCalib1V2FGeneralNewL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */
	
	FloatL *numScalarsL,	/*  4 'F' (I) num scalars:
				 *        [1] distType
				 *        [2] generate base vol 
				 *        [3] generate swaption matrix */
	FloatL *nfParamsL,	/*  5 'F' (I) N-fact params */
	TDateL *nfDatesL,	/*  6 'D' (I) array of dates */

	CharBlockL *inVTypeL,	/*  7 'L' (I) input vol type:LOG, NORM */
	FloatL *rateMatL,	/*  8 'F' (I) rate mat [0..nDates-1] */
	IntL *rateFreqL,	/*  9 'L' (I) rate freq [0..nDates-1] */
	FloatL *rateVolL,	/* 10 'F' (I) rate vol [0..nDates-1] */

	CharBlockL *outVTypeL,	/* 11 'L' (I) output vol type:LOG, NORM */
	TDateIntervalL *bvMatL,	/* 12 'F' (I) vol maturity */
	TDateL *bvDatesL,	/* 13 'D' (I) vol dates */
	IntL *swTypeL,		/* 14 'L' (I) [0]=type, [1]=freq */
	TDateIntervalL* swMatL,	/* 15 'F' (I) array of mat intervals */
	TDateIntervalL* swExpL,	/* 16 'F' (I) array of exp intervals */

	FloatL *spotVolsL,	/* 17 'F' (O) spot vol values */
	FloatL *bvRatesL,	/* 18 'F' (O) base vol values */
	FloatL *swVolL);	/* 19 'F' (O) swaption vol matrix */



/*
 * Calibrate spot volatilities and compute the volatility of simple average 
 * of overnight rates. 
 */
extern	DLL_EXPORT(int)
VnfmCalib1VNFAverageONL(
	TDateL *refDateL,       /*  1 'D' (I) reference date */
	TDateL *zcDateL,        /*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,        /*  3 'F' (I) zero coupon rates */
 
	FloatL *floatScalarsL,  /*  4 'F' (I) num scalars:
                                 *        [1] back bone q
                                 *        [2] minimum volatility rate mat */
	FloatL *nfParamsL,      /*  5 'F' (I) N-fact params */
 
	TDateL *rateDatesL,     /*  6 'D' (I) array of dates */
	FloatL *rateMatL,       /*  7 'F' (I) rate mat [0..nDates-1] */
	IntL   *rateFreqL,      /*  8 'L' (I) rate freq [0..nDates-1] */
	FloatL *rateVolL,       /*  9 'F' (I) rate vol [0..nDates-1] */
 
	TDateL *startDateL,     /* 10 'D' (I) start date of averaging period */
	TDateL *endDateL,       /* 11 'D' (I) end date of averaging period */
	TDateL *resetDateL,     /* 12 'D' (I) reset date of option */
 
	FloatL *avgONVolL);     /* 13 'F' (O) vol of average ON rates */



extern	DLL_EXPORT(int)	VnfmCalib1V2FVolCurveAndGenFutAdjL(
	TDateL *refDateL,		/*  1 'D' (I) reference date */
	TDateL *zcTDateL,		/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,		/*  3 'F' (I) zero coupon rates */
	
	FloatL *nfParamsL,		/*  4 'F' (I) N-fact params */
	TDateL *nfDatesL,		/*  5 'D' (I) array of dates */

	FloatL *volMatL,		/*  6 'F' (I) vol maturity (cms) */
	TDateL *volDatesL,		/*  7 'D' (I) vol dates */
	FloatL *volRatesL,		/*  8 'F' (I) vol values */

	TDateIntervalL *rateMatL,	/*  9 'F' (I) rate maturity */
	IntL *rateFreqL,		/* 10 'L' (I) rate frequency */
	CharBlockL *rateDayCountL,		/* 11 'C' (I) rate day count */
	TDateL *resetDatesL,		/* 12 'D' (I) adj dates */
	FloatL *fwdRateL,		/* 13 'F' (O) forward rate */
	FloatL *cvxAdjL,		/* 14 'F' (O) cvx adjustment */
	FloatL *mtmAdjL);		/* 15 'F' (O) mtm adjustment */



/*
 * Spot volatilities calibration (general form, N-factors)
 */

extern	DLL_EXPORT(int)	VnfmCalib1VVolCurveArbL(
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
	TDateL *rateResetL,     /* 10 'F' (I) array of rate reset[0..nDates-1]*/
	double *rateMatL,	/* 11 'F' (I) array of vol  mat [0..nDates-1] */
	IntL *rateFreqL,	/* 12 'L' (I) array of vol freq [0..nDates-1] */
	FloatL *rateVolL,	/* 13 'F' (I) array of vol [0..nDates-1] */
	FloatL *tStartL,	/* 14 'F' (I) start/end time for calibration */
				/*        Output Calibrated Model: */
	FloatL *oSigmaL);	/* 15 'F' (O) volatility vector (1-D) */


extern	DLL_EXPORT(int)	VnfmCalib1VVolCurveNewL(
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
	FloatL *oSigmaL);	/* 14 'F' (O) volatility vector (1-D) */


extern	DLL_EXPORT(int)	VnfmCalib1VVolCurveOldL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDatesL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRatesL,	/*  3 'F' (I) zero coupon rates */
				/*        Input Model: */
 	FloatL *backboneqL,	/*  4 'F' (I) 0=lognormal, 0.5=normal */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) volatility arrays */
	FloatL *rhoL,		/*  9 'F' (I) correlation arrays */
				/*        Input base vol curve: */
	IntL *calTypeL,		/* 10 'L' (I) calibration type [0] */
	TDateIntervalL *volMatL,/* 11 'F' (I) vol maturity [0] */
	IntL *volFreqL,		/* 12 'L' (I) volatility freqency [0] */
	TDateL *volDatesL,	/* 13 'D' (I) vol dates */
	FloatL *volRates1L,	/* 14 'F' (I) vol values index # 1*/
	FloatL *tStartL,	/* 15 'F' (I) start/end time for calibration */
				/*        Output Calibrated Model: */
	FloatL *oSigmaL);	/* 16 'F' (O) volatility vector (1-D) */

/*
 * 2 spot vol calibration
 */

extern	DLL_EXPORT(int)	VnfmCalib2VVolCurvesL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDatesL,	/*  2 'D' (I) zero coupon dates */
	FloatL *zcRatesL,	/*  3 'F' (I) zero coupon rates */
				/*        Input Model: */
 	FloatL *backboneqL,	/*  4 'F' (I) 0=lognormal, 0.5=normal */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array correlations */
				/*        Input base vol curve: */
	TDate *rateReset1L,	/* 10 'D' (I) rate 1 reset [0..nDates-1] */
	FloatL *rateMat1L,	/* 11 'F' (I) rate 1 mat [0..nDates-1] */
	IntL *rateFreq1L,	/* 12 'L' (I) rate 1 freq [0..nDates-1] */
	FloatL *rateVol1L,	/* 13 'F' (I) rate 1 vol [0..nDates-1] */
	TDate *rateReset2L,	/* 14 'D' (I) rate 2 reset [0..nDates-1] */
	FloatL *rateMat2L,	/* 15 'F' (I) rate 2 mat [0..nDates-1] */
	IntL *rateFreq2L,	/* 16 'L' (I) rate 2 freq [0..nDates-1] */
	FloatL *rateVol2L,	/* 17 'F' (I) rate 2 vol [0..nDates-1] */

	TDateL *calDatesL,	/* 18 'D' (I) array of calib start/end dates*/
				/*        [0] = time start 1-vol calib */
				/*        [1] = time start 2-vol calib */
				/*        [2] = time   end 2-vol calib */
				/*        [3] = time   end 1-vol calib */
				/*        Output Calibrated Model: */

	FloatL *outSigma1L,	/* 19 'F' (O) spot vol 1 array */
	FloatL *outSigma2L);	/* 20 'F' (O) spot vol 2 array */


/*
 * Parameter optimisation
 */

extern	DLL_EXPORT(int)	VnfmCalibParamShortTermL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDatesL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRatesL,	/*  3 'F' (I) zero coupon rates */

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

	FloatL *floatScalarsL,	/* 17 'F' (I) [1] = backbone q */

	FloatL *paOptL,		/* 18 'F' (O) optimal parameters */
	FloatL *optModL,	/* 19 'F' (O) optim. model value */
	FloatL *constModL);	/* 20 'F' (O) constr. model value */

extern	DLL_EXPORT(int)	VnfmCalibParamSquareL(
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
				/*            [2] = backbone q */
	FloatL *rateMatL,       /* 18 'F' (I) rate mat [0..nDates-1] */
	IntL *rateFreqL,        /* 19 'L' (I) rate freq [0..nDates-1] */
	FloatL *rateVolL,       /* 20 'F' (I) rate vol [0..nDates-1] */

	FloatL *paOptL,		/* 21 'F' (O) optimal param */
	FloatL *optModL,	/* 22 'F' (O) optim. benchmarks optim value */
	FloatL *constModL);	/* 23 'F' (O) constr. benchmarks optim value */

/*
 * 
 */

extern	DLL_EXPORT(int)	VnfmCheckSwaptionCalibrationL(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates	*/
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */
	
	FloatL *nfParamsL,	/*  4 'F' (I) N-F params (b1,b2,r,a) */
	FloatL *floatScalarsL,	/*  5 'F' (I) numerical scalars */
				/*	  [1] distribution type (0=N, 1=LN) */
				/*	  [2] min spot vol  */
				/*        [3] min vol calibration maturity  */
	IntL *swTypeL,		/*  6 'L' (I) mattype, freq  */
	TDateIntervalL* swMatL,	/*  7 'F' (I) array of mat intervals */
	TDateIntervalL* swExpL,	/*  8 'F' (I) array of exp intervals */
	FloatL *swVolL,		/*  9 'F' (I) swaption volatilities */

	IntL   *volCheckFlagL,  /* 10 'L' (I) Vol check flags  */
				/*        [1] check spot vol ratio flag */
				/*            1=Yes, 0=No               */
				/*        [2] vol adjutment flag      */
				/*            1=Yes, 0=No        */
				/*        [3] 1 = output adjusted swap vol*/
				/*            0 = output vol correction */
	FloatL *numScalarsL,    /* 11 'L' (I) numerical scalars  */
				/*        [1] max spot vol raio >1.0     */
				/*        [2] min vol adjust amount(optional)*/

	FloatL *tExpToCheckL,	/* 12 'F' (I) array of exp to check */
	IntL   *finalMatL,      /* 13 'L' (I) array of final flags  */
                                /*            TRUE = final; FALSE = CMS  */

	FloatL *tExpFailedL,	/* 14 'F' (O) output swaption vols */
	FloatL *tMatFailedL,	/* 15 'F' (O) output swaption vols  */
	FloatL *outputSWVolL);  /* 16 'F' (O) output modified swaption vols  */


/*
 *
 */
extern	DLL_EXPORT(int)	VnfmSmoothSwaptionMatrixL(
	TDateL *refDateL,	/* 01 'D' (I) reference date */
	TDateL *zcDateL,	/* 02 'D' (I) zero coupondates */
	FloatL *zcRateL,	/* 03 'F' (I) zero coupon rates */

	IntL *swTypeL,		/* 04 'L' (I) array of matrix param [2]: */
				/*        [0] matr type (0=vertical, 1=diag) */
				/*        [1] vol frequency (1,2,4,12) */
	double *swMatL,		/* 05 'F' (I) array of mat intervals */
	double *swExpL,		/* 06 'F' (I) array of exp intervals */
	double *midMktL,	/* 07 'F' (I) mid market matrix */
	double *bidToMidL,	/* 08 'F' (I) bid to mid matrix (>=0) */
	double *matWeigthL,	/* 09 'F' (I) maturity weights */
	double *expWeigthL,	/* 10 'F' (I) expiration weights */

	long *integerScalarsL,	/* 11 'L' (I) numeric scalars */
	double *doubleScalarsL,	/* 12 'F' (I) float scalars */
	double *nfParamsInL,	/* 13 'F' (I) */

	double *nfParamsOutL,	/* 14 'F' (O) output parameters */
	double *outputMktL,	/* 15 'F' (O) output swaption matrix */
	double *spvolMatL);	/* 16 'F' (O) output spot vol matrix */


/*
 *
 */
extern  DLL_EXPORT(int) VnfmDataCreateO(
	TDateL *refDateL,       /* 1 'D' (I) reference date */
	TDateL *zcDatesL,       /* 2 'D' (I) array of zero coupon dates */
	FloatL *zcRatesL,       /* 3 'F' (I) array of zero coupon rates */
 
	FloatL *backboneqL,     /* 4 'F' (I) back bone Q */
	FloatL *betaL,          /* 5 'F' (I) array of mr coeff */
	FloatL *alphaL,         /* 6 'F' (I) array of weight coeff */
	TDateL *dateL,          /* 7 'D' (I) array of dates */
	FloatL *sigmaL,         /* 8 'F' (I) volatility arrays */
	FloatL *rhoL,           /* 9 'F' (I) correlation arrays */
 
	VnfmData **thatO);      /* 10 '' (O) VnfmData structure */



/*
 *
 */
extern	DLL_EXPORT(int)	VnfmAvgQBVolL(
	VnfmData *that,           /* 01 (I) VnfmData structure */
	double *S1L,              /* 02 (I) start */
	double *S2L,              /* 03 (I) expiration */
	double *tResetL,          /* 04 (I) rate reset */
	double *tMatL,            /* 05 (I) rate maturity */
	int    *freqL,            /* 06 (I) rate frequency */
	double *retValL);         /* 07 (O) volatility */


/*
 *
 */
extern	DLL_EXPORT(int)	VnfmAvgQBCorrL(
	VnfmData *that,           /* 01 (I) VnfmData structure */
	double *SL,               /* 02 (I) Expiration */
	double *tReset1L,         /* 03 (I) rate reset 1 */
	double *tMat1L,           /* 04 (I) rate maturity 1 */
	int    *freq1L,           /* 05 (I) rate frequency 1 */
	double *tReset2L,         /* 06 (I) rate reset 2 */
	double *tMat2L,           /* 07 (I) rate maturity 2 */
	int    *freq2L,           /* 08 (I) rate frequency 2 */
	double *retValL);         /* 09 (O) correlation */

/*
 *
 */
extern	DLL_EXPORT(int)	VnfmAvgQBCorr2L(
	VnfmData *that,           /* 01 (I) VnfmData structure */
	double *ObsL,             /* 02 (I) Expiration */
	double *SL,               /* 03 (I) Expiration */
	double *tReset1L,         /* 04 (I) rate reset 1 */
	double *tMat1L,           /* 05 (I) rate maturity 1 */
	int    *freq1L,           /* 06 (I) rate frequency 1 */
	double *tReset2L,         /* 07 (I) rate reset 2 */
	double *tMat2L,           /* 08 (I) rate maturity 2 */
	int    *freq2L,           /* 09 (I) rate frequency 2 */
	double *retValL);         /* 10 (O) correlation */

/*
 *
 */
extern	DLL_EXPORT(int)	VnfmAvgQBSprdVolL(
	VnfmData *that,           /* 01 (I) VnfmData structure */
	double *SL,               /* 02 (I) Expiration */
	double *tReset1L,         /* 03 (I) rate reset 1 */
	double *tMat1L,           /* 04 (I) rate maturity 1 */
	int    *freq1L,           /* 05 (I) rate frequency 1 */
	double *tReset2L,         /* 06 (I) rate reset 2 */
	double *tMat2L,           /* 07 (I) rate maturity 2 */
	int    *freq2L,           /* 08 (I) rate frequency 2 */
	double *retValL);         /* 09 (O) spread vol */

/*
 *
 */
extern	DLL_EXPORT(int) VnfmAvgQBOrthFactorsL(
        VnfmData *that,	  	  /* (I) VnfmData structure */
	double   *tStartL,        /* (I) Start of observation period */
	double   *tEndL,          /* (I) End of observation period */
	double   *rateExpYrsL,    /* (I) Rate expirations in years */
	long     *rateFreqL,      /* (I) Rate frequencies in years */
	double   *rateMatYrsL,    /* (I) Rate maturities in years */
	TMatrix2D **eigVectL,     /* (O) Array of eigenvectors */
	TMatrix2D **eigValL);     /* (O) Array of eigenvalues */


/*
 * Logfile for wrappers
 */
#ifdef	_vnfm_SOURCE

#if defined(CDDEV)

#endif	/*__DEBUG__*/

#if defined(UNIX)
# define VNFM_LOGFNAM	"tfcal5.log"
#elif  defined(__MSC__)
# define VNFM_LOGFNAM	"c:\\tfcal5.log"
#else
# define VNFM_LOGFNAM	"tfcal5.log"
#endif
#define	VNFM_VERSION	"TFCAL5 - "__DATE__ " " __TIME__

#if defined(CDDEV)
#define	VNFM_LOGOPEN	{if (DrlStdoutFileSet(VNFM_LOGFNAM, "w") < 0) \
			return(FAILURE); vnfmFpLog = stdout;\
			fprintf(stdout, "\t%s: %s\n", VNFM_VERSION, routine);}
#define	VNFM_LOGCLOSE	{DrlStdoutClose();}
#else
#define	VNFM_LOGOPEN
#define	VNFM_LOGCLOSE
#endif


#endif	/*_vnfm_SOURCE*/

#endif	/*_vnfmwrap_H*/
