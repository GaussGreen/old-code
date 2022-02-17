/****************************************************************
 * Module:	VNFM
 * Submodule:	WRAP
 * File:	wrpobj.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vnfmwrpobj_H
#define	_vnfmwrpobj_H
#include "drlstd.h"		/* DLL_EXPORT() */
#include "volcurvo.h"
#include "zcurveo.h"

#ifdef	__DEBUG__
#include "drlio.h"		/* redirect stdout to log */
#endif	/*__DEBUG__*/

#ifdef	_vnfm_SOURCE
#include "drlproc.h"		/* DrlGlobFlag */
#endif	/*_vnfm_SOURCE*/

extern  DLL_EXPORT(int) VnfmGenerateFwdVolO(
        TDateL *refDateL,       /*  1 'D' (I) reference date */
        TDateL *zcDateL,        /*  2 'D' (I) zero coupondates */
        FloatL *zcRateL,        /*  3 'F' (I) zero coupon rates */
 
        FloatL *backboneqL,     /*  4 'L' (I) dist type (0=LN, 0.5=N) */
        FloatL *betaL,          /*  5 'F' (I) array of mr coeff */
        FloatL *alphaL,         /*  6 'F' (I) array of weight coeff */
        TDateL *dateL,          /*  7 'D' (I) array of dates */
        FloatL *sigmaL,         /*  8 'F' (I) array of spot volatilities */
        FloatL *rhoL,           /*  9 'F' (I) array of correlation */
 
        TDateL *obsStartDatesL, /* 10 'D' (I) observation start dates */
        TDateL *obsEndDatesL,   /* 11 'D' (I) observation end dates */
        TDateL *resetDatesL,    /* 12 'D' (I) rate reset dates */
        double *maturitiesL,    /* 13 'F' (I) rate maturities */
        long *frequenciesL,     /* 14 'L' (I) rate frequencies */
 
        TVolCurve **fwdVolO);   /* 15 'F' (O) Vols output */


extern	DLL_EXPORT(int)	VnfmGenerateVolCurveO(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcTDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

	FloatL *backboneqL,	/*  4 'L' (I) (0=LN, 0.5=N) */
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
	TVolCurve **bVolO);	/* 15 'F' (O) vol values */

/*
 *
 */
extern	DLL_EXPORT(int)	VnfmSmoothSwaptionMatrixO(
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

	TMatrix2D **nfParamsOutO,/* 14 'F' (O) output parameters */
	TMatrix2D **outputMktO,	/* 15 'F' (O) output swaption matrix */
	TMatrix2D **spvolMatO);	/* 16 'F' (O) output spot vol matrix */



extern	DLL_EXPORT(int)	VnfmCalib1V2FGeneralO(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */
	
	FloatL *backboneqL,	/*  4 'F' (I) num scalars:
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

	TVolCurve **spotVolO,	/* 15 'F' (O) spot vol values */
	TVolCurve **bVolO,	/* 16 'F' (O) base vol values */
	TMatrix2D **swVolO);	/* 17 'F' (O) swaption vol matrix */


extern	DLL_EXPORT(int)	VnfmCalib1V2FGeneralNewO(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */
	
	FloatL *floatScalarsL,  /*  4 'F' (I) num scalars:
				 *        [1] back bone q
				 *        [2] qL (1=normal,0=lognormal)
				 *        [3] qR (1=normal,0=lognormal)
				 *        [4] Fsh
				 *        [5] generate base vol
				 *        [6] generate swaption matrix
				 *        [7] minimum volatility rate mat */

	FloatL *nfParamsL,	/*  5 'F' (I) N-fact params */
	TDateL *nfDatesL,	/*  6 'D' (I) array of dates */

	CharBlockL *inVTypeL,   /*  7 'L' (I) input vol type:LOG, NORM */
	FloatL *rateMatL,	/*  8 'F' (I) rate mat [0..nDates-1] */
	IntL *rateFreqL,	/*  9 'L' (I) rate freq [0..nDates-1] */
	FloatL *rateVolL,	/* 10 'F' (I) rate vol [0..nDates-1] */

	CharBlockL *outVTypeL,  /* 11 'L' (I) output vol type:LOG, NORM */
	TDateIntervalL *bvMatL,	/* 12 'F' (I) vol maturity */
	TDateL *bvDatesL,	/* 13 'D' (I) vol dates */
	IntL *swTypeL,		/* 14 'L' (I) [0]=type, [1]=freq */
	TDateIntervalL* swMatL,	/* 15 'F' (I) array of mat intervals */
	TDateIntervalL* swExpL,	/* 16 'F' (I) array of exp intervals */

	TVolCurve **spotVolO,	/* 17 'F' (O) spot vol values */
	TVolCurve **bVolO,	/* 18 'F' (O) base vol values */
	TMatrix2D **swVolO);	/* 19 'F' (O) swaption vol matrix */


DLL_EXPORT(int)
VnfmGenerateSwaptionMatrixO(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

	FloatL *backboneqL,	/*  4 'L' (I) [1] dist type (0=LN, 0,5=N) */
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
	TDateIntervalL* swMatL,	/* 12 'F' (I) array of mat intervals */
	TDateIntervalL* swExpL,	/* 13 'F' (I) array of exp intervals */
	TMatrix2D  **swVolO);	/* 14 'F' (O) matrix of swaption vols */


extern  DLL_EXPORT(int) VnfmCalib1VVolCurveOldO(
        TDateL *refDateL,       /*  1 'D' (I) reference date */
        TDateL *zcDatesL,       /*  2 'D' (I) zero coupondates */
        FloatL *zcRatesL,       /*  3 'F' (I) zero coupon rates */
                                /*        Input Model: */
        FloatL *backboneqL,     /*  4 'L' (I) (0=LN, 0,5=N) */
        FloatL *betaL,          /*  5 'F' (I) array of mr coeff */
        FloatL *alphaL,         /*  6 'F' (I) array of weight coeff */
        TDateL *dateL,          /*  7 'D' (I) array of dates */
        FloatL *sigmaL,         /*  8 'F' (I) volatility arrays */
        FloatL *rhoL,           /*  9 'F' (I) correlation arrays */
                                /*        Input base vol curve: */
        IntL *calTypeL,         /* 10 'L' (I) calibration type [0] */
        TDateIntervalL *volMatL,/* 11 'F' (I) vol maturity [0] */
        IntL *volFreqL,         /* 12 'L' (I) volatility freqency [0] */
        TDateL *volDatesL,      /* 13 'D' (I) vol dates */
        FloatL *volRates1L,  /* 14 'F' (I) vol values index # 1*/
        FloatL *tStartL,     /* 15 'F' (I) start/end time for calibration */
                                /*        Output Calibrated Model: */
        TVolCurve **oSigmaO);    /* 16 'F' (O) volatility vector (1-D) */
 

extern  DLL_EXPORT(int) VnfmCalib1VVolCurveNewO(
        TDateL *refDateL,       /*  1 'D' (I) reference date */
        TDateL *zcDatesL,       /*  2 'D' (I) zero coupondates */
        FloatL *zcRatesL,       /*  3 'F' (I) zero coupon rates */
                                /*        Input Model: */
        FloatL *backboneqL,     /*  4 'L' (I) (0=LN, 0,5=N) */
        FloatL *betaL,          /*  5 'F' (I) array of mr coeff */
        FloatL *alphaL,         /*  6 'F' (I) array of weight coeff */
        TDateL *dateL,          /*  7 'D' (I) array of dates */
        FloatL *sigmaL,         /*  8 'F' (I) volatility arrays */
        FloatL *rhoL,           /*  9 'F' (I) correlation arrays */
                                /*        Input base vol curve: */
	double *rateMatL,       /* 10,'F' <I> array of vol mat [0..nDates-1] */
        IntL   *rateFreqL,      /* 11,'L' <I> array of vol freq [0..nDates-1] */
        FloatL *rateVolL,       /* 12,'F' <I> array of vol [0..nDates-1] */
        FloatL *tStartL,        /* 13,'F' (I) start/end time for calibration */
                                /*        Output Calibrated Model: */
        TVolCurve **oSigmaL);   /* 14 'F' (O) volatility vector (1-D) */
 

extern  DLL_EXPORT(int) VnfmCalib1VVolCurveArbO(
        TDateL *refDateL,       /*  1 'D' (I) reference date */
        TDateL *zcDatesL,       /*  2 'D' (I) zero coupondates */
        FloatL *zcRatesL,       /*  3 'F' (I) zero coupon rates */
                                /*        Input Model: */
        FloatL *backboneqL,     /*  4 'L' (I) (0=LN, 0,5=N) */
        FloatL *betaL,          /*  5 'F' (I) array of mr coeff */
        FloatL *alphaL,         /*  6 'F' (I) array of weight coeff */
        TDateL *dateL,          /*  7 'D' (I) array of dates */
        FloatL *sigmaL,         /*  8 'F' (I) volatility arrays */
        FloatL *rhoL,           /*  9 'F' (I) correlation arrays */
                                /*        Input benchmark vol curve: */
	TDateL *rateResetL,     /* 10,'F' <I> array of rate reset[0..nDates-1]*/
	double *rateMatL,       /* 11,'F' <I> array of vol mat [0..nDates-1] */
        IntL   *rateFreqL,      /* 12,'L' <I> array of vol freq [0..nDates-1] */
        FloatL *rateVolL,       /* 13,'F' <I> array of vol [0..nDates-1] */
        FloatL *tStartL,        /* 14,'F' (I) start/end time for calibration */
                                /*        Output Calibrated Model: */
        TVolCurve **oSigmaL);   /* 15 'F' (O) volatility vector (1-D) */


extern  DLL_EXPORT(int) VnfmCalibParamShortTermO(
        TDateL *refDateL,       /*  1 'D' (I) reference date */
        TDateL *zcDatesL,       /*  2 'D' (I) zero coupondates */
        FloatL *zcRatesL,       /*  3 'F' (I) zero coupon rates */
 
        TDateL *datesL,         /*  4 'D' (I) array of dates */
                                /*        Input Model: */
        IntL *paOptFlagL,       /*  5 'L' (I) optimize flag */
        FloatL *paMinL,         /*  6 'F' (I) min value */
        FloatL *paMaxL,         /*  7 'F' (I) max value */
        FloatL *paMidL,         /*  8 'F' (I) initial value */
 
        IntL *nOptBenchL,       /*  9 'L' (I) # of optimized benchmarks */
        char *optBenchL,     /* 10 'C' (I) optim. benchmarks */
        FloatL *optWeiL,        /* 11 'F' (I) optim. market weight */
        FloatL *optMidL,        /* 12 'F' (I) optim. market value */
 
        IntL *nConstBenchL,     /* 13 'L' (I) # of constraints benchmarks */
        char *constBenchL,   /* 14 'C' (I) constr. benchmark */
        FloatL *constMinL,      /* 15 'F' (I) constr. min value */
        FloatL *constMaxL,      /* 16 'F' (I) constr. max value */
 
        FloatL *floatScalarsL,  /* 17 'F' (I) [1] = dist (0=LN, 0.5=N) */
 
        TMatrix2D **paOptL,         /* 18 'F' (O) optimal parameters */
        TMatrix2D **optModL,        /* 19 'F' (O) optim. model value */
        TMatrix2D **constModL);     /* 20 'F' (O) constr. model value */

extern  DLL_EXPORT(int) VnfmCalibParamSquareO(
        TDateL *refDateL,       /*  1 'D' (I) reference date */
        TDateL *zcDatesL,       /*  2 'D' (I) zero coupondates */
        FloatL *zcRatesL,       /*  3 'F' (I) zero coupon rates */
 
        TDateL *datesL,         /*  4 'D' (I) array of dates */
 
        IntL *paOptFlagL,       /*  5 'L' (I) optimize param flags */
        FloatL *paMinL,         /*  6 'F' (I) optimize param min value */
        FloatL *paMaxL,         /*  7 'F' (I) optimize param max value */
        FloatL *paMidL,         /*  8 'F' (I) optimize param init value */
 
        IntL *nOptBenchL,       /*  9 'L' (I) # of optim. benchmarks */
        char *optBenchL,     /* 10 'C' (I) optim. benchmarks */
        FloatL *optWeightL,     /* 11 'F' (I) optim. benchmarks weights */
        FloatL *optMidL,        /* 12 'F' (I) optim. benchmarks market value */
 
        IntL *nConstBenchL,     /* 13 'L' (I) # of constr. benchmarks */
        char *constBenchL,   /* 14 'C' (I) constr. benchmarks */
        FloatL *constMinL,      /* 15 'F' (I) constr. benchmarks min value */
        FloatL *constMaxL,      /* 16 'F' (I) constr. benchmarks max value */
 
        FloatL *floatScalarsL,  /* 17 'F' (I) [1] = min spot vol */
                                /*            [2] = dist (0=LN, 0.5=N) */
 
	FloatL *rateMatL,       /* 18 'F' (I) rate mat [0..nDates-1] */
	IntL *rateFreqL,        /* 19 'L' (I) rate freq [0..nDates-1] */
	FloatL *rateVolL,       /* 20 'F' (I) rate vol [0..nDates-1] */
 
        TMatrix2D **paOptL,      /* 20 'F' (O) optimal param */
        TMatrix2D **optModL,     /* 21 'F' (O) optim. benchmarks optim value */
        TMatrix2D **constModL);  /* 22 'F' (O) constr. benchmarks optim value */


extern  DLL_EXPORT(int) VnfmCalib2VVolCurvesO(
        TDateL *refDateL,       /*  1 'D' (I) reference date */
        TDateL *zcDatesL,       /*  2 'D' (I) zero coupon dates */
        FloatL *zcRatesL,       /*  3 'F' (I) zero coupon rates */
                                /*        Input Model: */
        FloatL *backboneqL,            /*  4 'L' (I) number of dimensions */
        FloatL *betaL,          /*  5 'F' (I) array of mr coeff */
        FloatL *alphaL,         /*  6 'F' (I) array of weight coeff */
        TDateL *dateL,          /*  7 'D' (I) array of dates */
        FloatL *sigmaL,         /*  8 'F' (I) array of volatilities */
        FloatL *rhoL,           /*  9 'F' (I) array correlations */
                                /*        Input base vol curve: */
        TDate *rateReset1L,     /* 10 'D' (I) rate 1 reset [0..nDates-1] */
        FloatL *rateMat1L,      /* 11 'F' (I) rate 1 mat [0..nDates-1] */
        IntL *rateFreq1L,       /* 12 'L' (I) rate 1 freq [0..nDates-1] */
        FloatL *rateVol1L,      /* 13 'F' (I) rate 1 vol [0..nDates-1] */
        TDate *rateReset2L,     /* 14 'D' (I) rate 2 reset [0..nDates-1] */
        FloatL *rateMat2L,      /* 15 'F' (I) rate 2 mat [0..nDates-1] */
        IntL *rateFreq2L,       /* 16 'L' (I) rate 2 freq [0..nDates-1] */
        FloatL *rateVol2L,      /* 17 'F' (I) rate 2 vol [0..nDates-1] */
 
        TDateL *calDatesL,      /* 18 'D' (I) array of calib start/end dates*/
                                /*        [0] = time start 1-vol calib */
                                /*        [1] = time start 2-vol calib */
                                /*        [2] = time   end 2-vol calib */
                                /*        [3] = time   end 1-vol calib */
                                /*        Output Calibrated Model: */
 
        TMatrix2D **outSigma1L,     /* 19 'F' (O) spot vol 1 array */
        TMatrix2D **outSigma2L);    /* 20 'F' (O) spot vol 2 array */
 
extern  DLL_EXPORT(int) VnfmGenerateCorrelationMatrixO(
        TDateL *refDateL,       /*  1 'D' (I) reference date */
        TDateL *zcDateL,        /*  2 'D' (I) zero coupondates */
        FloatL *zcRateL,        /*  3 'F' (I) zero coupon rates */

        FloatL *backboneqL,     /*  4 'F' (I) [1] backbone q (0=LN, 0.5=N) */
        FloatL *betaL,          /*  5 'F' (I) array of mr coeff */
        FloatL *alphaL,         /*  6 'F' (I) array of weight coeff */
        TDateL *dateL,          /*  7 'D' (I) array of dates */
        FloatL *sigmaL,         /*  8 'F' (I) array of spot volatilities */
        FloatL *rhoL,           /*  9 'F' (I) array of correlation */

        TDateIntervalL* tExpL,  /* 10 'F' (I) expiration time */
        TDateIntervalL* tMatL,  /* 11 'F' (I) array of mat times */
        TMatrix2D  **corrO);    /* 12 'F' (O) matrix of correlations */

extern  DLL_EXPORT(int) VnfmCheckSwaptionCalibrationO(
        TDateL *refDateL,       /*  1 'D' (I) reference date */
        TDateL *zcDateL,        /*  2 'D' (I) zero coupondates  */
        FloatL *zcRateL,        /*  3 'F' (I) zero coupon rates */

	FloatL *nfParamsL,      /*  4 'F' (I) N-F params (b1,b2,r,a) */
        FloatL *floatScalarsL,  /*  5 'F' (I) numerical scalars */
                                /*        [1] vol twk size */
                                /*        [2] distribution type (0=LN, 0.5=N) */
                                /*        [3] min spot vol  */
        IntL *swTypeL,          /*  6 'L' (I) mattype, freq  */
        TDateIntervalL *swMatL, /*  7 'F' (I) array of mat intervals */
        TDateIntervalL *swExpL, /*  8 'F' (I) array of exp intervals */
        FloatL *swVolL,         /*  9 'F' (I) swaption volatilities */

	IntL   *volCheckFlagL,  /* 10 'L' (I) Vol check flags  */
				/*        [1] check spot vol ratio flag */
				/*            1=Yes, 0=No               */
				/*        [2] vol adjutment flag      */
				/*            1=Yes, 0=No        */
				/*        [3] 1 = output adjusted swap vol*/
				/*            0 = output vol correction */
	FloatL *numScalarsL,    /* 11 'L' (I) numerical scalars  */
				/*        [1] min vol adjust amount  */
				/*        [2] max spot vol raio >1.0     */

        FloatL *tExpToCheckL,   /* 12 'F' (I) array of exp to check */
	IntL   *finalFagL,      /* 13 'L' (I) array of final flags  */
                                /*            TRUE = final; FALSE = CMS  */

        TMatrix2D **tExpFailedO, /* 14 'F' (O) output expiration list */
        TMatrix2D **tMatFailedO, /* 15 'F' (O) output maturiry list  */
	TMatrix2D **outputSWVolO); /* 16 'F' (O) output mod swaption vols  */

#endif	/*_vnfmwrpobj_H*/
