/****************************************************************
 * Module:	VNFM
 * Submodule:	CALI
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vnfmcali_H
#define	_vnfmcali_H
#include <stdio.h>
#include <math.h>

#include "vnfmanly.h"



/*
 * Single spot volatility calibration routines.
 */

extern	DLL_EXPORT(int)	VnfmCalib1VVolCurve(
	VnfmData *that,		/* (I/O) model parameters */
	int calType1,		/* (I) 0=cms, 1=final */
	double volMat1,		/* (I) swaption maturity */
	int volFreq1,		/* (I) 0=simple, 1, 2, 4, 12 */
	TCurve *volCurve1,	/* (I) input base volatility curve */
	double tStart,		/* (I) time to start vol calib */
	double tEnd,		/* (I) time to end vol calib */
	int *bootstrapError,	/* (O) reports bootstrap failure iff NULL */
	double *resValue);	/* (O) <0 if calib OK, >0 otherwise */

extern	DLL_EXPORT(int)	VnfmCalib1VCMSFinalVolTCurve(
	VnfmData *that,		/* (I/O) Model parameters */
	int idxStart,		/* (I) TL idx of 1st expiry date calibrated */
	int idxEnd,		/* (I) TL idx of last expiry date calibrated */
	TBoolean finalMaturity,	/* (I) F=cms(use intval), T=final(use date) */
	TDate finalMatDate,	/* (I) Final mat date (used if finalMat=T)*/
	TDateInterval volMat1,	/* (I) Mat interval (used if finalMat=F)*/
	TCurve *volCurve1);	/* (I) Input base volatility curve */


/*
 * Two spot volatility calibration routines.
 */
extern	DLL_EXPORT(int)	VnfmCalib2VVolCurves(
	VnfmData *that,		/* (I/O) model parameters */
	int calType1,		/* (I) 0=sw cms, 1=sw final */
	double volMat1,		/* (I) swaption maturity */
	int volFreq1,		/* (I) 0=simple, 1, 2, 4, 12 */
	TCurve *volCurve1,	/* (I) input base volatility curve */
	int calType2,		/* (I) 0=sw cms, 1=sw final */
	double volMat2,		/* (I) swaption maturity */
	int volFreq2,		/* (I) 0=bvol, 1=swaption */
	TCurve *volCurve2,	/* (I) input base volatility curve */
	double tStart1,		/* (I) time to start 1-vol calib */
	double tStart2,		/* (I) time to start 2-vol calib */
	double tEnd2,		/* (I) time to   end 2-vol calib */
	double tEnd1,		/* (I) time to   end 1-vol calib */
	int errMsgFlag,		/* (I) reports bootstrap failure iff TRUE */
	double *resValue);	/* (O) <0 if calib OK, >0 otherwise */


/*
 * Convenience routines for 2F volatility bootstrapping calibration.
 */

extern	DLL_EXPORT(int)	VnfmCalib1V2FGeneral(
	TCurve *zcCurve,	/* (I) zero coupon term structure */
 	double backboneq,	/* (I) backbone */
	double beta1,		/* (I) 2F parameters */
	double beta2,		/* (I) idem  */
	double alpha,		/* (I) idem  */
	double rho,		/* (I) idem  */
	TDate refDate,		/* (I) today's date */
	int nDates,		/* (I) # dates */
	TDate *dates,		/* (I) timeline dates [0..nDates-1] */
	double *rateMat,	/* (I) arrays of rate mat [0..nDates-1] */
	int *rateFreq,		/* (I) arrays of rate freq [0..nDates-1] */
	double *rateVol,	/* (I) arrays of rate vol [0..nDates-1] */

	int nSpotVols,		/* (I) # of spot vols (can be 0) */
	double *spotVols,	/* (O) array of spot vols (can be NULL) */

	int nBvDates,		/* (I) # of base vol dates */
	TDate *bvDates,		/* (I) array of base vol dates */
	TDateInterval bvMat,	/* (I) base vol maturity */
	TCurve **bvCurve,	/* (O) model base volatility curve (or NULL) */
	TSwaptionMatrix2D *swMat);/* (O) model swaption matrix (or NULL) */

/*
 * Calibration test
 */

extern	DLL_EXPORT(int)	VnfmCheckSwaptionCalibration(
	TCurve *zcCurve,		/*  (I) zero curve */
	TSwaptionMatrix2D *swoptMat,	/*  (I) input matrix */

	double *nfParamsL,		/*  (I) mr, alphas, etc. */
	double backboneQ,		/*  (I) 0.5=N, 0=LN, etc. */
	double volTwkSize,		/* (I) vol tweak size */
	double spotVolMin,		/*  (I) minimum spot volatility */
	double tMatMin,         	/*  (I) Minimum calibration maturity */

	int     volBoundFlag,           /*  (I) vol adj type: 1, -1, 0  */
	int     spotVolRatioFlag,       /*  (I) spot vol ratio check flag  */
	double  spotVolRatio,           /*  (I) Max spot vol ratio   */

	int 	nMat,			/*  (I) # of maturities  */
	double *tMat,			/*  (I) array of maturities */
	long   *finalFlag,		/*  (I) TRUE=final, FALSE=fwd mat */

	double *tExpFailed,		/*  (O) failed exp times (or NULL) */
	double *tMatFailed,		/*  (O) failed mat times (or NULL) */
	double *volAdj,                 /*  (O) vol adjustment needed  */
	int    *nFailed);		/*  (O) # of failures (or NULL) */



extern	DLL_EXPORT(int)	VnfmCheckAdjustSwaptionCalibration(
	TCurve *zcCurve,		/*  (I) zero curve */
	TSwaptionMatrix2D *swoptMat,	/*  (B) input matrix */

	double *nfParamsL,		/*  (I) mr, alphas, etc. */
	double backboneQ,		/*  (I) 0.5=N, 0=LN, etc. */
	double volTwkSize,		/* (I) vol tweak size */
	double spotVolMin,		/*  (I) minimum spot volatility */
	double tMatMin,         	/*  (I) Minimum calibration maturity */

	int     volAdjFlag,             /*  (I) vol adjustment flag   */
        double  volAdjAmount,           /*  (I) vol adjustment amount */

	int     volBoundFlag,           /*  (I) vol adj type: 1, -1, 0  */
	int     spotVolRatioFlag,       /*  (I) spot vol ratio check flag  */
	double  spotVolRatio,           /*  (I) Max spot vol ratio  */

	int 	nMat,			/*  (I) # of maturities  */
	double *tMat,			/*  (I) array of maturities */
	long   *finalFlag,		/*  (I) TRUE=final, FALSE=fwd mat */

	double **tExpFailedTotal,	/*  (O) failed exp times (or NULL) */
	double **tMatFailedTotal,	/*  (O) failed mat times (or NULL) */
	int    *nFailedTotal);		/*  (O) Total # of failures  */

#endif	/* _vnfmcali_H */

