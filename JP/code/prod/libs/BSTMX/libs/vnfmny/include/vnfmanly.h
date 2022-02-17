/****************************************************************
 * Module:	VNFM
 * Submodule:	ANLY
 * File:	
 * Function:	
 * Author:	Christian Daher & David Fung
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vnfm_H
#define	_vnfm_H
#include "drlstd.h"			/* platform compatibility */
#include <stdio.h>

#include "cerrsup.h"
#include "tcurve.h"	/* term structure library */

/*#define	VNFM_V5X	*/
#define	VNFM_NDIMMAX	32	/* maximum # of dimemsions */

/* NormDist & LognormalDist identifier */
#define VNFM_NORMAL_DIST	(0.5e0)
#define VNFM_LOGNORMAL_DIST	(0e0)

/*t-@CDOC(idxn=VnfmData,catn=structdef)--------------------------
 * VnfmData structure.
 * 
 * <br><br>
 * The structure <i>VnfmData</i> contains the data for
 * an N-factor model with exponential factors.
 *
 * <H3>Conventions for timeline</H3>
 *
 * The time line  is defined by
 * <blockquote>
 * 0=t_0 < t_1 < ... t_{M-1}.
 * </blockquote>
 * If x is any piecewise constant time-dependent quantity,
 * (for example sigma(t), rho(t)),
 * then x_i designates the (constant) value of x(t) of the
 * interval t_i &lt;= t < t_{i+1}.
 * (or t_i &lt; t < +infty if i=N-1).
 * 
 * <h3>Conventions for storage of matrices as vectors</h3>

 * The matrix of correlations (rho_{ij}(t))
 * with 1 &lt;=i &lt; j &lt;= N is
 * stored as a vector <i> rho[k]</i> with the matrix index
 * (i+1,j+1) with 1 &lt;= i &lt; j &lt;= N
 * correponding to the vector index
 * <blockquote>
 * k = {N(N-1)/ 2} - {(N-i)(N-i-1)/ 2} + (j-i-1)
 * </blockquote>
 * with  0&lt;= k &lt;= {N(N-1)/ 2}-1.
 * <br><br>
 * Matrices of covariances J_{ij}(S)
 * for 1&lt;= i&lt;= j &lt;= N
 * are stored as a vectors <i> J[k]</i> with the matrix index
 * (i+1,j+1) with 1&lt;= i&lt;= j &lt;= N
 * correponding to the vector index
 * <blockquote>
 * k = {N(N+1)/ 2} - {(N-i)(N-i+1)/ 2} + (j-i)
 * </blockquote>
 * with  0&lt;= k &lt;= {N(N+1)/ 2}-1.
 * 
 */

typedef	struct {
	int	fNDates;	/* number time steps in timeline */
	double	*fTime,		/* years (ACT/365F) at each node */
		*fDt;		/* step size (dt[i] = t[i+1]-t[i]) */
	TDate	*fDate;		/* dates corresponding to time nodes */
	TCurve	*fZcCurve;	/* zero curve */
#ifndef	VNFM_V5X
	int	fNZDates;	/* number of elements in fZTime */
	double	*fZTime;	/* years (ACT/365F) at each Zero Date in
				 * fZcCurve. Note: 
				 * fZTime[0] = 0.0 (i.e. REFDATE, or fDate[0]),
				 * fZTime[i] = fZcCurve->fArray[i-1].fDate
				 * - REFDATE, for i=1..fZcCurve->fNumItems */
#endif	/*VNFM_V5X*/
				/* Model parameters: */
	int	fNf;		/* number of factors */
	double	*fBeta,		/* array of MR [0..nF-1] */
		*fAlpha,	/* array of factor weight [0..nF-1] */
		**fSigma,	/* sigma[0..nF-1][0..nDates-1] */
		**fRho;		/* correlation[0..nF*(nF-1)/2-1][0..nDates-1] */

	double	*fRate,		/* storage of 1-period cc fwd rate */
		*fZero,		/* storage of zero coupon bonds */
		**fJ;		/* storage of the J integrals */

	double	fBackBoneQ;	/* Backbone coef: 0=normal, 1=lognormal */
	double	fVolNorm;	/* Normal reference volatility */
	double	fVolLogn;	/* Lognormal reference volatility */

	double	fQLeft;		/* Mapping coefficient: left Q */
	double	fQRight;	/* Mapping coefficient: right Q */
	double	fQFShift;	/* Mapping coefficient: fwd shift  */

} VnfmData;

/*e*/


/*
 * Memory Management and Constructors/Destructors
 */

extern	VnfmData*	VnfmNew(int nDates, int nDim, int nZDates);
extern	VnfmData*	VnfmNewTimeLine(
	int nDim,		/* (I) # of factors */
 	double backboneq,	/* (I) backbone */
	TDate refDate,		/* (I) reference date */
	int nDates,		/* (I) # dates (can be -1) */
	TDate *dates,		/* (I) array of dates (can be NULL) */
	TCurve *volCurve,	/* (I) used for timeline only (can be NULL) */
	TCurve *zcCurve);	/* (I) zero curve */

extern	VnfmData*	VnfmNew2SpotVols(
	TDate refDate,		/* (I) vol reference date */
 	double backboneq,	/* (I) backbone */
	double beta1,		/* (I) 2F parameter */
	double beta2,		/* (I) 2F parameter */
	double alpha1,		/* (I) 2F parameter */
	double alpha2,		/* (I) 2F parameter */
	int nDates,		/* (I) # dates */
	TDate *dates,		/* (I) array of dates [0..nDates-1] */
	double *sigma1,		/* (I) 2F parameter [0..nDates-1] */
	double *sigma2,		/* (I) 2F parameter [0..nDates-1] */
	double *rho,		/* (I) 2F parameter [0..nDates-1] */
	TCurve *zcCurve);	/* (I) zero curve */

extern	VnfmData*	VnfmNew2FactSimple(
	TDate refDate,		/* (I) reference date */
	int nDates,		/* (I) # dates (can be -1) */
	TDate *date,		/* (I) array of dates (can be NULL) */
	TCurve *volCurve,	/* (I) used for timeline only (can be NULL) */
	TCurve *zcCurve,	/* (I) zero curve */
 	double backboneq,	/* (I) backbone */
	double beta1,		/* (I) 2F parameter */
	double beta2,		/* (I) 2F parameter */
	double alpha1,		/* (I) 2F parameter */
	double alpha2,		/* (I) 2F parameter */
	double sigma1,		/* (I) 2F parameter */
	double sigma2,		/* (I) 2F parameter */
	double rho);		/* (I) 2F parameter */
extern	int		VnfmFree(VnfmData *);


extern	DLL_EXPORT(int)	VnfmSetSmile(
	VnfmData *that,		/* (B) vnfm data */
	double bbq,		/* (I) backbone q (0=normal, 1=lognormal) */
	double qLeft,		/* (I) mapping coefficient */
	double qRight,		/* (I) mapping coefficient */
	double qFShift);	/* (I) mapping coefficient */




extern	VnfmData*	VnfmCopy(VnfmData *to, VnfmData *from);
extern	DLL_EXPORT(int) VnfmCheckValid(VnfmData *that);

extern	DLL_EXPORT(int)	VnfmCorrMinorDeterm(
	VnfmData *that,		/* (I) model parameters */
	int tpIdx,		/* (I) time point index */
	int order, 		/* (I) order */
	double *determ);	/* (O) determinant */

/*
 * I/O on file and FILE pointer
 */
extern	int	VnfmFileRead(VnfmData **, char *) ;
extern	int	VnfmFileWrite(VnfmData *, char *) ;
extern	int	VnfmFpRead(VnfmData **, FILE *) ;
extern	int	VnfmFpWrite(VnfmData *, FILE *) ;
extern	int	VnfmWrapRead(VnfmData **that, FloatL *backboneqL, FloatL *betaL,
			FloatL *alphaL, TDateL *dateL, FloatL *sigmaL,
			FloatL *rhoL, TCurve *zcCurve);
extern	int	VnfmWrapWrite(VnfmData *that, FloatL *backboneqL, FloatL *betaL,
			FloatL *alphaL, TDateL *dateL, FloatL *sigmaL,
			FloatL *rhoL);
extern	int	VnfmWrapReadSimple(VnfmData **thatp,
			FloatL *backboneqL, FloatL *paramsL,
			TDateL *datesL, TCurve *zcCurve);

extern	DLL_EXPORT(int)	VnfmToGtoParameters(
	VnfmData *that,		/* (I) model parameters */
 	double *backboneq,	/* (I) backbone */
	double *meanRevers,	/* (O) array of mean revers [0..numFact-1] */
	TDate *today,		/* (O) where volatility starts */
	TDate *dates,		/* (O) dates up to which spot vols apply */
	double **spotVols,	/* (O) sp vols[0..numFact-1][0..numDates-1] */
	double **correlations);	/* (O) corr[0..nF*(nF-1)/2-1][0..nDates-1] */


/*
 * Analytical Coefficients Computation
 */

extern	DLL_EXPORT(int)	VnfmComputeCoeff(VnfmData *that);
extern	int	VnfmComputeFwdRates(VnfmData *that);
extern	int	VnfmPrintCoeff(VnfmData *that, FILE *fp);
extern	int	VnfmComputeJ(VnfmData *that);
extern	int	VnfmComputeJPartial(VnfmData *that, int idxStart, int idxEnd);


extern	int	VnfmA(VnfmData *that, double S, double T, double *A);
extern	int	VnfmQ(VnfmData *that, double S, double T, double *y, double *Q);
extern	int	VnfmB(VnfmData *that, double S,
			double tMat, int freq, double *y, double *B);
extern	int	VnfmG(VnfmData *that, double S, double T, double *G);

extern	int	VnfmQBCoeff(VnfmData *that,	
			double S, double rateReset, double rateMat,
			int rateFreq, double *yield, double *QB);


extern	int	VnfmJ(VnfmData *that, double S1, double S2,
			double T, double *J);

extern	int	VnfmI(VnfmData *that, double S1, double S2,
			double T, double *I);

/*
 * Volatilities and Correlations Computation (low-level)
 */

extern	DLL_EXPORT(int)	VnfmAvgQBVol(VnfmData *that, double S1, double S2,
				double tReset, double tMat, int freq,
				KVolType vType,
				double *retVal);
extern	DLL_EXPORT(int)	VnfmAvgQBCorr(VnfmData *that, double S,
				double tReset1, double tMat1, int freq1,
				double tReset2, double tMat2, int freq2,
				double *retVal);
extern	DLL_EXPORT(int)	VnfmAvgQBCorr2(VnfmData *that, double Obs, double S,
				double tReset1, double tMat1, int freq1,
				double tReset2, double tMat2, int freq2,
				double *retVal);
extern	DLL_EXPORT(int)	VnfmAvgQBSprdVol(VnfmData *that, double S,
				double tReset1, double tMat1, int freq1,
				double tReset2, double tMat2, int freq2,
				double *retVal);


extern	int	VnfmVolCurve(VnfmData *that,
				TDateInterval lMat, int freq, int volType,
				TDate bvRefDate, int nBvDates, TDate *bvDates,
				KVolType vType,
				double *bvRates);
extern	int	VnfmSwaptionVolMatrix(VnfmData *that,
				int type, int freq, int adjFlag,
				int nExp, TDateInterval *lExp,
				int nMat, TDateInterval *lMat, 
				KVolType vType,
				double **vol);
extern	int	VnfmAvgCorrMatrix(VnfmData *that,
				TDateInterval lExp, int nMat,
				TDateInterval *lMat, double **corr);

extern int 	VnfmAvgVolONRate(
			VnfmData *that,		/* (I) model parameters */
			TDate   startDate,	/* (I) start date of period */
			TDate   endDate,	/* (I) end date of period */
			TDate   resetDate,	/* (I) reset date of option */
			double *avgVol);	/* (O) vol of avg ON rate */


extern	int	VnfmGenerateVolTCurve(
	VnfmData *that,		/* (I) model parameters */
	TDateInterval lMat,	/* (I) maturity of rate */
	int freq,		/* (I) rate frequency */
	int volType,		/* (I) 0=base, 1=fwd */
	TDate bvRefDate,	/* (I) fwd ref date */
	int nDates,		/* (I) # base vol dates */
	TDate *dates,		/* (I) vol dates (or NULL) */
	TCurve **bvCurve);	/* (O) output base vol curve */

extern	int	VnfmAvgQBOrthFactors(
	VnfmData *that,
	double tStart,		/* (I) start of observation period (yrs) */
	double tEnd,		/* (I) end of observation period (yrs) */
	int numRates,		/* (I) # of maturities  */
	double *rateExpYrs,	/* (I) array of expirations (yrs from tEnd) */
	int *rateFreq,		/* (I) array of rate frequencies */
	double *rateMatYrs,	/* (I) array of maturities (yrs) */
	double **eigVect,	/* (O) array of eigen vectors */
	double *eigVal);	/* (O) array of eigen values */


/*
 * Forward/Futures Adjustment
 */
extern	int	VnfmFutAdjIntegralsCompute(
	VnfmData *that,		/* (I) model parameters */
	double **faInteg);	/* (O) should be allocated */
extern	int	VnfmFutAdjIntegralsInterp(
	VnfmData *that,		/* (I) model parameters */
	double **faInteg,	/* (I) contains integrals to be interp */
	double T,		/* (I) time to perform the interp */
	double *faIntegInterp);	/* (O) should be allocated */
extern	int	VnfmFutAdjRateComputeAdjFactor(
	VnfmData *that,		/* (I) model parameters */
	double **faInteg,	/* (I) contains integrals to be interp */
	double tReset,		/* (I) rate reset */
	double tMat,		/* (I) rate maturity */
	int freq,		/* (I) rate frequency */
	double *retVal);	/* (O) adjustment factor */
extern	int	VnfmFwdFutAdjustment(
	VnfmData *that,		/* (I) model parameters */
	TCurve *zcCurve,	/* (I) zero curve */
	TDateInterval rateMat,	/* (I) rate maturity */
	int rateFreq,		/* (I) rate frequency */
	long rateDayCount,	/* (I) day count convention */
	int nResetDates,	/* (I) # of desired dates */
	TDate *resetDates,	/* (I) array of desired dates */
	double *fwdRate,	/* (O) array of fwd rates (or NULL) */
	double *cvxAdj,		/* (O) array of cvx adj (or NULL) */
	double *mtmAdj);	/* (O) array of mtm adj (or NULL) */


/*
 * Low level spot volatility bootstrapping routines.
 */

extern	DLL_EXPORT(int)	VnfmVolCalib1VArbitrary(
	VnfmData *that,		/* (I/O) model parameters */
	int idxStart,		/* (I) first timeline idx calibrated */
	int idxEnd,		/* (I) last timeline idx calibrated */
	double *tMat1,		/* (I) array of swaption mat [0..idxEnd] */
	int *freq1,		/* (I) array of swaption freq [0..idxEnd] */
	double *vol1,		/* (I) array of swaption vol [0..idxEnd] */
	KVolType vType,		/* (I) LOGVOL, NORMVOL */
	int errMsgFlag,		/* (I) write to errlog iff TRUE */
	double *resValue);	/* (O) <0 if calib OK, >0 otherwise */

extern	DLL_EXPORT(int)	VnfmVolCalib1VArbitraryNew(
	VnfmData *that,		/* (I/O) model parameters */
	int idxStart,		/* (I) first timeline idx calibrated */
	int idxEnd,		/* (I) last timeline idx calibrated */
	double *reset1,         /* (I) array of swaption resets [0..idxEnd] */
	double *tMat1,		/* (I) array of swaption mat [0..idxEnd] */
	int *freq1,		/* (I) array of swaption freq [0..idxEnd] */
	double *vol1,		/* (I) array of swaption vol [0..idxEnd] */
	KVolType vType,		/* (I) LOGVOL, NORMVOL */
	int errMsgFlag,		/* (I) write to errlog iff TRUE */
	double *resValue);	/* (O) <0 if calib OK, >0 otherwise */

extern	DLL_EXPORT(int)	VnfmVolCalib1VFailureId(
	VnfmData *that,		/* (I) Model parameters */
	double *failedExp,	/* (O) Expiration which failed (in years) */
	double *failedMat);	/* (O) Fwd maturity which failed (in years) */

extern	DLL_EXPORT(int)	VnfmVolCalib1VArbitraryMinSpotVol(
	VnfmData *that,		/* (I) model parameters */
	int idxStart,		/* (I) first timeline idx calibrated */
	int idxEnd,		/* (I) last timeline idx calibrated */
	double *tMat1,		/* (I) array of swaption mat [0..idxEnd] */
	int *freq1,		/* (I) array of swaption freq [0..idxEnd] */
	double *vol1,		/* (I) array of swaption vol [0..idxEnd] */
	int errMsgFlag,		/* (I) write to errlog iff TRUE */
	int    spotVolRatioFlag,/* (I) spot vol ratio check flag  */
	double spotVolRatio,    /* (I) spot vol ratio check flag  */
	double minSpotVol);	/* (I) minimum spot volatility */

extern  DLL_EXPORT(int) VnfmSpreadVolCalib1VArbitrary(
	VnfmData *that,		/* (I/O) model parameters */
	int idxStart,		/* (I) first timeline idx calibrated */
	int idxEnd,		/* (I) last timeline idx calibrated */
	double *vol1,		/* (I) array of swaption vol [0..idxEnd] */
	KVolType vType,		/* (I) LOGVOL, NORMVOL */
	int errMsgFlag,		/* (I) write to errlog iff TRUE */
	double *resValue);	/* (O) <0 if calib OK, >0 otherwise */

extern	DLL_EXPORT(int)	VnfmVolCalib2VArbitrary(
	VnfmData *that,		/* (I/O) model parameters */
	int idxStart,		/* (I) first timeline idx calibrated */
	int idxEnd,		/* (I) last timeline idx calibrated */
	TDate *dReset1,		/* (I) 1st array of swapt reset [0..idxEnd] */
	double *tMat1,		/* (I) 1st array of swapt mat [0..idxEnd] */
	int *freq1,		/* (I) 1st array of swapt freq [0..idxEnd] */
	double *vol1,		/* (I) 1st array of swapt vol [0..idxEnd] */
	TDate *dReset2,		/* (I) 2st array of swapt reset [0..idxEnd] */
	double *tMat2,		/* (I) 2nd array of swapt mat [0..idxEnd] */
	int *freq2,		/* (I) 2nd array of swapt freq [0..idxEnd] */
	double *vol2,		/* (I) 2nd array of swaption vol [0..idxEnd] */
	int errMsgFlag,		/* (I) write to errlog iff TRUE */
	double *resValue);	/* (O) <0 if calib OK, >0 otherwise */



/*--------------------------------------------------------------
 * Macros used for some simple integrals:
 * <br>
 *	K(a,b,\alpha) &=& integral_a^b e^{-\alpha*t} dt \\
 *	L(a,\alpha)   &=& integral_0^a e^{-\alpha*t} dt
 * <br>
 */

#ifdef	_vnfm_SOURCE
#undef	K
#define K(a,b,alpha)    (ABS((alpha)*((b)-(a))) < 1e-6 ? \
                        (exp(-(alpha)*(a))*(  \
                         (b) - (a) - 0.5*(alpha)*((b)-(a))*((b)-(a)) \
			+ (alpha)*(alpha)*((b)-(a))*((b)-(a))*((b)-(a))/6) ) : \
                        (exp(-(alpha)*(a)) - exp(-(alpha)*(b))) / (alpha))
#undef	L
#define	L(a,alpha) 	((1 - exp(-(alpha)*(a))) / (alpha))

#undef	NDIM
#define	NDIM		(that->fNf)

#undef	NDATES
#define	NDATES		(that->fNDates)

#undef	NZDATES
#ifndef	VNFM_V5X
#define	NZDATES		(that->fNZDates)
#else
#define	NZDATES		(that->fNDates)
#endif	/*VNFM_V5X*/

#undef	REFDATE
#define	REFDATE		(that->fDate[0])

#undef	VNFMIDX
#define	VNFMIDX(texp, idx)	DrlDoubleArrayFloorIdx(TT, NDATES, (texp)+\
				45.66210046e-6, (idx))

#undef	BETA
#define	BETA	that->fBeta

#undef	ALPHA
#define	ALPHA	that->fAlpha

#undef	SIGMA
#define	SIGMA	that->fSigma

#undef	RHO
#define	RHO		that->fRho

#undef	RATE
#define	RATE	that->fRate

#undef	JJ
#define	JJ		that->fJ

#undef	TT
#define	TT		that->fTime

#undef	ZTT
#ifndef	VNFM_V5X
#define	ZTT		that->fZTime
#else
#define	ZTT		that->fTime
#endif	/*VNFM_V5X*/


/*
 * Function that transforms the rate
 */
#undef	FN
#ifdef	__OLDMAP__
#define	FN(yield, q) (IS_ALMOST_ZERO((q)) ?  (yield) :                     \
                     (IS_ALMOST_ZERO((q)-0.5) ? (1.0) :                    \
                            pow((double)(yield), (double)(1.0-2.0*(q)))))
#else
#define	FN(yield, q) \
		((q) * that->fVolLogn * (yield) + \
			(1e0 - (q)) * that->fVolNorm)
#endif


/*
 * Index mapping
 */

#undef	RHOIDX
/*#define	RHOIDX(i,j)	(nDim*(nDim-1)/2 - (nDim-(i))*(nDim-(i)-1)/2 +\
			 ((j)-(i)-1))*/
#define	RHOIDX(i,j)	((i) <= (j) ? \
		(nDim*(nDim-1)/2 - (nDim-(i))*(nDim-(i)-1)/2 + ((j)-(i)-1)) : \
		(nDim*(nDim-1)/2 - (nDim-(j))*(nDim-(j)-1)/2 + ((i)-(j)-1)))


#undef	JIDX
/*#define	JIDX(i,j)	(nDim*(nDim+1)/2 - (nDim-(i))*(nDim-(i)+1)/2 +\
			 ((j)-(i)))*/
#define	JIDX(i,j)	((i) <= (j) ? \
	(nDim*(nDim+1)/2 - (nDim-(i))*(nDim-(i)+1)/2 + ((j)-(i))) : \
	(nDim*(nDim+1)/2 - (nDim-(j))*(nDim-(j)+1)/2 + ((i)-(j))) )

#undef	HIDX
#define	HIDX(i,j)	(nDim*(i)+(j))

#undef	SQR
#define	SQR(x)	((x)*(x))

#undef	SQRT_SIGN
#define	SQRT_SIGN(x)	((x) < 0e0 ? sqrt(-(x)) : sqrt (x))

/*
 * Macros 
 */
#define	VNFM_PRINT_LEVEL	1	/* default logging level */

#if defined(CDDEV) && (defined(UNIX) || defined(WIN32) || defined(_WIN32))
extern  DLL_EXPORT(long)    DrlGlobFlagGet(int);
#undef	GTO_IF_LOGGING
#define GTO_IF_LOGGING(printStatement)    \
if (DrlGlobFlagGet(18L) >= VNFM_PRINT_LEVEL) \
{                                         \
    printStatement;                       \
}
#else
/* Using */
#undef	GTO_IF_LOGGING
#define GTO_IF_LOGGING(printStatement)    \
if (GtoLoggingGet() >= VNFM_PRINT_LEVEL) \
{                                         \
    printStatement;                       \
}

#endif


#endif	/*_vnfm_SOURCE*/

/*
 * Global File Pointer for Logging
 */
extern	FILE	*vnfmFpLog;


#endif	/* _vnfm_H */

