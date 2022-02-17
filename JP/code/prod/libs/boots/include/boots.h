/****************************************************************
 * Module:	Boots
 * File:	
 * Function:	
 * Author:	CA
 * Revision:
 *****************************************************************/
#ifndef	_Boots_H
#define	_Boots_H
#include "drlstd.h"         /* platform compatibility */
#include "drlts.h"			/* DCurve */
#include "drlsmat.h"		/* DSwoptMat */

#include <stdio.h>

/*t--------------------------------------------------------------
 * BootsData structure.
 *
 * See Notes for notation.
 *
 *
 *
 *
 */
typedef	struct {
    /*Actual data*/
	DCurve	*fZcCurve;    /* Zero curve*/
    double *newTimeLine,  /* Storage of new time line*/
           *newZeroBond;  /* and zeroBond line*/
                          
	int freq;             /* Swaption frequency 
                             annual/semi-annual*/
	int newNbDates;       /* Maximum size for the timeline*/
	int nbMatMax;		  /* Maximum number of maturities*/
	int nbExpMax;		  /* Maximum number of expiries*/
	int fCorr;            /* if 1: Forward correls as targets (optimizer)*/
						  /* if 0: Swaption correls as targets (optimizer)*/

	double dt;            /* Time Step*/
	int expIndex;		  /* Keep track of the bootstrapping index
                           * (corresponding expiry index)*/

    double **lambdaMat;   /* Matrix copy of the Lambda Tensor
                             for a special expiry*/
    double ***lambda;	  /* Tensor : link coeffs volFRA/volSWAP
                           * at each expiry date*/

    /*Bootstrapping results*/
    double **kovMat;	  /* Covariance Matrix at a given expiry*/
    double ***kovTens;	  /* Covariance tensor, stores the values
                           * of the successive K matrices*/
    
    /*Temporay arrays*/ 
	double *alpha;	       /* Swaption vols minus covariance 
                            * already bootstrapped
                              (the diagonal of the corresponding matrix)*/
	double **invLambdaMat; /* Storage of LambdaMat inverse
                            * for efficiency*/
	double **invLambdaMatT;/* Storage of LambdaMat
                            * inverse transposed for efficiency*/
	double **alphaMat;	   /* Swaptions vols minus covariance
                            * already bootstrapped
                            * at each time-step of the boostrapping*/
	double **corrMat;	   /* Correlation matrix as implied by the Low-Specs.*/
	double **productMat;   /* Intermediary container for matrices products*/
    double **productMat2;  /* Intermediary container for matrices products*/
    double **productMat3;  /* Intermediary container for matrices products*/
	int    nbSwpMat;       /* Number of maturities to fit at a given expiry*/
	double *swpMat;        /* Points on a vector of swaption maturities*/

} BootsData;

/*t--------------------------------------------------------------
 * Target Correlations structure.
 * Wraps every target correlation for every expiry. 
 *
 */
typedef	struct {
    double *val;           /* Contains the value of a correlation
                            * or of a volatility. [volatilitity if x == y] */ 
	double *x;			   /* Contains the index of the forward1 */
	double *y;			   /* Contains the index of the forward2 */
	int nb;				   /* Contains the number of target correlations*/
	int method;			   /* Method type, see written doc for full description*/
} Target;

/*t--------------------------------------------------------------
 * Additional Input Params structure.
 * Anything that is not Swaption Matrix or Target Correlation. 
 * Input by the user.
 */
typedef	struct {
	int viewIndex;			 
	double *smoothingCoeffs; /* Smoothing coefficient for lambdaInverse*/
	double *initialThetas;   /* Initial values for the Thetas vector*/
} AddInput;

/*t--------------------------------------------------------------
 * SwMat structure.
 * Wraps market data related to swaption matrix. 
 *
 */
typedef	struct {
	double **swaptionMatrix; /* Initial Swaption Matrix*/
	double *expiry;          /* Expiry vector*/
	double *maturity;		 /* Maturity vector*/
	int nbExp;				 /* Total number of expiries*/
	int nbMat;				 /* Total number of maturities*/
	int freq;				 /* Frequency (1=A, 2=S, etc.) */
} SwMat;

/*
 * Memory Management and Constructors/Destructors 
 * BootsData.c
 */
extern  BootsData*  
BootsNew(
    int nbMaxMat,     /*(I)*/
    int nbMaxExp,     /*(I)*/
    int freq);        /*(I)*/

extern BootsData*
BootsNewCopy(
    BootsData* that); /*(I/O) Model parameters*/

extern  int 
BootsFree(
    BootsData *);     /*(I)*/
     
extern int
BootsPrintCoeff(
    BootsData *that, /*(I) model parameters*/
    double expiry,   /*(I) expiry in years*/
    FILE *fp);       /*(O) ouptut file*/

extern	int
BootsWrapRead(
    BootsData **that, /*(O)*/
    DCurve *zcCurve,  /*(I)*/
    int nbMaxMat,     /*(I)*/
    int nbMaxExp,     /*(I)*/
    int freq);        /*(I)*/

extern int
BootsComputeLambda(
    BootsData *that); /* (I/O) model parameters*/

extern int
BootsProcessSwaptionMatrix(
    SwMat    *swaptionMatrix,    /*(I/O) Swaption matrix 
                                  * before and after interpolation*/
    int     nbExpiryToMatch,     /*(I) Number of expiries to match*/
    double  *expiryToMatch,      /*(I) Values of the points
                                  * to match (in expiry)*/
    int     nbMaturityToMatch,   /*(I) Number of maturities to match*/
    double     *maturityToMatch, /*(I) Values of the of the points
                                  * to match (in maturity)*/
    double  terminalVol,         /*(I) Total decrease in vol*/        
    double  terminalSlope,       /*(I) Interpolation type*/
    int     freq);               /*(I) Frequency*/

/*
 * Methods, updating routines.
 * bootsmeth.c
 */
extern int 
BootsGetNewCorrMatrix(
    BootsData* that,        /*(I/O) model parameters */
    Target* targetCorr,     /*(I) Target correlations input*/
    double* Thetas);        /*(I) Thetas parameters*/

extern int 
BootsGetNewKovMatrix(
    BootsData* that,        /*(I/O) Model parameters */
    int numSwap);           /*(I) Number of maturities to be matched
                             *    in the swaption line*/
extern int 
BootsPrecomputeNewKov(
    BootsData* that,        /*(I/O) Model parameters */
    SwMat* sMat,            /*(I) Input SwaptionMatrix*/
    AddInput* addI,         /*(I) Additional Input Parameters*/
    int expiryIndex,        /*(I) ExpiryIndex for update*/
    int nbSwp,              /*(I) Number of points to use from SwaptionMatrix*/
    double* swpMat);        /*(I) Swaption maturities to be used 
                             *    (maturity indexes starting at 0)*/


extern	int 
BootsCopyResults(
    BootsData* that,   /*(I/O) Model parameters*/
    int expiryIndex);  /*(I) Expiry Index for storage in tensor*/


extern int 
BootsOptimizeCorrelations(
    BootsData* that,        /*(I) Model parameters */
    Target* target,         /*(I) Target correlations*/
    double* initialThetas); /*(I) Initial thetas*/

/*
* Main.
* bootsmain.c
*/
extern int 
BootsMain(
    BootsData* that,  /*(I/O) Model parameters*/
    SwMat* sMat,      /*(I) Initial SwaptionMatrix*/
    Target* target,   /*(I) Target correlations*/
    AddInput* addI,   /*(I) Additional Input Parameters*/
	int nbExpTM,      /*(I) Number of expiries to match*/
	double* exp,      /*(I) Vector of expiries*/
	int nbMatTM,      /*(I) Number of maturities to match*/
	double* mat);     /*(I) Vector of maturities*/

extern int 
BootsCalibNoOpt(
    BootsData* that,   /*(I/O) Model parameters */
    SwMat* sMat,       /*(I) Initial SwaptionMatrix*/
	double** corr,     /*(I) The correlation matrix (numMat*numMat)*/
    AddInput* addI,    /*(I) Additional Input Parameters*/
	int numExp,        /*(I) Number of expiries to match*/
	double* swExp,     /*(I) Where "BootsPrecomputeNewKov" will be used
	                    *    First expiry should be first time step*/
	int numMat,        /*(I) Number of tenors to match*/
	double* swMat);    /*(I) The values of the tenors*/

/*
 * Output routines.
 * bootsout.c
 */
extern int 
BootsRecomputeSwaptionVol(
    BootsData* that,       /*(I/O) Model parameters */
    SwMat* sMat,           /*(I) Initial SwaptionMatrix*/
    double** outSwpMat,    /*(O) Swaption Matrix from model*/
    double** errSwpMat);   /*(O) Errors on output SwaptionMatrix*/

extern int 
BootsDataRateVol(
    BootsData* that,        /*(I)*/ 
    double    rateObs,      /*(I)*/
    double    rateReset,    /*(I)*/
    double    rateStart,    /*(I)*/
    double    rateEnd,      /*(I)*/
    int    rateFreq,        /*(I), specific frequency of the rate, 
                             * different from timeline freq*/  
    int    outputType,      /*(I), 0 = %vol, 1 = bp vol*/
    double *outputVol);     /*(O)*/
      
extern int 
BootsDataRateCorrel(
    BootsData* that,        /*(I)*/ 
    double    rateObs,      /*(I)*/
    double    rateReset,    /*(I)*/
    double    rateStart1,   /*(I)*/
    double    rateEnd1,     /*(I)*/
    double    rateStart2,   /*(I)*/
    double    rateEnd2,     /*(I)*/
    int    rateFreq,        /*(I), specific frequency of the rate, 
                             * different from timeline freq*/  
    int    outputType,      /*(I), 0 = %vol, 1 = bp vol*/
    double *outputCorrels); /*(O)*/

extern int
BootsPrintMat(
    double*** matrix,   /*(I) matrix to be printed */
    int nl,             /*(I) number of cols */
    int nc,             /*(I) number of cols */
    FILE *fp);          /*(O) ouptut file*/

extern int 
BootsImpForCorr(
    BootsData* that,
    int index, 
    double** outMat);

extern int 
BootsGenSwapCorrs(
    BootsData* that,
    double** inputForCorr,
	double expiry,
	int nbFor,
    double*** outputForCorr
);



/*
 * Calibration routines.
 * calibtool.c
 */

extern int
frprmn1(BootsData *that,    /*(I) BootsData parameters*/
        Target* target,     /*(I) Target correlations input*/
        double* p,          /*(I) Starting point*/
        int n,              /*(I) Dimension*/
        double  ftol,       /*(I) Tolerance on convergence*/
        int *iter,          /*(I) Maximum number of iterations*/
        double  *fret);     /*(O) Minimum value of the function*/

/*
 * Objective functions definition.
 * objfunc.c
 */

extern int
BootsObjFunc(
    BootsData* that,    /*(I) Model parameters */
    Target* target,     /*(I) Target correlations */
    double* Thetas,     /*(I) Values of Input parameters*/
    double* res);       /*(O) Output value of the function*/

extern int
BootsdObjFunc(
    BootsData* that,    /*(I) Model parameters */
    Target* target,     /*(I) Target correlations */
    double* thetas,     /*(I) Values of Input parameters*/
    double* derivative);/*(O) Value of the derivative*/

/***************************************************************
 *
 *	Macros used in the source code
 *
 ***************************************************************/
#ifdef	_Boots_SOURCE



/*gives index in a ZeroTimeline */

#undef	LAMBDA
#define	LAMBDA	that->lambda

#undef	LAMBDAM
#define	LAMBDAM	that->lambdaMat

#undef	KOVM
#define	KOVM that->kovMat

#undef	KOVT
#define	KOVT that->kovTens


#endif


#endif	/* _Boots_H */

