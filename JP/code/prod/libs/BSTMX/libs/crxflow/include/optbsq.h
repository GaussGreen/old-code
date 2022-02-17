/*********************************************************************************
 * OPTBSQ.H 
 * read IR and credit info from files
 *
 ********************************************************************************/

#ifndef __CRX_OPTBSQ_H__
#define __CRX_OPTBSQ_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <common/include/drmacros.h>
#include <common/include/drtypes.h>
#include <common/include/drutils.h>

#include "crxutil.h"    
    
#define  OPT_MIN_FWD             1E-6
#define  OPT_MIN_VOL             1E-8
#define  OPT_MAX_VOL             2.50
#define  OPT_Q_SHIFT             1E-6
#define  OPT_MQ_RESN             1E-14      /* min df/dq resolution in NR  */
    
/*f----------------------------------------------------------------------------
 *      BS price formula
 *      no discounting!
 */
int BSPricer(
    double     *result,                   /* (O) Price, Vega, or Delta          */
    KOptType   optType,                   /* (I) Option type                    */
	double     Y,                         /* (I) Fwd yield                      */
	double     K,                         /* (I) Strike                         */
	double     T,                         /* (I) Option expiration in years     */
	double     s,                         /* (I) Annualized volatility          */
    KOptResult optResult);                /* (I) Return result type             */

/*f----------------------------------------------------------------------------
 *      Implied vol of a call or put using Black-Scholes.
 */
int     BSImpVol   (
    double     *impVol,                   /* (O) Implied BS vol                 */
    double     yield,                     /* (I) Fwd yield                      */
    double     strike,                    /* (I) Strike                         */
    double     expiry,                    /* (I) Option expiration              */
    double     price,                     /* (I) Price of option                */
    KOptType   optType,                   /* (I) Option type                    */
    double     volGuess);                 /* (I) Initial vol guess              */

    
/*f----------------------------------------------------------------------------
 * Single Q Black Pricer
 *
 * Price of a Call/Put with Vega using Q version of Black&Scholes.
 */
int BSQPricer(
    double     *result,                   /* (O) Price & Vega                   */
    KOptType   optType,                   /* (I) Option type                    */
    double     Y,                         /* (I) Fwd yield                      */
    double     K,                         /* (I) Strike                         */
    double     T,                         /* (I) Option expiration in years     */
    double     s,                         /* (I) Annualized volatility          */
    double     Q,                         /* (I) Q weight                       */
    KOptResult optResult);                /* Return result type                 */

/*f----------------------------------------------------------------------------
 *      Implied vol of a call or put using single Q Black-Scholes.
 */
int     BSQImpVol   (
    double     *impVol,                   /* (O) Implied BS vol                 */
    double     yield,                     /* (I) Fwd yield                      */
    double     strike,                    /* (I) Strike                         */
    double     expiry,                    /* (I) Option expiration              */
    double     Q,                         /* (I) Q parameter                    */
    double     price,                     /* (I) Price of option                */
    KOptType   optType,                   /* (I)Option type                     */
    double     volGuess);                 /* (I) Initial vol guess              */


/**
 * Normal density function
 */
double  NormDens (double x); 

/*f----------------------------------------------------------------------------
 * Double precision cumulative normal function
 */
double NormCum (double x);

/*f----------------------------------------------------------------------------
 * Inverse cumulative normal distribution
 *
 * Based on Risk Magazine
 */
double NormCumInv (double prob);


/*f----------------------------------------------------------------------------
 */
double  BSQIntCR (double x1, double x2,double a,double b,double q,double x0,double C, double sig);


/*f----------------------------------------------------------------------------
 * Exp(a^2) times difference of cumulative normal distribution
 * Computes exp(b+a^2/2)[N(x2-a)-N(x1-a)] (via erfc)
 */
double  BSQInt (double  a, double b, double x1, double x2);
    
/*f----------------------------------------------------------------------------
 * Cumulative error function weighted by exponential factor. Computes
 * exp(b+a^2/2)erfc(x-a)
 *
 * Alib comment:  The routine has a relative accuracy no worse than 1.0E-14, 
 * where relative accuracy is defined as (computed - truth)/truth, and truth
 * comes from a continued fraction calculation.  This is essentially 
 * accurate to the next to last decimal digit of machine accuracy on the Sun.
 */
double ExpCErrFcn (double a, double b, double x);

/******************************************************************************
 * Crx2QCalibrate
 * Given the Q-Smile details, calculates the 2-Q 'A' proportionality constant
 * and the 2Q total volatility sigmaQ
 *****************************************************************************/
int Crx2QCalibrate(
    double  forward,        /**<The forward that you want to match           */
    double  atmOptionPrice, /**<The ATM option price you want to match       */
    double  qL,             /**<Left q parameter                             */
    double  qR,             /**<Right q parameter                            */
    double  muQ,            /**<Q forward shift parameter                    */
    double* A,              /**<(O) Calibrated Q A proprtionality parameter  */
    double* sigmaQ          /**<(O) Calibrated Q total volatility parameter  */
    );

/******************************************************************************
 * Crx2QForward
 * Calculates the forward price under a 2-Q model
 * The underlying price, P, is assumed to be a function of a normal R Var, X:
 * P = A(1 + (exp(X qL)-1))/qL) if X<=0
 *   = A(1 + (exp(X qL)-1))/qR) if X>=0
 * where X ~ Normal(muQ, sigmaQ^2)
 * Returns FAILURE if anything goes wrong, otherwise SUCCESS
 *****************************************************************************/
int Crx2QForward(
    double A,         /**<(I) Proportionality constant                       */
    double muQ,       /**<(I) Expectation of driver variable, X.             */
    double sigmaQ,    /**<(I) Total vol (sigma sqrt(T)) of driver X.         */
    double qL,        /**<(I) Left side Q                                    */
    double qR,        /**<(I) Right side Q                                   */
    double* forward   /**<(O) Calculated forward                             */
    );

/******************************************************************************
 * Crx2QOptionPrice
 * Calculates the price of call and put options under a 2-Q model
 * The underlying price, P, is assumed to be a function of a normal RV, X:
 * P = A(1 + (exp(X qL)-1))/qL) if X<=0
 *   = A(1 + (exp(X qL)-1))/qR) if X>=0
 * where X ~ Normal(muQ, sigmaQ^2)
 * Returns FAILURE if anything goes wrong, otherwise SUCCESS
 *****************************************************************************/
int Crx2QOptionPrice(
    int optType,     /**<GtoOPTION_CALL or GtoOPTION_PUT                     */
    double   K,      /**<(I) Strike                                          */
    double   A,      /**<(I) Forward proportionality constant                */
    double   muQ,    /**<(I) Forward shift                                   */
    double   sigmaQ, /**<(I) Total volatility of driver, X                   */
    double   qL,     /**<(I) Left side Q                                     */
    double   qR,     /**<(I) Right side Q                                    */
    double*  price   /**<(O) Option price                                    */
    );

/******************************************************************************
 * Crx2QImpliedVol
 * Calculates the implied total vol (sigma sqrt(T)!) of call and put options 
 * under a 2-Q model
 * The underlying price, P, is assumed to be a function of a normal RV, X:
 * P = A(1 + (exp(X qL)-1))/qL) if X<=0
 *   = A(1 + (exp(X qL)-1))/qR) if X>=0
 * where X ~ Normal(muQ, sigmaQ^2)
 * This routine calculates sigmaQ so that the price matches that given
 * Returns FAILURE if anything goes wrong, otherwise SUCCESS
 *****************************************************************************/
int Crx2QImpliedVol(
    int optType,    /**<GtoOPTION_CALL or GtoOPTION_PUT                      */
    double   K,     /**<(I) Strike                                           */
    double   A,     /**<(I) Proprtional calibration parameter                */
    double   muQ,   /**<(I) Forward shift = expectation of Gaussian driver   */
    double   price, /**<(I) Option price to match                            */
    double   qL,    /**<(I) Left side Q (0 is normal)                        */
    double   qR,    /**<(I) Right side Q (0 is normal)                       */
    double*  vol    /**<(O) Resulting 2-Q volatility                         */
    );


/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif
