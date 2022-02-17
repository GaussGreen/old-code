#ifndef ESL_BAS_DOT_H
#define ESL_BAS_DOT_H

/** NOTE: This file should be only included through 'esl_bas.c'
 */


#ifdef  __cplusplus
extern "C" {
#endif

#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"
	
extern const double QCUTOFF;               /* Normal model for |q|<QCUTOFF  */

#define     INVSQRT2PI      0.3989422804   /* 1/sqrt(2*pi)                  */

#ifndef    BS2Q_FAILURE_VALUE
#define    BS2Q_FAILURE_VALUE  -1E10
#endif

/*****  NormalH  ************************************************************/
/**
*      Normal cumulative distribution accrding to J.Hull
*/
double  NormalH (double  x);

/*****  Normal_InvH  ********************************************************/
/**
*      Normal cumulative inverse distribution Risk Magazine
*/
double Normal_InvH (double prob);

/*****  NDensity  ************************************************************/
/**
*       Normal density.
*/
double  NDensity (double x);

/*****  Binary_BS  **********************************************************/
/**
       Binary price using Black&Scholes.
*/
double  Binary_BS ( double  S,    /**< Price of the underlying         */
                    double  K,    /**< Strike                          */
                    double  T,    /**< Binary expiration in years      */
                    double  r,    /**< Risk free rate                  */
                    double  s,    /**< Annualized volatility           */
                    char    CoP); /**< Call or put binary ('C' or 'P') */

/*****  Call_BS  ************************************************************/
/**
       Price of a Call using Black&Scholes.
*/
double  Call_BS (   double S,  /**< Price of the underlying    */
                    double K,  /**< Strike                     */
                    double T,  /**< Option expiration in years */
                    double r,  /**< Risk free rate             */
                    double s); /**< Annualized volatility      */

/*****  Put_BS  *************************************************************/
/**
       Price of a Put using Black&Scholes.
*/
double  Put_BS (double S,  /**< Price of the underlying    */
                double K,  /**< Strike                     */
                double T,  /**< Option expiration in years */
                double r,  /**< Risk free rate             */
                double s); /**< Annualized volatility      */

/*****  D_Call_BS  **********************************************************/
/**
       Delta of a Call using Black&Scholes.
*/
double  D_Call_BS ( double S,   /**< Price of the underlying    */
                    double K,   /**< Strike                     */
                    double T,   /**< Option expiration in years */
                    double r,   /**< Risk free rate             */
                    double s);  /**< Annualized volatility      */

/*****  D_Put_BS  ***********************************************************/
/**
       Delta of a Put using Black&Scholes.
*/
double  D_Put_BS (  double S,  /**< Price of the underlying    */
                    double K,  /**< Strike                     */
                    double T,  /**< Option expiration in years */
                    double r,  /**< Risk free rate             */
                    double s); /**< Annualized volatility      */

/*****  T_Call_BS  **********************************************************/
/**
       Theta of a Call using Black&Scholes.
*/
double  T_Call_BS ( double S, /**< Price of the underlying    */
                    double K, /**< Strike                     */
                    double T, /**< Option expiration in years */
                    double r, /**< Risk free rate             */
                    double s);/**< Annualized volatility      */

/*****  T_Put_BS  ***********************************************************/
/**
       Theta of a Put using Black&Scholes.
*/
double  T_Put_BS (  double S,   /**< Price of the underlying    */
                    double K,   /**< Strike                     */
                    double T,   /**< Option expiration in years */
                    double r,   /**< Risk free rate             */
                    double s);  /**< Annualized volatility      */

/*****  R_Call_BS  **********************************************************/
/**
       Rho of a Call using Black&Scholes.
*/
double  R_Call_BS ( double S, /**< Price of the underlying    */
                    double K, /**< Strike                     */
                    double T, /**< Option expiration in years */
                    double r, /**< Risk free rate             */
                    double s);/**< Annualized volatility      */

/*****  R_Put_BS  ***********************************************************/
/**
       Rho of a Put using Black&Scholes.
*/
double  R_Put_BS (  double  S,  /**< Price of the underlying    */
                    double  K,  /**< Strike                     */
                    double  T,  /**< Option expiration in years */
                    double  r,  /**< Risk free rate             */
                    double  s); /**< Annualized volatility      */

/*****  pdf  ****************************************************************/
/**
*      Density function of the normal distribution.
*/
double  pdf (   double  x);

/*****  BS_Density  **********************************************************/
/**
*      Density function in a Black & Scholes environment.
*/
double  BS_Density (	double S,  /**< (I) Price of the underlying       */
                        double S0, /**< (I) Initial price or expectation  */
                        double T,  /**< (I) Option expiration in years    */
                        double s   /**< (I) Annualized volatility         */
		);

/*****  Call_BSQ  ***********************************************************/
/**
*      Price of a Call using Q version of Black&Scholes.
*/
double  Call_BSQ (double Y,  /**< (I) Fwd yield                     */
                  double K,  /**< (I) Strike                        */
                  double T,  /**< (I) Option expiration in years    */
                  double s,  /**< (I) Annualized volatility         */
                  double Q); /**< (I) Q weight                      */

/*****  Put_BSQ  ************************************************************/
/**
*      Price of a Put using Black&Scholes.
*/
double  Put_BSQ (double Y,   /**< (I) Fwd yield                     */
                 double K,   /**< (I) Strike                        */
                 double T,   /**< (I) Option expiration in years    */
                 double s,   /**< (I) Annualized volatility         */
                 double Q);  /**< (I) Q weight                      */

/*****  Vega_BSQ  ***********************************************************/
/**
*      Vega of put/call Black&Scholes.
*/
double  Vega_BSQ (double Y,  /**< (I) Fwd yield                     */
                  double K,  /**< (I) Strike                        */
                  double T,  /**< (I) Option expiration in years    */
                  double s,  /**< (I) Annualized volatility         */
                  double Q); /**< (I) Q weight                      */
	
/*****  Gamma_BS  **********************************************************/
/**
*      Gamma of an option using Black&Scholes.
*      Remember that the Gamma of a Call is equal to the Gamma of a Put.
*/
double  Gamma_BS (  double        S,   /**< (I) Price of the underlying    */
                    double        K,   /**< (I) Strike                     */
                    double        T,   /**< (I) Option expiration in years */
                    double        r,   /**< (I) Risk free rate             */
                    double        s);  /**< (I) Annualized volatility      */

/*****  Vega_BS  ************************************************************/
/**
*      Vega of an option using Black&Scholes.
*      Remember that the Vega of a Call is equal to the Vega of a Put.
*/
double  Vega_BS (   double         S,   /**< (I) Price of the underlying    */
                    double         K,   /**< (I) Strike                     */
                    double         T,   /**< (I) Option expiration in years */
                    double         r,   /**< (I) Risk free rate             */
                    double         s);  /**< (I) Annualized volatility      */

/*****  CCEquation_BS2Q  ****************************************************/
/*
*       Evalute calibration equation for calibration constant.
*/
double  CCEquation_BS2Q (double    S,      /**< (I) Annualized volatility  */
                         double    QLeft,  /**< (I) Q left                 */
                         double    QRight, /**< (I) Q right                */
                         double    FwdSh,  /**< (I) Fwd shift              */
                         double    CC);    /**< (I) Constant               */

/*****  ConvexityC_BS2Q  ****************************************************/
/**
*      Calibrate constant in 2q formulas. 
*/
int     ConvexityC_BS2Q (double    S,      /**< (I) Annualized volatility  */
                         double    QLeft,  /**< (I) Q left                 */
                         double    QRight, /**< (I) Q right                */
                         double    FwdSh,  /**< (I) Fwd shift              */
                         double    *CC);   /**< (O) Constant               */

/*****  Option_BS2Q  ********************************************************/
/**
*      Price of a call or put using 2Q version of Black&Scholes.
*      Note: vol is defined by Y = Y/(1+fsh) * F((1+fsh)/(1+q*fsh)*vol)
*/
double  Option_BS2Q (double Y,      /**< (I) Fwd yield              */
                     double K,      /**< (I) Strike                 */
                     double T,      /**< (I) Option expiration      */
                     double S,      /**< (I) Annualized volatility  */
                     char   CoP,    /**< (I) Call or put            */
                     double QLeft,  /**< (I) Q left                 */
                     double QRight, /**< (I) Q right                */
                     double FwdSh); /**< (I) Fwd shift              */

/*****  ImpVol_BS2QD  ********************************************************/
/**
        Derive the Implied Vol from a 2Q option price using the dichotomy method.
  		
        Performs the same function as ImpVol_BS2Q, except it is more robust than
        the Newton-Raphson used there.
 */

int     ImpVol_BS2QD (double Y,       /**< (I) Fwd yield          */
                     double K,        /**< (I) Strike             */
                     double T,        /**< (I) Option expiration  */
                     double P,        /**< (I) Price of option    */
                     char   CoP,      /**< (I) Call or put        */
                     double QLeft,    /**< (I) Q left             */
                     double QRight,   /**< (I) Q right            */
                     double FwdSh,    /**< (I) Fwd shift          */
                     double VolGuess, /**< (I) Initial vol guess  */
		     double *ImpVol); /**< (I) Initial vol guess  */

/*****  ImpVol_BS2Q  ********************************************************/
/**
*      Implied vol of a call or put using 2Q version of Black&Scholes.
*/
double  ImpVol_BS2Q (double Y,         /**< (I) Fwd yield              */
                     double K,         /**< (I) Strike                 */
                     double T,         /**< (I) Option expiration      */
                     double P,         /**< (I) Price of option        */
                     char   CoP,       /**< (I) Call or put            */
                     double QLeft,     /**< (I) Q left                 */
                     double QRight,    /**< (I) Q right                */
                     double FwdSh,     /**< (I) Fwd shift              */
                     double VolGuess); /**< (I) Initial vol guess      */

/*****  NormalAS  *************************************************************/
/**
*       Cumulative function of the normal distribution.
*       See Handbook of Mathematical Functions 
*       by Milton Abramowitz and Irene A. Stegun
*       26.2.19 page 932                       
*/

double  NormalAS (double  x);

/*****  KO_ExpAS  **************************************************************/
/**
*	Get the conditional knock-out probability, expectation and variance.
*/
void	KO_ExpAS (
		double 	*P,  /**< Output: probability of not knocking-out   */
                double 	*ES, /**< Output: expected value cond on not K-out  */
                double	*ES2,/**< Output: expected sqr value cond on not KO */
                double 	S,   /**< Underlying expectation                    */
                double 	S0,  /**< Initial price                             */
                double 	H,   /**< Barrier                                   */
                char   	UoD, /**< Up & out or down & out ('U' or 'D')       */
                double 	T,   /**< Option expiration in years                */
                double 	s);  /**< Annualized volatility                     */

/*****  Binary_BS  ***********************************************************/
/**
*      Binary price using Black&Scholes.
*/
double  Binary_BSAS (
		   double       S,   /**< Price of the underlying         */
                   double       K,   /**< Strike                          */
                   double       T,   /**< Binary expiration in years      */
                   double       r,   /**< Risk free rate                  */
                   double       s,   /**< Annualized volatility           */
                   char         CoP);/**< Call or put binary ('C' or 'P') */

/*****  Call_BSAS ************************************************************/
/**
*      Price of a Call using Black&Scholes.
*/
double  Call_BSAS (
		 double S, /**< Price of the underlying    */
                 double K, /**< Strike                     */
                 double T, /**< Option expiration in years */
                 double r, /**< Risk free rate             */
                 double s);/**< Annualized volatility      */

/*****  Put_BSAS  ************************************************************/
/**
*      Price of a Put using Black&Scholes.
*/
double  Put_BSAS (
		double S, /**< Price of the underlying    */
                double K, /**< Strike                     */
                double T, /**< Option expiration in years */
                double r, /**< Risk free rate             */
                double s);/**< Annualized volatility      */

/*****  D_Call_BSAS  *********************************************************/
/**
*      Delta of a Call using Black&Scholes.
*/
double  D_Call_BSAS (
		   double  S, /**< Price of the underlying    */
                   double  K, /**< Strike                     */
                   double  T, /**< Option expiration in years */
                   double  r, /**< Risk free rate             */
                   double  s);/**< Annualized volatility      */

/*****  D_Put_BS  ************************************************************/
/**
*      Delta of a Put using Black&Scholes.
*/
double  D_Put_BSAS (
		  double  S, /**< Price of the underlying    */
                  double  K, /**< Strike                     */
                  double  T, /**< Option expiration in years */
                  double  r, /**< Risk free rate             */
                  double  s);/**< Annualized volatility      */

/*****  T_Call_BSAS  *********************************************************/
/**
*      Theta of a Call using Black&Scholes.
*/
double  T_Call_BSAS (
		   double S, /**< Price of the underlying    */
                   double K, /**< Strike                     */
                   double T, /**< Option expiration in years */
                   double r, /**< Risk free rate             */
                   double s);/**< Annualized volatility      */

/*****  T_Put_BSAS  **********************************************************/
/**
*      Theta of a Put using Black&Scholes.
*/
double  T_Put_BSAS (
		  double S, /**< Price of the underlying    */
                  double K, /**< Strike                     */
                  double T, /**< Option expiration in years */
                  double r, /**< Risk free rate             */
                  double s);/**< Annualized volatility      */

/*****  R_Call_BSAS  *********************************************************/
/**
*      Rho of a Call using Black&Scholes.
*/
double  R_Call_BSAS (
	  	   double S, /**< Price of the underlying    */
                   double K, /**< Strike                     */
                   double T, /**< Option expiration in years */
                   double r, /**< Risk free rate             */
                   double s);/**< Annualized volatility      */	

/*****  R_Put_BSAS  **********************************************************/
/**
*      Rho of a Put using Black&Scholes.
*/
double  R_Put_BSAS (
		  double S,  /**< Price of the underlying    */
                  double K,  /**< Strike                     */
                  double T,  /**< Option expiration in years */
                  double r,  /**< Risk free rate             */
                  double s); /**< Annualized volatility      */


#ifdef  __cplusplus
}
#endif


#endif



