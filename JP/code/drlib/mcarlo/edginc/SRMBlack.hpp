
//----------------------------------------------------------------------------
//
//   Group       : QR cross asset
//
//   Filename    : SRMBlack.hpp
//
//   Description : 2Q Black-Scholes utility functions
//
//   Author      : Henrik Rasmussen
//
//   Date        : 2 May 2006
//
//
//----------------------------------------------------------------------------

#ifndef SRM_BLACK_HPP
#define SRM_BLACK_HPP

#include "edginc/config.hpp"

DRLIB_BEGIN_NAMESPACE

struct SRMBlack // namespace would be better, but beware of VC6
{

/*****  Normal_InvH  ********************************************************/
/*
*      Normal cumulative inverse distribution Risk Magazine
*/
static double Normal_InvH(double prob);

/*****  NDensity  ************************************************************/
/*
*       Normal density. From bas::NDensity
*/
static double NDensity(double x);

/*  Normal cumulative distribution accrding to J.Hull
    This is here to make the numbers match. From bas::NormalH */
static double NormalH(double  x);

/*****  Call_BS  ************************************************************/
/*
*      Price of a Call using Black&Scholes. From bas::Call_BS
*/
static double Call_BS(double S,  /* Price of the underlying    */
                      double K,  /* Strike                     */
			      	  double T,  /* Option expiration in years */
				      double r,  /* Risk free rate             */
				      double s); /* Annualized volatility      */

/*****  Put_BS  *************************************************************/
/*
*      Price of a Put using Black&Scholes. From bas::Put_BS
*/
static double Put_BS(double S,  /* Price of the underlying    */
                     double K,  /* Strike                     */
					 double T,  /* Option expiration in years */
					 double r,  /* Risk free rate             */
					 double s); /* Annualized volatility      */

/*****  CCEquation_BS2Q  ****************************************************
 *       Evalute calibration equation for calibration constant.
 *       From bas::CCEquation_BS2Q
 */
static double CCEquation_BS2Q(
                              double    S,         /* (I) Annualized volatility  */
							  double    QLeft,     /* (I) Q left                 */
							  double    QRight,    /* (I) Q right                */
							  double    FwdSh,     /* (I) Fwd shift              */
							  double    CC);       /* (I) Constant               */

/*****  ConvexityC_BS2Q  ****************************************************/
/*
*      Calibrate constant in 2q formulas. From bas::ConvexityC_BS2Q
*/
static double ConvexityC_BS2Q(
					   double    S,     /* (I) Annualized volatility  */
                       double    QLeft, /* (I) Q left                 */
                       double    QRight,/* (I) Q right                */
                       double    FwdSh);/* (I) Fwd shift              */

/*****  Option_BS2Q  ********************************************************/
/*
*      Price of a call or put using 2Q version of Black&Scholes.
*      Note: vol is defined by Y = Y/(1+fsh) * F((1+fsh)/(1+q*fsh)*vol)
*      From bas::Option_BS2Q
*/
static double Option_BS2Q(
    double Y,                /* (I) Fwd yield              */
    double K,                /* (I) Strike                 */
    double T,                /* (I) Option expiration      */
    double S,                /* (I) Annualized volatility  */
    char   CoP,              /* (I) Call or put            */
    double QLeft,            /* (I) Q left                 */
    double QRight,           /* (I) Q right                */
    double FwdSh);           /* (I) Fwd shift              */

/*****  ImpVol_BS2Q  ********************************************************/
/*
*      Price of a call or put using 2Q version of Black&Scholes.
*      From bas::ImpVol_BS2Q
*/
static double ImpVol_BS2Q(
    double Y,                /* (I) Fwd yield              */
    double K,                /* (I) Strike                 */
    double T,                /* (I) Option expiration      */
    double P,                /* (I) Price of option        */
    char   CoP,              /* (I) Call or put            */
    double QLeft,            /* (I) Q left                 */
    double QRight,           /* (I) Q right                */
    double FwdSh,            /* (I) Fwd shift              */
    double VolGuess);         /* (I) Initial vol guess      */

};

DRLIB_END_NAMESPACE

#endif

