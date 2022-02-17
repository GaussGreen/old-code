///////////////////////////////////////////////////////////////////////////////
// File:             ito33/numeric/gammadistribution.h
// Purpose:          useful functions related to gamma distributions
// Author:           Ito33 Canada
// Created:          April 24, 2006
// RCS-ID:           $Id: gammadistribution.h,v 1.1 2006/04/28 16:19:04 yann Exp $
// Copyright:        (c) 2006 - Trilemma LLP all rights reserved
///////////////////////////////////////////////////////////////////////////////

/**

  @file ito33/numeric/gammadistribution.h

  @brief Defines useful functions related to gamma distributions. See
  numerical recipes. See Chapter 6.

 */

#ifndef _ITO33_NUMERIC_GAMMADISTRIBUTION_H_
#define _ITO33_NUMERIC_GAMMADISTRIBUTION_H_

namespace ito33
{

namespace numeric
{


/**
    Computes the natural log of the gamma function.  
    Use the numerical recipes handbook.

    @param dX strictly positive value

    @return value at dX
 */
double GammaLn(double dX);


/**
    Computes the chi square cumulative distribution based 
    on the GammaP distribution. See numerical recipes
    handbook Chapter 6.
    
    @param dA degree of freedom
    @param dX value at which we want to compute the distribution
    
    @return chi square distribution with degree A at x
 */
double Chi2Cdf(double dA, double dX);

/**
    Computes the inverse chi square cumulative distribution using
    newton bisection method.
       
    @param dA degree of freedom.
    @param dX value at which we want to compute the distribution.

    @return inverse chi square distribution with degree A at x
 */
double InvChi2Cdf(double dA, double dX);

/**
    Computes the cumulative gamma distribution using a 
    serie decomposition. See numerical recipes Chapter 6
    for details.
   
    @param dA degree of freedom
    @param dX value at which we want to compute the distribution

    @return value for degree dA at dX
 */
double GammaSeries(double dA, double dX);


/**
    Computes the incomplete Gamma P function. See numerical
    recipes Chapter 6 for details.

    @param dA degree of freedom
    @param dX value at which we want to compute the distribution
    
    @return incomplete gamma P function for degree dA at dX
 */
double IncompleteGammaP(double dA, double dX);

/**
    Computes the incomplete Gamma Q function. See numerical
    recipes Chapter 6 for details.

    @param dA degree of freedom
    @param dX value at which we want to compute the distribution

    @return incomplete gamma Q function for degree dA at dX
    
 */
double IncompleteGammaQ(double dA, double dX);

/**
    Computes the cumulative gamma distribution using a 
    continuous fraction. See numerical recipes Chapter 6
    for details.
   
    @param dA degree of freedom
    @param dX value at which we want to compute the distribution

    @return value for degree dA at dX
 */
double GammaContinuousFraction(double dA, double dX);

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_GAMMADISTRIBUTION_H_
