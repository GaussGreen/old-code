/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parameterizedhr.h
// Purpose:     Paramaterized hazard rate surface
// Author:      David
// Created:     04/01/26
// RCS-ID:      $Id: parameterizedhr.h,v 1.12 2004/10/04 18:04:04 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/parameterizedhr.h
   @brief Paramaterized hazard rate surface
 */


#ifndef _IHG_PARAMETERIZEDHR_H_
#define _IHG_PARAMETERIZEDHR_H_

#include "ihg/parameterizedsurface.h"
#include "ito33/ihg/hazardrate.h"


namespace ito33
{

namespace ihg
{

/** 
  Paramaterized hazard rate surface

  Based on the ideas in Tavella, Klopfer, Implying Local Volatility,
  Wilmott Magazine, Aug. 2001. It is expected that this class will be used
  for calibration.

  This is NOT a general hazard rate surface class.  It implements a hazard
  rate surface that is parameterized in space at discrete times. In other 
  words, a parameterized hr curve is specified in the space direction, such
  as hr = c_1 + c_2 * S + c_3 * S^2.  The c_i parameters can change at 
  specified times, so should be written as c_i(t).  Between the discrete
  times, linear interpolation is used to construct the parameters.

  The derived classes must implement a specific formula

  If a different formula is desired, a new class needs to be made.  A general
  base hazard rate surface class could be defined, but seems like overkill 

  See also the paramterized volatility class(es)
*/
  class ParameterizedHR : public HazardRate, public ParameterizedSurface
{
public:

  /** 
    Constructor

    Initialize the surface class
  */
  ParameterizedHR(size_t nNbParams) 
    : HazardRate(), ParameterizedSurface(nNbParams) 
  { 
  }

  /**
    Destructor

    Nothing to do
  */
  virtual ~ParameterizedHR() { }


  /** 
    Get the hazard rates at the specified points and time

    Pure virtual in base hazard rate class.

    @param dTime the time to get the hazard rate values at
    @param pdSpots the asset values at which to get vol values
    @param pdValues is overwritten with the hazard rate values
    @param nNumber is the number of asset values
  */
  virtual void GetHazardRates(double dTime,
                              const double *pdSpots, 
                              double *pdValues, 
                              size_t nNumber) const = 0;

  /**
    Get the parameter values at a specific time

    Used for calibration.  Cannot simply average (interpolate) the values
    at adjacent times, since the type of parameterization is 
    determined by the derived class

    @param dTime the time at which the parameters are requested
    @param pdParamValues is filled in with the actual paramater values
    @param dLeftBound is a left bound for which the returned params should be valid
    @param dRightBound is a right bound for which the returned params should be valid
  */
  virtual void GetParams(double dTime, double* pdParamValues,
    double dLeftBound, double dRightBound) const = 0;

  /** 
    Get default values for the parameters

    These default values can be used to start calibration, if a better guess
    is not known.  Space is allocated by the caller.

    @param pdDefaultParams is filled in with the default parameter values
  */
  virtual void GetParamDefaults(double* pdDefaultParams) const = 0;

  /**
    Get the lower and upper bounds to the paramaters

    Used for calibration

    @param pdLowerBounds is filled in with the lower param bounds
    @param pdUpperBounds is filled in with the upper param bounds
  */
  virtual void GetParamBounds(double* pdLowerBounds, double* pdUpperBounds) const = 0;

  /**
    Get an indication of the surface smoothness

    Used for calibration as a regularization term. Part of the objective
    function is to maximize the smoothness (minimize the derivative) of 
    the surface.  A smoothness of 0 means flat. Typically, this
    measure will be related to partial derivatives of the surface,
    especially the time derivative.

    @param dTime is the time at which to compute the smoothness
    @param dS is the point at which to compute the smoothness
    @return the smoothness measure (should be >= 0)
  */
  virtual double GetSmoothnessMeasure(double dTime, double dS) = 0;

  /**
    Output in gnuplot data file format

    Output in a format to be used by gnuplot 'splot' command.
    gnuplot> set view 50,50,1,1
    gnuplot> splot "surface.dat" with lines

    @param out the stream to which the vol surface is written
  */
  void WriteGnuPlot(std::ostream &out);

};

} // namespace ihg

} // namespace ito33


#endif // #ifndef _IHG_PARAMETERIZEDHR_H_

