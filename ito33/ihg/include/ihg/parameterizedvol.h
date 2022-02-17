/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parameterizedvol.h
// Purpose:     Paramterized volatility surface class
// Author:      David
// Created:     04.01.24
// RCS-ID:      $Id: parameterizedvol.h,v 1.14 2004/10/04 18:04:04 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/parameterizedvol.h
 */

#ifndef _IHG_PARAMETERIZEDVOL_H_
#define _IHG_PARAMETERIZEDVOL_H_

#include "ito33/beforestd.h"
#include <list>
#include "ito33/afterstd.h"

#include "ito33/array.h"

#include "ito33/ihg/volatility.h"
#include "ihg/parameterizedsurface.h"

namespace ito33
{

namespace ihg
{

  class VolatilityVisitor;

/** 
  Paramaterized volatility surface class

  Based on the equations given in Tavella, Klopfer, Implying Local Volatility,
  Wilmott Magazine, Aug. 2001. It is expected that this class will be used
  for calibration.

  This is NOT a general vol surface class.  It implements a vol surface
  that is parameterized in space at discrete times. In other words,
  a parameterized vol curve is specified in the space direction, such
  as vol = c_1 + c_2 * S + c_3 * S^2.  The c_i parameters can change at 
  specified times, so should be written as c_i(t).  Between the discrete
  times, linear interpolation is used to construct the parameters.

  The derived class must specify the actual equations used.
*/
  class ParameterizedVol : public Volatility, public ParameterizedSurface
{
public:

  /** 
    Constructor

    Initialize the surface class
  */
  ParameterizedVol(size_t nNbParams) 
    : Volatility(), ParameterizedSurface(nNbParams) { }

  /**
    Destructor

    Nothing to do
  */
  virtual ~ParameterizedVol() {}
  
  /** 
    Average volatility value

    Used to determine grid sizes.  Pure virtual in base class

    @param dSpotInitial is the place at which to compute an average value
  */
  virtual double AverageValue(double dSpotInitial) = 0;


  /** 
    Get the volatilities at the specified points and time

    Pure virtual in base class.

    @param dTime the time to get the vol values at
    @param Spots the asset values at which to get vol values
    @param Vols is overwritten with the vol values
    @param nNbSpots is the number of asset values
  */
  virtual void GetVols(double dTime, const double *Spots,
                       double *Vols, size_t nNbSpots) const = 0;


  /** 
    Get the volatility squared at the specified points and time

    Pure virtual in base class.

    @param dTime the time to get the vol values at
    @param Spots the asset values at which to get vol values
    @param Vols is overwritten with the vol values
    @param nNbSpots is the number of asset values
  */  
  virtual void GetVolsSquared(double dTime, const double *Spots,
                              double *Vols, size_t nNbSpots) const;

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

  /**
      Dump all data of this instrument in XML format.

      @param tagParent the parent tag under which our tag(s) should be created
      @noexport
   */
  virtual void Dump(XML::Tag& tagParent) const = 0;

    /**
      This method is part of the implementation of the visitor pattern.

      Visitor pattern allows to do different things for different kinds of
      volatility without adding virtual functions to do all of them in the
      base volatility class itself nor doing the dreaded switchs on typeid of
      the concrete type. Instead, define a new class deriving from
      VolatilityVisitor and do whatever is required in its methods.

      @noexport
   */
  virtual void Visit(VolatilityVisitor& visitor) const = 0;
};

} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_PARAMETERIZEDVOL_H_

