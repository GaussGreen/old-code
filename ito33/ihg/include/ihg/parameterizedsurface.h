/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parameterizedsurface.h
// Purpose:     Paramterized surface class
// Author:      David
// Created:     04.01.26
// RCS-ID:      $Id: parameterizedsurface.h,v 1.9 2004/10/04 18:04:04 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/parameterizedsurface.h
 */

#ifndef _IHG_PARAMETERIZEDSURFACE_H_
#define _IHG_PARAMETERIZEDSURFACE_H_

#include "ito33/beforestd.h"
#include <list>
#include "ito33/afterstd.h"

#include "ito33/array.h"

namespace ito33
{

namespace ihg
{

/** 
  Paramaterized surface class

  Based on the ideas in Tavella, Klopfer, Implying Local Volatility,
  Wilmott Magazine, Aug. 2001.  This is to be used for any local
  surface (volatility, hazard rate) in the ihg model.  A local
  surface is assumed to be parameterized.  Specifically, a 
  parameterized function is used in the spave direction, such
  as f(S) = c_1 + c_2 * S + c_3 * S^2. The c_i parameters can change at 
  specified times, so should be written as c_i(t).  Between the discrete
  times, linear interpolation is used to construct the parameters.

  Derived classes will specify the actual equations.

*/
class ParameterizedSurface
{
public:

  /** 
    Constructor

    Parameters must be set with the SetParams function. The number of params
    is set here.
  */
  ParameterizedSurface(size_t nNbParams) : m_nNbParams(nNbParams) { }

  /**
    Destructor

    Nothing to do
  */
  virtual ~ParameterizedSurface() {}
  
  /** 
    Write the surface to the specified data stream.

    Added as a helper function for the write() functions in
    the exisiting volatility and hazard rate classes.

    @param out the stream to which the vol surface is written
  */
  void Write(std::ostream &out) const;

  /**
    Return the number of underlying paramaters

    Assumed to be constant at all time points

    @return the number of underlying values
  */
  size_t GetNbParams() const
  {
    return m_nNbParams;
  }

  /**
    Set the parameter values at a specific time

    Used for calibration

    @param dTime the time at which these parameters apply
    @pdParamValues the actual paramater values
  */
  void SetParams(double dTime, double* pdParamValues);

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
    @param dTime is the point at which to compute the smoothness
    @return the smoothness measure (should be >= 0)
  */
  virtual double GetSmoothnessMeasure(double dTime, double dS) = 0;


protected:

  /// Find the parameter data interval
  void FindInterval(double dTime) const;

  /// List of discrete times at which the params are defined
  mutable std::list<double> m_listOfTimes;

  /// Parameter data, corresponding to the above times
  mutable std::list< Array<double> > m_listOfParams;

  /// Number of params
  size_t m_nNbParams;

  /// Pointers (iterators) into the lists for finding current time interval
  mutable std::list< Array<double> >::iterator m_iterParamsLeft;
  mutable std::list< Array<double> >::iterator m_iterParamsRight;
  mutable std::list< double >::iterator m_iterTimesLeft;
  mutable std::list< double >::iterator m_iterTimesRight;

};

} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_PARAMETERIZEDSURFACE_H_

