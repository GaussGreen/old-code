/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/polynomialhr.h
// Purpose:     Paramaterized hazard rate surface using polynomials
// Author:      David
// Created:     04/01/26
// RCS-ID:      $Id: polynomialhr.h,v 1.12 2004/10/04 18:04:04 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ihg/polynomialhr.h
  @brief Paramaterized hazard rate surface using polynomials
  */


#ifndef _IHG_POLYNOMIALHR_H_
#define _IHG_POLYNOMIALHR_H_

#include "ihg/parameterizedhr.h"

namespace ito33
{

namespace ihg
{

/** 
  Paramaterized hazard rate surface

  Based on the ideas in Tavella, Klopfer, Implying Local Volatility,
  Wilmott Magazine, Aug. 2001. It is expected that this class will be used
  for calibration.

  The actual equations being used are:
  hazardrate(S,t) = c1 + to be filled in ...

  See also the paramterized volatility class(es)
*/
class PolynomialHR : public ParameterizedHR
{
public:

  /** 
    Constructor

    Initialize the surface class
  */
  PolynomialHR() : ParameterizedHR(2) { }

  /**
    Destructor

    Nothing to do
  */
  virtual ~PolynomialHR() { }


  /** 
    Get the hazard rates at the specified points and time

    Pure virtual in base hazard rate class.

    @param dTime the time to get the hazard rate values at
    @param pdSpots the asset values at which to get vol values
    @param pdValues is overwritten with the hazard rate values
    @param nNumber is the number of asset values
  */
  void GetHazardRates(double dTime, const double *pdSpots, double *pdValues, 
                      size_t nNumber) const;

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
  void GetParams(double dTime, double* pdParamValues,
    double dLeftBound, double dRightBound) const;


  /** 
    Get default values for the parameters

    These default values can be used to start calibration, if a better guess
    is not known.  Space is allocated by the caller.

    @param pdDefaultParams is filled in with the default parameter values
  */
  void GetParamDefaults(double* pdDefaultParams) const
  {
    pdDefaultParams[0] = 0.0;
    pdDefaultParams[1] = 0.0;
  }


  /**
    Get the lower and upper bounds to the paramaters

    Used for calibration

    @param pdLowerBounds is filled in with the lower param bounds
    @param pdUpperBounds is filled in with the upper param bounds
  */
  void GetParamBounds(double* pdLowerBounds, double* pdUpperBounds) const
  {
    pdLowerBounds[0] = 0.0; pdUpperBounds[0] = 1.0;
    pdLowerBounds[1] = 0.0; pdUpperBounds[1] = 1.0;
  }

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
  double GetSmoothnessMeasure(double dTime, double dS);


   /**
      Dump all data of this instrument in XML format.

      Note that this method does @b  dump the contents of the associated
      volatility.

      @param tagParent the parent tag under which our tag(s) should be created

   */
 void Dump(XML::Tag& /* tagParent */) const
 {
   ASSERT_MSG(false, "Dump method not implemented for polynomialHR");
 }


};

} // namespace ihg

} // namespace ito33


#endif // #ifndef _IHG_POLYNOMIALHR_H_

