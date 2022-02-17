/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/polynomialvol.h
// Purpose:     Parameterized volatility surface class using polynomials
// Author:      David
// Created:     04.01.24
// RCS-ID:      $Id: polynomialvol.h,v 1.9 2004/10/04 18:04:04 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/polynomialvol.h
 */

#ifndef _IHG_POLYNOMIALVOL_H_
#define _IHG_POLYNOMIALVOL_H_

#include "ito33/beforestd.h"
#include <list>
#include "ito33/afterstd.h"

#include <nag.h>

#include "ito33/array.h"

#include "ihg/parameterizedvol.h"

namespace ito33
{

namespace ihg
{

/** 
  Paramaterized volatility surface class using polynomials

  Based on the equations given in Tavella, Klopfer, Implying Local Volatility,
  Wilmott Magazine, Aug. 2001. It is expected that this class will be used
  for calibration.

  The actual equations from Tavella are:
  vol(S,t) = c1 + (c2 / S^c5) + (c3 * S^C4)

*/
class PolynomialVol : public ParameterizedVol
{
public:

  /** 
    Constructor

    Initialize the surface class
  */
  PolynomialVol() : ParameterizedVol(5) { }

  /**
    Destructor

    Nothing to do
  */
  virtual ~PolynomialVol() {}
  
  /** 
    Average volatility value

    Used to determine grid sizes.  Pure virtual in base class

    @param dSpotInitial is the place at which to compute an average value
  */
  double AverageValue(double dSpotInitial);


  /** 
    Get the volatilities at the specified points and time

    Pure virtual in base class.

    @param dTime the time to get the vol values at
    @param Spots the asset values at which to get vol values
    @param Vols is overwritten with the vol values
    @param nNbSpots is the number of asset values
  */
  void GetVols(double dTime, const double *Spots,
                double *Vols, size_t nNbSpots);

  
  /**
    Get the parameter values at a specific time



    @param dTime the time at which the parameters are requested
    @pdParamValues is filled in with the actual paramater values
  */
  void GetParams(double dTime, double* pdParamValues) const;

  /**
    Get the parameter values at a specific time

    Used for calibration.  Cannot simply average (interpolate) the values
    at adjacent times, since polynomials are used.   Use NAG to do a least
    squares fit

    @param dTime the time at which the parameters are requested
    @param pdParamValues is filled in with the actual paramater values
    @param dLeftBound is a left bound for which the returned params should be valid
    @param dRightBound is a right bound for which the returned params should be valid
  */
  void GetParams(double dTime, double* pdParamValues,
    double dLeftBound, double dRightBound);


  /** 
    Get default values for the parameters

    These default values can be used to start calibration, if a better guess
    is not known.  Space is allocated by the caller.

    @param pdDefaultParams is filled in with the default parameter values
  */
  void GetParamDefaults(double* pdDefaultParams) const
  {
    pdDefaultParams[0] = 0.2;
    pdDefaultParams[1] = 0.0;
    pdDefaultParams[2] = 0.0;
    pdDefaultParams[3] = 0.0;
    pdDefaultParams[4] = 0.0;
  }


  /**
    Get the lower and upper bounds to the paramaters

    Used for calibration

    @param pdLowerBounds is filled in with the lower param bounds
    @param pdUpperBounds is filled in with the upper param bounds
  */
  void GetParamBounds(double* pdLowerBounds, double* pdUpperBounds) const
  {
    pdLowerBounds[0] = -5.0; pdUpperBounds[0] = 5.0;
    pdLowerBounds[1] = -100.0; pdUpperBounds[1] = 100.0;
    pdLowerBounds[2] = -100.0; pdUpperBounds[2] = 100.0;
    pdLowerBounds[3] = 0.0; pdUpperBounds[3] = 0.7;
    pdLowerBounds[4] = 0.0; pdUpperBounds[4] = 0.7;
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

  /// The fitting objective function
  void __stdcall ParamFitObjFun(Integer iN, double pdX[], double* objf, double pdG[], Nag_Comm* comm);

protected:
  /// The number of vol values used during fitting for GetParams
  size_t m_nNbPointsToFit;

  /// Vol values to fit during a GetParam() call
  Array<double> m_pdVolsToFit;

  /// The points for vol values to fit during a GetParam() call
  Array<double> m_pdPointsToFit;

};

} // namespace finance

} // namespace ito33

#endif // #ifndef _IHG_POLYNOMIALVOL_H_

