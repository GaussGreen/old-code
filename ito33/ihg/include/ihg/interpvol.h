/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/interpvol.h
// Purpose:     Parameterized volatility surface class using interpolation
// Author:      David
// Created:     04.02.02
// RCS-ID:      $Id: interpvol.h,v 1.13 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/interpvol.h
 */

#ifndef _IHG_INTERPVOL_H_
#define _IHG_INTERPVOL_H_

#include "ito33/beforestd.h"
#include <list>
#include "ito33/afterstd.h"

#include "ito33/array.h"

#include "ihg/parameterizedvol.h"
#include "ihg/xml/volatility_visitor.h"

namespace ito33
{

namespace ihg
{

/** 
  Paramaterized volatility surface class using interpolation

  Based on the ideas given in Tavella, Klopfer, Implying Local Volatility,
  Wilmott Magazine, Aug. 2001. It is expected that this class will be used
  for calibration.

  In this class, a series of discrete points are saved.  Vol values
  bewteen these points are found by interpolation.

*/
class InterpVol : public ParameterizedVol
{
public:

  /** 
    Constructor

    Initialize the surface class
  */
  InterpVol(double* pdPoints, size_t nNbPoints) 
    : ParameterizedVol(nNbPoints) 
  {
    m_nNbPoints = nNbPoints;
    m_pdPoints = Array<double>(m_nNbPoints);
    for (size_t nIdx = 0; nIdx < m_nNbPoints; nIdx++)
      m_pdPoints[nIdx] = pdPoints[nIdx];  
  }

  /**
    Destructor

    Nothing to do
  */
  virtual ~InterpVol() {}
  
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
                double *Vols, size_t nNbSpots) const;


  /**
    Get the parameter values at a specific time

    Used for calibration.  The parameter values are
    simply interpolated.

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
    for (size_t nIdx = 0; nIdx < m_nNbPoints; nIdx++)
      pdDefaultParams[nIdx] = 0.2;
  }


  /**
    Get the lower and upper bounds to the paramaters

    Used for calibration

    @param pdLowerBounds is filled in with the lower param bounds
    @param pdUpperBounds is filled in with the upper param bounds
  */
  void GetParamBounds(double* pdLowerBounds, double* pdUpperBounds) const
  {
    for (size_t nIdx = 0; nIdx < m_nNbPoints; nIdx++)
    {
      pdLowerBounds[nIdx] = 0.0; 
      pdUpperBounds[nIdx] = 2.0;
    }
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

  void PrintTo(std::ostream& /* out */) const {}

  ito33::shared_ptr<Volatility> Perturb(double dShift) const
  {
    Array<double> pdTmp(m_nNbPoints);
    for (size_t nIdx = 0; nIdx < m_nNbPoints; nIdx++)
      pdTmp[nIdx] = m_pdPoints[nIdx] + dShift;

    return ito33::shared_ptr<InterpVol>(new InterpVol(pdTmp.Get(), m_nNbPoints) );
  }

  /**
      Dump all data of this instrument in XML format.

      This method is usually called by the function doing the pricing,
      calibration &c for all instruments involved at once but can also be
      called "manually" if needed.

      Note that this method does @b not dump the contents of the associated
      SessionData because when we dump several instruments using the same session
      we don't have to duplicate the common session data many times.

      @param tagParent the parent tag under which our tag(s) should be created
      @return the tag we created (so that caller could add more stuff to it
              before closing it)
      @noexport
   */
  void Dump(ito33::XML::Tag& /* tagParent */) const
  {
    ASSERT_MSG(false, "Dump method not implemented for interpvol"); 
  }


    /**
      This method is part of the implementation of the visitor pattern.

      Visitor pattern allows to do different things for different kinds of
      volatility without adding virtual functions to do all of them in the
      base volatility class itself nor doing the dreaded switchs on typeid of
      the concrete type. Instead, define a new class deriving from
      VolatilityVisitor and do whatever is required in its methods.

      @noexport
   */
  void Visit(VolatilityVisitor& /*visitor*/) const 
  {
    ASSERT_MSG(false, "visitor pattern not implemented for interpvol"); 
  }

protected:

  /// The points at which vol values are stored
  Array<double> m_pdPoints;

  /// The number of point/vol values
  size_t m_nNbPoints;


};

} // namespace finance

} // namespace ito33

#endif // #ifndef _IHG_INTERPVOL_H_

