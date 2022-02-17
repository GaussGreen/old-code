/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/volatiliyflat.h
// Purpose:     implementation of flat volatility class
// Author:      Wang
// Created:     2004/03/17
// RCS-ID:      $Id: volatilityflat.h,v 1.31 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/volatilityflat.h
 */

#ifndef _ITO33_IHG_VOLATILITYFLAT_H_
#define _ITO33_IHG_VOLATILITYFLAT_H_

#include "ito33/ihg/volatility.h"

namespace ito33
{

namespace ihg
{

 class VolatilityVisitor;

/**
    Flat (constant) volatility class.
 */
class ITO33_IHG_DLLDECL VolatilityFlat : public Volatility
{
public:

  /**
     Ctor constructs a flat volatility from a given value

     @param dValue the value of the flat volatility
   */
  VolatilityFlat(double dValue);

  /** 
     Returns the value of the flat volatility

     @return the value of the flat volatility
   */
  double GetValue() const { return m_dValue; }

  void GetVols(double dTime,
               const double *pdS, double *pdVols, size_t nNbS) const;

  /**
     Gets the squares of the volatility at a given time for an array of spot

     It's more efficient than the default implementation in the base class.

     @param dTime the time to get the vol values at
     @param pdS Array of spot.
     @param pdVolsSquared Array of the squared values of (computed) volatilities
     @param nNbS The number of spots.
   */
  void GetVolsSquared(double dTime,
                      const double *pdS, 
                      double *pdVolsSquared, 
                      size_t nNbS) const;

  shared_ptr<Volatility> Perturb(double dShift) const
  {
    return shared_ptr<Volatility>( new VolatilityFlat(m_dValue + dShift) );
  }

  void Dump(ito33::XML::Tag& tagParent) const;
  
  void Visit(VolatilityVisitor& visitor) const;

  bool IsTimeOnly() const 
  {
    return true;
  } 

 virtual
   void GetModelParameters(finance::ModelParametersConsumer& visitor) const;


private:

  /// The value of the flat volatility
  double m_dValue;

}; // class FlatVolatility


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLATILITYFLAT_H_

