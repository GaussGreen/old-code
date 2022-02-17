/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/volatiliycombo.h
// Purpose:     Combination of two volatility classes
// Author:      David
// Created:     2004/06/02
// RCS-ID:      $Id: volatilitycombo.h,v 1.14 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/volatilitycombo.h

    Combination of time and spot components volatility using multiplication.  
 */

#ifndef _ITO33_IHG_VOLATILITYCOMBO_H_
#define _ITO33_IHG_VOLATILITYCOMBO_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/volatilitywithtimecomponent.h"

namespace ito33
{

namespace ihg
{

  class ITO33_IHG_DLLDECL SpotComponent;

  class VolatilityVisitor;

/**
    Volatility consisting of a time component and a spot component combined by 
    multiplication.

    @nocreate
 */
class ITO33_IHG_DLLDECL VolatilityCombo : public VolatilityWithTimeComponent
{
public:

  /**
     @internal
     @brief Ctor constructs a volatility using a spot component. 

     @param pSpotComponent the spot component of the volatility
     
     @noexport
   */
  VolatilityCombo(shared_ptr<SpotComponent> pSpotComponent);

  /**
     @internal
     @brief Ctor constructs a volatility using a spot component and a 
     piecewise constant time component.

     @param pSpotComponent the spot component of the volatility
     @param pDates the time points of the time component
     @param pdValues the values of the time component
     @param nNbTimes the number of values of the time component

     @noexport
   */
  VolatilityCombo(shared_ptr<SpotComponent> pSpotComponent,
                  const Date* pDates, const double* pdValues, size_t nNbTimes);

  /**
     Ctor constructs a volatility using a spot component and a 
     piecewise constant time component.

     @param pSpotComponent the spot component of the volatility
     @param dates the time points of the time component
     @param values the values of the time component
   */
  VolatilityCombo(shared_ptr<SpotComponent> pSpotComponent,
                  const std::vector<Date>& dates, 
                  const std::vector<double>& values);

  void GetVols(double dTime,
               const double *pdS, double *pdVols, size_t nNbS) const;

  void Dump(XML::Tag& tagParent) const;
  
  void Visit(VolatilityVisitor& visitor) const;


private:

  /// The spot component of the volatility
  shared_ptr<SpotComponent> m_pSpotComponent;

}; // class VolatilityCombo


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLATILITYCOMBO_H_
