/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardratecombo.h
// Purpose:     hazard rate combined by a spot and a time components
// Author:      David
// Created:     2004/06/02
// RCS-ID:      $Id: hazardratecombo.h,v 1.39 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/ihg/hazardratecombo.h
   @brief Hazard rate combied(Addition or multiplication) by a spot component
          and a piecewise constant time component.

   See also ito33/ihg/VolatilityCombo.h
 */

#ifndef _ITO33_IHG_HAZARDRATECOMBO_H_
#define _ITO33_IHG_HAZARDRATECOMBO_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/hazardratewithtimecomponent.h"

#include "ito33/ihg/common.h"


namespace ito33
{

namespace XML
{
  class Tag;
}

namespace ihg
{

  class HazardRateVisitor;
  class ITO33_IHG_DLLDECL SpotComponent;

/**
    Hazard rate given by multiplication of a spot component and a piecewise
    constant time component.
 */
class ITO33_IHG_DLLDECL HazardRateCombo : public HazardRateWithTimeComponent
{

public:

  /**
     @internal
     @brief Ctor constructs a hazard rate using a spot component.

     @param pSpotComponent the spot component of the hazard rate

     @noexport
   */
  HazardRateCombo(const shared_ptr<SpotComponent>& pSpotComponent)
                : HazardRateWithTimeComponent(),
                  m_pSpotComponent(pSpotComponent)
  {
  }

  /**
     @internal
     @brief Ctor constructs a hazard rate using a spot component and a
     piecewise constant time component.

     @param pSpotComponent the spot component of the hazard rate
     @param pDates the time points of the time component
     @param pdValues the values of the time component
     @param nNbTimes the number of values of the time component

     @noexport
   */
  HazardRateCombo(const shared_ptr<SpotComponent>& pSpotComponent,
                  const Date* pDates, const double* pdValues, size_t nNbTimes);

  /**
     Ctor constructs a hazard rate using a spot component and a
     piecewise constant time component.

     @param pSpotComponent the spot component of the hazard rate
     @param dates the time points of the time component
     @param values the values of the time component
   */
  HazardRateCombo(const shared_ptr<SpotComponent>& pSpotComponent,
                  const std::vector<Date>& dates,
                  const std::vector<double>& values);

  /// Gets the spot component
  shared_ptr<SpotComponent> GetSpotComponent() const
  {
    return m_pSpotComponent;
  }

  // Default dtor is ok

  void GetHazardRates(double dTime, const double *pdS, double *pdValues,
                      size_t nNbS) const;

  virtual void Dump(XML::Tag& tagParent) const;

  virtual void Visit(HazardRateVisitor& visitor) const;

  virtual void 
  GetModelParameters(finance::ModelParametersConsumer& visitor) const;

private:

  /// The spot component of the hazard rate
  shared_ptr<SpotComponent> m_pSpotComponent;

}; // class HazardRateCombo


} // namespace ihg

} // namespace ito33


#endif // #ifndef _ITO33_IHG_HAZARDRATECOMBO_H_

