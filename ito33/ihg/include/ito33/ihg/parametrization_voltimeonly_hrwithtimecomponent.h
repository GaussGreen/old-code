/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization_voltimeonly_hrwithtimecomponent.h
// Purpose:     parametrization using a vol time only and a hr with time component
// Author:      ITO 33
// Created:     2005/07/14
// RCS-ID:      $Id: parametrization_voltimeonly_hrwithtimecomponent.h,v 1.8 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization_voltimeonly_hrwithtimecomponent.h
    @brief parametrization using a vol time only and a hr with time component
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_VOLTIMEONLY_HRWITHTIMECOMPONENT_H_
#define _ITO33_IHG_PARAMETRIZATION_VOLTIMEONLY_HRWITHTIMECOMPONENT_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/parametrization.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivatives;
}


namespace ihg
{

  class ITO33_IHG_DLLDECL SpotComponent;
  class ITO33_IHG_DLLDECL VolatilityTimeOnly;
  class ITO33_IHG_DLLDECL HazardRateWithTimeComponent;
  class ParametrizationVisitor;

/**
   An IHG parametrization for a time only volatility and a hazard rate
   with time component.
 */
class ITO33_IHG_DLLDECL
ParametrizationVolTimeOnlyHRWithTimeComponent : public Parametrization
{
public:

  /// Default ctor assumes a time only hazard rate
  ParametrizationVolTimeOnlyHRWithTimeComponent() { }

  /// Default dtor is ok

  /**
     Sets the spot component of the hazard rate.
 
     @param pSpotComponent the spot component of the hazard rate
   */
  void SetSpotComponent(const shared_ptr<SpotComponent>& pSpotComponent);

  /**
     Calibrates with a general (but even number of) derivative list.

     @param derivatives the list of derivatives to calibrate
   */
  void CalibrateWithDerivatives(const finance::Derivatives& derivatives);

  virtual void Calibrate(const finance::BasketGoodType& basket);

  /**
     Gets the calibrated time only volatility.

     @return the calibrated time only volatility
   */
  shared_ptr<VolatilityTimeOnly> GetVolatility() const 
  { 
    return m_pVolatility; 
  }
  
  /**
     Gets the calibrated hazard rate with time component. It will be time 
     only if spot component is not set.

     @return the calibrated hazard rate with time component. 
   */
  shared_ptr<HazardRateWithTimeComponent> GetHazardRate() const
  {
    return m_pHazardRate;
  }

  /**
      Dumps all data stored in this object in XML format.

      @param tagParent the parent tag under which our tag(s) should be created
      @noexport
   */
  void Dump(ito33::XML::Tag& tagParent) const;

  /**
      Support for visitor pattern when reading from XML.

      @param visitor the parametrization visitor object
      @noexport
   */
  void Visit(ParametrizationVisitor& visitor) const;

  virtual shared_ptr<finance::TheoreticalModel> GetTheoreticalModel();

private:

  /// The user entered spot component
  shared_ptr<SpotComponent> m_pSpotComponent;

  /// The calibrated time only volatility (for output)
  shared_ptr<VolatilityTimeOnly> m_pVolatility;

  /// The calibrated hazard rate (for output)
  shared_ptr<HazardRateWithTimeComponent> m_pHazardRate;

}; // class ParametrizationVolTimeOnlyHRWithTimeComponent


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_PARAMETRIZATION_VOLTIMEONLY_HRWITHTIMECOMPONENT_H_
