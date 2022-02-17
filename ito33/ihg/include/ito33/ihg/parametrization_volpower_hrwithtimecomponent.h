/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization_volpower_hrwithtimecomponent.h
// Purpose:     parametrization using a vol power and a hr with time component
// Author:      ITO 33
// Created:     2005/01/05
// RCS-ID:      $Id: parametrization_volpower_hrwithtimecomponent.h,v 1.14 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization_volpower_hrwithtimecomponent.h
    @brief parametrization using a vol power and a hr with time component
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_VOLPOWER_HRWITHTIMECOMPONENT_H_
#define _ITO33_IHG_PARAMETRIZATION_VOLPOWER_HRWITHTIMECOMPONENT_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/parametrization.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Option;
  class ITO33_DLLDECL TermStructureCDS;
}


namespace ihg
{

  class ITO33_IHG_DLLDECL SpotComponent;
  class ITO33_IHG_DLLDECL VolatilityPower;
  class ITO33_IHG_DLLDECL HazardRateWithTimeComponent;
  class ParametrizationVisitor;

/**
   An IHG parametrization for a power volatility and a hazard rate
   with time component.
 */
class ITO33_IHG_DLLDECL
ParametrizationVolPowerHRWithTimeComponent : public Parametrization
{
public:

  /// Default ctor assumes a time only hazard rate
  ParametrizationVolPowerHRWithTimeComponent() { }

  /// Default dtor is ok

  /**
     Sets the spot component of the hazard rate.
 
     @param pSpotComponent the spot component of the hazard rate
   */
  void SetSpotComponent(const shared_ptr<SpotComponent>& pSpotComponent);

  /**
     Calibrates with two options and a cds term structure.

     @param option1 the first option to help calibrate a power vol
     @param option2 the second option to help calibrate a power vol
     @param tsCDS a cds term structure
   */
  void CalibrateWithOptionsAndCDSs(const finance::Option& option1, 
                                   const finance::Option& option2, 
                                   const finance::TermStructureCDS& tsCDS);

  virtual void Calibrate(const finance::BasketGoodType& basket);

  /**
     Gets the calibrated power volatility.

     @return the calibrated flat volatility
   */
  shared_ptr<VolatilityPower> GetVolatility() const { return m_pVolatility; }

  
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

  /// The calibrated power volatility (for output)
  shared_ptr<VolatilityPower> m_pVolatility;

  /// The calibrated hazard rate (for output)
  shared_ptr<HazardRateWithTimeComponent> m_pHazardRate;

}; // class ParametrizationVolPowerHRWithTimeComponent


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_PARAMETRIZATION_VOLPOWER_HRWITHTIMECOMPONENT_H_
