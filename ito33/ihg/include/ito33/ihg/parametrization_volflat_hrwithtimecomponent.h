/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization_volflat_hrwithtimecomponent.h
// Purpose:     parametrization using a vol flat and a hr with time component
// Author:      Wang
// Created:     2004/06/11
// RCS-ID:      $Id: parametrization_volflat_hrwithtimecomponent.h,v 1.36 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization_volflat_hrwithtimecomponent.h
    @brief parametrization using a vol flat and a hr with time component
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_VOLFLAT_HRWITHTIMECOMPONENT_H_
#define _ITO33_IHG_PARAMETRIZATION_VOLFLAT_HRWITHTIMECOMPONENT_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/termstructurederivative.h"
#include "ito33/ihg/parametrization.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Option;
  class ITO33_DLLDECL TermStructureParBond;
  class ITO33_DLLDECL TermStructureCDS;
  class ITO33_DLLDECL TermStructureOption;
}


namespace ihg
{

  class ITO33_IHG_DLLDECL SpotComponent;
  class ITO33_IHG_DLLDECL VolatilityFlat;
  class ITO33_IHG_DLLDECL HazardRateWithTimeComponent;
  class ParametrizationVisitor;

/**
   An IHG parametrization for a constant volatility and a hazard rate
   with time component.
 */
class ITO33_IHG_DLLDECL
ParametrizationVolFlatHRWithTimeComponent : public Parametrization
{
public:

  /// Default ctor assumes a time only hazard rate
  ParametrizationVolFlatHRWithTimeComponent() { }

  /// Default dtor is ok

  /**
     Sets the spot component of the hazard rate.
 
     @param pSpotComponent the spot component of the hazard rate
   */
  void SetSpotComponent(const shared_ptr<SpotComponent>& pSpotComponent);

  /**
     Calibrates with an option and a parbond term structure.

     @param option an option to help calibrate a vol flat
     @param tsParBond a parbond term structure
   */
  void 
  CalibrateWithOptionAndParBonds
  (const finance::Option& option, const finance::TermStructureParBond& tsParBond);

  /**
     Calibrates with an option and a cds term structure.

     @param option an option to help calibrate a vol flat
     @param tsCDS a cds term structure
   */
  void CalibrateWithOptionAndCDSs(const finance::Option& option, 
                                  const finance::TermStructureCDS& tsCDS);

  /**
     Calibrates with an option and an option term structure.

     @param option an option to help calibrate a vol flat
     @param tsOption an option term structure
   */
  void CalibrateWithOptionAndOptions
       (const finance::Option& option, 
        const finance::TermStructureOption& tsOption);
  
  /**
     Calibrates with a derivative and a term structure of
     derivative

     @param deriv a derivative to help calibrate a vol flat
     @param tsDerivs a derivative term structure

     @noexport
  */
  void Calibrate(const finance::Derivative& deriv, 
                 const finance::TermStructureDerivative& tsDerivs);
  
  virtual void Calibrate(const finance::BasketGoodType& basket);
  
  /**
     Gets the calibrated flat volatility.

     @return the calibrated flat volatility
   */
  shared_ptr<VolatilityFlat> GetVolatility() const 
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

  /// The calibrated flat volatility (for output)
  shared_ptr<VolatilityFlat> m_pVolatility;

  /// The calibrated hazard rate (for output)
  shared_ptr<HazardRateWithTimeComponent> m_pHazardRate;

}; // class ParametrizationVolFlatHRWithTimeComponent


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_PARAMETRIZATION_VOLFLAT_HRWITHTIMECOMPONENT_H_
