/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization_volwithtimecomponent.h
// Purpose:     parametrization using a vol with time component (hazard rate given)
// Created:     2005/03/04
// RCS-ID:      $Id: parametrization_volwithtimecomponent.h,v 1.7 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization_volwithtimecomponent.h
    @brief parametrization using a vol with time component, hazard rate given.
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_VOLWITHTIMECOMPONENT_H_
#define _ITO33_IHG_PARAMETRIZATION_VOLWITHTIMECOMPONENT_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/parametrization.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL TermStructureOption;
}

namespace XML
{
  class Tag;
}

namespace ihg
{
  class ITO33_IHG_DLLDECL SpotComponent;
  class ITO33_IHG_DLLDECL HazardRate;
  class ITO33_IHG_DLLDECL VolatilityWithTimeComponent;
  class ParametrizationVisitor;


/**
   An IHG parametrization using any hazard rate for a volatility
   with time component.
 */
class ITO33_IHG_DLLDECL 
ParametrizationVolWithTimeComponent : public Parametrization
{
public:

  /**
     Default constructors, the spot component and the time component is
     combined by multiplication.

     @param pHR the given hazard rate.
   */
  ParametrizationVolWithTimeComponent(const shared_ptr<HazardRate>& pHR);

  /// Default dtor is ok

  /**
     Get the user entered hazard rate.

     @return the user entered hazard rate
   */
  shared_ptr<HazardRate> GetHazardRate() const { return m_pHazardRate; }

  /**
     Sets the spot component of the volatility.
 
     @param pSpotComponent the spot component of the volatility
   */
  void SetSpotComponent(const shared_ptr<SpotComponent>& pSpotComponent);

  /**
     Calibrates with an option term structure.

     @param tsOption an option term structure

     @return the calibrated volatility with a time component.
   */
  shared_ptr<VolatilityWithTimeComponent>
  CalibrateWithOptions(const finance::TermStructureOption& tsOption);

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


private:

  /// The user entered spot component of volatility
  shared_ptr<SpotComponent> m_pSpotComponent;

  /// The user entered volatility
  shared_ptr<HazardRate> m_pHazardRate; 

}; // class ParametrizationVolWithTimeComponent


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_PARAMETRIZATION_VOLWITHTIMECOMPONENT_H_
