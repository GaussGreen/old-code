/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization_volpower_hrpower.h
// Purpose:     parametrization using a power volatility and spot component
//              power hazard rate
// Author:      Ito33
// Created:     2005/01/03
// RCS-ID:      $Id: parametrization_volpower_hrpower.h,v 1.11 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization_volpower_hrpower.h
    @brief parametrization using a power volatility and spot component
           power hazard rate
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_VOLPOWER_HRPOWER_H_
#define _ITO33_IHG_PARAMETRIZATION_VOLPOWER_HRPOWER_H_

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

  class ITO33_IHG_DLLDECL VolatilityPower;
  class ITO33_IHG_DLLDECL SpotComponent;
  class ITO33_IHG_DLLDECL HazardRateCombo;
  class ITO33_IHG_DLLDECL HRSpotComponentPower;
  class ParametrizationVisitor;


/**
   An IHG parametrization for a power volatility and spot component
   power hazard rate
 */
class ITO33_IHG_DLLDECL ParametrizationVolPowerHRPower : public Parametrization
{
public:

  /**
     Constructor.
   */
  ParametrizationVolPowerHRPower() {}

  /// Default dtor is ok


  /**
     Calibrates power volatility and spot component power hazard
     rate to two options and a CDS term structure

     @param option1 first option 
     @param option2 second option 
     @param tsCDS the CDS term structure
   */
  void CalibrateWithOptionsAndCDSs(const finance::Option &option1, 
                                   const finance::Option &option2,
                                   const finance::TermStructureCDS& tsCDS);

  virtual void Calibrate(const finance::BasketGoodType& basket);

  /**
     Gets the calibrated power volatility.

     @return the calibrated power volatility
   */
  shared_ptr<VolatilityPower> GetVolatility() const 
  { 
    return m_pVolatility; 
  }

  
  /**
     Gets the calibrated hazard rate with spot and time component.

     @return the calibrated hazard rate with spot and time component. 
   */
  shared_ptr<HazardRateCombo> GetHazardRate() const
  {
    return m_pHazardRate;
  }

  /*
    Gets the hazard spot component power calibrated to the first and 
    last cds of the term structure passed in

    @return the calibrated spot component to the 1st and last cds
  */
  shared_ptr<HRSpotComponentPower> GetHRSpotComponentPower() const
  {
    return m_pHRSpotComponentPower;
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

  /// The calibrated volatility (for output)
  shared_ptr<VolatilityPower> m_pVolatility;
  
  /// The calibrated hazard rate (for output)
  shared_ptr<HazardRateCombo> m_pHazardRate;
  
  /// The calibrated spot component (for output)
  shared_ptr<HRSpotComponentPower> m_pHRSpotComponentPower;

}; // class ParametrizationVolPower


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_PARAMETRIZATION_VOLPOWER_HRPOWER_H_

