/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization_voltanh_hrwithspotcomponentpower.h
// Purpose:     parametrization using a tanh volatility and spot component
//              power hazard rate
// Created:     2005/02/07
// RCS-ID:      $Id: parametrization_voltanh_hrwithspotcomponentpower.h,v 1.7 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization_voltanh_hrwithspotcomponentpower.h
    @brief parametrization using a tanh volatility and spot component
           power hazard rate
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_VOLTANH_HRWITHSPOTCOMPONENTPOWER_H_
#define _ITO33_IHG_PARAMETRIZATION_VOLTANH_HRWITHSPOTCOMPONENTPOWER_H_

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

  class ITO33_IHG_DLLDECL VolatilityTanh;
  class ITO33_IHG_DLLDECL SpotComponent;
  class ITO33_IHG_DLLDECL HazardRateCombo;
  class ITO33_IHG_DLLDECL HRSpotComponentPower;
  class ParametrizationVisitor;


/**
   An IHG parametrization for a tanh volatility and spot component
   power hazard rate
 */
class ITO33_IHG_DLLDECL
ParametrizationVolTanhHRWithSpotComponentPower : public Parametrization
{
public:

  /**
     Constructor.
   */
  ParametrizationVolTanhHRWithSpotComponentPower() {}

  /// Default dtor is ok


  /**
     Calibrates tanh volatility and spot component tanh hazard
     rate to two options and a CDS term structure

     @param option1 first option 
     @param option2 second option 
     @param tsCDS the CDS term structure
   */
  void CalibrateWithOptionsAndCDSs(const finance::Option &option1, 
                                   const finance::Option &option2,
                                   const finance::TermStructureCDS& tsCDS);

 /**
     Gets the calibrated tanh volatility.

     @return the calibrated tanh volatility
   */
  shared_ptr<VolatilityTanh> GetVolatility() const 
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
   
  void Dump(ito33::XML::Tag& tagParent) const;

  void Visit(ParametrizationVisitor& visitor) const;


private:

  /// The calibrated volatility (for output)
  shared_ptr<VolatilityTanh> m_pVolatility;
  
  /// The calibrated hazard rate (for output)
  shared_ptr<HazardRateCombo> m_pHazardRate;
  
  /// The calibrated spot component (for output)
  shared_ptr<HRSpotComponentPower> m_pHRSpotComponentPower;

}; // class ParametrizationVolTanhHRWithSpotComponentPower


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_PARAMETRIZATION_VOLTANH_HRWITHSPOTCOMPONENTPOWER_H_

