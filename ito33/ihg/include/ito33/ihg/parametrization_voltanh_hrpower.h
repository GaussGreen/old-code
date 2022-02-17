/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization_voltanh_hrpower.h
// Purpose:     parametrization using a tanh vol and a power hr
// Author:      ITO 33
// Created:     2005/07/29
// RCS-ID:      $Id: parametrization_voltanh_hrpower.h,v 1.7 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization_voltanh_hrpower.h
    @brief parametrization using a tanh vol and a power hr
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_VOLTANH_HRPOWER_H_
#define _ITO33_IHG_PARAMETRIZATION_VOLTANH_HRPOWER_H_

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
  
  class ITO33_IHG_DLLDECL VolatilityTanh;
  class ITO33_IHG_DLLDECL HazardRatePower;
  class ParametrizationVisitor;

/**
   An IHG parametrization for a flat volatility and power hazard rate.
 */
class ITO33_IHG_DLLDECL
ParametrizationVolTanhHRPower : public Parametrization
{
public:

  /// Empty ctor 
  ParametrizationVolTanhHRPower() { }

  /// Default dtor is ok


  /**
     Calibrates with a general derivative list.

     @param derivatives the list of derivatives to calibrate
   */
  void CalibrateWithDerivatives(const finance::Derivatives& derivatives);

  virtual void Calibrate(const finance::BasketGoodType& basket);

  /**
     Gets the calibrated tanh volatility.

     @return the calibrated tanh volatility
   */
  shared_ptr<VolatilityTanh> GetVolatility() const 
  { 
    return m_pVolatility; 
  }
  
  /**
     Gets the calibrated power hazard rate.

     @return the calibrated power hazard rate. 
   */
  shared_ptr<HazardRatePower> GetHazardRate() const
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

  /// The calibrated tanh volatility (for output)
  shared_ptr<VolatilityTanh> m_pVolatility;

  /// The calibrated power hazard rate (for output)
  shared_ptr<HazardRatePower> m_pHazardRate;

}; // class ParametrizationVolTanhHRPower


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_PARAMETRIZATION_VOLTANH_HRPOWER_H_
