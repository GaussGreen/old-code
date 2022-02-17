/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization_voltanh.h
// Purpose:     parametrization using a tanh volatility
// Created:     2005/02/04
// RCS-ID:      $Id: parametrization_voltanh.h,v 1.4 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization_voltanh.h
    @brief parametrization using a tanh volatility
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_VOLTANH_H_
#define _ITO33_IHG_PARAMETRIZATION_VOLTANH_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/parametrization.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Option;
}

namespace ihg
{

  class ITO33_IHG_DLLDECL VolatilityTanh;
  class ITO33_IHG_DLLDECL HazardRate;

/**
   An IHG parametrization for a tanh volatility 
 */
class ITO33_IHG_DLLDECL 
ParametrizationVolTanh : public Parametrization
{
public:

  /**
     Constructor.

     @param pHR The hazard rate for pricing
   */
  ParametrizationVolTanh(const shared_ptr<HazardRate>& pHR);

  /// Default dtor is ok

  /**
     Get the user entered hazard rate.

     @return the user entered hazard rate
   */
  shared_ptr<HazardRate> GetHazardRate() const
  {
    return m_pHazardRate;
  }

  /**
     Calibrates tanh volatility to two options.

     @param option1 first option to help calibrate an tanh vol
     @param option2 second option to help calibrate an tanh vol
   */
  shared_ptr<VolatilityTanh> 
  CalibrateWithOptions(const finance::Option& option1, 
                       const finance::Option& option2);
    
  void Dump(ito33::XML::Tag& tagParent) const;

  void Visit(ParametrizationVisitor& visitor) const;


private:

  /// The input hazard rate 
  shared_ptr<HazardRate> m_pHazardRate;

}; // class ParametrizationVolTanh


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_PARAMETRIZATION_VOLTANH_H_
