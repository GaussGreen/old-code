/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization_volpower.h
// Purpose:     parametrization using a power volatility
// Author:      Ito33
// Created:     2004/12/14
// RCS-ID:      $Id: parametrization_volpower.h,v 1.7 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization_volpower.h
    @brief parametrization using a power volatility
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_VOLPOWER_H_
#define _ITO33_IHG_PARAMETRIZATION_VOLPOWER_H_

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

  class ITO33_IHG_DLLDECL VolatilityPower;
  class ITO33_IHG_DLLDECL HazardRate;
  class ParametrizationVisitor;

/// An IHG parametrization for a power volatility
class ITO33_IHG_DLLDECL ParametrizationVolPower : public Parametrization
{
public:

  /**
      Constructor.

      @param pHR The hazard rate for pricing
   */
  ParametrizationVolPower(const shared_ptr<HazardRate>& pHR);

  // Default dtor is ok

  /**
      Gets the user entered hazard rate.

      @return the user entered hazard rate
   */
  shared_ptr<HazardRate> GetHazardRate() const
  {
    return m_pHazardRate;
  }

  /**
      Calibrates power volatility to two options.

      @param option1 first option to help calibrate a power vol
      @param option2 second option to help calibrate a power vol
   */
  shared_ptr<VolatilityPower> 
  CalibrateWithOptions(const finance::Option& option1, 
                       const finance::Option& option2);

    
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

  /// The input hazard rate 
  shared_ptr<HazardRate> m_pHazardRate;

}; // class ParametrizationVolPower


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_PARAMETRIZATION_VOLPOWER_H_
