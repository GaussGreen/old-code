/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization_hrwithspotcomponentpower.h
// Purpose:     parametrization with a hr spot component power (vol given)
// Author:      Ito33
// Created:     2004/11/23
// RCS-ID:      $Id: parametrization_hrwithspotcomponentpower.h,v 1.8 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization_hrwithspotcomponentpower.h

    @brief Parametrization using a hazard rate with spot component 
    power, volatility is given.
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_HRWITHSPOTCOMPONENTPOWER_H_
#define _ITO33_IHG_PARAMETRIZATION_HRWITHSPOTCOMPONENTPOWER_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/parametrization.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{
  class ITO33_DLLDECL TermStructureCDS;
}

namespace ihg
{
  class ITO33_IHG_DLLDECL Volatility;
  class ITO33_IHG_DLLDECL HazardRateCombo;
  class ITO33_IHG_DLLDECL HRSpotComponentPower;
  class ParametrizationVisitor;

/**
   An IHG parametrization using any volatility and a hazard rate
   with spot component power.
 */
class ITO33_IHG_DLLDECL 
ParametrizationHRWithSpotComponentPower : public Parametrization
{
public:

  /// Default ctor assumes a time only hazard rate
  ParametrizationHRWithSpotComponentPower(shared_ptr<Volatility> pVolatility);

  /**
     Get the user entered volatility.

     @return the user entered volatility
   */
  shared_ptr<Volatility> GetVolatility() const
  {
    return m_pVolatility;
  }

  /**
     Calibrates with a cds term structure.

     @param tsCDS a cds term structure

     @return the calibrated hazard rate 
   */
  shared_ptr<HazardRateCombo>
    CalibrateWithCDSs(const finance::TermStructureCDS& tsCDS);

  /**
    Get the hazard spot component power calibrated to the first 
    and last cds of the term structure passed in.

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


private:

  /**
      Check the validity of the volatility entered by the user.
   */
  void CheckVolatility();

  /// The user entered volatility.
  shared_ptr<Volatility> m_pVolatility;

  /// The spot component calibrated to the 1st and last cds.
  shared_ptr<HRSpotComponentPower> m_pHRSpotComponentPower;
  
}; // class ParametrizationHRWithSpotComponentPower


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_PARAMETRIZATION_HRWITHSPOTCOMPONENTPOWER_H_

