/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/parametrization_hrwithtimecomponent.h
// Purpose:     parametrization using a hr with time component (vol given)
// Author:      Wang
// Created:     2004/06/14
// RCS-ID:      $Id: parametrization_hrwithtimecomponent.h,v 1.31 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/parametrization_hrwithtimecomponent.h
    @brief parametrization using a hr with time component with a given vol
 */

#ifndef _ITO33_IHG_PARAMETRIZATION_HRWITHTIMECOMPONENT_H_
#define _ITO33_IHG_PARAMETRIZATION_HRWITHTIMECOMPONENT_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/parametrization.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL TermStructureCDS;
  class ITO33_DLLDECL TermStructureParBond;
  class ITO33_DLLDECL TermStructureEDS;
  class ITO33_DLLDECL TermStructureOption;

  class ITO33_DLLDECL Derivative;
  template<class T> class TermStructure;
  typedef TermStructure<Derivative> TermStructureDerivative;
}

namespace XML
{
  class Tag;
}

namespace ihg
{
  class ITO33_IHG_DLLDECL SpotComponent;
  class ITO33_IHG_DLLDECL Volatility;
  class ITO33_IHG_DLLDECL HazardRateWithTimeComponent;
  class ParametrizationVisitor;


/**
   An IHG parametrization using any volatility for a hazard rate
   with time component.
 */
class ITO33_IHG_DLLDECL 
ParametrizationHRWithTimeComponent : public Parametrization
{
public:

  /// Default ctor assumes a time only hazard rate
  ParametrizationHRWithTimeComponent() { }

  /// Default dtor is ok

  /**
     Set the user entered volatility.

     @param pVolatility the user entered volatility
   */
  void SetVolatility(shared_ptr<Volatility> pVolatility)
  {
    m_pVolatility = pVolatility;
  }

  /**
     Get the user entered volatility.

     @return the user entered volatility
   */
  shared_ptr<Volatility> GetVolatility() const
  {
    return m_pVolatility;
  }

  /**
     Sets the spot component of the hazard rate.
 
     @param pSpotComponent the spot component of the hazard rate
   */
  void SetSpotComponent(const shared_ptr<SpotComponent>& pSpotComponent);

  /**
     Calibrates with a cds term structure. The calibrated hazard rate 
     will be time only if spot component is not set.

     @param tsCDS a cds term structure

     @return the calibrated hazard rate with a time component. If spot component
             is not set, it will be time only
   */
  shared_ptr<HazardRateWithTimeComponent>
  CalibrateWithCDSs(const finance::TermStructureCDS& tsCDS);

  /**
     Calibrates with a parbond term structure. The calibrated hazard rate 
     will be time only if spot component is not set.

     @param tsParBond a parbond term structure

     @return the calibrated hazard rate with a time component. If spot component
             is not set, it will be time only
   */
  shared_ptr<HazardRateWithTimeComponent>
  CalibrateWithParBonds(const finance::TermStructureParBond& tsParBond);

  /**
     Calibrates with an eds term structure. The calibrated hazard rate 
     will be time only if spot component is not set.

     @param tsEDS an eds term structure

     @return the calibrated hazard rate with a time component. If spot component
             is not set, it will be time only
   */
  shared_ptr<HazardRateWithTimeComponent>
  CalibrateWithEDSs(const finance::TermStructureEDS& tsEDS);

  /**
     Calibrates with an option term structure, the volatility has to be set. 

     @param tsOption an option term structure

     @return the calibrated hazard rate with a time component.
   */
  shared_ptr<HazardRateWithTimeComponent>
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

  /// Valid for CDS and ParBond etc, but not for EDS
  shared_ptr<HazardRateWithTimeComponent>
    CalibrateWithTermStructurePossibleTimeOnly
    (const finance::TermStructureDerivative& tsDerivs);

private:
  /**
      Check the validity of the volatility entered by the user
   */
  void CheckVolatility();

  /// The user entered spot component
  shared_ptr<SpotComponent> m_pSpotComponent;

  /// The user entered volatility
  shared_ptr<Volatility> m_pVolatility;

}; // class ParametrizationHRWithTimeComponent


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_PARAMETRIZATION_HRWITHTIMECOMPONENT_H_
