/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/xml/parametrization_visitor_goodtype.h
// Purpose:     Visitor for parametrization classes
// Created:     2005/01/27
// RCS-ID:      $Id: parametrization_visitor_goodtype.h,v 1.9 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/xml/parametrization_visitor_goodtype.h
    @brief Visitor class accepting different parametrizations.

    A complete implementation of the visitor pattern for all the 
    parametrizations. 
 */

#ifndef _ITO33_IHG_XML_PARAMETRIZATION_VISITOR_GOODTYPE_H_
#define _ITO33_IHG_XML_PARAMETRIZATION_VISITOR_GOODTYPE_H_

#include "ito33/ihg/parametrization_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_hrwithspotcomponentpower.h"
#include "ito33/ihg/parametrization_volflat_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_volflat_hrwithspotcomponentpower.h"
#include "ito33/ihg/parametrization_volflat_hrpower.h"
#include "ito33/ihg/parametrization_volpower.h"
#include "ito33/ihg/parametrization_volpower_hrpower.h"
#include "ito33/ihg/parametrization_volpower_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_voltanh.h"
#include "ito33/ihg/parametrization_voltanh_hrwithspotcomponentpower.h"
#include "ito33/ihg/parametrization_voltanh_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_voltanh_hrpower.h"
#include "ito33/ihg/parametrization_voltimeonly_hrwithtimecomponent.h"

#include "ito33/ihg/parametrization_volwithtimecomponent.h"

#include "ihg/xml/parametrization_visitor.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace ihg
{


/// Parametrization visitor for all the parametrizations.
class ITO33_DLLDECL ParametrizationVisitorGoodType : public ParametrizationVisitor
{
public:

  ParametrizationVisitorGoodType() {}

  /// parametrization with hazard rate with time component
  void OnParametrizationHRWithTimeComponent
       (const ParametrizationHRWithTimeComponent& param)
  {
    m_pParametrizationHRWithTimeComponent = 
      shared_ptr<ParametrizationHRWithTimeComponent>
      (new ParametrizationHRWithTimeComponent(param));
  }

  shared_ptr<ParametrizationHRWithTimeComponent>
  GetParametrizationHRWithTimeComponent() 
  { 
    return m_pParametrizationHRWithTimeComponent; 
  }


  /// parametrization with hazard rate with spot component power
  void OnParametrizationHRWithSpotComponentPower
       (const ParametrizationHRWithSpotComponentPower& param)
  {
    m_pParametrizationHRWithSpotComponentPower = 
      shared_ptr<ParametrizationHRWithSpotComponentPower>
      (new ParametrizationHRWithSpotComponentPower(param));
  }

  shared_ptr<ParametrizationHRWithSpotComponentPower>
  GetParametrizationHRWithSpotComponentPower() 
  { 
    return m_pParametrizationHRWithSpotComponentPower; 
  }

  /// parametrization vol flat with hazard rate with time component
  void OnParametrizationVolFlatHRWithTimeComponent
       (const ParametrizationVolFlatHRWithTimeComponent& param)
  {
    m_pParametrizationVolFlatHRWithTimeComponent = 
      shared_ptr<ParametrizationVolFlatHRWithTimeComponent>
      (new ParametrizationVolFlatHRWithTimeComponent(param));
  }

  shared_ptr<ParametrizationVolFlatHRWithTimeComponent>
  GetParametrizationVolFlatHRWithTimeComponent() 
  { 
    return m_pParametrizationVolFlatHRWithTimeComponent; 
  }

  /// parametrization vol flat with hazard rate with spot component power
  void OnParametrizationVolFlatHRWithSpotComponentPower
       (const ParametrizationVolFlatHRWithSpotComponentPower& param)
  {
    m_pParametrizationVolFlatHRWithSpotComponentPower = 
      shared_ptr<ParametrizationVolFlatHRWithSpotComponentPower>
      (new ParametrizationVolFlatHRWithSpotComponentPower(param));
  }

  shared_ptr<ParametrizationVolFlatHRWithSpotComponentPower>
  GetParametrizationVolFlatHRWithSpotComponentPower() 
  { 
    return m_pParametrizationVolFlatHRWithSpotComponentPower; 
  }

  /// parametrization vol flat with hazard rate power
  void OnParametrizationVolFlatHRPower
       (const ParametrizationVolFlatHRPower& param)
  {
    m_pParametrizationVolFlatHRPower = 
      shared_ptr<ParametrizationVolFlatHRPower>
      (new ParametrizationVolFlatHRPower(param));
  }

  shared_ptr<ParametrizationVolFlatHRPower>
  GetParametrizationVolFlatHRPower() 
  { 
    return m_pParametrizationVolFlatHRPower; 
  }

  /// parametrization vol power
  void OnParametrizationVolPower(const ParametrizationVolPower& param)
  {
    m_pParametrizationVolPower = 
      shared_ptr<ParametrizationVolPower> (new ParametrizationVolPower(param) );
  }

  shared_ptr<ParametrizationVolPower> GetParametrizationVolPower() 
  { 
    return m_pParametrizationVolPower; 
  }

  /// parametrization vol power hr power
  void 
  OnParametrizationVolPowerHRPower(const ParametrizationVolPowerHRPower& param)
  {
    m_pParametrizationVolPowerHRPower = 
      shared_ptr<ParametrizationVolPowerHRPower> 
      (new ParametrizationVolPowerHRPower(param) );
  }

  shared_ptr<ParametrizationVolPowerHRPower> GetParametrizationVolPowerHRPower() 
  { 
    return m_pParametrizationVolPowerHRPower; 
  }

  /// parametrization vol power with hazard rate with time component
  void OnParametrizationVolPowerHRWithTimeComponent
       (const ParametrizationVolPowerHRWithTimeComponent& param)
  {
    m_pParametrizationVolPowerHRWithTimeComponent = 
      shared_ptr<ParametrizationVolPowerHRWithTimeComponent>
      (new ParametrizationVolPowerHRWithTimeComponent(param));
  }

  shared_ptr<ParametrizationVolPowerHRWithTimeComponent>
  GetParametrizationVolPowerHRWithTimeComponent() 
  { 
    return m_pParametrizationVolPowerHRWithTimeComponent; 
  }

  /// parametrization vol tanh
  void OnParametrizationVolTanh
       (const ParametrizationVolTanh& param)
  {
    m_pParametrizationVolTanh = 
      shared_ptr<ParametrizationVolTanh> (new ParametrizationVolTanh(param));
  }

  shared_ptr<ParametrizationVolTanh> GetParametrizationVolTanh() 
  { 
    return m_pParametrizationVolTanh; 
  }

  /// parametrization vol tanh with hazard rate with time component
  void OnParametrizationVolTanhHRWithTimeComponent
       (const ParametrizationVolTanhHRWithTimeComponent& param)
  {
    m_pParametrizationVolTanhHRWithTimeComponent = 
      shared_ptr<ParametrizationVolTanhHRWithTimeComponent>
      (new ParametrizationVolTanhHRWithTimeComponent(param));
  }

  shared_ptr<ParametrizationVolTanhHRWithTimeComponent>
  GetParametrizationVolTanhHRWithTimeComponent() 
  { 
    return m_pParametrizationVolTanhHRWithTimeComponent; 
  }

  /// parametrization vol tanh with hazard rate with spot component power
  void OnParametrizationVolTanhHRWithSpotComponentPower
       (const ParametrizationVolTanhHRWithSpotComponentPower& param)
  {
    m_pParametrizationVolTanhHRWithSpotComponentPower = 
      shared_ptr<ParametrizationVolTanhHRWithSpotComponentPower>
      (new ParametrizationVolTanhHRWithSpotComponentPower(param));
  }

  shared_ptr<ParametrizationVolTanhHRWithSpotComponentPower>
  GetParametrizationVolTanhHRWithSpotComponentPower() 
  { 
    return m_pParametrizationVolTanhHRWithSpotComponentPower; 
  }

  /// parametrization vol tanh with power hazard rate
  void OnParametrizationVolTanhHRPower
       (const ParametrizationVolTanhHRPower& param)
  {
    m_pParametrizationVolTanhHRPower = 
      shared_ptr<ParametrizationVolTanhHRPower>
      (new ParametrizationVolTanhHRPower(param));
  }

  shared_ptr<ParametrizationVolTanhHRPower>
  GetParametrizationVolTanhHRPower() 
  { 
    return m_pParametrizationVolTanhHRPower; 
  }

  /// parametrization vol with time component
  void OnParametrizationVolWithTimeComponent
       (const ParametrizationVolWithTimeComponent& param)
  {
    m_pParametrizationVolWithTimeComponent 
         = make_ptr( new ParametrizationVolWithTimeComponent(param) );
  }

  shared_ptr<ParametrizationVolWithTimeComponent> 
  GetParametrizationVolWithTimeComponent() 
  { 
    return m_pParametrizationVolWithTimeComponent; 
  }

  /// parametrization vol time only with hazard rate with time component
  void OnParametrizationVolTimeOnlyHRWithTimeComponent
       (const ParametrizationVolTimeOnlyHRWithTimeComponent& param)
  {
    m_pParametrizationVolTimeOnlyHRWithTimeComponent = 
      shared_ptr<ParametrizationVolTimeOnlyHRWithTimeComponent>
      (new ParametrizationVolTimeOnlyHRWithTimeComponent(param));
  }

  shared_ptr<ParametrizationVolTimeOnlyHRWithTimeComponent>
  GetParametrizationVolTimeOnlyHRWithTimeComponent() 
  { 
    return m_pParametrizationVolTimeOnlyHRWithTimeComponent; 
  }

private:

  shared_ptr<ParametrizationHRWithTimeComponent> 
    m_pParametrizationHRWithTimeComponent;

  shared_ptr<ParametrizationHRWithSpotComponentPower> 
    m_pParametrizationHRWithSpotComponentPower;

  shared_ptr<ParametrizationVolFlatHRWithTimeComponent> 
    m_pParametrizationVolFlatHRWithTimeComponent;

  shared_ptr<ParametrizationVolFlatHRWithSpotComponentPower> 
    m_pParametrizationVolFlatHRWithSpotComponentPower;

  shared_ptr<ParametrizationVolFlatHRPower> 
    m_pParametrizationVolFlatHRPower;

  shared_ptr<ParametrizationVolPower> 
    m_pParametrizationVolPower;

  shared_ptr<ParametrizationVolPowerHRPower> 
    m_pParametrizationVolPowerHRPower;

  shared_ptr<ParametrizationVolPowerHRWithTimeComponent> 
    m_pParametrizationVolPowerHRWithTimeComponent;

  shared_ptr<ParametrizationVolTanh> 
    m_pParametrizationVolTanh;

  shared_ptr<ParametrizationVolTanhHRWithTimeComponent> 
    m_pParametrizationVolTanhHRWithTimeComponent;

  shared_ptr<ParametrizationVolTanhHRWithSpotComponentPower> 
    m_pParametrizationVolTanhHRWithSpotComponentPower;

  shared_ptr<ParametrizationVolTanhHRPower> 
    m_pParametrizationVolTanhHRPower;

  shared_ptr<ParametrizationVolWithTimeComponent> 
    m_pParametrizationVolWithTimeComponent;

  shared_ptr<ParametrizationVolTimeOnlyHRWithTimeComponent> 
    m_pParametrizationVolTimeOnlyHRWithTimeComponent;
};

} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_XML_PARAMETRIZATION_VISITOR_GOODTYPE_H_
