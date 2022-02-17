/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/xml/parametrization_visitor.h
// Purpose:     Visitor for parametrization classes
// Author:      ITO33
// Created:     2004/12/02
// RCS-ID:      $Id: parametrization_visitor.h,v 1.8 2006/01/10 17:25:07 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/xml/parametrization_visitor.h
    @brief Visitor class accepting different parametrizations.

    This is the central part of the implementation of the visitor pattern for
    the parametrizations. To do something different depending on the exact 
    type of parametrization object you have to define a new class inheriting 
    from this one and do whatever is required in its methods.
 */

#ifndef _ITO33_IHG_XML_PARAMETRIZATION_VISITOR_H_
#define _ITO33_IHG_XML_PARAMETRIZATION_VISITOR_H_

#include "ito33/debug.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace ihg
{
  class ITO33_DLLDECL ParametrizationHRWithTimeComponent;
  class ITO33_DLLDECL ParametrizationHRWithSpotComponentPower;

  class ITO33_DLLDECL ParametrizationVolFlatHRWithTimeComponent;
  class ITO33_DLLDECL ParametrizationVolFlatHRWithSpotComponentPower;
  class ITO33_DLLDECL ParametrizationVolFlatHRPower;
  class ITO33_DLLDECL ParametrizationVolPower;
  class ITO33_DLLDECL ParametrizationVolPowerHRSpotComponentPower;
  class ITO33_DLLDECL ParametrizationVolPowerHRWithTimeComponent;
  class ITO33_DLLDECL ParametrizationVolTanh;
  class ITO33_DLLDECL ParametrizationVolTanhHRWithSpotComponentPower;
  class ITO33_DLLDECL ParametrizationVolTanhHRWithTimeComponent;
  class ITO33_DLLDECL ParametrizationVolTanhHRPower;
  class ITO33_DLLDECL ParametrizationVolTimeOnlyHRWithTimeComponent;

  class ParametrizationVolWithTimeComponent;


/**
    Parametrization visitor.

    For now, only have one parametrization visitor class. If specific
    visitors for single (or particular) parametrizations are needed,
    then this class can be split as was done for derivative_visitor,
    derivative_visitor_goodtype, and XXXVisitor.
    
    Using a visitor might be less usual than querying the data directly for
    but it allows us to not lose the type information about the hazard rate
    and keep the code maintainable and extensible.
 */
class ParametrizationVisitor
{
public:

  ParametrizationVisitor() {}

  /// Virtual dtor for any base class
  virtual ~ParametrizationVisitor() { }

  /// parametrization with hazard rate with time component
  virtual void OnParametrizationHRWithTimeComponent
               (const ParametrizationHRWithTimeComponent&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization with hazard rate with spot component power
  virtual void OnParametrizationHRWithSpotComponentPower
               (const ParametrizationHRWithSpotComponentPower&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization vol flat with hazard rate with time component
  virtual void OnParametrizationVolFlatHRWithTimeComponent
               (const ParametrizationVolFlatHRWithTimeComponent&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization vol flat with hazard rate with spot component power
  virtual void OnParametrizationVolFlatHRWithSpotComponentPower
               (const ParametrizationVolFlatHRWithSpotComponentPower&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization vol flat with hazard rate power
  virtual void OnParametrizationVolFlatHRPower
               (const ParametrizationVolFlatHRPower&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization vol power
  virtual void OnParametrizationVolPower(const ParametrizationVolPower&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization vol power hr power
  virtual void 
  OnParametrizationVolPowerHRPower(const ParametrizationVolPowerHRPower&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization vol power with hazard rate with time component
  virtual void OnParametrizationVolPowerHRWithTimeComponent
               (const ParametrizationVolPowerHRWithTimeComponent&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization vol tanh
  virtual void OnParametrizationVolTanh(const ParametrizationVolTanh&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization vol tanh hr with spot component power
  virtual void 
  OnParametrizationVolTanhHRWithSpotComponentPower
  (const ParametrizationVolTanhHRWithSpotComponentPower&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization vol tanh hazard rate with time component
  virtual void OnParametrizationVolTanhHRWithTimeComponent
               (const ParametrizationVolTanhHRWithTimeComponent&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization vol tanh hazard rate power
  virtual void OnParametrizationVolTanhHRPower
               (const ParametrizationVolTanhHRPower&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization with volatility with time component
  virtual void OnParametrizationVolWithTimeComponent
               (const ParametrizationVolWithTimeComponent&)
  { 
    FAIL("Not implemented."); 
  }

  /// parametrization with time only volatility and hr with time component
  virtual void OnParametrizationVolTimeOnlyHRWithTimeComponent
               (const ParametrizationVolTimeOnlyHRWithTimeComponent&)
  { 
    FAIL("Not implemented."); 
  }
};


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_XML_PARAMETRIZATION_VISITOR_H_
