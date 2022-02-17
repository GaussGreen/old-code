/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardrate_visitor.h
// Purpose:     Visitor for hazardrate-derived classes
// Author:      Yann d'halluin
// Created:     2004-18-07
// RCS-ID:      $Id: hazardrate_visitor.h,v 1.2 2006/01/10 17:25:07 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/hazardrate_visitor.h
    @brief Visitor class accepting different Hazard rate kinds.

    This is the central part of the implementation of the visitor pattern for
    the hazardrate. To do something different depending on the exact type of
    hazard rate object you have to define a new class inheriting from this one
    and do whatever is required in its methods.

    @TODO this file is not updated as it is not used till now.
 */

#ifndef _ITO33_IHG_HAZARDRATE_VISITOR_H_
#define _ITO33_IHG_HAZARDRATE_VISITOR_H_

#include "ito33/debug.h"

#include "ito33/ihg/common.h"

namespace ito33
{

namespace ihg
{

  class ITO33_IHG_DLLDECL HazardRateFlat;
  class ITO33_IHG_DLLDECL HazardRateCombo;
  class ITO33_IHG_DLLDECL HazardRateTimeOnly;
  class ITO33_IHG_DLLDECL HazardRatePower;
  class ITO33_IHG_DLLDECL HazardRateLinear;
  class ITO33_IHG_DLLDECL HazardRateDecay;

/**
   Hazard rate visitor.

    An object of a class derived from this one must be passed by the user code
    to any functions which may access heterogeneous collections of Hazard Rate,
    e.g. XML::Reader in IHG currently.

    Using a visitor might be less usual than querying the data directly for
    but it allows us to not lose the type information about the hazard rate
    and keep the code maintainable and extensible.
 */
class HazardRateVisitor
{
public:

  /*
    We have three kinds of hazard rate classes

    1. exported in public interface.
    2. for internal use.
    3. deprecated classes and other classes that even visitor has trouble to
       get out any information (e.g. CallBack).

    They are treated in different ways.
    1. The visit functions for the first kind of classes are declared purely 
       virtual in base visitor. 
    2. The visit functions for the second kind of classes are implemented in
       base visitor by asserting. So a special visitor needs to implement the
       visit function only when it needs to check the concerned class, of
       course, for internal use.
    3. Visit functions for these classes are not supported.
   */

  /// @name visit functions for public hazard rate classes
  //@{

  /// Visits hazard rate combo
  virtual void OnHazardRateCombo(const HazardRateCombo &hrCombo) = 0;
  
  /// Visits Time only hazard rate
  virtual void OnHazardRateTimeOnly(const HazardRateTimeOnly &hrTimeOnly) = 0;

  /// Visits power hazard rate
  virtual void OnHazardRatePower(const HazardRatePower &hrPower) = 0;

  //@}

  /// @name visit functions for internal hazard rate classes
  //@{

  /// Flat hazard rate
  virtual void OnHazardRateFlat(const HazardRateFlat&)
  {
    FAIL ( "OnHazardRateFlat() not implemented in the visitor." );
  }

  /// decay hazard rate
  virtual void OnHazardRateDecay(const HazardRateDecay&)
  {
    FAIL ( "OnHazardRateDecay() not implemented in the visitor." );
  }

  /// linear hazard rate
  virtual void OnHazardRateLinear(const HazardRateLinear&)
  {
    FAIL ( "OnHazardRateLinear() not implemented in the visitor." );
  }

  //@}

  /// Virtual dtor for any base class
  virtual ~HazardRateVisitor() { }
};

} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HAZARDRATE_VISITOR_H_
