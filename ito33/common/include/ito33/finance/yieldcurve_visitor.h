/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/yieldcurve_visitor.h
// Purpose:     Visitor for YieldCurve-derived classes
// Author:      Yann d'Halluin
// Created:     2004-06-14
// RCS-ID:      $Id: yieldcurve_visitor.h,v 1.7 2006/08/21 14:27:26 zhang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/yieldcurve_visitor.h
    @brief Visitor class accepting different yield curve kinds.

    This is the central part of the implementation of the visitor pattern for
    the yieldcurves. To do something different depending on the exact type of
    a yield curve object you have to define a new class inheriting from this one
    and do whatever is required in its methods.
 */

#ifndef _ITO33_FINANCE_YIELDCURVE_VISITOR_H_
#define _ITO33_FINANCE_YIELDCURVE_VISITOR_H_

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL YieldCurveFlat;
class ITO33_DLLDECL YieldCurveAnnuallyCompounded;
class ITO33_DLLDECL YieldCurveSwap;
class ITO33_DLLDECL YieldCurveZeroCoupon;


/**
    YieldCurve visitor.

    An object of a class derived from this one must be passed by the user code
    to any functions which may access heterogeneous collections of Yield curves,
    e.g. XML::Reader in IHG currently.

    Using a visitor might be less usual than querying the data directly for
    but it allows us to not lose the type information about the derivatives
    and keep the code maintainable and extensible.
 */
class ITO33_DLLDECL YieldCurveVisitor
{
public:

  /// Called for a flat yield curve
  virtual void OnYieldCurveFlat(const YieldCurveFlat& ycf) = 0;

  /// Called for an annually compounded yield curve
  virtual void OnYieldCurveAnnuallyCompounded(
                const YieldCurveAnnuallyCompounded& ycac) = 0;

  /// Called for a swap curve
  virtual void OnYieldCurveSwap(const YieldCurveSwap& ycs) = 0;

  /// Called for a zero coupon rate curve
  virtual void OnYieldCurveZeroCoupon(const YieldCurveZeroCoupon& ycs) = 0;

  /// Virtual dtor for any base class
  virtual ~YieldCurveVisitor() { }
};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_YIELDCURVE_VISITOR_H_
