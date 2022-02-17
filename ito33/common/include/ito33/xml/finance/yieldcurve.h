/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/yieldcurve.h
// Purpose:     Names of elements and attributes used in XML for a YieldCurve
// Author:      Vadim Zeitlin
// Created:     2004-05-04
// RCS-ID:      $Id: yieldcurve.h,v 1.16 2006/08/23 09:17:57 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/yieldcurve.h
    @brief Contains the names of the elements used in the XML description of a
           yield curve.

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_YIELDCURVE_H_
#define _ITO33_XML_FINANCE_YIELDCURVE_H_

#include "ito33/sharedptr.h"
#include "ito33/factory.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/string.h"

namespace xml 
{ 
  class node; 
}
//=============================================================================
//=              factory related                                              =
//=============================================================================
namespace ito33
{

/**
    YieldCurveFactory allows creation of yieldcurve by name.

    This factory allows creating objects of type finance::YieldCurve using
    string and creator functions taking xml::node as parameter.
 */
  typedef Factory<std::string, finance::YieldCurve, xml::node> YieldCurveFactory;

/**
    Macro to be used to define a factory for a specific YieldCurve-derived
    class.

    This macro supposes that a function ReadYieldCurve() taking an xml::node
    and returning a new YieldCurve is defined.

    @param name of the XML tag from which yieldcurve is read
    @param yieldcurve the concrete type to be created for this key
 */
#define ITO33_DEFINE_YIELDCURVE_READER(name, yieldcurve)                      \
    /* use unnamed namespace to hide the declarations inside it */            \
    namespace                                                                 \
    {                                                                         \
        static YieldCurveFactory                                              \
        FactoryFor ## yieldcurve (name, Read ## yieldcurve);                  \
    } struct Dummy /* just to force semicolon after the macro */


} // namespace ito33


//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{
namespace XML
{

/**
  Gets YieldCurve from given node

  @param node given node
  @return shared pointer to YieldCurve
  */
shared_ptr<finance::YieldCurve> GetYieldCurveFromNode(const xml::node& node);

} // namespace XML
} // namespace ito33


//=============================================================================
//=              tag name macros                                              =
//=============================================================================

/**
   @name Tag name macros
*/
//@{

// tags related to yieldcurve legs
#define XML_TAG_SWAPCURVELEG_RATE              "rate"
#define XML_TAG_SWAPCURVELEG_MATURITY_DURATION "maturity_duration_part"
#define XML_TAG_SWAPCURVELEG_MATURITY_UNIT     "maturity_unit_part"

#define XML_TAG_CASHRATE_ROOT "cash_rate"

#define XML_TAG_SWAPRATE_ROOT "swap_rate"

#define XML_TAG_ZEROCOUPONRATE_ROOT "zero_coupon_rate"

// unlike other elements, there is no root yield curve tag because we typically
// have several of them and each has its own name
#define XML_TAG_YIELDCURVEFLAT_ROOT               "yield_curve_flat"
#define XML_TAG_YIELDCURVEANNUALLYCOMPOUNDED_ROOT "yield_curve_annually_compounded"
#define XML_TAG_YIELDCURVESWAP_ROOT               "swap_curve"
#define XML_TAG_YIELDCURVEZEROCOUPON_ROOT         "zero_coupon_rate_curve"

// common tag
#define XML_TAG_YIELDCURVE_REFERENCE_DATE "reference_date"

// for annaully compounded yield curve
#define XML_TAG_YIELDCURVE_LEGS           "legs"
#define XML_TAG_YIELDCURVE_LEG            "leg"
#define XML_TAG_YIELDCURVE_LEG_DAY        "day"

// for swap curve
#define XML_TAG_YIELDCURVE_SWAP_CASHBASIS   "basis_for_cash_rate"
#define XML_TAG_YIELDCURVE_SWAP_SWAPBASIS   "basis_for_swap_rate"

//@}

#endif // _ITO33_XML_FINANCE_YIELDCURVE_H_

