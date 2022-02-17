/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/xml/volatility.h
// Purpose:     Names of elements and attributes used in XML for Volatility
// Author:      Yann d'Halluin
// Created:     2004-05-19
// RCS-ID:      $Id: volatility.h,v 1.13 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/xml/volatility.h
    @brief Contains the names of the elements used in the XML description of a
           volatility. Only flat volatility is currently supported.
 */

#ifndef _IHG_XML_VOLATILITY_H_
#define _IHG_XML_VOLATILITY_H_

#include "ito33/sharedptr.h"
#include "ito33/factory.h"
#include "ito33/string.h"

#include "ito33/xml/finance/common.h"

#define XML_TAG_VOLATILITYCOMBO_ROOT            "volatility_combo"

#define XML_TAG_VOLATILITY_TIMECOMPONENT        "time_component"

#define XML_TAG_VOLATILITYTIMEONLY_ROOT  "volatility_time_only"
#define XML_TAG_VOLATILITYFLAT_ROOT      "volatility_flat"
#define XML_TAG_VOLATILITY_FLAT          "flat"

#define XML_TAG_VOLATILITYPOWER_ROOT     "volatility_power"
#define XML_TAG_VOLATILITYPOWER_ALPHA    "alpha"
#define XML_TAG_VOLATILITYPOWER_BETA     "beta"
#define XML_TAG_VOLATILITYPOWER_S0       "S0"

#define XML_TAG_VOLATILITYTANH_ROOT     "volatility_tanh"
#define XML_TAG_VOLATILITYTANH_LEFT     "left"
#define XML_TAG_VOLATILITYTANH_RIGHT    "right"
#define XML_TAG_VOLATILITYTANH_SCALE    "scale"
#define XML_TAG_VOLATILITYTANH_S0       "S0"

#include "ito33/factory.h"
#include "ito33/string.h"

namespace xml { class node; }

namespace ito33
{

  namespace ihg
  {
    class Volatility;
  }
/**
    VolatilityFactory allows creation of volatility by name.

    This factory allows creating objects of type ihg::Volatility using
    string and creator functions taking xml::node as parameter.
 */
typedef Factory<std::string, ihg::Volatility, xml::node> VolatilityFactory;

/**
    Macro to be used to define a factory for a specific volatility-derived
    class.

    This macro supposes that a function ReadVolatility() taking an xml::node
    and returning a new Volatility is defined.

    @param name of the XML tag from which volatility is read
    @param derivative the concrete type to be created for this key
 */
#define ITO33_DEFINE_VOLATILITY_READER(name, volatility)                      \
    /* use unnamed namespace to hide the declarations inside it */            \
    namespace                                                                 \
    {                                                                         \
        static VolatilityFactory                                              \
        FactoryFor ## volatility (name, Read ## volatility);                  \
    } struct Dummy /* just to force semicolon after the macro */


namespace ihg
{

namespace XML
{

/**
    Restore a volatility object from XML. Otherwise, a NULL pointer
    is returned.

    @param node the root tag in DOM tree where we should find volatility object
    @return the new volatility object
 */
 shared_ptr<ihg::Volatility> ReadVolatility(const xml::node &node);

} // namespace XML

} // namespace ihg

} // namespace ito33

#endif // _ITO33_XML_FINANCE_VOLATILITY_H_

