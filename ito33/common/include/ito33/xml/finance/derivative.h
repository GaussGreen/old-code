/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/derivative.h
// Purpose:     DerivativeFactory class used for creating Derivatives from XML
// Author:      Vadim Zeitlin
// Created:     2004-05-04
// RCS-ID:      $Id: derivative.h,v 1.10 2006/07/28 21:01:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/derivative.h
    @brief 
 */

#ifndef _ITO33_XML_FINANCE_DERIVATIVE_H_
#define _ITO33_XML_FINANCE_DERIVATIVE_H_

#include "ito33/factory.h"
#include "ito33/string.h"

/// Market price
#define XML_TAG_FINANCE_MARKETPRICE "market_price"

namespace xml { class node; }

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivative;
}

/**
    DerivativeFactory allows creation of derivatives by name.

    This factory allows creating objects of type finance::Derivative using
    string and creator functions taking xml::node as parameter.

    Because the creator in Factory specifies that the third parameter is passed
    as pointer(so can be optional), this requires xml::node be passed as 
    pointer to the specific ReadDerivative function(weired). 
 */
typedef Factory<std::string, finance::Derivative, xml::node> DerivativeFactory;

/**
    Macro to be used to define a factory for a specific Derivative-derived
    class.

    This macro supposes that a function ReadDerivative() taking an xml::node
    and returning a new Derivative is defined.

    @param name of the XML tag from which derivative is read
    @param derivative the concrete type to be created for this key
 */
#define ITO33_DEFINE_DERIVATIVE_READER(name, derivative)                      \
    /* use unnamed namespace to hide the declarations inside it */            \
    namespace                                                                 \
    {                                                                         \
        static DerivativeFactory                                              \
        FactoryFor ## derivative (name, Read ## derivative);                  \
    } struct Dummy /* just to force semicolon after the macro */

namespace XML
{

/**
    Gets optional data of a derivative from xml node.

    @param node xml node
    @param derivative object whose optional data will be filled.
 */
void GetOptionalDerivativeDataFromNode
     (const xml::node& node, finance::Derivative& derivative);

/**
    Reads the market price, common to all derivative instrument.

    @param node of the xml tag to read
    @param derivative derivative contract 

 */
void GetMarketPrice(const xml::node& node, finance::Derivative& derivative);

} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_FINANCE_DERIVATIVE_H_
