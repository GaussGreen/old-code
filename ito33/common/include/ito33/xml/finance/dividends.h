/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/dividends.h
// Purpose:     Names of elements and attributes used in XML for Dividends
// Author:      Vadim Zeitlin
// Created:     2004-05-04
// RCS-ID:      $Id: dividends.h,v 1.12 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/dividends.h
    @brief Contains the names of the elements used in the XML description of a
           (discrete) dividends.

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_DIVIDENDS_H_
#define _ITO33_XML_FINANCE_DIVIDENDS_H_

#include "ito33/sharedptr.h"
/**
   @name Tag name macros
*/
//@{

#define XML_TAG_DIVIDEND              "dividend"

#define XML_VALUE_DIVIDEND_TYPE_YIELD "yield"

#define XML_VALUE_DIVIDEND_TYPE_CASH  "cash"

//@}

namespace xml 
{ 
  /// Forward declaration
  class node; 
}

namespace ito33
{

namespace finance
{
  /// Forward declaration
  class ITO33_DLLDECL Dividends;
}

namespace XML
{

/**
    Restore a Dividends object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new dividends object
 */
shared_ptr<finance::Dividends> ReadDividends(const xml::node& node);

} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_FINANCE_DIVIDENDS_H_

