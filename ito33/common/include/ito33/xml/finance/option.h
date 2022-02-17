/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/option.h
// Purpose:     Names of elements and attributes used in XML for an Option.
// Author:      Vadim Zeitlin
// Created:     2004-05-04
// RCS-ID:      $Id: option.h,v 1.12 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/option.h
    @brief Contains the names of the elements used in the XML description of a
           stock option.

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_OPTION_H_
#define _ITO33_XML_FINANCE_OPTION_H_

#include "ito33/sharedptr.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_OPTION_ROOT          "option"

#define XML_TAG_OPTION_IMPLIED_VOL   "implied_volatility"

//@}

namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Option;
}

namespace XML
{
  /// Restore option from xml node
  bool Restore(const xml::node& node, shared_ptr<finance::Option>& pOption);
}

} // namespace ito33


#endif // _ITO33_XML_FINANCE_OPTION_H_
