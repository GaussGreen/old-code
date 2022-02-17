/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/equity.h
// Purpose:     Names of elements and attributes used in XML for an Equity
// Created:     2006/03/16
// RCS-ID:      $Id: equity.h,v 1.4 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/equity.h
    @brief Contains the names of the elements used in the XML description of
           the Equity class.
 */

#ifndef _ITO33_XML_FINANCE_EQUITY_H_
#define _ITO33_XML_FINANCE_EQUITY_H_

#include "ito33/sharedptr.h"

/**
    @name Names of the elements and values in the XML description of the
          Equity class.
 */
//@{

/// The name of the root tag containing equity description
#define XML_TAG_EQUITY_ROOT             "equity"

#define XML_TAG_EQUITY_SPOTSHAREPRICE   "spot_share_price"

#define XML_TAG_EQUITY_BORROWCURVE      "borrow_curve"

#define XML_TAG_EQUITY_PREVIOUS_SHARE_PRICE    "previous_share_price"

//@}

namespace xml { class node; }

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Equity;
}

namespace XML
{

/**
    Get an Equity object from XML. Exception is thrown if format is wrong.

    @param node the Equity node in DOM tree
    @return the equity object constructed from information read from the node
 */
shared_ptr<finance::Equity> GetEquityFromNode(const xml::node& node);

} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_FINANCE_EQUITY_H_

