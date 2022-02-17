/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/moneymarket.h
// Purpose:     Names of elements and attributes used in XML for a MoneyMarket
// Author:      Vadim Zeitlin
// Created:     2004-05-04
// RCS-ID:      $Id: moneymarket.h,v 1.11 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/moneymarket.h
    @brief Contains the names of the elements used in the XML description of a
           money market.

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_MONEYMARKET_H_
#define _ITO33_XML_FINANCE_MONEYMARKET_H_

#include "ito33/sharedptr.h"

/**
    @name Names of the elements and values in the XML description of the
          MoneyMarket class.
 */
//@{

#define XML_TAG_MONEYMARKET_ROOT         "money_market"

#define XML_TAG_MONEYMARKET_YIELDCURVE   "yield_curve"

#define XML_TAG_MONEYMARKET_NUMERAIRE    "currency"

// @}


namespace xml { class node; }

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL MoneyMarket;
}

namespace XML
{

/**
  Gets MoneyMarket from given node if possible.

  @param node MoneyMarket node in XML tree
  @return completed MoneyMarket object with information we read from node
  */
shared_ptr<finance::MoneyMarket> GetMoneyMarketFromNode(const xml::node& node);

} // namespace XML

} // namespace ito33



#endif // _ITO33_XML_FINANCE_MONEYMARKET_H_

