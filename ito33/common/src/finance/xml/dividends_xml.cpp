/////////////////////////////////////////////////////////////////////////////
// Name:        dividends_xml.cpp
// Purpose:     Restore Derivatives object from XML document
// Author:      Vadim Zeitlin
// Created:     2004-05-04
// RCS-ID:      $Id: dividends_xml.cpp,v 1.8 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/dividends.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/dividends.h"
#include "ito33/xml/finance/common.h"

#include <stdio.h>

using namespace ito33;
using namespace ito33::XML;

shared_ptr<finance::Dividends>
ito33::XML::ReadDividends(const xml::node& node)
{
  shared_ptr<finance::Dividends> dividends(new finance::Dividends);

  for ( xml::node::const_iterator i = node.begin(); i != node.end(); ++i )
  {
    if ( strcmp(i->get_name(), XML_TAG_DIVIDEND) == 0 )
    {
      finance::Dividend::Type type;

      const std::string 
        strType = GetNodeByName(*i, XML_TAG_FINANCE_TYPE).get_content();

      if ( strcmp(strType.c_str(), XML_VALUE_DIVIDEND_TYPE_YIELD) == 0 )
      {
        type = finance::Dividend::Yield;
      }
      else if ( strcmp(strType.c_str(), XML_VALUE_DIVIDEND_TYPE_CASH) == 0 )
      {
        type = finance::Dividend::Cash;
      }
      else // skip dividend with unknown type
      {
        continue;
      }

      dividends->Add
                 (
                  type,
                  GetDateFromNode(GetNodeByName(*i, XML_TAG_FINANCE_DATE)),
                  GetDoubleFromNode(GetNodeByName(*i, XML_TAG_FINANCE_VALUE))
                 );
    }
    //else: ignore unknown tags
  }

  dividends->Validate();

  return dividends;
}

