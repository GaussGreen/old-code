/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/logcontract_xml.cpp
// Purpose:     Restore Option object from XML document
// Created:     2006/07/18
// RCS-ID:      $Id: logcontract_xml.cpp,v 1.3 2006/07/28 21:01:10 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/logcontract.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/logcontract.h"

using namespace ito33;
using namespace ito33::XML;
using finance::LogContract;

ITO33_FORCE_LINK_THIS_MODULE(logcontract_xml);

/**
    Restore a log contract object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new log contract object to be deleted by the caller
 */
static
finance::Derivative *ReadLogContract(const xml::node *pNode)
{
  const xml::node& node = *pNode;

  Date maturityDate = GetDateFromName(node, XML_TAG_FINANCE_MATURITY);
  Date startOfReturnPeriod = GetDateFromName(node, XML_TAG_LOGCONTRACT_T0);
  
  finance::LogContract *
    pLogContract = new LogContract(maturityDate, startOfReturnPeriod);

  xml::node::const_iterator 
    pNodeFound = node.find(XML_TAG_LOGCONTRACT_S0);

  if ( pNodeFound != node.end() )
    pLogContract->SetStartSharePrice( GetDoubleFromNode(*pNodeFound) );

  GetOptionalDerivativeDataFromNode(*pNode, *pLogContract);

  GetMarketPrice(node, *pLogContract);

  return pLogContract; 
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_LOGCONTRACT_ROOT, LogContract);
