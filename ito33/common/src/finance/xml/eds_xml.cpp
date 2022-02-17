/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/eds_xml.cpp
// Purpose:     Restore EDS object from XML document
// Created:     2005/01/26
// RCS-ID:      $Id: eds_xml.cpp,v 1.12 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/eds.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/barrier.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/cashflowstream_all.h"
#include "ito33/xml/finance/eds.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(eds_xml);

using namespace ito33;
using namespace ito33::XML;

/**
    Restore a EDS object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session data tag in DOM tree
   
    @return the new option object to be deleted by the caller
 */
static
finance::Derivative *ReadEDS(const xml::node *pNode)
{
  // Recovery rate
  double dRecoveryRate = GetDoubleFromName(*pNode,XML_TAG_FINANCE_RECOVERYRATE);
  
  shared_ptr<finance::CashFlowStreamUniform>
    pSpreadStream = GetCashFlowStreamUniformInNode(
                      GetNodeByName(*pNode, XML_TAG_EDS_SPREADSTREAM));

  // Barrier
  double dBarrier = GetDoubleFromName(*pNode, XML_TAG_BARRIER);

  // todo: market price? Will we use EDS for calibration? 
  finance::EDS*
    pEDS = new ito33::finance::EDS(dRecoveryRate, pSpreadStream, dBarrier);

  GetOptionalDerivativeDataFromNode(*pNode, *pEDS);

  GetMarketPrice(*pNode, *pEDS);

  return pEDS;
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_EDS_ROOT, EDS);

// This funtion is used when reading an EDS termstructure
bool ito33::XML::Restore(const xml::node& node, shared_ptr<finance::EDS>& pEDS)
{
  if( strcmp(node.get_name(), XML_TAG_EDS_ROOT) != 0)
    return false;

  pEDS = shared_ptr<finance::EDS>
              ( dynamic_cast<finance::EDS*>(ReadEDS(&node)) );

  xml::node::const_iterator pNodeFound = node.find(XML_TAG_FINANCE_MARKETPRICE);
  if(pNodeFound != node.end())
    pEDS->SetMarketPrice( GetDoubleFromNode(*pNodeFound) );

  return true;
}

