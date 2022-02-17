/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/cds_xml.cpp
// Purpose:     Restore cds object from XML document
// Author:      Yann d'halluin
// Created:     2004-06-25
// RCS-ID:      $Id: cds_xml.cpp,v 1.18 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"
#include "ito33/date.h"
#include "ito33/useexception.h"

#include "ito33/finance/frequency.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/cds.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/cds.h"
#include "ito33/xml/finance/cashflowstream_all.h"
#include "ito33/xml/finance/common.h"

using namespace ito33;
using namespace ito33::XML;

ITO33_FORCE_LINK_THIS_MODULE(cds_xml);


/**
    Restore a cds object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
finance::Derivative *ReadCDS(const xml::node *pNode)
{
    //Recovery rate
  double dRecoveryRate = GetDoubleFromName(*pNode,XML_TAG_FINANCE_RECOVERYRATE);
  
  shared_ptr<finance::CashFlowStreamUniform>
    pSpreadStream = GetCashFlowStreamUniformInNode(
                      GetNodeByName(*pNode, XML_TAG_CDS_SPREADSTREAM));

  finance::Derivative *deriv =
      new ito33::finance::CDS(dRecoveryRate, pSpreadStream);
 
  GetOptionalDerivativeDataFromNode(*pNode, *deriv);

  GetMarketPrice(*pNode, *deriv); 

  return deriv;

}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_CDS_ROOT, CDS);


bool ito33::XML::Restore(const xml::node& node, shared_ptr<finance::CDS>& pCDS)
{
  if( strcmp(node.get_name(), XML_TAG_CDS_ROOT) != 0)
    return false;

  pCDS = shared_ptr<finance::CDS>
              ( dynamic_cast<finance::CDS*>(ReadCDS(&node)) );

  xml::node::const_iterator pNodeFound = node.find(XML_TAG_FINANCE_MARKETPRICE);
  if(pNodeFound != node.end())
    pCDS->SetMarketPrice( GetDoubleFromNode(*pNodeFound) );


  return true;
}

