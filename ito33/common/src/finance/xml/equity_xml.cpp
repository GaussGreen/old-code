/////////////////////////////////////////////////////////////////////////////
// Name:        equitydata_xml.cpp
// Purpose:     Restore Equity object from XML document
// Created:     2006/03/23
// RCS-ID:      $Id: equity_xml.cpp,v 1.4 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/date.h"

#include "ito33/finance/equity.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/numeraire.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/equity.h"
#include "ito33/xml/finance/issuer.h"
#include "ito33/xml/finance/yieldcurve.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/dividends.h"

using namespace ito33;

namespace ito33
{
namespace XML
{


shared_ptr<finance::Equity> GetEquityFromNode(const xml::node& node)
{

  xml::node::const_iterator pNodeFound;

  // Numeraire
  shared_ptr<finance::Numeraire> pNumeraire(new finance::Numeraire(
    GetNodeByName(node, XML_TAG_CURRENCY).get_content() ) );

  // Spot
  double dSpot = GetDoubleFromName(node, XML_TAG_EQUITY_SPOTSHAREPRICE);

  // Construct and then look for optional data
  shared_ptr<finance::Equity> pEquity(new finance::Equity(dSpot, pNumeraire));

  // Issuer is optional
  pNodeFound = node.find(XML_TAG_ISSUER_ROOT);

  if ( pNodeFound != node.end() )
  {
    shared_ptr<finance::Issuer> pIssuer = GetIssuerFromNode(*pNodeFound);
    pEquity->SetIssuer(pIssuer);
  }

  // Borrow curve is optional
  pNodeFound = node.find(XML_TAG_EQUITY_BORROWCURVE);

  if ( pNodeFound != node.end() )
  {
    shared_ptr<finance::YieldCurve> pYC = GetYieldCurveFromNode(*pNodeFound);
    pEquity->SetBorrowCurve(pYC);
  }
  
  // Dividends are optional
  pNodeFound = node.find(XML_TAG_FINANCE_DIVIDENDS);

  if ( pNodeFound != node.end() )
  {
    shared_ptr<finance::Dividends> pDividends = ReadDividends(*pNodeFound);
    pEquity->SetDividends(pDividends);
  }

  //previous spot is optional
  pNodeFound = node.find(XML_TAG_EQUITY_PREVIOUS_SHARE_PRICE);
  
  if ( pNodeFound != node.end() )
  { 
    double dPreviousSharePrice = GetDoubleFromNode(*pNodeFound);
    pEquity->SetPreviousSharePrice(dPreviousSharePrice);
  }


  return pEquity;

}

} // namespace XML

} // namespace ito33

