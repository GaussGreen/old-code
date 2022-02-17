/////////////////////////////////////////////////////////////////////////////
// Name:        derivative_xml.cpp
// Purpose:     contains all derivative reader definitions
// Author:      Vadim Zeitlin
// Created:     2004-05-12
// RCS-ID:      $Id: derivative_xml.cpp,v 1.12 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/derivative.h"
#include "ito33/finance/numeraire.h"

#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/issuer.h"
#include "ito33/xml/finance/derivative.h"

#include "ito33/xml/read.h"

#include "ito33/constants.h"


// define the field of the Factory declared in ito33/xml/finance/derivative.h
ITO33_IMPLEMENT_THE_FACTORY(ito33::DerivativeFactory);


namespace ito33
{

namespace XML
{

void GetMarketPrice(const xml::node& node, finance::Derivative& derivative)
{
  // market price is optional
  xml::node::const_iterator i;

  if ( (i = node.find(XML_TAG_FINANCE_MARKETPRICE)) != node.end() ) 
  {
    double 
      dMarketPrice = XML::GetDoubleFromName(node, XML_TAG_FINANCE_MARKETPRICE); 
    
    derivative.SetMarketPrice( dMarketPrice );
  }
         
} // GetMarketPrice


void GetOptionalDerivativeDataFromNode
     (const xml::node& node, finance::Derivative& derivative)
{  
  xml::node::const_iterator i;

  // Look for issuer
  if ( ( i = node.find(XML_TAG_ISSUER_ROOT) ) != node.end() )
    derivative.SetIssuer( GetIssuerFromNode(*i) );

  // look for currency
  if ( (i = node.find(XML_TAG_CURRENCY) ) != node.end() ) 
  {
    shared_ptr<finance::Numeraire> 
      pNumeraire( new finance::Numeraire(i->get_content()) );
   
    derivative.SetNumeraire( pNumeraire );
  }
         
} 

} // namespace XML

} // namespace ito33
