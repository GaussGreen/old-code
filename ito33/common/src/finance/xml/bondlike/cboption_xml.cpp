/////////////////////////////////////////////////////////////////////////////
// Name:        cboption_xml.cpp
// Purpose:     Restore cb option object from XML document
// Author:      Nabil
// Created:     2005/07/12
// RCS-ID:      $Id: cboption_xml.cpp,v 1.5 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/bondlike/cboption.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/read_frequency.h"
#include "ito33/xml/finance/read_daycountconvention.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/bondlike/cboption.h"
#include "ito33/xml/finance/bondlike/convertiblebond.h"
#include "ito33/xml/finance/floatingrates.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::finance;

ITO33_FORCE_LINK_THIS_MODULE(cboption_xml);

/**
    Restore a CB Option object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new cb option object to be deleted by the caller
 */
finance::Derivative* ReadCBOption(const xml::node *pNode)
{
  shared_ptr<ConvertibleBond> 
    pcb;
  
  Date maturityDate = 
    GetDateFromName(*pNode, XML_TAG_MATURITYDATE); 
  
  xml::node::const_iterator pNodeFound;

  Restore(*pNode, pcb);
 
  pNodeFound = pNode->find(XML_TAG_FLOATINGRATES_ROOT);

  shared_ptr<FloatingRates> 
  pfloatingRates ( GetFloatingRatesFromNode(*pNodeFound) );
  CBOption 
    *pCBOption = new CBOption(pcb, pfloatingRates, maturityDate);

  pNodeFound = pNode->find(XML_TAG_ASWNOTIONALTYPE);
  
  if ( pNodeFound != pNode->end() )
  {
    ASWNotionalIs
      aswNotionalIs = GetEnumFromNode
                            (
                              *pNodeFound,
                              SIZEOF(g_ASWNotionalIs),
                              g_ASWNotionalIs
                            );
    
    pCBOption->SetASWNotionalIs(aswNotionalIs);
  }

  GetOptionalDerivativeDataFromNode(*pNode, *pCBOption);

  GetMarketPrice(*pNode, *pCBOption);

  return pCBOption;
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_CBOPTION_ROOT, CBOption);
