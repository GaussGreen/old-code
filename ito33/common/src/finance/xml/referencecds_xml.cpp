/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/referencecds_xml.cpp
// Purpose:     Restore reference cds object from XML document
// Created:     2006/05/17
// RCS-ID:      $Id: referencecds_xml.cpp,v 1.6 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"
#include "ito33/date.h"

#include "ito33/finance/frequency.h"
#include "ito33/finance/referencecds.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/referencecds.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/read_daycountconvention.h"
#include "ito33/xml/finance/read_frequency.h"

using namespace ito33;
using namespace ito33::XML;

ITO33_FORCE_LINK_THIS_MODULE(referencecds_xml);


/**
    Restore a reference cds object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new reference cds object to be deleted by the caller
 */
static
finance::Derivative *ReadReferenceCDS(const xml::node *pNode)
{
  size_t nMonthsToMaturity = 
    GetLongFromName(*pNode, XML_TAG_REFERENCECDS_MATURITY_MONTHS);  

  Date::DayCountConvention 
    dcc = GetDayCountConventionFromName(*pNode, XML_TAG_DAYCOUNTCONVENTION);

  finance::Frequency 
    freq = GetFrequencyFromName(*pNode, XML_TAG_PAYMENTFREQUENCY);

  double 
    dRecoveryRate = GetDoubleFromName(*pNode, XML_TAG_FINANCE_RECOVERYRATE);
  
  finance::ReferenceCDS* pRefCDS =
      new ito33::finance::ReferenceCDS(nMonthsToMaturity,
                                       freq, dcc, dRecoveryRate);
 
  // Get the (mostly) optional parameters
  GetOptionalDerivativeDataFromNode(*pNode, *pRefCDS);

  GetMarketPrice(*pNode, *pRefCDS); 

  xml::node::const_iterator pNodeSpread;
  pNodeSpread = pNode->find(XML_TAG_REFERENCECDS_SPREAD);
  if( pNodeSpread != pNode->end() )
    pRefCDS->SetSpread( GetDoubleFromNode(*pNodeSpread) );
  
  return pRefCDS;
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_REFERENCECDS_ROOT, ReferenceCDS);


bool ito33::XML::Restore(const xml::node& node, 
                         shared_ptr<finance::ReferenceCDS>& pRefCDS)
{
  if( strcmp(node.get_name(), XML_TAG_REFERENCECDS_ROOT) != 0)
    return false;

  pRefCDS = shared_ptr<finance::ReferenceCDS>
            ( dynamic_cast<finance::ReferenceCDS*>(ReadReferenceCDS(&node)) );

  xml::node::const_iterator pNodeFound = node.find(XML_TAG_FINANCE_MARKETPRICE);
  if (pNodeFound != node.end())
    pRefCDS->SetMarketPrice( GetDoubleFromNode(*pNodeFound) );

  return true;
}

