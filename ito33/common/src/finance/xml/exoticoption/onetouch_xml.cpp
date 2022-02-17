/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/exoticoption/onetouch_xml.cpp
// Purpose:     Restore OneTouch object from XML document
// Created:     2005/07/04
// RCS-ID:      $Id: onetouch_xml.cpp,v 1.7 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/exoticoption/onetouch.h"
#include "ito33/finance/exoticoption/fxonetouch.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/barrier.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/exoticoption/onetouch.h"

#include "ito33/xml/finance/read_barriertype.h"
#include "ito33/xml/finance/read_rebatetype.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(onetouch_xml);

using namespace ito33;
using namespace ito33::XML;
using finance::OneTouch;
using finance::FXOneTouch;

static
finance::Derivative* ReadOneTouch(const xml::node* pNode)
{
  const xml::node& node = *pNode;

  // maturity date
  Date maturityDate = GetDateFromName(node, XML_TAG_FINANCE_MATURITY);

  // Barrier
  double dBarrier = GetDoubleFromName(node, XML_TAG_BARRIER);

  finance::BarrierType 
    barrierType = GetBarrierTypeFromName(node, XML_TAG_BARRIER_TYPE);

  finance::RebateType 
    rebateType = GetRebateTypeFromName(node, XML_TAG_REBATE_TYPE);

  finance::OneTouch*
    pOneTouch = new finance::OneTouch
                    (maturityDate, dBarrier, barrierType, rebateType);

  GetOptionalDerivativeDataFromNode(node, *pOneTouch);

  GetMarketPrice(node, *pOneTouch);

  return pOneTouch;
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_ONETOUCH_ROOT, OneTouch);

// This funtion is used when reading termstructures
bool ito33::XML::Restore(const xml::node& node, shared_ptr<OneTouch>& pOneTouch)
{
  if ( strcmp(node.get_name(), XML_TAG_ONETOUCH_ROOT) != 0 )
    return false;

#ifndef NDEBUG
  pOneTouch = 
    shared_ptr<OneTouch>( dynamic_cast<OneTouch*>(ReadOneTouch(&node)) );
#else
  pOneTouch = 
    shared_ptr<OneTouch>( static_cast<OneTouch*>(ReadOneTouch(&node)) );
#endif

  xml::node::const_iterator 
    pNodeFound = node.find(XML_TAG_FINANCE_MARKETPRICE);

  if ( pNodeFound != node.end() )
    pOneTouch->SetMarketPrice( GetDoubleFromNode(*pNodeFound) );

  return true;
}


static
finance::Derivative* ReadFXOneTouch(const xml::node* pNode)
{
  const xml::node& node = *pNode;

  // maturity date
  Date maturityDate = GetDateFromName(node, XML_TAG_FINANCE_MATURITY);

  // Black-Scholes barrier
  double dBSBarrier = GetDoubleFromName(node, XML_TAG_FXONETOUCH_BSBARRIER);

  // barrier type
  finance::BarrierType 
    barrierType = GetBarrierTypeFromName(node, XML_TAG_BARRIER_TYPE);

  // reference volatility
  double dRefVol = GetDoubleFromName(node, XML_TAG_FXONETOUCH_REFVOL);

  finance::FXOneTouch*
    pFXOneTouch = new finance::FXOneTouch(maturityDate, dBSBarrier, 
                                          barrierType, dRefVol);

  // read and set the market quote, if defined
  xml::node::const_iterator 
    pNodeFound = pNode->find(XML_TAG_FXONETOUCH_QUOTE);

  GetOptionalDerivativeDataFromNode(node, *pFXOneTouch);

  if ( pNodeFound != pNode->end() )
    pFXOneTouch->SetMarketQuote( GetDoubleFromNode(*pNodeFound) );
  
  return pFXOneTouch;
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_FXONETOUCH_ROOT, FXOneTouch);

// This funtion is used when reading termstructures
bool ito33::XML::Restore(const xml::node& node, 
                         shared_ptr<FXOneTouch>& pFXOneTouch)
{
  if ( strcmp(node.get_name(), XML_TAG_FXONETOUCH_ROOT) != 0 )
    return false;

#ifndef NDEBUG
  pFXOneTouch = 
    shared_ptr<FXOneTouch>( dynamic_cast<FXOneTouch*>(ReadFXOneTouch(&node)) );
#else
  pFXOneTouch = 
    shared_ptr<FXOneTouch>( static_cast<FXOneTouch*>(ReadFXOneTouch(&node)) );
#endif

  return true;
}
