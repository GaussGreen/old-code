/////////////////////////////////////////////////////////////////////////////
// Name:        bond_xml.cpp
// Purpose:     Restore Bond object from XML document
// Author:      Vadim Zeitlin
// Created:     2004-05-08
// RCS-ID:      $Id: pepslike_xml.cpp,v 1.5 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/bondlike/pepslike.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/callfixedshare.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/bondlike/pepslike.h"
#include "ito33/xml/finance/bondlike/bondliketerms_all.h"
#include "ito33/xml/finance/bondlike/convertiblelike.h"
#include "ito33/xml/finance/bondlike/callschedule.h"
#include "ito33/xml/finance/bondlike/callfixedshare.h"
#include "ito33/xml/finance/bondlike/conversionschedule.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::finance;

ITO33_FORCE_LINK_THIS_MODULE(pepslike_xml);

/**
    Restore an PEPS-like object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new bond object to be deleted by the caller
 */
static
finance::Derivative *ReadPEPSLike(const xml::node *pNode)
{
  shared_ptr<BondLikeTerms> pBondLikeTerms = GetBondLikeTermsFromNode(*pNode);

  double dMinRatio = GetDoubleFromName
            (*pNode, XML_TAG_PEPSLIKE_MIN_CONVERSION_RATIO);
  double dMaxRatio = GetDoubleFromName
            (*pNode, XML_TAG_PEPSLIKE_MAX_CONVERSION_RATIO);

  PEPSLike *pPEPSLike = new PEPSLike(pBondLikeTerms, dMaxRatio, dMinRatio);

  GetOptionalConvertibleLikeDataFromNode(*pNode, *pPEPSLike);

  // xml::node::const_iterator *pNodeFound;

  // we can forget call_type,
  // pNodeFound = pNode->find(XML_TAG_PEPSLIKE_CALL_TYPE);

  // optional call provision
  shared_ptr<CallSchedule> pCall;
  
  if(Restore(*pNode, pCall))
    pPEPSLike->SetCallFixedCash(pCall);
  else
  {
    shared_ptr<CallFixedShare> pCallFixedShare;

    if(Restore(*pNode, pCallFixedShare))
      pPEPSLike->SetCallFixedShare(pCallFixedShare);
  }

  // conversion provision
  xml::node::const_iterator 
    pNodeFound = pNode->find
                 (XML_TAG_PEPSLIKE_HAS_OPTIONAL_CONVERSION);

  if ( pNodeFound != pNode->end() && GetBoolFromNode(*pNodeFound) )
    pPEPSLike->EnableOptionalConversion(); 

  return pPEPSLike;
}


ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_PEPSLIKE_ROOT, PEPSLike);
