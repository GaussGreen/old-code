/////////////////////////////////////////////////////////////////////////////
// Name:        percslike_xml.cpp
// Purpose:     Restore PERCSLike object from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-12-06
// RCS-ID:      $Id: percslike_xml.cpp,v 1.3 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/bondlike/percslike.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/callfixedshare.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/bondlike/percslike.h"
#include "ito33/xml/finance/bondlike/bondliketerms_all.h"
#include "ito33/xml/finance/bondlike/convertiblelike.h"
#include "ito33/xml/finance/bondlike/callschedule.h"
#include "ito33/xml/finance/bondlike/callfixedshare.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::finance;

ITO33_FORCE_LINK_THIS_MODULE(percslike_xml);

/**
    Restore an percs-like object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new bond object to be deleted by the caller
 */
static
finance::Derivative *ReadPERCSLike(const xml::node *pNode)
{
  shared_ptr<BondLikeTerms> pBondLikeTerms = GetBondLikeTermsFromNode(*pNode);

  double dRatio = GetDoubleFromName
            (*pNode, XML_TAG_PERCSLIKE_CONVERSION_RATIO_AT_MATURITY);
  double dCapPrice = GetDoubleFromName
            (*pNode, XML_TAG_PERCSLIKE_CAP_PRICE);

  PERCSLike *pPERCSLike = new PERCSLike(pBondLikeTerms, dCapPrice, dRatio);

  GetOptionalConvertibleLikeDataFromNode(*pNode, *pPERCSLike);

  // xml::node::const_iterator *pNodeFound;

  // we can forget call_type,
  // pNodeFound = pNode->find(XML_TAG_PERCSLIKE_CALL_TYPE);

  // optional call provision
  shared_ptr<CallSchedule> pCall;
  
  if(Restore(*pNode, pCall))
    pPERCSLike->SetCallSchedule(pCall);

  return pPERCSLike;
}


ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_PERCSLIKE_ROOT, PERCSLike);
