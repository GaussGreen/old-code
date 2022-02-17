/////////////////////////////////////////////////////////////////////////////
// Name:        bond_xml.cpp
// Purpose:     Restore Bond object from XML document
// Author:      Vadim Zeitlin
// Created:     2004-05-08
// RCS-ID:      $Id: bond_xml.cpp,v 1.6 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/bondlike/bond.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/bondlike/bond.h"
#include "ito33/xml/finance/bondlike/bondliketerms_all.h"

#include "ito33/xml/finance/bondlike/callschedule.h"
#include "ito33/xml/finance/bondlike/putschedule.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::finance;

ITO33_FORCE_LINK_THIS_MODULE(bond_xml);

void GetBondConstructorData(const xml::node& node, shared_ptr<BondTerms>& pBondTerms)
{
  pBondTerms = GetBondTermsFromNode(node);
}

void ito33::XML::GetOptionalBondDataFor(Bond& bond, const xml::node& node)
{
  GetOptionalDerivativeDataFromNode(node, bond);

  shared_ptr<finance::CallSchedule> pCalls;
  if ( Restore(node, pCalls) )
    bond.SetCallSchedule(pCalls);

  shared_ptr<finance::PutSchedule> pPuts;
  if ( Restore(node, pPuts) )
    bond.SetPutSchedule(pPuts);
}

/**
    Restore an Bond object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new bond object to be deleted by the caller
 */
static
finance::Derivative *ReadBond(const xml::node *pNode)
{
  shared_ptr<BondTerms> pBondTerms;
  GetBondConstructorData(*pNode, pBondTerms);

  Bond *pBond = new Bond(pBondTerms);

  GetOptionalBondDataFor(*pBond, *pNode);

  return pBond;
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_BOND_ROOT, Bond);
