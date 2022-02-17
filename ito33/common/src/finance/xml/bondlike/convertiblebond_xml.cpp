/////////////////////////////////////////////////////////////////////////////
// Name:        bond_xml.cpp
// Purpose:     Restore Bond object from XML document
// Author:      Vadim Zeitlin
// Created:     2004-05-08
// RCS-ID:      $Id: convertiblebond_xml.cpp,v 1.13 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/triggeraspercentageof.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/termstructure.h"
#include "ito33/xml/finance/bondlike/cb_base.h"
#include "ito33/xml/finance/bondlike/convertiblebond.h"
#include "ito33/xml/finance/bondlike/bondliketerms_all.h"
#include "ito33/xml/finance/bondlike/conversionschedule.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::finance;

ITO33_FORCE_LINK_THIS_MODULE(convertiblebond_xml);

void GetConvertibleBondConstructorData
        (
          const xml::node& node,
          shared_ptr<BondTerms>& pBondTerms,
          shared_ptr<ConversionSchedule>& pConversions
        )
{
  pBondTerms = GetBondTermsFromNode(node);

  if( !Restore(node, pConversions) )
  {
    typedef MissingNodeException Exception;
    throw EXCEPTION_MSG(node, "XML_TAG_BONDLIKE_CONVERSIONSCHEDULE_ROOT");
  }
}


/**
    Restore an Bond object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new bond object to be deleted by the caller
 */
static
finance::Derivative *ReadConvertibleBond(const xml::node *pNode)
{
  shared_ptr<BondTerms> pBondTerms;
  shared_ptr<ConversionSchedule> pConversions;
  GetConvertibleBondConstructorData(*pNode, pBondTerms, pConversions);

  ConvertibleBond *pBond = new ConvertibleBond(pBondTerms, pConversions);

  GetOptionalCBBaseDataFromNode(*pNode, *pBond);

  return pBond;
}
ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_CONVERTIBLEBOND_ROOT, ConvertibleBond);

// This funtion is used when reading a cb option
bool 
ito33::XML::Restore(const xml::node& node, 
                    shared_ptr<ConvertibleBond>& pcb)
{
  xml::node::const_iterator iter;
  if( (iter = node.find(XML_TAG_CONVERTIBLEBOND_ROOT) ) == node.end() )
    return false;

  pcb = shared_ptr<ConvertibleBond>( static_cast<ConvertibleBond*>
                                    (ReadConvertibleBond(&(*iter))) );

  return true;
}
