/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/bondlike/attachedwarrantcb_xml.cpp
// Purpose:     Restore attached warrant bond object from XML document
// Author:      ITO33 
// Created:     2005-03-04
// RCS-ID:      $Id: attachedwarrantcb_xml.cpp,v 1.4 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"
#include "ito33/finance/bondlike/sharedependentconversion.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/bondlike/cb_base.h"
#include "ito33/xml/finance/bondlike/bondliketerms_all.h"
#include "ito33/xml/finance/bondlike/sharedependentconversion.h"
#include "ito33/xml/finance/bondlike/attachedwarrantconvertiblebond.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::finance;

ITO33_FORCE_LINK_THIS_MODULE(attachedwarrantcb_xml);

void GetConvertibleBondConstructorData
        (
          const xml::node& node,
          shared_ptr<BondTerms>& pBondTerms,
          shared_ptr<ShareDependentConversion>& pShareDeConv
        )
{
  pBondTerms = GetBondTermsFromNode(node);

  if( !Restore(node, pShareDeConv) )
  {
    typedef MissingNodeException Exception;
    throw EXCEPTION_MSG(node, "XML_TAG_BONDLIKE_SHAREDEPENDENT_ROOT");
  }
}


/**
    Restore an attached warrant convertible bond object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new bond object to be deleted by the caller
 */
static
finance::Derivative *ReadAttachedWarrantCB(const xml::node *pNode)
{
  shared_ptr<BondTerms> pBondTerms;
  shared_ptr<ShareDependentConversion> pShareDeConv;
  GetConvertibleBondConstructorData(*pNode, pBondTerms, pShareDeConv);

  AttachedWarrantConvertibleBond *pWarrant 
    = new AttachedWarrantConvertibleBond(pBondTerms, pShareDeConv);

  GetOptionalCBBaseDataFromNode(*pNode, *pWarrant);

  return pWarrant;
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_ATTACHEDWARRANTCONVERTIBLEBOND_ROOT, 
                               AttachedWarrantCB);
