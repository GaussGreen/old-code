/////////////////////////////////////////////////////////////////////////////
// Name:        reset_xml.cpp
// Purpose:     Restore resettable convertible Bond object from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-10-21
// RCS-ID:      $Id: reset_xml.cpp,v 1.5 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/triggeraspercentageof.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/termstructure.h"
#include "ito33/xml/finance/bondlike/cb_base.h"
#include "ito33/xml/finance/bondlike/convertiblebond.h"
#include "ito33/xml/finance/bondlike/resetconversionschedule.h"
#include "ito33/xml/finance/bondlike/reset.h"
#include "ito33/xml/finance/bondlike/bondliketerms_all.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::finance;

ITO33_FORCE_LINK_THIS_MODULE(reset_xml);


void GetResetConstructorData
        (
          const xml::node& node,
          shared_ptr<BondTerms>& pBondTerms,
          shared_ptr<ResetConversionSchedule>& pResets
        )
{
  pBondTerms = GetBondTermsFromNode(node);

  pResets = GetResetConversionScheduleFromNode(node);
}


/**
    Restore an Bond object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new bond object to be deleted by the caller
 */
static
finance::Derivative *ReadReset(const xml::node *pNode)
{
  shared_ptr<BondTerms> pBondTerms;
  shared_ptr<ResetConversionSchedule> pResets;
  GetResetConstructorData(*pNode, pBondTerms, pResets);

  Reset *pBond = new Reset(pBondTerms, pResets);

  GetOptionalCBBaseDataFromNode(*pNode, *pBond);

  return pBond;
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_RESET_ROOT, Reset);
