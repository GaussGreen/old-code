/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/cds_xml.cpp
// Purpose:     Restore cds object from XML document
// Author:      Yann d'halluin
// Created:     2004-06-25
// RCS-ID:      $Id: parbond_xml.cpp,v 1.3 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"
#include "ito33/date.h"
#include "ito33/useexception.h"

#include "ito33/finance/frequency.h"
#include "ito33/finance/parbond.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/read_frequency.h"
#include "ito33/xml/finance/read_daycountconvention.h"
#include "ito33/xml/finance/parbond.h"
#include "ito33/xml/finance/common.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::finance;

ITO33_FORCE_LINK_THIS_MODULE(parbond_xml);


/**
    Restore a parbond object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
finance::Derivative *ReadParBond(const xml::node *pNode)
{
  double
    dRecoveryRate = GetDoubleFromName(*pNode,XML_TAG_FINANCE_RECOVERYRATE),
    dYTM = GetDoubleFromName(*pNode, XML_TAG_PARBOND_YTM),
    dSpread = GetDoubleFromName(*pNode,XML_TAG_PARBOND_SPREAD);

  Frequency freq = GetFrequencyFromName(*pNode, XML_TAG_PAYMENTFREQUENCY);

  Date::DayCountConvention dcc
    = GetDayCountConventionFromName(*pNode, XML_TAG_DAYCOUNTCONVENTION);

  size_t nMaturity = GetLongFromName(*pNode, XML_TAG_PARBOND_MATURITY);

  Date contractingDate = GetDateFromName(*pNode,
                                         XML_TAG_PARBOND_CONTRACTING_DATE);

  finance::ParBond*
    pParBond( new finance::ParBond(contractingDate, nMaturity, dYTM,
                                   dSpread, freq, dcc, dRecoveryRate) );

  GetOptionalDerivativeDataFromNode(*pNode, *pParBond);

  return pParBond;
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_PARBOND_ROOT, ParBond);


bool ito33::XML::Restore
               (const xml::node& node, shared_ptr<finance::ParBond>& pParBond)
{
  if( strcmp(node.get_name(), XML_TAG_PARBOND_ROOT) != 0)
    return false;

  pParBond = shared_ptr<finance::ParBond>
              ( dynamic_cast<finance::ParBond*>(ReadParBond(&node)) );

  return true;
}

