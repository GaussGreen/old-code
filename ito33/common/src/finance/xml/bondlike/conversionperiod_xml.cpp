/////////////////////////////////////////////////////////////////////////////
// Name:        conversionperiod_xml.cpp
// Purpose:     Restore conversion period objects from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: conversionperiod_xml.cpp,v 1.9 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/bondlike/conversionperiod.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/bondlike/conversionperiod.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/bondlike/cocotype.h"

using namespace ito33::finance;

namespace ito33
{

namespace XML
{ 

bool Restore(const xml::node& node,
             shared_ptr<ConversionPeriod>& pConversionPeriod)
{
  if(strcmp(node.get_name(), XML_TAG_BONDLIKE_CONVERSIONPERIOD_ROOT) != 0)
    return false;

  Date startDate
    = GetDateFromName(node, XML_TAG_BONDLIKE_STARTDATE);

  Date endDate
    = GetDateFromName(node, XML_TAG_BONDLIKE_ENDDATE);

  double dRatio
    = GetDoubleFromName(node, XML_TAG_BONDLIKE_RATIO);

  pConversionPeriod = make_ptr( new ConversionPeriod
                                    ( startDate, endDate, dRatio ) );

  double dCash
    = GetDoubleFromName(node, XML_TAG_BONDLIKE_CONVERSIONPERIOD_CASH);

  pConversionPeriod->SetCash(dCash);

  xml::node::const_iterator
    pnodeTrigger = node.find(XML_TAG_BONDLIKE_TRIGGERRATE);


  if(pnodeTrigger != node.end())
  {
    double dTrigger = GetDoubleFromNode(*pnodeTrigger);

    CoCoType cocoType = GetEnumFromName
                              (
                                node,
                                XML_TAG_BONDLIKE_COCOTYPE,
                                SIZEOF(g_coCoTypes),
                                g_coCoTypes
                              );
    
    double dChangeRate
      = GetDoubleFromName(node, XML_TAG_BONDLIKE_CHANGERATE);

    double dExtremeTrigger
      = GetDoubleFromName(node, XML_TAG_BONDLIKE_EXTREMETRIGGERRATE);

    pConversionPeriod->SetCoCo
                  ( dTrigger, cocoType, dChangeRate, dExtremeTrigger);
  }

  return true;

}

} // namespace XML

} // namespace ito33

