/////////////////////////////////////////////////////////////////////////////
// Name:        conversionschedule_xml.cpp
// Purpose:     Restore conversionschedule objects from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: conversionschedule_xml.cpp,v 1.9 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/bondlike/conversionschedule.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/bondlike/conversionperiod.h"
#include "ito33/xml/finance/bondlike/conversion.h"
#include "ito33/xml/finance/bondlike/conversionschedule.h"
#include "ito33/xml/finance/bondlike/cocotype.h"

using namespace ito33::finance;

namespace ito33
{

namespace XML
{ 
  

bool Restore(const xml::node& node,
             shared_ptr<ConversionSchedule>& pConversionSchedule)
{
  xml::node::const_iterator
    pnodeRoot = node.find(XML_TAG_BONDLIKE_CONVERSIONSCHEDULE_ROOT);
  if( pnodeRoot == node.end() )
    return false;

  const xml::node& nodeRoot = *pnodeRoot;

  pConversionSchedule 
    = shared_ptr<ConversionSchedule>(new ConversionSchedule);
  
  RestoreCommonConversionData(nodeRoot, *pConversionSchedule);

  xml::node::const_iterator pnodeFound;

  // islasttriggermet is optional, defaults to false
  bool bLastTrigger = false;
  pnodeFound = nodeRoot.find(XML_TAG_BONDLIKE_ISLASTTRIGGERMET);
  if ( pnodeFound != nodeRoot.end() )
    bLastTrigger = 
      GetBoolFromName(nodeRoot, XML_TAG_BONDLIKE_ISLASTTRIGGERMET);

  pConversionSchedule->SetIsLastTriggerConditionMet(bLastTrigger);

  pnodeFound
    = nodeRoot.find(XML_TAG_BONDLIKE_CONVERSIONSCHEDULE_CONVERSIONPERIODS);

  if(pnodeFound != node.end())
  {
    xml::node::const_iterator pnodePeriod;
    for(pnodePeriod = pnodeFound->begin();
        pnodePeriod != pnodeFound->end();
        pnodePeriod++)
    {
      // each node is a conversion period except the first and the last one
      shared_ptr<ConversionPeriod> pConversionPeriod;
      if(Restore(*pnodePeriod, pConversionPeriod))
        pConversionSchedule->AddConversionPeriod(pConversionPeriod);
    }
  }

  return true;
}

} // namespace XML

} // namespace ito33

