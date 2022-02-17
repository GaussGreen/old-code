/////////////////////////////////////////////////////////////////////////////
// Name:        callschedule_xml.cpp
// Purpose:     Restore bondliketerms objects from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: callschedule_xml.cpp,v 1.18 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/call.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/bondlike/callschedule.h"
#include "ito33/xml/finance/bondlike/call.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/common.h"

using namespace ito33::finance;

namespace ito33
{

namespace XML
{ 
  

bool Restore(const xml::node& node,
             shared_ptr<finance::CallSchedule>& pCallSchedule)
{
  xml::node::const_iterator
    pnodeRoot = node.find(XML_TAG_BONDLIKE_CALLSCHEDULE_ROOT);
  if( pnodeRoot == node.end() )
    return false;

  const xml::node& nodeRoot = *pnodeRoot;

  pCallSchedule = shared_ptr<finance::CallSchedule>(new CallSchedule);
  
  RestoreBondLikeCallData(nodeRoot, *pCallSchedule);

  xml::node::const_iterator pnodeFound;

  pnodeFound = nodeRoot.find(XML_TAG_BONDLIKE_CALLSCHEDULE_CALLPERIODS);

  if(pnodeFound != node.end())
  {
    xml::node::const_iterator pnodePeriod;
    for(pnodePeriod = pnodeFound->begin();
        pnodePeriod != pnodeFound->end();
        pnodePeriod++)
    {
      // each node is a call period except the first and the last one
      shared_ptr<CallPeriod> pCallPeriod;
      if(Restore(*pnodePeriod, pCallPeriod))
        pCallSchedule->AddCallPeriod(pCallPeriod);
    }
  }

  return true;
}


bool Restore(const xml::node& node,
             shared_ptr<finance::CallPeriod>& pCallPeriod)
{
  if (strcmp(node.get_name(), XML_TAG_BONDLIKE_CALLPERIOD_ROOT) != 0)
    return false;

  xml::node::const_iterator pSearchNode;

  Date startDate = GetDateFromName(node, XML_TAG_BONDLIKE_STARTDATE);

  Date endDate = GetDateFromName(node, XML_TAG_BONDLIKE_ENDDATE);

  // Check if strike was set. Default to 1.0.
  pSearchNode = node.find(XML_TAG_FINANCE_STRIKE);

  if (pSearchNode != node.end())
  {
    double dStrike = GetDoubleFromNode(*pSearchNode);

    pCallPeriod = CallPeriod::CreateWithStrike
                              (startDate, endDate, dStrike);
  }

  // Check if yield to call was set. This resets the strike
  pSearchNode = node.find(XML_TAG_FINANCE_YIELD);

  if (pSearchNode != node.end())
  {
    double dYield = GetDoubleFromNode(*pSearchNode);

    pCallPeriod = CallPeriod::CreateWithYield
                          (startDate, endDate, dYield);
  }
  
  // Check if trigger rate was set
  pSearchNode = node.find(XML_TAG_BONDLIKE_TRIGGERRATE);

  if (pSearchNode != node.end())
    pCallPeriod->SetTrigger( GetDoubleFromNode(*pSearchNode) );

  return true;
}

} // namespace XML

} // namespace ito33

