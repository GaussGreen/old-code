/////////////////////////////////////////////////////////////////////////////
// Name:        sharedependentconversion_xml.cpp
// Purpose:     Restore share dependent conversion object from XML document
// Author:      ITO 33
// Created:     2005-03-04
// RCS-ID:      $Id: sharedependentconversion_xml.cpp,v 1.5 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/bondlike/sharedependentconversion.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/bondlike/conversion.h"
#include "ito33/xml/finance/bondlike/sharedependentconversion.h"
#include "ito33/xml/finance/bondlike/cocotype.h"

using namespace ito33::finance;

namespace ito33
{

namespace XML
{ 

bool Restore(const xml::node& node,
        shared_ptr<finance::ShareDependentConversion>& pShareDepConv)
{
  xml::node::const_iterator
    pnodeRoot = node.find(XML_TAG_BONDLIKE_SHAREDEPENDENT_ROOT);
  if( pnodeRoot == node.end() )
    return false;

  const xml::node& nodeRoot = *pnodeRoot;

  Date startDate 
    = GetDateFromName(nodeRoot, XML_TAG_BONDLIKE_STARTDATE);

  Date endDate 
    = GetDateFromName(nodeRoot, XML_TAG_BONDLIKE_ENDDATE);

  double dBaseRatio
    = GetDoubleFromName(nodeRoot, XML_TAG_BONDLIKE_SHAREDEPENDENT_BASERATIO);

  double dIncrementalShareFactor
    = GetDoubleFromName(nodeRoot, XML_TAG_BONDLIKE_SHAREDEPENDENT_SHAREFACTOR);
 
  pShareDepConv = make_ptr
      ( new ShareDependentConversion
            ( startDate, endDate, dBaseRatio, dIncrementalShareFactor ) );

  xml::node::const_iterator
    pnodeCap = nodeRoot.find(XML_TAG_BONDLIKE_SHAREDEPENDENT_CAPRATIO);

  if ( pnodeCap != node.end() )
    pShareDepConv->SetCapRatio( GetDoubleFromNode(*pnodeCap) );

  xml::node::const_iterator pNodeResetDate =
    nodeRoot.find(XML_TAG_BONDLIKE_SHAREDEPENDENT_RESETDATE);

  if ( pNodeResetDate != node.end() )
  {
    Date resetDate = GetDateFromNode(*pNodeResetDate);
    pShareDepConv->SetResetDate(resetDate);

    xml::node::const_iterator pNodeCurrentConversionRatio =
    nodeRoot.find(XML_TAG_BONDLIKE_SHAREDEPENDENT_CURRENT_RATIO);

    if ( pNodeCurrentConversionRatio != node.end() )
    {
      double dCurrentConvRatio=GetDoubleFromNode(*pNodeCurrentConversionRatio);
      pShareDepConv->SetCurrentRatio(dCurrentConvRatio);

    }
  }

  xml::node::const_iterator pnodeFixedStrike = 
    nodeRoot.find(XML_TAG_BONDLIKE_SHAREDEPENDENT_FIXEDSTRIKE);

  if ( pnodeFixedStrike != node.end() )
  {  
    double dStrike = GetDoubleFromNode(*pnodeFixedStrike);
    pShareDepConv->SetFixedStrike(dStrike);
  }

  //read keep accrued and forfeit coupon
  RestoreCommonConversionData(nodeRoot, *pShareDepConv);

  xml::node::const_iterator
    pnodeTrigger = nodeRoot.find(XML_TAG_BONDLIKE_TRIGGERRATE);

  if ( pnodeTrigger != node.end() )
  {  
    double dTrigger = GetDoubleFromNode(*pnodeTrigger);

    CoCoType cocoType = GetEnumFromName
                              (
                                nodeRoot,
                                XML_TAG_BONDLIKE_COCOTYPE,
                                SIZEOF(g_coCoTypes),
                                g_coCoTypes
                              );
    
    double dChangeRate
      = GetDoubleFromName(nodeRoot, XML_TAG_BONDLIKE_CHANGERATE);

    double dExtremeTrigger
      = GetDoubleFromName(nodeRoot, XML_TAG_BONDLIKE_EXTREMETRIGGERRATE);

    // islasttriggermet is optional, defaults to false
    bool bLastTrigger = false;
    xml::node::const_iterator pSearchNode;
    pSearchNode = nodeRoot.find(XML_TAG_BONDLIKE_ISLASTTRIGGERMET);
    if ( pSearchNode != nodeRoot.end() )
      bLastTrigger = GetBoolFromName(nodeRoot, 
                                     XML_TAG_BONDLIKE_ISLASTTRIGGERMET);

    pShareDepConv->SetCoCo
        ( dTrigger, cocoType, dChangeRate, dExtremeTrigger, bLastTrigger);
  }

  return true;
}

} // namespace XML

} // namespace ito33

