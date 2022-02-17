/////////////////////////////////////////////////////////////////////////////
// Name:        bondliketerms_xml.cpp
// Purpose:     Restore bondliketerms objects from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: call_xml.cpp,v 1.12 2006/06/03 19:39:04 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#include "ito33/finance/bondlike/call.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/bondlike/call.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/bondlike/trigger.h"

namespace ito33
{

namespace XML
{ 
  
  using namespace finance;

void RestoreBondLikeCallData(const xml::node& node, Call& call)
{
  call.SetKeepAccrued
       ( GetBoolFromName(node, XML_TAG_BONDLIKE_KEEPACCRUED) );

  call.SetForfeitCoupon
       ( GetBoolFromName(node, XML_TAG_BONDLIKE_FORFEITCOUPON) );

  xml::node::const_iterator pNodeFound;

  pNodeFound = node.find(XML_TAG_BONDLIKE_NOTICEPERIOD);

  if ( pNodeFound != node.end() )
    call.SetNoticePeriod(GetLongFromNode(*pNodeFound));

  xml::node::const_iterator pNodeMakeWhole;
  if( (pNodeMakeWhole = node.find(XML_TAG_BONDLIKE_CALL_MAKEWHOLE_ROOT))
              != node.end())
  {
    static const EnumValuesNames<MakeWholeType> makeWholeTypes[] =
    {
      { XML_VALUE_BONDLIKE_CALL_MAKEWHOLE_TYPE_COUPON,
                                        MakeWholeType_Coupon },
      { XML_VALUE_BONDLIKE_CALL_MAKEWHOLE_TYPE_PREMIUM,
                                        MakeWholeType_Premium },
    };

    if(GetEnumFromName
          (
            *pNodeMakeWhole,
            XML_TAG_FINANCE_TYPE,
            SIZEOF(makeWholeTypes),
            makeWholeTypes
          )
        == MakeWholeType_Premium
      )
    {
      call.SetPremiumMakeWhole
          (
            GetDoubleFromName
              (
                *pNodeMakeWhole, 
                XML_TAG_BONDLIKE_CALL_MAKEWHOLE_PREMIUM
              )
          );
    }
    else
    {
      call.SetCouponMakeWhole
           (
             GetBoolFromName(*pNodeMakeWhole,
                             XML_TAG_BONDLIKE_CALL_MAKEWHOLE_PVCOUPON)
           );
    }
  }

  pNodeFound = node.find(XML_TAG_BONDLIKE_CALL_TRIGGERPERIOD);

  if ( pNodeFound != node.end() )
  {
    size_t nTriggerPeriod = GetLongFromNode(*pNodeFound);

    size_t nTriggerHistory = GetLongFromName
                             ( node, XML_TAG_BONDLIKE_CALL_TRIGGERHISTORY );

    call.SetTriggerCheckPeriod(nTriggerPeriod, nTriggerHistory);
  }

  // reading of TriggerAsPercentageOf is not symetric to writing.
  // the reason is that sometimes we don't want to write this tag
  // when it is not specified
  pNodeFound = node.find(XML_TAG_BONDLIKE_TAPO);

  if ( pNodeFound != node.end() )
  {
    TriggerAsPercentageOf
      triggerAs = GetEnumFromNode
                  (
                    *pNodeFound,
                    SIZEOF(g_triggerAsPercentageOfs),
                    g_triggerAsPercentageOfs
                  );
    call.SetTriggerAsPercentageOf(triggerAs);  
  
  }


}

} // namespace XML

} // namespace ito33

