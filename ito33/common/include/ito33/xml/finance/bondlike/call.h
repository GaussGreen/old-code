/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/call.h
// Purpose:     Names of elements and attributes used in XML for a call
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: call.h,v 1.15 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/call.h
    @brief Contains the names of the elements used in the XML description of a
           call

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_CALL_H_
#define _ITO33_XML_FINANCE_BONDLIKE_CALL_H_

namespace xml
{
  class node;
}

//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Call;
}

namespace XML
{

/**
    Restores bondlike::Call Data in an object of specific Call class from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the tag in DOM tree presents the specific Call
    @param call the reference to a call object
 */
void RestoreBondLikeCallData
          (
            const xml::node& node,
            finance::Call& call
          );

}

}

/**
   @name Tag name macros
*/
//@{
#define XML_TAG_BONDLIKE_CALL_MAKEWHOLE_ROOT              "makewhole"

#define XML_VALUE_BONDLIKE_CALL_MAKEWHOLE_TYPE_COUPON     "coupon"

#define XML_VALUE_BONDLIKE_CALL_MAKEWHOLE_TYPE_PREMIUM    "premium"

#define XML_TAG_BONDLIKE_CALL_MAKEWHOLE_PREMIUM           "premium_value"
  
#define XML_TAG_BONDLIKE_CALL_MAKEWHOLE_PVCOUPON          "PV_coupon_makewhole"

#define XML_TAG_BONDLIKE_CALL_TRIGGERPERIOD               "trigger_period"

#define XML_TAG_BONDLIKE_CALL_TRIGGERHISTORY              "trigger_history"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_CALL_H_

