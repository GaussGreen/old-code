/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/callschedule.h
// Purpose:     Names of elements and attributes used in XML for a callschedule
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: callschedule.h,v 1.18 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/callschedule.h
    @brief Contains the names of the elements used in the XML description of a
           callschedule

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_CALLSCHEDULE_H_
#define _ITO33_XML_FINANCE_BONDLIKE_CALLSCHEDULE_H_

#include "ito33/sharedptr.h"

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

  class ITO33_DLLDECL CallSchedule;
  class ITO33_DLLDECL CallPeriod;

}

namespace XML
{

/**
  Restore CallSchedule in given node, if possible.

  @param node given node, normally root node of ConvertibleBond
  @param pCallSchedule (output)
  @return true if the function succeeds, false if no.
  */
bool Restore(const xml::node& node,
             shared_ptr<finance::CallSchedule>& pCallSchedule);

/**
  Restore CallPeriod in given node, if possible.

  @param node given node which may be root node of call period
  @param pCallPeriod (output)
  @return true if the function succeeds, false if no.
  */
bool Restore(const xml::node& node,
             shared_ptr<finance::CallPeriod>& pCallPeriod);


}

}

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_BONDLIKE_CALLSCHEDULE_ROOT                  "call_schedule"

#define XML_TAG_BONDLIKE_CALLSCHEDULE_CALLPERIODS           "call_periods"

#define XML_TAG_BONDLIKE_CALLPERIOD_ROOT                    "call_period"
 
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_CALLSCHEDULE_H_

