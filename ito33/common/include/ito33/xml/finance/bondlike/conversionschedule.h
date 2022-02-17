/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/conversionschedule.h
// Purpose:     Names of elements and attributes used in XML for a conversionschedule
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: conversionschedule.h,v 1.16 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/conversionschedule.h
    @brief Contains the names of the elements used in the XML description of a
           conversionschedule

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_CONVERSIONSCHEDULE_H_
#define _ITO33_XML_FINANCE_BONDLIKE_CONVERSIONSCHEDULE_H_

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
  class ITO33_DLLDECL ConversionSchedule;
}

namespace XML
{

/**
  Restore ConversionSchedule in given node, if possible.

  @param node given node, normally root node of ConvertibleBond
  @param pConversionSchedule (output)
  @return true if the function succeeds, false if not.
  */
bool Restore(const xml::node& node,
        shared_ptr<finance::ConversionSchedule>& pConversionSchedule);
}
}

/**
   @name Tag name macros
*/
//@{
#define XML_TAG_BONDLIKE_CONVERSIONSCHEDULE_ROOT              "conversion_schedule"

#define XML_TAG_BONDLIKE_CONVERSIONSCHEDULE_CONVERSIONPERIODS "conversion_periods"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_CONVERSIONSCHEDULE_H_
