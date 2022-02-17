/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/conversionperiod.h
// Purpose:     Names of elements and attributes used in XML for a callschedule
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: conversionperiod.h,v 1.14 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/conversionperiod.h
    @brief Contains the names of the elements used in the XML description of a
           conversionperiod

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_CONVERSIONPERIOD_H_
#define _ITO33_XML_FINANCE_BONDLIKE_CONVERSIONPERIOD_H_

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
  class ITO33_DLLDECL ConversionPeriod;
}

namespace XML
{

/**
  Restore ConversionPeriod in given node, if possible.

  @param node given node which may be root node of call period
  @param pConversionPeriod (output)
  @return true if the function succeeds, false if not.
  */
bool Restore(const xml::node& node,
             shared_ptr<finance::ConversionPeriod>& pConversionPeriod);

} // namespace XML

}

/**
   @name Tag name macros
*/
//@{
#define XML_TAG_BONDLIKE_CONVERSIONPERIOD_ROOT             "conversion_period"

#define XML_TAG_BONDLIKE_CONVERSIONPERIOD_CASH             "cash_value"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_CONVERSIONPERIOD_H_

