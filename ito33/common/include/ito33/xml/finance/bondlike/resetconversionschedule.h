/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/resetconversionschedule.h
// Purpose:     Names of XML elements/attributes for a resetconversionschedule
// Author:      Yann and David
// Created:     2004/10/19
// RCS-ID:      $Id: resetconversionschedule.h,v 1.8 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/resetconversionschedule.h
    @brief Contains the names of the elements used in the XML description of a
           resetconversionschedule
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_RESETCONVERSIONSCHEDULE_H_
#define _ITO33_XML_FINANCE_BONDLIKE_RESETCONVERSIONSCHEDULE_H_

#include "ito33/sharedptr.h"
#include "ito33/enum_values_names.h"
#include "ito33/finance/bondlike/resetflooredby.h"

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
  class ITO33_DLLDECL ResetConversionSchedule;
}

namespace XML
{

/**
  Get ResetConversionSchedule in given node. If the node doesn't have 
  ResetConversionSchedule tag or is bad formatted, the function throws exception.

  @param node given node, normally root node of Reset
  @return shared ptr to ResetConversionSchedule object.
  */
shared_ptr<finance::ResetConversionSchedule>
GetResetConversionScheduleFromNode(const xml::node& node);
}
}


/**
   @name Tag name macros
*/
//@{

#define XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_ROOT        "reset_conversion_schedule"

#define XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_CONVPRICERESET "conversion_price_resets"

#define XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_INITIAL     "initial_conversion_price"

#define XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_CURRENT     "current_conversion_price"

#define XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_START       "conversion_start"

#define XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_END         "conversion_end"

#define XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_CASH        "cash_value"

#define XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_FLOOREDBY            "floored_by"
#define XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_FLOOREDBY_INITIAL    "initial_conversion_price"
#define XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_FLOOREDBY_PREVAILING "prevailing_conversion_price"
//@}

namespace ito33
{
/// Mapping between xml tag name and finance type
const EnumValuesNames<ito33::finance::ResetFlooredBy> g_resetFlooredBy[] =
{
  { XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_FLOOREDBY_INITIAL,
                                      ito33::finance::ResetFlooredBy_InitialConversionPrice},
  { XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_FLOOREDBY_PREVAILING,
                                      ito33::finance::ResetFlooredBy_PrevailingConversionPrice}
};

}
#endif // _ITO33_XML_FINANCE_BONDLIKE_RESETSCHEDULE_H_
