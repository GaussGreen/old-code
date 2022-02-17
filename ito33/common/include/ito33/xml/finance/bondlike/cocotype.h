/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/cocotype.h
// Purpose:     Names of elements and attributes for contingent conversion
// Author:      ITO 33
// Created:     22/02/2005
// RCS-ID:      $Id: cocotype.h,v 1.6 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/cocotype.h
    @brief Contains the names of the elements used in the XML description 
      for contingentconversion

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_COCOTYPE_H_
#define _ITO33_XML_FINANCE_BONDLIKE_COCOTYPE_H_

#include "ito33/enum_values_names.h"
#include "ito33/finance/bondlike/cocotype.h"

/**
   @name Tag name macros
*/
//@{
#define XML_TAG_BONDLIKE_COCOTYPE                          "coco_type"
#define XML_VALUE_BONDLIKE_COCOTYPE_QUARTER_FOREVER        "check_quarterly_and_convert_as_of_check_date"
#define XML_VALUE_BONDLIKE_COCOTYPE_QUARTER_QUARTER        "check_quarterly_and_convert_during_next_quarter"
#define XML_VALUE_BONDLIKE_COCOTYPE_ANYTIME_FOREVER        "check_any_time_and_convert_as_of_check_date"
#define XML_VALUE_BONDLIKE_COCOTYPE_ANYTIME_NOW            "check_any_time_and_convert_on_check_date"

#define XML_TAG_BONDLIKE_CHANGERATE           "change_rate_in_trigger"
#define XML_TAG_BONDLIKE_EXTREMETRIGGERRATE   "extreme_trigger_rate"
#define XML_TAG_BONDLIKE_ISLASTTRIGGERMET     "is_last_trigger_condition_met"
//@}

//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{
/// Declaration of mapping between finance types and macro names
const EnumValuesNames<finance::CoCoType> 
g_coCoTypes[] =
{
  { XML_VALUE_BONDLIKE_COCOTYPE_QUARTER_FOREVER, finance::CoCoType_CheckQuarterlyAndConvertAsOfCheckDate },
  { XML_VALUE_BONDLIKE_COCOTYPE_QUARTER_QUARTER, finance::CoCoType_CheckQuarterlyAndConvertDuringNextQuarter },
  { XML_VALUE_BONDLIKE_COCOTYPE_ANYTIME_FOREVER, finance::CoCoType_CheckAnyTimeAndConvertAsOfCheckDate },
  { XML_VALUE_BONDLIKE_COCOTYPE_ANYTIME_NOW, finance::CoCoType_CheckAnyTimeAndConvertOnCheckDate }
 };

}

#endif // _ITO33_XML_FINANCE_BONDLIKE_CONVERSIONSCHEDULE_H_

