/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/timeunit.h
// Purpose:     Names of elements and attributes used in XML for day 
//              counting convention.
// Author:      Yann d'Halluin
// Created:     2004-06-25
// RCS-ID:      $Id: timeunit.h,v 1.5 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/timeunit.h
    @brief Contains the names of the elements used in the XML description of
           TimeUnit as well as reader functions

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_TIMEUNIT_H_
#define _ITO33_XML_FINANCE_TIMEUNIT_H_

#include "ito33/finance/timeunit.h"
#include "ito33/enum_values_names.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_TIMEUNIT "time_unit"

#define XML_TAG_VALUE_TIMEUNIT_DAY           "Day"
#define XML_TAG_VALUE_TIMEUNIT_MONTH         "Month"
#define XML_TAG_VALUE_TIMEUNIT_YEAR          "Year"

//@}

namespace ito33
{

/// Mapping between macro name and finane type
const EnumValuesNames<finance::TimeUnit> 
g_timeUnits[] =
{
  { XML_TAG_VALUE_TIMEUNIT_DAY, finance::TimeUnit_Day },
  { XML_TAG_VALUE_TIMEUNIT_MONTH, finance::TimeUnit_Month},
  { XML_TAG_VALUE_TIMEUNIT_YEAR, finance::TimeUnit_Year}
};

}

namespace ito33
{

namespace XML
{

/// Dump helper
inline const char* GetNameOfTimeUnit(finance::TimeUnit unit)
{
  return GetNameFromEnumValue(unit, SIZEOF(g_timeUnits), g_timeUnits);
}

}

}
#endif // _ITO33_XML_FINANCE_TIMEUNIT_H_

