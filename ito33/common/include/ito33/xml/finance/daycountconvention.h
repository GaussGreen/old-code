/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/daycountconvention.h
// Purpose:     Names of elements and attributes used in XML for day 
//              counting convention.
// Author:      Yann d'Halluin
// Created:     2004-06-25
// RCS-ID:      $Id: daycountconvention.h,v 1.11 2006/06/08 09:45:35 nabil Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/daycountconvention.h
    @brief Contains the names of the elements used in the XML description of
           day count convention as well as reader functions

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_DAYCOUNTCONVENTION_H_
#define _ITO33_XML_FINANCE_DAYCOUNTCONVENTION_H_

#include "ito33/date.h"
#include "ito33/enum_values_names.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_DAYCOUNTCONVENTION "day_count_convention"

#define XML_TAG_VALUE_DCC_30360                   "30360"    
#define XML_TAG_VALUE_DCC_30E360                  "30360european" 
#define XML_TAG_VALUE_DCC_30U360                  "30U360" 
#define XML_TAG_VALUE_DCC_ActAct                  "actact" 
//#define XML_TAG_VALUE_DCC_ActAct_ISDA             "actact_ISDA" 
#define XML_TAG_VALUE_DCC_Act360                  "act360" 
#define XML_TAG_VALUE_DCC_Act365                  "act365" 
#define XML_TAG_VALUE_DCC_Act365L                 "act365L" 
// Day count conventions with NO EOM (Non End Of Month)
#define XML_TAG_VALUE_DCC_30360_NO_EOM            "30360_NOEOM"    
#define XML_TAG_VALUE_DCC_30E360_NO_EOM           "30360european_NOEOM" 
#define XML_TAG_VALUE_DCC_30U360_NO_EOM           "30U360_NOEOM" 
#define XML_TAG_VALUE_DCC_ActAct_NO_EOM           "actact_NOEOM" 
//#define XML_TAG_VALUE_DCC_ActAct_ISDA_NO_EOM      "actact_ISDA_NOEOM" 
#define XML_TAG_VALUE_DCC_Act360_NO_EOM           "act360_NOEOM" 
#define XML_TAG_VALUE_DCC_Act365_NO_EOM           "act365_NOEOM" 
#define XML_TAG_VALUE_DCC_Act365L_NO_EOM          "act365L_NOEOM" 

//@}

namespace ito33
{

/// Mapping between macro name and finance type
const EnumValuesNames<Date::DayCountConvention> 
g_dayCountConventions[] =
{
  { XML_TAG_VALUE_DCC_30360, Date::DayCountConvention_30360 },
  { XML_TAG_VALUE_DCC_30E360, Date::DayCountConvention_30E360 },
  { XML_TAG_VALUE_DCC_30U360, Date::DayCountConvention_30U360 },
  { XML_TAG_VALUE_DCC_ActAct, Date::DayCountConvention_ActAct },
  //{ XML_TAG_VALUE_DCC_ActAct_ISDA, Date::DayCountConvention_ActAct_ISDA },
  { XML_TAG_VALUE_DCC_Act360, Date::DayCountConvention_Act360 },
  { XML_TAG_VALUE_DCC_Act365, Date::DayCountConvention_Act365 },
  { XML_TAG_VALUE_DCC_Act365L, Date::DayCountConvention_Act365L },
  // Day count conventions with NO EOM (Non End Of Month)
  { XML_TAG_VALUE_DCC_30360_NO_EOM, Date::DayCountConvention_30360_NO_EOM },
  { XML_TAG_VALUE_DCC_30E360_NO_EOM, Date::DayCountConvention_30E360_NO_EOM },
  { XML_TAG_VALUE_DCC_30U360_NO_EOM, Date::DayCountConvention_30U360_NO_EOM },
  { XML_TAG_VALUE_DCC_ActAct_NO_EOM, Date::DayCountConvention_ActAct_NO_EOM },
  /*{ XML_TAG_VALUE_DCC_ActAct_ISDA_NO_EOM, 
      Date::DayCountConvention_ActAct_ISDA_NO_EOM },*/
  { XML_TAG_VALUE_DCC_Act360_NO_EOM, Date::DayCountConvention_Act360_NO_EOM },
  { XML_TAG_VALUE_DCC_Act365_NO_EOM, Date::DayCountConvention_Act365_NO_EOM },
  { XML_TAG_VALUE_DCC_Act365L_NO_EOM, Date::DayCountConvention_Act365L_NO_EOM }
};

}

namespace ito33
{

namespace XML
{


/// Dump helper 
inline const char* GetNameOfDayCountConvention(Date::DayCountConvention dcc)
{
  return GetNameFromEnumValue(dcc, SIZEOF(g_dayCountConventions),
                              g_dayCountConventions);
}


}

}

#endif // _ITO33_XML_FINANCE_DAYCOUNTCONVENTION_H_
