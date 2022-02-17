/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/frequency.h
// Purpose:     Names of elements and attributes used in XML for Frequency.
// Created:     2005/02/22
// RCS-ID:      $Id: frequency.h,v 1.10 2006/06/05 18:01:21 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/frequency.h
    @brief Contains the names of the elements used in the XML description of
           Frequency.
 */

#ifndef _ITO33_XML_FINANCE_FREQUENCY_H_
#define _ITO33_XML_FINANCE_FREQUENCY_H_

#include "ito33/finance/frequency.h"

#include "ito33/enum_values_names.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_PAYMENTFREQUENCY  "payment_frequency"

#define XML_TAG_VALUE_FREQUENCY_ANNUAL     "annual"
#define XML_TAG_VALUE_FREQUENCY_SEMIANNUAL "semiannual"
#define XML_TAG_VALUE_FREQUENCY_QUARTERLY  "quarterly"
#define XML_TAG_VALUE_FREQUENCY_BIMONTHLY  "bimonthly"
#define XML_TAG_VALUE_FREQUENCY_MONTHLY    "monthly"
#define XML_TAG_VALUE_FREQUENCY_DAILY      "daily"
#define XML_TAG_VALUE_FREQUENCY_WEEKLY     "weekly"
#define XML_TAG_VALUE_FREQUENCY_BIWEEKLY   "biweekly"

//@}

//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{

/// Mapping between macro name and finance type
const EnumValuesNames<finance::Frequency> 
g_frequencys[] =
{
  { XML_TAG_VALUE_FREQUENCY_ANNUAL, finance::Frequency_Annual },
  { XML_TAG_VALUE_FREQUENCY_SEMIANNUAL, finance::Frequency_SemiAnnual },
  { XML_TAG_VALUE_FREQUENCY_QUARTERLY, finance::Frequency_Quarterly },
  { XML_TAG_VALUE_FREQUENCY_BIMONTHLY, finance::Frequency_BiMonthly },
  { XML_TAG_VALUE_FREQUENCY_MONTHLY, finance::Frequency_Monthly }
};

}

namespace ito33
{

namespace XML
{

/// Dump helper for frequency
inline const char* GetNameOfFrequency(finance::Frequency unit)
{
  return GetNameFromEnumValue(unit, SIZEOF(g_frequencys), g_frequencys);
}


}

}

#endif // _ITO33_XML_FINANCE_FREQUENCY_H_
