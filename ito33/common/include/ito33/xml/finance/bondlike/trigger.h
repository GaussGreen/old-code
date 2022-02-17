/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/trigger.h
// Purpose:     Names of elements and attributes used in XML for trigger
// Created:     2004/09/13
// RCS-ID:      $Id: trigger.h,v 1.11 2006/03/24 10:18:27 pedro Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/trigger.h
    @brief Contains the names of the elements used in the XML description of a
           trigger
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_TRIGGER_H_
#define _ITO33_XML_FINANCE_BONDLIKE_TRIGGER_H_

#include "ito33/enum_values_names.h"

#include "ito33/finance/bondlike/triggerincurrencyof.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_BONDLIKE_TAPO               "trigger_as_percentage_of"
#define XML_VALUE_BONDLIKE_TAPO_PRINCIPAL   "principal"
#define XML_VALUE_BONDLIKE_TAPO_ISSUE_PRICE "issue_price"
#define XML_VALUE_BONDLIKE_TAPO_CLAIM       "claim"

#define XML_TAG_BONDLIKE_TRIGGERINCURRENCYOF               "trigger_defined_in_currency_of"
#define XML_VALUE_BONDLIKE_TRIGGERINCURRENCYOF_UNDERLYING  "underlying"
#define XML_VALUE_BONDLIKE_TRIGGERINCURRENCYOF_BOND        "bond"

//@}

namespace ito33
{

/// Mapping between xml tag name and finance variable type
const EnumValuesNames<finance::TriggerInCurrencyOf> 
  g_triggerInCurrencyOfs[] =
    {
      { XML_VALUE_BONDLIKE_TRIGGERINCURRENCYOF_UNDERLYING,
                ito33::finance::TriggerInCurrencyOf_Underlying},
      { XML_VALUE_BONDLIKE_TRIGGERINCURRENCYOF_BOND, 
                ito33::finance::TriggerInCurrencyOf_Derivative}
    };

/// Mapping between xml tag name and finance variable type
const EnumValuesNames<finance::TriggerAsPercentageOf>
  g_triggerAsPercentageOfs[] =
    {
      { XML_VALUE_BONDLIKE_TAPO_PRINCIPAL,
        finance::TriggerAsPercentageOf_Principal },
      { XML_VALUE_BONDLIKE_TAPO_ISSUE_PRICE,
        finance::TriggerAsPercentageOf_IssuePrice },
      { XML_VALUE_BONDLIKE_TAPO_CLAIM,
        finance::TriggerAsPercentageOf_Claim }
    };

}

#endif // _ITO33_XML_FINANCE_BONDLIKE_TRIGGER_H_
