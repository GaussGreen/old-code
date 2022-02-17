/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/swappayofftype.h
// Purpose:     Elements and attributes in XML for variance swap payoffs
// Created:     2006/07/12
// RCS-ID:      $Id: swappayofftype.h,v 1.2 2006/08/03 20:05:12 dave Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/swappayofftype.h
    @brief 
 */

#ifndef _ITO33_XML_FINANCE_SWAP_PAYOFF_TYPE_H_
#define _ITO33_XML_FINANCE_SWAP_PAYOFF_TYPE_H_

#include "ito33/finance/swappayofftype.h"
#include "ito33/enum_values_names.h"

//=============================================================================
//=              tag name macros                                              =
//=============================================================================

#define XML_TAG_VARIANCESWAP_PAYOFF_TYPE                "swap_payoff_type"

#define XML_VALUE_VARIANCESWAP_PAYOFF_TYPE_STANDARD     "standard"

#define XML_VALUE_VARIANCESWAP_PAYOFF_TYPE_GAMMA        "gamma"

#define XML_VALUE_VARIANCESWAP_PAYOFF_TYPE_CONDITIONAL  "conditional"

#define XML_VALUE_VARIANCESWAP_PAYOFF_TYPE_CALL         "call"

#define XML_VALUE_VARIANCESWAP_PAYOFF_TYPE_PUT          "put"


//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{

const EnumValuesNames<finance::SwapPayoffType> g_swapPayoffTypes[] =
{
  { XML_VALUE_VARIANCESWAP_PAYOFF_TYPE_STANDARD, finance::SwapPayoff_Standard },
  { XML_VALUE_VARIANCESWAP_PAYOFF_TYPE_GAMMA, finance::SwapPayoff_Gamma },
  { XML_VALUE_VARIANCESWAP_PAYOFF_TYPE_CONDITIONAL, finance::SwapPayoff_Conditional },
  { XML_VALUE_VARIANCESWAP_PAYOFF_TYPE_CALL, finance::SwapPayoff_Call },
  { XML_VALUE_VARIANCESWAP_PAYOFF_TYPE_PUT, finance::SwapPayoff_Put },
};

}

#endif // #define _ITO33_XML_FINANCE_SWAP_PAYOFF_TYPE_H_
