/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/swaptype.h
// Purpose:     Names of elements and attributes in XML for variance swap types
// Created:     2006/02/21
// RCS-ID:      $Id: swaptype.h,v 1.2 2006/04/10 10:55:51 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/swaptype.h
    @brief 
 */

#ifndef _ITO33_XML_FINANCE_SWAP_TYPE_H_
#define _ITO33_XML_FINANCE_SWAP_TYPE_H_

#include "ito33/finance/swaptype.h"
#include "ito33/enum_values_names.h"

//=============================================================================
//=              tag name macros                                              =
//=============================================================================

#define XML_TAG_VARIANCESWAP_SWAP_TYPE                   "swap_type"

#define XML_VALUE_VARIANCESWAP_SWAP_TYPE_VARIANCE        "variance"

#define XML_VALUE_VARIANCESWAP_SWAP_TYPE_VOLATILITY      "volatility"


//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{

const EnumValuesNames<finance::SwapType> g_swapTypes[] =
{
  { XML_VALUE_VARIANCESWAP_SWAP_TYPE_VARIANCE, finance::Swap_Variance },
  { XML_VALUE_VARIANCESWAP_SWAP_TYPE_VOLATILITY, finance::Swap_Volatility },
};

}

#endif // #define _ITO33_XML_FINANCE_SWAP_TYPE_H_
