/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/returntype.h
// Purpose:     Names of elements and attributes in XML for return calc method
// Created:     2006/02/21
// RCS-ID:      $Id: returntype.h,v 1.2 2006/04/10 10:55:51 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/returntype.h
    @brief 
 */

#ifndef _ITO33_XML_FINANCE_RETURN_TYPE_H_
#define _ITO33_XML_FINANCE_RETURN_TYPE_H_

#include "ito33/finance/returntype.h"
#include "ito33/enum_values_names.h"

//=============================================================================
//=              tag name macros                                              =
//=============================================================================

#define XML_TAG_VARIANCESWAP_RETURN_TYPE                 "return_method"

#define XML_VALUE_VARIANCESWAP_RETURN_TYPE_ACTUAL        "actual"

#define XML_VALUE_VARIANCESWAP_RETURN_TYPE_LOG           "log"


//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{

const EnumValuesNames<finance::ReturnType> g_returnTypes[] =
  {
    { XML_VALUE_VARIANCESWAP_RETURN_TYPE_ACTUAL, finance::Return_Actual},
    { XML_VALUE_VARIANCESWAP_RETURN_TYPE_LOG, finance::Return_Log},
  };

}

#endif // #define _ITO33_XML_VARIANCESWAP_RETURN_TYPE_H_
