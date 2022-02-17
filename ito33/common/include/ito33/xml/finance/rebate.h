/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/rebate.h
// Purpose:     Names of elements and attributes used in XML for rebate
// Created:     2005/07/05
// RCS-ID:      $Id: rebate.h,v 1.3 2006/03/24 10:18:27 pedro Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/rebate.h
    @brief Contains the names of the elements used in the XML description of
           rebate.
 */

#ifndef _ITO33_XML_FINANCE_REBATE_H_
#define _ITO33_XML_FINANCE_REBATE_H_

#include "ito33/finance/rebatetype.h"
#include "ito33/enum_values_names.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_REBATE_TYPE                       "rebate_type"

#define XML_REBATE_TYPE_UP                        "immediate"
#define XML_REBATE_TYPE_DOWN                      "at_maturity"

//@}

//=============================================================================
//=              helper for read/write                                                      =
//=============================================================================
namespace ito33
{

/// Mapping between macro name and finance type
const EnumValuesNames<finance::RebateType> 
g_rebateTypes[] =
{
  { XML_REBATE_TYPE_UP, finance::Rebate_Immediate },
  { XML_REBATE_TYPE_DOWN, finance::Rebate_AtMaturity }
};


} // namespace ito33

#endif // #ifndef _ITO33_XML_FINANCE_REBATE_H_
