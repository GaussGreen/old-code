/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/optiontype.h
// Purpose:     Names of elements and attributes used in XML for option type.
// Created:     March 29, 2005
// RCS-ID:      $Id: optiontype.h,v 1.4 2006/03/24 10:18:27 pedro Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/optiontype.h
    @brief Contains the names of the elements used in the XML description of
           OptionType.
 */

#ifndef _ITO33_XML_FINANCE_OPTIONTYPE_H_
#define _ITO33_XML_FINANCE_OPTIONTYPE_H_

#include "ito33/finance/optiontype.h"
#include "ito33/enum_values_names.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_OPTION_TYPE                        "type"

#define XML_VALUE_OPTION_TYPE_PUT                  "put"
#define XML_VALUE_OPTION_TYPE_CALL                 "call"

//@}

//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{

/// Mapping between macro name and finance type
const EnumValuesNames<finance::OptionType> 
g_optionTypes[] =
{
  { XML_VALUE_OPTION_TYPE_PUT, finance::Option_Put },
  { XML_VALUE_OPTION_TYPE_CALL, finance::Option_Call }
};


}

#endif // _ITO33_XML_FINANCE_OPTIONTYPE_H_
