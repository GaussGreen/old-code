/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/generalizedpepslikecall.h
// Purpose:     Names of elements and attributes used in XML for a 
//              generalized peps like call
// Author:      ZHANG Yunzhi
// Created:     2005-3-29
// RCS-ID:      $Id: generalizedpepslikecall.h,v 1.5 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/generalizedpepslikecall.h
    @brief Contains the names of the elements used in the XML description of a
           generalized peps like call
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_H_
#define _ITO33_XML_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_H_

#include "ito33/sharedptr.h"
#include "ito33/enum_values_names.h"

#include "ito33/finance/bondlike/generalizedpepslikecalltype.h"
//=============================================================================
//=              readers                                                      =
//=============================================================================

namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  /// Forward declaration
  class ITO33_DLLDECL GeneralizedPEPSLikeCall;
}

namespace XML
{

/**
  Restore GeneralizedPEPSLikeCall in given node, if possible.

  @param node given node, normally root node of GeneralizedPEPSLike
  @param pGeneralizedPEPSLikeCall (output)
  @return true if the function succeeds, false if no.
  */
bool Restore
      (
        const xml::node& node,
        shared_ptr<finance::GeneralizedPEPSLikeCall>& pGeneralizedPEPSLikeCall
      );

}

}

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_ROOT  \
                                                   "generalized_peps_like_call"

// call type
#define XML_TAG_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_TYPE             "call_type"
#define XML_VALUE_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_TYPE_FIXED   "fixed_share"
#define XML_VALUE_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_TYPE_VARIABLE \
                                                               "variable_share"
//@}

namespace ito33
{

/// Mapping between xml tag name and finance type
const EnumValuesNames<finance::GeneralizedPEPSLikeCallType> 
g_GeneralizedPEPSLikeCallTypes[] =
{
  { XML_VALUE_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_TYPE_FIXED, 
      finance::GeneralizedPEPSLikeCallType_FixedShare }, 
  { XML_VALUE_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_TYPE_VARIABLE, 
      finance::GeneralizedPEPSLikeCallType_VariableShare }
};

}

#endif // _ITO33_XML_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_H_
