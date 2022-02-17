/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/sharedependentconversion.h
// Purpose:     Names of elements and attributes used in XML for share
//              dependent conversion
// Author:      ITO 33
// Created:     2004/12/31
// RCS-ID:      $Id: sharedependentconversion.h,v 1.8 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/sharedependentconversion.h
    @brief Contains the names of the elements used in the XML description of 
           share dependent conversion

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_SHAREDEPENDENTCONVERSION_H_
#define _ITO33_XML_FINANCE_BONDLIKE_SHAREDEPENDENTCONVERSION_H_

#include "ito33/sharedptr.h"

namespace xml
{
  class node;
}

//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ShareDependentConversion;
}

namespace XML
{

/**
    Restore ShareDependentConversion in given node, if possible.

    @param node given node, normally root node of ConvertibleBond
    @param pShareDepConversion (output)

    @return true if the function succeeds, false if not.
 */
bool Restore(const xml::node& node,
        shared_ptr<finance::ShareDependentConversion>& pShareDepConversion);
}
}

/**
   @name Tag name macros
*/
//@{
#define XML_TAG_BONDLIKE_SHAREDEPENDENT_ROOT       "share_dependent_conversion"

#define XML_TAG_BONDLIKE_SHAREDEPENDENT_BASERATIO    "base_ratio"

#define XML_TAG_BONDLIKE_SHAREDEPENDENT_CAPRATIO     "cap_ratio"

#define XML_TAG_BONDLIKE_SHAREDEPENDENT_RESETDATE    "reset_date"

#define XML_TAG_BONDLIKE_SHAREDEPENDENT_SHAREFACTOR  "incremental_share_factor"

#define XML_TAG_BONDLIKE_SHAREDEPENDENT_FIXEDSTRIKE  "fixed_strike"

#define XML_TAG_BONDLIKE_SHAREDEPENDENT_CURRENT_RATIO "current_ratio"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_SHAREDEPENDENTCONVERSION_H_
