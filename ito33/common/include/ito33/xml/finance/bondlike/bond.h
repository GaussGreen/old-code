/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/bond.h
// Purpose:     Names of elements and attributes used in XML for a bond
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: bond.h,v 1.9 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/bond.h
    @brief Contains the names of the elements used in the XML description of a
           BOND

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_BOND_H_
#define _ITO33_XML_FINANCE_BONDLIKE_BOND_H_

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
  class ITO33_DLLDECL Bond;
}

namespace XML
{

/**
   Get from XML the optinal data, such as call schedule, for a specific bond

   @param bond Bond or ConvertibleBond
   @param node Bond tag or ConvertibleBond tag
 */
void GetOptionalBondDataFor(finance::Bond& bond, const xml::node& node);

}
}

/**
   @name Tag name macros
*/
//@{
#define XML_TAG_BOND_ROOT                      "bond"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_BOND_H_
