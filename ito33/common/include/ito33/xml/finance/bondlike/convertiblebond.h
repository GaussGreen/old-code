/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/convertiblebond.h
// Purpose:     Names of elements and attributes used in XML for 
//              a convertible bond
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: convertiblebond.h,v 1.17 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/convertiblebond.h
    @brief Contains the names of the elements used in the XML description of a
           convertible bond

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_CONVERTIBLEBOND_H_
#define _ITO33_XML_FINANCE_BONDLIKE_CONVERTIBLEBOND_H_

#include "ito33/sharedptr.h"

/**
   @name Tag name macros
*/
//@{
#define XML_TAG_CONVERTIBLEBOND_ROOT                "convertible_bond"
//@}

namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ConvertibleBond;
}

namespace XML
{

  /// Restore function for convertible bond
  bool Restore(const xml::node& node, 
               shared_ptr<finance::ConvertibleBond>& pcb);
}

}

#endif // _ITO33_XML_FINANCE_BONDLIKE_CONVERTIBLEBOND_H_
