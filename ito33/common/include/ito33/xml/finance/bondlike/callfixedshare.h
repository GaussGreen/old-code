/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/callfixedshare.h
// Purpose:     Names of elements and attributes used in XML for a call fixed
//                share
// Author:      Yann d'Halluin
// Created:     2004-10-1
// RCS-ID:      $Id: callfixedshare.h,v 1.8 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/callfixedshare.h
    @brief Contains the names of the elements used in the XML description of a
           call fixed share

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_CALLFIXEDSHARE_H_
#define _ITO33_XML_FINANCE_BONDLIKE_CALLFIXEDSHARE_H_

#include "ito33/sharedptr.h"

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
  class ITO33_DLLDECL CallFixedShare;
}

namespace XML
{

/**
  Restore CallFixedShare in given node, if possible.

  @param node given node, normally root node of PEPSLike
  @param pCallFixedShare (output)
  @return true if the function succeeds, false if no.
  */
bool Restore(const xml::node& node,
             shared_ptr<finance::CallFixedShare>& pCallFixedShare);

}

}

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_BONDLIKE_CALLFIXEDSHARE_ROOT               "call_fixed_share"

//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_CALLFIXEDSHARE_H_

