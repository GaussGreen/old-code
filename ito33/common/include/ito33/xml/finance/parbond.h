/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/parbond.h
// Purpose:     Names of elements and attributes used in XML for a parbond.
// Author:      Yann d'Halluin
// Created:     2004-06-25
// RCS-ID:      $Id: parbond.h,v 1.7 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/parbond.h
    @brief Contains the recovery rate of the parbond

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_PARBOND_H_
#define _ITO33_XML_FINANCE_PARBOND_H_

#include "ito33/sharedptr.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_PARBOND_ROOT         "par_bond"

#define XML_TAG_PARBOND_YTM          "risk_free_YTM"
#define XML_TAG_PARBOND_SPREAD       "spread"
#define XML_TAG_PARBOND_MATURITY     "maturity"
#define XML_TAG_PARBOND_CONTRACTING_DATE     "contracting_date"

//@}

namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ParBond;
}

namespace XML
{
  /// Restore parbond from xml node
  bool Restore(const xml::node& node, shared_ptr<finance::ParBond>& pParBond);
}
  
}
#endif // _ITO33_XML_FINANCE_PARBOND_H_
