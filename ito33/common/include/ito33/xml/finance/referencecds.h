/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/referencecds.h
// Purpose:     Names of elements and attributes in XML for a reference cds
// Created:     2006/05/17
// RCS-ID:      $Id: referencecds.h,v 1.2 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/referencecds.h
    @brief XML tags for reference cds

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_REFERENCECDS_H_
#define _ITO33_XML_FINANCE_REFERENCECDS_H_

#include "ito33/common.h"

/**
   @name name macros
*/
//@{

#define XML_TAG_REFERENCECDS_ROOT             "reference_cds"

#define XML_TAG_REFERENCECDS_ISSUEDATE        "issue_date"

#define XML_TAG_REFERENCECDS_MATURITY_MONTHS  "months_to_maturity"

#define XML_TAG_REFERENCECDS_SPREAD           "spread"
//@}

namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ReferenceCDS;
}

namespace XML
{
  /// Restore reference cds from an xml node
  bool Restore(const xml::node& node, 
               shared_ptr<finance::ReferenceCDS>& pReferenceCDS);
}
  
}
#endif // _ITO33_XML_FINANCE_REFERENCECDS_H_

