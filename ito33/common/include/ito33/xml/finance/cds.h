/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/cds.h
// Purpose:     Names of elements and attributes used in XML for a cds.
// Author:      Yann d'Halluin
// Created:     2004-06-25
// RCS-ID:      $Id: cds.h,v 1.12 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/cds.h
    @brief Contains the recovery rate of the cds

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_CDS_H_
#define _ITO33_XML_FINANCE_CDS_H_

#include "ito33/common.h"

/**
   @name name macros
*/
//@{

#define XML_TAG_CDS_ROOT         "cds"

#define XML_TAG_CDS_SPREADSTREAM "spread_stream"

//@}

namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL CDS;
}

namespace XML
{
  /// Restore cds from an xml node
  bool Restore(const xml::node& node, shared_ptr<finance::CDS>& pCDS);
}
  
}
#endif // _ITO33_XML_FINANCE_CDS_H_

