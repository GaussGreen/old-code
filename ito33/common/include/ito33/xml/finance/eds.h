/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/eds.h
// Purpose:     Names of elements and attributes used in XML for a EDS.
// Created:     2005/01/26
// RCS-ID:      $Id: eds.h,v 1.8 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/xml/finance/eds.h
   @brief Contains the elements of the EDS
 */

#ifndef _ITO33_XML_FINANCE_EDS_H_
#define _ITO33_XML_FINANCE_EDS_H_

#include "ito33/sharedptr.h"

#define XML_TAG_EDS_ROOT "eds"

#define XML_TAG_EDS_SPREADSTREAM "spread_stream"

namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL  EDS;
}


namespace XML
{
  /// Restore EDS from XML node
  bool Restore(const xml::node& node, shared_ptr<finance::EDS>& pEDS);
}
  

} // namespace ito33

#endif // #ifndef _ITO33_XML_FINANCE_EDS_H_
