/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/read_rebatetype.h
// Purpose:     Helper function to read the type of a rebate
// Created:     2005/07/05
// RCS-ID:      $Id: read_rebatetype.h,v 1.2 2006/03/23 09:27:15 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/read_rebatetype.h
    @brief Helper function to read the type of a rebate.
 */

#ifndef _ITO33_XML_FINANCE_READ_REBATETYPE_H_
#define _ITO33_XML_FINANCE_READ_REBATETYPE_H_

#include "ito33/xml/read.h"
#include "ito33/xml/finance/rebate.h"

namespace ito33
{

namespace XML
{

/// RebateType reader
inline finance::RebateType
GetRebateTypeFromName(const xml::node& node, const char* name)
{
  return GetEnumFromName(node, name, SIZEOF(g_rebateTypes), g_rebateTypes);
}


} // namespace XML

} // namespace ito33

#endif // #ifndef _ITO33_XML_FINANCE_READ_REBATETYPE_H_
