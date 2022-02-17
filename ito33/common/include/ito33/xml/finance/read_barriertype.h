/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/read_barriertype.h
// Purpose:     Helper function to read the type of a barrier
// Created:     2005/07/05
// RCS-ID:      $Id: read_barriertype.h,v 1.2 2006/03/23 09:27:15 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/read_barriertype.h
    @brief Helper function to read the type of a barrier.
 */

#ifndef _ITO33_XML_FINANCE_READ_BARRIERTYPE_H_
#define _ITO33_XML_FINANCE_READ_BARRIERTYPE_H_

#include "ito33/xml/read.h"
#include "ito33/xml/finance/barrier.h"

namespace ito33
{

namespace XML
{

/// Reader for barrier type
inline finance::BarrierType
GetBarrierTypeFromName(const xml::node& node, const char* name)
{
  return GetEnumFromName(node, name, SIZEOF(g_barrierTypes), g_barrierTypes);
}


} // namespace XML

} // namespace ito33

#endif // #ifndef _ITO33_XML_FINANCE_READ_BARRIERTYPE_H_
