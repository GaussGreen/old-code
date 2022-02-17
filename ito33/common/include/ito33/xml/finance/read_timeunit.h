/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/read_timeunit.h
// Purpose:     Helper function to read a time unit 
// Created:     2005/04/19
// RCS-ID:      $Id: read_timeunit.h,v 1.2 2006/03/23 09:27:15 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/read_timeunit.h
    @brief Contains the names of the elements used in the XML description of
           TimeUnit as well as reader functions
 */

#ifndef _ITO33_XML_FINANCE_READ_TIMEUNIT_H_
#define _ITO33_XML_FINANCE_READ_TIMEUNIT_H_

#include "ito33/xml/read.h"
#include "ito33/xml/finance/timeunit.h"

namespace ito33
{

namespace XML
{

/// Reader for time unit
inline finance::TimeUnit 
GetTimeUnitFromName(const xml::node& node, const char* name)
{
  return GetEnumFromName(node, name, SIZEOF(g_timeUnits), g_timeUnits);
}


}

}

#endif // _ITO33_XML_FINANCE_READ_TIMEUNIT_H_
