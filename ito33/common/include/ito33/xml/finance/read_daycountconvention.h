/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/read_daycountconvention.h
// Purpose:     Helper function to read a counting convention.
// Created:     2005/04/19
// RCS-ID:      $Id: read_daycountconvention.h,v 1.3 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/read_daycountconvention.h
    @brief Helper function to read a counting convention.

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_READ_DAYCOUNTCONVENTION_H_
#define _ITO33_XML_FINANCE_READ_DAYCOUNTCONVENTION_H_

#include "ito33/xml/read.h"
#include "ito33/xml/finance/daycountconvention.h"

namespace ito33
{

namespace XML
{


/// Reader for day count convention
inline Date::DayCountConvention GetDayCountConventionFromName
                              (const xml::node& node, const char* name)
{
  return GetEnumFromName(node, name,
                         SIZEOF(g_dayCountConventions), g_dayCountConventions);
}


}

}

#endif // _ITO33_XML_FINANCE_DAYCOUNTCONVENTION_H_
