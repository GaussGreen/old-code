/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/read_frequency.h
// Purpose:     Helper function to read a Frequency
// Created:     2005/04/19
// RCS-ID:      $Id: read_frequency.h,v 1.2 2006/03/23 09:27:15 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/read_frequency.h
    @brief Helper function to read a Frequency.
 */

#ifndef _ITO33_XML_FINANCE_READ_FREQUENCY_H_
#define _ITO33_XML_FINANCE_READ_FREQUENCY_H_

#include "ito33/xml/read.h"
#include "ito33/xml/finance/frequency.h"

namespace ito33
{

namespace XML
{


/// reader for frequency
inline finance::Frequency GetFrequencyFromName
                          (const xml::node& node, const char* name)
{
  return GetEnumFromName(node, name, SIZEOF(g_frequencys), g_frequencys);
}


}

}

#endif // _ITO33_XML_FINANCE_READ_FREQUENCY_H_
