/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/write_vector.h
// Purpose:     functions for writing vector in XML format
// Author:      ZHANG Yunzhi
// Created:     04.12.11
// RCS-ID:      $Id: write_vector.h,v 1.3 2006/03/23 09:27:14 yann Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/write_vector.h
    @brief functions for writing vector in XML format
 */

#ifndef _ITO33_XML_WRITE_VECTOR_H_
#define _ITO33_XML_WRITE_VECTOR_H_

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/xml/write.h"

namespace ito33
{

namespace XML
{
/// Dump vector
template<class T>
void DumpVector(Tag& tagParent,
                const std::vector<T>& pValues,
                const char* strArrayName, 
                const char* strElementName)
{
  Tag tagRoot = tagParent.Element(strArrayName);

  for(size_t n = 0; n < pValues.size(); n++)
    tagRoot.Element(strElementName)(pValues[n]);
}
} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_WRITE_VECTOR_H_


