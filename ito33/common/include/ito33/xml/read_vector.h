/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/read.h
// Purpose:     Helpers for parsing XML
// Author:      ZHANG yunzhi
// Created:     2004-12-10
// RCS-ID:      $Id: read_vector.h,v 1.4 2006/05/09 21:17:01 dave Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/read_vector.h
    @brief Helper functions used for parsing vectors.
 */

#ifndef _ITO33_XML_READ_VECTOR_H_
#define _ITO33_XML_READ_VECTOR_H_

#include "ito33/xml/read.h"

#include <xmlwrapp/node.h>

namespace ito33
{

namespace XML
{

/// Declaration of GetValueFromNode function
template<class T>
T GetValueFromNode(const xml::node& node);

/// Declaraion of GetValueFromNode for a date
template<> inline
Date GetValueFromNode<Date>(const xml::node& node)
{
  return GetDateFromNode(node);
}

/// Declaration of GetValueFromNode for a double
template<> inline
double GetValueFromNode<double>(const xml::node& node)
{
  return GetDoubleFromNode(node);
}

/// Declaration of GetVectorFromNode
template<class T>
  std::vector<T> GetVectorFromNode
        (
          const xml::node& node,
          const char* strArrayName, 
          const char* strElementName
        )
{
  xml::node::const_iterator pNodeRoot = node.find(strArrayName);

  if(pNodeRoot == node.end())
  {
    typedef MissingNodeException Exception;
    throw EXCEPTION_MSG(node, strArrayName);
  }

  std::vector<T> pValues;
  for (xml::node::const_iterator pNodeElement = pNodeRoot->begin();
        pNodeElement != pNodeRoot->end();
        pNodeElement++)
  {
    if  ( strcmp(pNodeElement->get_name(), strElementName)  == 0 )
      pValues.push_back(GetValueFromNode<T>(*pNodeElement));
  } //end loop over nodes
  return pValues;
}


} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_READ_VECTOR_H_

