/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/underlyingprocess.cpp
// Purpose:     Restore underlying process from xml document
// Created:     2006/06/02
// RCS-ID:      $Id: underlyingprocess_xml.cpp,v 1.1 2006/06/22 10:11:50 nabil Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/xml/read.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/underlyingprocess.h"

namespace ito33
{

namespace XML
{

double GetPostDefaultVolatilityFromNode(const xml::node &node)
{
  double dPostDefaultVolatility = 0.0;

  xml::node::const_iterator pNodeFound;

  // Volatility after default is optional 
  pNodeFound = node.find(XML_TAG_POST_DEFAULT_VOLATILITY);

  if ( pNodeFound != node.end() )
    dPostDefaultVolatility = GetDoubleFromNode(*pNodeFound);

  return dPostDefaultVolatility;
}


} // namespace XML

} // namespace ito33
