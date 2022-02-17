/////////////////////////////////////////////////////////////////////////////
// Name:        callfixedshare_xml.cpp
// Purpose:     Restore bondliketerms objects from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-12-03
// RCS-ID:      $Id: callfixedshare_xml.cpp,v 1.4 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/bondlike/callfixedshare.h"
#include "ito33/finance/bondlike/call.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/bondlike/callfixedshare.h"
#include "ito33/xml/finance/bondlike/call.h"
#include "ito33/xml/finance/bondlike/common.h"

using namespace ito33::finance;

namespace ito33
{

namespace XML
{ 
  

bool Restore(const xml::node& node,
             shared_ptr<finance::CallFixedShare>& pCallFixedShare)
{
  xml::node::const_iterator
    pnodeRoot = node.find(XML_TAG_BONDLIKE_CALLFIXEDSHARE_ROOT);

  if( pnodeRoot == node.end() )
    return false;

  const xml::node& nodeRoot = *pnodeRoot;

  pCallFixedShare
    = shared_ptr<finance::CallFixedShare>
         (new CallFixedShare
                (
                  GetDateFromName(nodeRoot, XML_TAG_BONDLIKE_STARTDATE),
                  GetDateFromName(nodeRoot, XML_TAG_BONDLIKE_ENDDATE),
                  GetDoubleFromName(nodeRoot, XML_TAG_BONDLIKE_RATIO)
                )
         );

  xml::node::const_iterator
    pnodeTrigger = nodeRoot.find(XML_TAG_BONDLIKE_TRIGGERRATE);

  if (pnodeTrigger != node.end())
    pCallFixedShare->SetTrigger( GetDoubleFromNode(*pnodeTrigger) );

  RestoreBondLikeCallData(nodeRoot, *pCallFixedShare);

  return true;
}


} // namespace XML

} // namespace ito33

