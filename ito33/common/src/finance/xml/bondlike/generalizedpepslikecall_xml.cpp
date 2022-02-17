/////////////////////////////////////////////////////////////////////////////
// Name:        generalizedpepslikecall_xml.cpp
// Purpose:     Restore GeneralizedPEPSLikeCall objects from XML document
// Author:      ZHANG Yunzhi
// Created:     2005-03-28
// RCS-ID:      $Id: generalizedpepslikecall_xml.cpp,v 1.2 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/bondlike/generalizedpepslikecall.h"
#include "ito33/finance/bondlike/call.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/bondlike/generalizedpepslikecall.h"
#include "ito33/xml/finance/bondlike/call.h"
#include "ito33/xml/finance/bondlike/common.h"

using namespace ito33::finance;

namespace ito33
{

namespace XML
{ 

bool Restore
        (
          const xml::node& node,
          shared_ptr<finance::GeneralizedPEPSLikeCall>& pGeneralizedPEPSLikeCall
        )
{
  xml::node::const_iterator
    pnodeRoot = node.find(XML_TAG_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_ROOT);

  if( pnodeRoot == node.end() )
    return false;

  const xml::node& nodeRoot = *pnodeRoot;

  pGeneralizedPEPSLikeCall
    = shared_ptr<finance::GeneralizedPEPSLikeCall>
         (new GeneralizedPEPSLikeCall
                (
                  GetDateFromName(nodeRoot, XML_TAG_BONDLIKE_STARTDATE),
                  GetDateFromName(nodeRoot, XML_TAG_BONDLIKE_ENDDATE),
                  GetDoubleFromName(nodeRoot, XML_TAG_BONDLIKE_TRIGGERRATE),
                  GetEnumFromName
                    (
                      nodeRoot,
                      XML_TAG_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_TYPE,
                      SIZEOF(g_GeneralizedPEPSLikeCallTypes),
                      g_GeneralizedPEPSLikeCallTypes
                    )
                )
         );

  RestoreBondLikeCallData(nodeRoot, *pGeneralizedPEPSLikeCall);

  return true;
}


} // namespace XML

} // namespace ito33

