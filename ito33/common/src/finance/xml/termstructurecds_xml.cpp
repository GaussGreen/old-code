/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/termstructurecds_xml.cpp
// Purpose:     Restore term structure cds object from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-09-24
// RCS-ID:      $Id: termstructurecds_xml.cpp,v 1.7 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"
#include "ito33/date.h"
#include "ito33/useexception.h"

#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/referencecds.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/termstructurederivative.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/cds.h"
#include "ito33/xml/finance/referencecds.h"
#include "ito33/xml/finance/cashflowstream_all.h"
#include "ito33/xml/finance/termstructure.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;

ITO33_FORCE_LINK_THIS_MODULE(termstructurecds_xml);
                             

namespace ito33
{

namespace XML
{


bool Restore(const xml::node& node, shared_ptr<TermStructureCDS>& pCDSCurve)
{
  xml::node::const_iterator pNodeRoot;

  if( (pNodeRoot = node.find(XML_TAG_TERMSTRUCTURE_CDS)) == node.end() )
    return false;

  // clear old curve
  pCDSCurve = shared_ptr<TermStructureCDS>(new TermStructureCDS);
  xml::node::const_iterator pNodeEle;
  for(pNodeEle = pNodeRoot->begin();
      pNodeEle != pNodeRoot->end();
      pNodeEle++)
  {
    // Try to read a normal CDS
    shared_ptr<finance::CDS> pCDS;
    if( Restore( *pNodeEle, pCDS) )
      pCDSCurve->Add( pCDS );

    // Try to read a reference CDS
    shared_ptr<finance::ReferenceCDS> pRefCDS;
    if( Restore( *pNodeEle, pRefCDS) )
      pCDSCurve->Add( pRefCDS );
  }

  return true;
}


shared_ptr<TermStructureCDS> GetTermStructureCDSInNode
        (
        const xml::node& node
        )
{
  shared_ptr<TermStructureCDS> pCDSCurve;

  if( ! Restore(node, pCDSCurve) )
  {
    typedef MissingNodeException Exception;
    throw EXCEPTION_MSG(node, XML_TAG_TERMSTRUCTURE_EDS);
  }
  return pCDSCurve;
}


}

}
