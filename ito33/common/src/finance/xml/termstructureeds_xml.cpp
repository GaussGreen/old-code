/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/termstructureeds_xml.cpp
// Purpose:     Restore term structure eds object from XML document
// Author:      ITO33
// Created:     2005/02/10
// RCS-ID:      $Id: termstructureeds_xml.cpp,v 1.3 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"
#include "ito33/date.h"
#include "ito33/useexception.h"

#include "ito33/finance/termstructureeds.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/termstructurederivative.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/eds.h"
#include "ito33/xml/finance/cashflowstream_all.h"
#include "ito33/xml/finance/termstructure.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;

ITO33_FORCE_LINK_THIS_MODULE(termstructureeds_xml);
                             

namespace ito33
{

namespace XML
{

bool Restore(const xml::node& node, 
             shared_ptr<finance::TermStructureEDS>& pEDSCurve)
{
  xml::node::const_iterator pNodeRoot;

  if( (pNodeRoot = node.find(XML_TAG_TERMSTRUCTURE_EDS)) == node.end() )
  {
    return false;
  }

  pEDSCurve = shared_ptr<TermStructureEDS>(new TermStructureEDS);
  xml::node::const_iterator pNodeEle;
  for(pNodeEle = pNodeRoot->begin();
      pNodeEle != pNodeRoot->end();
      pNodeEle++)
  {
    shared_ptr<finance::EDS> pEDS;
    if( Restore( *pNodeEle, pEDS) )
      pEDSCurve->Add( pEDS );
  }

  return true;
}

} // namespace XML

} // namespace ito33

