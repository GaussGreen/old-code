/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/termstructureparbond_xml.cpp
// Purpose:     Restore term structure parbond object from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-09-24
// RCS-ID:      $Id: termstructureparbond_xml.cpp,v 1.2 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"
#include "ito33/date.h"
#include "ito33/useexception.h"

#include "ito33/finance/termstructureparbond.h"
#include "ito33/finance/parbond.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/termstructurederivative.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/parbond.h"
#include "ito33/xml/finance/cashflowstream_all.h"
#include "ito33/xml/finance/termstructure.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;

ITO33_FORCE_LINK_THIS_MODULE(termstructureparbond_xml);
                             

namespace ito33
{

namespace XML
{

bool Restore(const xml::node& node, 
             shared_ptr<finance::TermStructureParBond>& pParBondCurve)
{
  xml::node::const_iterator pNodeRoot;

  if( (pNodeRoot = node.find(XML_TAG_TERMSTRUCTURE_PARBOND))
                     == node.end() )
    return false;

  pParBondCurve = 
    shared_ptr<TermStructureParBond>(new TermStructureParBond);
  xml::node::const_iterator pNodeEle;
  for(pNodeEle = pNodeRoot->begin();
      pNodeEle != pNodeRoot->end();
      pNodeEle++)
  {
    shared_ptr<finance::ParBond> pParBond;
    if( Restore( *pNodeEle, pParBond) )
      pParBondCurve->Add( pParBond );
  }

  return pParBondCurve;
}


}

}
