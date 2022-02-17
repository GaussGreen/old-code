/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/termstructureoption_xml.cpp
// Purpose:     Restore an option term structure from XML document
// Created:     2005/03/04
// RCS-ID:      $Id: termstructureoption_xml.cpp,v 1.5 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/termstructureoption.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/option.h"
#include "ito33/xml/finance/termstructure.h"

ITO33_FORCE_LINK_THIS_MODULE(termstructureoption_xml);
                             
namespace ito33
{

  using namespace finance;

namespace XML
{


bool Restore(const xml::node& node, 
             shared_ptr<finance::TermStructureOption>& pOptionCurve)
{
  xml::node::const_iterator pNodeRoot;

  if ( (pNodeRoot = node.find(XML_TAG_TERMSTRUCTURE_OPTION)) == node.end() )
    return false;

  pOptionCurve = shared_ptr<TermStructureOption>(new TermStructureOption);
  
  xml::node::const_iterator pNodeEle;
  for (pNodeEle = pNodeRoot->begin(); pNodeEle != pNodeRoot->end(); pNodeEle++)
  {
    shared_ptr<Option> pOption;
    
    if ( ito33::XML::Restore( *pNodeEle, pOption) )
      pOptionCurve->Add( pOption );
  }

  return true;
}


} // namespace XML

} // namespace ito33
