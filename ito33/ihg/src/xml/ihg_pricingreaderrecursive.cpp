/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/xml/ihg_pricingreaderrecursive.cpp
// Purpose:     code for reading pricing IHG XML files in 
//              ihg/tests/testsuite//bondlike/xmlfiles
// Author:      ITO 33 Canada
// Created:     April 4, 2005
// RCS-ID:      $Id: ihg_pricingreaderrecursive.cpp,v 1.7 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <string.h>

#include "ihg/xml/pricingreaderrecursive.h"
#include "ihg/xml/common.h"
#include "ihg/xml/modeloutput.h"
#include "ihg/xml/underlyingprocess.h"

#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/sessiondata.h"
#include "ito33/xml/finance/underlyingprocess.h"

#include "ito33/xml/read.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/finance/derivative_visitor.h"

#include <xmlwrapp/init.h>
#include <xmlwrapp/document.h>
#include <xmlwrapp/tree_parser.h>

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::ihg::XML;

using ito33::ihg::TheoreticalModel;

// ----------------------------------------------------------------------------
// reading theoretical model information
// ----------------------------------------------------------------------------
void PricingReaderRecursive::ReadTheoreticalModel(
                                    shared_ptr<TheoreticalModel>& pModel) const
{
  xml::node 
    nodeModel = GetNodeByNameRecursive(GetMainNode(), XML_TAG_MODEL);

  xml::node nodeUP = GetNodeByName(nodeModel, XML_TAG_UNDERLYING_PROCESS);

  shared_ptr<UnderlyingProcess> 
    pUnderlyingProcess = ReadUnderlyingProcess(nodeUP);

  pModel = make_ptr( new TheoreticalModel(pUnderlyingProcess) );

} //end ReadTheoreticalModel

// ----------------------------------------------------------------------------
// reading derivatives information
// ----------------------------------------------------------------------------

void PricingReaderRecursive::ReadDerivatives(
                                     finance::DerivativeVisitor& visitor) const
{
  xml::node nodeDerivs 
    = GetNodeByNameRecursive(GetRootNode(), XML_TAG_IHG_DERIVATIVES);

  const xml::node::const_iterator end = nodeDerivs.end();
  for ( xml::node::const_iterator i = nodeDerivs.begin(); i != end; ++i )
  {
    const xml::node& node = *i;
    const char *name = node.get_name();
    shared_ptr<finance::Derivative>
      derivative(DerivativeFactory::Create(name, &node));
    if ( derivative )
    {
      derivative->Visit(visitor);
    }
    //else: ignore unknown derivatives
  }

}

void PricingReaderRecursive::ReadDerivative(
      shared_ptr<finance::Derivative>& pDerivative) const
{
  xml::node nodeDerivs 
    = GetNodeByNameRecursive(GetRootNode(), XML_TAG_IHG_DERIVATIVES);

  const xml::node::const_iterator end = nodeDerivs.end();
  for ( xml::node::const_iterator i = nodeDerivs.begin(); i != end; ++i )
  {
    const xml::node& node = *i;
    const char *name = node.get_name();
    shared_ptr<finance::Derivative>
      derivative(DerivativeFactory::Create(name, &node));
    if ( derivative )
    {
      pDerivative = derivative;

      return;
    }
    //else: ignore unknown derivatives
  }

}

// ----------------------------------------------------------------------------
// reading session data information
// ----------------------------------------------------------------------------

shared_ptr<finance::SessionData> PricingReaderRecursive::ReadSessionData() const
{

 const xml::node rootNode = GetRootNode();
 xml::node tmpNode = GetNodeByNameRecursive(rootNode, XML_TAG_SESSIONDATA);

 return ito33::XML::GetSessionDataFromNode(tmpNode);
}
