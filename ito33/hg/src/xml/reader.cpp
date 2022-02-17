/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/xml/reader.cpp
// Purpose:     code for reading pricing HG XML files
// Created:     2005/04/15
// RCS-ID:      $Id: reader.cpp,v 1.7 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/finance/modeloutput.h"

#include "hg/xml/underlyingprocess.h"
#include "hg/xml/common.h"
#include "hg/xml/reader.h"

#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/modeloutput.h"
#include "ito33/xml/finance/underlyingprocess.h"

#include "ito33/xml/read.h"

#include "ito33/hg/theoreticalmodel.h"

#include "ito33/finance/derivative_visitor.h"

namespace ito33
{
  using namespace XML;

namespace hg
{


// ----------------------------------------------------------------------------
// reading theoretical model information
// ----------------------------------------------------------------------------

void
XML::Reader::ReadTheoreticalModel(shared_ptr<TheoreticalModel>& pModel) const
{
  xml::node nodeModel = GetNodeByName(GetRootNode(), XML_TAG_MODEL);

  xml::node nodeUP = GetNodeByName(nodeModel, XML_TAG_UNDERLYING_PROCESS);

  shared_ptr<UnderlyingProcess> 
    pUnderlyingProcess = ReadUnderlyingProcess(nodeUP);

  pModel = make_ptr( new TheoreticalModel(pUnderlyingProcess) );

  xml::node::iterator nodeSR = nodeModel.find(XML_TAG_HG_SR);

  if ( nodeSR != nodeModel.end() )
    pModel->SetSharpeRatio( GetDoubleFromNode(*nodeSR) );
}


// ----------------------------------------------------------------------------
// read output
// ----------------------------------------------------------------------------

bool XML::Reader::ReadOutput(finance::ModelOutput& output) const
{
  const xml::node node = GetRootNode();
  
  xml::node::const_iterator pNodeFound = node.find(XML_TAG_OUTPUT);

  if ( pNodeFound == node.end() )
    return false;

  Restore(*pNodeFound, output);

  return true;
}


} // namespace hg

} // namespace ito33
