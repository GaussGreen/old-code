/////////////////////////////////////////////////////////////////////////////
// Name:        xml/ihg_pricingreader.cpp
// Purpose:     code for reading pricing IHG XML files
// Author:      ITO33
// Created:     2004/12/01
// RCS-ID:      $Id: ihg_pricingreader.cpp,v 1.11 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ihg/xml/pricingreader.h"
#include "ihg/xml/common.h"
#include "ihg/xml/modeloutput.h"
#include "ihg/xml/underlyingprocess.h"

#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/modeloutput.h"
#include "ito33/xml/finance/underlyingprocess.h"

#include "ito33/xml/read.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/modeloutput.h"

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

void
PricingReader::ReadTheoreticalModel(shared_ptr<TheoreticalModel>& pModel) const
{
  xml::node nodeModel = GetNodeByName(GetMainNode(), XML_TAG_MODEL);

  xml::node nodeUP = GetNodeByName(nodeModel, XML_TAG_UNDERLYING_PROCESS);

  shared_ptr<UnderlyingProcess> 
    pUnderlyingProcess = ReadUnderlyingProcess(nodeUP);

  pModel = make_ptr( new TheoreticalModel(pUnderlyingProcess) );
}


// ----------------------------------------------------------------------------
// read output
// ----------------------------------------------------------------------------

bool PricingReader::ReadOutput(finance::ModelOutput& output) const
{
  const xml::node nodeRoot = GetMainNode();
  xml::node::const_iterator pNodeFound = nodeRoot.find(XML_TAG_OUTPUT);

  if ( pNodeFound == nodeRoot.end() )
    return false;

  ito33::XML::Restore(*pNodeFound, output);

  return true;
}
