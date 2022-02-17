/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/xml/calibrationreader.cpp
// Purpose:     code for reading calibration HG XML files
// Created:     2005/05/20
// RCS-ID:      $Id: calibrationreader.cpp,v 1.4 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/xml/read.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/derivative.h"

#include "hg/xml/common.h"
#include "hg/xml/underlyingprocess.h"
#include "hg/xml/parametrization.h"
#include "hg/xml/calibrationreader.h"

#include "ito33/finance/derivative.h"
#include "ito33/finance/derivatives.h"

#include "ito33/hg/parametrization.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::hg::XML;

namespace ito33
{

namespace hg
{

shared_ptr<Parametrization> 
XML::CalibrationReader::ReadParametrization() const
{
  xml::node 
    nodeCalibration = GetNodeByName
                      (GetRootNode(), XML_TAG_HG_CALIBRATION);

  // expected an underlying process as initial guess and as 
  // "parametrization type"
  xml::node 
    nodeUP = GetNodeByName
             (
               GetNodeByName(nodeCalibration, XML_TAG_PARAMETRIZATION_INIT),
               XML_TAG_UNDERLYING_PROCESS
             );

  shared_ptr<UnderlyingProcess> 
    pUnderlyingProcess = ReadUnderlyingProcess(nodeUP);

  shared_ptr<Parametrization>
    pParametrization( new Parametrization(pUnderlyingProcess) );

  return pParametrization;
}

shared_ptr<finance::Derivatives> 
XML::CalibrationReader::ReadDerivatives() const
{
  shared_ptr<finance::SessionData> pSessionData( ReadSessionData() );

  shared_ptr<finance::Derivatives> pDerivatives(new finance::Derivatives);

  xml::node 
    nodeCalibration = GetNodeByName
                      (GetRootNode(), XML_TAG_HG_CALIBRATION);

  // Read in weights, if any.  These are needed when adding to the 
  // derivatives container below
  xml::node::const_iterator 
    searchNode = nodeCalibration.find(XML_TAG_DERIVATIVEWEIGHTS);

  std::vector<double> pdWeights;
  if ( searchNode != nodeCalibration.end() )
  {
    xml::node nodeDerivWeights = *searchNode;

    xml::node::const_iterator iter;
    const xml::node::const_iterator endWeight = nodeDerivWeights.end();
    for ( iter = nodeDerivWeights.begin(); iter != endWeight; ++iter )
    {
      if ( strcmp(iter->get_name(), XML_TAG_DERIVATIVEWEIGHT) == 0 )
      {
        double dWeight = GetDoubleFromNode(*iter);
        pdWeights.push_back(dWeight);
      }
    }
  } // if derivative weights are specified

  // Now read in the derivatives
  size_t nWeightCounter = 0;
  size_t nNbWeights = pdWeights.size();
  xml::node nodeDerivs = GetNodeByName(nodeCalibration, XML_TAG_DERIVATIVES);

  const xml::node::const_iterator end = nodeDerivs.end();
  for ( xml::node::const_iterator i = nodeDerivs.begin(); i != end; ++i )
  {
    shared_ptr<finance::Derivative> 
      pDerivative( DerivativeFactory::Create(i->get_name(), &(*i) ) );

    if (pDerivative)
    {     
      pDerivative->SetSessionData(pSessionData);
      if (nNbWeights > 0 && nWeightCounter < nNbWeights)
        pDerivatives->AddWithWeight(pDerivative, pdWeights[nWeightCounter++]);
      else
        pDerivatives->Add(pDerivative);
    }
  }

  return pDerivatives;
}


} // namespace hg

} // namespace ito33
