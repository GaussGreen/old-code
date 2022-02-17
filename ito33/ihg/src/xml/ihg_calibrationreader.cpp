/////////////////////////////////////////////////////////////////////////////
// Name:        xml/ihg_calibrationreader.cpp
// Purpose:     code for reading calibration IHG XML files
// Author:      ITO33
// Created:     2004/12/01
// RCS-ID:      $Id: ihg_calibrationreader.cpp,v 1.14 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/common.h"

#include "ito33/finance/termstructure_enumerator.h"

#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/termstructureparbond.h"
#include "ito33/finance/termstructureeds.h"
#include "ito33/finance/termstructureoption.h"
#include "ito33/finance/derivatives.h"

#include "ito33/ihg/parametrization_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_hrwithspotcomponentpower.h"
#include "ito33/ihg/parametrization_volpower.h"
#include "ito33/ihg/parametrization_volpower_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_volpower_hrpower.h"
#include "ito33/ihg/parametrization_voltanh.h"
#include "ito33/ihg/parametrization_voltanh_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_voltanh_hrwithspotcomponentpower.h"

#include "ito33/ihg/parametrization_volwithtimecomponent.h"

#include "ito33/xml/finance/derivative.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/termstructure.h"
#include "ito33/xml/finance/option.h"
#include "ito33/xml/finance/common.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::ihg::XML;


// ----------------------------------------------------------------------------
// reading calibration information
// ----------------------------------------------------------------------------

void 
CalibrationReader::ReadCalibration
                  (
                    ihg::ParametrizationVisitor& param_visitor,
                    finance::TermStructureEnumerator& termstructures,
                    finance::DerivativeVisitor& deriv_visitor,
                    finance::Derivatives& derivatives) const
{
  xml::node nodeCalibration
                  = GetNodeByName(GetMainNode(), XML_TAG_IHG_CALIBRATION);

  // Try to read each of the different parametrization types
  xml::node::const_iterator pNode;

  // loop over nodes looking for parametrization
  for ( xml::node::const_iterator i = nodeCalibration.begin();
        i != nodeCalibration.end();
        ++i )
  {
    const xml::node& node = *i;
    const char *name = node.get_name();
    shared_ptr<ihg::Parametrization>
      parametrization(ParametrizationFactory::Create(name, &node));
    if ( parametrization )
    {
      parametrization->Visit(param_visitor);
    }
  } // loop over nodes looking for parametrization

  // Also read in the contracts to calibrate.  Try to read any basic
  // contract type, termstructures, and a general derivative list

  shared_ptr<finance::SessionData> pSessionData( ReadSessionData() );

  ReadTermStructures(nodeCalibration, termstructures);
  
  termstructures.SetSessionData(pSessionData);

  for ( xml::node::const_iterator i = nodeCalibration.begin();
        i != nodeCalibration.end();
        ++i )
  {
    const xml::node& node = *i;
    const char *name = node.get_name();
    shared_ptr<finance::Derivative>
      derivative(DerivativeFactory::Create(name, &node));
    if ( derivative )
    {
      derivative->SetSessionData(pSessionData);

      derivative->Visit(deriv_visitor);
    }
  } // loop over nodes looking for derivatives

  xml::node::const_iterator 
    searchNode = nodeCalibration.find(XML_TAG_DERIVATIVES);

  if ( searchNode != nodeCalibration.end() )
  {
    xml::node nodeDerivs = GetNodeByName(nodeCalibration, XML_TAG_DERIVATIVES);
    const xml::node::const_iterator end = nodeDerivs.end();
    for ( xml::node::const_iterator i = nodeDerivs.begin(); i != end; ++i )
    {
      shared_ptr<finance::Derivative> 
        pDerivative( DerivativeFactory::Create(i->get_name(), &(*i) ) );

      if (pDerivative)
      {
        pDerivative->SetSessionData(pSessionData);
        derivatives.Add(pDerivative);
      }
    }

  } // found derivatives node

}
