/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_volflat_hrwithtimecomponent.cpp
// Purpose:     implement for ParametrizationVolFlatHRWithTimeComponent class
// Author:      Wang
// Created:     03/11/04
// RCS-ID:      $Id: parametrization_volflat_hrwithspotcomponentpower.cpp,v 1.13 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/option.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/basket_goodtype.h"

#include "ihg/volflathrspotcomponentpowercalibrator.h"
#include "ihg/hrspotcomponentpowercalibrator.h"

#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/version.h"
#include "ito33/ihg/error.h"

#include "ito33/ihg/parametrization_volflat_hrwithspotcomponentpower.h"

#include "ito33/xml/write.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/parametrization_visitor.h"
#include "ihg/xml/common.h"

extern const ito33::finance::Error
  ITO33_CALIBRATION_FAIL,  
  ITO33_EMPTY_CDSCURVE;

extern const ito33::ihg::Error
  ITO33_IHG_CALIBRATION_DEGENERATE;

namespace ito33
{

namespace ihg
{


void ParametrizationVolFlatHRWithSpotComponentPower::Calibrate
     (const finance::Derivative& derivative,
      const finance::TermStructureCDS& tsCDS)
{
  CHECK_COND(!tsCDS.GetAll().empty(), ITO33_EMPTY_CDSCURVE);

  // Write the calibration parameters to XML file when debug requested
  if ( IsDebugOutputEnabled() )
  {
    std::ofstream ofs(GetDebugOutputFile().c_str());
    ito33::XML::RootTag tagRoot(XML_TAG_IHG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_IHG_ROOT_VERSION, ITO33_IHG_VERSION_DOT_STRING);

    {
      tsCDS.GetAll().front()->GetSessionData()->Dump(tagRoot);
    }    

    {
      ito33::XML::Tag tagParam(XML_TAG_IHG_CALIBRATION, tagRoot);
      Dump(tagParam);
      tsCDS.Dump(tagParam);
      derivative.Dump(tagParam);
    }
  }

  VolFlatHRSpotComponentPowerCalibrator calibrator;
  
  try
  {
    //check that the cds termstructure contains at least two elements
    CheckSize(tsCDS,2);

    calibrator.Calibrate(derivative, tsCDS);
  }
  catch(ito33::numeric::Exception)
  {
    // Even if the calibration fails, save the last calibration values
    m_pVolatility = calibrator.GetVolatility();

    m_pHazardRate = calibrator.GetHazardRate();

    m_pHRSpotComponentPower = calibrator.GetHRSpotComponentPower();

    // If beta is close to zero, the code likely failed because the data 
    // is degenerate, and a time-only parametrization should be used
    if ( fabs(m_pHRSpotComponentPower->GetBeta()) < 1.e-3)
      throw EXCEPTION(ITO33_IHG_CALIBRATION_DEGENERATE);
    else
      throw EXCEPTION(ITO33_CALIBRATION_FAIL);

  }

  // Save the calibrated values
  m_pVolatility = calibrator.GetVolatility();

  m_pHazardRate = calibrator.GetHazardRate();

  m_pHRSpotComponentPower = calibrator.GetHRSpotComponentPower();
}

void ParametrizationVolFlatHRWithSpotComponentPower::CalibrateWithOptionAndCDSs
     (const finance::Option& option, const finance::TermStructureCDS& tsCDS)
{
  Calibrate(option, tsCDS);
}

void ParametrizationVolFlatHRWithSpotComponentPower::Calibrate
     (const finance::BasketGoodType& basket)
{
  CalibrateWithOptionAndCDSs(*basket.GetOption(),
                             *basket.GetTermStructureCDS());
}

void ParametrizationVolFlatHRWithSpotComponentPower::Dump
(ito33::XML::Tag& tagParent) const
{
  {
  // create the root tag
  ito33::XML::Tag tagParametrization
    (XML_TAG_PARAMETRIZATION_VOLFLATHRWITHSPOTCOMPONENTPOWER_ROOT, tagParent);
  }

  // nothing else to do.  No input data in this class.
}


void ParametrizationVolFlatHRWithSpotComponentPower::Visit
(ParametrizationVisitor& visitor) const
{
  visitor.OnParametrizationVolFlatHRWithSpotComponentPower(*this);
}


} // namespace ihg

} // namespace ito33
