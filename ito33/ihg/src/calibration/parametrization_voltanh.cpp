/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_voltanh.cpp
// Purpose:     implement for parametrization for tanh volatiltiy
// Created:     2005/02/04
// RCS-ID:      $Id: parametrization_voltanh.cpp,v 1.7 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/option.h"
#include "ito33/finance/sessiondata.h"

#include "ihg/voltanhcalibrator_newton.h"

#include "ito33/ihg/volatilitytanh.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/version.h"

#include "ito33/ihg/parametrization_voltanh.h"

#include "ito33/xml/write.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/parametrization_visitor.h"
#include "ihg/xml/common.h"

extern const ito33::Error ITO33_NULL_PARAM;
extern const ito33::finance::Error ITO33_CALIBRATION_FAIL;

namespace ito33
{

namespace ihg
{


typedef VolTanhCalibratorNewton VolTanhCalibrator ;

ParametrizationVolTanh::ParametrizationVolTanh
                        ( const shared_ptr<HazardRate>& pHR )
{
  CHECK_COND(pHR, ITO33_NULL_PARAM);

  m_pHazardRate = pHR;
}

shared_ptr<VolatilityTanh> 
ParametrizationVolTanh::CalibrateWithOptions
(const finance::Option& option1, const finance::Option& option2)
{
  if ( IsDebugOutputEnabled() )
  {
    std::ofstream ofs(GetDebugOutputFile().c_str());
    ito33::XML::RootTag tagRoot(XML_TAG_IHG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_IHG_ROOT_VERSION, ITO33_IHG_VERSION_DOT_STRING);

    {
      option1.GetSessionData()->Dump(tagRoot);
    }    

    {
      ito33::XML::Tag tagParam(XML_TAG_IHG_CALIBRATION, tagRoot);
      Dump(tagParam);
      option1.Dump(tagParam);
      option2.Dump(tagParam);
    }
  }

  shared_ptr<VolatilityTanh> pVol;

  try
  {
    // The shift to the tanh function is determined by the average of the
    // two strikes
    double dStrike1 = option1.GetStrike();
    double dStrike2 = option2.GetStrike();
    
    double dScale, dS0 = 0.5 * (dStrike1 + dStrike2);;
    if ( dStrike1 < dStrike2 )
      dScale = 2. / (dStrike2 - dStrike1);
    else
      dScale = 2. / (dStrike1 - dStrike2);

    VolTanhCalibrator calibrator(dScale, dS0);

    if (dStrike1 < dStrike2)
      pVol = calibrator.Calibrate(option1, option2, m_pHazardRate);
    else
      pVol = calibrator.Calibrate(option2, option1, m_pHazardRate);
  }
  catch(ito33::numeric::Exception)
  {
    throw EXCEPTION(ITO33_CALIBRATION_FAIL);

  } // end catch

  return pVol;
}

void ParametrizationVolTanh::Dump(ito33::XML::Tag& tagParent) const
{
  // create the root tag
  ito33::XML::Tag 
    tagParametrization(XML_TAG_PARAMETRIZATION_VOLTANH_ROOT, tagParent);

  if (m_pHazardRate)
    m_pHazardRate->Dump(tagParametrization);
} 

void ParametrizationVolTanh::Visit(ParametrizationVisitor& visitor) const
{
  visitor.OnParametrizationVolTanh(*this);
}


} // namespace ihg

} // namespace ito33
