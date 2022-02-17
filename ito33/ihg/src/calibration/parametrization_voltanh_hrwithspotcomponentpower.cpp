/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_voltanh_hrwithspotcomponentpower.cpp
// Purpose:     implement for parametrization for tanh volatiltiy
// Author:      Ito33
// Created:     2004/12/14
// RCS-ID:      $Id: parametrization_voltanh_hrwithspotcomponentpower.cpp,v 1.8 2006/08/20 09:31:04 wang Exp $
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

#include "ito33/ihg/volatilitytanh.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/version.h"
#include "ito33/ihg/error.h"

#include "ito33/ihg/parametrization_voltanh_hrwithspotcomponentpower.h"

#include "ihg/voltanhhrpowercalibrator.h"

#include "ito33/xml/write.h"
#include "ihg/xml/common.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/parametrization_visitor.h"

extern const ito33::finance::Error ITO33_CALIBRATION_FAIL;

extern const ito33::ihg::Error ITO33_IHG_CALIBRATION_DEGENERATE;

namespace ito33
{

namespace ihg
{


void 
ParametrizationVolTanhHRWithSpotComponentPower::CalibrateWithOptionsAndCDSs
(
  const finance::Option& option1, 
  const finance::Option& option2,
  const finance::TermStructureCDS& tsCDS
)
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
      tsCDS.Dump(tagParam);
    }
  }

  VolTanhHRPowerCalibrator calibrator;

  try
  {
    // The two options need to be passed so that the strike of the first option
    // is at the left of the strike of the second option.
    double dStrike1 = option1.GetStrike();
    double dStrike2 = option2.GetStrike();

    double dScale, dS0 = 0.5 * (dStrike1 + dStrike2);;
    if ( dStrike1 < dStrike2 )
      dScale = 2. / (dStrike2 - dStrike1);
    else
      dScale = 2. / (dStrike1 - dStrike2);

    if ( dStrike1 < dStrike2 )
      calibrator.Calibrate(option1, option2, dScale, dS0, tsCDS);
    else
      calibrator.Calibrate(option2, option1, dScale, dS0, tsCDS);
  }
  catch(ito33::numeric::Exception)
  {
    // Save the results, even if it did not converge
    m_pVolatility = calibrator.GetVolatility();
    m_pHazardRate = calibrator.GetHazardRate();
    m_pHRSpotComponentPower = calibrator.GetHRSpotComponentPower();

    // If beta is close to zero, the code likely failed because the data 
    // is degenerate, and a time-only parametrization should be used
    if ( fabs(m_pHRSpotComponentPower->GetBeta()) < 1.e-3)
      throw EXCEPTION(ITO33_IHG_CALIBRATION_DEGENERATE);
    else
      throw EXCEPTION(ITO33_CALIBRATION_FAIL);

  }//end catch

  // Save the converged results
  m_pVolatility = calibrator.GetVolatility();
  m_pHazardRate = calibrator.GetHazardRate();
  m_pHRSpotComponentPower = calibrator.GetHRSpotComponentPower();
}

void 
ParametrizationVolTanhHRWithSpotComponentPower::Dump
(ito33::XML::Tag& tagParent) const
{
  // create the root tag
  ito33::XML::Tag 
    tagParametrization(XML_TAG_PARAMETRIZATION_VOLTANHHRWITHSPOTCOMPONENTPOWER_ROOT,
                       tagParent);
} 

void 
ParametrizationVolTanhHRWithSpotComponentPower::Visit
(ParametrizationVisitor& visitor) const
{
  visitor.OnParametrizationVolTanhHRWithSpotComponentPower(*this);
}


} // namespace ihg

} // namespace ito33
