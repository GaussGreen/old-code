/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_volpower.cpp
// Purpose:     implement for parametrization for power volatiltiy
// Author:      Ito33
// Created:     2004/12/14
// RCS-ID:      $Id: parametrization_volpower.cpp,v 1.8 2006/08/20 09:31:04 wang Exp $
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
#include "ito33/finance/sessiondata.h"

#include "ihg/volpowercalibrator.h"

#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/version.h"

#include "ito33/ihg/parametrization_volpower.h"

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

ParametrizationVolPower::ParametrizationVolPower
( const shared_ptr<HazardRate>& pHR )
{
  if ( !pHR )
    throw EXCEPTION_MSG
          (
            ITO33_NULL_PARAM,
            TRANS("Hazard rate has not been set properly.")
          );

  m_pHazardRate = pHR;

}//ParametrizationVolPower::ParametrizationVolPower


shared_ptr<VolatilityPower> 
ParametrizationVolPower::CalibrateWithOptions
( const finance::Option& option1, const finance::Option& option2)
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


  VolPowerCalibrator calibrator;
  shared_ptr<VolatilityPower> pVol;

  try
  {
    // The first option passed to the calibrator is used to obtain a
    // guess for alpha. The guess will be better if the option is
    // close to at-the-money
    double dStrike1 = option1.GetStrike();
    double dStrike2 = option2.GetStrike();
    double dS0 = option1.GetSessionData()->GetSpotSharePrice();
    if ( fabs(dStrike1-dS0) < fabs(dStrike2-dS0) )
      pVol = calibrator.Calibrate(option1, option2, m_pHazardRate);
    else
      pVol = calibrator.Calibrate(option2, option1, m_pHazardRate);
  }
  catch(ito33::Exception)
  {
    throw EXCEPTION(ITO33_CALIBRATION_FAIL);
  }//end catch

  return pVol;

}//ParametrizationVolPower::CalibrateWithOptions


void ParametrizationVolPower::Dump(ito33::XML::Tag& tagParent) const
{
  // create the root tag
  ito33::XML::Tag 
    tagParametrization(XML_TAG_PARAMETRIZATION_VOLPOWER_ROOT,
                        tagParent);

  if (m_pHazardRate)
    m_pHazardRate->Dump(tagParametrization);
} //ParametrizationVolPower::Dump


void ParametrizationVolPower::Visit(ParametrizationVisitor& visitor) const
{
  visitor.OnParametrizationVolPower(*this);
}//ParametrizationVolPower::Visit


} // namespace ihg

} // namespace ito33
