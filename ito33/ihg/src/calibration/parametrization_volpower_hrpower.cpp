/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_volpower.cpp
// Purpose:     implement for parametrization for power volatiltiy
// Author:      Ito33
// Created:     2004/12/14
// RCS-ID:      $Id: parametrization_volpower_hrpower.cpp,v 1.24 2006/08/20 09:31:04 wang Exp $
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
#include "ito33/finance/cdslike.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/basket_goodtype.h"

#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/version.h"
#include "ito33/ihg/error.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/ihg/parametrization_volpower_hrpower.h"
#include "ito33/ihg/hrspotcomponentpower.h"

#include "ihg/volpowerhrpowercalibrator.h"

#include "ito33/xml/write.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/parametrization_visitor.h"
#include "ihg/xml/common.h"

extern const ito33::finance::Error ITO33_CALIBRATION_FAIL;

extern const ito33::ihg::Error ITO33_IHG_CALIBRATION_DEGENERATE;

namespace ito33
{

namespace ihg
{

void ParametrizationVolPowerHRPower::CalibrateWithOptionsAndCDSs
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


  VolPowerHRPowerCalibrator calibrator;

  try
  {
    // The first option passed to the calibrator is used to obtain a
    // guess for volatility. The guess will be better if the option is
    // close to at-the-money
    double dStrike1 = option1.GetStrike();
    double dStrike2 = option2.GetStrike();
    double dS0 = option1.GetSessionData()->GetSpotSharePrice();
    if ( fabs(dStrike1-dS0) < fabs(dStrike2-dS0) )
      calibrator.Calibrate(option1, option2, tsCDS);
    else
      calibrator.Calibrate(option2, option1, tsCDS);
  }
  catch(const ito33::numeric::Exception&)
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

}//ParametrizationVolPowerHRPower::CalibrateWithOptionsAndCDSs


void ParametrizationVolPowerHRPower::Calibrate
     (const finance::BasketGoodType& basket)
{
  CalibrateWithOptionsAndCDSs(*basket.GetOption(),
                              *basket.GetOption2(),
                              *basket.GetTermStructureCDS());                             
}

void ParametrizationVolPowerHRPower::Dump(ito33::XML::Tag& tagParent) const
{
  // create the root tag
  ito33::XML::Tag 
    tagParametrization(XML_TAG_PARAMETRIZATION_VOLPOWERHRPOWER_ROOT,
                        tagParent);

  // nothing else to do!

} //ParametrizationVolPower::Dump


void ParametrizationVolPowerHRPower::Visit
  (ParametrizationVisitor& visitor) const
{
  visitor.OnParametrizationVolPowerHRPower(*this);
}//ParametrizationVolPower::Visit

shared_ptr<finance::TheoreticalModel> 
ParametrizationVolPowerHRPower::GetTheoreticalModel()
{
  shared_ptr<TheoreticalModel> 
    pTM( new TheoreticalModel(m_pVolatility, m_pHazardRate) );

  return pTM;
}

} // namespace ihg

} // namespace ito33
